/*****************************************************************************
 *
 *  psi_sor.c
 *
 *  A solution of the Poisson equation for the potential and
 *  charge densities stored in the psi_t object.
 *
 *  The simple Poisson equation looks like
 *
 *    nabla^2 \psi = - rho_elec / epsilon
 *
 *  where psi is the potential, rho_elec is the free charge density, and
 *  epsilon is a permeability.
 *
 *  $Id: psi_sor.c 2762 2015-10-21 12:48:43Z stratford $
 *
 *  Edinburgh Soft Matter and Statistical Physics Group and
 *  Edinburgh Parallel Computing Centre
 *
 *  (c) 2013-2015 The University of Edinburgh
 *
 *  Contributing Authors:
 *    Kevin Stratford (kevin@epcc.ed.ac.uk)
 *    Ignacio Pagonabarraga (ipagonabarraga@ub.edu)
 *
 *****************************************************************************/

#include <assert.h>
#include <float.h>
#include <math.h>
#include <mpi.h>

#include "pe.h"
#include "coords.h"
#include "physics.h"
#include "psi_s.h"
#include "psi_sor.h"
#include "psi.h"
#include "control.h"
#include "model.h"
#include "util.h"

/*****************************************************************************
 *
 *  psi_sor_solve
 *
 *  If the f_vare_t argument is NULL, the uniform epsilon solver is used.
 *  If the argument is present, the non-uniform solver is used.
 *
 *****************************************************************************/

int psi_sor_solve(psi_t * obj, f_vare_t fepsilon) {

  assert(obj);

  if (fepsilon == NULL) psi_sor_poisson(obj);
  if (fepsilon != NULL) psi_sor_vare_poisson(obj, fepsilon);

  return 0;
}

/*****************************************************************************
 *
 *  psi_sor_poisson
 *
 *  Uniform permittivity. The differencing is a seven
 *  point stencil for \nabla^2 \psi. So
 *
 *  epsilon [ psi(i+1,j,k) - 2 psi(i,j,k) + psi(i-1,j,k)
 *          + psi(i,j+1,k) - 2 psi(i,j,k) + psi(i,j-1,k)
 *          + psi(i,j,k+1) - 2 psi(i,j,k) + psi(i,j,k-1) ] = -rho_elec(i,j,k)
 *
 *  We use the asymptotic estimate of the spectral radius for
 *  the Jacobi iteration
 *      radius ~= 1 - (pi^2 / 2N^2)
 *  where N is the linear dimension of the problem. It's important
 *  to get this right to keep the number of iterations as small as
 *  possible.
 *
 *  If this is an initial solve, the initial norm of the residual
 *  may be quite large (e.g., psi(t = 0)  = 0; rhs \neq 0); in this
 *  case a relative tolerance would be appropriate to decide
 *  termination. On subsequent calls, we might expect the initial
 *  residual to be relatively small (psi not charged much since
 *  previous solve), and an absolute tolerance might be appropriate.
 *
 *  The actual residual is checked against both at every 'ncheck'
 *  iterations, and either condition met will result in termination
 *  of the iteration. If neither criterion is met, the iteration will
 *  finish after 'niteration' iterations.
 *
 *  See, e.g., Press et al. Chapter 19.
 *
 *****************************************************************************/

int psi_sor_poisson(psi_t * obj) {

  int niteration = 1000;       /* Maximum number of iterations */
  const int ncheck = 5;        /* Check global residual every n iterations */
  
  int ic, jc, kc, index;
  int nhalo;
  int n;                       /* Relaxation iterations */
  int pass;                    /* Red/black iteration */
  int kst;                     /* Start kc index for red/black iteration */
  int nlocal[3];
  int xs, ys, zs;              /* Memory strides */
  double rho_elec;             /* Right-hand side */
  double residual;             /* Residual at given point */
  double rnorm[2];             /* Initial and current norm of residual */
  double rnorm_local[2];       /* Local values */
  double epsilon;              /* Uniform permittivity */
  double dpsi;
  double omega;                /* Over-relaxation parameter 1 < omega < 2 */
  double radius;               /* Spectral radius of Jacobi iteration */
  double tol_rel;              /* Relative tolerance */
  double tol_abs;              /* Absolute tolerance */
  double eunit, beta;

  /* int index_nbr, coords_nbr[3];*/

  MPI_Comm comm;               /* Cartesian communicator */

  nhalo = coords_nhalo();
  coords_nlocal(nlocal);
  comm = cart_comm();

  assert(nhalo >= 1);

  /* The red/black operation needs to be tested for odd numbers
   * of points in parallel. */
  assert(nlocal[X] % 2 == 0);
  assert(nlocal[Y] % 2 == 0);
  assert(nlocal[Z] % 2 == 0);

  zs = 1;
  ys = zs*(nlocal[Z] + 2*nhalo);
  xs = ys*(nlocal[Y] + 2*nhalo);

  /* Compute initial norm of the residual */

  radius = 1.0 - 0.5*pow(4.0*atan(1.0)/imax(L(X),L(Z)), 2);

  psi_epsilon(obj, &epsilon);
  psi_reltol(obj, &tol_rel);
  psi_abstol(obj, &tol_abs);
  psi_maxits(obj, &niteration);

  psi_beta(obj, &beta);
  psi_unit_charge(obj, &eunit);

  rnorm_local[0] = 0.0;

  for (ic = 1; ic <= nlocal[X]; ic++) {
    for (jc = 1; jc <= nlocal[Y]; jc++) {
      for (kc = 1; kc <= nlocal[Z]; kc++) {

	index = coords_index(ic, jc, kc);

	psi_rho_elec(obj, index, &rho_elec);

        /* 6-point stencil of Laplacian */

	dpsi = obj->psi[index + xs] + obj->psi[index - xs]
	     + obj->psi[index + ys] + obj->psi[index - ys]
	     + obj->psi[index + zs] + obj->psi[index - zs]
	     - 6.0*obj->psi[index];

	/* Non-dimensional potential in Poisson eqn requires e/kT */
	rnorm_local[0] += fabs(epsilon*dpsi + eunit*beta*rho_elec);
      }
    }
  }

  /* Iterate to solution */

  omega = 1.0;

  for (n = 0; n < niteration; n++) {

    /* Compute current normal of the residual */

    rnorm_local[1] = 0.0;

    for (pass = 0; pass < 2; pass++) {

      for (ic = 1; ic <= nlocal[X]; ic++) {
	for (jc = 1; jc <= nlocal[Y]; jc++) {
	  kst = 1 + (ic + jc + pass) % 2;
	  for (kc = kst; kc <= nlocal[Z]; kc += 2) {

	    index = coords_index(ic, jc, kc);

	    psi_rho_elec(obj, index, &rho_elec);

	    /* 6-point stencil of Laplacian */

	    dpsi = obj->psi[index + xs] + obj->psi[index - xs]
	         + obj->psi[index + ys] + obj->psi[index - ys]
	         + obj->psi[index + zs] + obj->psi[index - zs]
	      - 6.0*obj->psi[index];

	    /* Non-dimensional potential in Poisson eqn requires e/kT */
	    residual = epsilon*dpsi + eunit*beta*rho_elec;
	    obj->psi[index] -= omega*residual / (-6.0*epsilon);
	    rnorm_local[1] += fabs(residual);
	  }
	}
      }

      /* Recompute relaxation parameter and next pass */

      if (n == 0 && pass == 0) {
	omega = 1.0 / (1.0 - 0.5*radius*radius);
      }
      else {
	omega = 1.0 / (1.0 - 0.25*radius*radius*omega);
      }
      assert(1.0 < omega);
      assert(omega < 2.0);

      psi_halo_psi(obj);
      psi_halo_psijump(obj);

    }

    if ((n % ncheck) == 0) {
      /* Compare residual and exit if small enough */

      MPI_Allreduce(rnorm_local, rnorm, 2, MPI_DOUBLE, MPI_SUM, comm);

      if (rnorm[1] < tol_abs) {

	if (is_psi_resid_step()) {
	  info("\n");
	  info("SOR solver converged to absolute tolerance\n");
	  info("SOR residual per site %14.7e at %d iterations\n",
	       rnorm[1]/(L(X)*L(Y)*L(Z)), n);
	}
	break;
      }

      if (rnorm[1] < tol_abs || rnorm[1] < tol_rel*rnorm[0]) {

	if (is_psi_resid_step()) {
	  info("\n");
	  info("SOR solver converged to relative tolerance\n");
	  info("SOR residual per site %14.7e at %d iterations\n",
	       rnorm[1]/(L(X)*L(Y)*L(Z)), n);
	}
	break;
      }
    }
 
    if (n == niteration-1) {
      info("\n");
      info("SOR solver exceeded %d iterations\n", n+1);
      info("SOR residual %le (initial) %le (final)\n\n", rnorm[0], rnorm[1]);
    }
  }

  return 0;
}

/*****************************************************************************
 *
 *  psi_sor_vare_poisson
 *
 *  This is essentially a copy of the above, but it allows a spatially
 *  varying permittivity epsilon:
 *
 *    div [epsilon(r) grad phi(r) ] = -rho(r)
 * 
 ****************************************************************************/

int psi_sor_vare_poisson(psi_t * obj, f_vare_t fepsilon) {

  int niteration = 2000;       /* Maximum number of iterations */
  const int ncheck = 1;        /* Check global residual every n iterations */
  
  int ic, jc, kc, index;
  int n;                       /* Relaxation iterations */
  int pass;                    /* Red/black iteration */
  int kst;                     /* Start kc index for red/black iteration */
  int nlocal[3];
  int xs, ys, zs;              /* Memory strides */

  double rho_elec;             /* Right-hand side */
  double residual;             /* Residual at given point */
  double rnorm[2];             /* Initial and current norm of residual */
  double rnorm_local[2];       /* Local values */

  double depsi;                /* Differenced left-hand side */
  double eps0, eps1;           /* Permittivity values */

  double omega;                /* Over-relaxation parameter 1 < omega < 2 */
  double radius;               /* Spectral radius of Jacobi iteration */

  double tol_rel;              /* Relative tolerance */
  double tol_abs;              /* Absolute tolerance */

  double eunit, beta;

  MPI_Comm comm;               /* Cartesian communicator */

  coords_nlocal(nlocal);
  comm = cart_comm();

  assert(coords_nhalo() >= 1);

  /* The red/black operation needs to be tested for odd numbers
   * of points in parallel. */

  assert(nlocal[X] % 2 == 0);
  assert(nlocal[Y] % 2 == 0);
  assert(nlocal[Z] % 2 == 0);

  coords_strides(&xs, &ys, &zs);

  /* Compute initial norm of the residual */

  radius = 1.0 - 0.5*pow(4.0*atan(1.0)/imax(L(X),L(Z)), 2);

  psi_reltol(obj, &tol_rel);
  psi_abstol(obj, &tol_abs);
  psi_maxits(obj, &niteration);
  psi_beta(obj, &beta);
  psi_unit_charge(obj, &eunit);

  rnorm_local[0] = 0.0;

  for (ic = 1; ic <= nlocal[X]; ic++) {
    for (jc = 1; jc <= nlocal[Y]; jc++) {
      for (kc = 1; kc <= nlocal[Z]; kc++) {

	depsi = 0.0;

	index = coords_index(ic, jc, kc);

	psi_rho_elec(obj, index, &rho_elec);
	fepsilon(index, &eps0);

	/* Laplacian part of operator */

        depsi += eps0*(-6.0*obj->psi[index]
	      + obj->psi[index + xs] + obj->psi[index - xs]
	      + obj->psi[index + ys] + obj->psi[index - ys]
	      + obj->psi[index + zs] + obj->psi[index - zs]);

	/* Additional terms in generalised Poisson equation */

	fepsilon(index + xs, &eps1);
	depsi += 0.25*eps1*(obj->psi[index + xs] - obj->psi[index - xs]);

	fepsilon(index - xs, &eps1);
	depsi -= 0.25*eps1*(obj->psi[index + xs] - obj->psi[index - xs]);

	fepsilon(index + ys, &eps1);
	depsi += 0.25*eps1*(obj->psi[index + ys] - obj->psi[index - ys]);

	fepsilon(index - ys, &eps1);
	depsi -= 0.25*eps1*(obj->psi[index + ys] - obj->psi[index - ys]);

	fepsilon(index + zs, &eps1);
	depsi += 0.25*eps1*(obj->psi[index + zs] - obj->psi[index - zs]);

	fepsilon(index - zs, &eps1);
	depsi -= 0.25*eps1*(obj->psi[index + zs] - obj->psi[index - zs]);

	/* Non-dimensional potential in Poisson eqn requires e/kT */
	rnorm_local[0] += fabs(depsi + eunit*beta*rho_elec);
      }
    }
  }

  /* Iterate to solution */

  omega = 1.0;

  for (n = 0; n < niteration; n++) {

    /* Compute current normal of the residual */

    rnorm_local[1] = 0.0;

    for (pass = 0; pass < 2; pass++) {

      for (ic = 1; ic <= nlocal[X]; ic++) {
	for (jc = 1; jc <= nlocal[Y]; jc++) {
	  kst = 1 + (ic + jc + pass) % 2;
	  for (kc = kst; kc <= nlocal[Z]; kc += 2) {

	    depsi  = 0.0;

	    index = coords_index(ic, jc, kc);

	    psi_rho_elec(obj, index, &rho_elec);
	    fepsilon(index, &eps0);

	    /* Laplacian part of operator */

	    depsi += eps0*(-6.0*obj->psi[index]
		  + obj->psi[index + xs] + obj->psi[index - xs]
		  + obj->psi[index + ys] + obj->psi[index - ys]
		  + obj->psi[index + zs] + obj->psi[index - zs]);

	    /* Additional terms in generalised Poisson equation */

	    fepsilon(index + xs, &eps1);
	    depsi += 0.25*eps1*(obj->psi[index + xs] - obj->psi[index - xs]);

	    fepsilon(index - xs, &eps1);
	    depsi -= 0.25*eps1*(obj->psi[index + xs] - obj->psi[index - xs]);

	    fepsilon(index + ys, &eps1);
	    depsi += 0.25*eps1*(obj->psi[index + ys] - obj->psi[index - ys]);

	    fepsilon(index - ys, &eps1);
	    depsi -= 0.25*eps1*(obj->psi[index + ys] - obj->psi[index - ys]);

	    fepsilon(index + zs, &eps1);
	    depsi += 0.25*eps1*(obj->psi[index + zs] - obj->psi[index - zs]);

	    fepsilon(index - zs, &eps1);
	    depsi -= 0.25*eps1*(obj->psi[index + zs] - obj->psi[index - zs]);

	    /* Non-dimensional potential in Poisson eqn requires e/kT */
	    residual = depsi + eunit*beta*rho_elec;
	    obj->psi[index] -= omega*residual / (-6.0*eps0);
	    rnorm_local[1] += fabs(residual);

	  }
	}
      }

      psi_halo_psi(obj);
      psi_halo_psijump(obj);

    }

    /* Recompute relation parameter */
    /* Note: The default Chebychev acceleration causes a convergence problem */ 
    omega = 1.0 / (1.0 - 0.25*radius*radius*omega);

    if ((n % ncheck) == 0) {

      /* Compare residual and exit if small enough */

      MPI_Allreduce(rnorm_local, rnorm, 2, MPI_DOUBLE, MPI_SUM, comm);

      if (rnorm[1] < tol_abs) {

	if (is_psi_resid_step()) {
	  info("\n");
	  info("SOR (heterogeneous) solver converged to absolute tolerance\n");
	  info("SOR residual per site %14.7e at %d iterations\n",
	       rnorm[1]/(L(X)*L(Y)*L(Z)), n);
	}
	break;
      }

      if (rnorm[1] < tol_rel*rnorm[0]) {

	if (is_psi_resid_step()) {
	  info("\n");
	  info("SOR (heterogeneous) solver converged to relative tolerance\n");
	  info("SOR residual per site %14.7e at %d iterations\n",
	       rnorm[1]/(L(X)*L(Y)*L(Z)), n);
	}
	break;
      }

      if (n == niteration-1) {
	info("\n");
	info("SOR solver (heterogeneous) exceeded %d iterations\n", n+1);
	info("SOR residual %le (initial) %le (final)\n\n", rnorm[0], rnorm[1]);
      }
    }
  }

  return 0;
}
