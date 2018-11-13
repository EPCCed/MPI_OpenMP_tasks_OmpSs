/*****************************************************************************
 *
 *  brownian.c
 *
 *  $Id: brownian.c,v 1.4 2010-10-15 12:40:02 kevin Exp $
 *
 *  Edinburgh Soft Matter and Statistical Physics Group and
 *  Edinburgh Parallel Computing Centre
 *
 *  Kevin Stratford (kevin@epcc.ed.ac.uk)
 *  (c) 2007 The University of Edinburgh
 *
 *****************************************************************************/

#include <assert.h>
#include <math.h>

#include "pe.h"
#include "ran.h"
#include "coords.h"
#include "colloids.h"
#include "colloids_halo.h"
#include "physics.h"
#include "util.h"

void brownian_step_ermak_buckholz(void);
void brownian_step_no_inertia(void);
void brownian_step_no_inertia_test(void);
void brownian_set_random(void);

static double dt_ = 1.0;  /* Time step; fixed here to be 1 LB step */
#define NTIME 16

/*****************************************************************************
 *
 *  do_brownian_dynamics
 *
 *  Perform a number of Brownian Dynamics steps.
 *
 *****************************************************************************/

void do_brownian_dynamics() {

  int     ia;
  int     ic, jc, kc;
  colloid_t * p_colloid;
  int i, n, nt, nmax = 100000;
  double u[3][NTIME];
  double r[3][NTIME];
  double s[3][NTIME];
  double correlator_uu[NTIME];
  double correlator_ur[NTIME];
  double correlator_rr[NTIME];
  double correlator_ss[NTIME];
  double mass, scaleu, scalex, beta;
  double dr[3];

  /* Brownian dynamics requires an update. */
  assert(0);

  brownian_step_no_inertia_test();
  if (1) return;

  mass  = (4.0/3.0)*pi_*pow(2.3, 3);
  beta  = 6.0*pi_*get_eta_shear()*2.3/mass;
  scaleu = sqrt(mass/(3.0*get_kT()));
  scalex = beta*scaleu;

  dt_ = 0.1/beta;
  dt_ = 1.0;

  info("beta = %f\n", beta);
  info("scale u is %f\n", scaleu);
  info("scale x is %f\n", scalex);

  for (nt = 0; nt < NTIME; nt++) {
    for (i = 0; i < 3; i++) {
      u[i][nt] = 0.0;
      r[i][nt] = 0.0;
      s[i][nt] = 0.0;
    }
    correlator_uu[nt] = 0.0;
    correlator_ur[nt] = 0.0;
    correlator_rr[nt] = 0.0;
    correlator_ss[nt] = 0.0;
  }

  for (n = 1; n <= nmax; n++) {

  /* Set random numbers for each particle */
    brownian_set_random();
    colloids_halo_state();
    /* brownian_step_ermak_buckholz();*/
    brownian_step_no_inertia();

    /* diagnostics */

    for (nt = NTIME-1; nt >= 1; nt--) {
      for (i = 0; i < 3; i++) {
	u[i][nt] = u[i][nt-1];
	r[i][nt] = r[i][nt-1];
	s[i][nt] = s[i][nt-1];
      }
    }

    for (ic = 1; ic <= Ncell(X); ic++) {
      for (jc = 1; jc <= Ncell(Y); jc++) {
	for (kc = 1; kc <= Ncell(Z); kc++) {

	  p_colloid = colloids_cell_list(ic, jc, kc);

	  while (p_colloid) {
	    for (ia = 0; ia < 3; ia++) {
	      r[ia][0] = p_colloid->s.r[ia];
	      u[ia][0] = p_colloid->s.v[ia];
	      s[ia][0] = p_colloid->s.s[ia];
	    }

	    p_colloid = p_colloid->next;
	  }
	}
      }
    }

    if ( n > NTIME) {
      for (nt = 0; nt < NTIME; nt++) {
	for (i = 0; i < 3; i++) {
	  correlator_uu[nt] += scaleu*scaleu*u[i][nt]*u[i][0];
	  /* -ve sign here because r[0] is current step and nt previous */
	  dr[i] = -(r[i][nt] - r[i][0]);
	  if (dr[i] > L(i)/2.0) dr[i] -= L(i);
	  if (dr[i] <-L(i)/2.0) dr[i] += L(i);
	  correlator_ur[nt] += scaleu*u[i][nt]*scalex*dr[i];
	  correlator_rr[nt] += dr[i]*dr[i];
	  correlator_ss[nt] += (s[i][nt] - s[i][0])*(s[i][nt] - s[i][0]);
	}
      }
    }

    colloids_cell_update();
  }

  info("\nResults Ermak and Buckholz \n");
  for (n = 0; n < NTIME; n++) {
    double tau = n*beta*dt_;
    info("%6.3f %6.3f %6.3f ", tau, exp(-tau), correlator_uu[n]/nmax);
    info("%6.3f %6.3f ", 1.0 -  exp(-tau), sqrt(3.0)*correlator_ur[n]/nmax);
    info("%6.3f %6.3f\n", 2.0*(tau - 1.0 + exp(-tau)),
	 scalex*scalex*correlator_rr[n]/nmax);
  }

  /* No inertia test */
  for (n = 0; n < NTIME; n++) {
    double tau = n*dt_;
    double difft = get_kT() / (6.0*pi_*get_eta_shear()*2.3);
    double diffr = get_kT() / (8.0*pi_*get_eta_shear()*2.3*2.3*2.3);
    info("%8.3f %8.3f %8.3f ",  tau, 2.0*tau, correlator_rr[n]/(3.0*difft*nmax));
    info("%8.3f %8.3f\n", 4.0*tau, correlator_ss[n]/(diffr*nmax));
  }

  return;
}

/*****************************************************************************
 *
 *  brownian_step_no_inertia
 *
 *  Brownian dynamics with no velocity; only position is updated.
 *
 *  Rotational and translation parts; see, for example, Meriget etal
 *  J. Chem Phys., 121 6078 (2004).
 *  
 *  r(t + dt) = r(t) +  rgamma_t dt Force + r_random
 *  s(t + dt) = s(t) + (rgamma_r dt Torque + t_random) x s(t)
 *
 *  The translational and rotational friction coefficients are
 *  gamma_t  = 6 pi eta a
 *  gamma_r  = 8.0 pi eta a^3
 *
 *  The variances of the random translational and rotational
 *  contributions are related to the diffusion constants (kT/gamma)
 *  <r_i . r_j> = 2 dt (kT / gamma_t) delta_ij
 *  <t_i . t_j> = 2 dt (kT / gamma_r) delta_ij
 * 
 *  This requires six random variates per time step.
 *
 *****************************************************************************/

void brownian_step_no_inertia() {

  int       ic, jc, kc;
  colloid_t * p_colloid;
  double    ran[3];
  double    eta, kT;
  double    sigma, rgamma;

  kT = get_kT();
  eta = get_eta_shear();

  /* Update each particle */

  for (ic = 0; ic <= Ncell(X) + 1; ic++) {
    for (jc = 0; jc <= Ncell(Y) + 1; jc++) {
      for (kc = 0; kc <= Ncell(Z) + 1; kc++) {

	p_colloid = colloids_cell_list(ic, jc, kc);

	while (p_colloid) {

	  /* Translational motion */

	  rgamma = 1.0 / (6.0*pi_*eta*p_colloid->s.ah);
	  sigma = sqrt(2.0*dt_*kT*rgamma);

	  ran[X] = sigma*p_colloid->random[0];
	  ran[Y] = sigma*p_colloid->random[1];
	  ran[Z] = sigma*p_colloid->random[2];

	  p_colloid->s.r[X] += dt_*rgamma*p_colloid->force[X] + ran[X];
	  p_colloid->s.r[Y] += dt_*rgamma*p_colloid->force[Y] + ran[Y];
	  p_colloid->s.r[Z] += dt_*rgamma*p_colloid->force[Z] + ran[Z];

	  /* Rotational motion */

	  rgamma = 1.0 / (8.0*pi_*eta*pow(p_colloid->s.ah, 3));
	  sigma = sqrt(2.0*dt_*kT*rgamma);

	  ran[X] = dt_*rgamma*p_colloid->torque[X]
	    + sigma*p_colloid->random[3+X];
	  ran[Y] = dt_*rgamma*p_colloid->torque[Y]
	    + sigma*p_colloid->random[3+Y];
	  ran[Z] = dt_*rgamma*p_colloid->torque[Z]
	    + sigma*p_colloid->random[3+Z];

	  rotate_vector(p_colloid->s.s, ran);

	  /* Next colloid */

	  p_colloid = p_colloid->next;
	}

	/* Next cell */
      }
    }
  }

  return;
}

/*****************************************************************************
 *
 *  brownian_step_ermak_buckholz
 *
 *  This is the method of Ermak and Buckholz J. Comp. Phys 35, 169 (1980).
 *  It is valid when the time step is of the same order as the decay time
 *  for the velocity autocorrelation function.
 *
 *  See also Allen and Tildesley.
 *
 *  Only the translational part is present at the moment. This requires
 *  6 random numbers per particle per timestep. The updates for the
 *  position and momentum must be correlated (see e.g., Allen & Tildesley
 *  Appendix G.3)   
 *
 *****************************************************************************/

void brownian_step_ermak_buckholz() {

  int     ic, jc, kc;
  colloid_t * pc;

  double ranx, rany, ranz;
  double cx, cy, cz;

  double c0, c1, c2;
  double xi, xidt;
  double sigma_r, sigma_v;
  double c12, c21;
  double rmass, kT;

  kT = get_kT();

  /* Update each particle */

  for (ic = 0; ic <= Ncell(X) + 1; ic++) {
    for (jc = 0; jc <= Ncell(Y) + 1; jc++) {
      for (kc = 0; kc <= Ncell(Z) + 1; kc++) {

	pc = colloids_cell_list(ic, jc, kc);

	while (pc) {

	  rmass = 1.0/((4.0/3.0)*pi_*pow(pc->s.ah, 3.0));

	  /* Friction coefficient is xi, and related quantities */

	  xi = 6.0*pi_*get_eta_shear()*pc->s.ah*rmass;
	  xidt = xi*dt_;

	  c0 = exp(-xidt);
	  c1 = (1.0 - c0) / xi;
	  c2 = (dt_ - c1) / xi;

	  /* Velocity update */

	  sigma_v = sqrt(rmass*kT*(1.0 - c0*c0));

	  ranx = pc->random[0];
	  rany = pc->random[1];
	  ranz = pc->random[2];

	  cx = c0*pc->s.v[X] + rmass*c1*pc->force[X] + sigma_v*ranx;
	  cy = c0*pc->s.v[Y] + rmass*c1*pc->force[Y] + sigma_v*rany;
	  cz = c0*pc->s.v[Z] + rmass*c1*pc->force[Z] + sigma_v*ranz;

	  pc->s.v[X] = cx;
	  pc->s.v[Y] = cy;
	  pc->s.v[Z] = cz;

	  /* Generate correlated random pairs */

	  sigma_r = sqrt(rmass*kT*(2.0*xidt - 3.0 + 4.0*c0 - c0*c0)/(xi*xi));
	  c12 = rmass*kT*(1.0 - c0)*(1.0 - c0) / (sigma_v*sigma_r*xi);
	  c21 = sqrt(1.0 - c12*c12);

	  ranx = (c12*ranx + c21*pc->random[3]);
	  rany = (c12*rany + c21*pc->random[4]);
	  ranz = (c12*ranz + c21*pc->random[5]);

	  /* Position update */

	  cx = c1*pc->s.v[X] + rmass*c2*pc->force[X] + sigma_r*ranx;
	  cy = c1*pc->s.v[Y] + rmass*c2*pc->force[Y] + sigma_r*rany;
	  cz = c1*pc->s.v[Z] + rmass*c2*pc->force[Z] + sigma_r*ranz;

	  pc->s.r[X] += cx;
	  pc->s.r[Y] += cy;
	  pc->s.r[Z] += cz;

	  pc = pc->next;
	}
      }
    }
  }

  return;
}

/*****************************************************************************
 *
 *  brownian_set_random
 *
 *  Set the table of random numbers for each particle.
 *
 *****************************************************************************/

void brownian_set_random() {

  int     ic, jc, kc;
  colloid_t * p_colloid;

  /* Set random numbers for each particle */

  for (ic = 1; ic <= Ncell(X); ic++) {
    for (jc = 1; jc <= Ncell(Y); jc++) {
      for (kc = 1; kc <= Ncell(Z); kc++) {

	p_colloid = colloids_cell_list(ic, jc, kc);

	while (p_colloid) {

	  p_colloid->random[0] = ran_parallel_gaussian();
	  p_colloid->random[1] = ran_parallel_gaussian();
	  p_colloid->random[2] = ran_parallel_gaussian();
	  p_colloid->random[3] = ran_parallel_gaussian();
	  p_colloid->random[4] = ran_parallel_gaussian();
	  p_colloid->random[5] = ran_parallel_gaussian();

	  p_colloid = p_colloid->next;
	}
      }
    }
  }

  return;
}

/*****************************************************************************
 *
 *  brownian_step_no_interia_test
 *
 *  We test the single-particle position and orientational
 *  correlators:
 *
 *  < [r(t) - r(0)][r(t) - r(0)] > = 2D_t t I     t >> m / gamma_t
 *  < |s(t) - s(0)|^2 >            = 4D_r t       2D_r t << 1
 *
 *  where I is the identity matrix.
 *
 *  See, for example, Jan Dhont 'An Introduction to dynamics of colloids'
 *  Chapter 2.
 *
 *  The correlators should be good to about 1% after 100,000 iterations.
 *
 *****************************************************************************/

void brownian_step_no_inertia_test() {

  const int nmax  = 10000;

  int     ia;
  int     ic, jc, kc;
  colloid_t * p_colloid;
  int     i, n, nt;
  int     ncell[3];
  double  diffr, difft;
  double  r[3][NTIME];
  double  s[3][NTIME];
  double  correlator_rr[3][NTIME];
  double  correlator_ss[NTIME];
  double  dr[3];


  difft = get_kT() / (6.0*pi_*get_eta_shear()*2.3);
  diffr = get_kT() / (8.0*pi_*get_eta_shear()*2.3*2.3*2.3);

  /* Initialise */

  for (i = 0; i < 3; i++) ncell[i] = Ncell(i);

  for (nt = 0; nt < NTIME; nt++) {
    for (i = 0; i < 3; i++) {
      correlator_rr[i][nt] = 0.0;
    }
    correlator_ss[nt] = 0.0;
  }

  /* Start the main iteration loop */

  for (n = 1; n <= nmax; n++) {

    brownian_set_random();

    colloids_halo_state();
    brownian_step_no_inertia();
    colloids_cell_update();

    /* Rotate the tables of position and orientation backwards
     * and store the current values */

    for (nt = NTIME-1; nt >= 1; nt--) {
      for (i = 0; i < 3; i++) {
	r[i][nt] = r[i][nt-1];
	s[i][nt] = s[i][nt-1];
      }
    }

    for (ic = 1; ic <= ncell[X]; ic++) {
      for (jc = 1; jc <= ncell[Y]; jc++) {
	for (kc = 1; kc <= ncell[Z]; kc++) {

	  p_colloid = colloids_cell_list(ic, jc, kc);

	  while (p_colloid) {

	    for (ia = 0; ia < 3; ia++) {
	      r[ia][0] = p_colloid->s.r[ia];
	      s[ia][0] = p_colloid->s.s[ia];
	    }

	    p_colloid = p_colloid->next;
	  }
	}
      }
    }

    /* Work out the contribution to the correlators */

    if ( n > NTIME) {
      for (nt = 0; nt < NTIME; nt++) {
	for (i = 0; i < 3; i++) {
	  /* correct for periodic boundaries */
	  dr[i] = (r[i][nt] - r[i][0]);
	  if (dr[i] > L(i)/2.0) dr[i] -= L(i);
	  if (dr[i] <-L(i)/2.0) dr[i] += L(i);
	  correlator_rr[i][nt] += dr[i]*dr[i];
	  correlator_ss[nt] += (s[i][nt] - s[i][0])*(s[i][nt] - s[i][0]);
	}
      }
    }

  }

  /* Results */

  info("\n\n");
  info("           Translational                       Orientational\n");
  info("    Time   theory     <rx>     <ry>     <rz>   ");
  info("theory     <ss>\n");
  
  for (n = 0; n < NTIME; n++) {
    double tau = n*dt_;
    info("%8.3f %8.3f ",  tau, 2.0*tau);
    info("%8.3f ", correlator_rr[X][n]/(difft*nmax));
    info("%8.3f ", correlator_rr[Y][n]/(difft*nmax));
    info("%8.3f ", correlator_rr[Z][n]/(difft*nmax));
    info("%8.3f %8.3f\n", 4.0*tau, correlator_ss[n]/(diffr*nmax));
  }

  return;
}
