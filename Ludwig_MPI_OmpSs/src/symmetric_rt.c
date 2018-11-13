/****************************************************************************
 *
 *  symmetric_rt.c
 *
 *  Run time initialisation for the symmetric phi^4 free energy.
 *
 *  $Id: symmetric_rt.c 2762 2015-10-21 12:48:43Z stratford $
 *
 *  Edinburgh Soft Matter and Statistical Physics Group
 *  and Edinburgh Parallel Computing Centre
 *
 *  Kevin Stratford (kevin@epcc.ed.ac.uk)
 *  (c) 2010 The University of Edinburgh
 *
 ****************************************************************************/

#include <assert.h>
#include <math.h>
#include <string.h>

#include "pe.h"
#include "coords.h"
#include "noise.h"
#include "runtime.h"
#include "free_energy.h"
#include "symmetric.h"

#include "physics.h"
#include "util.h"

#define  DEFAULT_NOISE      0.1
#define  DEFAULT_PATCH_SIZE 1
#define  DEFAULT_PATCH_VOL  0.5
#define  DEFAULT_RADIUS     8.0

static int symmetric_init_block(field_t * phi, double xi0);
static int symmetric_init_bath(field_t * phi);
int symmetric_init_spinodal(field_t * phi);
int symmetric_init_spinodal_patches(field_t * phi);
int symmetric_init_drop(field_t * fphi, double xi0, double radius);

/****************************************************************************
 *
 *  symmetric_run_time
 *
 ****************************************************************************/

void symmetric_run_time(void) {

  double a;
  double b;
  double kappa;

  info("Symmetric phi^4 free energy selected.\n");
  info("\n");

  /* Parameters */

  RUN_get_double_parameter("A", &a);
  RUN_get_double_parameter("B", &b);
  RUN_get_double_parameter("K", &kappa);

  info("Parameters:\n");
  info("Bulk parameter A      = %12.5e\n", a);
  info("Bulk parameter B      = %12.5e\n", b);
  info("Surface penalty kappa = %12.5e\n", kappa);

  symmetric_free_energy_parameters_set(a, b, kappa);

  info("Surface tension       = %12.5e\n", symmetric_interfacial_tension());
  info("Interfacial width     = %12.5e\n", symmetric_interfacial_width());

  /* Set free energy function pointers. */

  fe_density_set(symmetric_free_energy_density);
  fe_chemical_potential_set(symmetric_chemical_potential);
  fe_isotropic_pressure_set(symmetric_isotropic_pressure);
  fe_chemical_stress_set(symmetric_chemical_stress);

  return;
}

/*****************************************************************************
 *
 *  symmetric_rt_initial_conditions
 *
 *****************************************************************************/

int symmetric_rt_initial_conditions(field_t * phi) {

  int p;
  char value[BUFSIZ];
  char filestub[FILENAME_MAX];
  double radius;

  io_info_t * iohandler = NULL;

  assert(phi);

  p = RUN_get_string_parameter("phi_initialisation", value, BUFSIZ);

  /* Default is spinodal */
  if (p == 0 || strcmp(value, "spinodal") == 0) {
    info("Initialising phi for spinodal\n");
    symmetric_init_spinodal(phi);
  }

  if (p != 0 && strcmp(value, "patches") == 0) {
    info("Initialising phi in patches\n");
    symmetric_init_spinodal_patches(phi);
  }

  if (p != 0 && strcmp(value, "block") == 0) {
    info("Initialisng phi as block\n");
    symmetric_init_block(phi, symmetric_interfacial_width());
  }

  if (p != 0 && strcmp(value, "bath") == 0) {
    info("Initialising phi for bath\n");
    symmetric_init_bath(phi);
  }

  if (p != 0 && strcmp(value, "drop") == 0) {
    radius = DEFAULT_RADIUS;
    RUN_get_double_parameter("phi_init_drop_radius", &radius);
    info("Initialising droplet radius:     %14.7e\n", radius);
    symmetric_init_drop(phi, symmetric_interfacial_width(), radius);
  }

  if (p != 0 && strcmp(value, "from_file") == 0) {
    info("Initial order parameter requested from file\n");
    strcpy(filestub, "phi.init"); /* A default */
    RUN_get_string_parameter("phi_file_stub", filestub, FILENAME_MAX);
    info("Attempting to read phi from file: %s\n", filestub);

    field_io_info(phi, &iohandler);
    io_read_data(iohandler, filestub, phi);
  }

  return 0;
}

/*****************************************************************************
 *
 *  symmetric_init_drop
 *
 *****************************************************************************/

int symmetric_init_drop(field_t * fphi, double xi0, double radius) {

  int nlocal[3];
  int noffset[3];
  int index, ic, jc, kc;

  double position[3];
  double centre[3];
  double phi, r, rxi0;
  double phistar = 1.0;      /* "Amplitude", can be negative. */

  assert(fphi);

  coords_nlocal(nlocal);
  coords_nlocal_offset(noffset);

  rxi0 = 1.0/xi0;

  centre[X] = 0.5*L(X);
  centre[Y] = 0.5*L(Y);
  centre[Z] = 0.5*L(Z);

  RUN_get_double_parameter("phi_init_drop_amplitude", &phistar);
  info("Initialising droplet amplitude:  %14.7e\n", phistar);

  for (ic = 1; ic <= nlocal[X]; ic++) {
    for (jc = 1; jc <= nlocal[Y]; jc++) {
      for (kc = 1; kc <= nlocal[Z]; kc++) {

        index = coords_index(ic, jc, kc);
        position[X] = 1.0*(noffset[X] + ic) - centre[X];
        position[Y] = 1.0*(noffset[Y] + jc) - centre[Y];
        position[Z] = 1.0*(noffset[Z] + kc) - centre[Z];

        r = sqrt(dot_product(position, position));

        phi = phistar*tanh(rxi0*(r - radius));
	field_scalar_set(fphi, index, phi);
      }
    }
  }

  return 0;
}


/*****************************************************************************
 *
 * symmetric_init_block
 *
 *  Initialise two blocks with interfaces at z = Lz/4 and z = 3Lz/4.
 *
 *****************************************************************************/

static int symmetric_init_block(field_t * phi, double xi0) {

  int nlocal[3];
  int noffset[3];
  int ic, jc, kc, index;
  double z, z1, z2;
  double phi0;

  coords_nlocal(nlocal);
  coords_nlocal_offset(noffset);

  z1 = 0.25*L(Z);
  z2 = 0.75*L(Z);

  for (ic = 1; ic <= nlocal[X]; ic++) {
    for (jc = 1; jc <= nlocal[Y]; jc++) { 
      for (kc = 1; kc <= nlocal[Z]; kc++) {

	index = coords_index(ic, jc, kc);
	z = noffset[Z] + kc;

	if (z > 0.5*L(Z)) {
	  phi0 = tanh((z-z2)/xi0);
	}
	else {
	  phi0 = -tanh((z-z1)/xi0);
	}

	field_scalar_set(phi, index, phi0);
      }
    }
  }

  return 0;
}

/*****************************************************************************
 *
 *  symmetric_init_bath
 *
 *  Initialise one interface at z = Lz/8. This is inended for
 *  capillary rise in systems with z not periodic.
 *
 *****************************************************************************/

static int symmetric_init_bath(field_t * phi) {

  int nlocal[3];
  int noffset[3];
  int ic, jc, kc, index;
  double z, z0;
  double phi0, xi0;

  coords_nlocal(nlocal);
  coords_nlocal_offset(noffset);

  z0 = 0.25*L(Z);
  xi0 = 1.13;

  for (ic = 1; ic <= nlocal[X]; ic++) {
    for (jc = 1; jc <= nlocal[Y]; jc++) { 
      for (kc = 1; kc <= nlocal[Z]; kc++) {

	index = coords_index(ic, jc, kc);
	z = noffset[Z] + kc;
	phi0 = tanh((z-z0)/xi0);

	field_scalar_set(phi, index, phi0);
      }
    }
  }

  return 0;
}

/*****************************************************************************
 *
 *  symmetric_init_spinodal
 *
 *  Some random noise is required to initiate spinodal decomposition.
 *
 *****************************************************************************/

int symmetric_init_spinodal(field_t * phi) {

  int seed = 13;
  int ic, jc, kc, index;
  int nlocal[3];

  double ran;
  double phi0;                        /* Mean of phi */
  double phi1;                        /* Final phi value */
  double noise0 = DEFAULT_NOISE;      /* Amplitude of random noise */

  noise_t * rng = NULL;

  assert(phi);

  coords_nlocal(nlocal);
  physics_phi0(&phi0);

  RUN_get_int_parameter("random_seed", &seed);
  RUN_get_double_parameter("noise", &noise0);

  noise_create(&rng);
  noise_init(rng, seed);

  for (ic = 1; ic <= nlocal[X]; ic++) {
    for (jc = 1; jc <= nlocal[Y]; jc++) {
      for (kc = 1; kc <= nlocal[Z]; kc++) {

	index = coords_index(ic, jc, kc);

	noise_uniform_double_reap(rng, index, &ran);
	phi1 = phi0 + noise0*(ran - 0.5);
	field_scalar_set(phi, index, phi1);

      }
    }
  }

  noise_free(rng);

  return 0;
}

/*****************************************************************************
 *
 *  symmetric_init_spinodal_patches
 *
 *  Also for spinodal, but a slightly different strategy of patches,
 *  which is better for large composition ratios.
 *
 *  Generally, the further away from 50:50 one moves, the larger
 *  the patch size must be to prevent diffusion (via the mobility)
 *  washing out the spinodal decomposition.
 *
 *****************************************************************************/

int symmetric_init_spinodal_patches(field_t * phi) {

  int ic, jc, kc, index;
  int ip, jp, kp;
  int nlocal[3];
  int ipatch, jpatch, kpatch;
  int seed;
  int count = 0;
  int patch = DEFAULT_PATCH_SIZE;

  double volminus1 = DEFAULT_PATCH_VOL;
  double phi1;
  double ran_uniform;

  noise_t * rng = NULL;

  assert(phi);

  coords_nlocal(nlocal);

  RUN_get_int_parameter("random_seed", &seed);
  RUN_get_int_parameter("phi_init_patch_size", &patch);
  RUN_get_double_parameter("phi_init_patch_vol", &volminus1);

  noise_create(&rng);
  noise_init(rng, seed);

  for (ic = 1; ic <= nlocal[X]; ic += patch) {
    for (jc = 1; jc <= nlocal[Y]; jc += patch) {
      for (kc = 1; kc <= nlocal[Z]; kc += patch) {

	index = coords_index(ic, jc, kc);

	/* Uniform patch */
	phi1 = 1.0;
	noise_uniform_double_reap(rng, index, &ran_uniform);
	if (ran_uniform < volminus1) phi1 = -1.0;

	ipatch = dmin(nlocal[X], ic + patch - 1);
	jpatch = dmin(nlocal[Y], jc + patch - 1);
	kpatch = dmin(nlocal[Z], kc + patch - 1);

	for (ip = ic; ip <= ipatch; ip++) {
	  for (jp = jc; jp <= jpatch; jp++) {
	    for (kp = kc; kp <= kpatch; kp++) {

	      index = coords_index(ip, jp, kp);
	      field_scalar_set(phi, index, phi1);
	      count += 1;
	    }
	  }
	}

	/* Next patch */
      }
    }
  }

  noise_free(rng);

  assert(count == nlocal[X]*nlocal[Y]*nlocal[Z]);

  return 0;
}
