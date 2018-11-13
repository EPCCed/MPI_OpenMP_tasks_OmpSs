/*****************************************************************************
 *
 *  stats_surfactant.c
 *
 *  Some routines to perform analysis of the surfactant model.
 *
 *  $Id: stats_surfactant.c,v 1.2 2010-10-15 12:40:03 kevin Exp $
 *
 *  Edinburgh Soft Matter and Statistical Physics Group and
 *  Edinburgh Parallel Computing Centre
 *
 *  Kevin Stratford (kevin@epcc.ed.ac.uk)
 *  (c) 2009 The University of Edinburgh
 *
 *****************************************************************************/

#include <assert.h>
#include <math.h>

#include "pe.h"
#include "coords.h"
#include "field.h"
#include "control.h"
#include "util.h"
#include "surfactant.h"
#include "stats_surfactant.h"

/*****************************************************************************
 *
 *  stats_surfactant_1d
 *
 *  For a 1-d model in which we have a block initialisation of the
 *  composition, work out the profile of the free energy, and hence
 *  the current interfacial tension.
 *
 *****************************************************************************/

int stats_surfactant_1d(field_t * fphi) {

  int index;
  int ic = 1, jc = 1, kc;
  int nlocal[3];
  double e, e0;
  double psi_0, psi_b;
  double sigma, sigma0;
  double phi[2];

  /* This is not run in parallel, so assert it's serial.
   * We also require surfactant */

  assert(fphi);
  assert(0); /* Check nf = 2 in refactored version */
  assert(pe_size() == 1);

  coords_nlocal(nlocal);

  /* We assume z = 1 is a reasonable choice for the background
   * free energy level, which we need to subtract to find the
   * excess. */

  kc = 1;
  index = coords_index(ic, jc, kc);
  e0 = surfactant_free_energy_density(index);
  field_scalar_array(fphi, index, phi);

  psi_b = phi[1];

  /* To compute the surface tension, run through both interfaces
   * and divide the final excess free energy by 2. We also record
   * the maximum surfactant value (psi_0) along the way. */ 

  sigma = 0.0;
  psi_0 = 0.0;

  for (kc = 1; kc <= nlocal[Z]; kc++) {

    index = coords_index(ic, jc, kc);

    e = surfactant_free_energy_density(index);
    sigma += 0.5*(e - e0);
    field_scalar_array(fphi, index, phi);
    psi_0 = dmax(psi_0, phi[1]);
  }

  /* Compute the fractional reduction in the surface tension
   * below the bare surface value */

  sigma0 = surfactant_interfacial_tension();
  sigma = (sigma - sigma0)/sigma0;

  /* The sqrt(t) is the usual dependance for analysis of the
   * diffusion problem, so is included here. */

  info("Surfactant: %d %12.5e %12.5e %12.5e %12.5e\n", get_step(),
       sqrt(1.0*get_step()), psi_b, psi_0, sigma);

  return 0;
}
