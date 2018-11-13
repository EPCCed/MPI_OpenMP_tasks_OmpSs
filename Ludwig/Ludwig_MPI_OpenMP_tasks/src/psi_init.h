/*****************************************************************************
*
*  psi_init.h
*
*  Various initial states for electrokinetics.
*
*  $Id: psi_init.h 2479 2014-09-25 11:24:42Z stratford $
*
*  Edinburgh Soft Matter and Statistical Physics Group and
*  Edinburgh Parallel Computing Centre
*
*  Oliver Henrich (o.henrich@ucl.ac.uk) wrote most of these.
*  (c) 2012 The University of Edinburgh
*
*****************************************************************************/

#ifndef PSI_INIT_H
#define PSI_INIT_H

#include "psi.h"
#include "map.h"

int psi_init_uniform(psi_t * obj, double rho_el);
int psi_init_gouy_chapman_set(psi_t * obj, map_t * map, double rho_el,
			      double sigma);
int psi_init_liquid_junction_set(psi_t * obj, double rho_el, double delta_el);

#endif 
