/*****************************************************************************
 *
 *  phi_force_colloid.h
 *
 *  $Id: phi_force_colloid.h,v 1.2 2010-10-15 12:40:03 kevin Exp $
 *
 *  Edinburgh Soft Matter and Statistical Physics Group and
 *  Edinburgh Parallel Computing Centre
 *
 *  (c) The University of Edinburgh (2009)
 *  Kevin Stratford (kevin@epcc.ed.ac.uk)
 *
 *****************************************************************************/

#ifndef PHI_FORCE_COLLOID_H
#define PHI_FORCE_COLLOID_H

#include "colloids.h"
#include "hydro.h"
#include "map.h"

__targetHost__ int phi_force_colloid(colloids_info_t * cinfo, field_t* q, field_grad_t* q_grad, hydro_t * hydro, map_t * map);

#endif
