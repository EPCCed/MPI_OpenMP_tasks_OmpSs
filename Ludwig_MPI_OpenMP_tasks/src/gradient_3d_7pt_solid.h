/*****************************************************************************
 *
 *  gradient_3d_7pt_solid.h
 *
 *  $Id: gradient_3d_7pt_solid.h 2645 2015-05-14 13:20:32Z agray3 $
 *
 *  Edinburgh Soft Matter and Statistical Physics Group and
 *  Edinburgh Parallel Computing Centre
 *
 *  Kevin Stratford (kevin@epcc.ed.ac.uk)
 *  (c) 2011 The University of Edinburgh
 *
 *****************************************************************************/

#ifndef GRADIENT_3D_7PT_SOLID_H
#define GRADIENT_3D_7PT_SOLID_H

#include "map.h"

__targetHost__ int gradient_3d_7pt_solid_map_set(map_t * map);
__targetHost__ int gradient_3d_7pt_solid_d2(const int nop, const double * field,double * t_field,
			     double * grad,double * t_grad, double * delsq, double * t_delsq);
#endif
