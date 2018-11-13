/*****************************************************************************
 *
 *  gradient_3d_27pt_solid.h
 *
 *  $Id: gradient_3d_27pt_solid.h 2565 2014-12-14 14:06:16Z agray3 $
 *
 *  Edinburgh Soft Matter and Statistical Physics Group and
 *  Edinburgh Parallel Computing Centre
 *
 *  Kevin Stratford (kevin@epcc.ed.ac.uk)
 *  (c) 2010 The University of Edinburgh
 *
 *****************************************************************************/

#ifndef GRADIENT_3D_27PT_SOLID_H
#define GRADIENT_3D_27PT_SOLID_H

#include "map.h"

int gradient_3d_27pt_solid_map_set(map_t * map);
int gradient_3d_27pt_solid_d2(const int nop, const double * field,double * t_field,
			      double * grad,double * t_grad, double * delsq, double * t_delsq);
#endif
