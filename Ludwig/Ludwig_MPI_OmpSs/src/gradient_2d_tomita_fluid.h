/*****************************************************************************
 *
 *  gradient_2d_tomita_fluid.h
 *
 *  $Id: gradient_2d_tomita_fluid.h 2565 2014-12-14 14:06:16Z agray3 $
 *
 *  Edinburgh Soft Matter and Statistical Physics Group and
 *  Edinburgh Parallel Computing Centre
 *
 *  Kevin Stratford (kevin@epcc.ed.ac.uk)
 *  (c) 2011 The University of Edinburgh
 *
 *****************************************************************************/

#ifndef GRADIENT_2D_TOMITA_FLUID_H
#define GRADIENT_2D_TOMITA_FLUID_H

int gradient_2d_tomita_fluid_d2(const int nop, const double * field,double * t_field,
				double * grad,double * t_grad, double * delsq, double * t_delsq);

int gradient_2d_tomita_fluid_d4(const int nop, const double * field,double * t_field,
				double * grad,double * t_grad, double * delsq, double * t_delsq);


#endif
