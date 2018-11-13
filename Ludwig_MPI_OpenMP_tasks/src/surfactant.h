/****************************************************************************
 *
 *  surfactant.h
 *
 *  $Id: surfactant.h,v 1.2 2010-10-15 12:40:03 kevin Exp $
 *
 *  Edinburgh Soft Matter and Statistical Physics Group
 *  and Edinburgh Parallel Computing Centre
 *
 *  Kevin Stratford (kevin@epcc.ed.ac.uk)
 *  (c) The University of Edinburgh (2009)
 *
 ****************************************************************************/

#ifndef SURFACTANT_H
#define SURFACTANT_H

#include "field.h"
#include "field_grad.h"

int fe_surfactant_create(field_t * phi, field_grad_t * grad);
void fe_surfactant_free(void);

int surfactant_fluid_parameters_set(double a, double b, double kappa);
int surfactant_parameters_set(double kt, double eps, double beta, double w);

double surfactant_interfacial_tension(void);
double surfactant_interfacial_width(void);
double surfactant_langmuir_isotherm(void);
double surfactant_free_energy_density(const int index);
double surfactant_chemical_potential(const int index, const int nop);
double surfactant_isotropic_pressure(const int index);
void   surfactant_chemical_stress(const int index, double s[3][3]);

#endif

