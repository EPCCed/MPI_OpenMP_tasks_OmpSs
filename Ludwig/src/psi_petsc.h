#ifdef PETSC
/*****************************************************************************
 *
 *  psi_petsc.h
 *
 *  $Id: psi_petsc.h 2479 2014-09-25 11:24:42Z stratford $
 *
 *  Edinburgh Soft Matter and Statistical Physics Group and
 *  Edinburgh Parallel Computing Centre
 *
 *  (c) 2012-2013 The University of Edinburgh
 *
 *  Contributing Authors:
 *    Oliver Henrich (ohenrich@epcc.ed.ac.uk)
 *    Kevin Stratford (kevin@epcc.ed.ac.uk)
 *
 *****************************************************************************/

#ifndef PSI_PETSC_H
#define PSI_PETSC_H

#include "psi.h"

int psi_petsc_init(psi_t * obj, f_vare_t fepsilon);
int psi_petsc_copy_psi_to_da(psi_t * obj);
int psi_petsc_copy_da_to_psi(psi_t * obj);
int psi_petsc_compute_laplacian(psi_t * obj);
int psi_petsc_set_rhs(psi_t * obj);
int psi_petsc_compute_matrix(psi_t * obj, f_vare_t fepsilon);
int psi_petsc_set_rhs_vare(psi_t * obj, f_vare_t epsilon);
int psi_petsc_solve(psi_t * obj, f_vare_t fepsilon);
int psi_petsc_poisson(psi_t * obj);
int psi_petsc_finish();

#endif
#endif

