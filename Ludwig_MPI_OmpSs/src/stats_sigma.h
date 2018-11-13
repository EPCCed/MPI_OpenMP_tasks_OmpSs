/*****************************************************************************
 *
 *  stats_sigma.h
 *
 *  $Id: stats_sigma.h 2479 2014-09-25 11:24:42Z stratford $
 *
 *  Edinburgh Soft Matter and Statistical Physics Group and
 *  Edinburgh Parallel Computing Centre
 *
 *  Kevin Stratford (kevin@epcc.ed.ac.uk)
 *  (c) 2011 The University of Edinburgh
 *
 *****************************************************************************/

#ifndef STATS_SIGMA_H
#define STATS_SIGMA_H

#include "field.h"

int stats_sigma_init(field_t * phi, int nswtich);
int stats_sigma_measure(field_t * phi, int ntime);

#endif
