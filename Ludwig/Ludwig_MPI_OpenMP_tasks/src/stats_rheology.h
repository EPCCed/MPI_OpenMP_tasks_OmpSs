/*****************************************************************************
 *
 *  stats_rheology.h
 *
 *  Edinburgh Soft Matter and Statistcal Physics Group and
 *  Edinburgh Parallel Computing Centre
 *
 *  $Id: stats_rheology.h,v 1.3 2009-10-14 17:16:01 kevin Exp $
 *
 *  Kevin Stratford (kevin@epcc.ed.ac.uk)
 *  (c) The University of Edinburgh (2009)
 *
 *****************************************************************************/

#ifndef STATS_RHEOLOGY_H
#define STATS_RHEOLOGY_H

#include "model.h"
#include "hydro.h"

void stats_rheology_init(void);
void stats_rheology_finish(void);

int stats_rheology_stress_profile_accumulate(lb_t * lb, hydro_t * hydro);
int stats_rheology_mean_stress(lb_t * lb, const char * filename);

void stats_rheology_free_energy_density_profile(const char *);
void stats_rheology_stress_profile_zero(void);
void stats_rheology_stress_profile(const char *);
void stats_rheology_stress_section(const char *);

#endif
