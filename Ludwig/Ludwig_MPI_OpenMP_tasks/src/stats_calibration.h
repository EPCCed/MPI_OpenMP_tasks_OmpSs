/*****************************************************************************
 *
 *  stats_calibration.c
 *
 *  $Id: stats_calibration.h 2479 2014-09-25 11:24:42Z stratford $
 *
 *  Edinburgh Soft Matter and Statistical Physics Group and
 *  Edinburgh Parallel Computing Centre
 *
 *  Kevin Stratford (kevin@epcc.ed.ac.uk)
 *  (c) 2011 The University of Edinburgh
 *
 *****************************************************************************/

#ifndef STATS_CALIBRATION_H
#define STATS_CALIBRATION_H

#include "colloids.h"
#include "hydro.h"
#include "map.h"

int stats_calibration_init(colloids_info_t * cinfo, int nswitch);
int stats_calibration_accumulate(colloids_info_t * cinfo, int ntimestep,
				 hydro_t * hydro, map_t * map);
int stats_calibration_finish(void);

#endif
