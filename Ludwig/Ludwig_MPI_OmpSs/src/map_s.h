/*****************************************************************************
 *
 *  map_s.h
 *
 *  $Id: map_s.h 2644 2015-05-11 09:16:17Z agray3 $
 *
 *  Edinburgh Soft Matter and Statistical Physics Group and
 *  Edinburgh Parallel Computing Centre
 *
 *  Kevin Stratford (kevin@epcc.ed.ac.uk)
 *  (c) 2012 The University of Edinburgh
 *
 *****************************************************************************/

#ifndef MAP_S_H
#define MAP_S_H

#include <mpi.h>
#include "map.h"

struct map_s {
  int is_porous_media;        /* Flag for porous media */
  int ndata;                  /* Additional fields associated with map */
  char * status;              /* Status (one of enum_status) */
  char * t_status;              /* Status (one of enum_status) on target */
  double * data;              /* Additional site lattice property */
  MPI_Datatype halostatus[3]; /* Halo datatype for status */
  MPI_Datatype halodata[3];   /* Halo datatype for data */
  io_info_t * info;           /* I/O handler */

  map_t * tcopy;                /* Copy of this structure on target*/

};

#endif
