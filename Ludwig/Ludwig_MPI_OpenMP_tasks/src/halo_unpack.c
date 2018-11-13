
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "pe.h"
#include "coords.h"
#include "model.h"
#include "lb_model_s.h"
#include "targetDP.h"
#include "io_harness.h"
#include "control.h"
#include "unpack.h"



/*****************************************************************************
 *
 *  unpack_halo_buffers
 *
 *  Unpacks buffers from MPI messages into Halo values for next step of simulation
 *
 *****************************************************************************/

void unpack_halo_buffers_i(lb_t * lb, int id) {

  int i,ic, jc, kc;
  int n, p;
  int index, indexhalo;
  int count;
  int nlocal[3];

  coords_nlocal(nlocal);

  double* fptr; /*pointer to the distribution*/
  if (get_step()) /* we are in the timestep, so use target copy */
    fptr=lb->tcopy->f;
  else
    fptr=lb->f;  /* we are not in a timestep, use host copy */

  int nsendXZ, nsendYZ, nsendXY;
  nsendYZ = NVEL*lb->ndist*nlocal[Y]*nlocal[Z];
  nsendXZ = NVEL*lb->ndist*nlocal[X]*nlocal[Z];
  nsendXY = NVEL*lb->ndist*nlocal[X]*nlocal[Y];

  int nsendX, nsendY, nsendZ;
  nsendX = NVEL*lb->ndist*nlocal[X];
  nsendY = NVEL*lb->ndist*nlocal[Y];
  nsendZ = NVEL*lb->ndist*nlocal[Z];

  int nsend;
  nsend = NVEL*lb->ndist*1;

  /* Unpack Planes */
  if (id>=0 && id < 6){ 
    unpack_planes_i(lb, nlocal, fptr, id);
    free_planes_i(lb, id);
  }
  
  
  /* Unpack Edges */
  if(id>=6 && id <18){
    unpack_edges_i(lb, nlocal, fptr, id); 
    free_edges_i(lb, id);
  }
  

  /* Unpack corners buffers */
  if(id >=18 && id < 26){
    unpack_corners_i(lb, nlocal, fptr, id);
    free_corners_i(lb, id);
  }
  return;
}
