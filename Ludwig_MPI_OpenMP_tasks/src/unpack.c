
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


void unpack_planes(lb_t * lb, int* nlocal, double* fptr ){

  int count,p,n;
  int ic,jc,kc;
  int index, indexhalo;

  /* X */

  count = 0;
  for (p = 0; p < NVEL; p++) {
    for (n = 0; n < lb->ndist; n++) {
      for (jc = 1; jc <= nlocal[Y]; jc++) {
        for (kc = 1; kc <= nlocal[Z]; kc++) {

          index = coords_index(0, jc, kc);
          indexhalo = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
          fptr[indexhalo] = lb->hl.recvbackYZ[count];

          index = coords_index(nlocal[X] + 1, jc, kc);
          indexhalo = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
          fptr[indexhalo] = lb->hl.recvforwYZ[count];
          ++count;
        }
      }
    }
  }
  //assert(count == nsendYZ);

  /* Y */
  count = 0;
  for (p = 0; p < NVEL; p++) {
    for (n = 0; n < lb->ndist; n++) {
      for (ic = 1; ic <= nlocal[X]; ic++) {
        for (kc = 1; kc <= nlocal[Z]; kc++) {

          index = coords_index(ic, 0, kc);
          indexhalo = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
          fptr[indexhalo] = lb->hl.recvbackXZ[count];

          index = coords_index(ic, nlocal[Y] + 1, kc);
          indexhalo = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
          fptr[indexhalo] = lb->hl.recvforwXZ[count];
          ++count;
        }
      }
    }
  }
  //assert(count == nsendXZ);

  /* Z  */
  count = 0;
  for (p = 0; p < NVEL; p++) {
    for (n = 0; n < lb->ndist; n++) {
      for (ic = 1; ic <= nlocal[X]; ic++) {
        for (jc = 1; jc <= nlocal[Y]; jc++) {

          index = coords_index(ic, jc, 0);
          indexhalo = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
          fptr[indexhalo] = lb->hl.recvbackXY[count];

          index = coords_index(ic, jc, nlocal[Z] + 1);
          indexhalo = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
          fptr[indexhalo] = lb->hl.recvforwXY[count];
          ++count;
        }
      }
    }
  }
  //assert(count == nsendXY);


}


void unpack_planes_i(lb_t * lb, int* nlocal, double* fptr, int id ){

  int count,p,n;
  int ic,jc,kc;
  int index, indexhalo;

  /* X */
  if (id == 0 || id ==1){
    count = 0;
    for (p = 0; p < NVEL; p++) {
      for (n = 0; n < lb->ndist; n++) {
        for (jc = 1; jc <= nlocal[Y]; jc++) {
          for (kc = 1; kc <= nlocal[Z]; kc++) {
            if (id == 1 ){
              index = coords_index(0, jc, kc);
              indexhalo = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
              fptr[indexhalo] = lb->hl.recvbackYZ[count];

            }
            if ( id ==0 ){
              index = coords_index(nlocal[X] + 1, jc, kc);
              indexhalo = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
              fptr[indexhalo] = lb->hl.recvforwYZ[count];
            }

            ++count;
          }
        }
      }
    }
  }
  //assert(count == nsendYZ);

  /* Y */
  if (id == 2 || id ==3){
    count = 0;
    for (p = 0; p < NVEL; p++) {
      for (n = 0; n < lb->ndist; n++) {
        for (ic = 1; ic <= nlocal[X]; ic++) {
          for (kc = 1; kc <= nlocal[Z]; kc++) {
            if(id==3) {
              index = coords_index(ic, 0, kc);
              indexhalo = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
              fptr[indexhalo] = lb->hl.recvbackXZ[count];
            }
            if(id==2) {
              index = coords_index(ic, nlocal[Y] + 1, kc);
              indexhalo = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
              fptr[indexhalo] = lb->hl.recvforwXZ[count];
            }
            ++count;
          }
        }
      }
    }
    //assert(count == nsendXZ);
  }

  /* Z  */

  if (id == 4  || id == 5){
    count = 0;
    for (p = 0; p < NVEL; p++) {
      for (n = 0; n < lb->ndist; n++) {
        for (ic = 1; ic <= nlocal[X]; ic++) {
          for (jc = 1; jc <= nlocal[Y]; jc++) {
            if (id ==5){
              index = coords_index(ic, jc, 0);
              indexhalo = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
              fptr[indexhalo] = lb->hl.recvbackXY[count];
            }
            if(id==4){
              index = coords_index(ic, jc, nlocal[Z] + 1);
              indexhalo = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
              fptr[indexhalo] = lb->hl.recvforwXY[count];
            }

            ++count;

          }
        }
      }
    }
    //assert(count == nsendXY);
  }

}

void unpack_edges(lb_t * lb, int* nlocal, double* fptr){
  int count,p,n;
  int ic,jc,kc;
  int index, indexhalo;
  /* X */
  count = 0;
  for (p = 0; p < NVEL; p++) {
    for (n = 0; n < lb->ndist; n++) {
      for (ic = 1; ic <= nlocal[X]; ic++) {

        index = coords_index(ic, nlocal[Y] + 1, nlocal[Z] + 1);
        indexhalo = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
        fptr[indexhalo] = lb->hl .recvXpp[count];

        index = coords_index(ic, nlocal[Y] + 1, 0);
        indexhalo = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
        fptr[indexhalo] = lb->hl .recvXpn[count];

        index = coords_index(ic, 0, nlocal[Z] + 1);
        indexhalo = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
        fptr[indexhalo] = lb->hl .recvXnp[count];

        index = coords_index(ic, 0, 0);
        indexhalo = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
        fptr[indexhalo] = lb->hl .recvXnn[count];
        ++count;
      }
    }
  }
  //assert(count == nsendX);
  /* Y */
  count = 0;
  for (p = 0; p < NVEL; p++) {
    for (n = 0; n < lb->ndist; n++) {
      for (jc = 1; jc <= nlocal[Y]; jc++) {

        index = coords_index(nlocal[X] + 1, jc, nlocal[Z] + 1);
        indexhalo = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
        fptr[indexhalo] = lb->hl .recvYpp[count];

        index = coords_index(nlocal[X] + 1, jc, 0);
        indexhalo = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
        fptr[indexhalo] = lb->hl .recvYpn[count];

        index = coords_index(0, jc, nlocal[Z] + 1);
        indexhalo = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
        fptr[indexhalo] = lb->hl .recvYnp[count];

        index = coords_index(0, jc, 0);
        indexhalo = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
        fptr[indexhalo] = lb->hl .recvYnn[count];
        ++count;
      }
    }
  }
  //assert(count == nsendY);
  /* Z */
  count = 0;
  for (p = 0; p < NVEL; p++) {
    for (n = 0; n < lb->ndist; n++) {
      for (kc = 1; kc <= nlocal[Z]; kc++) {

        index = coords_index(nlocal[X] + 1, nlocal[Y] + 1, kc);
        indexhalo = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
        fptr[indexhalo] = lb->hl .recvZpp[count];

        index = coords_index(nlocal[X] + 1, 0, kc);
        indexhalo = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
        fptr[indexhalo] = lb->hl .recvZpn[count];

        index = coords_index(0, nlocal[Y] + 1, kc);
        indexhalo = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
        fptr[indexhalo] = lb->hl .recvZnp[count];

        index = coords_index(0, 0, kc);
        indexhalo = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
        fptr[indexhalo] = lb->hl .recvZnn[count];
        ++count;
      }
    }
  }
  // assert(count == nsendZ);
}


void unpack_edges_i(lb_t * lb, int* nlocal, double* fptr, int id){
  int count,p,n;
  int ic,jc,kc;
  int index, indexhalo;

  /* X */
  if (id >= 6 && id < 10 ){

    count = 0;
    for (p = 0; p < NVEL; p++) {
      for (n = 0; n < lb->ndist; n++) {
        for (ic = 1; ic <= nlocal[X]; ic++) {

          if ( id == 6) {
            index = coords_index(ic, 0, 0);
            indexhalo = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
            fptr[indexhalo] = lb->hl .recvXnn[count];
          }
          if ( id == 7){
            index = coords_index(ic, 0, nlocal[Z] + 1);
            indexhalo = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
            fptr[indexhalo] = lb->hl .recvXnp[count];
          }
          if ( id == 8){
            index = coords_index(ic, nlocal[Y] + 1, 0);
            indexhalo = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
            fptr[indexhalo] = lb->hl .recvXpn[count];
          }
          if (id == 9){
            index = coords_index(ic, nlocal[Y] + 1, nlocal[Z] + 1);
            indexhalo = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
            fptr[indexhalo] = lb->hl .recvXpp[count];
          }

          ++count;
        }
      }
    }

  }


  //assert(count == nsendX);
  /* Y */
  if ( id >= 10 && id < 14){
    count = 0;
    for (p = 0; p < NVEL; p++) {
      for (n = 0; n < lb->ndist; n++) {
        for (jc = 1; jc <= nlocal[Y]; jc++) {

          if (id == 10){
            index = coords_index(0, jc, 0);
            indexhalo = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
            fptr[indexhalo] = lb->hl .recvYnn[count];
          }
          if (id == 11){
            index = coords_index(0, jc, nlocal[Z] + 1);
            indexhalo = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
            fptr[indexhalo] = lb->hl .recvYnp[count];
          }
          if (id == 12){
            index = coords_index(nlocal[X] + 1, jc, 0);
            indexhalo = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
            fptr[indexhalo] = lb->hl .recvYpn[count];
          }
          if (id == 13){
            index = coords_index(nlocal[X] + 1, jc, nlocal[Z] + 1);
            indexhalo = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
            fptr[indexhalo] = lb->hl .recvYpp[count];
          }

          ++count;
        }
      }
    }
  }
  //assert(count == nsendY);
  /* Z */
  if ( id >= 14 && id < 18 ){
    count = 0;
    for (p = 0; p < NVEL; p++) {
      for (n = 0; n < lb->ndist; n++) {
        for (kc = 1; kc <= nlocal[Z]; kc++) {

          if (id == 14){
            index = coords_index(0, 0, kc);
            indexhalo = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
            fptr[indexhalo] = lb->hl .recvZnn[count];
          }
          if (id == 15){
            index = coords_index(0, nlocal[Y] + 1, kc);
            indexhalo = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
            fptr[indexhalo] = lb->hl .recvZnp[count];
          }
          if (id == 16){
            index = coords_index(nlocal[X] + 1, 0, kc);
            indexhalo = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
            fptr[indexhalo] = lb->hl .recvZpn[count];
          }
          if (id == 17){
            index = coords_index(nlocal[X] + 1, nlocal[Y] + 1, kc);
            indexhalo = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
            fptr[indexhalo] = lb->hl .recvZpp[count];
          }

          ++count;
        }
      }
    }
    // assert(count == nsendZ);
  }

}

void unpack_corners(lb_t * lb, int* nlocal, double* fptr){

  int count,p,n;
  int ic,jc,kc;
  int index, indexhalo;

  count = 0;
  for (p = 0; p < NVEL; p++) {
    for (n = 0; n < lb->ndist; n++) {

      index = coords_index(nlocal[X] + 1, nlocal[Y] + 1, nlocal[Z] + 1);
      indexhalo = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
      fptr[indexhalo] = lb->hl .recvppp[count];

      index = coords_index(nlocal[X] + 1, nlocal[Y] + 1, 0);
      indexhalo = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
      fptr[indexhalo] = lb->hl .recvppn[count];

      index = coords_index(nlocal[X] + 1, 0, nlocal[Z] + 1);
      indexhalo = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
      fptr[indexhalo] = lb->hl .recvpnp[count];

      index = coords_index(nlocal[X] + 1, 0, 0);
      indexhalo = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
      fptr[indexhalo] = lb->hl .recvpnn[count];

      index = coords_index(0, nlocal[Y] + 1, nlocal[Z] + 1);
      indexhalo = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
      fptr[indexhalo] = lb->hl .recvnpp[count];

      index = coords_index(0, nlocal[Y] + 1, 0);
      indexhalo = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
      fptr[indexhalo] = lb->hl .recvnpn[count];

      index = coords_index(0, 0, nlocal[Z] + 1);
      indexhalo = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
      fptr[indexhalo] = lb->hl .recvnnp[count];

      index = coords_index(0, 0, 0);
      indexhalo = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
      fptr[indexhalo] = lb->hl .recvnnn[count];
      count++;
    }
  }
}



void unpack_corners_i(lb_t * lb, int* nlocal, double* fptr, int id){

  int count,p,n;
  int ic,jc,kc;
  int index, indexhalo;

  if ( id >= 18 && id < 26){

    count = 0;
    for (p = 0; p < NVEL; p++) {
      for (n = 0; n < lb->ndist; n++) {
        if(id == 18 ){
          index = coords_index(0, 0, 0);
          indexhalo = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
          fptr[indexhalo] = lb->hl .recvnnn[count];
        }
        if (id == 19){
          index = coords_index(0, 0, nlocal[Z] + 1);
          indexhalo = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
          fptr[indexhalo] = lb->hl .recvnnp[count];
        }
        if (id == 20){
          index = coords_index(0, nlocal[Y] + 1, 0);
          indexhalo = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
          fptr[indexhalo] = lb->hl .recvnpn[count];
        }
        if (id == 21){
          index = coords_index(0, nlocal[Y] + 1, nlocal[Z] + 1);
          indexhalo = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
          fptr[indexhalo] = lb->hl .recvnpp[count];
        }
        if (id == 22){
          index = coords_index(nlocal[X] + 1, 0, 0);
          indexhalo = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
          fptr[indexhalo] = lb->hl .recvpnn[count];
        }
        if ( id == 23){
          index = coords_index(nlocal[X] + 1, 0, nlocal[Z] + 1);
          indexhalo = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
          fptr[indexhalo] = lb->hl .recvpnp[count];
        }
        if (id == 24){
          index = coords_index(nlocal[X] + 1, nlocal[Y] + 1, 0);
          indexhalo = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
          fptr[indexhalo] = lb->hl .recvppn[count];
        }
        if (id == 25){
          index = coords_index(nlocal[X] + 1, nlocal[Y] + 1, nlocal[Z] + 1);
          indexhalo = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
          fptr[indexhalo] = lb->hl .recvppp[count];
        }
        count++;
      }
    }
  }

}

void free_planes(lb_t * lb){

  free(lb->hl .sendforwYZ);
  free(lb->hl .recvforwYZ);
  free(lb->hl .sendbackYZ);
  free(lb->hl .recvbackYZ);
  free(lb->hl .sendforwXZ);
  free(lb->hl .recvforwXZ);
  free(lb->hl .sendbackXZ);
  free(lb->hl .recvbackXZ);
  free(lb->hl .sendforwXY);
  free(lb->hl .recvforwXY);
  free(lb->hl .sendbackXY);
  free(lb->hl .recvbackXY);

}
void free_planes_i(lb_t * lb, int id){

  switch (id){
  case 0:
    free(lb->hl .recvforwYZ);
    free(lb->hl .sendbackYZ);
    break;
  case 1:
    free(lb->hl .recvbackYZ);
    free(lb->hl .sendforwYZ);
    break;
  case 2:
    free(lb->hl .recvforwXZ);
    free(lb->hl .sendbackXZ);
    break;
  case 3:
    free(lb->hl .recvbackXZ);
    free(lb->hl .sendforwXZ);
    break;
  case 4:
    free(lb->hl .recvforwXY);
    free(lb->hl .sendbackXY);
    break;
  case 5:
    free(lb->hl .recvbackXY);
    free(lb->hl .sendforwXY);
    break;
  default:
    exit(0);
  }

}

void free_edges(lb_t * lb){
  free(lb->hl .sendXnn);
  free(lb->hl .recvXnn);
  free(lb->hl .sendXnp);
  free(lb->hl .recvXnp);
  free(lb->hl .sendXpn);
  free(lb->hl .recvXpn);
  free(lb->hl .sendXpp);
  free(lb->hl .recvXpp);
  free(lb->hl .sendYnn);
  free(lb->hl .recvYnn);
  free(lb->hl .sendYnp);
  free(lb->hl .recvYnp);

  free(lb->hl .sendYpn);
  free(lb->hl .recvYpn);
  free(lb->hl .sendYpp);
  free(lb->hl .recvYpp);
  free(lb->hl .sendZnn);
  free(lb->hl .recvZnn);
  free(lb->hl .sendZnp);
  free(lb->hl .recvZnp);
  free(lb->hl .sendZpn);
  free(lb->hl .recvZpn);
  free(lb->hl .sendZpp);
  free(lb->hl .recvZpp);
}


void free_edges_i(lb_t * lb, int id){

  switch (id){
  case 6:
    free(lb->hl .recvXnn);
    free(lb->hl .sendXpp);
    break;
  case 7:
    free(lb->hl .recvXnp);
    free(lb->hl .sendXpn);
    break;
  case 8:
    free(lb->hl .recvXpn);
    free(lb->hl .sendXnp);
    break;
  case 9:
    free(lb->hl .recvXpp);
    free(lb->hl .sendXnn);
    break;
  case 10:
    free(lb->hl .recvYnn);
    free(lb->hl .sendYpp);
    break;
  case 11:
    free(lb->hl .recvYnp);
    free(lb->hl .sendYpn);
    break;
  case 12:
    free(lb->hl .recvYpn);
    free(lb->hl .sendYnp);
    break;
  case 13:
    free(lb->hl .recvYpp);
    free(lb->hl .sendYnn);
    break;
  case 14:
    free(lb->hl .recvZnn);
    free(lb->hl .sendZpp);
    break;
  case 15:
    free(lb->hl .recvZnp);
    free(lb->hl .sendZpn);
    break;
  case 16:
    free(lb->hl .recvZpn);
    free(lb->hl .sendZnp);
    break;
  case 17:
    free(lb->hl .recvZpp);
    free(lb->hl .sendZnn);
    break;

  default:
    exit(0);
  }
}


void free_corners(lb_t * lb){
  free(lb->hl .sendnnn);
  free(lb->hl .sendnnp);
  free(lb->hl .sendnpn);
  free(lb->hl .sendnpp);
  free(lb->hl .sendpnn);
  free(lb->hl .sendpnp);
  free(lb->hl .sendppn);
  free(lb->hl .sendppp);

  free(lb->hl .recvnnn);
  free(lb->hl .recvnnp);
  free(lb->hl .recvnpn);
  free(lb->hl .recvnpp);
  free(lb->hl .recvpnn);
  free(lb->hl .recvpnp);
  free(lb->hl .recvppn);
  free(lb->hl .recvppp);
}

void free_corners_i(lb_t * lb, int id){
  
  switch(id){
  case 18:
    free(lb->hl .recvnnn);
    free(lb->hl .sendppp);
    break;
  case 19:
    free(lb->hl .recvnnp);
    free(lb->hl .sendppn);
    break;
  case 20:
    free(lb->hl .recvnpn);
    free(lb->hl .sendpnp);
    break;
  case 21:
    free(lb->hl .recvnpp);
    free(lb->hl .sendpnn);
    break;
  case 22:
    free(lb->hl .recvpnn);
    free(lb->hl .sendnpp);
    break;
  case 23:
    free(lb->hl .recvpnp);
    free(lb->hl .sendnpn);
    break;
  case 24:
    free(lb->hl .recvppn);
    free(lb->hl .sendnnp);
    break;
  case 25:
    free(lb->hl .recvppp);
    free(lb->hl .sendnnn);
    break;
  default:
    exit(0);
  }
}
