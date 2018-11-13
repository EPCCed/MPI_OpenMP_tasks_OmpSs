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



/* Send edges parallel to x-direction */
int prepareBuffersX(lb_t * lb, double* fptr, int* nlocal, int id){

  int count,p,n,ic;
  int index,indexreal;

  count = 0;
  for (p = 0; p < NVEL; p++) {
    for (n = 0; n < lb->ndist; n++) {
      for (ic = 1; ic <= nlocal[X]; ic++) {

        if(id==6){
          index = coords_index(ic, nlocal[Y], nlocal[Z]);
          indexreal = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
          lb->hl .sendXpp[count] = fptr[indexreal];
        }
        if (id==7){
          index = coords_index(ic, nlocal[Y], 1);
          indexreal = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
          lb->hl .sendXpn[count] = fptr[indexreal];
        }
        if (id==8){
          index = coords_index(ic, 1, nlocal[Z]);
          indexreal = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
          lb->hl .sendXnp[count] = fptr[indexreal];
        }
        if (id==9){
          index = coords_index(ic, 1, 1);
          indexreal = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
          lb->hl .sendXnn[count] = fptr[indexreal];
        }

        ++count;
      }
    }
  }
  return count;

}

/* Send edges parallel to y-direction */
int prepareBuffersY(lb_t * lb, double* fptr, int* nlocal, int id){

  int count,p,n,jc;
  int index,indexreal;

  count = 0;
  for (p = 0; p < NVEL; p++) {
    for (n = 0; n < lb->ndist; n++) {
      for (jc = 1; jc <= nlocal[Y]; jc++) {

        if(id==10){
          index = coords_index(nlocal[X], jc, nlocal[Z]);
          indexreal = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
          lb->hl .sendYpp[count] = fptr[indexreal];
        }
        if(id==11){
          index = coords_index(nlocal[X], jc, 1);
          indexreal = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
          lb->hl .sendYpn[count] = fptr[indexreal];
        }
        if(id==12){
          index = coords_index(1, jc, nlocal[Z]);
          indexreal = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
          lb->hl .sendYnp[count] = fptr[indexreal];
        }

        if(id==13){
          index = coords_index(1, jc, 1);
          indexreal = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
          lb->hl .sendYnn[count] = fptr[indexreal];
        }

        ++count;
      }
    }
  }


  return count;

}

/* Send edges parallel to z-direction */
int prepareBuffersZ(lb_t * lb, double* fptr, int* nlocal, int id){
  
  int count,p,n,kc;
  int index,indexreal;
  
  count = 0;
  for (p = 0; p < NVEL; p++) {
    for (n = 0; n < lb->ndist; n++) {
      for (kc = 1; kc <= nlocal[Z]; kc++) {
	
	if (id==14){
	  index = coords_index(nlocal[X], nlocal[Y], kc);
	  indexreal = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
	  lb->hl .sendZpp[count] = fptr[indexreal];
	}
	if(id==15){
	  index = coords_index(nlocal[X], 1, kc);
	  indexreal = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
	  lb->hl .sendZpn[count] = fptr[indexreal];
	}
	
	if (id==16){
	  index = coords_index(1, nlocal[Y], kc);
	  indexreal = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
	  lb->hl .sendZnp[count] = fptr[indexreal];
	}
	if(id==17){
	  index = coords_index(1, 1, kc);
	  indexreal = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
	  lb->hl .sendZnn[count] = fptr[indexreal];
	}
	
	++count;
      }
    }
  }
  
  return count;
}


  

int edgeBuffers(lb_t * lb, double* fptr,int* nlocal, int id){

  int count,p,n;
  int index,indexreal;
  
  count = 0;
  for (p = 0; p < NVEL; p++) {
    for (n = 0; n < lb->ndist; n++) {

      if(id==18){
	index = coords_index(nlocal[X], nlocal[Y], nlocal[Z]);
	indexreal = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
	lb->hl .sendppp[count] = fptr[indexreal];
      }
      if(id==19){
	index = coords_index(nlocal[X], nlocal[Y], 1);
	indexreal = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
	lb->hl .sendppn[count] = fptr[indexreal];
      }
      if(id==20){
	index = coords_index(nlocal[X], 1, nlocal[Z]);
	indexreal = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
	lb->hl .sendpnp[count] = fptr[indexreal];
      }
      if(id==21){
	index = coords_index(nlocal[X], 1, 1);
	indexreal = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
	lb->hl .sendpnn[count] = fptr[indexreal];
      }
      if(id==22){
	index = coords_index(1, nlocal[Y], nlocal[Z]);
	indexreal = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
	lb->hl .sendnpp[count] = fptr[indexreal];
      }
      if(id==23){
	index = coords_index(1, nlocal[Y], 1);
	indexreal = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
	lb->hl .sendnpn[count] = fptr[indexreal];
      }
      if(id==24){
	index = coords_index(1, 1, nlocal[Z]);
	indexreal = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
	lb->hl .sendnnp[count] = fptr[indexreal];
      }
      if(id==25){
	index = coords_index(1, 1, 1);
	indexreal = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
	lb->hl .sendnnn[count] = fptr[indexreal];
      }
     
      count++;
    }
  }

  return count;

}
