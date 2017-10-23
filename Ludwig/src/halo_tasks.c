
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


/*****************************************************************************
 *
 *  lb_halo_via_copy_nonblocking_tasks
 *
 *  A version of the halo swap which uses a flat buffer to copy the
 *  relevant data rathar than MPI data types. This version is NON-
 *  BLOCKING and sends 26 messages instead of 6: 6 planes, 12 edges
 *  and 8 corners.
 *
 *  It works for both MODEL and MODEL_R, but the loop order favours
 *  MODEL_R.
 *
 *****************************************************************************/

int lb_halo_via_copy_nonblocking_tasks_start(lb_t * lb) {

  /* Send messages to 6 neighbouring planes
     12 neighbouring edges
     8 neighbouring corners */



  int id;
  /* Create 6 taks, one per plane */
  for (id = 0; id < 6; id++)
#pragma omp task default(none) shared(lb)firstprivate(id)       \
  depend(out:lb->hl.recvreq[id],lb->hl.sendreq[id])
    halo_planes_tasks_i(lb, id);


  for (id = 6; id < 18; id++)
#pragma omp task default(none) shared(lb)firstprivate(id)       \
  depend(out:lb->hl.recvreq[id],lb->hl.sendreq[id])
    halo_edges_tasks_i(lb, id);



  for (id = 18; id < 26; id++){
#pragma omp task  default(none) shared(lb)firstprivate(id)      \
  depend(out:lb->hl.recvreq[id],lb->hl.sendreq[id])
    halo_corners_tasks_i(lb,id);
  }

  //#pragma omp taskwait

    return 0;
}

int lb_halo_via_copy_nonblocking_tasks_end(lb_t * lb) {

  int n;
  int recvcount;
  int sendcount;

  for (n=0; n<26; n++){
#pragma omp task depend (inout:lb->hl.recvreq[n])
    MPI_Waitany(26, lb->hl.recvreq, &recvcount, lb->hl.status);
  }


  /* Copy data from MPI buffers */
  //TODO:Watch out! Not calling the _task version here!!
  //#pragma omp task depend (in: lb->hl.recvreq[0:25])
  for (n=0; n<26; n++){
#pragma omp task default(none) shared(lb)firstprivate(n)   \
  depend(out:lb->hl.recvreq[id])
    unpack_halo_buffers_i(lb,n);
  }
  
  //#pragma omp taskwait

  for (n=0; n<26; n++){
#pragma omp task depend (in:lb->hl.sendreq[n])
    MPI_Waitany(26, lb->hl.sendreq, &sendcount, lb->hl.status);
  }
#pragma omp taskwait

  return 0;
}


/*****************************************************************************
 *
 *  halo_planes_tasks_i
 *
 *  Sends/Recv 1 MPI messages per plane and direction
 *
 *****************************************************************************/

void halo_planes_tasks_i(lb_t * lb, int id) {

  int ic, jc, kc;
  int n, p;
  int index, indexhalo, indexreal;
  int count;
  int nlocal[3];

  const int tagf = 900;
  const int tagb = 901;

  /* The ranks of neighbouring planes */
  int pforwX, pbackX, pforwY, pbackY, pforwZ, pbackZ;

  MPI_Comm comm = cart_comm();

  assert(lb);

  coords_nlocal(nlocal);

  double* fptr; /*pointer to the distribution*/
  if (get_step()) /* we are in the timestep, so use target copy */
    fptr=lb->tcopy->f;
  else
    fptr=lb->f;  /* we are not in a timestep, use host copy */

  /* Allocate size of sendplanes and number of elements send with each plane */
  int nsendXZ, nsendYZ, nsendXY;

  nsendYZ = NVEL*lb->ndist*nlocal[Y]*nlocal[Z];
  nsendXZ = NVEL*lb->ndist*nlocal[X]*nlocal[Z];
  nsendXY = NVEL*lb->ndist*nlocal[X]*nlocal[Y];


  if (id == 0 || id ==1){
    int k;


    /* PPM, NMM, P=Positive, M=Middle, N=Negative for the XYZ directions respectively */
    pforwX = nonblocking_cart_neighb(PMM);
    pbackX = nonblocking_cart_neighb(NMM);


    /* Receive planes in the X-direction */
    if(id==0){
      lb->hl .recvforwYZ = (double *) malloc(nsendYZ*sizeof(double));
      lb->hl .sendbackYZ = (double *) malloc(nsendYZ*sizeof(double));

      MPI_Irecv(&lb->hl.recvforwYZ[0], nsendYZ, MPI_DOUBLE, pforwX, tagb+id+20+id, comm,
                &lb->hl.recvreq[0]);
    }  else{
      lb->hl .recvbackYZ = (double *) malloc(nsendYZ*sizeof(double));
      lb->hl .sendforwYZ = (double *) malloc(nsendYZ*sizeof(double));

      MPI_Irecv(&lb->hl .recvbackYZ[0], nsendYZ, MPI_DOUBLE, pbackX, tagf+id+20+id, comm,
                &lb->hl.recvreq[1]);
    }

    count = 0;
    for (p = 0; p < NVEL; p++) {
      for (n = 0; n < lb->ndist; n++) {
        for (jc = 1; jc <= nlocal[Y]; jc++) {
          for (kc = 1; kc <= nlocal[Z]; kc++) {
            if(id==1){
              index = coords_index(nlocal[X], jc, kc);
              indexreal = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
              lb->hl.sendforwYZ[count] = fptr[indexreal];
            }else if (id==0){
              index = coords_index(1, jc, kc);
              indexreal = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
              lb->hl.sendbackYZ[count] = fptr[indexreal];
            }
            count++;
          }
        }
      }
    }
    assert(count == nsendYZ);

    /* Send in the X-direction (YZ plane) */

    if(id==0)
      MPI_Issend(&lb->hl .sendbackYZ [0], nsendYZ, MPI_DOUBLE, pbackX, tagb+id+20+id, comm,
                 &lb->hl.sendreq[0]);
    else
      MPI_Issend(&lb->hl .sendforwYZ [0], nsendYZ, MPI_DOUBLE, pforwX, tagf+id+20+id, comm,
                 &lb->hl.sendreq[1]);

  }

  if (id == 2 || id ==3){

    /* Receive planes in the Y-direction */
    pforwY = nonblocking_cart_neighb(MPM);
    pbackY = nonblocking_cart_neighb(MNM);


    if(id==2) {

      lb->hl .sendbackXZ = (double *) malloc(nsendXZ*sizeof(double));
      lb->hl .recvforwXZ = (double *) malloc(nsendXZ*sizeof(double));
      MPI_Irecv(&lb->hl .recvforwXZ[0], nsendXZ, MPI_DOUBLE, pforwY, tagb+id+30, comm,
                &lb->hl.recvreq[2]);
    }else {

      lb->hl .sendforwXZ = (double *) malloc(nsendXZ*sizeof(double));
      lb->hl .recvbackXZ = (double *) malloc(nsendXZ*sizeof(double));
      MPI_Irecv(&lb->hl .recvbackXZ[0], nsendXZ, MPI_DOUBLE, pbackY, tagf+id+30, comm,
                &lb->hl.recvreq[3]);
    }


    /* Send in the Y-direction (XZ plane) */
    count = 0;
    for (p = 0; p < NVEL; p++) {
      for (n = 0; n < lb->ndist; n++) {
        for (ic = 1; ic <= nlocal[X]; ic++) {
          for (kc = 1; kc <= nlocal[Z]; kc++) {

            if(id==3) {
              index = coords_index(ic, nlocal[Y], kc);
              indexreal = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
              lb->hl .sendforwXZ[count] = fptr[indexreal];
            }else if(id==2){
              index = coords_index(ic, 1, kc);
              indexreal = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
              lb->hl .sendbackXZ[count] = fptr[indexreal];
            }

            ++count;
          }
        }
      }
    }
    assert(count == nsendXZ);
    if(id==2){
      MPI_Issend(&lb->hl .sendbackXZ[0], nsendXZ, MPI_DOUBLE, pbackY, tagb+id+30, comm,
                 &lb->hl.sendreq[2]);
    }else if(id==3){
      MPI_Issend(&lb->hl .sendforwXZ[0], nsendXZ, MPI_DOUBLE, pforwY, tagf+id+30, comm,
                 &lb->hl.sendreq[3]);
    }


  }
  if (id == 4  || id == 5){

    /* Receive planes in the Z-direction */
    pforwZ = nonblocking_cart_neighb(MMP);
    pbackZ = nonblocking_cart_neighb(MMN);

    if(id==4){
      lb->hl .sendbackXY = (double *) malloc(nsendXY*sizeof(double));
      lb->hl .recvforwXY = (double *) malloc(nsendXY*sizeof(double));

      MPI_Irecv(&lb->hl.recvforwXY[0], nsendXY, MPI_DOUBLE, pforwZ, tagb+id, comm,
                &lb->hl.recvreq[4]);
    }else{

      lb->hl .sendforwXY = (double *) malloc(nsendXY*sizeof(double));
      lb->hl .recvbackXY = (double *) malloc(nsendXY*sizeof(double));

      MPI_Irecv(&lb->hl.recvbackXY[0], nsendXY, MPI_DOUBLE, pbackZ, tagf+id, comm,
                &lb->hl.recvreq[5]);
    }

    /* Finally, Send in the Z-direction (XY plane) */
    count = 0;
    for (p = 0; p < NVEL; p++) {
      for (n = 0; n < lb->ndist; n++) {
        for (ic = 1; ic <= nlocal[X]; ic++) {
          for (jc = 1; jc <= nlocal[Y]; jc++) {
            if(id==5){
              index = coords_index(ic, jc, nlocal[Z]);
              indexreal = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
              lb->hl .sendforwXY[count] = fptr[indexreal];
            }else if (id==4){
              index = coords_index(ic, jc, 1);
              indexreal = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
              lb->hl .sendbackXY[count] = fptr[indexreal];
            }
            ++count;
          }
        }
      }
    }

    assert(count == nsendXY);
    if(id==4){
      MPI_Issend(&lb->hl .sendbackXY[0], nsendXY, MPI_DOUBLE, pbackZ, tagb+id, comm,
                 &lb->hl.sendreq[4]);
    }else{
      MPI_Issend(&lb->hl .sendforwXY[0], nsendXY, MPI_DOUBLE, pforwZ, tagf+id, comm,
                 &lb->hl.sendreq[5]);
    }


  }

  return;

}


/*****************************************************************************
 *
 *  halo_edges
 *
 *  Sends 12 MPI messages (the edges) for the Non-Blocking version.
 *
 *****************************************************************************/


void halo_edges_tasks_i(lb_t * lb, int id) {

  int ic, jc, kc;
  int n, p;
  int index, indexhalo, indexreal;
  int count;
  int nlocal[3];

  const int tagnn = 903;
  const int tagnp = 904;
  const int tagpn = 905;
  const int tagpp = 906;

  /* Ranks of neighbouring edges.
     Xpp, X is the direction, parallel to which edges are
     sent pp refers to the YZ directions respectively*/
  int Xpp, Xpn, Xnp, Xnn, Ypp, Ypn, Ynp, Ynn, Zpp, Zpn, Znp, Znn;

  MPI_Comm comm = cart_comm();

  assert(lb);

  coords_nlocal(nlocal);

  double* fptr; /*pointer to the distribution*/
  if (get_step()) /* we are in the timestep, so use target copy */
    fptr=lb->tcopy->f;
  else
    fptr=lb->f;  /* we are not in a timestep, use host copy */

  int nsendX, nsendY, nsendZ;
  nsendX = NVEL*lb->ndist*nlocal[X];
  nsendY = NVEL*lb->ndist*nlocal[Y];
  nsendZ = NVEL*lb->ndist*nlocal[Z];


  if (id >= 6 && id < 10 ){

    /* Receive edges parallel to x-direction */
    Xnn = nonblocking_cart_neighb(MNN);
    Xnp = nonblocking_cart_neighb(MNP);
    Xpn = nonblocking_cart_neighb(MPN);
    Xpp = nonblocking_cart_neighb(MPP);

    switch(id) {
    case 6 :
      lb->hl .sendXpp = (double *) malloc(nsendX*sizeof(double));
      lb->hl .recvXnn = (double *) malloc(nsendX*sizeof(double));
      MPI_Irecv(&lb->hl.recvXnn[0], nsendX, MPI_DOUBLE, Xnn, tagpp, comm, &lb->hl.recvreq[6]);
      count = prepareBuffersX(lb, fptr, nlocal, id);
      MPI_Issend(&lb->hl .sendXpp[0], nsendX, MPI_DOUBLE, Xpp, tagpp, comm, &lb->hl.sendreq[6]);
      break;
    case 7 :
      lb->hl .recvXnp = (double *) malloc(nsendX*sizeof(double));
      lb->hl .sendXpn = (double *) malloc(nsendX*sizeof(double));
      MPI_Irecv(&lb->hl.recvXnp[0], nsendX, MPI_DOUBLE, Xnp, tagpn, comm, &lb->hl.recvreq[7]);
      count = prepareBuffersX(lb, fptr, nlocal, id);
      MPI_Issend(&lb->hl .sendXpn[0], nsendX, MPI_DOUBLE, Xpn, tagpn, comm, &lb->hl.sendreq[7]);
      break;
    case 8:
      lb->hl .recvXpn = (double *) malloc(nsendX*sizeof(double));
      lb->hl .sendXnp = (double *) malloc(nsendX*sizeof(double));
      MPI_Irecv(&lb->hl.recvXpn[0], nsendX, MPI_DOUBLE, Xpn, tagnp, comm, &lb->hl.recvreq[8]);
      count = prepareBuffersX(lb, fptr, nlocal,id);
      MPI_Issend(&lb->hl .sendXnp[0], nsendX, MPI_DOUBLE, Xnp, tagnp, comm, &lb->hl.sendreq[8]);
      break;
    case 9:
      lb->hl .sendXnn = (double *) malloc(nsendX*sizeof(double));
      lb->hl .recvXpp = (double *) malloc(nsendX*sizeof(double));
      MPI_Irecv(&lb->hl.recvXpp[0], nsendX, MPI_DOUBLE, Xpp, tagnn, comm, &lb->hl.recvreq[9]);
      count = prepareBuffersX(lb, fptr, nlocal, id);
      MPI_Issend(&lb->hl .sendXnn[0], nsendX, MPI_DOUBLE, Xnn, tagnn, comm, &lb->hl.sendreq[9]);
      break;
    default:
      info("Error in Edge X\n");
      exit(1);
    }

    assert(count == nsendX);
  }


  if ( id >= 10 && id < 14){

    /* Receive edges parallel to y-direction */
    Ynn = nonblocking_cart_neighb(NMN);
    Ynp = nonblocking_cart_neighb(NMP);
    Ypn = nonblocking_cart_neighb(PMN);
    Ypp = nonblocking_cart_neighb(PMP);

    switch(id){
    case 10:
      lb->hl .recvYnn = (double *) malloc(nsendY*sizeof(double));
      lb->hl .sendYpp = (double *) malloc(nsendY*sizeof(double));
      MPI_Irecv(&lb->hl.recvYnn[0], nsendY, MPI_DOUBLE, Ynn, tagpp, comm, &lb->hl.recvreq[10]);
      count=prepareBuffersY(lb, fptr, nlocal, id);
      MPI_Issend(&lb->hl .sendYpp[0], nsendY, MPI_DOUBLE, Ypp, tagpp, comm, &lb->hl.sendreq[10]);
      break;

    case 11:
      lb->hl .recvYnp = (double *) malloc(nsendY*sizeof(double));
      lb->hl .sendYpn = (double *) malloc(nsendY*sizeof(double));
      MPI_Irecv(&lb->hl.recvYnp[0], nsendY, MPI_DOUBLE, Ynp, tagpn, comm, &lb->hl.recvreq[11]);
      count=prepareBuffersY(lb, fptr, nlocal, id);
      MPI_Issend(&lb->hl .sendYpn[0], nsendY, MPI_DOUBLE, Ypn, tagpn, comm, &lb->hl.sendreq[11]);
      break;

    case 12:
      lb->hl .recvYpn = (double *) malloc(nsendY*sizeof(double));
      lb->hl .sendYnp = (double *) malloc(nsendY*sizeof(double));
      MPI_Irecv(&lb->hl.recvYpn[0], nsendY, MPI_DOUBLE, Ypn, tagnp, comm, &lb->hl.recvreq[12]);
      count=prepareBuffersY(lb, fptr, nlocal, id);
      MPI_Issend(&lb->hl .sendYnp[0], nsendY, MPI_DOUBLE, Ynp, tagnp, comm, &lb->hl.sendreq[12]);
      break;

    case 13:
      lb->hl .recvYpp = (double *) malloc(nsendY*sizeof(double));
      lb->hl .sendYnn = (double *) malloc(nsendY*sizeof(double));
      MPI_Irecv(&lb->hl.recvYpp[0], nsendY, MPI_DOUBLE, Ypp, tagnn, comm, &lb->hl.recvreq[13]);
      count=prepareBuffersY(lb, fptr, nlocal, id);
      MPI_Issend(&lb->hl .sendYnn[0], nsendY, MPI_DOUBLE, Ynn, tagnn, comm, &lb->hl.sendreq[13]);
      break;

    default:
      info("ERROR in Edge Y\n");
      exit(1);
    }

    /* Send edges parallel to y-direction */
    assert(count == nsendY);
  }


  if ( id >= 14 && id < 18 ){

    /* Receive edges parallel to z-direction */
    Znn = nonblocking_cart_neighb(NNM);
    Znp = nonblocking_cart_neighb(NPM);
    Zpn = nonblocking_cart_neighb(PNM);
    Zpp = nonblocking_cart_neighb(PPM);

    switch(id){

    case 14:
      lb->hl .recvZnn = (double *) malloc(nsendZ*sizeof(double));
      lb->hl .sendZpp = (double *) malloc(nsendZ*sizeof(double));
      MPI_Irecv(&lb->hl.recvZnn[0], nsendZ, MPI_DOUBLE, Znn, tagpp, comm, &lb->hl.recvreq[14]);
      count=prepareBuffersZ(lb, fptr, nlocal, id);
      MPI_Issend(&lb->hl .sendZpp[0] , nsendZ, MPI_DOUBLE, Zpp, tagpp, comm, &lb->hl.sendreq[14]);
      break;
    case 15:
      lb->hl .recvZnp = (double *) malloc(nsendZ*sizeof(double));
      lb->hl .sendZpn = (double *) malloc(nsendZ*sizeof(double));
      MPI_Irecv(&lb->hl.recvZnp[0], nsendZ, MPI_DOUBLE, Znp, tagpn, comm, &lb->hl.recvreq[15]);
      count=prepareBuffersZ(lb, fptr, nlocal, id);
      MPI_Issend(&lb->hl .sendZpn[0], nsendZ, MPI_DOUBLE, Zpn, tagpn, comm, &lb->hl.sendreq[15]);
      break;
    case 16:
      lb->hl .recvZpn = (double *) malloc(nsendZ*sizeof(double));
      lb->hl .sendZnp = (double *) malloc(nsendZ*sizeof(double));
      MPI_Irecv(&lb->hl.recvZpn[0], nsendZ, MPI_DOUBLE, Zpn, tagnp, comm, &lb->hl.recvreq[16]);
      count=prepareBuffersZ(lb, fptr, nlocal, id);
      MPI_Issend(&lb->hl .sendZnp[0], nsendZ, MPI_DOUBLE, Znp, tagnp, comm, &lb->hl.sendreq[16]);
      break;
    case 17:
      lb->hl .recvZpp = (double *) malloc(nsendZ*sizeof(double));
      lb->hl .sendZnn = (double *) malloc(nsendZ*sizeof(double));
      MPI_Irecv(&lb->hl.recvZpp[0], nsendZ, MPI_DOUBLE, Zpp, tagnn, comm, &lb->hl.recvreq[17]);
      count=prepareBuffersZ(lb, fptr, nlocal, id);
      MPI_Issend(&lb->hl .sendZnn[0], nsendZ, MPI_DOUBLE, Znn, tagnn, comm, &lb->hl.sendreq[17]);
      break;

    default:
      info("Error in Z direction\n");
      exit(1);
    }

    assert(count == nsendZ);

  }



  return;
}


/*****************************************************************************
 *
 *  halo_corners
 *
 *  Sends 8 MPI messages (the corners) for the Non-Blocking version.
 *
 *****************************************************************************/

void halo_corners_tasks_i(lb_t * lb, int id) {

  int ic, jc, kc;
  int n, p;
  int count;
  int index, indexhalo, indexreal;
  int nlocal[3];

  const int tagnnn = 907;
  const int tagnnp = 908;
  const int tagnpn = 909;
  const int tagnpp = 910;
  const int tagpnn = 911;
  const int tagpnp = 912;
  const int tagppn = 913;
  const int tagppp = 914;

  /* Ranks of neighbouring corners XYZ direction*/
  int ppp, ppn, pnp, pnn, npp, npn, nnp, nnn;

  MPI_Comm comm = cart_comm();

  assert(lb);

  coords_nlocal(nlocal);

  double* fptr;
  if (get_step())
    fptr=lb->tcopy->f;
  else
    fptr=lb->f;

  int nsend;
  nsend = NVEL*lb->ndist*1;

  /* Allocate message sizes for plane send/receives */

  nnn = nonblocking_cart_neighb(NNN);
  nnp = nonblocking_cart_neighb(NNP);
  npn = nonblocking_cart_neighb(NPN);
  npp = nonblocking_cart_neighb(NPP);
  pnn = nonblocking_cart_neighb(PNN);
  pnp = nonblocking_cart_neighb(PNP);
  ppn = nonblocking_cart_neighb(PPN);
  ppp = nonblocking_cart_neighb(PPP);

  switch(id) {

  case 18:
    lb->hl .recvnnn = (double *) malloc(nsend*sizeof(double));
    lb->hl .sendppp = (double *) malloc(nsend*sizeof(double));
    MPI_Irecv(&lb->hl.recvnnn[0], nsend, MPI_DOUBLE, nnn, tagppp, comm, &lb->hl.recvreq[18]);
    count =  edgeBuffers(lb, fptr, nlocal, id);
    MPI_Issend(&lb->hl .sendppp[0], nsend, MPI_DOUBLE, ppp, tagppp, comm, &lb->hl.sendreq[18]);

    break;
  case 19:
    lb->hl .recvnnp = (double *) malloc(nsend*sizeof(double));
    lb->hl .sendppn = (double *) malloc(nsend*sizeof(double));

    MPI_Irecv(&lb->hl.recvnnp[0], nsend, MPI_DOUBLE, nnp, tagppn, comm, &lb->hl.recvreq[19]);
    count =  edgeBuffers(lb, fptr, nlocal, id);
    MPI_Issend(&lb->hl .sendppn[0], nsend, MPI_DOUBLE, ppn, tagppn, comm, &lb->hl.sendreq[19]);

    break;
  case 20:
    lb->hl .recvnpn = (double *) malloc(nsend*sizeof(double));
    lb->hl .sendpnp = (double *) malloc(nsend*sizeof(double));

    MPI_Irecv(&lb->hl.recvnpn[0], nsend, MPI_DOUBLE, npn, tagpnp, comm, &lb->hl.recvreq[20]);
    count =  edgeBuffers(lb, fptr, nlocal, id);
    MPI_Issend(&lb->hl .sendpnp[0], nsend, MPI_DOUBLE, pnp, tagpnp, comm, &lb->hl.sendreq[20]);


    break;
  case 21:
    lb->hl .recvnpp = (double *) malloc(nsend*sizeof(double));
    lb->hl .sendpnn = (double *) malloc(nsend*sizeof(double));

    MPI_Irecv(&lb->hl.recvnpp[0], nsend, MPI_DOUBLE, npp, tagpnn, comm, &lb->hl.recvreq[21]);
    count =  edgeBuffers(lb, fptr, nlocal, id);
    MPI_Issend(&lb->hl .sendpnn[0], nsend, MPI_DOUBLE, pnn, tagpnn, comm, &lb->hl.sendreq[21]);

    break;

  case 22:
    lb->hl .recvpnn = (double *) malloc(nsend*sizeof(double));
    lb->hl .sendnpp = (double *) malloc(nsend*sizeof(double));

    MPI_Irecv(&lb->hl.recvpnn[0], nsend, MPI_DOUBLE, pnn, tagnpp, comm, &lb->hl.recvreq[22]);
    count =  edgeBuffers(lb, fptr, nlocal, id);
    MPI_Issend(&lb->hl .sendnpp[0], nsend, MPI_DOUBLE, npp, tagnpp, comm, &lb->hl.sendreq[22]);

    break;

  case 23:
    lb->hl .recvpnp = (double *) malloc(nsend*sizeof(double));
    lb->hl .sendnpn = (double *) malloc(nsend*sizeof(double));

    MPI_Irecv(&lb->hl.recvpnp[0], nsend, MPI_DOUBLE, pnp, tagnpn, comm, &lb->hl.recvreq[23]);
    count =  edgeBuffers(lb, fptr, nlocal, id);
    MPI_Issend(&lb->hl .sendnpn[0], nsend, MPI_DOUBLE, npn, tagnpn, comm, &lb->hl.sendreq[23]);

    break;
  case 24:
    lb->hl .recvppn = (double *) malloc(nsend*sizeof(double));
    lb->hl .sendnnp = (double *) malloc(nsend*sizeof(double));

    MPI_Irecv(&lb->hl.recvppn[0], nsend, MPI_DOUBLE, ppn, tagnnp, comm, &lb->hl.recvreq[24]);
    count =  edgeBuffers(lb, fptr, nlocal, id);
    MPI_Issend(&lb->hl .sendnnp[0], nsend, MPI_DOUBLE, nnp, tagnnp, comm, &lb->hl.sendreq[24]);
    break;

  case 25:
    lb->hl .recvppp = (double *) malloc(nsend*sizeof(double));
    lb->hl .sendnnn = (double *) malloc(nsend*sizeof(double));
    MPI_Irecv(&lb->hl.recvppp[0], nsend, MPI_DOUBLE, ppp, tagnnn, comm, &lb->hl.recvreq[25]);
    count =  edgeBuffers(lb, fptr, nlocal, id);
    MPI_Issend(&lb->hl .sendnnn[0], nsend, MPI_DOUBLE, nnn, tagnnn, comm, &lb->hl.sendreq[25]);
    break;

  default:
    exit(1);
  }

  return;
}

