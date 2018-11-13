#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "pe.h"
#include "coords.h"
#include "coords_field.h"

/*********************************************************************/
/*
  /*    Halo exchange via copy to a tmp buffer
  /*
  /*********************************************************************/

int coords_field_halo_edges(int nhcomm, int na, void * mbuf, MPI_Datatype mpidata) {
  
  int sz;
  int ic, jc, kc;
  int ia, index;
  int nh;
  int ireal, ihalo;
  int icount, nsend;
  int pforw, pback;
  int nlocal[3];

  unsigned char * buf;
  unsigned char * sendXnn;
  unsigned char * sendXnp;
  unsigned char * recvXnn;
  unsigned char * recvXnp;
  unsigned char * sendXpn;
  unsigned char * sendXpp;
  unsigned char * recvXpn;
  unsigned char * recvXpp;

  MPI_Comm comm;
  MPI_Request sendreq[12],recvreq[12];
  MPI_Status status[4];

  const int tagnn = 903;
  const int tagnp = 904;
  const int tagpn = 905;
  const int tagpp = 906;

  assert(mbuf);
  assert(mpidata == MPI_CHAR || mpidata == MPI_DOUBLE);

  buf = (unsigned char *) mbuf;

  /* Ranks of neighbouring edges.
     Xpp, X is the direction, parallel to which edges are
     sent pp refers to the YZ directions respectively*/
  int Xpp, Xpn, Xnp, Xnn, Ypp, Ypn, Ynp, Ynn, Zpp, Zpn, Znp, Znn;

  comm = cart_comm();

  if (mpidata == MPI_CHAR) sz = sizeof(char);
  if (mpidata == MPI_DOUBLE) sz = sizeof(double);

  coords_nlocal(nlocal);

  /* X-direction */

  nsend = nhcomm*na*nlocal[X];
  sendXnn = (unsigned char *) malloc(nsend*sz);
  sendXnp = (unsigned char *) malloc(nsend*sz);
  sendXpn = (unsigned char *) malloc(nsend*sz);
  sendXpp = (unsigned char *) malloc(nsend*sz);

  recvXnn = (unsigned char *) malloc(nsend*sz);
  recvXnp = (unsigned char *) malloc(nsend*sz);
  recvXpn = (unsigned char *) malloc(nsend*sz);
  recvXpp = (unsigned char *) malloc(nsend*sz);

  /* Receive edges parallel to x-direction*/
  Xnn = nonblocking_cart_neighb(MNN);
  Xnp = nonblocking_cart_neighb(MNP);
  Xpn = nonblocking_cart_neighb(MPN);
  Xpp = nonblocking_cart_neighb(MPP);

  /* Load send buffers */

  icount = 0;

  for (nh = 0; nh < nhcomm; nh++) {
    for (ic = 1; ic <= nlocal[X]; ic++) {
      for (ia = 0; ia < na; ia++) {

        index = coords_index(ic, 1, 1);
        ireal = na*index + ia;
        memcpy(sendXnn + icount*sz, buf + ireal*sz, sz);

        index = coords_index(ic, 1, nlocal[Z]);
        ireal = na*index + ia;
        memcpy(sendXnp + icount*sz, buf + ireal*sz, sz);

        index = coords_index(ic, nlocal[Y], 1);
        ireal = na*index + ia;
        memcpy(sendXpn + icount*sz, buf + ireal*sz, sz);


        index = coords_index(ic, nlocal[Y], nlocal[Z]);
        ireal = na*index + ia;
        memcpy(sendXpp + icount*sz, buf + ireal*sz, sz);

        count+=1;

      }
    }
  }

  assert(icount == nsend);

  if (cart_size(X) == 1) {
    //memcpy(recvback, sendforw, nsend*sz);
    //memcpy(recvforw, sendback, nsend*sz);
    //req[2] = MPI_REQUEST_NULL;
    //req[3] = MPI_REQUEST_NULL;
  }
  else {

    MPI_Irecv(recvXnn, nsend, mpidata, Xnn, tagpp, comm, recvreq + 0);
    MPI_Irecv(recvXnp, nsend, mpidata, Xnp, tagpn, comm, recvreq + 1);
    MPI_Irecv(recvXpn, nsend, mpidata, Xpn, tagnp, comm, recvreq + 2);
    MPI_Irecv(recvXpp, nsend, mpidata, Xpp, tagnn, comm, recvreq + 3);

    MPI_Issend(sendXpp, nsend, mpidata, Xpp, tagpp, comm, sendreq + 0);
    MPI_Issend(sendXpn, nsend, mpidata, Xpn, tagpn, comm, sendreq + 1);
    MPI_Issend(sendXnp, nsend, mpidata, Xnp, tagnp, comm, sendreq + 2);
    MPI_Issend(sendXnn, nsend, mpidata, Xnn, tagnn, comm, sendreq + 3);



    /* Wait for receives */
    MPI_Waitall(4, recvreq, status);
  }

  /* Unload */

  icount = 0;
  for (nh = 0; nh < nhcomm; nh++) {
    for (ic = 1; ic <= nlocal[X]; ic++) {
      for (ia = 0; ia < na; ia++) {

        index = coords_index(ic, nlocal[Y] + 1, nlocal[Z] + 1);
        ihalo = na*index +ia;
        memcpy(buf + ihalo*sz, recvXpp + icount*sz, sz);

        index = coords_index(ic, nlocal[Y] + 1, 0);
        ihalo = na*index +ia;
        memcpy(buf + ihalo*sz, recvXpp + icount*sz, sz);

        index = coords_index(ic, 0, nlocal[Z] + 1);
        ihalo = na*index + ia;
        memcpy(buf + ihalo*sz, recvXnp + icount*sz, sz);

        index = coords_index(ic, 0, 0);
        ihalo = na*index + ia;
        memcpy(buf + ihalo*sz, recvXnn + icount*sz, sz);
        icount += 1;

      }
    }
  }

  assert(icount == nsend);

  free(recvXpp);
  free(recvXnn);
  free(recvXnp);
  free(recvXpn);

  MPI_Waitall(4, sendreq , status);

  free(sendXpp);
  free(sendXnn);
  free(sendXnp);
  free(sendXpn);

  /* Y direction */

  /* Receive edges parallel to y-direction*/
  Ynn = nonblocking_cart_neighb(NMN);
  Ynp = nonblocking_cart_neighb(NMP);
  Ypn = nonblocking_cart_neighb(PMN);
  Ypp = nonblocking_cart_neighb(PMP);
  
  nsend = nhcomm*na*nlocal[Y];

  sendYnn = (unsigned char *) malloc(nsend*sz);
  sendYnp = (unsigned char *) malloc(nsend*sz);
  sendYpn = (unsigned char *) malloc(nsend*sz);
  sendYpp = (unsigned char *) malloc(nsend*sz);

  recvYnn = (unsigned char *) malloc(nsend*sz);
  recvYnp = (unsigned char *) malloc(nsend*sz);
  recvYpn = (unsigned char *) malloc(nsend*sz);
  recvYpp = (unsigned char *) malloc(nsend*sz);

  /* if (sendforw == NULL) fatal("malloc(sendforw) failed\n"); */
  /* if (sendback == NULL) fatal("malloc(sendback) failed\n"); */
  /* if (recvforw == NULL) fatal("malloc(recvforw) failed\n"); */
  /* if (recvback == NULL) fatal("malloc(recvback) failed\n"); */

  /* Load buffers */

  icount = 0;

  for (nh = 0; nh < nhcomm; nh++) {
    for (jc = 1 ; jc <= nlocal[Y]; ic++) {
      for (ia = 0; ia < na; ia++) {

        index = coords_index(1, jc, 1);
        ireal = na*index + ia;
        memcpy(sendYnn + icount*sz, buf + ireal*sz, sz);

        index = coords_index(1, jc, nlocal[Z]);
        ireal = na*index + ia;
        memcpy(sendYnp + icount*sz, buf + ireal*sz, sz);

        index = coords_index(nlocal[X], jc, 1);
        ireal = na*index + ia;
        memcpy(sendYpn + icount*sz, buf + ireal*sz, sz);

        index = coords_index(nlocal[X], jc, nlocal[Z]);
        ireal = na*index + ia;
        memcpy(sendYpp + icount*sz, buf + ireal*sz, sz);

        icount += 1;

      }
    }
  }

  assert(icount == nsend);

  if (cart_size(Y) == 1) {
    //memcpy(recvback, sendforw, nsend*sz);
    //memcpy(recvforw, sendback, nsend*sz);
    //req[6] = MPI_REQUEST_NULL;
    //req[7] = MPI_REQUEST_NULL;

  }
  else {

    MPI_Irecv(recvYnn, nsend, mpidata, Ynn, tagpp, comm, recvreq + 4);
    MPI_Irecv(recvYnp, nsend, mpidata, Ynp, tagpn, comm, recvreq + 5);
    MPI_Irecv(recvYpn, nsend, mpidata, Ypn, tagnp, comm, recvreq + 6);
    MPI_Irecv(recvYpp, nsend, mpidata, Ypp, tagnn, comm, recvreq + 7);

    MPI_Issend(sendYpp, nsend, mpidata, Ypp, tagpp, comm, sendreq + 4);
    MPI_Issend(sendYpn, nsend, mpidata, Ypn, tagpn, comm, sendreq + 5);
    MPI_Issend(sendYnp, nsend, mpidata, Ynp, tagnp, comm, sendreq + 6);
    MPI_Issend(sendYnn, nsend, mpidata, Ynn, tagnn, comm, sendreq + 7);

    /* Wait for receives */
    MPI_Waitall(4, recvreq+4, status);
  }

  /* Unload */

  icount = 0;

  for (nh = 0; nh < nhcomm; nh++) {
    for (jc = 1; jc <= nlocal[Y]; jc++) {
      for (ia = 0; ia < na; ia++) {

        index = coords_index(nlocal[X] + nhcomm, jc + nh, nlocal[Z] + nhcomm);
        ihalo = na*index + ia;
        memcpy(buf + ihalo*sz, recvYpp + icount*sz, sz);

        index = coords_index(nlocal[X] + nhcomm, jc + nh, 0);
        ihalo = na*index + ia;
        memcpy(buf + ihalo*sz, recvYpn + icount*sz, sz);

        index = coords_index(0, jc + nh, nlocal[Z] + nhcomm);
        ihalo = na*index + ia;
        memcpy(buf + ihalo*sz, recvYnp + icount*sz, sz);

        index = coords_index(0, jc + nh, 0);
        ihalo = na*index + ia;
        memcpy(buf + ihalo*sz, recvYnn + icount*sz, sz);

        icount += 1;
      }
    }
  }

  assert(icount == nsend);

  free(recvYnn);
  free(recvYpp);
  free(recvYpn);
  free(recvYnp);

  MPI_Waitall(4, sendreq + 4, status);

  free(sendYpp);
  free(sendYnn);
  free(sendYnp);
  free(sendYpn);

  /* Z direction */

  /* Receive edges parallel to z-direction*/
  Znn = nonblocking_cart_neighb(NNM);
  Znp = nonblocking_cart_neighb(NPM);
  Zpn = nonblocking_cart_neighb(PNM);
  Zpp = nonblocking_cart_neighb(PPM);

  nsend = nhcomm*na*nlocal[Z];
  sendZnn = (unsigned char *) malloc(nsend*sz);
  sendZnp = (unsigned char *) malloc(nsend*sz);
  sendZpn = (unsigned char *) malloc(nsend*sz);
  sendZpp = (unsigned char *) malloc(nsend*sz);

  recvZnn = (unsigned char *) malloc(nsend*sz);
  recvZnp = (unsigned char *) malloc(nsend*sz);
  recvZpn = (unsigned char *) malloc(nsend*sz);
  recvZpp = (unsigned char *) malloc(nsend*sz);

  /* if (sendforw == NULL) fatal("malloc(sendforw) failed\n"); */
  /* if (sendback == NULL) fatal("malloc(sendback) failed\n"); */
  /* if (recvforw == NULL) fatal("malloc(recvforw) failed\n"); */
  /* if (recvback == NULL) fatal("malloc(recvback) failed\n"); */

  /* Load */
  /* Some adjustment in the load required for 2d systems (X-Y) */

  icount = 0;

  for (nh = 0; nh < nhcomm; nh++) {
    for (kc = 1; kc <= nlocal[Z]; kc++) {
      for (ia = 0; ia < na; ia++) {

        index = coords_index(1, 1, kc);
        ireal = na*index + ia;
        memcpy(sendZnn + icount*sz, buf + ireal*sz, sz);

        index = coords_index(1, nlocal[Y], kc);
        ireal = na*index + ia;
        memcpy(sendZnp + icount*sz, buf + ireal*sz, sz);

        index = coords_index(nlocal[X], 1, kc);
        ireal = na*index + ia;
        memcpy(sendZpn + icount*sz, buf + ireal*sz, sz);

        index = coords_index(nlocal[X], nlocal[Y], jc);
        ireal = na*index + ia;
        memcpy(sendZpp + icount*sz, buf + ireal*sz, sz);

        icount += 1;
      }
    }
  }

  assert(icount == nsend);

  if (cart_size(Z) == 1) {
    /* memcpy(recvback, sendforw, nsend*sz); */
    /* memcpy(recvforw, sendback, nsend*sz); */
    /* req[2] = MPI_REQUEST_NULL; */
    /* req[3] = MPI_REQUEST_NULL; */
  }
  else {


    MPI_Irecv(recvZnn, nsend, mpidata, Znn, tagpp, comm, recvreq + 8);
    MPI_Irecv(recvZnp, nsend, mpidata, Znp, tagpn, comm, recvreq + 9);
    MPI_Irecv(recvZpn, nsend, mpidata, Zpn, tagnp, comm, recvreq + 10);
    MPI_Irecv(recvZpp, nsend, mpidata, Zpp, tagnn, comm, recvreq + 11);

    MPI_Issend(sendZnn, nsend, mpidata, Zpp, tagpp, comm, sendreq + 8);
    MPI_Issend(sendZpn, nsend, mpidata, Zpn, tagpn, comm, sendreq + 9);
    MPI_Issend(sendZnp, nsend, mpidata, Znp, tagnp, comm, sendreq + 10);
    MPI_Issend(sendZnn, nsend, mpidata, Znn, tagnn, comm, sendreq + 11);

    /* Wait before unloading */
    MPI_Waitall(4, recvreq + 8, status);
  }

  /* Unload */

  icount = 0;

  for (nh = 0; nh < nhcomm; nh++) {
    for (kc = 1; kc <= nlocal[Z]; kc++) {
      for (ia = 0; ia < na; ia++) {

        index = coords_index(nlocal[X]+nhcomm, nlocal[Y] + nhcomm, kc);
        ihalo = na * index +ia;
        memcpy(buf + ihalo*sz, recvZpp + icount*sz, sz);

        index = coords_index(nlocal[X] + nhcomm, 0, kc);
        ihalo = na*index + ia;
        memcpy(buf + ihalo*sz, recvZpn + icount*sz, sz);

        index = coords_index(0, nlocal[Y] + nhcomm , kc);
        ihalo = na * index +ia;
        memcpy(buf + ihalo*sz, recvZnp + icount*sz, sz);

        index = coords_index(0, 0, kc);
        ihalo = na*index + ia;
        memcpy(buf + ihalo*sz, recvZnn + icount*sz, sz);

        icount += 1;

      }
    }
  }

  assert(icount == nsend);

  free(recvZpp);
  free(recvZnn);
  free(recvZnp);
  free(recvZpn);

  MPI_Waitall(4, sendreq + 8, status);

  free(sendYpp);
  free(sendYnn);
  free(sendYnp);
  free(sendYpn);

  return 0;
}
