
/*********************************************************************/
/*
/*    Halo exchange via copy to a tmp buffer
/*
/*********************************************************************/

int coords_field_halo_planes(int nhcomm, int na, void * mbuf,
			    MPI_Datatype mpidata) {
  int sz;
  int ic, jc, kc;
  int ia, index;
  int nh;
  int ireal, ihalo;
  int icount, nsend;
  int pforw, pback;
  int nlocal[3];

  unsigned char * buf;
  unsigned char * sendforw;
  unsigned char * sendback;
  unsigned char * recvforw;
  unsigned char * recvback;

  MPI_Comm comm;
  MPI_Request req[4];
  MPI_Status status[2];

  const int tagf = 2015;
  const int tagb = 2016;

  assert(mbuf);
  assert(mpidata == MPI_CHAR || mpidata == MPI_DOUBLE);

  buf = (unsigned char *) mbuf;

  comm = cart_comm();
  if (mpidata == MPI_CHAR) sz = sizeof(char);
  if (mpidata == MPI_DOUBLE) sz = sizeof(double);

  coords_nlocal(nlocal);

  /* X-direction */

  nsend = nhcomm*na*nlocal[Y]*nlocal[Z];
  sendforw = (unsigned char *) malloc(nsend*sz);
  sendback = (unsigned char *) malloc(nsend*sz);
  recvforw = (unsigned char *) malloc(nsend*sz);
  recvback = (unsigned char *) malloc(nsend*sz);
  if (sendforw == NULL) fatal("malloc(sendforw) failed\n");
  if (sendback == NULL) fatal("malloc(sendback) failed\n");
  if (recvforw == NULL) fatal("malloc(recvforw) failed\n");
  if (recvback == NULL) fatal("malloc(recvback) failed\n");

  /* Load send buffers */

  icount = 0;

  for (nh = 0; nh < nhcomm; nh++) {
    for (jc = 1; jc <= nlocal[Y]; jc++) {
      for (kc = 1; kc <= nlocal[Z]; kc++) {
	for (ia = 0; ia < na; ia++) {
	  /* Backward going... */
	  index = coords_index(1 + nh, jc, kc);
	  ireal = na*index + ia;
	  memcpy(sendback + icount*sz, buf + ireal*sz, sz);
	  /* ...and forward going. */
	  index = coords_index(nlocal[X] - nh, jc, kc);
	  ireal = na*index +ia;
	  memcpy(sendforw + icount*sz, buf + ireal*sz, sz);
	  icount += 1;
	}
      }
    }
  }

  assert(icount == nsend);

  if (cart_size(X) == 1) {
    memcpy(recvback, sendforw, nsend*sz);
    memcpy(recvforw, sendback, nsend*sz);
    req[2] = MPI_REQUEST_NULL;
    req[3] = MPI_REQUEST_NULL;
  }
  else {
     /* Receive planes in the X-direction */
    /* PPM, NMM, P=Positive, M=Middle, N=Negative for the XYZ directions respectively */
    int pforwX = nonblocking_cart_neighb(PMM);
    int pbackX = nonblocking_cart_neighb(NMM);
    
    pforw = cart_neighb(FORWARD, X);
    pback = cart_neighb(BACKWARD, X);
    assert(pforwX == pforw);
    assert(pbackX == pback);
    MPI_Irecv(recvforw, nsend, mpidata, pforw, tagb, comm, req);
    MPI_Irecv(recvback, nsend, mpidata, pback, tagf, comm, req + 1);
    MPI_Issend(sendback, nsend, mpidata, pback, tagb, comm, req + 2);
    MPI_Issend(sendforw, nsend, mpidata, pforw, tagf, comm, req + 3);
    /* Wait for receives */
    MPI_Waitall(2, req, status);
  }

  /* Unload */

  icount = 0;

  for (nh = 0; nh < nhcomm; nh++) {
    for (jc = 1; jc <= nlocal[Y]; jc++) {
      for (kc = 1; kc <= nlocal[Z]; kc++) {
	for (ia = 0; ia < na; ia++) {
	  index = coords_index(nlocal[X] + 1 + nh, jc, kc);
	  ihalo = na*index +ia;
	  memcpy(buf + ihalo*sz, recvforw + icount*sz, sz);
	  index = coords_index(0 - nh, jc, kc);
	  ihalo = na*index + ia;
	  memcpy(buf + ihalo*sz, recvback + icount*sz, sz);
	  icount += 1;
	}
      }
    }
  }

  assert(icount == nsend);

  free(recvback);
  free(recvforw);

  MPI_Waitall(2, req + 2, status);

  free(sendback);
  free(sendforw);

  /* Y direction */

  nsend = nhcomm*na*(nlocal[X] + 2*nhcomm)*nlocal[Z];
  sendforw = (unsigned char *) malloc(nsend*sz);
  sendback = (unsigned char *) malloc(nsend*sz);
  recvforw = (unsigned char *) malloc(nsend*sz);
  recvback = (unsigned char *) malloc(nsend*sz);
  if (sendforw == NULL) fatal("malloc(sendforw) failed\n");
  if (sendback == NULL) fatal("malloc(sendback) failed\n");
  if (recvforw == NULL) fatal("malloc(recvforw) failed\n");
  if (recvback == NULL) fatal("malloc(recvback) failed\n");

  /* Load buffers */

  icount = 0;
  
  for (nh = 0; nh < nhcomm; nh++) {
    for (ic = 1 ; ic <= nlocal[X]; ic++) {
      for (kc = 1; kc <= nlocal[Z]; kc++) {
	for (ia = 0; ia < na; ia++) {
	  index = coords_index(ic, 1 + nh, kc);
	  ireal = na*index + ia;
	  memcpy(sendback + icount*sz, buf + ireal*sz, sz);
	  index = coords_index(ic, nlocal[Y] - nh, kc);
	  ireal = na*index + ia;
	  memcpy(sendforw + icount*sz, buf + ireal*sz, sz);
	  icount += 1;
	}
      }
    }
  }

  assert(icount == nsend);

  if (cart_size(Y) == 1) {
    memcpy(recvback, sendforw, nsend*sz);
    memcpy(recvforw, sendback, nsend*sz);
    req[2] = MPI_REQUEST_NULL;
    req[3] = MPI_REQUEST_NULL;
  }
  else {
    pforw = cart_neighb(FORWARD, Y);
    pback = cart_neighb(BACKWARD, Y);
    MPI_Irecv(recvforw, nsend, mpidata, pforw, tagb, comm, req);
    MPI_Irecv(recvback, nsend, mpidata, pback, tagf, comm, req + 1);
    MPI_Issend(sendback, nsend, mpidata, pback, tagb, comm, req + 2);
    MPI_Issend(sendforw, nsend, mpidata, pforw, tagf, comm, req + 3);
    /* Wait for receives */
    MPI_Waitall(2, req, status);
  }

  /* Unload */

  icount = 0;

  for (nh = 0; nh < nhcomm; nh++) {
    for (ic = 1; ic <= nlocal[X]; ic++) {
      for (kc = 1; kc <= nlocal[Z]; kc++) {
	for (ia = 0; ia < na; ia++) {
	  index = coords_index(ic, 0 - nh, kc);
	  ihalo = na*index + ia;
	  memcpy(buf + ihalo*sz, recvback + icount*sz, sz);
	  index = coords_index(ic, nlocal[Y] + 1 + nh, kc);
	  ihalo = na*index + ia;
	  memcpy(buf + ihalo*sz, recvforw + icount*sz, sz);
	  icount += 1;
	}
      }
    }
  }

  assert(icount == nsend);

  free(recvback);
  free(recvforw);

  MPI_Waitall(2, req + 2, status);

  free(sendback);
  free(sendforw);

  /* Z direction */

  nsend = nhcomm*na*(nlocal[X] + 2*nhcomm)*(nlocal[Y] + 2*nhcomm);
  sendforw = (unsigned char *) malloc(nsend*sz);
  sendback = (unsigned char *) malloc(nsend*sz);
  recvforw = (unsigned char *) malloc(nsend*sz);
  recvback = (unsigned char *) malloc(nsend*sz);
  if (sendforw == NULL) fatal("malloc(sendforw) failed\n");
  if (sendback == NULL) fatal("malloc(sendback) failed\n");
  if (recvforw == NULL) fatal("malloc(recvforw) failed\n");
  if (recvback == NULL) fatal("malloc(recvback) failed\n");

  /* Load */
  /* Some adjustment in the load required for 2d systems (X-Y) */

  icount = 0;

  for (nh = 0; nh < nhcomm; nh++) {
    for (ic = 1; ic <= nlocal[X]; ic++) {
      for (jc = 1; jc <= nlocal[Y]; jc++) {
	for (ia = 0; ia < na; ia++) {
	  //kc = imin(1 + nh, nlocal[Z]);
	  index = coords_index(ic, jc, 1+nh);
	  ireal = na*index + ia;
	  memcpy(sendback + icount*sz, buf + ireal*sz, sz);
	  //kc = imax(nlocal[Z] - nh, 1);
	  index = coords_index(ic, jc, nlocal[Z] - nh);
	  ireal = na*index + ia;
	  memcpy(sendforw + icount*sz, buf + ireal*sz, sz);
	  icount += 1;
	}
      }
    }
  }

  assert(icount == nsend);

  if (cart_size(Z) == 1) {
    memcpy(recvback, sendforw, nsend*sz);
    memcpy(recvforw, sendback, nsend*sz);
    req[2] = MPI_REQUEST_NULL;
    req[3] = MPI_REQUEST_NULL;
  }
  else {
    pforw = cart_neighb(FORWARD, Z);
    pback = cart_neighb(BACKWARD, Z);
    MPI_Irecv(recvforw, nsend, mpidata, pforw, tagb, comm, req);
    MPI_Irecv(recvback, nsend, mpidata, pback, tagf, comm, req + 1);
    MPI_Issend(sendback, nsend, mpidata, pback, tagb, comm, req + 2);
    MPI_Issend(sendforw, nsend, mpidata, pforw, tagf, comm, req + 3);
    /* Wait before unloading */
    MPI_Waitall(2, req, status);
  }

  /* Unload */

  icount = 0;

  for (nh = 0; nh < nhcomm; nh++) {
    for (ic = 1; ic <= nlocal[X]; ic++) {
      for (jc = 1; jc <= nlocal[Y]; jc++) {
	for (ia = 0; ia < na; ia++) {
	  index = coords_index(ic, jc, 0 - nh);
	  ihalo = na * index +ia;
	  memcpy(buf + ihalo*sz, recvback + icount*sz, sz);
	  index = coords_index(ic, jc, nlocal[Z] + 1 + nh);
	  ihalo = na*index + ia;
	  memcpy(buf + ihalo*sz, recvforw + icount*sz, sz);
	  icount += 1;
	}
      }
    }
  }

  assert(icount == nsend);

  free(recvback);
  free(recvforw);

  MPI_Waitall(2, req + 2, status);

  free(sendback);
  free(sendforw);

  return 0;
}

