/****************************************************************************
 *
 *  advection_bcs.c
 *
 *  Advection boundary conditions.
 *
 *  $Id: advection_bcs.c 2854 2016-02-09 12:48:10Z agray3 $
 *
 *  Edinburgh Soft Matter and Statistical Physics Group and
 *  Edinburgh Parallel Computing Centre
 *
 *  Kevin Stratford (kevin@epcc.ed.ac.uk)
 *  (c) 2009 The University of Edinburgh
 *
 ****************************************************************************/

#include <assert.h>
#include <stdlib.h>

#include "pe.h"
#include "wall.h"
#include "coords.h"
#include "coords_field.h"
#include "advection_s.h"
#include "advection_bcs.h"
#include "psi_gradients.h"
#include "map_s.h"
#include "timer.h"

/*****************************************************************************
 *
 *  advection_bcs_no_normal_fluxes
 *
 *  Set normal fluxes at solid fluid interfaces to zero.
 *
 *****************************************************************************/

__targetEntry__ void advection_bcs_no_normal_flux_lattice(int nf, advflux_t * flux, map_t * map) {


  int baseIndex;
  __targetTLP__(baseIndex,tc_nSites){

    int iv=0;
    int i;

    int n;
    int indexf[VVL];
    
    double mask[VVL], maskw[VVL], maske[VVL], masky[VVL], maskz[VVL];
    
    int coords[3];
    targetCoords3D(coords,tc_Nall,baseIndex);
    
    // if not a halo site:


#if VVL == 1    
    /*restrict operation to the interior lattice sites*/ 
    if (coords[0] >= (tc_nhalo) &&
    	coords[1] >= (tc_nhalo-1) &&
    	coords[2] >= (tc_nhalo-1) &&
    	coords[0] < (tc_Nall[X]-tc_nhalo) &&
	coords[1] < (tc_Nall[Y]-tc_nhalo)  &&
	coords[2] < (tc_Nall[Z]-tc_nhalo)) 
#endif
{


	/* work out which sites in this chunk should be included */
	int includeSite[VVL];
	__targetILP__(iv) includeSite[iv]=0;
	
	int coordschunk[3][VVL];
		
	__targetILP__(iv){
	  for(i=0;i<3;i++){
	    targetCoords3D(coords,tc_Nall,baseIndex+iv);
	    coordschunk[i][iv]=coords[i];
	  }
	}

	__targetILP__(iv){
	  
	  if ((coordschunk[0][iv] >= (tc_nhalo) &&
	       coordschunk[1][iv] >= (tc_nhalo-1) &&
	       coordschunk[2][iv] >= (tc_nhalo-1) &&
	       coordschunk[0][iv] < tc_Nall[X]-(tc_nhalo) &&
	       coordschunk[1][iv] < tc_Nall[Y]-(tc_nhalo)  &&
	       coordschunk[2][iv] < tc_Nall[Z]-(tc_nhalo)))
	    
	    includeSite[iv]=1;
	}

	

      int index1[VVL];
      
      
      __targetILP__(iv) index1[iv]=targetIndex3D(coordschunk[0][iv]-1,coordschunk[1][iv],coordschunk[2][iv],tc_Nall);
      __targetILP__(iv){
	if (includeSite[iv])
	  maskw[iv] = (map->status[index1[iv]]==MAP_FLUID);    
      }
      
      __targetILP__(iv) index1[iv]=targetIndex3D(coordschunk[0][iv]+1,coordschunk[1][iv],coordschunk[2][iv],tc_Nall);
      __targetILP__(iv){
	if (includeSite[iv])
	  maske[iv] = (map->status[index1[iv]]==MAP_FLUID);
      }    
      
      __targetILP__(iv) index1[iv]= targetIndex3D(coordschunk[0][iv],coordschunk[1][iv]+1,coordschunk[2][iv],tc_Nall);
      __targetILP__(iv){
	if (includeSite[iv])
	  masky[iv] = (map->status[index1[iv]]==MAP_FLUID);
      }    
      
      __targetILP__(iv) index1[iv]= targetIndex3D(coordschunk[0][iv],coordschunk[1][iv],coordschunk[2][iv]+1,tc_Nall);
      __targetILP__(iv){
	if (includeSite[iv]) 
	  maskz[iv] = (map->status[index1[iv]]==MAP_FLUID);    
      }      
      
      __targetILP__(iv) index1[iv]= targetIndex3D(coordschunk[0][iv],coordschunk[1][iv],coordschunk[2][iv],tc_Nall);
      __targetILP__(iv){
	if (includeSite[iv])
	  mask[iv] = (map->status[index1[iv]]==MAP_FLUID);  
      }  
      
      for (n = 0;  n < nf; n++) {
	
	__targetILP__(iv) indexf[iv] = ADVADR(tc_nSites,nf,baseIndex+iv,n);


	__targetILP__(iv){
	  if (includeSite[iv])
	    flux->fw[indexf[iv]] *= mask[iv]*maskw[iv];
	}
	__targetILP__(iv){
	  if (includeSite[iv])
	    flux->fe[indexf[iv]] *= mask[iv]*maske[iv];
	}
	__targetILP__(iv){
	  if (includeSite[iv])
	    flux->fy[indexf[iv]] *= mask[iv]*masky[iv];
	}
	__targetILP__(iv){
	  if (includeSite[iv])
	    flux->fz[indexf[iv]] *= mask[iv]*maskz[iv];
	}
      }
 }
  }
  

  return;
}


int advection_bcs_no_normal_flux(int nf, advflux_t * flux, map_t * map) {

  int n;
  int nlocal[3];
  int ic, jc, kc, index, indexf;
  int status;

  assert(nf > 0);
  assert(flux);
  assert(map);

  coords_nlocal(nlocal);

  int nhalo;
  nhalo = coords_nhalo();


  int Nall[3];
  Nall[X]=nlocal[X]+2*nhalo;  Nall[Y]=nlocal[Y]+2*nhalo;  Nall[Z]=nlocal[Z]+2*nhalo;

  int nSites=Nall[X]*Nall[Y]*Nall[Z];


  //copy lattice shape constants to target ahead of execution
  copyConstToTarget(&tc_nSites,&nSites, sizeof(int));
  copyConstToTarget(&tc_nhalo,&nhalo, sizeof(int));
  copyConstToTarget(tc_Nall,Nall, 3*sizeof(int));

  double* tmpptr;

#ifndef KEEPFIELDONTARGET
  map_t* t_map = map->tcopy; //target copy of field structure


  copyFromTarget(&tmpptr,&(t_map->status),sizeof(char*)); 
  copyToTarget(tmpptr,map->status,nSites*sizeof(char));


  advflux_t* t_flux = flux->tcopy; //target copy of flux structure

  copyFromTarget(&tmpptr,&(t_flux->fe),sizeof(double*));
  copyToTarget(tmpptr,flux->fe,nf*nSites*sizeof(double));

  copyFromTarget(&tmpptr,&(t_flux->fw),sizeof(double*));
  copyToTarget(tmpptr,flux->fw,nf*nSites*sizeof(double));

  copyFromTarget(&tmpptr,&(t_flux->fy),sizeof(double*));
  copyToTarget(tmpptr,flux->fy,nf*nSites*sizeof(double));

  copyFromTarget(&tmpptr,&(t_flux->fz),sizeof(double*));
  copyToTarget(tmpptr,flux->fz,nf*nSites*sizeof(double));
#endif

  TIMER_start(ADVECTION_BCS_KERNEL);
  advection_bcs_no_normal_flux_lattice  __targetLaunch__(nSites) (nf,flux->tcopy, map->tcopy);
  targetSynchronize();
  TIMER_stop(ADVECTION_BCS_KERNEL);

#ifndef KEEPFIELDONTARGET
  copyFromTarget(&tmpptr,&(t_flux->fe),sizeof(double*));
  copyFromTarget(flux->fe,tmpptr,nf*nSites*sizeof(double));

  copyFromTarget(&tmpptr,&(t_flux->fw),sizeof(double*));
  copyFromTarget(flux->fw,tmpptr,nf*nSites*sizeof(double));

  copyFromTarget(&tmpptr,&(t_flux->fy),sizeof(double*));
  copyFromTarget(flux->fy,tmpptr,nf*nSites*sizeof(double));

  copyFromTarget(&tmpptr,&(t_flux->fz),sizeof(double*));
  copyFromTarget(flux->fz,tmpptr,nf*nSites*sizeof(double));
#endif

  return 0;
}

/*****************************************************************************
 *
 *  advective_bcs_no_flux
 *
 *  Set normal fluxes at solid fluid interfaces to zero.
 *
 *****************************************************************************/

int advective_bcs_no_flux(int nf, double * fx, double * fy, double * fz,
			  map_t * map) {
  int n;
  int nlocal[3];
  int ic, jc, kc, index, indexf;
  int status;

  double mask, maskx, masky, maskz;

  assert(nf > 0);
  assert(fx);
  assert(fy);
  assert(fz);
  assert(map);

  coords_nlocal(nlocal);

  for (ic = 0; ic <= nlocal[X]; ic++) {
    for (jc = 0; jc <= nlocal[Y]; jc++) {
      for (kc = 0; kc <= nlocal[Z]; kc++) {

	index = coords_index(ic + 1, jc, kc);
	map_status(map, index, &status);
	maskx = (status == MAP_FLUID);

	index = coords_index(ic, jc + 1, kc);
	map_status(map, index, &status);
	masky = (status == MAP_FLUID);

	index = coords_index(ic, jc, kc + 1);
	map_status(map, index, &status);
	maskz = (status == MAP_FLUID);

	index = coords_index(ic, jc, kc);
	map_status(map, index, &status);
	mask = (status == MAP_FLUID);

	for (n = 0;  n < nf; n++) {

	  coords_field_index(index, n, nf, &indexf);
	  fx[indexf] *= mask*maskx;
	  fy[indexf] *= mask*masky;
	  fz[indexf] *= mask*maskz;

	}

      }
    }
  }

  return 0;
}

/*****************************************************************************
 *
 *  advective_bcs_no_flux_d3qx
 *
 *  Set normal fluxes at solid fluid interfaces to zero.
 *
 *****************************************************************************/

int advective_bcs_no_flux_d3qx(int nf, double ** flx, map_t * map) {

  int n;
  int nlocal[3];
  int ic, jc, kc, index0, index1;
  int status;
  int c;
  double * mask;

  assert(nf > 0);
  assert(flx);
  assert(map);

  mask = (double*) calloc(PSI_NGRAD, sizeof(double)); 

  coords_nlocal(nlocal);

  for (ic = 1; ic <= nlocal[X]; ic++) {
    for (jc = 1; jc <= nlocal[Y]; jc++) {
      for (kc = 1; kc <= nlocal[Z]; kc++) {

	index0 = coords_index(ic, jc, kc);
	map_status(map, index0, &status);
	mask[0] = (status == MAP_FLUID);

	for (c = 1; c < PSI_NGRAD; c++) {

	  index1 = coords_index(ic + psi_gr_cv[c][X], jc + psi_gr_cv[c][Y], kc + psi_gr_cv[c][Z]);
	  map_status(map, index1, &status);
	  mask[c] = (status == MAP_FLUID);

	  for (n = 0;  n < nf; n++) {
	    flx[nf*index0 + n][c - 1] *= mask[0]*mask[c];
	  }

	}

      }
    }
  }

  return 0;
}
/*****************************************************************************
 *
 *  advection_bcs_wall
 *
 *  For the case of flat walls, we kludge the order parameter advection
 *  by borrowing the adjacent fluid value.
 *
 *  The official explanation is this may be viewed as a no gradient
 *  condition on the order parameter near the wall.
 *
 *  This allows third and fourth order (x-direction) advective fluxes
 *  to be computed at interface one cell away from wall. Fluxes at
 *  the wall will always be zero.
 *
 ****************************************************************************/

int advection_bcs_wall(field_t * fphi) {

  int ic, jc, kc, index, index1;
  int nlocal[3];
  int nf;
  double q[NQAB];

  if (wall_at_edge(X) == 0) return 0;

  assert(fphi);

  field_nf(fphi, &nf);
  coords_nlocal(nlocal);
  assert(nf <= NQAB);

  if (cart_coords(X) == 0) {
    ic = 1;

    for (jc = 1; jc <= nlocal[Y]; jc++) {
      for (kc = 1; kc <= nlocal[Z]; kc++) {

	index  = coords_index(ic, jc, kc);
	index1 = coords_index(ic-1, jc, kc);

	field_scalar_array(fphi, index, q);
	field_scalar_array_set(fphi, index1, q);
      }
    }
  }

  if (cart_coords(X) == cart_size(X) - 1) {

    ic = nlocal[X];

    for (jc = 1; jc <= nlocal[Y]; jc++) {
      for (kc = 1; kc <= nlocal[Z]; kc++) {

	index = coords_index(ic, jc, kc);
	index1 = coords_index(ic+1, jc, kc);

	field_scalar_array(fphi, index, q);
	field_scalar_array_set(fphi, index1, q);

      }
    }
  }

  return 0;
}
