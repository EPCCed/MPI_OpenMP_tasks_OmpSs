/*
 * targetDP_C.c: API Implementation for targetDP: C version
 * Alan Gray
 *
 * Copyright 2015 The University of Edinburgh
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License. 
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "targetDP.h"



//The targetMalloc function allocates memory on the target.
void targetMalloc(void **address_of_ptr,size_t size){

  void* ptr;
  ptr = malloc(size);

  if(!ptr){
    printf("malloc failed\n");
    exit(1);
  }
    

  *address_of_ptr=ptr;



  return;
  
}

//The targetCalloc function allocates, and initialises to zero, memory on the target.
void targetCalloc(void **address_of_ptr,size_t size){

  void* ptr;
  ptr = calloc(1,size);

  if(!ptr){
    printf("calloc failed\n");
    exit(1);
  }
    

  *address_of_ptr=ptr;



  return;
  
}


// The targetMalloc function allocates unified memory that can be accessed
// on the host or the target.
void targetMallocUnified(void **address_of_ptr,size_t size){

  void* ptr;
  ptr = malloc(size);

  if(!ptr){
    printf("malloc failed\n");
    exit(1);
  }
    

  *address_of_ptr=ptr;



  return;
  
}

// The targetCalloc function allocates unified memory that can be accessed
// on the host or the target, and is initialised to 0
void targetCallocUnified(void **address_of_ptr,size_t size){

  void* ptr;
  ptr = calloc(1,size);

  if(!ptr){
    printf("calloc failed\n");
    exit(1);
  }
    

  *address_of_ptr=ptr;


  return;
  
}

void targetMallocHost(void **address_of_ptr,size_t size){

  targetMalloc(address_of_ptr,size);

  return;

}

//The targetFree function deallocates memory on the target.
void targetFree(void *ptr){
  
  free(ptr);
  return;
  
}


//The copyToTarget function copies data from the host to the target.
void copyToTarget(void *targetData,const void* data,const size_t size){

  memcpy(targetData,data,size);
  return;

}

//The copyFromTarget function copies data from the target to the host.
void copyFromTarget(void *data,const void* targetData,const size_t size){

  memcpy(data,targetData,size);
  return;

}

// The targetInit3D initialises the environment required to perform any of the
// “3D” operations defined below.
void targetInit3D(int extents[3], size_t nfieldsmax, int nhalo){
  return;
}

//deprecated
void targetInit(int extents[3], size_t nfieldsmax, int nhalo){
  return;
}



// The targetFinalize3D finalises the targetDP 3D environment.
void targetFinalize3D(){
  return;
}

// deprecated
void targetFinalize(){
  return;
}


void copyToTargetMasked(double *targetData,const double* data,const size_t nsites,
			const size_t nfields,char* siteMask){

  int i,j;
  for (i=0;i<nfields;i++){
    for (j=0;j<nsites;j++){
      if(siteMask[j]) targetData[i*nsites+j]=data[i*nsites+j];
    }
  }
  return;
  
}

//
void copyToTargetMaskedAoS(double *targetData,const double* data,const size_t nsites,
			const size_t nfields,char* siteMask){

  int i,j;
  for (j=0;j<nsites;j++){
    if(siteMask[j]){
      for (i=0;i<nfields;i++){
	targetData[j*nfields+i]=data[j*nfields+i];
      }
    }
  }

  return;
  
}

void copyFromTargetMasked(double *data,const double* targetData,const size_t nsites,
			const size_t nfields,char* siteMask){

  int i, j;
  for (i=0;i<nfields;i++){
    for (j=0;j<nsites;j++){
       if(siteMask[j]) data[i*nsites+j]=targetData[i*nsites+j];
    }
  }
  return;

}

//
void copyFromTargetMaskedAoS(double *data,const double* targetData,const size_t nsites,
			const size_t nfields,char* siteMask){

  int i, j;
  for (j=0;j<nsites;j++){
    if(siteMask[j]){
      for (i=0;i<nfields;i++){
	data[j*nfields+i]=targetData[j*nfields+i];
      }
    }
  }
  return;
    
}

// The targetSynchronize function is used to block until 
// the preceding __targetLaunch__ has completed.
void targetSynchronize(){
  return;
}


//
void targetZero(double* array,size_t size){

  int i;

#ifdef _OPENMP
#pragma omp parallel for 
#endif
  for(i=0;i<size;i++){
    
    array[i]=0.;
    
  }

}


void targetSetConstant(double* array,double value,size_t size){

  int i;

#ifdef _OPENMP
#pragma omp parallel for 
#endif
  for(i=0;i<size;i++){
    
    array[i]=value;
    
  }

}


__targetHost__ void targetAoS2SoA(double* array, size_t nsites, size_t nfields)
{
  
  int i,k;
  double* tmpbuf = (double*) malloc(nsites*nfields*sizeof(double));
  double* tmpbuf2 = (double*) malloc(nsites*nfields*sizeof(double));
  copyFromTarget(tmpbuf,array,nsites*nfields*sizeof(double));	  
  for(i=0;i<nsites;i++){
    for(k=0;k<nfields;k++){
      //tmpbuf2[i*nfields+k]=tmpbuf[k*nsites+i];
      tmpbuf2[k*nsites+i]=tmpbuf[i*nfields+k];
    }
  }
  copyToTarget(array,tmpbuf2,nsites*nfields*sizeof(double));	  
  free(tmpbuf);
  free(tmpbuf2);
}

__targetHost__ void targetSoA2AoS(double* array, size_t nsites, size_t nfields)
{
  
  int i,k;
  double* tmpbuf = (double*) malloc(nsites*nfields*sizeof(double));
  double* tmpbuf2 = (double*) malloc(nsites*nfields*sizeof(double));
  copyFromTarget(tmpbuf,array,nsites*nfields*sizeof(double));	  
  for(i=0;i<nsites;i++){
    for(k=0;k<nfields;k++){
      tmpbuf2[i*nfields+k]=tmpbuf[k*nsites+i];
      //tmpbuf2[k*nsites+i]=tmpbuf[i*nfields+k];
    }
  }
  copyToTarget(array,tmpbuf2,nsites*nfields*sizeof(double));	  
  free(tmpbuf);
  free(tmpbuf2);
}

//
void checkTargetError(const char *msg){

  return;

}


__targetHost__ void copyDeepDoubleArrayToTarget(void* targetObjectAddress,void* hostObjectAddress,void* hostComponentAddress,int size)
{

  //calculate the offset between the object start and the component of interest
  long int offset=((char*)hostComponentAddress)-((char*)hostObjectAddress);

  //calculate the position of this component on the target
  char* targetComponentAddress=((char*)targetObjectAddress)+offset;

  //get a copy of this component pointer on the target
  void* tmpptr;
  copyFromTarget(&tmpptr,targetComponentAddress,sizeof(double*)); 

  //copy the deep array 
  void** ptrToHostComponent= (void**) hostComponentAddress;
  copyToTarget(tmpptr,*ptrToHostComponent,size*sizeof(double));


}


__targetHost__ void copyDeepDoubleArrayFromTarget(void* hostObjectAddress,void* targetObjectAddress,void* hostComponentAddress,int size)
{

  //calculate the offset between the object start and the component of interest
  long int offset=((char*)hostComponentAddress)-((char*)hostObjectAddress);

  //calculate the position of this component on the target
  char* targetComponentAddress=((char*)targetObjectAddress)+offset;

  //get a copy of this component pointer on the target
  void* tmpptr;
  copyFromTarget(&tmpptr,targetComponentAddress,sizeof(double*)); 

  //copy the deep array 
  void** ptrToHostComponent= (void**) hostComponentAddress;
  copyFromTarget(*ptrToHostComponent,tmpptr,size*sizeof(double));


}

double targetDoubleSum(double* array, size_t size){

  int i;

  double result=0.;

#ifdef _OPENMP
#pragma omp parallel for reduction(+:result)
#endif
  for (i=0;i<size;i++)
    result+=array[i];

  return result;
}
