/* This release was prepared by Dana Akhmetova <danaak@kth.se>/<danieka@gmail.com> on behalf of the INTERTWinE European Exascale Project <http://www.intertwine-project.eu> */
/* iPIC3D was originally developed by Stefano Markidis and Giovanni Lapenta. 
 * This release was contributed by Alec Johnson and Ivy Bo Peng.
 * Publications that use results from iPIC3D need to properly cite  
 * 'S. Markidis, G. Lapenta, and Rizwan-uddin. "Multi-scale simulations of 
 * plasma with iPIC3D." Mathematics and Computers in Simulation 80.7 (2010): 1509-1519.'
 *
 *        Copyright 2015 KTH Royal Institute of Technology
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at 
 *
 *         http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */


#include <mpi.h>
#include "phdf5.h"
#include "ipicdefs.h"
#include "errors.h"
#include "debug.h"
#include "Alloc.h"
#include "MPIdata.h"

#ifdef PHDF5

PHDF5fileClass::PHDF5fileClass(string filestr, int nd, const int *coord, MPI_Comm mpicomm){

  SetDefaultGroups();

  filename    = filestr;
  ndim        = nd;
  comm        = mpicomm;
  for (int i=0; i<ndim; i++)
    mpicoord[i]  = coord[i];

}

PHDF5fileClass::PHDF5fileClass(string filestr){

  SetDefaultGroups();

  filename    = filestr;
}

void PHDF5fileClass::SetDefaultGroups(void){
  grpnames[0] = "Fields";
  grpnames[1] = "Particles";
  grpnames[2] = "Parameters";
}

void PHDF5fileClass::OpenPHDF5file(){

  file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

  ReadPHDF5param();

}

void PHDF5fileClass::CreatePHDF5file(double *L, int *dglob, int *dlocl, bool bp){

  hid_t   acc_t;
  hid_t   Ldataspace;
  hid_t   Ldataset;
  hid_t   ndataspace;
  hid_t   ndataset;
  herr_t  status;
  hsize_t d[1];

  /* ---------------------------------------------- */
  /* 0- Initialize the some of the class parameters */
  /* ---------------------------------------------- */

  bparticles  = bp;
  for (int i=0; i<ndim; i++){
    LxLyLz[i]    = L[i];
    dim[i]       = (hsize_t)dglob[i];
    chdim[i]     = (hsize_t)dlocl[i];
  }

  /* ----------------------------------------------------- */
  /* 1- Set the access template for the parallel HDF5 file */
  /* ----------------------------------------------------- */

  acc_t = H5Pcreate(H5P_FILE_ACCESS);

  /* --------------------------------------- */
  /* 2- Tell HDF5 that we want to use MPI-IO */
  /* --------------------------------------- */

  #ifdef USING_PARALLEL_HDF5
  MPI_Info info;
  MPI_Info_create(&info);
  const int stripe_size = 1024*256;
  const int cb_buffer_size = stripe_size*8;
  const int stripe_count = 16;
  // hint values must be ascii strings,
  // so convert using sprintf
  const int max_chars=32;
  char str_cb_buffer_size[max_chars];
  char str_stripe_size[max_chars];
  char str_stripe_count[max_chars];
  snprintf(str_cb_buffer_size,max_chars,"%d",cb_buffer_size);
  snprintf(str_stripe_size,max_chars,"%d",stripe_size);
  snprintf(str_stripe_count,max_chars,"%d",stripe_count);
  //char* CB_BUFFER_SIZE="1048576"; // 1MB
  //char* STRIPE_SIZE="131072"; // 128k
  //char* STRIPE_COUNT="16"; /* must be an ascii string */
  //char* CB_NODES="8"; /* number of aggregators */
  if(!MPIdata::get_rank()) dprint(str_stripe_count);
  if(!MPIdata::get_rank()) dprint(str_stripe_size);
  if(!MPIdata::get_rank()) dprint(str_cb_buffer_size);
  MPI_Info_set(info, (char*)"striping_factor", str_stripe_count);
  MPI_Info_set(info, (char*)"striping_unit", str_stripe_size);
  MPI_Info_set(info, (char*)"cb_buffer_size", str_cb_buffer_size);
  MPI_Info_set(info, (char*)"romio_cb_write", (char*)"enable");
  //MPI_Info_set(info, "cb_nodes", CB_NODES);
  H5Pset_fapl_mpio(acc_t, comm, info);
  //H5Pset_fapl_mpio(acc_t, comm, MPI_INFO_NULL);

  #else
  eprintf("WriteMethod==Parallel in input file "
          "requires setting USING_PARALLEL_HDF5 in ipicdefs.h");
  #endif

  /* ------------------------------------------------------- */
  /* 3- Load file identifier and release the access template */
  /* ------------------------------------------------------- */

  file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, acc_t);
  H5Pclose(acc_t);

  /* -------------------- */
  /* 4- Set up the groups */
  /* -------------------- */

  for (int i=0; i<ngrp; i++){
    if (!(grpnames[i].c_str()=="Particles"&&!bparticles)){
      hid_t grp;
      grp = H5Gcreate2(file_id, grpnames[i].c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Gclose(grp);
    }
  }

  /* -------------------------------------- */
  /* Create and fill the Parameters dataset */
  /* -------------------------------------- */

  d[0] = 3;

  status = H5LTmake_dataset(file_id, "/Parameters/LxLyLz", 1, d, H5T_NATIVE_DOUBLE, LxLyLz);
  status = H5LTmake_dataset(file_id, "/Parameters/ncell" , 1, d, H5T_NATIVE_INT   , dglob);

}

void PHDF5fileClass::ClosePHDF5file(){

  H5Fclose(file_id);

}

int PHDF5fileClass::WritePHDF5dataset(string grpname, string datasetname, const_arr3_double data, int nx, int ny, int nz)
{

  /* -------------------------- */
  /* Local variables and arrays */
  /* -------------------------- */

  string dname;
  double *buffer;

  hid_t const h5type = H5T_NATIVE_DOUBLE;

  hid_t glob_dspace;
  hid_t locl_dspace;
  hid_t dataset_prop;
  hid_t dataset;
  hid_t dataspace;
  hid_t dataset_xfer;

  /* --------------------------------- */
  /* Check that dimensions are correct */
  /* --------------------------------- */

  if (bparticles && grpname.c_str()=="Particles"){
    warning_printf("Particle data is not going to be written,"
      " because the 'bparticles' flag is currently turn to FALSE");
    return (2);
  }

  for (int i=0; i<ndim; i++){
    if (dim[i]%chdim[i]!=0){
      eprintf("Grid size is not a multiple of the chunk size in the %d dimension,"
        "\tGlob: %d %d %d\n"
        "\tLocl: %d %d %d\n",
        i, dim[0],dim[1],dim[2],
	chdim[0],chdim[1],chdim[2]);
      return 1;
    }
  }

  /* ----------------------- */
  /* Copy raw data to buffer */
  /* ----------------------- */

  if (nx!=chdim[0] || ny!=chdim[1] || nz!=chdim[2]){
    eprintf("data size is not equal to HDF5 chunk size ");
    return 1;
  }

  buffer = new double[nx*ny*nz];
  int l = 0;
  for (int i = 0; i < nx; i++)
    for (int j = 0; j < ny; j++)
      for (int k = 0; k < nz; k++)
        buffer[l++] = data[i][j][k];

  /* -------------------------------------------------------- */
  /* 5- Set the stride, count and block values for each chunk */
  /*    And set the offset for each chunk                     */
  /* -------------------------------------------------------- */

  hsize_t *stride = new hsize_t[ndim];
  hsize_t *count  = new hsize_t[ndim];
  hsize_t *block  = new hsize_t[ndim];
  hsize_t *offset = new hsize_t[ndim];

  for (int i=0; i<ndim; i++){
    stride[i] = 1;
    count[i]  = 1;
    block[i]  = chdim[i];
    offset[i] = mpicoord[i]*chdim[i];
  }

  /* ---------------------------------- */
  /* 6- Create data spaces for our data */
  /* ---------------------------------- */

  glob_dspace = H5Screate_simple(ndim, dim,   NULL);
  locl_dspace = H5Screate_simple(ndim, chdim, NULL);

  /* --------------------------------------- */
  /* 7- Create the dataset for the HDF5 file */
  /* --------------------------------------- */

  dataset_prop = H5Pcreate(H5P_DATASET_CREATE);

  H5Pset_chunk(dataset_prop, ndim, chdim);

  dname   = "/"+grpname+"/"+datasetname;
  dataset = H5Dcreate2(file_id, dname.c_str(), h5type, glob_dspace, H5P_DEFAULT, dataset_prop, H5P_DEFAULT);

  H5Pclose(dataset_prop);

  dataspace = H5Dget_space(dataset);
  H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, stride, count, block);

  /* --------------------------------- */
  /* 8- Set the parallel transfer mode */
  /* --------------------------------- */

  dataset_xfer = H5Pcreate(H5P_DATASET_XFER);
  #ifdef USING_PARALLEL_HDF5
  H5Pset_dxpl_mpio(dataset_xfer, H5FD_MPIO_COLLECTIVE);
  #else
  eprintf("WriteMethod==Parallel in input file "
          "requires setting USING_PARALLEL_HDF5 in ipicdefs.h");
  #endif

  /* ---------------------------- */
  /* 9- Write data to the dataset */
  /* ---------------------------- */

  H5Dwrite(dataset, h5type, locl_dspace, dataspace, dataset_xfer, buffer);
  //
  if(!MPIdata::get_rank())
  {
    // show the output mode that is actually being used
    H5D_mpio_actual_io_mode_t actual_io_mode;
    H5Pget_mpio_actual_io_mode(dataset_xfer, &actual_io_mode);
    if(actual_io_mode == H5D_MPIO_NO_COLLECTIVE)
      dprintf("H5D_MPIO_NO_COLLECTIVE");
    else if(actual_io_mode == H5D_MPIO_CHUNK_INDEPENDENT)
      dprintf("H5D_MPIO_CHUNK_INDEPENDENT");
    else if(actual_io_mode == H5D_MPIO_CHUNK_COLLECTIVE)
      dprintf("H5D_MPIO_CHUNK_COLLECTIVE");
    else if(actual_io_mode == H5D_MPIO_CONTIGUOUS_COLLECTIVE)
      dprintf("H5D_MPIO_CONTIGUOUS_COLLECTIVE");
    else
      dprintf("unrecognized output method");

    // show the chunking that is actually used
    H5D_mpio_actual_chunk_opt_mode_t actual_chunk_opt_mode;
    H5Pget_mpio_actual_chunk_opt_mode(dataset_xfer, &actual_chunk_opt_mode);
    if(actual_chunk_opt_mode == H5D_MPIO_NO_CHUNK_OPTIMIZATION)
      dprintf("H5D_MPIO_NO_CHUNK_OPTIMIZATION");
    else if(actual_chunk_opt_mode == H5D_MPIO_MULTI_CHUNK)
      dprintf("H5D_MPIO_MULTI_CHUNK");
    //else if(actual_chunk_opt_mode == H5D_MPIO_MULTI_CHUNK_NO_OPT)
    //  dprintf("H5D_MPIO_MULTI_CHUNK_NO_OPT");
    else if(actual_chunk_opt_mode == H5D_MPIO_LINK_CHUNK)
      dprintf("H5D_MPIO_LINK_CHUNK");
    else
      dprintf("unrecognized chunking method");
  }

  delete [] buffer;

  /* ------------------------------------------------------ */
  /* Close dataset related variables created with H5*create */
  /* ------------------------------------------------------ */

  H5Pclose(dataset_xfer);
  H5Dclose(dataset);
  H5Sclose(locl_dspace);

  /* ---------------------------------------------------- */
  /* Close the remaining variables created with H5*create */
  /* ---------------------------------------------------- */

  delete [] stride;
  delete [] count;
  delete [] block;
  delete [] offset;

  return 0;
}

void PHDF5fileClass::ReadPHDF5param(){

  herr_t  status;
  string  dname;
  int     datadims[3];
  double  L[3];

  dname   = "/Parameters/ncell";
  status = H5LTread_dataset_int(file_id, dname.c_str(), datadims);

  dname   = "/Parameters/LxLyLz";
  status = H5LTread_dataset_double(file_id, dname.c_str(), L);

  ndim = 3;
  if (datadims[0]<=1 || datadims[1]<=1 || datadims[2]<=1) ndim = 2;

  for (int i=0; i<ndim; i++){
    dim[i]    = datadims[i];
    LxLyLz[i] = L[i];
  }

}

void PHDF5fileClass::ReadPHDF5dataset_double(string datasetname, arr3_double data){

  herr_t  status;
  double *filedata;

  filedata = new double[dim[0]*dim[1]*dim[2]];

  status = H5LTread_dataset_double(file_id, datasetname.c_str(), filedata);

  for (int i=0; i<dim[0]; i++)
    for (int j=0; j<dim[1]; j++)
      for (int k=0; k<dim[2]; k++)
        data[i][j][k]=filedata[i+j*dim[2]+k*dim[1]*dim[0]];

}

int PHDF5fileClass::getPHDF5ncx(){
  return (int)dim[0];
}

int PHDF5fileClass::getPHDF5ncy(){
  if (ndim<2) return 1;
  return (int)dim[1];
}

int PHDF5fileClass::getPHDF5ncz(){
  if (ndim<3) return 1;
  return (int)dim[2];
}

int PHDF5fileClass::getPHDF5ndim(){
  return ndim;
}

#endif
