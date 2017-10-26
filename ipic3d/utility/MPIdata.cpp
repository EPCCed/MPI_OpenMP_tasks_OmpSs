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

#ifndef NO_MPI
  #include <mpi.h>
#endif
//#include <stdio.h>
#include <assert.h>
#include "MPIdata.h"
#ifndef NO_MPI
#include <mpi.h>
#endif
#include "ompdefs.h" // for omp_get_max_threads
#include <iostream>

using std::cout;
// code to check that init() is called before instance()
//
// no need for this to have more than file scope
int MPIdata::rank=-1;
int MPIdata::nprocs=-1;
static bool MPIdata_is_initialized=false;
bool MPIdata_assert_initialized()
{
  assert(MPIdata_is_initialized);
  return true;
}

MPIdata& MPIdata::instance()
{
  // This is executed on the first call to check that
  // MPIdata has first been initialized.
  static bool check = MPIdata_assert_initialized();
  static MPIdata* instance = new MPIdata;
  // After the first call, this is the only line
  // that is actually executed.
  return *instance;
}

void MPIdata::init(int *argc, char ***argv) {
  assert(!MPIdata_is_initialized);
  int provided;
 #ifdef NO_MPI
  rank = 0;
  nprocs = 1;
 #else // NO_MPI
  /* Initialize the MPI API */
  //  MPI_Init(argc, argv);
  MPI_Init_thread( 0, 0, MPI_THREAD_MULTIPLE, &provided );

  /* Set rank */
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  if (rank == 0)
    {
      if(provided == MPI_THREAD_SINGLE){
	cout << "The provided mode for the thread support is " << provided << " or MPI_THREAD_SINGLE" << "\n";
      }
      if(provided == MPI_THREAD_FUNNELED){
	cout << "The provided mode for the thread support is " << provided << " or MPI_THREAD_FUNNELED" << "\n";
      } 
      if(provided == MPI_THREAD_SERIALIZED){
	cout << "The provided mode for the thread support is " << provided << " or MPI_THREAD_SERIALIZED" << "\n";
      }
      if(provided == MPI_THREAD_MULTIPLE){
	cout << "The provided mode for the thread support is " << provided << " or MPI_THREAD_MULTIPLE" << "\n";
      }
    }

  /* Set total number of processors */
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
 #endif // NO_MPI

  MPIdata_is_initialized = true;
}

void MPIdata::exit(int code) {
  finalize_mpi();
  ::exit(code);
}

void MPIdata::finalize_mpi() {
 #ifndef NO_MPI
  MPI_Finalize();
 #endif
}

void MPIdata::Print(void) {
  printf("\n"
    "Number of processes = %d\n"
    "-------------------------\n"
    "Number of threads = %d\n"
    "-------------------------\n",
     get_nprocs(),
     omp_get_max_threads());
}

// extern MPIdata *mpi; // instantiated in iPIC3D.cpp

