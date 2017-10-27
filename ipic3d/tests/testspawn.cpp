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

#include "mpi.h"
#include <stdio.h>
#include <iostream>
#include "../utility/debug.cpp"
#include "../utility/asserts.cpp"
#include "../utility/errors.cpp"
#include "../utility/MPIdata.cpp"
// not sure on what systems this will be available
//#include <ext/stdio_filebuf.h>
using namespace std;

// to run this program on deep:
//
//   module use $IPIC_HOME/env/deep
//   module load ipic-offload
//
//   cd $IPIC_HOME/tests
//   mkdir build; cd build
//   cmake ..
//   make
//   isession
//   mpiexec -np 2 -u 4 ./testspawn
//

//int MPIrank;
//int MPInumprocs;

#define debugout cout << "(" << MPIrank << ")" \
  << "DEBUG " << __func__ << ", " << __FILE__ << ":" << __LINE__ << ": "

// hacked mechanisms to print output one process at a time
//
// for some reason I still get threads overwriting one another,
// although not as badly as without this.
//
void output_barrier()
{
  fflush(stdout);
  cout << flush;
  MPI_Barrier(MPI_COMM_WORLD);
}
int barrier_ret_nprocs()
{
  output_barrier();
  return MPIdata::get_nprocs();
}
int barrier_ret_1()
{
  output_barrier();
  return 1;
}
//#define criticalout \
//  for(int rank=0;rank<barrier_ret_nprocs();rank++) \
//      if(rank==MPIdata::get_rank()) 
#define criticalout
#define masterout \
  for(int i=0;i<barrier_ret_1();i++) \
    if(!MPIdata::get_rank())
// Creating MPI barriers to atomize output is a poor solution;
#define dout criticalout debugout
// a better solution is to write output to separate files:
//ofstream outfile;
//#define dout outfile;

void test_spawn(char**argv)
{
  const int nprocs = 2;
  // define communicator
  MPI_Comm comm;
  const int reorder = 1;
  const int ndims = 1;
  /*const*/ int dims[ndims] = {2};
  /*const*/ int periods[ndims] = {1};
  const int XDIR = 0;
  const int UPWARD = 1;
  const int DNWARD = -1;
  // create a communicator for upward communication
  MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, reorder, &comm);

  int cart_rank;
  MPI_Comm_rank(comm, &cart_rank);
  criticalout dprint(cart_rank);

  int coords[ndims];
  MPI_Cart_coords(comm, cart_rank, ndims, coords);
  criticalout dprint(coords[0]);

  //MPI_Comm_dup(comm, &dn_comm);
  //MPI_Comm_dup(MPI_COMM_WORLD, &up_comm);
  //MPI_Comm_dup(MPI_COMM_WORLD, &dn_comm);

  // This does not actually do a shift; rather, it returns the
  // ranks of the neighbors that would be used to do this shift.
  // We use this as a mechanism to identify neighbors.
  // Shifting upward means that rank_source will be
  // lower_neighbor and rank_dest will be upper_neighbor.
  int lower_neighbor;
  int upper_neighbor;
  MPI_Cart_shift(comm, XDIR, UPWARD, &lower_neighbor, &upper_neighbor);
  criticalout dprint(lower_neighbor); // from lower
  criticalout dprint(upper_neighbor); // to upper

  // showing that we can propagate a message upward
  if(1)
  {
    const int count=16;
    char recvbuff[count];
    char sendbuff[count];
    const int sendcount = 1+ // because of terminating '\0'
      snprintf(sendbuff,count,"hello from %d", MPIdata::get_rank());
    criticalout dprint(sendbuff);
    // receive message from lower
    MPI_Request recv_request;
    MPI_Irecv(recvbuff, count, MPI_CHAR, lower_neighbor, 0, comm, &recv_request);
    // send a message to upper
    MPI_Request send_request;
    MPI_Isend(sendbuff, sendcount, MPI_CHAR, upper_neighbor, 0, comm, &send_request);
    // wait for message to arrive
    MPI_Status recv_status;
    MPI_Wait(&recv_request, &recv_status);
    // print the message
    criticalout dprint(recvbuff);
  }

  // spawn more processes
  if(1)
  {
    MPI_Comm intercomm;
    MPI_Comm_get_parent(&intercomm);
    const bool was_spawned =  (MPI_COMM_NULL == intercomm) ? false : true;
    if(!was_spawned)
    {
      dprintf("hey, let's make some babies!");
      MPI_Comm_spawn(
        "./testspawn",
        &argv[1],
        nprocs,
        MPI_INFO_NULL,
        0,
        MPI_COMM_WORLD,
        &intercomm,
        MPI_ERRCODES_IGNORE);
    }
    else
    {
      dprintf("waaaaa!");
    }
  }

  // probably MPI_Finalize takes care of this anyway...
  MPI_Comm_free(&comm);
}

int main(int argc, char **argv)
{
  MPIdata::init(&argc, &argv);
  //MPI_Init(&argc, &argv);
  //MPI_Comm_rank(MPI_COMM_WORLD, &MPIrank);
  //MPI_Comm_size(MPI_COMM_WORLD, &MPInumprocs);

  test_spawn(argv);

  MPIdata::finalize_mpi();
  //MPI_Finalize();
}
