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

#include "Com3DNonblk.h"
#include "ompdefs.h"
//#include "HelpersCom3DNonblk.h"


//isCenterFlag: 1 communicateCenter; 0 communicateNode
void NBDerivedHaloComm(int nx, int ny, int nz, double ***vector,const VirtualTopology3D * vct, EMfields3D *EMf,bool isCenterFlag, bool isFaceOnlyFlag, bool needInterp, bool isParticle)
{

  MPI_Errhandler_set(MPI_COMM_WORLD,MPI_ERRORS_RETURN);
  MPI_Status  stat[12];
  MPI_Request reqList[12];          //at most 6 requests x 2 (send recv)
  int communicationCnt[6] = {0,0,0,0,0,0};  //1 if there communication on that dir
  int recvcnt = 0,sendcnt = 0;              //request counter

  const int tag_XL=1,tag_YL=2,tag_ZL=3,tag_XR=4,tag_YR=5,tag_ZR=6;//To address same rank as left and right neighbour in periodic case
  const int myrank        = vct->getCartesian_rank();
  const MPI_Comm comm       = isParticle ?vct->getParticleComm()      :vct->getFieldComm();
  const int right_neighborX = isParticle ?vct->getXright_neighbor_P() :vct->getXright_neighbor();
  const int left_neighborX  = isParticle ?vct->getXleft_neighbor_P()  :vct->getXleft_neighbor();
  const int right_neighborY = isParticle ?vct->getYright_neighbor_P() :vct->getYright_neighbor();
  const int left_neighborY  = isParticle ?vct->getYleft_neighbor_P()  :vct->getYleft_neighbor();
  const int right_neighborZ = isParticle ?vct->getZright_neighbor_P() :vct->getZright_neighbor();
  const int left_neighborZ  = isParticle ?vct->getZleft_neighbor_P()  :vct->getZleft_neighbor();
  bool isCenterDim = (isCenterFlag&&!needInterp);//This is to address moment interp use nodes dimension but exchange center face
  const MPI_Datatype yzFacetype = EMf->getYZFacetype(isCenterDim);
  const MPI_Datatype xzFacetype = EMf->getXZFacetype(isCenterDim);
  const MPI_Datatype xyFacetype = EMf->getXYFacetype(isCenterDim);
  const MPI_Datatype xEdgetype = EMf->getXEdgetype(isCenterDim);
  const MPI_Datatype yEdgetype = EMf->getYEdgetype(isCenterDim);
  const MPI_Datatype zEdgetype = EMf->getZEdgetype(isCenterDim);
  const MPI_Datatype xEdgetype2 = EMf->getXEdgetype2(isCenterDim);
  const MPI_Datatype yEdgetype2 = EMf->getYEdgetype2(isCenterDim);
  const MPI_Datatype zEdgetype2 = EMf->getZEdgetype2(isCenterDim);
  const MPI_Datatype cornertype = EMf->getCornertype(isCenterDim);
  int offset;




  //Face Exchange as long as neighbor exits on that direction
  //Tag is based on the sender: XL for sender = XR for receiver
  //Post Recv before Send
  int id_thread, nthreads; 

#pragma omp parallel default(shared) private(id_thread,nthreads) //shared(recvcnt,sendcnt)
{

    id_thread = omp_get_thread_num();
    nthreads = omp_get_num_threads();

    	if (nthreads==2)
	{
		if(id_thread==0)
			if(left_neighborX != MPI_PROC_NULL && left_neighborX != myrank )
			{
				#pragma omp critical
			{
				MPI_Irecv(&vector[0][1][1], 1, yzFacetype, left_neighborX,tag_XR, comm, &reqList[recvcnt]);
				recvcnt++;
				communicationCnt[0] = 1;
			}
			}
		if(id_thread==1)
			if(right_neighborX != MPI_PROC_NULL && right_neighborX != myrank )
			{
				#pragma omp critical
			{ 
				MPI_Irecv(&vector[nx-1][1][1], 1, yzFacetype, right_neighborX,tag_XL, comm, &reqList[recvcnt]);
				recvcnt++;
				communicationCnt[1] = 1;
				}
			}
		if(id_thread==0)
			if(left_neighborY != MPI_PROC_NULL && left_neighborY != myrank )
			{
				#pragma omp critical 
			{
				MPI_Irecv(&vector[1][0][1], 1, xzFacetype, left_neighborY,tag_YR, comm, &reqList[recvcnt]);
				recvcnt++;
				communicationCnt[2] = 1;
			}
			}
		if(id_thread==1)
			if(right_neighborY != MPI_PROC_NULL && right_neighborY != myrank )
			{
				#pragma omp critical 
			{ 
				MPI_Irecv(&vector[1][ny-1][1], 1, xzFacetype, right_neighborY,tag_YL, comm, &reqList[recvcnt]);
					recvcnt++;
					communicationCnt[3] = 1;
			}
			}
		if(id_thread==0)
			if(left_neighborZ != MPI_PROC_NULL && left_neighborZ != myrank )
			{
				#pragma omp critical
			{ 
				MPI_Irecv(&vector[1][1][0],1, xyFacetype, left_neighborZ,tag_ZR, comm, &reqList[recvcnt]);
				recvcnt++; 
				communicationCnt[4] = 1;
				}
			}
		if(id_thread==1)
			if(right_neighborZ != MPI_PROC_NULL&& right_neighborZ != myrank )
			{
				#pragma omp critical 
			{
				MPI_Irecv(&vector[1][1][nz-1], 1, xyFacetype, right_neighborZ,tag_ZL, comm, &reqList[recvcnt]);
				recvcnt++;
				communicationCnt[5] = 1;
				}
			}
		#pragma omp barrier
		if(id_thread==0)
		{
			sendcnt = recvcnt;
			offset = (isCenterFlag ?0:1);
		}
		#pragma omp barrier
		if(id_thread==0)
			if(communicationCnt[0] == 1)
			{
				#pragma omp critical
				{ 
					MPI_Isend(&vector[1+offset][1][1],1, yzFacetype, left_neighborX, tag_XL, comm, &reqList[sendcnt]);
					sendcnt++;
				}
			}
		if(id_thread==1)
			if(communicationCnt[1] == 1)
			{
				#pragma omp critical 
				{
					MPI_Isend(&vector[nx-2-offset][1][1], 1, yzFacetype, right_neighborX,tag_XR, comm, &reqList[sendcnt]);
					sendcnt++;
				}
			}
		if(id_thread==0)
			if(communicationCnt[2] == 1)
			{
				#pragma omp critical 
				{
					MPI_Isend(&vector[1][1+offset][1],1, xzFacetype, left_neighborY, tag_YL, comm, &reqList[sendcnt]);
					sendcnt++;
				}
			}
		if(id_thread==1)
			if(communicationCnt[3] == 1)
			{
				#pragma omp critical 
				{
					MPI_Isend(&vector[1][ny-2-offset][1], 1, xzFacetype, right_neighborY,tag_YR, comm, &reqList[sendcnt]);
					sendcnt++;
				}
			}
		if(id_thread==0)
			if(communicationCnt[4] == 1)
			{
				#pragma omp critical
				{
					MPI_Isend(&vector[1][1][1+offset],1, xyFacetype, left_neighborZ, tag_ZL, comm, &reqList[sendcnt]);
					sendcnt++;
				}
			}
		if(id_thread==1)
			if(communicationCnt[5] == 1)
			{
				#pragma omp critical 
				{
					MPI_Isend(&vector[1][1][nz-2-offset], 1, xyFacetype, right_neighborZ,tag_ZR, comm, &reqList[sendcnt]);
					sendcnt++;
				}
			}
		#pragma omp barrier
		if(id_thread==0)
		{
			assert_eq(recvcnt,sendcnt-recvcnt);
		}
		#pragma omp barrier
		//Buffer swap if any (done before waiting for receiving msg done to delay sync)
		if(id_thread==0)
		{
			if (right_neighborX == myrank &&left_neighborX== myrank)
			{
				for (int iy = 1; iy < ny-1; iy++)
					for (int iz = 1; iz < nz-1; iz++)
					{
						vector[0][iy][iz] = vector[nx-2][iy][iz];
						vector[nx-1][iy][iz] = vector[1][iy][iz];
					}
			}
			if (right_neighborY == myrank &&left_neighborY == myrank)
			{
				for (int ix = 1; ix < nx-1; ix++)
					for (int iz = 1; iz < nz-1; iz++)
					{
						vector[ix][0][iz] = vector[ix][ny-2][iz];
						vector[ix][ny-1][iz] = vector[ix][1][iz];
					}
			}
			if (right_neighborZ == myrank &&left_neighborZ == myrank)
			{
				for (int ix = 1; ix < nx-1; ix++)
					for (int iy = 1; iy < ny-1; iy++)
					{
						vector[ix][iy][0]= vector[ix][iy][nz-2];
						vector[ix][iy][nz-1] = vector[ix][iy][1];
					}
			}
		}
		#pragma omp barrier
		//Need to finish receiving + sending before edge exchange
		if(id_thread==0)
		{
			if(sendcnt>0)
			{
				MPI_Waitall(sendcnt,&reqList[0], &stat[0]);
				bool stopFlag = false;
				for(int si=0;si< sendcnt;si++)
				{
					int error_code = stat[si].MPI_ERROR;
					if (error_code != MPI_SUCCESS)
					{
						stopFlag = true;
						char error_string[100];
						int length_of_error_string, error_class;
						MPI_Error_class(error_code, &error_class);
						MPI_Error_string(error_class, error_string, &length_of_error_string);
						dprintf("MPI_Waitall error at %d%s\n",si, error_string);
					}
				}
				if(stopFlag)
					exit (EXIT_FAILURE);
			}
		}
		#pragma omp barrier
		if( !isFaceOnlyFlag )
		{
			//Exchange yEdge only when Z X neighbours exist
			//if Zleft + Zright, use merged yEdgeType2 to send two edges in one msg
			//Otherwise, only send one yEdge
			#pragma omp barrier
			if(id_thread==0)
			{
				recvcnt = 0, sendcnt = 0;
			}
			#pragma omp barrier
			if(id_thread==0)
				if(communicationCnt[0] == 1)
				{
					if(communicationCnt[4] == 1 && communicationCnt[5] == 1)
					{
						#pragma omp critical
						MPI_Irecv(&vector[0][1][0], 1,yEdgetype2, left_neighborX, tag_XR, comm, &reqList[recvcnt++]);
					}
					else if(communicationCnt[4] == 1)
					{
						#pragma omp critical
						MPI_Irecv(&vector[0][1][0], 1,yEdgetype, left_neighborX, tag_XR, comm, &reqList[recvcnt++]);
					}
					else if(communicationCnt[5] == 1)
					{
						#pragma omp critical
						MPI_Irecv(&vector[0][1][nz-1],1,yEdgetype, left_neighborX, tag_XR, comm, &reqList[recvcnt++]);
					}
				}
			if(id_thread==1)
				if(communicationCnt[1] == 1)
				{
					if(communicationCnt[4] == 1 && communicationCnt[5] == 1)
					{
						#pragma omp critical
						MPI_Irecv(&vector[nx-1][1][0],1,yEdgetype2,right_neighborX,tag_XL, comm, &reqList[recvcnt++]);
					}
					else if(communicationCnt[4] == 1)
					{
						#pragma omp critical
						MPI_Irecv(&vector[nx-1][1][0], 1, yEdgetype, right_neighborX, tag_XL, comm, &reqList[recvcnt++]);
					}
					else if(communicationCnt[5] == 1)
					{
						#pragma omp critical
						MPI_Irecv(&vector[nx-1][1][nz-1],1, yEdgetype, right_neighborX, tag_XL, comm, &reqList[recvcnt++]);
					}
				}
			//Exchange zEdge only when X Y neighbours exist
			//if Xleft + Xright, use merged zEdgeType2
			//Otherwise, only send one zEdge
			if(id_thread==0)
				if(communicationCnt[2] == 1)
				{
					if(communicationCnt[0] == 1 && communicationCnt[1] == 1)
					{
						#pragma omp critical
						MPI_Irecv(&vector[0][0][1], 1, zEdgetype2,left_neighborY, tag_YR,comm, &reqList[recvcnt++]);
					}
					else if(communicationCnt[0] == 1)
					{
						#pragma omp critical
						MPI_Irecv(&vector[0][0][1], 1, zEdgetype, left_neighborY, tag_YR,comm, &reqList[recvcnt++]);
					}
					else if(communicationCnt[1] == 1)
					{
						#pragma omp critical
						MPI_Irecv(&vector[nx-1][0][1],1, zEdgetype, left_neighborY, tag_YR,comm, &reqList[recvcnt++]);
					}
				}
			if(id_thread==1)
				if(communicationCnt[3] == 1)
				{
					if(communicationCnt[0] == 1 && communicationCnt[1] == 1)
					{
						#pragma omp critical
						MPI_Irecv(&vector[0][ny-1][1], 1, zEdgetype2,right_neighborY,tag_YL,comm, &reqList[recvcnt++]);
					}
					else if(communicationCnt[0] == 1)
					{
						#pragma omp critical
						MPI_Irecv(&vector[0][ny-1][1], 1, zEdgetype, right_neighborY,tag_YL,comm, &reqList[recvcnt++]);
					}
					else if(communicationCnt[1] == 1)
					{
						#pragma omp critical
						MPI_Irecv(&vector[nx-1][ny-1][1],1, zEdgetype, right_neighborY,tag_YL,comm, &reqList[recvcnt++]);
					}
				}
			//Exchange xEdge only when Y Z neighbours exist
			//if Yleft + Yright exist, use merged xEdgeType2
			//Otherwise, only send one xEdge
			if(id_thread==0)
				if(communicationCnt[4] == 1)
				{
					if(communicationCnt[2] == 1 && communicationCnt[3] == 1)
					{
						#pragma omp critical
						MPI_Irecv(&vector[1][0][0],1,xEdgetype2,left_neighborZ, tag_ZR, comm, &reqList[recvcnt++]);
					}
					else if(communicationCnt[2] == 1)
					{
						#pragma omp critical
						MPI_Irecv(&vector[1][0][0],1,xEdgetype, left_neighborZ, tag_ZR,comm, &reqList[recvcnt++]);
					}
					else if(communicationCnt[3] == 1)
					{
						#pragma omp critical
						MPI_Irecv(&vector[1][ny-1][0],1, xEdgetype, left_neighborZ,tag_ZR,comm, &reqList[recvcnt++]);
					}
				}
			if(id_thread==1)
				if(communicationCnt[5] == 1)
				{
					if(communicationCnt[2] == 1 && communicationCnt[3] == 1)
					{
						#pragma omp critical
						MPI_Irecv(&vector[1][0][nz-1],1,xEdgetype2,right_neighborZ, tag_ZL,comm, &reqList[recvcnt++]);
					}
					else if(communicationCnt[2] == 1)
					{
						#pragma omp critical
						MPI_Irecv(&vector[1][0][nz-1],1, xEdgetype,right_neighborZ, tag_ZL,comm, &reqList[recvcnt++]);
					}
					else if(communicationCnt[3] == 1)
					{
						#pragma omp critical
						MPI_Irecv(&vector[1][ny-1][nz-1],1,xEdgetype,right_neighborZ, tag_ZL,comm, &reqList[recvcnt++]);
					}
				}
			#pragma omp barrier
			if(id_thread==0)
				sendcnt = recvcnt;
			#pragma omp barrier
			if(id_thread==0)
				if(communicationCnt[0] == 1)
				{
					if(communicationCnt[4] == 1 && communicationCnt[5] == 1)
					{
						#pragma omp critical
						MPI_Isend(&vector[1][1][0], 1,yEdgetype2,left_neighborX, tag_XL, comm, &reqList[sendcnt++]);
					}
					else if(communicationCnt[4] == 1)
					{
						#pragma omp critical
						MPI_Isend(&vector[1][1][0], 1,yEdgetype, left_neighborX, tag_XL, comm, &reqList[sendcnt++]);
					}
					else if(communicationCnt[5] == 1)
					{
						#pragma omp critical
						MPI_Isend(&vector[1][1][nz-1],1,yEdgetype, left_neighborX, tag_XL, comm, &reqList[sendcnt++]);
					}
				}
			if(id_thread==1)
				if(communicationCnt[1] == 1)
				{
					if(communicationCnt[4] == 1 && communicationCnt[5] == 1)
					{
						#pragma omp critical
						MPI_Isend(&vector[nx-2][1][0],1, yEdgetype2,right_neighborX, tag_XR, comm,&reqList[sendcnt++]);
					}
					else if(communicationCnt[4] == 1)
					{
						#pragma omp critical
						MPI_Isend(&vector[nx-2][1][0],1, yEdgetype,right_neighborX, tag_XR, comm, &reqList[sendcnt++]);
					}
					else if(communicationCnt[5] == 1)
					{
						#pragma omp critical
						MPI_Isend(&vector[nx-2][1][nz-1],1,yEdgetype,right_neighborX, tag_XR,comm,&reqList[sendcnt++]);
					}
				}
			if(id_thread==0)
				if(communicationCnt[2] == 1)
				{
					if(communicationCnt[0] == 1 && communicationCnt[1] == 1)
					{
						#pragma omp critical
						MPI_Isend(&vector[0][1][1],1, zEdgetype2, left_neighborY, tag_YL, comm, &reqList[sendcnt++]);
					}
					else if(communicationCnt[0] == 1)
					{
						#pragma omp critical
						MPI_Isend(&vector[0][1][1],1, zEdgetype,left_neighborY, tag_YL, comm, &reqList[sendcnt++]);
					}
					else if(communicationCnt[1] == 1)
					{
						#pragma omp critical
						MPI_Isend(&vector[nx-1][1][1],1,zEdgetype,left_neighborY, tag_YL, comm, &reqList[sendcnt++]);
					}
				}
			if(id_thread==1)
				if(communicationCnt[3] == 1)
				{
					if(communicationCnt[0] == 1 && communicationCnt[1] == 1)
					{
						#pragma omp critical
						MPI_Isend(&vector[0][ny-2][1],1, zEdgetype2, right_neighborY, tag_YR, comm, &reqList[sendcnt++]);
					}
					else if(communicationCnt[0] == 1)
					{
						#pragma omp critical
						MPI_Isend(&vector[0][ny-2][1],1, zEdgetype,right_neighborY, tag_YR, comm, &reqList[sendcnt++]);
					}
					else if(communicationCnt[1] == 1)
					{
						#pragma omp critical
						MPI_Isend(&vector[nx-1][ny-2][1],1,zEdgetype,right_neighborY, tag_YR, comm, &reqList[sendcnt++]);
					}
				}
			if(id_thread==0)
				if(communicationCnt[4] == 1)
				{
					if(communicationCnt[2] == 1 && communicationCnt[3] == 1)
					{
						#pragma omp critical
						MPI_Isend(&vector[1][0][1],1, xEdgetype2, left_neighborZ, tag_ZL, comm, &reqList[sendcnt++]);
					}
					else if(communicationCnt[2] == 1)
					{
						#pragma omp critical
						MPI_Isend(&vector[1][0][1],1, xEdgetype,left_neighborZ, tag_ZL, comm, &reqList[sendcnt++]);
					}
					else if(communicationCnt[3] == 1)
					{
						#pragma omp critical
						MPI_Isend(&vector[1][ny-1][1],1,xEdgetype,left_neighborZ, tag_ZL, comm, &reqList[sendcnt++]);
					}
				}
			if(id_thread==1)
				if(communicationCnt[5] == 1)
				{
					if(communicationCnt[2] == 1 && communicationCnt[3] == 1)
					{
						#pragma omp critical
						MPI_Isend(&vector[1][0][nz-2],1,xEdgetype2, right_neighborZ, tag_ZR, comm, &reqList[sendcnt++]);
					}
					else if(communicationCnt[2] == 1)
					{
						#pragma omp critical
						MPI_Isend(&vector[1][0][nz-2],1,xEdgetype,right_neighborZ, tag_ZR, comm, &reqList[sendcnt++]);
					}
					else if(communicationCnt[3] == 1)
					{
						#pragma omp critical
						MPI_Isend(&vector[1][ny-1][nz-2],1, xEdgetype,right_neighborZ, tag_ZR, comm, &reqList[sendcnt++]);
					}
				}
			#pragma omp barrier
			if(id_thread==0)
				assert_eq(recvcnt,sendcnt-recvcnt);
			#pragma omp barrier
			//Swap Local Edges
			if(id_thread==0)
			{
				if(right_neighborX == myrank &&left_neighborX== myrank)
				{
					if(right_neighborZ != MPI_PROC_NULL )
					{
						for (int iy = 1; iy < ny-1; iy++)
						{
							vector[0][iy][nz-1]= vector[nx-2][iy][nz-1];
							vector[nx-1][iy][nz-1] = vector[1][iy][nz-1];
						}
					}
					if(left_neighborZ != MPI_PROC_NULL )
					{
						for (int iy = 1; iy < ny-1; iy++)
						{
							vector[0][iy][0]= vector[nx-2][iy][0];
							vector[nx-1][iy][0] = vector[1][iy][0];
						}
					}
					if(right_neighborY != MPI_PROC_NULL )
					{
						for (int iz = 1; iz < nz-1; iz++)
						{
							vector[0][ny-1][iz]= vector[nx-2][ny-1][iz];
							vector[nx-1][ny-1][iz] = vector[1][ny-1][iz];
						}
					}
					if(left_neighborY != MPI_PROC_NULL )
					{
						for (int iz = 1; iz < nz-1; iz++)
						{
							vector[0][0][iz] = vector[nx-2][0][iz];
							vector[nx-1][0][iz]= vector[1][0][iz];
						}
					}
				}
				if(right_neighborY == myrank &&left_neighborY == myrank)
				{
					if(right_neighborX != MPI_PROC_NULL )
					{
						for (int iz = 1; iz < nz-1; iz++)
						{
							vector[nx-1][0][iz]= vector[nx-1][ny-2][iz];
							vector[nx-1][ny-1][iz] = vector[nx-1][1][iz];
						}
					}
					if(left_neighborX != MPI_PROC_NULL)
					{
						for (int iz = 1; iz < nz-1; iz++)
						{
							vector[0][0][iz] = vector[0][ny-2][iz];
							vector[0][ny-1][iz]= vector[0][1][iz];
						}
					}
					if(right_neighborZ != MPI_PROC_NULL)
					{
						for (int ix = 1; ix < nx-1; ix++)
						{
							vector[ix][0][nz-1]= vector[ix][ny-2][nz-1];
							vector[ix][ny-1][nz-1] = vector[ix][1][nz-1];
						}
					}
					if(left_neighborZ != MPI_PROC_NULL)
					{
						for (int ix = 1; ix < nx-1; ix++)
						{
							vector[ix][0][0] = vector[ix][ny-2][0];
							vector[ix][ny-1][0]= vector[ix][1][0];
						}
					}
				}
				if(right_neighborZ == myrank &&left_neighborZ == myrank)
				{
					if(right_neighborY != MPI_PROC_NULL )
					{
						for (int ix = 1; ix < nx-1; ix++)
						{
							vector[ix][ny-1][0]= vector[ix][ny-1][nz-2];
							vector[ix][ny-1][nz-1] = vector[ix][ny-1][1];
						}
					}
					if(left_neighborY != MPI_PROC_NULL )
					{
						for (int ix = 1; ix < nx-1; ix++)
						{
							vector[ix][0][0]= vector[ix][0][nz-2];
							vector[ix][0][nz-1] = vector[ix][0][1];
						}
					}
					if(right_neighborX != MPI_PROC_NULL )
					{
						for (int iy = 1; iy < ny-1; iy++)
						{
							vector[nx-1][iy][0]= vector[nx-1][iy][nz-2];
							vector[nx-1][iy][nz-1] = vector[nx-1][iy][1];
						}
					}
					if(left_neighborX != MPI_PROC_NULL )
					{
						for (int iy = 1; iy < ny-1; iy++)
						{
							vector[0][iy][0]= vector[0][iy][nz-2];
							vector[0][iy][nz-1] = vector[0][iy][1];
						}
					}
				}
			}
			//Need to finish receiving edges for corner exchange
			#pragma omp barrier
			if(id_thread==0)
				if(sendcnt>0)
				{
					MPI_Waitall(sendcnt,&reqList[0],&stat[0]);
					bool stopFlag = false;
					for(int si=0;si< sendcnt;si++)
					{
						int error_code = stat[si].MPI_ERROR;
						if (error_code != MPI_SUCCESS)
						{
							stopFlag = true;
							char error_string[100];
							int length_of_error_string, error_class;

							MPI_Error_class(error_code, &error_class);
							MPI_Error_string(error_class, error_string, &length_of_error_string);
							dprintf("MPI_Waitall error at %d%s\n",si, error_string);
						}
					}
					if(stopFlag)
						exit (EXIT_FAILURE);;
				}
			#pragma omp barrier
			//Corner Exchange only needed if XYZ neighbours all exist
			//4 corners communicated in one message
			//Assume Non-periodic will be handled in BC
			//Define corner types for X communication
			if(id_thread==0)
				recvcnt = 0,sendcnt = 0;
			#pragma omp barrier
			if((communicationCnt[2] == 1 || communicationCnt[3] == 1) && (communicationCnt[4] == 1 || communicationCnt[5] == 1))
			{
				//if XLeft exists, send 4 corners to XLeft
				if(id_thread==0)
					if(communicationCnt[0] == 1)
					{
						#pragma omp critical
						MPI_Irecv(&vector[0][0][0],1, cornertype, left_neighborX, tag_XR,comm, &reqList[recvcnt++]);
					}
				//if XRight exist
				if(id_thread==1)
					if(communicationCnt[1] == 1)
					{
						#pragma omp critical
						MPI_Irecv(&vector[nx-1][0][0],1,cornertype, right_neighborX, tag_XL,comm, &reqList[recvcnt++]);
					}
				#pragma omp barrier
				if(id_thread==0)
					sendcnt=recvcnt;
				#pragma omp barrier
				if(id_thread==0)
					if(communicationCnt[0] == 1)
					{
						#pragma omp critical
						MPI_Isend(&vector[1][0][0], 1,cornertype,left_neighborX, tag_XL, comm, &reqList[sendcnt++]);
					}
				if(id_thread==1)
					if(communicationCnt[1] == 1)
					{
						#pragma omp critical
						MPI_Isend(&vector[nx-2][0][0],1,cornertype,right_neighborX, tag_XR,comm, &reqList[sendcnt++]);
					}
			}
			#pragma omp barrier
			if(id_thread==0)
				assert_eq(recvcnt,sendcnt-recvcnt);
			#pragma omp barrier
			//Delay local data copy
			if(id_thread==0)
			{
				if (left_neighborX== myrank && right_neighborX == myrank)
				{
					if( (left_neighborY != MPI_PROC_NULL) && (left_neighborZ != MPI_PROC_NULL))
					{
						vector[0][0][0] = vector[nx-2][0][0];
						vector[nx-1][0][0]= vector[1][0][0];
					}
					if( (left_neighborY != MPI_PROC_NULL) && (right_neighborZ != MPI_PROC_NULL))
					{
						vector[0][0][nz-1]= vector[nx-2][0][nz-1];
						vector[nx-1][0][nz-1] = vector[1][0][nz-1];
					}
					if( (right_neighborY != MPI_PROC_NULL) && (left_neighborZ != MPI_PROC_NULL))
					{
						vector[0][ny-1][0]= vector[nx-2][ny-1][0];
						vector[nx-1][ny-1][0] = vector[1][ny-1][0];
					}
					if( (right_neighborY != MPI_PROC_NULL) && (right_neighborZ != MPI_PROC_NULL))
					{
						vector[0][ny-1][nz-1]= vector[nx-2][ny-1][nz-1];
						vector[nx-1][ny-1][nz-1] = vector[1][ny-1][nz-1];
					}
				}
				else if (left_neighborY== myrank && right_neighborY == myrank)
				{
					if( (left_neighborX != MPI_PROC_NULL) && (left_neighborZ != MPI_PROC_NULL))
					{
						vector[0][0][0] = vector[0][ny-2][0];
						vector[0][ny-1][0]= vector[0][1][0];
					}
					if( (left_neighborX != MPI_PROC_NULL) && (right_neighborZ != MPI_PROC_NULL))
					{
						vector[0][0][nz-1]= vector[0][ny-2][nz-1];
						vector[0][ny-1][nz-1]= vector[0][1][nz-1];
					}
					if( (right_neighborX != MPI_PROC_NULL) && (left_neighborZ != MPI_PROC_NULL))
					{
						vector[nx-1][0][0]= vector[nx-1][ny-2][0];
						vector[nx-1][ny-1][0] = vector[nx-1][1][0];
					}
					if( (right_neighborX != MPI_PROC_NULL) && (right_neighborZ != MPI_PROC_NULL))
					{
						vector[nx-1][0][nz-1]= vector[nx-1][ny-2][nz-1];
						vector[nx-1][ny-1][nz-1] = vector[nx-1][1][nz-1];
					}
				}
				else if (left_neighborZ== myrank && right_neighborZ == myrank)
				{
					if( (left_neighborY != MPI_PROC_NULL) && (left_neighborX != MPI_PROC_NULL))
					{
						vector[0][0][0] = vector[0][0][nz-2];
						vector[0][0][nz-1]= vector[0][0][1];
					}
					if( (left_neighborY != MPI_PROC_NULL) && (right_neighborX != MPI_PROC_NULL))
					{
						vector[nx-1][0][0]= vector[nx-1][0][nz-2];
						vector[nx-1][0][nz-1] = vector[nx-1][0][1];
					}
					if( (right_neighborY != MPI_PROC_NULL) && (left_neighborX != MPI_PROC_NULL))
					{
						vector[0][ny-1][0]= vector[0][ny-1][nz-2];
						vector[0][ny-1][nz-1] = vector[0][ny-1][1];
					}
					if( (right_neighborY != MPI_PROC_NULL) && (right_neighborX != MPI_PROC_NULL))
					{
						vector[nx-1][ny-1][0]= vector[nx-1][ny-1][nz-2];
						vector[nx-1][ny-1][nz-1] = vector[nx-1][ny-1][1];
					}
				}
			}
			if(id_thread==0)
				if(sendcnt>0)
				{
					MPI_Waitall(sendcnt,&reqList[0], &stat[0]);
					bool stopFlag = false;
					for(int si=0;si< sendcnt;si++)
					{
						int error_code = stat[si].MPI_ERROR;
						if (error_code != MPI_SUCCESS)
						{
							stopFlag = true;
							char error_string[100];
							int length_of_error_string, error_class;
							MPI_Error_class(error_code, &error_class);
							MPI_Error_string(error_class, error_string, &length_of_error_string);
							dprintf("MPI_Waitall error at %d%s\n",si, error_string);
						}
					}
					if(stopFlag)
						exit (EXIT_FAILURE);
				}
			#pragma omp barrier
		}
	}
    else if (nthreads==4)
	{
		if(id_thread==0)
			if(left_neighborX != MPI_PROC_NULL && left_neighborX != myrank )
			{
				#pragma omp critical
			{
				MPI_Irecv(&vector[0][1][1], 1, yzFacetype, left_neighborX,tag_XR, comm, &reqList[recvcnt]);
				recvcnt++;
				communicationCnt[0] = 1;
			}
			}
		if(id_thread==1)
			if(right_neighborX != MPI_PROC_NULL && right_neighborX != myrank )
			{
				#pragma omp critical
			{ 
				MPI_Irecv(&vector[nx-1][1][1], 1, yzFacetype, right_neighborX,tag_XL, comm, &reqList[recvcnt]);
				recvcnt++;
				communicationCnt[1] = 1;
				}
			}
		if(id_thread==2)
			if(left_neighborY != MPI_PROC_NULL && left_neighborY != myrank )
			{
				#pragma omp critical 
			{
				MPI_Irecv(&vector[1][0][1], 1, xzFacetype, left_neighborY,tag_YR, comm, &reqList[recvcnt]);
				recvcnt++;
				communicationCnt[2] = 1;
			}
			}
		if(id_thread==3)
			if(right_neighborY != MPI_PROC_NULL && right_neighborY != myrank )
			{
				#pragma omp critical 
			{ 
				MPI_Irecv(&vector[1][ny-1][1], 1, xzFacetype, right_neighborY,tag_YL, comm, &reqList[recvcnt]);
					recvcnt++;
					communicationCnt[3] = 1;
			}
			}
		if(id_thread==0)
			if(left_neighborZ != MPI_PROC_NULL && left_neighborZ != myrank )
			{
				#pragma omp critical
			{ 
				MPI_Irecv(&vector[1][1][0],1, xyFacetype, left_neighborZ,tag_ZR, comm, &reqList[recvcnt]);
				recvcnt++; 
				communicationCnt[4] = 1;
				}
			}
		if(id_thread==1)
			if(right_neighborZ != MPI_PROC_NULL&& right_neighborZ != myrank )
			{
				#pragma omp critical 
			{
				MPI_Irecv(&vector[1][1][nz-1], 1, xyFacetype, right_neighborZ,tag_ZL, comm, &reqList[recvcnt]);
				recvcnt++;
				communicationCnt[5] = 1;
				}
			}
		#pragma omp barrier
		if(id_thread==0)
		{
			sendcnt = recvcnt;
			offset = (isCenterFlag ?0:1);
		}
		#pragma omp barrier
		if(id_thread==0)
			if(communicationCnt[0] == 1)
			{
				#pragma omp critical
				{ 
					MPI_Isend(&vector[1+offset][1][1],1, yzFacetype, left_neighborX, tag_XL, comm, &reqList[sendcnt]);
					sendcnt++;
				}
			}
		if(id_thread==1)
			if(communicationCnt[1] == 1)
			{
				#pragma omp critical 
				{
					MPI_Isend(&vector[nx-2-offset][1][1], 1, yzFacetype, right_neighborX,tag_XR, comm, &reqList[sendcnt]);
					sendcnt++;
				}
			}
		if(id_thread==2)
			if(communicationCnt[2] == 1)
			{
				#pragma omp critical 
				{
					MPI_Isend(&vector[1][1+offset][1],1, xzFacetype, left_neighborY, tag_YL, comm, &reqList[sendcnt]);
					sendcnt++;
				}
			}
		if(id_thread==3)
			if(communicationCnt[3] == 1)
			{
				#pragma omp critical 
				{
					MPI_Isend(&vector[1][ny-2-offset][1], 1, xzFacetype, right_neighborY,tag_YR, comm, &reqList[sendcnt]);
					sendcnt++;
				}
			}
		if(id_thread==0)
			if(communicationCnt[4] == 1)
			{
				#pragma omp critical
				{
					MPI_Isend(&vector[1][1][1+offset],1, xyFacetype, left_neighborZ, tag_ZL, comm, &reqList[sendcnt]);
					sendcnt++;
				}
			}
		if(id_thread==1)
			if(communicationCnt[5] == 1)
			{
				#pragma omp critical 
				{
					MPI_Isend(&vector[1][1][nz-2-offset], 1, xyFacetype, right_neighborZ,tag_ZR, comm, &reqList[sendcnt]);
					sendcnt++;
				}
			}
		#pragma omp barrier
		if(id_thread==0)
		{
			assert_eq(recvcnt,sendcnt-recvcnt);
		}
		#pragma omp barrier
		//Buffer swap if any (done before waiting for receiving msg done to delay sync)
		if(id_thread==0)
		{
			if (right_neighborX == myrank &&left_neighborX== myrank)
			{
				for (int iy = 1; iy < ny-1; iy++)
					for (int iz = 1; iz < nz-1; iz++)
					{
						vector[0][iy][iz] = vector[nx-2][iy][iz];
						vector[nx-1][iy][iz] = vector[1][iy][iz];
					}
			}
			if (right_neighborY == myrank &&left_neighborY == myrank)
			{
				for (int ix = 1; ix < nx-1; ix++)
					for (int iz = 1; iz < nz-1; iz++)
					{
						vector[ix][0][iz] = vector[ix][ny-2][iz];
						vector[ix][ny-1][iz] = vector[ix][1][iz];
					}
			}
			if (right_neighborZ == myrank &&left_neighborZ == myrank)
			{
				for (int ix = 1; ix < nx-1; ix++)
					for (int iy = 1; iy < ny-1; iy++)
					{
						vector[ix][iy][0]= vector[ix][iy][nz-2];
						vector[ix][iy][nz-1] = vector[ix][iy][1];
					}
			}
		}
		#pragma omp barrier
		//Need to finish receiving + sending before edge exchange
		if(id_thread==0)
		{
			if(sendcnt>0)
			{
				MPI_Waitall(sendcnt,&reqList[0], &stat[0]);
				bool stopFlag = false;
				for(int si=0;si< sendcnt;si++)
				{
					int error_code = stat[si].MPI_ERROR;
					if (error_code != MPI_SUCCESS)
					{
						stopFlag = true;
						char error_string[100];
						int length_of_error_string, error_class;
						MPI_Error_class(error_code, &error_class);
						MPI_Error_string(error_class, error_string, &length_of_error_string);
						dprintf("MPI_Waitall error at %d%s\n",si, error_string);
					}
				}
				if(stopFlag)
					exit (EXIT_FAILURE);
			}
		}
		#pragma omp barrier
		if( !isFaceOnlyFlag )
		{
			//Exchange yEdge only when Z X neighbours exist
			//if Zleft + Zright, use merged yEdgeType2 to send two edges in one msg
			//Otherwise, only send one yEdge
			#pragma omp barrier
			if(id_thread==0)
			{
				recvcnt = 0, sendcnt = 0;
			}
			#pragma omp barrier
			if(id_thread==0)
				if(communicationCnt[0] == 1)
				{
					if(communicationCnt[4] == 1 && communicationCnt[5] == 1)
					{
						#pragma omp critical
						MPI_Irecv(&vector[0][1][0], 1,yEdgetype2, left_neighborX, tag_XR, comm, &reqList[recvcnt++]);
					}
					else if(communicationCnt[4] == 1)
					{
						#pragma omp critical
						MPI_Irecv(&vector[0][1][0], 1,yEdgetype, left_neighborX, tag_XR, comm, &reqList[recvcnt++]);
					}
					else if(communicationCnt[5] == 1)
					{
						#pragma omp critical
						MPI_Irecv(&vector[0][1][nz-1],1,yEdgetype, left_neighborX, tag_XR, comm, &reqList[recvcnt++]);
					}
				}
			if(id_thread==1)
				if(communicationCnt[1] == 1)
				{
					if(communicationCnt[4] == 1 && communicationCnt[5] == 1)
					{
						#pragma omp critical
						MPI_Irecv(&vector[nx-1][1][0],1,yEdgetype2,right_neighborX,tag_XL, comm, &reqList[recvcnt++]);
					}
					else if(communicationCnt[4] == 1)
					{
						#pragma omp critical
						MPI_Irecv(&vector[nx-1][1][0], 1, yEdgetype, right_neighborX, tag_XL, comm, &reqList[recvcnt++]);
					}
					else if(communicationCnt[5] == 1)
					{
						#pragma omp critical
						MPI_Irecv(&vector[nx-1][1][nz-1],1, yEdgetype, right_neighborX, tag_XL, comm, &reqList[recvcnt++]);
					}
				}
			//Exchange zEdge only when X Y neighbours exist
			//if Xleft + Xright, use merged zEdgeType2
			//Otherwise, only send one zEdge
			if(id_thread==2)
				if(communicationCnt[2] == 1)
				{
					if(communicationCnt[0] == 1 && communicationCnt[1] == 1)
					{
						#pragma omp critical
						MPI_Irecv(&vector[0][0][1], 1, zEdgetype2,left_neighborY, tag_YR,comm, &reqList[recvcnt++]);
					}
					else if(communicationCnt[0] == 1)
					{
						#pragma omp critical
						MPI_Irecv(&vector[0][0][1], 1, zEdgetype, left_neighborY, tag_YR,comm, &reqList[recvcnt++]);
					}
					else if(communicationCnt[1] == 1)
					{
						#pragma omp critical
						MPI_Irecv(&vector[nx-1][0][1],1, zEdgetype, left_neighborY, tag_YR,comm, &reqList[recvcnt++]);
					}
				}
			if(id_thread==3)
				if(communicationCnt[3] == 1)
				{
					if(communicationCnt[0] == 1 && communicationCnt[1] == 1)
					{
						#pragma omp critical
						MPI_Irecv(&vector[0][ny-1][1], 1, zEdgetype2,right_neighborY,tag_YL,comm, &reqList[recvcnt++]);
					}
					else if(communicationCnt[0] == 1)
					{
						#pragma omp critical
						MPI_Irecv(&vector[0][ny-1][1], 1, zEdgetype, right_neighborY,tag_YL,comm, &reqList[recvcnt++]);
					}
					else if(communicationCnt[1] == 1)
					{
						#pragma omp critical
						MPI_Irecv(&vector[nx-1][ny-1][1],1, zEdgetype, right_neighborY,tag_YL,comm, &reqList[recvcnt++]);
					}
				}
			//Exchange xEdge only when Y Z neighbours exist
			//if Yleft + Yright exist, use merged xEdgeType2
			//Otherwise, only send one xEdge
			if(id_thread==0)
				if(communicationCnt[4] == 1)
				{
					if(communicationCnt[2] == 1 && communicationCnt[3] == 1)
					{
						#pragma omp critical
						MPI_Irecv(&vector[1][0][0],1,xEdgetype2,left_neighborZ, tag_ZR, comm, &reqList[recvcnt++]);
					}
					else if(communicationCnt[2] == 1)
					{
						#pragma omp critical
						MPI_Irecv(&vector[1][0][0],1,xEdgetype, left_neighborZ, tag_ZR,comm, &reqList[recvcnt++]);
					}
					else if(communicationCnt[3] == 1)
					{
						#pragma omp critical
						MPI_Irecv(&vector[1][ny-1][0],1, xEdgetype, left_neighborZ,tag_ZR,comm, &reqList[recvcnt++]);
					}
				}
			if(id_thread==1)
				if(communicationCnt[5] == 1)
				{
					if(communicationCnt[2] == 1 && communicationCnt[3] == 1)
					{
						#pragma omp critical
						MPI_Irecv(&vector[1][0][nz-1],1,xEdgetype2,right_neighborZ, tag_ZL,comm, &reqList[recvcnt++]);
					}
					else if(communicationCnt[2] == 1)
					{
						#pragma omp critical
						MPI_Irecv(&vector[1][0][nz-1],1, xEdgetype,right_neighborZ, tag_ZL,comm, &reqList[recvcnt++]);
					}
					else if(communicationCnt[3] == 1)
					{
						#pragma omp critical
						MPI_Irecv(&vector[1][ny-1][nz-1],1,xEdgetype,right_neighborZ, tag_ZL,comm, &reqList[recvcnt++]);
					}
				}
			#pragma omp barrier
			if(id_thread==0)
				sendcnt = recvcnt;
			#pragma omp barrier
			if(id_thread==0)
				if(communicationCnt[0] == 1)
				{
					if(communicationCnt[4] == 1 && communicationCnt[5] == 1)
					{
						#pragma omp critical
						MPI_Isend(&vector[1][1][0], 1,yEdgetype2,left_neighborX, tag_XL, comm, &reqList[sendcnt++]);
					}
					else if(communicationCnt[4] == 1)
					{
						#pragma omp critical
						MPI_Isend(&vector[1][1][0], 1,yEdgetype, left_neighborX, tag_XL, comm, &reqList[sendcnt++]);
					}
					else if(communicationCnt[5] == 1)
					{
						#pragma omp critical
						MPI_Isend(&vector[1][1][nz-1],1,yEdgetype, left_neighborX, tag_XL, comm, &reqList[sendcnt++]);
					}
				}
			if(id_thread==1)
				if(communicationCnt[1] == 1)
				{
					if(communicationCnt[4] == 1 && communicationCnt[5] == 1)
					{
						#pragma omp critical
						MPI_Isend(&vector[nx-2][1][0],1, yEdgetype2,right_neighborX, tag_XR, comm,&reqList[sendcnt++]);
					}
					else if(communicationCnt[4] == 1)
					{
						#pragma omp critical
						MPI_Isend(&vector[nx-2][1][0],1, yEdgetype,right_neighborX, tag_XR, comm, &reqList[sendcnt++]);
					}
					else if(communicationCnt[5] == 1)
					{
						#pragma omp critical
						MPI_Isend(&vector[nx-2][1][nz-1],1,yEdgetype,right_neighborX, tag_XR,comm,&reqList[sendcnt++]);
					}
				}
			if(id_thread==2)
				if(communicationCnt[2] == 1)
				{
					if(communicationCnt[0] == 1 && communicationCnt[1] == 1)
					{
						#pragma omp critical
						MPI_Isend(&vector[0][1][1],1, zEdgetype2, left_neighborY, tag_YL, comm, &reqList[sendcnt++]);
					}
					else if(communicationCnt[0] == 1)
					{
						#pragma omp critical
						MPI_Isend(&vector[0][1][1],1, zEdgetype,left_neighborY, tag_YL, comm, &reqList[sendcnt++]);
					}
					else if(communicationCnt[1] == 1)
					{
						#pragma omp critical
						MPI_Isend(&vector[nx-1][1][1],1,zEdgetype,left_neighborY, tag_YL, comm, &reqList[sendcnt++]);
					}
				}
			if(id_thread==3)
				if(communicationCnt[3] == 1)
				{
					if(communicationCnt[0] == 1 && communicationCnt[1] == 1)
					{
						#pragma omp critical
						MPI_Isend(&vector[0][ny-2][1],1, zEdgetype2, right_neighborY, tag_YR, comm, &reqList[sendcnt++]);
					}
					else if(communicationCnt[0] == 1)
					{
						#pragma omp critical
						MPI_Isend(&vector[0][ny-2][1],1, zEdgetype,right_neighborY, tag_YR, comm, &reqList[sendcnt++]);
					}
					else if(communicationCnt[1] == 1)
					{
						#pragma omp critical
						MPI_Isend(&vector[nx-1][ny-2][1],1,zEdgetype,right_neighborY, tag_YR, comm, &reqList[sendcnt++]);
					}
				}
			if(id_thread==0)
				if(communicationCnt[4] == 1)
				{
					if(communicationCnt[2] == 1 && communicationCnt[3] == 1)
					{
						#pragma omp critical
						MPI_Isend(&vector[1][0][1],1, xEdgetype2, left_neighborZ, tag_ZL, comm, &reqList[sendcnt++]);
					}
					else if(communicationCnt[2] == 1)
					{
						#pragma omp critical
						MPI_Isend(&vector[1][0][1],1, xEdgetype,left_neighborZ, tag_ZL, comm, &reqList[sendcnt++]);
					}
					else if(communicationCnt[3] == 1)
					{
						#pragma omp critical
						MPI_Isend(&vector[1][ny-1][1],1,xEdgetype,left_neighborZ, tag_ZL, comm, &reqList[sendcnt++]);
					}
				}
			if(id_thread==1)
				if(communicationCnt[5] == 1)
				{
					if(communicationCnt[2] == 1 && communicationCnt[3] == 1)
					{
						#pragma omp critical
						MPI_Isend(&vector[1][0][nz-2],1,xEdgetype2, right_neighborZ, tag_ZR, comm, &reqList[sendcnt++]);
					}
					else if(communicationCnt[2] == 1)
					{
						#pragma omp critical
						MPI_Isend(&vector[1][0][nz-2],1,xEdgetype,right_neighborZ, tag_ZR, comm, &reqList[sendcnt++]);
					}
					else if(communicationCnt[3] == 1)
					{
						#pragma omp critical
						MPI_Isend(&vector[1][ny-1][nz-2],1, xEdgetype,right_neighborZ, tag_ZR, comm, &reqList[sendcnt++]);
					}
				}
			#pragma omp barrier
			if(id_thread==0)
				assert_eq(recvcnt,sendcnt-recvcnt);
			#pragma omp barrier
			//Swap Local Edges
			if(id_thread==0)
			{
				if(right_neighborX == myrank &&left_neighborX== myrank)
				{
					if(right_neighborZ != MPI_PROC_NULL )
					{
						for (int iy = 1; iy < ny-1; iy++)
						{
							vector[0][iy][nz-1]= vector[nx-2][iy][nz-1];
							vector[nx-1][iy][nz-1] = vector[1][iy][nz-1];
						}
					}
					if(left_neighborZ != MPI_PROC_NULL )
					{
						for (int iy = 1; iy < ny-1; iy++)
						{
							vector[0][iy][0]= vector[nx-2][iy][0];
							vector[nx-1][iy][0] = vector[1][iy][0];
						}
					}
					if(right_neighborY != MPI_PROC_NULL )
					{
						for (int iz = 1; iz < nz-1; iz++)
						{
							vector[0][ny-1][iz]= vector[nx-2][ny-1][iz];
							vector[nx-1][ny-1][iz] = vector[1][ny-1][iz];
						}
					}
					if(left_neighborY != MPI_PROC_NULL )
					{
						for (int iz = 1; iz < nz-1; iz++)
						{
							vector[0][0][iz] = vector[nx-2][0][iz];
							vector[nx-1][0][iz]= vector[1][0][iz];
						}
					}
				}
				if(right_neighborY == myrank &&left_neighborY == myrank)
				{
					if(right_neighborX != MPI_PROC_NULL )
					{
						for (int iz = 1; iz < nz-1; iz++)
						{
							vector[nx-1][0][iz]= vector[nx-1][ny-2][iz];
							vector[nx-1][ny-1][iz] = vector[nx-1][1][iz];
						}
					}
					if(left_neighborX != MPI_PROC_NULL)
					{
						for (int iz = 1; iz < nz-1; iz++)
						{
							vector[0][0][iz] = vector[0][ny-2][iz];
							vector[0][ny-1][iz]= vector[0][1][iz];
						}
					}
					if(right_neighborZ != MPI_PROC_NULL)
					{
						for (int ix = 1; ix < nx-1; ix++)
						{
							vector[ix][0][nz-1]= vector[ix][ny-2][nz-1];
							vector[ix][ny-1][nz-1] = vector[ix][1][nz-1];
						}
					}
					if(left_neighborZ != MPI_PROC_NULL)
					{
						for (int ix = 1; ix < nx-1; ix++)
						{
							vector[ix][0][0] = vector[ix][ny-2][0];
							vector[ix][ny-1][0]= vector[ix][1][0];
						}
					}
				}
				if(right_neighborZ == myrank &&left_neighborZ == myrank)
				{
					if(right_neighborY != MPI_PROC_NULL )
					{
						for (int ix = 1; ix < nx-1; ix++)
						{
							vector[ix][ny-1][0]= vector[ix][ny-1][nz-2];
							vector[ix][ny-1][nz-1] = vector[ix][ny-1][1];
						}
					}
					if(left_neighborY != MPI_PROC_NULL )
					{
						for (int ix = 1; ix < nx-1; ix++)
						{
							vector[ix][0][0]= vector[ix][0][nz-2];
							vector[ix][0][nz-1] = vector[ix][0][1];
						}
					}
					if(right_neighborX != MPI_PROC_NULL )
					{
						for (int iy = 1; iy < ny-1; iy++)
						{
							vector[nx-1][iy][0]= vector[nx-1][iy][nz-2];
							vector[nx-1][iy][nz-1] = vector[nx-1][iy][1];
						}
					}
					if(left_neighborX != MPI_PROC_NULL )
					{
						for (int iy = 1; iy < ny-1; iy++)
						{
							vector[0][iy][0]= vector[0][iy][nz-2];
							vector[0][iy][nz-1] = vector[0][iy][1];
						}
					}
				}
			}
			//Need to finish receiving edges for corner exchange
			#pragma omp barrier
			if(id_thread==0)
				if(sendcnt>0)
				{
					MPI_Waitall(sendcnt,&reqList[0],&stat[0]);
					bool stopFlag = false;
					for(int si=0;si< sendcnt;si++)
					{
						int error_code = stat[si].MPI_ERROR;
						if (error_code != MPI_SUCCESS)
						{
							stopFlag = true;
							char error_string[100];
							int length_of_error_string, error_class;

							MPI_Error_class(error_code, &error_class);
							MPI_Error_string(error_class, error_string, &length_of_error_string);
							dprintf("MPI_Waitall error at %d%s\n",si, error_string);
						}
					}
					if(stopFlag)
						exit (EXIT_FAILURE);;
				}
			#pragma omp barrier
			//Corner Exchange only needed if XYZ neighbours all exist
			//4 corners communicated in one message
			//Assume Non-periodic will be handled in BC
			//Define corner types for X communication
			if(id_thread==0)
				recvcnt = 0,sendcnt = 0;
			#pragma omp barrier
			if((communicationCnt[2] == 1 || communicationCnt[3] == 1) && (communicationCnt[4] == 1 || communicationCnt[5] == 1))
			{
				//if XLeft exists, send 4 corners to XLeft
				if(id_thread==0)
					if(communicationCnt[0] == 1)
					{
						#pragma omp critical
						MPI_Irecv(&vector[0][0][0],1, cornertype, left_neighborX, tag_XR,comm, &reqList[recvcnt++]);
					}
				//if XRight exist
				if(id_thread==1)
					if(communicationCnt[1] == 1)
					{
						#pragma omp critical
						MPI_Irecv(&vector[nx-1][0][0],1,cornertype, right_neighborX, tag_XL,comm, &reqList[recvcnt++]);
					}
				#pragma omp barrier
				if(id_thread==0)
					sendcnt=recvcnt;
				#pragma omp barrier
				if(id_thread==0)
					if(communicationCnt[0] == 1)
					{
						#pragma omp critical
						MPI_Isend(&vector[1][0][0], 1,cornertype,left_neighborX, tag_XL, comm, &reqList[sendcnt++]);
					}
				if(id_thread==1)
					if(communicationCnt[1] == 1)
					{
						#pragma omp critical
						MPI_Isend(&vector[nx-2][0][0],1,cornertype,right_neighborX, tag_XR,comm, &reqList[sendcnt++]);
					}
			}
			#pragma omp barrier
			if(id_thread==0)
				assert_eq(recvcnt,sendcnt-recvcnt);
			#pragma omp barrier
			//Delay local data copy
			if(id_thread==0)
			{
				if (left_neighborX== myrank && right_neighborX == myrank)
				{
					if( (left_neighborY != MPI_PROC_NULL) && (left_neighborZ != MPI_PROC_NULL))
					{
						vector[0][0][0] = vector[nx-2][0][0];
						vector[nx-1][0][0]= vector[1][0][0];
					}
					if( (left_neighborY != MPI_PROC_NULL) && (right_neighborZ != MPI_PROC_NULL))
					{
						vector[0][0][nz-1]= vector[nx-2][0][nz-1];
						vector[nx-1][0][nz-1] = vector[1][0][nz-1];
					}
					if( (right_neighborY != MPI_PROC_NULL) && (left_neighborZ != MPI_PROC_NULL))
					{
						vector[0][ny-1][0]= vector[nx-2][ny-1][0];
						vector[nx-1][ny-1][0] = vector[1][ny-1][0];
					}
					if( (right_neighborY != MPI_PROC_NULL) && (right_neighborZ != MPI_PROC_NULL))
					{
						vector[0][ny-1][nz-1]= vector[nx-2][ny-1][nz-1];
						vector[nx-1][ny-1][nz-1] = vector[1][ny-1][nz-1];
					}
				}
				else if (left_neighborY== myrank && right_neighborY == myrank)
				{
					if( (left_neighborX != MPI_PROC_NULL) && (left_neighborZ != MPI_PROC_NULL))
					{
						vector[0][0][0] = vector[0][ny-2][0];
						vector[0][ny-1][0]= vector[0][1][0];
					}
					if( (left_neighborX != MPI_PROC_NULL) && (right_neighborZ != MPI_PROC_NULL))
					{
						vector[0][0][nz-1]= vector[0][ny-2][nz-1];
						vector[0][ny-1][nz-1]= vector[0][1][nz-1];
					}
					if( (right_neighborX != MPI_PROC_NULL) && (left_neighborZ != MPI_PROC_NULL))
					{
						vector[nx-1][0][0]= vector[nx-1][ny-2][0];
						vector[nx-1][ny-1][0] = vector[nx-1][1][0];
					}
					if( (right_neighborX != MPI_PROC_NULL) && (right_neighborZ != MPI_PROC_NULL))
					{
						vector[nx-1][0][nz-1]= vector[nx-1][ny-2][nz-1];
						vector[nx-1][ny-1][nz-1] = vector[nx-1][1][nz-1];
					}
				}
				else if (left_neighborZ== myrank && right_neighborZ == myrank)
				{
					if( (left_neighborY != MPI_PROC_NULL) && (left_neighborX != MPI_PROC_NULL))
					{
						vector[0][0][0] = vector[0][0][nz-2];
						vector[0][0][nz-1]= vector[0][0][1];
					}
					if( (left_neighborY != MPI_PROC_NULL) && (right_neighborX != MPI_PROC_NULL))
					{
						vector[nx-1][0][0]= vector[nx-1][0][nz-2];
						vector[nx-1][0][nz-1] = vector[nx-1][0][1];
					}
					if( (right_neighborY != MPI_PROC_NULL) && (left_neighborX != MPI_PROC_NULL))
					{
						vector[0][ny-1][0]= vector[0][ny-1][nz-2];
						vector[0][ny-1][nz-1] = vector[0][ny-1][1];
					}
					if( (right_neighborY != MPI_PROC_NULL) && (right_neighborX != MPI_PROC_NULL))
					{
						vector[nx-1][ny-1][0]= vector[nx-1][ny-1][nz-2];
						vector[nx-1][ny-1][nz-1] = vector[nx-1][ny-1][1];
					}
				}
			}
			if(id_thread==0)
				if(sendcnt>0)
				{
					MPI_Waitall(sendcnt,&reqList[0], &stat[0]);
					bool stopFlag = false;
					for(int si=0;si< sendcnt;si++)
					{
						int error_code = stat[si].MPI_ERROR;
						if (error_code != MPI_SUCCESS)
						{
							stopFlag = true;
							char error_string[100];
							int length_of_error_string, error_class;
							MPI_Error_class(error_code, &error_class);
							MPI_Error_string(error_class, error_string, &length_of_error_string);
							dprintf("MPI_Waitall error at %d%s\n",si, error_string);
						}
					}
					if(stopFlag)
						exit (EXIT_FAILURE);
				}
			#pragma omp barrier
		}
	}
    else if (nthreads==8)
	{
		if(id_thread==0)
			if(left_neighborX != MPI_PROC_NULL && left_neighborX != myrank )
			{
				#pragma omp critical
			{
				MPI_Irecv(&vector[0][1][1], 1, yzFacetype, left_neighborX,tag_XR, comm, &reqList[recvcnt]);
				recvcnt++;
				communicationCnt[0] = 1;
			}
			}
		if(id_thread==1)
			if(right_neighborX != MPI_PROC_NULL && right_neighborX != myrank )
			{
				#pragma omp critical
			{ 
				MPI_Irecv(&vector[nx-1][1][1], 1, yzFacetype, right_neighborX,tag_XL, comm, &reqList[recvcnt]);
				recvcnt++;
				communicationCnt[1] = 1;
				}
			}
		if(id_thread==2)
			if(left_neighborY != MPI_PROC_NULL && left_neighborY != myrank )
			{
				#pragma omp critical 
			{
				MPI_Irecv(&vector[1][0][1], 1, xzFacetype, left_neighborY,tag_YR, comm, &reqList[recvcnt]);
				recvcnt++;
				communicationCnt[2] = 1;
			}
			}
		if(id_thread==3)
			if(right_neighborY != MPI_PROC_NULL && right_neighborY != myrank )
			{
				#pragma omp critical 
			{ 
				MPI_Irecv(&vector[1][ny-1][1], 1, xzFacetype, right_neighborY,tag_YL, comm, &reqList[recvcnt]);
					recvcnt++;
					communicationCnt[3] = 1;
			}
			}
		if(id_thread==4)
			if(left_neighborZ != MPI_PROC_NULL && left_neighborZ != myrank )
			{
				#pragma omp critical
			{ 
				MPI_Irecv(&vector[1][1][0],1, xyFacetype, left_neighborZ,tag_ZR, comm, &reqList[recvcnt]);
				recvcnt++; 
				communicationCnt[4] = 1;
				}
			}
		if(id_thread==5)
			if(right_neighborZ != MPI_PROC_NULL&& right_neighborZ != myrank )
			{
				#pragma omp critical 
			{
				MPI_Irecv(&vector[1][1][nz-1], 1, xyFacetype, right_neighborZ,tag_ZL, comm, &reqList[recvcnt]);
				recvcnt++;
				communicationCnt[5] = 1;
				}
			}
		#pragma omp barrier
		if(id_thread==0)
		{
			sendcnt = recvcnt;
			offset = (isCenterFlag ?0:1);
		}
		#pragma omp barrier
		if(id_thread==0)
			if(communicationCnt[0] == 1)
			{
				#pragma omp critical
				{ 
					MPI_Isend(&vector[1+offset][1][1],1, yzFacetype, left_neighborX, tag_XL, comm, &reqList[sendcnt]);
					sendcnt++;
				}
			}
		if(id_thread==1)
			if(communicationCnt[1] == 1)
			{
				#pragma omp critical 
				{
					MPI_Isend(&vector[nx-2-offset][1][1], 1, yzFacetype, right_neighborX,tag_XR, comm, &reqList[sendcnt]);
					sendcnt++;
				}
			}
		if(id_thread==2)
			if(communicationCnt[2] == 1)
			{
				#pragma omp critical 
				{
					MPI_Isend(&vector[1][1+offset][1],1, xzFacetype, left_neighborY, tag_YL, comm, &reqList[sendcnt]);
					sendcnt++;
				}
			}
		if(id_thread==3)
			if(communicationCnt[3] == 1)
			{
				#pragma omp critical 
				{
					MPI_Isend(&vector[1][ny-2-offset][1], 1, xzFacetype, right_neighborY,tag_YR, comm, &reqList[sendcnt]);
					sendcnt++;
				}
			}
		if(id_thread==4)
			if(communicationCnt[4] == 1)
			{
				#pragma omp critical
				{
					MPI_Isend(&vector[1][1][1+offset],1, xyFacetype, left_neighborZ, tag_ZL, comm, &reqList[sendcnt]);
					sendcnt++;
				}
			}
		if(id_thread==5)
			if(communicationCnt[5] == 1)
			{
				#pragma omp critical 
				{
					MPI_Isend(&vector[1][1][nz-2-offset], 1, xyFacetype, right_neighborZ,tag_ZR, comm, &reqList[sendcnt]);
					sendcnt++;
				}
			}
		#pragma omp barrier
		if(id_thread==0)
		{
			assert_eq(recvcnt,sendcnt-recvcnt);
		}
		#pragma omp barrier
		//Buffer swap if any (done before waiting for receiving msg done to delay sync)
		if(id_thread==0)
		{
			if (right_neighborX == myrank &&left_neighborX== myrank)
			{
				for (int iy = 1; iy < ny-1; iy++)
					for (int iz = 1; iz < nz-1; iz++)
					{
						vector[0][iy][iz] = vector[nx-2][iy][iz];
						vector[nx-1][iy][iz] = vector[1][iy][iz];
					}
			}
			if (right_neighborY == myrank &&left_neighborY == myrank)
			{
				for (int ix = 1; ix < nx-1; ix++)
					for (int iz = 1; iz < nz-1; iz++)
					{
						vector[ix][0][iz] = vector[ix][ny-2][iz];
						vector[ix][ny-1][iz] = vector[ix][1][iz];
					}
			}
			if (right_neighborZ == myrank &&left_neighborZ == myrank)
			{
				for (int ix = 1; ix < nx-1; ix++)
					for (int iy = 1; iy < ny-1; iy++)
					{
						vector[ix][iy][0]= vector[ix][iy][nz-2];
						vector[ix][iy][nz-1] = vector[ix][iy][1];
					}
			}
		}
		#pragma omp barrier
		//Need to finish receiving + sending before edge exchange
		if(id_thread==0)
		{
			if(sendcnt>0)
			{
				MPI_Waitall(sendcnt,&reqList[0], &stat[0]);
				bool stopFlag = false;
				for(int si=0;si< sendcnt;si++)
				{
					int error_code = stat[si].MPI_ERROR;
					if (error_code != MPI_SUCCESS)
					{
						stopFlag = true;
						char error_string[100];
						int length_of_error_string, error_class;
						MPI_Error_class(error_code, &error_class);
						MPI_Error_string(error_class, error_string, &length_of_error_string);
						dprintf("MPI_Waitall error at %d%s\n",si, error_string);
					}
				}
				if(stopFlag)
					exit (EXIT_FAILURE);
			}
		}
		#pragma omp barrier
		if( !isFaceOnlyFlag )
		{
			//Exchange yEdge only when Z X neighbours exist
			//if Zleft + Zright, use merged yEdgeType2 to send two edges in one msg
			//Otherwise, only send one yEdge
			#pragma omp barrier
			if(id_thread==0)
			{
				recvcnt = 0, sendcnt = 0;
			}
			#pragma omp barrier
			if(id_thread==0)
				if(communicationCnt[0] == 1)
				{
					if(communicationCnt[4] == 1 && communicationCnt[5] == 1)
					{
						#pragma omp critical
						MPI_Irecv(&vector[0][1][0], 1,yEdgetype2, left_neighborX, tag_XR, comm, &reqList[recvcnt++]);
					}
					else if(communicationCnt[4] == 1)
					{
						#pragma omp critical
						MPI_Irecv(&vector[0][1][0], 1,yEdgetype, left_neighborX, tag_XR, comm, &reqList[recvcnt++]);
					}
					else if(communicationCnt[5] == 1)
					{
						#pragma omp critical
						MPI_Irecv(&vector[0][1][nz-1],1,yEdgetype, left_neighborX, tag_XR, comm, &reqList[recvcnt++]);
					}
				}
			if(id_thread==1)
				if(communicationCnt[1] == 1)
				{
					if(communicationCnt[4] == 1 && communicationCnt[5] == 1)
					{
						#pragma omp critical
						MPI_Irecv(&vector[nx-1][1][0],1,yEdgetype2,right_neighborX,tag_XL, comm, &reqList[recvcnt++]);
					}
					else if(communicationCnt[4] == 1)
					{
						#pragma omp critical
						MPI_Irecv(&vector[nx-1][1][0], 1, yEdgetype, right_neighborX, tag_XL, comm, &reqList[recvcnt++]);
					}
					else if(communicationCnt[5] == 1)
					{
						#pragma omp critical
						MPI_Irecv(&vector[nx-1][1][nz-1],1, yEdgetype, right_neighborX, tag_XL, comm, &reqList[recvcnt++]);
					}
				}
			//Exchange zEdge only when X Y neighbours exist
			//if Xleft + Xright, use merged zEdgeType2
			//Otherwise, only send one zEdge
			if(id_thread==2)
				if(communicationCnt[2] == 1)
				{
					if(communicationCnt[0] == 1 && communicationCnt[1] == 1)
					{
						#pragma omp critical
						MPI_Irecv(&vector[0][0][1], 1, zEdgetype2,left_neighborY, tag_YR,comm, &reqList[recvcnt++]);
					}
					else if(communicationCnt[0] == 1)
					{
						#pragma omp critical
						MPI_Irecv(&vector[0][0][1], 1, zEdgetype, left_neighborY, tag_YR,comm, &reqList[recvcnt++]);
					}
					else if(communicationCnt[1] == 1)
					{
						#pragma omp critical
						MPI_Irecv(&vector[nx-1][0][1],1, zEdgetype, left_neighborY, tag_YR,comm, &reqList[recvcnt++]);
					}
				}
			if(id_thread==3)
				if(communicationCnt[3] == 1)
				{
					if(communicationCnt[0] == 1 && communicationCnt[1] == 1)
					{
						#pragma omp critical
						MPI_Irecv(&vector[0][ny-1][1], 1, zEdgetype2,right_neighborY,tag_YL,comm, &reqList[recvcnt++]);
					}
					else if(communicationCnt[0] == 1)
					{
						#pragma omp critical
						MPI_Irecv(&vector[0][ny-1][1], 1, zEdgetype, right_neighborY,tag_YL,comm, &reqList[recvcnt++]);
					}
					else if(communicationCnt[1] == 1)
					{
						#pragma omp critical
						MPI_Irecv(&vector[nx-1][ny-1][1],1, zEdgetype, right_neighborY,tag_YL,comm, &reqList[recvcnt++]);
					}
				}
			//Exchange xEdge only when Y Z neighbours exist
			//if Yleft + Yright exist, use merged xEdgeType2
			//Otherwise, only send one xEdge
			if(id_thread==4)
				if(communicationCnt[4] == 1)
				{
					if(communicationCnt[2] == 1 && communicationCnt[3] == 1)
					{
						#pragma omp critical
						MPI_Irecv(&vector[1][0][0],1,xEdgetype2,left_neighborZ, tag_ZR, comm, &reqList[recvcnt++]);
					}
					else if(communicationCnt[2] == 1)
					{
						#pragma omp critical
						MPI_Irecv(&vector[1][0][0],1,xEdgetype, left_neighborZ, tag_ZR,comm, &reqList[recvcnt++]);
					}
					else if(communicationCnt[3] == 1)
					{
						#pragma omp critical
						MPI_Irecv(&vector[1][ny-1][0],1, xEdgetype, left_neighborZ,tag_ZR,comm, &reqList[recvcnt++]);
					}
				}
			if(id_thread==5)
				if(communicationCnt[5] == 1)
				{
					if(communicationCnt[2] == 1 && communicationCnt[3] == 1)
					{
						#pragma omp critical
						MPI_Irecv(&vector[1][0][nz-1],1,xEdgetype2,right_neighborZ, tag_ZL,comm, &reqList[recvcnt++]);
					}
					else if(communicationCnt[2] == 1)
					{
						#pragma omp critical
						MPI_Irecv(&vector[1][0][nz-1],1, xEdgetype,right_neighborZ, tag_ZL,comm, &reqList[recvcnt++]);
					}
					else if(communicationCnt[3] == 1)
					{
						#pragma omp critical
						MPI_Irecv(&vector[1][ny-1][nz-1],1,xEdgetype,right_neighborZ, tag_ZL,comm, &reqList[recvcnt++]);
					}
				}
			#pragma omp barrier
			if(id_thread==0)
				sendcnt = recvcnt;
			#pragma omp barrier
			if(id_thread==0)
				if(communicationCnt[0] == 1)
				{
					if(communicationCnt[4] == 1 && communicationCnt[5] == 1)
					{
						#pragma omp critical
						MPI_Isend(&vector[1][1][0], 1,yEdgetype2,left_neighborX, tag_XL, comm, &reqList[sendcnt++]);
					}
					else if(communicationCnt[4] == 1)
					{
						#pragma omp critical
						MPI_Isend(&vector[1][1][0], 1,yEdgetype, left_neighborX, tag_XL, comm, &reqList[sendcnt++]);
					}
					else if(communicationCnt[5] == 1)
					{
						#pragma omp critical
						MPI_Isend(&vector[1][1][nz-1],1,yEdgetype, left_neighborX, tag_XL, comm, &reqList[sendcnt++]);
					}
				}
			if(id_thread==1)
				if(communicationCnt[1] == 1)
				{
					if(communicationCnt[4] == 1 && communicationCnt[5] == 1)
					{
						#pragma omp critical
						MPI_Isend(&vector[nx-2][1][0],1, yEdgetype2,right_neighborX, tag_XR, comm,&reqList[sendcnt++]);
					}
					else if(communicationCnt[4] == 1)
					{
						#pragma omp critical
						MPI_Isend(&vector[nx-2][1][0],1, yEdgetype,right_neighborX, tag_XR, comm, &reqList[sendcnt++]);
					}
					else if(communicationCnt[5] == 1)
					{
						#pragma omp critical
						MPI_Isend(&vector[nx-2][1][nz-1],1,yEdgetype,right_neighborX, tag_XR,comm,&reqList[sendcnt++]);
					}
				}
			if(id_thread==2)
				if(communicationCnt[2] == 1)
				{
					if(communicationCnt[0] == 1 && communicationCnt[1] == 1)
					{
						#pragma omp critical
						MPI_Isend(&vector[0][1][1],1, zEdgetype2, left_neighborY, tag_YL, comm, &reqList[sendcnt++]);
					}
					else if(communicationCnt[0] == 1)
					{
						#pragma omp critical
						MPI_Isend(&vector[0][1][1],1, zEdgetype,left_neighborY, tag_YL, comm, &reqList[sendcnt++]);
					}
					else if(communicationCnt[1] == 1)
					{
						#pragma omp critical
						MPI_Isend(&vector[nx-1][1][1],1,zEdgetype,left_neighborY, tag_YL, comm, &reqList[sendcnt++]);
					}
				}
			if(id_thread==3)
				if(communicationCnt[3] == 1)
				{
					if(communicationCnt[0] == 1 && communicationCnt[1] == 1)
					{
						#pragma omp critical
						MPI_Isend(&vector[0][ny-2][1],1, zEdgetype2, right_neighborY, tag_YR, comm, &reqList[sendcnt++]);
					}
					else if(communicationCnt[0] == 1)
					{
						#pragma omp critical
						MPI_Isend(&vector[0][ny-2][1],1, zEdgetype,right_neighborY, tag_YR, comm, &reqList[sendcnt++]);
					}
					else if(communicationCnt[1] == 1)
					{
						#pragma omp critical
						MPI_Isend(&vector[nx-1][ny-2][1],1,zEdgetype,right_neighborY, tag_YR, comm, &reqList[sendcnt++]);
					}
				}
			if(id_thread==4)
				if(communicationCnt[4] == 1)
				{
					if(communicationCnt[2] == 1 && communicationCnt[3] == 1)
					{
						#pragma omp critical
						MPI_Isend(&vector[1][0][1],1, xEdgetype2, left_neighborZ, tag_ZL, comm, &reqList[sendcnt++]);
					}
					else if(communicationCnt[2] == 1)
					{
						#pragma omp critical
						MPI_Isend(&vector[1][0][1],1, xEdgetype,left_neighborZ, tag_ZL, comm, &reqList[sendcnt++]);
					}
					else if(communicationCnt[3] == 1)
					{
						#pragma omp critical
						MPI_Isend(&vector[1][ny-1][1],1,xEdgetype,left_neighborZ, tag_ZL, comm, &reqList[sendcnt++]);
					}
				}
			if(id_thread==5)
				if(communicationCnt[5] == 1)
				{
					if(communicationCnt[2] == 1 && communicationCnt[3] == 1)
					{
						#pragma omp critical
						MPI_Isend(&vector[1][0][nz-2],1,xEdgetype2, right_neighborZ, tag_ZR, comm, &reqList[sendcnt++]);
					}
					else if(communicationCnt[2] == 1)
					{
						#pragma omp critical
						MPI_Isend(&vector[1][0][nz-2],1,xEdgetype,right_neighborZ, tag_ZR, comm, &reqList[sendcnt++]);
					}
					else if(communicationCnt[3] == 1)
					{
						#pragma omp critical
						MPI_Isend(&vector[1][ny-1][nz-2],1, xEdgetype,right_neighborZ, tag_ZR, comm, &reqList[sendcnt++]);
					}
				}
			#pragma omp barrier
			if(id_thread==0)
				assert_eq(recvcnt,sendcnt-recvcnt);
			#pragma omp barrier
			//Swap Local Edges
			if(id_thread==0)
			{
				if(right_neighborX == myrank &&left_neighborX== myrank)
				{
					if(right_neighborZ != MPI_PROC_NULL )
					{
						for (int iy = 1; iy < ny-1; iy++)
						{
							vector[0][iy][nz-1]= vector[nx-2][iy][nz-1];
							vector[nx-1][iy][nz-1] = vector[1][iy][nz-1];
						}
					}
					if(left_neighborZ != MPI_PROC_NULL )
					{
						for (int iy = 1; iy < ny-1; iy++)
						{
							vector[0][iy][0]= vector[nx-2][iy][0];
							vector[nx-1][iy][0] = vector[1][iy][0];
						}
					}
					if(right_neighborY != MPI_PROC_NULL )
					{
						for (int iz = 1; iz < nz-1; iz++)
						{
							vector[0][ny-1][iz]= vector[nx-2][ny-1][iz];
							vector[nx-1][ny-1][iz] = vector[1][ny-1][iz];
						}
					}
					if(left_neighborY != MPI_PROC_NULL )
					{
						for (int iz = 1; iz < nz-1; iz++)
						{
							vector[0][0][iz] = vector[nx-2][0][iz];
							vector[nx-1][0][iz]= vector[1][0][iz];
						}
					}
				}
				if(right_neighborY == myrank &&left_neighborY == myrank)
				{
					if(right_neighborX != MPI_PROC_NULL )
					{
						for (int iz = 1; iz < nz-1; iz++)
						{
							vector[nx-1][0][iz]= vector[nx-1][ny-2][iz];
							vector[nx-1][ny-1][iz] = vector[nx-1][1][iz];
						}
					}
					if(left_neighborX != MPI_PROC_NULL)
					{
						for (int iz = 1; iz < nz-1; iz++)
						{
							vector[0][0][iz] = vector[0][ny-2][iz];
							vector[0][ny-1][iz]= vector[0][1][iz];
						}
					}
					if(right_neighborZ != MPI_PROC_NULL)
					{
						for (int ix = 1; ix < nx-1; ix++)
						{
							vector[ix][0][nz-1]= vector[ix][ny-2][nz-1];
							vector[ix][ny-1][nz-1] = vector[ix][1][nz-1];
						}
					}
					if(left_neighborZ != MPI_PROC_NULL)
					{
						for (int ix = 1; ix < nx-1; ix++)
						{
							vector[ix][0][0] = vector[ix][ny-2][0];
							vector[ix][ny-1][0]= vector[ix][1][0];
						}
					}
				}
				if(right_neighborZ == myrank &&left_neighborZ == myrank)
				{
					if(right_neighborY != MPI_PROC_NULL )
					{
						for (int ix = 1; ix < nx-1; ix++)
						{
							vector[ix][ny-1][0]= vector[ix][ny-1][nz-2];
							vector[ix][ny-1][nz-1] = vector[ix][ny-1][1];
						}
					}
					if(left_neighborY != MPI_PROC_NULL )
					{
						for (int ix = 1; ix < nx-1; ix++)
						{
							vector[ix][0][0]= vector[ix][0][nz-2];
							vector[ix][0][nz-1] = vector[ix][0][1];
						}
					}
					if(right_neighborX != MPI_PROC_NULL )
					{
						for (int iy = 1; iy < ny-1; iy++)
						{
							vector[nx-1][iy][0]= vector[nx-1][iy][nz-2];
							vector[nx-1][iy][nz-1] = vector[nx-1][iy][1];
						}
					}
					if(left_neighborX != MPI_PROC_NULL )
					{
						for (int iy = 1; iy < ny-1; iy++)
						{
							vector[0][iy][0]= vector[0][iy][nz-2];
							vector[0][iy][nz-1] = vector[0][iy][1];
						}
					}
				}
			}
			//Need to finish receiving edges for corner exchange
			#pragma omp barrier
			if(id_thread==0)
				if(sendcnt>0)
				{
					MPI_Waitall(sendcnt,&reqList[0],&stat[0]);
					bool stopFlag = false;
					for(int si=0;si< sendcnt;si++)
					{
						int error_code = stat[si].MPI_ERROR;
						if (error_code != MPI_SUCCESS)
						{
							stopFlag = true;
							char error_string[100];
							int length_of_error_string, error_class;

							MPI_Error_class(error_code, &error_class);
							MPI_Error_string(error_class, error_string, &length_of_error_string);
							dprintf("MPI_Waitall error at %d%s\n",si, error_string);
						}
					}
					if(stopFlag)
						exit (EXIT_FAILURE);;
				}
			#pragma omp barrier
			//Corner Exchange only needed if XYZ neighbours all exist
			//4 corners communicated in one message
			//Assume Non-periodic will be handled in BC
			//Define corner types for X communication
			if(id_thread==0)
				recvcnt = 0,sendcnt = 0;
			#pragma omp barrier
			if((communicationCnt[2] == 1 || communicationCnt[3] == 1) && (communicationCnt[4] == 1 || communicationCnt[5] == 1))
			{
				//if XLeft exists, send 4 corners to XLeft
				if(id_thread==0)
					if(communicationCnt[0] == 1)
					{
						#pragma omp critical
						MPI_Irecv(&vector[0][0][0],1, cornertype, left_neighborX, tag_XR,comm, &reqList[recvcnt++]);
					}
				//if XRight exist
				if(id_thread==1)
					if(communicationCnt[1] == 1)
					{
						#pragma omp critical
						MPI_Irecv(&vector[nx-1][0][0],1,cornertype, right_neighborX, tag_XL,comm, &reqList[recvcnt++]);
					}
				#pragma omp barrier
				if(id_thread==0)
					sendcnt=recvcnt;
				#pragma omp barrier
				if(id_thread==0)
					if(communicationCnt[0] == 1)
					{
						#pragma omp critical
						MPI_Isend(&vector[1][0][0], 1,cornertype,left_neighborX, tag_XL, comm, &reqList[sendcnt++]);
					}
				if(id_thread==1)
					if(communicationCnt[1] == 1)
					{
						#pragma omp critical
						MPI_Isend(&vector[nx-2][0][0],1,cornertype,right_neighborX, tag_XR,comm, &reqList[sendcnt++]);
					}
			}
			#pragma omp barrier
			if(id_thread==0)
				assert_eq(recvcnt,sendcnt-recvcnt);
			#pragma omp barrier
			//Delay local data copy
			if(id_thread==0)
			{
				if (left_neighborX== myrank && right_neighborX == myrank)
				{
					if( (left_neighborY != MPI_PROC_NULL) && (left_neighborZ != MPI_PROC_NULL))
					{
						vector[0][0][0] = vector[nx-2][0][0];
						vector[nx-1][0][0]= vector[1][0][0];
					}
					if( (left_neighborY != MPI_PROC_NULL) && (right_neighborZ != MPI_PROC_NULL))
					{
						vector[0][0][nz-1]= vector[nx-2][0][nz-1];
						vector[nx-1][0][nz-1] = vector[1][0][nz-1];
					}
					if( (right_neighborY != MPI_PROC_NULL) && (left_neighborZ != MPI_PROC_NULL))
					{
						vector[0][ny-1][0]= vector[nx-2][ny-1][0];
						vector[nx-1][ny-1][0] = vector[1][ny-1][0];
					}
					if( (right_neighborY != MPI_PROC_NULL) && (right_neighborZ != MPI_PROC_NULL))
					{
						vector[0][ny-1][nz-1]= vector[nx-2][ny-1][nz-1];
						vector[nx-1][ny-1][nz-1] = vector[1][ny-1][nz-1];
					}
				}
				else if (left_neighborY== myrank && right_neighborY == myrank)
				{
					if( (left_neighborX != MPI_PROC_NULL) && (left_neighborZ != MPI_PROC_NULL))
					{
						vector[0][0][0] = vector[0][ny-2][0];
						vector[0][ny-1][0]= vector[0][1][0];
					}
					if( (left_neighborX != MPI_PROC_NULL) && (right_neighborZ != MPI_PROC_NULL))
					{
						vector[0][0][nz-1]= vector[0][ny-2][nz-1];
						vector[0][ny-1][nz-1]= vector[0][1][nz-1];
					}
					if( (right_neighborX != MPI_PROC_NULL) && (left_neighborZ != MPI_PROC_NULL))
					{
						vector[nx-1][0][0]= vector[nx-1][ny-2][0];
						vector[nx-1][ny-1][0] = vector[nx-1][1][0];
					}
					if( (right_neighborX != MPI_PROC_NULL) && (right_neighborZ != MPI_PROC_NULL))
					{
						vector[nx-1][0][nz-1]= vector[nx-1][ny-2][nz-1];
						vector[nx-1][ny-1][nz-1] = vector[nx-1][1][nz-1];
					}
				}
				else if (left_neighborZ== myrank && right_neighborZ == myrank)
				{
					if( (left_neighborY != MPI_PROC_NULL) && (left_neighborX != MPI_PROC_NULL))
					{
						vector[0][0][0] = vector[0][0][nz-2];
						vector[0][0][nz-1]= vector[0][0][1];
					}
					if( (left_neighborY != MPI_PROC_NULL) && (right_neighborX != MPI_PROC_NULL))
					{
						vector[nx-1][0][0]= vector[nx-1][0][nz-2];
						vector[nx-1][0][nz-1] = vector[nx-1][0][1];
					}
					if( (right_neighborY != MPI_PROC_NULL) && (left_neighborX != MPI_PROC_NULL))
					{
						vector[0][ny-1][0]= vector[0][ny-1][nz-2];
						vector[0][ny-1][nz-1] = vector[0][ny-1][1];
					}
					if( (right_neighborY != MPI_PROC_NULL) && (right_neighborX != MPI_PROC_NULL))
					{
						vector[nx-1][ny-1][0]= vector[nx-1][ny-1][nz-2];
						vector[nx-1][ny-1][nz-1] = vector[nx-1][ny-1][1];
					}
				}
			}
			if(id_thread==0)
				if(sendcnt>0)
				{
					MPI_Waitall(sendcnt,&reqList[0], &stat[0]);
					bool stopFlag = false;
					for(int si=0;si< sendcnt;si++)
					{
						int error_code = stat[si].MPI_ERROR;
						if (error_code != MPI_SUCCESS)
						{
							stopFlag = true;
							char error_string[100];
							int length_of_error_string, error_class;
							MPI_Error_class(error_code, &error_class);
							MPI_Error_string(error_class, error_string, &length_of_error_string);
							dprintf("MPI_Waitall error at %d%s\n",si, error_string);
						}
					}
					if(stopFlag)
						exit (EXIT_FAILURE);
				}
			#pragma omp barrier
		}
	}
}

//if this is NodeInterpolation operation
	if(needInterp)
	{
		addFace(nx, ny, nz, vector,  vct);
		addEdgeZ(nx, ny, nz, vector, vct);
		//      if(id_thread==0)
		addEdgeY(nx, ny, nz, vector, vct);
		//      if(id_thread==0)
		addEdgeX(nx, ny, nz, vector, vct);
		//      if(id_thread==0)
		addCorner(nx, ny, nz, vector, vct);
	}
}


void communicateNodeBC(int nx, int ny, int nz, arr3_double _vector,
		       int bcFaceXrght, int bcFaceXleft,
		       int bcFaceYrght, int bcFaceYleft,
		       int bcFaceZrght, int bcFaceZleft,
		       const VirtualTopology3D * vct, EMfields3D *EMf)
{
  double ***vector=_vector.fetch_arr3();
  NBDerivedHaloComm(nx, ny, nz, vector, vct, EMf, false,false,false,false);

  // ////////////////////////////////////////////////////////////////////////
  // ///////////////// APPLY the boundary conditions ////////////////////////
  // ////////////////////////////////////////////////////////////////////////
  BCface(nx, ny, nz, vector, bcFaceXrght, bcFaceXleft, bcFaceYrght, bcFaceYleft, bcFaceZrght, bcFaceZleft, vct);
}

void communicateNodeBoxStencilBC( int nx, int ny, int nz, arr3_double _vector,
				  int bcFaceXrght, int bcFaceXleft,
				  int bcFaceYrght, int bcFaceYleft,
				  int bcFaceZrght, int bcFaceZleft,
				  const VirtualTopology3D * vct, EMfields3D *EMf)
{

  double ***vector=_vector.fetch_arr3();

  NBDerivedHaloComm(nx, ny, nz, vector, vct, EMf, false,true,false,false);

  // ////////////////////////////////////////////////////////////////////////
  // ///////////////// APPLY the boundary conditions ////////////////////////
  // ////////////////////////////////////////////////////////////////////////
  BCface(nx, ny, nz, vector, bcFaceXrght, bcFaceXleft, bcFaceYrght, bcFaceYleft, bcFaceZrght, bcFaceZleft, vct);

}

void communicateNodeBoxStencilBC_P(int nx, int ny, int nz, arr3_double _vector,
				   int bcFaceXrght, int bcFaceXleft,
				   int bcFaceYrght, int bcFaceYleft,
				   int bcFaceZrght, int bcFaceZleft,
				   const VirtualTopology3D * vct, EMfields3D *EMf)
{

  double ***vector=_vector.fetch_arr3();

  NBDerivedHaloComm(nx, ny, nz, vector, vct, EMf, false,true,false,true);

  // ////////////////////////////////////////////////////////////////////////
  // ///////////////// APPLY the boundary conditions ////////////////////////
  // ////////////////////////////////////////////////////////////////////////
  BCface_P(nx, ny, nz, vector, bcFaceXrght, bcFaceXleft, bcFaceYrght, bcFaceYleft, bcFaceZrght, bcFaceZleft, vct);

}

void communicateNodeBC_P(int nx, int ny, int nz, arr3_double _vector,
			 int bcFaceXrght, int bcFaceXleft,
			 int bcFaceYrght, int bcFaceYleft,
			 int bcFaceZrght, int bcFaceZleft,
			 const VirtualTopology3D * vct, EMfields3D *EMf)
{
  double ***vector=_vector.fetch_arr3();
  NBDerivedHaloComm(nx, ny, nz, vector, vct, EMf, false,false,false,true);

  // ////////////////////////////////////////////////////////////////////////
  // ///////////////// APPLY the boundary conditions ////////////////////////
  // ////////////////////////////////////////////////////////////////////////
  BCface_P(nx, ny, nz, vector, bcFaceXrght, bcFaceXleft, bcFaceYrght, bcFaceYleft, bcFaceZrght, bcFaceZleft, vct);
}

void communicateNode_P(int nx, int ny, int nz, double*** vector,
		       const VirtualTopology3D * vct, EMfields3D *EMf)
{
  NBDerivedHaloComm(nx, ny, nz, vector, vct, EMf, false,false,false,true);
}


void communicateCenterBC(int nx, int ny, int nz, arr3_double _vector,
			 int bcFaceXrght, int bcFaceXleft,
			 int bcFaceYrght, int bcFaceYleft,
			 int bcFaceZrght, int bcFaceZleft,
			 const VirtualTopology3D * vct, EMfields3D *EMf)
{
  double ***vector=_vector.fetch_arr3();
  NBDerivedHaloComm(nx, ny, nz, vector, vct, EMf, true, false,false,false);

  // ////////////////////////////////////////////////////////////////////////
  // ///////////////// APPLY the boundary conditions ////////////////////////
  // ////////////////////////////////////////////////////////////////////////
  BCface(nx, ny, nz, vector, bcFaceXrght, bcFaceXleft, bcFaceYrght, bcFaceYleft, bcFaceZrght, bcFaceZleft, vct);

}

void communicateCenterBC_P(   int nx, int ny, int nz, arr3_double _vector,
			      int bcFaceXrght, int bcFaceXleft,
			      int bcFaceYrght, int bcFaceYleft,
			      int bcFaceZrght, int bcFaceZleft,
			      const VirtualTopology3D * vct, EMfields3D *EMf)
{

  double ***vector=_vector.fetch_arr3();
  NBDerivedHaloComm(nx, ny, nz, vector, vct, EMf, true, false,false,true);

  // ////////////////////////////////////////////////////////////////////////
  // ///////////////// APPLY the boundary conditions ////////////////////////
  // ////////////////////////////////////////////////////////////////////////
  BCface_P(nx, ny, nz, vector, bcFaceXrght, bcFaceXleft, bcFaceYrght, bcFaceYleft, bcFaceZrght, bcFaceZleft, vct);
}

void communicateCenterBoxStencilBC(   int nx, int ny, int nz, arr3_double _vector,
				      int bcFaceXrght, int bcFaceXleft,
				      int bcFaceYrght, int bcFaceYleft,
				      int bcFaceZrght, int bcFaceZleft,
				      const VirtualTopology3D * vct, EMfields3D *EMf)
{
  double ***vector=_vector.fetch_arr3();
  NBDerivedHaloComm(nx, ny, nz, vector, vct, EMf, true,true,false,false);

  // ////////////////////////////////////////////////////////////////////////
  // ///////////////// APPLY the boundary conditions ////////////////////////
  // ////////////////////////////////////////////////////////////////////////
  BCface(nx, ny, nz, vector, bcFaceXrght, bcFaceXleft, bcFaceYrght, bcFaceYleft, bcFaceZrght, bcFaceZleft, vct);

}

void communicateCenterBoxStencilBC_P( int nx, int ny, int nz, arr3_double _vector,
				      int bcFaceXrght, int bcFaceXleft,
				      int bcFaceYrght, int bcFaceYleft,
				      int bcFaceZrght, int bcFaceZleft,
				      const VirtualTopology3D * vct, EMfields3D *EMf)
{

  double ***vector=_vector.fetch_arr3();
  NBDerivedHaloComm(nx, ny, nz, vector, vct, EMf, true,true,false,true);

  // ////////////////////////////////////////////////////////////////////////
  // ///////////////// APPLY the boundary conditions ////////////////////////
  // ////////////////////////////////////////////////////////////////////////
  BCface_P(nx, ny, nz, vector, bcFaceXrght, bcFaceXleft, bcFaceYrght, bcFaceYleft, bcFaceZrght, bcFaceZleft, vct);

}

/** add the values of ghost cells faces to the 3D physical vector */
void addFace(int nx, int ny, int nz, double ***vector, const VirtualTopology3D * vct)
{
  const int nxr = nx-2;
  const int nyr = ny-2;
  const int nzr = nz-2;

  // Xright
  if (vct->hasXrghtNeighbor_P())
    {
      //#pragma omp for
      for (int j = 1; j <= nyr; j++)
	for (int k = 1; k <= nzr; k++)
	  vector[nx - 2][j][k] += vector[nx - 1][j][k];
    }
  // XLEFT
  if (vct->hasXleftNeighbor_P())
    {
      //#pragma omp for
      for (int j = 1; j <= nyr; j++)
	for (int k = 1; k <= nzr; k++)
	  vector[1][j][k] += vector[0][j][k];
    }

  // Yright
  if (vct->hasYrghtNeighbor_P())
    {
      //#pragma omp for
      for (int i = 1; i <= nxr; i++)
	for (int k = 1; k <= nzr; k++)
	  vector[i][ny - 2][k] += vector[i][ny - 1][k];
    }
  // Yleft
  if (vct->hasYleftNeighbor_P())
    {
      //#pragma omp for
      for (int i = 1; i <= nxr; i++)
	for (int k = 1; k <= nzr; k++)
	  vector[i][1][k] += vector[i][0][k];
    }
  // Zright
  if (vct->hasZrghtNeighbor_P())
    {
      //    #pragma omp for
      for (int i = 1; i <= nxr; i++)
	for (int j = 1; j <= nyr; j++)
	  vector[i][j][nz - 2] += vector[i][j][nz - 1];
    }
  // ZLEFT
  if (vct->hasZleftNeighbor_P())
    {
      //    #pragma omp for
      for (int i = 1; i <= nxr; i++)
	for (int j = 1; j <= nyr; j++)
	  vector[i][j][1] += vector[i][j][0];
    }
}
/** insert the ghost cells Edge Z in the 3D physical vector */
void addEdgeZ(int nx, int ny, int nz, double ***vector, const VirtualTopology3D * vct)
{
  if (vct->hasXrghtNeighbor_P() && vct->hasYrghtNeighbor_P()) {
    for (int i = 1; i < (nz - 1); i++)
      vector[nx - 2][ny - 2][i] += vector[nx - 1][ny - 1][i];
  }
  if (vct->hasXleftNeighbor_P() && vct->hasYleftNeighbor_P()) {
    for (int i = 1; i < (nz - 1); i++)
      vector[1][1][i] += vector[0][0][i];
  }
  if (vct->hasXrghtNeighbor_P() && vct->hasYleftNeighbor_P()) {
    for (int i = 1; i < (nz - 1); i++)
      vector[nx - 2][1][i] += vector[nx - 1][0][i];
  }
  if (vct->hasXleftNeighbor_P() && vct->hasYrghtNeighbor_P()) {
    for (int i = 1; i < (nz - 1); i++)
      vector[1][ny - 2][i] += vector[0][ny - 1][i];
  }
}
/** add the ghost cell values Edge Y to the 3D physical vector */
void addEdgeY(int nx, int ny, int nz, double ***vector, const VirtualTopology3D * vct)
{
  if (vct->hasXrghtNeighbor_P() && vct->hasZrghtNeighbor_P()) {
    for (int i = 1; i < (ny - 1); i++)
      vector[nx - 2][i][nz - 2] += vector[nx - 1][i][nz - 1];
  }
  if (vct->hasXleftNeighbor_P() && vct->hasZleftNeighbor_P()) {
    for (int i = 1; i < (ny - 1); i++)
      vector[1][i][1] += vector[0][i][0];
  }
  if (vct->hasXleftNeighbor_P() && vct->hasZrghtNeighbor_P()) {
    for (int i = 1; i < (ny - 1); i++)
      vector[1][i][nz - 2] += vector[0][i][nz - 1];
  }
  if (vct->hasXrghtNeighbor_P() && vct->hasZleftNeighbor_P()) {
    for (int i = 1; i < (ny - 1); i++)
      vector[nx - 2][i][1] += vector[nx - 1][i][0];
  }
}
/** add the ghost values Edge X to the 3D physical vector */
void addEdgeX(int nx, int ny, int nz, double ***vector, const VirtualTopology3D * vct)
{
  if (vct->hasYrghtNeighbor_P() && vct->hasZrghtNeighbor_P()) {
    for (int i = 1; i < (nx - 1); i++)
      vector[i][ny - 2][nz - 2] += vector[i][ny - 1][nz - 1];
  }
  if (vct->hasYleftNeighbor_P() && vct->hasZleftNeighbor_P()) {
    for (int i = 1; i < (nx - 1); i++)
      vector[i][1][1] += vector[i][0][0];
  }
  if (vct->hasYleftNeighbor_P() && vct->hasZrghtNeighbor_P()) {
    for (int i = 1; i < (nx - 1); i++)
      vector[i][1][nz - 2] += vector[i][0][nz - 1];
  }
  if (vct->hasYrghtNeighbor_P() && vct->hasZleftNeighbor_P()) {
    for (int i = 1; i < (nx - 1); i++)
      vector[i][ny - 2][1] += vector[i][ny - 1][0];
  }
}
/** add ghost cells values Corners in the 3D physical vector */
void addCorner(int nx, int ny, int nz, double ***vector,const VirtualTopology3D * vct)
{
  if (vct->hasXrghtNeighbor_P() && vct->hasYrghtNeighbor_P() && vct->hasZrghtNeighbor_P())
    vector[nx - 2][ny - 2][nz - 2] += vector[nx - 1][ny - 1][nz - 1];
  if (vct->hasXleftNeighbor_P() && vct->hasYrghtNeighbor_P() && vct->hasZrghtNeighbor_P())
    vector[1][ny - 2][nz - 2] += vector[0][ny - 1][nz - 1];
  if (vct->hasXrghtNeighbor_P() && vct->hasYleftNeighbor_P() && vct->hasZrghtNeighbor_P())
    vector[nx - 2][1][nz - 2] += vector[nx - 1][0][nz - 1];
  if (vct->hasXleftNeighbor_P() && vct->hasYleftNeighbor_P() && vct->hasZrghtNeighbor_P())
    vector[1][1][nz - 2] += vector[0][0][nz - 1];
  if (vct->hasXrghtNeighbor_P() && vct->hasYrghtNeighbor_P() && vct->hasZleftNeighbor_P())
    vector[nx - 2][ny - 2][1] += vector[nx - 1][ny - 1][0];
  if (vct->hasXleftNeighbor_P() && vct->hasYrghtNeighbor_P() && vct->hasZleftNeighbor_P())
    vector[1][ny - 2][1] += vector[0][ny - 1][0] ;
  if (vct->hasXrghtNeighbor_P() && vct->hasYleftNeighbor_P() && vct->hasZleftNeighbor_P())
    vector[nx - 2][1][1] += vector[nx - 1][0][0] ;
  if (vct->hasXleftNeighbor_P() && vct->hasYleftNeighbor_P() && vct->hasZleftNeighbor_P())
    vector[1][1][1] += vector[0][0][0];

}
/** communicate and sum shared ghost cells */
void communicateInterp(int nx, int ny, int nz, double*** vector, const VirtualTopology3D * vct, EMfields3D *EMf)
{
  NBDerivedHaloComm(nx, ny, nz, vector, vct, EMf, true,false,true,true);
}