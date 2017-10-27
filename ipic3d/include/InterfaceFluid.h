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

/*******************************************************
InterfaceFluid_H handels all interaction between
a fluid model (BATSRUS/SWMF) and the particle model
in iPic3D

Writen by Lars Daldorff (daldorff@umich.edu) 15 Jan 2013

********************************************************/


#ifndef InterfaceFluid_H
#define InterfaceFluid_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <stdlib.h>
#include <string>
#include <unistd.h>


#include "Alloc.h"
#include "VCtopology3D.h"

using namespace std;

class InterfaceFluid
{ 
private:
  
  // BATSRUS output header file name 
  string FilenameRoot;
  // variables to read data from header file
  string TempString;
  string *INVarName_V;
  
  int INnProc;      // number of processors in the fluid code
  int INnStep;      // iteration number
  int INnPlotvar;   // number of variables in the fluid file
  int INnEqpar;     // number of equation parameters
  int INnRoot_D[3]; // number of root blocks
  int INnIJK_D[3];  // number of cells per block
  double INtime,INxRange_I[6],INplotDx_D[3],INgridDx_D[3];

  int nIJK_D[3];    // number of cells in the passed uniform grid
  double XYZ_D[3];  // physical size of the passed domain

  double ****State_GV; // cell centered state variables

  string FluidModel;  // Type of fluid model

  int iRho,iUx,iUy,iUz,iBx,iBy,iBz,iP; // named indexes for State_GV

  int StartIdx_D[3];  // storage for starting grid indexes of this processor

  static const int NG = 1;

  // Eventually we will interpolate from FLUID grid to iPIC3D grid based on coordinates
  inline void getGlobalIndex(const int il, const int jl, const int kl,
                             int *ig, int *jg, int *kg)
  {
    // convert from local cpu index il to global domain index ig
    *ig =StartIdx_D[0] + il;

    // solution is constant across ignored dimensions indicated by nIJK_D == 1
    if(nIJK_D[1] == 1) *jg = 0; else *jg =StartIdx_D[1] + jl;
    if(nIJK_D[2] == 1) *kg = 0; else *kg =StartIdx_D[2] + kl;

  };

  void InitDataStorage()
  {
    for(int i =0; i<3; i++){
      // set domain size for dimension i
      XYZ_D[i] = (INxRange_I[2*i+1] - INxRange_I[2*i]);
      // set number of cells in dimension i
      nIJK_D[i] = floor(XYZ_D[i]/INgridDx_D[i] + 0.5);
      cout<<"XYZ_D[i] = "<<XYZ_D[i]<<",  nIJK_D[i] = "<<nIJK_D[i]<<endl;	
      if(nIJK_D[i] == 0) 
	{
	  nIJK_D[i] = 1;
          XYZ_D[i] = INgridDx_D[1];
	}
    }
    // Allocate state variable array [nI,nJ,nK,nVar], last index fastest
    State_GV = newArr4(double,nIJK_D[0],nIJK_D[1],nIJK_D[2],INnPlotvar);
    //cout.flush();
  };		
  
  void ReadVariables()
  {
    // read fluid data from files

    // construct file name
    string Filename = FilenameRoot;
    string ProcName; // we will add "pe0000" for the filename
    ostringstream ss; // helping variable for constuctiong "pe0000"
    ifstream myReadFile;
    
    // variables for geting index from position
    double dx,xyz_D[3];
    int ijk_D[3];

    // loop through files generated by the processors used by the fluid code
    for(int iProc=0;iProc<INnProc;iProc++)
      {
	ss<< setfill('0') << setw(4) << iProc; // use 4 spaces and fill with "0"
	Filename = FilenameRoot + "_pe" + ss.str() + ".idl";
	ss.str(""); // reset string
	ss.clear(); // clear any error messagnes
	cout<<" opening file : "<<Filename<<endl;
	
        myReadFile.open(Filename.c_str(),ifstream::in|ifstream::out);
        if(myReadFile.is_open())
	  {
	    while (!myReadFile.eof())
	      {
		// read cell size and cell center coordinates
		myReadFile >> dx >> xyz_D[0] >> xyz_D[1] >> xyz_D[2];
		// get cell index in the passed fluid grid
                for(int iDim=0;iDim<3;iDim++)
		  {
		    ijk_D[iDim]=0;	
		    if(nIJK_D[iDim] > 1 ) ijk_D[iDim] = floor((xyz_D[iDim] - INxRange_I[iDim*2]-0.5*INgridDx_D[iDim])/INgridDx_D[iDim] + 0.5);
		  }
		// read cell state variables
		for(int iVar=0; iVar<INnPlotvar; iVar++)
		  { 
		    myReadFile >> State_GV[ijk_D[0]][ijk_D[1]][ijk_D[2]][iVar]; 
		  }
	      }
	  }
	else
	  {
	    cout<<"Can not open "<<Filename<<", will abort \n"<<flush;
	    abort();
	  }
	cout.flush();
	myReadFile.close();
      }
    //PrintStateVar();
  };

  // for debugging
  void PrintStateVar()
  {
    cout.setf(ios::fixed,ios::floatfield);
    cout.precision(4);
    cout.width(10);
    for(int iVar=0;iVar<INnPlotvar;iVar++)
      {
	cout<<endl;
	cout<<" Variable : "<<INVarName_V[iVar]<<endl;
	for(int k=0;k<nIJK_D[2];k++)
	  {
	    for(int j=0;j<nIJK_D[1];j++)
	      { 
		for(int i=0;i<nIJK_D[0];i++)
		  { 
		    cout<<setw(10)<<State_GV[i][j][k][iVar]<<", ";
		  }
		cout<<endl;
	      }	
	    cout<<endl;
	  }
      }
  };
  
 public:
   
  // constructor
  InterfaceFluid(){
    InitFluid();
  };

  void InitFluid()
  {
    cout<<"Constructor::InterfaceFluid"<<endl;cout.flush();
    FilenameRoot="./TestFluid/Input/cut_var_2_t0000.000_n0000000";
    FluidModel = "mhd";
    iRho = 0;
    iUx  = 1;
    iUy  = 2;
    iUz  = 3;
    iBx  = 4;
    iBy  = 5;
    iBz  = 6;
    iP   = 7;

    // read the BATSRUS header file that contains the grid information
    ReadIdlAscii();
  };

  void setGlobalStartIndex(VCtopology3D *vct)
  {
         StartIdx_D[0] = vct->getCoordinates(0)*getFluidNxc()/(double)vct->getXLEN(); 
         StartIdx_D[1] = vct->getCoordinates(1)*getFluidNyc()/(double)vct->getYLEN();
         StartIdx_D[2] = vct->getCoordinates(2)*getFluidNzc()/(double)vct->getZLEN(); 

         cout<<" setGlobalStartIndex :: "<<StartIdx_D[0]<<", "<<StartIdx_D[1]<<", "<<StartIdx_D[2]<<endl;
         cout<<" setGlobalStartIndex :: "<<StartIdx_D[0]<<", "<<StartIdx_D[1]<<", "<<StartIdx_D[2]<<endl;
  };	

  // destructor
  ~InterfaceFluid(){
    delete INVarName_V;
    delArr4(State_GV,nIJK_D[0],nIJK_D[1],nIJK_D[2]);
  };

  // nIJK_D includes 1 guard/ghost cell layer...
  inline double getFluidNxc()
  {
    return(nIJK_D[0]-2*NG);
  };
  inline double getFluidNyc()
  {
    if(nIJK_D[1] > 2*NG)
      return(nIJK_D[1]-2*NG);
    else
      return(1);
  };
  inline double getFluidNzc()
  {
    if(nIJK_D[2] > 2*NG)
      return(nIJK_D[2]-2*NG);
    else
      return(1);
  };

  inline double getFluidLx(){ return(XYZ_D[0]); };
  inline double getFluidLy(){ return(XYZ_D[1]); };
  inline double getFluidLz(){ return(XYZ_D[2]); };

  // Get the bulk fluid velocities 
  inline double getFluidUx(const int i, const int j,const int k, const int is)
  {
    int ig,jg,kg;
    getGlobalIndex(i,j,k,&ig,&jg,&kg);
    return(State_GV[ig][jg][kg][iUx]);
  };
  inline double getFluidUy(const int i, const int j,const int k, const int is)
  {
    int ig,jg,kg;
    getGlobalIndex(i,j,k,&ig,&jg,&kg);
    return(State_GV[ig][jg][kg][iUy]);
  };
  inline double getFluidUz(const int i, const int j,const int k, const int is)
  {
    int ig,jg,kg;
    getGlobalIndex(i,j,k,&ig,&jg,&kg);
    return(State_GV[ig][jg][kg][iUz]);
  };

  // Get the fluid thermal velocities 
  inline double getFluidUthx(const int i, const int j,const int k, const int is)
  {
    int ig,jg,kg;
    getGlobalIndex(i,j,k,&ig,&jg,&kg);
    return(sqrt(State_GV[ig][jg][kg][iP]/State_GV[ig][jg][kg][iRho]));
    //return(1.0);
  };
  inline double getFluidUthy(const int i, const int j,const int k, const int is)
  {
    int ig,jg,kg;
    getGlobalIndex(i,j,k,&ig,&jg,&kg);
    return(sqrt(State_GV[ig][jg][kg][iP]/State_GV[ig][jg][kg][iRho]));
    //return(1.0);
  };
  inline double getFluidUthz(const int i, const int j,const int k, const int is)
  {
    int ig,jg,kg;
    getGlobalIndex(i,j,k,&ig,&jg,&kg);
    //return(sqrt(State_GV[ig][jg][kg][iP]/State_GV[ig][jg][kg][iRho]));
    return(1.0);
  };

  // this should return NUMBER DENSITY
  inline double getFluidRhoCenter(const int ii,const int jj, const int kk, const int is)
  {
    // is : index specis
    int i,j,k;
    getGlobalIndex(ii,jj,kk,&i,&j,&k);

    //cout<<ii<<", "<<jj<<", "<<kk<<", "<<i<<", "<<j<<", "<<k<<endl;

    if(FluidModel == "mhd")
        {
          // Electrons and iones have same density, ignoring is
          //cout<<"mhd: ne=ni"<<endl;
          return(State_GV[i][j][k][iRho]);
	  //return(1.0);
        } 
     else if (FluidModel == "hallmhd")
        { 
          cout<<"hallmhd not implementing, aborting! \n"<<flush;
          abort();
        }
     else
        {
          cout<<" Do not know FluidModel : "<<FluidModel<<", will abort!\n"<<flush;
          abort();
        }
  };

  // Get the Electic field as from the fluid description
  inline void setFluidFieldsCenter(
      double *Ex, double *Ey, double *Ez, 
      double *Bx, double *By, double *Bz, 
      const int ii,const int jj, const int kk)
  {
    int i,j,k;
    // solution is constant across ignored dimentions
    // 2 cell domain is handeld as a ignord dimention
    getGlobalIndex(ii,jj,kk,&i,&j,&k);
    if(FluidModel == "mhd")
        {
          //cout<<"mhd: E = - UxB"<<endl;
          (*Ex) = State_GV[i][j][k][iUz]*State_GV[i][j][k][iBy] - State_GV[i][j][k][iUy]*State_GV[i][j][k][iBz];
          (*Ey) = State_GV[i][j][k][iUx]*State_GV[i][j][k][iBz] - State_GV[i][j][k][iUz]*State_GV[i][j][k][iBx];
          (*Ez) = State_GV[i][j][k][iUy]*State_GV[i][j][k][iBx] - State_GV[i][j][k][iUx]*State_GV[i][j][k][iBy];

          (*Bx) = State_GV[i][j][k][iBx];
          (*By) = State_GV[i][j][k][iBy];
          (*Bz) = State_GV[i][j][k][iBz];
        } 
     else if (FluidModel == "hallmhd")
        { 
          cout<<"hallmhd not implementing, aborting! \n"<<flush;
          abort();
        }
     else
        {
          cout<<" Do not know FluidModel : "<<FluidModel<<", will abort!\n"<<flush;
          abort();
        }
  };


  // read the header file for the idl files
  void ReadIdlAscii()
  {
    
    string Filename = FilenameRoot+".h";
    ifstream myReadFile;
    
    // Helping variable to set a max waiting time when looking for the data header file.
    const int nMaxWait = 1000;
    int iWait = 0;
    
    // open with rw so it hopefully fails if the file write is not finished
    // writing to.
    myReadFile.open(Filename.c_str(),ifstream::in|ifstream::out);
    char output[100];
    while (!myReadFile.good()) 
      {
      // if file not ready we will wait for it
      sleep(0.1); 
      iWait++;
      cout<<iWait<<" Waiting......"<<endl;
      cout.flush();
      // file never came, we abort
      if(iWait > nMaxWait)
	{
	cout<<"we have waited too long for the file, aborting"<<endl;
	cout.flush();
	abort();
	}
      myReadFile.open(Filename.c_str(),ifstream::in|ifstream::out);	
      } 
    if(myReadFile.is_open()) 
      {	
	myReadFile >> TempString;
	myReadFile >> INnProc >> TempString;
	myReadFile >> INnStep >> TempString;
	myReadFile >> INtime  >> TempString;
	for(int i =0;i<6;i++){myReadFile >> INxRange_I[i];}
	myReadFile >> TempString;
	for(int i =0;i<3;i++){myReadFile >> INplotDx_D[i];}
	for(int i =0;i<3;i++){myReadFile >> INgridDx_D[i];}
        for(int i =0;i<4;i++){myReadFile >> TempString;}
	myReadFile >> INnPlotvar >> TempString;
	myReadFile >> INnEqpar   >> TempString;
        for(int i =0;i<2;i++){myReadFile >> TempString;}
        INVarName_V = new string[INnPlotvar+INnEqpar];
        for(int i =0;i<(INnPlotvar+INnEqpar);i++){myReadFile >> INVarName_V[i];}
        // text without interest
        for(int i =0;i<4;i++){myReadFile >> TempString;}
	for(int i =0;i<3;i++){myReadFile >> INnRoot_D[i];}
	myReadFile >> TempString;
	for(int i =0;i<3;i++){myReadFile >> INnIJK_D[i];}
        /*
	  myReadFile >> TempString;
	  while (!myReadFile.eof()) {
	  myReadFile >> output;
	  cout<<output<<endl;
	  cout.flush();
          }
        */
      }
    myReadFile.close();
    
    InitDataStorage();
    ReadVariables();
  };
  
};

#endif