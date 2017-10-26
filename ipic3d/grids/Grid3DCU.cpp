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
#include "Grid3DCU.h"
#include "Collective.h"
#include "MPIdata.h"
#include "Alloc.h"
#include "Com3DNonblk.h"
#include "debug.h"
#include "ipicmath.h"

/*! constructor */
Grid3DCU::Grid3DCU(CollectiveIO * col, VirtualTopology3D * vct):
  _vct(*vct)
{
  assert_le(col->getXLEN(),col->getNxc());
  assert_le(col->getYLEN(),col->getNyc());
  assert_le(col->getZLEN(),col->getNzc());

  // get number of cells restricted to regular (untruncated) subdomain
  const int nxc_rr = ceiling_of_ratio(col->getNxc(), col->getXLEN());
  const int nyc_rr = ceiling_of_ratio(col->getNyc(), col->getYLEN());
  const int nzc_rr = ceiling_of_ratio(col->getNzc(), col->getZLEN());

  nxc_r = nxc_rr;
  nyc_r = nyc_rr;
  nzc_r = nzc_rr;

  num_cells_rr = nxc_rr*nyc_rr*nzc_rr;

  // Truncate the number of cells appropriately in the upper process.
  // (We truncate rather than extend to avoid causing any load imbalance.)
  //
  if(vct->isXupper()) nxc_r = col->getNxc()-nxc_r*(col->getXLEN()-1);
  if(vct->isYupper()) nyc_r = col->getNyc()-nyc_r*(col->getYLEN()-1);
  if(vct->isZupper()) nzc_r = col->getNzc()-nzc_r*(col->getZLEN()-1);
  //
  assert_lt(0,nxc_r); assert_le(nxc_r,nxc_rr);
  assert_lt(0,nyc_r); assert_le(nyc_r,nyc_rr);
  assert_lt(0,nzc_r); assert_le(nzc_r,nzc_rr);

  // These restrictions should be removed.
  //
  // An objection to removing them is that the user can always
  // increase the number of mesh cells appropriately to be a
  // multiple of XLEN, but as a rejoinder, it is desirable to be
  // able to change the number of processors without changing
  // the discretized problem in any way and without divisibility
  // restrictions.
  //
  //assert_divides(vct->getXLEN(), col->getNxc());
  //assert_divides(vct->getYLEN(), col->getNyc());
  //assert_divides(vct->getZLEN(), col->getNzc());
  //
  //assert_eq(nxc_r,nxc_rr);
  //assert_eq(nyc_r,nyc_rr);
  //assert_eq(nzc_r,nzc_rr);

  // add two for ghost cells
  nxc = nxc_r + 2;
  nyc = nyc_r + 2;
  nzc = nzc_r + 2;

  dx = col->getDx();
  dy = col->getDy();
  dz = col->getDz();
  assert(dx == col->getLx() / col->getNxc());
  assert(dy == col->getLy() / col->getNyc());
  assert(dz == col->getLz() / col->getNzc());

  // local grid dimensions and boundaries of active nodes
  //
  // width of an ordinary subdomain
  //
  const double xWidth = dx*nxc_rr;
  const double yWidth = dy*nyc_rr;
  const double zWidth = dz*nzc_rr;
  //const double xWidth = (col->getLx() / (double) vct->getXLEN());
  //const double yWidth = (col->getLy() / (double) vct->getYLEN());
  //const double zWidth = (col->getLz() / (double) vct->getZLEN());
  //assert_almost_eq(dx*(nxc_r),xWidth,dx*1e-8);
  //assert_almost_eq(dy*(nyc_r),yWidth,dy*1e-8);
  //assert_almost_eq(dz*(nzc_r),zWidth,dz*1e-8);
  //
  xStart = vct->getCoordinates(0) * xWidth;
  yStart = vct->getCoordinates(1) * yWidth;
  zStart = vct->getCoordinates(2) * zWidth;
  //
  xEnd = xStart + xWidth - (nxc_rr-nxc_r)*dx;
  yEnd = yStart + yWidth - (nyc_rr-nyc_r)*dy;
  zEnd = zStart + zWidth - (nzc_rr-nzc_r)*dz;

  init_derived_parameters();
}

//Grid3DCU::Grid3DCU(
//  int nxc_, int nyc_, int nzc_,
//  double dx_, double dy_, double dz_,
//  double xStart_, double yStart_, double zStart_)
//: nxc(nxc_), nyc(nyc_), nzc(nzc_),
//  dx(dx_), dy(dy_), dz(dz_),
//  xStart(xStart_), yStart(yStart_), zStart(zStart_)
//{
//  const double xWidth = dx*(nxc-2);
//  const double yWidth = dy*(nyc-2);
//  const double zWidth = dz*(nzc-2);
//
//  xEnd = xStart + xWidth;
//  yEnd = yStart + yWidth;
//  zEnd = zStart + zWidth;
//
//  init_derived_parameters();
//}

// set derived convenience parameters
void Grid3DCU::init_derived_parameters()
{
  epsilon = (nxc+nyc+nzc)*1e-15;
  nxc_minus_epsilon = nxc-epsilon;
  nyc_minus_epsilon = nyc-epsilon;
  nzc_minus_epsilon = nzc-epsilon;
  assert_lt(int(floor(nxc_minus_epsilon)),nxc);
  assert_lt(int(floor(nyc_minus_epsilon)),nyc);
  assert_lt(int(floor(nzc_minus_epsilon)),nzc);
  xStart_g = xStart - dx;
  yStart_g = yStart - dy;
  zStart_g = zStart - dz;
  // calculation conveniences
  //
  VOL = dx * dy * dz;
  invVOL = 1.0 / VOL;
  invdx = 1.0 / dx;
  invdy = 1.0 / dy;
  invdz = 1.0 / dz;
  //
  nxn = nxc + 1;
  nyn = nyc + 1;
  nzn = nzc + 1;
  cxlast = nxc-1;
  cylast = nyc-1;
  czlast = nzc-1;

  // arrays allocation: nodes ---> the first node has index 1, the last has index nxn-2!
  pfloat_node_xcoord = new pfloat[nxn];
  pfloat_node_ycoord = new pfloat[nyn];
  pfloat_node_zcoord = new pfloat[nzn];
  node_xcoord = new double[nxn];
  node_ycoord = new double[nyn];
  node_zcoord = new double[nzn];
  for (int i=0; i<nxn; i++) node_xcoord[i] = xStart + (i - 1) * dx;
  for (int j=0; j<nyn; j++) node_ycoord[j] = yStart + (j - 1) * dy;
  for (int k=0; k<nzn; k++) node_zcoord[k] = zStart + (k - 1) * dz;
  for (int i=0; i<nxn; i++) pfloat_node_xcoord[i] = node_xcoord[i];
  for (int j=0; j<nyn; j++) pfloat_node_ycoord[j] = node_ycoord[j];
  for (int k=0; k<nzn; k++) pfloat_node_zcoord[k] = node_zcoord[k];
  // arrays allocation: cells ---> the first cell has index 1, the last has index ncn-2!
  center_xcoord = new double[nxc];
  center_ycoord = new double[nyc];
  center_zcoord = new double[nzc];
  for(int i=0; i<nxc; i++) center_xcoord[i] = .5*(node_xcoord[i]+node_xcoord[i+1]);
  for(int j=0; j<nyc; j++) center_ycoord[j] = .5*(node_ycoord[j]+node_ycoord[j+1]);
  for(int k=0; k<nzc; k++) center_zcoord[k] = .5*(node_zcoord[k]+node_zcoord[k+1]);
  //num_cells_r = nxc_r*nyc_r*nzc_r;
  //num_cells = nxc*nyc*nzc;
}

/** deallocate the local grid */
Grid3DCU::~Grid3DCU() {
  delete [] node_xcoord;
  delete [] node_ycoord;
  delete [] node_zcoord;
  delete [] center_xcoord;
  delete [] center_ycoord;
  delete [] center_zcoord;
}

/** print the local grid info */
void Grid3DCU::print()const
{
  const VirtualTopology3D *ptVCT = &get_vct();
  printf("\nSubgrid (%d,%d,%d)\n",
    ptVCT->getCoordinates(0),
    ptVCT->getCoordinates(1),
    ptVCT->getCoordinates(2));
  printf("Number of cells: X:%d, Y:%d, Z:%d\n",
    nxc - 2,
    nyc - 2,
    nzc - 2);
  printf(
    "Xin = %g; Xfin = %g\n"
    "Yin = %g; Yfin = %g\n"
    "Zin = %g; Zfin = %g\n\n",
    node_xcoord[1], node_xcoord[nxn - 2],
    node_ycoord[1], node_ycoord[nyn - 2],
    node_zcoord[1], node_zcoord[nzn - 2]);
}

/** calculate gradient on nodes, given a scalar field defined on central points  */
void Grid3DCU::gradC2N(arr3_double gradXN, arr3_double gradYN, arr3_double gradZN, const_arr3_double scFieldC)const
{
  for (register int i = 1; i < nxn - 1; i++)
    for (register int j = 1; j < nyn - 1; j++)
      for (register int k = 1; k < nzn - 1; k++) {
        gradXN[i][j][k] = .25 * (scFieldC[i][j][k] - scFieldC[i - 1][j][k]) * invdx + .25 * (scFieldC[i][j][k - 1] - scFieldC[i - 1][j][k - 1]) * invdx + .25 * (scFieldC[i][j - 1][k] - scFieldC[i - 1][j - 1][k]) * invdx + .25 * (scFieldC[i][j - 1][k - 1] - scFieldC[i - 1][j - 1][k - 1]) * invdx;
        gradYN[i][j][k] = .25 * (scFieldC[i][j][k] - scFieldC[i][j - 1][k]) * invdy + .25 * (scFieldC[i][j][k - 1] - scFieldC[i][j - 1][k - 1]) * invdy + .25 * (scFieldC[i - 1][j][k] - scFieldC[i - 1][j - 1][k]) * invdy + .25 * (scFieldC[i - 1][j][k - 1] - scFieldC[i - 1][j - 1][k - 1]) * invdy;
        gradZN[i][j][k] = .25 * (scFieldC[i][j][k] - scFieldC[i][j][k - 1]) * invdz + .25 * (scFieldC[i - 1][j][k] - scFieldC[i - 1][j][k - 1]) * invdz + .25 * (scFieldC[i][j - 1][k] - scFieldC[i][j - 1][k - 1]) * invdz + .25 * (scFieldC[i - 1][j - 1][k] - scFieldC[i - 1][j - 1][k - 1]) * invdz;
      }
}

/** calculate gradient on nodes, given a scalar field defined on central points  */
void Grid3DCU::gradN2C(arr3_double gradXC, arr3_double gradYC, arr3_double gradZC, const_arr3_double scFieldN)const
{
  for (register int i = 1; i < nxc - 1; i++)
    for (register int j = 1; j < nyc - 1; j++)
      for (register int k = 1; k < nzc - 1; k++) {
        gradXC[i][j][k] = .25 * (scFieldN[i + 1][j][k] - scFieldN[i][j][k]) * invdx + .25 * (scFieldN[i + 1][j][k + 1] - scFieldN[i][j][k + 1]) * invdx + .25 * (scFieldN[i + 1][j + 1][k] - scFieldN[i][j + 1][k]) * invdx + .25 * (scFieldN[i + 1][j + 1][k + 1] - scFieldN[i][j + 1][k + 1]) * invdx;
        gradYC[i][j][k] = .25 * (scFieldN[i][j + 1][k] - scFieldN[i][j][k]) * invdy + .25 * (scFieldN[i][j + 1][k + 1] - scFieldN[i][j][k + 1]) * invdy + .25 * (scFieldN[i + 1][j + 1][k] - scFieldN[i + 1][j][k]) * invdy + .25 * (scFieldN[i + 1][j + 1][k + 1] - scFieldN[i + 1][j][k + 1]) * invdy;
        gradZC[i][j][k] = .25 * (scFieldN[i][j][k + 1] - scFieldN[i][j][k]) * invdz + .25 * (scFieldN[i + 1][j][k + 1] - scFieldN[i + 1][j][k]) * invdz + .25 * (scFieldN[i][j + 1][k + 1] - scFieldN[i][j + 1][k]) * invdz + .25 * (scFieldN[i + 1][j + 1][k + 1] - scFieldN[i + 1][j + 1][k]) * invdz;
      }
}

/** calculate divergence on central points, given a vector field defined on nodes  */
void Grid3DCU::divN2C(arr3_double divC, const_arr3_double vecFieldXN, const_arr3_double vecFieldYN, const_arr3_double vecFieldZN)const
{
  double compX;
  double compY;
  double compZ;
  for (register int i = 1; i < nxc - 1; i++)
    for (register int j = 1; j < nyc - 1; j++)
      for (register int k = 1; k < nzc - 1; k++) {
        compX = .25 * (vecFieldXN[i + 1][j][k] - vecFieldXN[i][j][k]) * invdx + .25 * (vecFieldXN[i + 1][j][k + 1] - vecFieldXN[i][j][k + 1]) * invdx + .25 * (vecFieldXN[i + 1][j + 1][k] - vecFieldXN[i][j + 1][k]) * invdx + .25 * (vecFieldXN[i + 1][j + 1][k + 1] - vecFieldXN[i][j + 1][k + 1]) * invdx;
        compY = .25 * (vecFieldYN[i][j + 1][k] - vecFieldYN[i][j][k]) * invdy + .25 * (vecFieldYN[i][j + 1][k + 1] - vecFieldYN[i][j][k + 1]) * invdy + .25 * (vecFieldYN[i + 1][j + 1][k] - vecFieldYN[i + 1][j][k]) * invdy + .25 * (vecFieldYN[i + 1][j + 1][k + 1] - vecFieldYN[i + 1][j][k + 1]) * invdy;
        compZ = .25 * (vecFieldZN[i][j][k + 1] - vecFieldZN[i][j][k]) * invdz + .25 * (vecFieldZN[i + 1][j][k + 1] - vecFieldZN[i + 1][j][k]) * invdz + .25 * (vecFieldZN[i][j + 1][k + 1] - vecFieldZN[i][j + 1][k]) * invdz + .25 * (vecFieldZN[i + 1][j + 1][k + 1] - vecFieldZN[i + 1][j + 1][k]) * invdz;
        divC[i][j][k] = compX + compY + compZ;
      }
}

/** calculate divergence on central points, given a Tensor field defined on nodes  */
void Grid3DCU::divSymmTensorN2C(arr3_double divCX, arr3_double divCY, arr3_double divCZ, const_arr4_double pXX, const_arr4_double pXY, const_arr4_double pXZ, const_arr4_double pYY, const_arr4_double pYZ, const_arr4_double pZZ, int ns)const
{
  double comp1X, comp2X, comp3X;
  double comp1Y, comp2Y, comp3Y;
  double comp1Z, comp2Z, comp3Z;
  for (register int i = 1; i < nxc - 1; i++)
    for (register int j = 1; j < nyc - 1; j++)
      for (register int k = 1; k < nzc - 1; k++) {
        comp1X = .25 * (pXX[ns][i + 1][j][k] - pXX[ns][i][j][k]) * invdx + .25 * (pXX[ns][i + 1][j][k + 1] - pXX[ns][i][j][k + 1]) * invdx + .25 * (pXX[ns][i + 1][j + 1][k] - pXX[ns][i][j + 1][k]) * invdx + .25 * (pXX[ns][i + 1][j + 1][k + 1] - pXX[ns][i][j + 1][k + 1]) * invdx;
        comp2X = .25 * (pXY[ns][i + 1][j][k] - pXY[ns][i][j][k]) * invdx + .25 * (pXY[ns][i + 1][j][k + 1] - pXY[ns][i][j][k + 1]) * invdx + .25 * (pXY[ns][i + 1][j + 1][k] - pXY[ns][i][j + 1][k]) * invdx + .25 * (pXY[ns][i + 1][j + 1][k + 1] - pXY[ns][i][j + 1][k + 1]) * invdx;
        comp3X = .25 * (pXZ[ns][i + 1][j][k] - pXZ[ns][i][j][k]) * invdx + .25 * (pXZ[ns][i + 1][j][k + 1] - pXZ[ns][i][j][k + 1]) * invdx + .25 * (pXZ[ns][i + 1][j + 1][k] - pXZ[ns][i][j + 1][k]) * invdx + .25 * (pXZ[ns][i + 1][j + 1][k + 1] - pXZ[ns][i][j + 1][k + 1]) * invdx;
        comp1Y = .25 * (pXY[ns][i][j + 1][k] - pXY[ns][i][j][k]) * invdy + .25 * (pXY[ns][i][j + 1][k + 1] - pXY[ns][i][j][k + 1]) * invdy + .25 * (pXY[ns][i + 1][j + 1][k] - pXY[ns][i + 1][j][k]) * invdy + .25 * (pXY[ns][i + 1][j + 1][k + 1] - pXY[ns][i + 1][j][k + 1]) * invdy;
        comp2Y = .25 * (pYY[ns][i][j + 1][k] - pYY[ns][i][j][k]) * invdy + .25 * (pYY[ns][i][j + 1][k + 1] - pYY[ns][i][j][k + 1]) * invdy + .25 * (pYY[ns][i + 1][j + 1][k] - pYY[ns][i + 1][j][k]) * invdy + .25 * (pYY[ns][i + 1][j + 1][k + 1] - pYY[ns][i + 1][j][k + 1]) * invdy;
        comp3Y = .25 * (pYZ[ns][i][j + 1][k] - pYZ[ns][i][j][k]) * invdy + .25 * (pYZ[ns][i][j + 1][k + 1] - pYZ[ns][i][j][k + 1]) * invdy + .25 * (pYZ[ns][i + 1][j + 1][k] - pYZ[ns][i + 1][j][k]) * invdy + .25 * (pYZ[ns][i + 1][j + 1][k + 1] - pYZ[ns][i + 1][j][k + 1]) * invdy;
        comp1Z = .25 * (pXZ[ns][i][j][k + 1] - pXZ[ns][i][j][k]) * invdz + .25 * (pXZ[ns][i + 1][j][k + 1] - pXZ[ns][i + 1][j][k]) * invdz + .25 * (pXZ[ns][i][j + 1][k + 1] - pXZ[ns][i][j + 1][k]) * invdz + .25 * (pXZ[ns][i + 1][j + 1][k + 1] - pXZ[ns][i + 1][j + 1][k]) * invdz;
        comp2Z = .25 * (pYZ[ns][i][j][k + 1] - pYZ[ns][i][j][k]) * invdz + .25 * (pYZ[ns][i + 1][j][k + 1] - pYZ[ns][i + 1][j][k]) * invdz + .25 * (pYZ[ns][i][j + 1][k + 1] - pYZ[ns][i][j + 1][k]) * invdz + .25 * (pYZ[ns][i + 1][j + 1][k + 1] - pYZ[ns][i + 1][j + 1][k]) * invdz;
        comp3Z = .25 * (pZZ[ns][i][j][k + 1] - pZZ[ns][i][j][k]) * invdz + .25 * (pZZ[ns][i + 1][j][k + 1] - pZZ[ns][i + 1][j][k]) * invdz + .25 * (pZZ[ns][i][j + 1][k + 1] - pZZ[ns][i][j + 1][k]) * invdz + .25 * (pZZ[ns][i + 1][j + 1][k + 1] - pZZ[ns][i + 1][j + 1][k]) * invdz;
        divCX[i][j][k] = comp1X + comp2X + comp3X;
        divCY[i][j][k] = comp1Y + comp2Y + comp3Y;
        divCZ[i][j][k] = comp1Z + comp2Z + comp3Z;
      }
}

/** calculate divergence on nodes, given a vector field defined on central points  */
void Grid3DCU::divC2N(arr3_double divN, const_arr3_double vecFieldXC, const_arr3_double vecFieldYC, const_arr3_double vecFieldZC)const
{
  double compX;
  double compY;
  double compZ;
  for (register int i = 1; i < nxn - 1; i++)
    for (register int j = 1; j < nyn - 1; j++)
      for (register int k = 1; k < nzn - 1; k++) {
        compX = .25 * (vecFieldXC[i][j][k] - vecFieldXC[i - 1][j][k]) * invdx + .25 * (vecFieldXC[i][j][k - 1] - vecFieldXC[i - 1][j][k - 1]) * invdx + .25 * (vecFieldXC[i][j - 1][k] - vecFieldXC[i - 1][j - 1][k]) * invdx + .25 * (vecFieldXC[i][j - 1][k - 1] - vecFieldXC[i - 1][j - 1][k - 1]) * invdx;
        compY = .25 * (vecFieldYC[i][j][k] - vecFieldYC[i][j - 1][k]) * invdy + .25 * (vecFieldYC[i][j][k - 1] - vecFieldYC[i][j - 1][k - 1]) * invdy + .25 * (vecFieldYC[i - 1][j][k] - vecFieldYC[i - 1][j - 1][k]) * invdy + .25 * (vecFieldYC[i - 1][j][k - 1] - vecFieldYC[i - 1][j - 1][k - 1]) * invdy;
        compZ = .25 * (vecFieldZC[i][j][k] - vecFieldZC[i][j][k - 1]) * invdz + .25 * (vecFieldZC[i - 1][j][k] - vecFieldZC[i - 1][j][k - 1]) * invdz + .25 * (vecFieldZC[i][j - 1][k] - vecFieldZC[i][j - 1][k - 1]) * invdz + .25 * (vecFieldZC[i - 1][j - 1][k] - vecFieldZC[i - 1][j - 1][k - 1]) * invdz;
        divN[i][j][k] = compX + compY + compZ;
      }
}

/** calculate curl on nodes, given a vector field defined on central points  */
void Grid3DCU::curlC2N(arr3_double curlXN, arr3_double curlYN, arr3_double curlZN, const_arr3_double vecFieldXC, const_arr3_double vecFieldYC, const_arr3_double vecFieldZC)const
{
  double compZDY, compYDZ;
  double compXDZ, compZDX;
  double compYDX, compXDY;
  for (register int i = 1; i < nxn - 1; i++)
    for (register int j = 1; j < nyn - 1; j++)
      for (register int k = 1; k < nzn - 1; k++) {
        // curl - X
        compZDY = .25 * (vecFieldZC[i][j][k] - vecFieldZC[i][j - 1][k]) * invdy + .25 * (vecFieldZC[i][j][k - 1] - vecFieldZC[i][j - 1][k - 1]) * invdy + .25 * (vecFieldZC[i - 1][j][k] - vecFieldZC[i - 1][j - 1][k]) * invdy + .25 * (vecFieldZC[i - 1][j][k - 1] - vecFieldZC[i - 1][j - 1][k - 1]) * invdy;
        compYDZ = .25 * (vecFieldYC[i][j][k] - vecFieldYC[i][j][k - 1]) * invdz + .25 * (vecFieldYC[i - 1][j][k] - vecFieldYC[i - 1][j][k - 1]) * invdz + .25 * (vecFieldYC[i][j - 1][k] - vecFieldYC[i][j - 1][k - 1]) * invdz + .25 * (vecFieldYC[i - 1][j - 1][k] - vecFieldYC[i - 1][j - 1][k - 1]) * invdz;
        // curl - Y
        compXDZ = .25 * (vecFieldXC[i][j][k] - vecFieldXC[i][j][k - 1]) * invdz + .25 * (vecFieldXC[i - 1][j][k] - vecFieldXC[i - 1][j][k - 1]) * invdz + .25 * (vecFieldXC[i][j - 1][k] - vecFieldXC[i][j - 1][k - 1]) * invdz + .25 * (vecFieldXC[i - 1][j - 1][k] - vecFieldXC[i - 1][j - 1][k - 1]) * invdz;
        compZDX = .25 * (vecFieldZC[i][j][k] - vecFieldZC[i - 1][j][k]) * invdx + .25 * (vecFieldZC[i][j][k - 1] - vecFieldZC[i - 1][j][k - 1]) * invdx + .25 * (vecFieldZC[i][j - 1][k] - vecFieldZC[i - 1][j - 1][k]) * invdx + .25 * (vecFieldZC[i][j - 1][k - 1] - vecFieldZC[i - 1][j - 1][k - 1]) * invdx;
        // curl - Z
        compYDX = .25 * (vecFieldYC[i][j][k] - vecFieldYC[i - 1][j][k]) * invdx + .25 * (vecFieldYC[i][j][k - 1] - vecFieldYC[i - 1][j][k - 1]) * invdx + .25 * (vecFieldYC[i][j - 1][k] - vecFieldYC[i - 1][j - 1][k]) * invdx + .25 * (vecFieldYC[i][j - 1][k - 1] - vecFieldYC[i - 1][j - 1][k - 1]) * invdx;
        compXDY = .25 * (vecFieldXC[i][j][k] - vecFieldXC[i][j - 1][k]) * invdy + .25 * (vecFieldXC[i][j][k - 1] - vecFieldXC[i][j - 1][k - 1]) * invdy + .25 * (vecFieldXC[i - 1][j][k] - vecFieldXC[i - 1][j - 1][k]) * invdy + .25 * (vecFieldXC[i - 1][j][k - 1] - vecFieldXC[i - 1][j - 1][k - 1]) * invdy;

        curlXN[i][j][k] = compZDY - compYDZ;
        curlYN[i][j][k] = compXDZ - compZDX;
        curlZN[i][j][k] = compYDX - compXDY;
      }
}

/** calculate curl on central points, given a vector field defined on nodes  */
void Grid3DCU::curlN2C(arr3_double curlXC, arr3_double curlYC, arr3_double curlZC,
  const_arr3_double vecFieldXN, const_arr3_double vecFieldYN, const_arr3_double vecFieldZN)const
{
  double compZDY, compYDZ;
  double compXDZ, compZDX;
  double compYDX, compXDY;
  for (register int i = 1; i < nxc - 1; i++)
    for (register int j = 1; j < nyc - 1; j++)
      for (register int k = 1; k < nzc - 1; k++) {
        // curl - X
        compZDY = .25 * (vecFieldZN[i][j + 1][k] - vecFieldZN[i][j][k]) * invdy + .25 * (vecFieldZN[i][j + 1][k + 1] - vecFieldZN[i][j][k + 1]) * invdy + .25 * (vecFieldZN[i + 1][j + 1][k] - vecFieldZN[i + 1][j][k]) * invdy + .25 * (vecFieldZN[i + 1][j + 1][k + 1] - vecFieldZN[i + 1][j][k + 1]) * invdy;
        compYDZ = .25 * (vecFieldYN[i][j][k + 1] - vecFieldYN[i][j][k]) * invdz + .25 * (vecFieldYN[i + 1][j][k + 1] - vecFieldYN[i + 1][j][k]) * invdz + .25 * (vecFieldYN[i][j + 1][k + 1] - vecFieldYN[i][j + 1][k]) * invdz + .25 * (vecFieldYN[i + 1][j + 1][k + 1] - vecFieldYN[i + 1][j + 1][k]) * invdz;
        // curl - Y
        compXDZ = .25 * (vecFieldXN[i][j][k + 1] - vecFieldXN[i][j][k]) * invdz + .25 * (vecFieldXN[i + 1][j][k + 1] - vecFieldXN[i + 1][j][k]) * invdz + .25 * (vecFieldXN[i][j + 1][k + 1] - vecFieldXN[i][j + 1][k]) * invdz + .25 * (vecFieldXN[i + 1][j + 1][k + 1] - vecFieldXN[i + 1][j + 1][k]) * invdz;
        compZDX = .25 * (vecFieldZN[i + 1][j][k] - vecFieldZN[i][j][k]) * invdx + .25 * (vecFieldZN[i + 1][j][k + 1] - vecFieldZN[i][j][k + 1]) * invdx + .25 * (vecFieldZN[i + 1][j + 1][k] - vecFieldZN[i][j + 1][k]) * invdx + .25 * (vecFieldZN[i + 1][j + 1][k + 1] - vecFieldZN[i][j + 1][k + 1]) * invdx;
        // curl - Z
        compYDX = .25 * (vecFieldYN[i + 1][j][k] - vecFieldYN[i][j][k]) * invdx + .25 * (vecFieldYN[i + 1][j][k + 1] - vecFieldYN[i][j][k + 1]) * invdx + .25 * (vecFieldYN[i + 1][j + 1][k] - vecFieldYN[i][j + 1][k]) * invdx + .25 * (vecFieldYN[i + 1][j + 1][k + 1] - vecFieldYN[i][j + 1][k + 1]) * invdx;
        compXDY = .25 * (vecFieldXN[i][j + 1][k] - vecFieldXN[i][j][k]) * invdy + .25 * (vecFieldXN[i][j + 1][k + 1] - vecFieldXN[i][j][k + 1]) * invdy + .25 * (vecFieldXN[i + 1][j + 1][k] - vecFieldXN[i + 1][j][k]) * invdy + .25 * (vecFieldXN[i + 1][j + 1][k + 1] - vecFieldXN[i + 1][j][k + 1]) * invdy;


        curlXC[i][j][k] = compZDY - compYDZ;
        curlYC[i][j][k] = compXDZ - compZDX;
        curlZC[i][j][k] = compYDX - compXDY;
      }



}

/** calculate laplacian on nodes, given a scalar field defined on nodes */
void Grid3DCU::lapN2N(arr3_double lapN, const_arr3_double scFieldN,EMfields3D *EMf)const
{
  const VirtualTopology3D *vct = &get_vct();
  // calculate laplacian as divercence of gradient
  // allocate 3 gradients: defined on central points
  array3_double gradXC(nxc, nyc, nzc);
  array3_double gradYC(nxc, nyc, nzc);
  array3_double gradZC(nxc, nyc, nzc);

  gradN2C(gradXC, gradYC, gradZC, scFieldN);
  // communicate with BC
  communicateCenterBC(nxc, nyc, nzc, gradXC, 1, 1, 1, 1, 1, 1, vct, EMf);
  communicateCenterBC(nxc, nyc, nzc, gradYC, 1, 1, 1, 1, 1, 1, vct, EMf);
  communicateCenterBC(nxc, nyc, nzc, gradZC, 1, 1, 1, 1, 1, 1, vct, EMf);
  divC2N(lapN, gradXC, gradYC, gradZC);
}

/** calculate laplacian on central points, given a scalar field defined on central points */
void Grid3DCU::lapC2C(arr3_double lapC, const_arr3_double scFieldC)const
{
  const VirtualTopology3D *vct = &get_vct();
  // calculate laplacian as divercence of gradient
  // allocate 3 gradients: defined on nodes
  array3_double gradXN(nxn, nyn, nzn);
  array3_double gradYN(nxn, nyn, nzn);
  array3_double gradZN(nxn, nyn, nzn);

  gradC2N(gradXN, gradYN, gradZN, scFieldC);
  if (vct->getYleft_neighbor() == MPI_PROC_NULL) {
    for (int ii = 0; ii < nxn; ii++)
      for (int kk = 0; kk < nzn; kk++) {
        gradXN[ii][0][kk] = 0.0;
        gradXN[ii][1][kk] = 0.0;
        gradXN[ii][2][kk] = 0.0;
        gradZN[ii][0][kk] = 0.0;
        gradZN[ii][1][kk] = 0.0;
        gradZN[ii][2][kk] = 0.0;
        // gradYN[ii][1][kk] = gradYN[ii][2][kk];
      }
  }
  if (vct->getYright_neighbor() == MPI_PROC_NULL) {
    for (int ii = 0; ii < nxn; ii++)
      for (int kk = 0; kk < nzn; kk++) {
        gradXN[ii][nyn - 1][kk] = 0.0;
        gradXN[ii][nyn - 2][kk] = 0.0;
        gradXN[ii][nyn - 3][kk] = 0.0;
        gradZN[ii][nyn - 1][kk] = 0.0;
        gradZN[ii][nyn - 2][kk] = 0.0;
        gradZN[ii][nyn - 3][kk] = 0.0;
        // gradYN[ii][nyn-2][kk] = gradYN[ii][nyc-3][kk];
      }
  }
  divN2C(lapC, gradXN, gradYN, gradZN);
}

/** calculate laplacian on central points, given a scalar field defined on central points for Poisson */
void Grid3DCU::lapC2Cpoisson(arr3_double lapC, arr3_double scFieldC,EMfields3D *EMf)const
{
  const VirtualTopology3D *vct = &get_vct();
  // communicate first the scFieldC
  communicateCenterBoxStencilBC(nxc, nyc, nzc, scFieldC, 1, 1, 1, 1, 1, 1, vct, EMf);
  for (register int i = 1; i < nxc - 1; i++)
    for (register int j = 1; j < nyc - 1; j++)
      for (register int k = 1; k < nzc - 1; k++)
        lapC[i][j][k] = (scFieldC[i - 1][j][k] - 2 * scFieldC[i][j][k] + scFieldC[i + 1][j][k]) * invdx * invdx + (scFieldC[i][j - 1][k] - 2 * scFieldC[i][j][k] + scFieldC[i][j + 1][k]) * invdy * invdy + (scFieldC[i][j][k - 1] - 2 * scFieldC[i][j][k] + scFieldC[i][j][k + 1]) * invdz * invdz;
}

/** calculate divergence on  boundaries */
void Grid3DCU::divBCleft(arr3_double divBC, const_arr3_double vectorX, const_arr3_double vectorY, const_arr3_double vectorZ, int leftActiveNode, int dirDER)const
{
  double compX, compY, compZ;
  switch (dirDER) {
    case 0:                    // DIVERGENCE DIRECTION X
      for (register int j = 1; j < nyn - 2; j++)
        for (register int k = 1; k < nzn - 2; k++) {
          compX = .25 * (vectorX[leftActiveNode + 1][j][k] - vectorX[leftActiveNode][j][k]) * invdx + .25 * (vectorX[leftActiveNode + 1][j][k + 1] - vectorX[leftActiveNode][j][k + 1]) * invdx + .25 * (vectorX[leftActiveNode + 1][j + 1][k] - vectorX[leftActiveNode][j + 1][k]) * invdx + .25 * (vectorX[leftActiveNode + 1][j + 1][k + 1] - vectorX[leftActiveNode][j + 1][k + 1]) * invdx;
          compY = .25 * (vectorY[leftActiveNode][j + 1][k] - vectorY[leftActiveNode][j][k]) * invdy + .25 * (vectorY[leftActiveNode][j + 1][k + 1] - vectorY[leftActiveNode][j][k + 1]) * invdy + .25 * (vectorY[leftActiveNode + 1][j + 1][k] - vectorY[leftActiveNode + 1][j][k]) * invdy + .25 * (vectorY[leftActiveNode + 1][j + 1][k + 1] - vectorY[leftActiveNode + 1][j][k + 1]) * invdy;
          compZ = .25 * (vectorZ[leftActiveNode][j][k + 1] - vectorZ[leftActiveNode][j][k]) * invdz + .25 * (vectorZ[leftActiveNode + 1][j][k + 1] - vectorZ[leftActiveNode + 1][j][k]) * invdz + .25 * (vectorZ[leftActiveNode][j + 1][k + 1] - vectorZ[leftActiveNode][j + 1][k]) * invdz + .25 * (vectorZ[leftActiveNode + 1][j + 1][k + 1] - vectorZ[leftActiveNode + 1][j + 1][k]) * invdz;
          divBC[leftActiveNode][j][k] = compX + compY + compZ;
        }
      break;
    case 1:                    // DIVERGENCE DIRECTION Y
      for (register int i = 1; i < nxn - 2; i++)
        for (register int k = 1; k < nzn - 2; k++) {
          compX = .25 * (vectorX[i + 1][leftActiveNode][k] - vectorX[i][leftActiveNode][k]) * invdx + .25 * (vectorX[i + 1][leftActiveNode][k + 1] - vectorX[i][leftActiveNode][k + 1]) * invdx + .25 * (vectorX[i + 1][leftActiveNode + 1][k] - vectorX[i][leftActiveNode + 1][k]) * invdx + .25 * (vectorX[i + 1][leftActiveNode + 1][k + 1] - vectorX[i][leftActiveNode + 1][k + 1]) * invdx;
          compY = .25 * (vectorY[i][leftActiveNode + 1][k] - vectorY[i][leftActiveNode][k]) * invdy + .25 * (vectorY[i][leftActiveNode + 1][k + 1] - vectorY[i][leftActiveNode][k + 1]) * invdy + .25 * (vectorY[i + 1][leftActiveNode + 1][k] - vectorY[i + 1][leftActiveNode][k]) * invdy + .25 * (vectorY[i + 1][leftActiveNode + 1][k + 1] - vectorY[i + 1][leftActiveNode][k + 1]) * invdy;
          compZ = .25 * (vectorZ[i][leftActiveNode][k + 1] - vectorZ[i][leftActiveNode][k]) * invdz + .25 * (vectorZ[i + 1][leftActiveNode][k + 1] - vectorZ[i + 1][leftActiveNode][k]) * invdz + .25 * (vectorZ[i][leftActiveNode + 1][k + 1] - vectorZ[i][leftActiveNode + 1][k]) * invdz + .25 * (vectorZ[i + 1][leftActiveNode + 1][k + 1] - vectorZ[i + 1][leftActiveNode + 1][k]) * invdz;
          divBC[i][leftActiveNode][k] = compX + compY + compZ;
        }
      break;
    case 2:                    // DIVERGENCE DIRECTION Z
      for (register int i = 1; i < nxn - 2; i++)
        for (register int j = 1; j < nyn - 2; j++) {
          compX = .25 * (vectorX[i + 1][j][leftActiveNode] - vectorX[i][j][leftActiveNode]) * invdx + .25 * (vectorX[i + 1][j][leftActiveNode + 1] - vectorX[i][j][leftActiveNode + 1]) * invdx + .25 * (vectorX[i + 1][j + 1][leftActiveNode] - vectorX[i][j + 1][leftActiveNode]) * invdx + .25 * (vectorX[i + 1][j + 1][leftActiveNode + 1] - vectorX[i][j + 1][leftActiveNode + 1]) * invdx;
          compY = .25 * (vectorY[i][j + 1][leftActiveNode] - vectorY[i][j][leftActiveNode]) * invdy + .25 * (vectorY[i][j + 1][leftActiveNode + 1] - vectorY[i][j][leftActiveNode + 1]) * invdy + .25 * (vectorY[i + 1][j + 1][leftActiveNode] - vectorY[i + 1][j][leftActiveNode]) * invdy + .25 * (vectorY[i + 1][j + 1][leftActiveNode + 1] - vectorY[i + 1][j][leftActiveNode + 1]) * invdy;
          compZ = .25 * (vectorZ[i][j][leftActiveNode + 1] - vectorZ[i][j][leftActiveNode]) * invdz + .25 * (vectorZ[i + 1][j][leftActiveNode + 1] - vectorZ[i + 1][j][leftActiveNode]) * invdz + .25 * (vectorZ[i][j + 1][leftActiveNode + 1] - vectorZ[i][j + 1][leftActiveNode]) * invdz + .25 * (vectorZ[i + 1][j + 1][leftActiveNode + 1] - vectorZ[i + 1][j + 1][leftActiveNode]) * invdz;
          divBC[i][j][leftActiveNode] = compX + compY + compZ;
        }
      break;

  }


}

/** calculate divergence on  boundaries */
void Grid3DCU::divBCright(arr3_double divBC, const_arr3_double vectorX, const_arr3_double vectorY, const_arr3_double vectorZ, int rightActiveNode, int dirDER)const
{
  double compX, compY, compZ;


  switch (dirDER) {
    case 0:                    // DIVERGENCE DIRECTION X
      for (register int j = 1; j < nxn - 2; j++)
        for (register int k = 1; k < nzn - 2; k++) {
          compX = .25 * (vectorX[rightActiveNode][j][k] - vectorX[rightActiveNode - 1][j][k]) * invdx + .25 * (vectorX[rightActiveNode][j][k + 1] - vectorX[rightActiveNode - 1][j][k + 1]) * invdx + .25 * (vectorX[rightActiveNode][j + 1][k] - vectorX[rightActiveNode - 1][j + 1][k]) * invdx + .25 * (vectorX[rightActiveNode][j + 1][k + 1] - vectorX[rightActiveNode - 1][j + 1][k + 1]) * invdx;
          compY = .25 * (vectorY[rightActiveNode][j + 1][k] - vectorY[rightActiveNode][j][k]) * invdy + .25 * (vectorY[rightActiveNode][j + 1][k + 1] - vectorY[rightActiveNode][j][k + 1]) * invdy + .25 * (vectorY[rightActiveNode - 1][j + 1][k] - vectorY[rightActiveNode - 1][j][k]) * invdy + .25 * (vectorY[rightActiveNode - 1][j + 1][k + 1] - vectorY[rightActiveNode - 1][j][k + 1]) * invdy;
          compZ = .25 * (vectorZ[rightActiveNode][j][k + 1] - vectorZ[rightActiveNode][j][k]) * invdz + .25 * (vectorZ[rightActiveNode - 1][j][k + 1] - vectorZ[rightActiveNode - 1][j][k]) * invdz + .25 * (vectorZ[rightActiveNode][j + 1][k + 1] - vectorZ[rightActiveNode][j + 1][k]) * invdz + .25 * (vectorZ[rightActiveNode - 1][j + 1][k + 1] - vectorZ[rightActiveNode - 1][j + 1][k]) * invdz;
          divBC[rightActiveNode][j][k] = compX + compY + compZ;
        }
      break;
    case 1:                    // DIVERGENCE DIRECTION Y
      for (register int i = 1; i < nxn - 2; i++)
        for (register int k = 1; k < nzn - 2; k++) {
          compX = .25 * (vectorX[i + 1][rightActiveNode][k] - vectorX[i][rightActiveNode][k]) * invdx + .25 * (vectorX[i + 1][rightActiveNode][k + 1] - vectorX[i][rightActiveNode][k + 1]) * invdx + .25 * (vectorX[i + 1][rightActiveNode - 1][k] - vectorX[i][rightActiveNode - 1][k]) * invdx + .25 * (vectorX[i + 1][rightActiveNode - 1][k + 1] - vectorX[i][rightActiveNode - 1][k + 1]) * invdx;
          compY = .25 * (vectorY[i][rightActiveNode][k] - vectorY[i][rightActiveNode - 1][k]) * invdy + .25 * (vectorY[i][rightActiveNode][k + 1] - vectorY[i][rightActiveNode - 1][k + 1]) * invdy + .25 * (vectorY[i + 1][rightActiveNode][k] - vectorY[i + 1][rightActiveNode - 1][k]) * invdy + .25 * (vectorY[i + 1][rightActiveNode + 1][k + 1] - vectorY[i + 1][rightActiveNode][k + 1]) * invdy;
          compZ = .25 * (vectorZ[i][rightActiveNode][k + 1] - vectorZ[i][rightActiveNode][k]) * invdz + .25 * (vectorZ[i + 1][rightActiveNode][k + 1] - vectorZ[i + 1][rightActiveNode][k]) * invdz + .25 * (vectorZ[i][rightActiveNode - 1][k + 1] - vectorZ[i][rightActiveNode - 1][k]) * invdz + .25 * (vectorZ[i + 1][rightActiveNode - 1][k + 1] - vectorZ[i + 1][rightActiveNode - 1][k]) * invdz;
          divBC[i][rightActiveNode][k] = compX + compY + compZ;
        }
      break;
    case 2:                    // DIVERGENCE DIRECTION Z
      for (register int i = 1; i < nxn - 2; i++)
        for (register int j = 1; j < nyn - 2; j++) {
          compX = .25 * (vectorX[i + 1][j][rightActiveNode] - vectorX[i][j][rightActiveNode]) * invdx + .25 * (vectorX[i + 1][j][rightActiveNode - 1] - vectorX[i][j][rightActiveNode - 1]) * invdx + .25 * (vectorX[i + 1][j + 1][rightActiveNode] - vectorX[i][j + 1][rightActiveNode]) * invdx + .25 * (vectorX[i + 1][j + 1][rightActiveNode - 1] - vectorX[i][j + 1][rightActiveNode - 1]) * invdx;
          compY = .25 * (vectorY[i][j + 1][rightActiveNode] - vectorY[i][j][rightActiveNode]) * invdy + .25 * (vectorY[i][j + 1][rightActiveNode - 1] - vectorY[i][j][rightActiveNode - 1]) * invdy + .25 * (vectorY[i + 1][j + 1][rightActiveNode] - vectorY[i + 1][j][rightActiveNode]) * invdy + .25 * (vectorY[i + 1][j + 1][rightActiveNode - 1] - vectorY[i + 1][j][rightActiveNode - 1]) * invdy;
          compZ = .25 * (vectorZ[i][j][rightActiveNode] - vectorZ[i][j][rightActiveNode - 1]) * invdz + .25 * (vectorZ[i + 1][j][rightActiveNode] - vectorZ[i + 1][j][rightActiveNode - 1]) * invdz + .25 * (vectorZ[i][j + 1][rightActiveNode] - vectorZ[i][j + 1][rightActiveNode - 1]) * invdz + .25 * (vectorZ[i + 1][j + 1][rightActiveNode] - vectorZ[i + 1][j + 1][rightActiveNode - 1]) * invdz;
          divBC[i][j][rightActiveNode] = compX + compY + compZ;
        }
      break;

  }


}

/** calculate derivative on left boundary */
void Grid3DCU::derBC(arr3_double derBC, const_arr3_double vector, int leftActiveNode, int dirDER)const
{
  switch (dirDER) {
    case 0:                    // DERIVATIVE DIRECTION X
      for (register int j = 1; j < nyc - 1; j++)
        for (register int k = 1; k < nzc - 1; k++)
          derBC[leftActiveNode][j][k] = .25 * (vector[leftActiveNode + 1][j][k] - vector[leftActiveNode][j][k]) * invdx + .25 * (vector[leftActiveNode + 1][j][k + 1] - vector[leftActiveNode][j][k + 1]) * invdx + .25 * (vector[leftActiveNode + 1][j + 1][k] - vector[leftActiveNode][j + 1][k]) * invdx + .25 * (vector[leftActiveNode + 1][j + 1][k + 1] - vector[leftActiveNode][j + 1][k + 1]) * invdx;;

      break;
    case 1:                    // DIVERGENCE DIRECTION Y
      for (register int i = 1; i < nxc - 1; i++)
        for (register int k = 1; k < nzc - 1; k++)
          derBC[i][leftActiveNode][k] = .25 * (vector[i][leftActiveNode + 1][k] - vector[i][leftActiveNode][k]) * invdy + .25 * (vector[i][leftActiveNode + 1][k + 1] - vector[i][leftActiveNode][k + 1]) * invdy + .25 * (vector[i + 1][leftActiveNode + 1][k] - vector[i + 1][leftActiveNode][k]) * invdy + .25 * (vector[i + 1][leftActiveNode + 1][k + 1] - vector[i + 1][leftActiveNode][k + 1]) * invdy;
      break;
    case 2:                    // DIVERGENCE DIRECTION Z
      for (register int i = 1; i < nxc - 1; i++)
        for (register int j = 1; j < nyc - 1; j++)
          derBC[i][j][leftActiveNode] = .25 * (vector[i][j][leftActiveNode + 1] - vector[i][j][leftActiveNode]) * invdz + .25 * (vector[i + 1][j][leftActiveNode + 1] - vector[i + 1][j][leftActiveNode]) * invdz + .25 * (vector[i][j + 1][leftActiveNode + 1] - vector[i][j + 1][leftActiveNode]) * invdz + .25 * (vector[i + 1][j + 1][leftActiveNode + 1] - vector[i + 1][j + 1][leftActiveNode]) * invdz;
      break;

  }
}

/** interpolate on nodes from central points: do this for the magnetic field*/
void Grid3DCU::interpC2N(arr3_double vecFieldN, const_arr3_double vecFieldC)const
{
  for (register int i = 1; i < nxn - 1; i++)
    for (register int j = 1; j < nyn - 1; j++)
      for (register int k = 1; k < nzn - 1; k++)
        vecFieldN[i][j][k] = .125 * (vecFieldC[i][j][k] + vecFieldC[i - 1][j][k] + vecFieldC[i][j - 1][k] + vecFieldC[i][j][k - 1] + vecFieldC[i - 1][j - 1][k] + vecFieldC[i - 1][j][k - 1] + vecFieldC[i][j - 1][k - 1] + vecFieldC[i - 1][j - 1][k - 1]);
}

/** interpolate on central points from nodes */
void Grid3DCU::interpN2C(arr3_double vecFieldC, const_arr3_double vecFieldN)const
{
  for (register int i = 1; i < nxc - 1; i++)
    for (register int j = 1; j < nyc - 1; j++)
      for (register int k = 1; k < nzc - 1; k++)
        vecFieldC[i][j][k] = .125 * (vecFieldN[i][j][k] + vecFieldN[i + 1][j][k] + vecFieldN[i][j + 1][k] + vecFieldN[i][j][k + 1] + vecFieldN[i + 1][j + 1][k] + vecFieldN[i + 1][j][k + 1] + vecFieldN[i][j + 1][k + 1] + vecFieldN[i + 1][j + 1][k + 1]);
}

/** interpolate on central points from nodes */
void Grid3DCU::interpN2C(arr4_double vecFieldC, int ns, const_arr4_double vecFieldN)const
{
  for (register int i = 1; i < nxc - 1; i++)
    for (register int j = 1; j < nyc - 1; j++)
      for (register int k = 1; k < nzc - 1; k++)
        vecFieldC[ns][i][j][k] = .125 * (vecFieldN[ns][i][j][k] + vecFieldN[ns][i + 1][j][k] + vecFieldN[ns][i][j + 1][k] + vecFieldN[ns][i][j][k + 1] + vecFieldN[ns][i + 1][j + 1][k] + vecFieldN[ns][i + 1][j][k + 1] + vecFieldN[ns][i][j + 1][k + 1] + vecFieldN[ns][i + 1][j + 1][k + 1]);
}


