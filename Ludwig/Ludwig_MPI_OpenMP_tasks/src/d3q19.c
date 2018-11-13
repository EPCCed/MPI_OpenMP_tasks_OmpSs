/*****************************************************************************
 *
 *  d3q19.c
 *
 *  D3Q19 definitions.
 *
 *  Edinburgh Soft Matter and Statistical Physics Group and
 *  Edinburgh Parallel Computing Centre
 *
 *  (c) 2008-2014 The University of Edinburgh
 *  Contributing authors:
 *    Kevin Stratford (kevin@epcc.ed.ac.uk)
 *    Ronojoy Adhikari computed this D3Q19 basis.
 *
 *****************************************************************************/

#include "pe.h"
#include "d3q19.h"

/*****************************************************************************
 *
 *  There are 19 eigenvectors:
 *
 *  rho             (eigenvector implicitly {1})
 *  rho cv[p][X]    (x-component of velocity)
 *  rho cv[p][Y]    (y-component of velocity) 
 *  rho cv[p][Z]    (z-component of velocity)
 *  q[p][X][X]      (xx component of deviatoric stress)
 *  q[p][X][Y]      (xy component of deviatoric stress)
 *  q[p][X][Z]      (xz ...
 *  q[p][Y][Y]      (yy ...
 *  q[p][Y][Z]      (yz ...
 *  q[p][Z][Z]      (zz ...
 *  chi1[p]         (1st ghost mode)
 *  chi2[p]         (2nd ghost mode)
 *  jchi1[p][X]     (x-component of ghost current chi1[p]*rho*cv[p][X])
 *  jchi1[p][Y]     (y-component of ghost current chi1[p]*rho*cv[p][Y])
 *  jchi1[p][Z]     (z-component of ghost current chi1[p]*rho*cv[p][Z])
 *  jchi2[p][X]     (x-component of ghost current chi2[p]*rho*cv[p][X])
 *  jchi2[p][Y]     (y-component of ghost current chi2[p]*rho*cv[p][Y])
 *  jchi2[p][Z]     (z-component of ghost current chi2[p]*rho*cv[p][Z])
 *  chi3[p]         (3rd ghost mode)
 *
 *
 *  We define the following:
 *
 *  cv[NVEL][3]      lattice velocities (integers)
 *  q_[NVEL][3][3]   kinetic projector c[p][i]*c[p][j] - c_s^2 d_[i][j]
 *  wv[NVEL]         quadrature weight for each velocity
 *  norm_[NVEL]      normaliser for each mode
 *
 *  ma_[NVEL][NVEL]  full matrix of eigenvectors (doubles)
 *  mi_[NVEL][NVEL]  inverse of ma_[][]
 *
 *  blocklens        reduced halo datatype blocklengths
 *  disp_fwd         reduced halo datatype displacements
 *  disp_bwd         reduced halo datatype displacements
 *
 *  Reduced halo swap information
 *  CVXBLOCK         number of separate blocks to send in x-direction
 *  CVYBLOCK         ditto                            ... y-direction
 *  CVZBLOCK         ditto                            ... z-direction
 *
 *  For each direction there is then an array of ...
 *
 *  blocklen         block lengths
 *  disp_fwd         displacements of block from start (forward direction)
 *  disp_bwd         displacements of block from start (backward direction)
 *
 *
 *****************************************************************************/


const int cv[NVEL][3] = {{ 0,  0,  0},
			 { 1,  1,  0}, { 1,  0,  1}, { 1,  0,  0},
			 { 1,  0, -1}, { 1, -1,  0}, { 0,  1,  1},
			 { 0,  1,  0}, { 0,  1, -1}, { 0,  0,  1},
			 { 0,  0, -1}, { 0, -1,  1}, { 0, -1,  0},
			 { 0, -1, -1}, {-1,  1,  0}, {-1,  0,  1},
			 {-1,  0,  0}, {-1,  0, -1}, {-1, -1,  0}};

const  double q_[NVEL][3][3] = {
  {{-1.0/3.0, 0.0, 0.0},{ 0.0,-1.0/3.0, 0.0},{ 0.0, 0.0,-1.0/3.0}},
  {{ 2.0/3.0, 1.0, 0.0},{ 1.0, 2.0/3.0, 0.0},{ 0.0, 0.0,-1.0/3.0}},
  {{ 2.0/3.0, 0.0, 1.0},{ 0.0,-1.0/3.0, 0.0},{ 1.0, 0.0, 2.0/3.0}},
  {{ 2.0/3.0, 0.0, 0.0},{ 0.0,-1.0/3.0, 0.0},{ 0.0, 0.0,-1.0/3.0}},
  {{ 2.0/3.0, 0.0,-1.0},{ 0.0,-1.0/3.0, 0.0},{-1.0, 0.0, 2.0/3.0}},
  {{ 2.0/3.0,-1.0, 0.0},{-1.0, 2.0/3.0, 0.0},{ 0.0, 0.0,-1.0/3.0}},
  {{-1.0/3.0, 0.0, 0.0},{ 0.0, 2.0/3.0, 1.0},{ 0.0, 1.0, 2.0/3.0}},
  {{-1.0/3.0, 0.0, 0.0},{ 0.0, 2.0/3.0, 0.0},{ 0.0, 0.0,-1.0/3.0}},
  {{-1.0/3.0, 0.0, 0.0},{ 0.0, 2.0/3.0,-1.0},{ 0.0,-1.0, 2.0/3.0}},
  {{-1.0/3.0, 0.0, 0.0},{ 0.0,-1.0/3.0, 0.0},{ 0.0, 0.0, 2.0/3.0}},
  {{-1.0/3.0, 0.0, 0.0},{ 0.0,-1.0/3.0, 0.0},{ 0.0, 0.0, 2.0/3.0}},
  {{-1.0/3.0, 0.0, 0.0},{ 0.0, 2.0/3.0,-1.0},{ 0.0,-1.0, 2.0/3.0}},
  {{-1.0/3.0, 0.0, 0.0},{ 0.0, 2.0/3.0, 0.0},{ 0.0, 0.0,-1.0/3.0}},
  {{-1.0/3.0, 0.0, 0.0},{ 0.0, 2.0/3.0, 1.0},{ 0.0, 1.0, 2.0/3.0}},
  {{ 2.0/3.0,-1.0, 0.0},{-1.0, 2.0/3.0, 0.0},{ 0.0, 0.0,-1.0/3.0}},
  {{ 2.0/3.0, 0.0,-1.0},{ 0.0,-1.0/3.0, 0.0},{-1.0, 0.0, 2.0/3.0}},
  {{ 2.0/3.0, 0.0, 0.0},{ 0.0,-1.0/3.0, 0.0},{ 0.0, 0.0,-1.0/3.0}},
  {{ 2.0/3.0, 0.0, 1.0},{ 0.0,-1.0/3.0, 0.0},{ 1.0, 0.0, 2.0/3.0}},
  {{ 2.0/3.0, 1.0, 0.0},{ 1.0, 2.0/3.0, 0.0},{ 0.0, 0.0,-1.0/3.0}}};


const double chi1[NVEL] = {0.0, -2.0,  1.0,  1.0, 
                                 1.0, -2.0,  1.0, 
                                 1.0,  1.0, -2.0,
                                -2.0,  1.0,  1.0,
                                 1.0, -2.0,  1.0,
                                 1.0,  1.0, -2.0};

const double chi2[NVEL] = {0.0,  0.0, -1.0,  1.0,
                                -1.0,  0.0,  1.0,
                                -1.0,  1.0,  0.0,
                                 0.0,  1.0, -1.0,
                                 1.0,  0.0, -1.0,
                                 1.0, -1.0,  0.0};

const double chi3[NVEL] = {1.0,  1.0,  1.0, -2.0,
                                 1.0,  1.0,  1.0,
                                -2.0,  1.0, -2.0,
                                -2.0,  1.0, -2.0,
                                 1.0,  1.0,  1.0,
                                -2.0,  1.0,  1.0};

#define w0 (12.0/36.0)
#define w1  (2.0/36.0)
#define w2  (1.0/36.0)

#define c0        0.0
#define c1        1.0
#define c2        2.0
#define r3   (1.0/3.0)
#define r6   (1.0/6.0)
#define t3   (2.0/3.0)
#define r2   (1.0/2.0)
#define r4   (1.0/4.0)
#define r8   (1.0/8.0)

#define wc ( 1.0/72.0)
#define wb ( 3.0/72.0)
#define wa ( 6.0/72.0)
#define t4 (16.0/72.0)
#define wd ( 1.0/48.0)
#define we ( 3.0/48.0)

const double wv[NVEL] = {w0,
			 w2, w2, w1, w2, w2, w2, w1, w2, w1,
			 w1, w2, w1, w2, w2, w2, w1, w2, w2}; 
const double norm_[NVEL] = {c1,
      3.0, 3.0, 3.0, 9.0/2.0, 9.0, 9.0, 9.0/2.0, 9.0, 9.0/2.0,
      3.0/4.0, 3.0/2.0,3.0/2.0,3.0/2.0, 9.0/4.0,9.0/2.0,9.0/2.0,9.0/2.0, 0.5};

const double ma_[NVEL][NVEL] = {
{ c1, c1, c1, c1, c1, c1, c1, c1, c1, c1, c1, c1, c1, c1, c1, c1, c1, c1, c1},
{ c0, c1, c1, c1, c1, c1, c0, c0, c0, c0, c0, c0, c0, c0,-c1,-c1,-c1,-c1,-c1},
{ c0, c1, c0, c0, c0,-c1, c1, c1, c1, c0, c0,-c1,-c1,-c1, c1, c0, c0, c0,-c1},
{ c0, c0, c1, c0,-c1, c0, c1, c0,-c1, c1,-c1, c1, c0,-c1, c0, c1, c0,-c1, c0},
{-r3, t3, t3, t3, t3, t3,-r3,-r3,-r3,-r3,-r3,-r3,-r3,-r3, t3, t3, t3, t3, t3},
{ c0, c1, c0, c0, c0,-c1, c0, c0, c0, c0, c0, c0, c0, c0,-c1, c0, c0, c0, c1},
{ c0, c0, c1, c0,-c1, c0, c0, c0, c0, c0, c0, c0, c0, c0, c0,-c1, c0, c1, c0},
{-r3, t3,-r3,-r3,-r3, t3, t3, t3, t3,-r3,-r3, t3, t3, t3, t3,-r3,-r3,-r3, t3},
{ c0, c0, c0, c0, c0, c0, c1, c0,-c1, c0, c0,-c1, c0, c1, c0, c0, c0, c0, c0},
{-r3,-r3, t3,-r3, t3,-r3, t3,-r3, t3, t3, t3, t3,-r3, t3,-r3, t3,-r3, t3,-r3},
{ c0,-c2, c1, c1, c1,-c2, c1, c1, c1,-c2,-c2, c1, c1, c1,-c2, c1, c1, c1,-c2},
{ c0,-c2, c1, c1, c1,-c2, c0, c0, c0, c0, c0, c0, c0, c0, c2,-c1,-c1,-c1, c2},
{ c0,-c2, c0, c0, c0, c2, c1, c1, c1, c0, c0,-c1,-c1,-c1,-c2, c0, c0, c0, c2},
{ c0, c0, c1, c0,-c1, c0, c1, c0,-c1,-c2, c2, c1, c0,-c1, c0, c1, c0,-c1, c0},
{ c0, c0,-c1, c1,-c1, c0, c1,-c1, c1, c0, c0, c1,-c1, c1, c0,-c1, c1,-c1, c0},
{ c0, c0,-c1, c1,-c1, c0, c0, c0, c0, c0, c0, c0, c0, c0, c0, c1,-c1, c1, c0},
{ c0, c0, c0, c0, c0, c0, c1,-c1, c1, c0, c0,-c1, c1,-c1, c0, c0, c0, c0, c0},
{ c0, c0,-c1, c0, c1, c0, c1, c0,-c1, c0, c0, c1, c0,-c1, c0,-c1, c0, c1, c0},
{ c1, c1, c1,-c2, c1, c1, c1,-c2, c1,-c2,-c2, c1,-c2, c1, c1, c1,-c2, c1, c1}
};

const double mi_[NVEL][NVEL] = {
{w0, c0, c0, c0,-r2, c0, c0,-r2, c0,-r2, c0, c0, c0, c0, c0, c0, c0, c0, r6},
{w2, wa, wa, c0, wa, r4, c0, wa, c0,-wb,-wb,-wa,-wa, c0, c0, c0, c0, c0, wc},
{w2, wa, c0, wa, wa, c0, r4,-wb, c0, wa, wd, wb, c0, wb,-we,-r8, c0,-r8, wc},
{w1, r6, c0, c0, r6, c0, c0,-wa, c0,-wa, wb, wa, c0, c0, r8, r4, c0, c0,-w1},
{w2, wa, c0,-wa, wa, c0,-r4,-wb, c0, wa, wd, wb, c0,-wb,-we,-r8, c0, r8, wc},
{w2, wa,-wa, c0, wa,-r4, c0, wa, c0,-wb,-wb,-wa, wa, c0, c0, c0, c0, c0, wc},
{w2, c0, wa, wa,-wb, c0, c0, wa, r4, wa, wd, c0, wb, wb, we, c0, r8, r8, wc},
{w1, c0, r6, c0,-wa, c0, c0, r6, c0,-wa, wb, c0, wa, c0,-r8, c0,-r4, c0,-w1},
{w2, c0, wa,-wa,-wb, c0, c0, wa,-r4, wa, wd, c0, wb,-wb, we, c0, r8,-r8, wc},
{w1, c0, c0, r6,-wa, c0, c0,-wa, c0, r6,-wa, c0, c0,-r6, c0, c0, c0, c0,-w1},
{w1, c0, c0,-r6,-wa, c0, c0,-wa, c0, r6,-wa, c0, c0, r6, c0, c0, c0, c0,-w1},
{w2, c0,-wa, wa,-wb, c0, c0, wa,-r4, wa, wd, c0,-wb, wb, we, c0,-r8, r8, wc},
{w1, c0,-r6, c0,-wa, c0, c0, r6, c0,-wa, wb, c0,-wa, c0,-r8, c0, r4, c0,-w1},
{w2, c0,-wa,-wa,-wb, c0, c0, wa, r4, wa, wd, c0,-wb,-wb, we, c0,-r8,-r8, wc},
{w2,-wa, wa, c0, wa,-r4, c0, wa, c0,-wb,-wb, wa,-wa, c0, c0, c0, c0, c0, wc},
{w2,-wa, c0, wa, wa, c0,-r4,-wb, c0, wa, wd,-wb, c0, wb,-we, r8, c0,-r8, wc},
{w1,-r6, c0, c0, r6, c0, c0,-wa, c0,-wa, wb,-wa, c0, c0, r8,-r4, c0, c0,-w1},
{w2,-wa, c0,-wa, wa, c0, r4,-wb, c0, wa, wd,-wb, c0,-wb,-we, r8, c0, r8, wc},
{w2,-wa,-wa, c0, wa, r4, c0, wa, c0,-wb,-wb, wa, wa, c0, c0, c0, c0, c0, wc}
};


const int xblocklen_cv[CVXBLOCK] = {5};
const int xdisp_fwd_cv[CVXBLOCK] = {1};
const int xdisp_bwd_cv[CVXBLOCK] = {14};

const int yblocklen_cv[CVYBLOCK] = {1, 3, 1};
const int ydisp_fwd_cv[CVYBLOCK] = {1, 6, 14};
const int ydisp_bwd_cv[CVYBLOCK] = {5, 11, 18};

const int zblocklen_cv[CVZBLOCK] = {1, 1, 1, 1, 1};
const int zdisp_fwd_cv[CVZBLOCK] = {2, 6, 9, 11, 15};
const int zdisp_bwd_cv[CVZBLOCK] = {4, 8, 10, 13, 17};
