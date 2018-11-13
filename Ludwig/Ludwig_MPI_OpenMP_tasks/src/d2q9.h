/*****************************************************************************
 *
 *  d2q9.h
 *
 *  D2Q9 definitions.
 *
 *  Edinburgh Soft Matter and Statistical Physics Group and
 *  Edinburgh Parallel Computing Centre
 *
 *  Kevin Stratford (kevin@epcc.ed.ac.uk)
 *  (c) 2010-2014 The University of Edinburgh
 *
 *****************************************************************************/

#ifndef D2Q9_MODEL_H
#define D2Q9_MODEL_H

enum {NDIM = 2};
enum {NVEL = 9};
enum {CVXBLOCK = 1};
enum {CVYBLOCK = 3};
enum {CVZBLOCK = 1};

extern const    int cv[NVEL][3];
extern const double wv[NVEL];
extern const double q_[NVEL][3][3];
extern const double norm_[NVEL];
extern const double ma_[NVEL][NVEL];
extern const double mi_[NVEL][NVEL];

extern const int xblocklen_cv[CVXBLOCK];
extern const int xdisp_fwd_cv[CVXBLOCK];
extern const int xdisp_bwd_cv[CVXBLOCK];

extern const int yblocklen_cv[CVYBLOCK];
extern const int ydisp_fwd_cv[CVYBLOCK];
extern const int ydisp_bwd_cv[CVYBLOCK];

extern const int zblocklen_cv[CVZBLOCK];
extern const int zdisp_fwd_cv[CVZBLOCK];
extern const int zdisp_bwd_cv[CVZBLOCK];

#endif
