/*****************************************************************************
 *
 *  runtime.h
 *
 *  Runtime input interface.
 *
 *  $Id: runtime.h,v 1.3 2010-10-15 12:40:03 kevin Exp $
 *
 *  Edinburgh Soft Matter and Statistical Physics Group and
 *  Edinburgh Parallel Computing Centre
 *
 *  Kevin Stratford (kevin@epcc.ed.ac.uk)
 *  (c) 2010 The University of Edinburgh
 *
 *****************************************************************************/

#ifndef RUNTIME_H
#define RUNTIME_H

__targetHost__ void RUN_read_input_file(const char *);
__targetHost__ int  RUN_get_int_parameter(const char *, int *);
__targetHost__ int  RUN_get_int_parameter_vector(const char *, int []);
__targetHost__ int  RUN_get_double_parameter(const char *, double *);
__targetHost__ int  RUN_get_double_parameter_vector(const char *, double []);
__targetHost__ int  RUN_get_string_parameter(const char *, char *, const int);
__targetHost__ int  RUN_get_active_keys(void);

#endif
