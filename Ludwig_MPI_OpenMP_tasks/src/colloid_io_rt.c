/*****************************************************************************
 *
 *  colloid_io_rt.c
 *
 *  Run time colloid I/O settings.
 *
 *  $Id: colloid_io_rt.c 2479 2014-09-25 11:24:42Z stratford $
 *
 *  Edinburgh Soft Matter and Statistical Physics Group and
 *  Edinburgh Parallel Computing Centre
 *
 *  Kevin Stratford (kevin@epcc.ed.ac.uk)
 *  (c) 2010 The University of Edinburgh
 *
 *****************************************************************************/

#include <assert.h>
#include <string.h>

#include "pe.h"
#include "runtime.h"
#include "colloid_io_rt.h"


/*****************************************************************************
 *
 *  colloid_io_run_time
 *
 *****************************************************************************/

int colloid_io_run_time(colloids_info_t * cinfo, colloid_io_t ** pcio) {

  int io_grid[3] = {1, 1, 1};
  char tmp[BUFSIZ];

  colloid_io_t * cio = NULL;

  assert(cinfo);

  RUN_get_int_parameter_vector("default_io_grid", io_grid);
  RUN_get_int_parameter_vector("colloid_io_grid", io_grid);

  colloid_io_create(io_grid, cinfo, &cio);
  assert(cio);

  /* Default format to ascii, parallel; then check user input */

  colloid_io_format_input_ascii_set(cio);
  colloid_io_format_output_ascii_set(cio);

  RUN_get_string_parameter("colloid_io_format", tmp, BUFSIZ);

  if (strncmp("BINARY", tmp, 5) == 0 || strncmp("binary", tmp, 5) == 0) {
    colloid_io_format_input_binary_set(cio);
    colloid_io_format_output_binary_set(cio);
  }

  RUN_get_string_parameter("colloid_io_format_input", tmp, BUFSIZ);

  if (strncmp("ASCII",  tmp, 5) == 0 || strncmp("ascii", tmp, 5) == 0) {
    colloid_io_format_input_ascii_set(cio);
  }

  if (strncmp("ASCII_SERIAL", tmp, 12) == 0 ||
      strncmp("ascii_serial", tmp, 12) == 0) {
    colloid_io_format_input_ascii_set(cio);
    colloid_io_format_input_serial_set(cio);
  }

  if (strncmp("BINARY", tmp, 6) == 0 || strncmp("binary", tmp, 6) == 0) {
    colloid_io_format_input_binary_set(cio);
  }

  if (strncmp("BINARY_SERIAL", tmp, 13) == 0 ||
      strncmp("binary_serial", tmp, 13) == 0) {
    colloid_io_format_input_binary_set(cio);
    colloid_io_format_input_serial_set(cio);
  }

  RUN_get_string_parameter("colloid_io_format_output", tmp, BUFSIZ);

  if (strncmp("ASCII",  tmp, 5) == 0 || strncmp("ascii", tmp, 5) == 0) {
    colloid_io_format_output_ascii_set(cio);
  }

  if (strncmp("BINARY", tmp, 6) == 0 || strncmp("binary", tmp, 6) == 0) {
    colloid_io_format_output_binary_set(cio);
  }

  colloid_io_info(cio);

  *pcio = cio;

  return 0;
}

