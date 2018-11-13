/******************************************************************************
 *
 *  runtime.c
 *
 *  Routines which read the input file at run time and provide access
 *  to the parameters specified therein.
 *
 *  The input file is parsed to create a list of key value pairs,
 *  which are stored as strings. The routines here remain completely
 *  agnostic about the meaning of the strings. The key / value pairs
 *  should be space separated, e.g., 
 *
 *  # The temperature of the fluid is
 *  temperature 0.01
 *
 *  Lines starting with #, and blank lines, are disregarded as comments.
 *
 *  The main code can get hold of the appropriate values by querying
 *  the "database" of keys to get their corresponding value. The
 *  relevant keys and types must clearly be known, e.g.,
 *
 *  RUN_get_double_parameter("temperature", &value);
 *
 *  We demand that the keys are unique, i.e., they only appear
 *  once in the input file.
 *
 *  In parallel, the root process in pe_comm() is responsible
 *  for reading the input file, and the key value pair list is then
 *  broadcast to all other processes.
 *
 *  $Id: runtime.c,v 1.4 2010-10-15 12:40:03 kevin Exp $
 *
 *  Edinburgh Soft Matter and Statistical Physics Group and
 *  Edinburgh Parallel Computing Centre
 *
 *  Kevin Stratford (kevin@epcc.ed.ac.uk)
 *  (c) 2010-2015 The University of Edinburgh
 *
 *****************************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "pe.h"
#include "runtime.h"

#define NKEY_LENGTH 128           /* Maximum key / value string length */

static void add_key_pair(const char *, int lineno);
static void key_broadcast(int);
static int  is_valid_key_pair(const char *, int lineno);
static int  look_up_key(const char *, char *);

struct key_pair {
  char key[NKEY_LENGTH];
  int  is_active;
  int  input_line_no;
  struct key_pair * next;
};

static struct key_pair * p_keylist = NULL;


/*****************************************************************************
 *
 *  RUN_read_user_input
 *
 *  Read the input file, and construct the list of key value pairs.
 *
 *****************************************************************************/

__targetHost__ void RUN_read_input_file(const char * input_file_name) {

  FILE * fp_input;
  int    nkeys = 0;
  int    nline = 0;
  char   line[NKEY_LENGTH];

  /* Read the file and work out number of valid key lines */

  if (pe_rank() == 0) {

    fp_input = fopen(input_file_name, "r");

    if (fp_input == NULL) {
      fatal("Input file '%s' not found.\n", input_file_name);
    }
    else {

      while (fgets(line, NKEY_LENGTH, fp_input)) {
	nline += 1;
	/* Look at the line and add it if it's a key. */
	if (is_valid_key_pair(line, nline)) {
	  add_key_pair(line, nline);
	  ++nkeys;
	}
      }
    }

    fclose(fp_input);
  }

  info("Read %d user parameters from %s\n", nkeys, input_file_name);

  key_broadcast(nkeys);

  return;
}

/*****************************************************************************
 *
 *  key_broadcast
 *
 *  Make the keys available to all MPI processes. As the number of
 *  keys could be quite large, it's worth restricting this to one
 *  MPI_Bcast().
 *
 *****************************************************************************/

static void key_broadcast(int nkeys) {

  char * packed_keys;
  int n = 0;

  /* Broacdcast the number of keys and set up the message. */

  MPI_Bcast(&nkeys, 1, MPI_INT, 0, pe_comm());

  packed_keys = (char *) malloc(nkeys*NKEY_LENGTH*sizeof(char));
  if (packed_keys == NULL) fatal("malloc(packed_keys) failed\n");

  /* Pack message */

  if (pe_rank() == 0) {
    struct key_pair * p_key = p_keylist;

    while (p_key) {
      strncpy(packed_keys + n*NKEY_LENGTH, p_key->key, NKEY_LENGTH);
      ++n;
      p_key = p_key->next;
    }
  }

  MPI_Bcast(packed_keys, nkeys*NKEY_LENGTH, MPI_CHAR, 0, pe_comm());

  /* Unpack message and set up the list */

  if (pe_rank() != 0) {
    for (n = 0; n < nkeys; n++) {
      add_key_pair(packed_keys + n*NKEY_LENGTH, 0);
    }
  }

  free(packed_keys);

  return;
}

/*****************************************************************************
 *
 *  RUN_get_double_parameter
 *
 *  Query the keys for a scalar double matching te given key.
 *
 *****************************************************************************/

__targetHost__ int RUN_get_double_parameter(const char * key, double * value) {

  int key_present = 0;
  char str_value[NKEY_LENGTH];

  key_present = look_up_key(key, str_value);

  if (key_present) {
    /* Parse the value as a double */
    *value = atof(str_value);
  }

  return key_present;
}

/*****************************************************************************
 *
 *  run_get_int_parameter
 *
 *  Query the keys for a scalr int.
 *
 *****************************************************************************/

__targetHost__ int RUN_get_int_parameter(const char * key, int * value) {

  int key_present = 0;
  char str_value[NKEY_LENGTH];

  key_present = look_up_key(key, str_value);

  if (key_present) {
    /* Parse the value as integer */
    *value = atoi(str_value);
  }

  return key_present;
}

/*****************************************************************************
 *
 *  run_get_double_parameter_vector
 *
 *  Query keys for a 3-vector of double.
 *
 *****************************************************************************/

__targetHost__
int RUN_get_double_parameter_vector(const char * key, double v[]) {

  int key_present = 0;
  char str_value[NKEY_LENGTH];

  key_present = look_up_key(key, str_value);

  if (key_present) {
    /* Parse the value as a 3-vector of double */
    if (sscanf(str_value, "%lf_%lf_%lf", &v[0], &v[1], &v[2]) != 3) {
      fatal("Could not parse input key %s as double[3]\n", key);
    }
  }

  return key_present;
}

/*****************************************************************************
 *
 *  run_get_int_parameter_vector
 *
 *  Query keys for a 3-vector of int.
 *
 *****************************************************************************/

__targetHost__
int RUN_get_int_parameter_vector(const char * key, int v[]) {

  int key_present = 0;
  char str_value[NKEY_LENGTH];

  key_present = look_up_key(key, str_value);

  if (key_present) {
    /* Parse the value as a 3-vector of ints */
    if (sscanf(str_value, "%d_%d_%d", &v[0], &v[1], &v[2]) != 3) {
      fatal("Could not parse input key %s as int[3]\n", key);
    }
  }

  return key_present;
}

/*****************************************************************************
 *
 *  run_get_string_parameter
 *
 *  Query the key list for a string. Any truncation is treated as
 *  fatal to prevent problems down the line.
 *
 *****************************************************************************/

__targetHost__
int RUN_get_string_parameter(const char * key, char * value, const int len) {

  int key_present = 0;
  char str_value[NKEY_LENGTH];

  key_present = look_up_key(key, str_value);

  if (key_present) {
    /* Just copy the string across */
    if (strlen(str_value) >= len) fatal("truncated input string %s\n", key);
    strncpy(value, str_value, len);
  }

  return key_present;
}

/*****************************************************************************
 *
 *  RUN_get_active_keys
 *
 *  Count up the number of active key / value pairs.
 *
 *****************************************************************************/

__targetHost__
int RUN_get_active_keys() {

  int nkeys = 0;
  struct key_pair * p_key = p_keylist;

  while (p_key) {
    if (p_key->is_active) ++nkeys;
    p_key = p_key->next;
  }

  return nkeys;
}

/*****************************************************************************
 *
 *  is_valid_key
 *
 *  Checks a line of the input file (one string) is a valid key.
 *
 *  Invalid strings, along with comments introduced via #, and blank
 *  lines return 0.
 *
 *****************************************************************************/

static int is_valid_key_pair(const char * line, int lineno) {

  char a[NKEY_LENGTH];
  char b[NKEY_LENGTH];

  if (strncmp("#",  line, 1) == 0) return 0;
  if (strncmp("\n", line, 1) == 0) return 0;

  /* Minimal syntax checks. The user will need to sort these
   * out. */

  if (sscanf(line, "%s %s", a, b) != 2) {
    /* This does not look like a key value pair... */
    fatal("Please check input file syntax at line %d:\n %s\n", lineno, line);
  }
  else {
    /* Check against existing keys for duplicate definitions. */

    struct key_pair * p_key = p_keylist;

    while (p_key) {

      /* We must compare for exact equality against existing key. */
      sscanf(p_key->key, "%s ", b);

      if (strcmp(b, a) == 0) {
	info("At line %d: %s\n", lineno, line); 
	fatal("Duplication of parameters in input file: %s %s\n", a, b);
      }

      p_key = p_key->next;
    }
  }

  return 1;
}

/*****************************************************************************
 *
 *  add_key_pair
 *
 *  Put a new key on the list.
 *
 *****************************************************************************/

static void add_key_pair(const char * key, int lineno) {

  struct key_pair * p_new;

  p_new = (struct key_pair *) malloc(sizeof(struct key_pair));

  if (p_new == NULL) {
    fatal("malloc(key_pair) failed\n");
  }
  else {
    /* Put the new key at the head of the list. */

    strncpy(p_new->key, key, NKEY_LENGTH);
    p_new->is_active = 1;
    p_new->input_line_no = lineno;
    p_new->next = p_keylist;

    p_keylist = p_new;
  }

  return;
}

/*****************************************************************************
 *
 *  look_up_key
 *
 *  Look through the list of keys to find one matching "key"
 *  and return the corrsponding value string.
 *
 *****************************************************************************/

static int look_up_key(const char * key, char * value) {

  int key_present = 0;
  struct key_pair * p_key = p_keylist;
  char a[NKEY_LENGTH];
  char b[NKEY_LENGTH];

  while (p_key) {

    sscanf(p_key->key, "%s %s", a, b);

    if (strcmp(a, key) == 0) {
      p_key->is_active = 0;
      key_present = 1;
      strncpy(value, b, NKEY_LENGTH);
      return key_present;
    }

    p_key = p_key->next;
  }

  return key_present;
}
