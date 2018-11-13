/*****************************************************************************
 *
 *  io_harness.c
 *
 *  Drivers for serial/parallel IO for lattice quantities.
 *
 *  Each quantity (e.g., distributions, order parameter) stored on
 *  the lattice should set up an io_info struct which tells the
 *  io_harness how to actually do the read/write.
 *
 *  Actual read and writes can be initiated with a call io_read() or
 *  io_write() with the appropriate io_info struct.
 *
 *  Parallel IO takes place by taking a Cartesian decomposition of
 *  the system which can be the same, or coarser than that of the
 *  lattice Cartesian communicator. Each IO communicator group so
 *  defined then deals with its own file.
 *
 *  $Id: io_harness.c 2799 2015-12-23 15:32:58Z stratford $
 *
 *  Edinburgh Soft Matter and Statistical Physics Group and
 *  Edinburgh Parallel Computing Centre
 *
 *  Kevin Stratford (kevin@epcc.ed.ac.uk)
 *  (c) 2007-2014 The University of Edinburgh
 *
 *****************************************************************************/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "pe.h"
#include "util.h"
#include "coords.h"
#include "leesedwards.h"
#include "io_harness.h"

struct io_decomposition_t {
  int n_io;         /* Total number of I/O groups (files) in decomposition */
  int index;        /* Index of this I/O group {0, 1, ...} */
  MPI_Comm comm;    /* MPI communicator for this group */
  int rank;         /* Rank of this process in communicator */
  int size;         /* Size of this group in processes */
  int ngroup[3];    /* Global I/O group topology XYZ */
  int coords[3];    /* Coordinates of this group in I/O topology XYZ */
  int nsite[3];     /* Size of file in lattice sites */
  int offset[3];    /* Offset of the file on the lattice */
};

struct io_info_s {
  struct io_decomposition_t * io_comm;
  size_t bytesize;
  int metadata_written;
  int processor_independent;
  int single_file_read;
  char metadata_stub[FILENAME_MAX];
  char name[FILENAME_MAX];

  int (* write_function)   (FILE *, const int, const int, const int);
  int (* read_function)    (FILE *, const int, const int, const int);
  int (* write_function_a) (FILE *, const int, const int, const int);
  int (* read_function_a)  (FILE *, const int, const int, const int);
  int (* write_function_b) (FILE *, const int, const int, const int);
  int (* read_function_b)  (FILE *, const int, const int, const int);
  io_rw_cb_ft write_data;
  io_rw_cb_ft write_ascii;
  io_rw_cb_ft write_binary;
  io_rw_cb_ft read_data;
  io_rw_cb_ft read_ascii;
  io_rw_cb_ft read_binary;
};

static io_info_t * io_info_allocate(void);
static void io_set_group_filename(char *, const char *, io_info_t *);
static long int io_file_offset(int, int, io_info_t *);
static struct io_decomposition_t * io_decomposition_allocate(void);
static struct io_decomposition_t * io_decomposition_create(const int grid[3]);
static void io_decomposition_destroy(struct io_decomposition_t *);

/*****************************************************************************
 *
 *  io_info_create
 *
 *  Return a pointer to an io_info object with default decomposition.
 *
 *****************************************************************************/

io_info_t * io_info_create() {

  int io_grid[3] = {1, 1, 1}; /* Default i/o grid */
  io_info_t * p_info;

  p_info = io_info_create_with_grid(io_grid);

  return p_info;
}

/*****************************************************************************
 *
 *  io_info_create_with_grid
 *
 *  Retrun a pointer to a new io_info object with specified decomposition.
 *
 *****************************************************************************/

io_info_t * io_info_create_with_grid(const int grid[3]) {

  io_info_t * p_info;
  struct io_decomposition_t * p_decomp;

  p_info = io_info_allocate();
  p_decomp = io_decomposition_create(grid);
  p_info->io_comm = p_decomp;
  io_info_set_processor_dependent(p_info);
  p_info->single_file_read = 0;

  return p_info;
}

/*****************************************************************************
 *
 *  io_decomposition_create
 *
 *  Set up an io_decomposition object using its size.
 *
 *****************************************************************************/

static struct io_decomposition_t * io_decomposition_create(const int grid[3]) {

  int i, colour;
  int noffset[3];
  struct io_decomposition_t * p = NULL;
  MPI_Comm comm = cart_comm();

  assert(comm != MPI_COMM_NULL);
  coords_nlocal_offset(noffset);

  p = io_decomposition_allocate();
  p->n_io = 1;

  for (i = 0; i < 3; i++) {
    if (cart_size(i) % grid[i] != 0) fatal("Bad I/O grid (dim %d)\n", i);
    p->ngroup[i] = grid[i];
    p->n_io *= grid[i];
    p->coords[i] = grid[i]*cart_coords(i)/cart_size(i);
    p->nsite[i] = N_total(i)/grid[i];
    p->offset[i] = noffset[i] - p->coords[i]*p->nsite[i];
  }

  colour = p->coords[X]
         + p->coords[Y]*grid[X]
         + p->coords[Z]*grid[X]*grid[Y];

  p->index = colour;

  MPI_Comm_split(comm, colour, cart_rank(), &p->comm);
  MPI_Comm_rank(p->comm, &p->rank);
  MPI_Comm_size(p->comm, &p->size);

  return p;
}

/*****************************************************************************
 *
 *  io_set_group_filename
 *
 *  Build the file name for this I/O group from the stub provided.
 *
 *****************************************************************************/

static void io_set_group_filename(char * filename_io, const char * stub,
				  io_info_t * info) {

  assert(stub);
  assert(strlen(stub) < FILENAME_MAX/2);  /* stub should not be too long */
  assert(info);
  assert(info->io_comm);
  assert(info->io_comm->n_io < 1000);     /* format restriction ... */


  sprintf(filename_io, "%s.%3.3d-%3.3d", stub, info->io_comm->n_io,
	  info->io_comm->index + 1);

  if (info->single_file_read) {
    sprintf(filename_io, "%s.%3.3d-%3.3d", stub, 1, 1);
  }

  return;
}

/*****************************************************************************
 *
 *  io_decomposition_allocate
 *
 *  Allocate an io_decomposition_t object or fail gracefully.
 *
 *****************************************************************************/

static struct io_decomposition_t * io_decomposition_allocate() {

  struct io_decomposition_t * p = NULL;

  p = (struct io_decomposition_t*) calloc(1, sizeof(struct io_decomposition_t));
  if (p == NULL) fatal("Failed to allocate io_decomposition_t\n");

  return p;
}

/*****************************************************************************
 *
 *  io_decomposition_destroy
 *
 *****************************************************************************/

static void io_decomposition_destroy(struct io_decomposition_t * p) {

  assert(p);
  MPI_Comm_free(&p->comm);
  free(p);
 
  return;
}

/*****************************************************************************
 *
 *  io_info_allocate
 *
 *  Return a pointer to newly allocated io_info object.
 *
 *****************************************************************************/

io_info_t * io_info_allocate() {

  io_info_t * p = NULL;

  p = (io_info_t*) calloc(1, sizeof(io_info_t));
  if (p == NULL) fatal("Failed to allocate io_info_t struct\n");

  return p;
}

/*****************************************************************************
 *
 *  io_info_destroy
 *
 *  Deallocate io_info_t struct.
 *
 *****************************************************************************/

void io_info_destroy(io_info_t * p) {

  assert(p != (io_info_t *) NULL);
  io_decomposition_destroy(p->io_comm);
  free(p);

  return;
}

/*****************************************************************************
 *
 *  io_info_set_write
 *
 *****************************************************************************/

void io_info_set_write(io_info_t * p,
		       int (* writer) (FILE *, int, int, int)) {

  assert(p != (io_info_t *) NULL);
  p->write_function = writer;

  return;
}

/*****************************************************************************
 *
 *  io_info_set_read
 *
 *****************************************************************************/

void io_info_set_read(io_info_t * p,
		      int (* reader) (FILE *, int, int, int)) {

  assert(p);
  p->read_function = reader;

  return;
}

/*****************************************************************************
 *
 *  io_info_set_name
 *
 *****************************************************************************/

void io_info_set_name(io_info_t * p, const char * name) {

  assert(p);
  assert(strlen(name) < FILENAME_MAX);
  strcpy(p->name, name);

  return;
}

/*****************************************************************************
 *
 *  io_info_set_processor_dependent
 *
 *****************************************************************************/

void io_info_set_processor_dependent(io_info_t * p) {

  assert(p);
  p->processor_independent = 0;

  return;
}

/*****************************************************************************
 *
 *  io_info_set_processor_independent
 *
 *****************************************************************************/

void io_info_set_processor_independent(io_info_t * p) {

  assert(p);
  p->processor_independent = 1;

  return;
}

/*****************************************************************************
 *
 *  io_info_set_bytesize
 *
 *****************************************************************************/

void io_info_set_bytesize(io_info_t * p, size_t size) {

  assert(p);
  p->bytesize = size;

  return;
}

/*****************************************************************************
 *
 *  io_file_offset
 *
 *  Compute the file offset required for processor decomposition
 *  indepenedent files.
 *
 *****************************************************************************/

static long int io_file_offset(int ic, int jc, io_info_t * info) {

  long int offset;
  int ifo, jfo, kfo;
  int noffset[3];

  assert(info);

  /* Work out the offset of local lattice site (ic, jc, kc=1) in the file */
  ifo = info->io_comm->offset[X] + ic - 1;
  jfo = info->io_comm->offset[Y] + jc - 1;
  kfo = info->io_comm->offset[Z];

  offset = (ifo*info->io_comm->nsite[Y]*info->io_comm->nsite[Z]
	  + jfo*info->io_comm->nsite[Z]
	  + kfo)*info->bytesize;

  /* Single file offset */

  if (info->single_file_read) {
    coords_nlocal_offset(noffset);
    ifo = noffset[X] + ic - 1;
    jfo = noffset[Y] + jc - 1;
    kfo = noffset[Z];
    offset = info->bytesize*(ifo*N_total(Y)*N_total(Z) + jfo*N_total(Z) + kfo);
  }

  return offset;
}

/*****************************************************************************
 *
 *  io_write_metadata
 *
 *****************************************************************************/

int io_write_metadata(io_info_t * info) {

  assert(info);

  io_write_metadata_file(info, info->metadata_stub);

  return 0;
}

/*****************************************************************************
 *
 *  io_write_metadata_file
 *
 *  This describes, in human-readable form, the contents of the set
 *  of 1 or more files produced by a call to io_write.
 *
 *****************************************************************************/

int io_write_metadata_file(io_info_t * info, char * filename_stub) {

  FILE * fp_meta;
  char filename_io[FILENAME_MAX];
  char subdirectory[FILENAME_MAX];
  char filename[FILENAME_MAX];
  int  nx, ny, nz;
  int n[3], noff[3];

  int token = 0;
  const int tag = 1293;
  MPI_Status status;

  /* Every group writes a file, ie., the information stub and
   * the details of the local group which allow the output to
   * be unmangled. */

  assert(info);
  coords_nlocal_offset(noff);

  pe_subdirectory(subdirectory);

  io_set_group_filename(filename, filename_stub, info);
  sprintf(filename_io, "%s%s.meta", subdirectory, filename);

  if (info->io_comm->rank == 0) {
    /* Write the information stub */

    nx = info->io_comm->ngroup[X];
    ny = info->io_comm->ngroup[Y];
    nz = info->io_comm->ngroup[Z];

    fp_meta = fopen(filename_io, "w");
    if (fp_meta == NULL) fatal("fopen(%s) failed\n", filename_io);

    fprintf(fp_meta, "Metadata for file set prefix:    %s\n", filename_stub);
    fprintf(fp_meta, "Data description:                %s\n", info->name);
    fprintf(fp_meta, "Data size per site (bytes):      %d\n",
	    (int) info->bytesize);
    fprintf(fp_meta, "is_bigendian():                  %d\n", is_bigendian());
    fprintf(fp_meta, "Number of processors:            %d\n", pe_size());
    fprintf(fp_meta, "Cartesian communicator topology: %d %d %d\n",
	    cart_size(X), cart_size(Y), cart_size(Z));
    fprintf(fp_meta, "Total system size:               %d %d %d\n",
	    N_total(X), N_total(Y), N_total(Z));
    /* Lees Edwards hardwired until refactor LE code dependencies */
    fprintf(fp_meta, "Lees-Edwards planes:             %d\n",
	    le_get_nplane_total());
    fprintf(fp_meta, "Lees-Edwards plane speed         %16.14f\n",
	    le_plane_uy_max());
    fprintf(fp_meta, "Number of I/O groups (files):    %d\n", nx*ny*nz);
    fprintf(fp_meta, "I/O communicator topology:       %d %d %d\n",
	    nx, ny, nz);
    fprintf(fp_meta, "Write order:\n");

  }
  else {
    MPI_Recv(&token, 1, MPI_INT, info->io_comm->rank - 1, tag,
	     info->io_comm->comm, &status);
    fp_meta = fopen(filename_io, "a");
    if (fp_meta == NULL) fatal("fopen(%s) failed\n", filename_io);
  }

  /* Local decomposition information */

  coords_nlocal(n);
  fprintf(fp_meta, "%3d %3d %3d %3d %d %d %d %d %d %d\n", info->io_comm->rank,
          cart_coords(X), cart_coords(Y), cart_coords(Z),
          n[X], n[Y], n[Z], noff[X], noff[Y], noff[Z]);

  if (ferror(fp_meta)) {
    perror("perror: ");
    fatal("File error on writing %s\n", filename_io);
  }
  fclose(fp_meta);

 if(info->io_comm->rank < info->io_comm->size - 1) {
   MPI_Ssend(&token, 1, MPI_INT, info->io_comm->rank + 1, tag,
	     info->io_comm->comm);
 }

 info->metadata_written = 1;

  return 0;
}

/*****************************************************************************
 *
 *  io_remove_metadata
 *
 *  Largely to clean up after automated tests; usually want to keep!
 *
 *****************************************************************************/

int io_remove_metadata(io_info_t * obj, const char * file_stub) {

  char subdirectory[FILENAME_MAX];
  char filename[FILENAME_MAX];
  char filename_io[FILENAME_MAX];

  assert(obj);
  assert(file_stub);

  if (obj->io_comm->rank == 0) {
    pe_subdirectory(subdirectory);
    io_set_group_filename(filename, file_stub, obj);
    sprintf(filename_io, "%s%s.meta", subdirectory, filename);
    remove(filename_io);
  }

  return 0;
}

/*****************************************************************************
 *
 *  io_remove
 *
 *  Remove filename on each IO root.
 *
 *****************************************************************************/

int io_remove(char * filename_stub, io_info_t * obj) {

  char subdirectory[FILENAME_MAX];
  char filename[FILENAME_MAX];

  assert(filename_stub);
  assert(obj);

  if (obj->io_comm->rank == 0) {
    pe_subdirectory(subdirectory);
    io_set_group_filename(filename, filename_stub, obj);
    remove(filename);
  }

  return 0;
}

/*****************************************************************************
 *
 *  io_info_format_set
 *
 *****************************************************************************/

int io_info_format_set(io_info_t * obj, int form_in, int form_out) {

  assert(obj);
  assert(form_in >= 0);
  assert(form_in <= IO_FORMAT_DEFAULT);
  assert(form_out >= 0);
  assert(form_out <= IO_FORMAT_DEFAULT);

  io_info_format_in_set(obj, form_in);
  io_info_format_out_set(obj, form_out);

  return 0;
}

/*****************************************************************************
 *
 *  io_info_format_in_set
 *
 *  Set input format.
 *
 *****************************************************************************/

int io_info_format_in_set(io_info_t * obj, int form_in) {

  assert(obj);
  assert(form_in >= 0);
  assert(form_in <= IO_FORMAT_DEFAULT);

  if (form_in == IO_FORMAT_NULL) return 0;

  switch (form_in) {
  case IO_FORMAT_ASCII_SERIAL:
    obj->read_data = obj->read_ascii;
    obj->processor_independent = 1;
    break;
  case IO_FORMAT_BINARY_SERIAL:
    obj->read_data = obj->read_binary;
    obj->processor_independent = 1;
    break;
  case IO_FORMAT_ASCII:
    obj->read_data = obj->read_ascii;
    obj->processor_independent = 0;
    break;
  case IO_FORMAT_BINARY:
  case IO_FORMAT_DEFAULT:
    obj->read_data = obj->read_binary;
    obj->processor_independent = 0;
    break;
  default:
    fatal("Bad i/o input format\n");
  }

  return 0;
}

/*****************************************************************************
 *
 *  io_info_format_out_set
 *
 *  No serial output at the moment (unless one MPI task).
 *
 *****************************************************************************/

int io_info_format_out_set(io_info_t * obj, int form_out) {

  assert(obj);
  assert(form_out >= 0);
  assert(form_out <= IO_FORMAT_DEFAULT);

  if (form_out == IO_FORMAT_NULL) return 0;

  switch (form_out) {
  case IO_FORMAT_ASCII:
    obj->write_data = obj->write_ascii;
    obj->processor_independent = 0;
    break;
  case IO_FORMAT_BINARY:
  case IO_FORMAT_DEFAULT:
    obj->write_data = obj->write_binary;
    obj->processor_independent = 0;
    break;
  default:
    fatal("Bad i/o output format\n");
  }

  return 0;
}

/*****************************************************************************
 *
 *  io_info_read_set
 *
 *****************************************************************************/

int io_info_read_set(io_info_t * obj, int format, io_rw_cb_ft f) {

  assert(obj);
  assert(format == IO_FORMAT_ASCII || format == IO_FORMAT_BINARY);
  assert(f);

  if (format == IO_FORMAT_ASCII) obj->read_ascii = f;
  if (format == IO_FORMAT_BINARY) obj->read_binary = f;

  return 0;
}

/*****************************************************************************
 *
 *  io_info_write_set
 *
 *****************************************************************************/

int io_info_write_set(io_info_t * obj, int format, io_rw_cb_ft f) {

  assert(obj);
  assert(format == IO_FORMAT_ASCII || format == IO_FORMAT_BINARY);
  assert(f);

  if (format == IO_FORMAT_ASCII) obj->write_ascii = f;
  if (format == IO_FORMAT_BINARY) obj->write_binary = f;

  return 0;
}

/*****************************************************************************
 *
 *  io_write
 *
 *  This is the driver to write lattice quantities on the lattice.
 *  The arguments are the filename stub and the io_info struct
 *  describing which quantity we are dealing with.
 *
 *  The third argument is an opaque pointer to the data object,
 *  which will be passed to the callback which does the write.
 *
 *  All writes are processor decomposition dependent at the moment.
 *
 *****************************************************************************/

int io_write_data(io_info_t * obj, const char * filename_stub, void * data) {

  FILE *    fp_state;
  char      filename_io[FILENAME_MAX];
  int       token = 0;
  int       ic, jc, kc, index;
  int       nlocal[3];
  const int io_tag = 140;

  MPI_Status status;

  assert(obj);
  assert(data);
  assert(obj->write_data);

  if (obj->metadata_written == 0) io_write_metadata(obj);

  coords_nlocal(nlocal);
  io_set_group_filename(filename_io, filename_stub, obj);

  if (obj->io_comm->rank == 0) {
    /* Open the file anew */
    fp_state = fopen(filename_io, "wb");
  }
  else {

    /* Non-io-root process. Block until we receive the token from the
     * previous process, after which we can re-open the file and write
     * our own data. */

    MPI_Recv(&token, 1, MPI_INT, obj->io_comm->rank - 1, io_tag,
	     obj->io_comm->comm, &status);
    fp_state = fopen(filename_io, "ab");
  }

  if (fp_state == NULL) fatal("Failed to open %s\n", filename_io);

  for (ic = 1; ic <= nlocal[X]; ic++) {
    for (jc = 1; jc <= nlocal[Y]; jc++) {
      for (kc = 1; kc <= nlocal[Z]; kc++) {
	index = coords_index(ic, jc, kc);
	obj->write_data(fp_state, index, data);
      }
    }
  }

  /* Check the error indicator on the stream and close */

  if (ferror(fp_state)) {
    perror("perror: ");
    fatal("File error on writing %s\n", filename_io);
  }
  fclose(fp_state);

  /* Pass the token to the next process to write */

  if (obj->io_comm->rank < obj->io_comm->size - 1) {
    MPI_Ssend(&token, 1, MPI_INT, obj->io_comm->rank + 1, io_tag,
	      obj->io_comm->comm);
  }

  return 0;
}

/*****************************************************************************
 *
 *  io_read_data
 *
 *  Driver for reads.
 *
 *****************************************************************************/

int io_read_data(io_info_t * obj, const char * filename_stub, void * data) {

  FILE *    fp_state;
  char      filename_io[FILENAME_MAX];
  long int  token = 0;
  int       ic, jc, kc, index;
  int       nlocal[3];
  long int  offset;
  const int io_tag = 141;

  MPI_Status status;

  assert(obj);
  assert(filename_stub);
  assert(data);

  coords_nlocal(nlocal);

  io_set_group_filename(filename_io, filename_stub, obj);

  if (obj->io_comm->rank == 0) {

    fp_state = fopen(filename_io, "r");
  }
  else {

    /* Non-io-root process. Block until we receive the token from the
     * previous process, after which we can re-open the file and read. */

    MPI_Recv(&token, 1, MPI_LONG, obj->io_comm->rank - 1, io_tag,
	     obj->io_comm->comm, &status);
    fp_state = fopen(filename_io, "r");
  }

  if (fp_state == NULL) fatal("Failed to open %s\n", filename_io);
  fseek(fp_state, token, SEEK_SET);

  for (ic = 1; ic <= nlocal[X]; ic++) {
    for (jc = 1; jc <= nlocal[Y]; jc++) {

      /* Work out where the read comes from if required */
      offset = io_file_offset(ic, jc, obj);
      if (obj->processor_independent) fseek(fp_state, offset, SEEK_SET);

      for (kc = 1; kc <= nlocal[Z]; kc++) {
	index = coords_index(ic, jc, kc);
	obj->read_data(fp_state, index, data);
      }
    }
  }

  /* The token is the current offset for processor-dependent output */

  token = ftell(fp_state);

  /* Check the error indicator on the stream and close */

  if (ferror(fp_state)) {
    perror("perror: ");
    fatal("File error on reading %s\n", filename_io);
  }
  fclose(fp_state);

  /* Pass the token to the next process to read */

  if (obj->io_comm->rank < obj->io_comm->size - 1) {
    MPI_Ssend(&token, 1, MPI_LONG, obj->io_comm->rank + 1, io_tag,
	      obj->io_comm->comm);
  }

  return 0;
}

/*****************************************************************************
 *
 *  io_info_single_file_set
 *
 *****************************************************************************/

void io_info_single_file_set(io_info_t * info) {

  assert(info);

  info->single_file_read = 1;

  return;
}

/*****************************************************************************
 *
 *  io_info_metadata_filestub_set
 *
 *****************************************************************************/

int io_info_metadata_filestub_set(io_info_t * info, const char * stub) {

  assert(info);
  assert(stub);
  assert(strlen(stub) < FILENAME_MAX);

  strncpy(info->metadata_stub, stub, FILENAME_MAX);

  return 0;
}
