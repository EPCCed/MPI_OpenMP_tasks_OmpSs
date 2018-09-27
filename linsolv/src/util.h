/*
 * util.h
 */

#ifndef UTIL_H
#define UTIL_H

#include <stddef.h> /* size_t */
#include <stdio.h> /* printf */
#include <stdlib.h> /* EXIT_FAILURE */
#include <float.h> /* DBL_EPSILON */

#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))

#if defined(DEBUG)
#  define DBG_MSG printf
#  define DBG_PRINT 1
#else
#  define DBG_MSG(...)
#  define DBG_PRINT 0
#endif

#define STR(token) #token
#define EXPSTR(macro) STR(macro)

#if !defined(USE_TASKLOOP) || USE_TASKLOOP < 0 || USE_TASKLOOP > 1
# undef USE_TASKLOOP
# define USE_TASKLOOP 0
#endif

#if !defined(SET_CHUNKSIZE) || SET_CHUNKSIZE < 1
# undef SET_CHUNKSIZE
# define SET_CHUNKSIZE 200
# define GRAINSIZECLAUSE
#else
# define GRAINSIZECLAUSE grainsize(SET_CHUNKSIZE)
# define GRAINSIZESTRING EXPSTR(GRAINSIZECLAUSE)
#endif

#if USE_TASKLOOP == 0
# define TASKVERSION "tiled"
#elif defined(GRAINSIZESTRING)
# define TASKVERSION "taskloop "GRAINSIZESTRING
#else
# define TASKVERSION "taskloop"
#endif

#define CHECK(expr) \
  if(!(expr)) \
  { \
    printf("Error: '%s' [%s:%i]\n", #expr, __FILE__, __LINE__); \
    exit(EXIT_FAILURE); \
  }

/*******************************************************************************
*
*******************************************************************************/
void  check_free(void *ptr);
void *check_malloc(size_t bytes);
void *check_calloc(size_t number, size_t bytes);
void *check_realloc(void *old, size_t bytes);

#endif /* UTIL_H */
