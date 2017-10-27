/*
 * util.c
 */

#include <stdlib.h>
#include "util.h"

/*******************************************************************************
*
*******************************************************************************/
void check_free(void *ptr)
{
  if(ptr != NULL)
    free(ptr);

  return;
} /* check_free() */

/*******************************************************************************
*
*******************************************************************************/
void *check_malloc(size_t bytes)
{
  if(bytes <= (size_t) 0)
    return NULL;

  void *tmp = malloc(bytes);
  CHECK(tmp != NULL);

  return tmp;
} /* check_malloc() */

/*******************************************************************************
*
*******************************************************************************/
void *check_calloc(size_t number, size_t bytes)
{
  if(number <= (size_t) 0 || bytes <= (size_t) 0)
	return NULL;

  void *tmp = calloc(number, bytes);
  CHECK(tmp != NULL);

  return tmp;
} /* check_calloc() */

/*******************************************************************************
*
*******************************************************************************/
void *check_realloc(void *old, size_t bytes)
{

  if(bytes <= (size_t) 0)
  {
	check_free(old);
    return NULL;
  }

  void *tmp = realloc(old, bytes);
  CHECK(tmp != NULL);

  return tmp;
} /* check_realloc() */
