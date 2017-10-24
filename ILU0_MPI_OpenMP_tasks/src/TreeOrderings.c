//
//  main.c
//  TreeOrderings
//
//  Created by José Ignacio Aliaga Estellés on 18/5/15.
//  Copyright (c) 2015 Universitat Jaume I. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>

void LevelOrdering(int *vec, int dim, int lev) {
  int i, j, k;
  int size, acum;
  
  k = 0; acum = 0;
  for (i=0; i<lev-1; i++) {
    size = 1 << (lev-i-1);
    acum += size;
    for (j=0; j<size; j++) {
      vec[k] = acum + (j/2);
      k++;
    }
  }
  vec[k] = -1;
}

void PostOrderOrdering(int *prm, int *vec, int dim, int lev) { //vec es treetab
  int i, j, k;
  int *izq_der = NULL, *pos = NULL;
  
  izq_der = (int *) malloc(dim * sizeof (int));
  for (i=0; i<dim; i++) izq_der[i] = 0;
  pos = (int *) malloc(lev * sizeof (int));
  pos[0] = 0; for (i=1; i<lev; i++) pos[i] = pos[i-1] + (1<<(lev-i));
  k = 0;
  for (i=0; i<(1<<(lev-1)); i++) {
    j = -1;
    do {
      j++;
      if (izq_der[j] == 0)
        vec[k] = k + (1 << (j+1));
      else
        vec[k] = k + 1;
      prm[pos[j]] = k;
      izq_der[j] = 1-izq_der[j];
      pos[j]++; k++;
    } while ((j<(lev-2)) && (izq_der[j] == 0));
      
  }
  vec[k] = -1; prm[k] = pos[lev-1];
  free (pos) ; pos = NULL;
  free (izq_der) ; izq_der = NULL;
}

void InvPerm (int *prm, int *inv, int dim) {
  int i;
  
  for (i=0; i<dim; i++)
    inv[prm[i]] = i;
}
