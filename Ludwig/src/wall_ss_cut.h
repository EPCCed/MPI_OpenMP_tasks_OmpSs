/*****************************************************************************
 *
 *  wall_ss_cut.h
 *
 *  Edinburgh Soft Matter and Statistical Physics Group and
 *  Edinburgh Parallel Computing Centre
 *
 *  (c) The University of Edinburgh (2014)
 *  Juho Lintuvuori (juho.lintuvuori@u-psud.fr)
 *
 *****************************************************************************/

#ifndef WALL_SS_CUT_H
#define WALL_SS_CUT_H

#include "colloids.h"
#include "interaction.h"

typedef struct wall_ss_cut_s wall_ss_cut_t;

int wall_ss_cut_create(wall_ss_cut_t ** pobj);
void wall_ss_cut_free(wall_ss_cut_t * obj);
int wall_ss_cut_info(wall_ss_cut_t * obj);
int wall_ss_cut_register(wall_ss_cut_t * obj, interact_t * parent);
int wall_ss_cut_compute(colloids_info_t * cinfo, void * self);
int wall_ss_cut_stats(void * self, double * stats);
int wall_ss_cut_single(wall_ss_cut_t * obj, double h, double * f, double * v);
int wall_ss_cut_param_set(wall_ss_cut_t * obj, double epsilon, double sigma,
			  int nu, double hc);

#endif
