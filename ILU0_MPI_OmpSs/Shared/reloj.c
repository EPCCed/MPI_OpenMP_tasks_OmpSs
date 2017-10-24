/*
 * Copyright (C) 2016 UJI
 *
 * This library is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation; either version 2.1 of the License, or (at
 * your option) any later version.
 *
 * This library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * See the GNU Lesser General Public License in COPYING.LGPL for more details.
 */
#include <sys/time.h>
#include <sys/types.h>
#include <sys/times.h>
#include <unistd.h>

static double  timetick;
static double  tstart = 0.0;
static double  ucpustart = 0.0;
static int     first = 1;

void reloj (double *elapsed, double *ucpu)
{
  struct tms      cpu;
  struct timeval  tp;
  
  if(first) {

    /* Initialize clock */
    timetick = 1.0 / (double)(sysconf(_SC_CLK_TCK));
    first = 0;
    gettimeofday(&tp, NULL); 
    tstart = (double)tp.tv_sec + (double)tp.tv_usec * 1.0e-6;

    /* Initialize CPU time */
    times(&cpu);
    ucpustart = (double)(cpu.tms_utime + cpu.tms_cutime) * timetick;

    /* Return values */
    *elapsed = 0.0e0;
    *ucpu = 0.0e0;

  }
  else  {

    /* Get clock time */
    gettimeofday(&tp, NULL); 
    *elapsed = (double)tp.tv_sec + (double)tp.tv_usec * 1.0e-6 - tstart;

    /* Get CPU time */
    times(&cpu);
    *ucpu = (double)(cpu.tms_utime + cpu.tms_cutime) * timetick - ucpustart;
  }
  
  return;

}

