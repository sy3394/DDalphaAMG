/*
 * Copyright (C) 2016, Matthias Rottmann, Artur Strebel, Simon Heybrock, Simone Bacchio, Bjoern Leder.
 * 
 * This file is part of the DDalphaAMG solver library.
 * 
 * The DDalphaAMG solver library is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * The DDalphaAMG solver library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * 
 * You should have received a copy of the GNU General Public License
 * along with the DDalphaAMG solver library. If not, see http://www.gnu.org/licenses/.
 * copied:11/30/2019
 * changed from sbacchio: take out SSE
 * glanced over: 12/08/2019
 */

/*
  Hyperthreading is not impremented.
 */

#include "main.h"
#include <omp.h>


void no_barrier(int id)
{
}

void no_hyperthread_barrier(void *barrier, int id)
{
}

void core_barrier(int core)
{
#ifdef OPENMP
#pragma omp barrier
#endif
}

void hyperthread_barrier(void *barrier, int hyperthead)
{
    // no hyperthreads for now => do nothing
}


void init_common_thread_data(struct common_thread_data *common)
{
    common->barrier = &core_barrier;
    common->thread_barrier = &hyperthread_barrier;
    common->workspace = NULL;
    MALLOC( common->workspace, char, 4*128*sizeof(double)); //4=#cores???? 128 doubles = 64 complex_double
}


void setup_threading(struct Thread *threading, struct common_thread_data *common, struct level_struct *l)
{
    // no hyperthreading for now
    setup_threading_external(threading, common, l, omp_get_num_threads(), 1, omp_get_thread_num(), 0);
}


void setup_threading_external(struct Thread *threading, struct common_thread_data *common, struct level_struct *l,
        int n_core, int n_thread, int core, int thread)
{
    threading->n_core = n_core;
    threading->n_thread = n_thread;
    threading->core = core;
    threading->thread = thread;
    threading->thread_barrier_data = 0;

    update_threading(threading, l);

    threading->barrier = common->barrier;
    threading->thread_barrier = common->thread_barrier;

    threading->workspace = common->workspace;
}

void setup_no_threading(struct Thread *no_threading, struct level_struct *l) {
  /********************
   * Thread struct for no threading case
   *******************/

  no_threading->core = 0;
  no_threading->n_core = 1;
  // no hyperthreading for now
  no_threading->thread = 0;
  no_threading->n_thread = 1;
  no_threading->workspace = NULL;
  
  update_threading(no_threading, l);
  
  no_threading->barrier = &no_barrier;
  no_threading->thread_barrier = &no_hyperthread_barrier;
  
  MALLOC( no_threading->workspace, char, 4*1024*sizeof(double));//4=#cores???? 1024 doubles = 512 complex_double 
}

void update_threading(struct Thread *threading, struct level_struct *l)
{
    struct level_struct *current = l;

    while(1)
    {
        compute_core_start_end(0, current->num_inner_lattice_sites,
			       threading->start_site+current->depth, threading->end_site+current->depth, 
			       current, threading);
        threading->n_site[current->depth]      = threading->end_site[current->depth] - threading->start_site[current->depth];
        threading->start_index[current->depth] = threading->start_site[current->depth]*current->num_lattice_site_var;
        threading->end_index[current->depth]   = threading->end_site[current->depth]*current->num_lattice_site_var;
        threading->n_index[current->depth]     = threading->n_site[current->depth]*current->num_lattice_site_var;

        if(current->next_level == NULL)
            break;
        current = current->next_level;
    }
}


void compute_core_start_end(int start, int end, int *core_start, int *core_end,
        struct level_struct *l, struct Thread *threading)
{
    // due to loop unrolling in low level functions
    int min_per_core = 3*40;
//     printf0("min_per_core = %d\n", min_per_core );
    compute_core_start_end_custom(start, end, core_start, core_end, l, threading, min_per_core);
}


void compute_core_start_end_custom(int start, int end, int *core_start, int *core_end,
				   struct level_struct *l, struct Thread *threading, int granularity) {
  /*******************
   * distribute vector entries from start to end among cores
   *******************/

  if(threading->thread != 0) //????
    {
        *core_start = 0;
        *core_end = 0;
        return;
    }
    // compute start and end indices for vector functions depending on thread
    *core_start = start;
    *core_end = end;

    // custom defined minimum size per core
    int min_per_core = granularity;

    int length = end-start;
    int per_core = ceil(((double)length/min_per_core)/threading->n_core)*min_per_core;
    int cores;
    if(per_core != 0)
       cores = length/per_core;
    else
        cores = 0;
    // there might be a remainder for one extra core
    int remainder = length - cores*per_core;

    *core_start += per_core*threading->core;
    *core_end = *core_start;
    if(threading->core < cores)
        *core_end += per_core;
    else if(threading->core == cores)
        *core_end += remainder;
}


void finalize_common_thread_data( struct common_thread_data *common ) {
  FREE( common->workspace, char, 4*128*sizeof(double));
}


void finalize_no_threading( struct Thread *no_threading ) {
  FREE( no_threading->workspace, char, 4*1024*sizeof(double));
}
