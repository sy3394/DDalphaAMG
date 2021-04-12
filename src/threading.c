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
  Nested OpenMP along hyperthreading is not impremented.
  In what follows:
    * core:   thread
    * thread: hyperthread
 */

#include "main.h"
#include <omp.h>


/********  WRAPPER FUNCTIONS FOR BARRIERS  **********************************/

//----- for no threading
void no_barrier(int id)
{
}

void no_hyperthread_barrier(void *barrier, int id)
{
}

//----- for threading
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


/*******  INITIALIZATION AND SETUP OF THREADS  *****************************/

// all threads except for no_threading share some baisc info
void init_common_thread_data(struct common_thread_data *common)
{
    common->barrier        = &core_barrier;
    common->thread_barrier = &hyperthread_barrier;
    common->workspace      = NULL;
    MALLOC( common->workspace, char, 50*4*128*sizeof(double)); //4=#cores???? 50*128 doubles = 50*64 complex_double;50=max#eigvec
}

// no hyperthreading for now (=> n_thread=1&&thread=0 hardcoded): need to add the funtion
void setup_threading(struct Thread *threading, struct common_thread_data *common, struct level_struct *l)
{
    setup_threading_external(threading, common, l, omp_get_num_threads(), 1, omp_get_thread_num(), 0);
}

// use info in common_thread_data
void setup_threading_external(struct Thread *threading, struct common_thread_data *common, struct level_struct *l,
        int n_core, int n_thread, int core, int thread)
{
    threading->n_core              = n_core;
    threading->n_thread            = n_thread;
    threading->core                = core;
    threading->thread              = thread;
    threading->thread_barrier_data = 0;

    update_threading(threading, l);

    threading->barrier        = common->barrier;
    threading->thread_barrier = common->thread_barrier;

    threading->workspace      = common->workspace;
}

void setup_no_threading(struct Thread *no_threading, struct level_struct *l) {
  /********************
   * Thread struct for no threading case
   * no_threading is a global variable
   *******************/

  no_threading->core      = 0;
  no_threading->n_core    = 1;
  // no hyperthreading for now
  no_threading->thread    = 0;
  no_threading->n_thread  = 1;
  no_threading->workspace = NULL;
  
  update_threading(no_threading, l);
  
  no_threading->barrier        = &no_barrier;
  no_threading->thread_barrier = &no_hyperthread_barrier;
  
  MALLOC( no_threading->workspace, char, 4*1024*sizeof(double));//4=#cores???? 1024 doubles = 512 complex_double 
}

// Compute indicies for work distribution over threads for each level
void update_threading(struct Thread *threading, struct level_struct *l)
{
    struct level_struct *current = l;

    while(1)
    {
      // split inner lattice sites over cores; store the indicies in start_site, end_site for the given depth
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

void finalize_common_thread_data( struct common_thread_data *common ) {
  FREE( common->workspace, char, 50*4*128*sizeof(double));
}


void finalize_no_threading( struct Thread *no_threading ) {
  FREE( no_threading->workspace, char, 4*1024*sizeof(double));
}

/*********  WORK DISTRIBUTION  ************************************************/

void compute_core_start_end(int start, int end, int *core_start, int *core_end,
        struct level_struct *l, struct Thread *threading) {
  // due to loop unrolling in low level functions, we set minimum entries per core
  int min_per_core = 1;//this may be obsolete; Check the usage of FOR?? and see what's needed!!!!!!
  compute_core_start_end_custom(start, end, core_start, core_end, l, threading, min_per_core);
}


void compute_core_start_end_custom(int start, int end, int *core_start, int *core_end,
				   struct level_struct *l, struct Thread *threading, int granularity) {
  /*******************
   * distribute vector entries from start to end among cores
   *******************/

  if(threading->thread != 0) 
    {
        *core_start = 0;
        *core_end = 0;
        return;
    }
#ifdef DEBUG
  if ( (end-start)%granularity != 0 )//could rise seg fault as error'0' is used.
    error0("compute_core_start_end_custom: each core needs multiple of %d entries\n", granularity);
#endif
  
  // compute start and end indices for vector functions depending on thread
  *core_start = start;
  *core_end   = end;

  // custom defined minimum size per core: For now, this is also a unit for the number of items per core.
  int min_per_core = granularity; // why did you redefine????
  
  int length   = end-start;
  int base_per_core  = floor(((double)length/min_per_core)/threading->n_core)*min_per_core; // base items per core
  int reminder  = length-base_per_core*threading->n_core;
  int ext_cores = reminder/min_per_core;

  // TODO: allow for nonzero value for reminder%min_per_core
  // As of now, each core should recieve multiple of min_per_core many elements.
  // To generalize, we need to make a special case in loops such as for ( phi_pt=phi->vector_buffer+start*nvec_phi, end_pt=phi->vector_buffer+end*nvec_phi, D_pt = op->D+(start*3),
  // nb_pt=neighbor+((start/12)*4); phi_pt<end_pt; phi_pt+=12*nvec_phi) in d_plus_clover_PRECISION
  // Note: phi_pt is jumped by (inner d.o.f.) x nvec_phi assuming end-start%12==0 while compute_core_start_end_custom(0, nv*n, &start, &end, l, threading, nv )
#ifdef DEBUG
  if ( reminder%min_per_core != 0 )
    error0("compute_core_start_end_custom: end-start=%d needs to be divisible by min_per_core (%d)\n", length, min_per_core );
#endif
  
  if( (threading->core+1)*min_per_core <= reminder ) {// we assign extra min_per_core items on the first ext_cores many cores
    *core_start += (base_per_core+min_per_core)*threading->core;
    *core_end = *core_start + base_per_core+min_per_core;
  } else if ( threading->core == threading->n_core-1 ) {// remaining (reminder%min_per_core many) elements put in the last thread
    *core_start += base_per_core*threading->core + min_per_core*ext_cores;
    *core_end = *core_start + base_per_core + reminder%min_per_core;
  } else { // on other cores, we assign only per_core many items
    *core_start += base_per_core*threading->core + min_per_core*ext_cores;
    *core_end = *core_start+base_per_core;
  }
}

