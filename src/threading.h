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
 * checked:11/30/2019
 * not changed from sbacchio
 * glanced over: 12/08/2019
 */

#ifndef THREADING_H
   #define THREADING_H
/**********
 * Note:
 *   thread: hyperthread; core:thread
 *   hyperthreading is not impremented
 *   dynamic memory allocation is done only in the master, and 
 *     when using PUBLIC_MALLOC no_threading should be used 
 *********/
struct level_struct;


struct common_thread_data
{
    void (*barrier)(int);                // barrier among cores
    void (*thread_barrier)(void *, int); // barrier among hyperthreads on a core
    // *common* workspace for *all* threads
    // sometimes threads need to exchange data, they can use this field
    char *workspace;
};

// holds information relevant for specific core/thread
typedef struct Thread
{
  int core;     // core id
  int n_core;   // total # cores
  // for SMT/hyperthreading: threads per core (1-4 on KNC)
  int thread;   // thread id
  int n_thread; // total # threads

  /* level_struct.num_inner_lattice_sites is split among cores
     These variables define start and end site for this specific *core* (not thread)
     but num_inner_lattice_sites depends on the level.
     Use level_struct.depth as index 
     If g.num_levels > 4, this will fail.
  */
  int start_site[4];
  int end_site[4];
  int n_site[4];
  // index = site*num_lattice_site_var = inner_vector_size
  int start_index[4];
  int end_index[4];
  int n_index[4];

  void (*barrier)(int);                // barrier among cores
  void (*thread_barrier)(void *, int); // barrier among hyperthreads on a core
  void *thread_barrier_data;
  
  // *common* workspace for *all* threads
  // sometimes threads need to exchange data, they can use this
  char *workspace;
} Thread;


/* flat omp: does not distinguish between threads over cores and hyperthreads within a core */
#ifdef FLAT_OMP 

#define CORE_BARRIER(threading) \
    do { \
    threading->thread_barrier(threading->thread_barrier_data, threading->core/60); \
    if(threading->core/60 == 0) \
        threading->barrier(threading->core); \
    threading->thread_barrier(threading->thread_barrier_data, threading->core/60); \
    } while(0)
#define HYPERTHREAD_BARRIER(threading) \
    do { \
    threading->thread_barrier(threading->thread_barrier_data, threading->core/60); \
    } while(0)

/* nested omp: first splits into cores.  Then, each core splits into hyperthreads (like DD preconditioner) */
#else

#define CORE_BARRIER(threading) \
    do { \
        threading->barrier(threading->core); \
    } while(0)
#define HYPERTHREAD_BARRIER(threading) \
    do { \
    threading->thread_barrier(threading->thread_barrier_data, threading->thread); \
    } while(0)

#endif


// used within a function; returns to the calling pt if not master hyperthread after passing through the barrier
// master execute the section of the code and bumps into the barrier
#define START_UNTHREADED_FUNCTION(threading) \
    if(threading->thread != 0) \
        return; \
    CORE_BARRIER(threading); \
    if(threading->core != 0) \
    { \
        CORE_BARRIER(threading); \
        return; \
    }
#define END_UNTHREADED_FUNCTION(threading) \
    CORE_BARRIER(threading);

/* Only one thread (master) will execute the code section between
   START_LOCKED_MASTER and END_LOCKED_MASTER, and it is protected by barriers
   among cores to prevent data races */
// not syncing hyperthread?????
#define START_LOCKED_MASTER(threading) \
    if(threading->thread == 0) \
        CORE_BARRIER(threading); \
    if(threading->core + threading->thread == 0) {
#define END_LOCKED_MASTER(threading) \
    } \
    if(threading->thread == 0) \
        CORE_BARRIER(threading);

#define MASTER(threading) \
    if(threading->core + threading->thread == 0)
#define START_MASTER(threading) \
     MASTER(threading) {
#define END_MASTER(threading) \
    }

#define SYNC_MASTER_TO_ALL(threading) \
    if(threading->thread == 0) \
        CORE_BARRIER(threading); \
    HYPERTHREAD_BARRIER(threading);
    
#define SYNC_CORES(threading) \
    if(threading->thread == 0) \
        CORE_BARRIER(threading);
#define SYNC_HYPERTHREADS(threading) \
    HYPERTHREAD_BARRIER(threading);

#define START_NO_HYPERTHREADS(threading) \
    if(threading->thread == 0) {
#define END_NO_HYPERTHREADS(threading) \
    }


#ifdef OPENMP
#include <omp.h>
#define DO_PRAGMA(EXP) _Pragma (#EXP )
#define THREADED(EXP) DO_PRAGMA ( omp parallel num_threads( EXP ) )
#else
#define THREADED(EXP)
static inline int omp_get_thread_num( void ) {
  return 0;
}
static inline int omp_get_num_threads( void ) {
  return 1;
}
#endif

void init_common_thread_data(struct common_thread_data *common);

void setup_threading(struct Thread *threading, struct common_thread_data *common, struct level_struct *l);
/* external means the caller gives us all info about threads, and is responsible to later set a proper barrier */
void setup_threading_external(struct Thread *threading, struct common_thread_data *common, struct level_struct *l,
        int n_core, int n_thread, int core, int thread);
void update_threading(struct Thread *threading, struct level_struct *l);
void setup_no_threading(struct Thread *no_threading, struct level_struct *l);

/* computes start and end indices for a core inside an array
   puts zero for other hyperthreads */
void compute_core_start_end(int start, int end, int *core_start, int *core_end,
        struct level_struct *l, struct Thread *threading);
void compute_core_start_end_custom(int start, int end, int *core_start, int *core_end,
        struct level_struct *l, struct Thread *threading, int granularity);

void finalize_common_thread_data( struct common_thread_data *common );

void finalize_no_threading( struct Thread *no_threading );

extern struct Thread *no_threading;

extern int threaded;

#endif // THREADING_H
