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
 * checked:11/29/2019
 * changed from sbacchio
 * glanced over: 12/18/2019
 */

#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <complex.h>
#include <math.h>
#include <time.h>
#include <stdarg.h>
#include <fabulous.h>

#ifndef MAIN_HEADER
  #define MAIN_HEADER

  #define STRINGLENGTH 500
  #define _FILE_OFFSET_BITS 64
  #define EPS_float 1E-6
  #define EPS_double 1E-14

  #define HAVE_TM       // flag for enable twisted mass
#define HAVE_MULT_TM
  //#define HAVE_TM1p1    // flag for enable doublet for twisted mass; unless g.n_flavours==2, Dirac matrix is degenerate, and each part is inverted individually, although the size of the memoery is doubled
  #define INIT_ONE_PREC // flag undef for enabling additional features in the lib

  #define num_loop 4
  #define SIMD_byte 32
  #define MUN_C
  #if num_loop == 1
    #define VECTOR_LOOP(j, jmax, jj, instructions) for( j=0; j<jmax; j++) { jj=0; instructions; }
  #else
    #define VECTOR_LOOP(j, jmax, jj, instructions) for( j=0; j<jmax; j+=num_loop) {_Pragma("unroll") _Pragma("vector aligned") _Pragma("ivdep") for( jj=0; jj<num_loop; jj++) { instructions; }} 
  #endif

  // These explicit repetitions are faster than for-loops.
//  #define DO_UNROLL(EXP) _Pragma (#EXP )
//  #define FORN( N, e ) DO_UNROLL(unroll (N)) for(int i=0; i < N; i++){ e }; // suggestion: This is relevant only when compiled with -O3 -time
  #define FOR2( e )  { e e }
  #define FOR6( e )  { e e e  e e e }
  #define FOR12( e ) { e e e  e e e  e e e  e e e }       // used only in the master thread
  #define FOR36( e ) { FOR12( e ) FOR12( e ) FOR12( e ) } // used only in the master thread
  #define FOR42( e ) { FOR36( e ) FOR6( e ) }             // used only in the master thread
  #define VECTOR_FOR( start, end, expression, update, l ) do{ \
            if ( l->depth == 0 ) {				      \
	      for ( start; end; )				      \
		FOR12( expression; update; )			      \
		  } else {					      \
	      for ( start; end; )				      \
		FOR2( expression; update; )			      \
		  }						      \
          } while(0)

  #define SQUARE( e ) (e)*(e)
  #define NORM_SQUARE_float( e ) SQUARE(crealf( e ))+SQUARE(cimagf( e ))
  #define NORM_SQUARE_double( e ) SQUARE(creal( e ))+SQUARE(cimag( e ))
  #define CSPLIT( e ) creal(e), cimag(e)

  #define MPI_double MPI_DOUBLE
  #define MPI_float MPI_FLOAT
  #define MPI_COMPLEX_double MPI_DOUBLE_COMPLEX
  #define MPI_COMPLEX_float MPI_COMPLEX
  #define FABULOUS_REAL_float FABULOUS_REAL_FLOAT 
  #define FABULOUS_REAL_double FABULOUS_REAL_DOUBLE 
  #define FABULOUS_COMPLEX_float FABULOUS_COMPLEX_FLOAT
  #define FABULOUS_COMPLEX_double FABULOUS_COMPLEX_DOUBLE
  #define I _Complex_I
  #define conj_double conj
  #define conj_float conjf
  #define cabs_double cabs
  #define cabs_float cabsf
  #define creal_double creal
  #define creal_float crealf
  #define cimag_double cimag
  #define cimag_float cimagf
  #define csqrt_double csqrt
  #define csqrt_float csqrtf
  #define sqrt_double sqrt
  #define sqrt_float sqrtf
  #define cpow_double cpow
  #define cpow_float cpowf
  #define pow_double pow
  #define pow_float powf
  #define abs_double fabs
  #define abs_float fabsf
  //printf("%d malloc of %s (%s:%d) kind: %s length: %d\n",g.my_rank,#variable, __FILE__, (int)__LINE__, #kind, (int) length);}
  //printf("%d free of %s (%s:%d) kind: %s length: %d\n",g.my_rank,#variable, __FILE__, (int)__LINE__, #kind, (int) length);

  // I temporaly replaced malloc by _mm_malloc; I might need to use this function only in some selected MALLOC statement!!!!
  #define MALLOC( variable, kind, length ) do{ if ( variable != NULL ) { \
  printf0("malloc of \"%s\" failed: pointer is not NULL (%s:%d).\n", #variable, __FILE__, __LINE__ ); } \
  if ( (length) > 0 ) { variable = (kind*) malloc( sizeof(kind) * (length) ); }	\
  if ( variable == NULL && (length) > 0 ) { \
  error0("malloc of \"%s\" failed: no memory allocated (%s:%d), current memory used: %lf GB.\n", \
  #variable, __FILE__, __LINE__, g.cur_storage/1024.0 ); } \
  g.cur_storage += (sizeof(kind) * (length))/(1024.0*1024.0); \
  if ( g.cur_storage > g.max_storage ) g.max_storage = g.cur_storage; }while(0)

  #define FREE( variable, kind, length ) do{ if ( variable != NULL ) { \
  free( variable ); variable = NULL; g.cur_storage -= (sizeof(kind) * (length))/(1024.0*1024.0); \
  } else {								\
  printf0("multiple free of \"%s\"? pointer is already NULL (%s:%d).\n", #variable, __FILE__, __LINE__ ); } }while(0)

  // if -std=c11 has been chosen, can use aligned_alloc; the following is Intel compiler specific
  #define ALI_MALLOC( variable, kind, length ) do{ if ( variable != NULL ) { \
    printf0("malloc of \"%s\" failed: pointer is not NULL (%s:%d).\n", #variable, __FILE__, __LINE__ ); } \
    if ( (length) > 0 ) { variable = (kind*) _mm_malloc( sizeof(kind) * (length), SIMD_byte ); } \
    if ( variable == NULL && (length) > 0 ) { \
    error0("malloc of \"%s\" failed: no memory allocated (%s:%d), current memory used: %lf GB.\n", \
	   #variable, __FILE__, __LINE__, g.cur_storage/1024.0 ); }	\
    g.cur_storage += (sizeof(kind) * (length))/(1024.0*1024.0); \
    if ( g.cur_storage > g.max_storage ) g.max_storage = g.cur_storage; }while(0)

  #define ALI_FREE( variable, kind, length ) do{ if ( variable != NULL ) { \
    _mm_free( variable ); variable = NULL; g.cur_storage -= (sizeof(kind) * (length))/(1024.0*1024.0); } else { \
    printf0("multiple free of \"%s\"? pointer is already NULL (%s:%d).\n", #variable, __FILE__, __LINE__ ); } }while(0)

  // allocate and deallocate macros (hugepages, aligned)
  #include <fcntl.h>
  #include <sys/mman.h>
  #define HUGE_PAGE_SIZE (2 * 1024 * 1024)
  #define ROUND_UP_TO_FULL_PAGE(x) \
    (((x) + HUGE_PAGE_SIZE - 1) / HUGE_PAGE_SIZE * HUGE_PAGE_SIZE)
  
  #define MALLOC_HUGEPAGES( variable, kind, length, alignment ) do { if ( variable != NULL ) { \
  printf0("malloc of \"%s\" failed: pointer is not NULL (%s:%d).\n", #variable, __FILE__, __LINE__ ); } \
  if ( (length) > 0 ) { \
  variable = (kind*)memalign( alignment, sizeof(kind)*((size_t)length)); } \
  if ( variable == NULL && (length) > 0 ) { \
  error0("malloc of \"%s\" failed: no memory allocated (%s:%d), current memory used: %lf GB.\n", \
  #variable, __FILE__, __LINE__, g.cur_storage/1024.0 ); } \
  g.cur_storage += (sizeof(kind) * (length))/(1024.0*1024.0); \
  if ( g.cur_storage > g.max_storage ) g.max_storage = g.cur_storage; }while(0)
  
  #define FREE_HUGEPAGES( addr, kind, length ) FREE( addr, kind, length )
  
  #ifdef DEBUG
    #define DPRINTF0 printf0
  #else
    #define DPRINTF0( ARGS, ... )
  #endif
    
  // Allocate memory in the master and assign the pt to variable in each thread (variable is assumed to be local to thread and thus private)
  #define PUBLIC_MALLOC( variable, kind, size ) do{ START_MASTER(threading) MALLOC( variable, kind, size ); \
  ((kind**)threading->workspace)[0] = variable; END_MASTER(threading) SYNC_MASTER_TO_ALL(threading) \
  variable = ((kind**)threading->workspace)[0]; SYNC_MASTER_TO_ALL(threading) }while(0)
  
  #define PUBLIC_FREE( variable, kind, size ) do{ SYNC_MASTER_TO_ALL(threading) \
  START_MASTER(threading) FREE( variable, kind, size ); END_MASTER(threading) SYNC_MASTER_TO_ALL(threading) variable = NULL; }while(0)
  
  #define ALI_PUBLIC_MALLOC( variable, kind, size ) do{ START_MASTER(threading) ALI_MALLOC( variable, kind, size ); \
  ((kind**)threading->workspace)[0] = variable; END_MASTER(threading) SYNC_MASTER_TO_ALL(threading) \
  variable = ((kind**)threading->workspace)[0]; SYNC_MASTER_TO_ALL(threading) }while(0)

  #define ALI_PUBLIC_FREE( variable, kind, size ) do{ SYNC_MASTER_TO_ALL(threading) \
  START_MASTER(threading) ALI_FREE( variable, kind, size ); END_MASTER(threading) SYNC_MASTER_TO_ALL(threading) variable = NULL; }while(0)

  #define PUBLIC_MALLOC2( variable, kind, size, thrdng ) do{ START_MASTER(thrdng) MALLOC( variable, kind, size ); \
  ((kind**)thrdng->workspace)[0] = variable; END_MASTER(thrdng) SYNC_MASTER_TO_ALL(thrdng) \
  variable = ((kind**)thrdng->workspace)[0]; SYNC_MASTER_TO_ALL(thrdng) }while(0)
  
  #define PUBLIC_FREE2( variable, kind, size, thrdng ) do{ SYNC_MASTER_TO_ALL(thrdng) \
  START_MASTER(thrdng) FREE( variable, kind, size ); END_MASTER(thrdng) SYNC_MASTER_TO_ALL(thrdng) variable = NULL; }while(0)
  
  
  #define ASSERT( expression ) do{ if ( !(expression) ) { \
  error0("assertion \"%s\" failed (%s:%d)\n       bad choice of input parameters (please read the user manual in /doc).\n", \
  #expression, __FILE__, __LINE__ ); } }while(0)
  
  #define IMPLIES( A, B ) !( A ) || ( B )
  #define XOR( A, B ) (( A ) && !( B )) || (!( A ) && ( B ))
  #define NAND( A, B ) !( (A) && (B) )
  #define DIVIDES( A, B ) A == 0 || ((double)(B)/(double)(A) - (double)((int)(B)/(int)(A))) == 0 
  #define ASCENDING( A, B, C ) ( (A)<=(B) ) && ( (B)<=(C) )
  #define MAX( A, B ) ( (A > B) ? A : B )
  #define MIN( A, B ) ( (A < B) ? A : B )
  
  #ifdef DEBUG
  #define DEBUGOUTPUT_ARRAY( A, FORMAT, INDEX ) do{ \
  char TMPSTRING[100]; sprintf( TMPSTRING, "%s[%d] = %s\n", #A, INDEX, FORMAT ); \
  printf0( TMPSTRING, A[INDEX] ); }while(0)
  #else
  #define DEBUGOUTPUT_ARRAY( A, FORMAT, INDEX )
  #endif
  
  #ifdef DEBUG
  #define DEBUGOUTPUT( A, FORMAT ) do{ \
  char TMPSTRING[100]; sprintf( TMPSTRING, "%s = %s\n", #A, FORMAT ); \
  printf0( TMPSTRING, A ); }while(0)
  #else
  #define DEBUGOUTPUT( A, FORMAT )
  #endif
  
  #include "threading.h"

  // enumerations
  enum { _EVEN, _ODD };
  enum { _NO_DEFAULT_SET, _DEFAULT_SET };
  enum { _NO_REORDERING, _REORDER };
  enum { _ADD, _COPY };
  enum { _ORDINARY, _SCHWARZ, _ODDEVEN, _INNER, _EVEN_INNER };
  enum { _RES, _NO_RES };// _RES: Use residual r as the right-side; _NO_RES: Use b as the right-hand side, which occurs when x=0 as r=b-Ax_0 = b.  
  enum { _STANDARD, _LIME }; //formats
  enum { _READ, _WRITE };
  enum { _NO_SHIFT };
  enum { _BTWN_ORTH = 20 };
  enum { _GLOBAL_FSOLVER, _K_CYCLE, _COARSE_SOLVER, _SMOOTHER };
  enum { _COARSE_GLOBAL };
  enum { _FULL_SYSTEM, _EVEN_SITES, _ODD_SITES };
  enum { _LEFT, _RIGHT, _NOTHING };
  enum { _PERIODIC, _ANTIPERIODIC, _TWISTED, _DIRICHLET };
  enum { _GIP, _PIP, _LA2, _LA6, _LA8, _LA, _CPY, _SET, _PR, _SC, _NC, _SM, _OP_COMM, _OP_IDLE, _ALLR, _GD_COMM, _GD_IDLE, _GRAM_SCHMIDT, _GRAM_SCHMIDT_ON_AGGREGATES,
	 _SM1, _SM2, _SM3, _SM4, _SMALL1, _SMALL2, _RS, _RL, _FIP, _FMVP, _FAB_COPY, _NUM_PROF }; // _NUM_PROF has always to be the last constant!
  enum { _VTS = 20 };
  enum { _TRCKD_VAL, _STP_TIME, _SLV_ITER, _SLV_TIME, _CRS_ITER, _CRS_TIME, _SLV_ERR, _CGNR_ERR, _NUM_OPTB };
  enum { _NVEC_OUTER, _NVEC_INNER }; //vector layout: spin first; vector first
  enum { _FGMRES, _BGMRES, _GCR, _IB, _DR, _IBDR, _QR, _QRIB, _QRDR, _QRIBDR };

  // structures
  typedef struct block_struct {
    int start, color, no_comm, *bt;
  } block_struct;
  
  #include "main_pre_def_float.h"
  #include "main_pre_def_double.h"
  
  extern complex_double _COMPLEX_double_ONE;
  extern complex_double _COMPLEX_double_ZERO;
  extern complex_double _COMPLEX_double_MINUS_ONE;
  extern complex_float  _COMPLEX_float_ONE;
  extern complex_float  _COMPLEX_float_ZERO;
  extern complex_float  _COMPLEX_float_MINUS_ONE;
  
  typedef struct plot_table_line {
    
    double values[_NUM_OPTB];
    struct plot_table_line *next;
    
  } plot_table_line;
  
  typedef struct var_table_entry { 
    
    void *pt;
    char name[STRINGLENGTH];
    char datatype[20];
    struct var_table_entry *next;
    
  } var_table_entry;
  
  typedef struct var_table {// top level struct
    
    int evaluation, multiplicative, shift_update, re_setup,
        track_error, track_cgn_error, average_over;
    char scan_var[STRINGLENGTH];
    double start_val, end_val, step_size, *output_table[6];
    var_table_entry *entry, *iterator; // iterator is temp field to iterate over chain of var_table_entry's; entry pts to the actual var_table_entry w/ data
    plot_table_line *p, *p_end;
    
  } var_table;
  
  typedef struct confbuffer_struct {
    
    double *data;
    struct confbuffer_struct *next;
    
  } confbuffer_struct;
  
  typedef struct {
    
    gmres_float_struct sp;
    gmres_double_struct dp;
    
  } gmres_MP_struct; // MP=Mixed Precision

  typedef struct level_struct {    
    
    // distributed: non-idling processes of previos level
    // gathered: non-idling processes of current level
    
    // distributed
    operator_double_struct op_double;
    operator_float_struct op_float;
    // odd_even
    operator_double_struct oe_op_double;
    operator_float_struct oe_op_float;
    // gathered / schwarz
    schwarz_double_struct s_double;
    schwarz_float_struct s_float;
    // interpolation / aggregation
    interpolation_double_struct is_double;
    interpolation_float_struct is_float;
    // gathering parameters and buffers
    gathering_double_struct gs_double;
    gathering_float_struct gs_float;
    // gmres used in V/K-cycle
    gmres_float_struct p_float;
    gmres_double_struct p_double;
    // gmres as a smoother: used when g.method >= 4
    gmres_float_struct sp_float;
    gmres_double_struct sp_double;
    // dummy gmres struct
    gmres_float_struct dummy_p_float;
    gmres_double_struct dummy_p_double;
    //profiling
    profiling_float_struct prof_float;
    profiling_double_struct prof_double;
    
    // communication
    MPI_Request *reqs;
    int parent_rank;          // rank that takes care of the sites on the idle processes in the family (a group of a parent and idle children)
    int idle, neighbor_rank[8];
    int num_processes;        // total #active processes at the given level
    int num_processes_dir[4]; // #processes in a given dir when mapped to the Cartesian process topology
    int comm_offset[4];       // =(#processes in the given dir on the top level)/(#processes in the given dir on the current level); some processes become idle
    // lattice
    int *global_lattice; // dims of global lattice at the given depth
    int *local_lattice;  // dims of local lattice on each process at the given depth
    int *block_lattice;  // dims of a block within a local lattice
    int num_eig_vect;        // #eigenvectors used for interpolation/restriction on a given level
    int num_parent_eig_vect; // #eigenvectors used for interpolation/restriction on a level one up
    int coarsening[4];       // dims of an aggregate, corresponding to a lattice site on a coarsened lattice.
    int global_splitting[4]; // #processes in the given dir in the Cartesian topology of processes on the global lattice
    int periodic_bc[4];      // specify whether the grid is periodic (true) or not (false) in each dimension
    // degrees of freedom on a site:
    //   fine lattice: 12 (i.e., complex d.o.f. of spin and color)
    //   coarser lattices: 2*num_eig_vect
    int num_lattice_site_var;    // # d.o.f. on a given site
    int level;                   // counts level of coarsening from the bottom, i.e., coarsest (=0) to the top (finest)
    int depth;                   // counts level of coarsening from the top, i.e., finest (=0) to the bottom (coarsest)
    int num_inner_lattice_sites; // number of sites in local volume
    int num_lattice_sites;       // number of sites in local volume + ghost shell (either fw or bw: one-sided boundry sites)
    int num_boundary_sites[4];
    long int inner_vector_size;  // complex d.o.f. in local volume = num_inner_lattice_sites * num_lattice_site_var
    long int vector_size;        // complex d.o.f. in local volume + ghost shell (either fw or bw) = num_lattice_sites * num_lattice_site_var
    long int schwarz_vector_size;// 2*vector_size - inner_vector_size = #local lattice sites + full ghost shell
    int D_size;
    int clover_size;
    int block_size;
    // buffer vectors: vbuf has num_eig_vect many vectors
    vector_float vbuf_float[5], sbuf_float[2];
    vector_double vbuf_double[5], sbuf_double[2];
    // storage + daggered-operator bufferes
    vector_double x; // used in io.c
    // local solver parameters
    double tol, relax_fac;
    int n_cy, post_smooth_iter, block_iter, setup_iter;
    
    // next coarser level
    struct level_struct *next_level;
    
  } level_struct;


  typedef struct global_struct {
    
    FILE *logfile;
    
    gmres_double_struct p;
    gmres_MP_struct p_MP;
    operator_double_struct op_double;
    operator_float_struct op_float;

    // parameters related the "fabulous" library
    int *solver;
    fabulous_orthoscheme *f_orthoscheme;
    fabulous_orthotype *f_orthotype;
    int *ortho_iter;           // #iteration for Iterated Schemas (IMGS and ICGS); Must be positive integer >= 2
    int *max_kept_direction;   // max #kept direction per iteration
    int *real_residual;        // if 1, compute X and R at each iteration such that the user can access them in fabulous CallBack funciton
    int logger_user_data_size; // #slots for user data when calling fabulous_set_iteration_user_data()
    int quiet;                 // if 1, no output to stdout when running
    int *k;                    // #deflating eigevecs at each level
    int *max_mvp;              
    int use_only_fgrmes_at_setup;
    
    // communication
    MPI_Comm comm_cart;
    MPI_Group global_comm_group;
    MPI_Request sreqs[8], rreqs[8];
    int num_processes, my_rank, my_coords[4], tv_io_single_file, num_openmp_processes;
    // string buffers
    char in[STRINGLENGTH], in_clov[STRINGLENGTH], source_list[STRINGLENGTH], tv_io_file_name[STRINGLENGTH];
    
    // geometry, method parameters
    //   process_grid: #processes in the mu dir in the Cartesian topology of processes
    int num_levels, num_desired_levels, process_grid[4], in_format,
        **global_lattice, **local_lattice, **block_lattice, 
        *post_smooth_iter, *block_iter, *setup_iter, *ncycle,
        method, odd_even, rhs, propagator_coords[4],
        interpolation, randomize, *num_eig_vect, num_coarse_eig_vect, kcycle, mixed_precision,
        restart, max_restart, kcycle_restart, kcycle_max_restart, coarse_iter, coarse_restart;
    double tol, coarse_tol, kcycle_tol, csw, rho, *relax_fac;

    // profiling, analysis, output
    int coarse_iter_count, *iter_counts, iter_count, iterator, print, conf_flag, setup_flag, in_setup;
    double coarse_time, *iter_times, prec_time, resids[num_loop], *output_table[8], cur_storage, max_storage, total_time,
      plaq_hopp, plaq_clov, norm_res, plaq, bicgstab_tol, twisted_bc[4], test;

    double m0, setup_m0;

#ifdef HAVE_TM
    // twisted mass parameters
    //   setup_mu should be smaller than or equal to mu
    //   As long as this is true, the performance of the solver part does not depend on setup_mu
    //   Otherwise, the solver becomes slower.
    int downprop, is_even_shifted_mu_nonzero, no_shift, even_shifted, n_chunk;
    double mu, setup_mu, mu_odd_shift, *mu_even_shift, *mu_factor;
#endif

#ifdef HAVE_TM1p1           
    int n_flavours;
    double epsbar, epsbar_ig5_odd_shift, epsbar_ig5_even_shift, *epsbar_factor;
#endif

    // index functions for external usage
    int (*conf_index_fct)(), (*vector_index_fct)();
    int *odd_even_table;
    int (*Cart_rank)(MPI_Comm comm, const int coords[], int *rank);
    int (*Cart_coords)(MPI_Comm comm, int rank, int maxdims, int coords[]);
    
    // bc: 0 dirichlet, 1 periodic, 2 anti-periodic
    int bc; 
    
    // number of rhs vectors (b) to be solved at the same time (hopefully)
    int num_rhs_vect, num_vect_now, num_vect_pass1, num_vect_pass2;// num_vect_now unnecessary!!!
    

    complex_double **gamma;
    var_table vt;
    
  } global_struct;

  extern global_struct g;
  
  // inline functions
  static inline void printf0( char* format, ... ) {
    START_MASTER(no_threading)
    if ( g.my_rank == 0 && g.print >= 0 ) {
      va_list argpt;
      va_start(argpt,format);
      vprintf(format,argpt);
#ifdef WRITE_LOGFILE
      vfprintf(g.logfile,format,argpt);
      fflush(g.logfile);
#endif
      va_end(argpt);
      fflush(0);
      //fflush(1)//my suggestion
    }
    END_MASTER(no_threading)
  }

  static inline void printf00( char* format, ... ) {
    if ( g.my_rank == 0 && g.print >= 0 ) {
      va_list argpt;
      va_start(argpt,format);
      vprintf(format,argpt);
#ifdef WRITE_LOGFILE
      vfprintf(g.logfile,format,argpt);
      fflush(g.logfile);
#endif
      va_end(argpt);
      fflush(0);
    }
  }

static inline void printfv_float ( vector_float *v ) {
  int vjj, vj;
  for ( int vi=0; vi<v->size; vi++ )
    VECTOR_LOOP(vj, v->num_vect_now, vjj, printf("v_%d[%d,%d]=%g %g ",vj+vjj,vi/v->l->num_lattice_site_var,vi%v->l->num_lattice_site_var,creal_float((float)v->vector_buffer[vi*v->num_vect+vj+vjj]),cimag_float((float)v->vector_buffer[vi*v->num_vect+vj+vjj]));)
}
static inline void printfv_double ( vector_double *v ) {
  int vjj, vj;
  for ( int vi=0; vi<v->size; vi++ )
    VECTOR_LOOP(vj, v->num_vect_now, vjj, printf("v_%d[%d]=%g %g ",vj+vjj,vi,creal_float((double)v->vector_buffer[vi*v->num_vect+vj+vjj]),cimag_double((float)v->vector_buffer[vi*v->num_vect+vj+vjj]));)
      }
  static inline void warning( char* format, ... ) {
    printf("\x1b[31mwarning, rank %d: ", g.my_rank);
    va_list argpt;
    va_start(argpt,format);
    vprintf(format,argpt);
#ifdef WRITE_LOGFILE
    vfprintf(g.logfile,format,argpt);
    fflush(g.logfile);
#endif
    va_end(argpt);
    printf("\x1b[0m");
    fflush(0);
  }

  static inline void warning0( char* format, ... ) {
    if ( g.my_rank == 0 && g.print >= 0 ) {
      printf("\x1b[31mwarning: ");
      va_list argpt;
      va_start(argpt,format);
      vprintf(format,argpt);
#ifdef WRITE_LOGFILE
      vfprintf(g.logfile,format,argpt);
      fflush(g.logfile);
#endif
      va_end(argpt);
      printf("\x1b[0m");
      fflush(0);
    }
  }
  
  static inline void error0( char* format, ... ) {
    if ( g.my_rank == 0 ) {
      printf("\x1b[31merror: ");
      va_list argpt;
      va_start(argpt,format);
      vprintf(format,argpt);
#ifdef WRITE_LOGFILE
      vfprintf(g.logfile,format,argpt);
      fflush(g.logfile);
#endif
      va_end(argpt);
      printf("\x1b[0m");
      fflush(0);
      MPI_Abort( MPI_COMM_WORLD, 0 );
    }
  }
  
#endif

// includes
#include "clifford.h"

#include "interpolation_float.h"
#include "interpolation_double.h"

#include "data_float.h"
#include "data_double.h"
#include "data_layout.h"
#include "io.h"
#include "init.h"
#include "readin.h"
#include "operator_float.h"
#include "operator_double.h"
#include "dirac.h"
#include "dirac_float.h"
#include "dirac_double.h"
#include "oddeven_float.h"
#include "oddeven_double.h"
#include "block_oddeven_float.h"
#include "block_oddeven_double.h"
#include "linalg.h"
#include "linalg_float.h"
#include "linalg_double.h"
#include "ghost_float.h"
#include "ghost_double.h"
#include "linsolve_float.h"
#include "linsolve_double.h"
#include "linsolve.h"
#include "preconditioner.h"
#include "vcycle_float.h"
#include "vcycle_double.h"
#include "solver_analysis.h"
#include "top_level.h"
#include "ghost.h"
#include "init_float.h"
#include "init_double.h"
#include "schwarz_double.h"
#include "schwarz_float.h"
#include "setup_float.h"
#include "setup_double.h"
#include "coarsening_float.h"
#include "coarsening_double.h"
#include "gathering_float.h"
#include "gathering_double.h"
#include "coarse_operator_float.h"
#include "coarse_operator_double.h"
#include "coarse_oddeven_float.h"
#include "coarse_oddeven_double.h"
#include "var_table.h"
#include "main_post_def_float.h"
#include "main_post_def_double.h"
#include "vector_float.h"
#include "vector_double.h"
#include "fabulous_float.h"
#include "fabulous_double.h"
#ifdef HAVE_LIME
#include <lime.h>
#include <lime_config.h>
#include <dcap-overload.h>
#include <lime_defs.h>
#include <lime_header.h>
#include <lime_writer.h>
#include <lime_reader.h>
#endif
#include "lime_io.h"
