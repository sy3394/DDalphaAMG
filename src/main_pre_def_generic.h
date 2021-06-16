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
 * checked: 11/29/2019
 * changed from sbacchio
 */

#ifndef MAIN_PRE_DEF_PRECISION_HEADER
  #define MAIN_PRE_DEF_PRECISION_HEADER
  

  typedef PRECISION complex complex_PRECISION;
  typedef PRECISION complex *config_PRECISION;
  typedef PRECISION complex *buffer_PRECISION;

  typedef struct {
    buffer_PRECISION vector_buffer;
    int num_vect;
    int num_vect_now;
    int layout;
    int type;
    int size;
    struct level_struct *l;
  } vector_PRECISION;

  typedef struct {
    int length[8], *boundary_table[8], max_length[4],
      comm_start[8], in_use[8], offset, comm, num_vect,
        num_even_boundary_sites[8], num_odd_boundary_sites[8],
        num_boundary_sites[8];
    buffer_PRECISION buffer[8];
    MPI_Request sreqs[8], rreqs[8];
  } comm_PRECISION_struct;
  
  typedef struct {

    int dist_local_lattice[4];    // dims of a local coarsened lattice?????
    int dist_inner_lattice_sites; // #local lattice inner sites of a coarsened lattice
    int  *permutation, *gather_list, gather_list_length;
    vector_PRECISION buffer, transfer_buffer;// buffer is used in gather routines
    MPI_Request *reqs;
    MPI_Group level_comm_group;
    MPI_Comm level_comm;
  } gathering_PRECISION_struct;
  
  typedef struct {//operators are row-major????
    double m0;
    config_PRECISION D, clover, clover_oo_inv;
    config_PRECISION odd_proj; //identity on the odd sites
    int oe_offset, self_coupling, num_even_sites, num_odd_sites,
      *index_table, pr_num_vect,
        *neighbor_table, *translation_table, table_dim[4],
        *backward_neighbor_table,
        table_mod_dim[4], *config_boundary_table[4];
    vector_PRECISION *buffer; // used in computing schur_complement
    buffer_PRECISION prnT, prnZ, prnY, prnX, prpT, prpZ, prpY, prpX;
    comm_PRECISION_struct c;
#ifdef HAVE_TM
    int is_even_shifted_mu_nonzero;
    double factor, mu, mu_odd_shift, *mu_even_shift, odd_shifted_mu, *diff_mu_eo, even_shift_avg;
    config_PRECISION tm_term;
#endif
#ifdef HAVE_TM1p1
    double epsbar, epsbar_ig5_odd_shift, epsbar_ig5_even_shift;
    config_PRECISION epsbar_term, clover_doublet_oo_inv;
#endif
  } operator_PRECISION_struct;

#ifdef HAVE_FABULOUS
  typedef struct {
    // max k-value = (max def space size) = (Max Krylov space size) x 2 - nrhs (<- Max Krylov space size = nrhs*iter_before_restart)
    int nrhs, dim, ldb, ldx, ldu, k, dsize, max_iter, mvp;
    vector_PRECISION B, X, B0, X0, C0;
    void *eigvals, *U;  
    fabulous_handle handle;
    struct level_struct *l;
    struct Thread *threading;
  } fabulous_PRECISION_struct;
#endif

  typedef struct {
    vector_PRECISION x, b, r, w, *V, *Z;
    complex_PRECISION **H, *y, *gamma, *c, *s;
    config_PRECISION *D, *clover;
    operator_PRECISION_struct *op;
#ifdef HAVE_FABULOUS
    fabulous_PRECISION_struct fab;
#endif
    PRECISION tol;
    int num_restart, restart_length, timing, print, kind,
      initial_guess_zero, layout, v_start, v_end, num_vect;
    long int total_storage;
    void (*preconditioner)();
    void (*eval_operator)();
  } gmres_PRECISION_struct;
  
  typedef struct {
    operator_PRECISION_struct op;
    vector_PRECISION buf[5];
    vector_PRECISION oe_buf[4];
    vector_PRECISION local_minres_buffer[3];
    int block_oe_offset, *index[4],  dir_length[4], num_blocks, //# blocks within a local lattice
      num_colors, // # colors used to color/partition the local lattice in SAP
        dir_length_even[4], dir_length_odd[4], *oe_index[4],
        num_block_even_sites, num_block_odd_sites, num_aggregates,
      block_vector_size, num_block_sites, // # sites in each block in a local lattice
      block_boundary_length[9],
        **block_list, *block_list_length;
    block_struct *block;
  } schwarz_PRECISION_struct;

  typedef struct {
    int num_agg, *agg_index[4], agg_length[4], *agg_boundary_index[4],
        *agg_boundary_neighbor[4], agg_boundary_length[4], num_bootstrap_vect;
    vector_PRECISION test_vector_vec, interpolation_vec;//*test_vector, *interpolatio *bootstrap_vector, tmp,
    complex_PRECISION *operator, *eigenvalues; //, *bootstrap_eigenvalues;
  } interpolation_PRECISION_struct;

  typedef struct {
    double time[_NUM_PROF];
    double flop[_NUM_PROF];
    double count[_NUM_PROF];
    char name[_NUM_PROF][50];
  } profiling_PRECISION_struct;
  
  #ifdef PROFILING
    #define PROF_PRECISION_START_UNTHREADED( TYPE ) do{ l->prof_PRECISION.time[TYPE] -= MPI_Wtime(); }while(0)
    #define PROF_PRECISION_START_THREADED( TYPE, threading ) do{ if(threading->core + threading->thread == 0) l->prof_PRECISION.time[TYPE] -= MPI_Wtime(); }while(0)
  #else
    #define PROF_PRECISION_START_UNTHREADED( TYPE )
    #define PROF_PRECISION_START_THREADED( TYPE, threading )
  #endif
  
  #ifdef PROFILING
    #define PROF_PRECISION_STOP_UNTHREADED( TYPE, COUNT ) do{ l->prof_PRECISION.time[TYPE] += MPI_Wtime(); \
    l->prof_PRECISION.count[TYPE] += COUNT; }while(0)
    #define PROF_PRECISION_STOP_THREADED( TYPE, COUNT, threading ) do{ if(threading->core + threading->thread == 0) { l->prof_PRECISION.time[TYPE] += MPI_Wtime(); \
    l->prof_PRECISION.count[TYPE] += COUNT; } }while(0)
  #else
    #define PROF_PRECISION_STOP_UNTHREADED( TYPE, COUNT )
    #define PROF_PRECISION_STOP_THREADED( TYPE, COUNT, threading )
  #endif

  #define GET_MACRO2(_1,_2,NAME,...) NAME
  #define GET_MACRO3(_1,_2,_3,NAME,...) NAME
  #define PROF_PRECISION_START(...) GET_MACRO2(__VA_ARGS__, PROF_PRECISION_START_THREADED, PROF_PRECISION_START_UNTHREADED, padding)(__VA_ARGS__)
  #define PROF_PRECISION_STOP(...) GET_MACRO3(__VA_ARGS__, PROF_PRECISION_STOP_THREADED, PROF_PRECISION_STOP_UNTHREADED, padding)(__VA_ARGS__)
  
#endif
