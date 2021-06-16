/*
 * Copyright (C) 2016, Simone Bacchio.
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
 * changed from sbacchio
 * NOTE: decide how lime_read_vector is used; that's not updated for multiple vectors
 */
 
#include "main.h"
#include "DDalphaAMG.h"

global_struct g;
static level_struct l;
static int (*conf_index_fct)(int t, int z, int y, int x, int mu); 
static int (*vector_index_fct)(int t, int z, int y, int x); 
struct common_thread_data *commonthreaddata;
struct Thread **threading;
struct Thread *no_threading;

/********************  Initialization and Setup  ****************************************************************/
/*
 * Mirror of 
 *   method_init( &argc, &argv, &l );
 *   setup_threading( ... );
 */
void DDalphaAMG_initialize( DDalphaAMG_init *mg_init, DDalphaAMG_parameters *mg_params, DDalphaAMG_status *mg_status ) {
  
  double t0, t1;
  t0 = MPI_Wtime();

  /*
   * BEGIN: method_init( &argc, &argv, &l );
   */
  MPI_Comm_rank( mg_init->comm_cart, &(g.my_rank) );

  g_init();
  l_init( &l );

  set_DDalphaAMG_parameters( mg_init, &l );
    
  conf_index_fct = NULL;
  vector_index_fct = NULL;

  int topo_type;
  MPI_Topo_test(mg_init->comm_cart, &topo_type);

  if(mg_init->Cart_coords==NULL)
    g.Cart_coords = MPI_Cart_coords;
  else
    g.Cart_coords = mg_init->Cart_coords;

  if(mg_init->Cart_rank==NULL)
    g.Cart_rank = MPI_Cart_rank;
  else
    g.Cart_rank = mg_init->Cart_rank;

  if( (g.Cart_rank==MPI_Cart_rank ||
       g.Cart_coords==MPI_Cart_coords) && topo_type!=MPI_CART ){
    warning0("Defining a new communicator\n");
    cart_define( mg_init->comm_cart, &l );
  } else
    cart_validate( mg_init->comm_cart, &l );
  
  if (mg_init->rnd_seeds == NULL)
    srand( time( 0 ) + 1000*g.my_rank );
  else
    srand( mg_init->rnd_seeds[g.my_rank] );

  data_layout_init( &l );
  operator_double_alloc( &(g.op_double), _ORDINARY, &l );
  operator_double_define( &(g.op_double), &l );
  MALLOC( g.odd_even_table, int, l.num_inner_lattice_sites );
  define_odd_even_table( &l );

  for ( int i=0; i<g.num_levels; i++ ) g.iter_counts[i] = 0;
  for ( int i=0; i<g.num_levels; i++ ) g.iter_times[i]  = 0;
  
  /*
   * BEGIN: setup_threading( ... );
   */
  no_threading = (struct Thread *)malloc(sizeof(struct Thread));
  setup_no_threading(no_threading, &l);

  commonthreaddata = NULL;
  MALLOC( commonthreaddata, struct common_thread_data, 1);
  init_common_thread_data(commonthreaddata);

  threading = NULL;
  MALLOC( threading, struct Thread *, g.num_openmp_processes);
  for(int i=0; i<g.num_openmp_processes; i++) {
    threading[i] = NULL;
    MALLOC( threading[i], struct Thread, 1);
  }
  THREADED(g.num_openmp_processes)
    setup_threading(threading[omp_get_thread_num()], commonthreaddata, &l);

  g.conf_flag = 0;
  g.setup_flag = 0;

  DDalphaAMG_get_parameters( mg_params );

  t1 = MPI_Wtime();

  mg_status->success = g.num_levels;
  mg_status->time = t1-t0;
  
}

void DDalphaAMG_get_parameters( DDalphaAMG_parameters *mg_params ){

  int i, j;
   
  mg_params->method = g.method;
  mg_params->interpolation = g.interpolation;
  mg_params->mixed_precision = g.mixed_precision;
  mg_params->kcycle_tolerance = g.kcycle_tol;
  mg_params->coarse_tolerance = g.tol[g.num_levels-1];
  mg_params->smoother_iterations = g.post_smooth_iter[0];
  mg_params->conf_index_fct = conf_index_fct;
  mg_params->vector_index_fct = vector_index_fct;
  mg_params->kappa = 0.5/(g.m0 + 4.);
#ifdef HAVE_TM
  mg_params->mu = g.mu;
  mg_params->mu_odd_shift = g.mu_odd_shift;
#else
  mg_params->mu = 0;
  mg_params->mu_odd_shift = 0;
#endif
#ifdef HAVE_TM1p1
  mg_params->epsbar = g.epsbar;
  mg_params->epsbar_ig5_odd_shift = g.epsbar_ig5_odd_shift;
  mg_params->epsbar_ig5_even_shift = g.epsbar_ig5_even_shift;
#else
  mg_params->epsbar = 0;
  mg_params->epsbar_ig5_odd_shift = 0;
  mg_params->epsbar_ig5_even_shift = 0;
#endif
  mg_params->print = g.print;
  
  for( i=0; i<g.num_levels; i++ ) {
    for( j=0; j<4; j++ )
      mg_params->block_lattice[i][j] = g.block_lattice[i][j];
    if( i<g.num_levels-1 )
      mg_params->mg_basis_vectors[i] = g.num_eig_vect[i];
    mg_params->setup_iterations[i] = g.setup_iter[i];
    mg_params->solver[i] = g.solver[i];
#ifdef HAVE_FABULOUS
    mg_params->fab_orthoscheme[i] = g.f_orthoscheme[i];
    mg_params->fab_orthotype[i] = g.f_orthotype[i];
    mg_params->fab_ortho_iter[i] = g.ortho_iter[i];
    mg_params->fab_max_kept_direction[i] = g.max_kept_direction[i];
    mg_params->fab_num_deflating_eig[i] = g.k[i];
#else
    mg_params->fab_orthoscheme[i] = 0;
    mg_params->fab_orthotype[i] = 0;
    mg_params->fab_ortho_iter[i] = 2;
    mg_params->fab_max_kept_direction[i] = -1;
    mg_params->fab_num_deflating_eig[i] = 0;
#endif
#ifdef HAVE_TM
    mg_params->mu_even_shift[i] = g.mu_even_shift[i];
    mg_params->mu_factor[i] = g.mu_factor[i];
#else
    mg_params->mu_even_shift[i] = 0;
    mg_params->mu_factor[i] = 1;
#endif
#ifdef HAVE_TM1p1
    mg_params->epsbar_factor[i] = g.epsbar_factor[i];
#else
    mg_params->epsbar_factor[i] = 1;
#endif
  }  
}

void DDalphaAMG_update_parameters( DDalphaAMG_parameters *mg_params, DDalphaAMG_status * mg_status ) {
  
  int i, j, re_setup=0, re_projs=0, re_dirac=0;
  level_struct *l_tmp;
  double t0, t1, m0;
  
  t0 = MPI_Wtime();
  
  for ( i=0; i<g.num_levels; i++ ) g.iter_counts[i] = 0;
  for ( i=0; i<g.num_levels; i++ ) g.iter_times[i]  = 0;
  
  // int method;
  if ( mg_params->method != g.method ) {
    g.method = mg_params->method;
    if( g.setup_flag ) {
      //TODO: test which cases work and what to do for making the other working
      warning0("Change of method parameter after setup not guaranteed to work\n");
    }
  } 

  // int interpolation;
  if ( g.interpolation != mg_params->interpolation ) {
    g.interpolation = mg_params->interpolation;
    if( g.setup_flag ) {
      //TODO: test which cases work and what to do for making the other working
      warning0("Change of interpolation parameter after setup not guaranteed to work\n");
    }
  } 
  
  // int mixed_precision;
  if ( mg_params->mixed_precision != g.mixed_precision ) {
    g.mixed_precision = mg_params->mixed_precision;
#ifndef INIT_ONE_PREC //change between 1 and 2 allowed
    if( g.setup_flag && mg_params->mixed_precision * g.mixed_precision == 0 ) {
      warning0("Change from mixed_precision==0 to !=0 (or viceversa) needs a new setup.\n");
      re_setup++;
    }
#else
    warning0("Change of mixed_precision needs a new setup.\n");
    re_setup++;
#endif
  }

  // int solver[MAX_MG_LEVELS]
  for ( i=0; i<g.num_levels; i++ )
    if ( mg_params->solver[i] != g.solver[i] ) {
      g.solver[i] = mg_params->solver[i];
      if( g.setup_flag ) {
	warning0("Changing the solver needs a new setup.\n");
	re_setup++;
      }
    }
      
  // int block_lattice[MAX_MG_LEVELS][4];
  for ( i=0; i<g.num_levels; i++ )
    for ( j=0; j<4; j++ )
      if (g.block_lattice[i][j] != mg_params->block_lattice[i][j]) {
        g.block_lattice[i][j] = mg_params->block_lattice[i][j];
        parameter_update(&l);
        if (g.setup_flag) {
          warning0("Change of block_lattice needs a new setup.\n");
          re_setup++;
        }
      }
  
  // int mg_basis_vectors[MAX_MG_LEVELS-1];
  l_tmp=&l;
  for ( i=0; i<g.num_levels-1; i++ ){
    if ( mg_params->mg_basis_vectors[i] != g.num_eig_vect[i] ) {
      g.num_eig_vect[i] = mg_params->mg_basis_vectors[i];
      if( i==0 )
        parameter_update(&l);
      if( g.setup_flag ) {
        if ( mg_params->mg_basis_vectors[i] < g.num_eig_vect[i] )
          re_projs++; //TODO: check if this works
        else { //TODO just compute the extra vectors
          warning0("Increasing mg_basis_vectors needs a new setup.\n");
          re_setup++;
        }
      }
    }
    if( g.setup_flag )
      l_tmp = l_tmp->next_level;
  }
  
  // int setup_iterations[MAX_MG_LEVELS];
  l_tmp=&l;
  for ( i=0; i<g.num_levels; i++ ){
    if ( mg_params->setup_iterations[i] != g.setup_iter[i] ) {
      g.setup_iter[i] = mg_params->setup_iterations[i];
      if( (g.setup_flag && i>0) || (!g.setup_flag && i==0) ) 
        //after setup, l.setup_iter[i] is used as a counter for total number of setup iters
        l_tmp->setup_iter = mg_params->setup_iterations[i];
    }
    if( g.setup_flag )
      l_tmp = l_tmp->next_level;
  }

  // int smoother_iterations;
  l_tmp=&l;
  for ( i=0; i<g.num_levels; i++ ){
    if ( mg_params->smoother_iterations != g.post_smooth_iter[i] ) {
      l_tmp->post_smooth_iter = mg_params->smoother_iterations;
      g.post_smooth_iter[i] = mg_params->smoother_iterations;
    }
    if( g.setup_flag )
      l_tmp = l_tmp->next_level;
  }

  l_tmp=&l;
  for ( i=0; i<g.num_levels; i++ ){
    if (l_tmp->level > 0) {
      // double kcycle_tolerance;
      if ( mg_params->kcycle_tolerance != g.kcycle_tol ) {
        g.kcycle_tol = mg_params->kcycle_tolerance;
        if( g.setup_flag || i==0 ) {  
          if ( g.mixed_precision )
            l_tmp->p_float.tol = g.kcycle_tol;
          else
            l_tmp->p_float.tol = g.kcycle_tol;
        }
      }
    } else {
      // double coarse_tolerance;
      if ( mg_params->coarse_tolerance != g.tol[g.num_levels-1] ){
        g.tol[g.num_levels-1] = mg_params->coarse_tolerance;
        if (g.setup_flag && g.mixed_precision )
          l_tmp->p_float.tol = g.tol[g.num_levels-1];
        else if(g.setup_flag)
          l_tmp->p_float.tol = g.tol[g.num_levels-1];
      }
    }
    
    if( g.setup_flag )
      l_tmp = l_tmp->next_level;
    else
      break;
  }
  
  // double kappa;
  m0 = 1./(2.*mg_params->kappa)-4.; 
  if( g.m0 != m0 ){
    g.m0 = m0;
    THREADED(threading[0]->n_core)
      if ( g.setup_flag )
        m0_update( g.m0, &l, threading[omp_get_thread_num()] );
      else if ( g.conf_flag )
        m0_update_double( g.m0, &(g.op_double), &l, threading[omp_get_thread_num()] );
    re_dirac++;
  }
  
  // double mu;
  // double mu_odd_shift;
  // double* mu_even_shift;
  // double mu_factor[MAX_MG_LEVELS];
#ifdef HAVE_TM
  int update_mu = 0, update_mu_odd = 0, update_mu_even = 0;
  for ( i=0; i<g.num_rhs_vect; i++ ) {
    update_mu_odd  += mg_params->mu_odd_shift != g.mu_odd_shift?1:0;
    for ( i=0; i<g.num_rhs_vect; i++ ) update_mu_even += mg_params->mu_even_shift[i] != g.mu_even_shift[i]?1:0;
  }
  for ( i=0; i<g.num_levels; i++ )
    if (mg_params->mu_factor[i] != g.mu_factor[i] ) {
      g.mu_factor[i] = mg_params->mu_factor[i];
      update_mu = 1;
    }

  if( update_mu || mg_params->mu != g.mu || update_mu_odd || update_mu_even ) {
    g.setup_mu = mg_params->mu;
    g.mu = mg_params->mu;
    for ( i=0; i<g.num_rhs_vect; i++ ) g.mu_even_shift[i] = mg_params->mu_even_shift[i];
    g.mu_odd_shift = mg_params->mu_odd_shift;
    THREADED(threading[0]->n_core)
      if ( g.setup_flag )
        tm_term_update( g.mu, &l, threading[omp_get_thread_num()] );
      else if ( g.conf_flag )
        tm_term_double_setup( g.mu, g.mu_even_shift, g.mu_odd_shift, g.mu_factor[0], &(g.op_double), &l, threading[omp_get_thread_num()] ); 
    re_dirac++;
  }
  
#else
  if ( mg_params->mu != 0 || mg_params->mu_odd_shift != 0 || mg_params->mu_even_shift != 0 )
    warning0("Parameters mu, mu_odd_shift, mu_even_shift not supported without defining HAVE_TM flag.");
#endif

  // int fab_orthoscheme[MAX_MG_LEVELS;
  // int fab_orthotype[MAX_MG_LEVELS;
  // int fab_ortho_iter[MAX_MG_LEVELS;
  // int fab_max_kept_direction[MAX_MG_LEVELS;
  // int fab_num_deflating_eig[MAX_MG_LEVELS;
#ifdef HAVE_FABULOUS
  int update_orth = 0, update_advanced = 0;
  for ( i=0; i<g.num_levels; i++ ) {
    if ( mg_params->fab_orthoscheme[i] != g.f_orthoscheme[i] ) {
      update_orth += 1;
      g.f_orthoscheme[i] = mg_params->fab_orthoscheme[i];
    }
    if ( mg_params->fab_orthotype[i] != g.f_orthotype[i] ) {
      update_orth += 1;
      g.f_orthotype[i] = mg_params->fab_orthotype[i];
    }
    if ( mg_params->fab_ortho_iter[i] != g.ortho_iter[i] ) {
      update_orth += 1;
      g.ortho_iter[i] = mg_params->fab_ortho_iter[i];
    }
    if ( mg_params->fab_max_kept_direction[i] != g.max_kept_direction[i] ){
      update_advanced += 1;
      g.max_kept_direction[i] = mg_params->fab_max_kept_direction[i];
    }
    if ( mg_params->fab_num_deflating_eig[i] != g.k[i] ) {
      update_advanced += 1;
      g.k[i] = mg_params->fab_num_deflating_eig[i];
    }
  }
  // g.mixed_precision could have been changed so that re-setup is not done here
  if ( g.setup_flag && (update_orth || update_advanced ) ) {
    warning0("Changing FABULOUS solver parameters needs a new setup.\n");
    re_setup ++;
  }
#else
  int update_fab = 0;
  for ( i=0; i<g.num_levels; i++ ) {
    if ( mg_params->fab_orthoscheme[i] != 0 ||
	 mg_params->fab_orthotype[i] != 0 ||
	 mg_params->fab_ortho_iter[i] != 2 ||
	 mg_params->fab_max_kept_direction[i] != -1 ||
	 mg_params->fab_num_deflating_eig[i] != 0 )
      update_fab ++;
  }
  if ( update_fab )
    warning0("FABULOUS solver parameters are not supported without defining HAVE_FABULOUS flag.");
#endif
  
#ifdef HAVE_TM1p1  
  int update_eps = 0;

  for ( i=0; i<g.num_levels; i++ )
    if (mg_params->epsbar_factor[i] != g.epsbar_factor[i] ) {
      g.epsbar_factor[i] = mg_params->epsbar_factor[i];
      update_eps = 1;
     }

  if( update_eps || mg_params->epsbar != g.epsbar || mg_params->epsbar_ig5_odd_shift != g.epsbar_ig5_odd_shift || mg_params->epsbar_ig5_even_shift != g.epsbar_ig5_even_shift ){
    g.epsbar = mg_params->epsbar;
    g.epsbar_ig5_even_shift = mg_params->epsbar_ig5_even_shift;
    g.epsbar_ig5_odd_shift = mg_params->epsbar_ig5_odd_shift;
    THREADED(threading[0]->n_core)
      if ( g.setup_flag )
        epsbar_term_update( &l, threading[omp_get_thread_num()] );
      else if ( g.conf_flag )
        epsbar_term_double_setup( g.epsbar, g.epsbar_ig5_even_shift, g.epsbar_ig5_odd_shift, &(g.op_double), &l, threading[omp_get_thread_num()] ); 
    re_dirac++;
  }
#else
  if ( mg_params->epsbar != 0 || mg_params->epsbar_ig5_odd_shift != 0 || mg_params->epsbar_ig5_even_shift != 0 )
    warning0("Parameters epsbar, epsbar_odd_shift, epsbar_even_shift not supported without defining HAVE_TM1p1 flag.");
#endif

  // int (*conf_index_fct)(int t, int z, int y, int x, int mu);
  // int (*vector_index_fct)(int t, int z, int y, int x );
  conf_index_fct = mg_params->conf_index_fct;
  vector_index_fct = mg_params->vector_index_fct;
  
  // int print;
  g.print = mg_params->print;
  
  // UPDATING
  if ( re_setup && g.setup_flag ){ // destroy and repeate setup
    DDalphaAMG_setup( mg_status ); // TODO handle status

  } else if ( re_projs && g.setup_flag ) { //project again the operators
    if ( g.mixed_precision )
      THREADED(threading[0]->n_core)
        re_setup_float( &l, threading[omp_get_thread_num()] ); 
    else
      THREADED(threading[0]->n_core)
        re_setup_double( &l, threading[omp_get_thread_num()] );
    
  } else if ( re_dirac && g.setup_flag ) { //update just the oddeven and vecorized operators
    THREADED(threading[0]->n_core)
      finalize_operator_update( &l, threading[omp_get_thread_num()]);
  }
  
  
  DDalphaAMG_get_parameters( mg_params );
  
  t1 = MPI_Wtime();
  
  mg_status->success = 1+re_setup;// 1: OK, 2>=: re_setup needs to be done
  mg_status->time = t1-t0;
  //why here; no setup is actually done here???????
  for ( i=0; i<g.num_levels; i++ ) mg_status->iter_times[i] = g.iter_times[i];
  for ( i=0; i<g.num_levels; i++ ) mg_status->iter_counts[i] = g.iter_counts[i];
}

void DDalphaAMG_change_mu_sign( DDalphaAMG_status *mg_status ) {
  
  double t0, t1;
  t0 = MPI_Wtime();
  for ( int i=0; i<g.num_levels; i++ ) g.iter_counts[i] = 0;
  for ( int i=0; i<g.num_levels; i++ ) g.iter_times[i]  = 0;
  mg_status->success = 0;
  mg_status->info = 0;  
  
#ifdef HAVE_TM
  g.mu *= -1;
  for ( int i=0; i<g.num_rhs_vect; i++ ) g.mu_even_shift[i] *= -1;
  g.mu_odd_shift *= -1;

  if (g.conf_flag && !g.setup_flag ) {
    
    THREADED(threading[0]->n_core)
    tm_term_double_setup( g.mu, g.mu_even_shift, g.mu_odd_shift, g.mu_factor[0], &(g.op_double), &l, threading[omp_get_thread_num()]);
    
  } else if (g.conf_flag && g.setup_flag )
    THREADED(threading[0]->n_core) {
      tm_term_update( g.mu, &l, threading[omp_get_thread_num()] );
      finalize_operator_update( &l, threading[omp_get_thread_num()] );
    }
  mg_status->info = g.mu;
#else
  warning0("DDalphaAMG_change_mu_sign called without the flag HAVE_TM enabled. Doing nothing.\n");
  mg_status->info = 0;
#endif

  t1 = MPI_Wtime();
  
  mg_status->success = 1;
  mg_status->time = t1-t0;
  //why here; no setup is actually done here???????
  for ( int i=0; i<g.num_levels; i++ ) mg_status->iter_times[i] = g.iter_times[i];
  
}

void DDalphaAMG_set_configuration( double *gauge_field, DDalphaAMG_status *mg_status ) {
  
  int t, z, y, x, mu, i, j, k;
  double t0, t1;
  SU3_storage U=NULL;
  complex_double phase[4] = { _COMPLEX_double_ZERO, _COMPLEX_double_ZERO, _COMPLEX_double_ZERO, _COMPLEX_double_ZERO};
  int *ll=l.local_lattice, onb[4];
  for (i=0; i<4; i++)
    onb[i] = (g.my_coords[i]==g.process_grid[i]-1)?1:0;

  t0 = MPI_Wtime();
  mg_status->success = 0;
  mg_status->info = 0;  

  // START: dirac_setup
  if ( g.print > 0 ) printf0("%s\n", CLIFFORD_BASIS );
  if ( g.bc == _ANTIPERIODIC ) printf0("antiperiodic in time");
  else if ( g.bc == _TWISTED ) printf0("twisted (%.2f, %.2f, %.2f, %.2f)", g.twisted_bc[0], 
               g.twisted_bc[1], g.twisted_bc[2], g.twisted_bc[3]);
  else printf0("periodic in time");
  printf0(" boundary conditions \n");

  SU3_storage_alloc( &U, &l );
  
  if(g.bc == _ANTIPERIODIC && onb[T] ) {
    phase[Z] = 1; phase[Y] = 1; phase[X] = 1;
    for ( t=1, i=0, k=0; t<ll[T]+1; t++ ) {
      if (t<ll[T]) phase[T] = 1; 
      else phase[T] = -1;
      for ( z=1; z<ll[Z]+1; z++ )
        for ( y=1; y<ll[Y]+1; y++ )
          for ( x=1; x<ll[X]+1; x++ )
            for ( mu=0; mu<4; mu++ ) {
              if ( conf_index_fct != NULL )
                k = conf_index_fct( t-1, z-1, y-1, x-1, mu );
              for (j=0; j<9; j++, i++, k+=2) {
                g.op_double.D[i] = 0.5*phase[mu]*(gauge_field[k]+I*gauge_field[k+1]);
                U[t][z][y][x][mu][j] = phase[mu]*(gauge_field[k]+I*gauge_field[k+1]);
              }
            }     
    }
  }
  else if(g.bc == _TWISTED && ( onb[T] || onb[Z] || onb[Y] || onb[X] ))
    for ( t=1, i=0, k=0; t<ll[T]+1; t++ ) {
      if ( !onb[T] || t<ll[T] || g.twisted_bc[T]==0) phase[T] = 1; 
      else phase[T] = cexp(I*g.twisted_bc[T]);
      for ( z=1; z<ll[Z]+1; z++ ) {
        if ( !onb[Z] || z<ll[Z] || g.twisted_bc[Z]==0) phase[Z] = 1; 
        else phase[Z] = cexp(I*g.twisted_bc[Z]);
        for ( y=1; y<ll[Y]+1; y++ ) {
          if ( !onb[Y] || y<ll[Y] || g.twisted_bc[Y]==0) phase[Y] = 1; 
          else phase[Y] = cexp(I*g.twisted_bc[Y]);
          for ( x=1; x<ll[X]+1; x++ ) {
            if ( !onb[X] || x<ll[X] || g.twisted_bc[X]==0) phase[X] = 1; 
            else phase[X] = cexp(I*g.twisted_bc[X]);
            for ( mu=0; mu<4; mu++ ) {
              if ( conf_index_fct != NULL )
                k = conf_index_fct( t-1, z-1, y-1, x-1, mu );
              for (j=0; j<9; j++, i++, k+=2) {
                g.op_double.D[i] = 0.5*phase[mu]*(gauge_field[k]+I*gauge_field[k+1]);
                U[t][z][y][x][mu][j] = phase[mu]*(gauge_field[k]+I*gauge_field[k+1]);
              }
            }
          }
        }
      }
    }
  
  else
    for ( t=1, i=0, k=0; t<ll[T]+1; t++ )
      for ( z=1; z<ll[Z]+1; z++ )
        for ( y=1; y<ll[Y]+1; y++ )
          for ( x=1; x<ll[X]+1; x++ )
            for ( mu=0; mu<4; mu++ ) {
              if ( conf_index_fct != NULL )
                k = conf_index_fct( t-1, z-1, y-1, x-1, mu );
              for (j=0; j<9; j++, i++, k+=2) {
                g.op_double.D[i] = 0.5*(gauge_field[k]+I*gauge_field[k+1]);
                U[t][z][y][x][mu][j] = (gauge_field[k]+I*gauge_field[k+1]);
              }
            }
  
  SU3_ghost_update( &U, &l );
  if ( g.print > 0 ) printf0("Configuration stored...\n");
  
  compute_clover_term( U, &l );
  
  // calculate the plaquette
  g.plaq_clov = calc_plaq( U, &l );
  if(g.print > 0) printf0("average plaquette: %.13lf\n", g.plaq_clov);
    
  SU3_storage_free( &U, &l );
  //END: dirac_setup
  
  mg_status->success = 1;
  g.conf_flag = 1;
  mg_status->info = g.plaq;
  
  if ( g.setup_flag == 1 ) {
    if(l.s_double.op.clover != NULL)
      schwarz_double_setup( &(l.s_double), &(g.op_double), &l );
    if(l.s_float.op.clover != NULL)
      schwarz_float_setup( &(l.s_float), &(g.op_double), &l );
    THREADED(threading[0]->n_core)
    if ( g.mixed_precision ) 
      operator_updates_float( &l, threading[omp_get_thread_num()] );
    else
      operator_updates_double( &l, threading[omp_get_thread_num()] );
    
    mg_status->success++; //0: error, 1: OK, 2 re_setup done
  }    

  t1 = MPI_Wtime();
  mg_status->time = t1-t0;
  for ( i=0; i<g.num_levels; i++ ) mg_status->iter_counts[i] = 0;
  for ( i=0; i<g.num_levels; i++ ) mg_status->iter_times[i]  = 0;
  
}

void DDalphaAMG_setup( DDalphaAMG_status * mg_status ) {

  int i;
  double t0, t1;
  t0 = MPI_Wtime();
  for ( i=0; i<g.num_levels; i++ ) g.iter_counts[i] = 0;
  for ( i=0; i<g.num_levels; i++ ) g.iter_times[i]  = 0;
  mg_status->success = 0;
  mg_status->info = 0;
  
  if(g.conf_flag == 1) {
    if ( g.setup_flag )
      method_free( &l );
    THREADED(threading[0]->n_core)
    {
      method_setup( NULL, &l, threading[omp_get_thread_num()] );
      method_iterative_setup( g.setup_iter[0], &l, threading[omp_get_thread_num()] );
    }
    g.setup_flag = 1;
        
    l.setup_iter = g.setup_iter[0];
    mg_status->success = l.setup_iter;
    
  }
  
  t1 = MPI_Wtime();
  mg_status->time = t1-t0;
  for ( i=0; i<g.num_levels; i++ ) mg_status->iter_counts[i] = g.iter_counts[i];
  for ( i=0; i<g.num_levels; i++ ) mg_status->iter_times[i]  = g.iter_times[i];
}

void DDalphaAMG_update_setup( int iterations, DDalphaAMG_status * mg_status ) {

  if(g.setup_flag) {
    double t0, t1;
    t0 = MPI_Wtime();
    for ( int i=0; i<g.num_levels; i++ ) g.iter_counts[i] = 0;
    for ( int i=0; i<g.num_levels; i++ ) g.iter_times[i]  = 0;
    mg_status->success = 0;
    mg_status->info = 0;

    THREADED(threading[0]->n_core) 
    method_iterative_setup( iterations, &l, threading[omp_get_thread_num()] );
    
    l.setup_iter += iterations;
    mg_status->success = l.setup_iter;
    t1 = MPI_Wtime();
    mg_status->time = t1-t0;
    for ( int i=0; i<g.num_levels; i++ ) mg_status->iter_counts[i] = g.iter_counts[i];
    for ( int i=0; i<g.num_levels; i++ ) mg_status->iter_times[i]  = g.iter_times[i];
  }
  else {
    g.setup_iter[0] = iterations;
    DDalphaAMG_setup( mg_status );
  }
}

/*************  Solver Utilities  **************************************************/
// Given input, solve each num_loop chunk, and put them in output
static inline void solver( vector_double *vector_out, vector_double *vector_in )
{
  int nvec = vector_in->num_vect, vfac = nvec/g.num_rhs_vect;
  vector_double rhs = (g.mixed_precision==2 && g.method >= 0)?g.p_MP.dp.b:g.p.b;
  vector_double sol = (g.mixed_precision==2 && g.method >= 0)?g.p_MP.dp.x:g.p.x;
  
  for ( int i=0; i<g.num_rhs_vect; i+=num_loop ) {
    g.n_chunk = i/num_loop;
    THREADED(threading[0]->n_core) {
      int start = threading[omp_get_thread_num()]->start_index[0];
      int end = threading[omp_get_thread_num()]->end_index[0];
      vector_double_copy2( &rhs, vector_in, vfac*i, nvec, 1, start, end, &l );
      if ( g.method == -1 ) {
	cgn_double( &(g.p), &l, threading[omp_get_thread_num()] );
      } else if ( g.mixed_precision == 2 ) {
	fgmres_MP( &(g.p_MP), &l, threading[omp_get_thread_num()] );
      } else {
	solver_double( &(g.p), &l, threading[omp_get_thread_num()] );
      }
      vector_double_copy2( vector_out, &sol, vfac*i, nvec, -1, start, end, &l );
    }
  }
}

static inline void apply_preconditioner( vector_double *vector_out, vector_double *vector_in )
{
  int nvec = vector_in->num_vect, vfac = nvec/g.num_rhs_vect;
  vector_double rhs = (g.mixed_precision==2 && g.method >= 0)?g.p_MP.dp.b:g.p.b;
  vector_double sol = (g.mixed_precision==2 && g.method >= 0)?g.p_MP.dp.x:g.p.x;
  
  for ( int i=0; i<g.num_rhs_vect; i+=num_loop ) {
    g.n_chunk = i/num_loop;
    THREADED(threading[0]->n_core) {
      int start = threading[omp_get_thread_num()]->start_index[0];
      int end = threading[omp_get_thread_num()]->end_index[0];
      vector_double_copy2( &rhs, vector_in, vfac*i, nvec, 1, start, end, &l );
      preconditioner( &sol, NULL, &rhs, _NO_RES, &l, threading[omp_get_thread_num()] );
      vector_double_copy2( vector_out, &sol, vfac*i, nvec, -1, start, end, &l );
    }
  }
}

static inline void apply_operator( vector_double *vector_out, vector_double *vector_in )
{
  int nvec = vector_in->num_vect, vfac = nvec/g.num_rhs_vect;
  vector_double rhs = (g.mixed_precision==2 && g.method >= 0)?g.p_MP.dp.b:g.p.b;
  vector_double sol = (g.mixed_precision==2 && g.method >= 0)?g.p_MP.dp.x:g.p.x;
  
  for ( int i=0; i<g.num_rhs_vect; i+=num_loop ) {
    g.n_chunk = i/num_loop;
    THREADED(threading[0]->n_core) {
      int start = threading[omp_get_thread_num()]->start_index[0];
      int end = threading[omp_get_thread_num()]->end_index[0];
      vector_double_copy2( &rhs, vector_in, vfac*i, nvec, 1, start, end, &l );
      if ( g.mixed_precision == 2 ) {
	apply_operator_double( &sol, &rhs, &(g.p_MP.dp), &l, threading[omp_get_thread_num()] );
      } else {
	apply_operator_double( &sol, &rhs, &(g.p), &l, threading[omp_get_thread_num()] );
      }
      vector_double_copy2( vector_out, &sol, vfac*i, nvec, -1, start, end, &l );
    }
  }
}


static inline void change_epsbar_shift_sign( ) {
  
#ifdef HAVE_TM1p1  
  if ( g.epsbar_ig5_even_shift !=0 || g.epsbar_ig5_odd_shift !=0 ) {
    g.epsbar_ig5_even_shift *= -1;
    g.epsbar_ig5_odd_shift *= -1;

    if (g.conf_flag && !g.setup_flag ) {
      
      THREADED(threading[0]->n_core) {
        epsbar_term_double_setup( g.epsbar, g.epsbar_ig5_even_shift, g.epsbar_ig5_odd_shift, &(g.op_double),
                                  &l, threading[omp_get_thread_num()]);
      }
    } else if (g.conf_flag && g.setup_flag )
      THREADED(threading[0]->n_core) {
        epsbar_term_update( &l, threading[omp_get_thread_num()] );
        finalize_operator_update( &l, threading[omp_get_thread_num()] );
      }
  }
#else
  warning0("change_epsbar_shift_sign called without the flag HAVE_TM1p1 enabled. Doing nothing.\n");
#endif
}

enum {_SOLVE, _SOLVE_SQ, _SOLVE_SQ_ODD, _SOLVE_SQ_EVEN, _PRECOND, _OPERATOR, _RESTRICT, _PROLONGATE};

// NOTE RESIDUAL
//
// The _SOLVE_SQ invert the squared operator in two inversion.
// One has to be careful to return a solution with the right residual.
//
// We have:
//   D^2 = Dd D
//
//   Dd D x = b   direct solution
//   Dd  y  = b   first step
//      D x = y   second step
//
//   r1 = Dd   y - b
//   r2 =    D x - y
//   r  = Dd D x - b = Dd r2 + r1
//
//  |r| < tol -->  |r| < |Dd| |r2| + |r1| < tol
//
// For the residual we do
//  |r1| < tol/2
//  |r2| < (tol - |r1|)/|Dd| using |Dd|=8 since is |Dd|<8
//
// With relative residual we have
//  |r1|/|b| < tol/2
//  |r2|/|y| < (tol - |r1|)/|Dd|*|b|/|y|

static inline void DDalphaAMG_driver( double *vector1_out, double *vector1_in, double *vector2_out, double *vector2_in, double tol, DDalphaAMG_status *mg_status, int _TYPE ) {
  /*
   * Assume: in&out vectors are at the top level
   * vector2_out & vector2_in are used when HAVE_TM1p1
   */
  
  int t, z, y, x, i, j, k, n, mu;
  int *ll = l.local_lattice, *gl=l.global_lattice, sl[4], precision_changed, nrhs = g.num_rhs_vect, vsize = l.inner_vector_size, vfac = 1, nrhs_sf = nrhs;
  complex_double twisted_bc, tmp;
  double phase[4] = {_COMPLEX_double_ZERO, _COMPLEX_double_ZERO, _COMPLEX_double_ZERO, _COMPLEX_double_ZERO}, vmin = 1, vmax = EPS_float, vtmp;
  gmres_double_struct *p = g.mixed_precision==2?&(g.p_MP.dp):&(g.p);
  vector_double rhs, sol;
#ifdef HAVE_TM1p1
  nrhs_sf *=  g.n_flavours/g.num_indep_flav;
#endif
  double norm[nrhs_sf], norm2[nrhs_sf];
  DDalphaAMG_status tmp_status;

  double t0, t1;
  t0 = MPI_Wtime();
  for ( i=0; i<g.num_levels; i++ ) g.iter_counts[i] = 0;
  for ( i=0; i<g.num_levels; i++ ) g.iter_times[i]  = 0;
  mg_status->success = 0;
  mg_status->info = 0;

  ASSERT(vector1_out!=NULL);
  ASSERT(vector1_in!=NULL);
#ifdef HAVE_TM1p1
  if(g.n_flavours==2) {
    vfac *= 2;
    ASSERT(vector2_out!=NULL);
    ASSERT(vector2_in!=NULL);
  }
#endif

  vector_double_init( &rhs );
  vector_double_init( &sol );
  vector_double_alloc( &rhs, _INNER, vfac*nrhs, &l, no_threading );
  vector_double_alloc( &sol, _INNER, vfac*nrhs, &l, no_threading );
  rhs.num_vect_now = vfac*nrhs;
  sol.num_vect_now = vfac*nrhs;
  
  if(g.mixed_precision!=2)
    g.p.tol = tol;
  else
    g.p_MP.dp.tol = tol;

  for (i=0; i<4; i++)
    sl[i] = ll[i]*g.my_coords[i];

  // i: the starting index of the chunk (internal d.o.f. x real/imag) for each site in the given vector
  // j: index for sits in lexicographic order & internal d.o.f. (color faster than spin)
  buffer_double rhs_pt = rhs.vector_buffer, sol_pt = sol.vector_buffer;
  for (t=0, j=0; t<ll[T]; t++) {
    if (g.bc==_TWISTED) phase[T] = g.twisted_bc[T]*((double)sl[T]+t)/(double)gl[T];
    for (z=0; z<ll[Z]; z++) {
      if (g.bc==_TWISTED) phase[Z] = phase[T] + g.twisted_bc[Z]*((double)sl[Z]+z)/(double)gl[Z];
      for (y=0; y<ll[Y]; y++) {
        if (g.bc==_TWISTED) phase[Y] = phase[Z] + g.twisted_bc[Y]*((double)sl[Y]+y)/(double)gl[Y];
        for (x=0; x<ll[X]; x++) {
          if (g.bc==_TWISTED) {
            phase[X] = phase[Y] + g.twisted_bc[X]*((double)sl[X]+x)/(double)gl[X];
            twisted_bc = cexp(I*phase[X]);
          } else
            twisted_bc = 1.;
          if(vector_index_fct!=NULL )
            i = vector_index_fct( t, z, y, x );
          else 
            i = 2*j;
          
	  for ( n=0; n<nrhs; n++ ) 
            for ( mu=0; mu<4; mu++ )
              for ( k=0; k<3; k++, j++ ) {
		/* TODO:
		   to allow the possibility of input and output in _NVEC_INNER, just compute index first, then do the assignment
		 */
#ifndef BASIS4 
		rhs_pt[j*vfac*nrhs+n] = ((complex_double)vector1_in[i+2*(k+3*mu)+n*vsize] + I*(complex_double)vector1_in[i+2*(k+3*mu)+1+n*vsize]) * twisted_bc;
#else
		rhs_pt[j*vfac*nrhs+n] = ((complex_double)vector1_in[i+2*(k+3*(3-mu)+n*vsize)] + I*(complex_double)vector1_in[i+2*(k+3*(3-mu))+1+n*vsize]) * twisted_bc;
#endif
#ifdef HAVE_TM1p1
		if(g.n_flavours==2) {
#ifndef BASIS4
		  rhs_pt[j*vfac*nrhs+n+nrhs] = ((complex_double)vector2_in[i+2*(k+3*mu)+n*vsize] + I*(complex_double)vector2_in[i+2*(k+3*mu)+1+n*vsize]) * twisted_bc;
#else
		  rhs_pt[j*vfac*nrhs+n+nrhs] = ((complex_double)vector2_in[i+2*(k+3*(3-mu)+n*vsize)] + I*(complex_double)vector2_in[i+2*(k+3*(3-mu))+1+n*vsize]) * twisted_bc;
#endif
		}
#endif
		if(p->initial_guess_zero == 0) {
#ifndef BASIS4
                  sol_pt[j*vfac*nrhs+n] = ((complex_double)vector1_out[i+2*(k+3*mu)+n*vsize] + I*(complex_double)vector1_out[i+2*(k+3*mu)+1+n*vsize]) * twisted_bc;
#else
                  sol_pt[j*vfac*nrhs+n] = ((complex_double)vector1_out[i+2*(k+3*(3-mu)+n*vsize)] + I*(complex_double)vector1_out[i+2*(k+3*(3-mu))+1+n*vsize]) * twisted_bc;
#endif
#ifdef HAVE_TM1p1
		  if(g.n_flavours==2) {
#ifndef BASIS4
		    sol_pt[j*vfac*nrhs+n+nrhs] = ((complex_double)vector2_out[i+2*(k+3*mu)+n*vsize] + I*(complex_double)vector2_out[i+2*(k+3*mu)+1+n*vsize]) * twisted_bc;
#else
		    sol_pt[j*vfac*nrhs+n+nrhs] = ((complex_double)vector2_out[i+2*(k+3*(3-mu)+n*vsize)] + I*(complex_double)vector2_out[i+2*(k+3*(3-mu))+1+n*vsize]) * twisted_bc;
#endif
		  }
#endif
		}
		
#ifndef INIT_ONE_PREC
		if(g.mixed_precision==2) {
		  vtmp=cabs(rhs_pt[j*vfac*nrhs+n]);
		  if(vtmp > vmax)
		    vmax=vtmp;
		  if( vtmp > EPS_double && vtmp < vmin )
		    vmin=vtmp;
#ifdef HAVE_TM1p1
		  if(g.n_flavours==2) {
		    vtmp=cabs(rhs_pt[j*vfac*nrhs+n+nrhs]);
		    if(vtmp > vmax)
			vmax=vtmp;
		    if( vtmp > EPS_double && vtmp < vmin )
		      vmin=vtmp;
		  }
#endif
                }
              
#endif
              }
        }
      }
    }
  }
   
#ifndef INIT_ONE_PREC
    
  double gvmin, gvmax;
  if(g.mixed_precision==2) {
    MPI_Allreduce(&vmin, &gvmin, 1, MPI_DOUBLE, MPI_MIN, g.comm_cart);
    MPI_Allreduce(&vmax, &gvmax, 1, MPI_DOUBLE, MPI_MAX, g.comm_cart);
  }
  
  //switching to double precision on the fine level
  if(g.mixed_precision==2 && gvmin/gvmax<EPS_float) {
    warning0("Changing solver precision on fine level due to rhs elements (min/max=%e)\n", vmin/vmax);
    precision_changed=1;
    g.mixed_precision=1;
    p = &(g.p);
    p->tol = g.p_MP.dp.tol;
  } else precision_changed = 0;
#endif

  switch(_TYPE) {
    
  case _SOLVE :
    solver( &sol, &rhs );
    break;

  case _SOLVE_SQ :
    THREADED(threading[0]->n_core)
#ifdef HAVE_TM1p1
      if(g.n_flavours==2) 
        // sol = (D_h^{-1})*g5*tau1*(D_h^{-1})*g5*tau1*rhs
        tau1_gamma5_double( &rhs, &rhs, &l, threading[omp_get_thread_num()] );
      else
#endif
        // sol = (D_d^{-1})*g5*(D_u^{-1})*g5*rhs
        gamma5_double( &rhs, &rhs, &l, threading[omp_get_thread_num()] );
    
    // read NOTE RESIDUAL
    THREADED(threading[0]->n_core)
      global_norm_double( norm, &rhs, p->v_start, p->v_end, &l, threading[omp_get_thread_num()] );
    p->tol = tol/2.;
    solver( &sol, &rhs );
      
    THREADED(threading[0]->n_core)
#ifdef HAVE_TM1p1
      if(g.n_flavours==2) 
        tau1_gamma5_double(&rhs, &sol, &l, threading[omp_get_thread_num()] );
      else
#endif
        gamma5_double(&rhs, &sol, &l, threading[omp_get_thread_num()] );
 
#ifdef HAVE_TM1p1
    if(g.n_flavours==2) 
      change_epsbar_shift_sign( );
    else
#endif
      DDalphaAMG_change_mu_sign( &tmp_status );

    // read NOTE RESIDUAL
    THREADED(threading[0]->n_core)
      global_norm_double( norm2, &rhs, p->v_start, p->v_end, &l, threading[omp_get_thread_num()] );
    p->tol = 0;
    for ( i=0, p->tol=0; i<nrhs_sf; i++ ) p->tol += (tol-g.resids[i])*norm[i]/norm2[i]/8./nrhs_sf; 
    solver( &sol, &rhs );

    // DDalphaAMG_change_mu_sign( &tmp_status );
    warning0("sign of mu changed during the inversion of squared operator\n");
    break;
    
  case _SOLVE_SQ_ODD :    
    THREADED(threading[0]->n_core)
#ifdef HAVE_TM1p1
      if(g.n_flavours==2) 
        // sol = (D_h^{-1})*g5*tau1*(D_h^{-1})*g5*tau1*rhs
        tau1_gamma5_set_even_to_zero_double(&rhs, &rhs, &l, threading[omp_get_thread_num()]);
      else
#endif
        // sol = (D_d^{-1})*g5*(D_u^{-1})*g5*rhs
        gamma5_set_even_to_zero_double(&rhs, &rhs, &l, threading[omp_get_thread_num()]);

    // read NOTE RESIDUAL
    THREADED(threading[0]->n_core)
      global_norm_double( norm, &rhs, p->v_start, p->v_end, &l, threading[omp_get_thread_num()] );
    p->tol = tol/2.;
    solver( &sol, &rhs );

    THREADED(threading[0]->n_core)
#ifdef HAVE_TM1p1
      if(g.n_flavours==2) 
        tau1_gamma5_set_even_to_zero_double(&rhs, &sol, &l, threading[omp_get_thread_num()]);
      else
#endif
        gamma5_set_even_to_zero_double(&rhs, &sol, &l, threading[omp_get_thread_num()]);
 
#ifdef HAVE_TM1p1
    if(g.n_flavours==2) 
      change_epsbar_shift_sign( );
    else
#endif
      DDalphaAMG_change_mu_sign( &tmp_status );

    // read NOTE RESIDUAL
    THREADED(threading[0]->n_core)
      global_norm_double( norm2, &rhs, p->v_start, p->v_end, &l, threading[omp_get_thread_num()] );
    p->tol = 0;
    for ( i=0, p->tol=0; i<nrhs_sf; i++ ) p->tol += (tol-g.resids[i])*norm[i]/norm2[i]/8./nrhs_sf;
    solver( &sol, &rhs );

    // DDalphaAMG_change_mu_sign( &tmp_status );
    warning0("sign of mu changed during the inversion of squared operator\n");
    break;
    
  case _SOLVE_SQ_EVEN :    
    THREADED(threading[0]->n_core)
#ifdef HAVE_TM1p1
      if(g.n_flavours==2) 
        // sol = (D_h^{-1})*g5*tau1*(D_h^{-1})*g5*tau1*rhs
        tau1_gamma5_set_odd_to_zero_double(&rhs, &rhs, &l, threading[omp_get_thread_num()]);
      else
#endif
        // sol = (D_d^{-1})*g5*(D_u^{-1})*g5*rhs
        gamma5_set_odd_to_zero_double(&rhs, &rhs, &l, threading[omp_get_thread_num()]);

    // read NOTE RESIDUAL
    THREADED(threading[0]->n_core)
    global_norm_double( norm, &rhs, p->v_start, p->v_end, &l, threading[omp_get_thread_num()] );
    p->tol = tol/2.;
    solver( &sol, &rhs );

    THREADED(threading[0]->n_core)
#ifdef HAVE_TM1p1
      if(g.n_flavours==2) 
        tau1_gamma5_set_odd_to_zero_double(&rhs, &sol, &l, threading[omp_get_thread_num()]);
      else
#endif
        gamma5_set_odd_to_zero_double(&rhs, &sol, &l, threading[omp_get_thread_num()]);

#ifdef HAVE_TM1p1
    if(g.n_flavours==2) 
      change_epsbar_shift_sign( );
    else
#endif
      DDalphaAMG_change_mu_sign( &tmp_status );

    // read NOTE RESIDUAL
    THREADED(threading[0]->n_core)
      global_norm_double( norm2, &rhs, p->v_start, p->v_end, &l, threading[omp_get_thread_num()] );
    p->tol = 0;
    for ( i=0, p->tol=0; i<nrhs_sf; i++ ) p->tol += (tol-g.resids[i])*norm[i]/norm2[i]/8./nrhs_sf;
    solver( &sol, &rhs );

    // DDalphaAMG_change_mu_sign( &tmp_status );
    warning0("sign of mu changed during the inversion of squared operator\n");
    break;

  case _PRECOND :
    apply_preconditioner( &sol, &rhs );
    break;

  case _OPERATOR :
    apply_operator( &sol, &rhs );
    break;

  default :
    warning0("_TYPE not found in DDalphaAMG_driver. Returing vector in as vector out.");
    sol=rhs;
    break;
  }

  sol_pt = sol.vector_buffer;
  for (t=0, j=0; t<ll[T]; t++) {
    if (g.bc==_TWISTED) phase[T] = g.twisted_bc[T]*((double)sl[T]+t)/(double)gl[T];
    for (z=0; z<ll[Z]; z++) {
      if (g.bc==_TWISTED) phase[Z] = phase[T] + g.twisted_bc[Z]*((double)sl[Z]+z)/(double)gl[Z];
      for (y=0; y<ll[Y]; y++) {
        if (g.bc==_TWISTED) phase[Y] = phase[Z] + g.twisted_bc[Y]*((double)sl[Y]+y)/(double)gl[Y];
        for (x=0; x<ll[X]; x++) {
          if (g.bc==_TWISTED) {
            phase[X] = phase[Y] + g.twisted_bc[X]*((double)sl[X]+x)/(double)gl[X];
            twisted_bc = cexp(-I*phase[X]);
          } else
            twisted_bc = 1.;
          if(vector_index_fct!=NULL )
            i = vector_index_fct( t, z, y, x );
          else 
            i = 2*j;

	  for ( n=0; n<nrhs; n++ )
            for ( mu=0; mu<4; mu++ )
              for ( k=0; k<3; k++, j++ ) {
                tmp = sol_pt[j*vfac*nrhs+n] * twisted_bc;
#ifndef BASIS4 
                vector1_out[i+2*(k+3*mu)+n*vsize]   = creal(tmp);
                vector1_out[i+2*(k+3*mu)+1+n*vsize] = cimag(tmp);
#else
                vector1_out[i+2*(k+3*(3-mu))+n*vsize]   = creal(tmp);
                vector1_out[i+2*(k+3*(3-mu))+1+n*vsize] = cimag(tmp);
#endif
#ifdef HAVE_TM1p1
		if(g.n_flavours==2) {
		  tmp = sol_pt[j*vfac*nrhs+n+nrhs] * twisted_bc;
#ifndef BASIS4
		  vector2_out[i+2*(k+3*mu)+n*vsize]   = creal(tmp);
		  vector2_out[i+2*(k+3*mu)+1+n*vsize] = cimag(tmp);
#else
		  vector2_out[i+2*(k+3*(3-mu))+n*vsize]   = creal(tmp);
		  vector2_out[i+2*(k+3*(3-mu))+1+n*vsize] = cimag(tmp);
#endif
		}
#endif
              }
        }
      }
    }
  }
    
#ifndef INIT_ONE_PREC
  if (precision_changed)
    g.mixed_precision=2;
#endif

  vector_double_free( &rhs, &l, no_threading );
  vector_double_free( &sol, &l, no_threading );
  
  if ( g.max_rel_res_norm <= tol || _TYPE == _OPERATOR || _TYPE == _PRECOND )
    mg_status->success = 1;
  mg_status->info = g.max_rel_res_norm;
  t1 = MPI_Wtime();
  mg_status->time = t1-t0;
  for ( i=0; i<g.num_levels; i++ ) mg_status->iter_counts[i] = g.iter_counts[i];
  for ( i=0; i<g.num_levels; i++ ) mg_status->iter_times[i] = g.iter_times[i];
  
}

// this is mainly intended for coarse level operations or interpolation/restriction
static inline void DDalphaAMG_proj_driver( float *vector_out, float *vector_in, int level, DDalphaAMG_status *mg_status, int _TYPE ) {
  /*
   * Assume: in & out vectors are single flavor vectors
   */
  
  int t, z, y, x, i, j, k, n, mu, *ll, vsize = l.inner_vector_size;

  level_struct *ltmp=&l, *from, *to;
  for(i=0; i<level; i++)
    ltmp=ltmp->next_level;

  if(_TYPE==_RESTRICT) {
    from=ltmp;
    to=ltmp->next_level;
  } else if(_TYPE==_PROLONGATE) {
    from=ltmp->next_level;
    to=ltmp;    
  } else {
    from=ltmp;
    to=ltmp;    
  }
  vector_float *rhs = &(from->p_float.b);
  vector_float *sol = &(to->p_float.x);

  double t0, t1;
  t0 = MPI_Wtime();
  for ( i=0; i<g.num_levels; i++ ) g.iter_counts[i] = 0;
  for ( i=0; i<g.num_levels; i++ ) g.iter_times[i]  = 0;
  mg_status->success = 0;
  mg_status->info = 0;
  
  ASSERT(g.mixed_precision);
  ASSERT(vector_out!=NULL);
  ASSERT(vector_in!=NULL);

  // Repeat the application of the chosen operator for each chunk of num_loop
  for ( k=0; k<g.num_rhs_vect; k+=num_loop ) {
    g.n_chunk = k/num_loop;
    ll = from->local_lattice;
    for (t=0, j=0; t<ll[T]; t++) 
      for (z=0; z<ll[Z]; z++) 
	for (y=0; y<ll[Y]; y++) 
	  for (x=0; x<ll[X]; x++) {
	    if(vector_index_fct!=NULL )
	      i = vector_index_fct( t, z, y, x );
	    else // for the moment, I assume in:_NVEC_OUTER and rearrange to _NVEC_INNER
	      i = 2*j;
	    
	    for ( mu=0; mu<from->num_lattice_site_var; mu++, j++ )
	      for ( n=0; n<num_loop; n++ )
		rhs->vector_buffer[j*rhs->num_vect+n] = (vector_in[i+2*mu+(k+n)*vsize] + I*vector_in[i+2*mu+1+(k+n)*vsize]);
	  }
    
    switch(_TYPE) {
      
    case _RESTRICT :
      THREADED(threading[0]->n_core)
	restrict_float( sol, rhs, from, threading[omp_get_thread_num()] );
      break;
      
    case _PROLONGATE :
      THREADED(threading[0]->n_core)
	interpolate3_float( sol, rhs, to, threading[omp_get_thread_num()] );
      break;
      
    case _OPERATOR :
      THREADED(threading[0]->n_core)
	if ( from->level==0 && g.odd_even )
	  coarse_odd_even_float_test( sol, rhs, from, threading[omp_get_thread_num()] );
	else
	  apply_operator_float( sol, rhs, &(from->p_float), from, threading[omp_get_thread_num()] );
      break;

    default :
      warning0("_TYPE not found in DDalphaAMG_driver. Returing vector in as vector out.");
      sol=rhs;
      break;
    }
  
    ll = to->local_lattice;
    for (t=0, j=0; t<ll[T]; t++) 
      for (z=0; z<ll[Z]; z++) 
	for (y=0; y<ll[Y]; y++) 
	  for (x=0; x<ll[X]; x++) {
	    if(vector_index_fct!=NULL )
	      i = vector_index_fct( t, z, y, x );
	    else 
	      i = 2*j;

	    for ( mu=0; mu<to->num_lattice_site_var; mu++, j++ )
	      for ( n=0; n<num_loop; n++ ) {
		vector_out[i+2*mu+(k+n)*vsize]   = creal(sol->vector_buffer[j*sol->num_vect+n]);
		vector_out[i+2*mu+1+(k+n)*vsize] = cimag(sol->vector_buffer[j*sol->num_vect+n]);
	      }
	  }
  }

    
  mg_status->success = 1;
  mg_status->info = g.max_rel_res_norm;
  t1 = MPI_Wtime();
  mg_status->time = t1-t0;
  for ( i=0; i<g.num_levels; i++ ) mg_status->iter_counts[i] = g.iter_counts[i];
  for ( i=0; i<g.num_levels; i++ ) mg_status->iter_times[i] = g.iter_times[i];
  
}


static inline void set_n_flavours( int n) {

#ifdef HAVE_TM1p1
  if( n==1 )
    g.force_2flavours = 0;
  else
    g.force_2flavours = 1;
  THREADED(threading[0]->n_core)
    data_layout_n_flavours( n, &l, threading[omp_get_thread_num()] );
#else
  if( n==2 )
      error0("For DDalphaAMG_solve_doublet_*, HAVE_TM1p1 flag required\n");
#endif
    
}

void DDalphaAMG_solve( double *vector_out, double *vector_in, double tol, DDalphaAMG_status *mg_status )
{
  DDalphaAMG_driver( vector_out, vector_in, NULL, NULL, tol, mg_status, _SOLVE );
}

void DDalphaAMG_solve_doublet( double *vector1_out, double *vector1_in,
                               double *vector2_out, double *vector2_in,
                               double tol, DDalphaAMG_status *mg_status )
{
  set_n_flavours( 2 );
  DDalphaAMG_driver( vector1_out, vector1_in, vector2_out, vector2_in, tol, mg_status, _SOLVE );
  set_n_flavours( 1 );
}

void DDalphaAMG_solve_doublet_with_guess( double *vector1_out, double *vector1_in,
                                          double *vector2_out, double *vector2_in,
                                          double tol, DDalphaAMG_status *mg_status )
{
  gmres_double_struct *p = g.mixed_precision==2?&(g.p_MP.dp):&(g.p);
  set_n_flavours( 2 );
  p->initial_guess_zero = 0;
  DDalphaAMG_driver( vector1_out, vector1_in, vector2_out, vector2_in, tol, mg_status, _SOLVE );
  p->initial_guess_zero = 1;
  set_n_flavours( 1 );
}

void DDalphaAMG_solve_squared( double *vector_out, double *vector_in, double tol, DDalphaAMG_status *mg_status )
{
  DDalphaAMG_driver( vector_out, vector_in, NULL, NULL, tol, mg_status, _SOLVE_SQ );
}

void DDalphaAMG_solve_doublet_squared( double *vector1_out, double *vector1_in,
                                       double *vector2_out, double *vector2_in,
                                       double tol, DDalphaAMG_status *mg_status )
{
  set_n_flavours( 2 );
  DDalphaAMG_driver( vector1_out, vector1_in, vector2_out, vector2_in, tol, mg_status, _SOLVE_SQ );
  set_n_flavours( 1 );
}

void DDalphaAMG_solve_squared_odd( double *vector_out, double *vector_in, double tol, DDalphaAMG_status *mg_status )
{
  DDalphaAMG_driver( vector_out, vector_in, NULL, NULL, tol, mg_status, _SOLVE_SQ_ODD );
}

void DDalphaAMG_solve_doublet_squared_odd( double *vector1_out, double *vector1_in,
                                           double *vector2_out, double *vector2_in,
                                           double tol, DDalphaAMG_status *mg_status )
{
  set_n_flavours( 2 );
  DDalphaAMG_driver( vector1_out, vector1_in, vector2_out, vector2_in, tol, mg_status, _SOLVE_SQ_ODD );
  set_n_flavours( 1 );
}

void DDalphaAMG_solve_squared_even( double *vector_out, double *vector_in, double tol, DDalphaAMG_status *mg_status )
{
  DDalphaAMG_driver( vector_out, vector_in, NULL, NULL, tol, mg_status, _SOLVE_SQ_EVEN );
}

void DDalphaAMG_solve_doublet_squared_even( double *vector1_out, double *vector1_in,
                                            double *vector2_out, double *vector2_in,
                                            double tol, DDalphaAMG_status *mg_status )
{
  set_n_flavours( 2 );
  DDalphaAMG_driver( vector1_out, vector1_in, vector2_out, vector2_in, tol, mg_status, _SOLVE_SQ_EVEN );
  set_n_flavours( 1 );
}

void DDalphaAMG_apply_operator( double *vector_out, double *vector_in, DDalphaAMG_status *mg_status ) {
  DDalphaAMG_driver( vector_out, vector_in, NULL, NULL, 0, mg_status, _OPERATOR );
}

void DDalphaAMG_apply_operator_doublet( double *vector1_out, double *vector1_in,
                                        double *vector2_out, double *vector2_in, DDalphaAMG_status *mg_status )
{
  set_n_flavours( 2 );
  DDalphaAMG_driver( vector1_out, vector1_in, vector2_out, vector2_in, 0, mg_status, _OPERATOR );
  set_n_flavours( 1 );
}

void DDalphaAMG_preconditioner( double *vector_out, double *vector_in, DDalphaAMG_status * mg_status ) {
  DDalphaAMG_driver( vector_out, vector_in, NULL, NULL, 0, mg_status, _PRECOND );
}

void DDalphaAMG_preconditioner_doublet( double *vector1_out, double *vector1_in,
                                        double *vector2_out, double *vector2_in, DDalphaAMG_status *mg_status )
{
  set_n_flavours( 2 );
  DDalphaAMG_driver( vector1_out, vector1_in, vector2_out, vector2_in, 0, mg_status, _PRECOND );
  set_n_flavours( 1 );
}

void DDalphaAMG_restrict( float *vector_out, float *vector_in, int level, DDalphaAMG_status *mg_status ) {
  DDalphaAMG_proj_driver( vector_out, vector_in, level, mg_status, _RESTRICT );
}

void DDalphaAMG_prolongate( float *vector_out, float *vector_in, int level, DDalphaAMG_status *mg_status ) {
  DDalphaAMG_proj_driver( vector_out, vector_in, level, mg_status, _PROLONGATE );
}

void DDalphaAMG_apply_coarse_operator( float *vector_out, float *vector_in, int level, DDalphaAMG_status *mg_status ) {
  DDalphaAMG_proj_driver( vector_out, vector_in, level, mg_status, _OPERATOR );
}

void DDalphaAMG_free( void ) {
  method_free( &l );
  g.setup_flag = 0;
}


void DDalphaAMG_finalize( void ) {

  if (g.setup_flag)
    method_free( &l );

  finalize_common_thread_data(commonthreaddata);
  finalize_no_threading(no_threading);
  for(int i=0; i<g.num_openmp_processes; i++) {
    FREE( threading[i], struct Thread, 1);
  }
  FREE( threading, struct Thread *, g.num_openmp_processes);
  FREE( commonthreaddata, struct common_thread_data, 1);

  method_finalize( &l );

  free( no_threading );
  
}

/***************  Extra Functions  ***********************************************************************************/

MPI_Comm DDalphaAMG_get_communicator( void ){
  MPI_Comm comm;
  MPI_Comm_dup( g.comm_cart, &comm);
  return comm;
}

void DDalphaAMG_read_configuration( double *gauge_field, char *filename, int format, DDalphaAMG_status *mg_status ) {

  double plaq;

  if(format==1)
    lime_read_conf( gauge_field, filename, &plaq );
  else
    read_conf( gauge_field, filename, &plaq, &l );

}

/* File I/O */
void DDalphaAMG_read_vector( double *vector_in, int nvec, char *filename, int format, DDalphaAMG_status *mg_status ){

  if(format==1)
    lime_read_vector( vector_in, nvec, filename );
  else
    vector_io( vector_in, NULL, nvec, _NVEC_OUTER, filename, _READ, "vector", &l );

}

void DDalphaAMG_write_vector( double *vector_out, int nvec, char *filename, int format, DDalphaAMG_status *mg_status ){

  if(format==1)
    lime_write_vector( vector_out, nvec, filename );
  else
    vector_io( vector_out, NULL, nvec, _NVEC_OUTER, filename, _WRITE, "vector", &l );

}

/* for top-level vectors */
void DDalphaAMG_define_vector_const( double *vector, int nvec, double re, double im ) {

  THREADED(threading[0]->n_core)
  if(vector!=NULL){
    int start, end;
    compute_core_start_end_custom( 0, l.inner_vector_size*nvec, &start, &end, &l, threading[omp_get_thread_num()], l.num_lattice_site_var );
    buffer_double_define( (buffer_double) vector, re+I*im, start, end, &l );
  }
  else {
    warning0("Vector NULL when calling DDalphaAMG_define_vector_const!");
  }
}

void DDalphaAMG_define_vector_rand( double *vector, int nvec ) {

  if(vector!=NULL){
    int i, j;
    THREADED(threading[0]->n_core) {
      int start, end;
      compute_core_start_end_custom( 0, l.inner_vector_size, &start, &end, &l, threading[omp_get_thread_num()], l.num_lattice_site_var);
      for ( i=start*nvec; i<end*nvec; i++ )
	vector[i] = ((double)rand()/(double)RAND_MAX)-0.5 + (((double)rand()/(double)RAND_MAX)-0.5)*_Complex_I;
    }
  }
  else {
    warning0("Vector NULL when calling DDalphaAMG_define_vector_const!");
  }

}

void DDalphaAMG_define_pt_src( double *vector, int nvec, int *global_pos, int *spin_color_ind ) {

  if(vector!=NULL){
    int i, x, y, z, t, n, desired_rank, *ll = l.local_lattice, size = l.inner_vector_size;
    for ( n=0; n<nvec; n++ ) {
      t = global_pos[4*n + T]; z = global_pos[4*n + Z]; y = global_pos[4*n + Y]; x = global_pos[4*n + X];
      desired_rank = process_index( t, z, y, x, ll );
      if ( g.my_rank == desired_rank ) {
	t = t%ll[T]; z = z%ll[Z]; y = y%ll[Y]; x = x%ll[X];
	if(vector_index_fct!=NULL )
	  i = vector_index_fct( t, z, y, x );
	else
	  i = 2*lex_index(t, z, y, x, ll );
	vector[size*n+i+spin_color_ind[n]] = 1.0;
      }
    }
  }
  else {
    warning0("Vector NULL when calling DDalphaAMG_define_pt_src!");
  }
  
}

/* wrapper functions for some of internal routines */
void DDalphaAMG_vector_norm( double* norm, double *vector, int nvec ) {

  if(vector!=NULL){
#ifdef HAVE_TM1p1 
    ASSERT( g.n_flavours == 1 ); // if not, we might need to make some assumptions on the order of elements in vector w/r/t flavour
#endif
    ASSERT( nvec%num_loop == 0 );
    
    vector_double vec;
    vec.vector_buffer = (buffer_double) vector;
    vec.num_vect = nvec;
    vec.num_vect_now = nvec;
    vec.size = l.inner_vector_size;
    vec.layout = _NVEC_OUTER;
    vector_double_change_layout( &vec, &vec, _NVEC_INNER, no_threading );
    
    THREADED(threading[0]->n_core) 
      global_norm_double( norm, &vec, 0, l.inner_vector_size, &l, threading[omp_get_thread_num()] );

    vector_double_change_layout( &vec, &vec, _NVEC_OUTER, no_threading );
   }
  else {
    warning0("Vector NULL when calling DDalphaAMG_define_vector_const!");
  }
}

void DDalphaAMG_vector_saxpy( double *vector_out, int nvec, double *a, double *x, double *y ) {

  if(vector_out!=NULL && x!=NULL && y!=NULL){
#ifdef HAVE_TM1p1
    ASSERT( g.n_flavours == 1 );// if not, we might need to make some assumptions on the order of elements in vector w/r/t flavour  
#endif
    ASSERT( nvec%num_loop == 0 );
    
    int size = l.inner_vector_size;
    vector_double vec_out, xx, yy;
    vec_out.vector_buffer= (buffer_double) vector_out; xx.vector_buffer= (buffer_double) x; yy.vector_buffer= (buffer_double) y;
    vec_out.num_vect = nvec; xx.num_vect = nvec; yy.num_vect = nvec;
    vec_out.num_vect_now = nvec; xx.num_vect_now = nvec; yy.num_vect_now = nvec;
    vec_out.size = size; xx.size = size; yy.size = size;
    vec_out.layout = _NVEC_OUTER; xx.layout = _NVEC_OUTER; yy.layout = _NVEC_OUTER;
    vec_out.l = &l; xx.l = &l; yy.l = &l;

    complex_double alpha[nvec];
    for( int i=0; i< nvec; i++ ) alpha[i] = a[i];

    vector_double_change_layout( &vec_out, &vec_out, _NVEC_INNER, no_threading );
    vector_double_change_layout( &xx, &xx, _NVEC_INNER, no_threading );
    vector_double_change_layout( &yy, &yy, _NVEC_INNER, no_threading );
    THREADED(threading[0]->n_core) {
      int start, end;
      compute_core_start_end_custom( 0, l.inner_vector_size, &start, &end, &l, threading[omp_get_thread_num()], l.num_lattice_site_var );
      vector_double_saxpy( &vec_out, &xx, &yy, alpha, 0, 1, start, end, &l );
    }
    vector_double_change_layout( &vec_out, &vec_out, _NVEC_OUTER, no_threading );
    vector_double_change_layout( &xx, &xx, _NVEC_OUTER, no_threading );
    vector_double_change_layout( &yy, &yy, _NVEC_OUTER, no_threading );
  }
  else {
    warning0("Vector NULL when calling DDalphaAMG_define_vector_const!");
  }

}

/* intended for coarse-level vectors */
float* DDalphaAMG_coarse_vector_alloc( int level, int nvec ) {
  level_struct *ltmp=&l;
  for(int i=0; i<level; i++)
    ltmp=ltmp->next_level;
  buffer_float vector = NULL;
  MALLOC(vector, complex_float, nvec*ltmp->vector_size );
  return (float*) vector;
}

void DDalphaAMG_coarse_vector_free( float *vector, int level, int nvec ) {
  level_struct *ltmp=&l;
  for(int i=0; i<level; i++)
    ltmp=ltmp->next_level;
  FREE(vector, float, 2*nvec*ltmp->vector_size );
}

void DDalphaAMG_coarse_vector_rand( float *vector, int level, int nvec ) {
  level_struct *ltmp=&l;
  int i;
  for(i=0; i<level; i++)
    ltmp=ltmp->next_level;

  THREADED(threading[0]->n_core)
  if(vector!=NULL){
    int start, end;
    compute_core_start_end( 0, ltmp->inner_vector_size, &start, &end, ltmp, threading[omp_get_thread_num()]);
    for ( int i=start*nvec; i<end*nvec; i++ )
      vector[i] = ((double)rand()/(double)RAND_MAX)-0.5 + (((double)rand()/(double)RAND_MAX)-0.5)*_Complex_I;
  }
  else {
    warning0("Vector NULL when calling DDalphaAMG_define_vector_const!");
  }
}

void DDalphaAMG_coarse_vector_residual( float *resid, float *vector_out, float *vector_in, int level, int nvec ) {

#ifdef HAVE_TM1p1
    ASSERT( g.n_flavours == 1 );
#endif
  if ( nvec%num_loop != 0 )
    error0("The number of rhs in the input and output vector needs to be multiple of num_loop\n");
  
  level_struct *ltmp=&l;
  for(int i=0; i<level; i++)
    ltmp=ltmp->next_level;

  int j, jj, vsize = ltmp->inner_vector_size;
  float norm1[nvec], norm2[nvec];
  vector_float vec[2];
  vec[0].vector_buffer = (complex_float *) vector_in;
  vec[1].vector_buffer = (complex_float *) vector_out;
  for(int i=0; i<2; i++) {
    vec[i].num_vect = nvec;
    vec[i].num_vect_now = nvec;
    vec[i].size = ltmp->inner_vector_size;
    vec[i].layout = _NVEC_OUTER;
    vec[i].l = ltmp;
  }
  VECTOR_LOOP(j, nvec, jj, norm1[j+jj] = 1.; norm2[j+jj] = 1.;)
    
  if(vector_in!=NULL){
    vector_float_change_layout( &(vec[0]), &(vec[0]), _NVEC_INNER, no_threading );
    THREADED(threading[0]->n_core)
      global_norm_float( norm1, &(vec[0]), 0, ltmp->inner_vector_size, ltmp, threading[omp_get_thread_num()] );
    if(vector_out!=NULL){
      vector_float_change_layout( &(vec[1]), &(vec[1]), _NVEC_INNER, no_threading );
      int start, end;
      compute_core_start_end( 0, ltmp->inner_vector_size, &start, &end, ltmp, threading[omp_get_thread_num()]);
      vector_float_minus( &(vec[1]), &(vec[1]), &(vec[0]), start, end, ltmp );
      vector_float_change_layout( &(vec[1]), &(vec[1]), _NVEC_OUTER, no_threading );
    }
    vector_float_change_layout( &(vec[0]), &(vec[0]), _NVEC_OUTER, no_threading );
  }
  else {
    warning0("VectorIn NULL when calling DDalphaAMG_coarse_vector_residual!\nProceed to compute the norm of VectorOut");
  }
  if(vector_out!=NULL){
    vector_float_change_layout( &(vec[1]), &(vec[1]), _NVEC_INNER, no_threading );
    THREADED(threading[0]->n_core)
      global_norm_float( norm2, &(vec[1]), 0, ltmp->inner_vector_size, ltmp, threading[omp_get_thread_num()] );
    vector_float_change_layout( &(vec[0]), &(vec[0]), _NVEC_OUTER, no_threading );
  }
  else {
    warning0("VectorOut NULL when calling DDalphaAMG_coarse_vector_residual!");
    VECTOR_LOOP(j, nvec, jj, resid[j+jj] = norm1[j+jj];)
    return;
  }
  
  VECTOR_LOOP(j, nvec, jj, resid[j+jj] = norm2[j+jj]/norm1[j+jj];)
}

void DDalphaAMG_test_routine( DDalphaAMG_status *mg_status ) {
  
  double t0, t1;
  t0 = MPI_Wtime();

  printf00("\n");
  THREADED(threading[0]->n_core)
  test_routine( &l, threading[omp_get_thread_num()]);

  if (g.test < 1e-5)
    mg_status->success = 1;
  else
    mg_status->success = 0;    
  mg_status->info = g.test; //highest error
  t1 = MPI_Wtime();
  mg_status->time = t1-t0;
  for ( int i=0; i<g.num_levels; i++ ) mg_status->iter_counts[i] = 0;
  for ( int i=0; i<g.num_levels; i++ ) mg_status->iter_times[i] = 0;
}
