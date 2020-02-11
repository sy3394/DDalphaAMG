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
 *
 * Main Entries:
 *  method_init
 *  method_setup
 *  method_re_setup
 *  next_level_setup
 *  method_iterative_setup
 *
 * changed from sbacchio
 * checked: 12/09/2019
 * 1st cleanup:12/18/2019
 */

#include "main.h"
#include "DDalphaAMG.h"

//struct Thread *no_threading;//my proposal

complex_double _COMPLEX_double_ONE       = (complex_double) 1.0;
complex_double _COMPLEX_double_MINUS_ONE = (complex_double) (-1.0);
complex_double _COMPLEX_double_ZERO      = (complex_double) 0.0;
complex_float  _COMPLEX_float_ONE        = (complex_float)  1.0;
complex_float  _COMPLEX_float_MINUS_ONE  = (complex_float)  (-1.0);
complex_float  _COMPLEX_float_ZERO       = (complex_float)  0.0;


// initialize global structure
void g_init(){

  var_table_init( &(g.vt) );
  operator_double_init( &(g.op_double) );
  operator_float_init( &(g.op_float) );
  fgmres_double_struct_init( &(g.p) );
  fgmres_MP_struct_init( &(g.p_MP) );
  g.global_lattice = NULL;
  g.local_lattice = NULL;
  g.block_lattice = NULL;
  g.post_smooth_iter = NULL;
  g.block_iter = NULL;
  g.setup_iter = NULL;
  g.relax_fac = NULL;
#ifdef HAVE_TM
  g.mu_factor = NULL;
#endif
#ifdef HAVE_TM1p1
  g.epsbar_factor = NULL;
  g.n_flavours = 1;
#endif
  g.gamma = NULL;
  g.odd_even_table = NULL;
  g.cur_storage = 0;
  g.max_storage = 0;
  g.in_setup = 0;
  g.num_rhs_vect = 0;
  g.num_vect_now=0;
}

// initialize level structure
void l_init( level_struct *l ) {

  level_double_init( l );
  level_float_init( l );

  vector_double_init(&(l->x));
  l->next_level = NULL;
  l->reqs = NULL;
}

void method_init( int *argc, char ***argv, level_struct *l ) {
  
  /********************************************************************************* 
   * Sets up the global and level struct for the method and assignes values 
   * according to the inputfile. 
   * Used outside of the OpenMP parallel region
   * - int *argc: Argument count of main function. Determines if inputfile is 
   *   provided.
   * - char ***argv: In case inputfile is provided, contains name of this file.
   * 
   * CAUTION: changes in this function have no influence on the interface         
   * library, cf. DDalphaAMG_init.                                                    
   *********************************************************************************/

  //-------- open files for input and output ifdef WRITE_LOGFILE
  char inputfile[STRINGLENGTH];
  if ( *argc > 1 ) {
    strcpy( inputfile, (*argv)[1] );
  } else {
    strcpy( inputfile, "sample.ini" );
  }
#ifdef WRITE_LOGFILE
  g.logfile = fopen( "output.log", "w" );
  fprintf(g.logfile,"---------- log file -----------\n\n");
  fflush(g.logfile);
#endif
  
  //------- initialize and set global and level structure according to the input file
  g_init();
  l_init( l );
  lg_in( inputfile, l );
  data_layout_init( l );

  g.Cart_rank   = MPI_Cart_rank;
  g.Cart_coords = MPI_Cart_coords;
  cart_define( MPI_COMM_WORLD, l );

  //no_threading = (struct Thread *)malloc(sizeof(struct Thread));//my porosal
  //setup_no_threading(no_threading, l);//my proposal

  // The following fields are used in dirac_setup
  operator_double_alloc( &(g.op_double), _ORDINARY, l ); 
  operator_double_define( &(g.op_double), l );
  MALLOC( g.odd_even_table, int, l->num_inner_lattice_sites );
  define_odd_even_table( l );
}

void method_setup( vector_double *V, level_struct *l, struct Thread *threading ) {
  
  int nvec = num_loop;
  double t0=0, t1=0;

  ASSERT(l->depth == 0);
  
  START_LOCKED_MASTER(threading)
  g.in_setup = 1;
  if ( g.vt.evaluation ) {//?????
    l->level = g.num_levels-1;
  }
  
  if ( l->depth==0 )
    prof_init( l );
  
  //-------------- Setup g.p(_MP) for the solver part (gmres_PRECISION_struct in global_struct: used only in the solver part)
  //                  if g.mixed_precision == 2, prec = vcycle_float; otherwise, prec = preconditioner
#ifdef HAVE_TM1p1
  nvec *= 2;
#endif
  if ( g.method > 0 ) {//------------- FGMRES + alpha
#ifdef INIT_ONE_PREC
    if ( g.mixed_precision == 2 ) {// INIT_ONE_PREC is defiend in main.h
#endif
      fgmres_MP_struct_alloc( g.restart, g.max_restart, _INNER,
                              g.tol, _RIGHT, vcycle_float_new, &(g.p_MP), l );
      g.p.op            = &(g.op_double);
      g.p.eval_operator = d_plus_clover_double_new;
#if defined(INIT_ONE_PREC) && (defined (DEBUG) || defined (TEST_VECTOR_ANALYSIS))
      vector_double_alloc( &(g.p.b), _INNER, nvec, l, no_threading );
      vector_double_alloc( &(g.p.x), _INNER, nvec, l, no_threading );
#endif
#ifdef INIT_ONE_PREC
    } else {
#endif
      fgmres_double_struct_alloc( g.restart, g.max_restart, _INNER, g.tol,
                                  _GLOBAL_FGMRES, _RIGHT, preconditioner_new,
				  d_plus_clover_double_new, &(g.p), l );
    }
#ifdef INIT_ONE_PREC
  }
#endif
  else if ( g.method == 0 ) {//------ pure GMRES (no AMG)
#ifdef INIT_ONE_PREC
    if ( g.mixed_precision == 2 ) {
#endif
      fgmres_MP_struct_alloc( g.restart, g.max_restart, _INNER,
                              g.tol, _NOTHING, NULL, &(g.p_MP), l );
      g.p.op = &(g.op_double);//moved from fgmres_MP_struct_alloc: why only double? double at the top; single below????
      g.p.eval_operator = d_plus_clover_double_new;//moved from fgmres_MP_struct_alloc: why only double????
#if defined(INIT_ONE_PREC) && (defined (DEBUG) || defined (TEST_VECTOR_ANALYSIS))
      vector_double_alloc( &(g.p.b), _INNER, nvec, l, no_threading );
      vector_double_alloc( &(g.p.x), _INNER, nvec, l, no_threading );
#endif
#ifdef INIT_ONE_PREC
    } else {
#endif
      fgmres_double_struct_alloc( g.restart, g.max_restart, _INNER, g.tol,
                                  _GLOBAL_FGMRES, _NOTHING, NULL, d_plus_clover_double_new,
                                  &(g.p), l );
#ifdef INIT_ONE_PREC
    }
#endif
  } 
  else if ( g.method == -1 ) {//----- pure CGN (no AMG) 
    fgmres_double_struct_alloc( 4, g.restart*g.max_restart, _INNER, g.tol,
                                _GLOBAL_FGMRES, _NOTHING, NULL, d_plus_clover_double_new, &(g.p), l );
    fine_level_double_alloc( l );
  }
  END_LOCKED_MASTER(threading)
  SYNC_MASTER_TO_ALL(threading)

  //------------------ Set the level structure recursively
  //                     l->p_PRECISION at the top level is not allocated as it is not used  
  if ( g.method >= 0 ) {
    START_LOCKED_MASTER(threading)
    t0 = MPI_Wtime();
    if ( g.mixed_precision ) {
      smoother_float_def( l ); // define s_PRECISION, and p_PRECISION(gmres:k-cycle) if g.method=4 or 5 why here????
      if ( g.method >= 4 && g.odd_even )
        oddeven_setup_float( &(g.op_double), l );
    } else {
      smoother_double_def( l );
      if ( g.method >= 4 && g.odd_even )
        oddeven_setup_double( &(g.op_double), l );
    }
    END_LOCKED_MASTER(threading)
    SYNC_MASTER_TO_ALL(threading)
    if ( g.method > 0 )
      if ( g.interpolation && g.num_levels > 1 ) {
	// allocate memory for interpolation and define interpolation op at depth == 0
	if ( g.mixed_precision ){
	  START_LOCKED_MASTER(threading)
	  fine_level_float_alloc( l ); //p_PRECISION.b/x are allocated here.
	  // The initial step of setting up the interpolation operator out of test vectors using only smoothing 
	  interpolation_float_alloc( l );
	  END_LOCKED_MASTER(threading)
	  interpolation_float_define_new( V, l, threading );
	} else {
	  START_LOCKED_MASTER(threading)
	  fine_level_double_alloc( l ); //p_PRECISION.b/x are allocated here.
	  // The initial step of setting up the interpolation operator out of test vectors using only smoothing
	  interpolation_double_alloc( l );
	  END_LOCKED_MASTER(threading)
	  interpolation_double_define_new( V, l, threading );  //!!!!!
	}
        next_level_setup_new( V, l, threading );//recursive function: dictated by the chice of g.num_vect_now!!!!!!
      }
    START_LOCKED_MASTER(threading)
    t1 = MPI_Wtime();
    g.total_time = t1-t0;
    printf0("elapsed time: %lf seconds\n", t1-t0 );
    END_LOCKED_MASTER(threading)
  }

  //------------------- Report the setup result
  START_LOCKED_MASTER(threading)
#ifdef PARAMOUTPUT  
  if ( g.method >= -1 && g.print > 0 && !( g.vt.evaluation && g.vt.re_setup ) ) {
    if(threading->n_core>1) {
      printf0("\nrunning with %d openmp threads per core", threading->n_core );
    }
    printf0("\n+----------------------------------------------------------+\n");
    if ( g.method > 0 ) {
      printf0("| %d-level method                                           |\n", l->level+1);
      if ( l->level > 0 )
        printf0("| postsmoothing %s-cycle                                    |\n", g.kcycle?"K":(l->n_cy>1?"W":"V") );
    }
    switch( g.method ) {
      case -1: printf0("| pure CGN                                                 |\n"); break;
      case  0: printf0("| pure GMRES                                               |\n"); break;
      case  1: printf0("| FGMRES + additive Schwarz                                |\n"); break;
      case  2: printf0("| FGMRES + red-black multiplicative Schwarz                |\n"); break;
      case  3: printf0("| FGMRES + sixteen color multiplicative Schwarz            |\n"); break;
      default: printf0("| FGMRES + GMRES                                           |\n"); break;
    }
    if ( g.method >=0  )
      printf0("|          restart length: %-3d                             |\n", g.restart );
    printf0("|                      m0: %+9.6lf                       |\n", g.m0 );
    if(g.setup_m0!=g.m0)
      printf0("|                setup m0: %+9.6lf                       |\n", g.setup_m0 );
    printf0("|                     csw: %+9.6lf                       |\n", g.csw );
#ifdef HAVE_TM
    printf0("|                      mu: %+9.6lf                       |\n", g.mu);
    if(g.setup_mu!=g.mu)
      printf0("|                setup mu: %+9.6lf                       |\n", g.setup_mu );
    if(g.mu_odd_shift!=0.)
      printf0("|         mu on odd sites: %+9.6lf                       |\n", g.mu + g.mu_odd_shift );
    if(g.mu_even_shift!=0.)
      printf0("|        mu on even sites: %+9.6lf                       |\n", g.mu + g.mu_even_shift );
#endif
#ifdef HAVE_TM1p1
    if(g.epsbar)
      printf0("|                  epsbar: %+9.6lf                       |\n", g.epsbar);
    if(g.epsbar_ig5_odd_shift!=0.)
      printf0("|    ig5 epsbar odd sites: %+9.6lf                       |\n", g.epsbar_ig5_odd_shift );
    if(g.epsbar_ig5_even_shift!=0.)
      printf0("|   ig5 epsbar even sites: %+9.6lf                       |\n", g.epsbar_ig5_even_shift );
#endif
    if ( g.method > 0 ) {
      printf0("+----------------------------------------------------------+\n");
      printf0("|%17s cycles: %-6d                          |\n", "preconditioner", l->n_cy );
      printf0("|            inner solver: %-26s      |\n", g.method==4?"GMRES":"minimal residual iteration" );
      printf0("|               precision: %6s                          |\n", g.mixed_precision?"single":"double" );
    }
    for ( int i=0; i<g.num_levels; i++ ) {
      int *gl = g.global_lattice[i], *ll = g.local_lattice[i], *bl = g.block_lattice[i];
      printf0("+---------------------- depth %2d --------------------------+\n", i );
      printf0("|          global lattice: %-3d %-3d %-3d %-3d                 |\n", gl[0], gl[1], gl[2], gl[3] );
      printf0("|           local lattice: %-3d %-3d %-3d %-3d                 |\n", ll[0], ll[1], ll[2], ll[3] );
      if ( g.method > 0 ) {
        printf0("|           block lattice: %-3d %-3d %-3d %-3d                 |\n", bl[0], bl[1], bl[2], bl[3] );
        if ( i+1 < g.num_levels ) {
          printf0("|        post smooth iter: %-3d                             |\n", g.post_smooth_iter[i] );
          printf0("|     smoother inner iter: %-3d                             |\n", g.block_iter[i] );
          printf0("|              setup iter: %-3d                             |\n", g.setup_iter[i] );
          printf0("|            test vectors: %-3d                             |\n", g.num_eig_vect[i] );
        } else {
          printf0("|      coarge grid solver: %-30s  |\n", g.odd_even?"odd even GMRES":"GMRES" );
          printf0("|              iterations: %-6d                          |\n", g.coarse_iter );
          printf0("|                  cycles: %-6d                          |\n", g.coarse_restart );
          printf0("|               tolerance: %-5.0le                           |\n", g.coarse_tol );
        }
#ifdef HAVE_TM
        if( g.mu!=0. && g.mu_factor[i]!=1 )
          printf0("|                      mu: %+9.6lf                       |\n", g.mu * g.mu_factor[i] );
        if( g.mu_odd_shift!=0. && g.mu_factor[i]!=1 )
          printf0("|         mu on odd sites: %+9.6lf                       |\n", (g.mu + g.mu_odd_shift) * g.mu_factor[i] );
        if( g.mu_even_shift!=0. && g.mu_factor[i]!=1 )
          printf0("|        mu on even sites: %+9.6lf                       |\n", (g.mu + g.mu_even_shift) * g.mu_factor[i] );
#endif
#ifdef HAVE_TM1p1
        if( g.epsbar!=0. && g.epsbar_factor[i]!=1 )
          printf0("|                 epsbar: %+9.6lf                       |\n", g.epsbar * g.epsbar_factor[i] );
        if(g.epsbar_ig5_odd_shift!=0. && g.epsbar_factor[i]!=1)
          printf0("|  ig5 epsbar on odd sites: %+9.6lf                     |\n", (g.epsbar + g.epsbar_ig5_odd_shift) * g.epsbar_factor[i] );
        if(g.epsbar_ig5_even_shift!=0. && g.epsbar_factor[i]!=1)
          printf0("| ig5 epsbar on even sites: %+9.6lf                     |\n", (g.epsbar + g.epsbar_ig5_even_shift) *g.epsbar_factor[i] );
#endif
      }
    }
    if ( g.method > 0 && g.kcycle > 0 ) {
      printf0("+----------------------------------------------------------+\n");
      printf0("|          K-cycle length: %-6d                          |\n", g.kcycle_restart );
      printf0("|        K-cycle restarts: %-6d                          |\n", g.kcycle_max_restart );
      printf0("|       K-cycle tolerance: %-5.0le                           |\n", g.kcycle_tol );
    }
    printf0("+----------------------------------------------------------+\n");
    printf0("\n");
  }
#endif
  g.in_setup = 0;

  if ( l->depth==0 && g.method >=0 )
    prof_print( l );
  END_LOCKED_MASTER(threading)
  
#ifdef DEBUG
  test_routine_new( l, threading );
#endif
}

void method_re_setup( level_struct *l, struct Thread *threading ) {
  START_LOCKED_MASTER(threading)
  method_free( l );
  END_LOCKED_MASTER(threading)
  SYNC_MASTER_TO_ALL(threading)
  method_setup( NULL, l, threading );
}

// initialize top-level if depth = 0 and next_level in level_struct *l
void next_level_setup_new( vector_double *V, level_struct *l, struct Thread *threading ) {

  if ( l->level > 0 ) {//if not at the bottom: safety net
    int mu;

    START_LOCKED_MASTER(threading)
    // allocate storage for next level parameters and initialize them
    MALLOC( l->next_level, level_struct, 1 );
    l_init( l->next_level );
    
    // define next level parameters
    l->next_level->level                = l->level-1;
    l->next_level->depth                = l->depth+1;
    l->next_level->tol                  = l->tol;
    l->next_level->post_smooth_iter     = g.post_smooth_iter[l->depth+1];
    l->next_level->relax_fac            = g.relax_fac[l->depth+1];
    l->next_level->block_iter           = g.block_iter[l->depth+1];
    l->next_level->setup_iter           = g.setup_iter[l->depth+1];
    l->next_level->num_eig_vect         = l->level==1?l->num_eig_vect:g.num_eig_vect[l->depth+1];//???????
    l->next_level->num_parent_eig_vect  = l->num_eig_vect;
    l->next_level->num_lattice_site_var = 2 * l->num_eig_vect;
    l->next_level->n_cy                 = g.ncycle[l->depth+1];
    l->next_level->global_lattice       = g.global_lattice[l->depth+1];
    l->next_level->local_lattice        = g.local_lattice[l->depth+1];
    l->next_level->block_lattice        = g.block_lattice[l->depth+1];
    l->next_level->num_processes        = 1;
    for (mu=0; mu<4; mu++) { // set dims of aggregates at the next level
      if ( l->depth+2 < g.num_levels ) // if l->next_level->level > 0 , i.e., if the next level is not the bottom
        l->next_level->coarsening[mu] = g.global_lattice[l->depth+1][mu]/g.global_lattice[l->depth+2][mu];
      else                             // if the next level is the bottom, dims of aggregates at the next level = dims of local lattice
        l->next_level->coarsening[mu] = g.local_lattice[l->depth+1][mu];
      
      l->next_level->num_processes_dir[mu] = l->next_level->global_lattice[mu]/l->next_level->local_lattice[mu];
      l->next_level->comm_offset[mu]       = (l->num_processes_dir[mu]/l->next_level->num_processes_dir[mu])*l->comm_offset[mu];
      l->next_level->global_splitting[mu]  = l->next_level->global_lattice[mu] / l->next_level->local_lattice[mu];
      l->next_level->periodic_bc[mu]       = l->periodic_bc[mu];
      l->next_level->num_processes        *= l->next_level->num_processes_dir[mu];
    }    

    data_layout_init( l->next_level ); // distribution of lattice sites over ranks
    neighbor_define( l->next_level );  // rank neighbors

    // update threading struct for no threading with size info of the next level
    update_threading(no_threading, l);

    END_LOCKED_MASTER(threading)
    SYNC_MASTER_TO_ALL(threading)

    // update threading struct with size info of the next level 
    update_threading(threading, l); // each thread has its own threading struc

    //------ entrance to recursion: coarse_grid_correction_PRECISION_setup = recursive function (next_level_setup_new is called inside of this)
    if ( g.mixed_precision ) {
      START_LOCKED_MASTER(threading)
      level_double_init( l->next_level );
      next_level_float_setup( l );// gmres_float_struct p_float (k-cycle) is defined here.
      END_LOCKED_MASTER(threading)
      if ( l->depth == 0 )
        coarse_grid_correction_float_setup_new( l, threading );//recursive function
    } else {
      START_LOCKED_MASTER(threading)
      level_double_init( l->next_level );
      next_level_double_setup( l );// gmres_double_struct p_double (k-cycle) is defined here. 
      END_LOCKED_MASTER(threading)
      if ( l->depth == 0 )
        coarse_grid_correction_double_setup_new( l, threading );//recursive function
    }
  }
  
  START_LOCKED_MASTER(threading)
  if ( l->depth == 0 ) printf0("\ninitial coarse grid correction is defined\n");
  END_LOCKED_MASTER(threading)
}

void method_iterative_setup( int setup_iter, level_struct *l, struct Thread *threading ) {
  
  if ( g.method > 0 && g.interpolation && g.num_levels > 1 && setup_iter > 0 ) {
    
    double t0=0, t1=0;
    START_LOCKED_MASTER(threading)
    g.in_setup = 1;
    if ( l->depth==0 )//is this here to reset profiler????
      prof_init( l );
    END_LOCKED_MASTER(threading)

    MASTER(threading)
      t0 = MPI_Wtime();
    
    if ( g.setup_m0 != g.m0 ) {
      m0_update( (complex_double)g.setup_m0, l, threading );
#ifdef HAVE_TM
    }
    if ( g.setup_mu != g.mu ) {
      tm_term_update( (complex_double)g.setup_mu, l, threading );
      finalize_operator_update( l, threading );
    } else if (g.setup_m0 != g.m0) {
#endif
      finalize_operator_update( l, threading );//test routine
    }
    
    if ( g.mixed_precision )
      iterative_float_setup_new( setup_iter, l, threading );
    else
      iterative_double_setup_new( setup_iter, l, threading );

    
    if ( g.setup_m0 != g.m0 ) {
      m0_update( (complex_double)g.m0, l, threading );
#ifdef HAVE_TM
    }
    if ( g.setup_mu != g.mu ) {
      tm_term_update( (complex_double)g.mu, l, threading );
      finalize_operator_update( l, threading );
    } else if (g.setup_m0 != g.m0) {
#endif
      finalize_operator_update( l, threading );
    }

    MASTER(threading) {
      t1 = MPI_Wtime();
      g.total_time = t1-t0;
      printf0("\nperformed %d iterative setup steps\n", setup_iter );
      printf0("elapsed time: %lf seconds (%lf seconds on coarse grid)\n\n", t1-t0, g.coarse_time );
    }
    
    START_LOCKED_MASTER(threading)
    g.in_setup = 0;
    if ( l->depth==0 )
      prof_print( l );
    END_LOCKED_MASTER(threading)
      
#ifdef DEBUG
    test_routine_new( l, threading );
#endif
  }
}

// next_level_free -> next_level_PRECISION_free -> coarse_grid_correction_PRECISION_free -> next_level_free
void next_level_free( level_struct *l ) {

  if ( l->level > 0 ) { 
    if ( g.mixed_precision )
      next_level_float_free( l );
    else
      next_level_double_free( l );
    FREE( l->next_level, level_struct, 1 );
  }
}

void method_free( level_struct *l ) {
  
  if ( g.method>=0 ) {
    if ( g.mixed_precision ) {
      if ( g.method >= 4 && g.odd_even )
        oddeven_free_float( l );
      smoother_float_free( l );
    } else {
      if ( g.method >= 4 && g.odd_even )
        oddeven_free_double( l );
      smoother_double_free( l );
    }
    if ( g.method > 0 )
      if ( g.interpolation ) {
	if ( g.mixed_precision )
	  fine_level_float_free( l );
	else
	  fine_level_double_free( l );
        next_level_free( l );
      }
  } else if ( g.method == -1 ) {
    fine_level_double_free( l );
  }

#ifdef INIT_ONE_PREC
  if ( g.mixed_precision == 2 && g.method >= 0 ) {
#endif
    fgmres_MP_struct_free( &(g.p_MP), l );
#if defined (INIT_ONE_PREC) && (defined (DEBUG) || defined (TEST_VECTOR_ANALYSIS))
    vector_double_free( &(g.p.b), l, no_threading );
    vector_double_free( &(g.p.x), l, no_threading );
#endif
#ifdef INIT_ONE_PREC
  } else {
#endif
    fgmres_double_struct_free( &(g.p), l );
#ifdef INIT_ONE_PREC
  }
#endif

}

void method_finalize( level_struct *l ) {
  
  int ls = MAX(g.num_desired_levels,2);
  
  operator_double_free( &(g.op_double), _ORDINARY, l );
  FREE( g.odd_even_table, int, l->num_inner_lattice_sites );
  FREE( g.global_lattice[0], int, 4*ls );
  FREE( g.local_lattice[0], int, 4*ls );
  FREE( g.block_lattice[0], int, 4*ls );
  FREE( g.global_lattice, int*, ls );
  FREE( g.local_lattice, int*, ls );
  FREE( g.block_lattice, int*, ls );
  FREE( g.post_smooth_iter, int, ls );
  FREE( g.ncycle, int, ls );
  FREE( g.relax_fac, double, ls );
#ifdef HAVE_TM
  FREE( g.mu_factor, double, ls );
#endif
#ifdef HAVE_TM1p1
  FREE( g.epsbar_factor, double, ls );
#endif
  FREE( g.block_iter, int, ls );
  FREE( g.setup_iter, int, ls );
  FREE( g.num_eig_vect, int, ls );
  cart_free( l );
  var_table_free( &(g.vt) );
  
  if ( g.cur_storage )
    warning0("amount of not freed memory: %lf MB\n", g.cur_storage );
  
#ifdef WRITE_LOGFILE
  fprintf(g.logfile,"---------- end of log file -----------\n\n");
  fflush(g.logfile);
  fclose(g.logfile);
#endif
}
