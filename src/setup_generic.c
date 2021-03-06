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
 * checked: 12/09/2019
 * 1st cleanup:12/18/2019
 */

#include "main.h"

//shouldn't we gather many setup routines here?????? or put some in other files

static void inv_iter_2lvl_extension_setup_PRECISION( int setup_iter, level_struct *l, struct Thread *threading );//default
static void inv_iter_inv_fcycle_PRECISION( int setup_iter, level_struct *l, struct Thread *threading );
static void testvector_analysis_PRECISION( vector_PRECISION *test_vectors, level_struct *l, struct Thread *threading );
static void read_tv_from_file_PRECISION( level_struct *l, struct Thread *threading );


void interpolation_PRECISION_define( vector_double *V, level_struct *l, struct Thread *threading ) {
  /**********************************
   * The initial step of setting up the interpolation operator out of test vectors using only smoothing at the level = l->level  
   * Set up l->is_PRECISION.operator (interpolation_PRECISION_struct) in level_struct *l
   * The interpolation/restriction operators constructed out of test vectors on the level=l->level can be used to restrict vectors 
   *  on the l->level and interpolate vectors on the l->level->next_level
   * When V != NULL, test vectors are provided from a file, whose precision is assumed to be double
   **********************************/
  
  int k, i, j, jj, n = l->num_eig_vect, idof = l->num_lattice_site_var;
#ifdef DEBUG
  int pc =0, pi = 1, pn = n*6;
#endif
  PRECISION beta_r[n];
  complex_PRECISION beta_c[n]; 
  vector_PRECISION *buffer = NULL;
  int start = threading->start_index[l->depth];
  int end   = threading->end_index[l->depth];

  if ( V == NULL && g.interpolation < 3 ) { // if the initial test vectors are not in V nor supplied from a file later, compute them using smoothing
    if ( l->depth == 0) {
      tm_term_update( (complex_double)(g.setup_mu), l, threading ); // set (first num_loop) tm_term using a single shift for all rhs
      finalize_operator_update( l, threading ); // to update clover_oo_inv
    }
    
    PUBLIC_MALLOC( buffer, vector_PRECISION, 2 );
    for( i=0; i<2; i++){
      vector_PRECISION_init( &buffer[i] );
      vector_PRECISION_alloc( &buffer[i], (i==0)?_INNER:_ORDINARY, num_loop, l, threading );
    }
    START_MASTER(threading)
    if ( g.print > 0 ) printf0("initial definition --- depth: %d\n", l->depth );
#ifdef DEBUG
    if ( g.print > 0 ) { printf0("\033[0;42m\033[1;37m|"); fflush(0); }
#endif
    g.n_chunk = 0;
    END_MASTER(threading)

    for ( i=0; i<n; i+=num_loop ) {      
      START_LOCKED_MASTER(threading)
      // this would yield different results if we threaded it, so we don't (rand() is not thread-safe) 
      vector_PRECISION_define_random( &buffer[0], 0, l->inner_vector_size, l );
      END_LOCKED_MASTER(threading)

      for ( j = 0; j < g.post_smooth_iter[l->depth]-1; j++ ) {
	smoother_PRECISION( &buffer[1], NULL, &buffer[0], g.method >= 4?1:j+1, _NO_RES, l, threading );
	vector_PRECISION_copy( &buffer[0], &buffer[1], start, end, l );
      }	
      smoother_PRECISION( &buffer[1], NULL, &buffer[0], g.method >= 4?1:g.post_smooth_iter[l->depth], _NO_RES, l, threading );
      vector_PRECISION_copy2( &(l->is_PRECISION.test_vector_vec), &buffer[1], i, num_loop, -1, start, end, l );

#ifdef DEBUG
      pc += 6*num_loop;
      START_MASTER(threading)
      if ( pc >= 0.2*pi*pn ) { if ( g.print > 0 ) printf0("%4d%% |", 20*pi); if ( g.my_rank == 0 ) fflush(0); pi++; }
      END_MASTER(threading)
#endif
    }
    for( i=0; i<2; i++)
      vector_PRECISION_free( &buffer[i], l, threading );
    PUBLIC_FREE( buffer, vector_PRECISION, 2 );

    global_norm_PRECISION( beta_r, &(l->is_PRECISION.test_vector_vec), 0, l->inner_vector_size, l, threading );
    VECTOR_LOOP(j, n, jj, beta_c[j+jj] = (complex_PRECISION) beta_r[j+jj];)
    vector_PRECISION_real_scale( &(l->is_PRECISION.test_vector_vec), &(l->is_PRECISION.test_vector_vec), beta_c, 0, 1, start, end, l );
#ifdef DEBUG
    START_MASTER(threading)
    if ( g.print > 0 ) printf0("\033[0m\n");
    END_MASTER(threading)
#endif
    
  } else if ( g.interpolation < 3 ){
    // if the initial test vectors are provided, reorder them to make the interpolation op from them.
    // test vectors are read from a file later when g.interpolation == 3
    trans_PRECISION( &(l->is_PRECISION.test_vector_vec), V, l->s_PRECISION.op.translation_table, l, threading );
  }

  vector_PRECISION_copy( &(l->is_PRECISION.interpolation_vec), &(l->is_PRECISION.test_vector_vec), start, end, l );

  testvector_analysis_PRECISION( &(l->is_PRECISION.test_vector_vec), l, threading );

  gram_schmidt_on_aggregates_PRECISION( &(l->is_PRECISION.interpolation_vec), n, l, threading );
  define_interpolation_PRECISION_operator( &(l->is_PRECISION.interpolation_vec), l, threading );
}

void coarse_grid_correction_PRECISION_setup( level_struct *l, struct Thread *threading ) {

  if ( !l->idle ) {
    /**********************************************************************************
     * Operators are used in smoothing and coarse-grid correction as wel as at the top.
     * First construct coarse op in l (l->next_level->op_PRECISION).
     * Then, use this to define op in other struct's
     * Each struct contains op. struct but with different layout tables (Schwarz respects oddeven but op in g&l is of regular layout)
     *   l->next_level->s_PRECISION.op is a gathered version of l->next_level->op_PRECISION with Schwarz layout
     *   l->next_level->oe_op_PRECISION is a rearranged version of l->next_level->s_PRECISION.op
     * Op's in the below is just alias
     *   &(l->next_level->p_PRECISION.op) = &(l->next_level->s_PRECISION.op)
     *   This is true also for l->p_PRECISION.op at the top
     **********************************************************************************/
    
    //--- set up coarse op's (l->next_level->op_PRECISION.[clover/odd_proj/D]) using l->is_PRECISION.interpolation_vec in interpolation_PRECISION_struc
    START_LOCKED_MASTER(threading)
    coarse_operator_PRECISION_alloc( l );
    coarse_operator_PRECISION_setup( &(l->is_PRECISION.interpolation_vec), l );
    END_LOCKED_MASTER(threading)

    //--- set up schwarz_PRECISION_struct (next_level->s_PRECISION)
    START_LOCKED_MASTER(threading)
    if ( !l->next_level->idle ) {
      if ( l->next_level->level > 0 ) {
        schwarz_PRECISION_alloc( &(l->next_level->s_PRECISION), l->next_level );
        schwarz_layout_PRECISION_define( &(l->next_level->s_PRECISION), l->next_level );
      } else {// if not g.odd_even, this is used as op
        operator_PRECISION_alloc( &(l->next_level->s_PRECISION.op), _ORDINARY, l->next_level );
        operator_PRECISION_define( &(l->next_level->s_PRECISION.op), l->next_level );
        interpolation_PRECISION_alloc( l->next_level );
      }
    } else {
      interpolation_PRECISION_dummy_alloc( l->next_level );
    }
    // (l->next_level->s_PRECISION.op) <- (l->next_level->op_PRECISION)
    conf_PRECISION_gather( &(l->next_level->s_PRECISION.op), &(l->next_level->op_PRECISION), l->next_level );
    END_LOCKED_MASTER(threading)
      
    //--- set l->next_level->[s/p]_PRECISION.op; also l->next_level->oe_op_PRECISION if g.method=4,5 and g.odd_even
    if ( !l->next_level->idle && l->next_level->level > 0 ) {
      START_LOCKED_MASTER(threading)
      schwarz_PRECISION_boundary_update( &(l->next_level->s_PRECISION), l->next_level );
      END_LOCKED_MASTER(threading)
      if ( g.method >= 4 && g.odd_even ) {//???? g.method >= 4 => no AMG
	// alloc l->next_level->oe_op_PRECISION
        START_LOCKED_MASTER(threading)
	coarse_oddeven_alloc_PRECISION( l->next_level );
        END_LOCKED_MASTER(threading)
        coarse_oddeven_setup_PRECISION( &(l->next_level->s_PRECISION.op), _REORDER, l->next_level, threading );
      }
      // set l->next_level->p_PRECISION.op
      START_LOCKED_MASTER(threading)
      l->next_level->p_PRECISION.op = &(l->next_level->s_PRECISION.op);
      END_LOCKED_MASTER(threading)
    }
    
    //--- set oe_op_PRECISION if the next level is the bottom and in even-odd preconditioning
    if ( !l->next_level->idle && l->next_level->level == 0 && g.odd_even ) {
      START_LOCKED_MASTER(threading)
      coarse_oddeven_alloc_PRECISION( l->next_level );
      END_LOCKED_MASTER(threading)
      //set l->next_level->oe_op_PRECISION
      coarse_oddeven_setup_PRECISION( &(l->next_level->s_PRECISION.op), _NO_REORDERING, l->next_level, threading );
    }
  }
  
  // Call next_level_setup and then set up the interpolation op on the next level.
  // If l->next_level->level ==1, these have been done for the bottom.
  // Note: intp op constructed on the given level connects the given level and the level one below
  // recursive part: call itself
  if ( l->next_level->level > 0 ) {
    next_level_setup( NULL, l->next_level, threading );
    // set up interpolation_PRECISION_struct *is_PRECISION in l->next_level
    START_LOCKED_MASTER(threading)
    if ( !l->next_level->idle )
      interpolation_PRECISION_alloc( l->next_level );
    END_LOCKED_MASTER(threading)
    if ( !l->idle ) {
      START_LOCKED_MASTER(threading)
      //first fill the allocated memory with the random vectors
      if ( !l->next_level->idle )// we can optmize here!!!!!
        vector_PRECISION_define_random( &(l->next_level->is_PRECISION.test_vector_vec), 0,                                         
					    l->next_level->inner_vector_size, l->next_level );
      END_LOCKED_MASTER(threading)
      // then, replace those random vectors with restriction of available vectors on a finer grid
      vector_PRECISION vbuf, vbuf_next;
      vector_PRECISION_init(&vbuf);
      vector_PRECISION_init(&vbuf_next);
      vector_PRECISION_alloc( &vbuf, _INNER, num_loop, l, threading );
      vector_PRECISION_alloc( &vbuf_next, _INNER, num_loop, l, threading );
      vbuf.num_vect_now = num_loop;
      vbuf.num_vect_now = num_loop;
      for ( int i=0; i<MIN(l->next_level->num_eig_vect,l->num_eig_vect); i+=num_loop ) {
	vector_PRECISION_copy2( &(vbuf), &(l->is_PRECISION.test_vector_vec), i, num_loop, 1, 
				threading->start_index[l->depth], threading->end_index[l->depth], l );
	restrict_PRECISION( &(vbuf_next), &(vbuf), l, threading );
	vector_PRECISION_copy2( &(l->next_level->is_PRECISION.test_vector_vec), &(vbuf_next), i, num_loop, -1, 
				threading->start_index[l->depth+1], threading->end_index[l->depth+1], l );
      }
      vector_PRECISION_free( &vbuf, l, threading );
      vector_PRECISION_free( &vbuf_next, l, threading );
    }
    if ( !l->next_level->idle )
      interpolation_PRECISION_define( NULL, l->next_level, threading );
    
    coarse_grid_correction_PRECISION_setup( l->next_level, threading ); // recursive step
  }
}

void coarse_grid_correction_PRECISION_free( level_struct *l ) {
  
  next_level_free( l->next_level );

  if ( !l->idle ) {
    if ( !l->next_level->idle ) {
      if ( l->next_level->level > 0 ) {
        schwarz_PRECISION_free( &(l->next_level->s_PRECISION), l->next_level );
        if ( g.method >= 4 && g.odd_even ) {
          coarse_oddeven_free_PRECISION( l->next_level );
        }
      } else {
        operator_PRECISION_free( &(l->next_level->s_PRECISION.op), _ORDINARY, l->next_level );
        interpolation_PRECISION_free( l->next_level );
        if ( g.odd_even )
          coarse_oddeven_free_PRECISION( l->next_level );
      }
    } else {
      interpolation_PRECISION_dummy_free( l->next_level );
    }
    interpolation_PRECISION_free( l );
    coarse_operator_PRECISION_free( l );
  }  
}

/*********** ADAPTIVE/ITERATIVE SETUP ROUTINES  ****************************************************/

// called only in the iterative setup phase: method_iterative_setup
void iterative_PRECISION_setup( int setup_iter, level_struct *l, struct Thread *threading ) {

  if ( l->depth == 0 ) {
    switch ( g.interpolation ) {
      case 2: inv_iter_inv_fcycle_PRECISION( setup_iter, l, threading ); break;
      case 3: read_tv_from_file_PRECISION( l, threading ); break;
      default: inv_iter_2lvl_extension_setup_PRECISION( setup_iter, l, threading ); break;
    }
  }
  // test routine
  level_struct *lp = l;
  while( lp->level > 0 ) {
    testvector_analysis_PRECISION( &(lp->is_PRECISION.test_vector_vec), lp, threading );
    lp = lp->next_level;
    if ( lp == NULL )
      break;
  }
}

// recursive
void re_setup_PRECISION( level_struct *l, struct Thread *threading ) {
  
  if ( l->level > 0 ) {
    if ( !l->idle ) {
      vector_PRECISION_copy( &(l->is_PRECISION.interpolation_vec), &(l->is_PRECISION.test_vector_vec),
                               threading->start_index[l->depth], threading->end_index[l->depth], l );
      
      gram_schmidt_on_aggregates_PRECISION( &(l->is_PRECISION.interpolation_vec), l->num_eig_vect, l, threading );
      //if ( l->depth > -1 )
      gram_schmidt_on_aggregates_PRECISION( &(l->is_PRECISION.interpolation_vec), l->num_eig_vect, l, threading );
      define_interpolation_PRECISION_operator( &(l->is_PRECISION.interpolation_vec), l, threading );
      START_LOCKED_MASTER(threading)
      coarse_operator_PRECISION_setup( &(l->is_PRECISION.interpolation_vec), l );
      conf_PRECISION_gather( &(l->next_level->s_PRECISION.op), &(l->next_level->op_PRECISION), l->next_level );
      END_LOCKED_MASTER(threading)
      if ( !l->next_level->idle && l->next_level->level > 0 ) {
        START_LOCKED_MASTER(threading)
        schwarz_PRECISION_boundary_update( &(l->next_level->s_PRECISION), l->next_level );
        END_LOCKED_MASTER(threading)
        if ( g.method >= 4 && g.odd_even ) {
          coarse_oddeven_setup_PRECISION( &(l->next_level->s_PRECISION.op), _REORDER, l->next_level, threading );
        }
      }
      if ( !l->next_level->idle && l->next_level->level == 0 && g.odd_even ) {
        coarse_oddeven_setup_PRECISION( &(l->next_level->s_PRECISION.op), _NO_REORDERING, l->next_level, threading );
      }
#ifdef HAVE_FABULOUS
      if ( g.solver[l->depth+1] == _GCRO ) {
	/* P_l updated => D_(l+1) updated => its deflation space needs to be updated */
	if ( l->next_level->p_PRECISION.fab.U  != NULL ) free(l->next_level->p_PRECISION.fab.U);
	l->next_level->p_PRECISION.fab.U       = NULL;
	l->next_level->p_PRECISION.fab.ldu     = 0;
	l->next_level->p_PRECISION.fab.dsize   = 0;
	for ( int i=0; i<l->next_level->p_PRECISION.fab.k; i++ ) ((complex_PRECISION *) l->next_level->p_PRECISION.fab.eigvals)[i] = 0;
      }
#endif
      re_setup_PRECISION( l->next_level, threading );
    }
  }  
}

static void test_vector_PRECISION_update( int i, level_struct *l, struct Thread *threading ) {
  int j, jj, n = num_loop; // test vectors are processed in a chunk of num_loop
  PRECISION beta[n];
  complex_PRECISION beta_c[n];

  if ( l->level > 1 )
    test_vector_PRECISION_update( i, l->next_level, threading );
  
  if ( !l->idle ){
    if ( i < l->num_eig_vect ) {
      global_norm_PRECISION( beta,  &(l->p_PRECISION.x), 0, l->inner_vector_size, l, threading );
      VECTOR_LOOP(j, num_loop, jj, beta_c[j+jj] = beta[j+jj];)
      vector_PRECISION_real_scale( &(l->p_PRECISION.x), &(l->p_PRECISION.x), beta_c, 0, 1,
				   threading->start_index[l->depth], threading->end_index[l->depth], l );
      vector_PRECISION_copy2( &(l->is_PRECISION.test_vector_vec), &(l->p_PRECISION.x), i, num_loop, -1,
			      threading->start_index[l->depth], threading->end_index[l->depth], l );
    }
  }
}

static void set_kcycle_tol_PRECISION( PRECISION tol, level_struct *l ) {
  
  if ( !l->idle ) {
    l->p_PRECISION.tol = tol;
#ifdef HAVE_FABULOUS
    if ( g.solver[l->depth] && !g.use_only_fgrmes_at_setup && l->depth > 0 ) {
      PRECISION tolerance[1] = { tol };
      fabulous_PRECISION_struct *fab = &(l->p_PRECISION.fab);
      fabulous_set_parameters( fab->mvp, fab->max_iter, tolerance, 1, fab->handle );
    }
#endif
  }
  
  if ( l->level > 1 )
    set_kcycle_tol_PRECISION( tol, l->next_level );
}

static void inv_iter_2lvl_extension_setup_PRECISION( int setup_iter, level_struct *l, struct Thread *threading ) {
  /* Description: At each level, apply two-level AMG prec. on test vectors with post-smoothing and update
   *               coarse op at one-level below with updated test vectors
   */
  
  if ( !l->idle ) {
    vector_PRECISION buf1;
    gmres_PRECISION_struct gmres;
    int n_vect = l->num_eig_vect, j, jj, n = num_loop;
    PRECISION beta[n];
    complex_PRECISION beta_c[n];

    vector_PRECISION_init( &buf1 );
    vector_PRECISION_alloc( &buf1, _ORDINARY, num_loop, l, threading );
    
    START_LOCKED_MASTER(threading)
    fgmres_PRECISION_struct_init( &gmres );
    fgmres_PRECISION_struct_alloc( g.max_iter[g.num_levels-1], g.max_restart[g.num_levels-1], _ORDINARY, g.tol[g.num_levels-1], 
                                   _COARSE_SOLVER, _NOTHING, NULL, apply_coarse_operator_PRECISION, &gmres, l->next_level );
    if ( g.odd_even && l->next_level->level == 0 )
      gmres.v_end = l->next_level->oe_op_PRECISION.num_even_sites*l->next_level->num_lattice_site_var;
    END_LOCKED_MASTER(threading)
    
    for ( int k=0; k<setup_iter; k++ ) {
      int pc = 0;
#ifdef DEBUG
      int pi = 1, pn = l->num_eig_vect*l->post_smooth_iter;
#endif
      START_MASTER(threading)
      printf0("depth: %d, 2lvl correction step number %d...\n", l->depth, k+1 ); 
#ifdef DEBUG
      printf0("\033[0;42m\033[1;37m|"); fflush(0); 
#endif
      g.n_chunk = 0;
      END_MASTER(threading)

      //--- update test vectors using newly defined coarse op
      for ( int i=0; i<n_vect; i += num_loop ) {
	vector_PRECISION_copy2( &buf1, &(l->is_PRECISION.test_vector_vec), i, num_loop, 1, 
				threading->start_index[l->depth], threading->end_index[l->depth], l );
	restrict_PRECISION( &(gmres.b), &buf1, l, threading );
	if ( !l->next_level->idle ) {
	  if ( g.odd_even && l->next_level->level == 0 ) {
	    coarse_solve_odd_even_PRECISION( &gmres, &(l->next_level->oe_op_PRECISION), l->next_level, threading );
	  } else {
	    fgmres_PRECISION( &gmres, l->next_level, threading );
	  }
	}
	interpolate3_PRECISION( &buf1, &(gmres.x), l, threading );
	smoother_PRECISION( &buf1, NULL, &(l->is_PRECISION.test_vector_vec), l->post_smooth_iter, _RES, l, threading );
      
	global_norm_PRECISION( beta, &buf1, 0, l->inner_vector_size, l, threading );
	VECTOR_LOOP(j, n, jj, beta_c[j+jj] = beta[j+jj];)
	vector_PRECISION_real_scale( &buf1, &buf1, beta_c, 0, 1,
				     threading->start_index[l->depth], threading->end_index[l->depth], l );
	vector_PRECISION_copy2( &(l->is_PRECISION.test_vector_vec), &buf1, i, num_loop, -1,
				threading->start_index[l->depth], threading->end_index[l->depth], l );
      }
      pc += l->post_smooth_iter;
#ifdef DEBUG
      START_MASTER(threading)
      if ( pc >= 0.2*pi*pn ) { printf0("%4d%% |", 20*pi); fflush(0); pi++; }
      END_MASTER(threading)
#endif

#ifdef DEBUG
      START_MASTER(threading)
      printf0("\033[0m\n");
      END_MASTER(threading)
#endif
      //--- orthonormalize the test vectors aggregate-wise to prepare for defining interpolation op
      vector_PRECISION_copy( &(l->is_PRECISION.interpolation_vec), &(l->is_PRECISION.test_vector_vec),
			     threading->start_index[l->depth], threading->end_index[l->depth], l );
      gram_schmidt_on_aggregates_PRECISION( &(l->is_PRECISION.interpolation_vec), l->num_eig_vect, l, threading );
      //if ( l->depth > -1 )
      gram_schmidt_on_aggregates_PRECISION( &(l->is_PRECISION.interpolation_vec), l->num_eig_vect, l, threading );
      define_interpolation_PRECISION_operator( &(l->is_PRECISION.interpolation_vec), l, threading );

      //--- update the coarse op
      START_LOCKED_MASTER(threading)
      coarse_operator_PRECISION_setup( &(l->is_PRECISION.interpolation_vec), l );
      conf_PRECISION_gather( &(l->next_level->s_PRECISION.op), &(l->next_level->op_PRECISION), l->next_level );
      END_LOCKED_MASTER(threading)
      if ( !l->next_level->idle && l->next_level->level > 0 ) {
        START_LOCKED_MASTER(threading)
        schwarz_PRECISION_boundary_update( &(l->next_level->s_PRECISION), l->next_level );
        END_LOCKED_MASTER(threading)
        if ( g.method >= 4 && g.odd_even )
          coarse_oddeven_setup_PRECISION( &(l->next_level->s_PRECISION.op), _REORDER, l->next_level, threading );
      }
      if ( !l->next_level->idle && l->next_level->level == 0 && g.odd_even )
        coarse_oddeven_setup_PRECISION( &(l->next_level->s_PRECISION.op), _NO_REORDERING, l->next_level, threading );
    }//END:for ( int k=0; k<setup_iter; k++ 
    
    //--- dive one level down
    if ( l->level > 1 )
      inv_iter_2lvl_extension_setup_PRECISION( setup_iter, l->next_level, threading );

    vector_PRECISION_free( &buf1, l, threading );
    START_LOCKED_MASTER(threading)
    fgmres_PRECISION_struct_free( &gmres, l );
    END_LOCKED_MASTER(threading)
  }
}

static void inv_iter_inv_fcycle_PRECISION( int setup_iter, level_struct *l, struct Thread *threading ) {
  
  START_LOCKED_MASTER(threading)
  if ( l->depth == 0 )
    set_kcycle_tol_PRECISION( g.setup_tol, l );
  END_LOCKED_MASTER(threading)

  int nvec = l->num_eig_vect;
  
  if ( !l->idle ) {
    for ( int j=0; j<setup_iter; j++ ) {
      int pc = 0;
#ifdef DEBUG
      int pi = 1, pn = l->num_eig_vect*l->post_smooth_iter;
#endif
      
      START_LOCKED_MASTER(threading)
      if ( g.print > 0 ) printf0("depth: %d, bootstrap step number %d...\n", l->depth, j+1 );
#ifdef DEBUG
      if ( g.print > 0 ) { printf0("\033[0;42m\033[1;37m|"); if ( g.my_rank == 0 ) fflush(0); }
#endif
      g.n_chunk = 0;
      END_LOCKED_MASTER(threading)

      gram_schmidt_PRECISION( &(l->is_PRECISION.test_vector_vec), l->num_eig_vect, l, threading );
      
      int start = threading->start_index[l->depth];
      int end   = threading->end_index[l->depth];
      vector_PRECISION buf;// might be better to use buffer in l?????
      vector_PRECISION_init(&buf);
      vector_PRECISION_alloc( &buf, _INNER, num_loop, l, threading );
      buf.num_vect_now = num_loop;
      for ( int i=0; i<nvec; i += num_loop ) {
        vector_PRECISION_copy2( &buf, &(l->is_PRECISION.test_vector_vec), i, num_loop, 1, start, end, l );
	// apply K-cycle
	vcycle_PRECISION( &(l->p_PRECISION.x), NULL, &buf, _NO_RES, l, threading );
	test_vector_PRECISION_update( i, l, threading );
	
	pc += l->post_smooth_iter*num_loop;
#ifdef DEBUG
	START_MASTER(threading)
	if ( pc >= (int)((0.2*pi)*pn) ) { if ( g.print > 0 ) { printf0("%4d%% |", 20*pi); if ( g.my_rank == 0 ) fflush(0); } pi++; }
	END_MASTER(threading)
#endif
      }
      vector_PRECISION_free( &buf, l, threading );
#ifdef DEBUG
      START_MASTER(threading)
      if ( g.print > 0 ) printf0("\033[0m\n");
      END_MASTER(threading)
#endif      
      // redefine coarse Dirac op using new test vectors at the levels below l->level
      re_setup_PRECISION( l, threading );

      // go to the next level for each iteration if at the top 
      if ( l->depth == 0 && l->next_level->level > 0 ) {
        inv_iter_inv_fcycle_PRECISION( MAX(1,round( ((double)(j+1)*l->next_level->setup_iter)/
        ((double)setup_iter) )), l->next_level, threading );
      }
    }//END: for ( int j=0; j<setup_iter; j++ )
    if ( l->depth > 0 && l->next_level->level > 0 ) {
      inv_iter_inv_fcycle_PRECISION( MAX(1, round((double)(l->next_level->setup_iter*setup_iter)/
      ((double)l->setup_iter))), l->next_level, threading );
    }
  }
 
  if ( l->depth == 0 ) {
    START_LOCKED_MASTER(threading)
    set_kcycle_tol_PRECISION( g.kcycle_tol, l );
    END_LOCKED_MASTER(threading)
  }
}

// read test vectors from a file
// Note: test vectors should be in a single flavour
static  void read_tv_from_file_PRECISION( level_struct *l, struct Thread *threading ) {
  
  if ( l->depth == 0 ) {
    if ( g.tv_io_single_file ) {
      START_LOCKED_MASTER(threading)
	vector_io( NULL, NULL, l->num_eig_vect, _NVEC_INNER, g.tv_io_file_name, _READ, "test vectors", l );
      END_LOCKED_MASTER(threading)
      re_setup_PRECISION( l, threading );
      
    } else {
      START_LOCKED_MASTER(threading)

      int n = l->num_eig_vect, i;
      char filename[STRINGLENGTH+1];
      vector_double tmp;
      vector_double_init(&tmp);
      vector_double_alloc( &tmp, _INNER, n, l, no_threading );
      tmp.num_vect_now = n;
      tmp.layout = _NVEC_OUTER;
      
      for ( i=0; i<n; i++ ) {
        sprintf( filename, "%s.%02d", g.tv_io_file_name, i );
        printf0("%s.%02d\n", g.tv_io_file_name, i );
        vector_io( (double*)(tmp.vector_buffer+tmp.size*i), NULL, 1, _NVEC_OUTER, filename, _READ, "test vector", l );
      }
      vector_double_change_layout(&tmp, &tmp, _NVEC_INNER, no_threading );
      trans_PRECISION( &(l->is_PRECISION.test_vector_vec), &tmp, l->s_PRECISION.op.translation_table, l, no_threading );
      
      vector_double_free( &tmp, l, no_threading );

      END_LOCKED_MASTER(threading)

      re_setup_PRECISION( l, threading );
    }
  }
}

/*************  TEST ROUTINES ************************************************************************************/

void testvector_analysis_PRECISION( vector_PRECISION *test_vectors, level_struct *l, struct Thread *threading ) {
#ifdef TESTVECTOR_ANALYSIS
  START_UNTHREADED_FUNCTION(threading)
  if ( l->depth == 0 ) {
    int j, jj; 
    complex_PRECISION lambda1[l->num_eig_vect], lambda2[l->num_eig_vect], lambda[l->num_eig_vect];
    PRECISION mu1[l->num_eig_vect], mu2[l->num_eig_vect];
    vector_PRECISION vbuf[3];
    for ( int i=0; i<3; i ++ ) {
      vector_PRECISION_init( &(vbuf[i]) );
      vector_PRECISION_alloc( &(vbuf[i]), _ORDINARY, l->num_eig_vect, l, no_threading );
      vbuf[i].num_vect_now = l->num_eig_vect;
    }
    printf0("--------------------------------------- depth: %d ----------------------------------------\n", l->depth );
    apply_operator_PRECISION( &(vbuf[3]), test_vectors, &(l->p_PRECISION), l, no_threading ); // vbuf[3] = D*test_vectors
    coarse_gamma5_PRECISION( &(vbuf[0]), &(vbuf[3]), 0, l->inner_vector_size, l );            // vbuf[0] = g5*vbuf[3] = g5*D*test_vectors (<- D is g5-hermitian)
    global_inner_product_PRECISION( lambda1, test_vectors, &(vbuf[0]), 0, l->inner_vector_size, l, no_threading );   // lambda1 = test_vectors*g5*D*test_vectors
    global_inner_product_PRECISION( lambda2, test_vectors, test_vectors, 0, l->inner_vector_size, l, no_threading ); // lambda2 = test_vectors*test_vectors
    VECTOR_LOOP(j, l->num_eig_vect, jj, lambda[j+jj]=lambda1[j+jj]/lambda2[j+jj]; )                                  // lambda = labda1/lambda2
    vector_PRECISION_saxpy( &(vbuf[1]), &(vbuf[0]), test_vectors, lambda, 0, -1, 0, l->inner_vector_size, l ); // vbuf[1] = g5*D*test_vectors - lambda*test_vector
    global_norm_PRECISION( mu1, &(vbuf[1]), 0, l->inner_vector_size, l, no_threading );   // mu1 = ||g5*D*test_vectors - lambda*test_vector||
    global_norm_PRECISION( mu2, test_vectors, 0, l->inner_vector_size, l, no_threading ); // mu2 = ||test_vectors||
    for ( int i=0; i<l->num_eig_vect; i++ ) {
      printf0("vector #%02d: ", i+1 );
      printf0("singular value: %+lf%+lfi, singular vector precision: %le\n", (double)creal(lambda[i]), (double)cimag(lambda[i]), (double)(mu1[i]/mu2[i]) );
    }
    printf0("--------------------------------------- depth: %d ----------------------------------------\n", l->depth );
    for ( int i=0; i<3; i ++ ) vector_PRECISION_free( &(vbuf[i]), l, no_threading );
  }
  END_UNTHREADED_FUNCTION(threading)
#endif
}
