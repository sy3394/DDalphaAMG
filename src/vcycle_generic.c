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
 */

#include "main.h"
#include "vcycle_PRECISION.h"


// phi <- V/K-cycle(eta); compute corse-grid correction
void vcycle_PRECISION( vector_PRECISION *phi, vector_PRECISION *Dphi, vector_PRECISION *eta,
		       int res, level_struct *l, struct Thread *threading ) {

  if ( g.interpolation && l->level>0 ) {
    int nvec = num_loop;//eta->num_vect_now;
    for ( int i=0; i<l->n_cy; i++ ) {
      //--- compute the residual //presmoohting????
      if ( i==0 && res == _NO_RES ) {
	// if the initial guess is 0, the resid on the next level is restircted rhs
	l->next_level->p_PRECISION.b.num_vect_now = nvec;
        restrict_PRECISION( &(l->next_level->p_PRECISION.b), eta, l, threading );
      } else {
	// otherwise, compute b_(l+1) = R(b_l - D_l eta_l)
        int start = threading->start_index[l->depth];
        int end   = threading->end_index[l->depth];
	l->vbuf_PRECISION[0].num_vect_now = nvec; l->vbuf_PRECISION[1].num_vect_now = nvec; l->next_level->p_PRECISION.b.num_vect_now = nvec;
        apply_operator_PRECISION( &(l->vbuf_PRECISION[0]), phi, &(l->p_PRECISION), l, threading );
        vector_PRECISION_minus( &(l->vbuf_PRECISION[1]), eta, &(l->vbuf_PRECISION[0]), start, end, l );
        restrict_PRECISION( &(l->next_level->p_PRECISION.b), &(l->vbuf_PRECISION[1]), l, threading );
      }

      //--- recursive inversion part: base case = solve
      if ( !l->next_level->idle ) {
        START_MASTER(threading)
	if ( l->depth == 0 )
	  g.coarse_time -= MPI_Wtime();
        END_MASTER(threading) 
	if ( l->level > 1 ) {
	  // If the next level is not at the bottom
	  if ( g.kcycle )//default
	    solver_PRECISION( &(l->next_level->p_PRECISION), l->next_level, threading );//in this function, call vcycle as p_PRECISION use it as prec
	  else {
	    // dive one level down 
	    l->next_level->p_PRECISION.x.num_vect_now = nvec;
	    vcycle_PRECISION( &(l->next_level->p_PRECISION.x), NULL, &(l->next_level->p_PRECISION.b), _NO_RES, l->next_level, threading );
	  }
	} else { 
	  // if the next level is the bottom (coarse_apply_schur_complement_PRECISION is used only at the bottom if g.odd_even),
	    if ( g.odd_even ) {
	      coarse_solve_odd_even_PRECISION( &(l->next_level->p_PRECISION), &(l->next_level->oe_op_PRECISION), l->next_level, threading );
	    } else {
	      solver_PRECISION( &(l->next_level->p_PRECISION), l->next_level, threading );
	    }
	}
        START_MASTER(threading)
	if ( l->depth == 0 )
	  g.coarse_time += MPI_Wtime();
        END_MASTER(threading)
      }
      //--- interpolate the correction      
      l->next_level->p_PRECISION.x.num_vect_now = num_loop;
      if( i == 0 && res == _NO_RES )
        interpolate3_PRECISION( phi, &(l->next_level->p_PRECISION.x), l, threading );
      else
        interpolate_PRECISION( phi, &(l->next_level->p_PRECISION.x), l, threading );
      //--- Perform post smoothing: phi <- Smooth(phi;eta)
      smoother_PRECISION( phi, Dphi, eta, l->post_smooth_iter, _RES, l, threading );

      if ( !g.kcycle )//my addition
	vector_PRECISION_copy( &(l->p_PRECISION.x), phi, threading->start_index[l->depth], threading->end_index[l->depth], l );
      res = _RES;
    } //END: for ( int i=0; i<l->n_cy; i++ )
  } else { // if AMG is not chosen or at the bottom
    smoother_PRECISION( phi, Dphi, eta, (l->depth==0)?l->n_cy:l->post_smooth_iter, res, l, threading );
  }
}

void vcycle_timing_PRECISION( int n, level_struct *l, struct Thread *threading ) {
  
  ASSERT( g.mixed_precision );
  vector_PRECISION v1, v2;
  vector_PRECISION_init(&v1);
  vector_PRECISION_init(&v2);
  
  double t0=0, t1=0;
  vector_PRECISION_alloc(&v1, _INNER, 1, l, threading);
  vector_PRECISION_alloc(&v2, _INNER, 1, l, threading);

  START_LOCKED_MASTER(threading)
  vector_PRECISION_define_random( &v2, 0, l->inner_vector_size, l );
  END_LOCKED_MASTER(threading)
  
  START_MASTER(threading)
  t0 = MPI_Wtime();
  END_MASTER(threading)
  for ( int i=0; i<n; i++ ) {
    vcycle_PRECISION( &v1, NULL, &v2, _NO_RES, l, threading );
  }
  START_MASTER(threading)
  t1 = MPI_Wtime();
  printf0("100 v-cycles: %le seconds\n", t1-t0 );
  END_MASTER(threading)

  START_LOCKED_MASTER(threading)
  vector_PRECISION_free(&v1, l, threading);
  vector_PRECISION_free(&v2, l, threading);
  END_LOCKED_MASTER(threading)
}
