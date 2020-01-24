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
 * changed from sbacchio
 * checked:12/05/2019:some work undone in smoother
 * 1st cleanup: 12/18/2019
 * confirmed:not much diff from sbacchio:01/02/2020
 */

#include "main.h"
#include "vcycle_PRECISION.h"


// assume num_vect_now for phi and eta are defined
// phi <- V/K-cycle(eta)
void vcycle_PRECISION_new( vector_PRECISION *phi, vector_PRECISION *Dphi, vector_PRECISION *eta,
			   int res, level_struct *l, struct Thread *threading ) {

  if ( g.interpolation && l->level>0 ) {
    //int nvec = eta->num_vect_now;// we might just use this #????
    for ( int i=0; i<l->n_cy; i++ ) {
      // Pre-somooth and compute residual //presmoohting????
      if ( i==0 && res == _NO_RES ) {//printf0("vcy:first if %d\n",eta->num_vect_now);fflush(stdout);
	// if the initial guess is 0, the resid on the next level is restircted rhs
	//SYNC_MASTER_TO_ALL(threading)
	l->next_level->p_PRECISION.b.num_vect_now = eta->num_vect_now;//!!!!!!!
        restrict_PRECISION_new( &(l->next_level->p_PRECISION.b), eta, l, threading );
      } else {//printf0("vcycle_PRECISION:first else %d %d\n",phi->num_vect_now,eta->num_vect_now);fflush(stdout);
	// otherwise, compute R(b_l - D_l eta_l)
        int start = threading->start_index[l->depth];
        int end   = threading->end_index[l->depth];
	//SYNC_MASTER_TO_ALL(threading)
	l->vbuf_PRECISION[2].num_vect_now = phi->num_vect_now; l->vbuf_PRECISION[3].num_vect_now = eta->num_vect_now; l->next_level->p_PRECISION.b.num_vect_now = eta->num_vect_now;
        apply_operator_PRECISION( &(l->vbuf_PRECISION[2]), phi, &(l->p_PRECISION), l, threading );
        vector_PRECISION_minus_new( &(l->vbuf_PRECISION[3]), eta, &(l->vbuf_PRECISION[2]), start, end, l );
        restrict_PRECISION_new( &(l->next_level->p_PRECISION.b), &(l->vbuf_PRECISION[3]), l, threading );
      }

      // recursive part: compute corse-grid correction
      if ( !l->next_level->idle ) {
        START_MASTER(threading)
	if ( l->depth == 0 )
	  g.coarse_time -= MPI_Wtime();
        END_MASTER(threading)
	if ( l->level > 1 ) {//printf0("vcy: level>1 %d %d\n",l->depth,g.num_vect_now);fflush(stdout);
	  // If the next level is not at the bottom
	  if ( g.kcycle )//default
	    fgmres_PRECISION( &(l->next_level->p_PRECISION), l->next_level, threading );//in this function, call vcycle as p_PRECISION use it as prec; also use g.num_vect_now!!!!!
	  else {
	    //printf0("vcy: level>1 and not k-cycle %d %d\n", l->next_level->p_PRECISION.x.num_vect_now,l->next_level->p_PRECISION.b.num_vect_now);fflush(stdout);
	    // dive one level down 
	    //SYNC_MASTER_TO_ALL(threading)
	    l->next_level->p_PRECISION.x.num_vect_now = l->next_level->p_PRECISION.b.num_vect_now;//!!!!!!????????
	    vcycle_PRECISION_new( &(l->next_level->p_PRECISION.x), NULL, &(l->next_level->p_PRECISION.b), _NO_RES, l->next_level, threading );
	  }
	} else { // if the next level is the bottom
	  if ( g.odd_even ) {//printf0("vy:else odd even\n");fflush(stdout);
	    coarse_solve_odd_even_PRECISION_new( &(l->next_level->p_PRECISION), &(l->next_level->oe_op_PRECISION), l->next_level, threading );
	  } else {//printf0("vy:else\n");fflush(stdout);
	    fgmres_PRECISION( &(l->next_level->p_PRECISION), l->next_level, threading );
	  }
	}
        START_MASTER(threading)
	if ( l->depth == 0 )
	  g.coarse_time += MPI_Wtime();
        END_MASTER(threading)
      }
      // interpolate the correction      
      if( i == 0 && res == _NO_RES ) {
	SYNC_MASTER_TO_ALL(threading)
	l->next_level->p_PRECISION.x.num_vect_now = phi->num_vect_now;
        interpolate3_PRECISION_new( phi, &(l->next_level->p_PRECISION.x), l, threading );
      }
      else {
	//SYNC_MASTER_TO_ALL(threading)
	l->next_level->p_PRECISION.x.num_vect_now = phi->num_vect_now;
        interpolate_PRECISION_new( phi, &(l->next_level->p_PRECISION.x), l, threading );
      }
      // Perform post smoothing
      //printf0("vcy:post smooth\n");fflush(stdout);
      smoother_PRECISION_new( phi, Dphi, eta, l->post_smooth_iter, _RES, l, threading );
      res = _RES;
    }
  } else { // if AMG is not chosen or at the bottom
    smoother_PRECISION_new( phi, Dphi, eta, (l->depth==0)?l->n_cy:l->post_smooth_iter, res, l, threading );
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
  vector_PRECISION_define_random_new( &v2, 0, l->inner_vector_size, l );
  END_LOCKED_MASTER(threading)
  
  START_MASTER(threading)
  t0 = MPI_Wtime();
  END_MASTER(threading)
  for ( int i=0; i<n; i++ ) {
    vcycle_PRECISION_new( &v1, NULL, &v2, _NO_RES, l, threading );
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
