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
 * checked:11/29/2019: need work on handling of different #vectors
 * changed from sbacchio
 * checked:12/05/2019
 * 1st cleanup:12/18/2019
 */

#include "main.h"

void interpolation_PRECISION_alloc( level_struct *l ) {
  
  int k, n = l->num_eig_vect;
  
  MALLOC( l->is_PRECISION.eigenvalues, complex_PRECISION, n );
  MALLOC( l->is_PRECISION.operator, complex_PRECISION, n*l->inner_vector_size );

  vector_PRECISION_alloc(&(l->is_PRECISION.interpolation_vec), _ORDINARY, n, l, no_threading );
  vector_PRECISION_alloc(&(l->is_PRECISION.test_vector_vec), _INNER, n, l, no_threading );
  l->is_PRECISION.interpolation_vec.num_vect_now = n;
  l->is_PRECISION.test_vector_vec.num_vect_now = n;

}

void interpolation_PRECISION_free( level_struct *l ) {
  
  int n = l->num_eig_vect;

  FREE( l->is_PRECISION.eigenvalues, complex_PRECISION, n );
  FREE( l->is_PRECISION.operator, complex_PRECISION, n*l->inner_vector_size );

  vector_PRECISION_free(&(l->is_PRECISION.interpolation_vec), l, no_threading );
  vector_PRECISION_free(&(l->is_PRECISION.test_vector_vec), l, no_threading );
}

void interpolation_PRECISION_dummy_alloc( level_struct *l ) {
}

void interpolation_PRECISION_dummy_free( level_struct *l ) {
}

void define_interpolation_PRECISION_operator( vector_PRECISION *interpolation, level_struct *l, struct Thread *threading ) {
  /****************************************
   *INPUT:
   *  vector_PRECISION *interpolation: test vectors
   *OUTPUT
   *  l->is_PRECISION.operator: interpolation operator constructed out of test vectors
   *                            is in raw-major order and of inner_vector_size x num_eig_vect
   ***************************************/
  int i, j, jj, num_eig_vect = l->num_eig_vect;
  complex_PRECISION *operator = l->is_PRECISION.operator;
  int start, end;
  compute_core_start_end_custom(0, l->inner_vector_size, &start, &end, l, threading, l->num_lattice_site_var );

  if (interpolation->num_vect != num_eig_vect )
    error0("define_interpolation_PRECISION_operator: num_eig_vect and #vectors in interpolation are different\n");

  SYNC_CORES(threading)
  operator += start*num_eig_vect;
  for ( i = start; i<end; i++ )
    VECTOR_LOOP( j, num_eig_vect, jj, *operator = interpolation->vector_buffer[i*interpolation->num_vect+j+jj]; operator++;)
  SYNC_CORES(threading)
}

//work not distributed among threads
void interpolate_PRECISION( vector_PRECISION *phi, vector_PRECISION *phi_c, level_struct *l, struct Thread *threading ) {
  /**********************************************
   * Assume: phi.num_vect == phi_c.num_vect && l->level is the level where phi is defined
   * Input:
   *  vector_PRECISION *phi_c
   *  level_struct *l
   *  struct Thread *threading
   * Output
   *  vector_PRECISION *phi (+=interpolate_PRECISION*phi_c)
   * Description:
   *  interpolate phi_c based on eigenvectors in l to get phi
   ***********************************************/

  PROF_PRECISION_START( _PR, threading );
  int i, j, k, k1, k2, jj, jjj, sign = 1;
  int num_aggregates      = l->is_PRECISION.num_agg;
  int num_eig_vect        = l->num_eig_vect;
  int num_parent_eig_vect = l->num_parent_eig_vect;
  int aggregate_sites     = l->num_inner_lattice_sites / num_aggregates;
  int n_vect = MIN(phi->num_vect_now, phi_c->num_vect_now),  n_vect_phi = phi->num_vect, n_vect_phic = l->next_level->gs_PRECISION.transfer_buffer.num_vect;
  complex_PRECISION *operator, *phi_pt, *phi_c_pt;

  if ( n_vect == 0 )
    error0("interpolate_PRECISION: assumptions are not met\n");

  l->next_level->gs_PRECISION.transfer_buffer.num_vect_now = n_vect;
  START_LOCKED_MASTER(threading)
  vector_PRECISION_distribute( &(l->next_level->gs_PRECISION.transfer_buffer), phi_c, l->next_level );
  END_LOCKED_MASTER(threading)
  SYNC_HYPERTHREADS(threading)

  for ( i=threading->n_thread*threading->core + threading->thread; i<num_aggregates; i+=threading->n_core*threading->n_thread ) {
    phi_pt   = phi->vector_buffer + i*2*num_parent_eig_vect*aggregate_sites*n_vect_phi;
    phi_c_pt = l->next_level->gs_PRECISION.transfer_buffer.vector_buffer + i*2*num_eig_vect*n_vect_phic;
    operator = l->is_PRECISION.operator + i*2*num_eig_vect*num_parent_eig_vect*aggregate_sites;
    for ( k=0; k<aggregate_sites; k++ ) {
      for ( k1=0; k1<2; k1++ ) {
	for ( k2=0; k2<num_parent_eig_vect; k2++ ) {
	  for ( j=0; j<num_eig_vect; j++ ) {
	    VECTOR_LOOP(jj, n_vect, jjj, *phi_pt += (*operator) * phi_c_pt[j*n_vect_phic+jj+jjj]; phi_pt++;)
	    operator++;
	    phi_pt -= n_vect;
	  }
	  phi_pt += n_vect_phi;
	}
	phi_c_pt += sign*num_eig_vect*n_vect_phic; sign*=-1;
      }
    }
  }
  
  PROF_PRECISION_STOP( _PR, 1, threading );

  SYNC_HYPERTHREADS(threading)
}

void interpolate3_PRECISION( vector_PRECISION *phi, vector_PRECISION *phi_c, level_struct *l, struct Thread *threading ) {
  /**********************************************
   * Assume: phi.num_vect == phi_c.num_vect && l->level is the level where phi is defined
   * Input:
   *  vector_PRECISION *phi_c
   *  level_struct *l
   *  struct Thread *threading
   * Output
   *  vector_PRECISION *phi = interpolate(phi_c) (phi defined anew unlike interpolate_PRECISION above) 
   * Description:
   *  interpolate phi_c based on eigenvectors in l to get phi
   ***********************************************/

  PROF_PRECISION_START( _PR, threading );
  int i, j, k, k1, k2, jj, jjj, sign = 1;
  int num_aggregates      = l->is_PRECISION.num_agg;
  int num_eig_vect        = l->num_eig_vect;
  int num_parent_eig_vect = l->num_parent_eig_vect; 
  int aggregate_sites     = l->num_inner_lattice_sites / num_aggregates;
  int n_vect = MIN( phi->num_vect_now, phi_c->num_vect_now), n_vect_phi = phi->num_vect, n_vect_phic = l->next_level->gs_PRECISION.transfer_buffer.num_vect;
  complex_PRECISION *operator, *phi_pt, *phi_c_pt;

  if ( n_vect == 0 )
    error0("interpolate3_PRECISION: assumptions are not met\n");

  l->next_level->gs_PRECISION.transfer_buffer.num_vect_now = n_vect;
  START_LOCKED_MASTER(threading)
  vector_PRECISION_distribute( &(l->next_level->gs_PRECISION.transfer_buffer), phi_c, l->next_level );
  END_LOCKED_MASTER(threading)
  SYNC_HYPERTHREADS(threading)
  for ( i=threading->n_thread*threading->core + threading->thread; i<num_aggregates; i+=threading->n_core*threading->n_thread ) {
    phi_pt   = phi->vector_buffer + i*2*num_parent_eig_vect*aggregate_sites*n_vect_phi;
    phi_c_pt = l->next_level->gs_PRECISION.transfer_buffer.vector_buffer + i*2*num_eig_vect*n_vect_phic;
    operator = l->is_PRECISION.operator + i*2*num_eig_vect*num_parent_eig_vect*aggregate_sites;
    for ( k=0; k<aggregate_sites; k++ ) {
      for ( k1=0; k1<2; k1++ ) {
	for ( k2=0; k2<num_parent_eig_vect; k2++ ) {//SET phi_pt instead of set them=0 and accumulate  
	  VECTOR_LOOP(jj, n_vect, jjj, *phi_pt = phi_c_pt[0*n_vect_phic+jj+jjj] * (*operator); phi_pt++;)
	  operator++;
	  phi_pt -= n_vect;
	  for ( j=1; j<num_eig_vect; j++ ) {
	    VECTOR_LOOP(jj, n_vect, jjj, *phi_pt += phi_c_pt[j*n_vect_phic+jj+jjj] * (*operator); phi_pt++;)
	    operator++;
	    phi_pt -= n_vect;
	  }
	  phi_pt += n_vect_phi;
	}
	phi_c_pt += sign*num_eig_vect*n_vect_phic; sign*=-1;
      }
    }
  }
  PROF_PRECISION_STOP( _PR, 1, threading );

  SYNC_HYPERTHREADS(threading)
}

void restrict_PRECISION( vector_PRECISION *phi_c, vector_PRECISION *phi, level_struct *l, struct Thread *threading ) {
  /************************************************
   * Assume: phi.num_vect == phi_c.num_vect
   * Input:
   *  vector_PRECISION *phi
   *  level_struct *l
   *  struct Thread *threading
   * Output
   *  vector_PRECISION *phi_c <- restrict(phi)
   * Description:
   *  restrict phi based on eigenvectors in l to get phi_c
   ***********************************************/

  SYNC_CORES(threading)
  SYNC_HYPERTHREADS(threading)

  PROF_PRECISION_START( _PR, threading );
  int i, j, k, k1, k2, jj, jjj, sign = 1;
  int num_aggregates      = l->is_PRECISION.num_agg;
  int num_eig_vect        = l->num_eig_vect;
  int num_parent_eig_vect = l->num_parent_eig_vect;
  int aggregate_sites     = l->num_inner_lattice_sites / num_aggregates;
  int n_vect = MIN(phi->num_vect_now, phi_c->num_vect_now), n_vect_phi = phi->num_vect, n_vect_phic = l->next_level->gs_PRECISION.transfer_buffer.num_vect;
  complex_PRECISION *operator, *phi_pt, *phi_c_pt;

  if ( n_vect == 0 )
    error0("restrict_PRECISION: assumptions are not met\n");

  l->next_level->gs_PRECISION.transfer_buffer.num_vect_now = n_vect;
  for ( i=threading->n_thread*threading->core + threading->thread; i<num_aggregates; i+=threading->n_core*threading->n_thread ) {
    phi_pt   = phi->vector_buffer + i*aggregate_sites*2*num_parent_eig_vect*n_vect_phi;
    phi_c_pt = l->next_level->gs_PRECISION.transfer_buffer.vector_buffer + i*2*num_eig_vect*n_vect_phic;
    operator = l->is_PRECISION.operator + i*2*num_eig_vect*num_parent_eig_vect*aggregate_sites;

    for ( j=0; j<2*num_eig_vect; j++ )
      VECTOR_LOOP(jj, n_vect_phic, jjj, phi_c_pt[j*n_vect_phic+jj+jjj] = 0; )

    for ( k=0; k<aggregate_sites; k++ ) {
      for ( k1=0; k1<2; k1++ ) {
	for ( k2=0; k2<num_parent_eig_vect; k2++ ) {
	  for ( j=0; j<num_eig_vect; j++ ) {
	    VECTOR_LOOP(jj, n_vect, jjj, phi_c_pt[j*n_vect_phic+jj+jjj] += conj_PRECISION(*operator) * (*phi_pt); phi_pt++;)//complex conjugation on phi???
	    phi_pt -= n_vect;//move back to the head of the vector bundle.  note: In case of single rhs, *phi_pt stays fixed.
	    operator++;
	  }
	  phi_pt += n_vect_phi;
	}
	phi_c_pt += sign*num_eig_vect*n_vect_phic; sign*=-1;//switch back and forth btwn 1st spinor half and 2nd spinor half. Recall: sites is the slowest index
      }
    }
  }

  SYNC_HYPERTHREADS(threading)
  START_LOCKED_MASTER(threading)
  vector_PRECISION_gather( phi_c, &(l->next_level->gs_PRECISION.transfer_buffer), l->next_level );
  END_LOCKED_MASTER(threading)
  PROF_PRECISION_STOP( _PR, 1, threading );
}
