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
 */

#include "main.h"

void vector_PRECISION_init( vector_PRECISION *vec ) {
  
  vec->vector_buffer = NULL;
}


void vector_PRECISION_alloc( vector_PRECISION *vec, const int type, int num_vect, level_struct *l, Thread *threading ) {

  switch (type){
  case _ORDINARY : PUBLIC_MALLOC( vec->vector_buffer, complex_PRECISION, l->vector_size*num_vect );
    break;
  case _SCHWARZ : PUBLIC_MALLOC( vec->vector_buffer, complex_PRECISION, l->schwarz_vector_size*num_vect );
    break;
  case _INNER: PUBLIC_MALLOC( vec->vector_buffer, complex_PRECISION, l->inner_vector_size*num_vect );
    break;
  }

  vec->type = type;
  vec->num_vect = num_vect;
  vec->layout = _NV_LV_SP_CL_RI;
  vec->l = l; 
}


void vector_PRECISION_free( vector_PRECISION *vec, level_struct *l, Thread *threading ) {
  
  switch (vec->type){
  case _ORDINARY : PUBLIC_FREE( vec->vector_buffer, complex_PRECISION, l->vector_size*vec->num_vect );
    break;
  case _SCHWARZ : PUBLIC_FREE( vec->vector_buffer, complex_PRECISION, l->schwarz_vector_size*vec->num_vect );
    break;
  case _INNER: PUBLIC_FREE( vec->vector_buffer, complex_PRECISION, l->inner_vector_size*vec->num_vect );
    break;
  }
}


// vector storage for PRECISION precision
void vector_PRECISION_define( vector_PRECISION *phi, complex_PRECISION value, int start, int end, level_struct *l ) {
  
  int thread = omp_get_thread_num();
  if(thread == 0 && start != end)
    PROF_PRECISION_START( _SET );
  if ( phi->vector_buffer != NULL ) {
    int i;
    for ( i=start; i<end; i++ )
      phi->vector_buffer[i] = value;
  } else {
    error0("Error in \"vector_PRECISION_define\": pointer is null\n");
  }
  if(thread == 0 && start != end)
    PROF_PRECISION_STOP( _SET, 1 );
}


void vector_PRECISION_real_scale( vector_PRECISION *z, vector_PRECISION *x, complex_PRECISION alpha, 
                                  int start, int end, level_struct *l ) { 
   
  PRECISION *r_z = (PRECISION*)z->vector_buffer, *r_x = (PRECISION*)x->vector_buffer, r_alpha = creal_PRECISION(alpha); 
  int r_start = 2*start, r_end = 2*end; 
   
  int thread = omp_get_thread_num(); 
  if(thread == 0 && start != end) 
  PROF_PRECISION_START( _LA2 ); 
   
  REAL_VECTOR_FOR( int i=r_start, i<r_end, r_z[i] = r_alpha*r_x[i], i++, l ); 
   
  if(thread == 0 && start != end) 
  PROF_PRECISION_STOP( _LA2, (double)(end-start)/(double)l->inner_vector_size ); 
}


void vector_PRECISION_copy( vector_PRECISION *z, vector_PRECISION *x, int start, int end, level_struct *l ) {

  int thread = omp_get_thread_num();
  if(thread == 0 && start != end)
  PROF_PRECISION_START( _CPY );

  VECTOR_FOR( int i=start, i<end, z->vector_buffer[i] = x->vector_buffer[i], i++, l );

  if(thread == 0 && start != end)
  PROF_PRECISION_STOP( _CPY, (double)(end-start)/(double)l->inner_vector_size );
}


void vector_PRECISION_check_compatibility( vector_PRECISION *vec1, vector_PRECISION *vec2) {

   

}


void vector_PRECISION_change_layout( vector_PRECISION *vec_out, vector_PRECISION *vec_in, const int layout, Thread *threading ) {
  
  if(vec_in->layout==layout) return;

  int n, i, s, c, lv = 0, num_s, num_c;
  vector_PRECISION vec_tmp;
  if( vec_in->vector_buffer == vec_out->vector_buffer ){
    vector_PRECISION_init( &vec_tmp );
    vector_PRECISION_alloc( &vec_tmp, vec_in->type, vec_in->num_vect, vec_in->l, no_threading );
  } else {
    vec_tmp = *vec_out;
  }

  if(vec_in->l->depth == 0){
    num_s = 4;
    num_c = 3;
  } else {
    num_s = 2;
    num_c = vec_in->l->num_parent_eig_vect;
  }

  switch (vec_in->type){
  case _ORDINARY :
    lv = vec_in->l->num_lattice_sites;
    break;
  case _SCHWARZ :
    lv = 2*vec_in->l->num_lattice_sites - vec_in->l->num_inner_lattice_sites;
    break;
  case _INNER:
    lv = vec_in->l->num_inner_lattice_sites;
    break;
  }

  switch (layout){
  case _NV_LV_SP_CL_RI :
    for( n=0; n<vec_in->num_vect; n++ )
      for( i=0; i<lv; i++ )
        for( s=0; s<num_s; s++ )
	  for( c=0; c<num_c; c++ )
	    vec_tmp.vector_buffer[INDEX_NV_LV_SP_CL( n, vec_in->num_vect, i, lv, s, num_s, c, num_c )] = vec_in->vector_buffer[INDEX_LV_SP_CL_NV( n, vec_in->num_vect, i, lv, s, num_s, c, num_c )];

    vec_out->layout = _NV_LV_SP_CL_RI;
    break;
  case _LV_SP_CL_RI_NV : 
    for( i=0; i<lv; i++ )
      for( s=0; s<num_s; s++ )
        for( c=0; c<num_c; c++ )
          for( n=0; n<vec_in->num_vect; n++ )
            vec_tmp.vector_buffer[INDEX_LV_SP_CL_NV( n, vec_in->num_vect, i, lv, s, num_s, c, num_c )] = vec_in->vector_buffer[INDEX_NV_LV_SP_CL( n, vec_in->num_vect, i, lv, s, num_s, c, num_c )];
    
    vec_out->layout = _LV_SP_CL_RI_NV;
    break;
  }
 
  if( vec_in->vector_buffer == vec_out->vector_buffer ){
     vector_PRECISION_copy( vec_out, &vec_tmp, 0, lv*num_s*num_c*vec_out->num_vect, vec_out->l );
     vector_PRECISION_free( &vec_tmp, vec_in->l, no_threading ); 
  }

}

void vector_PRECISION_test_routine( level_struct *l, struct Thread *threading ) {

  PRECISION diff = 0;
  
  vector_PRECISION vp[3];

  for(int i=0; i<3; i++){
    vector_PRECISION_init( &vp[i] );
    vector_PRECISION_alloc( &vp[i], _ORDINARY, 4, l, threading );
  }
  
  START_LOCKED_MASTER(threading)

  vector_PRECISION_define_random( &vp[0], 0, 4*l->vector_size, l );
  vector_PRECISION_copy( &vp[1], &vp[0], 0, 4*l->vector_size, l );
  vector_PRECISION_change_layout( &vp[1], &vp[1], _LV_SP_CL_RI_NV, no_threading );
  vector_PRECISION_change_layout( &vp[1], &vp[1], _NV_LV_SP_CL_RI, no_threading ); 
  vector_PRECISION_minus( &vp[2], &vp[1], &vp[0], 0, 4*l->vector_size, l );
  diff = global_norm_PRECISION( &vp[2], 0, 4*l->vector_size, l, no_threading )/
    global_norm_PRECISION( &vp[0], 0, 4*l->vector_size, l, no_threading );
  
  test0_PRECISION("depth: %d, correctness of vector PRECISION change layout: %le\n", l->depth, diff );
 
  END_LOCKED_MASTER(threading)
  for(int i=0; i<3; i++){
    vector_PRECISION_free( &vp[i], l, threading );
  }
  if ( l->level == 0 )
    return;
  else
    vector_PRECISION_test_routine(l->next_level, threading);
}
