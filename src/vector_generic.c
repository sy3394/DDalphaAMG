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
  case _ORDINARY : 
    PUBLIC_MALLOC( vec->vector_buffer, complex_PRECISION, l->vector_size*num_vect );
    vec->size = l->vector_size;
    break;
  case _SCHWARZ : 
    PUBLIC_MALLOC( vec->vector_buffer, complex_PRECISION, l->schwarz_vector_size*num_vect );
    vec->size = l->schwarz_vector_size;
    break;
  case _INNER: 
    PUBLIC_MALLOC( vec->vector_buffer, complex_PRECISION, l->inner_vector_size*num_vect );
    vec->size = l->inner_vector_size;
    break;
  }

  vec->type = type;
  vec->num_vect = num_vect;
  vec->layout = _NV_LV_SV;
  vec->l = l; 
}


void vector_PRECISION_free( vector_PRECISION *vec, level_struct *l, struct Thread *threading ) {
  
  switch (vec->type){
  case _ORDINARY : PUBLIC_FREE( vec->vector_buffer, complex_PRECISION, l->vector_size*vec->num_vect );
    break;
  case _SCHWARZ : PUBLIC_FREE( vec->vector_buffer, complex_PRECISION, l->schwarz_vector_size*vec->num_vect );
    break;
  case _INNER : PUBLIC_FREE( vec->vector_buffer, complex_PRECISION, l->inner_vector_size*vec->num_vect );
    break;
  }
}


// vector storage for PRECISION precision
void vector_PRECISION_define( vector_PRECISION *phi, complex_PRECISION value, int start, int end, level_struct *l ) {
  
  int thread = omp_get_thread_num();
  if(thread == 0 && start != end)
    PROF_PRECISION_START( _SET );
  if ( phi->vector_buffer != NULL ) {
    //int i;
    //for ( i=start; i<end; i++ )
    //  phi->vector_buffer[i] = value;
    for(int i=start; i<end; i+=PRECISION_LENGTH)
#ifdef OPTIMIZE
      #pragma unroll
      #pragma vector aligned
#endif      
      for(int j=0; j<PRECISION_LENGTH; j++)
        phi->vector_buffer[i+j] = value;
  } else {
    error0("Error in \"vector_PRECISION_define\": pointer is null\n");
  }
  if(thread == 0 && start != end)
    PROF_PRECISION_STOP( _SET, 1 );
}


void vector_PRECISION_define_new( vector_PRECISION *phi, complex_PRECISION value, level_struct *l, struct Thread *threading ) {

  int start, end;
  compute_core_start_end(0, (phi->size)*(phi->num_vect), &start, &end, l, threading);
  int thread = omp_get_thread_num();
  if(thread == 0)
    PROF_PRECISION_START( _SET );

  if ( phi->vector_buffer != NULL ) {
    int i;
    for ( i=start; i<end; i++ )
      phi->vector_buffer[i] = value;
  } else {
    error0("Error in \"vector_PRECISION_define\": pointer is null\n");
  }
  if(thread == 0)
    PROF_PRECISION_STOP( _SET, 1 );
}


void vector_PRECISION_real_scale( vector_PRECISION *z, vector_PRECISION *x, complex_PRECISION alpha, 
                                  int start, int end, level_struct *l ) { 
   
  vector_PRECISION_check_comp( z, x );
  //z->layout = x->layout;

  int thread = omp_get_thread_num(); 
  if(thread == 0 && start != end) 
  PROF_PRECISION_START( _RS ); 
  
  if(z == x){
#ifdef OPTIMIZE
    PRECISION * restrict r_x = (PRECISION*)x->vector_buffer, r_alpha = creal_PRECISION(alpha);
#else    
    PRECISION *r_x = (PRECISION*)x->vector_buffer, r_alpha = creal_PRECISION(alpha); 
#endif    
    int r_start = 2*start, r_end = 2*end;
 
    //REAL_VECTOR_FOR( int i=r_start, i<r_end, r_x[i] = r_alpha*r_x[i], i++, l );
    for(int i=r_start; i<r_end; i+=PRECISION_LENGTH)
#ifdef OPTIMIZE
      #pragma unroll
      #pragma vector aligned
#endif
      for(int j=0; j<PRECISION_LENGTH; j++)
        r_x[i+j] = r_alpha*r_x[i+j];
  } else {
#ifdef OPTIMIZE
    PRECISION * restrict r_z = (PRECISION*)z->vector_buffer, * restrict r_x = (PRECISION*)x->vector_buffer, r_alpha = creal_PRECISION(alpha);
#else    
    PRECISION *r_z = (PRECISION*)z->vector_buffer, *r_x = (PRECISION*)x->vector_buffer, r_alpha = creal_PRECISION(alpha);
#endif    
    int r_start = 2*start, r_end = 2*end;
 
    //REAL_VECTOR_FOR( int i=r_start, i<r_end, r_z[i] = r_alpha*r_x[i], i++, l );
    for(int i=r_start; i<r_end; i+=PRECISION_LENGTH)
#ifdef OPTIMIZE
      #pragma unroll
      #pragma vector aligned
#endif  
      for(int j=0; j<PRECISION_LENGTH; j++)
        r_z[i+j] = r_alpha*r_x[i+j]; 
  }
  
  if(thread == 0 && start != end) 
  PROF_PRECISION_STOP( _RS, (double)(end-start)/(double)l->inner_vector_size ); 
}


void vector_PRECISION_copy( vector_PRECISION *z, vector_PRECISION *x, int start, int end, level_struct *l ) {
 
  if(z == x) return;

  //vector_PRECISION_check_comp( z, x );
  //z->layout = x->layout;
#ifdef OPTIMIZE
  buffer_PRECISION restrict z_pt=z->vector_buffer, restrict x_pt=x->vector_buffer;
#else
  buffer_PRECISION z_pt=z->vector_buffer, x_pt=x->vector_buffer;
#endif
  int thread = omp_get_thread_num();
  if(thread == 0 && start != end)
    PROF_PRECISION_START( _CPY );
    //VECTOR_FOR( int i=start, i<end, z_pt[i] = x_pt[i], i++, l );
  for(int i=start; i<end; i+=PRECISION_LENGTH)
#ifdef OPTIMIZE
    #pragma unroll
    #pragma vector aligned
#endif
    for(int j=0; j<PRECISION_LENGTH; j++)
      z_pt[i+j] = x_pt[i+j];

  if(thread == 0 && start != end)
    PROF_PRECISION_STOP( _CPY, (double)(end-start)/(double)l->inner_vector_size );
}


void vector_PRECISION_copy_new( vector_PRECISION *z, vector_PRECISION *x, level_struct *l, struct Thread *threading ) {

  if(z == x) return;

  int start, end;
  compute_core_start_end(0, (x->size)*(x->num_vect), &start, &end, l, threading);
  buffer_PRECISION z_pt=z->vector_buffer, x_pt=x->vector_buffer;
  int thread = omp_get_thread_num();
  if(thread == 0)
    PROF_PRECISION_START( _CPY );
    
  VECTOR_FOR( int i=start, i<end, z_pt[i] = x_pt[i], i++, l );

  if(thread == 0)
    PROF_PRECISION_STOP( _CPY, (double)(end-start)/(double)l->inner_vector_size );
}


void vector_PRECISION_check_comp( vector_PRECISION *vec1, vector_PRECISION *vec2) {

  if(vec1->num_vect != vec2->num_vect)
    error0("Error: The number of vectors have to be the same in both vectors\n");

  if(vec1->l->level != vec2->l->level)
    error0("Error: The level of multigrid must be the same in both vectors\n");

  if(vec1->type != vec2->type)
    error0("Error: The type must be the same in both vectors\n");

}


void vector_PRECISION_change_layout( vector_PRECISION *vec_out, vector_PRECISION *vec_in, const int layout, struct Thread *threading ) {
  
  if(vec_in->layout==layout) return;
 
  vector_PRECISION_check_comp( vec_out, vec_in );

  int n, i, sv, lv = 0, num_sv = vec_in->l->num_lattice_site_var;
  vector_PRECISION vec_tmp;

  if( vec_in->vector_buffer == vec_out->vector_buffer ){
    vector_PRECISION_init( &vec_tmp );
    vector_PRECISION_alloc( &vec_tmp, vec_in->type, vec_in->num_vect, vec_in->l, no_threading );
  } else {
    vec_tmp = *vec_out;
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
  case _NV_LV_SV :
    for( n=0; n<vec_in->num_vect; n++ )
      for( i=0; i<lv; i++ )
        for( sv=0; sv<num_sv; sv++ )
	  vec_tmp.vector_buffer[INDEX_NV_LV_SV( n, vec_in->num_vect, i, lv, sv, num_sv )] = vec_in->vector_buffer[INDEX_LV_SV_NV( n, vec_in->num_vect, i, lv, sv, num_sv )];

    vec_out->layout = _NV_LV_SV;
    break;
  case _LV_SV_NV : 
    for( i=0; i<lv; i++ )
      for( sv=0; sv<num_sv; sv++ )
        for( n=0; n<vec_in->num_vect; n++ )
          vec_tmp.vector_buffer[INDEX_LV_SV_NV( n, vec_in->num_vect, i, lv, sv, num_sv )] = vec_in->vector_buffer[INDEX_NV_LV_SV( n, vec_in->num_vect, i, lv, sv, num_sv )];
    
    vec_out->layout = _LV_SV_NV;
    break;
  }
 
  if( vec_in->vector_buffer == vec_out->vector_buffer ){
     vector_PRECISION_copy( vec_out, &vec_tmp, 0, lv*num_sv*vec_out->num_vect, vec_out->l );
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
  vector_PRECISION_change_layout( &vp[1], &vp[1], _LV_SV_NV, no_threading );
  vector_PRECISION_change_layout( &vp[1], &vp[1], _NV_LV_SV, no_threading ); 
  vector_PRECISION_minus( &vp[2], &vp[1], &vp[0], 0, 4*l->vector_size, l );
  diff = global_norm_PRECISION( &vp[2], 0, 4*l->vector_size, l, no_threading )/
    global_norm_PRECISION( &vp[0], 0, 4*l->vector_size, l, no_threading );
  
  test0_PRECISION("depth: %d, correctness of vector PRECISION change layout: %le\n", l->depth, diff );
 
  END_LOCKED_MASTER(threading)
  for(int i=0; i<3; i++){
    vector_PRECISION_free( &vp[i], l, threading );
  }
  if ( l->level == 0 && g.method == 0)
    return;
  else
    vector_PRECISION_test_routine(l->next_level, threading);
}
