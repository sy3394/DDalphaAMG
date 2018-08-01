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

/*void vector_PRECISION_alloc( vector_PRECISION *vec, const int type, int num_vect, level_struct *l ) {
  
  MALLOC( vec->vector_buffer, complex_PRECISION, num_vect );
}

void vector_PRECISION_free( vector_PRECISION *vec, const int type, int num_vect, level_struct *l ) {
  
  FREE( vec->vector_buffer, complex_PRECISION, num_vect );
}
*/

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
/*
void vector_PRECISION_test_routine( vector_PRECISION *vec, level_struct *l, struct Thread *threading ) {

}*/
