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
 * new file from sbacchio
 * checked: 12/03/2019
 * 1st cleanup:12/18/2019
 */

#include "main.h"

void vector_PRECISION_init( vector_PRECISION *vec ) {
  
  vec->vector_buffer = NULL;
}


void vector_PRECISION_alloc( vector_PRECISION *vec, const int type, int num_vect, level_struct *l, Thread *threading ) {

  if ( num_vect % num_loop != 0 )
    error0("vector_PRECISION_alloc: #vectors should be multiple of num_loop for vectorization\n");

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
  case _EVEN_INNER:
    PUBLIC_MALLOC( vec->vector_buffer, complex_PRECISION, (l->num_inner_lattice_sites/2)*l->num_lattice_site_var*num_vect );
    vec->size = (l->num_inner_lattice_sites/2)*l->num_lattice_site_var;
    break;
  }

  vec->type         = type;
  vec->num_vect     = num_vect;
  vec->num_vect_now = num_loop;
  vec->start        = 0;
  vec->end          = vec->size*num_vect;
  vec->layout       = _NVEC_INNER;//default layout is vector runs fastest
  vec->l            = l; 
}


void vector_PRECISION_free( vector_PRECISION *vec, level_struct *l, struct Thread *threading ) {//make sure threading (used in PUBLIC_FREE) is not freed yet!!!!
  // PUBLIC_FREE requires threading to exist (in most cases no_threading)!!!!!!!!!!!
  switch (vec->type){
  case _ORDINARY :   PUBLIC_FREE( vec->vector_buffer, complex_PRECISION, l->vector_size*vec->num_vect );
    break;
  case _SCHWARZ :    PUBLIC_FREE( vec->vector_buffer, complex_PRECISION, l->schwarz_vector_size*vec->num_vect );
    break;
  case _INNER :      PUBLIC_FREE( vec->vector_buffer, complex_PRECISION, l->inner_vector_size*vec->num_vect );
    break;
  case _EVEN_INNER : PUBLIC_FREE( vec->vector_buffer, complex_PRECISION, (l->num_inner_lattice_sites/2)*l->num_lattice_site_var*vec->num_vect );
    break;
  }
}


// set all entries of a slice of phi from start to end to complex_PRECISION value
void vector_PRECISION_define( vector_PRECISION *phi, complex_PRECISION value, int start, int end, level_struct *l ) {

  int thread = omp_get_thread_num();
  if(thread == 0 && start != end)
    PROF_PRECISION_START( _SET );
  if ( phi->vector_buffer != NULL ) {
    int i, j, jj;
    for ( i=start; i<end; i++ )
      VECTOR_LOOP(j, phi->num_vect, jj, phi->vector_buffer[i*phi->num_vect+j+jj] = value;)
  } else {
    error0("Error in \"vector_PRECISION_define\": pointer is null\n");
  }
  if(thread == 0 && start != end)
    PROF_PRECISION_STOP( _SET, 1 );
}


// set every entry of all vectors in phi to alpha: similar to vector_PRECISION_define
void set_boundary_PRECISION( vector_PRECISION *phi, complex_PRECISION alpha, level_struct *l, struct Thread *threading ) {
  
  PROF_PRECISION_START( _SET, threading );
  int i, j, jj, nvec = phi->num_vect;
  int start, end;
  compute_core_start_end_custom( l->inner_vector_size, l->vector_size , &start, &end, l, threading, l->num_lattice_site_var);

  //SYNC_CORES(threading)<-this was commented out last!!!!
  for( i=start; i<end; i++ )
    VECTOR_LOOP( j, nvec, jj, phi->vector_buffer[i*nvec+j+jj] = alpha; ) 

  SYNC_CORES(threading)
  PROF_PRECISION_STOP( _SET, (double)(l->vector_size-l->inner_vector_size)/(double)l->inner_vector_size, threading );
}


void vector_PRECISION_define_random( vector_PRECISION *phi, int start, int end, level_struct *l ) {
  
  int thread = omp_get_thread_num();
  if(thread == 0 && start != end)
    PROF_PRECISION_START( _SET );
  if ( phi != NULL ) {
    int i;
    int n_vect = phi->num_vect;
    for ( i=start*n_vect; i<end*n_vect; i++ )
      phi->vector_buffer[i] = (PRECISION)(((double)rand()/(double)RAND_MAX))-0.5 + ( (PRECISION)((double)rand()/(double)RAND_MAX)-0.5)*_Complex_I;
  } else {
    error0("Error in \"vector_PRECISION_define_random\": pointer is null\n");
  }
  if(thread == 0 && start != end)
    PROF_PRECISION_STOP( _SET, 1 );
}

// z <- x with z->vector_buffer starting at start pos of x (Caution: z and x share the pointer)
void vector_PRECISION_duplicate( vector_PRECISION *z, vector_PRECISION *x, int start, level_struct *l ) {
  /**************
   * start: site index
   *************/
  
  z->vector_buffer = x->vector_buffer + start*l->num_lattice_site_var*x->num_vect;
  z->size          = x->size - start*l->num_lattice_site_var;
  z->type          = x->type;
  z->num_vect      = x->num_vect;
  z->num_vect_now  = x->num_vect_now;
  z->start         = start*l->num_lattice_site_var*x->num_vect;
  z->end           = x->end;
  z->layout        = x->layout;
  z->l             = x->l;
}

// dir == 1: z <-x[i]; dir==-1: z[i] <- x
void vector_PRECISION_copy2( vector_PRECISION *z, vector_PRECISION *x, int vec_ind, int length, int dir, int start, int end,  level_struct *l ) {

  //if(z == x) return;

  if ( !(length % num_loop ==0 || length == 1) )
    error0("vector_PRECISION_copy2: assumptions are not met\n");

  int i, j, jj;
  int thread = omp_get_thread_num();
  if(thread == 0)
    PROF_PRECISION_START( _CPY );

  if ( length == 1) {
    if ( dir == 1 )
      for( i=start; i<end; i++)
	z->vector_buffer[i*z->num_vect] = x->vector_buffer[i*x->num_vect+vec_ind];
    else 
      for( i=start; i<end; i++)
	z->vector_buffer[i*z->num_vect+vec_ind] = x->vector_buffer[i*x->num_vect];
  } else{
    if ( dir == 1 )
      for( i=start; i<end; i++)
	VECTOR_LOOP(j, length, jj, z->vector_buffer[i*z->num_vect+j+jj] = x->vector_buffer[i*x->num_vect+vec_ind+j+jj];)
    else
      for( i=start; i<end; i++)
	VECTOR_LOOP(j, length, jj, z->vector_buffer[i*z->num_vect+vec_ind+j+jj] = x->vector_buffer[i*x->num_vect+j+jj];)
  }

  if(thread == 0)
    PROF_PRECISION_STOP( _CPY, (double)(z->size)/(double)l->inner_vector_size );
}

// z <- x
void vector_PRECISION_copy( vector_PRECISION *z, vector_PRECISION *x, int start, int end, level_struct *l ) {

  if(z == x) return;
  if( z->num_vect < x->num_vect_now ) //or introduce an arg n to specify how many vector in x to be copied and impose n <= z->num_vect;or use g.num_vect_now
    error0("Error: The number of vectors to be copied should be smaller than the number of vectors receiving buffers can hold\n");

  int i, j, jj;
  int thread = omp_get_thread_num();
  if(thread == 0 && start != end)
    PROF_PRECISION_START( _CPY );

  for( i=start; i<end; i++)
    VECTOR_LOOP(j, x->num_vect_now, jj, z->vector_buffer[i*z->num_vect+j+jj] = x->vector_buffer[i*x->num_vect+j+jj];)
  
  if(thread == 0 && start != end)
    PROF_PRECISION_STOP( _CPY, (double)(end-start)/(double)l->inner_vector_size );
}

// z <- alpha*x
void vector_PRECISION_scale( vector_PRECISION *z, vector_PRECISION *x, complex_PRECISION *alpha, int k, int start, int end, level_struct *l ) {

  int i, j, jj;
  int thread = omp_get_thread_num();
  if(thread == 0 && start != end)
  PROF_PRECISION_START( _LA6 );

  for( i=start; i<end; i++)
    VECTOR_LOOP(j, x->num_vect_now, jj, z->vector_buffer[i*z->num_vect+j+jj] = alpha[k*x->num_vect_now+j+jj]*x->vector_buffer[i*x->num_vect+j+jj];)

  if(thread == 0 && start != end)
  PROF_PRECISION_STOP( _LA6, (double)(end-start)/(double)l->inner_vector_size );
}

void vector_PRECISION_real_scale( vector_PRECISION *z, vector_PRECISION *x, complex_PRECISION *alpha, 
				      int n, int opt, int start, int end, level_struct *l ) { 
  /********************************
   * Input:
   *  int n: #sets of vectors used (this option is used in linsolve_generic.c)
   *  complex_PRECISION *alpha: contains n sets of x->num_vect_now entries
   *  opt == 0: z = alpha*x
   *      == 1: z = (1/alpha)*x    
   *
   *******************************/
  
  if ( z->num_vect_now != x->num_vect_now ) //can also involve g.num_vect_now!!!!!!
    error0("vector_PRECISION_real_scale: Expect in and out bandles of vectors contain the same # of vectors in use");
  int i, j, jj, nvec = x->num_vect_now;//MIN(z->num_vect,x->num_vect);//!!!!!!!!
  PRECISION r_alpha[nvec];

  if(opt){
    VECTOR_LOOP(j, nvec, jj, r_alpha[j+jj]=1.0/creal_PRECISION(alpha[n*nvec+j+jj]);)
  }else{
    VECTOR_LOOP(j, nvec, jj, r_alpha[j+jj]=creal_PRECISION(alpha[n*nvec+j+jj]);)
  }

  int thread = omp_get_thread_num(); 
  if(thread == 0 && start != end) 
    PROF_PRECISION_START( _RS ); 

  if(z == x)
    for( i=start; i<end; i++)// does this make differece in speed?????
      VECTOR_LOOP(j, nvec, jj, z->vector_buffer[i*z->num_vect+j+jj] *= r_alpha[j+jj];)
  else
    for( i=start; i<end; i++)
      VECTOR_LOOP(j, nvec, jj, z->vector_buffer[i*z->num_vect+j+jj] = r_alpha[j+jj]*x->vector_buffer[i*x->num_vect+j+jj];)

  if(thread == 0 && start != end) 
    PROF_PRECISION_STOP( _RS, (double)(end-start)/(double)l->inner_vector_size ); 
}

// z=x+y
void vector_PRECISION_plus( vector_PRECISION *z, vector_PRECISION *x, vector_PRECISION *y, int start, int end, level_struct *l ) {

  int i, j, jj, nvec = x->num_vect_now;
  int thread = omp_get_thread_num();

  if(thread == 0 && start != end)
    PROF_PRECISION_START( _LA2 );

  for( i=start; i<end; i++)
    VECTOR_LOOP(j, nvec, jj, z->vector_buffer[i*z->num_vect+j+jj] = x->vector_buffer[i*x->num_vect+j+jj] + y->vector_buffer[i*y->num_vect+j+jj];)

  if(thread == 0 && start != end)
    PROF_PRECISION_STOP( _LA2, (double)(end-start)/(double)l->inner_vector_size );
}

// z <- x - y
void vector_PRECISION_minus( vector_PRECISION *z, vector_PRECISION *x, vector_PRECISION *y, int start, int end, level_struct *l ) {

  int i, j, jj, nvec = x->num_vect_now;
  int thread = omp_get_thread_num();

  if(thread == 0 && start != end)
    PROF_PRECISION_START( _LA2 );

  for( i=start; i<end; i++)
    VECTOR_LOOP(j, nvec, jj, z->vector_buffer[i*z->num_vect+j+jj] = x->vector_buffer[i*x->num_vect+j+jj] - y->vector_buffer[i*y->num_vect+j+jj];)

  if(thread == 0 && start != end)
    PROF_PRECISION_STOP( _LA2, (double)(end-start)/(double)l->inner_vector_size );
}

// z = x + sign*alpha*y // why?????? we can simply negate alpha before app of this function
// Assume: sign= \pm 1
void vector_PRECISION_saxpy( vector_PRECISION *z, vector_PRECISION *x, vector_PRECISION *y, complex_PRECISION *alpha, int k, int sign, int start, int end, level_struct *l ) {
  /*********
   * k: alpha's index for each bundle of vectors
   *********/
  if ( z->num_vect_now != x->num_vect_now || z->num_vect_now != y->num_vect_now )
    error0("vector_PRECISION_saxpy: z->num_vect_now != x-<num_vect_now || z->num_vect_now !=y->num_vect_now\n");

  int i, j, jj, nvec = x->num_vect_now;
  int thread = omp_get_thread_num();

  if (thread == 0 && start != end )
    PROF_PRECISION_START( _LA8 );

  //int sgn = ( sign == 1 )? 1:-1;//we can shorten!!!!!!!!!
  if ( sign == 1)
    for( i=start; i<end; i++) 
      VECTOR_LOOP(j, nvec, jj, z->vector_buffer[i*z->num_vect+j+jj] = x->vector_buffer[i*x->num_vect+j+jj] + alpha[k*nvec+j+jj]*y->vector_buffer[i*y->num_vect+j+jj];)
  else
    for( i=start; i<end; i++)
      VECTOR_LOOP(j, nvec, jj, z->vector_buffer[i*z->num_vect+j+jj] = x->vector_buffer[i*x->num_vect+j+jj] - alpha[k*nvec+j+jj]*y->vector_buffer[i*y->num_vect+j+jj];)

  if( thread == 0 && start != end )
    PROF_PRECISION_STOP( _LA8, (double)(end-start)/(double)l->inner_vector_size );
}

// float is used in MP setting!!!!!!!??????
void vector_PRECISION_multi_saxpy( vector_PRECISION *z, vector_PRECISION *V, complex_PRECISION *alpha,
				       int sign, int count, int start, int end, level_struct *l, struct Thread *threading ) {
  /********************************
   * z <- z + sign*alpha*V
   * z: a vector
   * V: a set of vectors
   * count: #vectors in V
   * complex_PRECISION *alpha: contains count many sets of z->num_vect_now many entries
   ********************************/
  
  int c, i, j, jj, nvec = z->num_vect_now;
  int thread = omp_get_thread_num();
  if (thread == 0 && start != end )
    PROF_PRECISION_START( _LA8 );
  /*
  complex_PRECISION alpha_signed[count*z->num_vect];
  for ( c=0; c<count; c++ )
    VECTOR_LOOP(j, z->num_vect, jj, alpha_signed[c*z->num_vect+j+jj] = sign*alpha[c*z->num_vect+j+jj];)
  */
  for ( c=0; c<count; c++ )
    for ( i=start; i<end; i++)
      VECTOR_LOOP(j, nvec, jj, z->vector_buffer[i*z->num_vect+j+jj] += sign*alpha[c*nvec+j+jj]*V[c].vector_buffer[i*V[c].num_vect+j+jj];)

  if( thread == 0 && start != end )
    PROF_PRECISION_STOP( _LA8, (PRECISION)(count) );
}

// check to see if vec1 and vec1 are of the same structure
void vector_PRECISION_check_comp( vector_PRECISION *vec1, vector_PRECISION *vec2) {

  if(vec1->num_vect != vec2->num_vect)
    error0("Error: The current number of vectors used have to be the same in both vectors\n");

  if(vec1->l->level != vec2->l->level)
    error0("Error: The level of multigrid must be the same in both vectors\n");

  if(vec1->type != vec2->type)
    error0("Error: The type must be the same in both vectors\n");

}

// change the order of a bundle of vectors
// to be used only by master
void vector_PRECISION_change_layout( vector_PRECISION *vec_out, vector_PRECISION *vec_in, const int layout, struct Thread *threading ) {
  
  if(vec_in->layout==layout) return;
  vector_PRECISION_check_comp( vec_out, vec_in );

  level_struct *l = vec_in->l;
  
  PROF_PRECISION_START( _RL );
  int n, i, nvec = vec_in->num_vect, size = vec_in->size;
  vector_PRECISION vec_tmp;
  if ( vec_in->vector_buffer == vec_out->vector_buffer ) {
    vector_PRECISION_init( &vec_tmp );
    vector_PRECISION_alloc( &vec_tmp, vec_in->type, vec_in->num_vect, vec_in->l, no_threading );
  } else {
    vec_tmp = *vec_out;
  }

  switch (layout){
    case _NVEC_OUTER : // from vectors->spins->site (fastest->slowest) to spins->sites->vectors 
      for( n=0; n<nvec; n++ )
	for( i=0; i<size; i++ ){
	  vec_tmp.vector_buffer[n*size+i] = vec_in->vector_buffer[i*nvec+n];}
      vec_out->layout = _NVEC_OUTER;
      break;
    case _NVEC_INNER : // from spins->sites->vectors (fastest->slowest) to vectors->spins->site
      for( i=0; i<size; i++ )
	for( n=0; n<nvec; n++ )
	    vec_tmp.vector_buffer[i*nvec+n] = vec_in->vector_buffer[n*size+i];
      vec_out->layout = _NVEC_INNER;
      break;
  }

  if( vec_in->vector_buffer == vec_out->vector_buffer ){
    vector_PRECISION_copy( vec_out, &vec_tmp, 0, size, vec_out->l );
    vector_PRECISION_free(&(vec_tmp), vec_in->l, no_threading);
  }
  PROF_PRECISION_STOP( _RL, 1 );
}

// should be intra-process permutation!!!!!!!!!!
void trans_PRECISION( vector_PRECISION *out, vector_double *in, int *tt, level_struct *l, struct Thread *threading ) {
  /*********************************************************************************************************************
   * Assume: in->num_vect == out->num_vect && l->depth == 0 && *in and *out have separate memory allocations
   * Description: Permute the order of *in according to the translation table *tt and store the result in *out
   *              *out = P[*in]
   *********************************************************************************************************************/

  if ( l->depth != 0 || (double *) out->vector_buffer == (double *) in->vector_buffer )
    error0("trans_PRECISION: assumptions are not met0\n");
  if ( in->num_vect_now > out->num_vect )
    error0("trans_PRECISION: assumptions are not met1\n");

  int i, j, jj, k, index;
  buffer_PRECISION out_pt = out->vector_buffer; buffer_double in_pt = in->vector_buffer;
  int start = threading->start_site[l->depth];
  int end   = threading->end_site[l->depth];

  SYNC_CORES(threading)
  START_NO_HYPERTHREADS(threading)
#ifdef HAVE_TM1p1
  if( g.n_flavours == 2 )
    for ( i=start; i<end; i++ ) {
      index = tt[i];
      out_pt = out->vector_buffer + 24*index;
      in_pt  = in->vector_buffer + 24*i;
      VECTOR_LOOP( *out_pt = (complex_PRECISION) *in_pt; out_pt++; in_pt++; )//?????
    }
  else
#endif
  for ( i=start; i<end; i++ ) {
    index = tt[i];
    out_pt = out->vector_buffer + 12*index*out->num_vect;
    in_pt  = in->vector_buffer + 12*i*in->num_vect;
    for( k=0; k<12; k++){
      VECTOR_LOOP(j, in->num_vect_now, jj, *out_pt = (complex_PRECISION) *in_pt;
                                            out_pt++;
                                            in_pt++;)
      out_pt += out->num_vect-in->num_vect_now;
      in_pt  += in->num_vect-in->num_vect_now;
    }
  }
  END_NO_HYPERTHREADS(threading)
  SYNC_CORES(threading)
}

// should be intra-process permutation!!!!!!!!!!
void trans_back_PRECISION( vector_double *out, vector_PRECISION *in, int *tt, level_struct *l, struct Thread *threading ) {
  /*********************************************************************************************************************
   * Assume: in->num_vect == out->num_vect && l->depth == 0 && *in and *out have separate memory allocations
   * Description: Permute back the order of *in according to the translation table *tt and store the result in *out
   *              *out = P[*in] 
   *********************************************************************************************************************/

  if ( l->depth != 0 || (double *) out->vector_buffer == (double *) in->vector_buffer )
    error0("trans_PRECISION: assumptions are not met\n");
  if ( in->num_vect_now > out->num_vect )
    error0("trans_PRECISION: assumptions are not met\n");
  
  int i, j, jj, k, index;
  buffer_double out_pt = out->vector_buffer; buffer_PRECISION in_pt = in->vector_buffer;
  int start = threading->start_site[l->depth];
  int end   = threading->end_site[l->depth];
  
  SYNC_CORES(threading)
  START_NO_HYPERTHREADS(threading)
#ifdef HAVE_TM1p1
  if( g.n_flavours == 2 )
    for ( i=start; i<end; i++ ) {
      index = tt[i];
      in_pt = in->vector_buffer + 24*index;
      out_pt = out->vector_buffer + 24*i;
      FOR24( *out_pt = (complex_double) *in_pt; out_pt++; in_pt++; )
    }
  else
#endif
  for ( i=start; i<end; i++ ) {
    index = tt[i];
    in_pt = in->vector_buffer + 12*index*in->num_vect;
    out_pt = out->vector_buffer + 12*i*out->num_vect;
    for( k=0; k<12; k++){
      VECTOR_LOOP(j, in->num_vect_now, jj, *out_pt = (complex_double) *in_pt;
                                            out_pt++;
                                            in_pt++;)
      out_pt += out->num_vect-in->num_vect_now;
      in_pt  += in->num_vect-in->num_vect_now;
    }
  }
  END_NO_HYPERTHREADS(threading)
  SYNC_CORES(threading)
}

/*******************   TM1p1  ***********************************************************************************************/

void two_flavours_to_serial_PRECISION( vector_PRECISION *flav1, vector_PRECISION *flav2, vector_PRECISION *serial, level_struct *l, struct Thread *threading ) {

#ifdef HAVE_TM1p1

  /*
   * Order: spin0and1 of flav1
   *        spin0and1 of flav2
   *        spin2and3 of flav1
   *        spin2and3 of flav2
   */
  buffer_PRECISION serial_end;
  buffer_PRECISION serial_pt = serial->vector_buffer, flav1_pt = flav1->vector_buffer, flav2_pt = flav2->vector_buffer;
  
  if( g.n_flavours == 2 ) {
    serial_end = serial->vector_buffer + threading->end_index[l->depth];
    serial_pt += threading->start_index[l->depth];
    flav1_pt += threading->start_index[l->depth]/2;
    flav2_pt += threading->start_index[l->depth]/2;
  }
  else {
    serial_end = serial->vector_buffer + threading->end_index[l->depth]*2;
    serial_pt += threading->start_index[l->depth]*2;
    flav1_pt += threading->start_index[l->depth];
    flav2_pt += threading->start_index[l->depth];
  }

  while ( serial_pt < serial_end ) {
    FOR6( *serial_pt = (*flav1_pt); serial_pt++; flav1_pt++; )
    FOR6( *serial_pt = (*flav2_pt); serial_pt++; flav2_pt++; )
    FOR6( *serial_pt = (*flav1_pt); serial_pt++; flav1_pt++; )
    FOR6( *serial_pt = (*flav2_pt); serial_pt++; flav2_pt++; )
  }
#else
  START_MASTER(threading)
  warning0("two_flavours_to_serial_PRECISION called without HAVE_TM1p1 defined\n");
  END_MASTER(threading)
#endif
    
}

void serial_to_two_flavours_PRECISION( vector_PRECISION *flav1, vector_PRECISION *flav2, vector_PRECISION *serial, level_struct *l, struct Thread *threading ) {

#ifdef HAVE_TM1p1
  buffer_PRECISION serial_end;
  buffer_PRECISION serial_pt = serial->vector_buffer, flav1_pt = flav1->vector_buffer, flav2_pt = flav2->vector_buffer;

  if( g.n_flavours == 2 ) {
    serial_end = serial->vector_buffer + threading->end_index[l->depth];
    serial_pt += threading->start_index[l->depth];
    flav1_pt += threading->start_index[l->depth]/2;
    flav2_pt += threading->start_index[l->depth]/2;
  }
  else {
    serial_end = serial->vector_buffer + threading->end_index[l->depth]*2;
    serial_pt += threading->start_index[l->depth]*2;
    flav1_pt += threading->start_index[l->depth];
    flav2_pt += threading->start_index[l->depth];
  }

  while ( serial_pt < serial_end ) {
    FOR6( *flav1_pt = (*serial_pt); serial_pt++; flav1_pt++; )
    FOR6( *flav2_pt = (*serial_pt); serial_pt++; flav2_pt++; )
    FOR6( *flav1_pt = (*serial_pt); serial_pt++; flav1_pt++; )
    FOR6( *flav2_pt = (*serial_pt); serial_pt++; flav2_pt++; )
  }
#else
  START_MASTER(threading)
  warning0("two_flavours_to_serial_PRECISION called without HAVE_TM1p1 defined\n");
  END_MASTER(threading)
#endif
    
}

/*******************  TEST ROUTINES  ********************************************/

void free_alloc_PRECISION( level_struct *l, int n_v_old, int n_v ) {
/*
 * n_v_old : old number of vectors
 * n_v : new number of vectors
 */
  int i;
//This might cause an issue regarding # vectors !!!!!!!!
  schwarz_PRECISION_struct *s = &(l->s_PRECISION);
  operator_PRECISION_struct *op = &(l->op_PRECISION);
  error0("free_alloc_PRECISIO\n");
//*****
  int tm1p1 = 1;
#ifdef HAVE_TM1p1
  tm1p1 = 2;
#endif
  
  int svs = l->schwarz_vector_size*n_v_old; 

  if ( l->depth == 0 )
    for (int i=0; i<4; i++ )
      vector_PRECISION_free( &(s->oe_buf[i]), l, no_threading );
  
  for ( i=0; i<4; i++ )
    vector_PRECISION_free( &(s->buf[i]), l, no_threading );
  if ( g.method == 1 )
    vector_PRECISION_free( &(s->buf[4]), l, no_threading );
  for ( i=0; i<3; i++ )
    FREE( s->local_minres_buffer[i], complex_PRECISION, svs*tm1p1 );

  for (int i=0; i<2; i++ )
    vector_PRECISION_free( &(l->sbuf_PRECISION[i]), l, no_threading );
//
  if ( l->depth == 0 ) { 
    int its = (l->num_lattice_site_var/2)*l->num_lattice_sites*n_v_old;
    FREE( op->prnT, complex_PRECISION, its*tm1p1*8 );
  }

  for (int k=0; k<2; k++ )
    vector_PRECISION_free( &(op->buffer[k]), l, no_threading );
  FREE( op->buffer, vector_PRECISION, 2 );
//****
  svs = l->schwarz_vector_size*n_v;

  if ( l->depth == 0 ) {
    for (int i=0; i<4; i++ ) {
      vector_PRECISION_init( &(s->oe_buf[i]) );
      vector_PRECISION_alloc( &(s->oe_buf[i]), _INNER, n_v*tm1p1, l, no_threading );
    }
  }

  for ( i=0; i<4; i++ ) {
    vector_PRECISION_init( &(s->buf[i]) );
    vector_PRECISION_alloc( &(s->buf[i]), (i==0)?((l->depth==0)?_INNER:_ORDINARY):_SCHWARZ, n_v*tm1p1, l, no_threading );
  }
  
  if ( g.method == 1 ){
    vector_PRECISION_init( &(s->buf[4]) );
    vector_PRECISION_alloc( &(s->buf[4]), _SCHWARZ, n_v*tm1p1, l, no_threading );
  }

  for ( i=0; i<3; i++ ) {
    s->local_minres_buffer[i] = NULL;
    MALLOC( s->local_minres_buffer[i], complex_PRECISION, svs*tm1p1 );
  }

  for ( i=0; i<2; i++ ) {
    vector_PRECISION_init( &(l->sbuf_PRECISION[i]) );
    vector_PRECISION_alloc( &(l->sbuf_PRECISION[i]), (l->depth==0)?_INNER:_ORDINARY, n_v*tm1p1, l, no_threading );
  }
//
  op->prnT = NULL;
  if ( l->depth == 0 ) { 
    int its = (l->num_lattice_site_var/2)*l->num_lattice_sites*n_v;
    MALLOC( op->prnT, complex_PRECISION, its*tm1p1*8 );
    op->prnZ = op->prnT + its*tm1p1; op->prnY = op->prnZ + its*tm1p1; op->prnX = op->prnY + its*tm1p1;
    op->prpT = op->prnX + its*tm1p1; op->prpZ = op->prpT + its*tm1p1; op->prpY = op->prpZ + its*tm1p1; op->prpX = op->prpY + its*tm1p1;
  }

  op->buffer = NULL;
  MALLOC( op->buffer, vector_PRECISION, 2 );
  for (int k=0; k<2; k++ ){
    vector_PRECISION_init( &(op->buffer[k]) );
    vector_PRECISION_alloc( &(op->buffer[k]), _ORDINARY, n_v*tm1p1, l, no_threading );
  }

}
