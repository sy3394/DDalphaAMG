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
 * copied:11/30/2019
 * new file from sbacchio
 * checked:12/06/2019
 * 1sr cleanup: 12/18/2019
 */

#ifndef VECTOR_PRECISION_HEADER
  #define VECTOR_PRECISION_HEADER

  struct Thread;
  
  void vector_PRECISION_init( vector_PRECISION *vec );
  void vector_PRECISION_alloc( vector_PRECISION *vec, const int type, int num_vect, level_struct *l, struct Thread *threading );
  void vector_PRECISION_free( vector_PRECISION *vec, level_struct *l, Thread *threading);

  void vector_PRECISION_define( vector_PRECISION *phi, complex_PRECISION value, int start, int end, level_struct *l );
  void set_boundary_PRECISION( vector_PRECISION *phi, complex_PRECISION alpha, level_struct *l, struct Thread *threading );
  void vector_PRECISION_define_random( vector_PRECISION *phi, int start, int end, level_struct *l ); 

void vector_PRECISION_duplicate( vector_PRECISION *z, vector_PRECISION *x, int start, level_struct *l );
void vector_PRECISION_copy2( vector_PRECISION *z, vector_PRECISION *x, int start_ind, int length, int dir, int start, int end,  level_struct *l );
  void vector_PRECISION_copy( vector_PRECISION *z, vector_PRECISION *x, int start, int end, level_struct *l );
  void vector_PRECISION_scale( vector_PRECISION *z, vector_PRECISION *x, complex_PRECISION *alpha, int k, int start, int end, level_struct *l );
  void vector_PRECISION_real_scale( vector_PRECISION *z, vector_PRECISION *x, complex_PRECISION *alpha,
					int n, int opt, int start, int end, level_struct *l );

  void vector_PRECISION_plus( vector_PRECISION *z, vector_PRECISION *x, vector_PRECISION *y, int start, int end, level_struct *l );
  void vector_PRECISION_minus( vector_PRECISION *z, vector_PRECISION *x, vector_PRECISION *y, int start, int end, level_struct *l );
  void vector_PRECISION_saxpy( vector_PRECISION *z, vector_PRECISION *x, vector_PRECISION *y, complex_PRECISION *alpha, int k, int sign, int start, int end, level_struct *l );
  void vector_PRECISION_multi_saxpy( vector_PRECISION *z, vector_PRECISION *V, complex_PRECISION *alpha,
					 int sign, int count, int start, int end, level_struct *l, struct Thread *threading );

  void vector_PRECISION_check_comp( vector_PRECISION *vec1, vector_PRECISION *vec2 );
  void vector_PRECISION_change_layout( vector_PRECISION *vec_out, vector_PRECISION *vec_in, const int layout, struct Thread *threading );
  void trans_PRECISION( vector_PRECISION *out, vector_double *in, int *tt, level_struct *l, struct Thread *threading );
  void trans_back_PRECISION( vector_double *out, vector_PRECISION *in, int *tt, level_struct *l, struct Thread *threading );
  void two_flavours_to_serial_PRECISION( vector_PRECISION *flav1, vector_PRECISION *flav2, vector_PRECISION *serial, level_struct *l, struct Thread *threading );
  void serial_to_two_flavours_PRECISION( vector_PRECISION *flav1, vector_PRECISION *flav2, vector_PRECISION *serial, level_struct *l, struct Thread *threading );


  void free_alloc_PRECISION( level_struct *l, int n_v_i, int n_v_f ); 

#endif
