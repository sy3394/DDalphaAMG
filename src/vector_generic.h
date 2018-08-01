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

#ifndef VECTOR_PRECISION_HEADER
  #define VECTOR_PRECISION_HEADER

  struct Thread;
  
  void vector_PRECISION_init( vector_PRECISION *vec );
 // void vector_PRECISION_alloc( vector_PRECISION *vec, const int type, int num_vect, level_struct *l );
 void vector_PRECISION_define( vector_PRECISION *phi, complex_PRECISION value, int start, int end, level_struct *l );
 void vector_PRECISION_real_scale( vector_PRECISION *z, vector_PRECISION *x, complex_PRECISION alpha,
                                    int start, int end, level_struct *l );
 void vector_PRECISION_copy( vector_PRECISION *z, vector_PRECISION *x, int start, int end, level_struct *l ); // z := x
 // void vector_PRECISION_free( vector_PRECISION *vec, const int type, int num_vect, level_struct *l );
  
 // void vector_PRECISION_test_routine( vector_PRECISION *vec, level_struct *l, struct Thread *threading );
  
#endif
