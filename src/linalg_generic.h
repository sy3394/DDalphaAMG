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

#ifndef LINALG_PRECISION_HEADER
  #define LINALG_PRECISION_HEADER
  
  struct Thread;

  void global_norm_PRECISION( PRECISION *res, vector_PRECISION *x, int start, int end, level_struct *l, struct Thread *threading );
  void global_inner_product_PRECISION( complex_PRECISION *results, vector_PRECISION *phi, vector_PRECISION *psi, 
					   int start, int end, level_struct *l, struct Thread *threading );
  void process_multi_inner_product_PRECISION( int count, complex_PRECISION *results, vector_PRECISION *phi, vector_PRECISION *psi,
						  int start, int end, level_struct *l, struct Thread *threading );

  void local_xy_over_xx_PRECISION( complex_PRECISION *res, vector_PRECISION *phi, vector_PRECISION *psi, int start, int end, level_struct *l  );

  void gram_schmidt_PRECISION( vector_PRECISION *V, const int nvec, level_struct *l, struct Thread *threading );
  void gram_schmidt_on_aggregates_PRECISION( vector_PRECISION *phi, const int num_vec, level_struct *l, struct Thread *threading );

#endif
