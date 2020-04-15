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
 */

#ifndef FABULOUS_OP_H
   #define THREADING_H

  void fabulous_PRECISION_init( fabulous_PRECISION_struct *fab );
  void setup_fabulous_PRECISION( gmres_PRECISION_struct *p, int v_type, level_struct *l, struct Thread *threading );

  int64_t mvp_PRECISION(  void *user_env, int N,
			  const void *alpha, const void *XX, int ldx,
			  const void *beta, void *BB, int ldb);

  int64_t dot_product_PRECISION(void *user_env,
				int M, int N,
				const void *A_, int lda,
				const void *B_, int ldb,
				void *C_, int ldc);

  int64_t fabulous_rightprecond_PRECISION(void *user_env, int N,
					  const void *XX, int ldx,
					  void *BB, int ldb);

  void fabulous_PRECISION_free( fabulous_PRECISION_struct *fab, level_struct *l, struct Thread *threading );

#endif
