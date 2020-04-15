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
 * changed from sbacchio
 * checked: 12/08/2019: some reordering might be needed
 * 1st cleanup: 12/22/2019
 */

#ifndef COARSE_ODDEVEN_PRECISION_HEADER
  #define COARSE_ODDEVEN_PRECISION_HEADER

  struct Thread;

  void coarse_oddeven_alloc_PRECISION( level_struct *l );
  void coarse_oddeven_free_PRECISION( level_struct *l );
  void coarse_oddeven_setup_PRECISION( operator_PRECISION_struct *in, int reorder, level_struct *l,
                                       struct Thread *threading );
  void coarse_oddeven_PRECISION_set_self_couplings( level_struct *l, struct Thread *threading );


  void coarse_diag_ee_PRECISION_new( vector_PRECISION *y, vector_PRECISION *x, operator_PRECISION_struct *op, level_struct *l, struct Thread *threading );
  void coarse_diag_oo_PRECISION_new( vector_PRECISION *y, vector_PRECISION *x, operator_PRECISION_struct *op, level_struct *l, struct Thread *threading );
  int coarse_solve_odd_even_PRECISION_new( gmres_PRECISION_struct *p, operator_PRECISION_struct *op, level_struct *l, 
					    struct Thread *threading );
  void coarse_fabulous_solve_odd_even_PRECISION( fabulous_PRECISION_struct *fab, gmres_PRECISION_struct *p, struct Thread *threading );
  
  void coarse_apply_schur_complement_PRECISION_new( vector_PRECISION *out, vector_PRECISION *in,
						    operator_PRECISION_struct *op, level_struct *l, struct Thread *threading );
  void coarse_hopping_term_PRECISION_new( vector_PRECISION *out, vector_PRECISION *in, operator_PRECISION_struct *op,
                                      const int amount, level_struct *l, struct Thread *threading );
  void coarse_n_hopping_term_PRECISION_new( vector_PRECISION *out, vector_PRECISION *in, operator_PRECISION_struct *op,
                                        const int amount, level_struct *l, struct Thread *threading );

  void coarse_odd_even_PRECISION_test_new( vector_PRECISION *c4, vector_PRECISION *c1,
                                       level_struct *l, struct Thread *threading );


static inline void coarse_perform_fwd_bwd_subs_PRECISION_new( vector_PRECISION *x, vector_PRECISION *b, config_PRECISION A, int start, int end, level_struct *l ) {

  register int s, i, j, jj, jjj;
  int n2 = l->num_lattice_site_var;
  int oo_inv_size = SQUARE(n2);
  int nvec = b->num_vect_now, nvec_x = x->num_vect, nvec_b = b->num_vect;
  buffer_PRECISION x_pt = x->vector_buffer, b_pt = b->vector_buffer;
  
  if( nvec_x < nvec )
    error0("coarse_perform_fwd_bwd_subs_PRECISION: assumptions are not met\n");

  for ( s=start; s<end; s++ ) {  
    // solve x = U^(-1) L^(-1) b
    // forward substitution with L
    for ( i=0; i<n2; i++ ) {
      VECTOR_LOOP( jj, nvec, jjj, x_pt[i*nvec_x+jj+jjj] = b_pt[i*nvec_b+jj+jjj];)
	for ( j=0; j<i; j++ ) 
	  VECTOR_LOOP(jj, nvec, jjj, x_pt[i*nvec_x+jj+jjj] -= A[i*n2+j]*x_pt[j*nvec_x+jj+jjj];)//printf("%g ",creal_PRECISION(A[i*n2+j]));)
    }
    // backward substitution with U
    for ( i=n2-1; i>=0; i-- ) {
      for ( j=i+1; j<n2; j++ ) 
	VECTOR_LOOP(jj, nvec, jjj, x_pt[i*nvec_x+jj+jjj] -= A[i*n2+j]*x_pt[j*nvec_x+jj+jjj];)//printf("%g ",creal_PRECISION(A[i*n2+j]));)
	  VECTOR_LOOP( jj, nvec, jjj, x_pt[i*nvec_x+jj+jjj] /= A[i*(n2+1)];)//printf("%g ",creal_PRECISION(A[i*n2+j]));)
    }
    x_pt += n2*nvec_x;
    b_pt += n2*nvec_b;
    A += oo_inv_size;
  }
}

// used only in test routines
static inline void coarse_LU_multiply_PRECISION_new( vector_PRECISION *y, vector_PRECISION *x, config_PRECISION A, int start, int  end, level_struct *l ) {
  
  register int s, i, j, n2 = l->num_lattice_site_var, oo_inv_size = SQUARE(n2), jj, jjj;
  int  nvec = x->num_vect_now, nvec_x  = x->num_vect, nvec_y = y->num_vect;
  buffer_PRECISION x_pt = x->vector_buffer, y_pt = y->vector_buffer;

  if ( nvec_y < nvec )
    error0("coarse_LU_multiply_PRECISION: assumptions are not met\n");

  for ( s=start; s<end; s++ ) { 
    // y = Ax
    // multiplication with U
    for ( i=0; i<n2; i++ ) {
      VECTOR_LOOP( jj, nvec, jjj, y_pt[i*nvec_y+jj+jjj] = A[i*(n2+1)]*x_pt[i*nvec_x+jj+jjj];)
	for ( j=i+1; j<n2; j++ )
	  VECTOR_LOOP( jj, nvec, jjj, y_pt[i*nvec_y+jj+jjj] += A[i*n2+j]*x_pt[j*nvec_x+jj+jjj];)
    }
    // multiplication with L
    for ( i=n2-1; i>0; i-- )
      for ( j=0; j<i; j++ )
	VECTOR_LOOP(jj, nvec, jjj, y_pt[i*nvec_y+jj+jjj] += A[i*n2+j]*y_pt[j*nvec_y+jj+jjj];)
    x_pt += n2*nvec_x;
    y_pt += n2*nvec_y;
    A += oo_inv_size;
  }
}
  
#endif
