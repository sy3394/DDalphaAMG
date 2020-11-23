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


  void coarse_diag_ee_PRECISION( vector_PRECISION *y, vector_PRECISION *x, operator_PRECISION_struct *op, level_struct *l, struct Thread *threading );
  void coarse_diag_oo_PRECISION( vector_PRECISION *y, vector_PRECISION *x, operator_PRECISION_struct *op, level_struct *l, struct Thread *threading );
  int coarse_solve_odd_even_PRECISION( gmres_PRECISION_struct *p, operator_PRECISION_struct *op, level_struct *l, 
					    struct Thread *threading );
  void coarse_apply_schur_complement_PRECISION( vector_PRECISION *out, vector_PRECISION *in,
						    operator_PRECISION_struct *op, level_struct *l, struct Thread *threading );
  void coarse_hopping_term_PRECISION( vector_PRECISION *out, vector_PRECISION *in, operator_PRECISION_struct *op,
                                      const int amount, level_struct *l, struct Thread *threading );
  void coarse_n_hopping_term_PRECISION( vector_PRECISION *out, vector_PRECISION *in, operator_PRECISION_struct *op,
                                        const int amount, level_struct *l, struct Thread *threading );

  void coarse_odd_even_PRECISION_test( vector_PRECISION *c4, vector_PRECISION *c1,
                                       level_struct *l, struct Thread *threading );


//static inline void coarse_perform_fwd_bwd_subs_PRECISION( vector_PRECISION *x, vector_PRECISION *b, operator_PRECISION_struct *op, int start, int end, level_struct *l ) {
//static inline void coarse_perform_fwd_bwd_subs_PRECISION( vector_PRECISION *x, vector_PRECISION *b, config_PRECISION A, int start, int end, level_struct *l ) {
static inline void coarse_perform_fwd_bwd_subs_PRECISION( const buffer_PRECISION x_pt, const buffer_PRECISION b_pt, config_PRECISION A, 
							  const int nvec, const int nvec_b, const int nvec_x, level_struct *l ) {
  register int s, i, j, jj, jjj;
  int n2 = l->num_lattice_site_var;
  int oo_inv_size = SQUARE(n2);
#if 0
  int nvec = b->num_vect_now, nvec_x = x->num_vect, nvec_b = b->num_vect;
  config_PRECISION A = op->clover_oo_inv;
  buffer_PRECISION x_pt = x->vector_buffer+(op->num_even_sites+start)*n2*nvec_x, b_pt = b->vector_buffer+(op->num_even_sites+start)*n2*nvec_b;
#endif
#ifdef DEBUG
  if( nvec_x < nvec && nvec==num_loop )
    error0("coarse_perform_fwd_bwd_subs_PRECISION: assumptions are not met\n");
#endif

#if defined(HAVE_MULT_TM)
  //  A += oo_inv_size*start*num_loop;
  int nc = g.n_chunk, nv = oo_inv_size*(l->num_inner_lattice_sites/2)*num_loop;
  //  for ( s=start; s<end; s++ ) {  
    // solve x = U^(-1) L^(-1) b
    // forward substitution with L
    for ( i=0; i<n2; i++ ) {
      VECTOR_LOOP( jj, nvec, jjj, x_pt[i*nvec_x+jj+jjj] = b_pt[i*nvec_b+jj+jjj];)
	for ( j=0; j<i; j++ ) 
	  VECTOR_LOOP(jj, nvec, jjj, x_pt[i*nvec_x+jj+jjj] -= A[(i*n2+j)*num_loop+nc*nv+jj+jjj]*x_pt[j*nvec_x+jj+jjj];)
    }
    // backward substitution with U
    for ( i=n2-1; i>=0; i-- ) {
      for ( j=i+1; j<n2; j++ ) 
	VECTOR_LOOP(jj, nvec, jjj, x_pt[i*nvec_x+jj+jjj] -= A[(i*n2+j)*num_loop+nc*nv+jj+jjj]*x_pt[j*nvec_x+jj+jjj];)//printf("%g ",creal_PRECISION(A[i*n2+j]));)
      VECTOR_LOOP( jj, nvec, jjj, x_pt[i*nvec_x+jj+jjj] /= A[i*(n2+1)*num_loop+nc*nv+jj+jjj];)//printf("%g ",creal_PRECISION(A[i*n2+j]));)
    }
    /*
    x_pt += n2*nvec_x;
    b_pt += n2*nvec_b;
    A += oo_inv_size*num_loop;
  }*/
#else
  //  A += oo_inv_size*start;
    //  for ( s=start; s<end; s++ ) {  
    // solve x = U^(-1) L^(-1) b
    // forward substitution with L
    for ( i=0; i<n2; i++ ) {
      VECTOR_LOOP( jj, nvec, jjj, x_pt[i*nvec_x+jj+jjj] = b_pt[i*nvec_b+jj+jjj];)
	for ( j=0; j<i; j++ ) 
	  VECTOR_LOOP(jj, nvec, jjj, x_pt[i*nvec_x+jj+jjj] -= A[i*n2+j]*x_pt[j*nvec_x+jj+jjj];)
    }
    // backward substitution with U
    for ( i=n2-1; i>=0; i-- ) {
      for ( j=i+1; j<n2; j++ ) 
	VECTOR_LOOP(jj, nvec, jjj, x_pt[i*nvec_x+jj+jjj] -= A[i*n2+j]*x_pt[j*nvec_x+jj+jjj];)//printf("%g ",creal_PRECISION(A[i*n2+j]));)
      VECTOR_LOOP( jj, nvec, jjj, x_pt[i*nvec_x+jj+jjj] /= A[i*(n2+1)];)//printf("%g ",creal_PRECISION(A[i*n2+j]));)
    }
    /*
    x_pt += n2*nvec_x;
    b_pt += n2*nvec_b;
    A += oo_inv_size;
    }*/
#endif
}

// used only in test routines
//static inline void coarse_LU_multiply_PRECISION( vector_PRECISION *y, vector_PRECISION *x, operator_PRECISION_struct *op, int start, int  end, level_struct *l ) {
//static inline void coarse_LU_multiply_PRECISION( vector_PRECISION *y, vector_PRECISION *x, config_PRECISION A, int start, int  end, level_struct *l ) {
static inline void coarse_LU_multiply_PRECISION( const buffer_PRECISION y_pt, const buffer_PRECISION x_pt, const config_PRECISION A,
						 const int nvec, const int nvec_y, const int nvec_x, level_struct *l ) {
  register int i, j, jj, jjj;
  int n2 = l->num_lattice_site_var, oo_inv_size = SQUARE(n2);

#if 0
  int s, nvec = x->num_vect_now, nvec_x  = x->num_vect, nvec_y = y->num_vect;
  config_PRECISION A = op->clover_oo_inv;
  buffer_PRECISION x_pt = x->vector_buffer+(op->num_even_sites+start)*n2*nvec_x, y_pt = y->vector_buffer+(op->num_even_sites+start)*n2*nvec_y;
#endif
#ifdef DEBUG
  if ( nvec_y < nvec )
    error0("coarse_LU_multiply_PRECISION: assumptions are not met\n");
#endif

#if defined(HAVE_MULT_TM)
  //  A += oo_inv_size*start*num_loop;
    int nc = g.n_chunk, nv = oo_inv_size*(l->num_inner_lattice_sites/2)*num_loop;//#odd sits!!!!
  //  for ( s=start; s<end; s++ ) { 
    // y = Ax
    // multiplication with U
    for ( i=0; i<n2; i++ ) {
      VECTOR_LOOP( jj, nvec, jjj, y_pt[i*nvec_y+jj+jjj] = A[i*(n2+1)*num_loop+nc*nv+jj+jjj]*x_pt[i*nvec_x+jj+jjj];)
	for ( j=i+1; j<n2; j++ )
	  VECTOR_LOOP( jj, nvec, jjj, y_pt[i*nvec_y+jj+jjj] += A[(i*n2+j)*num_loop+nc*nv+jj+jjj]*x_pt[j*nvec_x+jj+jjj];)
    }
    // multiplication with L
    for ( i=n2-1; i>0; i-- )
      for ( j=0; j<i; j++ )
	VECTOR_LOOP(jj, nvec, jjj, y_pt[i*nvec_y+jj+jjj] += A[(i*n2+j)*num_loop+nc*nv+jj+jjj]*y_pt[j*nvec_y+jj+jjj];)
	  /*
    x_pt += n2*nvec_x;
    y_pt += n2*nvec_y;
    A += oo_inv_size*num_loop;
      }*/
#else
  //  A += oo_inv_size*start;
    //  for ( s=start; s<end; s++ ) { 
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
	  /*
    x_pt += n2*nvec_x;
    y_pt += n2*nvec_y;
    A += oo_inv_size;
    }*/
#endif
}
  
#endif
