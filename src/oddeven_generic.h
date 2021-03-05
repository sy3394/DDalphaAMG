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
 * 12/08/2019: some reordering might be needed
 * 1st cleanup: 12/22/2019
 */
 
#ifndef ODDEVEN_PRECISION_HEADER
  #define ODDEVEN_PRECISION_HEADER

  struct Thread;

  void oddeven_setup_PRECISION( operator_double_struct *in, level_struct *l );
  void oddeven_free_PRECISION( level_struct *l );

  void selfcoupling_cholesky_decomposition_PRECISION( const config_PRECISION output, config_double input );
#ifdef HAVE_TM
  void selfcoupling_LU_decomposition_PRECISION( const config_PRECISION output, config_double input, level_struct *l );
#endif
#ifdef HAVE_TM1p1
void selfcoupling_LU_doublet_decomposition_PRECISION( const config_PRECISION output, config_double input, level_struct *l );
#endif

  void hopping_term_PRECISION( vector_PRECISION *eta, vector_PRECISION *phi, operator_PRECISION_struct *op,
                               const int amount, level_struct *l, struct Thread *threading );
  void apply_schur_complement_PRECISION( vector_PRECISION *out, vector_PRECISION *in, operator_PRECISION_struct *op, level_struct *l, struct Thread *threading );

  void solve_oddeven_PRECISION( gmres_PRECISION_struct *p, operator_PRECISION_struct *op, level_struct *l, struct Thread *threading );
  
  void oddeven_to_serial_PRECISION( vector_double *out, vector_PRECISION *in, level_struct *l, struct Thread *threading );
  void serial_to_oddeven_PRECISION( vector_PRECISION *out, vector_double *in, level_struct *l, struct Thread *threading );
  
  void oddeven_PRECISION_test( level_struct *l );


  static inline void LU_perform_fwd_bwd_subs_PRECISION( vector_PRECISION *x, vector_PRECISION *b, config_PRECISION LU,
							int start, int end, level_struct *l ) {

    /*********************************************************************************
     * Solves L*U*x = b for x, i.e., the clover coupling for a single lattice 
     * site.
     * - vector_PRECISION *b: Right hand side.
     * - vector_PRECISION *x: Solution.
     * - config_PRECISION L: Lower matrix from modified LU decomposition
     * Note: U is given by u_{ii}=1, u_{ij}=l_{ji}* / l_{ii} 
     *********************************************************************************/
    
    register int s, i, j, n, jj, jjj, nvec = x->num_vect_now, nvec_x = x->num_vect, nvec_b = b->num_vect;
    buffer_PRECISION x_pt = x->vector_buffer, b_pt = b->vector_buffer;
    
    // TODO: when g.method>3 and HAVE_MULT_TM, nv should be changed.
#ifdef DEBUG
    if ( b->num_vect_now != nvec )
      error0("LU_perform_fwd_bwd_subs_PRECISION: assumptions are not met\n");
#endif
    
#ifdef HAVE_TM1p1
    if( g.n_flavours == 2) {
#ifdef DEBUG
      if ( nvec != 2*num_loop )
	error0("LU_perform_fwd_bwd_subs_PRECISION (doublet):assumptions are not met\n");
#endif
      int i2, j2;
#ifdef HAVE_MULT_TM
      int nc = g.n_chunk, nv = (l->num_inner_lattice_sites/2+1)*288*num_loop;
      //LU += start*288*num_loop;
      //    for ( id=start; id<end; id+=24 ) {
      for ( s=start; s<end; s++ ) {
	for ( n=0; n<2; n++ ) { // upper diagonal part and then lower diaognal part of self-coupling, both LU decomposed
	  // solve x = U^(-1) L^(-1) b
	  // forward substitution with L
	  /*
	    for ( i=0; i<12; i++ ) {
	    VECTOR_LOOP(jj, nvec, jjj, x_pt[i*nvec_x+jj+jjj] = b_pt[i*nvec_b+jj+jjj];)
	    for ( j=0; j<i; j++ ) {
            VECTOR_LOOP(jj, nvec/2, jjj, x_pt[(i%6)*nvec_x+(i/6)*nvex_x/2+jj+jjj] -= LU[(i*12+j)*num_loop+nc*nv+jj+jjj]*x_pt[(j%6)*nvec_x+(j/6)*nvec_x/2+jj+jjj];)
	    }
	    }*/
	  for ( i=0; i<6; i++ ) {
	    VECTOR_LOOP(jj, nvec, jjj, x_pt[i*nvec_x+jj+jjj] = b_pt[i*nvec_b+jj+jjj];)
	    for ( j=0; j<i; j++ ) 
	      VECTOR_LOOP(jj, nvec/2, jjj, x_pt[i*nvec_x+jj+jjj] -= LU[(i*12+j)*num_loop+nc*nv+jj+jjj]*x_pt[j*nvec_x+jj+jjj];)
	  }
	  for ( i=6; i<12; i++ ) {
	    i2 = i-6;
	    //VECTOR_LOOP(jj, nvec/2, jjj, x_pt[i2*nvec_x+nvec_x/2+jj+jjj] = b_pt[i2*nvec_b+nvec_b/2+jj+jjj];)
	    for ( j=0; j<6; j++ )
	      VECTOR_LOOP(jj, nvec/2, jjj, x_pt[i2*nvec_x+nvec_x/2+jj+jjj] -= LU[(i*12+j)*num_loop+nc*nv+jj+jjj]*x_pt[j*nvec_x+jj+jjj];)
	    for ( j=6; j<i; j++ ) {
	      j2 = j-6;
	      VECTOR_LOOP(jj, nvec/2, jjj, x_pt[i2*nvec_x+nvec_x/2+jj+jjj] -= LU[(i*12+j)*num_loop+nc*nv+jj+jjj]*x_pt[j2*nvec_x+nvec_x/2+jj+jjj];)
	    }
	  }
	  // backward substitution with U
	  for ( i=12-1; i>=6; i-- ) {
	    i2 = i-6;
	    for ( j=i+1; j<12; j++ ) {
	      j2 = j-6;
	      VECTOR_LOOP(jj, nvec/2, jjj, x_pt[i2*nvec_x+nvec_x/2+jj+jjj] -= LU[(i*12+j)*num_loop+nc*nv+jj+jjj]*x_pt[j2*nvec_x+nvec_x/2+jj+jjj];)
	    }
	    VECTOR_LOOP(jj, nvec/2, jjj, x_pt[i2*nvec_x+nvec_x/2+jj+jjj] /= LU[(i*(12+1))*num_loop+nc*nv+jj+jjj];)
	  }
	  for ( i=5; i>=0; i-- ) {
	    for ( j=i+1; j<6; j++ ) 
	      VECTOR_LOOP(jj, nvec/2, jjj, x_pt[i*nvec_x+jj+jjj] -= LU[(i*12+j)*num_loop+nc*nv+jj+jjj]*x_pt[j*nvec_x+jj+jjj];)
	    for ( j=6; j<12; j++ ) {
	      j2 = j-6;
	      VECTOR_LOOP(jj, nvec/2, jjj, x_pt[i*nvec_x+jj+jjj] -= LU[(i*12+j)*num_loop+nc*nv+jj+jjj]*x_pt[j2*nvec_x+nvec_x/2+jj+jjj];)
	    }
	    VECTOR_LOOP(jj, nvec/2, jjj, x_pt[i*nvec_x+jj+jjj] /= LU[(i*(12+1))*num_loop+nc*nv+jj+jjj];)
	  }
	  x_pt+=6*nvec_x;
	  b_pt+=6*nvec_b;
	  LU+=12*12*num_loop;
	}
      }
#else
      //LU += start*288;
      for ( s=start; s<end; s++ ) {
	for ( n=0; n<2; n++ ) { // upper diagonal part and then lower diaognal part of self-coupling, both LU decomposed
	  // solve x = U^(-1) L^(-1) b
	  
	  // forward substitution with L
	  for ( i=0; i<6; i++ ) {
	    VECTOR_LOOP(jj, nvec, jjj, x_pt[i*nvec_x+jj+jjj] = b_pt[i*nvec_b+jj+jjj];)
	   for ( j=0; j<i; j++ ) 
	     VECTOR_LOOP(jj, nvec/2, jjj, x_pt[i*nvec_x+jj+jjj] -= LU[i*12+j]*x_pt[j*nvec_x+jj+jjj];)
	   }
	  for ( i=6; i<12; i++ ) {
	    i2 = i-6;
	    //VECTOR_LOOP(jj, nvec, jjj, x_pt[i2*nvec_x+jj+jjj] = b_pt[i2*nvec_b+jj+jjj];)
	    for ( j=0; j<6; j++ ) 
	      VECTOR_LOOP(jj, nvec/2, jjj, x_pt[i2*nvec_x+nvec_x/2+jj+jjj] -= LU[i*12+j]*x_pt[j*nvec_x+jj+jjj];)
	    for ( j=6; j<i; j++ ) {
	      j2 = j-6;
	      VECTOR_LOOP(jj, nvec/2, jjj, x_pt[i2*nvec_x+nvec_x/2+jj+jjj] -= LU[i*12+j]*x_pt[j2*nvec_x+nvec_x/2+jj+jjj];)
	    }
	  }
	  // backward substitution with U
	  for ( i=12-1; i>=6; i-- ) {
	    i2 = i-6;
	    for ( j=i+1; j<12; j++ ) {
	      j2 = j-6;
	      VECTOR_LOOP(jj, nvec/2, jjj, x_pt[i2*nvec_x+nvec_x/2+jj+jjj] -= LU[i*12+j]*x_pt[j2*nvec_x+nvec_x/2+jj+jjj];)
	    }
	    VECTOR_LOOP(jj, nvec/2, jjj, x_pt[i2*nvec_x+nvec_x/2+jj+jjj] /= LU[i*(12+1)];)
	  }
	  for ( i=5; i>=0; i-- ) {
	    for ( j=i+1; j<6; j++ ) 
	      VECTOR_LOOP(jj, nvec/2, jjj, x_pt[i*nvec_x+jj+jjj] -= LU[i*12+j]*x_pt[j*nvec_x+jj+jjj];)
	    for ( j=6; j<12; j++ ) {
	      j2 = j-6;
	      VECTOR_LOOP(jj, nvec/2, jjj, x_pt[i*nvec_x+jj+jjj] -= LU[i*12+j]*x_pt[j2*nvec_x+nvec_x/2+jj+jjj];)
	    }
	    VECTOR_LOOP(jj, nvec/2, jjj, x_pt[i*nvec_x+jj+jjj] /= LU[i*(12+1)];)
	  }
	  x_pt+=6*nvec_x;
	  b_pt+=6*nvec_b;
	  LU+=12*12;
	}
      }
#endif
    } else
#endif
      {
#ifdef HAVE_MULT_TM
	int nc = g.n_chunk, nv = (l->num_inner_lattice_sites/2+1)*72;
	for ( s=start; s<end; s++ ) {//id+=12 ) {
	  for ( n=0; n<2; n++ ) {
	    // solve x = U^(-1) L^(-1) b
	    // forward substitution with L
	    for ( i=0; i<6; i++ ) {
	      VECTOR_LOOP(jj, nvec, jjj, x_pt[i*nvec_x+jj+jjj] = b_pt[i*nvec_b+jj+jjj];)
	      for ( j=0; j<i; j++ ) 
		VECTOR_LOOP(jj, nvec, jjj, x_pt[i*nvec_x+jj+jjj] -= LU[(i*6+j)*num_loop+nc*nv+jj+jjj]*x_pt[j*nvec_x+jj+jjj];)
	    }
	    // backward substitution with U
	    for ( i=6-1; i>=0; i-- ) {
	      for ( j=i+1; j<6; j++ )
		VECTOR_LOOP(jj, nvec, jjj, x_pt[i*nvec_x+jj+jjj] -= LU[(i*6+j)*num_loop+nc*nv+jj+jjj]*x_pt[j*nvec_x+jj+jjj];)
	      VECTOR_LOOP(jj, nvec, jjj,  x_pt[i*nvec_x+jj+jjj] /=LU[i*(6+1)*num_loop+nc*nv+jj+jjj];)
	    }
	    x_pt+=6*nvec_x;
	    b_pt+=6*nvec_b;
	    LU+=6*6*num_loop;
	  }
	}
#else
	for ( s=start; s<end; s++ ) {//id+=12 ) {
	  for ( n=0; n<2; n++ ) {
	    // solve x = U^(-1) L^(-1) b
	    // forward substitution with L
	    for ( i=0; i<6; i++ ) {
	      VECTOR_LOOP(jj, nvec, jjj, x_pt[i*nvec_x+jj+jjj] = b_pt[i*nvec_b+jj+jjj];)
	    for ( j=0; j<i; j++ ) 
	      VECTOR_LOOP(jj, nvec, jjj, x_pt[i*nvec_x+jj+jjj] -= LU[i*6+j]*x_pt[j*nvec_x+jj+jjj];)
	    }
	    // backward substitution with U
	    for ( i=6-1; i>=0; i-- ) {
	      for ( j=i+1; j<6; j++ ) 
		VECTOR_LOOP(jj, nvec, jjj, x_pt[i*nvec_x+jj+jjj] -= LU[i*6+j]*x_pt[j*nvec_x+jj+jjj];)
	      VECTOR_LOOP(jj, nvec, jjj,  x_pt[i*nvec_x+jj+jjj] /=LU[i*(6+1)];)
	    }
	    x_pt+=6*nvec_x;
	    b_pt+=6*nvec_b;
	    LU+=6*6;
	  }
	}
#endif
      }
  }

  static inline void LU_multiply_PRECISION( vector_PRECISION *y, vector_PRECISION *x, config_PRECISION LU,
					    int start, int end, level_struct *l ) {

    /*********************************************************************************
     * Applies the clover coupling term to a vector, by multiplying L^H 
     * and then L. 
     * - vector_PRECISION *x: Input vector.
     * - vector_PRECISION *y: Output vector.
     * - config_PRECISION LU: LU decomposition
     *********************************************************************************/

    register int s, i, j, n, jj, jjj, nvec = x->num_vect_now, nvec_x = x->num_vect, nvec_y = y->num_vect;
    buffer_PRECISION x_pt = x->vector_buffer, y_pt = y->vector_buffer;
    
    // TODO: when g.method>3 and HAVE_MULT_TM, nv should be changed.  
#ifdef DEBUG
    if ( nvec_y < nvec )
      error0("LU_multiply_PRECISION: assumptions are not met\n");
#endif
    
#ifdef HAVE_TM1p1
    if( g.n_flavours == 2) {
      int i2, j2;
#ifdef HAVE_MULT_TM
      //printf0("LU m %d %d %d\n",nvec,nvec_y,nvec_x);
      int nc = g.n_chunk, nv = (l->num_inner_lattice_sites/2+1)*288*num_loop;
      for ( s=start; s<end; s++ ) { 
	for ( n=0; n<2; n++ ) {// upper diagonal part and then lower diaognal part of self-coupling, both LU decomposed  
	  for ( i=0; i<6; i++ ) {
	    VECTOR_LOOP( jj, nvec/2, jjj, y_pt[i*nvec_y+jj+jjj] = LU[(i*(12+1))*num_loop+nc*nv+jj+jjj]*x_pt[i*nvec_x+jj+jjj];)
	    for ( j=i+1; j<6; j++ )
	      VECTOR_LOOP( jj, nvec/2, jjj, y_pt[i*nvec_y+jj+jjj] += LU[(i*12+j)*num_loop+nc*nv+jj+jjj]*x_pt[j*nvec_x+jj+jjj];)
	    for ( j=6; j<12; j++ ) {
	      j2 = j-6;
	      VECTOR_LOOP( jj, nvec/2, jjj, y_pt[i*nvec_y+jj+jjj] += LU[(i*12+j)*num_loop+nc*nv+jj+jjj]*x_pt[j2*nvec_x+nvec_x/2+jj+jjj];)
	    }
	  }
	  for ( i=6; i<12; i++ ) {
	    i2 = i-6;
	    VECTOR_LOOP( jj, nvec/2, jjj, y_pt[i2*nvec_y+nvec_y/2+jj+jjj] = LU[(i*(12+1))*num_loop+nc*nv+jj+jjj]*x_pt[i2*nvec_x+nvec_x/2+jj+jjj];)
	    for ( j=i+1; j<12; j++ ) {
	      j2 = j-6;
	      VECTOR_LOOP( jj, nvec/2, jjj, y_pt[i2*nvec_y+nvec_y/2+jj+jjj] += LU[(i*12+j)*num_loop+nc*nv+jj+jjj]*x_pt[j2*nvec_x+nvec_x/2+jj+jjj];)
	    }
	  }
	  // multiplication with L
	  for ( i=12-1; i>5; i-- ) {
	    i2 = i-6;
	    for ( j=0; j<6; j++ )
	      VECTOR_LOOP( jj, nvec/2, jjj, y_pt[i2*nvec_y+nvec_y/2+jj+jjj] += LU[(i*12+j)*num_loop+nc*nv+jj+jjj]*y_pt[j*nvec_y+jj+jjj];)
	    for ( j=6; j<i; j++ ) {
	      j2 = j-6;
	      VECTOR_LOOP( jj, nvec/2, jjj, y_pt[i2*nvec_y+nvec_y/2+jj+jjj] += LU[(i*12+j)*num_loop+nc*nv+jj+jjj]*y_pt[j2*nvec_y+nvec_y/2+jj+jjj];)
	    }
	  }
	  for ( i=5; i>0; i-- )
	    for ( j=0; j<i; j++ )
	      VECTOR_LOOP( jj, nvec/2, jjj, y_pt[i*nvec_y+jj+jjj] += LU[(i*12+j)*num_loop+nc*nv+jj+jjj]*y_pt[j*nvec_y+jj+jjj];)
	  x_pt+=6*nvec_x;
	  y_pt+=6*nvec_y;
	  LU+=12*12*num_loop;
	}
      }
#else
      for ( s=start; s<end; s++ ) { 
	for ( n=0; n<2; n++ ) {// upper diagonal part and then lower diaognal part of self-coupling, both LU decomposed  
	  for ( i=0; i<6; i++ ) {
	    VECTOR_LOOP( jj, nvec/2, jjj, y_pt[i*nvec_y+jj+jjj] = LU[i*(12+1)]*x_pt[i*nvec_x+jj+jjj];)
	    for ( j=i+1; j<6; j++ )
	      VECTOR_LOOP( jj, nvec/2, jjj, y_pt[i*nvec_y+jj+jjj] += LU[i*12+j]*x_pt[j*nvec_x+jj+jjj];)
	    for ( j=6; j<12; j++ ) {
	      j2 = j-6;
	      VECTOR_LOOP( jj, nvec/2, jjj, y_pt[i*nvec_y+jj+jjj] += LU[i*12+j]*x_pt[j2*nvec_x+nvec_x/2+jj+jjj];)
	    }
	  }
	  for ( i=6; i<12; i++ ) {
	    i2 = i-6;
	    VECTOR_LOOP( jj, nvec/2, jjj, y_pt[i2*nvec_y+nvec_y/2+jj+jjj] = LU[i*(12+1)]*x_pt[i2*nvec_x+nvec_x/2+jj+jjj];)
	    for ( j=i+1; j<12; j++ ) {
	      j2 = j-6;
	      VECTOR_LOOP( jj, nvec/2, jjj, y_pt[i2*nvec_y+nvec_y/2+jj+jjj] += LU[i*12+j]*x_pt[j2*nvec_x+nvec_x/2+jj+jjj];)
	    }
	  }
	  // multiplication with L
	  for ( i=12-1; i>5; i-- ) {
	    i2 = i-6;
	    for ( j=0; j<6; j++ )
	      VECTOR_LOOP( jj, nvec/2, jjj, y_pt[i2*nvec_y+nvec_y/2+jj+jjj] += LU[i*12+j]*y_pt[j*nvec_y+jj+jjj];)
	    for ( j=6; j<i; j++ ) {
	      j2 = j-6;
	      VECTOR_LOOP( jj, nvec/2, jjj, y_pt[i2*nvec_y+nvec_y/2+jj+jjj] += LU[i*12+j]*y_pt[j2*nvec_y+nvec_y/2+jj+jjj];)
	    }
	  for ( i=5; i>0; i-- )
	    for ( j=0; j<i; j++ )
	      VECTOR_LOOP( jj, nvec/2, jjj, y_pt[i*nvec_y+jj+jjj] += LU[i*12+j]*y_pt[j*nvec_y+jj+jjj];)
	  }
	  x_pt+=6*nvec_x;
	  y_pt+=6*nvec_y;
	  LU+=12*12;
	}
      }
#endif
    } else
#endif
      {
#ifdef HAVE_MULT_TM
      int nc = g.n_chunk, nv = (l->num_inner_lattice_sites/2+1)*72;
      for ( s=start; s<end; s++ ) {//id+=12 ) {
	for ( n=0; n<2; n++ ) {
	  for ( i=0; i<6; i++ ) {
            VECTOR_LOOP( jj, nvec, jjj, y_pt[i*nvec_y+jj+jjj] = LU[i*(6+1)*num_loop+nc*nv+jj+jjj]*x_pt[i*nvec_x+jj+jjj];)
	    for ( j=i+1; j<6; j++ )
              VECTOR_LOOP( jj, nvec, jjj, y_pt[i*nvec_y+jj+jjj] += LU[(i*6+j)*num_loop+nc*nv+jj+jjj]*x_pt[j*nvec_x+jj+jjj];)
	  }
	  // multiplication with L
	  for ( i=6-1; i>0; i-- )
	    for ( j=0; j<i; j++ )
	      VECTOR_LOOP(jj, nvec, jjj, y_pt[i*nvec_y+jj+jjj] += LU[(i*6+j)*num_loop+nc*nv+jj+jjj]*y_pt[j*nvec_y+jj+jjj];)

	  x_pt+=6*nvec_x;
	  y_pt+=6*nvec_y;
	  LU+=6*6*num_loop;
	}
      }
#else
      for ( s=start; s<end; s++ ) {//id+=12 ) {
	for ( n=0; n<2; n++ ) {
	  for ( i=0; i<6; i++ ) {
            VECTOR_LOOP( jj, nvec, jjj, y_pt[i*nvec_y+jj+jjj] = LU[i*(6+1)]*x_pt[i*nvec_x+jj+jjj];)
	    for ( j=i+1; j<6; j++ )
              VECTOR_LOOP( jj, nvec, jjj, y_pt[i*nvec_y+jj+jjj] += LU[i*6+j]*x_pt[j*nvec_x+jj+jjj];)
	  }
	  // multiplication with L
	  for ( i=6-1; i>0; i-- )
	    for ( j=0; j<i; j++ )
	      VECTOR_LOOP(jj, nvec, jjj, y_pt[i*nvec_y+jj+jjj] += LU[i*6+j]*y_pt[j*nvec_y+jj+jjj];)

	  x_pt+=6*nvec_x;
	  y_pt+=6*nvec_y;
	  LU+=6*6;
	}
      }
#endif
    }
}

  static inline void LLH_perform_fwd_bwd_subs_PRECISION( vector_PRECISION *x, vector_PRECISION *b, config_PRECISION L,
							 int start, int end ) {

    /*********************************************************************************
     * Solves L*(L^H)*x = b for x, i.e., the clover coupling for a single lattice 
     * site.
     * - vector_PRECISION *b: Right hand side.
     * - vector_PRECISION *x: Solution.
     * - config_PRECISION L: Cholesky factor ( lower triangular matrix )
     *********************************************************************************/

    register int s, i, j;
    int n, jj, jjj, nvec = x->num_vect_now, nvec_x = x->num_vect, nvec_b = b->num_vect;
    buffer_PRECISION x_pt = x->vector_buffer, b_pt = b->vector_buffer;

#ifdef DEBUG
    if ( b->num_vect_now != nvec )
      error0("LLH_perform_fwd_bwd_subs_PRECISION: assumptions are not met\n");
#endif
  
    for ( s=start; s<end; s++ ) {//id+=12 ) {
      for ( n=0; n<2; n++ ) {
	// forward substitution with L
	for ( i=0; i<6; i++ ) {
	  VECTOR_LOOP(jj, nvec, jjj, x_pt[i*nvec_x+jj+jjj] = b_pt[i*nvec_b+jj+jjj];)
	    for ( j=0; j<i; j++ ) {
	      VECTOR_LOOP(jj, nvec, jjj, x_pt[i*nvec_x+jj+jjj] -= *L * x_pt[j*nvec_x+jj+jjj];)
	      L++;
	    }
	  VECTOR_LOOP(jj, nvec, jjj, x_pt[i*nvec_x+jj+jjj] /= *L;)
	  L++;
	}
	L -= 21;
	// backward substitution with L^H
	for ( i=5; i>=0; i-- ) {
	  for ( j=i+1; j<6; j++ ) 
	    VECTOR_LOOP(jj, nvec, jjj, x_pt[i*nvec_x+jj+jjj] -= conj_PRECISION(L[(j*(j+1))/2 + i]) * x_pt[j*nvec_x+jj+jjj];)
	    VECTOR_LOOP(jj, nvec, jjj, x_pt[i*nvec_x+jj+jjj] /= conj_PRECISION(L[(i*(i+1))/2 + i]);)
	  }
	x_pt+=6*nvec_x;
	b_pt+=6*nvec_b;
	L+=21;
      }
    }
  }


  static inline void LLH_multiply_PRECISION( vector_PRECISION *y, vector_PRECISION *x, config_PRECISION L,
					     int start, int end ) {
    
    /*********************************************************************************
     * Applies the clover coupling term to a vector, by multiplying L^H 
     * and then L. 
     * - vector_PRECISION *x: Input vector.
     * - vector_PRECISION *y: Output vector.
     * - config_PRECISION L: Cholesky factor ( lower triangular matrix )
     *********************************************************************************/
    
    register int s, i, j, jj, jjj;
    int n, nvec = x->num_vect_now, nvec_x = x->num_vect, nvec_y =y->num_vect;
    complex_PRECISION z[6*nvec];
    buffer_PRECISION x_pt = x->vector_buffer, y_pt = y->vector_buffer;

    if ( nvec_y < nvec )
      error0("LLH_multiply_PRECISION: assumptions are not met\n");
    
    for ( s=start; s<end; s++ ) { //id+=12
      for ( n=0; n<2; n++ ) {
	// z = L^H x
	for ( j=0; j<6; j++ ) { // columns
	  for ( i=0; i<j; i++ ) { // rows
	    VECTOR_LOOP(jj, nvec, jjj, z[i*nvec+jj+jjj] += conj_PRECISION(*L)*x_pt[j*nvec_x+jj+jjj];)
	    L++;
	  }
	  VECTOR_LOOP(jj, nvec, jjj, z[j*nvec+jj+jjj] = conj_PRECISION(*L)*x_pt[j*nvec_x+jj+jjj];)
	  L++;
	}
	L-=21;
	// y = L*z;
	for ( i=0; i<6; i++ ) { // rows
	  VECTOR_LOOP(jj, nvec, jjj, y_pt[i*nvec_y+jj+jjj] = *L * z[0*nvec+jj+jjj];)
	  L++;
	  for ( j=1; j<=i; j++ ) { // columns
	    VECTOR_LOOP(jj, nvec, jjj, y_pt[i*nvec_y+jj+jjj] += *L * z[j*nvec+jj+jjj];)
	    L++;
	  }
	}
	x_pt+=6*nvec_x;
	y_pt+=6*nvec_y;
      }
    }
  }

#endif
