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
 * some work undone: 12/06/2019
 * some inline functions are defined here instead of the header file???
 * 12/08/2019/ some reordering might be needed
 * 1st cleanup:12/22/2019
 */
 
#include "main.h"

static void diag_ee_PRECISION( vector_PRECISION *y, vector_PRECISION *x, operator_PRECISION_struct *op, level_struct *l, int start, int end );
static void diag_oo_PRECISION( vector_PRECISION *y, vector_PRECISION *x, operator_PRECISION_struct *op, level_struct *l, struct Thread *threading );
static void diag_oo_inv_PRECISION( vector_PRECISION *y, vector_PRECISION *x, operator_PRECISION_struct *op, level_struct *l, int start, int end );

// TODO: In this function, clover contains all the self-couling term including tm term and LU decomposed
//       For now, if HAVE_MULT_TM, we compute LU decomposition for each shift.
//                Otherwise, tm_term on even sites are computed using the avg shift, which is not good at the top level.
void oddeven_setup_PRECISION( operator_double_struct *in, level_struct *l ) {// called only when g.method==4,5 at the top in init.c to setup g.op_PRECISION

/*********************************************************************************
* Reorder data layouts and index tables to allow for odd even preconditioning.
*********************************************************************************/ 
  error0("not tested!\n");
  int j, k, k_e, k_o, n=l->num_inner_lattice_sites, oe_offset=0, mu, nu, nvec,
    sc_size = g.csw ? 42:12, lu_dec_size = 42, bs, **bt = NULL,
    *eot = NULL, *nt = NULL, *tt = NULL, t, z, y, x, le[4], N[4];
  config_double sc_in = in->clover, nc_in = in->D;
  config_PRECISION Aee = NULL, Aoo = NULL;
  operator_PRECISION_struct *op = &(l->oe_op_PRECISION);

  op->m0 = in->m0;

#ifdef HAVE_TM
  int jj, jjj, nv = l->num_inner_lattice_sites*12;

  op->mu = in->mu;
  op->mu_odd_shift = in->mu_odd_shift;
  op->odd_shifted_mu = in->odd_shifted_mu;
  op->even_shift_avg = in->even_shift_avg;
  for ( j=0; j<g.num_rhs_vect; j++ ) {
    op->mu_even_shift[j] = in->mu_even_shift[j];
    op->diff_mu_eo[j] = in->diff_mu_eo[j];
    op->is_even_shifted_mu_nonzero = in->is_even_shifted_mu_nonzero;
  }
#ifdef HAVE_MULT_TM
  lu_dec_size = 72*g.num_rhs_vect;
  config_double tm_term_in = in->tm_term;
#else
  lu_dec_size = 72;
  complex_double even = I*(in->mu+in->even_shift_avg), odd = I*(in->mu+in->mu_odd_shift);
#endif
#endif
  
  for ( mu=0; mu<4; mu++ ) {
    le[mu] = l->local_lattice[mu];
    N[mu] = le[mu]+1;
    op->table_dim[mu] = N[mu];
  }

  // locally, the first site on the process is always an even site as it is at (0,0,0,0)
  // However, the parity of the first local site is not necessarily even when put back to the global lattice
  for ( mu=0; mu<4; mu++ )
    oe_offset += (l->local_lattice[mu]*(g.my_coords[mu]/l->comm_offset[mu]))%2;
  oe_offset = oe_offset%2;
  
  // calculate site numbers
  op->num_even_sites = 0;
  op->num_odd_sites = 0;
  op->oe_offset = oe_offset;
  
  for ( t=0; t<le[T]; t++ )
    for ( z=0; z<le[Z]; z++ )
      for ( y=0; y<le[Y]; y++ )
        for ( x=0; x<le[X]; x++ ) {
          if ( (t+z+y+x+oe_offset)%2 == 1 ) {
            op->num_odd_sites++;
          } else {
            op->num_even_sites++;
          }
        }
  
  // re-order clover term (i.e., self coupling)
  if ( g.csw ) {
    MALLOC( op->clover, complex_PRECISION, lu_dec_size*n );
    Aee = op->clover;
    Aoo = op->clover + op->num_even_sites*lu_dec_size;
    for ( t=0; t<le[T]; t++ )
      for ( z=0; z<le[Z]; z++ )
        for ( y=0; y<le[Y]; y++ )
          for ( x=0; x<le[X]; x++ ) {
            if ( (t+z+y+x+oe_offset)%2 == 1 ) {
              // odd sites
#ifndef HAVE_TM
              selfcoupling_cholesky_decomposition_PRECISION( Aoo, sc_in );
#else
#ifdef HAVE_MULT_TM
              complex_double buffer[12*g.num_rhs_vect+30];
              for(int i=0; i<12; i++)
                VECTOR_LOOP(jj, g.num_rhs_vect, jjj, buffer[i*g.num_rhs_vect+jj+jjj] = sc_in[i];)
	      for(int i=12; i<42; i++)
		buffer[12*g.num_rhs_vect+i] = sc_in[i];
	      for(int i=0; i<12; i++ )
		VECTOR_LOOP(jj, g.num_rhs_vect, jjj, buffer[i*g.num_rhs_vect+jj+jjj] += (complex_double) tm_term_in[i*num_loop+jj*nv+jjj];)
		selfcoupling_LU_decomposition_PRECISION( Aoo, buffer, l );
#else
	      complex_double buffer[42];
              for(int i=0; i<42; i++)
                buffer[i] = sc_in[i];
              for(int i=0; i<6; i++ )
                buffer[i] -= odd;
	      for(int i=6; i<12; i++ )
                buffer[i] += odd;
              selfcoupling_LU_decomposition_PRECISION( Aoo, buffer, l );
#endif
#endif
              Aoo += lu_dec_size;
            } else {
              // even sites
#ifndef HAVE_TM
              selfcoupling_cholesky_decomposition_PRECISION( Aee, sc_in );
#else
#ifdef HAVE_MULT_TM 
              complex_double buffer[12*g.num_rhs_vect+30];
              for(int i=0; i<12; i++)
                VECTOR_LOOP(jj, g.num_rhs_vect, jjj, buffer[i*g.num_rhs_vect+jj+jjj] = sc_in[i];)
              for(int i=12; i<42; i++)
                buffer[12*g.num_rhs_vect+i] = sc_in[i];
              for(int i=0; i<12; i++ )
                VECTOR_LOOP(jj, g.num_rhs_vect, jjj, buffer[i*g.num_rhs_vect+jj+jjj] += (complex_double) tm_term_in[i*num_loop+jj*nv+jjj];)
		selfcoupling_LU_decomposition_PRECISION( Aee, buffer, l );
#else
              complex_double buffer[42];
              for(int i=0; i<42; i++)
                buffer[i] = sc_in[i];
              for(int i=0; i<6; i++ )
                buffer[i] -= even;
	      for(int i=6; i<12; i++ )
                buffer[i] += even;
              selfcoupling_LU_decomposition_PRECISION( Aoo, buffer, l );
#endif
#endif
              Aee += lu_dec_size;
            }
            sc_in += sc_size;
#ifdef HAVE_MULT_TM
            tm_term_in += 12*num_loop;
#endif
          }
  } else {//no idea why tm_term is not taken into account?????
    MALLOC( op->clover, complex_PRECISION, 12*n );
    Aee = op->clover;
    Aoo = op->clover + op->num_even_sites*12;
    
    for ( t=0; t<le[T]; t++ )
      for ( z=0; z<le[Z]; z++ )
        for ( y=0; y<le[Y]; y++ )
          for ( x=0; x<le[X]; x++ ) {
            if ( (t+z+y+x+oe_offset)%2 == 1 ) {
              // odd sites
              FOR12( *Aoo = *sc_in; Aoo++; sc_in++; )
            } else {
              // even sites
              FOR12( *Aee = *sc_in; Aee++; sc_in++; )
            }
          }
  }

#ifdef HAVE_TM1p1

  int doublet_nvec = 1;
  int lu_doublet_dec_size = 288;
  config_double eps_term_in = in->epsbar_term;
  sc_in = in->clover; 
#ifdef HAVE_MULT_TM
  tm_term_in = in->tm_term;
  lu_doublet_dec_size *= num_loop;
#endif
  op->epsbar = in->epsbar;
  op->epsbar_ig5_even_shift = in->epsbar_ig5_even_shift;
  op->epsbar_ig5_odd_shift = in->epsbar_ig5_odd_shift;
    
  // re-order clover term (i.e., self coupling)
#ifdef HAVE_MULT_TM
  doublet_nvec *= g.num_rhs_vect;
#endif
  MALLOC( op->clover_doublet_oo_inv, complex_PRECISION, lu_doublet_dec_size*n*doublet_nvec );
  
  Aee = op->clover_doublet_oo_inv;
  Aoo = op->clover_doublet_oo_inv + op->num_even_sites*lu_doublet_dec_size*doublet_nvec;
  for ( t=0; t<le[T]; t++ )
    for ( z=0; z<le[Z]; z++ )
      for ( y=0; y<le[Y]; y++ )
        for ( x=0; x<le[X]; x++ ) {
          if ( (t+z+y+x+oe_offset)%2 == 1 ) {
            // odd sites
#ifdef HAVE_MULT_TM
	    complex_double buffer[12*g.num_rhs_vect+54];
	    for(int i=0; i<12; i++) //0-23
	      VECTOR_LOOP(jj, g.num_rhs_vect, jjj,buffer[(i+12)*g.num_rhs_vect+jj+jjj] = buffer[i*g.num_rhs_vect+jj+jjj] = (complex_double) sc_in[i];)
	    if ( g.csw )
              for(int i=0; i<30; i++) //24-53
                buffer[i+24*g.num_rhs_vect] = (complex_double) sc_in[i+12];
            else
              for(int i=0; i<30; i++) //24-53
                buffer[i+24*g.num_rhs_vect] = _COMPLEX_double_ZERO;
	    for(int i=0; i<12; i++) //54-65
              buffer[i+24*g.num_rhs_vect+30] = (complex_double) eps_term_in[i];
	    if (g.mu + g.mu_odd_shift != 0.0 || g.is_even_shifted_mu_nonzero )
              for(int i=0; i<12; i++) {
                VECTOR_LOOP(jj, g.num_rhs_vect, jjj,buffer[i*g.num_rhs_vect+jj+jjj] += (complex_double) tm_term_in[i*g.num_rhs_vect+jj+jjj];)
		VECTOR_LOOP(jj, g.num_rhs_vect, jjj,buffer[(i+12)*g.num_rhs_vect+jj+jjj] -= (complex_double) tm_term_in[i*g.num_rhs_vect+jj+jjj];)
              }
#else
            complex_double buffer[66]; // 42+12*2 = 66
            if ( g.csw ) {
              for(int i=0; i<12; i++) //0-23
                buffer[i+12] = buffer[i] = (complex_double) sc_in[i];
              for(int i=12; i<42; i++) //24-53
                buffer[i+12] = (complex_double) sc_in[i];
            } else {
              for(int i=0; i<12; i++) //0-23
                buffer[i+12] = buffer[i] = (complex_double) sc_in[i];
              for(int i=12; i<42; i++) //24-53
                buffer[i+12] = _COMPLEX_double_ZERO;
            }              
            for(int i=0; i<12; i++) //54-65
              buffer[i+54] = (complex_double) eps_term_in[i];
#ifdef HAVE_TM
            if (g.mu + g.mu_odd_shift != 0.0 || g.is_even_shifted_mu_nonzero )
              for(int i=0; i<12; i++) {
                buffer[i] += odd;
                buffer[i+12] -= odd;
              }
#endif
#endif
            selfcoupling_LU_doublet_decomposition_PRECISION( Aoo, buffer, l );
            Aoo += lu_doublet_dec_size;
          } else {
            // even sites
#ifdef HAVE_MULT_TM
	    complex_double buffer[12*g.num_rhs_vect+54];
	    for(int i=0; i<12; i++) //0-23
	      VECTOR_LOOP(jj, g.num_rhs_vect, jjj,buffer[(i+12)*g.num_rhs_vect+jj+jjj] = buffer[i*g.num_rhs_vect+jj+jjj] = (complex_double) sc_in[i];)
	    if ( g.csw )
              for(int i=0; i<30; i++) //24-53
                buffer[i+24*g.num_rhs_vect] = (complex_double) sc_in[i];
            else
              for(int i=0; i<30; i++) //24-53
                buffer[i+24*g.num_rhs_vect] = _COMPLEX_double_ZERO;
	    for(int i=0; i<12; i++) //54-65
              buffer[i+24*g.num_rhs_vect+30] = (complex_double) eps_term_in[i];
	    if (g.mu + g.mu_odd_shift != 0.0 || g.is_even_shifted_mu_nonzero )
              for(int i=0; i<12; i++) {
                VECTOR_LOOP(jj, g.num_rhs_vect, jjj,buffer[i*g.num_rhs_vect+jj+jjj] += (complex_double) tm_term_in[i*g.num_rhs_vect+jj+jjj];)
		VECTOR_LOOP(jj, g.num_rhs_vect, jjj,buffer[(i+12)*g.num_rhs_vect+jj+jjj] -= (complex_double) tm_term_in[i*g.num_rhs_vect+jj+jjj];)
              }
#else
            complex_double buffer[66]; // 42+12*2 = 66
            if ( g.csw ) {
              for(int i=0; i<12; i++) //0-23
                buffer[i+12] = buffer[i] = (complex_double) sc_in[i];
              for(int i=12; i<42; i++) //24-53
                buffer[i+12] = (complex_double) sc_in[i];
            } else {
              for(int i=0; i<12; i++) //0-23
                buffer[i+12] = buffer[i] = (complex_double) sc_in[i];
              for(int i=12; i<42; i++) //24-53
                buffer[i+12] = _COMPLEX_double_ZERO;
            }              
            for(int i=0; i<12; i++) //54-65
              buffer[i+54] = (complex_double) eps_term_in[i];
#ifdef HAVE_TM
            if (g.mu + g.mu_odd_shift != 0.0 || g.is_even_shifted_mu_nonzero )
              for(int i=0; i<12; i++) {
                buffer[i] += even;
                buffer[i+12] -= even;
              }
#endif
#endif
            selfcoupling_LU_doublet_decomposition_PRECISION( Aee, buffer, l );
            Aee += lu_doublet_dec_size;
          }
          sc_in += sc_size;
          eps_term_in += 12;
#ifdef HAVE_MULT_TM
          tm_term_in += 12;
#endif
        }
#endif
  
  // re-order hopping term (i.e., nearest neighbor coupling)
  MALLOC( op->D, complex_PRECISION, 36*n );
  
  k=0; k_e=0; k_o=0;
  for ( t=0; t<le[T]; t++ )
    for ( z=0; z<le[Z]; z++ )
      for ( y=0; y<le[Y]; y++ )
        for ( x=0; x<le[X]; x++ ) {
          if ( (t+z+y+x+oe_offset)%2 == 1 ) {
            for ( j=0; j<36; j++)
              (op->D+(k_e+op->num_even_sites)*36)[j] = (complex_PRECISION) (nc_in+k*36)[j];
            k_e++;
          } else {
            for ( j=0; j<36; j++)
              (op->D+k_o*36)[j] = (complex_PRECISION) (nc_in+k*36)[j];
            k_o++;
          }
          k++;
        }
       
  // define data layout
  MALLOC( op->index_table, int, N[T]*N[Z]*N[Y]*N[X] );
  eot = op->index_table;
  
  define_eot( eot, N, l );
  
  // neighbor table, translation table
  MALLOC( op->neighbor_table, int, 5*N[T]*N[Z]*N[Y]*N[X] );
  MALLOC( op->backward_neighbor_table, int, 5*N[T]*N[Z]*N[Y]*N[X] );
  MALLOC( op->translation_table, int, le[T]*le[Z]*le[Y]*le[X] );
  nt = op->neighbor_table;
  tt = op->translation_table;
  
  define_nt_bt_tt( nt, op->backward_neighbor_table, NULL, tt, eot, N, l );
  
  // boundary table
  for ( mu=0; mu<4; mu++ ) {
    bs = 1;
    le[mu] = 1;
    for ( nu=0; nu<4; nu++ )
      bs *= le[nu];
    
    MALLOC( op->c.boundary_table[2*mu], int, bs );
    op->c.boundary_table[2*mu+1] = op->c.boundary_table[2*mu];
    
    le[mu] = l->local_lattice[mu];
  }
  
  bt = op->c.boundary_table;
  define_eo_bt( bt, eot, op->c.num_even_boundary_sites, op->c.num_odd_boundary_sites, op->c.num_boundary_sites, N, l );
  
  nvec = num_loop;
#ifdef HAVE_TM1p1
  nvec *= 2;
#endif
  op->pr_num_vect = nvec;
  j = (l->num_lattice_site_var/2)*l->num_lattice_sites*nvec;
  MALLOC( op->prnT, complex_PRECISION, j*8 );
  op->prnZ = op->prnT + j; op->prnY = op->prnZ + j; op->prnX = op->prnY + j;
  op->prpT = op->prnX + j; op->prpZ = op->prpT + j; op->prpY = op->prpZ + j; op->prpX = op->prpY + j;  
  MALLOC( op->buffer, vector_PRECISION, 2 );
  for(int i=0; i<2; i++ ){
    vector_PRECISION_init( &(op->buffer[i]) );
    vector_PRECISION_alloc( &(op->buffer[i]), _ORDINARY, nvec, l, no_threading );
  }
  ghost_alloc_PRECISION( 0, &(op->c), l );
  ghost_sendrecv_init_PRECISION( _COARSE_GLOBAL, &(op->c), l ) ;
  l->sp_PRECISION.v_end = op->num_even_sites*l->num_lattice_site_var;
}


void oddeven_free_PRECISION( level_struct *l ) {

  warning0("oddeven_free_PRECISION: not tested\n");
  
  int mu, nu, nc_size = 9, lu_dec_size = 42, nvec = l->oe_op_PRECISION.pr_num_vect,
      *ll = l->local_lattice, n = l->num_inner_lattice_sites, bs;
#ifdef HAVE_TM
#ifdef HAVE_MULT_TM
  lu_dec_size = 72*g.num_rhs_vect;
#else
  lu_dec_size = 72;
#endif
#endif
      
  ghost_free_PRECISION( &(l->oe_op_PRECISION.c), l );
  FREE( l->oe_op_PRECISION.D, complex_PRECISION, 4*nc_size*n );
  if ( g.csw )
    FREE( l->oe_op_PRECISION.clover, complex_PRECISION, lu_dec_size*n );
  else
    FREE( l->oe_op_PRECISION.clover, complex_PRECISION, 12*n );
  FREE( l->oe_op_PRECISION.index_table, int, (ll[T]+1)*(ll[Z]+1)*(ll[Y]+1)*(ll[X]+1) );
  FREE( l->oe_op_PRECISION.neighbor_table, int, 5*(ll[T]+1)*(ll[Z]+1)*(ll[Y]+1)*(ll[X]+1) );
  FREE( l->oe_op_PRECISION.backward_neighbor_table, int, 5*(ll[T]+1)*(ll[Z]+1)*(ll[Y]+1)*(ll[X]+1) );
  FREE( l->oe_op_PRECISION.translation_table, int, ll[T]*ll[Z]*ll[Y]*ll[X] );
  
  for ( mu=0; mu<4; mu++ ) {
    bs = 1;
    for ( nu=0; nu<4; nu++ )
      if ( mu != nu )
        bs *= ll[nu];
    
    FREE( l->oe_op_PRECISION.c.boundary_table[2*mu], int, bs );
    l->oe_op_PRECISION.c.boundary_table[2*mu+1] = NULL;
  }

  for(int i=0; i<2; i++ )
    vector_PRECISION_free( &(l->oe_op_PRECISION.buffer[i]), l, no_threading );
  FREE( l->oe_op_PRECISION.buffer, vector_PRECISION, 2 );
  
  FREE( l->oe_op_PRECISION.prnT, complex_PRECISION, (l->num_lattice_site_var/2)*l->num_lattice_sites*8*nvec ); //elimminated the factor of 2 for HAVE_TM1p1 case as included in nvec
#ifdef HAVE_TM1p1
#ifdef HAVE_MULT_TM
  FREE( l->oe_op_PRECISION.clover_doublet_oo_inv, complex_PRECISION, 288*n*g.num_rhs_vect );
#else
  FREE( l->oe_op_PRECISION.clover_doublet_oo_inv, complex_PRECISION, 288*n );
#endif
#endif
}

void selfcoupling_cholesky_decomposition_PRECISION( const config_PRECISION output, config_double input ) {
  
/*********************************************************************************
* Performs a cholesky decomposition for a selfcoupling term.
* Input = [ A 0 ]   , A=A*, B=B*
*         [ 0 B ]
* Output ordering: diag(A), diag(B), triu(A,1) row major, triu(B,1) row major 
*                  (matlab notation). 
* Note: clover D_W is gamma_5 symmetric => off-diagonal parts are zero
*       clover term is Hermitian => only upper triangular part is necessary
*        => 3*7*2 = 42 d.o.f.
*       Output is of the same nunmber of d.o.f. ordered differently
*********************************************************************************/  
  
  register int i, j, k;
  int n, offset[4] = {0,12,6,27};
  config_double in_pt;
  config_PRECISION out_pt = output;
  complex_PRECISION L[6][6], s;
  
  for ( n=0; n<2; n++ ) {
    // construct initial L = tril(A) for n=0, L = tril(B) for n=1, L row major
    in_pt = input+offset[2*n];
    for ( j=0; j<6; j++ ) {
      L[j][j] = (complex_PRECISION) *in_pt; in_pt++;
    }
    
    in_pt = input+offset[2*n+1];
    for ( j=0; j<5; j++ ) {
      for ( i=j+1; i<6; i++ ) {
        L[i][j] = (complex_PRECISION) conj_double(*in_pt); in_pt++;
      }
    }
    
    // calculate cholesky factor L
    for ( i=0; i<6; i++ ) {
      for ( j=0; j<=i; j++ ) {
        s = L[i][j];
        for ( k=0; k<j; k++ )
          s = s - L[i][k]*conj_PRECISION(L[j][k]);
        if ( i > j )
          L[i][j] = s / L[j][j];
        else if ( cabs_PRECISION(L[i][i]) > EPS_PRECISION ) 
          L[i][i] = csqrt_PRECISION(s);
        else {
          L[i][i] = csqrt_PRECISION(s);
        }
      }  
    }
    
    // output = tril(L) row major 
    for ( i=0; i<6; i++ ) {
      for ( j=0; j<=i; j++ ) {
        *out_pt = L[i][j]; out_pt++;
      }
    }
  }
}

#ifdef HAVE_TM
void selfcoupling_LU_decomposition_PRECISION( const config_PRECISION output, config_double input, level_struct *l ) {

  /*********************************************************************************   
   * Performs a LU decomposition for a selfcoupling term.
   * Input = [ A 0 ]   , A=A*, B=B* (diagonals excluded)
   *         [ 0 B ]               
   * Input ordering: diag(A), diag(B), triu(A,1) row major, triu(B,1) row major
   *                  (matlab notation).
   * Output ordering: triu(L,1) + tril(U,0), i.e., output contains L and U without 
   *                  the diagonal of L which is equal to 1 for A&B => (15+21)*2
   *********************************************************************************/

  register int i, j, k;
  int n;
  config_double in_pt;
  config_PRECISION out_pt;

#ifdef HAVE_MULT_TM
  // TODO: to correctly set up output when called from oddeven_setup_PRECISION not from its block version, might need to take additional input
  //       nv must be prop. to num_inner_lattice_sites in the former case
  int jj,jjj, nv = 72*(l->num_inner_lattice_sites/2+1), nrt = (g.in_setup)?num_loop:g.num_rhs_vect;//!!!!!!
  int offset[4] = {0,12*nrt,6*nrt,12*nrt+15};
  // construct initial L = A for n=0, L = B for n=1, L row major
  for ( n=0; n<2; n++ ) {
    
    out_pt = output + n*36*num_loop;
    
    in_pt = input + offset[2*n];
    for ( j=0; j<6; j++ ) {
      VECTOR_LOOP(jj, nrt, jjj, out_pt[(6*j+j)*num_loop+jj*nv+jjj] = (complex_PRECISION) in_pt[jj+jjj];) in_pt += nrt;
    }
    
    in_pt = input+offset[2*n+1];
    for ( j=0; j<5; j++ )
      for ( i=j+1; i<6; i++ ) {
	VECTOR_LOOP(jj, nrt, jjj, out_pt[(6*j+i)*num_loop+jj*nv+jjj] = (complex_PRECISION) *in_pt;)
	VECTOR_LOOP(jj, nrt, jjj, out_pt[(6*i+j)*num_loop+jj*nv+jjj] = (complex_PRECISION) conj_double(*in_pt);) in_pt++;
      }
    
    // calculate LU
    for ( k=0; k<5; k++ ) {
      for ( i=k+1; i<6; i++ ) {// L: out(i,k) = out(i,k)/out(k,k)
        VECTOR_LOOP(jj, nrt, jjj, out_pt[(6*i+k)*num_loop+jj*nv+jjj] = out_pt[(6*i+k)*num_loop+jj*nv+jjj]/out_pt[(6*k+k)*num_loop+jj*nv+jjj];) 
        for ( j=k+1; j<6; j++ )  // U: out(i,j) = out(i,j)-out(i,k)*out(k,j)
          VECTOR_LOOP(jj, nrt, jjj,
		      out_pt[(6*i+j)*num_loop+jj*nv+jjj] = out_pt[(6*i+j)*num_loop+jj*nv+jjj]-out_pt[(6*i+k)*num_loop+jj*nv+jjj]*out_pt[(6*k+j)*num_loop+jj*nv+jjj];)
      }
    }
  }
#else
  int offset[4] = {0,12,6,27};
  // construct initial L = A for n=0, L = B for n=1, L row major
  for ( n=0; n<2; n++ ) {

    out_pt = output + n*36;
    
    in_pt = input + offset[2*n];
    for ( j=0; j<6; j++ ) {
      out_pt[6*j+j] = (complex_PRECISION) *in_pt; in_pt++;
    }
    
    in_pt = input+offset[2*n+1];
    for ( j=0; j<5; j++ ) {
      for ( i=j+1; i<6; i++ ) {
        out_pt[6*j+i] = (complex_PRECISION) *in_pt;
        out_pt[6*i+j] = (complex_PRECISION) conj_double(*in_pt); in_pt++;
      }
    }
    
    // calculate LU
    for ( k=0; k<5; k++ ) {
      for ( i=k+1; i<6; i++ ) {
        out_pt[6*i+k] = out_pt[6*i+k]/out_pt[6*k+k]; // L: out(i,k) = out(i,k)/out(k,k)
        for ( j=k+1; j<6; j++ )
          out_pt[6*i+j] = out_pt[6*i+j]-out_pt[6*i+k]*out_pt[6*k+j]; // U: out(i,j) = out(i,j)-out(i,k)*out(k,j)
      }
    }
  }
#endif
}
#endif

#ifdef HAVE_TM1p1
void selfcoupling_LU_doublet_decomposition_PRECISION( const config_PRECISION output, config_double input, level_struct *l ) {

  /*********************************************************************************   
   * Performs a LU decomposition for a selfcoupling term U(x).
   * U = [ u(A)  0   ]   , A=A*, B=B* (diagonals excluded)
   *     [  0   u(B) ]
   * u = [ A/B(mu) -eps*I    ]
   *     [ -eps*I  A/B(-mu) ]
   * Input ordering: diag(A), diag(B), triu(A,1) row major, triu(B,1) row major
   *                  (matlab notation), epsbar_term
   * Output ordering: triu(L,1) + tril(U,0), i.e., output contains L and U without 
   *                  the diagonal of L which is equal to 1 for u(mu) and then for u(-mu)
   *                  (11*12/2+12*13/2)*2 = 288
   *********************************************************************************/

  register int i, j, k;
  int n;
  config_double in_pt;
  config_PRECISION out_pt = output;


#ifdef HAVE_MULT_TM
  // TODO: to correctly set up output when called from oddeven_setup_PRECISION not from its block version, might need to take additional input
  //       nv must be prop. to num_inner_lattice_sites in the former case
  // nrt = num_loop during setup as the smallest supplied mu value is used in constructing interpolation op
  int jj, jjj, nv = 288*(l->num_inner_lattice_sites/2+1), nrt = (g.in_setup)?num_loop:g.num_rhs_vect;
  int offset[8] = {0,12*nrt,24*nrt,30+24*nrt,6*nrt,18*nrt,15+24*nrt,36+24*nrt}; //diagonal terms in the input exist for each rhs and flavor
  
  // construct initial L = A for n=0, L = B for n=1, L row major
  for ( n=0; n<2; n++ ) {

    out_pt = output + n*144*num_loop;

    // diagonal part
    in_pt = input+offset[4*n];
    for ( j=0; j<6; j++ ) {
      VECTOR_LOOP(jj, nrt, jjj, out_pt[(12*j+j)*num_loop+jj*nv+jjj] = (complex_PRECISION) in_pt[jj+jjj];) in_pt+=nrt;
    }
    in_pt = input+offset[4*n+1];
    for ( j=0; j<6; j++ ) {
      VECTOR_LOOP(jj, nrt, jjj, out_pt[(12*(j+6)+(j+6))*num_loop+jj*nv+jjj] = (complex_PRECISION) in_pt[jj+jjj];) in_pt+=nrt;
    }

    // triu(A,1) and  triu(B,1)
    in_pt = input+offset[4*n+2];
    for ( j=0; j<5; j++ ) {
      for ( i=j+1; i<6; i++ ) {
        VECTOR_LOOP(jj, nrt, jjj, out_pt[(12*(j+6)+(i+6))*num_loop+jj*nv+jjj] = out_pt[(12*j+i)*num_loop+jj*nv+jjj] = (complex_PRECISION) *in_pt;)
	VECTOR_LOOP(jj, nrt, jjj,out_pt[(12*(i+6)+(j+6))*num_loop+jj*nv+jjj] = out_pt[(12*i+j)*num_loop+jj*nv+jjj] = (complex_PRECISION) conj_double(*in_pt);)
	in_pt++;
      }
    }

    // epsbar_term
    in_pt = input+offset[4*n+3];
    for ( j=0; j<6; j++ ) {
      for ( i=0; i<6; i++ ) {
        VECTOR_LOOP(jj, nrt, jjj, out_pt[(12*(j+6)+i)*num_loop+jj*nv+jjj] = out_pt[(12*j+(i+6))*num_loop+jj*nv+jjj] = _COMPLEX_PRECISION_ZERO;)
      }
    }
    for ( j=0; j<6; j++ ) {
      VECTOR_LOOP(jj, nrt, jjj, out_pt[(12*(j+6)+j)*num_loop+jj*nv+jjj] = out_pt[(12*j+(j+6))*num_loop+jj*nv+jjj] = (complex_PRECISION) *in_pt;) in_pt++;
    }
    
    // calculate LU
    for ( k=0; k<11; k++ ) {
      for ( i=k+1; i<12; i++ ) {
        VECTOR_LOOP(jj, nrt, jjj, out_pt[(12*i+k)*num_loop+jj*nv+jjj] = out_pt[(12*i+k)*num_loop+jj*nv+jjj]/out_pt[(12*k+k)*num_loop+jj*nv+jjj];) // L: out(i,k) = out(i,k)/out(k,k)
        for ( j=k+1; j<12; j++ )
          VECTOR_LOOP(jj, nrt, jjj, out_pt[(12*i+j)*num_loop+jj*nv+jjj] =
		      out_pt[(12*i+j)*num_loop+jj*nv+jjj]-out_pt[(12*i+k)*num_loop+jj*nv+jjj]*out_pt[(12*k+j)*num_loop+jj*nv+jjj];) // U: out(i,j) = out(i,j)-out(i,k)*out(k,j)
      }
    }
  }
#else
  int offset[8] = {0,12,24,54,6,18,39,60};
  
  // construct initial L = A for n=0, L = B for n=1, L row major
  for ( n=0; n<2; n++ ) {

    out_pt = output + n*144;

    in_pt = input+offset[4*n];
    for ( j=0; j<6; j++ ) {
      out_pt[12*j+j] = (complex_PRECISION) *in_pt; in_pt++;
    }
    in_pt = input+offset[4*n+1];
    for ( j=0; j<6; j++ ) {
      out_pt[12*(j+6)+(j+6)] = (complex_PRECISION) *in_pt; in_pt++;
    }
    
    in_pt = input+offset[4*n+2];
    for ( j=0; j<5; j++ ) {
      for ( i=j+1; i<6; i++ ) {
        out_pt[12*(j+6)+(i+6)] = out_pt[12*j+i] = (complex_PRECISION) *in_pt;
        out_pt[12*(i+6)+(j+6)] = out_pt[12*i+j] = (complex_PRECISION) conj_double(*in_pt); in_pt++;
      }
    }
    
    in_pt = input+offset[4*n+3];
    for ( j=0; j<6; j++ ) {
      for ( i=0; i<6; i++ ) {
        out_pt[12*(j+6)+i] = out_pt[12*j+(i+6)] = _COMPLEX_PRECISION_ZERO;
      }
    }
    for ( j=0; j<6; j++ ) {
      out_pt[12*(j+6)+j] = out_pt[12*j+(j+6)] = (complex_PRECISION) *in_pt; in_pt++;
    }
    
    // calculate LU
    for ( k=0; k<11; k++ ) {
      for ( i=k+1; i<12; i++ ) {
        out_pt[12*i+k] = out_pt[12*i+k]/out_pt[12*k+k]; // L: out(i,k) = out(i,k)/out(k,k)
        for ( j=k+1; j<12; j++ )
          out_pt[12*i+j] = out_pt[12*i+j]-out_pt[12*i+k]*out_pt[12*k+j]; // U: out(i,j) = out(i,j)-out(i,k)*out(k,j)
      }
    }
  }
#endif
}
#endif

// out_o/e += D_oe/eo*in_e/o if amount==_ODD_SITES/_EVEN_SITES; works on both parts if amount==_FULL_SYSTEM
// in and out needs to be at least of size _ORDINARY 
void hopping_term_PRECISION( vector_PRECISION *eta, vector_PRECISION *phi, operator_PRECISION_struct *op,
                             const int amount, level_struct *l, struct Thread *threading ) {

  int start_even, end_even, start_odd, end_odd;
  int start=0, plus_dir_param=_FULL_SYSTEM, minus_dir_param=_FULL_SYSTEM;
  int n         = l->num_inner_lattice_sites;
  int *neighbor = op->neighbor_table;
  int nvec = phi->num_vect_now, nvec_phi = phi->num_vect, nvec_eta = eta->num_vect, nvec_op = op->pr_num_vect;

  SYNC_CORES(threading)  

#ifdef DEBUG
    if ( nvec != eta->num_vect_now || nvec_op < nvec_phi )
      error0("hopping_term_PRECISION: assumptions are not met\n");
#ifdef HAVE_TM1p1
  if ( nvec != 2*num_loop || nvec != nvec_phi )
    error0("hopping_term_PRECISION: assumptions are not met for TM1p1\n");
#endif
#endif
    
  if ( amount == _EVEN_SITES || amount == _ODD_SITES ) {
    compute_core_start_end_custom(0, op->num_even_sites, &start_even, &end_even, l, threading, 1 );
    compute_core_start_end_custom(op->num_even_sites, op->num_even_sites+op->num_odd_sites, &start_odd, &end_odd, l, threading, 1 );
  } else {
    compute_core_start_end_custom(0, l->num_inner_lattice_sites, &start, &n, l, threading, 1 );
  }
  
  SYNC_CORES(threading)  
  
  if ( amount == _EVEN_SITES ) {
    start = start_odd, n = end_odd;
    minus_dir_param = _ODD_SITES;
    plus_dir_param = _EVEN_SITES;
  } else if ( amount == _ODD_SITES ) {
    start = start_even, n = end_even;
    minus_dir_param = _EVEN_SITES;
    plus_dir_param = _ODD_SITES;
  }

  int i, *nb_pt;
  buffer_PRECISION phi_pt, eta_pt, end_pt;
  config_PRECISION D_pt;
  g.num_vect_pass1 = nvec; 
  
  // project in negative directions
  complex_PRECISION pbuf[6*nvec];
  for ( i=6*start, phi_pt=phi->vector_buffer+12*start*nvec_phi; i<6*n; i+=6, phi_pt+=12*nvec_phi ) {
    prp_T_PRECISION( op->prnT+i*nvec_op, phi_pt, nvec, nvec_op, nvec_phi );
    prp_Z_PRECISION( op->prnZ+i*nvec_op, phi_pt, nvec, nvec_op, nvec_phi );
    prp_Y_PRECISION( op->prnY+i*nvec_op, phi_pt, nvec, nvec_op, nvec_phi );
    prp_X_PRECISION( op->prnX+i*nvec_op, phi_pt, nvec, nvec_op, nvec_phi );
  }
  // start communication in negative direction
  START_LOCKED_MASTER(threading)
  g.num_vect_pass2 = nvec_op;
  ghost_sendrecv_PRECISION( op->prnT, T, -1, &(op->c), minus_dir_param, l );
  ghost_sendrecv_PRECISION( op->prnZ, Z, -1, &(op->c), minus_dir_param, l );
  ghost_sendrecv_PRECISION( op->prnY, Y, -1, &(op->c), minus_dir_param, l );
  ghost_sendrecv_PRECISION( op->prnX, X, -1, &(op->c), minus_dir_param, l );
  END_LOCKED_MASTER(threading) 
  // project plus dir and multiply with U dagger
  for ( phi_pt=phi->vector_buffer+12*start*nvec_phi, end_pt=phi->vector_buffer+12*n*nvec_phi, D_pt = op->D+36*start, nb_pt=neighbor+4*start; phi_pt<end_pt; phi_pt+=12*nvec_phi ) {
    // T dir
    i = 6*(*nb_pt); nb_pt++;
    prn_T_PRECISION( pbuf, phi_pt, nvec, nvec, nvec_phi );
    mvmh_PRECISION( op->prpT+i*nvec_op, D_pt, pbuf, nvec, nvec_op, nvec );
    mvmh_PRECISION( op->prpT+(i+3)*nvec_op, D_pt, pbuf+3*nvec, nvec, nvec_op, nvec ); D_pt += 9;
    // Z dir
    i = 6*(*nb_pt); nb_pt++;
    prn_Z_PRECISION( pbuf, phi_pt, nvec, nvec, nvec_phi );
    mvmh_PRECISION( op->prpZ+i*nvec_op, D_pt, pbuf, nvec, nvec_op, nvec );
    mvmh_PRECISION( op->prpZ+(i+3)*nvec_op, D_pt, pbuf+3*nvec, nvec, nvec_op, nvec ); D_pt += 9;
    // Y dir
    i = 6*(*nb_pt); nb_pt++;
    prn_Y_PRECISION( pbuf, phi_pt, nvec, nvec, nvec_phi );
    mvmh_PRECISION( op->prpY+i*nvec_op, D_pt, pbuf, nvec, nvec_op, nvec );
    mvmh_PRECISION( op->prpY+(i+3)*nvec_op, D_pt, pbuf+3*nvec, nvec, nvec_op, nvec ); D_pt += 9;
    // X dir
    i = 6*(*nb_pt); nb_pt++;
    prn_X_PRECISION( pbuf, phi_pt, nvec, nvec, nvec_phi );
    mvmh_PRECISION( op->prpX+i*nvec_op, D_pt, pbuf, nvec, nvec_op, nvec );
    mvmh_PRECISION( op->prpX+(i+3)*nvec_op, D_pt, pbuf+3*nvec, nvec, nvec_op, nvec ); D_pt += 9;
  }
  if ( amount == _EVEN_SITES ) {
    start = start_even, n = end_even;
  } else if ( amount == _ODD_SITES ) {
    start = start_odd, n = end_odd;
  }  
  // start communication in positive direction
  START_LOCKED_MASTER(threading)
  ghost_sendrecv_PRECISION( op->prpT, T, +1, &(op->c), plus_dir_param, l );
  ghost_sendrecv_PRECISION( op->prpZ, Z, +1, &(op->c), plus_dir_param, l );
  ghost_sendrecv_PRECISION( op->prpY, Y, +1, &(op->c), plus_dir_param, l );
  ghost_sendrecv_PRECISION( op->prpX, X, +1, &(op->c), plus_dir_param, l );
  // wait for communication in negative direction
  ghost_wait_PRECISION( op->prnT, T, -1, &(op->c), minus_dir_param, l );
  ghost_wait_PRECISION( op->prnZ, Z, -1, &(op->c), minus_dir_param, l );
  ghost_wait_PRECISION( op->prnY, Y, -1, &(op->c), minus_dir_param, l );
  ghost_wait_PRECISION( op->prnX, X, -1, &(op->c), minus_dir_param, l );
  END_LOCKED_MASTER(threading) 
  // multiply with U and lift up minus dir
  for ( eta_pt=eta->vector_buffer+12*start*nvec_eta, end_pt=eta->vector_buffer+12*n*nvec_eta, D_pt = op->D+36*start, nb_pt=neighbor+4*start; eta_pt<end_pt; eta_pt+=12*nvec_eta ) {
    // T dir
    i = 6*(*nb_pt); nb_pt++;
    mvm_PRECISION( pbuf, D_pt, op->prnT+i*nvec_op, nvec, nvec, nvec_op );
    mvm_PRECISION( pbuf+3*nvec, D_pt, op->prnT+(i+3)*nvec_op, nvec, nvec, nvec_op );
    pbp_su3_T_PRECISION( pbuf, eta_pt, nvec, nvec, nvec_eta ); D_pt += 9;
    // Z dir
    i = 6*(*nb_pt); nb_pt++;
    mvm_PRECISION( pbuf, D_pt, op->prnZ+i*nvec_op, nvec, nvec, nvec_op );
    mvm_PRECISION( pbuf+3*nvec, D_pt, op->prnZ+(i+3)*nvec_op, nvec, nvec, nvec_op );
    pbp_su3_Z_PRECISION( pbuf, eta_pt, nvec, nvec, nvec_eta ); D_pt += 9;
    // Y dir
    i = 6*(*nb_pt); nb_pt++;
    mvm_PRECISION( pbuf, D_pt, op->prnY+i*nvec_op, nvec, nvec, nvec_op );
    mvm_PRECISION( pbuf+3*nvec, D_pt, op->prnY+(i+3)*nvec_op, nvec, nvec, nvec_op );
    pbp_su3_Y_PRECISION( pbuf, eta_pt, nvec, nvec, nvec_eta ); D_pt += 9;
    // X dir
    i = 6*(*nb_pt); nb_pt++;
    mvm_PRECISION( pbuf, D_pt, op->prnX+i*nvec_op, nvec, nvec, nvec_op );
    mvm_PRECISION( pbuf+3*nvec, D_pt, op->prnX+(i+3)*nvec_op, nvec, nvec, nvec_op );
    pbp_su3_X_PRECISION( pbuf, eta_pt, nvec, nvec, nvec_eta ); D_pt += 9;
  }
  // wait for communication in positive direction
  START_LOCKED_MASTER(threading)
  ghost_wait_PRECISION( op->prpT, T, +1, &(op->c), plus_dir_param, l );
  ghost_wait_PRECISION( op->prpZ, Z, +1, &(op->c), plus_dir_param, l );
  ghost_wait_PRECISION( op->prpY, Y, +1, &(op->c), plus_dir_param, l );
  ghost_wait_PRECISION( op->prpX, X, +1, &(op->c), plus_dir_param, l );
  END_LOCKED_MASTER(threading) 
  // lift up plus dir
  for ( i=6*start, eta_pt=eta->vector_buffer+12*start*nvec_eta; i<6*n; i+=6, eta_pt+=12*nvec_eta ) {
    pbn_su3_T_PRECISION( op->prpT+i*nvec_op, eta_pt, nvec, nvec_op, nvec_eta );
    pbn_su3_Z_PRECISION( op->prpZ+i*nvec_op, eta_pt, nvec, nvec_op, nvec_eta );
    pbn_su3_Y_PRECISION( op->prpY+i*nvec_op, eta_pt, nvec, nvec_op, nvec_eta );
    pbn_su3_X_PRECISION( op->prpX+i*nvec_op, eta_pt, nvec, nvec_op, nvec_eta );
  }

  SYNC_CORES(threading)
}

static void diag_ee_PRECISION( vector_PRECISION *y, vector_PRECISION *x, operator_PRECISION_struct *op, 
			       level_struct *l, int start, int end ) {

/*********************************************************************************
* Applies the even-even block of the odd even decomposition to a vector.
* - vector_PRECISION *x: Input vector.
* - vector_PRECISION *y: Output vector.
* used when g.method > 3 in which case clover term is LU decomposed and clover_doublet_oo_inv is defined both on even and odd sites
* clover_oo_inv is not defined for non-block oe_op
*********************************************************************************/

  int i, j, jj, nvec = x->num_vect_now, nvec_x = x->num_vect, nvec_y = y->num_vect;
#ifdef HAVE_TM
  int factor = 1;
#ifdef HAVE_MULT_TM
  factor *= num_loop;
#endif
#endif
  
#ifdef DEBUG
  if ( nvec_y < nvec )
    error0("diag_ee_PRECISIO: assumptions are not met\n");
#endif
  
#ifdef HAVE_TM1p1
  if( g.n_flavours == 2) {
    vector_PRECISION x_pt, y_pt;
    vector_PRECISION_duplicate( &x_pt, x, start/12, l );
    vector_PRECISION_duplicate( &y_pt, y, start/12, l );
    config_PRECISION sc = op->clover_doublet_oo_inv + (start/12)*288*factor;
    LU_multiply_PRECISION( &y_pt, &x_pt, sc, start/12, end/12, l);
  } else {
#endif
    if ( g.csw ) {
      vector_PRECISION x_pt, y_pt;
      vector_PRECISION_duplicate( &x_pt, x, start/12, l );
      vector_PRECISION_duplicate( &y_pt, y, start/12, l );
#if defined(HAVE_TM) 
      config_PRECISION sc = op->clover + (start/12)*72*factor;
      LU_multiply_PRECISION( &y_pt, &x_pt, sc, start/12, end/12, l );
#else
      config_PRECISION sc = op->clover + (start/12)*42;
      LLH_multiply_PRECISION( &y_pt, &x_pt, sc, start/12, end/12 );
#endif
    } else {//here op->clover contains all the self-couling terms which is diff from other functions????
      buffer_PRECISION x_pt = x->vector_buffer+start*nvec_x, y_pt = y->vector_buffer+start*nvec_y;
      config_PRECISION sc = op->clover + start;
      for ( i=start/12; i<end/12; i++ )
	FOR12( VECTOR_LOOP( j, nvec, jj, y_pt[j+jj] = x_pt[j+jj]*(*sc); ) sc++; x_pt += nvec_x; y_pt += nvec_y;)
    }
#ifdef HAVE_TM1p1
  }
#endif
}

static void diag_oo_PRECISION( vector_PRECISION *y, vector_PRECISION *x, operator_PRECISION_struct *op,
				   level_struct *l, struct Thread *threading ) {

/*********************************************************************************
* Applies the odd-odd block of the odd even decomposition to a vector.
* - vector_PRECISION *x: Input vector.
* - vector_PRECISION *y: Output vector.
* used when g.method > 3 in which case clover term is LU decomposed and clover_doublet_oo_inv is defined both on even and odd sites
* clover_oo_inv is not defined for non-block oe_op
*********************************************************************************/

  START_UNTHREADED_FUNCTION(threading)

  int i, ne = op->num_even_sites, no = op->num_odd_sites;
  int j, jj, nvec = x->num_vect_now, nvec_x = x->num_vect, nvec_y = y->num_vect;
#ifdef HAVE_TM
  int factor = 1;
#ifdef HAVE_MULT_TM
  factor *= num_loop;
#endif
#endif
  
#ifdef DEBUG
  if ( nvec_y < nvec )
    error0("diag_oo_PRECISION: assumptions are not met\n");
#endif
    
  config_PRECISION sc = op->clover;
#ifdef HAVE_TM1p1
  if( g.n_flavours == 2) {
    vector_PRECISION x_pt, y_pt;
    vector_PRECISION_duplicate( &x_pt, x, ne, l );
    vector_PRECISION_duplicate( &y_pt, y, ne, l );
    config_PRECISION sc = op->clover_doublet_oo_inv + ne*288*factor;
    LU_multiply_PRECISION( &y_pt, &x_pt, sc, 0, no, l );
  } else {
#endif
    // diagonal blocks applied to the odd sites
    if ( g.csw ) {
      vector_PRECISION x_pt, y_pt;
      vector_PRECISION_duplicate( &x_pt, x, ne, l );
      vector_PRECISION_duplicate( &y_pt, y, ne, l );
#ifndef HAVE_TM
      sc += ne*42;
      LLH_multiply_PRECISION( &y_pt, &x_pt, sc, 0, no );
#else
      sc += ne*72*factor;
      LU_multiply_PRECISION( &y_pt, &x_pt, sc, 0, no, l );
#endif
    } else {//here op->clover contains all the self-couling terms which is diff from other functions
      buffer_PRECISION x_pt = x->vector_buffer+12*ne*nvec_x, y_pt = y->vector_buffer+12*ne*nvec_y;
      sc += ne*12;
      for ( i=0; i<no; i++ )
        FOR12(VECTOR_LOOP( j, nvec, jj, y_pt[j+jj] = x_pt[j+jj]*(*sc); ) sc++; x_pt += nvec_x; y_pt += nvec_y;)
    }
#ifdef HAVE_TM1p1
  }
#endif

  END_UNTHREADED_FUNCTION(threading)
}

// used when g.method > 3 in which case clover term is LU decomposed and clover_doublet_oo_inv is defined both on even and odd sites
// NOTE: clover_oo_inv is not defined for non-block oe_op
static void diag_oo_inv_PRECISION( vector_PRECISION *y, vector_PRECISION *x, operator_PRECISION_struct *op,
                            level_struct *l, int start, int end ) {

  int i, j, jj, nvec = x->num_vect_now, nvec_x = x->num_vect, nvec_y = y->num_vect;
#ifdef HAVE_TM
  int factor = 1;
#ifdef HAVE_MULT_TM
  factor *= num_loop;
#endif
#endif

#ifdef DEBUG
  if ( nvec_y < nvec )
    error0("diag_oo_inv_PRECISION: assumptions are not met\n");
#endif
  
#ifdef HAVE_TM1p1
  if( g.n_flavours == 2) {
    vector_PRECISION y_pt, x_pt;
    vector_PRECISION_duplicate( &y_pt, y, start/12, l );
    vector_PRECISION_duplicate( &x_pt, x, start/12, l );
    // inverse diagonal blocks applied to the odd sites
    config_PRECISION sc = op->clover_doublet_oo_inv + (start/12)*288*factor;
    LU_perform_fwd_bwd_subs_PRECISION( &y_pt, &x_pt, sc, start/12, end/12, l );
  } else {
#endif
  config_PRECISION sc = op->clover;
  if ( g.csw ) {
    vector_PRECISION y_pt, x_pt;
    vector_PRECISION_duplicate( &y_pt, y, start/12, l );
    vector_PRECISION_duplicate( &x_pt, x, start/12, l );
#if defined(HAVE_TM)
    sc += (start/12)*72*factor;
    LU_perform_fwd_bwd_subs_PRECISION( &y_pt, &x_pt, sc, start/12, end/12, l );
#else
    sc += (start/12)*42;
    LLH_perform_fwd_bwd_subs_PRECISION( &y_pt, &x_pt, sc, start/12, end/12 );
#endif
  } else {
    buffer_PRECISION x_pt = x->vector_buffer+start*nvec_x, y_pt = y->vector_buffer+start*nvec_y;
    sc += start;
    for ( i=start/12; i<end/12; i++ )
      FOR12( VECTOR_LOOP( j, nvec, jj, y_pt[j+jj] = x_pt[j+jj]/(*sc);) sc++; x_pt += nvec_x; y_pt += nvec_y;)
  }
#ifdef HAVE_TM1p1
  }
#endif
}

// used only when g.method==4, 5  
void apply_schur_complement_PRECISION( vector_PRECISION *out, vector_PRECISION *in, operator_PRECISION_struct *op,
    level_struct *l, struct Thread *threading ) {

  /*********************************************************************************
   * Applies the Schur complement to a vector.
   *********************************************************************************/

  int start_even, end_even, start_odd, end_odd;

  compute_core_start_end_custom(0, op->num_even_sites*l->num_lattice_site_var, &start_even, &end_even, l, threading, l->num_lattice_site_var );
  compute_core_start_end_custom(op->num_even_sites*l->num_lattice_site_var, l->inner_vector_size, &start_odd, &end_odd, l, threading, l->num_lattice_site_var );
#ifdef DEBUG
  if ( out->num_vect < in->num_vect_now )
    error0("apply_schur_complement_PRECISION: assumptions are not met\n");
#endif
  
  vector_PRECISION *tmp = op->buffer;
  tmp[0].num_vect_now = in->num_vect_now; tmp[1].num_vect_now = in->num_vect_now;

  SYNC_CORES(threading)
  vector_PRECISION_define( &tmp[0], 0, start_odd, end_odd, l );
  vector_PRECISION_define( &tmp[0], 0, start_even, end_even, l );
  
  SYNC_CORES(threading)
  PROF_PRECISION_START( _NC, threading );
  PROF_PRECISION_START( _SC, threading );
  diag_ee_PRECISION( out, in, op, l, start_even, end_even );
  SYNC_CORES(threading)
  PROF_PRECISION_STOP( _SC, 1, threading );
  
  hopping_term_PRECISION( &tmp[0], in, op, _ODD_SITES, l, threading );
  PROF_PRECISION_STOP( _NC, 0, threading );
  PROF_PRECISION_START( _SC, threading );
  diag_oo_inv_PRECISION( &tmp[1], &tmp[0], op, l, start_odd, end_odd );
  SYNC_CORES(threading)
  PROF_PRECISION_STOP( _SC, 0, threading );
  PROF_PRECISION_START( _NC, threading );
  hopping_term_PRECISION( &tmp[0], &tmp[1], op, _EVEN_SITES, l, threading );
  PROF_PRECISION_STOP( _NC, 1, threading );
  vector_PRECISION_minus( out, out, &tmp[0], start_even, end_even, l );
}

// used only in void preconditioner() when g.method==4, 5
void solve_oddeven_PRECISION( gmres_PRECISION_struct *p, operator_PRECISION_struct *op, level_struct *l, struct Thread *threading ) {
  /****************************************************************************************  
   * Descripton: Solve D*x = b in the even-odd preconditioning using Schur complement D_sc
   *   D_sc = D_ee - D_eo D_oo ^{-1} D_oe
   *   x_e = D_sc^{-1} (b_e - D_eo*D_oo^{-1}*b_o)
   *   x_o = D_oo^{-1} (b_o - D_oe*x_e)  
   ****************************************************************************************/

  int start;
  int end;
  compute_core_start_end_custom(op->num_even_sites*l->num_lattice_site_var, l->inner_vector_size, &start, &end, l, threading, l->num_lattice_site_var );

  int j, jj, nvec = num_loop;
  complex_PRECISION factor[nvec];
  VECTOR_LOOP(j, nvec, jj, factor[j+jj] = -1;)

  vector_PRECISION tmp = op->buffer[0]; // op->buffer[0] is of size _ORDINARY
  tmp.num_vect_now = nvec; p->x.num_vect_now = nvec; p->b.num_vect_now = nvec;
  
  // solve for x_e
  PROF_PRECISION_START( _SC, threading );
  diag_oo_inv_PRECISION( &tmp, &(p->b), op, l, start, end ); // tmp_o = D_oo^{-1}*b_o
  PROF_PRECISION_STOP( _SC, 0, threading );
  SYNC_CORES(threading)
  vector_PRECISION_scale( &tmp, &tmp, factor, 0, start, end, l ); // tmp_o = -D_oo^{-1}*b_o
  SYNC_CORES(threading)
  PROF_PRECISION_START( _NC, threading );
  hopping_term_PRECISION( &(p->b), &tmp, op, _EVEN_SITES, l, threading ); // b_e = b_e - D_eo*D_oo^{-1}*b_o
  PROF_PRECISION_STOP( _NC, 0, threading );

  // x_e = D_sc^{-1} (b_e - D_eo*D_oo^{-1}*b_o); op in p is set to apply_schur_complement_PRECISION
  if ( g.method == 4 )
    fgmres_PRECISION( p, l, threading );
  else if ( g.method == 5 )
    bicgstab_PRECISION( p, l, threading );
  
  // construct x_o from x_e 
  diag_oo_inv_PRECISION( &(p->x), &(p->b), op, l, start, end ); // x_o = D_oo^{-1}*b_o 
  SYNC_CORES(threading)
  vector_PRECISION_define( &tmp, 0, start, end, l ); // tmp_o = 0
  SYNC_CORES(threading)
  PROF_PRECISION_START( _NC, threading );
  hopping_term_PRECISION( &tmp, &(p->x), op, _ODD_SITES, l, threading ); // tmp_o = D_oe*x_e
  PROF_PRECISION_STOP( _NC, 1, threading );
  PROF_PRECISION_START( _SC, threading );
  diag_oo_inv_PRECISION( &(p->b), &tmp, op, l, start, end ); // b_o = D_oo^{-1}*D_oe*x_e
  PROF_PRECISION_STOP( _SC, 1, threading );
  SYNC_CORES(threading)
  vector_PRECISION_minus( &(p->x), &(p->x), &(p->b), start, end, l ); // x_o = D_oo^{-1}*(b_o - D_oe*x_e)
  SYNC_CORES(threading)
}

void oddeven_to_serial_PRECISION( vector_double *out, vector_PRECISION *in, level_struct *l, struct Thread *threading ) {

/*********************************************************************************
* Translates a vector from an odd even PRECISION precision layout to a serial 
* double precision layout.
*********************************************************************************/ 

  int i, j, k, jj, jjj;
  int nvec = in->num_vect_now, nvec_in = in->num_vect, nvec_out = out->num_vect;
  int nsv = l->num_lattice_site_var, *tt = l->oe_op_PRECISION.translation_table;
  int start = threading->start_site[l->depth];
  int end   = threading->end_site[l->depth];

  if ( nvec_out < nvec )
    error0("oddeven_to_serial_PRECISION: assumptions are not met\n");

  // this function seems to do some data reordering, barriers ensure that everything is in sync
  SYNC_CORES(threading)
  START_NO_HYPERTHREADS(threading)  
  for ( i=start; i<end; i++ ) {
    k = tt[i];
    for ( j=0; j<nsv; j++ ) {
      VECTOR_LOOP( jj, nvec, jjj, out->vector_buffer[(i*nsv+j)*nvec_out+jj+jjj] = (complex_double) in->vector_buffer[(k*nsv+j)*nvec_in+jj+jjj];)
    }
  }
  END_NO_HYPERTHREADS(threading)
  SYNC_CORES(threading)  
}

void serial_to_oddeven_PRECISION( vector_PRECISION *out, vector_double *in, level_struct *l, struct Thread *threading ) {

/*********************************************************************************
* Translates a vector from a serial double precision layout to an odd even
* PRECISION precision layout.
*********************************************************************************/ 

  int i, j, k, jj, jjj;
  int nvec = in->num_vect_now, nvec_in = in->num_vect, nvec_out = out->num_vect;
  int nsv = l->num_lattice_site_var, *tt = l->oe_op_PRECISION.translation_table;
  int start = threading->start_site[l->depth];
  int end   = threading->end_site[l->depth];

  if ( nvec_out < nvec )
    error0("oddeven_to_serial_PRECISION: assumptions are not met\n");

  // this function seems to do some data reordering, barriers ensure that everything is in sync
  SYNC_CORES(threading)
  START_NO_HYPERTHREADS(threading)  
  for ( i=start; i<end; i++ ) {
    k = tt[i];
    for ( j=0; j<nsv; j++ ) {
      VECTOR_LOOP( jj, nvec, jjj, out->vector_buffer[(k*nsv+j)*nvec_out+jj+jjj] = (complex_PRECISION) in->vector_buffer[(i*nsv+j)*nvec_in+jj+jjj];)
    }
  }
  END_NO_HYPERTHREADS(threading)
  SYNC_CORES(threading)  
}

/*******************  TEST ROUTINES  *****************************************************/

void oddeven_PRECISION_test( level_struct *l ) {

/*********************************************************************************
* In the first step this function checks the correctness of the odd even layout. 
* This is done by:
* - Applying D_W in odd-even structure to a vector.
* - Applying D_W directly to the same vector.
* - Compare solutions ( Difference should be close to 0 ).
*********************************************************************************/  

  int nvec = num_loop, nvecsf = num_loop;
#ifdef HAVE_TM1p1
  if ( g.n_flavours > 1 ) {
      nvec *= g.n_flavours;
      if( ! (g.epsbar != 0 || g.epsbar_ig5_odd_shift != 0 || g.epsbar_ig5_odd_shift != 0))
        nvecsf *=  g.n_flavours;
    }
#endif
  vector_double d[3];
  vector_PRECISION f[5];
  double diff1_d[nvecsf], diff2_d[nvecsf];
  PRECISION diff1[nvecsf], diff2[nvecsf];

  for(int i=0; i<3; i++){
    vector_double_init( &d[i] );
    vector_double_alloc( &d[i], _INNER, nvec, l, no_threading );
    d[i].num_vect_now = nvec;
  }

  for(int i=0; i<5; i++){                                                                 
    vector_PRECISION_init( &f[i] );                                                          
    vector_PRECISION_alloc( &f[i], _INNER, nvec, l, no_threading );                             
    f[i].num_vect_now = nvec;
  } 

  vector_double_define_random( &d[0], 0, l->inner_vector_size, l ); 
  serial_to_oddeven_PRECISION( &f[0], &d[0], l, no_threading );
   
  diag_ee_PRECISION( &f[1], &f[0], &(l->oe_op_PRECISION), l, 0, l->oe_op_PRECISION.num_even_sites*l->num_lattice_site_var );
  diag_oo_PRECISION( &f[1], &f[0], &(l->oe_op_PRECISION), l, no_threading );
  
  hopping_term_PRECISION( &f[1], &f[0], &(l->oe_op_PRECISION), _FULL_SYSTEM, l, no_threading );
  
  d_plus_clover_double( &d[1], &d[0], &(g.op_double), l, no_threading );
  oddeven_to_serial_PRECISION( &d[0], &f[1], l, no_threading );
  
  vector_double_minus( &d[2], &d[0], &d[1], 0, l->num_inner_lattice_sites, l );
  global_norm_double( diff1_d, &d[2], 0, l->num_inner_lattice_sites, l, no_threading );
  global_norm_double( diff2_d, &d[0], 0, l->num_inner_lattice_sites, l, no_threading );
  
  for( int i=0; i<nvecsf; i++ )
    test0_PRECISION("depth: %d, correctness of odd even layout: %le\n", l->depth, diff1_d[i]/diff2_d[i] );
    
  // -------------- //
  
  vector_PRECISION_copy( &f[3], &f[0], 0, l->inner_vector_size, l );
  diag_oo_PRECISION( &f[2], &f[3], &(l->oe_op_PRECISION), l, no_threading );
  diag_oo_inv_PRECISION( &f[3], &f[2], &(l->oe_op_PRECISION), l, l->oe_op_PRECISION.num_even_sites*l->num_lattice_site_var, l->inner_vector_size );
  vector_PRECISION_minus( &f[3], &f[3], &f[0], 0, l->inner_vector_size, l );

  global_norm_PRECISION( diff1, &f[3], 0, l->inner_vector_size, l, no_threading );
  global_norm_PRECISION( diff2, &f[0], 0, l->inner_vector_size, l, no_threading );
  
  for( int i=0; i<nvecsf; i++ )
    test0_PRECISION("depth: %d, correctness of odd even diagonal term: %le\n", l->depth, diff1[i]/diff2[i] );

  // -------------- //
    
  // transformation part
  vector_PRECISION_copy( &f[3], &f[0], 0, l->inner_vector_size, l );
  // even to odd
  // set odd part of f3 to 0. 
  vector_PRECISION_define( &f[2], 0, l->oe_op_PRECISION.num_even_sites*l->num_lattice_site_var, l->inner_vector_size, l );
  
  hopping_term_PRECISION( &f[2], &f[3], &(l->oe_op_PRECISION), _ODD_SITES, l, no_threading );
  diag_oo_inv_PRECISION( &f[4], &f[2], &(l->oe_op_PRECISION), l, l->oe_op_PRECISION.num_even_sites*l->num_lattice_site_var, l->inner_vector_size );
  vector_PRECISION_plus( &f[3], &f[3], &f[4], l->oe_op_PRECISION.num_even_sites*l->num_lattice_site_var, l->inner_vector_size, l );
  
  // block diagonal part
  apply_schur_complement_PRECISION( &f[2], &f[3], &(l->oe_op_PRECISION), l, no_threading );
  diag_oo_PRECISION( &f[2], &f[3], &(l->oe_op_PRECISION), l, no_threading );
  // back transformation part
  diag_oo_inv_PRECISION( &f[4], &f[3], &(l->oe_op_PRECISION), l, l->oe_op_PRECISION.num_even_sites*l->num_lattice_site_var, l->inner_vector_size );
  hopping_term_PRECISION( &f[2], &f[4], &(l->oe_op_PRECISION), _EVEN_SITES, l, no_threading );
  
  vector_PRECISION_minus( &f[0], &f[1], &f[2], 0, l->inner_vector_size, l );
  global_norm_PRECISION( diff1, &f[0], 0, l->inner_vector_size, l, no_threading );
  global_norm_PRECISION( diff2, &f[1], 0, l->inner_vector_size, l, no_threading );
  
  for( int i=0; i<nvecsf; i++ )
    test0_PRECISION("depth: %d, correctness of odd even schur complement: %le\n", l->depth, diff1[i]/diff2[i] );

  // -------------- //

  for(int i=0; i<3; i++)
    vector_double_free( &d[i], l, no_threading );

  for(int i=0; i<5; i++)
    vector_PRECISION_free( &f[i], l, no_threading );
}
