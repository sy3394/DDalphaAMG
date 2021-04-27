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
 * glanced over: 12/08/2019: some reordering needed to be done
 * 1st cleanup:12/22/2019
 */

#include "main.h"

static void coarse_selfcoupling_LU_decomposition_PRECISION( config_PRECISION output, operator_PRECISION_struct *op, int index, level_struct *l );
#ifdef HAVE_TM1p1
static void coarse_selfcoupling_LU_doublet_decomposition_PRECISION( config_PRECISION output, operator_PRECISION_struct *op, int index, level_struct *l );
#endif
static void coarse_diag_oo_inv_PRECISION( vector_PRECISION *y, vector_PRECISION *x, operator_PRECISION_struct *op, level_struct *l, struct Thread *threading );

void coarse_oddeven_alloc_PRECISION( level_struct *l ) {

  int nv = l->num_parent_eig_vect, oe_offset=0, mu, **bt = NULL,
    *eot = NULL, *nt = NULL, *tt = NULL, t, z, y, x, le[4], N[4];
  int nvec = num_loop;
  operator_PRECISION_struct *op = &(l->oe_op_PRECISION);

#ifdef HAVE_TM1p1
  nvec *= 2;
#endif
  operator_PRECISION_alloc( op, _ODDEVEN, l );

  // allocate buffers in l->oe_op_PRECISION
  MALLOC( op->buffer, vector_PRECISION, 2 );
  for (int k=0; k<2; k++ ){
    vector_PRECISION_init( &(op->buffer[k]) );
    vector_PRECISION_alloc( &(op->buffer[k]), _ORDINARY, nvec, l, no_threading );
  }
  for ( mu=0; mu<4; mu++ ) {
    le[mu] = l->local_lattice[mu];
    N[mu] = le[mu]+1;
    op->table_dim[mu] = N[mu];
  }

  for ( mu=0; mu<4; mu++ )
    oe_offset += (l->local_lattice[mu]*(g.my_coords[mu]/l->comm_offset[mu]))%2;
  oe_offset = oe_offset%2;

  // estimate site numbers
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

  // allocate memory for decomposed matrix for total self-coupling term, which has off-diagonal part too.
#ifdef HAVE_MULT_TM
  MALLOC( op->clover_oo_inv, complex_PRECISION, g.num_rhs_vect*SQUARE(2*nv)*op->num_odd_sites );
#else
  MALLOC( op->clover_oo_inv, complex_PRECISION, SQUARE(2*nv)*op->num_odd_sites );
#endif
#ifdef HAVE_TM1p1
#ifdef HAVE_MULT_TM
  MALLOC( op->clover_doublet_oo_inv, complex_PRECISION, g.num_rhs_vect*SQUARE(4*nv)*op->num_odd_sites );
#else
  MALLOC( op->clover_doublet_oo_inv, complex_PRECISION, SQUARE(4*nv)*op->num_odd_sites );
#endif
#endif
  // define data layout
  eot = op->index_table;
  define_eot( eot, N, l );

  // define neighbor table, translation table
  nt = op->neighbor_table;
  tt = op->translation_table;
  define_nt_bt_tt( nt, op->backward_neighbor_table, NULL, tt, eot, N, l );

  // define boundary table
  bt = op->c.boundary_table;
  define_eo_bt( bt, eot, op->c.num_even_boundary_sites, op->c.num_odd_boundary_sites, op->c.num_boundary_sites, N, l );

  // setup shell communication
  ghost_sendrecv_init_PRECISION( _COARSE_GLOBAL, &(op->c), l ) ;

  // set solver parameters
  if ( l->level == 0 )
    l->p_PRECISION.v_end = op->num_even_sites*l->num_lattice_site_var;
  else
    l->sp_PRECISION.v_end = op->num_even_sites*l->num_lattice_site_var;

}

void coarse_oddeven_free_PRECISION( level_struct *l ) {
  
  int nv = l->num_parent_eig_vect;
  operator_PRECISION_struct *op = &(l->oe_op_PRECISION);

  operator_PRECISION_free( op, _ODDEVEN, l );

#ifdef HAVE_MULT_TM
  FREE( op->clover_oo_inv, complex_PRECISION, g.num_rhs_vect*SQUARE(2*nv)*op->num_odd_sites );
#else
  FREE( op->clover_oo_inv, complex_PRECISION, SQUARE(2*nv)*op->num_odd_sites );
#endif
#ifdef HAVE_TM1p1
#ifdef HAVE_MULT_TM
  FREE( op->clover_doublet_oo_inv, complex_PRECISION, g.num_rhs_vect*SQUARE(4*nv)*op->num_odd_sites );
#else
  FREE( op->clover_doublet_oo_inv, complex_PRECISION, SQUARE(4*nv)*op->num_odd_sites );
#endif
#endif
  for (int k=0; k<2; k++ )
    vector_PRECISION_free( &(op->buffer[k]), l, no_threading );
  FREE( op->buffer, vector_PRECISION, 2 );
}

void coarse_oddeven_setup_PRECISION( operator_PRECISION_struct *in, int reorder, level_struct *l, 
                                     struct Thread *threading ) {

  operator_PRECISION_struct *op = &(l->oe_op_PRECISION);

  START_LOCKED_MASTER(threading)
  int ns=l->num_inner_lattice_sites, nv = l->num_parent_eig_vect, i,
    D_size = 4*SQUARE(2*nv),
    clover_size = (nv)*(nv*2+1),
    block_size = (nv)*(nv+1);
  config_PRECISION D_in = in->D,
    clover_in = in->clover,
    odd_proj_in = in->odd_proj;

  // setup clover term and neighbor-coupling term given *in w/ or w/t reordering
  if ( reorder ) {
    int t, z, y, x, index, *le = l->local_lattice, oe_offset = op->oe_offset,
      *it = in->index_table, *dt = in->table_dim;
    config_PRECISION D_oe = op->D, 
      D_eo = (op->D)+D_size*op->num_even_sites,
      clover_ee = op->clover,
      clover_oo = (op->clover)+clover_size*op->num_even_sites,
      odd_proj_ee = op->odd_proj,
      odd_proj_oo = op->odd_proj+block_size*op->num_even_sites;

    for ( t=0; t<le[T]; t++ )
      for ( z=0; z<le[Z]; z++ )
        for ( y=0; y<le[Y]; y++ )
          for ( x=0; x<le[X]; x++ ) {
            index = site_index( t, z, y, x, dt, it );
            if ( (t+z+y+x+oe_offset)%2 == 1 ) {
              for ( i=0; i<D_size; i++ ) 
                D_eo[i] = D_in[ index*D_size+i ];
              for ( i=0; i<clover_size; i++ )
                clover_oo[i] = clover_in[ index*clover_size+i ];
              for ( i=0; i<block_size; i++ )
                odd_proj_oo[i] = odd_proj_in[ index*block_size+i ];
              D_eo += D_size;
              clover_oo += clover_size;
              odd_proj_oo += block_size;
            } else {
              for ( i=0; i<D_size; i++ )
                D_oe[i] = D_in[ index*D_size+i ];
              for ( i=0; i<clover_size; i++ )
                clover_ee[i] = clover_in[ index*clover_size+i ];
              for ( i=0; i<block_size; i++ )
                odd_proj_ee[i] = odd_proj_in[ index*block_size+i ];
              D_oe += D_size;
              clover_ee += clover_size;
              odd_proj_ee += block_size;
            }
          }
    
  } else {
    for ( i=0; i<D_size*ns; i++ )
      op->D[i] = D_in[i];
    for ( i=0; i<clover_size*ns; i++ )
      op->clover[i] = clover_in[i];
    for ( i=0; i<block_size*ns; i++ ) {
      op->odd_proj[i] = odd_proj_in[i];
    }
    
  }
  END_LOCKED_MASTER(threading)
  
  op->m0 = in->m0;

#ifdef HAVE_TM
  tm_term_PRECISION_setup( in->mu, in->mu_even_shift, in->mu_odd_shift, 1, op, l, threading );
#endif  
#ifdef HAVE_TM1p1
  epsbar_term_PRECISION_setup( in->epsbar, in->epsbar_ig5_even_shift, in->epsbar_ig5_odd_shift, op, l, threading );
#endif

  coarse_oddeven_PRECISION_set_self_couplings( l, threading );

}

void coarse_oddeven_PRECISION_set_self_couplings( level_struct *l, struct Thread *threading ) {

  operator_PRECISION_struct *op = &(l->oe_op_PRECISION);
  int nv = l->num_parent_eig_vect, start, end;

  compute_core_start_end_custom( 0, op->num_odd_sites, &start, &end, l, threading, 1);

  SYNC_CORES(threading)//debug!!!!; added because its absence was the src of OpenMP issue
#ifdef HAVE_MULT_TM
  int size = SQUARE(2*nv)*num_loop;
#else
  int size = SQUARE(2*nv);
#endif
  for( int i=start; i<end; i++ )
    coarse_selfcoupling_LU_decomposition_PRECISION( op->clover_oo_inv+i*size, op, op->num_even_sites+i, l );

#ifdef HAVE_TM1p1
#ifdef HAVE_MULT_TM
  int size_doublet = SQUARE(4*nv)*num_loop;
#else
  int size_doublet = SQUARE(4*nv);
#endif
  for( int i=start; i<end; i++ )
    coarse_selfcoupling_LU_doublet_decomposition_PRECISION( op->clover_doublet_oo_inv+i*size_doublet, op, 
                                                            op->num_even_sites+i, l );
#endif

}


static void coarse_selfcoupling_LU_decomposition_PRECISION( config_PRECISION output, operator_PRECISION_struct *op, int index, level_struct *l ) {

  // clover = [ A B      , A=A*, D=D*, C = -B*
  //            C D ]
  //
  // order: upper triangle of A, upper triangle of D, B, each column major
  //
  // tm_term = [ E 0      , E=-E*, F=-F* diag. excluded
  //             0 F ]
  //
  // order: upper triangle of E, upper triangle of F
  //
  // output = [ A+E  B   
  //             C  D+F ] LU decomposed
  // Output ordering: triu(L,1) + tril(U,0), i.e., output contains L and U without
  //                  the diagonal of L which is equal to 1

  register int i, j, k, n = l->num_parent_eig_vect, n2 = 2*n;
  config_PRECISION clover = op->clover + n*(n2+1)*index;
#ifndef HAVE_MULT_TM
  // A
  for ( j=0; j<n; j++ ) {
    for ( i=0; i<j; i++ ) {
      output[n2*i+j] = *clover;
      output[i+n2*j] = conj_PRECISION(*clover);
      clover++;      
    }
    output[(n2+1)*j] = *clover;
    clover++; // diagonal entry
  }
  // D
  for ( j=n; j<n2; j++ ) {
    for ( i=n; i<j; i++ ) {
      output[n2*i+j] = *clover;
      output[i+n2*j] = conj_PRECISION(*clover);
      clover++;      
    }
    output[(n2+1)*j] = *clover;
    clover++; // diagonal entry
  }
  // B and C
  for ( j=n; j<n2; j++ ) {
    for ( i=0; i<n; i++ ) {
      output[n2*i+j] = *clover;
      output[i+n2*j] = -conj_PRECISION(*clover);
      clover++;      
    }
  }

#ifdef HAVE_TM // In this case, we compute tm_term using avg even shifts at coarsest & smallest shift at intermediate levels
  double even = (l->level == 0 && 0)?op->even_shift_avg: op->mu_even_shift[0];
#if 1 // In this case, we compute tm_term using avg even shifts at coarsest & smallest shift at intermediate levels
  complex_PRECISION mu_even = (complex_PRECISION) I*(op->mu+even);
  complex_PRECISION odd_factor = (complex_PRECISION) I*(op->mu_odd_shift - even);
  config_PRECISION odd_proj = op->odd_proj+n*(n+1)*index;
  if (op->mu + op->mu_odd_shift != 0.0 || op->is_even_shifted_mu_nonzero ) {
    // E
    for ( j=0; j<n; j++ ) {
      for ( i=0; i<j; i++ ) {
        output[n2*i+j] += -odd_factor*(*odd_proj);
        output[i+n2*j] += -conj_PRECISION(-odd_factor*(*odd_proj));
        odd_proj++;      
      }
      output[(n2+1)*j] += -1.* ( mu_even + odd_factor * (*odd_proj));
      odd_proj++; // diagonal entry
    }
    // F
    for ( j=n; j<n2; j++ ) {
      for ( i=n; i<j; i++ ) {
        output[n2*i+j] += odd_factor * (*odd_proj);
        output[i+n2*j] += -conj_PRECISION(odd_factor * (*odd_proj));
        odd_proj++;      
      }
      output[(n2+1)*j] += mu_even + odd_factor * (*odd_proj);
      odd_proj++; // diagonal entry
    }
  }
#else
  complex_PRECISION odd = (complex_PRECISION) I*op->odd_shifted_mu, odd_factor = (complex_PRECISION) I*(op->mu_odd_shift - even);
  config_PRECISION odd_proj = op->odd_proj+n*(n+1)*index;
  if (op->mu + op->mu_odd_shift != 0.0 || op->is_even_shifted_mu_nonzero ) {
    // E
    for ( j=0; j<n; j++ ) {
      for ( i=0; i<j; i++ ) {
        output[n2*i+j] += -odd_factor*(*odd_proj);
        output[i+n2*j] += -conj_PRECISION(-odd_factor*(*odd_proj));
        odd_proj++;      
      }
      output[(n2+1)*j] += -1.* ( odd - odd_factor * (1 - *odd_proj));
      odd_proj++; // diagonal entry
    }
    // F
    for ( j=n; j<n2; j++ ) {
      for ( i=n; i<j; i++ ) {
        output[n2*i+j] += odd_factor * (*odd_proj);
        output[i+n2*j] += -conj_PRECISION(odd_factor * (*odd_proj));
        odd_proj++;      
      }
      output[(n2+1)*j] += odd - odd_factor * (1 - *odd_proj);
      odd_proj++; // diagonal entry
    }
  }
#endif
#endif
    
  // compute LU decomposition
  // output = triu(L,1) + tril(U,0)
  // i.e., output contains L and U without the diagonal of L which is equal to 1
  // order: row major
  for ( k=0; k<n2; k++ ) {
    for ( i=k+1; i<n2; i++ ) {
      output[n2*i+k] = output[n2*i+k]/output[(n2+1)*k]; // output(i,k) = output(i,k)/output(k,k)
      for ( j=k+1; j<n2; j++ )
        output[n2*i+j] = output[n2*i+j]-output[n2*i+k]*output[n2*k+j]; // output(i,j) = output(i,j)-output(i,k)*output(k,j)
    }
  }
  
#else
  register int jj, jjj, block_size = n*(n+1), nv = l->num_inner_lattice_sites*block_size, nc = SQUARE(n2)*op->num_odd_sites, nrt = (g.in_setup)?num_loop:g.num_rhs_vect;
  
  // A
  for ( j=0; j<n; j++ ) {
    for ( i=0; i<j; i++ ) {
      VECTOR_LOOP(jj, nrt, jjj, output[(n2*i+j)*num_loop+jj*nc+jjj] = *clover; )
      VECTOR_LOOP(jj, nrt, jjj, output[(i+n2*j)*num_loop+jj*nc+jjj] = conj_PRECISION(*clover);)
      clover++;      
    }
    VECTOR_LOOP(jj, nrt, jjj, output[(n2+1)*j*num_loop+jj*nc+jjj] = *clover;)
    clover++; // diagonal entry
  }
  // D
  for ( j=n; j<n2; j++ ) {
    for ( i=n; i<j; i++ ) {
      VECTOR_LOOP(jj, nrt, jjj, output[(n2*i+j)*num_loop+jj*nc+jjj] = *clover;)
      VECTOR_LOOP(jj, nrt, jjj, output[(i+n2*j)*num_loop+jj*nc+jjj] = conj_PRECISION(*clover);)
      clover++;      
    }
    VECTOR_LOOP(jj, nrt, jjj, output[(n2+1)*j*num_loop+jj*nc+jjj] = *clover; )
    clover++; // diagonal entry
  }
  // B and C
  for ( j=n; j<n2; j++ ) {
    for ( i=0; i<n; i++ ) {
      VECTOR_LOOP(jj, nrt, jjj, output[(n2*i+j)*num_loop+jj*nc+jjj] = *clover;)
      VECTOR_LOOP(jj, nrt, jjj, output[(i+n2*j)*num_loop+jj*nc+jjj] = -conj_PRECISION(*clover);)
      clover++;      
    }
  }
  
  config_PRECISION tm_term = op->tm_term + n*(n+1)*index*num_loop;
  if (op->mu + op->mu_odd_shift != 0.0 || op->is_even_shifted_mu_nonzero ) {
    // E
    for ( j=0; j<n; j++ ) {
      for ( i=0; i<j; i++ ) {
        VECTOR_LOOP(jj, nrt, jjj, output[(n2*i+j)*num_loop+jj*nc+jjj] += tm_term[jj*nv+jjj];)
	VECTOR_LOOP(jj, nrt, jjj, output[(i+n2*j)*num_loop+jj*nc+jjj] += -conj_PRECISION(tm_term[jj*nv+jjj]); )
        tm_term+=num_loop;      
      }
      VECTOR_LOOP(jj, nrt, jjj, output[(n2+1)*j*num_loop+jj*nc+jjj] += tm_term[jj*nv+jjj]; );
      tm_term+=num_loop; // diagonal entry
    }
    // F
    for ( j=n; j<n2; j++ ) {
      for ( i=n; i<j; i++ ) {
        VECTOR_LOOP(jj, nrt, jjj, output[(n2*i+j)*num_loop+jj*nc+jjj] += tm_term[jj*nv+jjj]; )
	VECTOR_LOOP(jj, nrt, jjj, output[(i+n2*j)*num_loop+jj*nc+jjj] += -conj_PRECISION(tm_term[jj*nv+jjj]); )
        tm_term+=num_loop; 
      }
      VECTOR_LOOP(jj, nrt, jjj, output[(n2+1)*j*num_loop+jj*nc+jjj] += tm_term[jj*nv+jjj]; )
      tm_term+=num_loop; // diagonal entry
    }
  }

  // compute LU decomposition
  // output = triu(L,1) + tril(U,0)
  // i.e., output contains L and U without the diagonal of L which is equal to 1
  // order: row major
  for ( k=0; k<n2; k++ ) {
    for ( i=k+1; i<n2; i++ ) {
      // output(i,k) = output(i,k)/output(k,k)
      VECTOR_LOOP(jj, nrt, jjj, output[(n2*i+k)*num_loop+jj*nc+jjj] = output[(n2*i+k)*num_loop+jj*nc+jjj]/output[(n2+1)*k*num_loop+jj*nc+jjj];) 
      for ( j=k+1; j<n2; j++ )
	// output(i,j) = output(i,j)-output(i,k)*output(k,j)
        VECTOR_LOOP(jj, nrt, jjj,
		    output[(n2*i+j)*num_loop+jj*nc+jjj] = output[(n2*i+j)*num_loop+jj*nc+jjj]-output[(n2*i+k)*num_loop+jj*nc+jjj]*output[(n2*k+j)*num_loop+jj*nc+jjj]; )
    }
  }
#endif
}

#ifdef HAVE_TM1p1
static void coarse_selfcoupling_LU_doublet_decomposition_PRECISION( config_PRECISION output, operator_PRECISION_struct *op, int index, level_struct *l ) {

  // clover = [ A B      , A=A*, D=D*, C = -B*
  //            C D ]
  //
  // order: upper triangle of A, upper triangle of D, B, each column major
  //
  // tm_term = [ E 0      , E = -E*, F=-F* diag. excluded
  //             0 F ]    
  //
  // order: upper triangle of E, upper triangle of F
  //
  // epsbar_term = [ G 0      , G=-G*, H=-H* diag. excluded
  //                 0 H ]    , block diagonal in spinor index & block off-diagonal in flavor index
  //
  // order: upper triangle of G, upper triangle of H
  //
  // Input = [ 
  // output = [ A+E  G   B   0
  //             G  A-E  0   B
  //             C   0  D+F  H
  //             0   C   H  D-F ]  LU decomposed

  register int i, j, k, n = l->num_parent_eig_vect, n2 = 2*n, n3 = 3*n, n4 = 4*n;
#ifdef  HAVE_MULT_TM
  register int jj, jjj, block_size = n*(n+1), nv = l->num_inner_lattice_sites*block_size, nc = SQUARE(n4)*op->num_odd_sites, nrt = (g.in_setup)?num_loop:g.num_rhs_vect;

  // initialize the matrix to 0
  for ( j=0; j<n4; j++ )
    for ( i=0; i<n4; i++ ) 
      VECTOR_LOOP(jj, nrt, jjj, output[(i*n4+j)*num_loop+jj*nc+jjj] = _COMPLEX_PRECISION_ZERO;)

  config_PRECISION clover = op->clover + n*(n2+1)*index;
  // A
  for ( j=0; j<n; j++ ) {
    for ( i=0; i<j; i++ ) {
      VECTOR_LOOP(jj, nrt, jjj, output[2*(n4*i+j)*num_loop+jj*nc+jjj] = output[(n4*(2*i+1)+2*j+1)*num_loop+jj*nc+jjj] = *clover;)
      VECTOR_LOOP(jj, nrt, jjj, output[2*(n4*j+i)*num_loop+jj*nc+jjj] = output[(n4*(2*j+1)+2*i+1)*num_loop+jj*nc+jjj] = conj_PRECISION(*clover);)
      clover++;      
    }
    // diagonal entry
    VECTOR_LOOP(jj, nrt, jjj, output[((n4+1)*2*j)*num_loop+jj*nc+jjj] = output[((n4+1)*(2*j+1))*num_loop+jj*nc+jjj] = *clover;)
    clover++; 
  }
  // D
  for ( j=n; j<n2; j++ ) {
    for ( i=n; i<j; i++ ) {
      VECTOR_LOOP(jj, nrt, jjj, output[2*(n4*i+j)*num_loop+jj*nc+jjj] = output[(n4*(2*i+1)+2*j+1)*num_loop+jj*nc+jjj] = *clover;)
      VECTOR_LOOP(jj, nrt, jjj, output[2*(n4*j+i)*num_loop+jj*nc+jjj] = output[(n4*(2*j+1)+2*i+1)*num_loop+jj*nc+jjj] = conj_PRECISION(*clover);)
      clover++;      
    }
    // diagonal entry
    VECTOR_LOOP(jj, nrt, jjj, output[((n4+1)*2*j)*num_loop+jj*nc+jjj] = output[((n4+1)*(2*j+1))*num_loop+jj*nc+jjj] = *clover;)
    clover++; 
  }
  // B and C
  for ( j=n; j<n2; j++ ) {
    for ( i=0; i<n; i++ ) {
      VECTOR_LOOP(jj, nrt, jjj, output[2*(n4*i+j)*num_loop+jj*nc+jjj] = output[(n4*(2*i+1)+2*j+1)*num_loop+jj*nc+jjj] = *clover;)
      VECTOR_LOOP(jj, nrt, jjj, output[2*(n4*j+i)*num_loop+jj*nc+jjj] = output[(2*i+1+n4*(2*j+1))*num_loop+jj*nc+jjj] = -conj_PRECISION(*clover);)
      clover++;      
    }
  }

  config_PRECISION tm_term = op->tm_term + n*(n+1)*index*num_loop;
  if (op->mu + op->mu_odd_shift != 0.0 || op->is_even_shifted_mu_nonzero ) {
    // E
    for ( j=0; j<n; j++ ) {
      for ( i=0; i<j; i++ ) {
	VECTOR_LOOP(jj, nrt, jjj, output[2*(n4*i+j)*num_loop+jj*nc+jjj] += tm_term[jj*nv+jjj];)
        VECTOR_LOOP(jj, nrt, jjj, output[(n4*(2*i+1)+2*j+1)*num_loop+jj*nc+jjj] -= tm_term[jj*nv+jjj];)
        VECTOR_LOOP(jj, nrt, jjj, output[2*(n4*j+i)*num_loop+jj*nc+jjj] += -conj_PRECISION(tm_term[jj*nv+jjj]);)
        VECTOR_LOOP(jj, nrt, jjj, output[(n4*(2*j+1)+2*i+1)*num_loop+jj*nc+jjj] -= -conj_PRECISION(tm_term[jj*nv+jjj]);)
        tm_term+=num_loop;
      }
      // diagonal entry
      VECTOR_LOOP(jj, nrt, jjj, output[((n4+1)*2*j)*num_loop+jj*nc+jjj] += tm_term[jj*nv+jjj];)
      VECTOR_LOOP(jj, nrt, jjj, output[((n4+1)*(2*j+1))*num_loop+jj*nc+jjj] -= tm_term[jj*nv+jjj];)
      tm_term+=num_loop; 
    }
    // F
    for ( j=n; j<n2; j++ ) {
      for ( i=n; i<j; i++ ) {
	VECTOR_LOOP(jj, nrt, jjj, output[2*(n4*i+j)*num_loop+jj*nc+jjj] += tm_term[jj*nv+jjj];)
	VECTOR_LOOP(jj, nrt, jjj, output[(n4*(2*i+1)+2*j+1)*num_loop+jj*nc+jjj] -= tm_term[jj*nv+jjj];)
	VECTOR_LOOP(jj, nrt, jjj, output[2*(n4*j+i)*num_loop+jj*nc+jjj] += -conj_PRECISION(tm_term[jj*nv+jjj]);)
	VECTOR_LOOP(jj, nrt, jjj, output[(n4*(2*j+1)+2*i+1)*num_loop+jj*nc+jjj] -= -conj_PRECISION(tm_term[jj*nv+jjj]);)
        tm_term+=num_loop;
      }
      // diagonal entry
      VECTOR_LOOP(jj, nrt, jjj, output[((n4+1)*2*j)*num_loop+jj*nc+jjj] += tm_term[jj*nv+jjj];)
      VECTOR_LOOP(jj, nrt, jjj, output[((n4+1)*(2*j+1))*num_loop+jj*nc+jjj] -= tm_term[jj*nv+jjj];)
      tm_term+=num_loop; 
    }
  }

  config_PRECISION epsbar_term = op->epsbar_term + n*(n+1)*index;
  // G
  for ( j=0; j<n; j++ ) {
    for ( i=0; i<j; i++ ) {
      VECTOR_LOOP(jj, nrt, jjj, output[(n4*2*i+2*j+1)*num_loop+jj*nc+jjj] = output[(n4*(2*i+1)+2*j)*num_loop+jj*nc+jjj] = (*epsbar_term);)
      VECTOR_LOOP(jj, nrt, jjj, output[(2*i+1+n4*2*j)*num_loop+jj*nc+jjj] = output[(n4*(2*j+1)+2*i)*num_loop+jj*nc+jjj] = -conj_PRECISION(*epsbar_term);)
      epsbar_term++;      
    }
    // diagonal entry
    VECTOR_LOOP(jj, nrt, jjj, output[((n4+1)*2*j+1)*num_loop+jj*nc+jjj] = output[((n4+1)*2*j+n4)*num_loop+jj*nc+jjj] = (*epsbar_term);)
    epsbar_term++;
  }
  // H
  for ( j=n; j<n2; j++ ) {
    for ( i=n; i<j; i++ ) {
      VECTOR_LOOP(jj, nrt, jjj, output[(n4*2*i+2*j+1)*num_loop+jj*nc+jjj] = output[(n4*(2*i+1)+2*j)*num_loop+jj*nc+jjj] = (*epsbar_term);)
      VECTOR_LOOP(jj, nrt, jjj, output[(2*i+1+n4*2*j)*num_loop+jj*nc+jjj] = output[(n4*(2*j+1)+2*i)*num_loop+jj*nc+jjj] = -conj_PRECISION(*epsbar_term);)
      epsbar_term++;      
    }
    // diagonal entry
    VECTOR_LOOP(jj, nrt, jjj, output[((n4+1)*2*j+1)*num_loop+jj*nc+jjj] = output[((n4+1)*2*j+n4)*num_loop+jj*nc+jjj] = (*epsbar_term);)
    epsbar_term++; 
  }
    
  // compute LU decomposition
  // output = triu(L,1) + tril(U,0), i.e., lower strictly-triangular part of output is set to triu(L,1), and upper triangular part to tril(U,0)
  // i.e., output contains L and U without the diagonal of L which is equal to 1
  // order: row major
  for ( k=0; k<n4; k++ ) {//for each column
    for ( i=k+1; i<n4; i++ ) { // for L
      VECTOR_LOOP(jj, nrt, jjj,
		  output[(n4*i+k)*num_loop+jj*nc+jjj] = output[(n4*i+k)*num_loop+jj*nc+jjj]/output[((n4+1)*k)*num_loop+jj*nc+jjj];) // output(i,k) = output(i,k)/output(k,k)
      // for U
      for ( j=k+1; j<n4; j++ )// output(i,j) = output(i,j)-output(i,k)*output(k,j)
	VECTOR_LOOP(jj, nrt, jjj,
		  output[(n4*i+j)*num_loop+jj*nc+jjj] = output[(n4*i+j)*num_loop+jj*nc+jjj]-output[(n4*i+k)*num_loop+jj*nc+jjj]*output[(n4*k+j)*num_loop+jj*nc+jjj];) 
    }
  }
#else
  // initialize the matrix to 0
  for ( j=0; j<n4; j++ )
    for ( i=0; i<n4; i++ ) 
      output[(i*n4+j)] = _COMPLEX_PRECISION_ZERO;

  config_PRECISION clover = op->clover + n*(n2+1)*index;
  // A
  for ( j=0; j<n; j++ ) {
    for ( i=0; i<j; i++ ) {
      output[2*(n4*i+j)*num_loop+jj*nc+jjj] = output[(n4*(i*2+1)+2*j+1)] = *clover;
      output[2*(i+n4*j)] = output[(2*i+1+n4*(2*j+1))] = conj_PRECISION(*clover);
      clover++;      
    }
    // diagonal entry
    VECTOR_LOOP(jj, nrt, jjj, output[((n4+1)*2*j)*num_loop+jj*nc+jjj] = output[((n4+1)*(2*j+1))*num_loop+jj*nc+jjj] = *clover;)
    clover++; 
  }
  // D
  for ( j=n; j<n2; j++ ) {
    for ( i=n; i<j; i++ ) {
      VECTOR_LOOP(jj, nrt, jjj, output[2*(n4*i+j)*num_loop+jj*nc+jjj] = output[(n4*(i*2+1)+2*j+1)*num_loop+jj*nc+jjj] = *clover;)
      VECTOR_LOOP(jj, nrt, jjj, output[2*(i+n4*j)*num_loop+jj*nc+jjj] = output[(2*i+1+n4*(2*j+1))*num_loop+jj*nc+jjj] = conj_PRECISION(*clover);)
      clover++;      
    }
    // diagonal entry
    VECTOR_LOOP(jj, nrt, jjj, output[((n4+1)*2*j)*num_loop+jj*nc+jjj] = output[((n4+1)*(2*j+1))*num_loop+jj*nc+jjj] = *clover;)
    clover++; 
  }
  // B and C
  for ( j=n; j<n2; j++ ) {
    for ( i=0; i<n; i++ ) {
      VECTOR_LOOP(jj, nrt, jjj, output[2*(n4*i+j)*num_loop+jj*nc+jjj] = output[(n4*(i*2+1)+2*j+1)*num_loop+jj*nc+jjj] = *clover;)
      VECTOR_LOOP(jj, nrt, jjj, output[2*(i+n4*j)*num_loop+jj*nc+jjj] = output[(2*i+1+n4*(2*j+1))*num_loop+jj*nc+jjj] = -conj_PRECISION(*clover);)
      clover++;      
    }
  }

  config_PRECISION tm_term = op->tm_term + n*(n+1)*index*num_loop;
  if (op->mu + op->mu_odd_shift != 0.0 || op->is_even_shifted_mu_nonzero ) {
    // E
    for ( j=0; j<n; j++ ) {
      for ( i=0; i<j; i++ ) {
	VECTOR_LOOP(jj, nrt, jjj, output[2*(n4*i+j)*num_loop+jj*nc+jjj] += tm_term[jj*nv+jjj];)
        VECTOR_LOOP(jj, nrt, jjj, output[(n4*(i*2+1)+2*j+1)*num_loop+jj*nc+jjj] -= tm_term[jj*nv+jjj];)
        VECTOR_LOOP(jj, nrt, jjj, output[2*(n4*j+i)*num_loop+jj*nc+jjj] += conj_PRECISION(tm_term[jj*nv+jjj]);)
        VECTOR_LOOP(jj, nrt, jjj, output[(n4*(j*2+1)+2*i+1)*num_loop+jj*nc+jjj] -= conj_PRECISION(tm_term[jj*nv+jjj]);)
        tm_term+=num_loop;
      }
      // diagonal entry
      VECTOR_LOOP(jj, nrt, jjj, output[((n4+1)*2*j)*num_loop+jj*nc+jjj] += tm_term[jj*nv+jjj];)
      VECTOR_LOOP(jj, nrt, jjj, output[((n4+1)*2*j+1+n4)*num_loop+jj*nc+jjj] -= tm_term[jj*nv+jjj];)
      tm_term+=num_loop; 
    }
    // F
    for ( j=n; j<n2; j++ ) {
      for ( i=n; i<j; i++ ) {
	VECTOR_LOOP(jj, nrt, jjj, output[2*(n4*i+j)*num_loop+jj*nc+jjj] += tm_term[jj*nv+jjj];)
	VECTOR_LOOP(jj, nrt, jjj, output[(n4*(i*2+1)+2*j+1)*num_loop+jj*nc+jjj] -= tm_term[jj*nv+jjj];)
	VECTOR_LOOP(jj, nrt, jjj, output[2*(n4*j+i)*num_loop+jj*nc+jjj] += conj_PRECISION(tm_term[jj*nv+jjj]);)
	VECTOR_LOOP(jj, nrt, jjj, output[(n4*(j*2+1)+2*i+1)*num_loop+jj*nc+jjj] -= conj_PRECISION(tm_term[jj*nv+jjj]);)
        tm_term+=num_loop;
      }
      // diagonal entry
      VECTOR_LOOP(jj, nrt, jjj, output[((n4+1)*2*j)*num_loop+jj*nc+jjj] += tm_term[jj*nv+jjj];)
      VECTOR_LOOP(jj, nrt, jjj, output[((n4+1)*2*j+1+n4)*num_loop+jj*nc+jjj] -= tm_term[jj*nv+jjj];)
      tm_term+=num_loop; 
    }
  }

  config_PRECISION epsbar_term = op->epsbar_term + n*(n+1)*index;
  // G
  for ( j=0; j<n; j++ ) {
    for ( i=0; i<j; i++ ) {
      VECTOR_LOOP(jj, nrt, jjj, output[(n4*2*i+2*j+1)*num_loop+jj*nc+jjj] = output[(n4*(2*i+1)+2*j)*num_loop+jj*nc+jjj] = (*epsbar_term);)
      VECTOR_LOOP(jj, nrt, jjj, output[(2*i+1+n4*2*j)*num_loop+jj*nc+jjj] = output[(n4*(2*j+1)+2*i)*num_loop+jj*nc+jjj] = -conj_PRECISION(*epsbar_term);)
      epsbar_term++;      
    }
    // diagonal entry
    VECTOR_LOOP(jj, nrt, jjj, output[((n4+1)*2*j+1)*num_loop+jj*nc+jjj] = output[((n4+1)*2*j+n4)*num_loop+jj*nc+jjj] (*epsbar_term);)
    epsbar_term++;
  }
  // H
  for ( j=n; j<n2; j++ ) {
    for ( i=n; i<j; i++ ) {
      VECTOR_LOOP(jj, nrt, jjj, output[(n4*2*i+2*j+1)*num_loop+jj*nc+jjj] = output[(n4*(2*i+1)+2*j)*num_loop+jj*nc+jjj] = (*epsbar_term);)
      VECTOR_LOOP(jj, nrt, jjj, output[(2*i+1+n4*2*j)*num_loop+jj*nc+jjj] = output[(n4*(2*j+1)+2*i)*num_loop+jj*nc+jjj] = -conj_PRECISION(*epsbar_term);)
      epsbar_term++;      
    }
    // diagonal entry
    VECTOR_LOOP(jj, nrt, jjj, output[((n4+1)*2*j+1)*num_loop+jj*nc+jjj] = output[((n4+1)*2*j+n4)*num_loop+jj*nc+jjj] (*epsbar_term);)
    epsbar_term++; 
  }
    
  // compute LU decomposition
  // output = triu(L,1) + tril(U,0), i.e., lower strictly-triangular part of output is set to triu(L,1), and upper triangular part to tril(U,0)
  // i.e., output contains L and U without the diagonal of L which is equal to 1
  // order: row major
  for ( k=0; k<n4; k++ ) {//for each column
    for ( i=k+1; i<n4; i++ ) { // for L
      VECTOR_LOOP(jj, nrt, jjj,
		  output[(n4*i+k)*num_loop+jj*nc+jjj] = output[(n4*i+k)*num_loop+jj*nc+jjj]/output[((n4+1)*k)*num_loop+jj*nc+jjj];) // output(i,k) = output(i,k)/output(k,k)
      // for U
      for ( j=k+1; j<n4; j++ )// output(i,j) = output(i,j)-output(i,k)*output(k,j)
	VECTOR_LOOP(jj, nrt, jjj,
		  output[(n4*i+j)*num_loop+jj*nc+jjj] = output[(n4*i+j)*num_loop+jj*nc+jjj]-output[(n4*i+k)*num_loop+jj*nc+jjj]*output[(n4*k+j)*num_loop+jj*nc+jjj];) 
    }
  }
  
  //////////0
  for ( j=0; j<n; j++ ) {
    for ( i=0; i<n; i++ ) {
      output[n4*(i+0 )+(j+n3)] = _COMPLEX_PRECISION_ZERO;
      output[n4*(i+n )+(j+n2)] = _COMPLEX_PRECISION_ZERO;
      output[n4*(i+n2)+(j+n )] = _COMPLEX_PRECISION_ZERO;
      output[n4*(i+n3)+(j+0 )] = _COMPLEX_PRECISION_ZERO;
    }
  }

  config_PRECISION clover = op->clover + n*(n2+1)*index;
  // A
  for ( j=0; j<n; j++ ) {
    for ( i=0; i<j; i++ ) {
      output[n4*i+j] = *clover;
      output[i+n4*j] = conj_PRECISION(*clover);
      output[n4*(i+n)+(j+n)] = *clover;
      output[(i+n)+n4*(j+n)] = conj_PRECISION(*clover);
      clover++;      
    }
    output[(n4+1)*j] = *clover;
    output[(n4+1)*(j+n)] = *clover;
    clover++; // diagonal entry
  }
  // D
  for ( j=n2; j<n3; j++ ) {
    for ( i=n2; i<j; i++ ) {
      output[n4*i+j] = *clover;
      output[i+n4*j] = conj_PRECISION(*clover);
      output[n4*(i+n)+(j+n)] = *clover;
      output[(i+n)+n4*(j+n)] = conj_PRECISION(*clover);
      clover++;      
    }
    output[(n4+1)*j] = *clover;
    output[(n4+1)*(j+n)] = *clover;
    clover++; // diagonal entry
  }
  // B and C
  for ( j=n2; j<n3; j++ ) {
    for ( i=0; i<n; i++ ) {
      output[n4*i+j] = *clover;
      output[i+n4*j] = -conj_PRECISION(*clover);
      output[n4*(i+n)+(j+n)] = *clover;
      output[(i+n)+n4*(j+n)] = -conj_PRECISION(*clover);
      clover++;      
    }
  }

#ifdef HAVE_TM
  double even = (l->level == 0 && 0)?op->even_shift_avg: op->mu_even_shift[0];
  complex_PRECISION mu_even = (complex_PRECISION) I*(op->mu+even);
  complex_PRECISION odd_factor = (complex_PRECISION) I*(op->mu_odd_shift - even);
  config_PRECISION odd_proj = op->odd_proj+n*(n+1)*index;

  if (op->mu + op->mu_odd_shift != 0.0 || op->is_even_shifted_mu_nonzero ) {
    // E
    for ( j=0; j<n; j++ ) {
      for ( i=0; i<j; i++ ) {
        output[n4*i+j] += -odd_factor * (*odd_proj);
        output[i+n4*j] += -conj_PRECISION(-odd_factor * (*odd_proj));
        output[n4*(i+n)+(j+n)] -= -odd_factor * (*odd_proj);
        output[(i+n)+n4*(j+n)] -= -conj_PRECISION(-odd_factor * (*odd_proj));
        odd_proj++;
      }
      output[(n4+1)*j] += -1.* ( mu_even + odd_factor * (*odd_proj));
      output[(n4+1)*(j+n)] -= -1.* ( mu_even + odd_factor * (*odd_proj));
      odd_proj++; // diagonal entry
    }
    // F
    for ( j=n2; j<n3; j++ ) {
      for ( i=n2; i<j; i++ ) {
        output[n4*i+j] += odd_factor * (*odd_proj);
        output[i+n4*j] += -conj_PRECISION(odd_factor * (*odd_proj));
        output[n4*(i+n)+(j+n)] -= odd_factor * (*odd_proj);
        output[(i+n)+n4*(j+n)] -= -conj_PRECISION(odd_factor * (*odd_proj));
        odd_proj++;      
      }
      output[(n4+1)*j] += mu_even + odd_factor * (*odd_proj);
      output[(n4+1)*(j+n)] -= mu_even + odd_factor * (*odd_proj);
      odd_proj++; // diagonal entry
    }
  }
#endif

  config_PRECISION epsbar_term = op->epsbar_term + n*(n+1)*index;
  // G
  for ( j=n; j<n2; j++ ) {
    for ( i=0; i<(j-n); i++ ) {
      output[n4*i+j] = (*epsbar_term);
      output[(i+n)+n4*(j-n)] = -conj_PRECISION(*epsbar_term);
      output[n4*(i+n)+(j-n)] = (*epsbar_term);
      output[i+n4*j] = -conj_PRECISION(*epsbar_term);
      epsbar_term++;      
    }
    output[(n4+1)*(j-n)+n] = (*epsbar_term);
    output[(n4+1)*j-n] = (*epsbar_term);
    epsbar_term++; // diagonal entry
  }
  // H
  for ( j=n3; j<n4; j++ ) {
    for ( i=n2; i<(j-n); i++ ) {
      output[n4*i+j] = (*epsbar_term);
      output[(i+n)+n4*(j-n)] = -conj_PRECISION(*epsbar_term);
      output[n4*(i+n)+(j-n)] = (*epsbar_term);
      output[i+n4*j] = -conj_PRECISION(*epsbar_term);
      epsbar_term++;      
    }
    output[(n4+1)*(j-n)+n] = (*epsbar_term);
    output[(n4+1)*j-n] = (*epsbar_term);
    epsbar_term++; // diagonal entry
  }
    
  // compute LU decomposition
  // output = triu(L,1) + tril(U,0)
  // i.e., output contains L and U without the diagonal of L which is equal to 1
  // order: row major
  for ( k=0; k<n4; k++ ) {
    for ( i=k+1; i<n4; i++ ) {
      output[n4*i+k] = output[n4*i+k]/output[(n4+1)*k]; // output(i,k) = output(i,k)/output(k,k)
      for ( j=k+1; j<n4; j++ )
        output[n4*i+j] = output[n4*i+j]-output[n4*i+k]*output[n4*k+j]; // output(i,j) = output(i,j)-output(i,k)*output(k,j)
    }
  }
  complex_PRECISION tmp;
  for ( i=0; i<n4; i++ ) {
    for ( j=0; j<n; j++ ) {
      tmp = output[n4*i+n+j];
      output[n4*i+n+j] = output[n4*i+n2+j];
      output[n4*i+n2+j] = tmp;
    }
  }
#endif
  
#if 0
#ifdef HAVE_MULT_TM
  register int jj, jjj, block_size = n*(n+1), nv = l->num_inner_lattice_sites*block_size, nc = SQUARE(n4)*op->num_odd_sites, nrt = (g.in_setup)?num_loop:g.num_rhs_vect;
  // off-diagonal zero parts
  // 0
  for ( j=0; j<n; j++ ) {
    for ( i=0; i<n; i++ ) {
      VECTOR_LOOP(jj, nrt, jjj, output[(n4*(i+0 )+(j+n3))*num_loop+jj*nc+jjj] = _COMPLEX_PRECISION_ZERO;)
      VECTOR_LOOP(jj, nrt, jjj, output[(n4*(i+n )+(j+n2))*num_loop+jj*nc+jjj] = _COMPLEX_PRECISION_ZERO;)
      VECTOR_LOOP(jj, nrt, jjj, output[(n4*(i+n2)+(j+n ))*num_loop+jj*nc+jjj] = _COMPLEX_PRECISION_ZERO;)
      VECTOR_LOOP(jj, nrt, jjj, output[(n4*(i+n3)+(j+0 ))*num_loop+jj*nc+jjj] = _COMPLEX_PRECISION_ZERO;)
    }
  }

  config_PRECISION clover = op->clover + n*(n2+1)*index;
  // A
  for ( j=0; j<n; j++ ) {
    for ( i=0; i<j; i++ ) {
      VECTOR_LOOP(jj, nrt, jjj, output[(n4*i+j)*num_loop+jj*nc+jjj] = *clover;)
      VECTOR_LOOP(jj, nrt, jjj, output[(i+n4*j)*num_loop+jj*nc+jjj] = conj_PRECISION(*clover);)
      VECTOR_LOOP(jj, nrt, jjj, output[(n4*(i+n)+(j+n))*num_loop+jj*nc+jjj] = *clover;)
      VECTOR_LOOP(jj, nrt, jjj, output[((i+n)+n4*(j+n))*num_loop+jj*nc+jjj] = conj_PRECISION(*clover);)
      clover++;      
    }
    VECTOR_LOOP(jj, nrt, jjj, output[((n4+1)*j)*num_loop+jj*nc+jjj] = *clover;)
    VECTOR_LOOP(jj, nrt, jjj, output[((n4+1)*(j+n))*num_loop+jj*nc+jjj] = *clover;)
    clover++; // diagonal entry
  }
  // D
  for ( j=n2; j<n3; j++ ) {
    for ( i=n2; i<j; i++ ) {
      VECTOR_LOOP(jj, nrt, jjj, output[(n4*i+j)*num_loop+jj*nc+jjj] = *clover;)
      VECTOR_LOOP(jj, nrt, jjj, output[(i+n4*j)*num_loop+jj*nc+jjj] = conj_PRECISION(*clover);)
      VECTOR_LOOP(jj, nrt, jjj, output[(n4*(i+n)+(j+n))*num_loop+jj*nc+jjj] = *clover;)
      VECTOR_LOOP(jj, nrt, jjj, output[((i+n)+n4*(j+n))*num_loop+jj*nc+jjj] = conj_PRECISION(*clover);)
      clover++;      
    }
    VECTOR_LOOP(jj, nrt, jjj, output[((n4+1)*j)*num_loop+jj*nc+jjj] = *clover;)
    VECTOR_LOOP(jj, nrt, jjj, output[((n4+1)*(j+n))*num_loop+jj*nc+jjj] = *clover;)
    clover++; // diagonal entry
  }
  // B and C
  for ( j=n2; j<n3; j++ ) {
    for ( i=0; i<n; i++ ) {
      VECTOR_LOOP(jj, nrt, jjj, output[(n4*i+j)*num_loop+jj*nc+jjj] = *clover;)
      VECTOR_LOOP(jj, nrt, jjj, output[(i+n4*j)*num_loop+jj*nc+jjj] = -conj_PRECISION(*clover);)
      VECTOR_LOOP(jj, nrt, jjj, output[(n4*(i+n)+(j+n))*num_loop+jj*nc+jjj] = *clover;)
      VECTOR_LOOP(jj, nrt, jjj, output[((i+n)+n4*(j+n))*num_loop+jj*nc+jjj] = -conj_PRECISION(*clover);)
      clover++;      
    }
  }

  config_PRECISION tm_term = op->tm_term + n*(n+1)*index*num_loop;
  if (op->mu + op->mu_odd_shift != 0.0 || op->is_even_shifted_mu_nonzero ) {
    // E
    for ( j=0; j<n; j++ ) {
      for ( i=0; i<j; i++ ) {
        VECTOR_LOOP(jj, nrt, jjj, output[(n4*i+j)*num_loop+jj*nc+jjj] += tm_term[jj*nv+jjj];)
	VECTOR_LOOP(jj, nrt, jjj, output[(i+n4*j)*num_loop+jj*nc+jjj] += -conj_PRECISION(tm_term[jj*nv+jjj]);)
	VECTOR_LOOP(jj, nrt, jjj, output[(n4*(i+n)+(j+n))*num_loop+jj*nc+jjj] -= tm_term[jj*nv+jjj];)
	VECTOR_LOOP(jj, nrt, jjj, output[((i+n)+n4*(j+n))*num_loop+jj*nc+jjj] -= -conj_PRECISION(tm_term[jj*nv+jjj]);)
        tm_term+=num_loop;
      }
      VECTOR_LOOP(jj, nrt, jjj, output[((n4+1)*j)*num_loop+jj*nc+jjj] += tm_term[jj*nv+jjj];)
      VECTOR_LOOP(jj, nrt, jjj, output[((n4+1)*(j+n))*num_loop+jj*nc+jjj] -= tm_term[jj*nv+jjj];)
      tm_term+=num_loop; // diagonal entry
    }
    // F
    for ( j=n2; j<n3; j++ ) {
      for ( i=n2; i<j; i++ ) {
        VECTOR_LOOP(jj, nrt, jjj, output[(n4*i+j)*num_loop+jj*nc+jjj] += tm_term[jj*nv+jjj];)
	VECTOR_LOOP(jj, nrt, jjj, output[(i+n4*j)*num_loop+jj*nc+jjj] += -conj_PRECISION(tm_term[jj*nv+jjj]);)
	VECTOR_LOOP(jj, nrt, jjj, output[(n4*(i+n)+(j+n))*num_loop+jj*nc+jjj] -= tm_term[jj*nv+jjj];)
	VECTOR_LOOP(jj, nrt, jjj, output[((i+n)+n4*(j+n))*num_loop+jj*nc+jjj] -= -conj_PRECISION(tm_term[jj*nv+jjj]);)
        tm_term+=num_loop;
      }
      VECTOR_LOOP(jj, nrt, jjj, output[((n4+1)*j)*num_loop+jj*nc+jjj] += tm_term[jj*nv+jjj];)
      VECTOR_LOOP(jj, nrt, jjj, output[((n4+1)*(j+n))*num_loop+jj*nc+jjj] -= tm_term[jj*nv+jjj];)
      tm_term+=num_loop; // diagonal entry
    }
  }

  config_PRECISION epsbar_term = op->epsbar_term + n*(n+1)*index;
  // G
  for ( j=n; j<n2; j++ ) {
    for ( i=0; i<(j-n); i++ ) {
      VECTOR_LOOP(jj, nrt, jjj, output[(n4*i+j)*num_loop+jj*nc+jjj] = (*epsbar_term);)
      VECTOR_LOOP(jj, nrt, jjj, output[((i+n)+n4*(j-n))*num_loop+jj*nc+jjj] = -conj_PRECISION(*epsbar_term);)
      VECTOR_LOOP(jj, nrt, jjj, output[(n4*(i+n)+(j-n))*num_loop+jj*nc+jjj] = (*epsbar_term);)
      VECTOR_LOOP(jj, nrt, jjj, output[(i+n4*j)*num_loop+jj*nc+jjj] = -conj_PRECISION(*epsbar_term);)
      epsbar_term++;      
    }
    VECTOR_LOOP(jj, nrt, jjj, output[((n4+1)*(j-n)+n)*num_loop+jj*nc+jjj] = (*epsbar_term);)
    VECTOR_LOOP(jj, nrt, jjj, output[((n4+1)*j-n)*num_loop+jj*nc+jjj] = (*epsbar_term);)
    epsbar_term++; // diagonal entry
  }
  // H
  for ( j=n3; j<n4; j++ ) {
    for ( i=n2; i<(j-n); i++ ) {
      VECTOR_LOOP(jj, nrt, jjj, output[(n4*i+j)*num_loop+jj*nc+jjj] = (*epsbar_term);)
      VECTOR_LOOP(jj, nrt, jjj, output[((i+n)+n4*(j-n))*num_loop+jj*nc+jjj] = -conj_PRECISION(*epsbar_term);)
      VECTOR_LOOP(jj, nrt, jjj, output[(n4*(i+n)+(j-n))*num_loop+jj*nc+jjj] = (*epsbar_term);)
      VECTOR_LOOP(jj, nrt, jjj, output[(i+n4*j)*num_loop+jj*nc+jjj] = -conj_PRECISION(*epsbar_term);)
      epsbar_term++;      
    }
    VECTOR_LOOP(jj, nrt, jjj, output[((n4+1)*(j-n)+n)*num_loop+jj*nc+jjj] = (*epsbar_term);)
    VECTOR_LOOP(jj, nrt, jjj, output[((n4+1)*j-n)*num_loop+jj*nc+jjj] = (*epsbar_term);)
    epsbar_term++; // diagonal entry
  }
    
  // compute LU decomposition
  // output = triu(L,1) + tril(U,0), i.e., lower strictly-triangular part of output is set to triu(L,1), and upper triangular part to tril(U,0)
  // i.e., output contains L and U without the diagonal of L which is equal to 1
  // order: row major
  for ( k=0; k<n4; k++ ) {//for each column
    for ( i=k+1; i<n4; i++ ) { // for L
      VECTOR_LOOP(jj, nrt, jjj,
		  output[(n4*i+k)*num_loop+jj*nc+jjj] = output[(n4*i+k)*num_loop+jj*nc+jjj]/output[((n4+1)*k)*num_loop+jj*nc+jjj];) // output(i,k) = output(i,k)/output(k,k)
      // for U
      for ( j=k+1; j<n4; j++ )// output(i,j) = output(i,j)-output(i,k)*output(k,j)
	VECTOR_LOOP(jj, nrt, jjj,
		  output[(n4*i+j)*num_loop+jj*nc+jjj] = output[(n4*i+j)*num_loop+jj*nc+jjj]-output[(n4*i+k)*num_loop+jj*nc+jjj]*output[(n4*k+j)*num_loop+jj*nc+jjj];) 
    }
  }
#else
  // set the matrix up
  // 0
  for ( j=0; j<n; j++ ) {
    for ( i=0; i<n; i++ ) {
      output[n4*(i+0 )+(j+n3)] = _COMPLEX_PRECISION_ZERO;
      output[n4*(i+n )+(j+n2)] = _COMPLEX_PRECISION_ZERO;
      output[n4*(i+n2)+(j+n )] = _COMPLEX_PRECISION_ZERO;
      output[n4*(i+n3)+(j+0 )] = _COMPLEX_PRECISION_ZERO;
    }
  }

  config_PRECISION clover = op->clover + n*(n2+1)*index;
  // A
  for ( j=0; j<n; j++ ) {
    for ( i=0; i<j; i++ ) {
      output[n4*i+j] = *clover;
      output[i+n4*j] = conj_PRECISION(*clover);
      output[n4*(i+n)+(j+n)] = *clover;
      output[(i+n)+n4*(j+n)] = conj_PRECISION(*clover);
      clover++;      
    }
    output[(n4+1)*j] = *clover;
    output[(n4+1)*(j+n)] = *clover;
    clover++; // diagonal entry
  }
  // D
  for ( j=n2; j<n3; j++ ) {
    for ( i=n2; i<j; i++ ) {
      output[n4*i+j] = *clover;
      output[i+n4*j] = conj_PRECISION(*clover);
      output[n4*(i+n)+(j+n)] = *clover;
      output[(i+n)+n4*(j+n)] = conj_PRECISION(*clover);
      clover++;      
    }
    output[(n4+1)*j] = *clover;
    output[(n4+1)*(j+n)] = *clover;
    clover++; // diagonal entry
  }
  // B and C
  for ( j=n2; j<n3; j++ ) {
    for ( i=0; i<n; i++ ) {
      output[n4*i+j] = *clover;
      output[i+n4*j] = -conj_PRECISION(*clover);
      output[n4*(i+n)+(j+n)] = *clover;
      output[(i+n)+n4*(j+n)] = -conj_PRECISION(*clover);
      clover++;      
    }
  }

#ifdef HAVE_TM
  double even = (l->level == 0 && 0)?op->even_shift_avg: op->mu_even_shift[0];
  complex_PRECISION mu_even = (complex_PRECISION) I*(op->mu+even);
  complex_PRECISION odd_factor = (complex_PRECISION) I*(op->mu_odd_shift - even);
  config_PRECISION odd_proj = op->odd_proj+n*(n+1)*index;

  if (op->mu + op->mu_odd_shift != 0.0 || op->is_even_shifted_mu_nonzero ) {
    // E
    for ( j=0; j<n; j++ ) {
      for ( i=0; i<j; i++ ) {
        output[n4*i+j] += -odd_factor * (*odd_proj);
        output[i+n4*j] += -conj_PRECISION(-odd_factor * (*odd_proj));
        output[n4*(i+n)+(j+n)] -= -odd_factor * (*odd_proj);
        output[(i+n)+n4*(j+n)] -= -conj_PRECISION(-odd_factor * (*odd_proj));
        odd_proj++;
      }
      output[(n4+1)*j] += -1.* ( mu_even + odd_factor * (*odd_proj));
      output[(n4+1)*(j+n)] -= -1.* ( mu_even + odd_factor * (*odd_proj));
      odd_proj++; // diagonal entry
    }
    // F
    for ( j=n2; j<n3; j++ ) {
      for ( i=n2; i<j; i++ ) {
        output[n4*i+j] += odd_factor * (*odd_proj);
        output[i+n4*j] += -conj_PRECISION(odd_factor * (*odd_proj));
        output[n4*(i+n)+(j+n)] -= odd_factor * (*odd_proj);
        output[(i+n)+n4*(j+n)] -= -conj_PRECISION(odd_factor * (*odd_proj));
        odd_proj++;      
      }
      output[(n4+1)*j] += mu_even + odd_factor * (*odd_proj);
      output[(n4+1)*(j+n)] -= mu_even + odd_factor * (*odd_proj);
      odd_proj++; // diagonal entry
    }
  }
#endif

  config_PRECISION epsbar_term = op->epsbar_term + n*(n+1)*index;
  // G
  for ( j=n; j<n2; j++ ) {
    for ( i=0; i<(j-n); i++ ) {
      output[n4*i+j] = (*epsbar_term);
      output[(i+n)+n4*(j-n)] = -conj_PRECISION(*epsbar_term);
      output[n4*(i+n)+(j-n)] = (*epsbar_term);
      output[i+n4*j] = -conj_PRECISION(*epsbar_term);
      epsbar_term++;      
    }
    output[(n4+1)*(j-n)+n] = (*epsbar_term);
    output[(n4+1)*j-n] = (*epsbar_term);
    epsbar_term++; // diagonal entry
  }
  // H
  for ( j=n3; j<n4; j++ ) {
    for ( i=n2; i<(j-n); i++ ) {
      output[n4*i+j] = (*epsbar_term);
      output[(i+n)+n4*(j-n)] = -conj_PRECISION(*epsbar_term);
      output[n4*(i+n)+(j-n)] = (*epsbar_term);
      output[i+n4*j] = -conj_PRECISION(*epsbar_term);
      epsbar_term++;      
    }
    output[(n4+1)*(j-n)+n] = (*epsbar_term);
    output[(n4+1)*j-n] = (*epsbar_term);
    epsbar_term++; // diagonal entry
  }
    
  // compute LU decomposition
  // output = triu(L,1) + tril(U,0)
  // i.e., output contains L and U without the diagonal of L which is equal to 1
  // order: row major
  for ( k=0; k<n4; k++ ) {
    for ( i=k+1; i<n4; i++ ) {
      output[n4*i+k] = output[n4*i+k]/output[(n4+1)*k]; // output(i,k) = output(i,k)/output(k,k)
      for ( j=k+1; j<n4; j++ )
        output[n4*i+j] = output[n4*i+j]-output[n4*i+k]*output[n4*k+j]; // output(i,j) = output(i,j)-output(i,k)*output(k,j)
    }
  }
  complex_PRECISION tmp;
  for ( i=0; i<n4; i++ ) {
    for ( j=0; j<n; j++ ) {
      tmp = output[n4*i+n+j];
      output[n4*i+n+j] = output[n4*i+n2+j];
      output[n4*i+n2+j] = tmp;
    }
  }
#endif
#endif
}
#endif

// out_o/e += D_oe/eo*in_e/o if amount==_ODD_SITES/_EVEN_SITES; works on both parts if amount==_FULL_SYSTEM 
// in and out needs to be at least of size _ORDINARY
void coarse_hopping_term_PRECISION( vector_PRECISION *out, vector_PRECISION *in, operator_PRECISION_struct *op,
                                    const int amount, level_struct *l, struct Thread *threading ) {

  START_NO_HYPERTHREADS(threading)

  int mu, i, index, j, jj;
  int num_site_var  = l->num_lattice_site_var;
  int num_4link_var = 4*4*l->num_parent_eig_vect*l->num_parent_eig_vect;
  int num_link_var  = 4*l->num_parent_eig_vect*l->num_parent_eig_vect;
  int start             = 0;
  int num_lattice_sites = l->num_inner_lattice_sites;
  int plus_dir_param =_FULL_SYSTEM, minus_dir_param =_FULL_SYSTEM;
  int nvec = in->num_vect_now, nvec_in = in->num_vect, nvec_out=out->num_vect;
  vector_PRECISION in_pt, out_pt;
  config_PRECISION D_pt;

#ifdef DEBUG
  if ( nvec_out < nvec_in )
    error0("coarse_hopping_term_PRECISION: assunmptions are not met\n");
#endif
  
  g.num_vect_pass1 = nvec;
  vector_PRECISION_duplicate( &in_pt, in, 0, l );
  vector_PRECISION_duplicate( &out_pt, out, 0, l );

  int core_start;
  int core_end;
  
  // assumptions (1) self coupling has already been performed
  //          OR (2) "out" is initialized with zeros
  set_boundary_PRECISION( out, 0, l, threading );
  
  if ( amount == _EVEN_SITES ) {
    minus_dir_param = _ODD_SITES;
    plus_dir_param = _EVEN_SITES;
  } else if ( amount == _ODD_SITES ) {
    minus_dir_param = _EVEN_SITES;
    plus_dir_param = _ODD_SITES;
  }
  
  START_MASTER(threading)
  if ( op->c.comm ) {
    g.num_vect_pass2 = nvec_in;
    for ( mu=0; mu<4; mu++ ) {
      // communicate in -mu direction
      ghost_sendrecv_PRECISION( in->vector_buffer, mu, -1, &(op->c), minus_dir_param, l );
    }
  }
  END_MASTER(threading)
  SYNC_CORES(threading)
  
  if ( amount == _EVEN_SITES ) { // partition the odd sites
    start = op->num_even_sites, num_lattice_sites = op->num_odd_sites;
  } else if ( amount == _ODD_SITES ) { // partition the even sites
    start = 0; num_lattice_sites = op->num_even_sites;
  }
  compute_core_start_end_custom(start, start+num_lattice_sites, &core_start, &core_end, l, threading, 1);
  
  // compute U_mu^dagger coupling
  for ( i=core_start; i<core_end; i++ ) {
    index = 5*i;
    in_pt.vector_buffer = in->vector_buffer + num_site_var*op->neighbor_table[index]*nvec_in;
    D_pt = op->D + num_4link_var*op->neighbor_table[index] + 0*num_link_var;
    index++;
    out_pt.vector_buffer = out->vector_buffer + num_site_var*op->neighbor_table[index+T]*nvec_out;
    coarse_daggered_hopp_PRECISION( &out_pt, &in_pt, D_pt, l );
  }
  SYNC_CORES(threading)
  for ( i=core_start; i<core_end; i++ ) {
    index = 5*i;
    in_pt.vector_buffer = in->vector_buffer + num_site_var*op->neighbor_table[index]*nvec_in;
    D_pt = op->D + num_4link_var*op->neighbor_table[index] + 1*num_link_var;
    index++;
    out_pt.vector_buffer = out->vector_buffer + num_site_var*op->neighbor_table[index+Z]*nvec_out;
    coarse_daggered_hopp_PRECISION( &out_pt, &in_pt, D_pt, l );
  }
  SYNC_CORES(threading)
  for ( i=core_start; i<core_end; i++ ) {
    index = 5*i;
    in_pt.vector_buffer = in->vector_buffer + num_site_var*op->neighbor_table[index]*nvec_in;
    D_pt = op->D + num_4link_var*op->neighbor_table[index] + 2*num_link_var;
    index++;
    out_pt.vector_buffer = out->vector_buffer + num_site_var*op->neighbor_table[index+Y]*nvec_out;
    coarse_daggered_hopp_PRECISION( &out_pt, &in_pt, D_pt, l );
  }
  SYNC_CORES(threading)
  for ( i=core_start; i<core_end; i++ ) {
    index = 5*i;
    in_pt.vector_buffer = in->vector_buffer + num_site_var*op->neighbor_table[index]*nvec_in;
    D_pt = op->D + num_4link_var*op->neighbor_table[index] + 3*num_link_var;
    index++;
    out_pt.vector_buffer = out->vector_buffer + num_site_var*op->neighbor_table[index+X]*nvec_out;
    coarse_daggered_hopp_PRECISION( &out_pt, &in_pt, D_pt, l );
  }
  
  START_LOCKED_MASTER(threading)
  if ( op->c.comm ) {
    g.num_vect_pass2 = nvec_out;
    for ( mu=0; mu<4; mu++ ) {
      // communicate in +mu direction
      ghost_sendrecv_PRECISION( out->vector_buffer, mu, +1, &(op->c), plus_dir_param, l );    
    }
    g.num_vect_pass2 = nvec_in;
    for ( mu=0; mu<4; mu++ ) {
      // wait for -mu direction
      ghost_wait_PRECISION( in->vector_buffer, mu, -1, &(op->c), minus_dir_param, l );    
    }
  }
  END_LOCKED_MASTER(threading)
  
  if ( amount == _EVEN_SITES ) {
    start = 0; num_lattice_sites = op->num_even_sites;
  } else if ( amount == _ODD_SITES ) {
    start = op->num_even_sites, num_lattice_sites = op->num_odd_sites;
  }
  compute_core_start_end_custom(start, start+num_lattice_sites, &core_start, &core_end, l, threading, 1);
  
  // compute U_mu couplings
  for ( i=core_start; i<core_end; i++ ) {
    index = 5*i;
    out_pt.vector_buffer = out->vector_buffer + num_site_var*op->neighbor_table[index]*nvec_out;
    D_pt = op->D + num_4link_var*op->neighbor_table[index];
    index++;
    in_pt.vector_buffer = in->vector_buffer + num_site_var*op->neighbor_table[index+T]*nvec_in;
    coarse_hopp_PRECISION( &out_pt, &in_pt, D_pt, l );
    
    D_pt += num_link_var;
    in_pt.vector_buffer = in->vector_buffer + num_site_var*op->neighbor_table[index+Z]*nvec_in;
    coarse_hopp_PRECISION( &out_pt, &in_pt, D_pt, l );
    
    D_pt += num_link_var;
    in_pt.vector_buffer = in->vector_buffer + num_site_var*op->neighbor_table[index+Y]*nvec_in;
    coarse_hopp_PRECISION( &out_pt, &in_pt, D_pt, l );
    
    D_pt += num_link_var;
    in_pt.vector_buffer = in->vector_buffer + num_site_var*op->neighbor_table[index+X]*nvec_in;
    coarse_hopp_PRECISION( &out_pt, &in_pt, D_pt, l );
  }
  
  START_LOCKED_MASTER(threading)
  if ( op->c.comm ) {
    g.num_vect_pass2 = nvec_out;
    for ( mu=0; mu<4; mu++ ) {
      // wait for +mu direction
      ghost_wait_PRECISION( out->vector_buffer, mu, +1, &(op->c), plus_dir_param, l );
    }
  }
  END_LOCKED_MASTER(threading)

  END_NO_HYPERTHREADS(threading)
}

// out_o/e -= D_oe/eo*in_e/o if amount==_ODD_SITES/_EVEN_SITES; works on both parts if amount==_FULL_SYSTEM
void coarse_n_hopping_term_PRECISION( vector_PRECISION *out, vector_PRECISION *in, operator_PRECISION_struct *op,
                                      const int amount, level_struct *l, struct Thread *threading ) {

  START_NO_HYPERTHREADS(threading)

  int mu, i, index, j, jj;
  int num_site_var  = l->num_lattice_site_var;
  int num_4link_var = 4*4*l->num_parent_eig_vect*l->num_parent_eig_vect;
  int num_link_var  = 4*l->num_parent_eig_vect*l->num_parent_eig_vect;
  int start             = 0;
  int num_lattice_sites = l->num_inner_lattice_sites;
  int plus_dir_param = _FULL_SYSTEM, minus_dir_param = _FULL_SYSTEM;
  int nvec = in->num_vect_now, nvec_in = in->num_vect, nvec_out = out->num_vect;
  vector_PRECISION in_pt, out_pt;
  config_PRECISION D_pt;

#ifdef DEBUG
  if ( nvec_out < nvec )
    error0("coarse_n_hopping_term_PRECISION: assumptions are not met\n");
#endif
  
  vector_PRECISION_duplicate( &in_pt, in, 0, l );
  vector_PRECISION_duplicate( &out_pt, out, 0, l );
  g.num_vect_pass1 = nvec;

  int core_start;
  int core_end;
  
  // assumptions (1) self coupling has already been performed
  //          OR (2) "out" is initialized with zeros
  set_boundary_PRECISION( out, 0, l, threading );
  
  if ( amount == _EVEN_SITES ) {
    minus_dir_param = _ODD_SITES;
    plus_dir_param = _EVEN_SITES;
  } else if ( amount == _ODD_SITES ) {
    minus_dir_param = _EVEN_SITES;
    plus_dir_param = _ODD_SITES;
  }
  
  START_MASTER(threading)
  if ( op->c.comm ) {
    for ( mu=0; mu<4; mu++ ) {
      // communicate in -mu direction
      g.num_vect_pass2 = nvec_in;
      ghost_sendrecv_PRECISION( in->vector_buffer, mu, -1, &(op->c), minus_dir_param, l );
    }
  }
  END_MASTER(threading)
  SYNC_CORES(threading)
  
  if ( amount == _EVEN_SITES ) {
    start = op->num_even_sites, num_lattice_sites = op->num_odd_sites;
  } else if ( amount == _ODD_SITES ) {
    start = 0; num_lattice_sites = op->num_even_sites;
  }
  compute_core_start_end_custom(start, start+num_lattice_sites, &core_start, &core_end, l, threading, 1);
  
  // compute U_mu^dagger coupling
  for ( i=core_start; i<core_end; i++ ) {
    index = 5*i;
    in_pt.vector_buffer = in->vector_buffer + num_site_var*op->neighbor_table[index]*nvec_in;
    D_pt = op->D + num_4link_var*op->neighbor_table[index] + 0*num_link_var;
    index++;
    out_pt.vector_buffer = out->vector_buffer + num_site_var*op->neighbor_table[index+T]*nvec_out;
    coarse_n_daggered_hopp_PRECISION( &out_pt, &in_pt, D_pt, l );
  }
  SYNC_CORES(threading)
  for ( i=core_start; i<core_end; i++ ) {
    index = 5*i;
    in_pt.vector_buffer = in->vector_buffer + num_site_var*op->neighbor_table[index]*nvec_in;
    D_pt = op->D + num_4link_var*op->neighbor_table[index] + 1*num_link_var;
    index++;
    out_pt.vector_buffer = out->vector_buffer + num_site_var*op->neighbor_table[index+Z]*nvec_out;
    coarse_n_daggered_hopp_PRECISION( &out_pt, &in_pt, D_pt, l );
  }
  SYNC_CORES(threading)
  for ( i=core_start; i<core_end; i++ ) {
    index = 5*i;
    in_pt.vector_buffer = in->vector_buffer + num_site_var*op->neighbor_table[index]*nvec_in;
    D_pt = op->D + num_4link_var*op->neighbor_table[index] + 2*num_link_var;
    index++;
    out_pt.vector_buffer = out->vector_buffer + num_site_var*op->neighbor_table[index+Y]*nvec_out;
    coarse_n_daggered_hopp_PRECISION( &out_pt, &in_pt, D_pt, l );
  }
  SYNC_CORES(threading)
  for ( i=core_start; i<core_end; i++ ) {
    index = 5*i;
    in_pt.vector_buffer = in->vector_buffer + num_site_var*op->neighbor_table[index]*nvec_in;
    D_pt = op->D + num_4link_var*op->neighbor_table[index] + 3*num_link_var;
    index++;
    out_pt.vector_buffer = out->vector_buffer + num_site_var*op->neighbor_table[index+X]*nvec_out;
    coarse_n_daggered_hopp_PRECISION( &out_pt, &in_pt, D_pt, l );
  }
  
  START_LOCKED_MASTER(threading)
  if ( op->c.comm ) {
    g.num_vect_pass2 = nvec_out;
    for ( mu=0; mu<4; mu++ ) {
      // communicate in +mu direction
      ghost_sendrecv_PRECISION( out->vector_buffer, mu, +1, &(op->c), plus_dir_param, l );    
    }
    g.num_vect_pass2 = nvec_in;
    for ( mu=0; mu<4; mu++ ) {
      // wait for -mu direction
      ghost_wait_PRECISION( in->vector_buffer, mu, -1, &(op->c), minus_dir_param, l );    
    }
  }
  END_LOCKED_MASTER(threading)
  
  if ( amount == _EVEN_SITES ) {
    start = 0; num_lattice_sites = op->num_even_sites;
  } else if ( amount == _ODD_SITES ) {
    start = op->num_even_sites, num_lattice_sites = op->num_odd_sites;
  }
  compute_core_start_end_custom(start, start+num_lattice_sites, &core_start, &core_end, l, threading, 1);
  
  // compute U_mu couplings
  for ( i=core_start; i<core_end; i++ ) {
    index = 5*i;
    out_pt.vector_buffer = out->vector_buffer + num_site_var*op->neighbor_table[index]*nvec_out;
    D_pt = op->D + num_4link_var*op->neighbor_table[index];
    index++;
    in_pt.vector_buffer = in->vector_buffer + num_site_var*op->neighbor_table[index+T]*nvec_in;
    coarse_n_hopp_PRECISION( &out_pt, &in_pt, D_pt, l );
    
    D_pt += num_link_var;
    in_pt.vector_buffer = in->vector_buffer + num_site_var*op->neighbor_table[index+Z]*nvec_in;
    coarse_n_hopp_PRECISION( &out_pt, &in_pt, D_pt, l );
    
    D_pt += num_link_var;
    in_pt.vector_buffer = in->vector_buffer + num_site_var*op->neighbor_table[index+Y]*nvec_in;
    coarse_n_hopp_PRECISION( &out_pt, &in_pt, D_pt, l );
    
    D_pt += num_link_var;
    in_pt.vector_buffer = in->vector_buffer + num_site_var*op->neighbor_table[index+X]*nvec_in;
    coarse_n_hopp_PRECISION( &out_pt, &in_pt, D_pt, l );
  }
  
  START_LOCKED_MASTER(threading)
  if ( op->c.comm ) {
    g.num_vect_pass2 = nvec_out;
    for ( mu=0; mu<4; mu++ ) {
      // wait for +mu direction
      ghost_wait_PRECISION( out->vector_buffer, mu, +1, &(op->c), plus_dir_param, l );
    }
  }
  END_LOCKED_MASTER(threading)

  END_NO_HYPERTHREADS(threading)
}

// D_sc = D_ee - D_eo D_oo ^{-1} D_oe
// used at the bottom; also used in preconditioner if g.method==4,5
void coarse_apply_schur_complement_PRECISION( vector_PRECISION *out, vector_PRECISION *in, operator_PRECISION_struct *op, level_struct *l, struct Thread *threading ) {

  int start;
  int end;
  // compute start and end indices for core
  // this puts zero for all other hyperthreads, so we can call functions below with all hyperthreads
  compute_core_start_end_custom(op->num_even_sites*l->num_lattice_site_var, l->inner_vector_size, &start, &end, l, threading, l->num_lattice_site_var );

  vector_PRECISION *tmp = op->buffer;
  for ( int i=0; i<2; i++ ) tmp[i].num_vect_now = in->num_vect_now;

  SYNC_CORES(threading)
  PROF_PRECISION_START( _SC, threading );
  coarse_diag_ee_PRECISION( out, in, op, l, threading );
  PROF_PRECISION_STOP( _SC, 0, threading );
  SYNC_CORES(threading)
  vector_PRECISION_define( &tmp[0], 0, start, end, l );
  SYNC_CORES(threading)
  PROF_PRECISION_START( _NC, threading );
  coarse_hopping_term_PRECISION( &tmp[0], in, op, _ODD_SITES, l, threading );
  PROF_PRECISION_STOP( _NC, 0, threading );
  PROF_PRECISION_START( _SC, threading );
  coarse_diag_oo_inv_PRECISION( &tmp[1], &tmp[0], op, l, threading );
  PROF_PRECISION_STOP( _SC, 1, threading );
  PROF_PRECISION_START( _NC, threading );
  coarse_n_hopping_term_PRECISION( out, &tmp[1], op, _EVEN_SITES, l, threading );
  PROF_PRECISION_STOP( _NC, 1, threading );
}

// simply apply even part of self-couling terms w/t making use of LU decomposition
// y_e = D_ee*x_e
void coarse_diag_ee_PRECISION( vector_PRECISION *y, vector_PRECISION *x, operator_PRECISION_struct *op, level_struct *l, struct Thread *threading ) {

  int start, end;
  compute_core_start_end_custom( 0, op->num_even_sites, &start, &end, l, threading, 1 );
  // even sites
  coarse_self_couplings_PRECISION( y, x, op, start, end, l );
}

// used only in test routines
void coarse_diag_oo_PRECISION( vector_PRECISION *y, vector_PRECISION *x, operator_PRECISION_struct *op, level_struct *l, struct Thread *threading ) {
  
  int start, end;
  
#ifdef DEBUG
  if ( y->num_vect < x->num_vect_now )
    error0("coarse_diag_oo_PRECISION: assumptions are not met\n");
#endif
  
#ifdef HAVE_TM1p1
  config_PRECISION sc = (g.n_flavours==2) ? op->clover_doublet_oo_inv:op->clover_oo_inv;
#else
  config_PRECISION sc = op->clover_oo_inv;
#endif
  
  compute_core_start_end_custom( 0, op->num_odd_sites, &start, &end, l, threading, 1 );

#if 1
  int s, fac = 1, nvec = x->num_vect_now, nvec_x = x->num_vect, nvec_y = y->num_vect;
  int num_site_var=l->num_lattice_site_var;
  int oo_inv_size = SQUARE(num_site_var);
  buffer_PRECISION x_pt = x->vector_buffer+(op->num_even_sites+start)*num_site_var*nvec_x, y_pt = y->vector_buffer+(op->num_even_sites+start)*num_site_var*nvec_y;
#ifdef HAVE_MULT_TM
  fac *= num_loop;
#endif
#ifdef HAVE_TM1p1
  if (g.n_flavours==2)
    oo_inv_size *= 4;
#endif
  sc += oo_inv_size*start*fac;
  for ( s=start; s<end; s++ ) {
    coarse_LU_multiply_PRECISION( y_pt, x_pt, sc, nvec, nvec_y, nvec_x, l );
    y_pt += num_site_var*nvec_y;
    x_pt += num_site_var*nvec_x;
    sc += oo_inv_size*fac;
  }
#else
  coarse_LU_multiply_PRECISION( y, x, op, start, end, l );
#endif

}

static void coarse_diag_oo_inv_PRECISION( vector_PRECISION *y, vector_PRECISION *x, operator_PRECISION_struct *op, 
					      level_struct *l, struct Thread *threading ) {

  int start, end;

#ifdef DEBUG
  if ( y->num_vect < x->num_vect_now )
    error0("coarse_diag_oo_inv_PRECISION: assumptions are not met\n");
#endif
  
  compute_core_start_end_custom( 0, op->num_odd_sites, &start, &end, l, threading, 1 );

#ifdef HAVE_TM1p1
  config_PRECISION sc = (g.n_flavours==2) ? op->clover_doublet_oo_inv:op->clover_oo_inv;
#else
  config_PRECISION sc = op->clover_oo_inv;
#endif
  
#if 1
  int s, fac = 1, nvec = x->num_vect_now, nvec_x = x->num_vect, nvec_y = y->num_vect;
  int num_site_var = l->num_lattice_site_var;
  int oo_inv_size  = SQUARE(num_site_var);
  buffer_PRECISION x_pt = x->vector_buffer+(op->num_even_sites+start)*num_site_var*nvec_x, y_pt = y->vector_buffer+(op->num_even_sites+start)*num_site_var*nvec_y;
  //  config_PRECISION sc = op->clover_oo_inv;
#ifdef HAVE_MULT_TM
  fac *= num_loop;
#endif
#ifdef HAVE_TM1p1
  if (g.n_flavours==2)
    oo_inv_size *= 4;
#endif
  
  sc += oo_inv_size*start*fac;
  for ( s=start; s<end; s++ ) {
    coarse_perform_fwd_bwd_subs_PRECISION( y_pt, x_pt, sc, nvec, nvec_y, nvec_x, l );
    y_pt += num_site_var*nvec_y;
    x_pt += num_site_var*nvec_x;
    sc += oo_inv_size*fac;
  }
#else
  coarse_perform_fwd_bwd_subs_PRECISION( y, x, op, start, end, l );
#endif
}

int coarse_solve_odd_even_PRECISION( gmres_PRECISION_struct *p, operator_PRECISION_struct *op, level_struct *l, struct Thread *threading ) {
  /****************************************************************************************
   * Descripton: Solve D*x = b in the even-odd preconditioning using Schur complement D_sc
   *   D_sc = D_ee - D_eo D_oo ^{-1} D_oe
   *   x_e = D_sc^{-1} (b_e - D_eo*D_oo^{-1}*b_o)
   *   x_o = D_oo^{-1} (b_o - D_oe*x_e)  
   ****************************************************************************************/

  //p->b.num_vect_now = g.num_vect_now; p->x.num_vect_now = g.num_vect_now;//!!!!!!
  //p->b.num_vect_now = num_loop; p->x.num_vect_now = num_loop;//!!!!!!

  // solve for x_e
  SYNC_CORES(threading)
  PROF_PRECISION_START( _SC, threading ); 
  coarse_diag_oo_inv_PRECISION( &p->x, &p->b, op, l, threading ); // x_o = D_oo^{-1}*b_o
  PROF_PRECISION_STOP( _SC, 0, threading );
  PROF_PRECISION_START( _NC, threading );
  coarse_n_hopping_term_PRECISION( &p->b, &p->x, op, _EVEN_SITES, l, threading ); // b_e = b_e - D_eo*D_oo^{-1}*b_o
  PROF_PRECISION_STOP( _NC, 0, threading );

  int iter = solver_PRECISION( p, l, threading ); // x_e = D_sc^{-1} (b_e - D_eo*D_oo^{-1}*b_o)
  
  // construct x_o from x_e
  PROF_PRECISION_START( _NC, threading );
  coarse_n_hopping_term_PRECISION( &p->b, &p->x, op, _ODD_SITES, l, threading ); // b_o = b_o - D_oe*x_e
  PROF_PRECISION_STOP( _NC, 1, threading );
  PROF_PRECISION_START( _SC, threading );
  coarse_diag_oo_inv_PRECISION( &p->x, &p->b, op, l, threading ); // x_o = D_oo^{-1} (b_o - D_oe*x_e)
  PROF_PRECISION_STOP( _SC, 1, threading );
  SYNC_CORES(threading)

  START_LOCKED_MASTER(threading)
  g.iter_counts[l->depth] += iter;
  END_LOCKED_MASTER(threading)
    
  return iter;
}

/************************  TEST ROUTINES  ********************************************/

void coarse_odd_even_PRECISION_test( vector_PRECISION *out, vector_PRECISION *in, level_struct *l, struct Thread *threading ) {
  
  if ( g.odd_even ) {
#ifdef DEBUG
    if ( out->num_vect != in->num_vect || out->num_vect_now != in->num_vect_now )
      error0("coarse_odd_even_PRECISION_test: assumptions are not met\n");
#endif
    vector_PRECISION buf[2];

    for(int i=0; i<2; i++){
      vector_PRECISION_init( &buf[i] );
      vector_PRECISION_alloc( &buf[i], _ORDINARY, in->num_vect, l, threading );
      buf[i].num_vect_now = in->num_vect_now;
    }
    
    START_LOCKED_MASTER(threading)
    // transformation part
      // buf[0] = in
    vector_PRECISION_copy( &buf[0], in, 0, l->inner_vector_size, l );
    // set out_odd = 0
    vector_PRECISION_define( out, 0, l->oe_op_PRECISION.num_even_sites*l->num_lattice_site_var, l->inner_vector_size, l );
    END_LOCKED_MASTER(threading)

    // out_o += D_oe*in_e => out_o = D_oe*in_e 
    coarse_hopping_term_PRECISION( out, &buf[0], &(l->oe_op_PRECISION), _ODD_SITES, l, threading );
    // buff[1]_o = D_oo^{-1}out_o = D_oo^{-1}D_oe*in_e
    coarse_diag_oo_inv_PRECISION( &buf[1], out, &(l->oe_op_PRECISION), l, threading );

    START_LOCKED_MASTER(threading)
    // buf[0]_o = in_o + buf[1]_o = in_o+D_oo^{-1}D_oe*in_e 
    vector_PRECISION_plus( &buf[0], &buf[0], &buf[1], l->oe_op_PRECISION.num_even_sites*l->num_lattice_site_var, l->inner_vector_size, l ); 
    END_LOCKED_MASTER(threading)
    
    // block diagonal part
    coarse_apply_schur_complement_PRECISION( out, &buf[0], &(l->oe_op_PRECISION), l, threading ); // out_e = D_sc*in_e
    
    coarse_diag_oo_PRECISION( out, &buf[0], &(l->oe_op_PRECISION), l, threading ); // out_o = D_oo(in_o+D_oo^{-1}D_oe*in_e) = D_oo*in_o+D_oe*in_e
    
    // back transformation part
    coarse_diag_oo_inv_PRECISION( &buf[1], out, &(l->oe_op_PRECISION), l, threading ); // buf[1]_o = D_oo^{-1}(D_oo*in_o+D_oe*in_e) = D_oo^{-1}D_oe*in_e+in_o
    //out_e += D_eo*buff[1]_o => out_e = D_sc*in_e + D_eoD_oo^{-1}D_oe*in_e+D_eo*in_o = D_ee*in_e+D_eo*in_o
    coarse_hopping_term_PRECISION( out, &buf[1], &(l->oe_op_PRECISION), _EVEN_SITES, l, threading ); 
 
    for(int i=0; i<2; i++)
      vector_PRECISION_free( &buf[i], l, threading );                      
  }
}
