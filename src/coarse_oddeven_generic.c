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
static void coarse_diag_oo_inv_PRECISION_new( vector_PRECISION *y, vector_PRECISION *x, operator_PRECISION_struct *op, level_struct *l, struct Thread *threading );

void coarse_oddeven_alloc_PRECISION( level_struct *l ) {

  int nv = l->num_parent_eig_vect, oe_offset=0, mu, **bt = NULL,
    *eot = NULL, *nt = NULL, *tt = NULL, t, z, y, x, le[4], N[4];
  int nvec = num_loop;//(g.num_rhs_vect < l->num_eig_vect)? l->num_eig_vect:g.num_rhs_vect;
  operator_PRECISION_struct *op = &(l->oe_op_PRECISION);

  operator_PRECISION_alloc( op, _ODDEVEN, l );

  // buffers
  MALLOC( op->buffer, vector_PRECISION, 2 );
  for (int k=0; k<2; k++ ){
    vector_PRECISION_init( &(op->buffer[k]) );
#ifdef HAVE_TM1p1
    vector_PRECISION_alloc( &(op->buffer[k]), _ORDINARY, 2*nvec, l, no_threading );
#else
    vector_PRECISION_alloc( &(op->buffer[k]), _ORDINARY, nvec, l, no_threading );
#endif
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
  
  MALLOC( op->clover_oo_inv, complex_PRECISION, SQUARE(2*nv)*op->num_odd_sites );
#ifdef HAVE_TM1p1
  MALLOC( op->clover_doublet_oo_inv, complex_PRECISION, SQUARE(4*nv)*op->num_odd_sites );
#endif
  // define data layout
  eot = op->index_table;
  define_eot( eot, N, l );

  // neighbor table, translation table
  nt = op->neighbor_table;
  tt = op->translation_table;
  define_nt_bt_tt( nt, op->backward_neighbor_table, NULL, tt, eot, N, l );

  // boundary table
  bt = op->c.boundary_table;
  define_eo_bt( bt, eot, op->c.num_even_boundary_sites, op->c.num_odd_boundary_sites, op->c.num_boundary_sites, N, l );

  // ghost
  ghost_sendrecv_init_PRECISION( _COARSE_GLOBAL, &(op->c), l ) ;

  // solver
  if ( l->level == 0 )
    l->p_PRECISION.v_end = op->num_even_sites*l->num_lattice_site_var;
  else
    l->sp_PRECISION.v_end = op->num_even_sites*l->num_lattice_site_var;

}

void coarse_oddeven_free_PRECISION( level_struct *l ) {
  
  int nv = l->num_parent_eig_vect;
  operator_PRECISION_struct *op = &(l->oe_op_PRECISION);

  operator_PRECISION_free( op, _ODDEVEN, l );

  FREE( op->clover_oo_inv, complex_PRECISION, SQUARE(2*nv)*op->num_odd_sites );
#ifdef HAVE_TM1p1
  FREE( op->clover_doublet_oo_inv, complex_PRECISION, SQUARE(4*nv)*op->num_odd_sites );
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

  // neighbor couplings
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
  tm_term_PRECISION_setup( in->mu, in->mu_even_shift, in->mu_odd_shift, op, l, threading );
#endif  
#ifdef HAVE_TM1p1
  epsbar_term_PRECISION_setup( in->epsbar, in->epsbar_ig5_even_shift, in->epsbar_ig5_odd_shift, op, l, threading );
#endif

  coarse_oddeven_PRECISION_set_self_couplings( l, threading );//this is the src of OpenMP issue

}

void coarse_oddeven_PRECISION_set_self_couplings( level_struct *l, struct Thread *threading ) {

  operator_PRECISION_struct *op = &(l->oe_op_PRECISION);
  int nv = l->num_parent_eig_vect, start, end;

  compute_core_start_end_custom( 0, op->num_odd_sites, &start, &end, l, threading, 1);

  SYNC_CORES(threading)//debug!!!!
  int size = SQUARE(2*nv);
  for( int i=start; i<end; i++ )
    coarse_selfcoupling_LU_decomposition_PRECISION( op->clover_oo_inv+i*size, op, op->num_even_sites+i, l );

#ifdef HAVE_TM1p1
  int size_doublet = SQUARE(4*nv);
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

  register int i, j, k, n = l->num_parent_eig_vect, n2 = 2*n;
  config_PRECISION clover = op->clover + n*(n2+1)*index;
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

#ifdef HAVE_TM
  config_PRECISION tm_term = op->tm_term + n*(n+1)*index;
  if (op->mu + op->mu_odd_shift != 0.0 || op->mu + op->mu_even_shift != 0.0 ) {
    // E
    for ( j=0; j<n; j++ ) {
      for ( i=0; i<j; i++ ) {
        output[n2*i+j] += *tm_term;
        output[i+n2*j] += -conj_PRECISION(*tm_term);
        tm_term++;      
      }
      output[(n2+1)*j] += *tm_term;
      tm_term++; // diagonal entry
    }
    // F
    for ( j=n; j<n2; j++ ) {
      for ( i=n; i<j; i++ ) {
        output[n2*i+j] += *tm_term;
        output[i+n2*j] += -conj_PRECISION(*tm_term);
        tm_term++;      
      }
      output[(n2+1)*j] += *tm_term;
      tm_term++; // diagonal entry
    }
  }
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
}

#ifdef HAVE_TM1p1
static void coarse_selfcoupling_LU_doublet_decomposition_PRECISION( config_PRECISION output, operator_PRECISION_struct *op, int index, level_struct *l ) {

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
  // epsbar_term = [ G 0      , G=-G*, H=-H* diag. excluded
  //                 0 H ]
  //
  // order: upper triangle of G, upper triangle of H
  //
  // output = [ A+E  G   B   0
  //             G  A-E  0   B
  //             C   0  D+F  H
  //             0   C   H  D-F ]  LU decomposed

  register int i, j, k, n = l->num_parent_eig_vect, n2 = 2*n, n3 = 3*n, n4 = 4*n;
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
  config_PRECISION tm_term = op->tm_term + n*(n+1)*index;
  if (op->mu + op->mu_odd_shift != 0.0 || op->mu + op->mu_even_shift != 0.0 ) {
    // E
    for ( j=0; j<n; j++ ) {
      for ( i=0; i<j; i++ ) {
        output[n4*i+j] += *tm_term;
        output[i+n4*j] += -conj_PRECISION(*tm_term);
        output[n4*(i+n)+(j+n)] -= *tm_term;
        output[(i+n)+n4*(j+n)] -= -conj_PRECISION(*tm_term);
        tm_term++;      
      }
      output[(n4+1)*j] += *tm_term;
      output[(n4+1)*(j+n)] -= *tm_term;
      tm_term++; // diagonal entry
    }
    // F
    for ( j=n2; j<n3; j++ ) {
      for ( i=n2; i<j; i++ ) {
        output[n4*i+j] += *tm_term;
        output[i+n4*j] += -conj_PRECISION(*tm_term);
        output[n4*(i+n)+(j+n)] -= *tm_term;
        output[(i+n)+n4*(j+n)] -= -conj_PRECISION(*tm_term);
        tm_term++;      
      }
      output[(n4+1)*j] += *tm_term;
      output[(n4+1)*(j+n)] -= *tm_term;
      tm_term++; // diagonal entry
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
}
#endif

// out_o/e += D_eo/oe*in_e/o
// in and out needs to be at least of size _ORDINARY
void coarse_hopping_term_PRECISION_new( vector_PRECISION *out, vector_PRECISION *in, operator_PRECISION_struct *op,
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

  if ( nvec_out < nvec_in )
    error0("coarse_hopping_term_PRECISION: assunmptions are not met\n");

  g.num_vect_pass1 = nvec;
  vector_PRECISION_duplicate( &in_pt, in, 0, l );
  vector_PRECISION_duplicate( &out_pt, out, 0, l );

  int core_start;
  int core_end;
  
  // assumptions (1) self coupling has already been performed
  //          OR (2) "out" is initialized with zeros
  set_boundary_PRECISION_new( out, 0, l, threading );
  
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
      ghost_sendrecv_PRECISION_new( in->vector_buffer, mu, -1, &(op->c), minus_dir_param, l );
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
    coarse_daggered_hopp_PRECISION_new( &out_pt, &in_pt, D_pt, l );
  }
  SYNC_CORES(threading)
  for ( i=core_start; i<core_end; i++ ) {
    index = 5*i;
    in_pt.vector_buffer = in->vector_buffer + num_site_var*op->neighbor_table[index]*nvec_in;
    D_pt = op->D + num_4link_var*op->neighbor_table[index] + 1*num_link_var;
    index++;
    out_pt.vector_buffer = out->vector_buffer + num_site_var*op->neighbor_table[index+Z]*nvec_out;
    coarse_daggered_hopp_PRECISION_new( &out_pt, &in_pt, D_pt, l );
  }
  SYNC_CORES(threading)
  for ( i=core_start; i<core_end; i++ ) {
    index = 5*i;
    in_pt.vector_buffer = in->vector_buffer + num_site_var*op->neighbor_table[index]*nvec_in;
    D_pt = op->D + num_4link_var*op->neighbor_table[index] + 2*num_link_var;
    index++;
    out_pt.vector_buffer = out->vector_buffer + num_site_var*op->neighbor_table[index+Y]*nvec_out;
    coarse_daggered_hopp_PRECISION_new( &out_pt, &in_pt, D_pt, l );
  }
  SYNC_CORES(threading)
  for ( i=core_start; i<core_end; i++ ) {
    index = 5*i;
    in_pt.vector_buffer = in->vector_buffer + num_site_var*op->neighbor_table[index]*nvec_in;
    D_pt = op->D + num_4link_var*op->neighbor_table[index] + 3*num_link_var;
    index++;
    out_pt.vector_buffer = out->vector_buffer + num_site_var*op->neighbor_table[index+X]*nvec_out;
    coarse_daggered_hopp_PRECISION_new( &out_pt, &in_pt, D_pt, l );
  }
  
  START_LOCKED_MASTER(threading)
  if ( op->c.comm ) {
    g.num_vect_pass2 = nvec_out;
    for ( mu=0; mu<4; mu++ ) {
      // communicate in +mu direction
      ghost_sendrecv_PRECISION_new( out->vector_buffer, mu, +1, &(op->c), plus_dir_param, l );    
    }
    g.num_vect_pass2 = nvec_in;
    for ( mu=0; mu<4; mu++ ) {
      // wait for -mu direction
      ghost_wait_PRECISION_new( in->vector_buffer, mu, -1, &(op->c), minus_dir_param, l );    
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
    coarse_hopp_PRECISION_new( &out_pt, &in_pt, D_pt, l );
    
    D_pt += num_link_var;
    in_pt.vector_buffer = in->vector_buffer + num_site_var*op->neighbor_table[index+Z]*nvec_in;
    coarse_hopp_PRECISION_new( &out_pt, &in_pt, D_pt, l );
    
    D_pt += num_link_var;
    in_pt.vector_buffer = in->vector_buffer + num_site_var*op->neighbor_table[index+Y]*nvec_in;
    coarse_hopp_PRECISION_new( &out_pt, &in_pt, D_pt, l );
    
    D_pt += num_link_var;
    in_pt.vector_buffer = in->vector_buffer + num_site_var*op->neighbor_table[index+X]*nvec_in;
    coarse_hopp_PRECISION_new( &out_pt, &in_pt, D_pt, l );
  }
  
  START_LOCKED_MASTER(threading)
  if ( op->c.comm ) {
    g.num_vect_pass2 = nvec_out;
    for ( mu=0; mu<4; mu++ ) {
      // wait for +mu direction
      ghost_wait_PRECISION_new( out->vector_buffer, mu, +1, &(op->c), plus_dir_param, l );
    }
  }
  END_LOCKED_MASTER(threading)

  END_NO_HYPERTHREADS(threading)
}

// out_o/e -= D_eo/oe*in_e/o
void coarse_n_hopping_term_PRECISION_new( vector_PRECISION *out, vector_PRECISION *in, operator_PRECISION_struct *op,
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
  
  if ( nvec_out < nvec )
    error0("coarse_n_hopping_term_PRECISION: assumptions are not met\n");

  vector_PRECISION_duplicate( &in_pt, in, 0, l );
  vector_PRECISION_duplicate( &out_pt, out, 0, l );
  g.num_vect_pass1 = nvec;

  int core_start;
  int core_end;
  
  // assumptions (1) self coupling has already been performed
  //          OR (2) "out" is initialized with zeros
  set_boundary_PRECISION_new( out, 0, l, threading );
  
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
      ghost_sendrecv_PRECISION_new( in->vector_buffer, mu, -1, &(op->c), minus_dir_param, l );
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
    coarse_n_daggered_hopp_PRECISION_new( &out_pt, &in_pt, D_pt, l );
  }
  SYNC_CORES(threading)
  for ( i=core_start; i<core_end; i++ ) {
    index = 5*i;
    in_pt.vector_buffer = in->vector_buffer + num_site_var*op->neighbor_table[index]*nvec_in;
    D_pt = op->D + num_4link_var*op->neighbor_table[index] + 1*num_link_var;
    index++;
    out_pt.vector_buffer = out->vector_buffer + num_site_var*op->neighbor_table[index+Z]*nvec_out;
    coarse_n_daggered_hopp_PRECISION_new( &out_pt, &in_pt, D_pt, l );
  }
  SYNC_CORES(threading)
  for ( i=core_start; i<core_end; i++ ) {
    index = 5*i;
    in_pt.vector_buffer = in->vector_buffer + num_site_var*op->neighbor_table[index]*nvec_in;
    D_pt = op->D + num_4link_var*op->neighbor_table[index] + 2*num_link_var;
    index++;
    out_pt.vector_buffer = out->vector_buffer + num_site_var*op->neighbor_table[index+Y]*nvec_out;
    coarse_n_daggered_hopp_PRECISION_new( &out_pt, &in_pt, D_pt, l );
  }
  SYNC_CORES(threading)
  for ( i=core_start; i<core_end; i++ ) {
    index = 5*i;
    in_pt.vector_buffer = in->vector_buffer + num_site_var*op->neighbor_table[index]*nvec_in;
    D_pt = op->D + num_4link_var*op->neighbor_table[index] + 3*num_link_var;
    index++;
    out_pt.vector_buffer = out->vector_buffer + num_site_var*op->neighbor_table[index+X]*nvec_out;
    coarse_n_daggered_hopp_PRECISION_new( &out_pt, &in_pt, D_pt, l );
  }
  
  START_LOCKED_MASTER(threading)
  if ( op->c.comm ) {
    g.num_vect_pass2 = nvec_out;
    for ( mu=0; mu<4; mu++ ) {
      // communicate in +mu direction
      ghost_sendrecv_PRECISION_new( out->vector_buffer, mu, +1, &(op->c), plus_dir_param, l );    
    }
    g.num_vect_pass2 = nvec_in;
    for ( mu=0; mu<4; mu++ ) {
      // wait for -mu direction
      ghost_wait_PRECISION_new( in->vector_buffer, mu, -1, &(op->c), minus_dir_param, l );    
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
    coarse_n_hopp_PRECISION_new( &out_pt, &in_pt, D_pt, l );
    
    D_pt += num_link_var;
    in_pt.vector_buffer = in->vector_buffer + num_site_var*op->neighbor_table[index+Z]*nvec_in;
    coarse_n_hopp_PRECISION_new( &out_pt, &in_pt, D_pt, l );
    
    D_pt += num_link_var;
    in_pt.vector_buffer = in->vector_buffer + num_site_var*op->neighbor_table[index+Y]*nvec_in;
    coarse_n_hopp_PRECISION_new( &out_pt, &in_pt, D_pt, l );
    
    D_pt += num_link_var;
    in_pt.vector_buffer = in->vector_buffer + num_site_var*op->neighbor_table[index+X]*nvec_in;
    coarse_n_hopp_PRECISION_new( &out_pt, &in_pt, D_pt, l );
  }
  
  START_LOCKED_MASTER(threading)
  if ( op->c.comm ) {
    g.num_vect_pass2 = nvec_out;
    for ( mu=0; mu<4; mu++ ) {
      // wait for +mu direction
      ghost_wait_PRECISION_new( out->vector_buffer, mu, +1, &(op->c), plus_dir_param, l );
    }
  }
  END_LOCKED_MASTER(threading)

  END_NO_HYPERTHREADS(threading)
}

void coarse_apply_schur_complement_PRECISION_new( vector_PRECISION *out, vector_PRECISION *in, operator_PRECISION_struct *op, level_struct *l, struct Thread *threading ) {

  // start and end indices for vector functions depending on thread
  int start;
  int end;
  // compute start and end indices for core
  // this puts zero for all other hyperthreads, so we can call functions below with all hyperthreads
  compute_core_start_end_custom(op->num_even_sites*l->num_lattice_site_var, l->inner_vector_size, &start, &end, l, threading, l->num_lattice_site_var );//1 );

  vector_PRECISION *tmp = op->buffer;
  for ( int i=0; i<2; i++ ) tmp[i].num_vect_now = in->num_vect_now;

  SYNC_CORES(threading)
  PROF_PRECISION_START( _SC, threading );
  coarse_diag_ee_PRECISION_new( out, in, op, l, threading );
  PROF_PRECISION_STOP( _SC, 0, threading );
  SYNC_CORES(threading)
  vector_PRECISION_define_new( &tmp[0], 0, start, end, l );
  SYNC_CORES(threading)
  PROF_PRECISION_START( _NC, threading );
  coarse_hopping_term_PRECISION_new( &tmp[0], in, op, _ODD_SITES, l, threading );
  PROF_PRECISION_STOP( _NC, 0, threading );
  PROF_PRECISION_START( _SC, threading );
  coarse_diag_oo_inv_PRECISION_new( &tmp[1], &tmp[0], op, l, threading );
  PROF_PRECISION_STOP( _SC, 1, threading );
  PROF_PRECISION_START( _NC, threading );
  coarse_n_hopping_term_PRECISION_new( out, &tmp[1], op, _EVEN_SITES, l, threading );
  PROF_PRECISION_STOP( _NC, 1, threading );
}

void coarse_diag_ee_PRECISION_new( vector_PRECISION *y, vector_PRECISION *x, operator_PRECISION_struct *op, level_struct *l, struct Thread *threading ) {

  int start, end;
  compute_core_start_end_custom( 0, op->num_even_sites, &start, &end, l, threading, 1 );
  // even sites
  coarse_self_couplings_PRECISION_new( y, x, op, start, end, l );
}

// used only in test routines
void coarse_diag_oo_PRECISION_new( vector_PRECISION *y, vector_PRECISION *x, operator_PRECISION_struct *op, level_struct *l, struct Thread *threading ) {
  
  int start, end;
  int nvec_x = x->num_vect, nvec_y = y->num_vect;

  vector_PRECISION x_pt, y_pt;
  int num_site_var=l->num_lattice_site_var;
  int oo_inv_size = SQUARE(num_site_var);
  
  if ( nvec_y < x->num_vect_now )
    error0("coarse_diag_oo_PRECISION: assumptions are not met\n");

/*#ifdef HAVE_TM1p1
  config_PRECISION sc = (g.n_flavours==2) ? op->clover_doublet_oo_inv:op->clover_oo_inv;
#else*/
  config_PRECISION sc = op->clover_oo_inv;
//#endif
  
  compute_core_start_end_custom( 0, op->num_odd_sites, &start, &end, l, threading, 1 );

  vector_PRECISION_duplicate( &x_pt, x, op->num_even_sites+start, l );
  vector_PRECISION_duplicate( &y_pt, y, op->num_even_sites+start, l );
  sc += oo_inv_size*start;

  coarse_LU_multiply_PRECISION_new( &y_pt, &x_pt, sc, start, end, l );

}

static void coarse_diag_oo_inv_PRECISION_new( vector_PRECISION *y, vector_PRECISION *x, operator_PRECISION_struct *op, 
					      level_struct *l, struct Thread *threading ) {

  int start, end, j, jj;
  int nvec_x = x->num_vect, nvec_y = y->num_vect;
  vector_PRECISION x_pt, y_pt;

  if ( nvec_y < x->num_vect_now )
    error0("coarse_diag_oo_inv_PRECISION: assumptions are not met\n");

  compute_core_start_end_custom( 0, op->num_odd_sites, &start, &end, l, threading, 1 );
  vector_PRECISION_duplicate( &x_pt, x, op->num_even_sites+start, l );
  vector_PRECISION_duplicate( &y_pt, y, op->num_even_sites+start, l );

  // odd sites
  int num_site_var = l->num_lattice_site_var;
  int oo_inv_size  = SQUARE(num_site_var);

/*#ifdef HAVE_TM1p1
  config_PRECISION sc = (g.n_flavours==2) ? op->clover_doublet_oo_inv:op->clover_oo_inv;
#else*/
  config_PRECISION sc = op->clover_oo_inv;
//#endif

  sc += oo_inv_size*start;//src were all 0
  coarse_perform_fwd_bwd_subs_PRECISION_new( &y_pt, &x_pt, sc, start, end, l );

}

void coarse_solve_odd_even_PRECISION_new( gmres_PRECISION_struct *p, operator_PRECISION_struct *op, level_struct *l, struct Thread *threading ) {
  /****************************************************************************************
   * Descripton: Solve D*x = b in the even-odd preconditioning using Schur complement D_sc
   *   D_sc = D_ee - D_eo D_oo ^{-1} D_oe
   *   x_e = D_sc^{-1} (b_e - D_eo*D_oo^{-1}*b_o)
   *   x_o = D_oo^{-1} (b_o - D_oe*x_e)  
   ****************************************************************************************/

  //p->b.num_vect_now = g.num_vect_now; p->x.num_vect_now = g.num_vect_now;//!!!!!!
  p->b.num_vect_now = num_loop; p->x.num_vect_now = num_loop;//!!!!!!

  // solve for x_e
  SYNC_CORES(threading)
  PROF_PRECISION_START( _SC, threading ); 
  coarse_diag_oo_inv_PRECISION_new( &p->x, &p->b, op, l, threading ); // x_o = D_oo^{-1}*b_o
  PROF_PRECISION_STOP( _SC, 0, threading );
  PROF_PRECISION_START( _NC, threading );
  coarse_n_hopping_term_PRECISION_new( &p->b, &p->x, op, _EVEN_SITES, l, threading ); // b_e = b_e - D_eo*D_oo^{-1}*b_o
  PROF_PRECISION_STOP( _NC, 0, threading );

  solver_PRECISION( p, l, threading ); // x_e = D_sc^{-1} (b_e - D_eo*D_oo^{-1}*b_o)
  
  // construct x_o from x_e
  PROF_PRECISION_START( _NC, threading );
  coarse_n_hopping_term_PRECISION_new( &p->b, &p->x, op, _ODD_SITES, l, threading ); // b_o = b_o - D_oe*x_e
  PROF_PRECISION_STOP( _NC, 1, threading );
  PROF_PRECISION_START( _SC, threading );
  coarse_diag_oo_inv_PRECISION_new( &p->x, &p->b, op, l, threading ); // x_o = D_oo^{-1} (b_o - D_oe*x_e)
  PROF_PRECISION_STOP( _SC, 1, threading );
  SYNC_CORES(threading)
}
#if 0
void coarse_fabulous_solve_odd_even_PRECISION( fabulous_PRECISION_struct *fab, gmres_PRECISION_struct *p, struct Thread *threading ) {
  /****************************************************************************************
   * Descripton: Solve D*x = b in the even-odd preconditioning using Schur complement D_sc
   *   D_sc = D_ee - D_eo D_oo ^{-1} D_oe
   *   x_e = D_sc^{-1} (b_e - D_eo*D_oo^{-1}*b_o)
   *   x_o = D_oo^{-1} (b_o - D_oe*x_e)  
   ****************************************************************************************/
  level_struct *l = fab->l;
  p->b.num_vect_now = num_loop; p->x.num_vect_now = num_loop;//!!!!!!

  // solve for x_e
  SYNC_CORES(threading)
  PROF_PRECISION_START( _SC, threading );
  coarse_diag_oo_inv_PRECISION_new( &p->x, &p->b, fab->op, l, threading ); // x_o = D_oo^{-1}*b_o
  PROF_PRECISION_STOP( _SC, 0, threading );
  PROF_PRECISION_START( _NC, threading );
  coarse_n_hopping_term_PRECISION_new( &p->b, &p->x, fab->op, _EVEN_SITES, l, threading ); // b_e = b_e - D_eo*D_oo^{-1}*b_o
  PROF_PRECISION_STOP( _NC, 0, threading );
  
  fabulous_PRECISION( fab, p, threading ); // x_e = D_sc^{-1} (b_e - D_eo*D_oo^{-1}*b_o)
  
  // construct x_o from x_e
  PROF_PRECISION_START( _NC, threading );
  coarse_n_hopping_term_PRECISION_new( &p->b, &p->x, fab->op, _ODD_SITES, l, threading ); // b_o = b_o - D_oe*x_e
  PROF_PRECISION_STOP( _NC, 1, threading );
  PROF_PRECISION_START( _SC, threading );
  coarse_diag_oo_inv_PRECISION_new( &p->x, &p->b, fab->op, l, threading ); // x_o = D_oo^{-1} (b_o - D_oe*x_e)
  PROF_PRECISION_STOP( _SC, 1, threading );
  SYNC_CORES(threading)
}
#endif
/************************  TEST ROUTINES  ********************************************/

void coarse_odd_even_PRECISION_test_new( vector_PRECISION *out, vector_PRECISION *in, level_struct *l, struct Thread *threading ) {
  
  if ( g.odd_even ) {
    vector_PRECISION buf[2];

    for(int i=0; i<2; i++){
      vector_PRECISION_init( &buf[i] );
      vector_PRECISION_alloc( &buf[i], _ORDINARY, in->num_vect, l, threading );
      buf[i].num_vect_now = in->num_vect_now;
    }
    
    START_LOCKED_MASTER(threading)
    // transformation part
    vector_PRECISION_copy_new( &buf[0], in, 0, l->inner_vector_size, l );
    // even to odd
    vector_PRECISION_define_new( out, 0, l->oe_op_PRECISION.num_even_sites*l->num_lattice_site_var, l->inner_vector_size, l );
    END_LOCKED_MASTER(threading)

    coarse_hopping_term_PRECISION_new( out, &buf[0], &(l->oe_op_PRECISION), _ODD_SITES, l, threading );
    coarse_diag_oo_inv_PRECISION_new( &buf[1], out, &(l->oe_op_PRECISION), l, threading );

    START_LOCKED_MASTER(threading)
    vector_PRECISION_plus_new( &buf[0], &buf[0], &buf[1], l->oe_op_PRECISION.num_even_sites*l->num_lattice_site_var, l->inner_vector_size, l );
    END_LOCKED_MASTER(threading)
    
    // block diagonal part
    coarse_apply_schur_complement_PRECISION_new( out, &buf[0], &(l->oe_op_PRECISION), l, threading );
    
    coarse_diag_oo_PRECISION_new( out, &buf[0], &(l->oe_op_PRECISION), l, threading );
    
    // back transformation part
    //printf("100\n");fflush(stdout);
    coarse_diag_oo_inv_PRECISION_new( &buf[1], out, &(l->oe_op_PRECISION), l, threading );
    //printf("111\n");fflush(stdout);
    coarse_hopping_term_PRECISION_new( out, &buf[1], &(l->oe_op_PRECISION), _EVEN_SITES, l, threading );
 
    for(int i=0; i<2; i++)
      vector_PRECISION_free( &buf[i], l, threading );                      
  }
}
