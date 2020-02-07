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
 * 1st cleanup: 12/22/2019
 */

#include "main.h"

static void block_diag_ee_PRECISION_new( vector_PRECISION *eta, vector_PRECISION *phi, int start, schwarz_PRECISION_struct *s, level_struct *l, struct Thread *threading );
static void block_diag_oo_PRECISION_new( vector_PRECISION *eta, vector_PRECISION *phi, int start, schwarz_PRECISION_struct *s, level_struct *l, struct Thread *threading );
static void block_diag_oo_inv_PRECISION_new( vector_PRECISION *eta, vector_PRECISION *phi, int start, schwarz_PRECISION_struct *s, level_struct *l, struct Thread *threading ) ;

void schwarz_PRECISION_oddeven_setup( schwarz_PRECISION_struct *s, level_struct *l ) {
  
  int mu, i, d0, c0, b0, a0, d1, c1, b1, a1, t, z, y, x, agg_split[4], block_split[4], block_size[4];
  operator_PRECISION_struct *op = &(s->op);
  int n1 = s->num_block_even_sites;
#ifdef HAVE_TM
  config_PRECISION tm_term_pt = op->tm_term;
#endif

  for ( mu=0; mu<4; mu++ ) {
    agg_split[mu] = l->local_lattice[mu]/l->coarsening[mu];
    block_split[mu] = l->coarsening[mu]/l->block_lattice[mu];
    block_size[mu] = l->block_lattice[mu];
  }
  
  if ( g.csw ) {
    config_PRECISION clover_pt = op->clover, clover_oo_inv_pt = op->clover_oo_inv;
    complex_double buffer[42];
    int cs = 42;
    for ( d0=0; d0<agg_split[T]; d0++ )
      for ( c0=0; c0<agg_split[Z]; c0++ )
        for ( b0=0; b0<agg_split[Y]; b0++ )
          for ( a0=0; a0<agg_split[X]; a0++ )
            
            for ( d1=d0*block_split[T]; d1<(d0+1)*block_split[T]; d1++ )
              for ( c1=c0*block_split[Z]; c1<(c0+1)*block_split[Z]; c1++ )
                for ( b1=b0*block_split[Y]; b1<(b0+1)*block_split[Y]; b1++ )
                  for ( a1=a0*block_split[X]; a1<(a0+1)*block_split[X]; a1++ ) {

                    // skipping even sites
                    clover_pt += n1*cs;
#ifdef HAVE_TM
                    tm_term_pt += n1*12;
#endif
                    for ( t=d1*block_size[T]; t<(d1+1)*block_size[T]; t++ )
                      for ( z=c1*block_size[Z]; z<(c1+1)*block_size[Z]; z++ )
                        for ( y=b1*block_size[Y]; y<(b1+1)*block_size[Y]; y++ )
                          for ( x=a1*block_size[X]; x<(a1+1)*block_size[X]; x++ ) {
                            if (((t-d1*block_size[T])+(z-c1*block_size[Z])+
                                (y-b1*block_size[Y])+(x-a1*block_size[X]))%2 == 1 ) {
                              for ( i=0; i<42; i++ )
                                buffer[i] = (complex_double)clover_pt[i];
#ifdef HAVE_TM
                              if (g.mu + g.mu_odd_shift != 0.0 || g.mu + g.mu_even_shift != 0.0 )
                                for ( i=0; i<12; i++ )
                                  buffer[i] += (complex_double)tm_term_pt[i];
                              tm_term_pt += 12;
                              selfcoupling_LU_decomposition_PRECISION( clover_oo_inv_pt, buffer );
                              clover_oo_inv_pt += 72;
#else  
                              selfcoupling_cholesky_decomposition_PRECISION( clover_oo_inv_pt, buffer );
                              clover_oo_inv_pt += 42;
#endif
                              clover_pt += cs;

                            }
                          }
                  }
  }

#ifdef HAVE_TM1p1
  complex_double buffer[66];
  config_PRECISION clover_oo_inv_pt = op->clover_doublet_oo_inv, clover_pt = op->clover;
  int cs = g.csw ? 42:12;
  config_PRECISION eps_term_pt = op->epsbar_term;
#ifdef HAVE_TM
  tm_term_pt = op->tm_term;
#endif
  
  for ( d0=0; d0<agg_split[T]; d0++ )
    for ( c0=0; c0<agg_split[Z]; c0++ )
      for ( b0=0; b0<agg_split[Y]; b0++ )
        for ( a0=0; a0<agg_split[X]; a0++ )
          
          for ( d1=d0*block_split[T]; d1<(d0+1)*block_split[T]; d1++ )
            for ( c1=c0*block_split[Z]; c1<(c0+1)*block_split[Z]; c1++ )
              for ( b1=b0*block_split[Y]; b1<(b0+1)*block_split[Y]; b1++ )
                for ( a1=a0*block_split[X]; a1<(a0+1)*block_split[X]; a1++ ) {

                  // skipping even sites
                  clover_pt += n1*cs;
                  eps_term_pt += n1*12;
#ifdef HAVE_TM
                  tm_term_pt += n1*12;
#endif
                  for ( t=d1*block_size[T]; t<(d1+1)*block_size[T]; t++ )
                    for ( z=c1*block_size[Z]; z<(c1+1)*block_size[Z]; z++ )
                      for ( y=b1*block_size[Y]; y<(b1+1)*block_size[Y]; y++ )
                        for ( x=a1*block_size[X]; x<(a1+1)*block_size[X]; x++ ) {
                          if (((t-d1*block_size[T])+(z-c1*block_size[Z])+
                               (y-b1*block_size[Y])+(x-a1*block_size[X]))%2 == 1 ) {
                            if ( g.csw ) {
                              for( i=0; i<12; i++ ) //0-23
                                buffer[i+12] = buffer[i] = (complex_double) clover_pt[i];
                              for( i=12; i<42; i++ ) //24-53
                                buffer[i+12] = (complex_double) clover_pt[i];
                            } else {
                              for( i=0; i<12; i++ ) //0-23
                                buffer[i+12] = buffer[i] = (complex_double) clover_pt[i];
                              for( i=12; i<42; i++ ) //24-53
                                buffer[i+12] = _COMPLEX_double_ZERO;
                            }              
                            for( i=0; i<12; i++ ) //54-65
                              buffer[i+54] = (complex_double) eps_term_pt[i];
#ifdef HAVE_TM
                            if (g.mu + g.mu_odd_shift != 0.0 || g.mu + g.mu_even_shift != 0.0 )
                              for(int i=0; i<12; i++) { //0-23
                                buffer[i] += (complex_double) tm_term_pt[i];
                                buffer[i+12] -= (complex_double) tm_term_pt[i];
                              }
                            tm_term_pt += 12;
#endif
                            eps_term_pt += 12;
                            clover_pt += cs;
                            selfcoupling_LU_doublet_decomposition_PRECISION( clover_oo_inv_pt, buffer );
                            clover_oo_inv_pt += 288;
                          }
                        }
                }
#endif  
}

// eta <- block_hopping_term_PRECISION*phi
void block_hopping_term_PRECISION_new( vector_PRECISION *eta, vector_PRECISION *phi,
				       int start, int amount, schwarz_PRECISION_struct *s, level_struct *l, struct Thread *threading ) {
  
  START_UNTHREADED_FUNCTION(threading)

  int a1, a2, n1, n2;
  int *length_even = s->dir_length_even;
  int*length_odd   = s->dir_length_odd;
  int **index      = s->oe_index;
  int *neighbor    = s->op.neighbor_table;
  int nv           = l->num_lattice_site_var;
  int nvec = phi->num_vect_now, nvec_phi = phi->num_vect, nvec_eta = eta->num_vect;
  config_PRECISION D = s->op.D + (start/nv)*36;
  int i, j, k, *ind;
  config_PRECISION D_pt; 
  buffer_PRECISION lphi = phi->vector_buffer+start*nvec_phi, leta = eta->vector_buffer+start*nvec_eta;
  //  printf("block_hopping_term_PRECISION: %d %d %d\n",nvec,nvec_phi,nvec_eta);

  if ( nvec_eta < nvec )
    error0("block_hopping_term_PRECISION: assumptions are not met\n");

/*#ifdef HAVE_TM1p1  
  if ( g.n_flavours == 2 ) {
    complex_PRECISION buf1[25] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, *buf2 = buf1+12;
    // T direction
    if ( amount == _EVEN_SITES ) {
      a1 = 0; n1 = length_even[T];
      a2 = n1; n2 = a2 + length_odd[T];
    } else if ( amount == _ODD_SITES ) {
      a1 = length_even[T]; n1 = a1 + length_odd[T];
      a2 = 0; n2 = a1;
    } else {
      a1 = 0; n1 = length_even[T]+length_odd[T];
      a2 = 0; n2 = n1;
    }
    // "amount" of a block, +T coupling
    ind = index[T];
    for ( i=a1; i<n1; i++ ) {
      k = ind[i]; j = neighbor[4*k+T]; D_pt = D + 36*k + 9*T;
      dprp_T_PRECISION( buf1, lphi+24*j );
      mvm_PRECISION( buf2, D_pt, buf1 );
      mvm_PRECISION( buf2+3, D_pt, buf1+3 );
      mvm_PRECISION( buf2+6, D_pt, buf1+6 );
      mvm_PRECISION( buf2+9, D_pt, buf1+9 );
      dpbp_su3_T_PRECISION( buf2, leta+24*k );
    }
    // "amount" of a block, -T coupling
    for ( i=a2; i<n2; i++ ) {
      k = ind[i]; j = neighbor[4*k+T]; D_pt = D + 36*k + 9*T;
      dprn_T_PRECISION( buf1, lphi+24*k );
      mvmh_PRECISION( buf2, D_pt, buf1 );
      mvmh_PRECISION( buf2+3, D_pt, buf1+3 );
      mvmh_PRECISION( buf2+6, D_pt, buf1+6 );
      mvmh_PRECISION( buf2+9, D_pt, buf1+9 );
      dpbn_su3_T_PRECISION( buf2, leta+24*j );
    }
  
    // Z direction
    if ( amount == _EVEN_SITES ) {
      a1 = 0; n1 = length_even[Z];
      a2 = n1; n2 = a2 + length_odd[Z];
    } else if ( amount == _ODD_SITES ) {
      a1 = length_even[Z]; n1 = a1 + length_odd[Z];
      a2 = 0; n2 = a1;
    } else {
      a1 = 0; n1 = length_even[Z]+length_odd[Z];
      a2 = 0; n2 = n1;
    }
    // "amount" of a block, +Z coupling
    ind = index[Z];
    for ( i=a1; i<n1; i++ ) {
      k = ind[i]; j = neighbor[4*k+Z]; D_pt = D + 36*k + 9*Z;
      dprp_Z_PRECISION( buf1, lphi+24*j );
      mvm_PRECISION( buf2, D_pt, buf1 );
      mvm_PRECISION( buf2+3, D_pt, buf1+3 );
      mvm_PRECISION( buf2+6, D_pt, buf1+6 );
      mvm_PRECISION( buf2+9, D_pt, buf1+9 );
      dpbp_su3_Z_PRECISION( buf2, leta+24*k );
    }
    // "amount" of a block, -Z coupling
    for ( i=a2; i<n2; i++ ) {
      k = ind[i]; j = neighbor[4*k+Z]; D_pt = D + 36*k + 9*Z;
      dprn_Z_PRECISION( buf1, lphi+24*k );
      mvmh_PRECISION( buf2, D_pt, buf1 );
      mvmh_PRECISION( buf2+3, D_pt, buf1+3 );
      mvmh_PRECISION( buf2+6, D_pt, buf1+6 );
      mvmh_PRECISION( buf2+9, D_pt, buf1+9 );
      dpbn_su3_Z_PRECISION( buf2, leta+24*j );
    }
  
    // Y direction
    if ( amount == _EVEN_SITES ) {
      a1 = 0; n1 = length_even[Y];
      a2 = n1; n2 = a2 + length_odd[Y];
    } else if ( amount == _ODD_SITES ) {
      a1 = length_even[Y]; n1 = a1 + length_odd[Y];
      a2 = 0; n2 = a1;
    } else {
      a1 = 0; n1 = length_even[Y]+length_odd[Y];
      a2 = 0; n2 = n1;
    }
    // "amount" of a block, +Y coupling
    ind = index[Y];
    for ( i=a1; i<n1; i++ ) {
      k = ind[i]; j = neighbor[4*k+Y]; D_pt = D + 36*k + 9*Y;
      dprp_Y_PRECISION( buf1, lphi+24*j );
      mvm_PRECISION( buf2, D_pt, buf1 );
      mvm_PRECISION( buf2+3, D_pt, buf1+3 );
      mvm_PRECISION( buf2+6, D_pt, buf1+6 );
      mvm_PRECISION( buf2+9, D_pt, buf1+9 );
      dpbp_su3_Y_PRECISION( buf2, leta+24*k );
    }
    // "amount" of a block, -Y coupling
    for ( i=a2; i<n2; i++ ) {
      k = ind[i]; j = neighbor[4*k+Y]; D_pt = D + 36*k + 9*Y;
      dprn_Y_PRECISION( buf1, lphi+24*k );
      mvmh_PRECISION( buf2, D_pt, buf1 );
      mvmh_PRECISION( buf2+3, D_pt, buf1+3 );
      mvmh_PRECISION( buf2+6, D_pt, buf1+6 );
      mvmh_PRECISION( buf2+9, D_pt, buf1+9 );
      dpbn_su3_Y_PRECISION( buf2, leta+24*j );
    }
    
    // X direction
    if ( amount == _EVEN_SITES ) {
      a1 = 0; n1 = length_even[X];
      a2 = n1; n2 = a2 + length_odd[X];
    } else if ( amount == _ODD_SITES ) {
      a1 = length_even[X]; n1 = a1 + length_odd[X];
      a2 = 0; n2 = a1;
    } else {
      a1 = 0; n1 = length_even[X]+length_odd[X];
      a2 = 0; n2 = n1;
    }
    // "amount" of a block, +X coupling
    ind = index[X];
    for ( i=a1; i<n1; i++ ) {
      k = ind[i]; j = neighbor[4*k+X]; D_pt = D + 36*k + 9*X;
      dprp_X_PRECISION( buf1, lphi+24*j );
      mvm_PRECISION( buf2, D_pt, buf1 );
      mvm_PRECISION( buf2+3, D_pt, buf1+3 );
      mvm_PRECISION( buf2+6, D_pt, buf1+6 );
      mvm_PRECISION( buf2+9, D_pt, buf1+9 );
      dpbp_su3_X_PRECISION( buf2, leta+24*k );
    }
    // "amount" of a block, -X coupling
    for ( i=a2; i<n2; i++ ) {
      k = ind[i]; j = neighbor[4*k+X]; D_pt = D + 36*k + 9*X;
      dprn_X_PRECISION( buf1, lphi+24*k );
      mvmh_PRECISION( buf2, D_pt, buf1 );
      mvmh_PRECISION( buf2+3, D_pt, buf1+3 );
      mvmh_PRECISION( buf2+6, D_pt, buf1+6 );
      mvmh_PRECISION( buf2+9, D_pt, buf1+9 );
      dpbn_su3_X_PRECISION( buf2, leta+24*j );
    }
  } else {
#endif*/
    //complex_PRECISION buf1[13] = {0,0,0,0,0,0,0,0,0,0,0,0,0}, *buf2 = buf1+6;
    complex_PRECISION buf1[13*nvec], *buf2 = buf1+6*nvec; 
    // T direction
    if ( amount == _EVEN_SITES ) {
      a1 = 0; n1 = length_even[T];
      a2 = n1; n2 = a2 + length_odd[T];
    } else if ( amount == _ODD_SITES ) {
      a1 = length_even[T]; n1 = a1 + length_odd[T];
      a2 = 0; n2 = a1;
    } else {
      a1 = 0; n1 = length_even[T]+length_odd[T];
      a2 = 0; n2 = n1;
    }
  // "amount" of a block, +T coupling
  ind = index[T];
  for ( i=a1; i<n1; i++ ) {
    k = ind[i]; j = neighbor[4*k+T]; D_pt = D + 36*k + 9*T;
    prp_T_PRECISION_new( buf1, lphi+12*j*nvec_phi, nvec, nvec, nvec_phi );
    mvm_PRECISION_new( buf2, D_pt, buf1, nvec, nvec, nvec );
    mvm_PRECISION_new( buf2+3*nvec, D_pt, buf1+3*nvec, nvec, nvec, nvec );
    pbp_su3_T_PRECISION_new( buf2, leta+12*k*nvec_eta, nvec, nvec, nvec_eta );
  }
  // "amount" of a block, -T coupling
  for ( i=a2; i<n2; i++ ) {
    k = ind[i]; j = neighbor[4*k+T]; D_pt = D + 36*k + 9*T;
    prn_T_PRECISION_new( buf1, lphi+12*k*nvec_phi, nvec, nvec, nvec_phi );
    mvmh_PRECISION_new( buf2, D_pt, buf1, nvec, nvec, nvec );
    mvmh_PRECISION_new( buf2+3*nvec, D_pt, buf1+3*nvec, nvec, nvec, nvec );
    pbn_su3_T_PRECISION_new( buf2, leta+12*j*nvec_eta, nvec, nvec, nvec_eta );
  }
  
  // Z direction
  if ( amount == _EVEN_SITES ) {
    a1 = 0; n1 = length_even[Z];
    a2 = n1; n2 = a2 + length_odd[Z];
  } else if ( amount == _ODD_SITES ) {
    a1 = length_even[Z]; n1 = a1 + length_odd[Z];
    a2 = 0; n2 = a1;
  } else {
    a1 = 0; n1 = length_even[Z]+length_odd[Z];
    a2 = 0; n2 = n1;
  }
  // "amount" of a block, +Z coupling
  ind = index[Z];
  for ( i=a1; i<n1; i++ ) {
    k = ind[i]; j = neighbor[4*k+Z]; D_pt = D + 36*k + 9*Z;
    prp_Z_PRECISION_new( buf1, lphi+12*j*nvec_phi, nvec, nvec, nvec_phi );
    mvm_PRECISION_new( buf2, D_pt, buf1, nvec, nvec, nvec );
    mvm_PRECISION_new( buf2+3*nvec, D_pt, buf1+3*nvec, nvec, nvec, nvec );
    pbp_su3_Z_PRECISION_new( buf2, leta+12*k*nvec_eta, nvec, nvec, nvec_eta );
  }
  // "amount" of a block, -Z coupling
  for ( i=a2; i<n2; i++ ) {
    k = ind[i]; j = neighbor[4*k+Z]; D_pt = D + 36*k + 9*Z;
    prn_Z_PRECISION_new( buf1, lphi+12*k*nvec_phi, nvec, nvec, nvec_phi );
    mvmh_PRECISION_new( buf2, D_pt, buf1, nvec, nvec, nvec );
    mvmh_PRECISION_new( buf2+3*nvec, D_pt, buf1+3*nvec, nvec, nvec, nvec );
    pbn_su3_Z_PRECISION_new( buf2, leta+12*j*nvec_eta, nvec, nvec, nvec_eta );
  }
  
  // Y direction
  if ( amount == _EVEN_SITES ) {
    a1 = 0; n1 = length_even[Y];
    a2 = n1; n2 = a2 + length_odd[Y];
  } else if ( amount == _ODD_SITES ) {
    a1 = length_even[Y]; n1 = a1 + length_odd[Y];
    a2 = 0; n2 = a1;
  } else {
    a1 = 0; n1 = length_even[Y]+length_odd[Y];
    a2 = 0; n2 = n1;
  }
  // "amount" of a block, +Y coupling
  ind = index[Y];
  for ( i=a1; i<n1; i++ ) {
    k = ind[i]; j = neighbor[4*k+Y]; D_pt = D + 36*k + 9*Y;
    prp_Y_PRECISION_new( buf1, lphi+12*j*nvec_phi, nvec, nvec, nvec_phi );
    mvm_PRECISION_new( buf2, D_pt, buf1, nvec, nvec, nvec );
    mvm_PRECISION_new( buf2+3*nvec, D_pt, buf1+3*nvec, nvec, nvec, nvec );
    pbp_su3_Y_PRECISION_new( buf2, leta+12*k*nvec_eta, nvec, nvec, nvec_eta );
  }
  // "amount" of a block, -Y coupling
  for ( i=a2; i<n2; i++ ) {
    k = ind[i]; j = neighbor[4*k+Y]; D_pt = D + 36*k + 9*Y;
    prn_Y_PRECISION_new( buf1, lphi+12*k*nvec_phi, nvec, nvec, nvec_phi );
    mvmh_PRECISION_new( buf2, D_pt, buf1, nvec, nvec, nvec );
    mvmh_PRECISION_new( buf2+3*nvec, D_pt, buf1+3*nvec, nvec, nvec, nvec );
    pbn_su3_Y_PRECISION_new( buf2, leta+12*j*nvec_eta, nvec, nvec, nvec_eta );
  }
  
  // X direction
  if ( amount == _EVEN_SITES ) {
    a1 = 0; n1 = length_even[X];
    a2 = n1; n2 = a2 + length_odd[X];
  } else if ( amount == _ODD_SITES ) {
    a1 = length_even[X]; n1 = a1 + length_odd[X];
    a2 = 0; n2 = a1;
  } else {
    a1 = 0; n1 = length_even[X]+length_odd[X];
    a2 = 0; n2 = n1;
  }
  // "amount" of a block, +X coupling
  ind = index[X];
  for ( i=a1; i<n1; i++ ) {
    k = ind[i]; j = neighbor[4*k+X]; D_pt = D + 36*k + 9*X;
    prp_X_PRECISION_new( buf1, lphi+12*j*nvec_phi, nvec, nvec, nvec_phi );
    mvm_PRECISION_new( buf2, D_pt, buf1, nvec, nvec, nvec );
    mvm_PRECISION_new( buf2+3*nvec, D_pt, buf1+3*nvec, nvec, nvec, nvec );
    pbp_su3_X_PRECISION_new( buf2, leta+12*k*nvec_eta, nvec, nvec, nvec_eta );
  }
  // "amount" of a block, -X coupling
  for ( i=a2; i<n2; i++ ) {
    k = ind[i]; j = neighbor[4*k+X]; D_pt = D + 36*k + 9*X;
    prn_X_PRECISION_new( buf1, lphi+12*k*nvec_phi, nvec, nvec, nvec_phi );
    mvmh_PRECISION_new( buf2, D_pt, buf1, nvec, nvec, nvec );
    mvmh_PRECISION_new( buf2+3*nvec, D_pt, buf1+3*nvec, nvec, nvec, nvec );
    pbn_su3_X_PRECISION_new( buf2, leta+12*j*nvec_eta, nvec, nvec, nvec_eta );
  }
/*#ifdef HAVE_TM1p1
  }
#endif*/
  END_UNTHREADED_FUNCTION(threading)
}

// eta <- block_n_hopping_term_PRECISION*phi
void block_n_hopping_term_PRECISION_new( vector_PRECISION *eta, vector_PRECISION *phi,
					 int start, int amount, schwarz_PRECISION_struct *s, level_struct *l, struct Thread *threading ) {
  
  START_UNTHREADED_FUNCTION(threading)

  int a1, a2, n1, n2;
  int *length_even = s->dir_length_even;
  int *length_odd  = s->dir_length_odd;
  int **index   = s->oe_index;
  int *neighbor = s->op.neighbor_table;
  int nv        = l->num_lattice_site_var;
  int nvec = phi->num_vect_now, nvec_phi = phi->num_vect, nvec_eta =eta->num_vect;
  int i, j, k, *ind;
  buffer_PRECISION lphi = phi->vector_buffer+start*nvec_phi, leta = eta->vector_buffer+start*nvec_eta;
  config_PRECISION D_pt, D = s->op.D + (start/nv)*36;

  if ( nvec_eta < nvec )
    error0("block_n_hopping_term_PRECISION: assumptions are not met\n");

/*#ifdef HAVE_TM1p1
  if ( g.n_flavours == 2 ) {
    complex_PRECISION buf1[24], *buf2 = buf1+12;
    // T direction
    if ( amount == _EVEN_SITES ) {
      a1 = 0; n1 = length_even[T];
      a2 = n1; n2 = a2 + length_odd[T];
    } else if ( amount == _ODD_SITES ) {
      a1 = length_even[T]; n1 = a1 + length_odd[T];
      a2 = 0; n2 = a1;
    } else {
      a1 = 0; n1 = length_even[T]+length_odd[T];
      a2 = 0; n2 = n1;
    }
    // "amount" of a block, +T coupling
    ind = index[T];
    for ( i=a1; i<n1; i++ ) {
      k = ind[i]; j = neighbor[4*k+T]; D_pt = D + 36*k + 9*T;
      dprp_T_PRECISION( buf1, lphi+24*j );
      nmvm_PRECISION( buf2, D_pt, buf1 );
      nmvm_PRECISION( buf2+3, D_pt, buf1+3 );
      nmvm_PRECISION( buf2+6, D_pt, buf1+6 );
      nmvm_PRECISION( buf2+9, D_pt, buf1+9 );
      dpbp_su3_T_PRECISION( buf2, leta+24*k );
    }
    // "amount" of a block, -T coupling
    for ( i=a2; i<n2; i++ ) {
      k = ind[i]; j = neighbor[4*k+T]; D_pt = D + 36*k + 9*T;
      dprn_T_PRECISION( buf1, lphi+24*k );
      nmvmh_PRECISION( buf2, D_pt, buf1 );
      nmvmh_PRECISION( buf2+3, D_pt, buf1+3 );
      nmvmh_PRECISION( buf2+6, D_pt, buf1+6 );
      nmvmh_PRECISION( buf2+9, D_pt, buf1+9 );
      dpbn_su3_T_PRECISION( buf2, leta+24*j );
    }
    
    // Z direction
    if ( amount == _EVEN_SITES ) {
      a1 = 0; n1 = length_even[Z];
      a2 = n1; n2 = a2 + length_odd[Z];
    } else if ( amount == _ODD_SITES ) {
      a1 = length_even[Z]; n1 = a1 + length_odd[Z];
      a2 = 0; n2 = a1;
    } else {
      a1 = 0; n1 = length_even[Z]+length_odd[Z];
      a2 = 0; n2 = n1;
    }
    // "amount" of a block, +Z coupling
    ind = index[Z];
    for ( i=a1; i<n1; i++ ) {
      k = ind[i]; j = neighbor[4*k+Z]; D_pt = D + 36*k + 9*Z;
      dprp_Z_PRECISION( buf1, lphi+24*j );
      nmvm_PRECISION( buf2, D_pt, buf1 );
      nmvm_PRECISION( buf2+3, D_pt, buf1+3 );
      nmvm_PRECISION( buf2+6, D_pt, buf1+6 );
      nmvm_PRECISION( buf2+9, D_pt, buf1+9 );
      dpbp_su3_Z_PRECISION( buf2, leta+24*k );
    }
    // "amount" of a block, -Z coupling
    for ( i=a2; i<n2; i++ ) {
      k = ind[i]; j = neighbor[4*k+Z]; D_pt = D + 36*k + 9*Z;
      dprn_Z_PRECISION( buf1, lphi+24*k );
      nmvmh_PRECISION( buf2, D_pt, buf1 );
      nmvmh_PRECISION( buf2+3, D_pt, buf1+3 );
      nmvmh_PRECISION( buf2+6, D_pt, buf1+6 );
      nmvmh_PRECISION( buf2+9, D_pt, buf1+9 );
      dpbn_su3_Z_PRECISION( buf2, leta+24*j );
    }
    
    // Y direction
    if ( amount == _EVEN_SITES ) {
      a1 = 0; n1 = length_even[Y];
      a2 = n1; n2 = a2 + length_odd[Y];
    } else if ( amount == _ODD_SITES ) {
      a1 = length_even[Y]; n1 = a1 + length_odd[Y];
      a2 = 0; n2 = a1;
    } else {
      a1 = 0; n1 = length_even[Y]+length_odd[Y];
      a2 = 0; n2 = n1;
    }
    // "amount" of a block, +Y coupling
    ind = index[Y];
    for ( i=a1; i<n1; i++ ) {
      k = ind[i]; j = neighbor[4*k+Y]; D_pt = D + 36*k + 9*Y;
      dprp_Y_PRECISION( buf1, lphi+24*j );
      nmvm_PRECISION( buf2, D_pt, buf1 );
      nmvm_PRECISION( buf2+3, D_pt, buf1+3 );
      nmvm_PRECISION( buf2+6, D_pt, buf1+6 );
      nmvm_PRECISION( buf2+9, D_pt, buf1+9 );
      dpbp_su3_Y_PRECISION( buf2, leta+24*k );
    }
    // "amount" of a block, -Y coupling
    for ( i=a2; i<n2; i++ ) {
      k = ind[i]; j = neighbor[4*k+Y]; D_pt = D + 36*k + 9*Y;
      dprn_Y_PRECISION( buf1, lphi+24*k );
      nmvmh_PRECISION( buf2, D_pt, buf1 );
      nmvmh_PRECISION( buf2+3, D_pt, buf1+3 );
      nmvmh_PRECISION( buf2+6, D_pt, buf1+6 );
      nmvmh_PRECISION( buf2+9, D_pt, buf1+9 );
      dpbn_su3_Y_PRECISION( buf2, leta+24*j );
    }
    
    // X direction
    if ( amount == _EVEN_SITES ) {
      a1 = 0; n1 = length_even[X];
      a2 = n1; n2 = a2 + length_odd[X];
    } else if ( amount == _ODD_SITES ) {
      a1 = length_even[X]; n1 = a1 + length_odd[X];
      a2 = 0; n2 = a1;
    } else {
      a1 = 0; n1 = length_even[X]+length_odd[X];
      a2 = 0; n2 = n1;
    }
    // "amount" of a block, +X coupling
    ind = index[X];
    for ( i=a1; i<n1; i++ ) {
      k = ind[i]; j = neighbor[4*k+X]; D_pt = D + 36*k + 9*X;
      dprp_X_PRECISION( buf1, lphi+24*j );
      nmvm_PRECISION( buf2, D_pt, buf1 );
      nmvm_PRECISION( buf2+3, D_pt, buf1+3 );
      nmvm_PRECISION( buf2+6, D_pt, buf1+6 );
      nmvm_PRECISION( buf2+9, D_pt, buf1+9 );
      dpbp_su3_X_PRECISION( buf2, leta+24*k );
    }
    // "amount" of a block, -X coupling
    for ( i=a2; i<n2; i++ ) {
      k = ind[i]; j = neighbor[4*k+X]; D_pt = D + 36*k + 9*X;
      dprn_X_PRECISION( buf1, lphi+24*k );
      nmvmh_PRECISION( buf2, D_pt, buf1 );
      nmvmh_PRECISION( buf2+3, D_pt, buf1+3 );
      nmvmh_PRECISION( buf2+6, D_pt, buf1+6 );
      nmvmh_PRECISION( buf2+9, D_pt, buf1+9 );
      dpbn_su3_X_PRECISION( buf2, leta+24*j );
    } 
  } else {
#endif*/
    complex_PRECISION buf1[12*nvec], *buf2 = buf1+6*nvec;
    
    // T direction
    if ( amount == _EVEN_SITES ) {
      a1 = 0; n1 = length_even[T];
      a2 = n1; n2 = a2 + length_odd[T];
    } else if ( amount == _ODD_SITES ) {
      a1 = length_even[T]; n1 = a1 + length_odd[T];
      a2 = 0; n2 = a1;
    } else {
      a1 = 0; n1 = length_even[T]+length_odd[T];
      a2 = 0; n2 = n1;
    }
    // "amount" of a block, +T coupling
    ind = index[T];
    for ( i=a1; i<n1; i++ ) {
      k = ind[i]; j = neighbor[4*k+T]; D_pt = D + 36*k + 9*T;
      prp_T_PRECISION_new( buf1, lphi+12*j*nvec_phi, nvec, nvec, nvec_phi );
      nmvm_PRECISION_new( buf2, D_pt, buf1, nvec, nvec, nvec );
      nmvm_PRECISION_new( buf2+3*nvec, D_pt, buf1+3*nvec, nvec, nvec, nvec );
      pbp_su3_T_PRECISION_new( buf2, leta+12*k*nvec_eta, nvec, nvec, nvec_eta );
    }
    // "amount" of a block, -T coupling
    for ( i=a2; i<n2; i++ ) {
      k = ind[i]; j = neighbor[4*k+T]; D_pt = D + 36*k + 9*T;
      prn_T_PRECISION_new( buf1, lphi+12*k*nvec_phi, nvec, nvec, nvec_phi );
      nmvmh_PRECISION_new( buf2, D_pt, buf1, nvec, nvec, nvec );
      nmvmh_PRECISION_new( buf2+3*nvec, D_pt, buf1+3*nvec, nvec, nvec, nvec );
      pbn_su3_T_PRECISION_new( buf2, leta+12*j*nvec_eta, nvec, nvec, nvec_eta );
    }
    
    // Z direction
    if ( amount == _EVEN_SITES ) {
      a1 = 0; n1 = length_even[Z];
      a2 = n1; n2 = a2 + length_odd[Z];
    } else if ( amount == _ODD_SITES ) {
      a1 = length_even[Z]; n1 = a1 + length_odd[Z];
      a2 = 0; n2 = a1;
    } else {
      a1 = 0; n1 = length_even[Z]+length_odd[Z];
      a2 = 0; n2 = n1;
    }
    // "amount" of a block, +Z coupling
    ind = index[Z];
    for ( i=a1; i<n1; i++ ) {
      k = ind[i]; j = neighbor[4*k+Z]; D_pt = D + 36*k + 9*Z;
      prp_Z_PRECISION_new( buf1, lphi+12*j*nvec_phi, nvec, nvec, nvec_phi );
      nmvm_PRECISION_new( buf2, D_pt, buf1, nvec, nvec, nvec );
      nmvm_PRECISION_new( buf2+3*nvec, D_pt, buf1+3*nvec, nvec, nvec, nvec );
      pbp_su3_Z_PRECISION_new( buf2, leta+12*k*nvec_eta, nvec, nvec, nvec_eta );
    }
    // "amount" of a block, -Z coupling
    for ( i=a2; i<n2; i++ ) {
      k = ind[i]; j = neighbor[4*k+Z]; D_pt = D + 36*k + 9*Z;
      prn_Z_PRECISION_new( buf1, lphi+12*k*nvec_phi, nvec, nvec, nvec_phi );
      nmvmh_PRECISION_new( buf2, D_pt, buf1, nvec, nvec, nvec );
      nmvmh_PRECISION_new( buf2+3*nvec, D_pt, buf1+3*nvec, nvec, nvec, nvec );
      pbn_su3_Z_PRECISION_new( buf2, leta+12*j*nvec_eta, nvec, nvec, nvec_eta );
    }
    
    // Y direction
    if ( amount == _EVEN_SITES ) {
      a1 = 0; n1 = length_even[Y];
      a2 = n1; n2 = a2 + length_odd[Y];
    } else if ( amount == _ODD_SITES ) {
      a1 = length_even[Y]; n1 = a1 + length_odd[Y];
      a2 = 0; n2 = a1;
    } else {
      a1 = 0; n1 = length_even[Y]+length_odd[Y];
      a2 = 0; n2 = n1;
    }
    // "amount" of a block, +Y coupling
    ind = index[Y];
    for ( i=a1; i<n1; i++ ) {
      k = ind[i]; j = neighbor[4*k+Y]; D_pt = D + 36*k + 9*Y;
      prp_Y_PRECISION_new( buf1, lphi+12*j*nvec_phi, nvec, nvec, nvec_phi );
      nmvm_PRECISION_new( buf2, D_pt, buf1, nvec, nvec, nvec );
      nmvm_PRECISION_new( buf2+3*nvec, D_pt, buf1+3*nvec, nvec, nvec, nvec );
      pbp_su3_Y_PRECISION_new( buf2, leta+12*k*nvec_eta, nvec, nvec, nvec_eta );
    }
    // "amount" of a block, -Y coupling
    for ( i=a2; i<n2; i++ ) {
      k = ind[i]; j = neighbor[4*k+Y]; D_pt = D + 36*k + 9*Y;
      prn_Y_PRECISION_new( buf1, lphi+12*k*nvec_phi, nvec, nvec, nvec_phi );
      nmvmh_PRECISION_new( buf2, D_pt, buf1, nvec, nvec, nvec );
      nmvmh_PRECISION_new( buf2+3*nvec, D_pt, buf1+3*nvec, nvec, nvec, nvec );
      pbn_su3_Y_PRECISION_new( buf2, leta+12*j*nvec_eta, nvec, nvec, nvec_eta );
    }
    
    // X direction
    if ( amount == _EVEN_SITES ) {
      a1 = 0; n1 = length_even[X];
      a2 = n1; n2 = a2 + length_odd[X];
    } else if ( amount == _ODD_SITES ) {
      a1 = length_even[X]; n1 = a1 + length_odd[X];
      a2 = 0; n2 = a1;
    } else {
      a1 = 0; n1 = length_even[X]+length_odd[X];
      a2 = 0; n2 = n1;
    }
    // "amount" of a block, +X coupling
    ind = index[X];
    for ( i=a1; i<n1; i++ ) {
      k = ind[i]; j = neighbor[4*k+X]; D_pt = D + 36*k + 9*X;
      prp_X_PRECISION_new( buf1, lphi+12*j*nvec_phi, nvec, nvec, nvec_phi );
      nmvm_PRECISION_new( buf2, D_pt, buf1, nvec, nvec, nvec );
      nmvm_PRECISION_new( buf2+3*nvec, D_pt, buf1+3*nvec, nvec, nvec, nvec );
      pbp_su3_X_PRECISION_new( buf2, leta+12*k*nvec_eta, nvec, nvec, nvec_eta );
    }
    // "amount" of a block, -X coupling
    for ( i=a2; i<n2; i++ ) {
      k = ind[i]; j = neighbor[4*k+X]; D_pt = D + 36*k + 9*X;
      prn_X_PRECISION_new( buf1, lphi+12*k*nvec_phi, nvec, nvec, nvec_phi );
      nmvmh_PRECISION_new( buf2, D_pt, buf1, nvec, nvec, nvec );
      nmvmh_PRECISION_new( buf2+3*nvec, D_pt, buf1+3*nvec, nvec, nvec, nvec );
      pbn_su3_X_PRECISION_new( buf2, leta+12*j*nvec_eta, nvec, nvec, nvec_eta );
    }
/*#ifdef HAVE_TM1p1
  }
#endif*/
  END_UNTHREADED_FUNCTION(threading)
}

// out <- block_schur_complement_PRECISION*in
void apply_block_schur_complement_PRECISION_new( vector_PRECISION *out, vector_PRECISION *in, int start,
						 schwarz_PRECISION_struct *s, level_struct *l, struct Thread *threading ) {
  
  vector_PRECISION *tmp = s->oe_buf;
  for (int i=0; i<2; i++ ) s->oe_buf[i].num_vect_now = in->num_vect_now;

  if ( out->num_vect < in->num_vect_now )
    error0("apply_block_schur_complement_PRECISION: assumptions are not met\n");

  block_diag_ee_PRECISION_new( out, in, start, s, l, threading );
  START_LOCKED_MASTER(threading)
  vector_PRECISION_define_new( &tmp[0], 0, start + l->num_lattice_site_var*s->num_block_even_sites, start + s->block_vector_size, l );
  END_LOCKED_MASTER(threading)
  block_hopping_term_PRECISION_new( &tmp[0], in, start, _ODD_SITES, s, l, threading );
  block_diag_oo_inv_PRECISION_new( &tmp[1], &tmp[0], start, s, l, threading );
  block_n_hopping_term_PRECISION_new( out, &tmp[1], start, _EVEN_SITES, s, l, threading );
}

static void block_diag_ee_PRECISION_new( vector_PRECISION *eta, vector_PRECISION *phi,
					 int start, schwarz_PRECISION_struct *s, level_struct *l, struct Thread *threading ) {
 
  START_UNTHREADED_FUNCTION(threading)  
    int n1 = s->num_block_even_sites, nv = l->num_lattice_site_var;//+s->num_block_odd_sites
  clover_PRECISION_new( eta, phi, &(s->op), start, start+nv*n1, l, no_threading ); //!!!! changed to no_threading

  END_UNTHREADED_FUNCTION(threading)
}

// used only in test routines
// diagonal blocks applied to the odd sites of a block
static void block_diag_oo_PRECISION_new( vector_PRECISION *eta, vector_PRECISION *phi,
					 int start, schwarz_PRECISION_struct *s, level_struct *l, struct Thread *threading ) {

  START_UNTHREADED_FUNCTION(threading)

  int i, n1 = s->num_block_even_sites, n2 = s->num_block_odd_sites;
  int j, jj, nvec = phi->num_vect_now, nvec_phi = phi->num_vect, nvec_eta = eta->num_vect;
  vector_PRECISION leta, lphi;
  vector_PRECISION_duplicate( &leta, eta, n1+start/12, l );
  vector_PRECISION_duplicate( &lphi, phi, n1+start/12, l );
  //  printf("block_diag_oo_: %d\n",l->num_lattice_site_var);
  if ( nvec_eta < nvec )
    error0("block_diag_oo_PRECISION: assunmptions are not met\n");

/*#ifdef HAVE_TM1p1
  if ( g.n_flavours == 2 ) {
    int block_num = start/24/(n1+n2);
    //    config_PRECISION clover = s->op.clover_doublet_oo_inv+n1*288+(start/24)*288;
    config_PRECISION clover = s->op.clover_doublet_oo_inv-(block_num+1)*n1*288;
    LU_multiply_PRECISION( eta, phi, clover, n1*24+start, (n1+n2)*24+start );
  } else {
#endif*/
    if ( g.csw ) {
      int block_num = start/12/(n1+n2);
#ifndef HAVE_TM
      config_PRECISION clover = s->op.clover_oo_inv+(start/12-block_num*n1)*42;
      LLH_multiply_PRECISION_new( &leta, &lphi, clover, 0, n2 );
#else
      config_PRECISION clover =  s->op.clover_oo_inv+(start/12-block_num*n1)*72;
      LU_multiply_PRECISION_new( &leta, &lphi, clover, 0, n2 );
#endif
    } else {
      config_PRECISION clover = s->op.clover+n1*12+start;
#ifndef HAVE_TM
      for ( i=0; i<12*n2; i++ )
        VECTOR_LOOP(j, nvec, jj, leta.vector_buffer[i*nvec_eta+j+jj] = lphi.vector_buffer[i*nvec_phi+j+jj]*(clover[i]);)
#else
      config_PRECISION tm_term = s->op.tm_term+n1*12+start;
      for ( i=0; i<12*n2; i++ )
        VECTOR_LOOP(j, nvec, jj, leta.vector_buffer[i*nvec_eta+j+jj] = lphi.vector_buffer[i*nvec_phi+j+jj]*(clover[i]+tm_term[i]);)
#endif
    }
/*#ifdef HAVE_TM1p1
  }
#endif*/  
  END_UNTHREADED_FUNCTION(threading)
}

// inverted diagonal blocks applied to the odd sites of a block
static void block_diag_oo_inv_PRECISION_new( vector_PRECISION *eta, vector_PRECISION *phi, int start, schwarz_PRECISION_struct *s,
					     level_struct *l, struct Thread *threading ) {

  START_UNTHREADED_FUNCTION(threading)
  
  int i, n1 = s->num_block_even_sites, n2 = s->num_block_odd_sites;
  int j, jj, nvec = phi->num_vect_now, nvec_phi = phi->num_vect, nvec_eta = eta->num_vect;
  vector_PRECISION leta, lphi;
  vector_PRECISION_duplicate( &leta, eta, n1+start/12, l );
  vector_PRECISION_duplicate( &lphi, phi, n1+start/12, l );

  if ( nvec_eta < nvec )
    error0("block_diag_oo_inv_PRECISION: assumptions are not met\n");

/*#ifdef HAVE_TM1p1
  if ( g.n_flavours == 2 ) {

    int block_num = start/24/(n1+n2);
    config_PRECISION clover = s->op.clover_doublet_oo_inv-(block_num+1)*n1*288;
    LU_perform_fwd_bwd_subs_PRECISION( eta, phi, clover, n1*24+start, (n1+n2)*24+start );
  } else {
#endif*/

    if ( g.csw ) {
      int block_num = start/12/(n1+n2);
#ifndef HAVE_TM
      config_PRECISION clover = s->op.clover_oo_inv+(start/12-(block_num)*n1)*42;
      //    LU_perform_fwd_bwd_subs_PRECISION_new( &leta, &lphi, clover, 0, n2 );
#else
      config_PRECISION clover = s->op.clover_oo_inv+(start/12-(block_num)*n1)*72;
      //    LU_perform_fwd_bwd_subs_PRECISION_new( &leta, &lphi, clover, 0, n2 );
#endif
      LU_perform_fwd_bwd_subs_PRECISION_new( &leta, &lphi, clover, 0, n2 );
    } else {
      config_PRECISION clover = s->op.clover+n1*12+start;
#ifndef HAVE_TM
      for ( i=0; i<12*n2; i++ )
        VECTOR_LOOP(j, nvec, jj, leta.vector_buffer[i*nvec_eta+j+jj] = lphi.vector_buffer[i*nvec_phi+j+jj]/(clover[i]);)
#else
      config_PRECISION tm_term = s->op.tm_term+n1*12+start;
      for ( i=0; i<12*n2; i++ )
        VECTOR_LOOP(j, nvec, jj, leta.vector_buffer[i*nvec_eta+j+jj] = lphi.vector_buffer[i*nvec_phi+j+jj]/(clover[i]+tm_term[i]);)
#endif
    }
/*#ifdef HAVE_TM1p1
  }
#endif*/
  
  END_UNTHREADED_FUNCTION(threading)
}

//phi <- D_block_oddeven^-1 r????
void block_solve_oddeven_PRECISION_new( vector_PRECISION *phi, vector_PRECISION *r, vector_PRECISION *latest_iter,
    int start, schwarz_PRECISION_struct *s, level_struct *l, struct Thread *threading ) {
  
  START_UNTHREADED_FUNCTION(threading)

  vector_PRECISION *tmp = s->oe_buf;
  int end = start+s->block_vector_size;

  for ( int i=2; i<4; i++ )
    s->oe_buf[i].num_vect_now = r->num_vect_now;

  if ( phi->num_vect < r->num_vect_now )
    error0("block_solve_oddeven_PRECISION: assumptions are not met\n");

  // odd to even
  vector_PRECISION_copy_new( &tmp[3], r, start, end, l );
  block_diag_oo_inv_PRECISION_new( &tmp[2], &tmp[3], start, s, l, no_threading );
  block_n_hopping_term_PRECISION_new( &tmp[3], &tmp[2], start, _EVEN_SITES, s, l, no_threading );
  local_minres_PRECISION_new( NULL, &tmp[3], &tmp[2], start, s, l, no_threading );  
  // even to odd
  block_n_hopping_term_PRECISION_new( &tmp[3], &tmp[2], start, _ODD_SITES, s, l, no_threading );
  block_diag_oo_inv_PRECISION_new( &tmp[2], &tmp[3], start, s, l, no_threading );
  //  printf0("block_solve_oddeven2\n");
  // update phi, latest_iter
  vector_PRECISION_copy_new( latest_iter, &tmp[2], start, end, l );
  vector_PRECISION_plus_new( phi, phi, &tmp[2], start, end, l );
  // update r
  vector_PRECISION_copy_new( r, &tmp[3], start, start+l->num_lattice_site_var*s->num_block_even_sites, l );
  vector_PRECISION_define_new( r, 0, start+l->num_lattice_site_var*s->num_block_even_sites, end, l );
  //  printf0("block_solve_oddeven3\n");
  END_UNTHREADED_FUNCTION(threading)
}

// layout: even->odd => block
void oddeven_to_block_PRECISION_new( vector_PRECISION *out, vector_PRECISION *in, level_struct *l, struct Thread *threading ) {

  int i, j, k, m, jj, jjj;
  int nvec = in->num_vect_now, nvec_in = in->num_vect, nvec_out = out->num_vect;
  int nsv = l->num_lattice_site_var, *tt_oe = l->oe_op_PRECISION.translation_table, *tt_b = l->s_PRECISION.op.translation_table;
  int start = threading->start_site[l->depth];
  int end   = threading->end_site[l->depth];

  if ( nvec_out < nvec)
    error0("oddeven_to_serial_PRECISION: assumptions are not met\n");

  // this function seems to do some data reordering, barriers ensure that everything is in sync
  SYNC_CORES(threading)
  START_NO_HYPERTHREADS(threading)  
  for ( i=start; i<end; i++ ) {
    k = tt_oe[i]; m = tt_b[i];
    for ( j=0; j<nsv; j++ ) {
      VECTOR_LOOP( jj, nvec, jjj, out->vector_buffer[(m*nsv+j)*nvec_out+jj+jjj] = in->vector_buffer[(k*nsv+j)*nvec_in+jj+jjj]; )
    }
  }
  END_NO_HYPERTHREADS(threading)
  SYNC_CORES(threading)  
}

void block_to_oddeven_PRECISION_new( vector_PRECISION *out, vector_PRECISION *in, level_struct *l, struct Thread *threading ) {

  int i, j, k, m, jj, jjj;
  int nsv = l->num_lattice_site_var, *tt_oe = l->oe_op_PRECISION.translation_table, *tt_b = l->s_PRECISION.op.translation_table;
  int nvec = in->num_vect_now, nvec_in = in->num_vect, nvec_out = out->num_vect;
  int start = threading->start_site[l->depth];
  int end   = threading->end_site[l->depth];

  if ( nvec_out < nvec )
    error0("oddeven_to_serial_PRECISION: assumptions are not met\n");

  // this function seems to do some data reordering, barriers ensure that everything is in sync
  SYNC_CORES(threading)
  START_NO_HYPERTHREADS(threading)  
  for ( i=start; i<end; i++ ) {
    k = tt_oe[i]; m = tt_b[i];
    for ( j=0; j<nsv; j++ ) {
      VECTOR_LOOP( jj, nvec, jjj, out->vector_buffer[(k*nsv+j)*nvec_out+jj+jjj] = in->vector_buffer[(m*nsv+j)*nvec_in+jj+jjj]; )
    }
  }
  END_NO_HYPERTHREADS(threading)
  SYNC_CORES(threading)  
}

/*******************  TEST ROUTINES  *****************************************************/

void block_oddeven_PRECISION_test_new( level_struct *l, struct Thread *threading ) {
  START_UNTHREADED_FUNCTION(threading)

  schwarz_PRECISION_struct *s = &(l->s_PRECISION);
  
  int n_vect = g.num_rhs_vect;
  vector_PRECISION b1,b2,b3,b4,b5;
  PRECISION diff1[n_vect], diff2[n_vect];
  
  vector_PRECISION_init(&b1);
  vector_PRECISION_init(&b2);
  vector_PRECISION_init(&b3);
  vector_PRECISION_init(&b4);
  vector_PRECISION_init(&b5);  
 
  int vs = s->block_vector_size * s->num_blocks;
 
  MALLOC( b1.vector_buffer, complex_PRECISION, vs*n_vect );
  b1.num_vect = n_vect;
  b1.num_vect_now = n_vect;

  MALLOC( b2.vector_buffer, complex_PRECISION, vs*n_vect );
  b2.num_vect = n_vect;
  b2.num_vect_now = n_vect;

  MALLOC( b3.vector_buffer, complex_PRECISION, vs*n_vect );
  b3.num_vect = n_vect;
  b3.num_vect_now = n_vect;

  MALLOC( b4.vector_buffer, complex_PRECISION, vs*n_vect );
  b4.num_vect = n_vect;
  b4.num_vect_now = n_vect;

  MALLOC( b5.vector_buffer, complex_PRECISION, vs*n_vect );
  b5.num_vect = n_vect;
  b5.num_vect_now = n_vect;

  vector_PRECISION_define_random_new( &b1, 0, vs, l );

  for (int i = 0; i< s->num_blocks; i++ ) {//printf("block: %d %d\n",s->block[i].start,s->num_block_even_sites);
    block_diag_ee_PRECISION_new( &b2, &b1, s->block[i].start*l->num_lattice_site_var, s, l, no_threading );
    block_diag_oo_PRECISION_new( &b2, &b1, s->block[i].start*l->num_lattice_site_var, s, l, no_threading );
    block_hopping_term_PRECISION_new( &b2, &b1, s->block[i].start*l->num_lattice_site_var, _FULL_SYSTEM, s, l, no_threading );
    
    block_d_plus_clover_PRECISION_new( &b3, &b1, s->block[i].start*l->num_lattice_site_var, s, l, no_threading );
  }
  //  for(int i=0;i<vs*n_vect;i++)printf0("%d %d %d %13g %g %g\n",i/n_vect/l->num_lattice_site_var,(i/n_vect)%l->num_lattice_site_var,i%n_vect,creal_PRECISION(b3.vector_buffer[i]-b2.vector_buffer[i]),creal_PRECISION(b3.vector_buffer[i]),creal_PRECISION(b2.vector_buffer[i]));   
  vector_PRECISION_minus_new( &b3, &b3, &b2, 0, vs, l );
  global_norm_PRECISION_new( diff1, &b3, 0, vs, l, no_threading );
  global_norm_PRECISION_new( diff2, &b2, 0, vs, l, no_threading );
  
  for (int i = 0; i<n_vect; i++ )
    test0_PRECISION("depth: %d, correctness of block odd even layout: %le\n", l->depth, diff1[i]/diff2[i] );
  
  vector_PRECISION_copy_new( &b4, &b1, 0, s->block_vector_size, l );
  vector_PRECISION_define_new( &b3, 0, l->num_lattice_site_var*s->num_block_even_sites, s->block_vector_size, l );
  
  block_hopping_term_PRECISION_new( &b3, &b4, 0, _ODD_SITES, s, l, no_threading );
  block_diag_oo_inv_PRECISION_new( &b5, &b3, 0, s, l, no_threading );
  vector_PRECISION_plus_new( &b4, &b4, &b5, l->num_lattice_site_var*s->num_block_even_sites, s->block_vector_size, l );
  
  apply_block_schur_complement_PRECISION_new( &b3, &b4, 0, s, l, no_threading );
  block_diag_oo_PRECISION_new( &b3, &b4, 0, s, l, no_threading );
  
  block_diag_oo_inv_PRECISION_new( &b5, &b3, 0, s, l, no_threading );
  block_hopping_term_PRECISION_new( &b3, &b5, 0, _EVEN_SITES, s, l, no_threading );
  
  vector_PRECISION_minus_new( &b3, &b2, &b3, 0, s->block_vector_size, l );
  global_norm_PRECISION_new( diff1, &b3, 0, s->block_vector_size, l, no_threading );
  global_norm_PRECISION_new( diff2, &b2, 0, s->block_vector_size, l, no_threading );
  
  for (int i = 0; i<n_vect; i++ )
    test0_PRECISION("depth: %d, correctness of block odd even schur complement: %le\n", l->depth, diff1[i]/diff2[i] );
  
  FREE( b1.vector_buffer, complex_PRECISION, vs*n_vect );
  FREE( b2.vector_buffer, complex_PRECISION, vs*n_vect );
  FREE( b3.vector_buffer, complex_PRECISION, vs*n_vect );
  FREE( b4.vector_buffer, complex_PRECISION, vs*n_vect );
  FREE( b5.vector_buffer, complex_PRECISION, vs*n_vect );
  //  error0("stop\n");
  END_UNTHREADED_FUNCTION(threading)
}
