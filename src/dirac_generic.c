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
 * copied: 11/29/2019 
 * changed from sbacchio
 * checked: 12/04/2019(some functions TMp1 and update are not corrected.
 * 1st cleanup: 12/18/2019
 */

#include "main.h"

// eta = op->clover(+op->tm)*phi
void clover_PRECISION( vector_PRECISION *eta, vector_PRECISION *phi, operator_PRECISION_struct *op, int start, int end,
                       level_struct *l, struct Thread *threading ) {

#ifdef DEBUG
  if ( eta->num_vect < phi->num_vect_now )
    error0("clover_PRECISION: assumptions are not met\n");
#endif

  int j, jj;
  int nv = l->num_lattice_site_var, nvec = phi->num_vect_now, nvec_phi = phi->num_vect, nvec_eta = eta->num_vect;
  
  buffer_PRECISION lphi = phi->vector_buffer+start*nvec_phi, leta = eta->vector_buffer+start*nvec_eta;
  buffer_PRECISION leta_end = eta->vector_buffer+end*nvec_eta;

#ifdef HAVE_TM
  int nd = nv;
#ifdef HAVE_MULT_TM
  nd *= num_loop;
  config_PRECISION diag = op->tm_term+start*num_loop+num_loop*l->inner_vector_size*g.n_chunk;
#else
  config_PRECISION diag = op->odd_proj+start;
#endif
#endif

#ifdef PROFILING
  START_MASTER(threading)
  PROF_PRECISION_START( _SC );
  END_MASTER(threading)
#endif
    
  if ( g.csw == 0.0 ) { // op->clover contains only the mass term
    config_PRECISION clover = op->clover+(start/nv)*12;
#ifdef HAVE_TM
      //if ( g.mu + g.mu_odd_shift != 0.0 || g.mu + g.mu_even_shift != 0.0 ) {
      if ( g.mu + g.mu_odd_shift != 0.0 || g.is_even_shifted_mu_nonzero ) {
	while ( leta < leta_end ) {
	  set_diag_PRECISION( leta, lphi, clover, diag, op, nvec, nvec_eta, nvec_phi, 1 );
	  leta += nv*nvec_eta; lphi += nv*nvec_phi; clover += nv; diag += nd;
	}
      }
#endif // else
      while ( leta < leta_end )
	FOR12( VECTOR_LOOP(j, nvec, jj, leta[j+jj] = lphi[j+jj]*(*clover);)
	       leta += nvec_eta; lphi += nvec_phi; clover++; )
  } else { // op->clover also contains clover contribution
    config_PRECISION clover = op->clover+(start/nv)*42;
#ifdef HAVE_TM
      //if ( g.mu + g.mu_odd_shift != 0.0 || g.mu + g.mu_even_shift != 0.0 ) {
      if ( g.mu + g.mu_odd_shift != 0.0 || g.is_even_shifted_mu_nonzero ) {
	while ( leta < leta_end ) {
	  site_clover_PRECISION( leta, lphi, clover, nvec, nvec_eta, nvec_phi );
	  set_diag_PRECISION( leta, lphi, NULL, diag, op, nvec, nvec_eta, nvec_phi, 1 );
          leta += nv*nvec_eta; lphi += nv*nvec_phi; clover += 42; diag += nd;
        }
      }
#endif // else
      while ( leta < leta_end ) {
	site_clover_PRECISION( leta, lphi, clover, nvec, nvec_eta, nvec_phi );
	leta += 12*nvec_eta; lphi += 12*nvec_phi;
	clover+=42;
      }
  }

#ifdef HAVE_TM1p1
  config_PRECISION eps_term = op->epsbar_term+(start/nv)*12;  
  //lphi = phi->vector_buffer+start*nvec_phi, leta = eta->vector_buffer+start*nvec_eta;
  if ( g.n_flavours == 2 &&
       ( g.epsbar != 0 || g.epsbar_ig5_odd_shift != 0 || g.epsbar_ig5_odd_shift != 0 ) )
    apply_doublet_coupling_PRECISION(eta, phi, eps_term, start, end);
#endif

#ifdef PROFILING
  START_MASTER(threading)
  PROF_PRECISION_STOP( _SC, 1 );
  END_MASTER(threading)
#endif
    
}

// eta <- d_plus_clover*phi
void d_plus_clover_PRECISION( vector_PRECISION *eta, vector_PRECISION *phi, operator_PRECISION_struct *op, level_struct *l, struct Thread *threading ) {
  
  int n = l->num_inner_lattice_sites, *neighbor = op->neighbor_table, start, end, nv = l->num_lattice_site_var;
  int nvec = phi->num_vect_now, nvec_phi = phi->num_vect, nvec_eta = eta->num_vect, nvec_op = op->pr_num_vect;
  int i, j, *nb_pt;
  buffer_PRECISION phi_pt, eta_pt, end_pt;
  config_PRECISION D_pt;

#ifdef DEBUG
  if ( eta->num_vect < phi->num_vect_now )
    error0("d_plus_clover_PRECISION: assumptions are not met\n");
#endif
  
  compute_core_start_end_custom(0, nv*n, &start, &end, l, threading, nv );

  SYNC_MASTER_TO_ALL(threading)
  clover_PRECISION( eta, phi, op, start, end, l, threading );
  START_MASTER(threading)
  PROF_PRECISION_START( _NC ); 
  END_MASTER(threading)

  complex_PRECISION pbuf[6*nvec];//pay attention to #vectors in op->prn* defined in operator_PRECISION_alloc_projection_buffers!!!!
  for ( i=start*nvec_op/2, phi_pt=phi->vector_buffer+start*nvec_phi; i<end*nvec_op/2; i+=6*nvec_op, phi_pt+=12*nvec_phi ) {
    prp_T_PRECISION( op->prnT+i, phi_pt, nvec, nvec_op, nvec_phi );
    prp_Z_PRECISION( op->prnZ+i, phi_pt, nvec, nvec_op, nvec_phi );
    prp_Y_PRECISION( op->prnY+i, phi_pt, nvec, nvec_op, nvec_phi );
    prp_X_PRECISION( op->prnX+i, phi_pt, nvec, nvec_op, nvec_phi );
  }
  // start communication in negative direction
  START_LOCKED_MASTER(threading) 
  g.num_vect_pass1 = nvec; g.num_vect_pass2 = nvec_op;//temp fix!!!!!!
  ghost_sendrecv_PRECISION( op->prnT, T, -1, &(op->c), _FULL_SYSTEM, l );
  ghost_sendrecv_PRECISION( op->prnZ, Z, -1, &(op->c), _FULL_SYSTEM, l );
  ghost_sendrecv_PRECISION( op->prnY, Y, -1, &(op->c), _FULL_SYSTEM, l );
  ghost_sendrecv_PRECISION( op->prnX, X, -1, &(op->c), _FULL_SYSTEM, l );
  END_LOCKED_MASTER(threading) 

  // project plus dir and multiply with U dagger
  for ( phi_pt=phi->vector_buffer+start*nvec_phi, end_pt=phi->vector_buffer+end*nvec_phi, D_pt = op->D+(start*3), nb_pt=neighbor+((start/12)*4); phi_pt<end_pt; phi_pt+=12*nvec_phi ) {
    // T dir
    j = 6*(*nb_pt)*nvec_op; nb_pt++;
    prn_T_PRECISION( pbuf, phi_pt, nvec, nvec, nvec_phi );
    mvmh_PRECISION( op->prpT+j, D_pt, pbuf, nvec, nvec_op, nvec );
    mvmh_PRECISION( op->prpT+j+3*nvec_op, D_pt, pbuf+3*nvec, nvec, nvec_op, nvec ); D_pt += 9;
    // Z dir
    j = 6*(*nb_pt)*nvec_op; nb_pt++;
    prn_Z_PRECISION( pbuf, phi_pt, nvec, nvec, nvec_phi );
    mvmh_PRECISION( op->prpZ+j, D_pt, pbuf, nvec, nvec_op, nvec );
    mvmh_PRECISION( op->prpZ+j+3*nvec_op, D_pt, pbuf+3*nvec, nvec, nvec_op, nvec ); D_pt += 9;
    // Y dir
    j = 6*(*nb_pt)*nvec_op; nb_pt++;
    prn_Y_PRECISION( pbuf, phi_pt, nvec, nvec, nvec_phi );
    mvmh_PRECISION( op->prpY+j, D_pt, pbuf, nvec, nvec_op, nvec );
    mvmh_PRECISION( op->prpY+j+3*nvec_op, D_pt, pbuf+3*nvec, nvec, nvec_op, nvec ); D_pt += 9;
    // X dir
    j = 6*(*nb_pt)*nvec_op; nb_pt++;
    prn_X_PRECISION( pbuf, phi_pt, nvec, nvec, nvec_phi );
    mvmh_PRECISION( op->prpX+j, D_pt, pbuf, nvec, nvec_op, nvec );
    mvmh_PRECISION( op->prpX+j+3*nvec_op, D_pt, pbuf+3*nvec, nvec, nvec_op, nvec ); D_pt += 9;
  }
  // start communication in positive direction
  START_LOCKED_MASTER(threading)
  ghost_sendrecv_PRECISION( op->prpT, T, +1, &(op->c), _FULL_SYSTEM, l );
  ghost_sendrecv_PRECISION( op->prpZ, Z, +1, &(op->c), _FULL_SYSTEM, l );
  ghost_sendrecv_PRECISION( op->prpY, Y, +1, &(op->c), _FULL_SYSTEM, l );
  ghost_sendrecv_PRECISION( op->prpX, X, +1, &(op->c), _FULL_SYSTEM, l );
  // wait for communication in negative direction
  ghost_wait_PRECISION( op->prnT, T, -1, &(op->c), _FULL_SYSTEM, l );
  ghost_wait_PRECISION( op->prnZ, Z, -1, &(op->c), _FULL_SYSTEM, l );
  ghost_wait_PRECISION( op->prnY, Y, -1, &(op->c), _FULL_SYSTEM, l );
  ghost_wait_PRECISION( op->prnX, X, -1, &(op->c), _FULL_SYSTEM, l );
  END_LOCKED_MASTER(threading)

  // multiply with U and lift up minus dir
  for ( eta_pt=eta->vector_buffer+start*nvec_eta, end_pt=eta->vector_buffer+end*nvec_eta, D_pt = op->D+start*3, nb_pt=neighbor+(start/12)*4; eta_pt<end_pt; eta_pt+=12*nvec_eta ) {
    // T dir
    j = 6*(*nb_pt)*nvec_op; nb_pt++;
    mvm_PRECISION( pbuf, D_pt, op->prnT+j, nvec, nvec, nvec_op );
    mvm_PRECISION( pbuf+3*nvec, D_pt, op->prnT+j+3*nvec_op, nvec, nvec, nvec_op );
    pbp_su3_T_PRECISION( pbuf, eta_pt, nvec, nvec, nvec_eta ); D_pt += 9;
    // Z dir
    j = 6*(*nb_pt)*nvec_op; nb_pt++;
    mvm_PRECISION( pbuf, D_pt, op->prnZ+j, nvec, nvec, nvec_op );
    mvm_PRECISION( pbuf+3*nvec, D_pt, op->prnZ+j+3*nvec_op, nvec, nvec, nvec_op );
    pbp_su3_Z_PRECISION( pbuf, eta_pt, nvec, nvec, nvec_eta ); D_pt += 9;
    // Y dir
    j = 6*(*nb_pt)*nvec_op; nb_pt++;
    mvm_PRECISION( pbuf, D_pt, op->prnY+j, nvec, nvec, nvec_op );
    mvm_PRECISION( pbuf+3*nvec, D_pt, op->prnY+j+3*nvec_op, nvec, nvec, nvec_op );
    pbp_su3_Y_PRECISION( pbuf, eta_pt, nvec, nvec, nvec_eta ); D_pt += 9;
    // X dir
    j = 6*(*nb_pt)*nvec_op; nb_pt++;
    mvm_PRECISION( pbuf, D_pt, op->prnX+j, nvec, nvec, nvec_op );
    mvm_PRECISION( pbuf+3*nvec, D_pt, op->prnX+j+3*nvec_op, nvec, nvec, nvec_op );
    pbp_su3_X_PRECISION( pbuf, eta_pt, nvec, nvec, nvec_eta ); D_pt += 9;
  }
  // wait for communication in positive direction
  START_LOCKED_MASTER(threading)
  ghost_wait_PRECISION( op->prpT, T, +1, &(op->c), _FULL_SYSTEM, l );
  ghost_wait_PRECISION( op->prpZ, Z, +1, &(op->c), _FULL_SYSTEM, l );
  ghost_wait_PRECISION( op->prpY, Y, +1, &(op->c), _FULL_SYSTEM, l );
  ghost_wait_PRECISION( op->prpX, X, +1, &(op->c), _FULL_SYSTEM, l );
  END_LOCKED_MASTER(threading)
  // lift up plus dir
  for ( i=start*nvec_op/2, eta_pt=eta->vector_buffer+start*nvec_eta; i<end*nvec_op/2; i+=6*nvec_op, eta_pt+=12*nvec_eta ) {
    pbn_su3_T_PRECISION( op->prpT+i, eta_pt, nvec, nvec_op, nvec_eta );
    pbn_su3_Z_PRECISION( op->prpZ+i, eta_pt, nvec, nvec_op, nvec_eta );
    pbn_su3_Y_PRECISION( op->prpY+i, eta_pt, nvec, nvec_op, nvec_eta );
    pbn_su3_X_PRECISION( op->prpX+i, eta_pt, nvec, nvec_op, nvec_eta );
  }
  
  START_MASTER(threading)
  PROF_PRECISION_STOP( _NC, 1 );
  END_MASTER(threading)
  SYNC_MASTER_TO_ALL(threading)
}

// eta = blockD*phi
// Assume: sites are order blockwise, i.e., in Schwarz layout
void block_d_plus_clover_PRECISION( vector_PRECISION *eta, vector_PRECISION *phi, int start, schwarz_PRECISION_struct *s, level_struct *l, struct Thread *threading ) {

  START_UNTHREADED_FUNCTION(threading)

#ifdef DEBUG
  if ( eta->num_vect < phi->num_vect_now )
    error0("block_d_plus_clover_PRECISION: assumptions are not met\n");
#endif
  
  int n = s->num_block_sites, *length = s->dir_length, **index = s->index, *neighbor = s->op.neighbor_table, nv = l->num_lattice_site_var;
  int nvec = phi->num_vect_now, nvec_phi = phi->num_vect, nvec_eta = eta->num_vect;
  buffer_PRECISION lphi = phi->vector_buffer+start*nvec_phi, leta = eta->vector_buffer+start*nvec_eta;

  clover_PRECISION(eta, phi, &(s->op), start, start+nv*n, l, no_threading );

  int i, j, k, *ind;
  config_PRECISION D_pt;
  config_PRECISION D = s->op.D + (start/nv)*36;
  //printf0("block d_plus %d %d %d\n",nvec, nvec_phi,nvec_eta);
  //complex_PRECISION buf1[25]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, *buf2=buf1+6, *buf3=buf2+6, *buf4=buf3+6;
  complex_PRECISION buf1[25*nvec], *buf2=buf1+6*nvec, *buf3=buf2+6*nvec, *buf4=buf3+6*nvec;//why 25?????? 4(buffers)x6(projected_spin*color)
  
  // inner block couplings
  ind = index[T]; // T direction
  for ( i=0; i<length[T]; i++ ) {
    k = ind[i]; j = neighbor[4*k+T]; D_pt = D + 36*k + 9*T;
    prn_T_PRECISION( buf1, lphi+12*k*nvec_phi, nvec, nvec, nvec_phi ); // (1+gamma_T) phi(x) + projection
    prp_T_PRECISION( buf2, lphi+12*j*nvec_phi, nvec, nvec, nvec_phi ); // (1-gamma_T) phi(x+hat{T}) + projection
    mvmh_PRECISION( buf3, D_pt, buf1, nvec, nvec, nvec );     // U_T^dagger(x) (1+gamma_T) phi(x) - projected
    mvmh_PRECISION( buf3+3*nvec, D_pt, buf1+3*nvec, nvec, nvec, nvec ); // U_T^dagger(x) (1+gamma_T) phi(x) - projected
    mvm_PRECISION( buf4, D_pt, buf2, nvec, nvec, nvec );     // U_T(x) (1-gamma_T) phi(x+hat{T}) - projected
    mvm_PRECISION( buf4+3*nvec, D_pt, buf2+3*nvec, nvec, nvec, nvec ); // U_T(x) (1-gamma_T) phi(x+hat{T}) - projected
    pbn_su3_T_PRECISION( buf3, leta+12*j*nvec_eta, nvec, nvec, nvec_eta); // eta(x+hat{T}) -= U_T(x)^dagger(x) (1+gamma_T) phi(x) + lift back
    pbp_su3_T_PRECISION( buf4, leta+12*k*nvec_eta, nvec, nvec, nvec_eta ); // eta(x) -= U_T(x) (1-gamma_T) phi(x+hat{T}) + lift back
  }
  ind = index[Z]; // Z direction
  for ( i=0; i<length[Z]; i++ ) {
    k = ind[i]; j = neighbor[4*k+Z]; D_pt = D + 36*k + 9*Z;
    prn_Z_PRECISION( buf1, lphi+12*k*nvec_phi,nvec, nvec, nvec_phi ); // (1+gamma_Z) phi(x) + projection
    prp_Z_PRECISION( buf2, lphi+12*j*nvec_phi,nvec, nvec, nvec_phi ); // (1-gamma_Z) phi(x+hat{Z}) + projection
    mvmh_PRECISION( buf3, D_pt, buf1, nvec, nvec, nvec );     // U_Z^dagger(x) (1+gamma_Z) phi(x) - projected
    mvmh_PRECISION( buf3+3*nvec, D_pt, buf1+3*nvec, nvec, nvec, nvec ); // U_Z^dagger(x) (1+gamma_Z) phi(x) - projected
    mvm_PRECISION( buf4, D_pt, buf2, nvec, nvec, nvec );     // U_Z(x) (1-gamma_Z) phi(x+hat{Z}) - projected
    mvm_PRECISION( buf4+3*nvec, D_pt, buf2+3*nvec, nvec, nvec, nvec ); // U_Z(x) (1-gamma_Z) phi(x+hat{Z}) - projected
    pbn_su3_Z_PRECISION( buf3, leta+12*j*nvec_eta, nvec, nvec, nvec_eta ); // eta(x+hat{Z}) -= U_Z(x)^dagger(x) (1+gamma_Z) phi(x) + lift back
    pbp_su3_Z_PRECISION( buf4, leta+12*k*nvec_eta, nvec, nvec, nvec_eta ); // eta(x) -= U_Z(x) (1-gamma_Z) phi(x+hat{Z}) + lift back
  }
  ind = index[Y]; // Y direction
  for ( i=0; i<length[Y]; i++ ) {
    k = ind[i]; j = neighbor[4*k+Y]; D_pt = D + 36*k + 9*Y;
    prn_Y_PRECISION( buf1, lphi+12*k*nvec_phi, nvec, nvec, nvec_phi ); // (1+gamma_Y) phi(x) + projection
    prp_Y_PRECISION( buf2, lphi+12*j*nvec_phi, nvec, nvec, nvec_phi ); // (1-gamma_Y) phi(x+hat{Y}) + projection
    mvmh_PRECISION( buf3, D_pt, buf1, nvec, nvec, nvec );     // U_Y^dagger(x) (1+gamma_Y) phi(x) - projected
    mvmh_PRECISION( buf3+3*nvec, D_pt, buf1+3*nvec, nvec, nvec, nvec ); // U_Y^dagger(x) (1+gamma_Y) phi(x) - projected
    mvm_PRECISION( buf4, D_pt, buf2, nvec, nvec, nvec );     // U_Y(x) (1-gamma_Y) phi(x+hat{Y}) - projected
    mvm_PRECISION( buf4+3*nvec, D_pt, buf2+3*nvec, nvec, nvec, nvec ); // U_Y(x) (1-gamma_Y) phi(x+hat{Y}) - projected
    pbn_su3_Y_PRECISION( buf3, leta+12*j*nvec_eta, nvec, nvec, nvec_eta ); // eta(x+hat{Y}) -= U_Y(x)^dagger(x) (1+gamma_Y) phi(x) + lift back
    pbp_su3_Y_PRECISION( buf4, leta+12*k*nvec_eta, nvec, nvec, nvec_eta ); // eta(x) -= U_Y(x) (1-gamma_Y) phi(x+hat{Y}) + lift back
  }
  ind = index[X]; // X direction
  for ( i=0; i<length[X]; i++ ) {
    k = ind[i]; j = neighbor[4*k+X]; D_pt = D + 36*k + 9*X;
    prn_X_PRECISION( buf1, lphi+12*k*nvec_phi, nvec, nvec, nvec_phi ); // (1+gamma_X) phi(x) + projection
    prp_X_PRECISION( buf2, lphi+12*j*nvec_phi, nvec, nvec, nvec_phi ); // (1-gamma_X) phi(x+hat{X}) + projection
    mvmh_PRECISION( buf3, D_pt, buf1, nvec, nvec, nvec );     // U_X^dagger(x) (1+gamma_X) phi(x) - projected
    mvmh_PRECISION( buf3+3*nvec, D_pt, buf1+3*nvec, nvec, nvec, nvec ); // U_X^dagger(x) (1+gamma_X) phi(x) - projected
    mvm_PRECISION( buf4, D_pt, buf2, nvec, nvec, nvec );     // U_mu(x) (1-gamma_X) phi(x+hat{X}) - projected
    mvm_PRECISION( buf4+3*nvec, D_pt, buf2+3*nvec, nvec, nvec, nvec ); // U_mu(x) (1-gamma_X) phi(x+hat{X}) - projected
    pbn_su3_X_PRECISION( buf3, leta+12*j*nvec_eta, nvec, nvec, nvec_eta ); // eta(x+hat{X}) -= U_X(x)^dagger(x) (1+gamma_X) phi(x) + lift back
    pbp_su3_X_PRECISION( buf4, leta+12*k*nvec_eta, nvec, nvec, nvec_eta ); // eta(x) -= U_X(x) (1-gamma_X) phi(x+hat{X}) + lift back
  }
  END_UNTHREADED_FUNCTION(threading)
}

/********************  aggreagate operators for coarse op. setup  ************************************************/

// eta <-spin0and1_clover*phi
static void spin0and1_clover_PRECISION( vector_PRECISION *eta, vector_PRECISION *phi, config_PRECISION clover, level_struct *l ) {

#ifdef DEBUG
  if ( eta->num_vect < phi->num_vect_now )
    error0("spin0and1_clover_PRECISION: assumptions are not met\n");
#endif
  
  int j, jj, nvec = phi->num_vect_now, nvec_phi = phi->num_vect, nvec_eta=eta->num_vect;
  buffer_PRECISION eta_end = eta->vector_buffer + l->inner_vector_size*nvec_eta, leta = eta->vector_buffer, lphi = phi->vector_buffer;

  if ( g.csw == 0.0 ) {
    while ( leta < eta_end ) {
      FOR6( VECTOR_LOOP(j, nvec, jj, leta[j+jj] = lphi[j+jj]*(*clover);)
	    leta += nvec_eta; lphi += nvec_phi; clover++; )
      FOR6( VECTOR_LOOP(j, nvec_eta, jj, *leta = _COMPLEX_PRECISION_ZERO; leta++;))
      lphi += 6*nvec_phi; clover+=6;
    }
  } else {
    while ( leta < eta_end ) {
      spin0and1_site_clover_PRECISION( leta, lphi, clover, nvec, nvec_eta, nvec_phi );
      leta += 12*nvec_eta; lphi += 12*nvec_phi; clover+=42;
    }
  }
}

// eta <- spin2and3_clover*phi
static void spin2and3_clover_PRECISION( vector_PRECISION *eta, vector_PRECISION *phi, config_PRECISION clover, level_struct *l ) {

#ifdef DEBUG
  if ( eta->num_vect < phi->num_vect_now )
    error0("spin2and3_clover_PRECISION: assumptions are not met\n");
#endif
  
  int j, jj, nvec = phi->num_vect_now, nvec_phi = phi->num_vect, nvec_eta = eta->num_vect;
  buffer_PRECISION eta_end = eta->vector_buffer + l->inner_vector_size*nvec_eta, leta = eta->vector_buffer, lphi = phi->vector_buffer;
  if ( g.csw == 0.0 ) {
    while ( leta < eta_end ) {
      lphi += 6*nvec_phi; clover+=6;
      FOR6( VECTOR_LOOP(j, nvec_eta, jj, *leta = _COMPLEX_PRECISION_ZERO; leta++;))
      FOR6( VECTOR_LOOP(j, nvec, jj, leta[j+jj] = lphi[j+jj]*(*clover);)
	    leta += nvec_eta; lphi += nvec_phi; clover++;)
    }
  } else {
    while ( leta < eta_end ) {
      spin2and3_site_clover_PRECISION( leta, lphi, clover, nvec, nvec_eta, nvec_phi );
      leta += 12*nvec_eta; lphi += 12*nvec_phi; clover+=42;
    }
  }
}

// used in coarse_operator_PRECISION_setup: no_threading  
void diagonal_aggregate_PRECISION( vector_PRECISION *eta1, vector_PRECISION *eta2, vector_PRECISION *phi, config_PRECISION diag, level_struct *l ) {

  int j, jj, nvec = phi->num_vect_now, nvec_eta1 = eta1->num_vect, nvec_eta2 = eta2->num_vect, nvec_phi = phi->num_vect;
  buffer_PRECISION eta_end = eta1->vector_buffer + l->inner_vector_size*nvec_eta1;
  buffer_PRECISION eta1_pt = eta1->vector_buffer, eta2_pt = eta2->vector_buffer, phi_pt = phi->vector_buffer;

#ifdef DEBUG
  if ( nvec_eta1 != nvec_eta2 )
    error0("diagonal_aggregate_PRECISION: assumptions are not met\n");
#endif
  
  while ( eta1_pt < eta_end ) {
    FOR6( VECTOR_LOOP(j, nvec, jj, eta1_pt[j+jj] = phi_pt[j+jj]*(*diag); eta2_pt[j+jj] = _COMPLEX_PRECISION_ZERO;)
	  eta1_pt += nvec_eta1; eta2_pt += nvec_eta2; phi_pt += nvec_phi; diag++; )
    FOR6( VECTOR_LOOP(j, nvec, jj, eta2_pt[j+jj] = phi_pt[j+jj]*(*diag); eta1_pt[j+jj] = _COMPLEX_PRECISION_ZERO;)
	  eta1_pt += nvec_eta1; eta2_pt += nvec_eta2; phi_pt += nvec_phi; diag++; )
  }
}

// eta <- d_plus_clover_aggregate*phi
// used in coarse_operator_PRECISION_setup: no_threading
void d_plus_clover_aggregate_PRECISION( vector_PRECISION *eta1, vector_PRECISION *eta2, vector_PRECISION *phi, schwarz_PRECISION_struct *s, level_struct *l ) {
  
  int i, length, index1, index2, *index_dir, *neighbor = s->op.neighbor_table;
  int nvec = phi->num_vect_now, nvec_eta1 = eta1->num_vect, nvec_eta2 = eta2->num_vect, nvec_phi = phi->num_vect;
  buffer_PRECISION eta1_pt, eta2_pt, phi_pt;
  complex_PRECISION buffer1[12*nvec], buffer2[12*nvec];
  config_PRECISION D_pt, D = s->op.D;

#ifdef DEBUG
  if ( nvec_eta1 != nvec_eta2 || nvec_eta1 < nvec_phi )
    error0("d_plus_clover_aggregate_PRECISION: assumptions are not met\n");
#endif
  
  // add clover term/shift
  spin0and1_clover_PRECISION( eta1, phi, s->op.clover, l );
  spin2and3_clover_PRECISION( eta2, phi, s->op.clover, l );

  // T dir
  length = l->is_PRECISION.agg_length[T]; index_dir = l->is_PRECISION.agg_index[T];
  for ( i=0; i<length; i++ ) {
    index1 = index_dir[i]; index2 = neighbor[4*index1 + T];
    phi_pt = phi->vector_buffer + 12*index2*nvec_phi; D_pt = D + 36*index1+9*T;
    mvm_PRECISION( buffer1, D_pt, phi_pt, nvec, nvec, nvec_phi );
    mvm_PRECISION( buffer1+3*nvec, D_pt, phi_pt+3*nvec_phi, nvec, nvec, nvec_phi );
    mvm_PRECISION( buffer1+6*nvec, D_pt, phi_pt+6*nvec_phi, nvec, nvec, nvec_phi );
    mvm_PRECISION( buffer1+9*nvec, D_pt, phi_pt+9*nvec_phi, nvec, nvec, nvec_phi );
    phi_pt = phi->vector_buffer + 12*index1*nvec_phi;
    mvmh_PRECISION( buffer2, D_pt, phi_pt, nvec, nvec, nvec_phi );
    mvmh_PRECISION( buffer2+3*nvec, D_pt, phi_pt+3*nvec_phi, nvec, nvec, nvec_phi );
    mvmh_PRECISION( buffer2+6*nvec, D_pt, phi_pt+6*nvec_phi, nvec, nvec, nvec_phi );
    mvmh_PRECISION( buffer2+9*nvec, D_pt, phi_pt+9*nvec_phi, nvec, nvec, nvec_phi );
    eta1_pt = eta1->vector_buffer + 12*index1*nvec_eta1; eta2_pt = eta2->vector_buffer + 12*index1*nvec_eta2;
    twospin_p_T_PRECISION( eta1_pt, eta2_pt, buffer1, nvec, nvec_eta1, nvec_eta2, nvec );
    eta1_pt = eta1->vector_buffer + 12*index2*nvec_eta1; eta2_pt = eta2->vector_buffer + 12*index2*nvec_eta2;
    twospin_n_T_PRECISION( eta1_pt, eta2_pt, buffer2, nvec, nvec_eta1, nvec_eta2, nvec );
  }
  // Z dir
  length = l->is_PRECISION.agg_length[Z]; index_dir = l->is_PRECISION.agg_index[Z];
  for ( i=0; i<length; i++ ) {
    index1 = index_dir[i]; index2 = neighbor[4*index1 + Z];
    phi_pt = phi->vector_buffer + 12*index2*nvec_phi; D_pt = D + 36*index1+9*Z;
    mvm_PRECISION( buffer1, D_pt, phi_pt, nvec, nvec, nvec_phi );
    mvm_PRECISION( buffer1+3*nvec, D_pt, phi_pt+3*nvec_phi, nvec, nvec, nvec_phi );
    mvm_PRECISION( buffer1+6*nvec, D_pt, phi_pt+6*nvec_phi, nvec, nvec, nvec_phi );
    mvm_PRECISION( buffer1+9*nvec, D_pt, phi_pt+9*nvec_phi, nvec, nvec, nvec_phi );
    phi_pt = phi->vector_buffer + 12*index1*nvec_phi;
    mvmh_PRECISION( buffer2, D_pt, phi_pt, nvec, nvec, nvec_phi );
    mvmh_PRECISION( buffer2+3*nvec, D_pt, phi_pt+3*nvec_phi, nvec, nvec, nvec_phi );
    mvmh_PRECISION( buffer2+6*nvec, D_pt, phi_pt+6*nvec_phi, nvec, nvec, nvec_phi );
    mvmh_PRECISION( buffer2+9*nvec, D_pt, phi_pt+9*nvec_phi, nvec, nvec, nvec_phi );
    eta1_pt = eta1->vector_buffer + 12*index1*nvec_eta1; eta2_pt = eta2->vector_buffer + 12*index1*nvec_eta2;
    twospin_p_Z_PRECISION( eta1_pt, eta2_pt, buffer1, nvec, nvec_eta1, nvec_eta2, nvec );
    eta1_pt = eta1->vector_buffer + 12*index2*nvec_eta1; eta2_pt = eta2->vector_buffer + 12*index2*nvec_eta2;
    twospin_n_Z_PRECISION( eta1_pt, eta2_pt, buffer2, nvec, nvec_eta1, nvec_eta2, nvec );
  }
  // Y dir
  length = l->is_PRECISION.agg_length[Y]; index_dir = l->is_PRECISION.agg_index[Y];
  for ( i=0; i<length; i++ ) {
    index1 = index_dir[i]; index2 = neighbor[4*index1 + Y];
    phi_pt = phi->vector_buffer + 12*index2*nvec_phi; D_pt = D + 36*index1+9*Y;
    mvm_PRECISION( buffer1, D_pt, phi_pt, nvec, nvec, nvec_phi );
    mvm_PRECISION( buffer1+3*nvec, D_pt, phi_pt+3*nvec_phi, nvec, nvec, nvec_phi );
    mvm_PRECISION( buffer1+6*nvec, D_pt, phi_pt+6*nvec_phi, nvec, nvec, nvec_phi );
    mvm_PRECISION( buffer1+9*nvec, D_pt, phi_pt+9*nvec_phi, nvec, nvec, nvec_phi );
    phi_pt = phi->vector_buffer + 12*index1*nvec_phi;
    mvmh_PRECISION( buffer2, D_pt, phi_pt, nvec, nvec, nvec_phi );
    mvmh_PRECISION( buffer2+3*nvec, D_pt, phi_pt+3*nvec_phi, nvec, nvec, nvec_phi );
    mvmh_PRECISION( buffer2+6*nvec, D_pt, phi_pt+6*nvec_phi, nvec, nvec, nvec_phi );
    mvmh_PRECISION( buffer2+9*nvec, D_pt, phi_pt+9*nvec_phi, nvec, nvec, nvec_phi );
    eta1_pt = eta1->vector_buffer + 12*index1*nvec_eta1; eta2_pt = eta2->vector_buffer + 12*index1*nvec_eta2;
    twospin_p_Y_PRECISION( eta1_pt, eta2_pt, buffer1, nvec, nvec_eta1, nvec_eta2, nvec );
    eta1_pt = eta1->vector_buffer + 12*index2*nvec_eta1; eta2_pt = eta2->vector_buffer + 12*index2*nvec_eta2;
    twospin_n_Y_PRECISION( eta1_pt, eta2_pt, buffer2, nvec, nvec_eta1, nvec_eta2 , nvec );
  }
  // X dir
  length = l->is_PRECISION.agg_length[X]; index_dir = l->is_PRECISION.agg_index[X];
  for ( i=0; i<length; i++ ) {
    index1 = index_dir[i]; index2 = neighbor[4*index1 + X];
    phi_pt = phi->vector_buffer + 12*index2*nvec_phi; D_pt = D + 36*index1+9*X;
    mvm_PRECISION( buffer1, D_pt, phi_pt, nvec, nvec, nvec_phi );
    mvm_PRECISION( buffer1+3*nvec, D_pt, phi_pt+3*nvec_phi, nvec, nvec, nvec_phi );
    mvm_PRECISION( buffer1+6*nvec, D_pt, phi_pt+6*nvec_phi, nvec, nvec, nvec_phi );
    mvm_PRECISION( buffer1+9*nvec, D_pt, phi_pt+9*nvec_phi, nvec, nvec, nvec_phi );
    phi_pt = phi->vector_buffer + 12*index1*nvec_phi;
    mvmh_PRECISION( buffer2, D_pt, phi_pt, nvec, nvec, nvec_phi );
    mvmh_PRECISION( buffer2+3*nvec, D_pt, phi_pt+3*nvec_phi, nvec, nvec, nvec_phi );
    mvmh_PRECISION( buffer2+6*nvec, D_pt, phi_pt+6*nvec_phi, nvec, nvec, nvec_phi );
    mvmh_PRECISION( buffer2+9*nvec, D_pt, phi_pt+9*nvec_phi, nvec, nvec, nvec_phi );
    eta1_pt = eta1->vector_buffer + 12*index1*nvec_eta1; eta2_pt = eta2->vector_buffer + 12*index1*nvec_eta2;
    twospin_p_X_PRECISION( eta1_pt, eta2_pt, buffer1, nvec, nvec_eta1, nvec_eta2 , nvec );
    eta1_pt = eta1->vector_buffer + 12*index2*nvec_eta1; eta2_pt = eta2->vector_buffer + 12*index2*nvec_eta2;
    twospin_n_X_PRECISION( eta1_pt, eta2_pt, buffer2, nvec, nvec_eta1, nvec_eta2 , nvec );
  }
}

// used in coarse_operator_PRECISION_setup: no_threading  
void d_neighbor_aggregate_PRECISION( vector_PRECISION *eta1, vector_PRECISION *eta2, vector_PRECISION *phi, const int mu, schwarz_PRECISION_struct *s, level_struct *l ) {
  
  int i, length, index1, index2, *index_dir, *neighbor;
  int nvec = phi->num_vect_now, nvec_eta1 = eta1->num_vect, nvec_eta2 = eta2->num_vect, nvec_phi = phi->num_vect;
  buffer_PRECISION eta1_pt, eta2_pt, phi_pt;
  complex_PRECISION buffer1[12*nvec];
  config_PRECISION D_pt, D = s->op.D;

#ifdef DEBUG
  if ( nvec_eta1 != nvec_eta2 || nvec_eta1 < nvec_phi )
    error0("d_neighbor_aggregate_PRECISION: assumptions are not met\n");
#endif
  
  length    = l->is_PRECISION.agg_boundary_length[mu];
  index_dir = l->is_PRECISION.agg_boundary_index[mu];
  neighbor  = l->is_PRECISION.agg_boundary_neighbor[mu];
  
  // requires the positive boundaries of phi to be communicated befor
  if ( mu == T ) {
    // T dir
    for ( i=0; i<length; i++ ) {
      index1 = index_dir[i]; index2 = neighbor[i];
      phi_pt = phi->vector_buffer + 12*index2*nvec_phi; D_pt = D + 36*index1 + 9*T;
      mvm_PRECISION( buffer1, D_pt, phi_pt, nvec, nvec, nvec_phi );
      mvm_PRECISION( buffer1+3*nvec, D_pt, phi_pt+3*nvec_phi, nvec, nvec, nvec_phi );
      mvm_PRECISION( buffer1+6*nvec, D_pt, phi_pt+6*nvec_phi, nvec, nvec, nvec_phi );
      mvm_PRECISION( buffer1+9*nvec, D_pt, phi_pt+9*nvec_phi, nvec, nvec, nvec_phi );
      eta1_pt = eta1->vector_buffer + 12*index1*nvec_eta1; eta2_pt = eta2->vector_buffer + 12*index1*nvec_eta2;   
      twospin2_p_T_PRECISION( eta1_pt, eta2_pt, buffer1, nvec, nvec_eta1, nvec_eta2, nvec );
    }
  } else if ( mu == Z ) {
    // Z dir
    for ( i=0; i<length; i++ ) {
      index1 = index_dir[i]; index2 = neighbor[i];
      phi_pt = phi->vector_buffer + 12*index2*nvec_phi; D_pt = D + 36*index1 + 9*Z;
      mvm_PRECISION( buffer1, D_pt, phi_pt, nvec, nvec, nvec_phi );
      mvm_PRECISION( buffer1+3*nvec, D_pt, phi_pt+3*nvec_phi, nvec, nvec, nvec_phi );
      mvm_PRECISION( buffer1+6*nvec, D_pt, phi_pt+6*nvec_phi, nvec, nvec, nvec_phi );
      mvm_PRECISION( buffer1+9*nvec, D_pt, phi_pt+9*nvec_phi, nvec, nvec, nvec_phi );
      eta1_pt = eta1->vector_buffer + 12*index1*nvec_eta1; eta2_pt = eta2->vector_buffer + 12*index1*nvec_eta2;
      twospin2_p_Z_PRECISION( eta1_pt, eta2_pt, buffer1, nvec, nvec_eta1, nvec_eta2, nvec );
    }
  } else if ( mu == Y ) {
    // Y dir
    for ( i=0; i<length; i++ ) {
      index1 = index_dir[i]; index2 = neighbor[i];
      phi_pt = phi->vector_buffer + 12*index2*nvec_phi; D_pt = D + 36*index1 + 9*Y;
      mvm_PRECISION( buffer1, D_pt, phi_pt, nvec, nvec, nvec_phi );
      mvm_PRECISION( buffer1+3*nvec, D_pt, phi_pt+3*nvec_phi, nvec, nvec, nvec_phi );
      mvm_PRECISION( buffer1+6*nvec, D_pt, phi_pt+6*nvec_phi, nvec, nvec, nvec_phi );
      mvm_PRECISION( buffer1+9*nvec, D_pt, phi_pt+9*nvec_phi, nvec, nvec, nvec_phi );
      eta1_pt = eta1->vector_buffer + 12*index1*nvec_eta1; eta2_pt = eta2->vector_buffer + 12*index1*nvec_eta2;
      twospin2_p_Y_PRECISION( eta1_pt, eta2_pt, buffer1, nvec, nvec_eta1, nvec_eta2, nvec );
    }
  } else if ( mu == X ) {
    // X dir
    for ( i=0; i<length; i++ ) {
      index1 = index_dir[i]; index2 = neighbor[i];
      phi_pt = phi->vector_buffer + 12*index2*nvec_phi; D_pt = D + 36*index1 + 9*X;
      mvm_PRECISION( buffer1, D_pt, phi_pt, nvec, nvec, nvec_phi );
      mvm_PRECISION( buffer1+3*nvec, D_pt, phi_pt+3*nvec_phi, nvec, nvec, nvec_phi );
      mvm_PRECISION( buffer1+6*nvec, D_pt, phi_pt+6*nvec_phi, nvec, nvec, nvec_phi );
      mvm_PRECISION( buffer1+9*nvec, D_pt, phi_pt+9*nvec_phi, nvec, nvec, nvec_phi );
      eta1_pt = eta1->vector_buffer + 12*index1*nvec_eta1; eta2_pt = eta2->vector_buffer + 12*index1*nvec_eta2;
      twospin2_p_X_PRECISION( eta1_pt, eta2_pt, buffer1, nvec, nvec_eta1, nvec_eta2, nvec );
    }
  }
}
 

/*************************** Setup and Update Functions ************************************/

// used only in DDalphaAMG_interface.c: so not changed!!!!!!!!!!!!
void operator_updates_PRECISION( level_struct *l, struct Thread *threading ) {

  if ( l->level > 0 ) {
    if ( !l->idle ) {
      START_LOCKED_MASTER(threading)
      coarse_operator_PRECISION_setup( &(l->is_PRECISION.interpolation_vec), l );
      conf_PRECISION_gather( &(l->next_level->s_PRECISION.op), &(l->next_level->op_PRECISION), l->next_level );
      END_LOCKED_MASTER(threading)
      if ( !l->next_level->idle && l->next_level->level > 0 ) {
        START_LOCKED_MASTER(threading)
        schwarz_PRECISION_boundary_update( &(l->next_level->s_PRECISION), l->next_level );
        END_LOCKED_MASTER(threading)
        if ( g.method >= 4 && g.odd_even ) 
          coarse_oddeven_setup_PRECISION( &(l->next_level->s_PRECISION.op), _REORDER, l->next_level, threading );
      }
      if ( !l->next_level->idle && l->next_level->level == 0 && g.odd_even )
        coarse_oddeven_setup_PRECISION( &(l->next_level->s_PRECISION.op), _NO_REORDERING, l->next_level, threading );
      operator_updates_PRECISION( l->next_level, threading );
    }
  }  
}


void m0_update_PRECISION( PRECISION m0, operator_PRECISION_struct *op, level_struct *l, struct Thread *threading ) {
  
  // no hyperthreading in this function
  if(threading->thread != 0)
    return;

  config_PRECISION clover = op->clover;
  
  if ( clover != NULL && op->m0 != m0 ) {
    int i, j;
    complex_PRECISION m0_diff = m0 - op->m0;

    START_MASTER(threading)
    op->m0 = m0;
    END_MASTER(threading)

    if( m0_diff != 0 ) {
      if ( l->depth == 0 ) {
        int start = threading->start_site[l->depth];
        int n     = threading->n_site[l->depth];
        clover += start*(g.csw?42:12);
        for ( i=0; i<n; i++ ) {
          for ( j=0; j<12; j++ ) {
            clover[j] += m0_diff;
          }
          // clover term diag also stored as complex, so size is 2*15+2*6 = 42
          clover += (g.csw?42:12);
        }
      } else {
        int start = threading->start_site[l->depth];
        int n     = threading->n_site[l->depth];
        int k = l->num_parent_eig_vect;
        int sc_size = (l->num_parent_eig_vect)*(l->num_parent_eig_vect*2+1);
        clover += start*sc_size;
        for ( i=0; i<n; i++ ) {
          for ( j=0; j<k; j++ ) {
            if ( j>0 ) clover += j+1; 
            *clover += m0_diff;
          }
          clover ++;
          for ( j=0; j<k; j++ ) {
            if ( j>0 ) clover += j+1;
            *clover += m0_diff;
          }
          clover += 1 + SQUARE(k);
        }
      }
    }
  }
}

void tm_term_PRECISION_setup( double mu, double *even, double odd, double factor, operator_PRECISION_struct *op,
                              level_struct *l, struct Thread *threading ) {
  /* Description:
   *   - Setup meta-data for tm term (and tm_term itself if HAVE_MULT_TM)
   *   - For now, we assume only multi-shifts on even sites but a single shift on odd sites
   *   - tm term is single shifted in setup phase and multi-shifted in the solver phase
   *   - If HAVE_MULT_TM is switched on, 
   *       - tm term is computed separately for each rhs and put into tm_term
   *       - clover_oo_inv is also computed separately for each rhs for g.odd_even
   *   - If not,
   *       - tm term is computed on the fly
   *       - a single clover_oo_inv at the coarsest level is computed using the smallest even shift
   *       - only a single clover_oo_inv is required at the top if g.method<4 as all rhs has the same shift on odd sites
   *       - when g.method>=4, usage of the smallest even shift to construct clover_oo_inv for all rhs might lead to a loss of precision and slowing down the solver
   * Note:
   *  clover_oo_inv in oe_op_PRECISION is used if g.method<4
   */
  
#ifdef HAVE_TM
  if(threading->thread != 0)
    return;
  
  if ( op->diff_mu_eo != NULL) {
    //--- Set meta fields related to tm term
    int nrt = (g.in_setup && g.method>0)?num_loop:g.num_rhs_vect;
    START_LOCKED_MASTER(threading)
    op->factor         = factor;
    op->mu             = factor*mu;
    op->mu_odd_shift   = factor*odd;
    op->odd_shifted_mu = factor*(mu + odd); //for on the fly
    op->even_shift_avg = 0; //for on the fly
    op->is_even_shifted_mu_nonzero = 0;
    for ( int i=0; i<nrt; i++ ) {
      op->mu_even_shift[i]            = (g.in_setup && g.method>0)?factor*even[0]:factor*even[i];
      op->even_shift_avg             += g.mu_factor[l->depth]*g.mu_even_shift[i]/g.num_rhs_vect;//need to check!!!!
      op->diff_mu_eo[i]               = op->mu_even_shift[i] - factor*odd;//for on the fly 
      op->is_even_shifted_mu_nonzero += (factor*mu + op->mu_even_shift[i]!=0.0)?1:0;
    }
    END_LOCKED_MASTER(threading)

#ifdef HAVE_MULT_TM
    //--- Compute tm_term
    /*      g.num_rhs_vect many rhs are grouped into chunks of the size num_loop.
     *      Each of these chunks is processed together simultaneously.
     *      To adapt to this, tm_term is ordered num_loop vectors in vector-fast index one after another.
     */
      
    int i, j;
    int start, end;
    compute_core_start_end_custom(0, l->num_inner_lattice_sites, &start, &end, l, threading, 1);
    int n = end-start;
    config_PRECISION tm_term = op->tm_term;

    int nr, jj, jjj;
    config_PRECISION  odd_proj   = op->odd_proj;
    complex_PRECISION shift      = (complex_PRECISION) I*factor*mu;
    complex_PRECISION even_shift[nrt];
    complex_PRECISION odd_shift  = (complex_PRECISION) I*factor*odd;
    complex_PRECISION odd_factor[nrt];

    for( nr=0; nr<nrt; nr++ ) {
      even_shift[nr] = (complex_PRECISION) (I*op->mu_even_shift[nr]);
      odd_factor[nr] = (complex_PRECISION) (odd_shift - even_shift[nr]);
    }

    if( l->depth == 0 ) {
      int nv = l->num_inner_lattice_sites*12;
      complex_PRECISION tm_shift[nrt];

      tm_term  += start*12*num_loop;
      odd_proj += start*12;
      
      for ( i=0; i<n; i++ ) {
	// compute a shifted mu for each rhs
	for( nr=0; nr<nrt; nr++ ) {
	  if( cimag(even_shift[nr]) == 0. && cimag(odd_shift) == 0. )
	    tm_shift[nr] = shift;
	  else
	    tm_shift[nr] = shift + even_shift[nr] + odd_proj[0]*odd_factor[nr]; // odd_proj[0]=0 at even sites and 1 otherwise
	}
	// fill up tm_term for each internal d.o.f. at the given site, jumping over chunks
	// ??? we might put loop over chunks outside
	FOR6( VECTOR_LOOP( jj, nrt, jjj, tm_term[jj*nv+jjj] = - tm_shift[jj+jjj];) tm_term+=num_loop; )//spin indicies orderd in T,Z,Y,X?
	FOR6( VECTOR_LOOP( jj, nrt, jjj, tm_term[jj*nv+jjj] =   tm_shift[jj+jjj];) tm_term+=num_loop; )
	odd_proj += 12;
      }
    } else {
      int k, m     = l->num_parent_eig_vect;
      int tm_size  = m*(m+1);// there are 2 m-by-m upper triangular block matrices for each aggregate
      int tm_vsize = l->num_inner_lattice_sites*tm_size;
      
      tm_term  += start*tm_size*num_loop;
      odd_proj += start*tm_size;
      
      if( factor*g.even_shifted == 0. && cimag(odd_shift) == 0. ) {
	// note: approximated eigenvectors consisting P is aggregate-wise orthonormalized.
	
        for ( i=0; i<n; i++ ) { // for each site
	  for ( j=0; j<m; j++ ) { // for each upper column
	    // off-diagonal entries
	    for ( k=0; k<j; k++ )
	      VECTOR_LOOP( jj, nrt, jjj, tm_term[k*num_loop+tm_vsize*jj+jjj] = _COMPLEX_PRECISION_ZERO;)// for off-diagonal entries
	    tm_term += j*num_loop;
	    // diagonal entry
	    VECTOR_LOOP( jj, nrt, jjj, tm_term[jj*tm_vsize+jjj]= -1.* shift;)
	    tm_term+=num_loop;
	  }
	  for ( j=0; j<m; j++ ) {
	    for ( k=0; k<j; k++ ) 
	      VECTOR_LOOP( jj, nrt, jjj, tm_term[k*num_loop+tm_vsize*jj+jjj] = _COMPLEX_PRECISION_ZERO;)
	    tm_term += j*num_loop;
	    // diagonal entry
	    VECTOR_LOOP( jj, nrt, jjj, tm_term[tm_vsize*jj+jjj]= shift;)
	    tm_term+=num_loop;
	  }
	}
      } else {
	/*
	  For off-diagonal entries, we need TM_jk = mu_even_shift*(v_j(e)*v_k(e))+mu_odd_shift*(v_j(o)*v_k(o))
	  Note: v_j(e)*v_k(e)+v_j(o)*v_k(o) = \delta_jk & odd_proj_jk = v_j(o)*v_k(o)
	  So: TM_jk(j!=k) = mu_even_shift*(v_j(e)*v_k(e)+v_j(o)*v_k(o))+odd_factor*odd_proj_jk
	  For diagonal entries, we use 1= v_j(e)*v_j(e)+v_j(o)*v_j(o)
	  Diagonal entries: mu + mu_even_shift*(v_j(e)*v_j(e))+mu_odd_shift*(v_j(o)*v_j(o)) 
	 */
        for ( i=0; i<n; i++ ) { // for each site
	  // for spin 0,1, i.e., T,Z
          for ( j=0; j<m; j++ ) { // for each upper column 
            for ( k=0; k<j; k++ ) // for off-diagonal entries
              VECTOR_LOOP( jj, nrt, jjj, tm_term[k*num_loop+jjj+jj*tm_vsize] = -1. * odd_factor [jj+jjj]* odd_proj[k]; )
	    // for diagonal entry
            tm_term += j*num_loop;
            odd_proj += j;
            VECTOR_LOOP( jj, nrt, jjj, tm_term[jjj+jj*tm_vsize] = -1.* ( shift + even_shift[jj+jjj] + odd_factor[jj+jjj] * (*odd_proj));)
            tm_term+=num_loop;
            odd_proj++;
          } 

	  // for spin 2,3, i.e., Y,X
          for ( j=0; j<m; j++ ) { // for each upper column 
            for ( k=0; k<j; k++ ) // for off-diagonal entries  
              VECTOR_LOOP( jj, nrt, jjj, tm_term[k*num_loop+jjj+jj*tm_vsize] = odd_factor [jj+jjj]* odd_proj[k]; )
	    // for diagonal entry
            tm_term += j*num_loop;
            odd_proj += j;
            VECTOR_LOOP( jj, nrt, jjj, tm_term[jjj+jj*tm_vsize] =  shift + even_shift[jj+jjj] + odd_factor[jj+jjj] * (*odd_proj);)
            tm_term+=num_loop;
            odd_proj++;
          } 
        }
      }
    }
#endif
  }
#endif
}

void epsbar_term_PRECISION_setup( PRECISION epsbar, PRECISION even, PRECISION odd, operator_PRECISION_struct *op,
                                  level_struct *l, struct Thread *threading ) {

  // epsbar_term is off-diagonal in flavour index and (block) diagonal in color-spin/num_eig_vect index, defined at each site
#ifdef HAVE_TM1p1
  if(threading->thread != 0)
    return;

  config_PRECISION eps_term = op->epsbar_term;
  if ( eps_term != NULL ) {
    config_PRECISION odd_proj = op->odd_proj;
    complex_PRECISION shift = -epsbar;
    complex_PRECISION even_shift = I*even;
    complex_PRECISION odd_shift = I*odd;

    START_MASTER(threading)
    op->epsbar = epsbar;
    op->epsbar_ig5_even_shift = even;
    op->epsbar_ig5_odd_shift = odd;
    END_MASTER(threading)

    int i, j;
    int start, end;
    compute_core_start_end_custom(0, l->num_inner_lattice_sites, &start, &end, l, threading, 1);
    int n = end-start;
          
    if ( l->depth == 0 ) {
      eps_term += start*12;
      odd_proj += start*12;

      if( cimag(even_shift) == 0. && cimag(odd_shift) == 0. )
        for ( i=0; i<n; i++ ) {
          FOR12( *eps_term = shift; eps_term++; );
        }
      else
        for ( i=0; i<n; i++ ) {
          complex_PRECISION ig5_shift = even_shift + (*odd_proj)*(odd_shift - even_shift);
          FOR6( *eps_term = shift-ig5_shift; eps_term++; );
          FOR6( *eps_term = shift+ig5_shift; eps_term++; );
          odd_proj += 12;
      }
    } else {
      int k, m  = l->num_parent_eig_vect;
      int eps_size = m*(m+1);
      
      eps_term += start*eps_size;
      odd_proj += start*eps_size;
      
      if( cimag(even_shift) == 0. && cimag(odd_shift) == 0. ) {
        for ( i=0; i<2*n; i++ ) { // no multiplication by gamma_5 so each block term is identical
          for ( j=0; j<m; j++ ) {
            for ( k=0; k<j; k++ )
              eps_term[k] = _COMPLEX_PRECISION_ZERO;
            eps_term += j;
            *eps_term = shift;
            eps_term++;
          }
        } 
      } else {
        complex_PRECISION odd_factor = odd_shift - even_shift;
         
        for ( i=0; i<n; i++ ) {
          for ( j=0; j<m; j++ ) {
            for ( k=0; k<j; k++ ) 
              eps_term[k] = -1.* odd_factor*odd_proj[k] ;
            eps_term += j;
            odd_proj += j;
            *eps_term = shift - (even_shift + odd_factor * (*odd_proj));
            eps_term++;
            odd_proj++;
          } 
          
          for ( j=0; j<m; j++ ) {
            for ( k=0; k<j; k++ ) 
              eps_term[k] = odd_factor*odd_proj[k] ;
            eps_term += j;
            odd_proj += j;
            *eps_term = shift + (even_shift + odd_factor * (*odd_proj));
            eps_term++;
            odd_proj++;
          } 
        }
      }
    }  
  }
#endif
}


/*************************** Helper Functions ************************************/

void apply_twisted_bc_to_vector_PRECISION( vector_PRECISION *eta, vector_PRECISION *phi, double *theta, level_struct *l) {
  int t, z, y, x, j, jj, nvec = phi->num_vect_now, nvec_eta = eta->num_vect, nvec_phi = phi->num_vect;
  int *gl=l->global_lattice, sl[4];
  double phase[4];
  buffer_PRECISION eta_pt = eta->vector_buffer, phi_pt = phi->vector_buffer;
  complex_double twisted_bc;
  for (int i=0; i<4; i++)
    sl[i] = l->local_lattice[i]*g.my_coords[i];

#ifdef DEBUG
  if ( nvec_eta < nvec )
    error0("apply_twisted_bc_to_vector_PRECISION: assumptions are not met\n");
#endif
  
  for (t=0; t<l->local_lattice[0]; t++) {
    phase[T] = theta[T]*((double)sl[T]+t)/(double)gl[T];
    for (z=0; z<l->local_lattice[1]; z++) {
      phase[Z] = phase[T] + theta[Z]*((double)sl[Z]+z)/(double)gl[Z];
      for (y=0; y<l->local_lattice[2]; y++) {
        phase[Y] = phase[Z] + theta[Y]*((double)sl[Y]+y)/(double)gl[Y];
        for (x=0; x<l->local_lattice[3]; x++) {
          phase[X] = phase[Y] + theta[X]*((double)sl[X]+x)/(double)gl[X];
          twisted_bc = exp(I*phase[X]);
#ifdef HAVE_TM1p1
          if( g.n_flavours == 2 ) {
	    FOR12( VECTOR_LOOP2(j, nvec, jj, eta_pt[j+jj] = phi_pt[j+jj]*twisted_bc;)
		   eta_pt += nvec_eta; phi_pt += nvec_phi;)
          } else
#endif
	    FOR12( VECTOR_LOOP(j, nvec, jj, eta_pt[j+jj] = phi_pt[j+jj]*twisted_bc;)
		   eta_pt += nvec_eta; phi_pt += nvec_phi;)
        }
      }
    }
  }
}


// eta <- gamma5*phi: used in DDalphaAMG_interface    
void gamma5_PRECISION( vector_PRECISION *eta, vector_PRECISION *phi, level_struct *l, struct Thread *threading ) {

#ifdef DEBUG
  ASSERT(l->depth == 0);
  if ( eta->num_vect < phi->num_vect_now )
    error0("gamma5_PRECISION: assumptions are not met\n");
#endif
  
  int j, jj, nvec = phi->num_vect_now, nvec_phi = phi->num_vect, nvec_eta = eta->num_vect;
  buffer_PRECISION eta_end = eta->vector_buffer + threading->end_index[l->depth]*nvec_eta;
  buffer_PRECISION leta = eta->vector_buffer + threading->start_index[l->depth]*nvec_eta, lphi = phi->vector_buffer + threading->start_index[l->depth]*nvec_phi;

  while ( leta < eta_end ) {
    FOR6( VECTOR_LOOP( j, nvec, jj, leta[j+jj] = -lphi[j+jj];)
	  leta += nvec_eta; lphi += nvec_phi; )
    FOR6( VECTOR_LOOP( j, nvec, jj, leta[j+jj] = lphi[j+jj];)
	  leta += nvec_eta; lphi += nvec_phi; )
  }
}

// used in DDalphaAMG_interface
// eta_odd <- gamma_5*phi_odd; eta_even <- 0
void gamma5_set_even_to_zero_PRECISION( vector_PRECISION *eta, vector_PRECISION *phi, level_struct *l, struct Thread *threading ) {

#ifdef DEBUG
  ASSERT(l->depth == 0);
  if ( eta->num_vect < phi->num_vect_now )
    error0("gamma5_set_even_to_zero_PRECISION: assumptions are not met\n");
#endif
  
  int i = threading->start_site[l->depth], j, jj, nvec = phi->num_vect_now, nvec_phi = phi->num_vect, nvec_eta = eta->num_vect;
  buffer_PRECISION eta_end = eta->vector_buffer + threading->end_index[l->depth]*nvec_eta;
  buffer_PRECISION leta = eta->vector_buffer + threading->start_index[l->depth]*nvec_eta, lphi = phi->vector_buffer + threading->start_index[l->depth]*nvec_phi;

  while ( leta < eta_end ) {
    if(g.odd_even_table[i]==_ODD){
      FOR6( VECTOR_LOOP( j, nvec, jj, leta[j+jj] = -lphi[j+jj];) leta += nvec_eta; lphi += nvec_phi;)
      FOR6( VECTOR_LOOP( j, nvec, jj, leta[j+jj] = lphi[j+jj];) leta += nvec_eta; lphi += nvec_phi;)
    }
    else if(g.odd_even_table[i]==_EVEN){
      FOR12( VECTOR_LOOP( j, nvec, jj, leta[j+jj] = _COMPLEX_PRECISION_ZERO;) leta += nvec_eta; lphi += nvec_phi; )
    }
    i++;
  }
}

// used in DDalphaAMG_interface
// eta_odd <- 0; eta_even <- gamma_5*phi_even 
void gamma5_set_odd_to_zero_PRECISION( vector_PRECISION *eta, vector_PRECISION *phi, level_struct *l, struct Thread *threading ) {

#ifdef DEBUG
  ASSERT(l->depth == 0);
  if ( eta->num_vect < phi->num_vect_now )
    error0("gamma5_set_odd_to_zero_PRECISION: assumptions are not met\n");
#endif

  int i = threading->start_site[l->depth], j, jj, nvec = phi->num_vect_now, nvec_phi = phi->num_vect, nvec_eta = eta->num_vect;
  buffer_PRECISION eta_end = eta->vector_buffer + threading->end_index[l->depth]*nvec_eta;
  buffer_PRECISION leta = eta->vector_buffer + threading->start_index[l->depth]*nvec_eta, lphi = phi->vector_buffer + threading->start_index[l->depth]*nvec_eta;

  while ( leta < eta_end ) {
    if(g.odd_even_table[i]==_EVEN){
      FOR6( VECTOR_LOOP( j, nvec, jj, leta[j+jj] = -lphi[j+jj];) leta += nvec_eta; lphi += nvec_phi;)
      FOR6( VECTOR_LOOP( j, nvec, jj, leta[j+jj] = lphi[j+jj];) leta += nvec_eta; lphi += nvec_phi;)
    }
    else if(g.odd_even_table[i]==_ODD){
      FOR12( VECTOR_LOOP( j, nvec, jj, leta[j+jj] = _COMPLEX_PRECISION_ZERO;) leta += nvec_eta; lphi += nvec_phi; )
    }
    i++;
  }
}

/****************** TM1p1 ******************************************************************/

void tau1_gamma5_PRECISION( vector_PRECISION *eta, vector_PRECISION *phi, level_struct *l, struct Thread *threading ) {

#if defined(DEBUG) && defined(HAVE_TM1p1)
  ASSERT(l->depth == 0);
  if ( g.n_flavours == 2 && (eta->num_vect != eta->num_vect_now || phi->num_vect != phi->num_vect_now || eta->num_vect_now != phi->num_vect_now) )
    error0("tau1_gamma5_PRECISION: assumptions are not met\n");
#endif
  
#ifdef HAVE_TM1p1
  if ( g.n_flavours == 2 ) {
    int j, jj, nvec = phi->num_vect_now/2, nvec_eta = eta->num_vect, nvec_phi = phi->num_vect;
    buffer_PRECISION eta_end = eta->vector_buffer + threading->end_index[l->depth]*nvec_eta;
    buffer_PRECISION leta = eta->vector_buffer + threading->start_index[l->depth]*nvec_eta, lphi = phi->vector_buffer + threading->start_index[l->depth]*nvec_phi;

    while ( leta < eta_end ) {
      FOR6( VECTOR_LOOP( j, nvec, jj, leta[j+jj] = -lphi[j+jj+nvec_phi/2];)
	    VECTOR_LOOP( j, nvec, jj, leta[j+jj+nvec_eta/2] = -lphi[j+jj];)
	    leta += nvec_eta; lphi += nvec_phi;)
      FOR6( VECTOR_LOOP( j, nvec, jj, leta[j+jj] = lphi[j+jj+nvec_phi/2];)
            VECTOR_LOOP( j, nvec, jj, leta[j+jj+nvec_eta/2] = lphi[j+jj];)
            leta += nvec_eta; lphi += nvec_phi;)
    }
  } else 
#endif
    {
      START_MASTER(threading)
      warning0("tau1_gamma5_PRECISION called with g.n_flavours != 2\n");
      END_MASTER(threading)
      gamma5_PRECISION( eta, phi, l, threading );
    }
}

// used in DDalphaAMG_interface    
void tau1_gamma5_set_even_to_zero_PRECISION( vector_PRECISION *eta, vector_PRECISION *phi, level_struct *l, struct Thread *threading ) {

#if defined(DEBUG) && defined(HAVE_TM1p1)
  ASSERT(l->depth == 0);
  if ( g.n_flavours == 2 && (eta->num_vect != eta->num_vect_now || phi->num_vect != phi->num_vect_now || eta->num_vect_now != phi->num_vect_now) )
    error0("tau1_gamma5_set_even_to_zero__PRECISION: assumptions are not met\n");
#endif
  
#ifdef HAVE_TM1p1
  if ( g.n_flavours == 2 ) {
    int i = threading->start_site[l->depth], j, jj;
    int nvec = phi->num_vect_now, nvec_eta = eta->num_vect, nvec_phi = phi->num_vect;
    buffer_PRECISION eta_end = eta->vector_buffer + threading->end_index[l->depth]*nvec_eta;
    buffer_PRECISION leta = eta->vector_buffer + threading->start_index[l->depth]*nvec_eta, lphi = phi->vector_buffer + threading->start_index[l->depth]*nvec_phi;

    while ( leta < eta_end ) {
      
      if(g.odd_even_table[i]==_ODD){
	FOR6( VECTOR_LOOP( j, nvec/2, jj, leta[j+jj] = -lphi[j+jj+nvec_phi/2];)
	      VECTOR_LOOP( j, nvec/2, jj, leta[j+jj+nvec_eta/2] = -lphi[j+jj];)
	      leta += nvec_eta; lphi += nvec_phi;)
	FOR6( VECTOR_LOOP( j, nvec/2, jj, leta[j+jj] = lphi[j+jj+nvec_phi/2];)
	      VECTOR_LOOP( j, nvec/2, jj, leta[j+jj+nvec_eta/2] = lphi[j+jj];)
	      leta += nvec_eta; lphi += nvec_phi;)
      } else if(g.odd_even_table[i]==_EVEN){
	FOR12( VECTOR_LOOP2( j, nvec, jj, leta[j+jj] = _COMPLEX_PRECISION_ZERO;)
	       lphi += nvec_phi; leta += nvec_eta; );
      }
      i++;
    }
  } else 
#endif
    {
      START_MASTER(threading)
      warning0("tau1_gamma5_set_even_to_zero_PRECISION called with g.n_flavours != 2\n");
      END_MASTER(threading)
      gamma5_set_even_to_zero_PRECISION( eta, phi, l, threading );
    }
}

// used in DDalphaAMG_interface    
void tau1_gamma5_set_odd_to_zero_PRECISION( vector_PRECISION *eta, vector_PRECISION *phi, level_struct *l, struct Thread *threading ) {

#if defined(DEBUG) && defined(HAVE_TM1p1)
  ASSERT(l->depth == 0);
  if ( g.n_flavours == 2 && (eta->num_vect != eta->num_vect_now || phi->num_vect != phi->num_vect_now || eta->num_vect_now != phi->num_vect_now) )
    error0("tau1_gamma5_set_odd_to_zero_PRECISION: assumptions are not met\n");
#endif
  
#ifdef HAVE_TM1p1
  if ( g.n_flavours == 2 ) {
    int i = threading->start_site[l->depth], j, jj;
    int nvec = phi->num_vect_now, nvec_eta = eta->num_vect, nvec_phi = phi->num_vect;
    buffer_PRECISION eta_end = eta->vector_buffer + threading->end_index[l->depth]*nvec_eta;
    buffer_PRECISION leta = eta->vector_buffer + threading->start_index[l->depth]*nvec_eta, lphi = phi->vector_buffer + threading->start_index[l->depth]*nvec_phi;

    int offset[2] = {1,-1};
    while ( leta < eta_end ) {
      if(g.odd_even_table[i]==_EVEN){
	FOR6( VECTOR_LOOP( j, nvec/2, jj, leta[j+jj] = -lphi[j+jj+nvec_phi/2];)
              VECTOR_LOOP( j, nvec/2, jj, leta[j+jj+nvec_eta/2] = -lphi[j+jj];)
              leta += nvec_eta; lphi += nvec_phi;)
        FOR6( VECTOR_LOOP( j, nvec/2, jj, leta[j+jj] = lphi[j+jj+nvec_phi/2];)
              VECTOR_LOOP( j, nvec/2, jj, leta[j+jj+nvec_eta/2] = lphi[j+jj];)
              leta += nvec_eta; lphi += nvec_phi;)
      } else if(g.odd_even_table[i]==_ODD){
        FOR12( VECTOR_LOOP2( j, nvec, jj, leta[j+jj] = _COMPLEX_PRECISION_ZERO;)
               lphi += nvec_phi; leta += nvec_eta; );
      }
      i++;
    }
  } else 
#endif
    {
      START_MASTER(threading)
      warning0("tau1_gamma5_set_odd_to_zero_PRECISION called with g.n_flavours != 2\n");
      END_MASTER(threading)
      gamma5_set_odd_to_zero_PRECISION( eta, phi, l, threading );
    }
}

/******************  TEST ROUTINES *********************************************************/
// new function not impremented!!!!!!!
void two_flavours_test_PRECISION( operator_PRECISION_struct *op, level_struct *l, struct Thread *threading ) {

#ifdef HAVE_TM1p1
  int n_vect = num_loop, j, jj;
  double diff1[n_vect], diff2[n_vect];
  vector_double vd[4], vdd[4];
  vector_PRECISION vpp[2];
  for(int i=0; i<4; i++){
    // storage for each flavor content
    vector_double_init( &vd[i] );                                                         
    vector_double_alloc( &vd[i], _INNER, num_loop, l, threading );
    // storage for a vector with two flavors
    vector_double_init( &vdd[i] );
    vector_double_alloc( &vdd[i], _INNER, 2*num_loop, l, threading );
    vdd[i].num_vect_now = 2*num_loop;
  }                                                                                       
                                                                                          
  for(int i=0; i<2; i++){                                                                 
    vector_PRECISION_init( &vpp[i] );                                                      
    vector_PRECISION_alloc( &vpp[i], _INNER, 2*num_loop, l, threading );
    vpp[i].num_vect_now = 2*num_loop;
  }

  ASSERT(g.n_flavours==2);

  /*******  n_flavour == 1 below ********/
  data_layout_n_flavours( 1, l, threading );
  
  ASSERT(g.n_flavours==1);
  
  START_LOCKED_MASTER(threading)
  vector_double_define_random( &vd[0], 0, l->inner_vector_size, l ); // vd[0]: up quark
  vector_double_define_random( &vd[1], 0, l->inner_vector_size, l ); // vd[1]: down quark
  apply_operator_double( &vd[2], &vd[0], &(g.p), l, no_threading ); // vd[2] = d_plus_clover_double(mu,vd[0]) 
#ifdef HAVE_MULT_TM
  buffer_double_real_scale( g.op_double.tm_term, g.op_double.tm_term, -1, 0, l->inner_vector_size*g.num_rhs_vect, l ); 
#endif
  apply_operator_double( &vd[3], &vd[1], &(g.p), l, no_threading ); // vd[3] = d_plus_clover_double(-mu,vd[1])
#ifdef HAVE_MULT_TM
  buffer_double_real_scale( g.op_double.tm_term, g.op_double.tm_term, -1, 0, l->inner_vector_size*g.num_rhs_vect, l ); 
#endif
  add_diagonal_double( &vd[2], &vd[1], g.op_double.epsbar_term, l->inner_vector_size ); // vd[2] += eps_term*vd[1]
  add_diagonal_double( &vd[3], &vd[0], g.op_double.epsbar_term, l->inner_vector_size ); // vd[3] += eps_term*vd[0]  

  two_flavours_to_serial_double( &vd[0], &vd[1], &vdd[0], l, no_threading ); // vdd[0] = vd[0] _ vd[1] : 2 flav orig
  two_flavours_to_serial_double( &vd[2], &vd[3], &vdd[1], l, no_threading ); // vdd[1] = vd[2] _ vd[3] : 2 flav after
  END_LOCKED_MASTER(threading)

  /*******  n_flavour == 2 below ********/
  data_layout_n_flavours( 2, l, threading );
  START_LOCKED_MASTER(threading)
  trans_PRECISION( &vpp[0], &vdd[0], op->translation_table, l, no_threading );      // vpp[0] = vdd[0] with reordering to schwarz_index
  apply_operator_PRECISION( &vpp[1], &vpp[0], &(l->p_PRECISION), l, no_threading ); // vpp[1] = d_plus_clover_double(vpp[0])
  trans_back_PRECISION( &vdd[2], &vpp[1], op->translation_table, l, no_threading ); // vdd[2] = vpp[1] with reordering back to local_lex_index
  vector_double_minus( &vdd[3], &vdd[2], &vdd[1], 0, l->inner_vector_size, l ); // vdd[3] = vdd[2] - vdd[1]
  //  vdd[3].num_vect_now = num_loop;vdd[2].num_vect_now  = num_loop;
  global_norm_double( diff1, &vdd[3], 0, l->inner_vector_size, l, no_threading );
  global_norm_double( diff2, &vdd[2], 0, l->inner_vector_size, l, no_threading );
  for(int i=0; i<n_vect; i++)
    test0_PRECISION("depth: %d, correctness of doublet Dirac operator PRECISION: %le\n", l->depth, diff1[i]/diff2[i] );
  END_LOCKED_MASTER(threading)

  if(threading->n_core > 1) {
    trans_PRECISION( &vpp[0], &vdd[0], op->translation_table, l, threading );
    apply_operator_PRECISION( &vpp[1], &vpp[0], &(l->p_PRECISION), l, threading );
    trans_back_PRECISION( &vdd[2], &vpp[1], op->translation_table, l, threading );
    
    SYNC_MASTER_TO_ALL(threading)
    SYNC_CORES(threading)

    START_LOCKED_MASTER(threading)
    vector_double_minus( &vdd[3], &vdd[2], &vdd[1], 0, l->inner_vector_size, l );
    global_norm_double( diff1, &vdd[3], 0, l->inner_vector_size, l, no_threading );
    global_norm_double( diff2, &vdd[2], 0, l->inner_vector_size, l, no_threading );
    for(int i=0; i<n_vect; i++)
      test0_PRECISION("depth: %d, correctness of doublet Dirac operator PRECISION with threading: %le\n", l->depth, diff1[i]/diff2[i] );
    END_LOCKED_MASTER(threading)
  }    
  
  for(int i=0; i<4; i++){
    vector_double_free( &vd[i], l, threading );
    vector_double_free( &vdd[i], l, threading );
  }

  for(int i=0; i<2; i++)
    vector_PRECISION_free( &vpp[i], l, threading );

  START_LOCKED_MASTER(threading)
  if ( g.method >= 4 && g.odd_even )
    oddeven_PRECISION_test( l );
  END_LOCKED_MASTER(threading) 
#endif
    
}
