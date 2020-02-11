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

// at depth = 0
// eta <- op->clover(+op->tm)*phi
void clover_PRECISION_new( vector_PRECISION *eta, vector_PRECISION *phi, operator_PRECISION_struct *op, int start, int end,
                       level_struct *l, struct Thread *threading ) {

  if ( eta->num_vect < phi->num_vect_now )
    error0("clover_PRECISION: assumptions are not met\n");

  int i, j, jj;
  int nv = l->num_lattice_site_var, nvec = phi->num_vect_now, nvec_phi = phi->num_vect, nvec_eta = eta->num_vect;
  
  buffer_PRECISION lphi = phi->vector_buffer+start*nvec_phi, leta = eta->vector_buffer+start*nvec_eta;
  buffer_PRECISION leta_end = eta->vector_buffer+end*nvec_eta;
#ifdef PROFILING
  START_MASTER(threading)
  PROF_PRECISION_START( _SC );
  END_MASTER(threading)
#endif

#ifdef HAVE_TM
  config_PRECISION tm_term = op->tm_term+(start/nv)*12;
#endif

  if ( g.csw == 0.0 ) {

    config_PRECISION clover = op->clover+(start/nv)*12;// why csw==0 => no clover????????
/*#ifdef HAVE_TM1p1
    if( g.n_flavours == 2 ) {
#ifdef HAVE_TM
      if ( g.mu + g.mu_odd_shift != 0.0 || g.mu + g.mu_even_shift != 0.0 ) 
        while ( leta < leta_end ) {
          FOR6( *leta = (*lphi)*((*clover)+(*tm_term)); leta++; lphi++; clover++; tm_term++; );
          clover -= 6;
          tm_term -= 6;
          FOR6( *leta = (*lphi)*((*clover)-(*tm_term)); leta++; lphi++; clover++; tm_term++; );
          FOR6( *leta = (*lphi)*((*clover)+(*tm_term)); leta++; lphi++; clover++; tm_term++; );
          clover -= 6;
          tm_term -= 6;
          FOR6( *leta = (*lphi)*((*clover)-(*tm_term)); leta++; lphi++; clover++; tm_term++; );
        }
      else
#endif
        while ( leta < leta_end ) {
          FOR6( *leta = (*lphi)*(*clover); leta++; lphi++; clover++; );
          clover -= 6;
          FOR12( *leta = (*lphi)*(*clover); leta++; lphi++; clover++; );
          clover -= 6;
          FOR6( *leta = (*lphi)*(*clover); leta++; lphi++; clover++; );
        }
    } else {
#endif*/
#ifdef HAVE_TM
      if ( g.mu + g.mu_odd_shift != 0.0 || g.mu + g.mu_even_shift != 0.0 ) {
        while ( leta < leta_end )
          for( i=0; i<12; i++ ) {
	    //VECTOR_LOOP(j, nvec, jj, eta->vector_buffer[(eta->start+i)*nvec_eta+j+jj] = 
            VECTOR_LOOP(j, nvec, jj, *leta = (*lphi)*((*clover)+(*tm_term));
                                     leta++;
                                     lphi++;) 
	    leta += nvec_eta-nvec; lphi += nvec_phi-nvec;
            clover++;
            tm_term++;
	  }
      }
#endif // else
      while ( leta < leta_end )
        for( i=0; i<12; i++ ){
          VECTOR_LOOP(j, nvec, jj, *leta = (*lphi)*(*clover);
                                   leta++;
                                   lphi++;)
	  leta += nvec_eta-nvec; lphi += nvec_phi-nvec;
          clover++;
	}
/*#ifdef HAVE_TM1p1
    }
#endif*/

  } else { // need clover contribution

    config_PRECISION clover = op->clover+(start/nv)*42;
/*#ifdef HAVE_TM1p1
    if( g.n_flavours == 2 ) {
#ifdef HAVE_TM
      if ( g.mu + g.mu_odd_shift != 0.0 || g.mu + g.mu_even_shift != 0.0 ) 
        while ( leta < leta_end ) {
          doublet_site_clover_PRECISION( leta, lphi, clover );
          clover+=42;
          FOR6( *leta += (*lphi)*(*tm_term); leta++; lphi++; tm_term++; );
          tm_term -= 6;
          FOR6( *leta -=(*lphi)*(*tm_term); leta++; lphi++; tm_term++; );
          FOR6( *leta += (*lphi)*(*tm_term); leta++; lphi++; tm_term++; );
          tm_term -= 6;
          FOR6( *leta -= (*lphi)*(*tm_term); leta++; lphi++; tm_term++; );
        }
      else
#endif
        while ( leta < leta_end ) {
          doublet_site_clover_PRECISION( leta, lphi, clover );
          leta+=24; lphi+=24;
          clover+=42;
        }
    } else {
#endif*/
#ifdef HAVE_TM
      if ( g.mu + g.mu_odd_shift != 0.0 || g.mu + g.mu_even_shift != 0.0 ) 
        while ( leta < leta_end ) {
          site_clover_PRECISION_new( leta, lphi, clover, nvec, nvec_eta, nvec_phi );
          for( i=0; i<12; i++ ){
            VECTOR_LOOP(j, nvec, jj, *leta += (*lphi)*(*tm_term);
			             leta++;
                                     lphi++;)
	    leta += nvec_eta-nvec; lphi+=nvec_phi-nvec;
            tm_term++;
          }
          clover+=42;
        }
#endif // else
      while ( leta < leta_end ) {
	site_clover_PRECISION_new( leta, lphi, clover, nvec, nvec_eta, nvec_phi );
	leta += 12*nvec_eta; lphi += 12*nvec_phi;
	clover+=42;
      }
/*#ifdef HAVE_TM1p1
    }
#endif  */
  }
/*
#ifdef HAVE_TM1p1
  config_PRECISION eps_term = op->epsbar_term+(start/nv)*12;  
  lphi = phi->vector_buffer+start+phi_shift, leta = eta->vector_buffer+start+eta_shift;
  if ( g.n_flavours == 2 &&
       ( g.epsbar != 0 || g.epsbar_ig5_odd_shift != 0 || g.epsbar_ig5_odd_shift != 0 ) )
    while ( leta < leta_end ) { 
      lphi += 6;
      FOR6( *leta += (*lphi)*(*eps_term); leta++; lphi++; eps_term++; )
      lphi -= 12;
      eps_term -= 6;
      FOR6( *leta += (*lphi)*(*eps_term); leta++; lphi++; eps_term++; )
      lphi += 6;
    }
#endif
*/
  
#ifdef PROFILING
  START_MASTER(threading)
  PROF_PRECISION_STOP( _SC, 1 );
  END_MASTER(threading)
#endif
    
}

// eta <-spin0and1_clover*phi
static void spin0and1_clover_PRECISION_new( vector_PRECISION *eta, vector_PRECISION *phi, config_PRECISION clover, level_struct *l ) {
  
  if ( eta->num_vect < phi->num_vect_now )
    error0("spin0and1_clover_PRECISION: assumptions are not met\n");

  int nvec = phi->num_vect_now, nvec_phi = phi->num_vect, nvec_eta=eta->num_vect, i, j, jj;
  buffer_PRECISION eta_end = eta->vector_buffer + l->inner_vector_size*nvec_eta, leta = eta->vector_buffer, lphi = phi->vector_buffer;

  if ( g.csw == 0.0 ) {
    while ( leta < eta_end ) {
      for( i=0; i<6; i++ ){
        VECTOR_LOOP(j, nvec, jj, *leta = (*lphi)*(*clover); leta++; lphi++;) 
	leta += nvec_eta-nvec; lphi += nvec_phi-nvec;
        clover++;
      }
      for( i=0; i<6; i++ ){
        VECTOR_LOOP(j, nvec_eta, jj, *leta = _COMPLEX_PRECISION_ZERO; leta++;)
      }
      lphi += 6*nvec_phi; clover+=6;
    }
  } else {
    while ( leta < eta_end ) {
      spin0and1_site_clover_PRECISION_new( leta, lphi, clover, nvec, nvec_eta, nvec_phi );
      leta += 12*nvec_eta; lphi += 12*nvec_phi; clover+=42;
    }
  }
}

// eta <- spin2and3_clove*phi
static void spin2and3_clover_PRECISION_new( vector_PRECISION *eta, vector_PRECISION *phi, config_PRECISION clover, level_struct *l ) {

  if ( eta->num_vect < phi->num_vect_now )
    error0("spin2and3_clover_PRECISION: assumptions are not met\n");

  int i, j, jj;
  int nvec = phi->num_vect_now, nvec_phi = phi->num_vect, nvec_eta = eta->num_vect;
  buffer_PRECISION eta_end = eta->vector_buffer + l->inner_vector_size*nvec_eta, leta = eta->vector_buffer, lphi = phi->vector_buffer;
  if ( g.csw == 0.0 ) {
    while ( leta < eta_end ) {
      lphi += 6*nvec_phi; clover+=6;
      for( i=0; i<6; i++ ){
        VECTOR_LOOP(j, nvec_eta, jj, *leta = _COMPLEX_PRECISION_ZERO; leta++;)
      }
      for( i=0; i<6; i++ ){
        VECTOR_LOOP(j, nvec, jj, *leta = (*lphi)*(*clover); leta++; lphi++;)//!!!!!!!!+= -> =
	leta += nvec_eta-nvec; lphi += nvec_phi-nvec;
        clover++;
      }
    }
  } else {
    while ( leta < eta_end ) {
      spin2and3_site_clover_PRECISION_new( leta, lphi, clover, nvec, nvec_eta, nvec_phi );
      leta += 12*nvec_eta; lphi += 12*nvec_phi; clover+=42;
    }
  }
}

// eta <- d_plus_clover*phi
void d_plus_clover_PRECISION_new( vector_PRECISION *eta, vector_PRECISION *phi, operator_PRECISION_struct *op, level_struct *l, struct Thread *threading ) {
  
  int n = l->num_inner_lattice_sites, *neighbor = op->neighbor_table, start, end, nv = l->num_lattice_site_var;
  int nvec = phi->num_vect_now, nvec_phi = phi->num_vect, nvec_eta = eta->num_vect, nvec_op = op->pr_num_vect;//!!!!!!
  int i, j, *nb_pt;
  buffer_PRECISION phi_pt, eta_pt, end_pt;
  config_PRECISION D_pt;

  if ( eta->num_vect < phi->num_vect_now )
    error0("d_plus_clover_PRECISION: assumptions are not met\n");

  compute_core_start_end(0, nv*n, &start, &end, l, threading );

  SYNC_MASTER_TO_ALL(threading)
  clover_PRECISION_new( eta, phi, op, start, end, l, threading );//eta=op*phi where clover term of op is applied
  START_MASTER(threading)
  PROF_PRECISION_START( _NC ); 
  END_MASTER(threading)
/*
#ifdef HAVE_TM1p1
  if( g.n_flavours == 2 ) {
    complex_PRECISION pbuf[12];  
    for ( i=start/2, phi_pt=phi->vector_buffer+start+phi_shift; i<end/2; i+=nv/2, phi_pt+=nv ) {
      dprp_T_PRECISION( op->prnT+i, phi_pt );
      dprp_Z_PRECISION( op->prnZ+i, phi_pt );
      dprp_Y_PRECISION( op->prnY+i, phi_pt );
      dprp_X_PRECISION( op->prnX+i, phi_pt );
    }
    // start communication in negative direction
    START_LOCKED_MASTER(threading)
    ghost_sendrecv_PRECISION( op->prnT, T, -1, &(op->c), _FULL_SYSTEM, l );
    ghost_sendrecv_PRECISION( op->prnZ, Z, -1, &(op->c), _FULL_SYSTEM, l );
    ghost_sendrecv_PRECISION( op->prnY, Y, -1, &(op->c), _FULL_SYSTEM, l );
    ghost_sendrecv_PRECISION( op->prnX, X, -1, &(op->c), _FULL_SYSTEM, l );
    END_LOCKED_MASTER(threading) 
    // project plus dir and multiply with U dagger
    for ( phi_pt=phi->vector_buffer+start+phi_shift,c end_pt=phi->vector_buffer+end+phi_shift, D_pt = op->D+((start/nv)*36), nb_pt=neighbor+((start/nv)*4); phi_pt<end_pt; phi_pt+=nv ) {
      // T dir
      j = nv/2*(*nb_pt); nb_pt++;
      dprn_T_PRECISION( pbuf, phi_pt );
      mvmh_PRECISION( op->prpT+j, D_pt, pbuf );
      mvmh_PRECISION( op->prpT+j+3, D_pt, pbuf+3 );
      mvmh_PRECISION( op->prpT+j+6, D_pt, pbuf+6 );
      mvmh_PRECISION( op->prpT+j+9, D_pt, pbuf+9 ); D_pt += 9;
      // Z dir
      j = nv/2*(*nb_pt); nb_pt++;
      dprn_Z_PRECISION( pbuf, phi_pt );
      mvmh_PRECISION( op->prpZ+j, D_pt, pbuf );
      mvmh_PRECISION( op->prpZ+j+3, D_pt, pbuf+3 );
      mvmh_PRECISION( op->prpZ+j+6, D_pt, pbuf+6 );
      mvmh_PRECISION( op->prpZ+j+9, D_pt, pbuf+9 ); D_pt += 9;
      // Y dir
      j = nv/2*(*nb_pt); nb_pt++;
      dprn_Y_PRECISION( pbuf, phi_pt );
      mvmh_PRECISION( op->prpY+j, D_pt, pbuf );
      mvmh_PRECISION( op->prpY+j+3, D_pt, pbuf+3 );
      mvmh_PRECISION( op->prpY+j+6, D_pt, pbuf+6 );
      mvmh_PRECISION( op->prpY+j+9, D_pt, pbuf+9 ); D_pt += 9;
      // X dir
      j = nv/2*(*nb_pt); nb_pt++;
      dprn_X_PRECISION( pbuf, phi_pt );
      mvmh_PRECISION( op->prpX+j, D_pt, pbuf );
      mvmh_PRECISION( op->prpX+j+3, D_pt, pbuf+3 );
      mvmh_PRECISION( op->prpX+j+6, D_pt, pbuf+6 );
      mvmh_PRECISION( op->prpX+j+9, D_pt, pbuf+9 ); D_pt += 9;
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
    for ( eta_pt=eta->vector_buffer+start+eta_shift, end_pt=eta->vector_buffer+end+eta_shift, D_pt = op->D+(start/nv)*36, nb_pt=neighbor+(start/nv)*4; eta_pt<end_pt; eta_pt+=nv ) {
      // T dir
      j = nv/2*(*nb_pt); nb_pt++;
      mvm_PRECISION( pbuf, D_pt, op->prnT+j );
      mvm_PRECISION( pbuf+3, D_pt, op->prnT+j+3 );
      mvm_PRECISION( pbuf+6, D_pt, op->prnT+j+6 );
      mvm_PRECISION( pbuf+9, D_pt, op->prnT+j+9 );
      dpbp_su3_T_PRECISION( pbuf, eta_pt ); D_pt += 9;
      // Z dir
      j = nv/2*(*nb_pt); nb_pt++;
      mvm_PRECISION( pbuf, D_pt, op->prnZ+j );
      mvm_PRECISION( pbuf+3, D_pt, op->prnZ+j+3 );
      mvm_PRECISION( pbuf+6, D_pt, op->prnZ+j+6 );
      mvm_PRECISION( pbuf+9, D_pt, op->prnZ+j+9 );
      dpbp_su3_Z_PRECISION( pbuf, eta_pt ); D_pt += 9;
      // Y dir
      j = nv/2*(*nb_pt); nb_pt++;
      mvm_PRECISION( pbuf, D_pt, op->prnY+j );
      mvm_PRECISION( pbuf+3, D_pt, op->prnY+j+3 );
      mvm_PRECISION( pbuf+6, D_pt, op->prnY+j+6 );
      mvm_PRECISION( pbuf+9, D_pt, op->prnY+j+9 );
      dpbp_su3_Y_PRECISION( pbuf, eta_pt ); D_pt += 9;
      // X dir
      j = nv/2*(*nb_pt); nb_pt++;
      mvm_PRECISION( pbuf, D_pt, op->prnX+j );
      mvm_PRECISION( pbuf+3, D_pt, op->prnX+j+3 );
      mvm_PRECISION( pbuf+6, D_pt, op->prnX+j+6 );
      mvm_PRECISION( pbuf+9, D_pt, op->prnX+j+9 );
      dpbp_su3_X_PRECISION( pbuf, eta_pt ); D_pt += 9;
    }

    // wait for communication in positive direction
    START_LOCKED_MASTER(threading)
    ghost_wait_PRECISION( op->prpT, T, +1, &(op->c), _FULL_SYSTEM, l );
    ghost_wait_PRECISION( op->prpZ, Z, +1, &(op->c), _FULL_SYSTEM, l );
    ghost_wait_PRECISION( op->prpY, Y, +1, &(op->c), _FULL_SYSTEM, l );
    ghost_wait_PRECISION( op->prpX, X, +1, &(op->c), _FULL_SYSTEM, l );
    END_LOCKED_MASTER(threading)
      
    // lift up plus dir
    for ( i=start/2, eta_pt=eta->vector_buffer+start+eta_shift; i<end/2; i+=12, eta_pt+=24 ) {
      dpbn_su3_T_PRECISION( op->prpT+i, eta_pt );
      dpbn_su3_Z_PRECISION( op->prpZ+i, eta_pt );
      dpbn_su3_Y_PRECISION( op->prpY+i, eta_pt );
      dpbn_su3_X_PRECISION( op->prpX+i, eta_pt );
    }
  } else {
#endif
*/
  complex_PRECISION pbuf[6*nvec];//pay attention to #vectors in op->prn* defined in operator_PRECISION_alloc_projection_buffers!!!!
  for ( i=start*nvec_op/2, phi_pt=phi->vector_buffer+start*nvec_phi; i<end*nvec_op/2; i+=6*nvec_op, phi_pt+=12*nvec_phi ) {
    prp_T_PRECISION_new( op->prnT+i, phi_pt, nvec, nvec_op, nvec_phi );
    prp_Z_PRECISION_new( op->prnZ+i, phi_pt, nvec, nvec_op, nvec_phi );
    prp_Y_PRECISION_new( op->prnY+i, phi_pt, nvec, nvec_op, nvec_phi );
    prp_X_PRECISION_new( op->prnX+i, phi_pt, nvec, nvec_op, nvec_phi );
  }
  // start communication in negative direction
  START_LOCKED_MASTER(threading)
  g.num_vect_pass1 = nvec; g.num_vect_pass2 = nvec_op;//temp fix!!!!!!
  ghost_sendrecv_PRECISION_new( op->prnT, T, -1, &(op->c), _FULL_SYSTEM, l );
  ghost_sendrecv_PRECISION_new( op->prnZ, Z, -1, &(op->c), _FULL_SYSTEM, l );
  ghost_sendrecv_PRECISION_new( op->prnY, Y, -1, &(op->c), _FULL_SYSTEM, l );
  ghost_sendrecv_PRECISION_new( op->prnX, X, -1, &(op->c), _FULL_SYSTEM, l );
  END_LOCKED_MASTER(threading) 
  
  // project plus dir and multiply with U dagger
  for ( phi_pt=phi->vector_buffer+start*nvec_phi, end_pt=phi->vector_buffer+end*nvec_phi, D_pt = op->D+(start*3), nb_pt=neighbor+((start/12)*4); phi_pt<end_pt; phi_pt+=12*nvec_phi ) {
    // T dir
    j = 6*(*nb_pt)*nvec_op; nb_pt++;
    prn_T_PRECISION_new( pbuf, phi_pt, nvec, nvec, nvec_phi );
    mvmh_PRECISION_new( op->prpT+j, D_pt, pbuf, nvec, nvec_op, nvec );
    mvmh_PRECISION_new( op->prpT+j+3*nvec_op, D_pt, pbuf+3*nvec, nvec, nvec_op, nvec ); D_pt += 9;
    // Z dir
    j = 6*(*nb_pt)*nvec_op; nb_pt++;
    prn_Z_PRECISION_new( pbuf, phi_pt, nvec, nvec, nvec_phi );
    mvmh_PRECISION_new( op->prpZ+j, D_pt, pbuf, nvec, nvec_op, nvec );
    mvmh_PRECISION_new( op->prpZ+j+3*nvec_op, D_pt, pbuf+3*nvec, nvec, nvec_op, nvec ); D_pt += 9;
    // Y dir
    j = 6*(*nb_pt)*nvec_op; nb_pt++;
    prn_Y_PRECISION_new( pbuf, phi_pt, nvec, nvec, nvec_phi );
    mvmh_PRECISION_new( op->prpY+j, D_pt, pbuf, nvec, nvec_op, nvec );
    mvmh_PRECISION_new( op->prpY+j+3*nvec_op, D_pt, pbuf+3*nvec, nvec, nvec_op, nvec ); D_pt += 9;
    // X dir
    j = 6*(*nb_pt)*nvec_op; nb_pt++;
    prn_X_PRECISION_new( pbuf, phi_pt, nvec, nvec, nvec_phi );
    mvmh_PRECISION_new( op->prpX+j, D_pt, pbuf, nvec, nvec_op, nvec );
    mvmh_PRECISION_new( op->prpX+j+3*nvec_op, D_pt, pbuf+3*nvec, nvec, nvec_op, nvec ); D_pt += 9;
  }
  // start communication in positive direction
  START_LOCKED_MASTER(threading)
  ghost_sendrecv_PRECISION_new( op->prpT, T, +1, &(op->c), _FULL_SYSTEM, l );
  ghost_sendrecv_PRECISION_new( op->prpZ, Z, +1, &(op->c), _FULL_SYSTEM, l );
  ghost_sendrecv_PRECISION_new( op->prpY, Y, +1, &(op->c), _FULL_SYSTEM, l );
  ghost_sendrecv_PRECISION_new( op->prpX, X, +1, &(op->c), _FULL_SYSTEM, l );
  // wait for communication in negative direction
  ghost_wait_PRECISION_new( op->prnT, T, -1, &(op->c), _FULL_SYSTEM, l );
  ghost_wait_PRECISION_new( op->prnZ, Z, -1, &(op->c), _FULL_SYSTEM, l );
  ghost_wait_PRECISION_new( op->prnY, Y, -1, &(op->c), _FULL_SYSTEM, l );
  ghost_wait_PRECISION_new( op->prnX, X, -1, &(op->c), _FULL_SYSTEM, l );
  END_LOCKED_MASTER(threading)
  // multiply with U and lift up minus dir
  for ( eta_pt=eta->vector_buffer+start*nvec_eta, end_pt=eta->vector_buffer+end*nvec_eta, D_pt = op->D+start*3, nb_pt=neighbor+(start/12)*4; eta_pt<end_pt; eta_pt+=12*nvec_eta ) {
    // T dir
    j = 6*(*nb_pt)*nvec_op; nb_pt++;
    mvm_PRECISION_new( pbuf, D_pt, op->prnT+j, nvec, nvec, nvec_op );
    mvm_PRECISION_new( pbuf+3*nvec, D_pt, op->prnT+j+3*nvec_op, nvec, nvec, nvec_op );
    pbp_su3_T_PRECISION_new( pbuf, eta_pt, nvec, nvec, nvec_eta ); D_pt += 9;
    // Z dir
    j = 6*(*nb_pt)*nvec_op; nb_pt++;
    mvm_PRECISION_new( pbuf, D_pt, op->prnZ+j, nvec, nvec, nvec_op );
    mvm_PRECISION_new( pbuf+3*nvec, D_pt, op->prnZ+j+3*nvec_op, nvec, nvec, nvec_op );
    pbp_su3_Z_PRECISION_new( pbuf, eta_pt, nvec, nvec, nvec_eta ); D_pt += 9;
    // Y dir
    j = 6*(*nb_pt)*nvec_op; nb_pt++;
    mvm_PRECISION_new( pbuf, D_pt, op->prnY+j, nvec, nvec, nvec_op );
    mvm_PRECISION_new( pbuf+3*nvec, D_pt, op->prnY+j+3*nvec_op, nvec, nvec, nvec_op );
    pbp_su3_Y_PRECISION_new( pbuf, eta_pt, nvec, nvec, nvec_eta ); D_pt += 9;
    // X dir
    j = 6*(*nb_pt)*nvec_op; nb_pt++;
    mvm_PRECISION_new( pbuf, D_pt, op->prnX+j, nvec, nvec, nvec_op );
    mvm_PRECISION_new( pbuf+3*nvec, D_pt, op->prnX+j+3*nvec_op, nvec, nvec, nvec_op );
    pbp_su3_X_PRECISION_new( pbuf, eta_pt, nvec, nvec, nvec_eta ); D_pt += 9;
  }
  
  // wait for communication in positive direction
  START_LOCKED_MASTER(threading)
  ghost_wait_PRECISION_new( op->prpT, T, +1, &(op->c), _FULL_SYSTEM, l );
  ghost_wait_PRECISION_new( op->prpZ, Z, +1, &(op->c), _FULL_SYSTEM, l );
  ghost_wait_PRECISION_new( op->prpY, Y, +1, &(op->c), _FULL_SYSTEM, l );
  ghost_wait_PRECISION_new( op->prpX, X, +1, &(op->c), _FULL_SYSTEM, l );
  END_LOCKED_MASTER(threading)
  // lift up plus dir
  for ( i=start*nvec_op/2, eta_pt=eta->vector_buffer+start*nvec_eta; i<end*nvec_op/2; i+=6*nvec_op, eta_pt+=12*nvec_eta ) {
    pbn_su3_T_PRECISION_new( op->prpT+i, eta_pt, nvec, nvec_op, nvec_eta );
    pbn_su3_Z_PRECISION_new( op->prpZ+i, eta_pt, nvec, nvec_op, nvec_eta );
    pbn_su3_Y_PRECISION_new( op->prpY+i, eta_pt, nvec, nvec_op, nvec_eta );
    pbn_su3_X_PRECISION_new( op->prpX+i, eta_pt, nvec, nvec_op, nvec_eta );
  }
/*#ifdef HAVE_TM1p1
  }
#endif*/
  
  START_MASTER(threading)
  PROF_PRECISION_STOP( _NC, 1 );
  END_MASTER(threading)

  SYNC_MASTER_TO_ALL(threading)
}

//eta = blockD*phi 
void block_d_plus_clover_PRECISION_new( vector_PRECISION *eta, vector_PRECISION *phi, int start, schwarz_PRECISION_struct *s, level_struct *l, struct Thread *threading ) {

  START_UNTHREADED_FUNCTION(threading)

  if ( eta->num_vect < phi->num_vect_now )
    error0("block_d_plus_clover_PRECISION: assumptions are not met\n");
  
  int n = s->num_block_sites, *length = s->dir_length, **index = s->index, *neighbor = s->op.neighbor_table, nv = l->num_lattice_site_var;
  int nvec = phi->num_vect_now, nvec_phi = phi->num_vect, nvec_eta = eta->num_vect;
  buffer_PRECISION lphi = phi->vector_buffer+start*nvec_phi, leta = eta->vector_buffer+start*nvec_eta;

  clover_PRECISION_new(eta, phi, &(s->op), start, start+nv*n, l, no_threading );

  int i, j, k, *ind;
  config_PRECISION D_pt;
  config_PRECISION D = s->op.D + (start/nv)*36;
/*#ifdef HAVE_TM1p1
  if( g.n_flavours == 2 ) {
    complex_PRECISION buf1[50]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, *buf2=buf1+12, *buf3=buf2+12, *buf4=buf3+12;
    // inner block couplings
    ind = index[T]; // T direction
    for ( i=0; i<length[T]; i++ ) {
      k = ind[i]; j = neighbor[4*k+T]; D_pt = D + 36*k + 9*T;
      dprn_T_PRECISION( buf1, lphi+24*k ); // (1+gamma_T) phi(x) + projection
      dprp_T_PRECISION( buf2, lphi+24*j ); // (1-gamma_T) phi(x+hat{T}) + projection
      mvmh_PRECISION( buf3, D_pt, buf1 );     // U_T^dagger(x) (1+gamma_T) phi(x) - projected
      mvmh_PRECISION( buf3+3, D_pt, buf1+3 ); // U_T^dagger(x) (1+gamma_T) phi(x) - projected
      mvmh_PRECISION( buf3+6, D_pt, buf1+6 ); // U_T^dagger(x) (1+gamma_T) phi(x) - projected
      mvmh_PRECISION( buf3+9, D_pt, buf1+9 ); // U_T^dagger(x) (1+gamma_T) phi(x) - projected
      mvm_PRECISION( buf4, D_pt, buf2 );      // U_T(x) (1-gamma_T) phi(x+hat{T}) - projected
      mvm_PRECISION( buf4+3, D_pt, buf2+3 );     // U_T(x) (1-gamma_T) phi(x+hat{T}) - projected
      mvm_PRECISION( buf4+6, D_pt, buf2+6 ); // U_T(x) (1-gamma_T) phi(x+hat{T}) - projected
      mvm_PRECISION( buf4+9, D_pt, buf2+9 ); // U_T(x) (1-gamma_T) phi(x+hat{T}) - projected
      dpbn_su3_T_PRECISION( buf3, leta+24*j ); // eta(x+hat{T}) -= U_T(x)^dagger(x) (1+gamma_T) phi(x) + lift back
      dpbp_su3_T_PRECISION( buf4, leta+24*k ); // eta(x) -= U_T(x) (1-gamma_T) phi(x+hat{T}) + lift back
    }
    ind = index[Z]; // Z direction
    for ( i=0; i<length[Z]; i++ ) {
      k = ind[i]; j = neighbor[4*k+Z]; D_pt = D + 36*k + 9*Z;
      dprn_Z_PRECISION( buf1, lphi+24*k ); // (1+gamma_Z) phi(x) + projection
      dprp_Z_PRECISION( buf2, lphi+24*j ); // (1-gamma_Z) phi(x+hat{Z}) + projection
      mvmh_PRECISION( buf3, D_pt, buf1 );     // U_Z^dagger(x) (1+gamma_Z) phi(x) - projected
      mvmh_PRECISION( buf3+3, D_pt, buf1+3 );     // U_Z^dagger(x) (1+gamma_Z) phi(x) - projected
      mvmh_PRECISION( buf3+6, D_pt, buf1+6 ); // U_Z^dagger(x) (1+gamma_Z) phi(x) - projected
      mvmh_PRECISION( buf3+9, D_pt, buf1+9 ); // U_Z^dagger(x) (1+gamma_Z) phi(x) - projected
      mvm_PRECISION( buf4, D_pt, buf2 );     // U_Z(x) (1-gamma_Z) phi(x+hat{Z}) - projected
      mvm_PRECISION( buf4+3, D_pt, buf2+3 );     // U_Z(x) (1-gamma_Z) phi(x+hat{Z}) - projected
      mvm_PRECISION( buf4+6, D_pt, buf2+6 ); // U_Z(x) (1-gamma_Z) phi(x+hat{Z}) - projected
      mvm_PRECISION( buf4+9, D_pt, buf2+9 ); // U_Z(x) (1-gamma_Z) phi(x+hat{Z}) - projected
      dpbn_su3_Z_PRECISION( buf3, leta+24*j ); // eta(x+hat{Z}) -= U_Z(x)^dagger(x) (1+gamma_Z) phi(x) + lift back
      dpbp_su3_Z_PRECISION( buf4, leta+24*k ); // eta(x) -= U_Z(x) (1-gamma_Z) phi(x+hat{Z}) + lift back
    }
    ind = index[Y]; // Y direction
    for ( i=0; i<length[Y]; i++ ) {
      k = ind[i]; j = neighbor[4*k+Y]; D_pt = D + 36*k + 9*Y;
      dprn_Y_PRECISION( buf1, lphi+24*k ); // (1+gamma_Y) phi(x) + projection
      dprp_Y_PRECISION( buf2, lphi+24*j ); // (1-gamma_Y) phi(x+hat{Y}) + projection
      mvmh_PRECISION( buf3, D_pt, buf1 );     // U_Y^dagger(x) (1+gamma_Y) phi(x) - projected
      mvmh_PRECISION( buf3+3, D_pt, buf1+3 );     // U_Y^dagger(x) (1+gamma_Y) phi(x) - projected
      mvmh_PRECISION( buf3+6, D_pt, buf1+6 ); // U_Y^dagger(x) (1+gamma_Y) phi(x) - projected
      mvmh_PRECISION( buf3+9, D_pt, buf1+9 ); // U_Y^dagger(x) (1+gamma_Y) phi(x) - projected
      mvm_PRECISION( buf4, D_pt, buf2 );     // U_Y(x) (1-gamma_Y) phi(x+hat{Y}) - projected
      mvm_PRECISION( buf4+3, D_pt, buf2+3 );     // U_Y(x) (1-gamma_Y) phi(x+hat{Y}) - projected
      mvm_PRECISION( buf4+6, D_pt, buf2+6 ); // U_Y(x) (1-gamma_Y) phi(x+hat{Y}) - projected
      mvm_PRECISION( buf4+9, D_pt, buf2+9 ); // U_Y(x) (1-gamma_Y) phi(x+hat{Y}) - projected
      dpbn_su3_Y_PRECISION( buf3, leta+24*j ); // eta(x+hat{Y}) -= U_Y(x)^dagger(x) (1+gamma_Y) phi(x) + lift back
      dpbp_su3_Y_PRECISION( buf4, leta+24*k ); // eta(x) -= U_Y(x) (1-gamma_Y) phi(x+hat{Y}) + lift back
    }
    ind = index[X]; // X direction
    for ( i=0; i<length[X]; i++ ) {
      k = ind[i]; j = neighbor[4*k+X]; D_pt = D + 36*k + 9*X;
      dprn_X_PRECISION( buf1, lphi+24*k ); // (1+gamma_X) phi(x) + projection
      dprp_X_PRECISION( buf2, lphi+24*j ); // (1-gamma_X) phi(x+hat{X}) + projection
      mvmh_PRECISION( buf3, D_pt, buf1 );     // U_X^dagger(x) (1+gamma_X) phi(x) - projected
      mvmh_PRECISION( buf3+3, D_pt, buf1+3 );     // U_X^dagger(x) (1+gamma_X) phi(x) - projected
      mvmh_PRECISION( buf3+6, D_pt, buf1+6 ); // U_X^dagger(x) (1+gamma_X) phi(x) - projected
      mvmh_PRECISION( buf3+9, D_pt, buf1+9 ); // U_X^dagger(x) (1+gamma_X) phi(x) - projected
      mvm_PRECISION( buf4, D_pt, buf2 );     // U_mu(x) (1-gamma_X) phi(x+hat{X}) - projected
      mvm_PRECISION( buf4+3, D_pt, buf2+3 );     // U_mu(x) (1-gamma_X) phi(x+hat{X}) - projected
      mvm_PRECISION( buf4+6, D_pt, buf2+6 ); // U_mu(x) (1-gamma_X) phi(x+hat{X}) - projected
      mvm_PRECISION( buf4+9, D_pt, buf2+9 ); // U_mu(x) (1-gamma_X) phi(x+hat{X}) - projected
      dpbn_su3_X_PRECISION( buf3, leta+24*j ); // eta(x+hat{X}) -= U_X(x)^dagger(x) (1+gamma_X) phi(x) + lift back
      dpbp_su3_X_PRECISION( buf4, leta+24*k ); // eta(x) -= U_X(x) (1-gamma_X) phi(x+hat{X}) + lift back
    }    
  } else {
#endif*/   
    //complex_PRECISION buf1[25]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, *buf2=buf1+6, *buf3=buf2+6, *buf4=buf3+6;
    //complex_PRECISION *buf1=NULL;
    //MALLOC( buf1, complex_PRECISION, 25*n_vect );
    //complex_PRECISION *buf2=buf1+6*n_vect, *buf3=buf2+6*n_vect, *buf4=buf3+6*n_vect;
  complex_PRECISION buf1[25*nvec], *buf2=buf1+6*nvec, *buf3=buf2+6*nvec, *buf4=buf3+6*nvec;//why 25?????? 4(buffers)x6(projected_spin*color)
  
  // inner block couplings
  ind = index[T]; // T direction
  for ( i=0; i<length[T]; i++ ) {
    k = ind[i]; j = neighbor[4*k+T]; D_pt = D + 36*k + 9*T;
    prn_T_PRECISION_new( buf1, lphi+12*k*nvec_phi, nvec, nvec, nvec_phi ); // (1+gamma_T) phi(x) + projection
    prp_T_PRECISION_new( buf2, lphi+12*j*nvec_phi, nvec, nvec, nvec_phi ); // (1-gamma_T) phi(x+hat{T}) + projection
    mvmh_PRECISION_new( buf3, D_pt, buf1, nvec, nvec, nvec );     // U_T^dagger(x) (1+gamma_T) phi(x) - projected
    mvmh_PRECISION_new( buf3+3*nvec, D_pt, buf1+3*nvec, nvec, nvec, nvec ); // U_T^dagger(x) (1+gamma_T) phi(x) - projected
    mvm_PRECISION_new( buf4, D_pt, buf2, nvec, nvec, nvec );     // U_T(x) (1-gamma_T) phi(x+hat{T}) - projected
    mvm_PRECISION_new( buf4+3*nvec, D_pt, buf2+3*nvec, nvec, nvec, nvec ); // U_T(x) (1-gamma_T) phi(x+hat{T}) - projected
    pbn_su3_T_PRECISION_new( buf3, leta+12*j*nvec_eta, nvec, nvec, nvec_eta); // eta(x+hat{T}) -= U_T(x)^dagger(x) (1+gamma_T) phi(x) + lift back
    pbp_su3_T_PRECISION_new( buf4, leta+12*k*nvec_eta, nvec, nvec, nvec_eta ); // eta(x) -= U_T(x) (1-gamma_T) phi(x+hat{T}) + lift back
  }
  ind = index[Z]; // Z direction
  for ( i=0; i<length[Z]; i++ ) {
    k = ind[i]; j = neighbor[4*k+Z]; D_pt = D + 36*k + 9*Z;
    prn_Z_PRECISION_new( buf1, lphi+12*k*nvec_phi,nvec, nvec, nvec_phi ); // (1+gamma_Z) phi(x) + projection
    prp_Z_PRECISION_new( buf2, lphi+12*j*nvec_phi,nvec, nvec, nvec_phi ); // (1-gamma_Z) phi(x+hat{Z}) + projection
    mvmh_PRECISION_new( buf3, D_pt, buf1, nvec, nvec, nvec );     // U_Z^dagger(x) (1+gamma_Z) phi(x) - projected
    mvmh_PRECISION_new( buf3+3*nvec, D_pt, buf1+3*nvec, nvec, nvec, nvec ); // U_Z^dagger(x) (1+gamma_Z) phi(x) - projected
    mvm_PRECISION_new( buf4, D_pt, buf2, nvec, nvec, nvec );     // U_Z(x) (1-gamma_Z) phi(x+hat{Z}) - projected
    mvm_PRECISION_new( buf4+3*nvec, D_pt, buf2+3*nvec, nvec, nvec, nvec ); // U_Z(x) (1-gamma_Z) phi(x+hat{Z}) - projected
    pbn_su3_Z_PRECISION_new( buf3, leta+12*j*nvec_eta, nvec, nvec, nvec_eta ); // eta(x+hat{Z}) -= U_Z(x)^dagger(x) (1+gamma_Z) phi(x) + lift back
    pbp_su3_Z_PRECISION_new( buf4, leta+12*k*nvec_eta, nvec, nvec, nvec_eta ); // eta(x) -= U_Z(x) (1-gamma_Z) phi(x+hat{Z}) + lift back
  }
  ind = index[Y]; // Y direction
  for ( i=0; i<length[Y]; i++ ) {
    k = ind[i]; j = neighbor[4*k+Y]; D_pt = D + 36*k + 9*Y;
    prn_Y_PRECISION_new( buf1, lphi+12*k*nvec_phi, nvec, nvec, nvec_phi ); // (1+gamma_Y) phi(x) + projection
    prp_Y_PRECISION_new( buf2, lphi+12*j*nvec_phi, nvec, nvec, nvec_phi ); // (1-gamma_Y) phi(x+hat{Y}) + projection
    mvmh_PRECISION_new( buf3, D_pt, buf1, nvec, nvec, nvec );     // U_Y^dagger(x) (1+gamma_Y) phi(x) - projected
    mvmh_PRECISION_new( buf3+3*nvec, D_pt, buf1+3*nvec, nvec, nvec, nvec ); // U_Y^dagger(x) (1+gamma_Y) phi(x) - projected
    mvm_PRECISION_new( buf4, D_pt, buf2, nvec, nvec, nvec );     // U_Y(x) (1-gamma_Y) phi(x+hat{Y}) - projected
    mvm_PRECISION_new( buf4+3*nvec, D_pt, buf2+3*nvec, nvec, nvec, nvec ); // U_Y(x) (1-gamma_Y) phi(x+hat{Y}) - projected
    pbn_su3_Y_PRECISION_new( buf3, leta+12*j*nvec_eta, nvec, nvec, nvec_eta ); // eta(x+hat{Y}) -= U_Y(x)^dagger(x) (1+gamma_Y) phi(x) + lift back
    pbp_su3_Y_PRECISION_new( buf4, leta+12*k*nvec_eta, nvec, nvec, nvec_eta ); // eta(x) -= U_Y(x) (1-gamma_Y) phi(x+hat{Y}) + lift back
  }
  ind = index[X]; // X direction
  for ( i=0; i<length[X]; i++ ) {
    k = ind[i]; j = neighbor[4*k+X]; D_pt = D + 36*k + 9*X;
    prn_X_PRECISION_new( buf1, lphi+12*k*nvec_phi, nvec, nvec, nvec_phi ); // (1+gamma_X) phi(x) + projection
    prp_X_PRECISION_new( buf2, lphi+12*j*nvec_phi, nvec, nvec, nvec_phi ); // (1-gamma_X) phi(x+hat{X}) + projection
    mvmh_PRECISION_new( buf3, D_pt, buf1, nvec, nvec, nvec );     // U_X^dagger(x) (1+gamma_X) phi(x) - projected
    mvmh_PRECISION_new( buf3+3*nvec, D_pt, buf1+3*nvec, nvec, nvec, nvec ); // U_X^dagger(x) (1+gamma_X) phi(x) - projected
    mvm_PRECISION_new( buf4, D_pt, buf2, nvec, nvec, nvec );     // U_mu(x) (1-gamma_X) phi(x+hat{X}) - projected
    mvm_PRECISION_new( buf4+3*nvec, D_pt, buf2+3*nvec, nvec, nvec, nvec ); // U_mu(x) (1-gamma_X) phi(x+hat{X}) - projected
    pbn_su3_X_PRECISION_new( buf3, leta+12*j*nvec_eta, nvec, nvec, nvec_eta ); // eta(x+hat{X}) -= U_X(x)^dagger(x) (1+gamma_X) phi(x) + lift back
    pbp_su3_X_PRECISION_new( buf4, leta+12*k*nvec_eta, nvec, nvec, nvec_eta ); // eta(x) -= U_X(x) (1-gamma_X) phi(x+hat{X}) + lift back
  }
/*#ifdef HAVE_TM1p1
  }
#endif*/
  END_UNTHREADED_FUNCTION(threading)
}

/********************  aggreagate operators ************************************************/

// used in coarse_operator_PRECISION_setup: no_threading  
void diagonal_aggregate_PRECISION_new( vector_PRECISION *eta1, vector_PRECISION *eta2, vector_PRECISION *phi, config_PRECISION diag, level_struct *l ) {

  int i, j, jj, nvec = phi->num_vect_now, nvec_eta1 = eta1->num_vect, nvec_eta2 = eta2->num_vect, nvec_phi = phi->num_vect;
  buffer_PRECISION eta_end = eta1->vector_buffer + l->inner_vector_size*nvec_eta1;
  buffer_PRECISION eta1_pt = eta1->vector_buffer, eta2_pt = eta2->vector_buffer, phi_pt = phi->vector_buffer;

  if ( nvec_eta1 != nvec_eta2 )
    error0("diagonal_aggregate_PRECISION: assumptions are not met\n");

  while ( eta1_pt < eta_end ) {
    for( i=0; i<6; i++ ){
      VECTOR_LOOP(j, nvec, jj, *eta1_pt = (*phi_pt)*(*diag); *eta2_pt = _COMPLEX_PRECISION_ZERO; eta1_pt++; eta2_pt++; phi_pt++;)
      eta1_pt += nvec_eta1-nvec; eta2_pt += nvec_eta2-nvec; phi_pt += nvec_phi - nvec;
      diag++;
    }
    for( i=0; i<6; i++ ){
      VECTOR_LOOP(j, nvec, jj, *eta2_pt = (*phi_pt)*(*diag); *eta1_pt = _COMPLEX_PRECISION_ZERO; eta1_pt++; eta2_pt++; phi_pt++;)
      eta1_pt += nvec_eta1-nvec; eta2_pt += nvec_eta2-nvec; phi_pt += nvec_phi - nvec;
      diag++;
    }
  }
}

// <- phi????
// used in coarse_operator_PRECISION_setup: no_threading
void d_plus_clover_aggregate_PRECISION_new( vector_PRECISION *eta1, vector_PRECISION *eta2, vector_PRECISION *phi, schwarz_PRECISION_struct *s, level_struct *l ) {
  
  int i, length, index1, index2, *index_dir, *neighbor = s->op.neighbor_table;
  int nvec = phi->num_vect_now, nvec_eta1 = eta1->num_vect, nvec_eta2 = eta2->num_vect, nvec_phi = phi->num_vect;
  buffer_PRECISION eta1_pt, eta2_pt, phi_pt;
  complex_PRECISION buffer1[12*nvec], buffer2[12*nvec];
  config_PRECISION D_pt, D = s->op.D;
    
  if ( nvec_eta1 != nvec_eta2 || nvec_eta1 < nvec_phi )
    error0("d_plus_clover_aggregate_PRECISION: assumptions are not met\n");

  // add clover term/shift
  spin0and1_clover_PRECISION_new( eta1, phi, s->op.clover, l );
  spin2and3_clover_PRECISION_new( eta2, phi, s->op.clover, l );

  // T dir
  length = l->is_PRECISION.agg_length[T]; index_dir = l->is_PRECISION.agg_index[T];
  for ( i=0; i<length; i++ ) {
    index1 = index_dir[i]; index2 = neighbor[4*index1 + T];
    phi_pt = phi->vector_buffer + 12*index2*nvec_phi; D_pt = D + 36*index1+9*T;
    mvm_PRECISION_new( buffer1, D_pt, phi_pt, nvec, nvec, nvec_phi );
    mvm_PRECISION_new( buffer1+3*nvec, D_pt, phi_pt+3*nvec_phi, nvec, nvec, nvec_phi );
    mvm_PRECISION_new( buffer1+6*nvec, D_pt, phi_pt+6*nvec_phi, nvec, nvec, nvec_phi );
    mvm_PRECISION_new( buffer1+9*nvec, D_pt, phi_pt+9*nvec_phi, nvec, nvec, nvec_phi );
    phi_pt = phi->vector_buffer + 12*index1*nvec_phi;
    mvmh_PRECISION_new( buffer2, D_pt, phi_pt, nvec, nvec, nvec_phi );
    mvmh_PRECISION_new( buffer2+3*nvec, D_pt, phi_pt+3*nvec_phi, nvec, nvec, nvec_phi );
    mvmh_PRECISION_new( buffer2+6*nvec, D_pt, phi_pt+6*nvec_phi, nvec, nvec, nvec_phi );
    mvmh_PRECISION_new( buffer2+9*nvec, D_pt, phi_pt+9*nvec_phi, nvec, nvec, nvec_phi );
    eta1_pt = eta1->vector_buffer + 12*index1*nvec_eta1; eta2_pt = eta2->vector_buffer + 12*index1*nvec_eta2;
    twospin_p_T_PRECISION_new( eta1_pt, eta2_pt, buffer1, nvec, nvec_eta1, nvec_eta2, nvec );
    eta1_pt = eta1->vector_buffer + 12*index2*nvec_eta1; eta2_pt = eta2->vector_buffer + 12*index2*nvec_eta2;
    twospin_n_T_PRECISION_new( eta1_pt, eta2_pt, buffer2, nvec, nvec_eta1, nvec_eta2, nvec );
  }
  // Z dir
  length = l->is_PRECISION.agg_length[Z]; index_dir = l->is_PRECISION.agg_index[Z];
  for ( i=0; i<length; i++ ) {
    index1 = index_dir[i]; index2 = neighbor[4*index1 + Z];
    phi_pt = phi->vector_buffer + 12*index2*nvec_phi; D_pt = D + 36*index1+9*Z;
    mvm_PRECISION_new( buffer1, D_pt, phi_pt, nvec, nvec, nvec_phi );
    mvm_PRECISION_new( buffer1+3*nvec, D_pt, phi_pt+3*nvec_phi, nvec, nvec, nvec_phi );
    mvm_PRECISION_new( buffer1+6*nvec, D_pt, phi_pt+6*nvec_phi, nvec, nvec, nvec_phi );
    mvm_PRECISION_new( buffer1+9*nvec, D_pt, phi_pt+9*nvec_phi, nvec, nvec, nvec_phi );
    phi_pt = phi->vector_buffer + 12*index1*nvec_phi;
    mvmh_PRECISION_new( buffer2, D_pt, phi_pt, nvec, nvec, nvec_phi );
    mvmh_PRECISION_new( buffer2+3*nvec, D_pt, phi_pt+3*nvec_phi, nvec, nvec, nvec_phi );
    mvmh_PRECISION_new( buffer2+6*nvec, D_pt, phi_pt+6*nvec_phi, nvec, nvec, nvec_phi );
    mvmh_PRECISION_new( buffer2+9*nvec, D_pt, phi_pt+9*nvec_phi, nvec, nvec, nvec_phi );
    eta1_pt = eta1->vector_buffer + 12*index1*nvec_eta1; eta2_pt = eta2->vector_buffer + 12*index1*nvec_eta2;
    twospin_p_Z_PRECISION_new( eta1_pt, eta2_pt, buffer1, nvec, nvec_eta1, nvec_eta2, nvec );
    eta1_pt = eta1->vector_buffer + 12*index2*nvec_eta1; eta2_pt = eta2->vector_buffer + 12*index2*nvec_eta2;
    twospin_n_Z_PRECISION_new( eta1_pt, eta2_pt, buffer2, nvec, nvec_eta1, nvec_eta2, nvec );
  }
  // Y dir
  length = l->is_PRECISION.agg_length[Y]; index_dir = l->is_PRECISION.agg_index[Y];
  for ( i=0; i<length; i++ ) {
    index1 = index_dir[i]; index2 = neighbor[4*index1 + Y];
    phi_pt = phi->vector_buffer + 12*index2*nvec_phi; D_pt = D + 36*index1+9*Y;
    mvm_PRECISION_new( buffer1, D_pt, phi_pt, nvec, nvec, nvec_phi );
    mvm_PRECISION_new( buffer1+3*nvec, D_pt, phi_pt+3*nvec_phi, nvec, nvec, nvec_phi );
    mvm_PRECISION_new( buffer1+6*nvec, D_pt, phi_pt+6*nvec_phi, nvec, nvec, nvec_phi );
    mvm_PRECISION_new( buffer1+9*nvec, D_pt, phi_pt+9*nvec_phi, nvec, nvec, nvec_phi );
    phi_pt = phi->vector_buffer + 12*index1*nvec_phi;
    mvmh_PRECISION_new( buffer2, D_pt, phi_pt, nvec, nvec, nvec_phi );
    mvmh_PRECISION_new( buffer2+3*nvec, D_pt, phi_pt+3*nvec_phi, nvec, nvec, nvec_phi );
    mvmh_PRECISION_new( buffer2+6*nvec, D_pt, phi_pt+6*nvec_phi, nvec, nvec, nvec_phi );
    mvmh_PRECISION_new( buffer2+9*nvec, D_pt, phi_pt+9*nvec_phi, nvec, nvec, nvec_phi );
    eta1_pt = eta1->vector_buffer + 12*index1*nvec_eta1; eta2_pt = eta2->vector_buffer + 12*index1*nvec_eta2;
    twospin_p_Y_PRECISION_new( eta1_pt, eta2_pt, buffer1, nvec, nvec_eta1, nvec_eta2, nvec );
    eta1_pt = eta1->vector_buffer + 12*index2*nvec_eta1; eta2_pt = eta2->vector_buffer + 12*index2*nvec_eta2;
    twospin_n_Y_PRECISION_new( eta1_pt, eta2_pt, buffer2, nvec, nvec_eta1, nvec_eta2 , nvec );
  }
  // X dir
  length = l->is_PRECISION.agg_length[X]; index_dir = l->is_PRECISION.agg_index[X];
  for ( i=0; i<length; i++ ) {
    index1 = index_dir[i]; index2 = neighbor[4*index1 + X];
    phi_pt = phi->vector_buffer + 12*index2*nvec_phi; D_pt = D + 36*index1+9*X;
    mvm_PRECISION_new( buffer1, D_pt, phi_pt, nvec, nvec, nvec_phi );
    mvm_PRECISION_new( buffer1+3*nvec, D_pt, phi_pt+3*nvec_phi, nvec, nvec, nvec_phi );
    mvm_PRECISION_new( buffer1+6*nvec, D_pt, phi_pt+6*nvec_phi, nvec, nvec, nvec_phi );
    mvm_PRECISION_new( buffer1+9*nvec, D_pt, phi_pt+9*nvec_phi, nvec, nvec, nvec_phi );
    phi_pt = phi->vector_buffer + 12*index1*nvec_phi;
    mvmh_PRECISION_new( buffer2, D_pt, phi_pt, nvec, nvec, nvec_phi );
    mvmh_PRECISION_new( buffer2+3*nvec, D_pt, phi_pt+3*nvec_phi, nvec, nvec, nvec_phi );
    mvmh_PRECISION_new( buffer2+6*nvec, D_pt, phi_pt+6*nvec_phi, nvec, nvec, nvec_phi );
    mvmh_PRECISION_new( buffer2+9*nvec, D_pt, phi_pt+9*nvec_phi, nvec, nvec, nvec_phi );
    eta1_pt = eta1->vector_buffer + 12*index1*nvec_eta1; eta2_pt = eta2->vector_buffer + 12*index1*nvec_eta2;
    twospin_p_X_PRECISION_new( eta1_pt, eta2_pt, buffer1, nvec, nvec_eta1, nvec_eta2 , nvec );
    eta1_pt = eta1->vector_buffer + 12*index2*nvec_eta1; eta2_pt = eta2->vector_buffer + 12*index2*nvec_eta2;
    twospin_n_X_PRECISION_new( eta1_pt, eta2_pt, buffer2, nvec, nvec_eta1, nvec_eta2 , nvec );
  }
}

// used in coarse_operator_PRECISION_setup: no_threading  
void d_neighbor_aggregate_PRECISION_new( vector_PRECISION *eta1, vector_PRECISION *eta2, vector_PRECISION *phi, const int mu, schwarz_PRECISION_struct *s, level_struct *l ) {
  
  int i, length, index1, index2, *index_dir, *neighbor;
  int nvec = phi->num_vect_now, nvec_eta1 = eta1->num_vect, nvec_eta2 = eta2->num_vect, nvec_phi = phi->num_vect;
  buffer_PRECISION eta1_pt, eta2_pt, phi_pt;
  complex_PRECISION buffer1[12*nvec];
  config_PRECISION D_pt, D = s->op.D;
  
  if ( nvec_eta1 != nvec_eta2 || nvec_eta1 < nvec_phi )
    error0("d_neighbor_aggregate_PRECISION: assumptions are not met\n");

  length    = l->is_PRECISION.agg_boundary_length[mu];
  index_dir = l->is_PRECISION.agg_boundary_index[mu];
  neighbor  = l->is_PRECISION.agg_boundary_neighbor[mu];
  
  // requires the positive boundaries of phi to be communicated befor
  if ( mu == T ) {
    // T dir
    for ( i=0; i<length; i++ ) {
      index1 = index_dir[i]; index2 = neighbor[i];
      phi_pt = phi->vector_buffer + 12*index2*nvec_phi; D_pt = D + 36*index1 + 9*T;
      mvm_PRECISION_new( buffer1, D_pt, phi_pt, nvec, nvec, nvec_phi );
      mvm_PRECISION_new( buffer1+3*nvec, D_pt, phi_pt+3*nvec_phi, nvec, nvec, nvec_phi );
      mvm_PRECISION_new( buffer1+6*nvec, D_pt, phi_pt+6*nvec_phi, nvec, nvec, nvec_phi );
      mvm_PRECISION_new( buffer1+9*nvec, D_pt, phi_pt+9*nvec_phi, nvec, nvec, nvec_phi );
      eta1_pt = eta1->vector_buffer + 12*index1*nvec_eta1; eta2_pt = eta2->vector_buffer + 12*index1*nvec_eta2;   
      twospin2_p_T_PRECISION_new( eta1_pt, eta2_pt, buffer1, nvec, nvec_eta1, nvec_eta2, nvec );
    }
  } else if ( mu == Z ) {
    // Z dir
    for ( i=0; i<length; i++ ) {
      index1 = index_dir[i]; index2 = neighbor[i];
      phi_pt = phi->vector_buffer + 12*index2*nvec_phi; D_pt = D + 36*index1 + 9*Z;
      mvm_PRECISION_new( buffer1, D_pt, phi_pt, nvec, nvec, nvec_phi );
      mvm_PRECISION_new( buffer1+3*nvec, D_pt, phi_pt+3*nvec_phi, nvec, nvec, nvec_phi );
      mvm_PRECISION_new( buffer1+6*nvec, D_pt, phi_pt+6*nvec_phi, nvec, nvec, nvec_phi );
      mvm_PRECISION_new( buffer1+9*nvec, D_pt, phi_pt+9*nvec_phi, nvec, nvec, nvec_phi );
      eta1_pt = eta1->vector_buffer + 12*index1*nvec_eta1; eta2_pt = eta2->vector_buffer + 12*index1*nvec_eta2;
      twospin2_p_Z_PRECISION_new( eta1_pt, eta2_pt, buffer1, nvec, nvec_eta1, nvec_eta2, nvec );
    }
  } else if ( mu == Y ) {
    // Y dir
    for ( i=0; i<length; i++ ) {
      index1 = index_dir[i]; index2 = neighbor[i];
      phi_pt = phi->vector_buffer + 12*index2*nvec_phi; D_pt = D + 36*index1 + 9*Y;
      mvm_PRECISION_new( buffer1, D_pt, phi_pt, nvec, nvec, nvec_phi );
      mvm_PRECISION_new( buffer1+3*nvec, D_pt, phi_pt+3*nvec_phi, nvec, nvec, nvec_phi );
      mvm_PRECISION_new( buffer1+6*nvec, D_pt, phi_pt+6*nvec_phi, nvec, nvec, nvec_phi );
      mvm_PRECISION_new( buffer1+9*nvec, D_pt, phi_pt+9*nvec_phi, nvec, nvec, nvec_phi );
      eta1_pt = eta1->vector_buffer + 12*index1*nvec_eta1; eta2_pt = eta2->vector_buffer + 12*index1*nvec_eta2;
      twospin2_p_Y_PRECISION_new( eta1_pt, eta2_pt, buffer1, nvec, nvec_eta1, nvec_eta2, nvec );
    }
  } else if ( mu == X ) {
    // X dir
    for ( i=0; i<length; i++ ) {
      index1 = index_dir[i]; index2 = neighbor[i];
      phi_pt = phi->vector_buffer + 12*index2*nvec_phi; D_pt = D + 36*index1 + 9*X;
      mvm_PRECISION_new( buffer1, D_pt, phi_pt, nvec, nvec, nvec_phi );
      mvm_PRECISION_new( buffer1+3*nvec, D_pt, phi_pt+3*nvec_phi, nvec, nvec, nvec_phi );
      mvm_PRECISION_new( buffer1+6*nvec, D_pt, phi_pt+6*nvec_phi, nvec, nvec, nvec_phi );
      mvm_PRECISION_new( buffer1+9*nvec, D_pt, phi_pt+9*nvec_phi, nvec, nvec, nvec_phi );
      eta1_pt = eta1->vector_buffer + 12*index1*nvec_eta1; eta2_pt = eta2->vector_buffer + 12*index1*nvec_eta2;
      twospin2_p_X_PRECISION_new( eta1_pt, eta2_pt, buffer1, nvec, nvec_eta1, nvec_eta2, nvec );
    }
  }
}
 
// does this not affect eta and phi after returning to the calling pt?
void apply_twisted_bc_to_vector_PRECISION_new( vector_PRECISION *eta, vector_PRECISION *phi, double *theta, level_struct *l) {
  int t, z, y, x, i, j, jj, nvec = phi->num_vect_now, nvec_eta = eta->num_vect, nvec_phi = phi->num_vect;
  int *gl=l->global_lattice, sl[4];
  double phase[4];
  buffer_PRECISION eta_pt = eta->vector_buffer, phi_pt = phi->vector_buffer;
  complex_double twisted_bc;
  for (i=0; i<4; i++)
    sl[i] = l->local_lattice[i]*g.my_coords[i];
  
  if ( nvec_eta < nvec )
    error0("apply_twisted_bc_to_vector_PRECISION: assumptions are not met\n");

  for (t=0; t<l->local_lattice[0]; t++) {
    phase[T] = theta[T]*((double)sl[T]+t)/(double)gl[T];
    for (z=0; z<l->local_lattice[1]; z++) {
      phase[Z] = phase[T] + theta[Z]*((double)sl[Z]+z)/(double)gl[Z];
      for (y=0; y<l->local_lattice[2]; y++) {
        phase[Y] = phase[Z] + theta[Y]*((double)sl[Y]+y)/(double)gl[Y];
        for (x=0; x<l->local_lattice[3]; x++) {
          phase[X] = phase[Y] + theta[X]*((double)sl[X]+x)/(double)gl[X];
          twisted_bc = exp(I*phase[X]);
/*#ifdef HAVE_TM1p1
          if( g.n_flavours == 2 ) {
            FOR24( *eta->vector_buffer = (*phi->vector_buffer)*twisted_bc; phi->vector_buffer++; eta->vector_buffer++; );
          } else
#endif*/
	  for (i=0; i<12; i++){
	    VECTOR_LOOP(j, nvec, jj, *eta_pt = (*phi_pt)*twisted_bc; phi_pt++;  eta_pt++;)
	    eta_pt += nvec_eta-nvec; phi_pt += nvec_phi-nvec;
          }
        }
      }
    }
  }
}


/*************************** setup and update ************************************/

// used only in DDalphaAMG_interface.c: so not changed!!!!!!!!!!!!
void operator_updates_PRECISION( level_struct *l, struct Thread *threading ) {

  if ( l->level > 0 ) {
    if ( !l->idle ) {
      START_LOCKED_MASTER(threading)
      coarse_operator_PRECISION_setup_new( &(l->is_PRECISION.interpolation_vec), l );
      conf_PRECISION_gather( &(l->next_level->s_PRECISION.op), &(l->next_level->op_PRECISION), l->next_level );
      END_LOCKED_MASTER(threading)
      if ( !l->next_level->idle && l->next_level->level > 0 ) {
        START_LOCKED_MASTER(threading)
        schwarz_PRECISION_boundary_update( &(l->next_level->s_PRECISION), l->next_level );
        END_LOCKED_MASTER(threading)
        if ( g.method >= 4 && g.odd_even ) {
          coarse_oddeven_setup_PRECISION( &(l->next_level->s_PRECISION.op), _REORDER, l->next_level, threading );
        } else {
          coarse_operator_PRECISION_set_couplings( &(l->next_level->s_PRECISION.op), l->next_level, threading );
        }
      }
      if ( !l->next_level->idle && l->next_level->level == 0 && g.odd_even ) {
        coarse_oddeven_setup_PRECISION( &(l->next_level->s_PRECISION.op), _NO_REORDERING, l->next_level, threading );
      } else if ( !l->next_level->idle && l->next_level->level == 0 ) {
        coarse_operator_PRECISION_set_couplings( &(l->next_level->s_PRECISION.op), l->next_level, threading );
      }
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

void tm_term_PRECISION_setup( PRECISION mu, PRECISION even, PRECISION odd, operator_PRECISION_struct *op,
                              level_struct *l, struct Thread *threading ) {
   
#ifdef HAVE_TM
  if(threading->thread != 0)
    return;

  config_PRECISION tm_term = op->tm_term;
  if ( tm_term != NULL ) {
    config_PRECISION odd_proj = op->odd_proj;
    complex_PRECISION shift = I*mu;
    complex_PRECISION even_shift = I*even;
    complex_PRECISION odd_shift = I*odd;

    START_MASTER(threading)
    op->mu = mu;
    op->mu_even_shift = even;
    op->mu_odd_shift = odd;
    END_MASTER(threading)

    int i, j;
    int start, end;
    compute_core_start_end(0, l->num_inner_lattice_sites, &start, &end, l, threading);
    int n = end-start;
          
    if ( l->depth == 0 ) {
      complex_PRECISION tm_shift;
      tm_term += start*12;
      odd_proj += start*12;
      
      for ( i=0; i<n; i++ ) {
        if( cimag(even_shift) == 0. && cimag(odd_shift) == 0. )
          tm_shift = shift;
        else
          tm_shift = shift + even_shift + odd_proj[0]*(odd_shift - even_shift);
        FOR6( *tm_term = - tm_shift; tm_term++; )
        FOR6( *tm_term = tm_shift; tm_term++; )
        odd_proj += 12;
      }
    } else {
      int k, m  = l->num_parent_eig_vect;
      int tm_size = m*(m+1);
      
      tm_term += start*tm_size;
      odd_proj += start*tm_size;
      
      if( cimag(even_shift) == 0. && cimag(odd_shift) == 0. ) {
        
        for ( i=0; i<n; i++ ) {
          for ( j=0; j<m; j++ ) {
            for ( k=0; k<j; k++ )
              tm_term[k] = _COMPLEX_PRECISION_ZERO;
            tm_term += j;
            *tm_term = -1.* shift;
            tm_term++;
          }
          
          for ( j=0; j<m; j++ ) {
            for ( k=0; k<j; k++ )
              tm_term[k] = _COMPLEX_PRECISION_ZERO;
            tm_term += j;
            *tm_term = shift;
            tm_term++;
          }
        }
      } else {
        complex_PRECISION odd_factor = odd_shift - even_shift;
        
        for ( i=0; i<n; i++ ) {
          for ( j=0; j<m; j++ ) {
            for ( k=0; k<j; k++ ) 
              tm_term[k] = -1. * odd_factor * odd_proj[k] ;
            tm_term += j;
            odd_proj += j;
            *tm_term = -1.* ( shift + even_shift + odd_factor * (*odd_proj));
            tm_term++;
            odd_proj++;
          } 
          
          for ( j=0; j<m; j++ ) {
            for ( k=0; k<j; k++ ) 
              tm_term[k] = odd_factor * odd_proj[k] ;
            tm_term += j;
            odd_proj += j;
            *tm_term = ( shift + even_shift + odd_factor * (*odd_proj));
            tm_term++;
            odd_proj++;
          } 
        }
      }
    }  
  }
#endif
}

void epsbar_term_PRECISION_setup( PRECISION epsbar, PRECISION even, PRECISION odd, operator_PRECISION_struct *op,
                                  level_struct *l, struct Thread *threading ) {
  
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
    compute_core_start_end(0, l->num_inner_lattice_sites, &start, &end, l, threading);
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
        for ( i=0; i<2*n; i++ ) {
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

// eta <- gamma5*phi: not used anywhere
void gamma5_PRECISION_new( vector_PRECISION *eta, vector_PRECISION *phi, level_struct *l, struct Thread *threading ) {
  
  ASSERT(l->depth == 0);

  if ( eta->num_vect < phi->num_vect_now )
    error0("gamma5_PRECISION: assumptions are not met\n");

  int i, j, jj, nvec = phi->num_vect_now, nvec_phi = phi->num_vect, nvec_eta = eta->num_vect;
  buffer_PRECISION eta_end = eta->vector_buffer + threading->end_index[l->depth]*nvec_eta;
  buffer_PRECISION leta = eta->vector_buffer + threading->start_index[l->depth]*nvec_eta, lphi = phi->vector_buffer + threading->start_index[l->depth]*nvec_phi;
/*#ifdef HAVE_TM1p1
  if ( g.n_flavours == 2 ) {
    while ( leta < eta_end ) {
      FOR12( *leta = -(*lphi); lphi++; leta++; )
      FOR12( *leta =  (*lphi); lphi++; leta++; )
    }
  } else
#endif*/
  while ( leta < eta_end ) {
    for ( i=0; i<6; i++ ) {
      VECTOR_LOOP( j, nvec, jj, *leta = -(*lphi); lphi++; leta++;)
      leta += nvec_eta - nvec; lphi += nvec_phi - nvec;
      }
    for ( i=0; i<6; i++ ) {
      VECTOR_LOOP( j, nvec, jj, *leta = (*lphi); lphi++; leta++;)
      leta += nvec_eta - nvec; lphi += nvec_phi - nvec;
    }
  }
}

// used in DDalphaAMG_interface 
void gamma5_set_even_to_zero_PRECISION( vector_PRECISION *eta, vector_PRECISION *phi, level_struct *l, struct Thread *threading ) {

  ASSERT(l->depth == 0);
  
  int i = threading->start_site[l->depth];
  buffer_PRECISION eta_end = eta->vector_buffer + threading->end_index[l->depth];
  buffer_PRECISION leta = eta->vector_buffer + threading->start_index[l->depth], lphi = phi->vector_buffer + threading->start_index[l->depth];

#ifdef HAVE_TM1p1
  if ( g.n_flavours == 2 )
    while ( leta < eta_end ) {
      if(g.odd_even_table[i]==_ODD){
        FOR12( *leta = -(*lphi); lphi++; leta++; );
        FOR12( *leta = (*lphi); lphi++; leta++; );
      }
      else if(g.odd_even_table[i]==_EVEN){
        FOR24( *leta = _COMPLEX_PRECISION_ZERO; lphi++; leta++; );
      }
      i++;
    }
  else
#endif
    while ( leta < eta_end ) {
      if(g.odd_even_table[i]==_ODD){
        FOR6( *leta = -(*lphi); lphi++; leta++; );
        FOR6( *leta = (*lphi); lphi++; leta++; );
      }
      else if(g.odd_even_table[i]==_EVEN){
        FOR12( *leta = _COMPLEX_PRECISION_ZERO; lphi++; leta++; );
      }
      i++;
    }
}

// used in DDalphaAMG_interface
void gamma5_set_odd_to_zero_PRECISION( vector_PRECISION *eta, vector_PRECISION *phi, level_struct *l, struct Thread *threading ) {
  
  int i = threading->start_site[l->depth];
  buffer_PRECISION eta_end = eta->vector_buffer + threading->end_index[l->depth];
  buffer_PRECISION leta = eta->vector_buffer + threading->start_index[l->depth], lphi = phi->vector_buffer + threading->start_index[l->depth];

#ifdef HAVE_TM1p1
  if ( g.n_flavours == 2 )
    while ( leta < eta_end ) {
      if(g.odd_even_table[i]==_EVEN){
        FOR12( *leta = -(*lphi); lphi++; leta++; );
        FOR12( *leta = (*lphi); lphi++; leta++; );
      }
      else if(g.odd_even_table[i]==_ODD){
        FOR24( *leta = 0; lphi++; leta++; );
      }
      i++;
    }
  else
#endif
    while ( leta < eta_end ) {
      if(g.odd_even_table[i]==_EVEN){
        FOR6( *leta = -(*lphi); lphi++; leta++; );
        FOR6( *leta = (*lphi); lphi++; leta++; );
      }
      else if(g.odd_even_table[i]==_ODD){
        FOR12( *leta = 0; lphi++; leta++; );
      }
      i++;
    }
}

/****************** TM1p1 ******************************************************************/

void tau1_gamma5_PRECISION( vector_PRECISION *eta, vector_PRECISION *phi, level_struct *l, struct Thread *threading ) {
  
  ASSERT(l->depth == 0);
  error0("NOT CORRECTED\n");
#ifdef HAVE_TM1p1
  if ( g.n_flavours == 2 ) {
    buffer_PRECISION eta_end = eta->vector_buffer + threading->end_index[l->depth];
    buffer_PRECISION leta = eta->vector_buffer + threading->start_index[l->depth], lphi = phi->vector_buffer + threading->start_index[l->depth];
    complex_PRECISION b[6];
    while ( leta < eta_end ) {
      int i = 0;
      FOR6( b[i] =  (*lphi); lphi++; i++;   );
      FOR6( *leta = -(*lphi); lphi++; leta++; );
      i = 0;
      FOR6( *leta = - b[i] ; leta++; i++;   );
      i = 0;
      FOR6( b[i] =  (*lphi); lphi++; i++;   );
      FOR6( *leta =  (*lphi); lphi++; leta++; );
      i = 0;
      FOR6( *leta =   b[i] ; leta++; i++;   );
    }
  } else 
#endif
    {
      START_MASTER(threading)
      warning0("tau1_gamma5_PRECISION called with g.n_flavours != 2\n");
      END_MASTER(threading)
      gamma5_PRECISION_new( eta, phi, l, threading );
    }
}

// used in DDalphaAMG_interface    
void tau1_gamma5_set_even_to_zero_PRECISION( vector_PRECISION *eta, vector_PRECISION *phi, level_struct *l, struct Thread *threading ) {

  ASSERT(l->depth == 0);
  
#ifdef HAVE_TM1p1
  if ( g.n_flavours == 2 ) {
    int i = threading->start_site[l->depth];
    buffer_PRECISION eta_end = eta->vector_buffer + threading->end_index[l->depth];
    buffer_PRECISION leta = eta->vector_buffer + threading->start_index[l->depth], lphi = phi->vector_buffer + threading->start_index[l->depth];

    complex_PRECISION b[6];
    while ( leta < eta_end ) {
      if(g.odd_even_table[i]==_ODD){
        int i = 0;
        FOR6( b[i] =  (*lphi); lphi++; i++;   );
        FOR6( *leta = -(*lphi); lphi++; leta++; );
        i = 0;
        FOR6( *leta = - b[i] ; leta++; i++;   );
        i = 0;
        FOR6( b[i] =  (*lphi); lphi++; i++;   );
        FOR6( *leta =  (*lphi); lphi++; leta++; );
        i = 0;
        FOR6( *leta =   b[i] ; leta++; i++;   );
      } else if(g.odd_even_table[i]==_EVEN){
        FOR24( *leta = _COMPLEX_PRECISION_ZERO; lphi++; leta++; );
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

  ASSERT(l->depth == 0);
  
#ifdef HAVE_TM1p1
  if ( g.n_flavours == 2 ) {
    int i = threading->start_site[l->depth];
    buffer_PRECISION eta_end = eta->vector_buffer + threading->end_index[l->depth];
    buffer_PRECISION leta = eta->vector_buffer + threading->start_index[l->depth], lphi = phi->vector_buffer + threading->start_index[l->depth];
    
    complex_PRECISION b[6];
    while ( leta < eta_end ) {
      if(g.odd_even_table[i]==_EVEN){
        int i = 0;
        FOR6( b[i] =  (*lphi); lphi++; i++;   );
        FOR6( *leta = -(*lphi); lphi++; leta++; );
        i = 0;
        FOR6( *leta = - b[i] ; leta++; i++;   );
        i = 0;
        FOR6( b[i] =  (*lphi); lphi++; i++;   );
        FOR6( *leta =  (*lphi); lphi++; leta++; );
        i = 0;
        FOR6( *leta =   b[i] ; leta++; i++;   );
      } else if(g.odd_even_table[i]==_ODD){
        FOR24( *leta = _COMPLEX_PRECISION_ZERO; lphi++; leta++; );
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

/*
// not used anywhere
void scale_even_odd_PRECISION( vector_PRECISION *eta, vector_PRECISION *phi, complex_double even, complex_double odd, 
                               level_struct *l, struct Thread *threading ) {
   
  int i = threading->start_site[l->depth];
  buffer_PRECISION eta_end = eta->vector_buffer + threading->end_index[l->depth];
  buffer_PRECISION leta = eta->vector_buffer + threading->start_index[l->depth], lphi = phi->vector_buffer + threading->start_index[l->depth];

#ifdef HAVE_TM1p1
  if ( g.n_flavours == 2 )
    while ( leta < eta_end ) {
      if(g.odd_even_table[i]==_EVEN){
        FOR24( *leta = even*(*lphi); lphi++; leta++; );
      }
      else if(g.odd_even_table[i]==_ODD){
        FOR24( *leta = odd*(*lphi); lphi++; leta++; );
      }
      i++;
    }
  else
#endif
    while ( leta < eta_end ) {
      if(g.odd_even_table[i]==_EVEN) {
        FOR12( *leta = even*(*lphi); lphi++; leta++; );
      }
      else if(g.odd_even_table[i]==_ODD) {
        FOR12( *leta = odd*(*lphi); lphi++; leta++; );
      }
      i++;
    }
}

// not used anywhere
void set_even_to_zero_PRECISION( vector_PRECISION *eta, vector_PRECISION *phi, level_struct *l, struct Thread *threading ) {

  ASSERT(l->depth == 0);
  
  int i = threading->start_site[l->depth];
  buffer_PRECISION eta_end = eta->vector_buffer + threading->end_index[l->depth];
  buffer_PRECISION leta = eta->vector_buffer + threading->start_index[l->depth], lphi = phi->vector_buffer + threading->start_index[l->depth];

#ifdef HAVE_TM1p1
  if ( g.n_flavours == 2 )
    while ( leta < eta_end ) {
      if(g.odd_even_table[i]==_ODD) {
        FOR24( *leta = (*lphi); lphi++; leta++; );
      }
      else if(g.odd_even_table[i]==_EVEN) {
        FOR24( *leta = _COMPLEX_PRECISION_ZERO; lphi++; leta++; );
      }
      i++;
    }
  else
#endif
    while ( leta < eta_end ) {
      if(g.odd_even_table[i]==_ODD) {
        FOR12( *leta = (*lphi); lphi++; leta++; );
      }
      else if(g.odd_even_table[i]==_EVEN) {
        FOR12( *leta = _COMPLEX_PRECISION_ZERO; lphi++; leta++; );
      }
      i++;
    }
}

// not used anywhere
void set_odd_to_zero_PRECISION( vector_PRECISION *eta, vector_PRECISION *phi, level_struct *l, struct Thread *threading ) {
   
  int i = threading->start_site[l->depth];
  buffer_PRECISION eta_end = eta->vector_buffer + threading->end_index[l->depth];
  buffer_PRECISION leta = eta->vector_buffer + threading->start_index[l->depth], lphi = phi->vector_buffer + threading->start_index[l->depth];

#ifdef HAVE_TM1p1
  if ( g.n_flavours == 2 )
    while ( leta < eta_end ) {
      if(g.odd_even_table[i]==_EVEN){
        FOR24( *leta = (*lphi); lphi++; leta++; );
      }
      else if(g.odd_even_table[i]==_ODD){
        FOR24( *leta = 0; lphi++; leta++; );
      }
      i++;
    }
  else
#endif
    while ( leta < eta_end ) {
      if(g.odd_even_table[i]==_EVEN) {
        FOR12( *leta = (*lphi); lphi++; leta++; );
      }
      else if(g.odd_even_table[i]==_ODD) {
        FOR12( *leta = 0; lphi++; leta++; );
      }
      i++;
    }
}

*/
///// tanslations
// for TM1p1
void two_flavours_to_serial_PRECISION( vector_PRECISION *flav1, vector_PRECISION *flav2, vector_PRECISION *serial, level_struct *l, struct Thread *threading ) {

#ifdef HAVE_TM1p1

  /*
   * Order: spin0and1 of flav1
   *        spin0and1 of flav2
   *        spin2and3 of flav1
   *        spin2and3 of flav2
   */
  buffer_PRECISION serial_end;
  buffer_PRECISION serial_pt = serial->vector_buffer, flav1_pt = flav1->vector_buffer, flav2_pt = flav2->vector_buffer;
  
  if( g.n_flavours == 2 ) {
    serial_end = serial->vector_buffer + threading->end_index[l->depth];
    serial_pt += threading->start_index[l->depth];
    flav1_pt += threading->start_index[l->depth]/2;
    flav2_pt += threading->start_index[l->depth]/2;
  }
  else {
    serial_end = serial->vector_buffer + threading->end_index[l->depth]*2;
    serial_pt += threading->start_index[l->depth]*2;
    flav1_pt += threading->start_index[l->depth];
    flav2_pt += threading->start_index[l->depth];
  }

  while ( serial_pt < serial_end ) {
    FOR6( *serial_pt = (*flav1_pt); serial_pt++; flav1_pt++; )
    FOR6( *serial_pt = (*flav2_pt); serial_pt++; flav2_pt++; )
    FOR6( *serial_pt = (*flav1_pt); serial_pt++; flav1_pt++; )
    FOR6( *serial_pt = (*flav2_pt); serial_pt++; flav2_pt++; )
  }
#else
  START_MASTER(threading)
  warning0("two_flavours_to_serial_PRECISION called without HAVE_TM1p1 defined\n");
  END_MASTER(threading)
#endif
    
}

void serial_to_two_flavours_PRECISION( vector_PRECISION *flav1, vector_PRECISION *flav2, vector_PRECISION *serial, level_struct *l, struct Thread *threading ) {

#ifdef HAVE_TM1p1
  buffer_PRECISION serial_end;
  buffer_PRECISION serial_pt = serial->vector_buffer, flav1_pt = flav1->vector_buffer, flav2_pt = flav2->vector_buffer;

  if( g.n_flavours == 2 ) {
    serial_end = serial->vector_buffer + threading->end_index[l->depth];
    serial_pt += threading->start_index[l->depth];
    flav1_pt += threading->start_index[l->depth]/2;
    flav2_pt += threading->start_index[l->depth]/2;
  }
  else {
    serial_end = serial->vector_buffer + threading->end_index[l->depth]*2;
    serial_pt += threading->start_index[l->depth]*2;
    flav1_pt += threading->start_index[l->depth];
    flav2_pt += threading->start_index[l->depth];
  }

  while ( serial_pt < serial_end ) {
    FOR6( *flav1_pt = (*serial_pt); serial_pt++; flav1_pt++; )
    FOR6( *flav2_pt = (*serial_pt); serial_pt++; flav2_pt++; )
    FOR6( *flav1_pt = (*serial_pt); serial_pt++; flav1_pt++; )
    FOR6( *flav2_pt = (*serial_pt); serial_pt++; flav2_pt++; )
  }
#else
  START_MASTER(threading)
  warning0("two_flavours_to_serial_PRECISION called without HAVE_TM1p1 defined\n");
  END_MASTER(threading)
#endif
    
}

/******************  TEST ROUTINES *********************************************************/
// new function not impremented!!!!!!!
void two_flavours_test_PRECISION( operator_PRECISION_struct *op, level_struct *l, struct Thread *threading ) {

#ifdef HAVE_TM1p1
  double diff;
  
  vector_double vd[4], vdd[4];
  vector_PRECISION vpp[2];

  for(int i=0; i<4; i++){                                                                 
    vector_double_init( &vd[i] );                                                         
    vector_double_alloc( &vd[i], _INNER, 1, l, threading );
    vector_double_init( &vdd[i] );
    vector_double_alloc( &vdd[i], _INNER, 2, l, threading );                               
  }                                                                                       
                                                                                          
  for(int i=0; i<2; i++){                                                                 
    vector_PRECISION_init( &vpp[i] );                                                      
    vector_PRECISION_alloc( &vpp[i], _INNER, 2, l, threading );                            
  }

  ASSERT(g.n_flavours==2);

  data_layout_n_flavours( 1, l, threading );

  START_LOCKED_MASTER(threading)
  vector_double_define_random( &vd[0], 0, l->inner_vector_size, l );
  vector_double_define_random( &vd[1], 0, l->inner_vector_size, l );
  apply_operator_double( &vd[2], &vd[0], &(g.p), l, no_threading );
#ifdef HAVE_TM
  buffer_double_real_scale( g.op_double.tm_term, g.op_double.tm_term, -1, 0, l->inner_vector_size, l ); 
#endif
  apply_operator_double( &vd[3], &vd[1], &(g.p), l, no_threading );
#ifdef HAVE_TM
  buffer_double_real_scale( g.op_double.tm_term, g.op_double.tm_term, -1, 0, l->inner_vector_size, l ); 
#endif
  add_diagonal_double( &vd[2], &vd[1], g.op_double.epsbar_term, l->inner_vector_size );
  add_diagonal_double( &vd[3], &vd[0], g.op_double.epsbar_term, l->inner_vector_size );

  two_flavours_to_serial_double( &vd[0], &vd[1], &vdd[0], l, no_threading );
  two_flavours_to_serial_double( &vd[2], &vd[3], &vdd[1], l, no_threading );
  END_LOCKED_MASTER(threading)

  data_layout_n_flavours( 2, l, threading );

  START_LOCKED_MASTER(threading)
  trans_PRECISION( &vpp[0], &vdd[0], op->translation_table, l, no_threading );
  apply_operator_PRECISION( &vpp[1], &vpp[0], &(l->p_PRECISION), l, no_threading );
  trans_back_PRECISION( &vdd[2], &vpp[1], op->translation_table, l, no_threading );
  
  vector_double_minus( &vdd[3], &vdd[2], &vdd[1], 0, l->inner_vector_size, l );
  diff = global_norm_double( &vdd[3], 0, l->inner_vector_size, l, no_threading ) /
    global_norm_double( &vdd[2], 0, l->inner_vector_size, l, no_threading );
  
  test0_PRECISION("depth: %d, correctness of doublet Dirac operator PRECISION: %le\n", l->depth, diff );
  END_LOCKED_MASTER(threading)

  if(threading->n_core > 1) {
    trans_PRECISION( &vpp[0], &vdd[0], op->translation_table, l, threading );
    apply_operator_PRECISION( &vpp[1], &vpp[0], &(l->p_PRECISION), l, threading );
    trans_back_PRECISION( &vdd[2], &vpp[1], op->translation_table, l, threading );
    
    SYNC_MASTER_TO_ALL(threading)
    SYNC_CORES(threading)

    START_LOCKED_MASTER(threading)
    vector_double_minus( &vdd[3], &vdd[2], &vdd[1], 0, l->inner_vector_size, l );
    diff = global_norm_double( &vdd[3], 0, l->inner_vector_size, l, no_threading ) /
      global_norm_double( &vdd[2], 0, l->inner_vector_size, l, no_threading );
    
    test0_PRECISION("depth: %d, correctness of doublet Dirac operator PRECISION with threading: %le\n", l->depth, diff );
    END_LOCKED_MASTER(threading)
  }    
  
  for(int i=0; i<4; i++){
    vector_double_free( &vd[i], l, threading );
    vector_double_free( &vdd[i], l, threading );
  }

  for(int i=0; i<2; i++)
    vector_PRECISION_free( &vpp[i], l, threading );

  START_LOCKED_MASTER(threading)
  if ( g.method >=4 && g.odd_even )
    oddeven_PRECISION_test( l );
  END_LOCKED_MASTER(threading) 
#endif
    
}
