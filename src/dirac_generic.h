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
 * checked:
`* 1st cleanup:12/18/2019
 */

#ifndef DIRAC_PRECISION_HEADER
  #define DIRAC_PRECISION_HEADER

  struct Thread;
  

  void clover_PRECISION_new( vector_PRECISION *eta, vector_PRECISION *phi, operator_PRECISION_struct *op, int start, int end, level_struct *l, struct Thread *threading );
  void d_plus_clover_PRECISION_new( vector_PRECISION *eta, vector_PRECISION *phi, operator_PRECISION_struct *op, level_struct *l, struct Thread *threading );
  void block_d_plus_clover_PRECISION_new( vector_PRECISION *eta, vector_PRECISION *phi, int start, schwarz_PRECISION_struct *s, level_struct *l, struct Thread *threading );

  void diagonal_aggregate_PRECISION_new( vector_PRECISION *eta1, vector_PRECISION *eta2, vector_PRECISION *phi, config_PRECISION diag, level_struct *l );
  void d_plus_clover_aggregate_PRECISION_new( vector_PRECISION *eta1, vector_PRECISION *eta2, vector_PRECISION *phi, schwarz_PRECISION_struct *s, level_struct *l );
  void d_neighbor_aggregate_PRECISION_new( vector_PRECISION *eta1, vector_PRECISION *eta2, vector_PRECISION *phi, const int mu, schwarz_PRECISION_struct *s, level_struct *l );

  void apply_twisted_bc_to_vector_PRECISION_new( vector_PRECISION *eta, vector_PRECISION *phi, double *theta, level_struct *l);

  void operator_updates_PRECISION( level_struct *l, struct Thread *threading );
  void m0_update_PRECISION( PRECISION m0,operator_PRECISION_struct *op, level_struct *l, struct Thread *threading );
  void tm_term_PRECISION_setup( PRECISION mu, PRECISION even, PRECISION odd, operator_PRECISION_struct *op,
                              level_struct *l, struct Thread *threading );
  void epsbar_term_PRECISION_setup( PRECISION epsbar, PRECISION even, PRECISION odd, operator_PRECISION_struct *op,
                                    level_struct *l, struct Thread *threading );

  void gamma5_PRECISION_new( vector_PRECISION *eta, vector_PRECISION *phi, level_struct *l, struct Thread *threading );

  void gamma5_set_even_to_zero_PRECISION( vector_PRECISION *eta, vector_PRECISION *phi, level_struct *l, struct Thread *threading );//used in DDalphaAMG_interface 
  void gamma5_set_odd_to_zero_PRECISION( vector_PRECISION *eta, vector_PRECISION *phi, level_struct *l, struct Thread *threading );//used in DDalphaAMG_interface

  void tau1_gamma5_PRECISION( vector_PRECISION *eta, vector_PRECISION *phi, level_struct *l, struct Thread *threading );//used when TM1p1
//  void set_even_to_zero_PRECISION( vector_PRECISION *eta, vector_PRECISION *phi, level_struct *l, struct Thread *threading );
//  void set_odd_to_zero_PRECISION( vector_PRECISION *eta, vector_PRECISION *phi, level_struct *l, struct Thread *threading );

  void tau1_gamma5_set_even_to_zero_PRECISION( vector_PRECISION *eta, vector_PRECISION *phi, level_struct *l, struct Thread *threading );//used when TM1p1 in DDalphaAMG_interface 
  void tau1_gamma5_set_odd_to_zero_PRECISION( vector_PRECISION *eta, vector_PRECISION *phi, level_struct *l, struct Thread *threading ); //used when TM1p1 in DDalphaAMG_interface
//  void scale_even_odd_PRECISION( vector_PRECISION *eta, vector_PRECISION *phi, complex_double even, complex_double odd,
//                                 level_struct *l, struct Thread *threading );

  void two_flavours_to_serial_PRECISION( vector_PRECISION *flav1, vector_PRECISION *flav2, vector_PRECISION *serial, level_struct *l, struct Thread *threading );
  void serial_to_two_flavours_PRECISION( vector_PRECISION *flav1, vector_PRECISION *flav2, vector_PRECISION *serial, level_struct *l, struct Thread *threading );

  void two_flavours_test_PRECISION( operator_PRECISION_struct *op, level_struct *l, struct Thread *threading );


//eta+= diag*phi
  static inline void add_diagonal_PRECISION_new( const vector_PRECISION *eta, const vector_PRECISION *phi,
             const config_PRECISION diag, const int length ) {
    int i, j, jj, n_vect_phi = phi->num_vect_now, n_vect_eta = eta->num_vect_now;
    config_PRECISION diag_pt = diag;
    buffer_PRECISION phi_pt = phi->vector_buffer, eta_pt = eta->vector_buffer, eta_end = eta->vector_buffer + length*eta->num_vect;
    //        printf("add_diagonal_PRECISION:%d %d %d %d\n",n_vect_phi,n_vect_eta ,phi->num_vect,eta->num_vect);
    if ( n_vect_phi != n_vect_eta )
      error0("add_diagonal_PRECISION: assumptions are not met\n");
/*#ifdef HAVE_TM1p1
    if(g.n_flavours == 2)
      while ( eta_pt < eta_end ) {
        FOR6( *eta_pt += (*phi_pt)*(*diag_pt); eta_pt++; phi_pt++; diag_pt++; )
        diag_pt -= 6;
        FOR6( *eta_pt -= (*phi_pt)*(*diag_pt); eta_pt++; phi_pt++; diag_pt++; )
        FOR6( *eta_pt += (*phi_pt)*(*diag_pt); eta_pt++; phi_pt++; diag_pt++; )
        diag_pt -= 6;
        FOR6( *eta_pt -= (*phi_pt)*(*diag_pt); eta_pt++; phi_pt++; diag_pt++; )
      }
       while ( eta_pt < eta_end )
	 FOR6( VECTOR_LOOP(j, n_vect, jj, *eta_pt += (*phi_pt)*(*diag_pt); eta_pt++; phi_pt++;) 
	       eta_pt += eta->num_vect - n_vect_eta; phi_pt += phi->num_vect -n_vect_phi;
	       diag_pt++;)
	   diag -= 6;
       ...
     else
#endif*/
       while ( eta_pt < eta_end )
         for ( i=0; i<12; i++ ) {
           VECTOR_LOOP(j, n_vect_eta, jj, *eta_pt += (*phi_pt)*(*diag_pt); eta_pt++; phi_pt++;)//printf("%g ",creal_PRECISION(*eta_pt ));)
	   eta_pt += eta->num_vect-n_vect_eta; phi_pt += phi->num_vect -n_vect_eta;
	   //	   	   printf("%g ",creal_PRECISION((*diag_pt)));
           diag_pt++;
         }
  }

#ifdef HAVE_TM1p1
  static inline void apply_doublet_coupling_PRECISION( const vector_PRECISION *eta, const vector_PRECISION *phi,
             const config_PRECISION diag, const int length ) {
    config_PRECISION diag_pt = diag;
    buffer_PRECISION phi_pt = phi->vector_buffer, eta_pt = eta->vector_buffer, eta_end = eta->vector_buffer + length;
    while ( eta_pt < eta_end ) { 
      phi_pt += 6;
      FOR6( *eta_pt += (*phi_pt)*(*diag_pt); eta_pt++; phi_pt++; diag_pt++; )
      phi_pt -= 12;
      diag_pt -= 6;
      FOR6( *eta_pt += (*phi_pt)*(*diag_pt); eta_pt++; phi_pt++; diag_pt++; )
      phi_pt += 6;
    }
  }
#endif

  // eta = D*phi
  static inline void mvm_PRECISION_new( const buffer_PRECISION eta, const complex_PRECISION *D, const buffer_PRECISION phi, const int n_vect, const int n_vect_eta, const int n_vect_phi ) {
    int j, jj;
    VECTOR_LOOP(j, n_vect, jj, eta[0*n_vect_eta+j+jj]  = D[0]*phi[0*n_vect_phi+j+jj];
                               eta[0*n_vect_eta+j+jj] += D[1]*phi[1*n_vect_phi+j+jj];
                               eta[0*n_vect_eta+j+jj] += D[2]*phi[2*n_vect_phi+j+jj];
                               eta[1*n_vect_eta+j+jj]  = D[3]*phi[0*n_vect_phi+j+jj];
                               eta[1*n_vect_eta+j+jj] += D[4]*phi[1*n_vect_phi+j+jj];
                               eta[1*n_vect_eta+j+jj] += D[5]*phi[2*n_vect_phi+j+jj];
                               eta[2*n_vect_eta+j+jj]  = D[6]*phi[0*n_vect_phi+j+jj];
                               eta[2*n_vect_eta+j+jj] += D[7]*phi[1*n_vect_phi+j+jj];
                               eta[2*n_vect_eta+j+jj] += D[8]*phi[2*n_vect_phi+j+jj];)
  }

  // eta = D**H*phi
  static inline void mvmh_PRECISION_new( const buffer_PRECISION eta, const complex_PRECISION *D, const buffer_PRECISION phi, const int n_vect, const int n_vect_eta, const int n_vect_phi ) {
    int j, jj;
    //  printf("d_plus_clover_aggrmvmv_PRECISION: %d %d %d %d\n",n_vect, n_vect_eta, n_vect_phi );
    VECTOR_LOOP(j, n_vect, jj, eta[0*n_vect_eta+j+jj]  = conj_PRECISION(D[0])*phi[0*n_vect_phi+j+jj];
                               eta[1*n_vect_eta+j+jj]  = conj_PRECISION(D[1])*phi[0*n_vect_phi+j+jj];
                               eta[2*n_vect_eta+j+jj]  = conj_PRECISION(D[2])*phi[0*n_vect_phi+j+jj];
                               eta[0*n_vect_eta+j+jj] += conj_PRECISION(D[3])*phi[1*n_vect_phi+j+jj];
                               eta[1*n_vect_eta+j+jj] += conj_PRECISION(D[4])*phi[1*n_vect_phi+j+jj];
                               eta[2*n_vect_eta+j+jj] += conj_PRECISION(D[5])*phi[1*n_vect_phi+j+jj];
                               eta[0*n_vect_eta+j+jj] += conj_PRECISION(D[6])*phi[2*n_vect_phi+j+jj];
                               eta[1*n_vect_eta+j+jj] += conj_PRECISION(D[7])*phi[2*n_vect_phi+j+jj];
                               eta[2*n_vect_eta+j+jj] += conj_PRECISION(D[8])*phi[2*n_vect_phi+j+jj];)
  }
  
  // eta = -D*phi
  static inline void nmvm_PRECISION_new( const buffer_PRECISION eta, const complex_PRECISION *D, const buffer_PRECISION phi, const int n_vect, const int n_vect_eta, const int n_vect_phi ) {
    int j, jj;
    VECTOR_LOOP(j, n_vect, jj, eta[0*n_vect_eta+j+jj]  = - D[0]*phi[0*n_vect_phi+j+jj];
                               eta[0*n_vect_eta+j+jj] -=   D[1]*phi[1*n_vect_phi+j+jj];
                               eta[0*n_vect_eta+j+jj] -=   D[2]*phi[2*n_vect_phi+j+jj];
                               eta[1*n_vect_eta+j+jj]  = - D[3]*phi[0*n_vect_phi+j+jj];
                               eta[1*n_vect_eta+j+jj] -=   D[4]*phi[1*n_vect_phi+j+jj];
                               eta[1*n_vect_eta+j+jj] -=   D[5]*phi[2*n_vect_phi+j+jj];
                               eta[2*n_vect_eta+j+jj]  = - D[6]*phi[0*n_vect_phi+j+jj];
                               eta[2*n_vect_eta+j+jj] -=   D[7]*phi[1*n_vect_phi+j+jj];
                               eta[2*n_vect_eta+j+jj] -=   D[8]*phi[2*n_vect_phi+j+jj];)
  }

  // eta = -D**H*phi
  static inline void nmvmh_PRECISION_new( const buffer_PRECISION eta, const complex_PRECISION *D, const buffer_PRECISION phi, const int n_vect, const int n_vect_eta, const int n_vect_phi ) {
    int j, jj;
    VECTOR_LOOP(j, n_vect, jj, eta[0*n_vect_eta+j+jj]  = - conj_PRECISION(D[0])*phi[0*n_vect_phi+j+jj];
                               eta[1*n_vect_eta+j+jj]  = - conj_PRECISION(D[1])*phi[0*n_vect_phi+j+jj];
                               eta[2*n_vect_eta+j+jj]  = - conj_PRECISION(D[2])*phi[0*n_vect_phi+j+jj];
                               eta[0*n_vect_eta+j+jj] -=   conj_PRECISION(D[3])*phi[1*n_vect_phi+j+jj];
                               eta[1*n_vect_eta+j+jj] -=   conj_PRECISION(D[4])*phi[1*n_vect_phi+j+jj];
                               eta[2*n_vect_eta+j+jj] -=   conj_PRECISION(D[5])*phi[1*n_vect_phi+j+jj];
                               eta[0*n_vect_eta+j+jj] -=   conj_PRECISION(D[6])*phi[2*n_vect_phi+j+jj];
                               eta[1*n_vect_eta+j+jj] -=   conj_PRECISION(D[7])*phi[2*n_vect_phi+j+jj];
                               eta[2*n_vect_eta+j+jj] -=   conj_PRECISION(D[8])*phi[2*n_vect_phi+j+jj];)
  }

  // 1 - gamma_T
  static inline void prp_T_PRECISION_new( const buffer_PRECISION prp_pt, const buffer_PRECISION l_pt, const int n_vect, const int n_vect_prp, const int n_vect_l  ) {
    int j, jj;
    VECTOR_LOOP(j, n_vect, jj, prp_pt[0*n_vect_prp+j+jj] = l_pt[0*n_vect_l+j+jj] -GAMMA_T_SPIN0_VAL*l_pt[ 3*GAMMA_T_SPIN0_CO   *n_vect_l+j+jj];
                               prp_pt[1*n_vect_prp+j+jj] = l_pt[1*n_vect_l+j+jj] -GAMMA_T_SPIN0_VAL*l_pt[(3*GAMMA_T_SPIN0_CO+1)*n_vect_l+j+jj];
                               prp_pt[2*n_vect_prp+j+jj] = l_pt[2*n_vect_l+j+jj] -GAMMA_T_SPIN0_VAL*l_pt[(3*GAMMA_T_SPIN0_CO+2)*n_vect_l+j+jj];
                               prp_pt[3*n_vect_prp+j+jj] = l_pt[3*n_vect_l+j+jj] -GAMMA_T_SPIN1_VAL*l_pt[ 3*GAMMA_T_SPIN1_CO   *n_vect_l+j+jj];
                               prp_pt[4*n_vect_prp+j+jj] = l_pt[4*n_vect_l+j+jj] -GAMMA_T_SPIN1_VAL*l_pt[(3*GAMMA_T_SPIN1_CO+1)*n_vect_l+j+jj];
                               prp_pt[5*n_vect_prp+j+jj] = l_pt[5*n_vect_l+j+jj] -GAMMA_T_SPIN1_VAL*l_pt[(3*GAMMA_T_SPIN1_CO+2)*n_vect_l+j+jj];)
  }

  // 1 + gamma_T
  static inline void prn_T_PRECISION_new( const buffer_PRECISION prn_pt, const buffer_PRECISION l_pt, const int n_vect, const int n_vect_prn, const int n_vect_l ) {//start?
    int j, jj;
    VECTOR_LOOP(j, n_vect, jj, prn_pt[0*n_vect_prn+j+jj] = l_pt[0*n_vect_l+j+jj] +GAMMA_T_SPIN0_VAL*l_pt[ 3*GAMMA_T_SPIN0_CO   *n_vect_l+j+jj];
                               prn_pt[1*n_vect_prn+j+jj] = l_pt[1*n_vect_l+j+jj] +GAMMA_T_SPIN0_VAL*l_pt[(3*GAMMA_T_SPIN0_CO+1)*n_vect_l+j+jj];
                               prn_pt[2*n_vect_prn+j+jj] = l_pt[2*n_vect_l+j+jj] +GAMMA_T_SPIN0_VAL*l_pt[(3*GAMMA_T_SPIN0_CO+2)*n_vect_l+j+jj];
                               prn_pt[3*n_vect_prn+j+jj] = l_pt[3*n_vect_l+j+jj] +GAMMA_T_SPIN1_VAL*l_pt[ 3*GAMMA_T_SPIN1_CO   *n_vect_l+j+jj];
                               prn_pt[4*n_vect_prn+j+jj] = l_pt[4*n_vect_l+j+jj] +GAMMA_T_SPIN1_VAL*l_pt[(3*GAMMA_T_SPIN1_CO+1)*n_vect_l+j+jj];
                               prn_pt[5*n_vect_prn+j+jj] = l_pt[5*n_vect_l+j+jj] +GAMMA_T_SPIN1_VAL*l_pt[(3*GAMMA_T_SPIN1_CO+2)*n_vect_l+j+jj];)
  }

  // - (1 - gamma_T)
  static inline void pbp_su3_T_PRECISION_new( const buffer_PRECISION prp_su3_pt, const buffer_PRECISION l_pt, const int n_vect, const int n_vect_pbp, const int n_vect_l  ) {
    int j, jj;
    VECTOR_LOOP(j, n_vect, jj, l_pt[ 0*n_vect_l+j+jj] -= prp_su3_pt[0*n_vect_pbp+j+jj];
                               l_pt[ 1*n_vect_l+j+jj] -= prp_su3_pt[1*n_vect_pbp+j+jj];
                               l_pt[ 2*n_vect_l+j+jj] -= prp_su3_pt[2*n_vect_pbp+j+jj];
                               l_pt[ 3*n_vect_l+j+jj] -= prp_su3_pt[3*n_vect_pbp+j+jj];
                               l_pt[ 4*n_vect_l+j+jj] -= prp_su3_pt[4*n_vect_pbp+j+jj];
                               l_pt[ 5*n_vect_l+j+jj] -= prp_su3_pt[5*n_vect_pbp+j+jj];
                               l_pt[ 6*n_vect_l+j+jj] += GAMMA_T_SPIN2_VAL*prp_su3_pt[ 3*GAMMA_T_SPIN2_CO   *n_vect_pbp+j+jj];
                               l_pt[ 7*n_vect_l+j+jj] += GAMMA_T_SPIN2_VAL*prp_su3_pt[(3*GAMMA_T_SPIN2_CO+1)*n_vect_pbp+j+jj];
                               l_pt[ 8*n_vect_l+j+jj] += GAMMA_T_SPIN2_VAL*prp_su3_pt[(3*GAMMA_T_SPIN2_CO+2)*n_vect_pbp+j+jj];
                               l_pt[ 9*n_vect_l+j+jj] += GAMMA_T_SPIN3_VAL*prp_su3_pt[ 3*GAMMA_T_SPIN3_CO   *n_vect_pbp+j+jj];
                               l_pt[10*n_vect_l+j+jj] += GAMMA_T_SPIN3_VAL*prp_su3_pt[(3*GAMMA_T_SPIN3_CO+1)*n_vect_pbp+j+jj];
                               l_pt[11*n_vect_l+j+jj] += GAMMA_T_SPIN3_VAL*prp_su3_pt[(3*GAMMA_T_SPIN3_CO+2)*n_vect_pbp+j+jj];)
  }

  // -(1 + gamma_T)
  static inline void pbn_su3_T_PRECISION_new( const buffer_PRECISION prn_su3_pt, const buffer_PRECISION l_pt, const int n_vect, const int n_vect_prn, const int n_vect_l  ) {
    int j, jj;
    VECTOR_LOOP(j, n_vect, jj, l_pt[ 0*n_vect_l+j+jj] -= prn_su3_pt[0*n_vect_prn+j+jj];
                               l_pt[ 1*n_vect_l+j+jj] -= prn_su3_pt[1*n_vect_prn+j+jj];
                               l_pt[ 2*n_vect_l+j+jj] -= prn_su3_pt[2*n_vect_prn+j+jj];
                               l_pt[ 3*n_vect_l+j+jj] -= prn_su3_pt[3*n_vect_prn+j+jj];
                               l_pt[ 4*n_vect_l+j+jj] -= prn_su3_pt[4*n_vect_prn+j+jj];
                               l_pt[ 5*n_vect_l+j+jj] -= prn_su3_pt[5*n_vect_prn+j+jj];
                               l_pt[ 6*n_vect_l+j+jj] -= GAMMA_T_SPIN2_VAL*prn_su3_pt[ 3*GAMMA_T_SPIN2_CO   *n_vect_prn+j+jj];
                               l_pt[ 7*n_vect_l+j+jj] -= GAMMA_T_SPIN2_VAL*prn_su3_pt[(3*GAMMA_T_SPIN2_CO+1)*n_vect_prn+j+jj];
                               l_pt[ 8*n_vect_l+j+jj] -= GAMMA_T_SPIN2_VAL*prn_su3_pt[(3*GAMMA_T_SPIN2_CO+2)*n_vect_prn+j+jj];
                               l_pt[ 9*n_vect_l+j+jj] -= GAMMA_T_SPIN3_VAL*prn_su3_pt[ 3*GAMMA_T_SPIN3_CO   *n_vect_prn+j+jj];
                               l_pt[10*n_vect_l+j+jj] -= GAMMA_T_SPIN3_VAL*prn_su3_pt[(3*GAMMA_T_SPIN3_CO+1)*n_vect_prn+j+jj];
                               l_pt[11*n_vect_l+j+jj] -= GAMMA_T_SPIN3_VAL*prn_su3_pt[(3*GAMMA_T_SPIN3_CO+2)*n_vect_prn+j+jj];)
  }

  // 1 - gamma_Z
  static inline void prp_Z_PRECISION_new( const buffer_PRECISION prp_pt, const buffer_PRECISION l_pt, const int n_vect, const int n_vect_prp, const int n_vect_l  ) {
    int j, jj;
    VECTOR_LOOP(j, n_vect, jj, prp_pt[0*n_vect_prp+j+jj] = l_pt[0*n_vect_l+j+jj] -GAMMA_Z_SPIN0_VAL*l_pt[ 3*GAMMA_Z_SPIN0_CO   *n_vect_l+j+jj];
                               prp_pt[1*n_vect_prp+j+jj] = l_pt[1*n_vect_l+j+jj] -GAMMA_Z_SPIN0_VAL*l_pt[(3*GAMMA_Z_SPIN0_CO+1)*n_vect_l+j+jj];
                               prp_pt[2*n_vect_prp+j+jj] = l_pt[2*n_vect_l+j+jj] -GAMMA_Z_SPIN0_VAL*l_pt[(3*GAMMA_Z_SPIN0_CO+2)*n_vect_l+j+jj];
                               prp_pt[3*n_vect_prp+j+jj] = l_pt[3*n_vect_l+j+jj] -GAMMA_Z_SPIN1_VAL*l_pt[ 3*GAMMA_Z_SPIN1_CO   *n_vect_l+j+jj];
                               prp_pt[4*n_vect_prp+j+jj] = l_pt[4*n_vect_l+j+jj] -GAMMA_Z_SPIN1_VAL*l_pt[(3*GAMMA_Z_SPIN1_CO+1)*n_vect_l+j+jj];
                               prp_pt[5*n_vect_prp+j+jj] = l_pt[5*n_vect_l+j+jj] -GAMMA_Z_SPIN1_VAL*l_pt[(3*GAMMA_Z_SPIN1_CO+2)*n_vect_l+j+jj];)
  }

  // 1 + gamma_Z
  static inline void prn_Z_PRECISION_new( const buffer_PRECISION prn_pt, const buffer_PRECISION l_pt, const int n_vect, const int n_vect_prn, const int n_vect_l  ) {
    int j, jj;
    VECTOR_LOOP(j, n_vect, jj, prn_pt[0*n_vect_prn+j+jj] = l_pt[0*n_vect_l+j+jj] +GAMMA_Z_SPIN0_VAL*l_pt[ 3*GAMMA_Z_SPIN0_CO   *n_vect_l+j+jj];
                               prn_pt[1*n_vect_prn+j+jj] = l_pt[1*n_vect_l+j+jj] +GAMMA_Z_SPIN0_VAL*l_pt[(3*GAMMA_Z_SPIN0_CO+1)*n_vect_l+j+jj];
                               prn_pt[2*n_vect_prn+j+jj] = l_pt[2*n_vect_l+j+jj] +GAMMA_Z_SPIN0_VAL*l_pt[(3*GAMMA_Z_SPIN0_CO+2)*n_vect_l+j+jj];
                               prn_pt[3*n_vect_prn+j+jj] = l_pt[3*n_vect_l+j+jj] +GAMMA_Z_SPIN1_VAL*l_pt[ 3*GAMMA_Z_SPIN1_CO   *n_vect_l+j+jj];
                               prn_pt[4*n_vect_prn+j+jj] = l_pt[4*n_vect_l+j+jj] +GAMMA_Z_SPIN1_VAL*l_pt[(3*GAMMA_Z_SPIN1_CO+1)*n_vect_l+j+jj];
                               prn_pt[5*n_vect_prn+j+jj] = l_pt[5*n_vect_l+j+jj] +GAMMA_Z_SPIN1_VAL*l_pt[(3*GAMMA_Z_SPIN1_CO+2)*n_vect_l+j+jj];)
  }

  // - (1 - gamma_Z)
  static inline void pbp_su3_Z_PRECISION_new( const buffer_PRECISION prp_su3_pt, const buffer_PRECISION l_pt, const int n_vect, const int n_vect_pbp, const int n_vect_l  ) {
    int j, jj;
    VECTOR_LOOP(j, n_vect, jj, l_pt[ 0*n_vect_l+j+jj] -= prp_su3_pt[0*n_vect_pbp+j+jj];
                               l_pt[ 1*n_vect_l+j+jj] -= prp_su3_pt[1*n_vect_pbp+j+jj];
                               l_pt[ 2*n_vect_l+j+jj] -= prp_su3_pt[2*n_vect_pbp+j+jj];
                               l_pt[ 3*n_vect_l+j+jj] -= prp_su3_pt[3*n_vect_pbp+j+jj];
                               l_pt[ 4*n_vect_l+j+jj] -= prp_su3_pt[4*n_vect_pbp+j+jj];
                               l_pt[ 5*n_vect_l+j+jj] -= prp_su3_pt[5*n_vect_pbp+j+jj];
                               l_pt[ 6*n_vect_l+j+jj] += GAMMA_Z_SPIN2_VAL*prp_su3_pt[ 3*GAMMA_Z_SPIN2_CO   *n_vect_pbp+j+jj];
                               l_pt[ 7*n_vect_l+j+jj] += GAMMA_Z_SPIN2_VAL*prp_su3_pt[(3*GAMMA_Z_SPIN2_CO+1)*n_vect_pbp+j+jj];
                               l_pt[ 8*n_vect_l+j+jj] += GAMMA_Z_SPIN2_VAL*prp_su3_pt[(3*GAMMA_Z_SPIN2_CO+2)*n_vect_pbp+j+jj];
                               l_pt[ 9*n_vect_l+j+jj] += GAMMA_Z_SPIN3_VAL*prp_su3_pt[ 3*GAMMA_Z_SPIN3_CO   *n_vect_pbp+j+jj];
                               l_pt[10*n_vect_l+j+jj] += GAMMA_Z_SPIN3_VAL*prp_su3_pt[(3*GAMMA_Z_SPIN3_CO+1)*n_vect_pbp+j+jj];
                               l_pt[11*n_vect_l+j+jj] += GAMMA_Z_SPIN3_VAL*prp_su3_pt[(3*GAMMA_Z_SPIN3_CO+2)*n_vect_pbp+j+jj];)
  }

  // -(1 + gamma_Z)
  static inline void pbn_su3_Z_PRECISION_new( const buffer_PRECISION pbn_su3_pt, const buffer_PRECISION l_pt, const int n_vect, const int n_vect_pbn, const int n_vect_l  ) {
    int j, jj;
    VECTOR_LOOP(j, n_vect, jj, l_pt[ 0*n_vect_l+j+jj] -= pbn_su3_pt[0*n_vect_pbn+j+jj];
                               l_pt[ 1*n_vect_l+j+jj] -= pbn_su3_pt[1*n_vect_pbn+j+jj];
                               l_pt[ 2*n_vect_l+j+jj] -= pbn_su3_pt[2*n_vect_pbn+j+jj];
                               l_pt[ 3*n_vect_l+j+jj] -= pbn_su3_pt[3*n_vect_pbn+j+jj];
                               l_pt[ 4*n_vect_l+j+jj] -= pbn_su3_pt[4*n_vect_pbn+j+jj];
                               l_pt[ 5*n_vect_l+j+jj] -= pbn_su3_pt[5*n_vect_pbn+j+jj];
                               l_pt[ 6*n_vect_l+j+jj] -= GAMMA_Z_SPIN2_VAL*pbn_su3_pt[ 3*GAMMA_Z_SPIN2_CO   *n_vect_pbn+j+jj];
                               l_pt[ 7*n_vect_l+j+jj] -= GAMMA_Z_SPIN2_VAL*pbn_su3_pt[(3*GAMMA_Z_SPIN2_CO+1)*n_vect_pbn+j+jj];
                               l_pt[ 8*n_vect_l+j+jj] -= GAMMA_Z_SPIN2_VAL*pbn_su3_pt[(3*GAMMA_Z_SPIN2_CO+2)*n_vect_pbn+j+jj];
                               l_pt[ 9*n_vect_l+j+jj] -= GAMMA_Z_SPIN3_VAL*pbn_su3_pt[ 3*GAMMA_Z_SPIN3_CO   *n_vect_pbn+j+jj];
                               l_pt[10*n_vect_l+j+jj] -= GAMMA_Z_SPIN3_VAL*pbn_su3_pt[(3*GAMMA_Z_SPIN3_CO+1)*n_vect_pbn+j+jj];
                               l_pt[11*n_vect_l+j+jj] -= GAMMA_Z_SPIN3_VAL*pbn_su3_pt[(3*GAMMA_Z_SPIN3_CO+2)*n_vect_pbn+j+jj];)
  }

  // 1 - gamma_Y
  static inline void prp_Y_PRECISION_new( const buffer_PRECISION prp_pt, const buffer_PRECISION l_pt, const int n_vect, const int n_vect_prp, const int n_vect_l  ) {
    int j, jj;
    VECTOR_LOOP(j, n_vect, jj, prp_pt[0*n_vect_prp+j+jj] = l_pt[0*n_vect_l+j+jj] -GAMMA_Y_SPIN0_VAL*l_pt[ 3*GAMMA_Y_SPIN0_CO   *n_vect_l+j+jj];
                               prp_pt[1*n_vect_prp+j+jj] = l_pt[1*n_vect_l+j+jj] -GAMMA_Y_SPIN0_VAL*l_pt[(3*GAMMA_Y_SPIN0_CO+1)*n_vect_l+j+jj];
                               prp_pt[2*n_vect_prp+j+jj] = l_pt[2*n_vect_l+j+jj] -GAMMA_Y_SPIN0_VAL*l_pt[(3*GAMMA_Y_SPIN0_CO+2)*n_vect_l+j+jj];
                               prp_pt[3*n_vect_prp+j+jj] = l_pt[3*n_vect_l+j+jj] -GAMMA_Y_SPIN1_VAL*l_pt[ 3*GAMMA_Y_SPIN1_CO   *n_vect_l+j+jj];
                               prp_pt[4*n_vect_prp+j+jj] = l_pt[4*n_vect_l+j+jj] -GAMMA_Y_SPIN1_VAL*l_pt[(3*GAMMA_Y_SPIN1_CO+1)*n_vect_l+j+jj];
                               prp_pt[5*n_vect_prp+j+jj] = l_pt[5*n_vect_l+j+jj] -GAMMA_Y_SPIN1_VAL*l_pt[(3*GAMMA_Y_SPIN1_CO+2)*n_vect_l+j+jj];)
  }

  // 1 + gamma_Y
  static inline void prn_Y_PRECISION_new( const buffer_PRECISION prn_pt, const buffer_PRECISION l_pt, const int n_vect, const int n_vect_prn, const int n_vect_l  ) {
    int j, jj;
    VECTOR_LOOP(j, n_vect, jj, prn_pt[0*n_vect_prn+j+jj] = l_pt[0*n_vect_l+j+jj] +GAMMA_Y_SPIN0_VAL*l_pt[ 3*GAMMA_Y_SPIN0_CO   *n_vect_l+j+jj];
                               prn_pt[1*n_vect_prn+j+jj] = l_pt[1*n_vect_l+j+jj] +GAMMA_Y_SPIN0_VAL*l_pt[(3*GAMMA_Y_SPIN0_CO+1)*n_vect_l+j+jj];
                               prn_pt[2*n_vect_prn+j+jj] = l_pt[2*n_vect_l+j+jj] +GAMMA_Y_SPIN0_VAL*l_pt[(3*GAMMA_Y_SPIN0_CO+2)*n_vect_l+j+jj];
                               prn_pt[3*n_vect_prn+j+jj] = l_pt[3*n_vect_l+j+jj] +GAMMA_Y_SPIN1_VAL*l_pt[ 3*GAMMA_Y_SPIN1_CO   *n_vect_l+j+jj];
                               prn_pt[4*n_vect_prn+j+jj] = l_pt[4*n_vect_l+j+jj] +GAMMA_Y_SPIN1_VAL*l_pt[(3*GAMMA_Y_SPIN1_CO+1)*n_vect_l+j+jj];
                               prn_pt[5*n_vect_prn+j+jj] = l_pt[5*n_vect_l+j+jj] +GAMMA_Y_SPIN1_VAL*l_pt[(3*GAMMA_Y_SPIN1_CO+2)*n_vect_l+j+jj];)
  }

  // - (1 - gamma_Y)
  static inline void pbp_su3_Y_PRECISION_new( const buffer_PRECISION prp_su3_pt, const buffer_PRECISION l_pt, const int n_vect, const int n_vect_prp, const int n_vect_l  ) {
    int j, jj;
    VECTOR_LOOP(j, n_vect, jj, l_pt[ 0*n_vect_l+j+jj] -= prp_su3_pt[0*n_vect_prp+j+jj];
                               l_pt[ 1*n_vect_l+j+jj] -= prp_su3_pt[1*n_vect_prp+j+jj];
                               l_pt[ 2*n_vect_l+j+jj] -= prp_su3_pt[2*n_vect_prp+j+jj];
                               l_pt[ 3*n_vect_l+j+jj] -= prp_su3_pt[3*n_vect_prp+j+jj];
                               l_pt[ 4*n_vect_l+j+jj] -= prp_su3_pt[4*n_vect_prp+j+jj];
                               l_pt[ 5*n_vect_l+j+jj] -= prp_su3_pt[5*n_vect_prp+j+jj];
                               l_pt[ 6*n_vect_l+j+jj] += GAMMA_Y_SPIN2_VAL*prp_su3_pt[ 3*GAMMA_Y_SPIN2_CO   *n_vect_prp+j+jj];
                               l_pt[ 7*n_vect_l+j+jj] += GAMMA_Y_SPIN2_VAL*prp_su3_pt[(3*GAMMA_Y_SPIN2_CO+1)*n_vect_prp+j+jj];
                               l_pt[ 8*n_vect_l+j+jj] += GAMMA_Y_SPIN2_VAL*prp_su3_pt[(3*GAMMA_Y_SPIN2_CO+2)*n_vect_prp+j+jj];
                               l_pt[ 9*n_vect_l+j+jj] += GAMMA_Y_SPIN3_VAL*prp_su3_pt[ 3*GAMMA_Y_SPIN3_CO   *n_vect_prp+j+jj];
                               l_pt[10*n_vect_l+j+jj] += GAMMA_Y_SPIN3_VAL*prp_su3_pt[(3*GAMMA_Y_SPIN3_CO+1)*n_vect_prp+j+jj];
                               l_pt[11*n_vect_l+j+jj] += GAMMA_Y_SPIN3_VAL*prp_su3_pt[(3*GAMMA_Y_SPIN3_CO+2)*n_vect_prp+j+jj];)
  }

  // -(1 + gamma_Y)
  static inline void pbn_su3_Y_PRECISION_new( const buffer_PRECISION prn_su3_pt, const buffer_PRECISION l_pt, const int n_vect, const int n_vect_prn, const int n_vect_l  ) {
    int j, jj;
    VECTOR_LOOP(j, n_vect, jj, l_pt[ 0*n_vect_l+j+jj] -= prn_su3_pt[0*n_vect_prn+j+jj];
                               l_pt[ 1*n_vect_l+j+jj] -= prn_su3_pt[1*n_vect_prn+j+jj];
                               l_pt[ 2*n_vect_l+j+jj] -= prn_su3_pt[2*n_vect_prn+j+jj];
		               l_pt[ 3*n_vect_l+j+jj] -= prn_su3_pt[3*n_vect_prn+j+jj];
                               l_pt[ 4*n_vect_l+j+jj] -= prn_su3_pt[4*n_vect_prn+j+jj];
                               l_pt[ 5*n_vect_l+j+jj] -= prn_su3_pt[5*n_vect_prn+j+jj];
                               l_pt[ 6*n_vect_l+j+jj] -= GAMMA_Y_SPIN2_VAL*prn_su3_pt[ 3*GAMMA_Y_SPIN2_CO   *n_vect_prn+j+jj];
                               l_pt[ 7*n_vect_l+j+jj] -= GAMMA_Y_SPIN2_VAL*prn_su3_pt[(3*GAMMA_Y_SPIN2_CO+1)*n_vect_prn+j+jj];
                               l_pt[ 8*n_vect_l+j+jj] -= GAMMA_Y_SPIN2_VAL*prn_su3_pt[(3*GAMMA_Y_SPIN2_CO+2)*n_vect_prn+j+jj];
                               l_pt[ 9*n_vect_l+j+jj] -= GAMMA_Y_SPIN3_VAL*prn_su3_pt[ 3*GAMMA_Y_SPIN3_CO   *n_vect_prn+j+jj];
                               l_pt[10*n_vect_l+j+jj] -= GAMMA_Y_SPIN3_VAL*prn_su3_pt[(3*GAMMA_Y_SPIN3_CO+1)*n_vect_prn+j+jj];
                               l_pt[11*n_vect_l+j+jj] -= GAMMA_Y_SPIN3_VAL*prn_su3_pt[(3*GAMMA_Y_SPIN3_CO+2)*n_vect_prn+j+jj];)
  }

  // 1 - gamma_X
  static inline void prp_X_PRECISION_new( const buffer_PRECISION prp_pt, const buffer_PRECISION l_pt, const int n_vect, const int n_vect_prp, const int n_vect_l  ) {
    int j, jj;
    VECTOR_LOOP(j, n_vect, jj, prp_pt[0*n_vect_prp+j+jj] = l_pt[0*n_vect_l+j+jj] -GAMMA_X_SPIN0_VAL*l_pt[ 3*GAMMA_X_SPIN0_CO   *n_vect_l+j+jj];
                               prp_pt[1*n_vect_prp+j+jj] = l_pt[1*n_vect_l+j+jj] -GAMMA_X_SPIN0_VAL*l_pt[(3*GAMMA_X_SPIN0_CO+1)*n_vect_l+j+jj];
                               prp_pt[2*n_vect_prp+j+jj] = l_pt[2*n_vect_l+j+jj] -GAMMA_X_SPIN0_VAL*l_pt[(3*GAMMA_X_SPIN0_CO+2)*n_vect_l+j+jj];
                               prp_pt[3*n_vect_prp+j+jj] = l_pt[3*n_vect_l+j+jj] -GAMMA_X_SPIN1_VAL*l_pt[ 3*GAMMA_X_SPIN1_CO   *n_vect_l+j+jj];
                               prp_pt[4*n_vect_prp+j+jj] = l_pt[4*n_vect_l+j+jj] -GAMMA_X_SPIN1_VAL*l_pt[(3*GAMMA_X_SPIN1_CO+1)*n_vect_l+j+jj];
                               prp_pt[5*n_vect_prp+j+jj] = l_pt[5*n_vect_l+j+jj] -GAMMA_X_SPIN1_VAL*l_pt[(3*GAMMA_X_SPIN1_CO+2)*n_vect_l+j+jj];)
  }

  // 1 + gamma_X
  static inline void prn_X_PRECISION_new( const buffer_PRECISION prn_pt, const buffer_PRECISION l_pt, const int n_vect, const int n_vect_prn, const int n_vect_l  ) {
    int j, jj;
    VECTOR_LOOP(j, n_vect, jj, prn_pt[0*n_vect_prn+j+jj] = l_pt[0*n_vect_l+j+jj] +GAMMA_X_SPIN0_VAL*l_pt[ 3*GAMMA_X_SPIN0_CO   *n_vect_l+j+jj];
                               prn_pt[1*n_vect_prn+j+jj] = l_pt[1*n_vect_l+j+jj] +GAMMA_X_SPIN0_VAL*l_pt[(3*GAMMA_X_SPIN0_CO+1)*n_vect_l+j+jj];
                               prn_pt[2*n_vect_prn+j+jj] = l_pt[2*n_vect_l+j+jj] +GAMMA_X_SPIN0_VAL*l_pt[(3*GAMMA_X_SPIN0_CO+2)*n_vect_l+j+jj];
                               prn_pt[3*n_vect_prn+j+jj] = l_pt[3*n_vect_l+j+jj] +GAMMA_X_SPIN1_VAL*l_pt[ 3*GAMMA_X_SPIN1_CO   *n_vect_l+j+jj];
                               prn_pt[4*n_vect_prn+j+jj] = l_pt[4*n_vect_l+j+jj] +GAMMA_X_SPIN1_VAL*l_pt[(3*GAMMA_X_SPIN1_CO+1)*n_vect_l+j+jj];
                               prn_pt[5*n_vect_prn+j+jj] = l_pt[5*n_vect_l+j+jj] +GAMMA_X_SPIN1_VAL*l_pt[(3*GAMMA_X_SPIN1_CO+2)*n_vect_l+j+jj];)
  }

  // - (1 - gamma_Y)
  static inline void pbp_su3_X_PRECISION_new( const buffer_PRECISION prp_su3_pt, const buffer_PRECISION l_pt, const int n_vect, const int n_vect_pbp, const int n_vect_l  ) {
    int j, jj;
    VECTOR_LOOP(j, n_vect, jj, l_pt[ 0*n_vect_l+j+jj] -= prp_su3_pt[0*n_vect_pbp+j+jj];
                               l_pt[ 1*n_vect_l+j+jj] -= prp_su3_pt[1*n_vect_pbp+j+jj];
                               l_pt[ 2*n_vect_l+j+jj] -= prp_su3_pt[2*n_vect_pbp+j+jj];
                               l_pt[ 3*n_vect_l+j+jj] -= prp_su3_pt[3*n_vect_pbp+j+jj];
                               l_pt[ 4*n_vect_l+j+jj] -= prp_su3_pt[4*n_vect_pbp+j+jj];
                               l_pt[ 5*n_vect_l+j+jj] -= prp_su3_pt[5*n_vect_pbp+j+jj];
                               l_pt[ 6*n_vect_l+j+jj] += GAMMA_X_SPIN2_VAL*prp_su3_pt[ 3*GAMMA_X_SPIN2_CO   *n_vect_pbp+j+jj];
                               l_pt[ 7*n_vect_l+j+jj] += GAMMA_X_SPIN2_VAL*prp_su3_pt[(3*GAMMA_X_SPIN2_CO+1)*n_vect_pbp+j+jj];
                               l_pt[ 8*n_vect_l+j+jj] += GAMMA_X_SPIN2_VAL*prp_su3_pt[(3*GAMMA_X_SPIN2_CO+2)*n_vect_pbp+j+jj];
                               l_pt[ 9*n_vect_l+j+jj] += GAMMA_X_SPIN3_VAL*prp_su3_pt[ 3*GAMMA_X_SPIN3_CO   *n_vect_pbp+j+jj];
                               l_pt[10*n_vect_l+j+jj] += GAMMA_X_SPIN3_VAL*prp_su3_pt[(3*GAMMA_X_SPIN3_CO+1)*n_vect_pbp+j+jj];
                               l_pt[11*n_vect_l+j+jj] += GAMMA_X_SPIN3_VAL*prp_su3_pt[(3*GAMMA_X_SPIN3_CO+2)*n_vect_pbp+j+jj];)
  }

  // -(1 + gamma_X)
  static inline void pbn_su3_X_PRECISION_new( const buffer_PRECISION prn_su3_pt, const buffer_PRECISION l_pt, const int n_vect, const int n_vect_pbn, const int n_vect_l  ) {
    int j, jj;
    VECTOR_LOOP(j, n_vect, jj, l_pt[ 0*n_vect_l+j+jj] -= prn_su3_pt[0*n_vect_pbn+j+jj];
                               l_pt[ 1*n_vect_l+j+jj] -= prn_su3_pt[1*n_vect_pbn+j+jj];
                               l_pt[ 2*n_vect_l+j+jj] -= prn_su3_pt[2*n_vect_pbn+j+jj];
                               l_pt[ 3*n_vect_l+j+jj] -= prn_su3_pt[3*n_vect_pbn+j+jj];
                               l_pt[ 4*n_vect_l+j+jj] -= prn_su3_pt[4*n_vect_pbn+j+jj];
                               l_pt[ 5*n_vect_l+j+jj] -= prn_su3_pt[5*n_vect_pbn+j+jj];
                               l_pt[ 6*n_vect_l+j+jj] -= GAMMA_X_SPIN2_VAL*prn_su3_pt[ 3*GAMMA_X_SPIN2_CO   *n_vect_pbn+j+jj];
                               l_pt[ 7*n_vect_l+j+jj] -= GAMMA_X_SPIN2_VAL*prn_su3_pt[(3*GAMMA_X_SPIN2_CO+1)*n_vect_pbn+j+jj];
                               l_pt[ 8*n_vect_l+j+jj] -= GAMMA_X_SPIN2_VAL*prn_su3_pt[(3*GAMMA_X_SPIN2_CO+2)*n_vect_pbn+j+jj];
                               l_pt[ 9*n_vect_l+j+jj] -= GAMMA_X_SPIN3_VAL*prn_su3_pt[ 3*GAMMA_X_SPIN3_CO   *n_vect_pbn+j+jj];
                               l_pt[10*n_vect_l+j+jj] -= GAMMA_X_SPIN3_VAL*prn_su3_pt[(3*GAMMA_X_SPIN3_CO+1)*n_vect_pbn+j+jj];
                               l_pt[11*n_vect_l+j+jj] -= GAMMA_X_SPIN3_VAL*prn_su3_pt[(3*GAMMA_X_SPIN3_CO+2)*n_vect_pbn+j+jj];)
  }

/***************************  TM1p1   ********************************************/
#ifdef HAVE_TM1p1

//#define flav_gamma(k) ((k)>1?((k)*3+6):((k)*3))
#define flav_gamma(k) (3*(k)+6*((k)/2))

  // 1 - gamma_T
  static inline void dprp_T_PRECISION( const buffer_PRECISION prp_pt, const buffer_PRECISION l_pt ) {
    prp_pt[ 0] = l_pt[ 0] -GAMMA_T_SPIN0_VAL*l_pt[flav_gamma(GAMMA_T_SPIN0_CO)];
    prp_pt[ 1] = l_pt[ 1] -GAMMA_T_SPIN0_VAL*l_pt[flav_gamma(GAMMA_T_SPIN0_CO)+1];
    prp_pt[ 2] = l_pt[ 2] -GAMMA_T_SPIN0_VAL*l_pt[flav_gamma(GAMMA_T_SPIN0_CO)+2];
    prp_pt[ 3] = l_pt[ 3] -GAMMA_T_SPIN1_VAL*l_pt[flav_gamma(GAMMA_T_SPIN1_CO)];
    prp_pt[ 4] = l_pt[ 4] -GAMMA_T_SPIN1_VAL*l_pt[flav_gamma(GAMMA_T_SPIN1_CO)+1];
    prp_pt[ 5] = l_pt[ 5] -GAMMA_T_SPIN1_VAL*l_pt[flav_gamma(GAMMA_T_SPIN1_CO)+2];
    prp_pt[ 6] = l_pt[ 6] -GAMMA_T_SPIN0_VAL*l_pt[flav_gamma(GAMMA_T_SPIN0_CO)+6];
    prp_pt[ 7] = l_pt[ 7] -GAMMA_T_SPIN0_VAL*l_pt[flav_gamma(GAMMA_T_SPIN0_CO)+7];
    prp_pt[ 8] = l_pt[ 8] -GAMMA_T_SPIN0_VAL*l_pt[flav_gamma(GAMMA_T_SPIN0_CO)+8];
    prp_pt[ 9] = l_pt[ 9] -GAMMA_T_SPIN1_VAL*l_pt[flav_gamma(GAMMA_T_SPIN1_CO)+6];
    prp_pt[10] = l_pt[10] -GAMMA_T_SPIN1_VAL*l_pt[flav_gamma(GAMMA_T_SPIN1_CO)+7];
    prp_pt[11] = l_pt[11] -GAMMA_T_SPIN1_VAL*l_pt[flav_gamma(GAMMA_T_SPIN1_CO)+8];
  }

  // 1 + gamma_T
  static inline void dprn_T_PRECISION( const buffer_PRECISION prn_pt, const buffer_PRECISION l_pt ) {
    prn_pt[ 0] = l_pt[ 0] +GAMMA_T_SPIN0_VAL*l_pt[flav_gamma(GAMMA_T_SPIN0_CO)];
    prn_pt[ 1] = l_pt[ 1] +GAMMA_T_SPIN0_VAL*l_pt[flav_gamma(GAMMA_T_SPIN0_CO)+1];
    prn_pt[ 2] = l_pt[ 2] +GAMMA_T_SPIN0_VAL*l_pt[flav_gamma(GAMMA_T_SPIN0_CO)+2];
    prn_pt[ 3] = l_pt[ 3] +GAMMA_T_SPIN1_VAL*l_pt[flav_gamma(GAMMA_T_SPIN1_CO)];
    prn_pt[ 4] = l_pt[ 4] +GAMMA_T_SPIN1_VAL*l_pt[flav_gamma(GAMMA_T_SPIN1_CO)+1];
    prn_pt[ 5] = l_pt[ 5] +GAMMA_T_SPIN1_VAL*l_pt[flav_gamma(GAMMA_T_SPIN1_CO)+2];
    prn_pt[ 6] = l_pt[ 6] +GAMMA_T_SPIN0_VAL*l_pt[flav_gamma(GAMMA_T_SPIN0_CO)+6];
    prn_pt[ 7] = l_pt[ 7] +GAMMA_T_SPIN0_VAL*l_pt[flav_gamma(GAMMA_T_SPIN0_CO)+7];
    prn_pt[ 8] = l_pt[ 8] +GAMMA_T_SPIN0_VAL*l_pt[flav_gamma(GAMMA_T_SPIN0_CO)+8];
    prn_pt[ 9] = l_pt[ 9] +GAMMA_T_SPIN1_VAL*l_pt[flav_gamma(GAMMA_T_SPIN1_CO)+6];
    prn_pt[10] = l_pt[10] +GAMMA_T_SPIN1_VAL*l_pt[flav_gamma(GAMMA_T_SPIN1_CO)+7];
    prn_pt[11] = l_pt[11] +GAMMA_T_SPIN1_VAL*l_pt[flav_gamma(GAMMA_T_SPIN1_CO)+8];
  }

  // - (1 - gamma_T)
  static inline void dpbp_su3_T_PRECISION( const buffer_PRECISION prp_su3_pt, const buffer_PRECISION l_pt ) {
    l_pt[ 0] -= prp_su3_pt[ 0];
    l_pt[ 1] -= prp_su3_pt[ 1];
    l_pt[ 2] -= prp_su3_pt[ 2];
    l_pt[ 3] -= prp_su3_pt[ 3];
    l_pt[ 4] -= prp_su3_pt[ 4];
    l_pt[ 5] -= prp_su3_pt[ 5];
    l_pt[ 6] -= prp_su3_pt[ 6];
    l_pt[ 7] -= prp_su3_pt[ 7];
    l_pt[ 8] -= prp_su3_pt[ 8];
    l_pt[ 9] -= prp_su3_pt[ 9];
    l_pt[10] -= prp_su3_pt[10];
    l_pt[11] -= prp_su3_pt[11];
    l_pt[12] += GAMMA_T_SPIN2_VAL*prp_su3_pt[flav_gamma(GAMMA_T_SPIN2_CO)];
    l_pt[13] += GAMMA_T_SPIN2_VAL*prp_su3_pt[flav_gamma(GAMMA_T_SPIN2_CO)+1];
    l_pt[14] += GAMMA_T_SPIN2_VAL*prp_su3_pt[flav_gamma(GAMMA_T_SPIN2_CO)+2];
    l_pt[15] += GAMMA_T_SPIN3_VAL*prp_su3_pt[flav_gamma(GAMMA_T_SPIN3_CO)];
    l_pt[16] += GAMMA_T_SPIN3_VAL*prp_su3_pt[flav_gamma(GAMMA_T_SPIN3_CO)+1];
    l_pt[17] += GAMMA_T_SPIN3_VAL*prp_su3_pt[flav_gamma(GAMMA_T_SPIN3_CO)+2];
    l_pt[18] += GAMMA_T_SPIN2_VAL*prp_su3_pt[flav_gamma(GAMMA_T_SPIN2_CO)+6];
    l_pt[19] += GAMMA_T_SPIN2_VAL*prp_su3_pt[flav_gamma(GAMMA_T_SPIN2_CO)+7];
    l_pt[20] += GAMMA_T_SPIN2_VAL*prp_su3_pt[flav_gamma(GAMMA_T_SPIN2_CO)+8];
    l_pt[21] += GAMMA_T_SPIN3_VAL*prp_su3_pt[flav_gamma(GAMMA_T_SPIN3_CO)+6];
    l_pt[22] += GAMMA_T_SPIN3_VAL*prp_su3_pt[flav_gamma(GAMMA_T_SPIN3_CO)+7];
    l_pt[23] += GAMMA_T_SPIN3_VAL*prp_su3_pt[flav_gamma(GAMMA_T_SPIN3_CO)+8];
  }

  // -(1 + gamma_T)
  static inline void dpbn_su3_T_PRECISION( const buffer_PRECISION prn_su3_pt, const buffer_PRECISION l_pt ) {
    l_pt[ 0] -= prn_su3_pt[ 0];
    l_pt[ 1] -= prn_su3_pt[ 1];
    l_pt[ 2] -= prn_su3_pt[ 2];
    l_pt[ 3] -= prn_su3_pt[ 3];
    l_pt[ 4] -= prn_su3_pt[ 4];
    l_pt[ 5] -= prn_su3_pt[ 5];
    l_pt[ 6] -= prn_su3_pt[ 6];
    l_pt[ 7] -= prn_su3_pt[ 7];
    l_pt[ 8] -= prn_su3_pt[ 8];
    l_pt[ 9] -= prn_su3_pt[ 9];
    l_pt[10] -= prn_su3_pt[10];
    l_pt[11] -= prn_su3_pt[11];
    l_pt[12] -= GAMMA_T_SPIN2_VAL*prn_su3_pt[flav_gamma(GAMMA_T_SPIN2_CO)];
    l_pt[13] -= GAMMA_T_SPIN2_VAL*prn_su3_pt[flav_gamma(GAMMA_T_SPIN2_CO)+1];
    l_pt[14] -= GAMMA_T_SPIN2_VAL*prn_su3_pt[flav_gamma(GAMMA_T_SPIN2_CO)+2];
    l_pt[15] -= GAMMA_T_SPIN3_VAL*prn_su3_pt[flav_gamma(GAMMA_T_SPIN3_CO)];
    l_pt[16] -= GAMMA_T_SPIN3_VAL*prn_su3_pt[flav_gamma(GAMMA_T_SPIN3_CO)+1];
    l_pt[17] -= GAMMA_T_SPIN3_VAL*prn_su3_pt[flav_gamma(GAMMA_T_SPIN3_CO)+2];
    l_pt[18] -= GAMMA_T_SPIN2_VAL*prn_su3_pt[flav_gamma(GAMMA_T_SPIN2_CO)+6];
    l_pt[19] -= GAMMA_T_SPIN2_VAL*prn_su3_pt[flav_gamma(GAMMA_T_SPIN2_CO)+7];
    l_pt[20] -= GAMMA_T_SPIN2_VAL*prn_su3_pt[flav_gamma(GAMMA_T_SPIN2_CO)+8];
    l_pt[21] -= GAMMA_T_SPIN3_VAL*prn_su3_pt[flav_gamma(GAMMA_T_SPIN3_CO)+6];
    l_pt[22] -= GAMMA_T_SPIN3_VAL*prn_su3_pt[flav_gamma(GAMMA_T_SPIN3_CO)+7];
    l_pt[23] -= GAMMA_T_SPIN3_VAL*prn_su3_pt[flav_gamma(GAMMA_T_SPIN3_CO)+8];
  }


  // 1 - gamma_Z
  static inline void dprp_Z_PRECISION( const buffer_PRECISION prp_pt, const buffer_PRECISION l_pt ) {
    prp_pt[ 0] = l_pt[ 0] -GAMMA_Z_SPIN0_VAL*l_pt[flav_gamma(GAMMA_Z_SPIN0_CO)];
    prp_pt[ 1] = l_pt[ 1] -GAMMA_Z_SPIN0_VAL*l_pt[flav_gamma(GAMMA_Z_SPIN0_CO)+1];
    prp_pt[ 2] = l_pt[ 2] -GAMMA_Z_SPIN0_VAL*l_pt[flav_gamma(GAMMA_Z_SPIN0_CO)+2];
    prp_pt[ 3] = l_pt[ 3] -GAMMA_Z_SPIN1_VAL*l_pt[flav_gamma(GAMMA_Z_SPIN1_CO)];
    prp_pt[ 4] = l_pt[ 4] -GAMMA_Z_SPIN1_VAL*l_pt[flav_gamma(GAMMA_Z_SPIN1_CO)+1];
    prp_pt[ 5] = l_pt[ 5] -GAMMA_Z_SPIN1_VAL*l_pt[flav_gamma(GAMMA_Z_SPIN1_CO)+2];
    prp_pt[ 6] = l_pt[ 6] -GAMMA_Z_SPIN0_VAL*l_pt[flav_gamma(GAMMA_Z_SPIN0_CO)+6];
    prp_pt[ 7] = l_pt[ 7] -GAMMA_Z_SPIN0_VAL*l_pt[flav_gamma(GAMMA_Z_SPIN0_CO)+7];
    prp_pt[ 8] = l_pt[ 8] -GAMMA_Z_SPIN0_VAL*l_pt[flav_gamma(GAMMA_Z_SPIN0_CO)+8];
    prp_pt[ 9] = l_pt[ 9] -GAMMA_Z_SPIN1_VAL*l_pt[flav_gamma(GAMMA_Z_SPIN1_CO)+6];
    prp_pt[10] = l_pt[10] -GAMMA_Z_SPIN1_VAL*l_pt[flav_gamma(GAMMA_Z_SPIN1_CO)+7];
    prp_pt[11] = l_pt[11] -GAMMA_Z_SPIN1_VAL*l_pt[flav_gamma(GAMMA_Z_SPIN1_CO)+8];
  }

  // 1 + gamma_Z
  static inline void dprn_Z_PRECISION( const buffer_PRECISION prn_pt, const buffer_PRECISION l_pt ) {
    prn_pt[ 0] = l_pt[ 0] +GAMMA_Z_SPIN0_VAL*l_pt[flav_gamma(GAMMA_Z_SPIN0_CO)];
    prn_pt[ 1] = l_pt[ 1] +GAMMA_Z_SPIN0_VAL*l_pt[flav_gamma(GAMMA_Z_SPIN0_CO)+1];
    prn_pt[ 2] = l_pt[ 2] +GAMMA_Z_SPIN0_VAL*l_pt[flav_gamma(GAMMA_Z_SPIN0_CO)+2];
    prn_pt[ 3] = l_pt[ 3] +GAMMA_Z_SPIN1_VAL*l_pt[flav_gamma(GAMMA_Z_SPIN1_CO)];
    prn_pt[ 4] = l_pt[ 4] +GAMMA_Z_SPIN1_VAL*l_pt[flav_gamma(GAMMA_Z_SPIN1_CO)+1];
    prn_pt[ 5] = l_pt[ 5] +GAMMA_Z_SPIN1_VAL*l_pt[flav_gamma(GAMMA_Z_SPIN1_CO)+2];
    prn_pt[ 6] = l_pt[ 6] +GAMMA_Z_SPIN0_VAL*l_pt[flav_gamma(GAMMA_Z_SPIN0_CO)+6];
    prn_pt[ 7] = l_pt[ 7] +GAMMA_Z_SPIN0_VAL*l_pt[flav_gamma(GAMMA_Z_SPIN0_CO)+7];
    prn_pt[ 8] = l_pt[ 8] +GAMMA_Z_SPIN0_VAL*l_pt[flav_gamma(GAMMA_Z_SPIN0_CO)+8];
    prn_pt[ 9] = l_pt[ 9] +GAMMA_Z_SPIN1_VAL*l_pt[flav_gamma(GAMMA_Z_SPIN1_CO)+6];
    prn_pt[10] = l_pt[10] +GAMMA_Z_SPIN1_VAL*l_pt[flav_gamma(GAMMA_Z_SPIN1_CO)+7];
    prn_pt[11] = l_pt[11] +GAMMA_Z_SPIN1_VAL*l_pt[flav_gamma(GAMMA_Z_SPIN1_CO)+8];
  }

  // - (1 - gamma_Z)
  static inline void dpbp_su3_Z_PRECISION( const buffer_PRECISION prp_su3_pt, const buffer_PRECISION l_pt ) {
    l_pt[ 0] -= prp_su3_pt[ 0];
    l_pt[ 1] -= prp_su3_pt[ 1];
    l_pt[ 2] -= prp_su3_pt[ 2];
    l_pt[ 3] -= prp_su3_pt[ 3];
    l_pt[ 4] -= prp_su3_pt[ 4];
    l_pt[ 5] -= prp_su3_pt[ 5];
    l_pt[ 6] -= prp_su3_pt[ 6];
    l_pt[ 7] -= prp_su3_pt[ 7];
    l_pt[ 8] -= prp_su3_pt[ 8];
    l_pt[ 9] -= prp_su3_pt[ 9];
    l_pt[10] -= prp_su3_pt[10];
    l_pt[11] -= prp_su3_pt[11];
    l_pt[12] += GAMMA_Z_SPIN2_VAL*prp_su3_pt[flav_gamma(GAMMA_Z_SPIN2_CO)];
    l_pt[13] += GAMMA_Z_SPIN2_VAL*prp_su3_pt[flav_gamma(GAMMA_Z_SPIN2_CO)+1];
    l_pt[14] += GAMMA_Z_SPIN2_VAL*prp_su3_pt[flav_gamma(GAMMA_Z_SPIN2_CO)+2];
    l_pt[15] += GAMMA_Z_SPIN3_VAL*prp_su3_pt[flav_gamma(GAMMA_Z_SPIN3_CO)];
    l_pt[16] += GAMMA_Z_SPIN3_VAL*prp_su3_pt[flav_gamma(GAMMA_Z_SPIN3_CO)+1];
    l_pt[17] += GAMMA_Z_SPIN3_VAL*prp_su3_pt[flav_gamma(GAMMA_Z_SPIN3_CO)+2];
    l_pt[18] += GAMMA_Z_SPIN2_VAL*prp_su3_pt[flav_gamma(GAMMA_Z_SPIN2_CO)+6];
    l_pt[19] += GAMMA_Z_SPIN2_VAL*prp_su3_pt[flav_gamma(GAMMA_Z_SPIN2_CO)+7];
    l_pt[20] += GAMMA_Z_SPIN2_VAL*prp_su3_pt[flav_gamma(GAMMA_Z_SPIN2_CO)+8];
    l_pt[21] += GAMMA_Z_SPIN3_VAL*prp_su3_pt[flav_gamma(GAMMA_Z_SPIN3_CO)+6];
    l_pt[22] += GAMMA_Z_SPIN3_VAL*prp_su3_pt[flav_gamma(GAMMA_Z_SPIN3_CO)+7];
    l_pt[23] += GAMMA_Z_SPIN3_VAL*prp_su3_pt[flav_gamma(GAMMA_Z_SPIN3_CO)+8];
  }

  // -(1 + gamma_Z)
  static inline void dpbn_su3_Z_PRECISION( const buffer_PRECISION prn_su3_pt, const buffer_PRECISION l_pt ) {
    l_pt[ 0] -= prn_su3_pt[ 0];
    l_pt[ 1] -= prn_su3_pt[ 1];
    l_pt[ 2] -= prn_su3_pt[ 2];
    l_pt[ 3] -= prn_su3_pt[ 3];
    l_pt[ 4] -= prn_su3_pt[ 4];
    l_pt[ 5] -= prn_su3_pt[ 5];
    l_pt[ 6] -= prn_su3_pt[ 6];
    l_pt[ 7] -= prn_su3_pt[ 7];
    l_pt[ 8] -= prn_su3_pt[ 8];
    l_pt[ 9] -= prn_su3_pt[ 9];
    l_pt[10] -= prn_su3_pt[10];
    l_pt[11] -= prn_su3_pt[11];
    l_pt[12] -= GAMMA_Z_SPIN2_VAL*prn_su3_pt[flav_gamma(GAMMA_Z_SPIN2_CO)];
    l_pt[13] -= GAMMA_Z_SPIN2_VAL*prn_su3_pt[flav_gamma(GAMMA_Z_SPIN2_CO)+1];
    l_pt[14] -= GAMMA_Z_SPIN2_VAL*prn_su3_pt[flav_gamma(GAMMA_Z_SPIN2_CO)+2];
    l_pt[15] -= GAMMA_Z_SPIN3_VAL*prn_su3_pt[flav_gamma(GAMMA_Z_SPIN3_CO)];
    l_pt[16] -= GAMMA_Z_SPIN3_VAL*prn_su3_pt[flav_gamma(GAMMA_Z_SPIN3_CO)+1];
    l_pt[17] -= GAMMA_Z_SPIN3_VAL*prn_su3_pt[flav_gamma(GAMMA_Z_SPIN3_CO)+2];
    l_pt[18] -= GAMMA_Z_SPIN2_VAL*prn_su3_pt[flav_gamma(GAMMA_Z_SPIN2_CO)+6];
    l_pt[19] -= GAMMA_Z_SPIN2_VAL*prn_su3_pt[flav_gamma(GAMMA_Z_SPIN2_CO)+7];
    l_pt[20] -= GAMMA_Z_SPIN2_VAL*prn_su3_pt[flav_gamma(GAMMA_Z_SPIN2_CO)+8];
    l_pt[21] -= GAMMA_Z_SPIN3_VAL*prn_su3_pt[flav_gamma(GAMMA_Z_SPIN3_CO)+6];
    l_pt[22] -= GAMMA_Z_SPIN3_VAL*prn_su3_pt[flav_gamma(GAMMA_Z_SPIN3_CO)+7];
    l_pt[23] -= GAMMA_Z_SPIN3_VAL*prn_su3_pt[flav_gamma(GAMMA_Z_SPIN3_CO)+8];
  }


  // 1 - gamma_Y
  static inline void dprp_Y_PRECISION( const buffer_PRECISION prp_pt, const buffer_PRECISION l_pt ) {
    prp_pt[ 0] = l_pt[ 0] -GAMMA_Y_SPIN0_VAL*l_pt[flav_gamma(GAMMA_Y_SPIN0_CO)];
    prp_pt[ 1] = l_pt[ 1] -GAMMA_Y_SPIN0_VAL*l_pt[flav_gamma(GAMMA_Y_SPIN0_CO)+1];
    prp_pt[ 2] = l_pt[ 2] -GAMMA_Y_SPIN0_VAL*l_pt[flav_gamma(GAMMA_Y_SPIN0_CO)+2];
    prp_pt[ 3] = l_pt[ 3] -GAMMA_Y_SPIN1_VAL*l_pt[flav_gamma(GAMMA_Y_SPIN1_CO)];
    prp_pt[ 4] = l_pt[ 4] -GAMMA_Y_SPIN1_VAL*l_pt[flav_gamma(GAMMA_Y_SPIN1_CO)+1];
    prp_pt[ 5] = l_pt[ 5] -GAMMA_Y_SPIN1_VAL*l_pt[flav_gamma(GAMMA_Y_SPIN1_CO)+2];
    prp_pt[ 6] = l_pt[ 6] -GAMMA_Y_SPIN0_VAL*l_pt[flav_gamma(GAMMA_Y_SPIN0_CO)+6];
    prp_pt[ 7] = l_pt[ 7] -GAMMA_Y_SPIN0_VAL*l_pt[flav_gamma(GAMMA_Y_SPIN0_CO)+7];
    prp_pt[ 8] = l_pt[ 8] -GAMMA_Y_SPIN0_VAL*l_pt[flav_gamma(GAMMA_Y_SPIN0_CO)+8];
    prp_pt[ 9] = l_pt[ 9] -GAMMA_Y_SPIN1_VAL*l_pt[flav_gamma(GAMMA_Y_SPIN1_CO)+6];
    prp_pt[10] = l_pt[10] -GAMMA_Y_SPIN1_VAL*l_pt[flav_gamma(GAMMA_Y_SPIN1_CO)+7];
    prp_pt[11] = l_pt[11] -GAMMA_Y_SPIN1_VAL*l_pt[flav_gamma(GAMMA_Y_SPIN1_CO)+8];
  }

  // 1 + gamma_Y
  static inline void dprn_Y_PRECISION( const buffer_PRECISION prn_pt, const buffer_PRECISION l_pt ) {
    prn_pt[ 0] = l_pt[ 0] +GAMMA_Y_SPIN0_VAL*l_pt[flav_gamma(GAMMA_Y_SPIN0_CO)];
    prn_pt[ 1] = l_pt[ 1] +GAMMA_Y_SPIN0_VAL*l_pt[flav_gamma(GAMMA_Y_SPIN0_CO)+1];
    prn_pt[ 2] = l_pt[ 2] +GAMMA_Y_SPIN0_VAL*l_pt[flav_gamma(GAMMA_Y_SPIN0_CO)+2];
    prn_pt[ 3] = l_pt[ 3] +GAMMA_Y_SPIN1_VAL*l_pt[flav_gamma(GAMMA_Y_SPIN1_CO)];
    prn_pt[ 4] = l_pt[ 4] +GAMMA_Y_SPIN1_VAL*l_pt[flav_gamma(GAMMA_Y_SPIN1_CO)+1];
    prn_pt[ 5] = l_pt[ 5] +GAMMA_Y_SPIN1_VAL*l_pt[flav_gamma(GAMMA_Y_SPIN1_CO)+2];
    prn_pt[ 6] = l_pt[ 6] +GAMMA_Y_SPIN0_VAL*l_pt[flav_gamma(GAMMA_Y_SPIN0_CO)+6];
    prn_pt[ 7] = l_pt[ 7] +GAMMA_Y_SPIN0_VAL*l_pt[flav_gamma(GAMMA_Y_SPIN0_CO)+7];
    prn_pt[ 8] = l_pt[ 8] +GAMMA_Y_SPIN0_VAL*l_pt[flav_gamma(GAMMA_Y_SPIN0_CO)+8];
    prn_pt[ 9] = l_pt[ 9] +GAMMA_Y_SPIN1_VAL*l_pt[flav_gamma(GAMMA_Y_SPIN1_CO)+6];
    prn_pt[10] = l_pt[10] +GAMMA_Y_SPIN1_VAL*l_pt[flav_gamma(GAMMA_Y_SPIN1_CO)+7];
    prn_pt[11] = l_pt[11] +GAMMA_Y_SPIN1_VAL*l_pt[flav_gamma(GAMMA_Y_SPIN1_CO)+8];
  }

  // - (1 - gamma_Y)
  static inline void dpbp_su3_Y_PRECISION( const buffer_PRECISION prp_su3_pt, const buffer_PRECISION l_pt ) {
    l_pt[ 0] -= prp_su3_pt[ 0];
    l_pt[ 1] -= prp_su3_pt[ 1];
    l_pt[ 2] -= prp_su3_pt[ 2];
    l_pt[ 3] -= prp_su3_pt[ 3];
    l_pt[ 4] -= prp_su3_pt[ 4];
    l_pt[ 5] -= prp_su3_pt[ 5];
    l_pt[ 6] -= prp_su3_pt[ 6];
    l_pt[ 7] -= prp_su3_pt[ 7];
    l_pt[ 8] -= prp_su3_pt[ 8];
    l_pt[ 9] -= prp_su3_pt[ 9];
    l_pt[10] -= prp_su3_pt[10];
    l_pt[11] -= prp_su3_pt[11];
    l_pt[12] += GAMMA_Y_SPIN2_VAL*prp_su3_pt[flav_gamma(GAMMA_Y_SPIN2_CO)];
    l_pt[13] += GAMMA_Y_SPIN2_VAL*prp_su3_pt[flav_gamma(GAMMA_Y_SPIN2_CO)+1];
    l_pt[14] += GAMMA_Y_SPIN2_VAL*prp_su3_pt[flav_gamma(GAMMA_Y_SPIN2_CO)+2];
    l_pt[15] += GAMMA_Y_SPIN3_VAL*prp_su3_pt[flav_gamma(GAMMA_Y_SPIN3_CO)];
    l_pt[16] += GAMMA_Y_SPIN3_VAL*prp_su3_pt[flav_gamma(GAMMA_Y_SPIN3_CO)+1];
    l_pt[17] += GAMMA_Y_SPIN3_VAL*prp_su3_pt[flav_gamma(GAMMA_Y_SPIN3_CO)+2];
    l_pt[18] += GAMMA_Y_SPIN2_VAL*prp_su3_pt[flav_gamma(GAMMA_Y_SPIN2_CO)+6];
    l_pt[19] += GAMMA_Y_SPIN2_VAL*prp_su3_pt[flav_gamma(GAMMA_Y_SPIN2_CO)+7];
    l_pt[20] += GAMMA_Y_SPIN2_VAL*prp_su3_pt[flav_gamma(GAMMA_Y_SPIN2_CO)+8];
    l_pt[21] += GAMMA_Y_SPIN3_VAL*prp_su3_pt[flav_gamma(GAMMA_Y_SPIN3_CO)+6];
    l_pt[22] += GAMMA_Y_SPIN3_VAL*prp_su3_pt[flav_gamma(GAMMA_Y_SPIN3_CO)+7];
    l_pt[23] += GAMMA_Y_SPIN3_VAL*prp_su3_pt[flav_gamma(GAMMA_Y_SPIN3_CO)+8];
  }

  // -(1 + gamma_Y)
  static inline void dpbn_su3_Y_PRECISION( const buffer_PRECISION prn_su3_pt, const buffer_PRECISION l_pt ) {
    l_pt[ 0] -= prn_su3_pt[ 0];
    l_pt[ 1] -= prn_su3_pt[ 1];
    l_pt[ 2] -= prn_su3_pt[ 2];
    l_pt[ 3] -= prn_su3_pt[ 3];
    l_pt[ 4] -= prn_su3_pt[ 4];
    l_pt[ 5] -= prn_su3_pt[ 5];
    l_pt[ 6] -= prn_su3_pt[ 6];
    l_pt[ 7] -= prn_su3_pt[ 7];
    l_pt[ 8] -= prn_su3_pt[ 8];
    l_pt[ 9] -= prn_su3_pt[ 9];
    l_pt[10] -= prn_su3_pt[10];
    l_pt[11] -= prn_su3_pt[11];
    l_pt[12] -= GAMMA_Y_SPIN2_VAL*prn_su3_pt[flav_gamma(GAMMA_Y_SPIN2_CO)];
    l_pt[13] -= GAMMA_Y_SPIN2_VAL*prn_su3_pt[flav_gamma(GAMMA_Y_SPIN2_CO)+1];
    l_pt[14] -= GAMMA_Y_SPIN2_VAL*prn_su3_pt[flav_gamma(GAMMA_Y_SPIN2_CO)+2];
    l_pt[15] -= GAMMA_Y_SPIN3_VAL*prn_su3_pt[flav_gamma(GAMMA_Y_SPIN3_CO)];
    l_pt[16] -= GAMMA_Y_SPIN3_VAL*prn_su3_pt[flav_gamma(GAMMA_Y_SPIN3_CO)+1];
    l_pt[17] -= GAMMA_Y_SPIN3_VAL*prn_su3_pt[flav_gamma(GAMMA_Y_SPIN3_CO)+2];
    l_pt[18] -= GAMMA_Y_SPIN2_VAL*prn_su3_pt[flav_gamma(GAMMA_Y_SPIN2_CO)+6];
    l_pt[19] -= GAMMA_Y_SPIN2_VAL*prn_su3_pt[flav_gamma(GAMMA_Y_SPIN2_CO)+7];
    l_pt[20] -= GAMMA_Y_SPIN2_VAL*prn_su3_pt[flav_gamma(GAMMA_Y_SPIN2_CO)+8];
    l_pt[21] -= GAMMA_Y_SPIN3_VAL*prn_su3_pt[flav_gamma(GAMMA_Y_SPIN3_CO)+6];
    l_pt[22] -= GAMMA_Y_SPIN3_VAL*prn_su3_pt[flav_gamma(GAMMA_Y_SPIN3_CO)+7];
    l_pt[23] -= GAMMA_Y_SPIN3_VAL*prn_su3_pt[flav_gamma(GAMMA_Y_SPIN3_CO)+8];
  }


  // 1 - gamma_X
  static inline void dprp_X_PRECISION( const buffer_PRECISION prp_pt, const buffer_PRECISION l_pt ) {
    prp_pt[ 0] = l_pt[ 0] -GAMMA_X_SPIN0_VAL*l_pt[flav_gamma(GAMMA_X_SPIN0_CO)];
    prp_pt[ 1] = l_pt[ 1] -GAMMA_X_SPIN0_VAL*l_pt[flav_gamma(GAMMA_X_SPIN0_CO)+1];
    prp_pt[ 2] = l_pt[ 2] -GAMMA_X_SPIN0_VAL*l_pt[flav_gamma(GAMMA_X_SPIN0_CO)+2];
    prp_pt[ 3] = l_pt[ 3] -GAMMA_X_SPIN1_VAL*l_pt[flav_gamma(GAMMA_X_SPIN1_CO)];
    prp_pt[ 4] = l_pt[ 4] -GAMMA_X_SPIN1_VAL*l_pt[flav_gamma(GAMMA_X_SPIN1_CO)+1];
    prp_pt[ 5] = l_pt[ 5] -GAMMA_X_SPIN1_VAL*l_pt[flav_gamma(GAMMA_X_SPIN1_CO)+2];
    prp_pt[ 6] = l_pt[ 6] -GAMMA_X_SPIN0_VAL*l_pt[flav_gamma(GAMMA_X_SPIN0_CO)+6];
    prp_pt[ 7] = l_pt[ 7] -GAMMA_X_SPIN0_VAL*l_pt[flav_gamma(GAMMA_X_SPIN0_CO)+7];
    prp_pt[ 8] = l_pt[ 8] -GAMMA_X_SPIN0_VAL*l_pt[flav_gamma(GAMMA_X_SPIN0_CO)+8];
    prp_pt[ 9] = l_pt[ 9] -GAMMA_X_SPIN1_VAL*l_pt[flav_gamma(GAMMA_X_SPIN1_CO)+6];
    prp_pt[10] = l_pt[10] -GAMMA_X_SPIN1_VAL*l_pt[flav_gamma(GAMMA_X_SPIN1_CO)+7];
    prp_pt[11] = l_pt[11] -GAMMA_X_SPIN1_VAL*l_pt[flav_gamma(GAMMA_X_SPIN1_CO)+8];
  }

  // 1 + gamma_X
  static inline void dprn_X_PRECISION( const buffer_PRECISION prn_pt, const buffer_PRECISION l_pt ) {
    prn_pt[ 0] = l_pt[ 0] +GAMMA_X_SPIN0_VAL*l_pt[flav_gamma(GAMMA_X_SPIN0_CO)];
    prn_pt[ 1] = l_pt[ 1] +GAMMA_X_SPIN0_VAL*l_pt[flav_gamma(GAMMA_X_SPIN0_CO)+1];
    prn_pt[ 2] = l_pt[ 2] +GAMMA_X_SPIN0_VAL*l_pt[flav_gamma(GAMMA_X_SPIN0_CO)+2];
    prn_pt[ 3] = l_pt[ 3] +GAMMA_X_SPIN1_VAL*l_pt[flav_gamma(GAMMA_X_SPIN1_CO)];
    prn_pt[ 4] = l_pt[ 4] +GAMMA_X_SPIN1_VAL*l_pt[flav_gamma(GAMMA_X_SPIN1_CO)+1];
    prn_pt[ 5] = l_pt[ 5] +GAMMA_X_SPIN1_VAL*l_pt[flav_gamma(GAMMA_X_SPIN1_CO)+2];
    prn_pt[ 6] = l_pt[ 6] +GAMMA_X_SPIN0_VAL*l_pt[flav_gamma(GAMMA_X_SPIN0_CO)+6];
    prn_pt[ 7] = l_pt[ 7] +GAMMA_X_SPIN0_VAL*l_pt[flav_gamma(GAMMA_X_SPIN0_CO)+7];
    prn_pt[ 8] = l_pt[ 8] +GAMMA_X_SPIN0_VAL*l_pt[flav_gamma(GAMMA_X_SPIN0_CO)+8];
    prn_pt[ 9] = l_pt[ 9] +GAMMA_X_SPIN1_VAL*l_pt[flav_gamma(GAMMA_X_SPIN1_CO)+6];
    prn_pt[10] = l_pt[10] +GAMMA_X_SPIN1_VAL*l_pt[flav_gamma(GAMMA_X_SPIN1_CO)+7];
    prn_pt[11] = l_pt[11] +GAMMA_X_SPIN1_VAL*l_pt[flav_gamma(GAMMA_X_SPIN1_CO)+8];
  }

  // - (1 - gamma_X)
  static inline void dpbp_su3_X_PRECISION( const buffer_PRECISION prp_su3_pt, const buffer_PRECISION l_pt ) {
    l_pt[ 0] -= prp_su3_pt[ 0];
    l_pt[ 1] -= prp_su3_pt[ 1];
    l_pt[ 2] -= prp_su3_pt[ 2];
    l_pt[ 3] -= prp_su3_pt[ 3];
    l_pt[ 4] -= prp_su3_pt[ 4];
    l_pt[ 5] -= prp_su3_pt[ 5];
    l_pt[ 6] -= prp_su3_pt[ 6];
    l_pt[ 7] -= prp_su3_pt[ 7];
    l_pt[ 8] -= prp_su3_pt[ 8];
    l_pt[ 9] -= prp_su3_pt[ 9];
    l_pt[10] -= prp_su3_pt[10];
    l_pt[11] -= prp_su3_pt[11];
    l_pt[12] += GAMMA_X_SPIN2_VAL*prp_su3_pt[flav_gamma(GAMMA_X_SPIN2_CO)];
    l_pt[13] += GAMMA_X_SPIN2_VAL*prp_su3_pt[flav_gamma(GAMMA_X_SPIN2_CO)+1];
    l_pt[14] += GAMMA_X_SPIN2_VAL*prp_su3_pt[flav_gamma(GAMMA_X_SPIN2_CO)+2];
    l_pt[15] += GAMMA_X_SPIN3_VAL*prp_su3_pt[flav_gamma(GAMMA_X_SPIN3_CO)];
    l_pt[16] += GAMMA_X_SPIN3_VAL*prp_su3_pt[flav_gamma(GAMMA_X_SPIN3_CO)+1];
    l_pt[17] += GAMMA_X_SPIN3_VAL*prp_su3_pt[flav_gamma(GAMMA_X_SPIN3_CO)+2];
    l_pt[18] += GAMMA_X_SPIN2_VAL*prp_su3_pt[flav_gamma(GAMMA_X_SPIN2_CO)+6];
    l_pt[19] += GAMMA_X_SPIN2_VAL*prp_su3_pt[flav_gamma(GAMMA_X_SPIN2_CO)+7];
    l_pt[20] += GAMMA_X_SPIN2_VAL*prp_su3_pt[flav_gamma(GAMMA_X_SPIN2_CO)+8];
    l_pt[21] += GAMMA_X_SPIN3_VAL*prp_su3_pt[flav_gamma(GAMMA_X_SPIN3_CO)+6];
    l_pt[22] += GAMMA_X_SPIN3_VAL*prp_su3_pt[flav_gamma(GAMMA_X_SPIN3_CO)+7];
    l_pt[23] += GAMMA_X_SPIN3_VAL*prp_su3_pt[flav_gamma(GAMMA_X_SPIN3_CO)+8];
  }

  // -(1 + gamma_X)
  static inline void dpbn_su3_X_PRECISION( const buffer_PRECISION prn_su3_pt, const buffer_PRECISION l_pt ) {
    l_pt[ 0] -= prn_su3_pt[ 0];
    l_pt[ 1] -= prn_su3_pt[ 1];
    l_pt[ 2] -= prn_su3_pt[ 2];
    l_pt[ 3] -= prn_su3_pt[ 3];
    l_pt[ 4] -= prn_su3_pt[ 4];
    l_pt[ 5] -= prn_su3_pt[ 5];
    l_pt[ 6] -= prn_su3_pt[ 6];
    l_pt[ 7] -= prn_su3_pt[ 7];
    l_pt[ 8] -= prn_su3_pt[ 8];
    l_pt[ 9] -= prn_su3_pt[ 9];
    l_pt[10] -= prn_su3_pt[10];
    l_pt[11] -= prn_su3_pt[11];
    l_pt[12] -= GAMMA_X_SPIN2_VAL*prn_su3_pt[flav_gamma(GAMMA_X_SPIN2_CO)];
    l_pt[13] -= GAMMA_X_SPIN2_VAL*prn_su3_pt[flav_gamma(GAMMA_X_SPIN2_CO)+1];
    l_pt[14] -= GAMMA_X_SPIN2_VAL*prn_su3_pt[flav_gamma(GAMMA_X_SPIN2_CO)+2];
    l_pt[15] -= GAMMA_X_SPIN3_VAL*prn_su3_pt[flav_gamma(GAMMA_X_SPIN3_CO)];
    l_pt[16] -= GAMMA_X_SPIN3_VAL*prn_su3_pt[flav_gamma(GAMMA_X_SPIN3_CO)+1];
    l_pt[17] -= GAMMA_X_SPIN3_VAL*prn_su3_pt[flav_gamma(GAMMA_X_SPIN3_CO)+2];
    l_pt[18] -= GAMMA_X_SPIN2_VAL*prn_su3_pt[flav_gamma(GAMMA_X_SPIN2_CO)+6];
    l_pt[19] -= GAMMA_X_SPIN2_VAL*prn_su3_pt[flav_gamma(GAMMA_X_SPIN2_CO)+7];
    l_pt[20] -= GAMMA_X_SPIN2_VAL*prn_su3_pt[flav_gamma(GAMMA_X_SPIN2_CO)+8];
    l_pt[21] -= GAMMA_X_SPIN3_VAL*prn_su3_pt[flav_gamma(GAMMA_X_SPIN3_CO)+6];
    l_pt[22] -= GAMMA_X_SPIN3_VAL*prn_su3_pt[flav_gamma(GAMMA_X_SPIN3_CO)+7];
    l_pt[23] -= GAMMA_X_SPIN3_VAL*prn_su3_pt[flav_gamma(GAMMA_X_SPIN3_CO)+8];
  }

#endif
/************************************* END  **********************************************/

  static inline void twospin_p_T_PRECISION_new( const buffer_PRECISION out_spin0and1, const buffer_PRECISION out_spin2and3, const buffer_PRECISION in, const int n_vect, const int n_vect01, const int n_vect23, const int n_vect_in ) {
    int j, jj;                                                          
    VECTOR_LOOP(j, n_vect, jj,  out_spin0and1[ 0*n_vect01+j+jj] -= in[ 0*n_vect_in+j+jj];
                                out_spin0and1[ 1*n_vect01+j+jj] -= in[ 1*n_vect_in+j+jj];
                                out_spin0and1[ 2*n_vect01+j+jj] -= in[ 2*n_vect_in+j+jj];
                                out_spin0and1[ 3*n_vect01+j+jj] -= in[ 3*n_vect_in+j+jj];
                                out_spin0and1[ 4*n_vect01+j+jj] -= in[ 4*n_vect_in+j+jj];
                                out_spin0and1[ 5*n_vect01+j+jj] -= in[ 5*n_vect_in+j+jj];
                                out_spin0and1[ 6*n_vect01+j+jj] += GAMMA_T_SPIN2_VAL * in[ 3*GAMMA_T_SPIN2_CO   *n_vect_in+j+jj];
                                out_spin0and1[ 7*n_vect01+j+jj] += GAMMA_T_SPIN2_VAL * in[(3*GAMMA_T_SPIN2_CO+1)*n_vect_in+j+jj];
                                out_spin0and1[ 8*n_vect01+j+jj] += GAMMA_T_SPIN2_VAL * in[(3*GAMMA_T_SPIN2_CO+2)*n_vect_in+j+jj];
                                out_spin0and1[ 9*n_vect01+j+jj] += GAMMA_T_SPIN3_VAL * in[ 3*GAMMA_T_SPIN3_CO   *n_vect_in+j+jj ];
                                out_spin0and1[10*n_vect01+j+jj] += GAMMA_T_SPIN3_VAL * in[(3*GAMMA_T_SPIN3_CO+1)*n_vect_in+j+jj];
                                out_spin0and1[11*n_vect01+j+jj] += GAMMA_T_SPIN3_VAL * in[(3*GAMMA_T_SPIN3_CO+2)*n_vect_in+j+jj];
                                out_spin2and3[ 0*n_vect23+j+jj] += GAMMA_T_SPIN0_VAL * in[ 3*GAMMA_T_SPIN0_CO   *n_vect_in+j+jj ];
                                out_spin2and3[ 1*n_vect23+j+jj] += GAMMA_T_SPIN0_VAL * in[(3*GAMMA_T_SPIN0_CO+1)*n_vect_in+j+jj];
                                out_spin2and3[ 2*n_vect23+j+jj] += GAMMA_T_SPIN0_VAL * in[(3*GAMMA_T_SPIN0_CO+2)*n_vect_in+j+jj];
                                out_spin2and3[ 3*n_vect23+j+jj] += GAMMA_T_SPIN1_VAL * in[ 3*GAMMA_T_SPIN1_CO   *n_vect_in+j+jj ];
                                out_spin2and3[ 4*n_vect23+j+jj] += GAMMA_T_SPIN1_VAL * in[(3*GAMMA_T_SPIN1_CO+1)*n_vect_in+j+jj];
                                out_spin2and3[ 5*n_vect23+j+jj] += GAMMA_T_SPIN1_VAL * in[(3*GAMMA_T_SPIN1_CO+2)*n_vect_in+j+jj];
                                out_spin2and3[ 6*n_vect23+j+jj] -= in[ 6*n_vect_in+j+jj];
                                out_spin2and3[ 7*n_vect23+j+jj] -= in[ 7*n_vect_in+j+jj];
                                out_spin2and3[ 8*n_vect23+j+jj] -= in[ 8*n_vect_in+j+jj];
                                out_spin2and3[ 9*n_vect23+j+jj] -= in[ 9*n_vect_in+j+jj];
                                out_spin2and3[10*n_vect23+j+jj] -= in[10*n_vect_in+j+jj];
                                out_spin2and3[11*n_vect23+j+jj] -= in[11*n_vect_in+j+jj];)
  }

  static inline void twospin2_p_T_PRECISION_new( const buffer_PRECISION out_spin0and1, const buffer_PRECISION out_spin2and3, const buffer_PRECISION in, const int n_vect, const int n_vect01, const int n_vect23, const int n_vect_in ) {
    int j, jj;
    VECTOR_LOOP(j, n_vect, jj, out_spin0and1[ 0*n_vect01+j+jj] = in[ 0*n_vect_in+j+jj];
                               out_spin0and1[ 1*n_vect01+j+jj] = in[ 1*n_vect_in+j+jj];
                               out_spin0and1[ 2*n_vect01+j+jj] = in[ 2*n_vect_in+j+jj];
                               out_spin0and1[ 3*n_vect01+j+jj] = in[ 3*n_vect_in+j+jj];
                               out_spin0and1[ 4*n_vect01+j+jj] = in[ 4*n_vect_in+j+jj];
                               out_spin0and1[ 5*n_vect01+j+jj] = in[ 5*n_vect_in+j+jj];
                               out_spin0and1[ 6*n_vect01+j+jj] = -GAMMA_T_SPIN2_VAL * in[ 3*GAMMA_T_SPIN2_CO   *n_vect_in+j+jj];
                               out_spin0and1[ 7*n_vect01+j+jj] = -GAMMA_T_SPIN2_VAL * in[(3*GAMMA_T_SPIN2_CO+1)*n_vect_in+j+jj];
                               out_spin0and1[ 8*n_vect01+j+jj] = -GAMMA_T_SPIN2_VAL * in[(3*GAMMA_T_SPIN2_CO+2)*n_vect_in+j+jj];
                               out_spin0and1[ 9*n_vect01+j+jj] = -GAMMA_T_SPIN3_VAL * in[ 3*GAMMA_T_SPIN3_CO   *n_vect_in+j+jj];
                               out_spin0and1[10*n_vect01+j+jj] = -GAMMA_T_SPIN3_VAL * in[(3*GAMMA_T_SPIN3_CO+1)*n_vect_in+j+jj];
                               out_spin0and1[11*n_vect01+j+jj] = -GAMMA_T_SPIN3_VAL * in[(3*GAMMA_T_SPIN3_CO+2)*n_vect_in+j+jj];
                               out_spin2and3[ 0*n_vect23+j+jj] = -GAMMA_T_SPIN0_VAL * in[ 3*GAMMA_T_SPIN0_CO   *n_vect_in+j+jj];
                               out_spin2and3[ 1*n_vect23+j+jj] = -GAMMA_T_SPIN0_VAL * in[(3*GAMMA_T_SPIN0_CO+1)*n_vect_in+j+jj];
                               out_spin2and3[ 2*n_vect23+j+jj] = -GAMMA_T_SPIN0_VAL * in[(3*GAMMA_T_SPIN0_CO+2)*n_vect_in+j+jj];
                               out_spin2and3[ 3*n_vect23+j+jj] = -GAMMA_T_SPIN1_VAL * in[ 3*GAMMA_T_SPIN1_CO   *n_vect_in+j+jj];
                               out_spin2and3[ 4*n_vect23+j+jj] = -GAMMA_T_SPIN1_VAL * in[(3*GAMMA_T_SPIN1_CO+1)*n_vect_in+j+jj];
                               out_spin2and3[ 5*n_vect23+j+jj] = -GAMMA_T_SPIN1_VAL * in[(3*GAMMA_T_SPIN1_CO+2)*n_vect_in+j+jj];
                               out_spin2and3[ 6*n_vect23+j+jj] = in[ 6*n_vect_in+j+jj];
                               out_spin2and3[ 7*n_vect23+j+jj] = in[ 7*n_vect_in+j+jj];
                               out_spin2and3[ 8*n_vect23+j+jj] = in[ 8*n_vect_in+j+jj];
                               out_spin2and3[ 9*n_vect23+j+jj] = in[ 9*n_vect_in+j+jj];
                               out_spin2and3[10*n_vect23+j+jj] = in[10*n_vect_in+j+jj];
                               out_spin2and3[11*n_vect23+j+jj] = in[11*n_vect_in+j+jj];)
  }

  static inline void twospin_n_T_PRECISION_new( const buffer_PRECISION out_spin0and1, const buffer_PRECISION out_spin2and3, const buffer_PRECISION in, const int n_vect, const int n_vect01, const int n_vect23, const int n_vect_in ) {
    int j, jj;
    VECTOR_LOOP(j, n_vect, jj,  out_spin0and1[ 0*n_vect01+j+jj] -= in[ 0*n_vect_in+j+jj];
                                out_spin0and1[ 1*n_vect01+j+jj] -= in[ 1*n_vect_in+j+jj];
                                out_spin0and1[ 2*n_vect01+j+jj] -= in[ 2*n_vect_in+j+jj];
                                out_spin0and1[ 3*n_vect01+j+jj] -= in[ 3*n_vect_in+j+jj];
                                out_spin0and1[ 4*n_vect01+j+jj] -= in[ 4*n_vect_in+j+jj];
                                out_spin0and1[ 5*n_vect01+j+jj] -= in[ 5*n_vect_in+j+jj];
                                out_spin0and1[ 6*n_vect01+j+jj] -= GAMMA_T_SPIN2_VAL * in[ 3*GAMMA_T_SPIN2_CO   *n_vect_in+j+jj];
                                out_spin0and1[ 7*n_vect01+j+jj] -= GAMMA_T_SPIN2_VAL * in[(3*GAMMA_T_SPIN2_CO+1)*n_vect_in+j+jj];
                                out_spin0and1[ 8*n_vect01+j+jj] -= GAMMA_T_SPIN2_VAL * in[(3*GAMMA_T_SPIN2_CO+2)*n_vect_in+j+jj];
                                out_spin0and1[ 9*n_vect01+j+jj] -= GAMMA_T_SPIN3_VAL * in[ 3*GAMMA_T_SPIN3_CO   *n_vect_in+j+jj];
                                out_spin0and1[10*n_vect01+j+jj] -= GAMMA_T_SPIN3_VAL * in[(3*GAMMA_T_SPIN3_CO+1)*n_vect_in+j+jj];
                                out_spin0and1[11*n_vect01+j+jj] -= GAMMA_T_SPIN3_VAL * in[(3*GAMMA_T_SPIN3_CO+2)*n_vect_in+j+jj];
                                out_spin2and3[ 0*n_vect23+j+jj] -= GAMMA_T_SPIN0_VAL * in[ 3*GAMMA_T_SPIN0_CO   *n_vect_in+j+jj];
                                out_spin2and3[ 1*n_vect23+j+jj] -= GAMMA_T_SPIN0_VAL * in[(3*GAMMA_T_SPIN0_CO+1)*n_vect_in+j+jj];
                                out_spin2and3[ 2*n_vect23+j+jj] -= GAMMA_T_SPIN0_VAL * in[(3*GAMMA_T_SPIN0_CO+2)*n_vect_in+j+jj];
                                out_spin2and3[ 3*n_vect23+j+jj] -= GAMMA_T_SPIN1_VAL * in[ 3*GAMMA_T_SPIN1_CO   *n_vect_in+j+jj];
                                out_spin2and3[ 4*n_vect23+j+jj] -= GAMMA_T_SPIN1_VAL * in[(3*GAMMA_T_SPIN1_CO+1)*n_vect_in+j+jj];
                                out_spin2and3[ 5*n_vect23+j+jj] -= GAMMA_T_SPIN1_VAL * in[(3*GAMMA_T_SPIN1_CO+2)*n_vect_in+j+jj];
                                out_spin2and3[ 6*n_vect23+j+jj] -= in[ 6*n_vect_in+j+jj];
                                out_spin2and3[ 7*n_vect23+j+jj] -= in[ 7*n_vect_in+j+jj];
                                out_spin2and3[ 8*n_vect23+j+jj] -= in[ 8*n_vect_in+j+jj];
                                out_spin2and3[ 9*n_vect23+j+jj] -= in[ 9*n_vect_in+j+jj];
                                out_spin2and3[10*n_vect23+j+jj] -= in[10*n_vect_in+j+jj];
                                out_spin2and3[11*n_vect23+j+jj] -= in[11*n_vect_in+j+jj];)
  }

  static inline void twospin_p_Z_PRECISION_new( const buffer_PRECISION out_spin0and1, const buffer_PRECISION out_spin2and3, const buffer_PRECISION in, const int n_vect, const int n_vect01, const int n_vect23, const int n_vect_in ) {
    int j, jj;
    VECTOR_LOOP(j, n_vect, jj,  out_spin0and1[ 0*n_vect01+j+jj] -= in[ 0*n_vect_in+j+jj];
                                out_spin0and1[ 1*n_vect01+j+jj] -= in[ 1*n_vect_in+j+jj];
                                out_spin0and1[ 2*n_vect01+j+jj] -= in[ 2*n_vect_in+j+jj];
                                out_spin0and1[ 3*n_vect01+j+jj] -= in[ 3*n_vect_in+j+jj];
                                out_spin0and1[ 4*n_vect01+j+jj] -= in[ 4*n_vect_in+j+jj];
                                out_spin0and1[ 5*n_vect01+j+jj] -= in[ 5*n_vect_in+j+jj];
                                out_spin0and1[ 6*n_vect01+j+jj] += GAMMA_Z_SPIN2_VAL * in[ 3*GAMMA_Z_SPIN2_CO   *n_vect_in+j+jj];
                                out_spin0and1[ 7*n_vect01+j+jj] += GAMMA_Z_SPIN2_VAL * in[(3*GAMMA_Z_SPIN2_CO+1)*n_vect_in+j+jj];
                                out_spin0and1[ 8*n_vect01+j+jj] += GAMMA_Z_SPIN2_VAL * in[(3*GAMMA_Z_SPIN2_CO+2)*n_vect_in+j+jj];
                                out_spin0and1[ 9*n_vect01+j+jj] += GAMMA_Z_SPIN3_VAL * in[ 3*GAMMA_Z_SPIN3_CO   *n_vect_in+j+jj];
                                out_spin0and1[10*n_vect01+j+jj] += GAMMA_Z_SPIN3_VAL * in[(3*GAMMA_Z_SPIN3_CO+1)*n_vect_in+j+jj];
                                out_spin0and1[11*n_vect01+j+jj] += GAMMA_Z_SPIN3_VAL * in[(3*GAMMA_Z_SPIN3_CO+2)*n_vect_in+j+jj];
                                out_spin2and3[ 0*n_vect23+j+jj] += GAMMA_Z_SPIN0_VAL * in[ 3*GAMMA_Z_SPIN0_CO   *n_vect_in+j+jj];
                                out_spin2and3[ 1*n_vect23+j+jj] += GAMMA_Z_SPIN0_VAL * in[(3*GAMMA_Z_SPIN0_CO+1)*n_vect_in+j+jj];
                                out_spin2and3[ 2*n_vect23+j+jj] += GAMMA_Z_SPIN0_VAL * in[(3*GAMMA_Z_SPIN0_CO+2)*n_vect_in+j+jj];
                                out_spin2and3[ 3*n_vect23+j+jj] += GAMMA_Z_SPIN1_VAL * in[ 3*GAMMA_Z_SPIN1_CO   *n_vect_in+j+jj];
                                out_spin2and3[ 4*n_vect23+j+jj] += GAMMA_Z_SPIN1_VAL * in[(3*GAMMA_Z_SPIN1_CO+1)*n_vect_in+j+jj];
                                out_spin2and3[ 5*n_vect23+j+jj] += GAMMA_Z_SPIN1_VAL * in[(3*GAMMA_Z_SPIN1_CO+2)*n_vect_in+j+jj];
                                out_spin2and3[ 6*n_vect23+j+jj] -= in[ 6*n_vect_in+j+jj];
                                out_spin2and3[ 7*n_vect23+j+jj] -= in[ 7*n_vect_in+j+jj];
                                out_spin2and3[ 8*n_vect23+j+jj] -= in[ 8*n_vect_in+j+jj];
                                out_spin2and3[ 9*n_vect23+j+jj] -= in[ 9*n_vect_in+j+jj];
                                out_spin2and3[10*n_vect23+j+jj] -= in[10*n_vect_in+j+jj];
                                out_spin2and3[11*n_vect23+j+jj] -= in[11*n_vect_in+j+jj];)
  }

  static inline void twospin2_p_Z_PRECISION_new( const buffer_PRECISION out_spin0and1, const buffer_PRECISION out_spin2and3, const buffer_PRECISION in, const int n_vect, const int n_vect01, const int n_vect23, const int n_vect_in ) {
    int j, jj;
    VECTOR_LOOP(j, n_vect, jj,  out_spin0and1[ 0*n_vect01+j+jj] = in[ 0*n_vect_in+j+jj];
                                out_spin0and1[ 1*n_vect01+j+jj] = in[ 1*n_vect_in+j+jj];
                                out_spin0and1[ 2*n_vect01+j+jj] = in[ 2*n_vect_in+j+jj];
                                out_spin0and1[ 3*n_vect01+j+jj] = in[ 3*n_vect_in+j+jj];
                                out_spin0and1[ 4*n_vect01+j+jj] = in[ 4*n_vect_in+j+jj];
                                out_spin0and1[ 5*n_vect01+j+jj] = in[ 5*n_vect_in+j+jj];
                                out_spin0and1[ 6*n_vect01+j+jj] = -GAMMA_Z_SPIN2_VAL * in[ 3*GAMMA_Z_SPIN2_CO   *n_vect_in+j+jj];
                                out_spin0and1[ 7*n_vect01+j+jj] = -GAMMA_Z_SPIN2_VAL * in[(3*GAMMA_Z_SPIN2_CO+1)*n_vect_in+j+jj];
                                out_spin0and1[ 8*n_vect01+j+jj] = -GAMMA_Z_SPIN2_VAL * in[(3*GAMMA_Z_SPIN2_CO+2)*n_vect_in+j+jj];
                                out_spin0and1[ 9*n_vect01+j+jj] = -GAMMA_Z_SPIN3_VAL * in[ 3*GAMMA_Z_SPIN3_CO   *n_vect_in+j+jj];
                                out_spin0and1[10*n_vect01+j+jj] = -GAMMA_Z_SPIN3_VAL * in[(3*GAMMA_Z_SPIN3_CO+1)*n_vect_in+j+jj];
                                out_spin0and1[11*n_vect01+j+jj] = -GAMMA_Z_SPIN3_VAL * in[(3*GAMMA_Z_SPIN3_CO+2)*n_vect_in+j+jj];
                                out_spin2and3[ 0*n_vect23+j+jj] = -GAMMA_Z_SPIN0_VAL * in[ 3*GAMMA_Z_SPIN0_CO   *n_vect_in+j+jj];
                                out_spin2and3[ 1*n_vect23+j+jj] = -GAMMA_Z_SPIN0_VAL * in[(3*GAMMA_Z_SPIN0_CO+1)*n_vect_in+j+jj];
                                out_spin2and3[ 2*n_vect23+j+jj] = -GAMMA_Z_SPIN0_VAL * in[(3*GAMMA_Z_SPIN0_CO+2)*n_vect_in+j+jj];
                                out_spin2and3[ 3*n_vect23+j+jj] = -GAMMA_Z_SPIN1_VAL * in[ 3*GAMMA_Z_SPIN1_CO   *n_vect_in+j+jj];
                                out_spin2and3[ 4*n_vect23+j+jj] = -GAMMA_Z_SPIN1_VAL * in[(3*GAMMA_Z_SPIN1_CO+1)*n_vect_in+j+jj];
                                out_spin2and3[ 5*n_vect23+j+jj] = -GAMMA_Z_SPIN1_VAL * in[(3*GAMMA_Z_SPIN1_CO+2)*n_vect_in+j+jj];
                                out_spin2and3[ 6*n_vect23+j+jj] = in[ 6*n_vect_in+j+jj];
                                out_spin2and3[ 7*n_vect23+j+jj] = in[ 7*n_vect_in+j+jj];
                                out_spin2and3[ 8*n_vect23+j+jj] = in[ 8*n_vect_in+j+jj];
                                out_spin2and3[ 9*n_vect23+j+jj] = in[ 9*n_vect_in+j+jj];
                                out_spin2and3[10*n_vect23+j+jj] = in[10*n_vect_in+j+jj];
                                out_spin2and3[11*n_vect23+j+jj] = in[11*n_vect_in+j+jj];)
  }

  static inline void twospin_n_Z_PRECISION_new( const buffer_PRECISION out_spin0and1, const buffer_PRECISION out_spin2and3, const buffer_PRECISION in, const int n_vect, const int n_vect01, const int n_vect23, const int n_vect_in ) {
    int j, jj;
    VECTOR_LOOP(j, n_vect, jj,  out_spin0and1[ 0*n_vect01+j+jj] -= in[ 0*n_vect_in+j+jj];
                                out_spin0and1[ 1*n_vect01+j+jj] -= in[ 1*n_vect_in+j+jj];
                                out_spin0and1[ 2*n_vect01+j+jj] -= in[ 2*n_vect_in+j+jj];
                                out_spin0and1[ 3*n_vect01+j+jj] -= in[ 3*n_vect_in+j+jj];
                                out_spin0and1[ 4*n_vect01+j+jj] -= in[ 4*n_vect_in+j+jj];
                                out_spin0and1[ 5*n_vect01+j+jj] -= in[ 5*n_vect_in+j+jj];
                                out_spin0and1[ 6*n_vect01+j+jj] -= GAMMA_Z_SPIN2_VAL * in[ 3*GAMMA_Z_SPIN2_CO   *n_vect_in+j+jj];
                                out_spin0and1[ 7*n_vect01+j+jj] -= GAMMA_Z_SPIN2_VAL * in[(3*GAMMA_Z_SPIN2_CO+1)*n_vect_in+j+jj];
                                out_spin0and1[ 8*n_vect01+j+jj] -= GAMMA_Z_SPIN2_VAL * in[(3*GAMMA_Z_SPIN2_CO+2)*n_vect_in+j+jj];
                                out_spin0and1[ 9*n_vect01+j+jj] -= GAMMA_Z_SPIN3_VAL * in[ 3*GAMMA_Z_SPIN3_CO   *n_vect_in+j+jj];
                                out_spin0and1[10*n_vect01+j+jj] -= GAMMA_Z_SPIN3_VAL * in[(3*GAMMA_Z_SPIN3_CO+1)*n_vect_in+j+jj];
                                out_spin0and1[11*n_vect01+j+jj] -= GAMMA_Z_SPIN3_VAL * in[(3*GAMMA_Z_SPIN3_CO+2)*n_vect_in+j+jj];
                                out_spin2and3[ 0*n_vect23+j+jj] -= GAMMA_Z_SPIN0_VAL * in[ 3*GAMMA_Z_SPIN0_CO   *n_vect_in+j+jj];
                                out_spin2and3[ 1*n_vect23+j+jj] -= GAMMA_Z_SPIN0_VAL * in[(3*GAMMA_Z_SPIN0_CO+1)*n_vect_in+j+jj];
                                out_spin2and3[ 2*n_vect23+j+jj] -= GAMMA_Z_SPIN0_VAL * in[(3*GAMMA_Z_SPIN0_CO+2)*n_vect_in+j+jj];
                                out_spin2and3[ 3*n_vect23+j+jj] -= GAMMA_Z_SPIN1_VAL * in[ 3*GAMMA_Z_SPIN1_CO   *n_vect_in+j+jj];
                                out_spin2and3[ 4*n_vect23+j+jj] -= GAMMA_Z_SPIN1_VAL * in[(3*GAMMA_Z_SPIN1_CO+1)*n_vect_in+j+jj];
                                out_spin2and3[ 5*n_vect23+j+jj] -= GAMMA_Z_SPIN1_VAL * in[(3*GAMMA_Z_SPIN1_CO+2)*n_vect_in+j+jj];
                                out_spin2and3[ 6*n_vect23+j+jj] -= in[ 6*n_vect_in+j+jj];
                                out_spin2and3[ 7*n_vect23+j+jj] -= in[ 7*n_vect_in+j+jj];
                                out_spin2and3[ 8*n_vect23+j+jj] -= in[ 8*n_vect_in+j+jj];
                                out_spin2and3[ 9*n_vect23+j+jj] -= in[ 9*n_vect_in+j+jj];
                                out_spin2and3[10*n_vect23+j+jj] -= in[10*n_vect_in+j+jj];
                                out_spin2and3[11*n_vect23+j+jj] -= in[11*n_vect_in+j+jj];)
  }

  static inline void twospin_p_Y_PRECISION_new( const buffer_PRECISION out_spin0and1, const buffer_PRECISION out_spin2and3, const buffer_PRECISION in, const int n_vect, const int n_vect01, const int n_vect23, const int n_vect_in ) {
    int j, jj;
    VECTOR_LOOP(j, n_vect, jj,  out_spin0and1[ 0*n_vect01+j+jj] -= in[ 0*n_vect_in+j+jj];
                                out_spin0and1[ 1*n_vect01+j+jj] -= in[ 1*n_vect_in+j+jj];
                                out_spin0and1[ 2*n_vect01+j+jj] -= in[ 2*n_vect_in+j+jj];
                                out_spin0and1[ 3*n_vect01+j+jj] -= in[ 3*n_vect_in+j+jj];
                                out_spin0and1[ 4*n_vect01+j+jj] -= in[ 4*n_vect_in+j+jj];
                                out_spin0and1[ 5*n_vect01+j+jj] -= in[ 5*n_vect_in+j+jj];
                                out_spin0and1[ 6*n_vect01+j+jj] += GAMMA_Y_SPIN2_VAL * in[ 3*GAMMA_Y_SPIN2_CO   *n_vect_in+j+jj];
                                out_spin0and1[ 7*n_vect01+j+jj] += GAMMA_Y_SPIN2_VAL * in[(3*GAMMA_Y_SPIN2_CO+1)*n_vect_in+j+jj];
                                out_spin0and1[ 8*n_vect01+j+jj] += GAMMA_Y_SPIN2_VAL * in[(3*GAMMA_Y_SPIN2_CO+2)*n_vect_in+j+jj];
                                out_spin0and1[ 9*n_vect01+j+jj] += GAMMA_Y_SPIN3_VAL * in[ 3*GAMMA_Y_SPIN3_CO   *n_vect_in+j+jj];
                                out_spin0and1[10*n_vect01+j+jj] += GAMMA_Y_SPIN3_VAL * in[(3*GAMMA_Y_SPIN3_CO+1)*n_vect_in+j+jj];
                                out_spin0and1[11*n_vect01+j+jj] += GAMMA_Y_SPIN3_VAL * in[(3*GAMMA_Y_SPIN3_CO+2)*n_vect_in+j+jj];
                                out_spin2and3[ 0*n_vect23+j+jj] += GAMMA_Y_SPIN0_VAL * in[ 3*GAMMA_Y_SPIN0_CO   *n_vect_in+j+jj];
                                out_spin2and3[ 1*n_vect23+j+jj] += GAMMA_Y_SPIN0_VAL * in[(3*GAMMA_Y_SPIN0_CO+1)*n_vect_in+j+jj];
                                out_spin2and3[ 2*n_vect23+j+jj] += GAMMA_Y_SPIN0_VAL * in[(3*GAMMA_Y_SPIN0_CO+2)*n_vect_in+j+jj];
                                out_spin2and3[ 3*n_vect23+j+jj] += GAMMA_Y_SPIN1_VAL * in[ 3*GAMMA_Y_SPIN1_CO   *n_vect_in+j+jj];
                                out_spin2and3[ 4*n_vect23+j+jj] += GAMMA_Y_SPIN1_VAL * in[(3*GAMMA_Y_SPIN1_CO+1)*n_vect_in+j+jj];
                                out_spin2and3[ 5*n_vect23+j+jj] += GAMMA_Y_SPIN1_VAL * in[(3*GAMMA_Y_SPIN1_CO+2)*n_vect_in+j+jj];
                                out_spin2and3[ 6*n_vect23+j+jj] -= in[ 6*n_vect_in+j+jj];
                                out_spin2and3[ 7*n_vect23+j+jj] -= in[ 7*n_vect_in+j+jj];
                                out_spin2and3[ 8*n_vect23+j+jj] -= in[ 8*n_vect_in+j+jj];
                                out_spin2and3[ 9*n_vect23+j+jj] -= in[ 9*n_vect_in+j+jj];
                                out_spin2and3[10*n_vect23+j+jj] -= in[10*n_vect_in+j+jj];
                                out_spin2and3[11*n_vect23+j+jj] -= in[11*n_vect_in+j+jj];)
  }

  static inline void twospin2_p_Y_PRECISION_new( const buffer_PRECISION out_spin0and1, const buffer_PRECISION out_spin2and3, const buffer_PRECISION in, const int n_vect, const int n_vect01, const int n_vect23, const int n_vect_in ) {
    int j, jj;
    VECTOR_LOOP(j, n_vect, jj,  out_spin0and1[ 0*n_vect01+j+jj] = in[ 0*n_vect_in+j+jj];
                                out_spin0and1[ 1*n_vect01+j+jj] = in[ 1*n_vect_in+j+jj];
                                out_spin0and1[ 2*n_vect01+j+jj] = in[ 2*n_vect_in+j+jj];
                                out_spin0and1[ 3*n_vect01+j+jj] = in[ 3*n_vect_in+j+jj];
                                out_spin0and1[ 4*n_vect01+j+jj] = in[ 4*n_vect_in+j+jj];
                                out_spin0and1[ 5*n_vect01+j+jj] = in[ 5*n_vect_in+j+jj];
                                out_spin0and1[ 6*n_vect01+j+jj] = -GAMMA_Y_SPIN2_VAL * in[ 3*GAMMA_Y_SPIN2_CO   *n_vect_in+j+jj];
                                out_spin0and1[ 7*n_vect01+j+jj] = -GAMMA_Y_SPIN2_VAL * in[(3*GAMMA_Y_SPIN2_CO+1)*n_vect_in+j+jj];
                                out_spin0and1[ 8*n_vect01+j+jj] = -GAMMA_Y_SPIN2_VAL * in[(3*GAMMA_Y_SPIN2_CO+2)*n_vect_in+j+jj];
                                out_spin0and1[ 9*n_vect01+j+jj] = -GAMMA_Y_SPIN3_VAL * in[ 3*GAMMA_Y_SPIN3_CO   *n_vect_in+j+jj];
                                out_spin0and1[10*n_vect01+j+jj] = -GAMMA_Y_SPIN3_VAL * in[(3*GAMMA_Y_SPIN3_CO+1)*n_vect_in+j+jj];
                                out_spin0and1[11*n_vect01+j+jj] = -GAMMA_Y_SPIN3_VAL * in[(3*GAMMA_Y_SPIN3_CO+2)*n_vect_in+j+jj];
                                out_spin2and3[ 0*n_vect23+j+jj] = -GAMMA_Y_SPIN0_VAL * in[ 3*GAMMA_Y_SPIN0_CO   *n_vect_in+j+jj];
                                out_spin2and3[ 1*n_vect23+j+jj] = -GAMMA_Y_SPIN0_VAL * in[(3*GAMMA_Y_SPIN0_CO+1)*n_vect_in+j+jj];
                                out_spin2and3[ 2*n_vect23+j+jj] = -GAMMA_Y_SPIN0_VAL * in[(3*GAMMA_Y_SPIN0_CO+2)*n_vect_in+j+jj];
                                out_spin2and3[ 3*n_vect23+j+jj] = -GAMMA_Y_SPIN1_VAL * in[ 3*GAMMA_Y_SPIN1_CO   *n_vect_in+j+jj];
                                out_spin2and3[ 4*n_vect23+j+jj] = -GAMMA_Y_SPIN1_VAL * in[(3*GAMMA_Y_SPIN1_CO+1)*n_vect_in+j+jj];
                                out_spin2and3[ 5*n_vect23+j+jj] = -GAMMA_Y_SPIN1_VAL * in[(3*GAMMA_Y_SPIN1_CO+2)*n_vect_in+j+jj];
                                out_spin2and3[ 6*n_vect23+j+jj] = in[ 6*n_vect_in+j+jj];
                                out_spin2and3[ 7*n_vect23+j+jj] = in[ 7*n_vect_in+j+jj];
                                out_spin2and3[ 8*n_vect23+j+jj] = in[ 8*n_vect_in+j+jj];
                                out_spin2and3[ 9*n_vect23+j+jj] = in[ 9*n_vect_in+j+jj];
                                out_spin2and3[10*n_vect23+j+jj] = in[10*n_vect_in+j+jj];
                                out_spin2and3[11*n_vect23+j+jj] = in[11*n_vect_in+j+jj];)
  }

  static inline void twospin_n_Y_PRECISION_new( const buffer_PRECISION out_spin0and1, const buffer_PRECISION out_spin2and3, const buffer_PRECISION in, const int n_vect, const int n_vect01, const int n_vect23, const int n_vect_in ) {
    int j, jj;
    VECTOR_LOOP(j, n_vect, jj,  out_spin0and1[ 0*n_vect01+j+jj] -= in[ 0*n_vect_in+j+jj];
                                out_spin0and1[ 1*n_vect01+j+jj] -= in[ 1*n_vect_in+j+jj];
                                out_spin0and1[ 2*n_vect01+j+jj] -= in[ 2*n_vect_in+j+jj];
                                out_spin0and1[ 3*n_vect01+j+jj] -= in[ 3*n_vect_in+j+jj];
                                out_spin0and1[ 4*n_vect01+j+jj] -= in[ 4*n_vect_in+j+jj];
                                out_spin0and1[ 5*n_vect01+j+jj] -= in[ 5*n_vect_in+j+jj];
                                out_spin0and1[ 6*n_vect01+j+jj] -= GAMMA_Y_SPIN2_VAL * in[ 3*GAMMA_Y_SPIN2_CO   *n_vect_in+j+jj];
                                out_spin0and1[ 7*n_vect01+j+jj] -= GAMMA_Y_SPIN2_VAL * in[(3*GAMMA_Y_SPIN2_CO+1)*n_vect_in+j+jj];
                                out_spin0and1[ 8*n_vect01+j+jj] -= GAMMA_Y_SPIN2_VAL * in[(3*GAMMA_Y_SPIN2_CO+2)*n_vect_in+j+jj];
                                out_spin0and1[ 9*n_vect01+j+jj] -= GAMMA_Y_SPIN3_VAL * in[ 3*GAMMA_Y_SPIN3_CO   *n_vect_in+j+jj];
                                out_spin0and1[10*n_vect01+j+jj] -= GAMMA_Y_SPIN3_VAL * in[(3*GAMMA_Y_SPIN3_CO+1)*n_vect_in+j+jj];
                                out_spin0and1[11*n_vect01+j+jj] -= GAMMA_Y_SPIN3_VAL * in[(3*GAMMA_Y_SPIN3_CO+2)*n_vect_in+j+jj];
                                out_spin2and3[ 0*n_vect23+j+jj] -= GAMMA_Y_SPIN0_VAL * in[ 3*GAMMA_Y_SPIN0_CO   *n_vect_in+j+jj];
                                out_spin2and3[ 1*n_vect23+j+jj] -= GAMMA_Y_SPIN0_VAL * in[(3*GAMMA_Y_SPIN0_CO+1)*n_vect_in+j+jj];
                                out_spin2and3[ 2*n_vect23+j+jj] -= GAMMA_Y_SPIN0_VAL * in[(3*GAMMA_Y_SPIN0_CO+2)*n_vect_in+j+jj];
                                out_spin2and3[ 3*n_vect23+j+jj] -= GAMMA_Y_SPIN1_VAL * in[ 3*GAMMA_Y_SPIN1_CO   *n_vect_in+j+jj];
                                out_spin2and3[ 4*n_vect23+j+jj] -= GAMMA_Y_SPIN1_VAL * in[(3*GAMMA_Y_SPIN1_CO+1)*n_vect_in+j+jj];
                                out_spin2and3[ 5*n_vect23+j+jj] -= GAMMA_Y_SPIN1_VAL * in[(3*GAMMA_Y_SPIN1_CO+2)*n_vect_in+j+jj];
                                out_spin2and3[ 6*n_vect23+j+jj] -= in[ 6*n_vect_in+j+jj];
                                out_spin2and3[ 7*n_vect23+j+jj] -= in[ 7*n_vect_in+j+jj];
                                out_spin2and3[ 8*n_vect23+j+jj] -= in[ 8*n_vect_in+j+jj];
                                out_spin2and3[ 9*n_vect23+j+jj] -= in[ 9*n_vect_in+j+jj];
                                out_spin2and3[10*n_vect23+j+jj] -= in[10*n_vect_in+j+jj];
                                out_spin2and3[11*n_vect23+j+jj] -= in[11*n_vect_in+j+jj];)
  }

  static inline void twospin_p_X_PRECISION_new( const buffer_PRECISION out_spin0and1, const buffer_PRECISION out_spin2and3, const buffer_PRECISION in, const int n_vect, const int n_vect01, const int n_vect23, const int n_vect_in ) {
    int j, jj;
    VECTOR_LOOP(j, n_vect, jj,  out_spin0and1[ 0*n_vect01+j+jj] -= in[ 0*n_vect_in+j+jj];
                                out_spin0and1[ 1*n_vect01+j+jj] -= in[ 1*n_vect_in+j+jj];
                                out_spin0and1[ 2*n_vect01+j+jj] -= in[ 2*n_vect_in+j+jj];
                                out_spin0and1[ 3*n_vect01+j+jj] -= in[ 3*n_vect_in+j+jj];
                                out_spin0and1[ 4*n_vect01+j+jj] -= in[ 4*n_vect_in+j+jj];
                                out_spin0and1[ 5*n_vect01+j+jj] -= in[ 5*n_vect_in+j+jj];
                                out_spin0and1[ 6*n_vect01+j+jj] += GAMMA_X_SPIN2_VAL * in[ 3*GAMMA_X_SPIN2_CO   *n_vect_in+j+jj];
                                out_spin0and1[ 7*n_vect01+j+jj] += GAMMA_X_SPIN2_VAL * in[(3*GAMMA_X_SPIN2_CO+1)*n_vect_in+j+jj];
                                out_spin0and1[ 8*n_vect01+j+jj] += GAMMA_X_SPIN2_VAL * in[(3*GAMMA_X_SPIN2_CO+2)*n_vect_in+j+jj];
                                out_spin0and1[ 9*n_vect01+j+jj] += GAMMA_X_SPIN3_VAL * in[ 3*GAMMA_X_SPIN3_CO   *n_vect_in+j+jj];
                                out_spin0and1[10*n_vect01+j+jj] += GAMMA_X_SPIN3_VAL * in[(3*GAMMA_X_SPIN3_CO+1)*n_vect_in+j+jj];
                                out_spin0and1[11*n_vect01+j+jj] += GAMMA_X_SPIN3_VAL * in[(3*GAMMA_X_SPIN3_CO+2)*n_vect_in+j+jj];
                                out_spin2and3[ 0*n_vect23+j+jj] += GAMMA_X_SPIN0_VAL * in[ 3*GAMMA_X_SPIN0_CO   *n_vect_in+j+jj];
                                out_spin2and3[ 1*n_vect23+j+jj] += GAMMA_X_SPIN0_VAL * in[(3*GAMMA_X_SPIN0_CO+1)*n_vect_in+j+jj];
                                out_spin2and3[ 2*n_vect23+j+jj] += GAMMA_X_SPIN0_VAL * in[(3*GAMMA_X_SPIN0_CO+2)*n_vect_in+j+jj];
                                out_spin2and3[ 3*n_vect23+j+jj] += GAMMA_X_SPIN1_VAL * in[ 3*GAMMA_X_SPIN1_CO   *n_vect_in+j+jj];
                                out_spin2and3[ 4*n_vect23+j+jj] += GAMMA_X_SPIN1_VAL * in[(3*GAMMA_X_SPIN1_CO+1)*n_vect_in+j+jj];
                                out_spin2and3[ 5*n_vect23+j+jj] += GAMMA_X_SPIN1_VAL * in[(3*GAMMA_X_SPIN1_CO+2)*n_vect_in+j+jj];
                                out_spin2and3[ 6*n_vect23+j+jj] -= in[ 6*n_vect_in+j+jj];
                                out_spin2and3[ 7*n_vect23+j+jj] -= in[ 7*n_vect_in+j+jj];
                                out_spin2and3[ 8*n_vect23+j+jj] -= in[ 8*n_vect_in+j+jj];
                                out_spin2and3[ 9*n_vect23+j+jj] -= in[ 9*n_vect_in+j+jj];
                                out_spin2and3[10*n_vect23+j+jj] -= in[10*n_vect_in+j+jj];
                                out_spin2and3[11*n_vect23+j+jj] -= in[11*n_vect_in+j+jj];)
  }

  static inline void twospin2_p_X_PRECISION_new( const buffer_PRECISION out_spin0and1, const buffer_PRECISION out_spin2and3, const buffer_PRECISION in, const int n_vect, const int n_vect01, const int n_vect23, const int n_vect_in ) {
    int j, jj;
    VECTOR_LOOP(j, n_vect, jj,  out_spin0and1[ 0*n_vect01+j+jj] = in[ 0*n_vect_in+j+jj];
                                out_spin0and1[ 1*n_vect01+j+jj] = in[ 1*n_vect_in+j+jj];
                                out_spin0and1[ 2*n_vect01+j+jj] = in[ 2*n_vect_in+j+jj];
                                out_spin0and1[ 3*n_vect01+j+jj] = in[ 3*n_vect_in+j+jj];
                                out_spin0and1[ 4*n_vect01+j+jj] = in[ 4*n_vect_in+j+jj];
                                out_spin0and1[ 5*n_vect01+j+jj] = in[ 5*n_vect_in+j+jj];
                                out_spin0and1[ 6*n_vect01+j+jj] = -GAMMA_X_SPIN2_VAL * in[ 3*GAMMA_X_SPIN2_CO   *n_vect_in+j+jj];
                                out_spin0and1[ 7*n_vect01+j+jj] = -GAMMA_X_SPIN2_VAL * in[(3*GAMMA_X_SPIN2_CO+1)*n_vect_in+j+jj];
                                out_spin0and1[ 8*n_vect01+j+jj] = -GAMMA_X_SPIN2_VAL * in[(3*GAMMA_X_SPIN2_CO+2)*n_vect_in+j+jj];
                                out_spin0and1[ 9*n_vect01+j+jj] = -GAMMA_X_SPIN3_VAL * in[ 3*GAMMA_X_SPIN3_CO   *n_vect_in+j+jj];
                                out_spin0and1[10*n_vect01+j+jj] = -GAMMA_X_SPIN3_VAL * in[(3*GAMMA_X_SPIN3_CO+1)*n_vect_in+j+jj];
                                out_spin0and1[11*n_vect01+j+jj] = -GAMMA_X_SPIN3_VAL * in[(3*GAMMA_X_SPIN3_CO+2)*n_vect_in+j+jj];
                                out_spin2and3[ 0*n_vect23+j+jj] = -GAMMA_X_SPIN0_VAL * in[ 3*GAMMA_X_SPIN0_CO   *n_vect_in+j+jj];
                                out_spin2and3[ 1*n_vect23+j+jj] = -GAMMA_X_SPIN0_VAL * in[(3*GAMMA_X_SPIN0_CO+1)*n_vect_in+j+jj];
                                out_spin2and3[ 2*n_vect23+j+jj] = -GAMMA_X_SPIN0_VAL * in[(3*GAMMA_X_SPIN0_CO+2)*n_vect_in+j+jj];
                                out_spin2and3[ 3*n_vect23+j+jj] = -GAMMA_X_SPIN1_VAL * in[ 3*GAMMA_X_SPIN1_CO   *n_vect_in+j+jj];
                                out_spin2and3[ 4*n_vect23+j+jj] = -GAMMA_X_SPIN1_VAL * in[(3*GAMMA_X_SPIN1_CO+1)*n_vect_in+j+jj];
                                out_spin2and3[ 5*n_vect23+j+jj] = -GAMMA_X_SPIN1_VAL * in[(3*GAMMA_X_SPIN1_CO+2)*n_vect_in+j+jj];
                                out_spin2and3[ 6*n_vect23+j+jj] = in[ 6*n_vect_in+j+jj];
                                out_spin2and3[ 7*n_vect23+j+jj] = in[ 7*n_vect_in+j+jj];
                                out_spin2and3[ 8*n_vect23+j+jj] = in[ 8*n_vect_in+j+jj];
                                out_spin2and3[ 9*n_vect23+j+jj] = in[ 9*n_vect_in+j+jj];
                                out_spin2and3[10*n_vect23+j+jj] = in[10*n_vect_in+j+jj];
                                out_spin2and3[11*n_vect23+j+jj] = in[11*n_vect_in+j+jj];)
  }

  static inline void twospin_n_X_PRECISION_new( const buffer_PRECISION out_spin0and1, const buffer_PRECISION out_spin2and3, const buffer_PRECISION in, const int n_vect, const int n_vect01, const int n_vect23, const int n_vect_in ) {
    int j, jj;
    VECTOR_LOOP(j, n_vect, jj,  out_spin0and1[ 0*n_vect01+j+jj] -= in[ 0*n_vect_in+j+jj];
                                out_spin0and1[ 1*n_vect01+j+jj] -= in[ 1*n_vect_in+j+jj];
                                out_spin0and1[ 2*n_vect01+j+jj] -= in[ 2*n_vect_in+j+jj];
                                out_spin0and1[ 3*n_vect01+j+jj] -= in[ 3*n_vect_in+j+jj];
                                out_spin0and1[ 4*n_vect01+j+jj] -= in[ 4*n_vect_in+j+jj];
                                out_spin0and1[ 5*n_vect01+j+jj] -= in[ 5*n_vect_in+j+jj];
                                out_spin0and1[ 6*n_vect01+j+jj] -= GAMMA_X_SPIN2_VAL * in[ 3*GAMMA_X_SPIN2_CO   *n_vect_in+j+jj];
                                out_spin0and1[ 7*n_vect01+j+jj] -= GAMMA_X_SPIN2_VAL * in[(3*GAMMA_X_SPIN2_CO+1)*n_vect_in+j+jj];
                                out_spin0and1[ 8*n_vect01+j+jj] -= GAMMA_X_SPIN2_VAL * in[(3*GAMMA_X_SPIN2_CO+2)*n_vect_in+j+jj];
                                out_spin0and1[ 9*n_vect01+j+jj] -= GAMMA_X_SPIN3_VAL * in[ 3*GAMMA_X_SPIN3_CO   *n_vect_in+j+jj];
                                out_spin0and1[10*n_vect01+j+jj] -= GAMMA_X_SPIN3_VAL * in[(3*GAMMA_X_SPIN3_CO+1)*n_vect_in+j+jj];
                                out_spin0and1[11*n_vect01+j+jj] -= GAMMA_X_SPIN3_VAL * in[(3*GAMMA_X_SPIN3_CO+2)*n_vect_in+j+jj];
                                out_spin2and3[ 0*n_vect23+j+jj] -= GAMMA_X_SPIN0_VAL * in[ 3*GAMMA_X_SPIN0_CO   *n_vect_in+j+jj];
                                out_spin2and3[ 1*n_vect23+j+jj] -= GAMMA_X_SPIN0_VAL * in[(3*GAMMA_X_SPIN0_CO+1)*n_vect_in+j+jj];
                                out_spin2and3[ 2*n_vect23+j+jj] -= GAMMA_X_SPIN0_VAL * in[(3*GAMMA_X_SPIN0_CO+2)*n_vect_in+j+jj];
                                out_spin2and3[ 3*n_vect23+j+jj] -= GAMMA_X_SPIN1_VAL * in[ 3*GAMMA_X_SPIN1_CO   *n_vect_in+j+jj];
                                out_spin2and3[ 4*n_vect23+j+jj] -= GAMMA_X_SPIN1_VAL * in[(3*GAMMA_X_SPIN1_CO+1)*n_vect_in+j+jj];
                                out_spin2and3[ 5*n_vect23+j+jj] -= GAMMA_X_SPIN1_VAL * in[(3*GAMMA_X_SPIN1_CO+2)*n_vect_in+j+jj];
                                out_spin2and3[ 6*n_vect23+j+jj] -= in[ 6*n_vect_in+j+jj];
                                out_spin2and3[ 7*n_vect23+j+jj] -= in[ 7*n_vect_in+j+jj];
                                out_spin2and3[ 8*n_vect23+j+jj] -= in[ 8*n_vect_in+j+jj];
                                out_spin2and3[ 9*n_vect23+j+jj] -= in[ 9*n_vect_in+j+jj];
                                out_spin2and3[10*n_vect23+j+jj] -= in[10*n_vect_in+j+jj];
                                out_spin2and3[11*n_vect23+j+jj] -= in[11*n_vect_in+j+jj];)
  }

// for TM1p1
  static inline void doublet_site_clover_PRECISION( const buffer_PRECISION eta, const buffer_PRECISION phi, const config_PRECISION clover ) {
    // diagonal
    eta[ 0] = clover[ 0]*phi[ 0];
    eta[ 1] = clover[ 1]*phi[ 1];
    eta[ 2] = clover[ 2]*phi[ 2];
    eta[ 3] = clover[ 3]*phi[ 3];
    eta[ 4] = clover[ 4]*phi[ 4];
    eta[ 5] = clover[ 5]*phi[ 5];
    eta[ 6] = clover[ 0]*phi[ 6];
    eta[ 7] = clover[ 1]*phi[ 7];
    eta[ 8] = clover[ 2]*phi[ 8];
    eta[ 9] = clover[ 3]*phi[ 9];
    eta[10] = clover[ 4]*phi[10];
    eta[11] = clover[ 5]*phi[11];
    eta[12] = clover[ 6]*phi[12];
    eta[13] = clover[ 7]*phi[13];
    eta[14] = clover[ 8]*phi[14];
    eta[15] = clover[ 9]*phi[15];
    eta[16] = clover[10]*phi[16];
    eta[17] = clover[11]*phi[17];
    eta[18] = clover[ 6]*phi[18];
    eta[19] = clover[ 7]*phi[19];
    eta[20] = clover[ 8]*phi[20];
    eta[21] = clover[ 9]*phi[21];
    eta[22] = clover[10]*phi[22];
    eta[23] = clover[11]*phi[23];
    // spin 0 and 1 flav 1
    eta[0] += clover[12]*phi[1];
    eta[0] += clover[13]*phi[2];
    eta[0] += clover[14]*phi[3];
    eta[0] += clover[15]*phi[4];
    eta[0] += clover[16]*phi[5];
    eta[1] += clover[17]*phi[2];
    eta[1] += clover[18]*phi[3];
    eta[1] += clover[19]*phi[4];
    eta[1] += clover[20]*phi[5];
    eta[2] += clover[21]*phi[3];
    eta[2] += clover[22]*phi[4];
    eta[2] += clover[23]*phi[5];
    eta[3] += clover[24]*phi[4];
    eta[3] += clover[25]*phi[5];
    eta[4] += clover[26]*phi[5];
    eta[1] += conj_PRECISION(clover[12])*phi[0];
    eta[2] += conj_PRECISION(clover[13])*phi[0];
    eta[3] += conj_PRECISION(clover[14])*phi[0];
    eta[4] += conj_PRECISION(clover[15])*phi[0];
    eta[5] += conj_PRECISION(clover[16])*phi[0];
    eta[2] += conj_PRECISION(clover[17])*phi[1];
    eta[3] += conj_PRECISION(clover[18])*phi[1];
    eta[4] += conj_PRECISION(clover[19])*phi[1];
    eta[5] += conj_PRECISION(clover[20])*phi[1];
    eta[3] += conj_PRECISION(clover[21])*phi[2];
    eta[4] += conj_PRECISION(clover[22])*phi[2];
    eta[5] += conj_PRECISION(clover[23])*phi[2];
    eta[4] += conj_PRECISION(clover[24])*phi[3];
    eta[5] += conj_PRECISION(clover[25])*phi[3];
    eta[5] += conj_PRECISION(clover[26])*phi[4];
    // spin 0 and 1 flav 2
    eta[ 6] += clover[12]*phi[ 7];
    eta[ 6] += clover[13]*phi[ 8];
    eta[ 6] += clover[14]*phi[ 9];
    eta[ 6] += clover[15]*phi[10];
    eta[ 6] += clover[16]*phi[11];
    eta[ 7] += clover[17]*phi[ 8];
    eta[ 7] += clover[18]*phi[ 9];
    eta[ 7] += clover[19]*phi[10];
    eta[ 7] += clover[20]*phi[11];
    eta[ 8] += clover[21]*phi[ 9];
    eta[ 8] += clover[22]*phi[10];
    eta[ 8] += clover[23]*phi[11];
    eta[ 9] += clover[24]*phi[10];
    eta[ 9] += clover[25]*phi[11];
    eta[10] += clover[26]*phi[11];
    eta[ 7] += conj_PRECISION(clover[12])*phi[ 6];
    eta[ 8] += conj_PRECISION(clover[13])*phi[ 6];
    eta[ 9] += conj_PRECISION(clover[14])*phi[ 6];
    eta[10] += conj_PRECISION(clover[15])*phi[ 6];
    eta[11] += conj_PRECISION(clover[16])*phi[ 6];
    eta[ 8] += conj_PRECISION(clover[17])*phi[ 7];
    eta[ 9] += conj_PRECISION(clover[18])*phi[ 7];
    eta[10] += conj_PRECISION(clover[19])*phi[ 7];
    eta[11] += conj_PRECISION(clover[20])*phi[ 7];
    eta[ 9] += conj_PRECISION(clover[21])*phi[ 8];
    eta[10] += conj_PRECISION(clover[22])*phi[ 8];
    eta[11] += conj_PRECISION(clover[23])*phi[ 8];
    eta[10] += conj_PRECISION(clover[24])*phi[ 9];
    eta[11] += conj_PRECISION(clover[25])*phi[ 9];
    eta[11] += conj_PRECISION(clover[26])*phi[10];
    // spin 2 and 3 flav 1
    eta[12] += clover[28]*phi[14];
    eta[12] += clover[27]*phi[13];
    eta[12] += clover[29]*phi[15];
    eta[12] += clover[30]*phi[16];
    eta[12] += clover[31]*phi[17];
    eta[13] += clover[32]*phi[14];
    eta[13] += clover[33]*phi[15];
    eta[13] += clover[34]*phi[16];
    eta[13] += clover[35]*phi[17];
    eta[14] += clover[36]*phi[15];
    eta[14] += clover[37]*phi[16];
    eta[14] += clover[38]*phi[17];
    eta[15] += clover[39]*phi[16];
    eta[15] += clover[40]*phi[17];
    eta[16] += clover[41]*phi[17];
    eta[13] += conj_PRECISION(clover[27])*phi[12];
    eta[14] += conj_PRECISION(clover[28])*phi[12];
    eta[15] += conj_PRECISION(clover[29])*phi[12];
    eta[16] += conj_PRECISION(clover[30])*phi[12];
    eta[17] += conj_PRECISION(clover[31])*phi[12];
    eta[14] += conj_PRECISION(clover[32])*phi[13];
    eta[15] += conj_PRECISION(clover[33])*phi[13];
    eta[16] += conj_PRECISION(clover[34])*phi[13];
    eta[17] += conj_PRECISION(clover[35])*phi[13];
    eta[15] += conj_PRECISION(clover[36])*phi[14];
    eta[16] += conj_PRECISION(clover[37])*phi[14];
    eta[17] += conj_PRECISION(clover[38])*phi[14];
    eta[16] += conj_PRECISION(clover[39])*phi[15];
    eta[17] += conj_PRECISION(clover[40])*phi[15];
    eta[17] += conj_PRECISION(clover[41])*phi[16];
    // spin 2 and 3 flav 2
    eta[18] += clover[28]*phi[20];
    eta[18] += clover[27]*phi[19];
    eta[18] += clover[29]*phi[21];
    eta[18] += clover[30]*phi[22];
    eta[18] += clover[31]*phi[23];
    eta[19] += clover[32]*phi[20];
    eta[19] += clover[33]*phi[21];
    eta[19] += clover[34]*phi[22];
    eta[19] += clover[35]*phi[23];
    eta[20] += clover[36]*phi[21];
    eta[20] += clover[37]*phi[22];
    eta[20] += clover[38]*phi[23];
    eta[21] += clover[39]*phi[22];
    eta[21] += clover[40]*phi[23];
    eta[22] += clover[41]*phi[23];
    eta[19] += conj_PRECISION(clover[27])*phi[18];
    eta[20] += conj_PRECISION(clover[28])*phi[18];
    eta[21] += conj_PRECISION(clover[29])*phi[18];
    eta[22] += conj_PRECISION(clover[30])*phi[18];
    eta[23] += conj_PRECISION(clover[31])*phi[18];
    eta[20] += conj_PRECISION(clover[32])*phi[19];
    eta[21] += conj_PRECISION(clover[33])*phi[19];
    eta[22] += conj_PRECISION(clover[34])*phi[19];
    eta[23] += conj_PRECISION(clover[35])*phi[19];
    eta[21] += conj_PRECISION(clover[36])*phi[20];
    eta[22] += conj_PRECISION(clover[37])*phi[20];
    eta[23] += conj_PRECISION(clover[38])*phi[20];
    eta[22] += conj_PRECISION(clover[39])*phi[21];
    eta[23] += conj_PRECISION(clover[40])*phi[21];
    eta[23] += conj_PRECISION(clover[41])*phi[22];
  }

  static inline void spin0and1_site_clover_PRECISION_new( const buffer_PRECISION eta, const buffer_PRECISION phi, const config_PRECISION clover, const int n_vect, const int n_vect_eta, const int n_vect_phi ) {
    int j, jj;
    // diagonal//????????
    VECTOR_LOOP( j, n_vect, jj, eta[ 0*n_vect_eta+j+jj] = clover[ 0]*phi[ 0*n_vect_phi+j+jj];
                                eta[ 1*n_vect_eta+j+jj] = clover[ 1]*phi[ 1*n_vect_phi+j+jj];
                                eta[ 2*n_vect_eta+j+jj] = clover[ 2]*phi[ 2*n_vect_phi+j+jj];
                                eta[ 3*n_vect_eta+j+jj] = clover[ 3]*phi[ 3*n_vect_phi+j+jj];
                                eta[ 4*n_vect_eta+j+jj] = clover[ 4]*phi[ 4*n_vect_phi+j+jj];
                                eta[ 5*n_vect_eta+j+jj] = clover[ 5]*phi[ 5*n_vect_phi+j+jj];
                                eta[ 6*n_vect_eta+j+jj] = _COMPLEX_PRECISION_ZERO;
                                eta[ 7*n_vect_eta+j+jj] = _COMPLEX_PRECISION_ZERO;
                                eta[ 8*n_vect_eta+j+jj] = _COMPLEX_PRECISION_ZERO;
                                eta[ 9*n_vect_eta+j+jj] = _COMPLEX_PRECISION_ZERO;
                                eta[10*n_vect_eta+j+jj] = _COMPLEX_PRECISION_ZERO;
                                eta[11*n_vect_eta+j+jj] = _COMPLEX_PRECISION_ZERO;
                                // spin 0 and 1
                                eta[0*n_vect_eta+j+jj] += clover[12]*phi[1*n_vect_phi+j+jj];
                                eta[0*n_vect_eta+j+jj] += clover[13]*phi[2*n_vect_phi+j+jj];
                                eta[0*n_vect_eta+j+jj] += clover[14]*phi[3*n_vect_phi+j+jj];
                                eta[0*n_vect_eta+j+jj] += clover[15]*phi[4*n_vect_phi+j+jj];
                                eta[0*n_vect_eta+j+jj] += clover[16]*phi[5*n_vect_phi+j+jj];
                                eta[1*n_vect_eta+j+jj] += clover[17]*phi[2*n_vect_phi+j+jj];
                                eta[1*n_vect_eta+j+jj] += clover[18]*phi[3*n_vect_phi+j+jj];
                                eta[1*n_vect_eta+j+jj] += clover[19]*phi[4*n_vect_phi+j+jj];
                                eta[1*n_vect_eta+j+jj] += clover[20]*phi[5*n_vect_phi+j+jj];
                                eta[2*n_vect_eta+j+jj] += clover[21]*phi[3*n_vect_phi+j+jj];
                                eta[2*n_vect_eta+j+jj] += clover[22]*phi[4*n_vect_phi+j+jj];
                                eta[2*n_vect_eta+j+jj] += clover[23]*phi[5*n_vect_phi+j+jj];
                                eta[3*n_vect_eta+j+jj] += clover[24]*phi[4*n_vect_phi+j+jj];
                                eta[3*n_vect_eta+j+jj] += clover[25]*phi[5*n_vect_phi+j+jj];
                                eta[4*n_vect_eta+j+jj] += clover[26]*phi[5*n_vect_phi+j+jj];
                                eta[1*n_vect_eta+j+jj] += conj_PRECISION(clover[12])*phi[0*n_vect_phi+j+jj];
                                eta[2*n_vect_eta+j+jj] += conj_PRECISION(clover[13])*phi[0*n_vect_phi+j+jj];
                                eta[3*n_vect_eta+j+jj] += conj_PRECISION(clover[14])*phi[0*n_vect_phi+j+jj];
                                eta[4*n_vect_eta+j+jj] += conj_PRECISION(clover[15])*phi[0*n_vect_phi+j+jj];
                                eta[5*n_vect_eta+j+jj] += conj_PRECISION(clover[16])*phi[0*n_vect_phi+j+jj];
                                eta[2*n_vect_eta+j+jj] += conj_PRECISION(clover[17])*phi[1*n_vect_phi+j+jj];
                                eta[3*n_vect_eta+j+jj] += conj_PRECISION(clover[18])*phi[1*n_vect_phi+j+jj];
                                eta[4*n_vect_eta+j+jj] += conj_PRECISION(clover[19])*phi[1*n_vect_phi+j+jj];
                                eta[5*n_vect_eta+j+jj] += conj_PRECISION(clover[20])*phi[1*n_vect_phi+j+jj];
                                eta[3*n_vect_eta+j+jj] += conj_PRECISION(clover[21])*phi[2*n_vect_phi+j+jj];
                                eta[4*n_vect_eta+j+jj] += conj_PRECISION(clover[22])*phi[2*n_vect_phi+j+jj];
                                eta[5*n_vect_eta+j+jj] += conj_PRECISION(clover[23])*phi[2*n_vect_phi+j+jj];
                                eta[4*n_vect_eta+j+jj] += conj_PRECISION(clover[24])*phi[3*n_vect_phi+j+jj];
                                eta[5*n_vect_eta+j+jj] += conj_PRECISION(clover[25])*phi[3*n_vect_phi+j+jj];
                                eta[5*n_vect_eta+j+jj] += conj_PRECISION(clover[26])*phi[4*n_vect_phi+j+jj];)
  }

  static inline void spin2and3_site_clover_PRECISION_new( const buffer_PRECISION eta, const buffer_PRECISION phi, const config_PRECISION clover, const int n_vect, const int n_vect_eta, const int n_vect_phi ) {
    int j, jj;
                                // diagonal
    VECTOR_LOOP( j, n_vect, jj, eta[ 0*n_vect_eta+j+jj] = _COMPLEX_PRECISION_ZERO;
                                eta[ 1*n_vect_eta+j+jj] = _COMPLEX_PRECISION_ZERO;
                                eta[ 2*n_vect_eta+j+jj] = _COMPLEX_PRECISION_ZERO;
                                eta[ 3*n_vect_eta+j+jj] = _COMPLEX_PRECISION_ZERO;
                                eta[ 4*n_vect_eta+j+jj] = _COMPLEX_PRECISION_ZERO;
                                eta[ 5*n_vect_eta+j+jj] = _COMPLEX_PRECISION_ZERO;
                                eta[ 6*n_vect_eta+j+jj] = clover[ 6]*phi[ 6*n_vect_phi+j+jj];
                                eta[ 7*n_vect_eta+j+jj] = clover[ 7]*phi[ 7*n_vect_phi+j+jj];
                                eta[ 8*n_vect_eta+j+jj] = clover[ 8]*phi[ 8*n_vect_phi+j+jj];
                                eta[ 9*n_vect_eta+j+jj] = clover[ 9]*phi[ 9*n_vect_phi+j+jj];
                                eta[10*n_vect_eta+j+jj] = clover[10]*phi[10*n_vect_phi+j+jj];
                                eta[11*n_vect_eta+j+jj] = clover[11]*phi[11*n_vect_phi+j+jj];
                                // spin 2 and 3
                                eta[ 6*n_vect_eta+j+jj] += clover[28]*phi[ 8*n_vect_phi+j+jj];
                                eta[ 6*n_vect_eta+j+jj] += clover[27]*phi[ 7*n_vect_phi+j+jj];
                                eta[ 6*n_vect_eta+j+jj] += clover[29]*phi[ 9*n_vect_phi+j+jj];
                                eta[ 6*n_vect_eta+j+jj] += clover[30]*phi[10*n_vect_phi+j+jj];
                                eta[ 6*n_vect_eta+j+jj] += clover[31]*phi[11*n_vect_phi+j+jj];
                                eta[ 7*n_vect_eta+j+jj] += clover[32]*phi[ 8*n_vect_phi+j+jj];
                                eta[ 7*n_vect_eta+j+jj] += clover[33]*phi[ 9*n_vect_phi+j+jj];
                                eta[ 7*n_vect_eta+j+jj] += clover[34]*phi[10*n_vect_phi+j+jj];
                                eta[ 7*n_vect_eta+j+jj] += clover[35]*phi[11*n_vect_phi+j+jj];
                                eta[ 8*n_vect_eta+j+jj] += clover[36]*phi[ 9*n_vect_phi+j+jj];
                                eta[ 8*n_vect_eta+j+jj] += clover[37]*phi[10*n_vect_phi+j+jj];
                                eta[ 8*n_vect_eta+j+jj] += clover[38]*phi[11*n_vect_phi+j+jj];
                                eta[ 9*n_vect_eta+j+jj] += clover[39]*phi[10*n_vect_phi+j+jj];
                                eta[ 9*n_vect_eta+j+jj] += clover[40]*phi[11*n_vect_phi+j+jj];
                                eta[10*n_vect_eta+j+jj] += clover[41]*phi[11*n_vect_phi+j+jj];
                                eta[ 7*n_vect_eta+j+jj] += conj_PRECISION(clover[27])*phi[ 6*n_vect_phi+j+jj];
                                eta[ 8*n_vect_eta+j+jj] += conj_PRECISION(clover[28])*phi[ 6*n_vect_phi+j+jj];
                                eta[ 9*n_vect_eta+j+jj] += conj_PRECISION(clover[29])*phi[ 6*n_vect_phi+j+jj];
                                eta[10*n_vect_eta+j+jj] += conj_PRECISION(clover[30])*phi[ 6*n_vect_phi+j+jj];
                                eta[11*n_vect_eta+j+jj] += conj_PRECISION(clover[31])*phi[ 6*n_vect_phi+j+jj];
                                eta[ 8*n_vect_eta+j+jj] += conj_PRECISION(clover[32])*phi[ 7*n_vect_phi+j+jj];
                                eta[ 9*n_vect_eta+j+jj] += conj_PRECISION(clover[33])*phi[ 7*n_vect_phi+j+jj];
                                eta[10*n_vect_eta+j+jj] += conj_PRECISION(clover[34])*phi[ 7*n_vect_phi+j+jj];
                                eta[11*n_vect_eta+j+jj] += conj_PRECISION(clover[35])*phi[ 7*n_vect_phi+j+jj];
                                eta[ 9*n_vect_eta+j+jj] += conj_PRECISION(clover[36])*phi[ 8*n_vect_phi+j+jj];
                                eta[10*n_vect_eta+j+jj] += conj_PRECISION(clover[37])*phi[ 8*n_vect_phi+j+jj];
                                eta[11*n_vect_eta+j+jj] += conj_PRECISION(clover[38])*phi[ 8*n_vect_phi+j+jj];
                                eta[10*n_vect_eta+j+jj] += conj_PRECISION(clover[39])*phi[ 9*n_vect_phi+j+jj];
                                eta[11*n_vect_eta+j+jj] += conj_PRECISION(clover[40])*phi[ 9*n_vect_phi+j+jj];
                                eta[11*n_vect_eta+j+jj] += conj_PRECISION(clover[41])*phi[10*n_vect_phi+j+jj];)
  }

  static inline void site_clover_PRECISION_new( const buffer_PRECISION eta, const buffer_PRECISION phi, const config_PRECISION clover, const int n_vect, const int n_vect_eta, const int n_vect_phi ) {
    int j, jj;
    VECTOR_LOOP(j, n_vect, jj, // diagonal
                               eta[ 0*n_vect_eta+j+jj] = clover[ 0]*phi[ 0*n_vect_phi+j+jj];
                               eta[ 1*n_vect_eta+j+jj] = clover[ 1]*phi[ 1*n_vect_phi+j+jj];
                               eta[ 2*n_vect_eta+j+jj] = clover[ 2]*phi[ 2*n_vect_phi+j+jj];
                               eta[ 3*n_vect_eta+j+jj] = clover[ 3]*phi[ 3*n_vect_phi+j+jj];
                               eta[ 4*n_vect_eta+j+jj] = clover[ 4]*phi[ 4*n_vect_phi+j+jj];
                               eta[ 5*n_vect_eta+j+jj] = clover[ 5]*phi[ 5*n_vect_phi+j+jj];
                               eta[ 6*n_vect_eta+j+jj] = clover[ 6]*phi[ 6*n_vect_phi+j+jj];
                               eta[ 7*n_vect_eta+j+jj] = clover[ 7]*phi[ 7*n_vect_phi+j+jj];
                               eta[ 8*n_vect_eta+j+jj] = clover[ 8]*phi[ 8*n_vect_phi+j+jj];
                               eta[ 9*n_vect_eta+j+jj] = clover[ 9]*phi[ 9*n_vect_phi+j+jj];
                               eta[10*n_vect_eta+j+jj] = clover[10]*phi[10*n_vect_phi+j+jj];
                               eta[11*n_vect_eta+j+jj] = clover[11]*phi[11*n_vect_phi+j+jj];
                               // spin 0 and 1, row major
                               eta[0*n_vect_eta+j+jj] += clover[12]*phi[1*n_vect_phi+j+jj];
                               eta[0*n_vect_eta+j+jj] += clover[13]*phi[2*n_vect_phi+j+jj];
                               eta[0*n_vect_eta+j+jj] += clover[14]*phi[3*n_vect_phi+j+jj];
                               eta[0*n_vect_eta+j+jj] += clover[15]*phi[4*n_vect_phi+j+jj];
                               eta[0*n_vect_eta+j+jj] += clover[16]*phi[5*n_vect_phi+j+jj];
                               eta[1*n_vect_eta+j+jj] += clover[17]*phi[2*n_vect_phi+j+jj];
                               eta[1*n_vect_eta+j+jj] += clover[18]*phi[3*n_vect_phi+j+jj];
                               eta[1*n_vect_eta+j+jj] += clover[19]*phi[4*n_vect_phi+j+jj];
                               eta[1*n_vect_eta+j+jj] += clover[20]*phi[5*n_vect_phi+j+jj];
                               eta[2*n_vect_eta+j+jj] += clover[21]*phi[3*n_vect_phi+j+jj];
                               eta[2*n_vect_eta+j+jj] += clover[22]*phi[4*n_vect_phi+j+jj];
                               eta[2*n_vect_eta+j+jj] += clover[23]*phi[5*n_vect_phi+j+jj];
                               eta[3*n_vect_eta+j+jj] += clover[24]*phi[4*n_vect_phi+j+jj];
                               eta[3*n_vect_eta+j+jj] += clover[25]*phi[5*n_vect_phi+j+jj];
                               eta[4*n_vect_eta+j+jj] += clover[26]*phi[5*n_vect_phi+j+jj];
                               eta[1*n_vect_eta+j+jj] += conj_PRECISION(clover[12])*phi[0*n_vect_phi+j+jj];
                               eta[2*n_vect_eta+j+jj] += conj_PRECISION(clover[13])*phi[0*n_vect_phi+j+jj];
                               eta[3*n_vect_eta+j+jj] += conj_PRECISION(clover[14])*phi[0*n_vect_phi+j+jj];
                               eta[4*n_vect_eta+j+jj] += conj_PRECISION(clover[15])*phi[0*n_vect_phi+j+jj];
                               eta[5*n_vect_eta+j+jj] += conj_PRECISION(clover[16])*phi[0*n_vect_phi+j+jj];
                               eta[2*n_vect_eta+j+jj] += conj_PRECISION(clover[17])*phi[1*n_vect_phi+j+jj];
                               eta[3*n_vect_eta+j+jj] += conj_PRECISION(clover[18])*phi[1*n_vect_phi+j+jj];
                               eta[4*n_vect_eta+j+jj] += conj_PRECISION(clover[19])*phi[1*n_vect_phi+j+jj];
                               eta[5*n_vect_eta+j+jj] += conj_PRECISION(clover[20])*phi[1*n_vect_phi+j+jj];
                               eta[3*n_vect_eta+j+jj] += conj_PRECISION(clover[21])*phi[2*n_vect_phi+j+jj];
                               eta[4*n_vect_eta+j+jj] += conj_PRECISION(clover[22])*phi[2*n_vect_phi+j+jj];
                               eta[5*n_vect_eta+j+jj] += conj_PRECISION(clover[23])*phi[2*n_vect_phi+j+jj];
                               eta[4*n_vect_eta+j+jj] += conj_PRECISION(clover[24])*phi[3*n_vect_phi+j+jj];
                               eta[5*n_vect_eta+j+jj] += conj_PRECISION(clover[25])*phi[3*n_vect_phi+j+jj];
                               eta[5*n_vect_eta+j+jj] += conj_PRECISION(clover[26])*phi[4*n_vect_phi+j+jj];
                               // spin 2 and 3, row major
                               eta[ 6*n_vect_eta+j+jj] += clover[27]*phi[ 7*n_vect_phi+j+jj];
                               eta[ 6*n_vect_eta+j+jj] += clover[28]*phi[ 8*n_vect_phi+j+jj];
                               eta[ 6*n_vect_eta+j+jj] += clover[29]*phi[ 9*n_vect_phi+j+jj];
                               eta[ 6*n_vect_eta+j+jj] += clover[30]*phi[10*n_vect_phi+j+jj];
                               eta[ 6*n_vect_eta+j+jj] += clover[31]*phi[11*n_vect_phi+j+jj];
                               eta[ 7*n_vect_eta+j+jj] += clover[32]*phi[ 8*n_vect_phi+j+jj];
                               eta[ 7*n_vect_eta+j+jj] += clover[33]*phi[ 9*n_vect_phi+j+jj];
                               eta[ 7*n_vect_eta+j+jj] += clover[34]*phi[10*n_vect_phi+j+jj];
                               eta[ 7*n_vect_eta+j+jj] += clover[35]*phi[11*n_vect_phi+j+jj];
                               eta[ 8*n_vect_eta+j+jj] += clover[36]*phi[ 9*n_vect_phi+j+jj];
                               eta[ 8*n_vect_eta+j+jj] += clover[37]*phi[10*n_vect_phi+j+jj];
                               eta[ 8*n_vect_eta+j+jj] += clover[38]*phi[11*n_vect_phi+j+jj];
                               eta[ 9*n_vect_eta+j+jj] += clover[39]*phi[10*n_vect_phi+j+jj];
                               eta[ 9*n_vect_eta+j+jj] += clover[40]*phi[11*n_vect_phi+j+jj];
                               eta[10*n_vect_eta+j+jj] += clover[41]*phi[11*n_vect_phi+j+jj];
                               eta[ 7*n_vect_eta+j+jj] += conj_PRECISION(clover[27])*phi[ 6*n_vect_phi+j+jj];
                               eta[ 8*n_vect_eta+j+jj] += conj_PRECISION(clover[28])*phi[ 6*n_vect_phi+j+jj];
                               eta[ 9*n_vect_eta+j+jj] += conj_PRECISION(clover[29])*phi[ 6*n_vect_phi+j+jj];
                               eta[10*n_vect_eta+j+jj] += conj_PRECISION(clover[30])*phi[ 6*n_vect_phi+j+jj];
                               eta[11*n_vect_eta+j+jj] += conj_PRECISION(clover[31])*phi[ 6*n_vect_phi+j+jj];
                               eta[ 8*n_vect_eta+j+jj] += conj_PRECISION(clover[32])*phi[ 7*n_vect_phi+j+jj];
                               eta[ 9*n_vect_eta+j+jj] += conj_PRECISION(clover[33])*phi[ 7*n_vect_phi+j+jj];
                               eta[10*n_vect_eta+j+jj] += conj_PRECISION(clover[34])*phi[ 7*n_vect_phi+j+jj];
                               eta[11*n_vect_eta+j+jj] += conj_PRECISION(clover[35])*phi[ 7*n_vect_phi+j+jj];
                               eta[ 9*n_vect_eta+j+jj] += conj_PRECISION(clover[36])*phi[ 8*n_vect_phi+j+jj];
                               eta[10*n_vect_eta+j+jj] += conj_PRECISION(clover[37])*phi[ 8*n_vect_phi+j+jj];
                               eta[11*n_vect_eta+j+jj] += conj_PRECISION(clover[38])*phi[ 8*n_vect_phi+j+jj];
                               eta[10*n_vect_eta+j+jj] += conj_PRECISION(clover[39])*phi[ 9*n_vect_phi+j+jj];
                               eta[11*n_vect_eta+j+jj] += conj_PRECISION(clover[40])*phi[ 9*n_vect_phi+j+jj];
                               eta[11*n_vect_eta+j+jj] += conj_PRECISION(clover[41])*phi[10*n_vect_phi+j+jj];)
  }
 
#endif 
