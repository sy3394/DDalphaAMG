/*
 * copyright (C) 2016, Matthias Rottmann, Artur Strebel, Simon Heybrock, Simone Bacchio, Bjoern Leder.
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
 *  copied: 11/27/2019 (US system)
 * changed from sbacchio
 * checked: 12/07/2019
 * 1st cleanup:12/18/2019
 */

#include "main.h"

static void set_coarse_self_coupling_PRECISION( vector_PRECISION *buffer1, vector_PRECISION *buffer2, vector_PRECISION *V, const int n, level_struct *l );
static void set_coarse_neighbor_coupling_PRECISION( vector_PRECISION *buffer1, vector_PRECISION *buffer2, vector_PRECISION *V, const int mu, const int n, level_struct *l );
static void set_block_diagonal_PRECISION( vector_PRECISION *spin_0_1, vector_PRECISION *spin_2_3, vector_PRECISION *V, const int n, config_PRECISION block, level_struct *l );
static void coarse_spinwise_self_couplings_PRECISION( vector_PRECISION *eta1, vector_PRECISION *eta2, vector_PRECISION *phi, config_PRECISION clover, 

							  int length, level_struct *l );
void coarse_operator_PRECISION_alloc( level_struct *l ) {

  int nd = l->next_level->num_inner_lattice_sites;
  int  k = l->next_level->num_parent_eig_vect*2;   // internal d.o.f. on the next level

  l->next_level->D_size      = k*k*4*nd;           // (internal d.o.f.)**2*(# nearest neighbor couplings)*(# lattice sites)
  l->next_level->clover_size = ((k*(k+1))/2)*nd;   // upper triangle???
  l->next_level->block_size  = ((k/2*(k/2+1)))*nd; // upper triangle???

  // allocates space for setting up an operator and compute neighbor tables.  
  operator_PRECISION_alloc( &(l->next_level->op_PRECISION), _ORDINARY, l->next_level );

}

void coarse_operator_PRECISION_free( level_struct *l ) {
  
  operator_PRECISION_free( &(l->next_level->op_PRECISION), _ORDINARY, l->next_level );
  
}

void coarse_operator_PRECISION_setup( vector_PRECISION *V, level_struct *l ) {
  /**************************************
   * Input
   *  vector_PRECISION *V: a orthonomal set of test vectors consisting of interpolation op.
   * Output
   *  Set l->next_level->op_PRECISION.D
   * note: see Matius thesis p. 55
   *************************************/
  double t0, t1, tc;
  t0 = MPI_Wtime();
  vector_PRECISION buffer[3];

  if ( V->num_vect_now != l->num_eig_vect || V->num_vect != l->num_eig_vect )
    error0("oarse_operator_PRECISION_setup: assumptions are not met\n");

  int mu, i, j;
  int n           = l->num_eig_vect;
  int D_size      = l->next_level->D_size;
  int clover_size = l->next_level->clover_size;
  int block_size  = l->next_level->block_size;
  void (*aggregate_self_coupling)()     = (l->depth==0)?d_plus_clover_aggregate_PRECISION:coarse_aggregate_self_couplings_PRECISION;
  void (*aggregate_neighbor_coupling)() = (l->depth==0)?d_neighbor_aggregate_PRECISION   :coarse_aggregate_neighbor_couplings_PRECISION;
  void (*aggregate_block)()             = (l->depth==0)?diagonal_aggregate_PRECISION     :coarse_aggregate_block_diagonal_PRECISION;

  for ( i=0; i<3; i++ ) {
    vector_PRECISION_init(&(buffer[i]));
    vector_PRECISION_alloc( &(buffer[i]), _ORDINARY, num_loop, l, no_threading );
    //buffer[i] = l->vbuf_PRECISION[i+4];
    buffer[i].num_vect_now = num_loop; 
  }

  // define index tables for op application
  operator_PRECISION_define( &(l->next_level->op_PRECISION), l->next_level );
  // initialize op's  
  for ( j=0; j<D_size; j++ )
    l->next_level->op_PRECISION.D[j]        = _COMPLEX_PRECISION_ZERO;
  for ( j=0; j<clover_size; j++ )
    l->next_level->op_PRECISION.clover[j]   = _COMPLEX_PRECISION_ZERO;
  for ( j=0; j<block_size; j++ )
    l->next_level->op_PRECISION.odd_proj[j] = _COMPLEX_PRECISION_ZERO;

  if ( V->layout != _NVEC_INNER ) 
    error0("coarse_operator_PRECISION_setup: V->layout is not set correctly\n");

  // for all test vectors V[i]:
  for ( i=0; i<n; i+=num_loop ) {
    tc-=MPI_Wtime();
    vector_PRECISION_copy2( &buffer[2], V, i, num_loop, 1, 0, l->vector_size, l );
    tc+=MPI_Wtime();
    for ( mu=0; mu<4; mu++ ) {
      // update ghost cells of V[i]
      negative_sendrecv_PRECISION( &buffer[2], mu, &(l->s_PRECISION.op.c), l );
    }
    // calculate selfcoupling entries of the coarse grid operator
    //       apply self coupling of block-and-2spin-restricted dirac operator for each aggregate
    aggregate_self_coupling( &buffer[0], &buffer[1], &buffer[2], &(l->s_PRECISION), l );
    //       set l->next_level->op_PRECISION.clover)
    set_coarse_self_coupling_PRECISION( &buffer[0], &buffer[1], V, i, l );
    //       set l->s_PRECISION.op.odd_proj
    aggregate_block( &buffer[0], &buffer[1], &buffer[2], l->s_PRECISION.op.odd_proj, l );
    //       set l->next_level->op_PRECISION.odd_proj
    set_block_diagonal_PRECISION( &buffer[0], &buffer[1], V, i, l->next_level->op_PRECISION.odd_proj, l ); 
    for ( mu=0; mu<4; mu++ ) {
      //      finish updating ghostcells of V[i]
      negative_wait_PRECISION( mu, &(l->s_PRECISION.op.c), l );      
      //      apply 2spin-restricted dirac operator for direction mu for all aggregates
      aggregate_neighbor_coupling( &buffer[0], &buffer[1], &buffer[2], mu, &(l->s_PRECISION), l );      
      //     finally set l->next_level->op_PRECISION.D 
      tc-=MPI_Wtime();
      vector_PRECISION_copy2( V, &buffer[2], i, num_loop, -1, 0, l->vector_size, l );
      tc+=MPI_Wtime();
      set_coarse_neighbor_coupling_PRECISION( &buffer[0], &buffer[1], V, mu, i, l );
    }
  }
  for ( i=0; i<3; i++ ) vector_PRECISION_free( &(buffer[i]), l, no_threading);
  // set tm_term (and eps_term if HAVE_TM1p1)
  coarse_operator_PRECISION_setup_finalize( l, no_threading );
  t1 = MPI_Wtime();
  if ( g.print > 0 ) printf0("depth: %d, time spent for setting up next coarser operator: %lf seconds (%lf sec for copying) \n", l->depth, t1-t0, tc );

}

void coarse_operator_PRECISION_setup_finalize( level_struct *l, struct Thread *threading ) {

  int block_size = l->next_level->block_size;
  
  l->next_level->op_PRECISION.m0 = l->s_PRECISION.op.m0;
#ifdef HAVE_TM    
  //tm_term
  PRECISION mf = (g.mu_factor[l->depth]) ? g.mu_factor[l->next_level->depth]/g.mu_factor[l->depth]:0;
  if ( mf*l->s_PRECISION.op.mu + mf*l->s_PRECISION.op.mu_even_shift == 0 &&
       mf*l->s_PRECISION.op.mu + mf*l->s_PRECISION.op.mu_odd_shift == 0 )
    buffer_PRECISION_define( l->next_level->op_PRECISION.tm_term, _COMPLEX_double_ZERO, 0, block_size, l->next_level ); //block_size here is the size of the vector; cf coarse_operator_PRECISION_alloc
  else
    tm_term_PRECISION_setup( mf*l->s_PRECISION.op.mu, mf*l->s_PRECISION.op.mu_even_shift,
                             mf*l->s_PRECISION.op.mu_odd_shift, &(l->next_level->op_PRECISION),
                             l->next_level, threading ); 
#endif
#ifdef HAVE_TM1p1
  //eps_term
  PRECISION ef = (g.epsbar_factor[l->depth]) ? g.epsbar_factor[l->next_level->depth]/g.epsbar_factor[l->depth]:0; 
  if ( ef*l->s_PRECISION.op.epsbar == 0 &&  ef*l->s_PRECISION.op.epsbar_ig5_even_shift == 0 &&
       ef*l->s_PRECISION.op.epsbar_ig5_odd_shift == 0 )
    buffer_PRECISION_define( l->next_level->op_PRECISION.epsbar_term, _COMPLEX_double_ZERO, 0, block_size, l->next_level );
  else
    epsbar_term_PRECISION_setup( ef*l->s_PRECISION.op.epsbar, ef*l->s_PRECISION.op.epsbar_ig5_even_shift,
                                 ef*l->s_PRECISION.op.epsbar_ig5_odd_shift, &(l->next_level->op_PRECISION),
                                 l->next_level, threading );
#endif

}

// vector layout has changed from vector fastest to vector slowest
static void set_coarse_self_coupling_PRECISION( vector_PRECISION *spin_0_1, vector_PRECISION *spin_2_3, 
					     vector_PRECISION *V, const int n0, level_struct *l ) {
  
  int i, j, k, m, k1, k2, jj, jjj, n, p;
  int num_aggregates   = l->is_PRECISION.num_agg;
  int num_eig_vect     = l->next_level->num_parent_eig_vect;
  int aggregate_size   = l->num_inner_lattice_sites*l->num_parent_eig_vect*2/num_aggregates;
  int offset           = l->num_parent_eig_vect; // half of internal d.o.f
  int clover_site_size = (num_eig_vect*(2*num_eig_vect+1));
  int nvec_V = V->num_vect, nvec_01 = spin_0_1->num_vect, nvec_23 = spin_2_3->num_vect;
  config_PRECISION clover_pt, clover = l->next_level->op_PRECISION.clover;  
  
  // U(x) = [ A B      , A=A*, D=D*, C = -B*
  //          C D ]
  // storage order: upper triangle of A, upper triangle of D, B, columnwise
  // diagonal coupling
  for ( p=0; p<num_loop; p++ ) {
    n  = n0+p; 
    k1 = (n*(n+1))/2; k2 = (n*(n+1))/2+(num_eig_vect*(num_eig_vect+1))/2;
  
    for ( j=0; j<num_aggregates; j++ ) { // for each aggregate
      clover_pt   = clover + j*clover_site_size;
      
      for ( i=0; i<aggregate_size; ) { // for each site in aggregate
	// A
	for ( m=0; m<offset; m++, i++ ) {
	  VECTOR_LOOP(jj, n0, jjj, clover_pt[ k1+jj+jjj ] += 
	  	      conj_PRECISION( V->vector_buffer[(i+j*aggregate_size)*nvec_V+jj+jjj] ) * spin_0_1->vector_buffer[(i+j*aggregate_size)*nvec_01+p]; )
	  for ( k=n0; k<=n; k++ ) clover_pt[ k1+k ] += conj_PRECISION( V->vector_buffer[(i+j*aggregate_size)*nvec_V+k] ) * spin_0_1->vector_buffer[(i+j*aggregate_size)*nvec_01+p];
	}
	// D
	for ( m=0; m<offset; m++, i++ ) {
	  VECTOR_LOOP(jj, n0, jjj,clover_pt[ k2+jj+jjj ] += 
	  	      conj_PRECISION( V->vector_buffer[(i+j*aggregate_size)*nvec_V+jj+jjj] ) * spin_2_3->vector_buffer[(i+j*aggregate_size)*nvec_23+p]; )
	  for ( k=n0; k<=n; k++ ) clover_pt[ k2+k ] += conj_PRECISION( V->vector_buffer[(i+j*aggregate_size)*nvec_V+k] ) * spin_2_3->vector_buffer[(i+j*aggregate_size)*nvec_23+p]; 
	}
      }
    }

    for ( j=0; j<num_aggregates; j++ ) {
      clover_pt = clover + j*clover_site_size + num_eig_vect*(num_eig_vect+1+n);
      
      for ( i=0; i<aggregate_size; ) {
	// B
	for ( m=0; m<offset; m++, i++ )
	  VECTOR_LOOP(jj, num_eig_vect, jjj, clover_pt[ jj+jjj ] += 
		      conj_PRECISION( V->vector_buffer[(i+j*aggregate_size)*nvec_V+jj+jjj] ) * spin_2_3->vector_buffer[(i+j*aggregate_size)*nvec_23+p]; )
	i += offset;
      }
    }
  }
}

// n = num_eig_vect
static void set_block_diagonal_PRECISION( vector_PRECISION *spin_0_1, vector_PRECISION *spin_2_3, 
				       vector_PRECISION *V, const int n0, config_PRECISION block, level_struct *l ) {
  
  // U(x) = [ A 0      , A=A*, D=D*
  //          0 D ]
  // storage order: upper triangle of A, upper triangle of D, columnwise
  // suitable for tm_term and odd_proj
  int i, j, k, m, k1, k2, n, p, jj, jjj;
  int num_aggregates  = l->is_PRECISION.num_agg;
  int num_eig_vect    = l->next_level->num_parent_eig_vect;
  int aggregate_size  = l->num_inner_lattice_sites*l->num_parent_eig_vect*2/num_aggregates;
  int offset          = l->num_parent_eig_vect;
  int block_site_size = (num_eig_vect*(num_eig_vect+1));
  int nvec_V = V->num_vect, nvec_01 = spin_0_1->num_vect, nvec_23 = spin_2_3->num_vect;
  config_PRECISION block_pt;

  for ( p=0; p<num_loop; p++ ) {
    n = n0+p;
    k1 = (n*(n+1))/2; k2 = (n*(n+1))/2+block_site_size/2;
  
    for ( j=0; j<num_aggregates; j++ ) {
      block_pt = block + j*block_site_size;
      
      for ( i=0; i<aggregate_size; ) {
	// A
	for ( m=0; m<offset; m++, i++ ) {//k=0 <- k=n0
	  VECTOR_LOOP(jj, n0, jjj, block_pt[ k1+jj+jjj ] += 
	  	      conj_PRECISION( V->vector_buffer[(i+j*aggregate_size)*nvec_V+jj+jjj] ) * spin_0_1->vector_buffer[(i+j*aggregate_size)*nvec_01+p];)
	  for ( k=n0; k<=n; k++ ) block_pt[ k1+k ] += conj_PRECISION( V->vector_buffer[(i+j*aggregate_size)*nvec_V+k] ) * spin_0_1->vector_buffer[(i+j*aggregate_size)*nvec_01+p];
	}
	// D
	for ( m=0; m<offset; m++, i++ ) {
	  VECTOR_LOOP(jj, n0, jjj, block_pt[ k2+jj+jjj ] += 
	  	      conj_PRECISION( V->vector_buffer[(i+j*aggregate_size)*nvec_V+jj+jjj] ) * spin_2_3->vector_buffer[(i+j*aggregate_size)*nvec_23+p];)
	  for ( k=n0; k<=n; k++ ) block_pt[ k2+k ] += conj_PRECISION( V->vector_buffer[(i+j*aggregate_size)*nvec_V+k] ) * spin_2_3->vector_buffer[(i+j*aggregate_size)*nvec_23+p];
	}
      }
    }
  }
}

static void set_coarse_neighbor_coupling_PRECISION( vector_PRECISION *spin_0_1, vector_PRECISION *spin_2_3, 
                                             vector_PRECISION *V, const int mu, const int n0, level_struct *l ) {
  
  int i, i1, j, k, k1, k2, m, jj, jjj, n, p;
  int num_aggregates = l->is_PRECISION.num_agg;
  int num_eig_vect   = l->next_level->num_parent_eig_vect;
  int offset         = l->num_parent_eig_vect;
  int nlsv           = l->num_parent_eig_vect*2;
  int D_link_size              = num_eig_vect*num_eig_vect*4;
  int *index_dir               = l->is_PRECISION.agg_boundary_index[mu];
  int aggregate_boundary_sites = l->is_PRECISION.agg_boundary_length[mu]/num_aggregates;
  int aggregate_size           = l->num_inner_lattice_sites*l->num_parent_eig_vect*2/num_aggregates;
  int nvec_V = V->num_vect, nvec_01 = spin_0_1->num_vect, nvec_23 = spin_2_3->num_vect;
  config_PRECISION D_pt, D = l->next_level->op_PRECISION.D;
  
  if ( nvec_V != num_eig_vect )
    error0("set_coarse_neighbor_coupling_PRECISION: assumptions are not met\n");

  // U_mu(x) = [ A B      , U_-mu(x+muhat) = [ A* -C*
  //             C D ]                        -B*  D* ]
  // storage order: A, C, B, D, each column wise
  for ( p=0; p<num_loop; p++ ) {
    n  = n0+p;
    k1 = n*num_eig_vect;
    k2 = (n+num_eig_vect)*num_eig_vect;
    i1 = 0;
    for ( j=0; j<num_aggregates; j++ ) {
      D_pt = D+(j*4+mu)*D_link_size;
      for ( i=0; i<aggregate_boundary_sites; i++ ) {
	// A
	for ( m=0; m<offset; m++ )
	  VECTOR_LOOP(jj, num_eig_vect, jjj, D_pt[ k1+jj+jjj ] += 
		      conj_PRECISION( V->vector_buffer[(m+nlsv*index_dir[i1])*nvec_V+jj+jjj] ) * spin_0_1->vector_buffer[(m+nlsv*index_dir[i1])*nvec_01+p]; )
	// C
	for ( ; m<2*offset; m++ )
	  VECTOR_LOOP(jj, num_eig_vect, jjj, D_pt[ k2+jj+jjj ] += 
		      conj_PRECISION( V->vector_buffer[(m+nlsv*index_dir[i1])*nvec_V+jj+jjj] ) * spin_0_1->vector_buffer[(m+nlsv*index_dir[i1])*nvec_01+p]; )
	i1++;
      }
    }
    
    k1 = (n+2*num_eig_vect)*num_eig_vect;
    k2 = (n+3*num_eig_vect)*num_eig_vect;
    i1 = 0;
    for ( j=0; j<num_aggregates; j++ ) {
      D_pt = D+(j*4+mu)*D_link_size;
      for ( i=0; i<aggregate_boundary_sites; i++ ) {
	// B
	for ( m=0; m<offset; m++ )
	  VECTOR_LOOP(jj, num_eig_vect, jjj, D_pt[ k1+jj+jjj ] += 
		      conj_PRECISION( V->vector_buffer[(m+nlsv*index_dir[i1])*nvec_V+jj+jjj] ) * spin_2_3->vector_buffer[(m+nlsv*index_dir[i1])*nvec_23+p]; )
	// D
        for ( ; m<2*offset; m++ )
	  VECTOR_LOOP(jj, num_eig_vect, jjj, D_pt[ k2+jj+jjj ] += 
		      conj_PRECISION( V->vector_buffer[(m+nlsv*index_dir[i1])*nvec_V+jj+jjj] ) * spin_2_3->vector_buffer[(m+nlsv*index_dir[i1])*nvec_23+p]; )
	i1++;
      }
    }
  }
}

// eta <- coarse_self_couplings*phi
void coarse_self_couplings_PRECISION( vector_PRECISION *eta, vector_PRECISION *phi,
						 operator_PRECISION_struct *op, int start, int end, level_struct *l ) {

  int num_eig_vect = l->num_parent_eig_vect, nvec = phi->num_vect_now, nvec_phi = phi->num_vect, nvec_eta = eta->num_vect;
  int site_size    = l->num_lattice_site_var;
  int clover_size  = (2*num_eig_vect*num_eig_vect+num_eig_vect);
  int block_size   = (num_eig_vect*num_eig_vect+num_eig_vect);
  vector_PRECISION eta_pt, phi_pt;
  vector_PRECISION_duplicate( &eta_pt, eta, start, l );
  vector_PRECISION_duplicate( &phi_pt, phi, start, l );
  
  if ( nvec > nvec_eta )
    error0("coarse_self_couplings_PRECISION: assumptions are not met\n");

  coarse_self_couplings_clover_PRECISION( &eta_pt, &phi_pt, op->clover+start*clover_size, (end-start)*site_size, l );
#ifdef HAVE_TM // tm_term
  if (op->mu + op->mu_odd_shift != 0.0 || op->mu + op->mu_even_shift != 0.0 )
    coarse_add_anti_block_diagonal_PRECISION( &eta_pt, &phi_pt, op->tm_term+start*block_size, (end-start)*site_size, l );
#endif
/*#ifdef HAVE_TM1p1 //eps_term
  if ( g.n_flavours == 2 &&
       ( op->epsbar != 0 || op->epsbar_ig5_odd_shift != 0 || op->epsbar_ig5_odd_shift != 0 ) )
    coarse_add_doublet_coupling_PRECISION( eta+start*vector_size, phi+start*vector_size, 
                                           op->epsbar_term+start*block_size, (end-start)*vector_size, l );
#endif*/

}

// eta <- coarse_block_operator*phi
// start should be equal to site_index*l->num_lattice_site_var!!!!
void coarse_block_operator_PRECISION( vector_PRECISION *eta, vector_PRECISION *phi, int start, schwarz_PRECISION_struct *s, level_struct *l, struct Thread *threading ) {
  
  START_UNTHREADED_FUNCTION(threading)

  int n = s->num_block_sites, *length = s->dir_length, **index = s->index;
  int nvec = phi->num_vect_now, nvec_phi = phi->num_vect, nvec_eta = eta->num_vect;
  int *ind, *neighbor = s->op.neighbor_table, m = l->num_lattice_site_var, num_eig_vect = l->num_parent_eig_vect;
  vector_PRECISION lphi, leta;
  vector_PRECISION_duplicate( &lphi, phi, 0, l );//start/m, l );
  vector_PRECISION_duplicate( &leta, eta, 0, l );//start/m, l );//which is better?????

  if ( nvec_eta < nvec_phi )
    error0("coarse_block_operator_PRECISION: assumptions are not met\n");
  
  // site-wise self coupling
  coarse_self_couplings_PRECISION( &leta, &lphi, &(s->op), (start/m), (start/m)+n, l);

  // inner block couplings
  int hopp_size = 4 * SQUARE( num_eig_vect*2 );
  config_PRECISION D_pt, D = s->op.D + (start/m)*hopp_size;

  for ( int mu=0; mu<4; mu++ ) {
    ind = index[mu]; // mu direction
    for ( int i=0; i<length[mu]; i++ ) {
      int k = ind[i]; int j = neighbor[5*k+mu+1];
      D_pt = D + hopp_size*k + (hopp_size/4)*mu;
      leta.vector_buffer = eta->vector_buffer+(start+m*k)*nvec_eta;
      lphi.vector_buffer = phi->vector_buffer+(start+m*j)*nvec_phi;
      coarse_hopp_PRECISION( &leta, &lphi, D_pt, l );

      leta.vector_buffer = eta->vector_buffer+(start+m*j)*nvec_eta;
      lphi.vector_buffer = phi->vector_buffer+(start+m*k)*nvec_phi;
      coarse_daggered_hopp_PRECISION( &leta, &lphi, D_pt, l );
    }
  }
  END_UNTHREADED_FUNCTION(threading)
}

// eta1,eta2 <- coarse_aggregate_self_couplings_PRECISION_*phi (eta1&2 are initialized to 0's)
void coarse_aggregate_self_couplings_PRECISION( vector_PRECISION *eta1, vector_PRECISION *eta2, vector_PRECISION *phi,
						    schwarz_PRECISION_struct *s, level_struct *l ) {
  
  int i, mu, index1, index2, length, *index_dir;
  int *neighbor = s->op.neighbor_table;
  int n         = l->num_lattice_site_var;
  int Dls       = n*n;
  int Dss       = 4*n*n;
  int nvec = phi->num_vect_now, nvec_phi = phi->num_vect, nvec_eta1 = eta1->num_vect, nvec_eta2 = eta2->num_vect;

  vector_PRECISION eta1_pt, eta2_pt, phi_pt;
  config_PRECISION D_pt, D = s->op.D;
  vector_PRECISION_duplicate( &eta1_pt, eta1, 0, l );
  vector_PRECISION_duplicate( &eta2_pt, eta2, 0, l );
  vector_PRECISION_duplicate( &phi_pt, phi, 0, l );
  
  if ( nvec_eta1 != nvec_eta2 || nvec_eta1 < nvec_phi )
    error0("coarse_aggregate_self_couplings_PRECISION: assumptions are not met\n");

  vector_PRECISION_define( eta1, 0, 0, l->vector_size, l );
  vector_PRECISION_define( eta2, 0, 0, l->vector_size, l );  
  coarse_spinwise_self_couplings_PRECISION( eta1, eta2, phi, s->op.clover, l->inner_vector_size, l );
  
  for ( mu=0; mu<4; mu++ ) { // direction mu
    length = l->is_PRECISION.agg_length[mu]; index_dir = l->is_PRECISION.agg_index[mu];
    for ( i=0; i<length; i++ ) {
      index1 = index_dir[i]; index2 = neighbor[5*index1+mu+1]; D_pt = D + Dss*index1 + Dls*mu;
      phi_pt.vector_buffer  = phi->vector_buffer  + n*index2*nvec_phi; 
      eta1_pt.vector_buffer = eta1->vector_buffer + n*index1*nvec_eta1; 
      eta2_pt.vector_buffer = eta2->vector_buffer + n*index1*nvec_eta2;
      coarse_spinwise_n_hopp_PRECISION( &eta1_pt, &eta2_pt, &phi_pt, D_pt, l );
      phi_pt.vector_buffer  = phi->vector_buffer  + n*index1*nvec_phi; 
      eta1_pt.vector_buffer = eta1->vector_buffer + n*index2*nvec_eta1; 
      eta2_pt.vector_buffer = eta2->vector_buffer + n*index2*nvec_eta2;
      coarse_spinwise_n_daggered_hopp_PRECISION( &eta1_pt, &eta2_pt, &phi_pt, D_pt, l );
    }
  }
}

static void coarse_spinwise_self_couplings_PRECISION( vector_PRECISION *eta1, vector_PRECISION *eta2, vector_PRECISION *phi,
							  config_PRECISION clover, int length, level_struct *l ) {
  
  int num_eig_vect      = l->num_parent_eig_vect;
  int clover_step_size1 = (num_eig_vect * (num_eig_vect+1))/2;
  int clover_step_size2 = SQUARE(num_eig_vect);
  int nvec = phi->num_vect_now, nvec_phi = phi->num_vect, nvec_eta1 = eta1->num_vect, nvec_eta2 = eta2->num_vect;

  config_PRECISION clover_pt = clover;
  buffer_PRECISION phi_pt = phi->vector_buffer, eta1_pt = eta1->vector_buffer, eta2_pt = eta2->vector_buffer+num_eig_vect*nvec_eta2, phi_end_pt = phi->vector_buffer+length*nvec_phi;
  if ( nvec_eta1 != nvec_eta2 || nvec_eta1 < nvec_phi )
    error0("coarse_spinwise_self_couplings_PRECISION: assumptions are not met\n");

  // U(x) = [ A B      , A=A*, D=D*, C = -B*
  //          C D ]
  // storage order: upper triangle of A, upper triangle of D, B, columnwise
  // diagonal coupling
  while ( phi_pt < phi_end_pt ) {
    // A
    mvp_PRECISION( eta1_pt, clover_pt, phi_pt, num_eig_vect, nvec, nvec_eta1, nvec_phi );
    clover_pt += clover_step_size1; 
    phi_pt    += num_eig_vect*nvec_phi;
    eta1_pt   += num_eig_vect*nvec_eta1;
    // D
    mvp_PRECISION( eta2_pt, clover_pt, phi_pt, num_eig_vect, nvec, nvec_eta2, nvec_phi );
    clover_pt += clover_step_size1;
    phi_pt    -= num_eig_vect*nvec_phi; 
    eta2_pt   -= num_eig_vect*nvec_eta2;
    // C = -B*
    nmvh_PRECISION( eta1_pt, clover_pt, phi_pt, num_eig_vect, nvec, nvec_eta1, nvec_phi );
    phi_pt  += num_eig_vect*nvec_phi; 
    eta1_pt += num_eig_vect*nvec_eta1;
    // B
    mv_PRECISION( eta2_pt, clover_pt, phi_pt, num_eig_vect, nvec, nvec_eta2, nvec_phi );
    clover_pt += clover_step_size2; 
    phi_pt    += num_eig_vect*nvec_phi; 
    eta2_pt   += 3*num_eig_vect*nvec_eta2;
  }
}

// eta1, eta2 <-coarse_aggregate_neighbor_couplings_PRECISION*phi (eta1&2 are initialized to 0's)
void coarse_aggregate_neighbor_couplings_PRECISION( vector_PRECISION *eta1, vector_PRECISION *eta2, vector_PRECISION *phi,
							const int mu, schwarz_PRECISION_struct *s, level_struct *l ) {
  
  int i, index1, index2;
  int length     = l->is_PRECISION.agg_boundary_length[mu];
  int *index_dir = l->is_PRECISION.agg_boundary_index[mu];
  int *neighbor  = l->is_PRECISION.agg_boundary_neighbor[mu];
  int n = l->num_lattice_site_var;
  int Dls = n*n, Dss = 4*n*n;
  int nvec = phi->num_vect_now, nvec_phi = phi->num_vect, nvec_eta1 = eta1->num_vect, nvec_eta2 = eta2->num_vect;

  vector_PRECISION eta1_pt, eta2_pt, phi_pt;
  config_PRECISION D_pt, D = s->op.D;
  
  if ( nvec_eta1 != nvec_eta2 || nvec_eta1 < nvec_phi )
    error0("coarse_aggregate_neighbor_couplings_PRECISION: assumptions are not met\n");

  vector_PRECISION_define( eta1, 0, 0, l->vector_size, l );
  vector_PRECISION_define( eta2, 0, 0, l->vector_size, l ); 

  vector_PRECISION_duplicate( &eta1_pt, eta1, 0, l);
  vector_PRECISION_duplicate( &eta2_pt, eta2, 0, l );
  vector_PRECISION_duplicate( &phi_pt, phi, 0, l );
  
  // requires the positive boundaries of phi to be communicated befor
  for ( i=0; i<length; i++ ) {
    index1 = index_dir[i];
    index2 = neighbor[i];
    D_pt = D + Dss*index1 + Dls*mu;
    phi_pt.vector_buffer  = phi->vector_buffer  + n*index2*nvec_phi; 
    eta1_pt.vector_buffer = eta1->vector_buffer + n*index1*nvec_eta1;
    eta2_pt.vector_buffer = eta2->vector_buffer + n*index1*nvec_eta2;
    coarse_spinwise_hopp_PRECISION( &eta1_pt, &eta2_pt, &phi_pt, D_pt, l );
  }
}

// 
void coarse_aggregate_block_diagonal_PRECISION( vector_PRECISION *eta1, vector_PRECISION *eta2, vector_PRECISION *phi,
						    config_PRECISION block, level_struct *l ) {
  int length          = l->inner_vector_size;
  int num_eig_vect    = l->num_parent_eig_vect;
  int block_step_size = (num_eig_vect * (num_eig_vect+1))/2;
  int nvec = phi->num_vect_now, nvec_phi = phi->num_vect, nvec_eta1 = eta1->num_vect, nvec_eta2 = eta2->num_vect;
  config_PRECISION block_pt = block;
  vector_PRECISION phi_pt, eta1_pt, eta2_pt, phi_end_pt;
  vector_PRECISION_duplicate( &eta1_pt, eta1, 0, l);
  vector_PRECISION_duplicate( &eta2_pt, eta2, 0, l);
  vector_PRECISION_duplicate( &phi_pt, phi, 0, l);
  phi_end_pt.vector_buffer = phi->vector_buffer+length*nvec_phi;

  if ( nvec_eta1 != nvec_eta2 || nvec_eta1 < nvec_phi )
    error0("coarse_aggregate_block_diagonal_PRECISION: assumptions are not met\n");

  // U(x) = [ A 0      , A=A*, D=D* 
  //          0 D ]
  // storage order: upper triangle of A, upper triangle of D, columnwise
  // diagonal coupling
  while ( phi_pt.vector_buffer < phi_end_pt.vector_buffer ) {
    // A
    mvp_PRECISION( eta1_pt.vector_buffer, block_pt, phi_pt.vector_buffer, num_eig_vect, nvec, nvec_eta1, nvec_phi );
    vector_PRECISION_define( &eta2_pt, _COMPLEX_PRECISION_ZERO, 0, num_eig_vect, l );
    block_pt += block_step_size; 
    eta1_pt.vector_buffer += num_eig_vect*nvec_eta1; 
    eta2_pt.vector_buffer += num_eig_vect*nvec_eta2; 
    phi_pt.vector_buffer  += num_eig_vect*nvec_phi;
    // D
    vector_PRECISION_define( &eta1_pt, _COMPLEX_PRECISION_ZERO, 0, num_eig_vect, l );
    mvp_PRECISION( eta2_pt.vector_buffer, block_pt, phi_pt.vector_buffer, num_eig_vect, nvec, nvec_eta2, nvec_phi );
    block_pt += block_step_size; 
    eta1_pt.vector_buffer += num_eig_vect*nvec_eta1; 
    eta2_pt.vector_buffer += num_eig_vect*nvec_eta2;
    phi_pt.vector_buffer  += num_eig_vect*nvec_phi;
  }
}

void coarse_gamma5_PRECISION( vector_PRECISION *eta, vector_PRECISION *phi, int start, int end, level_struct *l ) {
  
  int j, k = l->num_lattice_site_var/2;
  int i, jj;
  int nvec = phi->num_vect_now, nvec_phi = phi->num_vect, nvec_eta = eta->num_vect;
  buffer_PRECISION eta_end, eta_pt, phi_pt;
  eta_end = eta->vector_buffer + end*nvec_eta;
  phi_pt = phi->vector_buffer + start*nvec_phi;
  eta_pt = eta->vector_buffer + start*nvec_eta;

  if( nvec_phi != nvec_eta )
    error0("coarse_gamma5_PRECISION: #vec eta should be #vec phi\n");

  if ( eta_pt != phi_pt ) {
    while ( eta_pt < eta_end ) {
      for ( j=0; j<k; j++ ) {
        VECTOR_LOOP( i, nvec, jj, *eta_pt = -(*phi_pt); eta_pt++; phi_pt++;) eta_pt += nvec_eta-nvec; phi_pt += nvec_phi-nvec;
      }
      for ( j=0; j<k; j++ ) {
        VECTOR_LOOP( i, nvec, jj, *eta_pt = *phi_pt; eta_pt++; phi_pt++;) eta_pt += nvec_eta-nvec; phi_pt += nvec_phi-nvec;
      }
    }
  } else {
    while ( eta_pt < eta_end ) {
      for ( j=0; j<k; j++ ) {
        VECTOR_LOOP( i, nvec, jj, *eta_pt = -(*eta_pt); eta_pt++;) eta_pt += nvec_eta-nvec;
      }
      eta_pt+=k*nvec_eta;
    }
  }
}

void apply_coarse_operator_PRECISION( vector_PRECISION *eta, vector_PRECISION *phi, operator_PRECISION_struct *op,
					  level_struct *l, struct Thread *threading ) {
  PROF_PRECISION_START( _SC, threading );
  int start;
  int end;
  compute_core_start_end_custom(0, l->num_inner_lattice_sites, &start, &end, l, threading, 1);

  coarse_self_couplings_PRECISION( eta, phi, op, start, end, l);

  PROF_PRECISION_STOP( _SC, 1, threading );
  PROF_PRECISION_START( _NC, threading );

  coarse_hopping_term_PRECISION( eta, phi, op, _FULL_SYSTEM, l, threading );

  PROF_PRECISION_STOP( _NC, 1, threading );
}

/**************** TEST ROUTINES *****************************************************************************************/

// used only in the test routines??? for TM1p1
void coarse_tau1_gamma5_PRECISION( vector_PRECISION *eta, vector_PRECISION *phi, int start, int end, level_struct *l ) {
  
#ifdef HAVE_TM1p1
  if ( g.n_flavours == 2 ) {
    int j, k=l->num_lattice_site_var/4;
    buffer_PRECISION eta_end, phi_pt, eta_pt;
    
    eta_end = eta->vector_buffer + end;
    phi_pt = phi->vector_buffer + start;
    eta_pt = eta->vector_buffer + start;
    
    ASSERT( eta_pt != phi_pt );
    while ( eta_pt < eta_end ) {
      phi_pt += k;
      for ( j=0; j<k; j++ ) {
        *eta_pt = -(*phi_pt);
        eta_pt++; phi_pt++;
      }
      phi_pt -= 2*k;
      for ( j=0; j<k; j++ ) {
        *eta_pt = -(*phi_pt);
        eta_pt++; phi_pt++;
      }
      phi_pt += 2*k;
      for ( j=0; j<k; j++ ) {
        *eta_pt = *phi_pt;
        eta_pt++; phi_pt++;
      }
      phi_pt -= 2*k;
      for ( j=0; j<k; j++ ) {
        *eta_pt = *phi_pt;
        eta_pt++; phi_pt++;
      }
      phi_pt += k;
    }
  } else 
#endif
    {
      warning0("coarse_tau1_gamma5_PRECISION called with g.n_flavours != 2\n");
      coarse_gamma5_PRECISION( eta, phi, start, end, l );
    }
}

void coarse_operator_PRECISION_test_routine( level_struct *l, struct Thread *threading ) {
#if 1
  if ( !l->idle ) {
    int ivs = l->inner_vector_size, civs = l->next_level->inner_vector_size;
    int n_vect = num_loop, j, jj;
    PRECISION diff1[n_vect], diff2[n_vect];
    vector_PRECISION vp[4], vc[3];
    complex_PRECISION factor[n_vect];

    for(int i=0; i<4; i++){
      vector_PRECISION_init( &vp[i] );
      vector_PRECISION_alloc( &vp[i], _ORDINARY, n_vect, l, threading );
      vp[i].num_vect_now = n_vect;
    }

    for(int i=0; i<3; i++){
      vector_PRECISION_init( &vc[i] );
      vector_PRECISION_alloc( &vc[i], _ORDINARY, n_vect, l->next_level, threading );                       
      vc[i].num_vect_now = n_vect;
    } 

#if 1
    START_LOCKED_MASTER(threading)
/*#ifdef HAVE_TM1p1
    if(g.n_flavours == 1)
#endif*/
/*    {
      diff = global_inner_product_PRECISION( &(l->is_PRECISION.interpolation[0]), &(l->is_PRECISION.interpolation[1]), 0, ivs, l, no_threading )
        / global_norm_PRECISION( &(l->is_PRECISION.interpolation[0]), 0, ivs, l, no_threading );
      
      test0_PRECISION("depth: %d, correctness of block_gram_schmidt: %le\n", l->depth, cabs(diff) );
    }*/
    if ( !l->next_level->idle )
      vector_PRECISION_define_random( &vc[0], 0, civs, l->next_level );
    vector_PRECISION_distribute( &vc[1], &vc[0], l->next_level );
    vector_PRECISION_gather( &vc[2], &vc[1], l->next_level );
    if ( !l->next_level->idle ) {
      vector_PRECISION_minus( &vc[1], &vc[0], &vc[2], 0, civs, l->next_level );
      global_norm_PRECISION( diff1, &vc[1], 0, civs, l->next_level, no_threading );
      global_norm_PRECISION( diff2, &vc[0], 0, civs, l->next_level, no_threading );
    }
    for(int i=0; i<n_vect; i++) 
      test0_PRECISION("depth: %d, correctness of gather( distribute( phi_c ) ) : %le\n", l->depth, diff1[i]/diff2[i] );
    if ( !l->next_level->idle )
      vector_PRECISION_define_random( &vc[0], 0, civs, l->next_level );
    END_LOCKED_MASTER(threading)

    START_LOCKED_MASTER(threading)
    interpolate3_PRECISION( &vp[0], &vc[0], l, no_threading );
    restrict_PRECISION( &vc[1], &vp[0], l, no_threading );
    if ( !l->next_level->idle ) {
      vector_PRECISION_minus( &vc[2], &vc[0], &vc[1], 0, civs, l->next_level );
      global_norm_PRECISION( diff1, &vc[2], 0, civs, l->next_level, no_threading );
      global_norm_PRECISION( diff2, &vc[0], 0, civs, l->next_level, no_threading );
      for(int i=0; i<n_vect; i++)
        test0_PRECISION("depth: %d, correctness of ( P* P - 1 ) phi_c: %le\n", l->depth, abs_PRECISION(diff1[i]/diff2[i]) );
    }    
    END_LOCKED_MASTER(threading)

    if(threading->n_core>1) {
      interpolate3_PRECISION( &vp[0], &vc[0], l, threading );
      restrict_PRECISION( &vc[1], &vp[0], l, threading );
      START_LOCKED_MASTER(threading)
      if ( !l->next_level->idle ) {
        vector_PRECISION_minus( &vc[2], &vc[0], &vc[1], 0, civs, l->next_level );
        global_norm_PRECISION( diff1, &vc[2], 0, civs, l->next_level, no_threading );
        global_norm_PRECISION( diff2, &vc[0], 0, civs, l->next_level, no_threading );
        for(int i=0; i<n_vect; i++)
          test0_PRECISION("depth: %d, correctness of ( P* P - 1 ) phi_c with threading: %le\n", l->depth, diff1[i]/diff2[i] );
      }
      END_LOCKED_MASTER(threading)
    }

    START_LOCKED_MASTER(threading)
    if (l->depth==0) 
      gamma5_PRECISION( &vp[1], &vp[0], l, no_threading );
    else
      coarse_gamma5_PRECISION( &vp[1], &vp[0], 0, ivs, l );
    restrict_PRECISION( &vc[1], &vp[1], l, no_threading );
    coarse_gamma5_PRECISION( &vc[2], &vc[1], 0, civs, l->next_level );
    if ( !l->next_level->idle ) {
      vector_PRECISION_minus( &vc[1], &vc[0], &vc[2], 0, civs, l->next_level );
      global_norm_PRECISION( diff1, &vc[1], 0, civs, l->next_level, no_threading );
      global_norm_PRECISION( diff2, &vc[0], 0, civs, l->next_level, no_threading );
      for(int i=0; i<n_vect; i++)
        test0_PRECISION("depth: %d, correctness of ( g5_c P* g5 P - 1 ) phi_c: %le\n", l->depth, diff1[i]/diff2[i] );
    }    
/*#ifdef HAVE_TM1p1
    if(g.n_flavours == 2) {
      if (l->depth==0) 
        tau1_gamma5_PRECISION( &vp[1], &vp[0], l, no_threading );
      else
        coarse_tau1_gamma5_PRECISION( &vp[1], &vp[0], 0, ivs, l );
      restrict_PRECISION( &vc[1], &vp[1], l, no_threading );
      coarse_tau1_gamma5_PRECISION( &vc[2], &vc[1], 0, civs, l->next_level );
      if ( !l->next_level->idle ) {
        vector_PRECISION_minus( &vc[1], &vc[0], &vc[2], 0, civs, l->next_level );
        diff = global_norm_PRECISION( &vc[1], 0, civs, l->next_level, no_threading ) / global_norm_PRECISION( &vc[0], 0, civs, l->next_level, no_threading );
        test0_PRECISION("depth: %d, correctness of ( tau1 g5_c P* tau1 g5 P - 1 ) phi_c: %le\n", l->depth, diff );
      }    
    }
#endif*/
    END_LOCKED_MASTER(threading)

    START_LOCKED_MASTER(threading)
    vector_PRECISION_change_layout(&vp[2],&vp[1], _NVEC_OUTER, no_threading );
    vector_PRECISION_change_layout(&vp[2],&vp[2], _NVEC_INNER, no_threading );
    //for(int i=0;i<ivs*n_vect;i++)printf0("%d %d %d %13g %g %g\n",i/n_vect/l->num_lattice_site_var,(i/n_vect)%l->num_lattice_site_var,i%n_vect,creal_PRECISION(vp[1].vector_buffer[i]-vp[2].vector_buffer[i]),creal_PRECISION(vp[1].vector_buffer[i]),creal_PRECISION(vp[2].vector_buffer[i]));   
    vector_PRECISION_minus( &vp[2], &vp[2], &vp[1], 0, ivs, l );
    global_norm_PRECISION( diff1, &vp[2], 0, ivs, l, no_threading );
    global_norm_PRECISION( diff2, &vp[1], 0, ivs, l, no_threading );
    for(int i=0; i<n_vect; i++)
      test0_PRECISION("depth: %d, correctness of vector layout change: %le\n", l->depth, diff1[i]/diff2[i] );
    END_LOCKED_MASTER(threading)
      
    START_LOCKED_MASTER(threading)
      vector_PRECISION_define_random( &vc[0], 0, civs, l->next_level );//printf("1111 %d\n", g.my_rank);fflush(stdout);
    interpolate3_PRECISION( &vp[0], &vc[0], l, no_threading );//printf("11222 %d\n", g.my_rank);fflush(stdout);
    vector_PRECISION_define( &vp[1], 0, 0, ivs, l );//printf("333 %d\n", g.my_rank);fflush(stdout);
    if (l->depth==0) 
      add_diagonal_PRECISION( &vp[1], &vp[0], l->s_PRECISION.op.odd_proj, ivs );
    else
      coarse_add_block_diagonal_PRECISION( &vp[1], &vp[0], l->s_PRECISION.op.odd_proj, ivs, l );
    //    printf("1334444 %d\n", g.my_rank);fflush(stdout);
    restrict_PRECISION( &vc[1], &vp[1], l, no_threading );//printf("555 %d\n", g.my_rank);fflush(stdout);
      
    VECTOR_LOOP(j, n_vect, jj, factor[j+jj] = -1.0;)

    vector_PRECISION_scale( &vc[1], &vc[1], factor, 0, 0, civs, l->next_level );
    coarse_add_block_diagonal_PRECISION( &vc[1], &vc[0], l->next_level->s_PRECISION.op.odd_proj, civs, l->next_level );//printf("66666 %d\n", g.my_rank);fflush(stdout);
    //for(int i=0;i<civs*n_vect;i++)printf0("%d %d %13g %g %g\n",i/n_vect,i%n_vect,creal_PRECISION(vc[0].vector_buffer[i]-vc[1].vector_buffer[i]),creal_PRECISION(vc[0].vector_buffer[i]),creal_PRECISION(vc[1].vector_buffer[i]));
    //    MPI_Barrier(g.comm_cart);
    global_norm_PRECISION( diff1, &vc[1], 0, civs, l->next_level, no_threading );//printf("77 %d\n", g.my_rank);fflush(stdout);
    global_norm_PRECISION( diff2, &vc[0], 0, civs, l->next_level, no_threading );//printf("888 %d\n", g.my_rank);fflush(stdout);
    for(int i=0; i<n_vect; i++)
      test0_PRECISION("depth: %d, correctness of ( P* 1odd P - 1odd_c ) phi_c: %le\n", l->depth, diff1[i]/diff2[i] );
    END_LOCKED_MASTER(threading)  

#ifdef HAVE_TM
    START_LOCKED_MASTER(threading)
    if (g.mu + g.mu_odd_shift != 0.0 || g.mu + g.mu_even_shift != 0.0 ) {
      vector_PRECISION_define( &vp[1], 0, 0, ivs, l );
      if (l->depth==0) 
        add_diagonal_PRECISION( &vp[1], &vp[0], l->s_PRECISION.op.tm_term, ivs );
      else
        coarse_add_anti_block_diagonal_PRECISION( &vp[1], &vp[0], l->s_PRECISION.op.tm_term, ivs, l );
      restrict_PRECISION( &vc[1], &vp[1], l, no_threading );
      
      VECTOR_LOOP(j, n_vect, jj, factor[j+jj] = -g.mu_factor[l->next_level->depth]/g.mu_factor[l->depth];)

      vector_PRECISION_scale( &vc[1], &vc[1], factor, 0, 0, civs, l->next_level );
      coarse_add_anti_block_diagonal_PRECISION( &vc[1], &vc[0], l->next_level->s_PRECISION.op.tm_term, civs, l->next_level );
      global_norm_PRECISION( diff1, &vc[1], 0, civs, l->next_level, no_threading );
      global_norm_PRECISION( diff2, &vc[0], 0, civs, l->next_level, no_threading );
      for(int i=0; i<n_vect; i++)
        test0_PRECISION("depth: %d, correctness of ( P* tm P - tm_c ) phi_c: %le\n", l->depth, diff1[i]/diff2[i] );
    }
    END_LOCKED_MASTER(threading)  
#endif

/*#ifdef HAVE_TM1p1
    START_LOCKED_MASTER(threading)
    if ( g.n_flavours == 2 &&
         ( g.epsbar != 0 || g.epsbar_ig5_odd_shift != 0 || g.epsbar_ig5_odd_shift != 0 ) ) {
      vector_PRECISION_define( &vp[1], 0, 0, ivs, l );
      if (l->depth==0) 
        apply_doublet_coupling_PRECISION( &vp[1], &vp[0], l->s_PRECISION.op.epsbar_term, ivs );
      else
        coarse_add_doublet_coupling_PRECISION( &vp[1], &vp[0], l->s_PRECISION.op.epsbar_term, ivs, l );
      restrict_PRECISION( &vc[1], &vp[1], l, no_threading );
      
      vector_PRECISION_scale( &vc[1], &vc[1], -g.epsbar_factor[l->next_level->depth]/g.epsbar_factor[l->depth], 0, civs, l->next_level );
      coarse_add_doublet_coupling_PRECISION( &vc[1], &vc[0], l->next_level->s_PRECISION.op.epsbar_term, civs, l->next_level );
      diff = global_norm_PRECISION( &vc[1], 0, civs, l->next_level, no_threading ) / global_norm_PRECISION( &vc[0], 0, civs, l->next_level, no_threading );
      test0_PRECISION("depth: %d, correctness of ( P* eps P - eps_c ) phi_c: %le\n", l->depth, diff );
    }
    END_LOCKED_MASTER(threading)  
#endif*/

    if ( l->level > 0 ) {
      START_LOCKED_MASTER(threading)
      interpolate3_PRECISION( &vp[0], &vc[0], l, no_threading );
      apply_operator_PRECISION( &vp[1], &vp[0], &(l->p_PRECISION), l, no_threading );      //l==0=>d_plus_clover
      
#ifdef HAVE_TM
      if (g.mu + g.mu_odd_shift != 0.0 || g.mu + g.mu_even_shift != 0.0 )
        if (g.mu_factor[l->depth] != g.mu_factor[l->next_level->depth]) {  
          VECTOR_LOOP(j, n_vect, jj, factor[j+jj] = (g.mu_factor[l->next_level->depth]/g.mu_factor[l->depth])-1.;)
          vector_PRECISION_scale( &vp[2], &vp[0], factor, 0, 0, ivs, l );
          if(l->depth == 0)
            add_diagonal_PRECISION( &vp[1], &vp[2], l->p_PRECISION.op->tm_term, ivs );//tm==0!!!!!
          else
            coarse_add_anti_block_diagonal_PRECISION( &vp[1], &vp[2], l->p_PRECISION.op->tm_term, ivs, l );
        }
#endif
      //      printf0("eval me 4 depth %d\n",l->depth);
      /*
      if (g.method == 2) {
	l->next_level->p_PRECISION.eval_operator = g.odd_even?coarse_apply_schur_complement_PRECISION:apply_coarse_operator_PRECISION;
	if ( g.odd_even )
	  l->next_level->p_PRECISION.op = &(l->next_level->oe_op_PRECISION);
	else
	  l->next_level->p_PRECISION.op = &(l->next_level->s_PRECISION.op);
	  }*/
      restrict_PRECISION( &vc[1], &vp[1], l, no_threading );
      if ( !l->next_level->idle ) {
        if ( l->level==1 && g.odd_even ) {
	    coarse_odd_even_PRECISION_test( &vc[2], &vc[0], l->next_level, no_threading );
	}
        else
          apply_operator_PRECISION( &vc[2], &vc[0], &(l->next_level->p_PRECISION), l->next_level, no_threading );

	//for(int i=0;i<civs*n_vect;i++)printf0("%d %d %d %13g %g %g\n",i/n_vect/l->next_level->num_lattice_site_var,(i/n_vect)%l->next_level->num_lattice_site_var,i%n_vect,creal_PRECISION(vc[1].vector_buffer[i]-vc[2].vector_buffer[i]),creal_PRECISION(vc[1].vector_buffer[i]),creal_PRECISION(vc[2].vector_buffer[i]));
        
        vector_PRECISION_minus( &vc[2], &vc[1], &vc[2], 0, civs, l->next_level );
        global_norm_PRECISION( diff1, &vc[2], 0, civs, l->next_level, no_threading );
        global_norm_PRECISION( diff2, &vc[1], 0, civs, l->next_level, no_threading );
																						    
        if ( l->level==1 && g.odd_even ) {
          for(int i=0; i<n_vect; i++) 
            test0_PRECISION("depth: %d, correctness of odd even preconditioned ( P* D P - D_c ) phi_c: %le\n", l->depth, diff1[i]/diff2[i] );
        } else {
          for(int i=0; i<n_vect; i++)
            test0_PRECISION("depth: %d, correctness of ( P* D P - D_c ) phi_c: %le\n", l->depth, diff1[i]/diff2[i] );
        }
      }      
      END_LOCKED_MASTER(threading)

      if(threading->n_core>1) {
        if ( !l->next_level->idle ) {
          if ( l->level==1 && g.odd_even )
            coarse_odd_even_PRECISION_test( &vc[2], &vc[0], l->next_level, threading );
          else
            apply_operator_PRECISION( &vc[2], &vc[0], &(l->next_level->p_PRECISION), l->next_level, threading );
        }
        START_LOCKED_MASTER(threading)
        if ( !l->next_level->idle ) {
          vector_PRECISION_minus( &vc[2], &vc[1], &vc[2], 0, civs, l->next_level );
          global_norm_PRECISION( diff1, &vc[2], 0, civs, l->next_level, no_threading );
          global_norm_PRECISION( diff2, &vc[1], 0, civs, l->next_level, no_threading );
          if ( l->level==1 && g.odd_even ) { //TODO: this test doesn't work without SSE!!
            for(int i=0; i<n_vect; i++)
              test0_PRECISION("depth: %d, correctness of odd even preconditioned ( P* D P - D_c ) phi_c with D_c threaded: %le\n", l->depth, diff1[i]/diff2[i] );
          } else {
            for(int i=0; i<n_vect; i++)
              test0_PRECISION("depth: %d, correctness of ( P* D P - D_c ) phi_c with D_c threaded: %le\n", l->depth, diff1[i]/diff2[i] );
          }
        }
        END_LOCKED_MASTER(threading)
        }
    }
	
    START_LOCKED_MASTER(threading)
    if ( l->level > 0 && l->depth > 0 && g.method == 3 && g.odd_even ) {
      vector_PRECISION_define_random( &vp[0], 0, ivs, l );
      block_to_oddeven_PRECISION( &vp[3], &vp[0], l, no_threading );
      coarse_diag_ee_PRECISION( &vp[2], &vp[3], &(l->oe_op_PRECISION), l, no_threading );
      coarse_diag_oo_PRECISION( &vp[2], &vp[3], &(l->oe_op_PRECISION), l, no_threading );
      coarse_hopping_term_PRECISION( &vp[2], &vp[3], &(l->oe_op_PRECISION), _FULL_SYSTEM, l, no_threading );
      oddeven_to_block_PRECISION( &vp[3], &vp[2], l, no_threading );
      apply_operator_PRECISION( &vp[1], &vp[0], &(l->p_PRECISION), l, no_threading );
      vector_PRECISION_minus( &vp[3], &vp[3], &vp[1], 0, ivs, l );
      global_norm_PRECISION( diff1, &vp[3], 0, ivs, l, no_threading );
      global_norm_PRECISION( diff2, &vp[1], 0, ivs, l, no_threading );
      for(int i=0; i<n_vect; i++)
        test0_PRECISION("depth: %d, correctness of odd even layout (smoother): %le\n", l->depth, diff1[i]/diff2[i] );
     
      block_to_oddeven_PRECISION( &vp[3], &vp[0], l, no_threading );
      coarse_odd_even_PRECISION_test( &vp[2], &vp[3], l, no_threading );
      oddeven_to_block_PRECISION( &vp[3], &vp[2], l, no_threading );
      apply_operator_PRECISION( &vp[1], &vp[0], &(l->p_PRECISION), l, no_threading );
      vector_PRECISION_minus( &vp[3], &vp[3], &vp[1], 0, ivs, l );
      global_norm_PRECISION( diff1, &vp[3], 0, ivs, l, no_threading );
      global_norm_PRECISION( diff2, &vp[1], 0, ivs, l, no_threading );
      for(int i=0; i<n_vect; i++)   
        test0_PRECISION("depth: %d, correctness of odd even preconditioned operator (smoother): %le\n", l->depth, diff1[i]/diff2[i] );
    }

    END_LOCKED_MASTER(threading)    
#endif

    for(int i=0; i<4; i++)
      vector_PRECISION_free( &vp[i], l, threading );
    for(int i=0; i<3; i++)
      vector_PRECISION_free( &vc[i], l->next_level, threading );
#if 1    
    if ( l->next_level->level > 0  && !l->next_level->idle ) {
      schwarz_PRECISION_mvm_testfun( &(l->next_level->s_PRECISION), l->next_level, threading );
    }

    if ( l->next_level->level > 0 && !l->next_level->idle )
      coarse_operator_PRECISION_test_routine( l->next_level, threading );
#endif

  }
#endif
}

