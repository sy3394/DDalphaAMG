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
 * copied: 11/27/2019 (US system)
 * changed from sbacchio
 * checked: 12/03/2019
 * checked: 12/08/2019
 * glanced over: 12/18/2019
 * confirmed: not much different from sbacchio:01/02/2020
 */

#include "main.h"

static void operator_PRECISION_alloc_projection_buffers( operator_PRECISION_struct *op, level_struct *l );
static void operator_PRECISION_free_projection_buffers( operator_PRECISION_struct *op, level_struct *l );


void operator_PRECISION_init( operator_PRECISION_struct *op ) {
  
  op->prnT = NULL;
  op->pr_num_vect = 0;
  op->index_table = NULL;
  op->neighbor_table = NULL;
  op->backward_neighbor_table = NULL;
  op->translation_table = NULL;
  op->D = NULL;
  op->clover = NULL;
  op->clover_oo_inv = NULL;
  op->m0 = 0;
#ifdef HAVE_TM
  op->mu = 0;
  op->mu_even_shift = 0;
  op->mu_odd_shift = 0;
  op->odd_proj = NULL;
  op->tm_term = NULL;
#endif
#ifdef HAVE_TM1p1
  op->epsbar = 0;
  op->epsbar_ig5_even_shift = 0;
  op->epsbar_ig5_odd_shift = 0;
  op->epsbar_term = NULL;
  op->clover_doublet_oo_inv = NULL;
#endif
  
  for ( int mu=0; mu<4; mu++ )
    op->config_boundary_table[mu] = NULL;
  
  for ( int i=0; i<8; i++ ) {
    op->c.boundary_table[i] = NULL;
    op->c.buffer[i] = NULL;
    op->c.in_use[i] = 0;
  }
  op->c.comm = 1;
  op->buffer = NULL;
}

void operator_PRECISION_alloc( operator_PRECISION_struct *op, const int type, level_struct *l ) {
  
  /*********************************************************************************
   * Allocates space for the operator struc and define index tables and comminicator for its application
   * - operator_PRECISION_struct *op: operator struct for which space is allocated.
   * - const int type: Defines the data layout type of the operator.
   *                   Possible values are: { _ORDINARY, _SCHWARZ }
   *********************************************************************************/

  int mu, nu, its = 1, its_boundary, nls, clover_site_size, coupling_site_size;
  
  //------------- allocate memory for op->D, op->clove, op->odd_proj, (op->tm_term if HAVE_TM,  op->epsbar_term if HAVE_TM1p1)
  if ( l->depth == 0 ) {
    clover_site_size = 42;//??????????
    coupling_site_size = 4*9; // coupling term in D_w
  } else {
    clover_site_size = (l->num_lattice_site_var*(l->num_lattice_site_var+1))/2;
    coupling_site_size = 4*l->num_lattice_site_var*l->num_lattice_site_var;
  }

  nls = (type==_SCHWARZ) ? (2*l->num_lattice_sites-l->num_inner_lattice_sites):l->num_inner_lattice_sites;

  MALLOC( op->D, complex_PRECISION, coupling_site_size*nls );
  MALLOC( op->clover, complex_PRECISION, clover_site_size*l->num_inner_lattice_sites );

  int block_site_size = ( l->depth == 0 ) ? 12 : (l->num_lattice_site_var/2*(l->num_lattice_site_var/2+1));//??????
  MALLOC( op->odd_proj, complex_PRECISION, block_site_size*l->num_inner_lattice_sites );
#ifdef HAVE_TM
  MALLOC( op->tm_term, complex_PRECISION, block_site_size*l->num_inner_lattice_sites );
#endif
#ifdef HAVE_TM1p1
  MALLOC( op->epsbar_term, complex_PRECISION, block_site_size*l->num_inner_lattice_sites );
#endif

  if ( type == _SCHWARZ && l->depth == 0 && g.odd_even ) {
    if( g.csw ) {
#ifdef HAVE_TM //we use LU here
      MALLOC( op->clover_oo_inv, complex_PRECISION, 72*(l->num_inner_lattice_sites/2+1) );
#else
      MALLOC( op->clover_oo_inv, complex_PRECISION, clover_site_size*(l->num_inner_lattice_sites/2+1) );
#endif
    }
#ifdef HAVE_TM1p1
    MALLOC( op->clover_doublet_oo_inv, complex_PRECISION, 12*12*2*(l->num_inner_lattice_sites/2+1) );
#endif
  }  

  //---------------- allocate memory for index tables
  if ( type ==_SCHWARZ ) {
    its_boundary = 2; // postive and negative boundaries???
  } else {
    its_boundary = 1; // only one side boundary???
  }
  for ( mu=0; mu<4; mu++ ) {
    its *= (l->local_lattice[mu]+its_boundary); //its:index table size
  }

  MALLOC( op->index_table, int, its );
  if ( type ==_ODDEVEN ) {
    MALLOC( op->neighbor_table, int, 5*its );
    MALLOC( op->backward_neighbor_table, int, 5*its );
  } else {
    MALLOC( op->neighbor_table, int, (l->depth==0?4:5)*l->num_inner_lattice_sites );
    MALLOC( op->backward_neighbor_table, int, (l->depth==0?4:5)*l->num_inner_lattice_sites );
  }
  MALLOC( op->translation_table, int, l->num_inner_lattice_sites );

  //-------------- allocate memory for comm_PRECISION_struc (this is where comm_PRECISION_struc is allocated!!!!!)
  if ( type != _ODDEVEN )
    operator_PRECISION_alloc_projection_buffers( op, l );
  
  ghost_alloc_PRECISION( 0, &(op->c), l );

  for ( mu=0; mu<4; mu++ ) {
    // its: #boundary sites in the mu dir
    its = 1;
    for ( nu=0; nu<4; nu++ ) {
      if ( mu != nu ) {
        its *= l->local_lattice[nu];
      }
    }
    op->c.num_boundary_sites[2*mu] = its;
    op->c.num_boundary_sites[2*mu+1] = its;
    MALLOC( op->c.boundary_table[2*mu], int, its );
    if ( type == _SCHWARZ ) {
      MALLOC( op->c.boundary_table[2*mu+1], int, its );
      MALLOC( op->config_boundary_table[mu], int, its );
    } else {
      op->c.boundary_table[2*mu+1] = op->c.boundary_table[2*mu];
    }
  }
}

void operator_PRECISION_free( operator_PRECISION_struct *op, const int type, level_struct *l ) {
  
  int mu, nu, its = 1, clover_site_size, coupling_site_size;
  
  if ( l->depth == 0 ) {
    clover_site_size = 42;
    coupling_site_size = 4*9;
  } else {
    clover_site_size = (l->num_lattice_site_var*(l->num_lattice_site_var+1))/2;
    coupling_site_size = 4*l->num_lattice_site_var*l->num_lattice_site_var;
  }
  
  int its_boundary;
  if ( type ==_SCHWARZ ) {
    its_boundary = 2;
  } else {
    its_boundary = 1;
  }
  for ( mu=0; mu<4; mu++ ) {
    its *= (l->local_lattice[mu]+its_boundary);
  }
  
  int nls = (type==_SCHWARZ) ? (2*l->num_lattice_sites-l->num_inner_lattice_sites) : l->num_inner_lattice_sites;
  FREE( op->D, complex_PRECISION, coupling_site_size*nls );
  FREE( op->clover, complex_PRECISION, clover_site_size*l->num_inner_lattice_sites );

  int block_site_size = ( l->depth == 0 ) ? 12 : (l->num_lattice_site_var/2*(l->num_lattice_site_var/2+1));
  FREE( op->odd_proj, complex_PRECISION, block_site_size*l->num_inner_lattice_sites );
#ifdef HAVE_TM
  FREE( op->tm_term, complex_PRECISION, block_site_size*l->num_inner_lattice_sites );
#endif
  if ( type == _SCHWARZ && l->depth == 0 && g.odd_even ) {
    if( g.csw ) {
#ifdef HAVE_TM //we use LU here
      FREE( op->clover_oo_inv, complex_PRECISION, 72*(l->num_inner_lattice_sites/2+1) );
#else
      FREE( op->clover_oo_inv, complex_PRECISION, clover_site_size*(l->num_inner_lattice_sites/2+1) );
#endif
    }
#ifdef HAVE_TM1p1
    FREE( op->clover_doublet_oo_inv, complex_PRECISION, 12*12*2*(l->num_inner_lattice_sites/2+1) );
#endif
  }  

#ifdef HAVE_TM1p1
  FREE( op->epsbar_term, complex_PRECISION, block_site_size*l->num_inner_lattice_sites );
#endif
  FREE( op->index_table, int, its );
  if ( type ==_ODDEVEN ) {
    FREE( op->neighbor_table, int, 5*its );
    FREE( op->backward_neighbor_table, int, 5*its );
  } else {
    FREE( op->neighbor_table, int, (l->depth==0?4:5)*l->num_inner_lattice_sites );
    FREE( op->backward_neighbor_table, int, (l->depth==0?4:5)*l->num_inner_lattice_sites );
  }
  FREE( op->translation_table, int, l->num_inner_lattice_sites );
  
  if ( type != _ODDEVEN )
    operator_PRECISION_free_projection_buffers( op, l );
  
  ghost_free_PRECISION( &(op->c), l );
  
  for ( mu=0; mu<4; mu++ ) {
    its = 1;
    for ( nu=0; nu<4; nu++ ) {
      if ( mu != nu ) {
        its *= l->local_lattice[nu];
      }
    }
    
    FREE( op->c.boundary_table[2*mu], int, its );
    if ( type == _SCHWARZ ) {
      FREE( op->c.boundary_table[2*mu+1], int, its );
      FREE( op->config_boundary_table[mu], int, its );
    } else {
      op->c.boundary_table[2*mu+1] = NULL;
    }
  }
}

static void operator_PRECISION_alloc_projection_buffers( operator_PRECISION_struct *op, level_struct *l ) {//!!!!!!!!

  // when used as preconditioner we usually do not need the projection buffers, unless
  // g.method >= 5: then oddeven_setup_float() is called in init.c, method_setup(). ???????
  int n_vect = num_loop;//(g.num_rhs_vect < l->num_eig_vect)? l->num_eig_vect:g.num_rhs_vect;//g.num_vect_now;//!!!!!!!g.num_rhs_vect
  if ( l->depth == 0 ) {
    int its = (l->num_lattice_site_var/2)*l->num_lattice_sites*n_vect; // half of spinor d.o.f. projected out
#ifdef HAVE_TM1p1
    its *= 2;
#endif
    MALLOC( op->prnT, complex_PRECISION, its*8 );
    op->prnZ = op->prnT + its; op->prnY = op->prnZ + its; op->prnX = op->prnY + its;
    op->prpT = op->prnX + its; op->prpZ = op->prpT + its; op->prpY = op->prpZ + its; op->prpX = op->prpY + its;
  }
  op->pr_num_vect = n_vect;
}

static void operator_PRECISION_free_projection_buffers( operator_PRECISION_struct *op, level_struct *l ) {

  if ( l->depth == 0 ) {
    int n_vect = op->pr_num_vect;
    int its = (l->num_lattice_site_var/2)*l->num_lattice_sites*n_vect;
#ifdef HAVE_TM1p1
    its *= 2;
#endif
    FREE( op->prnT, complex_PRECISION, its*8 );
  }
}

void operator_PRECISION_define( operator_PRECISION_struct *op, level_struct *l ) {
  /**********************************
   * Define tables (table_dim, index_table, neighbor_table, backward_neighbor_table
   * c.boundary_table, translation_table in *op where c is comm_PRECISION_struct).
   *********************************/
  
  int i, mu, t, z, y, x, *it = op->index_table,
    ls[4], le[4], l_st[4], l_en[4], *dt = op->table_dim;
  
  for ( mu=0; mu<4; mu++ ) {
    dt[mu] = l->local_lattice[mu]+1; //only positive boundaries are considered
    ls[mu] = 0;//are these necessary?????
    le[mu] = ls[mu] + l->local_lattice[mu];
    l_st[mu] = ls[mu];
    l_en[mu] = le[mu];
  }
  
  // Define index table: it[lex(site)] = local_lex_index
  // In the op array, sites are ordered from x(fastest) to t(slowest)
  // with inner cuboid first and then +T,+Z,+Y,+X boundaries
  // lexicographic inner cuboid and //?????why inner???
  // lexicographic +T,+Z,+Y,+X boundaries
  i=0;
  // inner hyper cuboid
  for ( t=ls[T]; t<le[T]; t++ )
    for ( z=ls[Z]; z<le[Z]; z++ )
      for ( y=ls[Y]; y<le[Y]; y++ )
        for ( x=ls[X]; x<le[X]; x++ ) {
          it[ lex_index( t, z, y, x, dt ) ] = i; i++;
        }
  // positive boundaries (buffers)
  for ( mu=0; mu<4; mu++ ) {
    l_st[mu] = le[mu];
    l_en[mu] = le[mu]+1;
    
    for ( t=l_st[T]; t<l_en[T]; t++ )
      for ( z=l_st[Z]; z<l_en[Z]; z++ )
        for ( y=l_st[Y]; y<l_en[Y]; y++ )
          for ( x=l_st[X]; x<l_en[X]; x++ ) {
            it[ lex_index( t, z, y, x, dt ) ] = i; i++;
          }
          
    l_st[mu] = ls[mu];
    l_en[mu] = le[mu];
  }
  
  //--- define (for the application of the entire operator)
  //      s->op.neighbor_table: local_lex_index+mu -> local_lex_index of the neighbor in the pos mu dir  
  //      s->op.backward_neighbor_table: local_lex_index+mu -> local_lex_index of the neighbor in the neg mu dir
  //      s->op.c.boundary_table: iter_index of inner bondaray sites -> local_lex_index
  //      s->op.translation_table: local_lex_index -> local_lex_index 
  define_nt_bt_tt( op->neighbor_table, op->backward_neighbor_table, op->c.boundary_table, op->translation_table, it, dt, l );
}

/********************  TEST ROUTINES *****************************************************************************/

void operator_PRECISION_test_routine( operator_PRECISION_struct *op, level_struct *l, struct Thread *threading ) {

/*********************************************************************************
* Checks for correctness of operator data layout by doing:
* - Applying D_W in double precision to a double vector.
* - Translates the same vector into PRECISION and apply D_W in PRECISION to this
*   vector and translate it back to double
* - Compare solutions ( Difference should be close to 0 ).
* If enabled, also tests odd even preconditioning.
*********************************************************************************/ 

  int ivs = l->inner_vector_size, n_vect=num_loop;
  double diff, diff1[n_vect], diff2[n_vect];
  
  vector_double vd[4];
  vector_PRECISION vp[2];

  for(int i=0; i<4; i++){
    vector_double_init( &vd[i] );
    vector_double_alloc( &vd[i], _INNER, n_vect, l, threading );
    vd[i].num_vect_now = n_vect;
  }
  
  for(int i=0; i<2; i++){
    vector_PRECISION_init( &vp[i] );
    vector_PRECISION_alloc( &vp[i], _INNER, n_vect, l, threading );
    vp[i].num_vect_now = n_vect;
  }

  START_LOCKED_MASTER(threading)
  
  vector_double_define_random_new( &vd[0], 0, l->inner_vector_size, l ); 
  apply_operator_double( &vd[1], &vd[0], &(g.p), l, no_threading );

  trans_PRECISION_new( &vp[0], &vd[0], op->translation_table, l, no_threading );
  apply_operator_PRECISION( &vp[1], &vp[0], &(l->p_PRECISION), l, no_threading );
  trans_back_PRECISION_new( &vd[2], &vp[1], op->translation_table, l, no_threading );
  
  vector_double_minus_new( &vd[3], &vd[2], &vd[1], 0, l->inner_vector_size, l );
  global_norm_double_new( diff1, &vd[3], 0, ivs, l, no_threading );
  global_norm_double_new( diff2, &vd[2], 0, ivs, l, no_threading );
  
  for(int i=0; i<n_vect; i++)
    test0_PRECISION("depth: %d, correctness of schwarz PRECISION Dirac operator: %le\n", l->depth, diff1[i]/diff2[i] );
  END_LOCKED_MASTER(threading)

  if(threading->n_core > 1) {
    apply_operator_PRECISION( &vp[1], &vp[0], &(l->p_PRECISION), l, threading );
    //SYNC_MASTER_TO_ALL(threading)
    //SYNC_CORES(threading)

    START_LOCKED_MASTER(threading)
    trans_back_PRECISION_new( &vd[2], &vp[1], op->translation_table, l, no_threading );
    vector_double_minus_new( &vd[3], &vd[2], &vd[1], 0, l->inner_vector_size, l );
    global_norm_double_new( diff1, &vd[3], 0, ivs, l, no_threading );
    global_norm_double_new( diff2, &vd[2], 0, ivs, l, no_threading );

    if ( diff > EPS_PRECISION )
      printf0("\x1b[31m");
    for(int i=0; i<n_vect; i++)
      printf0("depth: %d, correctness of schwarz PRECISION Dirac operator with threading: %le\n", l->depth, diff1[i]/diff2[i] );
    if ( diff > EPS_PRECISION )
      printf0("\x1b[0m");
    if(diff > g.test) g.test = diff;

    END_LOCKED_MASTER(threading) 
  }    

  for(int i=0; i<4; i++){
    vector_double_free( &vd[i], l, threading );
  }

  for(int i=0; i<2; i++){
    vector_PRECISION_free( &vp[i], l, threading );
  }

  START_LOCKED_MASTER(threading)
  if ( g.method >= 5 && g.odd_even )
    oddeven_PRECISION_test_new( l );
  END_LOCKED_MASTER(threading) 
}
