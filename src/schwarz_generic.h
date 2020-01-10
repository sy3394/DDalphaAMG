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
 * copied:11/29/2019
 * changed from sbacchio
 * checked: 12/08/2019
 * 1st cleanup:12/18/2019
 * confirmed:not much diff from sbacchio:01/03/2020
 */

#ifndef SCHWARZ_PRECISION_HEADER
  #define SCHWARZ_PRECISION_HEADER

struct Thread;
  
  void smoother_PRECISION_def( level_struct *l );
  void smoother_PRECISION_free( level_struct *l );
  void smoother_PRECISION_new( vector_PRECISION *phi, vector_PRECISION *Dphi, vector_PRECISION *eta,
			       int n, const int res, level_struct *l, struct Thread *threading );    
  
  void schwarz_PRECISION_init( schwarz_PRECISION_struct *s, level_struct *l );
  void schwarz_PRECISION_alloc( schwarz_PRECISION_struct *s, level_struct *l );
  void schwarz_PRECISION_free( schwarz_PRECISION_struct *s, level_struct *l );
  void schwarz_PRECISION_def( schwarz_PRECISION_struct *s, operator_double_struct *op, level_struct *l );
  void schwarz_layout_PRECISION_define( schwarz_PRECISION_struct *s, level_struct *l );
  void schwarz_PRECISION_setup( schwarz_PRECISION_struct *s, operator_double_struct *op_in, level_struct *l );//used in DDalphaAMG_interface

  void schwarz_PRECISION_boundary_update( schwarz_PRECISION_struct *s, level_struct *l );
  void block_PRECISION_boundary_op_new( vector_PRECISION *eta, vector_PRECISION *phi, int k,
                                    schwarz_PRECISION_struct *s, level_struct *l );
  void n_block_PRECISION_boundary_op_new( vector_PRECISION *eta, vector_PRECISION *phi, int k,
                                      schwarz_PRECISION_struct *s, level_struct *l );
  void coarse_block_PRECISION_boundary_op_new( vector_PRECISION *eta, vector_PRECISION *phi,
                                           int k, schwarz_PRECISION_struct *s, level_struct *l );
  void n_coarse_block_PRECISION_boundary_op_new( vector_PRECISION *eta, vector_PRECISION *phi,
						 int k, schwarz_PRECISION_struct *s, level_struct *l );
  
  void schwarz_PRECISION_new( vector_PRECISION *phi, vector_PRECISION *D_phi, vector_PRECISION *eta, const int cycles, int res, 
			      schwarz_PRECISION_struct *s, level_struct *l, struct Thread *threading );
  void additive_schwarz_PRECISION_new( vector_PRECISION *phi, vector_PRECISION *D_phi, vector_PRECISION *eta, const int cycles, int res, 
                                   schwarz_PRECISION_struct *s, level_struct *l, struct Thread *threading );
  void red_black_schwarz_PRECISION_new( vector_PRECISION *phi, vector_PRECISION *D_phi, vector_PRECISION *eta, const int cycles, int res,
                                    schwarz_PRECISION_struct *s, level_struct *l, struct Thread *threading );
  void sixteen_color_schwarz_PRECISION_new( vector_PRECISION *phi, vector_PRECISION *D_phi, vector_PRECISION *eta, const int cycles, int res, 
                                        schwarz_PRECISION_struct *s, level_struct *l, struct Thread *threading );
  
  void schwarz_PRECISION_mvm_testfun_new( schwarz_PRECISION_struct *s, level_struct *l, struct Thread *threading );
  
  static inline int connect_link_PRECISION( int t, int z, int y, int x, int mu, int dir, int *dt, int *it, 
                                            schwarz_PRECISION_struct *s, level_struct *l ) {
    
    int coord[4];
    coord[T]=t; coord[Z]=z; coord[Y]=y; coord[X]=x;
    coord[mu]+=dir;
    if ( l->global_splitting[mu] > 1 ) {
      return site_mod_index( coord[T], coord[Z], coord[Y], coord[X], dt, it );
    } else {
      coord[mu] = (coord[mu]+l->local_lattice[mu])%l->local_lattice[mu];// This is mod as in site_mod_index
      return site_index( coord[T], coord[Z], coord[Y], coord[X], dt, it );
    }
  }

#endif 
