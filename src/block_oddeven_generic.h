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
 * 12/08/2019: some reordering might be needed
 */
 
#ifndef BLOCK_ODDEVEN_PRECISION_HEADER
  #define BLOCK_ODDEVEN_PRECISION_HEADER

  struct Thread;

  void schwarz_PRECISION_oddeven_setup( schwarz_PRECISION_struct *s, level_struct *l );
  
//  void apply_block_schur_complement_PRECISION( vector_PRECISION *out, vector_PRECISION *in, int start,
//                                               schwarz_PRECISION_struct *s, level_struct *l, struct Thread *threading );
  void apply_block_schur_complement_PRECISION( vector_PRECISION *out, vector_PRECISION *in, int start,
                                               schwarz_PRECISION_struct *s, level_struct *l, struct Thread *threading );

//  void block_hopping_term_PRECISION( vector_PRECISION *eta, vector_PRECISION *phi,
//                                     int start, int amount, schwarz_PRECISION_struct *s, level_struct *l, struct Thread *threading );
  void block_hopping_term_PRECISION( vector_PRECISION *eta, vector_PRECISION *phi,
                                     int start, int amount, schwarz_PRECISION_struct *s, level_struct *l, struct Thread *threading );
//  void block_n_hopping_term_PRECISION( vector_PRECISION *eta, vector_PRECISION *phi, 
//                                       int start, int amount, schwarz_PRECISION_struct *s, level_struct *l, struct Thread *threading );
  void block_n_hopping_term_PRECISION( vector_PRECISION *eta, vector_PRECISION *phi, 
                                       int start, int amount, schwarz_PRECISION_struct *s, level_struct *l, struct Thread *threading );


//  void block_solve_oddeven_PRECISION( vector_PRECISION *phi, vector_PRECISION *r, vector_PRECISION *latest_iter,
//                                      int start, schwarz_PRECISION_struct *s, level_struct *l, struct Thread *threading );
  void block_solve_oddeven_PRECISION( vector_PRECISION *phi, vector_PRECISION *r, vector_PRECISION *latest_iter,
                                      int start, schwarz_PRECISION_struct *s, level_struct *l, struct Thread *threading );

  void oddeven_to_block_PRECISION( vector_PRECISION *out, vector_PRECISION *in, level_struct *l, struct Thread *threading );
  void block_to_oddeven_PRECISION( vector_PRECISION *out, vector_PRECISION *in, level_struct *l, struct Thread *threading );

//  void block_oddeven_PRECISION_test( level_struct *l, struct Thread *threading );
  void block_oddeven_PRECISION_test( level_struct *l, struct Thread *threading );

#endif
