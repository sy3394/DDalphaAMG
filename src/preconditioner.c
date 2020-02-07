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
 * glanced over: 12/08/2019
 * 1st cleanup:12/18/2019
 */

#include "main.h"
#include "preconditioner.h"

// chosen as preconditioner only for solver part and called only at the top level
void preconditioner_new( vector_double *phi, vector_double *Dphi, vector_double *eta,
                      const int res, level_struct *l, struct Thread *threading ) {

  if ( g.method == 0 )
    vector_double_copy_new( phi, eta, threading->start_index[l->depth], threading->end_index[l->depth], l );
  else if ( g.method < 5 || !g.odd_even ) {
    if ( g.mixed_precision ) {
      START_LOCKED_MASTER(threading)// my addition
      for ( int i=0; i<2; i++ ) 
	l->sbuf_float[i].num_vect_now = eta->num_vect_now;//g.num_vect_now;//!!!!!!!!!!!
      END_LOCKED_MASTER(threading)
      trans_float_new( &(l->sbuf_float[0]), eta, l->s_float.op.translation_table, l, threading );
      vcycle_float_new( &(l->sbuf_float[1]), NULL, &(l->sbuf_float[0]), res, l, threading );
      trans_back_float_new( phi, &(l->sbuf_float[1]), l->s_float.op.translation_table, l, threading );
    } else {
      START_LOCKED_MASTER(threading)// my addition  
      for ( int i=0; i<2; i++ ) 
	l->sbuf_double[i].num_vect_now = eta->num_vect_now;//g.num_vect_now;//!!!!!!!!!!!
      END_LOCKED_MASTER(threading)
      trans_double_new( &(l->sbuf_double[0]), eta, l->s_double.op.translation_table, l, threading );
      vcycle_double_new( &(l->sbuf_double[1]), NULL, &(l->sbuf_double[0]), res, l, threading );
      trans_back_double_new( phi, &(l->sbuf_double[1]), l->s_double.op.translation_table, l, threading );
    }
  } else {
    if ( g.mixed_precision ) {
      START_LOCKED_MASTER(threading)
      l->sp_float.x.num_vect_now = eta->num_vect_now;//g.num_vect_now;//!!!!!!!!
      l->sp_float.b.num_vect_now = phi->num_vect_now;//g.num_vect_now;//!!!!!!!!
      l->sp_float.num_restart = l->n_cy;
      l->sp_float.initial_guess_zero = res;
      END_LOCKED_MASTER(threading)
      serial_to_oddeven_float_new( &(l->sp_float.b), eta, l, threading );
      solve_oddeven_float_new( &(l->sp_float), &(l->oe_op_float), l, threading );
      oddeven_to_serial_float_new( phi, &(l->sp_float.x), l, threading );
    } else {
      START_LOCKED_MASTER(threading)
      l->sp_double.x.num_vect_now = eta->num_vect_now;//g.num_vect_now;//!!!!!!!!
      l->sp_double.b.num_vect_now = phi->num_vect_now;//g.num_vect_now;//!!!!!!!!
      l->sp_double.num_restart = l->n_cy;
      l->sp_double.initial_guess_zero = res;
      END_LOCKED_MASTER(threading)
      serial_to_oddeven_double_new( &(l->sp_double.b), eta, l, threading );
      solve_oddeven_double_new( &(l->sp_double), &(l->oe_op_double), l, threading );
      oddeven_to_serial_double_new( phi, &(l->sp_double.x), l, threading );
    }
    
  }
  ASSERT( g.mixed_precision != 2 );
  ASSERT( Dphi == NULL );
}
