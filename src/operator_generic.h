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
 * not changed from sbacchio
 * checked: 12/08/2019
 * glanced over:12/18/2019
 */

#ifndef OPERATOR_PRECISION_HEADER
  #define OPERATOR_PRECISION_HEADER

  struct Thread;
  
  void operator_PRECISION_init( operator_PRECISION_struct *op );
  void operator_PRECISION_alloc( operator_PRECISION_struct *op, const int type, level_struct *l );
  void operator_PRECISION_free( operator_PRECISION_struct *op, const int type, level_struct *l );
  void operator_PRECISION_define( operator_PRECISION_struct *op, level_struct *l );

void operator_PRECISION_set_couplings( operator_PRECISION_struct *op, level_struct *l );// no effect
  void operator_PRECISION_set_self_couplings( operator_PRECISION_struct *op, level_struct *l );// no effect
  void operator_PRECISION_set_neighbor_couplings( operator_PRECISION_struct *op, level_struct *l );// no effect
  
  void operator_PRECISION_test_routine( operator_PRECISION_struct *op, level_struct *l, struct Thread *threading );
  
#endif
