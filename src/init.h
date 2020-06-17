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
 * minor change from sbacchio
 * 1st cleanup:12/18/2019
 */

#ifndef INIT_HEADER
  #define INIT_HEADER

  struct init;

  void l_init( level_struct *l );
  void g_init( ); //level_struct *l );//!!!!!!!!!!

  void method_init( int *argc, char ***argv, level_struct *l );
  void method_setup( vector_double *V, level_struct *l, struct Thread *threading );
  void method_re_setup( level_struct *l, struct Thread *threading );
  void next_level_setup( vector_double *V, level_struct *l, struct Thread *threading );
  void method_iterative_setup( int setup_iter, level_struct *l, struct Thread *threading );
  void method_free( level_struct *l );
  void method_finalize( level_struct *l );
  
  void next_level_setup( vector_double *V, level_struct *l, struct Thread *threading );
  void next_level_free( level_struct *l );
  
  
#endif
