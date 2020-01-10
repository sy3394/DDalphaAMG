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
 * copied: 11/30/2019
 * changed from sbacchio
 * checked: 12/08/2019
 * 1st cleanup:12/19/2019
 */

#ifndef DATA_PRECISION_HEADER
  #define DATA_PRECISION_HEADER
  
  void buffer_PRECISION_define( complex_PRECISION *phi, complex_PRECISION value, int start, int end, level_struct *l );
  void buffer_PRECISION_copy( complex_PRECISION *z, complex_PRECISION *x, int start, int end, level_struct *l ); // z := x

#endif
