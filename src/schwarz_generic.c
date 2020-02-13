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
 * checked: 11/29/2019
 * changed from sbacchio
 * checked:12/06/2019: some work remains to be done, also ordering and static
 * 1st cleanup:12/18/2019
 */

#include "main.h"

void smoother_PRECISION_def( level_struct *l ) {
  
  if ( g.method >= 0 )
    schwarz_PRECISION_def( &(l->s_PRECISION), &(g.op_double), l );// Copies the Dirac operator and the clover term from g.op_double into the Schwarz s
  
  l->p_PRECISION.op = &(l->s_PRECISION.op);
  l->p_PRECISION.v_start = 0;
  l->p_PRECISION.v_end = l->inner_vector_size;
  l->p_PRECISION.eval_operator = (l->depth > 0)?apply_coarse_operator_PRECISION_new:d_plus_clover_PRECISION_new;
}

void schwarz_PRECISION_def( schwarz_PRECISION_struct *s, operator_double_struct *op, level_struct *l ) {

  schwarz_PRECISION_alloc( s, l );
  schwarz_layout_PRECISION_define( s, l );
  schwarz_PRECISION_setup( s, op, l );
}

void smoother_PRECISION_free( level_struct *l ) {
  
  if ( g.method >= 0 )
    schwarz_PRECISION_free( &(l->s_PRECISION), l );
}

// moved from vcycle
void smoother_PRECISION_new( vector_PRECISION *phi, vector_PRECISION *Dphi, vector_PRECISION *eta,
			     int n, const int res, level_struct *l, struct Thread *threading ) {

  ASSERT( phi->vector_buffer != eta->vector_buffer );

  START_MASTER(threading);
  PROF_PRECISION_START( _SM );
  END_MASTER(threading);

  if ( g.method == 1 ) {
    additive_schwarz_PRECISION_new( phi, Dphi, eta, n, res, &(l->s_PRECISION), l, threading );
  } else if ( g.method == 2 ) {
    red_black_schwarz_PRECISION_new( phi, Dphi, eta, n, res, &(l->s_PRECISION), l, threading );
  } else if ( g.method == 3 ) {
    sixteen_color_schwarz_PRECISION_new( phi, Dphi, eta, n, res, &(l->s_PRECISION), l, threading );
  } else {
    // commented out in milla???????
  }
  START_MASTER(threading);
  PROF_PRECISION_STOP( _SM, n );
  END_MASTER(threading);
}

/**************  SCHWARZ  *******************/

void schwarz_PRECISION_init( schwarz_PRECISION_struct *s, level_struct *l ) {

  int i;

  operator_PRECISION_init( &(s->op) );
  
  s->index[T]    = NULL;
  s->oe_index[T] = NULL;
  s->block       = NULL;
  for( i=0; i<5; i++ )
    vector_PRECISION_init(&(s->buf[i]));
  for ( i=0; i<2; i++ )
    vector_PRECISION_init(&(l->sbuf_PRECISION[i]));
  for( i=0; i<4; i++ )
    vector_PRECISION_init(&(s->oe_buf[i]));
  for( i=0; i<3; i++ )
    s->local_minres_buffer[i] = NULL;
  s->block_list        = NULL;
  s->block_list_length = NULL;
  s->num_colors        = 0;
}

void schwarz_PRECISION_alloc( schwarz_PRECISION_struct *s, level_struct *l ) {
  
  int i, j, n, mu, nu, *bl = l->block_lattice;
  
  //why are they here??????
  if ( g.method == 4 ) {//FGMRES + GMRES
    fgmres_PRECISION_struct_alloc( l->block_iter, 1, (l->depth==0)?_INNER:_ORDINARY,
                                   EPS_PRECISION, _COARSE_GMRES, _NOTHING, NULL,
                                   (l->depth==0)?(g.odd_even?apply_schur_complement_PRECISION_new:d_plus_clover_PRECISION_new):
                                   (g.odd_even?coarse_apply_schur_complement_PRECISION_new:apply_coarse_operator_PRECISION_new),
                                   &(l->sp_PRECISION), l );
  } else if ( g.method == 5 ) {//FGMRES + biCGstab (no AMG)
    fgmres_PRECISION_struct_alloc( 5, 1, (l->depth==0)?_INNER:_ORDINARY,
                                   EPS_PRECISION, _COARSE_GMRES, _NOTHING, NULL,
                                   (l->depth==0)?(g.odd_even?apply_schur_complement_PRECISION_new:d_plus_clover_PRECISION_new):
                                   (g.odd_even?coarse_apply_schur_complement_PRECISION_new:apply_coarse_operator_PRECISION_new),
                                   &(l->sp_PRECISION), l );
  }

  //--------- allocates memory for setting up an operator and compute neighbor tables. 
  operator_PRECISION_alloc( &(s->op), _SCHWARZ, l );
  if ( l->level > 0 && l->depth > 0 ) l->p_PRECISION.op = &(s->op);//this is hard to find!!!!!!!!

  //--------- compute the size of index tables and allocate memory for them
  // initialize dir_length (#block lattice sites excluding inner boundary sites in mu dir)
 for( mu=0; mu<4; mu++ ) 
    s->dir_length[mu] = bl[T]*bl[Z]*bl[Y]*bl[X]/bl[mu]*(bl[mu]-1);
  // allocate memory for index[mu] accordingly  
  MALLOC( s->index[T], int, MAX(1,s->dir_length[T]+s->dir_length[Z]+s->dir_length[Y]+s->dir_length[X]) );
  for( mu=1; mu<4; mu++ )
    s->index[mu] = s->index[mu-1]+s->dir_length[mu-1];

  if ( l->depth == 0 && g.odd_even ) {
    // allocate memory for oe_index[mu]
    MALLOC( s->oe_index[T], int, MAX(1,s->dir_length[T]+s->dir_length[Z]+s->dir_length[Y]+s->dir_length[X]) );
    for( mu=1; mu<4; mu++ )
      s->oe_index[mu] = s->oe_index[mu-1]+s->dir_length[mu-1];
  }
  
  //------------ compute info associated with blocking and allocate memory for block struc
  // compute #blcoks and #sites in the block
  s->num_blocks = 1;
  s->num_block_sites = 1;
  for ( mu=0; mu<4; mu++ ) {
    s->num_block_sites *= bl[mu];
    s->num_blocks *= l->local_lattice[mu]/bl[mu];
  }
  s->block_vector_size = s->num_block_sites*l->num_lattice_site_var;
  
  if ( g.method == 3 ) {
    MALLOC( s->block_list, int*, 16 );
    s->block_list[0] = NULL;
    MALLOC( s->block_list[0], int, s->num_blocks );
    j = s->num_blocks/16;
    for ( i=1; i<16; i++ )
      s->block_list[i] = s->block_list[0]+i*j;
  } else if ( g.method == 2 ) {
    MALLOC( s->block_list_length, int, 8 );
    MALLOC( s->block_list, int*, 8 );
    for ( i=0; i<8; i++ ) {
      s->block_list[i] = NULL;
      MALLOC( s->block_list[i], int, s->num_blocks );
    }
  }
  
  // for each block, allocate space for block_struct
  MALLOC( s->block, block_struct, s->num_blocks );
  
  /* s->block_boundary_length[2*mu]   : starting index of inner and outer boundary sites in pos mu dir
     s->block_boundary_length[2*mu+1] : starting index of inner and outer boundary sites in neg mu dir */
  n = 0; 
  for ( mu=0; mu<4; mu++ ) {
    i = bl[T]*bl[Z]*bl[Y]*bl[X]/bl[mu];
    s->block_boundary_length[2*mu] = n;       // plus mu dir
    s->block_boundary_length[2*mu+1] = n+2*i; // minus mu dir: 2*2(extra 2 for neighber)
    n += 4*i;
  }
  s->block_boundary_length[8] = n;            // = (4:#dirs)x[(#inner boundary sites in mu dir)+(#ext. boundary sites in mu dir)]x(2:pos&neg)
  
  for ( i=0; i<s->num_blocks; i++ ) {
    s->block[i].bt = NULL;
    MALLOC( s->block[i].bt, int, n );
  }

  //---------- allocate memory for vector buffers
  int nvec = num_loop;//(g.num_rhs_vect<l->num_eig_vect)?l->num_eig_vect:g.num_rhs_vect;//!!!!!!!!
  int svs  = l->schwarz_vector_size;

#ifdef HAVE_TM1p1
  svs *= 2;
  nvec *= 2;
#endif

  if ( l->depth == 0 )
    for ( i=0; i<4; i++ ) 
      vector_PRECISION_alloc( &(s->oe_buf[i]), _INNER, nvec, l, no_threading );
  
  for ( i=0; i<4; i++ ) 
    vector_PRECISION_alloc( &(s->buf[i]), (i==0)?((l->depth==0)?_INNER:_ORDINARY):_SCHWARZ, nvec, l, no_threading );

  if ( g.method == 1 )
    vector_PRECISION_alloc( &(s->buf[4]), _SCHWARZ, nvec, l, no_threading );

  for ( i=0; i<2; i++ )
    vector_PRECISION_alloc( &(l->sbuf_PRECISION[i]), (l->depth==0)?_INNER:_ORDINARY, nvec, l, no_threading );
  
  // these buffers are introduced to make local_minres_PRECISION thread-safe
  for ( i=0; i<3; i++)
    MALLOC( s->local_minres_buffer[i], complex_PRECISION, svs*nvec );
}


void schwarz_PRECISION_free( schwarz_PRECISION_struct *s, level_struct *l ) {
  
  int i, n, mu, nu, *bl = l->block_lattice;
  
  if ( g.method == 4 || g.method == 5 )
    fgmres_PRECISION_struct_free( &(l->sp_PRECISION), l );

  //------------ free operator
  operator_PRECISION_free( &(s->op), _SCHWARZ, l );
  
  //------------ free memory for index tables
  FREE( s->index[T], int, MAX(1,s->dir_length[T]+s->dir_length[Z]+s->dir_length[Y]+s->dir_length[X]) );
  for ( mu=1; mu<4; mu++)
    s->index[mu] = NULL;
  if ( l->depth == 0 && g.odd_even ) {
    FREE( s->oe_index[T], int, MAX(1,s->dir_length[T]+s->dir_length[Z]+s->dir_length[Y]+s->dir_length[X]) );
    for (mu=1; mu<4; mu++)
      s->oe_index[mu] = NULL;
  }
  
  //---------------- free memory for block struc
  for ( i=0; i<s->num_blocks; i++ )
    FREE( s->block[i].bt, int, s->block_boundary_length[8] );
  if ( g.method == 3 ) {
    FREE( s->block_list[0], int, s->num_blocks );
    FREE( s->block_list, int*, 16 );
  } else if ( g.method == 2 ) {
    FREE( s->block_list_length, int, 8 );
    for ( i=0; i<8; i++ )
      FREE( s->block_list[i], int, s->num_blocks );
    FREE( s->block_list, int*, 8 );
  }
  
  FREE( s->block, block_struct, s->num_blocks );

  //------------- free buffer memories
  int nvec = s->buf[0].num_vect;
  int svs = l->schwarz_vector_size;

#ifdef HAVE_TM1p1
  svs *= 2;
#endif
  if ( l->depth == 0 )
    for ( i=0; i<4; i++ )
      vector_PRECISION_free( &(s->oe_buf[i]), l, no_threading );
  
  for ( i=0; i<4; i++ )
    vector_PRECISION_free( &(s->buf[i]), l, no_threading );

  if ( g.method == 1 )
    vector_PRECISION_free( &(s->buf[4]), l, no_threading );
  
  for ( i=0; i<2; i++ )
    vector_PRECISION_free( &(l->sbuf_PRECISION[i]), l, no_threading );

  for (i=0; i<3; i++) {
    FREE( s->local_minres_buffer[i], complex_PRECISION, svs*nvec );
    s->local_minres_buffer[i] = NULL;
  }
  
}

void schwarz_layout_PRECISION_define( schwarz_PRECISION_struct *s, level_struct *l ) {

  int a0, b0, c0, d0, a1, b1, c1, d1, block_split[4], block_size[4], agg_split[4], 
      i, j, k, mu, index, x, y, z, t, ls[4], le[4], l_st[4], l_en[4], *dt = s->op.table_dim,
      *dt_mod = s->op.table_mod_dim, *it = s->op.index_table, *count[4];

  // Define coloring    
  if ( g.method == 1 )        // Additive
    s->num_colors = 1;
  else if ( g.method == 2 )   // Red-Black
    s->num_colors = 2;
  else if ( g.method == 3 ) { // 16 Color
    int flag = 0;
    for ( mu=0; mu<4; mu++ ) {
      if ( (l->local_lattice[mu]/l->block_lattice[mu]) % 2 == 1 )
        flag = 1;
    }
    if ( flag == 0 )
      s->num_colors = 16;
    else {
      s->num_colors = 2;
      printf0("depth: %d, switching to red black schwarz as smoother\n", l->depth );
    }
  }
  
  const int sigma[16] = {0,1,3,2,6,4,5,7,15,14,12,13,9,11,10,8};
  const int color_to_comm[16][2] = { {T,-1}, {X,+1}, {Y,+1}, {X,-1}, {Z,+1}, {Y,-1}, {X,+1}, {Y,+1},
                               {T,+1}, {X,-1}, {Y,-1}, {X,+1}, {Z,-1}, {Y,+1}, {X,-1}, {Y,-1}  };
  int color_counter[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  
  s->num_block_sites = 1;
  s->block_oe_offset = 0;
  s->num_aggregates = 1;
  
  if ( g.method == 2 ) {
    for ( i=0; i<8; i++ )
      s->block_list_length[i]=0;
  }
  
  for ( mu=0; mu<4; mu++ ) {
    s->num_block_sites *= l->block_lattice[mu]; // #sites in a block, which is within an aggregate
    s->block_oe_offset += ((l->local_lattice[mu]/l->block_lattice[mu])*(g.my_coords[mu]/l->comm_offset[mu]))%2;//??????
    ls[mu] = 0;
    le[mu] = ls[mu] + l->local_lattice[mu];
    dt[mu] = l->local_lattice[mu]+2;            // dim of local lattice on a process with ghost cells
    dt_mod[mu] = l->local_lattice[mu]+2;        // dim of local lattice on a process with ghost cells
    l_st[mu] = ls[mu];
    l_en[mu] = le[mu];
    agg_split[mu] = l->local_lattice[mu]/l->coarsening[mu];   // # aggregates in mu dir
    block_split[mu] = l->coarsening[mu]/l->block_lattice[mu]; // # blocks within each aggregate
    block_size[mu] = l->block_lattice[mu];                    // dims of a block on a local lattice
    s->num_aggregates *= agg_split[mu];                       // # aggregates in a local lattice
  }
  s->block_oe_offset = s->block_oe_offset%2;//%2 is already applied????
  s->block_vector_size = s->num_block_sites*l->num_lattice_site_var;

  // Schwarz indexing: first counts blocks within an aggregate and then those in another aggregate both in lexicographical order: x->...->t
  i = 0; j = 0; // i counts lattice index; j counts blocks
  // inner hyper cuboid
  count[T] = &d1; count[Z] = &c1; count[Y] = &b1; count[X] = &a1;
  // for each aggregate
  for ( d0=0; d0<agg_split[T]; d0++ )
    for ( c0=0; c0<agg_split[Z]; c0++ )
      for ( b0=0; b0<agg_split[Y]; b0++ )
        for ( a0=0; a0<agg_split[X]; a0++ ) {
          // for each block in the aggregate
          for ( d1=d0*block_split[T]; d1<(d0+1)*block_split[T]; d1++ )
            for ( c1=c0*block_split[Z]; c1<(c0+1)*block_split[Z]; c1++ )
              for ( b1=b0*block_split[Y]; b1<(b0+1)*block_split[Y]; b1++ )
                for ( a1=a0*block_split[X]; a1<(a0+1)*block_split[X]; a1++ ) {
                  
                  s->block[j].start = i; // iteration index for the first site in the j^th block
                  s->block[j].no_comm = 1;//????
                  if ( s->num_colors == 1 ) {
                    s->block[j].color = 0;
                  } else if ( s->num_colors == 2 ) {
                    s->block[j].color = ( d1+c1+b1+a1+s->block_oe_offset )%2;//each block is either red or black
                  } else if ( s->num_colors == 16 ) {
                    for ( k=0; k<16; k++ )
                      if ( sigma[k] == 8*(d1%2)+4*(c1%2)+2*(b1%2)+1*(a1%2) ) {
                        s->block[j].color = k;
                        s->block_list[k][color_counter[k]] = j;
                        color_counter[k]++;
                        break;
                      }
                  }
                  
                  if ( s->num_colors == 1 || s->num_colors == 2 ) {
                    for ( mu=0; mu<4; mu++ ) {
                      if ( ( (*count[mu]) == 0 ) || ( (*count[mu]+1) == le[mu]/block_size[mu] ) )
                        s->block[j].no_comm = 0;//????????
                    }
                    
                    if ( s->num_colors == 2 ) {
                      // calculate boundary correspondence of the block
                      int count_plus=0, count_minus=0, count_inner=0, index;
                      for ( mu=0; mu<4; mu++ ) {
                        if ( (*count[mu]) == 0 )
                          count_minus++;
                        if ( (*count[mu]+1) == le[mu]/block_size[mu] )
                          count_plus++;
                        if ( (*count[mu]) != 0 && (*count[mu]+1) != le[mu]/block_size[mu] )
                          count_inner++;
                      }
                      
                      if ( count_inner == 4 ) {
                        index = 4*s->block[j].color;
                      } else if ( count_minus == 0 ) {
                        if ( s->block[j].color == 0 ) index = 1;
                        else index = 7;
                      } else if ( count_plus == 0 ) {
                        if ( s->block[j].color == 0 ) index = 3;
                        else index = 5;
                      } else {
                        index = 2 + 4*s->block[j].color;
                      }
                      
                      s->block_list[index][s->block_list_length[index]] = j;
                      s->block_list_length[index]++;
                    }
                    
                  } else if ( s->num_colors == 16 ) {
                    k = s->block[j].color;
                    if ( k == 0 ) {
                      for ( mu=0; mu<4; mu++ ) {
                        if ( (*count[mu]) == 0 )
                          s->block[j].no_comm = 0;
                      }
                    } else {
                      mu = color_to_comm[k][0];
                      if ( (color_to_comm[k][1] == +1 && (*count[mu]+1) == le[mu]/block_size[mu]) ||
                           (color_to_comm[k][1] == -1 && (*count[mu]) == 0 ) )
                        s->block[j].no_comm = 0;
                    }
                  }
                  
                  j++;//End of encoding the block information
                  
                  // --- set up index table (it=s->op.index_table: local_lex_index -> Schwarz_index)
		  // Here, a local lattice with its boundaries (including corners)
		  // are indexed in x->...->t(fast->slow) within a block
		  // The result is stored in it[lex_index]
		  // Start with inner sites
                  if ( l->depth == 0 && g.odd_even ) {
                    // odd even on the blocks: even sites first, then odd sites
                    // even sites
                    for ( t=d1*block_size[T]; t<(d1+1)*block_size[T]; t++ )
                      for ( z=c1*block_size[Z]; z<(c1+1)*block_size[Z]; z++ )
                        for ( y=b1*block_size[Y]; y<(b1+1)*block_size[Y]; y++ )
                          for ( x=a1*block_size[X]; x<(a1+1)*block_size[X]; x++ ) {
			    // if the given site is even
                            if (((t-d1*block_size[T])+(z-c1*block_size[Z])+
                                 (y-b1*block_size[Y])+(x-a1*block_size[X]))%2 == 0 ) {
                              index = lex_index( t, z, y, x, dt );
                              it[index] = i;
                              i++;
                            }
                          }
                    // odd sites
                    for ( t=d1*block_size[T]; t<(d1+1)*block_size[T]; t++ )
                      for ( z=c1*block_size[Z]; z<(c1+1)*block_size[Z]; z++ )
                        for ( y=b1*block_size[Y]; y<(b1+1)*block_size[Y]; y++ )
                          for ( x=a1*block_size[X]; x<(a1+1)*block_size[X]; x++ ) {
			    // if the given site is odd
                            if (((t-d1*block_size[T])+(z-c1*block_size[Z])+
                                  (y-b1*block_size[Y])+(x-a1*block_size[X]))%2 == 1 ) {
                              index = lex_index( t, z, y, x, dt );
                              it[index] = i;
                              i++;
                            }
                          }
                  } else {
                    // no odd even
                    for ( t=d1*block_size[T]; t<(d1+1)*block_size[T]; t++ )
                      for ( z=c1*block_size[Z]; z<(c1+1)*block_size[Z]; z++ )
                        for ( y=b1*block_size[Y]; y<(b1+1)*block_size[Y]; y++ )
                          for ( x=a1*block_size[X]; x<(a1+1)*block_size[X]; x++ ) {
                            index = lex_index( t, z, y, x, dt );
                            it[index] = i;
                            i++;
                          }
                  }
                }
        }// END: loop over blocks and aggregates
        
  

  // Then, index postive ghost shell (exterior boundary sites of the local lattice): +t(fastest)->...->+t(slowest)
  for ( mu=0; mu<4; mu++ ) {
    l_st[mu] = le[mu];
    l_en[mu] = le[mu]+1;
    for ( t=l_st[T]; t<l_en[T]; t++ )
      for ( z=l_st[Z]; z<l_en[Z]; z++ )
        for ( y=l_st[Y]; y<l_en[Y]; y++ )
          for ( x=l_st[X]; x<l_en[X]; x++ ) {
            index = lex_index( t, z, y, x, dt );
            it[index] = i;
            i++;
          }  
    l_st[mu] = ls[mu];
    l_en[mu] = le[mu];
  }
  
  // Index negative ghost shell: -x->...->-t
  for ( mu=0; mu<4; mu++ ) {
    l_st[mu] = ls[mu]-1;
    l_en[mu] = ls[mu];
    for ( t=l_st[T]; t<l_en[T]; t++ )
      for ( z=l_st[Z]; z<l_en[Z]; z++ )
        for ( y=l_st[Y]; y<l_en[Y]; y++ )
          for ( x=l_st[X]; x<l_en[X]; x++ ) {
            index = lex_mod_index( t, z, y, x, dt );
            it[index] = i;
            i++;
          }  
    l_st[mu] = ls[mu];
    l_en[mu] = le[mu];
  }

  //--- block boundary table (block[j].bt: iter_index (boundary) -> Schwarz_index)
  i = 0; j = 0;// i:conts sites (pos/neg:pos->neg inner boundary and its neighbor in pos dir first, mu:T->X, sites:x->t); j:counts blocks
  // for each aggregate
  for ( d0=0; d0<agg_split[T]; d0++ )
    for ( c0=0; c0<agg_split[Z]; c0++ )
      for ( b0=0; b0<agg_split[Y]; b0++ )
        for ( a0=0; a0<agg_split[X]; a0++ )
	  // for each block in the aggregate
          for ( d1=d0*block_split[T]; d1<(d0+1)*block_split[T]; d1++ )
            for ( c1=c0*block_split[Z]; c1<(c0+1)*block_split[Z]; c1++ )
              for ( b1=b0*block_split[Y]; b1<(b0+1)*block_split[Y]; b1++ )
                for ( a1=a0*block_split[X]; a1<(a0+1)*block_split[X]; a1++ ) {

                  i=0;
                  int block_start[4], block_end[4], tmp;
                  
                  block_start[T] = d1*block_size[T]; block_start[Z] = c1*block_size[Z];
                  block_start[Y] = b1*block_size[Y]; block_start[X] = a1*block_size[X];
                  block_end[T] = (d1+1)*block_size[T]; block_end[Z] = (c1+1)*block_size[Z];
                  block_end[Y] = (b1+1)*block_size[Y]; block_end[X] = (a1+1)*block_size[X];
                  
                  for ( mu=0; mu<4; mu++ ) {
                    tmp = block_start[mu];
                                        
                    // positive inner boundary in mu dir: neighbor in plus dir
                    block_start[mu] = block_end[mu]-1;
                    for ( t=block_start[T]; t<block_end[T]; t++ )
                      for ( z=block_start[Z]; z<block_end[Z]; z++ )
                        for ( y=block_start[Y]; y<block_end[Y]; y++ )
                          for ( x=block_start[X]; x<block_end[X]; x++ ) {
                            s->block[j].bt[i] = site_index( t, z, y, x, dt, it );
                            i++;
                            s->block[j].bt[i] = connect_link_PRECISION( t, z, y, x, mu, +1, dt, it, s, l );
                            i++;
                          }
                          
                    block_start[mu] = tmp;
                    tmp = block_end[mu];
                        
                    // negative inner boundary in mu dir: neighbor in minus dir
                    block_end[mu] = block_start[mu]+1;
                    for ( t=block_start[T]; t<block_end[T]; t++ )
                      for ( z=block_start[Z]; z<block_end[Z]; z++ )
                        for ( y=block_start[Y]; y<block_end[Y]; y++ )
                          for ( x=block_start[X]; x<block_end[X]; x++ ) {
                            s->block[j].bt[i] = site_index( t, z, y, x, dt, it );
                            i++;
                            s->block[j].bt[i] = connect_link_PRECISION( t, z, y, x, mu, -1, dt, it, s, l );
                            i++;
                          }
                          block_end[mu] = tmp;
                  }
                  j++;
                }
  
  //--- index tables for block dirac operator
  if ( l->depth == 0 && g.odd_even ) {
    i = 1;
    for ( mu=0; mu<4; mu++ ) i*= block_size[mu];
    s->num_block_even_sites = i/2 + i%2;
    s->num_block_odd_sites  = i/2;

    count[T] = &t; count[Z] = &z; count[Y] = &y; count[X] = &x;    
    for ( mu=0; mu<4; mu++ ) {
      // even sites, plus dir ( = odd sites, minus dir )
      i=0; j=0;//i:counts even part of a block; j:counts???
      for ( t=0; t<block_size[T]; t++ )
        for ( z=0; z<block_size[Z]; z++ )
          for ( y=0; y<block_size[Y]; y++ )
            for ( x=0; x<block_size[X]; x++ ) {
              if ( (t+z+y+x)%2 == 0 ) {
                if ( *(count[mu]) < block_size[mu]-1 ) {
                  s->oe_index[mu][j] = i; 
                  j++;
                }
                i++;
              }
            }
      s->dir_length_even[mu] = j;
      // odd sites, plus dir ( = even sites, minus dir )
      j=0;
      for ( t=0; t<block_size[T]; t++ )
        for ( z=0; z<block_size[Z]; z++ )
          for ( y=0; y<block_size[Y]; y++ )
            for ( x=0; x<block_size[X]; x++ ) {
              if ( (t+z+y+x)%2 == 1 ) {
                if ( *(count[mu]) < block_size[mu]-1 ) {//????????
                  s->oe_index[mu][s->dir_length_even[mu]+j] = i; 
                  j++;
                }
                i++;
              }
            }
      s->dir_length_odd[mu] = j;
    }
  }
  
  count[T] = &t; count[Z] = &z; count[Y] = &y; count[X] = &x;
  for ( mu=0; mu<4; mu++ ) {
    j=0;
    for ( t=0; t<block_size[T]; t++ )
      for ( z=0; z<block_size[Z]; z++ )
        for ( y=0; y<block_size[Y]; y++ )
          for ( x=0; x<block_size[X]; x++ ) {
            if ( *(count[mu]) < block_size[mu]-1 ) {//???????
              s->index[mu][j] =  site_index( t, z, y, x, dt, it ); j++;
            }
          }
  }
  
  //--- define neighbor table (for the application of the entire operator),
  //      s->op.neighbor_table: Schwarz index+mu -> Schwarz index of the neighbor in the pos mu dir
  //      s->op.backward_neighbor_table: Schwarz index+mu -> Schwarz index of the neighbor in the neg mu dir 
  //      s->op.c.boundary_table: iter_index of inner bondaray sites -> Schwarz index
  //      s->op.translation_table: local_lex_index -> Schwarz_index
  define_nt_bt_tt( s->op.neighbor_table, s->op.backward_neighbor_table, s->op.c.boundary_table, s->op.translation_table, it, dt, l );
}

void schwarz_PRECISION_setup( schwarz_PRECISION_struct *s, operator_double_struct *op_in, level_struct *l ) {

/*********************************************************************************  
* Copies the Dirac operator and the clover term from op_in into the Schwarz 
* struct (this function is depth 0 only).
* - operator_double_struct *op_in: Input operator.                                  
*********************************************************************************/

  int i, index, n = l->num_inner_lattice_sites, *tt = s->op.translation_table;
  config_PRECISION D_out_pt, clover_out_pt, odd_proj_out_pt;
  config_double D_in_pt = op_in->D, clover_in_pt = op_in->clover, odd_proj_in_pt = op_in->odd_proj;
  
  s->op.m0 = op_in->m0;

  for ( i=0; i<n; i++ ) {
    index = tt[i];
    D_out_pt = s->op.D + 36*index;
    FOR36( *D_out_pt = (complex_PRECISION) *D_in_pt; D_out_pt++; D_in_pt++; );
  }
  
  if ( g.csw != 0 ) {
    for ( i=0; i<n; i++ ) {
      index = tt[i];
      clover_out_pt = s->op.clover + 42*index;
      FOR42( *clover_out_pt = (complex_PRECISION) *clover_in_pt; clover_out_pt++; clover_in_pt++; );
    }
  } else {
    for ( i=0; i<n; i++ ) {
      index = tt[i];
      clover_out_pt = s->op.clover + 12*index;
      FOR12( *clover_out_pt = (complex_PRECISION) *clover_in_pt; clover_out_pt++; clover_in_pt++; );
    }
  }

  for ( i=0; i<n; i++ ) {
    index = tt[i];
    odd_proj_out_pt = s->op.odd_proj + 12*index;
    FOR12( *odd_proj_out_pt = (complex_PRECISION) *odd_proj_in_pt; odd_proj_out_pt++; odd_proj_in_pt++; );
  }

#ifdef HAVE_TM
  tm_term_PRECISION_setup( (PRECISION) (g.mu_factor[l->depth]*op_in->mu), (PRECISION) (g.mu_factor[l->depth]*op_in->mu_even_shift), (PRECISION) (g.mu_factor[l->depth]*op_in->mu_odd_shift), &(s->op), l, no_threading );
#endif  

#ifdef HAVE_TM1p1
  epsbar_term_PRECISION_setup( (PRECISION) (g.epsbar_factor[l->depth]*op_in->epsbar), (PRECISION) (g.epsbar_factor[l->depth]*op_in->epsbar_ig5_even_shift), (PRECISION) (g.epsbar_factor[l->depth]*op_in->epsbar_ig5_odd_shift), &(s->op), l, no_threading );
#endif

  schwarz_PRECISION_boundary_update( s, l );
  
  if ( g.odd_even )
    schwarz_PRECISION_oddeven_setup( s, l );

}
//intact
void schwarz_PRECISION_boundary_update( schwarz_PRECISION_struct *s, level_struct *l ) {

/*********************************************************************************
* Updates the current level hopping term in "s->op.D" on the process boundaries
* in all negative directions. This is necessary for enabling Schwarz to perform
* local block residual updates on demand.
* Updates hopping terms at the ghost sites on the negative boundaries
*********************************************************************************/   

  int i, t, z, y, x, mu, nu, index, *it = s->op.index_table, *dt = s->op.table_dim,
      ls[4], le[4], buf_length[4], link_size;
  buffer_PRECISION buf[4] = {NULL,NULL,NULL,NULL}, rbuf[4] = {NULL,NULL,NULL,NULL};
  config_PRECISION D=s->op.D;
  
  for ( mu=0; mu<4; mu++ ) {
    ls[mu] = 0;
    le[mu] = l->local_lattice[mu];
    buf_length[mu] = 0;
  }
  
  if ( l->depth == 0 )
    link_size = 4*9; // size of hopping term in D at each site (spin x color^2)?????there are only four non-zero spinor entries???
  else
    link_size = 4*SQUARE(l->num_lattice_site_var);
  
  // allocate buffers
  //if ( l->global_splitting[mu] > 1 ) {//my proposal!!!!!!!!!
  for ( mu=0; mu<4; mu++ ) {
    if ( l->global_splitting[mu] > 1 ) {
      buf_length[mu] = link_size;
      for ( nu=0; nu<4; nu++ ) {
        if ( nu != mu )
          buf_length[mu] *= le[nu];
      }
      MALLOC( buf[mu], complex_PRECISION, buf_length[mu] );
      MALLOC( rbuf[mu], complex_PRECISION, buf_length[mu] );
    }
  }
  
  // post recv for desired directions
  for ( mu=0; mu<4; mu++ ) {
    if ( l->global_splitting[mu] > 1 ) {
      MPI_Irecv( rbuf[mu], buf_length[mu], MPI_COMPLEX_PRECISION, l->neighbor_rank[2*mu+1],
                 2*mu+1, g.comm_cart, &(s->op.c.rreqs[2*mu+1]) );
    }
  }
  
  // buffer data for send and send it
  for ( mu=0; mu<4; mu++ ) {
    if ( l->global_splitting[mu] > 1 ) {
      ls[mu] = l->local_lattice[mu]-1;
      i=0;
      for ( t=ls[T]; t<le[T]; t++ )
        for ( z=ls[Z]; z<le[Z]; z++ )
          for ( y=ls[Y]; y<le[Y]; y++ )
            for ( x=ls[X]; x<le[X]; x++ ) {
              index = site_index( t, z, y, x, dt, it );
              buffer_PRECISION_copy( buf[mu]+i*link_size, D+index*link_size, 0, link_size, l );
              i++;
            }
      MPI_Isend( buf[mu], buf_length[mu], MPI_COMPLEX_PRECISION, l->neighbor_rank[2*mu],
                 2*mu+1, g.comm_cart, &(s->op.c.sreqs[2*mu+1]) );
      ls[mu] = 0;
    }
  }
  
  // store links in desired ordering after recv
  for ( mu=0; mu<4; mu++ ) {
    if ( l->global_splitting[mu] > 1 ) {
      MPI_Wait( &(s->op.c.rreqs[2*mu+1]), MPI_STATUS_IGNORE );
      ls[mu] = -1;
      le[mu] = 0;
      i=0;
      for ( t=ls[T]; t<le[T]; t++ )
        for ( z=ls[Z]; z<le[Z]; z++ )
          for ( y=ls[Y]; y<le[Y]; y++ )
            for ( x=ls[X]; x<le[X]; x++ ) {
              index = site_mod_index( t, z, y, x, dt, it );
              buffer_PRECISION_copy( D+index*link_size, rbuf[mu]+i*link_size, 0, link_size, l );
              i++;
            }
      ls[mu] = 0;
      le[mu] = l->local_lattice[mu];
    }
  }
  
  // free buffers
  for ( mu=0; mu<4; mu++ ) {
    if ( l->global_splitting[mu] > 1 ) {
      MPI_Wait( &(s->op.c.sreqs[2*mu+1]), MPI_STATUS_IGNORE );
      FREE( buf[mu], complex_PRECISION, buf_length[mu] );
      FREE( rbuf[mu], complex_PRECISION, buf_length[mu] );
    }
  }
}

// eta <- block_PRECISION_boundary_op*phi
void block_PRECISION_boundary_op_new( vector_PRECISION *eta, vector_PRECISION *phi, int k,
                                  schwarz_PRECISION_struct *s, level_struct *l ) {
  // k: number of current block
  int *bbl = s->block_boundary_length;
  int nvec = phi->num_vect_now, nvec_phi = phi->num_vect, nvec_eta = eta->num_vect;//!!!!!!!
  int i, mu, index, neighbor_index;
  config_PRECISION D_pt, D = s->op.D;
  buffer_PRECISION phi_pt, eta_pt;

/*#ifdef HAVE_TM1p1
  if( g.n_flavours == 2 ) {
    complex_PRECISION buf1[24], *buf2=buf1+12;
    mu=T;
    // plus mu direction
    for ( i=bbl[2*mu]; i<bbl[2*mu+1]; i+=2 ) {
      index = s->block[k].bt[i];
      neighbor_index = s->block[k].bt[i+1];
      D_pt = D + 36*index + 9*mu;
      phi_pt = phi->vector_buffer + 24*neighbor_index;
      eta_pt = eta->vector_buffer + 24*index;
      dprp_T_PRECISION( buf1, phi_pt );
      mvm_PRECISION( buf2, D_pt, buf1 );
      mvm_PRECISION( buf2+3, D_pt, buf1+3 );
      mvm_PRECISION( buf2+6, D_pt, buf1+6 );
      mvm_PRECISION( buf2+9, D_pt, buf1+9 );
      dpbp_su3_T_PRECISION( buf2, eta_pt );
    }
    // minus mu direction
    for ( i=bbl[2*mu+1]; i<bbl[2*mu+2]; i+=2 ) {
      index = s->block[k].bt[i];
      neighbor_index = s->block[k].bt[i+1];
      D_pt = D + 36*neighbor_index + 9*mu;
      phi_pt = phi->vector_buffer + 24*neighbor_index;
      eta_pt = eta->vector_buffer + 24*index;
      dprn_T_PRECISION( buf1, phi_pt );
      mvmh_PRECISION( buf2, D_pt, buf1 );
      mvmh_PRECISION( buf2+3, D_pt, buf1+3 );
      mvmh_PRECISION( buf2+6, D_pt, buf1+6 );
      mvmh_PRECISION( buf2+9, D_pt, buf1+9 );
      dpbn_su3_T_PRECISION( buf2, eta_pt );
    }
    
    mu=Z;
    // plus mu direction
    for ( i=bbl[2*mu]; i<bbl[2*mu+1]; i+=2 ) {
      index = s->block[k].bt[i];
      neighbor_index = s->block[k].bt[i+1];
      D_pt = D + 36*index + 9*mu;
      phi_pt = phi->vector_buffer + 24*neighbor_index;
      eta_pt = eta->vector_buffer + 24*index;
      dprp_Z_PRECISION( buf1, phi_pt );
      mvm_PRECISION( buf2, D_pt, buf1 );
      mvm_PRECISION( buf2+3, D_pt, buf1+3 );
      mvm_PRECISION( buf2+6, D_pt, buf1+6 );
      mvm_PRECISION( buf2+9, D_pt, buf1+9 );
      dpbp_su3_Z_PRECISION( buf2, eta_pt );
    }
    // minus mu direction
    for ( i=bbl[2*mu+1]; i<bbl[2*mu+2]; i+=2 ) {
      index = s->block[k].bt[i];
      neighbor_index = s->block[k].bt[i+1];
      D_pt = D + 36*neighbor_index + 9*mu;
      phi_pt = phi->vector_buffer + 24*neighbor_index;
      eta_pt = eta->vector_buffer + 24*index;
      dprn_Z_PRECISION( buf1, phi_pt );
      mvmh_PRECISION( buf2, D_pt, buf1 );
      mvmh_PRECISION( buf2+3, D_pt, buf1+3 );
      mvmh_PRECISION( buf2+6, D_pt, buf1+6 );
      mvmh_PRECISION( buf2+9, D_pt, buf1+9 );
      dpbn_su3_Z_PRECISION( buf2, eta_pt );
    }
    
    mu=Y;
    // plus mu direction
    for ( i=bbl[2*mu]; i<bbl[2*mu+1]; i+=2 ) {
      index = s->block[k].bt[i];
      neighbor_index = s->block[k].bt[i+1];
      D_pt = D + 36*index + 9*mu;
      phi_pt = phi->vector_buffer + 24*neighbor_index;
      eta_pt = eta->vector_buffer + 24*index;
      dprp_Y_PRECISION( buf1, phi_pt );
      mvm_PRECISION( buf2, D_pt, buf1 );
      mvm_PRECISION( buf2+3, D_pt, buf1+3 );
      mvm_PRECISION( buf2+6, D_pt, buf1+6 );
      mvm_PRECISION( buf2+9, D_pt, buf1+9 );
      dpbp_su3_Y_PRECISION( buf2, eta_pt );
    }
    // minus mu direction
    for ( i=bbl[2*mu+1]; i<bbl[2*mu+2]; i+=2 ) {
      index = s->block[k].bt[i];
      neighbor_index = s->block[k].bt[i+1];
      D_pt = D + 36*neighbor_index + 9*mu;
      phi_pt = phi->vector_buffer + 24*neighbor_index;
      eta_pt = eta->vector_buffer + 24*index;
      dprn_Y_PRECISION( buf1, phi_pt );
      mvmh_PRECISION( buf2, D_pt, buf1 );
      mvmh_PRECISION( buf2+3, D_pt, buf1+3 );
      mvmh_PRECISION( buf2+6, D_pt, buf1+6 );
      mvmh_PRECISION( buf2+9, D_pt, buf1+9 );
      dpbn_su3_Y_PRECISION( buf2, eta_pt );
    }
  
    mu=X;
    // plus mu direction
    for ( i=bbl[2*mu]; i<bbl[2*mu+1]; i+=2 ) {
      index = s->block[k].bt[i];
      neighbor_index = s->block[k].bt[i+1];
      D_pt = D + 36*index + 9*mu;
      phi_pt = phi->vector_buffer + 24*neighbor_index;
      eta_pt = eta->vector_buffer + 24*index;
      dprp_X_PRECISION( buf1, phi_pt );
      mvm_PRECISION( buf2, D_pt, buf1 );
      mvm_PRECISION( buf2+3, D_pt, buf1+3 );
      mvm_PRECISION( buf2+6, D_pt, buf1+6 );
      mvm_PRECISION( buf2+9, D_pt, buf1+9 );
      dpbp_su3_X_PRECISION( buf2, eta_pt );
    }
    // minus mu direction
    for ( i=bbl[2*mu+1]; i<bbl[2*mu+2]; i+=2 ) {
      index = s->block[k].bt[i];
      neighbor_index = s->block[k].bt[i+1];
      D_pt = D + 36*neighbor_index + 9*mu;
      phi_pt = phi->vector_buffer + 24*neighbor_index;
      eta_pt = eta->vector_buffer + 24*index;
      dprn_X_PRECISION( buf1, phi_pt );
      mvmh_PRECISION( buf2, D_pt, buf1 );
      mvmh_PRECISION( buf2+3, D_pt, buf1+3 );
      mvmh_PRECISION( buf2+6, D_pt, buf1+6 );
      mvmh_PRECISION( buf2+9, D_pt, buf1+9 );
      dpbn_su3_X_PRECISION( buf2, eta_pt );
    }  
  } else {
#endif*/
    complex_PRECISION buf1[12*nvec], *buf2=buf1+6*nvec;
    mu=T;
    // plus mu direction
    for ( i=bbl[2*mu]; i<bbl[2*mu+1]; i+=2 ) {
      index = s->block[k].bt[i];
      neighbor_index = s->block[k].bt[i+1];
      D_pt = D + 36*index + 9*mu;
      phi_pt = phi->vector_buffer + 12*neighbor_index*nvec_phi;
      eta_pt = eta->vector_buffer + 12*index*nvec_eta;
      prp_T_PRECISION_new( buf1, phi_pt, nvec, nvec, nvec_phi);
      mvm_PRECISION_new( buf2, D_pt, buf1, nvec, nvec, nvec );
      mvm_PRECISION_new( buf2+3*nvec, D_pt, buf1+3*nvec, nvec, nvec, nvec );
      pbp_su3_T_PRECISION_new( buf2, eta_pt, nvec, nvec, nvec_eta);
    }
    // minus mu direction
    for ( i=bbl[2*mu+1]; i<bbl[2*mu+2]; i+=2 ) {
      index = s->block[k].bt[i];
      neighbor_index = s->block[k].bt[i+1];
      D_pt = D + 36*neighbor_index + 9*mu;
      phi_pt = phi->vector_buffer + 12*neighbor_index*nvec_phi;
      eta_pt = eta->vector_buffer + 12*index*nvec_eta;
      prn_T_PRECISION_new( buf1, phi_pt, nvec, nvec, nvec_phi );
      mvmh_PRECISION_new( buf2, D_pt, buf1 , nvec, nvec, nvec);
      mvmh_PRECISION_new( buf2+3*nvec, D_pt, buf1+3*nvec, nvec, nvec, nvec );
      pbn_su3_T_PRECISION_new( buf2, eta_pt, nvec, nvec, nvec_eta );
    }
    
    mu=Z;
    // plus mu direction
    for ( i=bbl[2*mu]; i<bbl[2*mu+1]; i+=2 ) {
      index = s->block[k].bt[i];
      neighbor_index = s->block[k].bt[i+1];
      D_pt = D + 36*index + 9*mu;
      phi_pt = phi->vector_buffer + 12*neighbor_index*nvec_phi;
      eta_pt = eta->vector_buffer + 12*index*nvec_eta;
      prp_Z_PRECISION_new( buf1, phi_pt, nvec, nvec, nvec_phi );
      mvm_PRECISION_new( buf2, D_pt, buf1, nvec, nvec, nvec );
      mvm_PRECISION_new( buf2+3*nvec, D_pt, buf1+3*nvec, nvec, nvec, nvec );
      pbp_su3_Z_PRECISION_new( buf2, eta_pt, nvec, nvec, nvec_eta );
    }
    // minus mu direction
    for ( i=bbl[2*mu+1]; i<bbl[2*mu+2]; i+=2 ) {
      index = s->block[k].bt[i];
      neighbor_index = s->block[k].bt[i+1];
      D_pt = D + 36*neighbor_index + 9*mu;
      phi_pt = phi->vector_buffer + 12*neighbor_index*nvec_phi;
      eta_pt = eta->vector_buffer + 12*index*nvec_eta;
      prn_Z_PRECISION_new( buf1, phi_pt, nvec, nvec, nvec_phi );
      mvmh_PRECISION_new( buf2, D_pt, buf1, nvec, nvec, nvec );
      mvmh_PRECISION_new( buf2+3*nvec, D_pt, buf1+3*nvec, nvec, nvec, nvec );
      pbn_su3_Z_PRECISION_new( buf2, eta_pt, nvec, nvec, nvec_eta );
    }
    
    mu=Y;
    // plus mu direction
    for ( i=bbl[2*mu]; i<bbl[2*mu+1]; i+=2 ) {
      index = s->block[k].bt[i];
      neighbor_index = s->block[k].bt[i+1];
      D_pt = D + 36*index + 9*mu;
      phi_pt = phi->vector_buffer + 12*neighbor_index*nvec_phi;
      eta_pt = eta->vector_buffer + 12*index*nvec_eta;
      prp_Y_PRECISION_new( buf1, phi_pt, nvec, nvec, nvec_phi );
      mvm_PRECISION_new( buf2, D_pt, buf1, nvec, nvec, nvec );
      mvm_PRECISION_new( buf2+3*nvec, D_pt, buf1+3*nvec, nvec, nvec, nvec );
      pbp_su3_Y_PRECISION_new( buf2, eta_pt, nvec, nvec, nvec_eta );
    }
    // minus mu direction
    for ( i=bbl[2*mu+1]; i<bbl[2*mu+2]; i+=2 ) {
      index = s->block[k].bt[i];
      neighbor_index = s->block[k].bt[i+1];
      D_pt = D + 36*neighbor_index + 9*mu;
      phi_pt = phi->vector_buffer + 12*neighbor_index*nvec_phi;
      eta_pt = eta->vector_buffer + 12*index*nvec_eta;
      prn_Y_PRECISION_new( buf1, phi_pt, nvec, nvec, nvec_phi );
      mvmh_PRECISION_new( buf2, D_pt, buf1, nvec, nvec, nvec );
      mvmh_PRECISION_new( buf2+3*nvec, D_pt, buf1+3*nvec, nvec, nvec, nvec );
      pbn_su3_Y_PRECISION_new( buf2, eta_pt, nvec, nvec, nvec_eta );
    }
  
    mu=X;
    // plus mu direction
    for ( i=bbl[2*mu]; i<bbl[2*mu+1]; i+=2 ) {
      index = s->block[k].bt[i];
      neighbor_index = s->block[k].bt[i+1];
      D_pt = D + 36*index + 9*mu;
      phi_pt = phi->vector_buffer + 12*neighbor_index*nvec_phi;
      eta_pt = eta->vector_buffer + 12*index*nvec_eta;
      prp_X_PRECISION_new( buf1, phi_pt, nvec, nvec, nvec_phi );
      mvm_PRECISION_new( buf2, D_pt, buf1, nvec, nvec, nvec );
      mvm_PRECISION_new( buf2+3*nvec, D_pt, buf1+3*nvec, nvec, nvec, nvec );
      pbp_su3_X_PRECISION_new( buf2, eta_pt, nvec, nvec, nvec_eta );
    }
    // minus mu direction
    for ( i=bbl[2*mu+1]; i<bbl[2*mu+2]; i+=2 ) {
      index = s->block[k].bt[i];
      neighbor_index = s->block[k].bt[i+1];
      D_pt = D + 36*neighbor_index + 9*mu;
      phi_pt = phi->vector_buffer + 12*neighbor_index*nvec_phi;
      eta_pt = eta->vector_buffer + 12*index*nvec_eta;
      prn_X_PRECISION_new( buf1, phi_pt, nvec, nvec, nvec_phi );
      mvmh_PRECISION_new( buf2, D_pt, buf1, nvec, nvec, nvec );
      mvmh_PRECISION_new( buf2+3*nvec, D_pt, buf1+3*nvec, nvec, nvec, nvec );
      pbn_su3_X_PRECISION_new( buf2, eta_pt, nvec, nvec, nvec_eta );
    }  
/*#ifdef HAVE_TM1p1
  }
#endif*/
}

// eta <- n_block_PRECISION_boundary_op*phi
void n_block_PRECISION_boundary_op_new( vector_PRECISION *eta, vector_PRECISION *phi, int k,
                                    schwarz_PRECISION_struct *s, level_struct *l ) {
  // k: number of current block
  int *bbl = s->block_boundary_length;
  int nvec = phi->num_vect_now, nvec_phi = phi->num_vect, nvec_eta = eta->num_vect;//!!!!!!
  int i, mu, index, neighbor_index;
  config_PRECISION D_pt, D = s->op.D; //pt to complex: D_pt = 3x3 color matrix & there are 4 such matrices
  buffer_PRECISION phi_pt, eta_pt;
  complex_PRECISION buf1[12*nvec], *buf2=buf1+6*nvec;

  mu=T;
  // plus mu direction
  for ( i=bbl[2*mu]; i<bbl[2*mu+1]; i+=2 ) {
    index = s->block[k].bt[i];
    neighbor_index = s->block[k].bt[i+1];
    D_pt = D + 36*index + 9*mu;
    phi_pt = phi->vector_buffer + 12*neighbor_index*nvec_phi;
    eta_pt = eta->vector_buffer + 12*index*nvec_eta;
    prp_T_PRECISION_new( buf1, phi_pt, nvec, nvec, nvec_phi );
    nmvm_PRECISION_new( buf2, D_pt, buf1, nvec, nvec, nvec );
    nmvm_PRECISION_new( buf2+3*nvec, D_pt, buf1+3*nvec, nvec, nvec, nvec );
    pbp_su3_T_PRECISION_new( buf2, eta_pt, nvec, nvec, nvec_eta );
  }  
  // minus mu direction
    for ( i=bbl[2*mu+1]; i<bbl[2*mu+2]; i+=2 ) {
      index = s->block[k].bt[i];
      neighbor_index = s->block[k].bt[i+1];
      D_pt = D + 36*neighbor_index + 9*mu;
      phi_pt = phi->vector_buffer + 12*neighbor_index*nvec_phi;
      eta_pt = eta->vector_buffer + 12*index*nvec_eta;
      prn_T_PRECISION_new( buf1, phi_pt, nvec, nvec, nvec_phi );
      nmvmh_PRECISION_new( buf2, D_pt, buf1, nvec, nvec, nvec );
      nmvmh_PRECISION_new( buf2+3*nvec, D_pt, buf1+3*nvec, nvec, nvec, nvec );
      pbn_su3_T_PRECISION_new( buf2, eta_pt, nvec, nvec, nvec_eta );
    }
    
    mu=Z;
    // plus mu direction
    for ( i=bbl[2*mu]; i<bbl[2*mu+1]; i+=2 ) {
      index = s->block[k].bt[i];
      neighbor_index = s->block[k].bt[i+1];
      D_pt = D + 36*index + 9*mu;
      phi_pt = phi->vector_buffer + 12*neighbor_index*nvec_phi;
      eta_pt = eta->vector_buffer + 12*index*nvec_eta;
      prp_Z_PRECISION_new( buf1, phi_pt, nvec, nvec, nvec_phi );
      nmvm_PRECISION_new( buf2, D_pt, buf1, nvec, nvec, nvec );
      nmvm_PRECISION_new( buf2+3*nvec, D_pt, buf1+3*nvec, nvec, nvec, nvec );
      pbp_su3_Z_PRECISION_new( buf2, eta_pt, nvec, nvec, nvec_eta );
    }
    // minus mu direction
    for ( i=bbl[2*mu+1]; i<bbl[2*mu+2]; i+=2 ) {
      index = s->block[k].bt[i];
      neighbor_index = s->block[k].bt[i+1];
      D_pt = D + 36*neighbor_index + 9*mu;
      phi_pt = phi->vector_buffer + 12*neighbor_index*nvec_phi;
      eta_pt = eta->vector_buffer + 12*index*nvec_eta;
      prn_Z_PRECISION_new( buf1, phi_pt, nvec, nvec, nvec_phi );
      nmvmh_PRECISION_new( buf2, D_pt, buf1, nvec, nvec, nvec );
      nmvmh_PRECISION_new( buf2+3*nvec, D_pt, buf1+3*nvec, nvec, nvec, nvec );
      pbn_su3_Z_PRECISION_new( buf2, eta_pt, nvec, nvec, nvec_eta );
    }
    
    mu=Y;
    // plus mu direction
    for ( i=bbl[2*mu]; i<bbl[2*mu+1]; i+=2 ) {
      index = s->block[k].bt[i];
      neighbor_index = s->block[k].bt[i+1];
      D_pt = D + 36*index + 9*mu;
      phi_pt = phi->vector_buffer + 12*neighbor_index*nvec_phi;
      eta_pt = eta->vector_buffer + 12*index*nvec_eta;
      prp_Y_PRECISION_new( buf1, phi_pt, nvec, nvec, nvec_phi );
      nmvm_PRECISION_new( buf2, D_pt, buf1, nvec, nvec, nvec );
      nmvm_PRECISION_new( buf2+3*nvec, D_pt, buf1+3*nvec, nvec, nvec, nvec );
      pbp_su3_Y_PRECISION_new( buf2, eta_pt, nvec, nvec, nvec_eta );
    }
    // minus mu direction
    for ( i=bbl[2*mu+1]; i<bbl[2*mu+2]; i+=2 ) {
      index = s->block[k].bt[i];
      neighbor_index = s->block[k].bt[i+1];
      D_pt = D + 36*neighbor_index + 9*mu;
      phi_pt = phi->vector_buffer + 12*neighbor_index*nvec_phi;
      eta_pt = eta->vector_buffer + 12*index*nvec_eta;
      prn_Y_PRECISION_new( buf1, phi_pt, nvec, nvec, nvec_phi );
      nmvmh_PRECISION_new( buf2, D_pt, buf1, nvec, nvec, nvec );
      nmvmh_PRECISION_new( buf2+3*nvec, D_pt, buf1+3*nvec, nvec, nvec, nvec );
      pbn_su3_Y_PRECISION_new( buf2, eta_pt, nvec, nvec, nvec_eta );
    }
    
    mu=X;
    // plus mu direction
    for ( i=bbl[2*mu]; i<bbl[2*mu+1]; i+=2 ) {
      index = s->block[k].bt[i];
      neighbor_index = s->block[k].bt[i+1];
      D_pt = D + 36*index + 9*mu;
      phi_pt = phi->vector_buffer + 12*neighbor_index*nvec_phi;
      eta_pt = eta->vector_buffer + 12*index*nvec_eta;
      prp_X_PRECISION_new( buf1, phi_pt, nvec, nvec, nvec_phi );
      nmvm_PRECISION_new( buf2, D_pt, buf1, nvec, nvec, nvec );
      nmvm_PRECISION_new( buf2+3*nvec, D_pt, buf1+3*nvec, nvec, nvec, nvec );
      pbp_su3_X_PRECISION_new( buf2, eta_pt, nvec, nvec, nvec_eta );
    }
    // minus mu direction
    for ( i=bbl[2*mu+1]; i<bbl[2*mu+2]; i+=2 ) {
      index = s->block[k].bt[i];
      neighbor_index = s->block[k].bt[i+1];
      D_pt = D + 36*neighbor_index + 9*mu;
      phi_pt = phi->vector_buffer + 12*neighbor_index*nvec_phi;
      eta_pt = eta->vector_buffer + 12*index*nvec_eta;
      prn_X_PRECISION_new( buf1, phi_pt, nvec, nvec, nvec_phi );
      nmvmh_PRECISION_new( buf2, D_pt, buf1, nvec, nvec, nvec );
      nmvmh_PRECISION_new( buf2+3*nvec, D_pt, buf1+3*nvec, nvec, nvec, nvec );
      pbn_su3_X_PRECISION_new( buf2, eta_pt, nvec, nvec, nvec_eta );
    } 

/*#ifdef HAVE_TM1p1
  }
#endif*/ 
}

// eta <- coarse_block_PRECISION_boundary_op*phi
void coarse_block_PRECISION_boundary_op_new( vector_PRECISION *eta, vector_PRECISION *phi,
                                         int k, schwarz_PRECISION_struct *s, level_struct *l ) {
  // k: number of current block
  int i, mu, index, neighbor_index;
  int *bbl = s->block_boundary_length;
  int n    = l->num_lattice_site_var;
  int nvec = phi->num_vect_now, nvec_phi = phi->num_vect, nvec_eta = eta->num_vect;
  vector_PRECISION phi_pt, eta_pt;

  config_PRECISION D_pt, D = s->op.D;
  int link_size = SQUARE(2*l->num_parent_eig_vect), site_size=4*link_size;

  vector_PRECISION_duplicate( &phi_pt, phi, 0, l );
  vector_PRECISION_duplicate( &eta_pt, eta, 0, l );

  for ( mu=0; mu<4; mu++ ) {
    // plus mu direction
    for ( i=bbl[2*mu]; i<bbl[2*mu+1]; i+=2 ) {
      index = s->block[k].bt[i];
      neighbor_index = s->block[k].bt[i+1];
      phi_pt.vector_buffer = phi->vector_buffer + n*neighbor_index*nvec_phi;
      eta_pt.vector_buffer = eta->vector_buffer + n*index*nvec_eta;
      D_pt = D + site_size*index + link_size*mu;
      coarse_hopp_PRECISION_new( &eta_pt, &phi_pt, D_pt, l );
    }
    // minus mu direction
    for ( i=bbl[2*mu+1]; i<bbl[2*mu+2]; i+=2 ) {
      index = s->block[k].bt[i];
      neighbor_index = s->block[k].bt[i+1];
      phi_pt.vector_buffer = phi->vector_buffer + n*neighbor_index*nvec_phi;
      eta_pt.vector_buffer = eta->vector_buffer + n*index*nvec_eta;
      D_pt = D + site_size*neighbor_index + link_size*mu;
      coarse_daggered_hopp_PRECISION_new( &eta_pt, &phi_pt, D_pt, l );
    }
  }
}

// eta <- n_coarse_block_PRECISION_boundary_op*phi
void n_coarse_block_PRECISION_boundary_op_new( vector_PRECISION *eta, vector_PRECISION *phi,
					       int k, schwarz_PRECISION_struct *s, level_struct *l ) {
  // k: number of current block
  int i, mu, index, neighbor_index;
  int *bbl = s->block_boundary_length;
  int n    = l->num_lattice_site_var;
  int link_size = SQUARE(2*l->num_parent_eig_vect);
  int site_size = 4*link_size;
  int nvec_phi = phi->num_vect, nvec_eta = eta->num_vect;
  config_PRECISION D_pt, D = s->op.D;
  vector_PRECISION phi_pt, eta_pt; 

  vector_PRECISION_duplicate( &phi_pt, phi, 0, l );
  vector_PRECISION_duplicate( &eta_pt, eta, 0, l );

  for ( mu=0; mu<4; mu++ ) {
    // plus mu direction
    for ( i=bbl[2*mu]; i<bbl[2*mu+1]; i+=2 ) {
      index = s->block[k].bt[i];
      neighbor_index = s->block[k].bt[i+1];
      phi_pt.vector_buffer = phi->vector_buffer + n*neighbor_index*nvec_phi;
      eta_pt.vector_buffer = eta->vector_buffer + n*index*nvec_eta;
      D_pt = D + site_size*index + link_size*mu;
      coarse_n_hopp_PRECISION_new( &eta_pt, &phi_pt, D_pt, l );
    }
    // minus mu direction
    for ( i=bbl[2*mu+1]; i<bbl[2*mu+2]; i+=2 ) {
      index = s->block[k].bt[i];
      neighbor_index = s->block[k].bt[i+1];
      phi_pt.vector_buffer = phi->vector_buffer + n*neighbor_index*nvec_phi;
      eta_pt.vector_buffer = eta->vector_buffer + n*index*nvec_eta;
      D_pt = D + site_size*neighbor_index + link_size*mu;
      coarse_n_daggered_hopp_PRECISION_new( &eta_pt, &phi_pt, D_pt, l );
    }
  }
}

void schwarz_PRECISION_new( vector_PRECISION *phi, vector_PRECISION *D_phi, vector_PRECISION *eta, const int cycles, int res,
                        schwarz_PRECISION_struct *s, level_struct *l, struct Thread *threading ) {
  
  START_NO_HYPERTHREADS(threading)

  int color, k, mu, i, nb = s->num_blocks, init_res = res, nvec = phi->num_vect_now, jj, jjj;
  vector_PRECISION *r = &(s->buf[0]), *Dphi = &(s->buf[3]), *latest_iter = &(s->buf[1]), *x = &(s->buf[2]);
  void (*block_op)()      = (l->depth==0)?block_d_plus_clover_PRECISION_new:            coarse_block_operator_PRECISION_new;
  void (*boundary_op)()   = (l->depth==0)?block_PRECISION_boundary_op_new:              coarse_block_PRECISION_boundary_op_new,
       (*n_boundary_op)() = (l->depth==0)?n_block_PRECISION_boundary_op_new:            n_coarse_block_PRECISION_boundary_op_new,
       (*block_solve)()   = (l->depth==0&&g.odd_even)?block_solve_oddeven_PRECISION_new:local_minres_PRECISION_new;
  
  complex_PRECISION relax_factor[nvec];
  VECTOR_LOOP(jj, nvec, jjj, relax_factor[jj+jjj]=l->relax_fac;)

  SYNC_CORES(threading)
  
  int nb_thread_start;
  int nb_thread_end;
  compute_core_start_end_custom(0, nb, &nb_thread_start, &nb_thread_end, l, threading, 1);
  
  g.num_vect_pass1 = nvec;//temp fix!!!
  r->num_vect_now = nvec; Dphi->num_vect_now = nvec; latest_iter->num_vect_now; x->num_vect_now = nvec;
  if ( res == _NO_RES ) {
    vector_PRECISION_copy_new( r, eta, nb_thread_start*s->block_vector_size, nb_thread_end*s->block_vector_size, l );
    vector_PRECISION_define_new( x, 0, nb_thread_start*s->block_vector_size, nb_thread_end*s->block_vector_size, l );
  } else {
    vector_PRECISION_copy_new( x, phi, nb_thread_start*s->block_vector_size, nb_thread_end*s->block_vector_size, l );
  }
    
  START_MASTER(threading)
  if ( res == _NO_RES ) {
    vector_PRECISION_define_new( x, 0, l->inner_vector_size, l->schwarz_vector_size, l );
  }
  END_MASTER(threading)
  
  SYNC_CORES(threading)
  
  for ( k=0; k<cycles; k++ ) {
    
    for ( color=0; color<s->num_colors; color++ ) {
      if ( res == _RES ) {
        START_LOCKED_MASTER(threading)
        for ( mu=0; mu<4; mu++ ) {
	  g.num_vect_pass2 = (k==0 && init_res == _RES)?x->num_vect:latest_iter->num_vect;
          ghost_update_PRECISION_new( (k==0 && init_res == _RES)?x:latest_iter, mu, +1, &(s->op.c), l );
          ghost_update_PRECISION_new( (k==0 && init_res == _RES)?x:latest_iter, mu, -1, &(s->op.c), l );
        }
        END_LOCKED_MASTER(threading)
      } else {
        // we need a barrier between black and white blocks
        SYNC_CORES(threading)
      }
        
      for ( i=nb_thread_start; i<nb_thread_end; i++ ) {
        // for all blocks of current color NOT involved in communication
        if ( color == s->block[i].color && s->block[i].no_comm ) {
          // calculate block residual
          START_MASTER(threading)
          PROF_PRECISION_START( _SM1 );
          END_MASTER(threading)
          if ( res == _RES ) {
            if ( k==0 && init_res == _RES ) {
              block_op( Dphi, x, s->block[i].start*l->num_lattice_site_var, s, l, no_threading );
              boundary_op( Dphi, x, i, s, l, no_threading );
              vector_PRECISION_minus_new( r, eta, Dphi, s->block[i].start*l->num_lattice_site_var,
                                      s->block[i].start*l->num_lattice_site_var+s->block_vector_size, l );
            } else {
              n_boundary_op( r, latest_iter, i, s, l );
            }
          }
          START_MASTER(threading)
          PROF_PRECISION_STOP( _SM1, 1 );
          // local minres updates x, r and latest iter
          PROF_PRECISION_START( _SM2 );
          END_MASTER(threading)
          block_solve( x, r, latest_iter, s->block[i].start*l->num_lattice_site_var, s, l, no_threading );
          START_MASTER(threading)
          PROF_PRECISION_STOP( _SM2, 1 );
          END_MASTER(threading)
        }
      }
      
      if ( res == _RES ) {
        START_LOCKED_MASTER(threading)
	  g.num_vect_pass2 = (k==0 && init_res == _RES)?x->num_vect:latest_iter->num_vect;
        for ( mu=0; mu<4; mu++ ) {
          ghost_update_wait_PRECISION( (k==0 && init_res == _RES)?x:latest_iter, mu, +1, &(s->op.c), l );
          ghost_update_wait_PRECISION( (k==0 && init_res == _RES)?x:latest_iter, mu, -1, &(s->op.c), l );
        }
        END_LOCKED_MASTER(threading)
      } else {
        // we need a barrier between black and white blocks
        SYNC_CORES(threading)
      }
        
      for ( i=nb_thread_start; i<nb_thread_end; i++ ) {
        // for all blocks of current color involved in communication
        if ( color == s->block[i].color && !s->block[i].no_comm ) {
          // calculate block residual
          START_MASTER(threading)
          PROF_PRECISION_START( _SM3 );
          END_MASTER(threading)
          if ( res == _RES ) {
            if ( k==0 && init_res == _RES ) {
              block_op( Dphi, x, s->block[i].start*l->num_lattice_site_var, s, l, no_threading );
              boundary_op( Dphi, x, i, s, l, no_threading );
              vector_PRECISION_minus_new( r, eta, Dphi, s->block[i].start*l->num_lattice_site_var,
                                      s->block[i].start*l->num_lattice_site_var+s->block_vector_size, l );
            } else {
              n_boundary_op( r, latest_iter, i, s, l );
            }
          }
          START_MASTER(threading)
          PROF_PRECISION_STOP( _SM3, 1 );
          // local minres updates x, r and latest iter
          PROF_PRECISION_START( _SM4 );
          END_MASTER(threading)
	    block_solve( x, r, latest_iter, s->block[i].start*l->num_lattice_site_var, s, l, no_threading );//no_!!!!!!
          START_MASTER(threading)
          PROF_PRECISION_STOP( _SM4, 1 );
          END_MASTER(threading)
        }
      }
      res = _RES;
    }
  }
  
  for ( i=nb_thread_start; i<nb_thread_end; i++ ) {
    if ( l->relax_fac != 1.0 )
      vector_PRECISION_scale_new( phi, x, relax_factor, 0, s->block[i].start*l->num_lattice_site_var, s->block[i].start*l->num_lattice_site_var+s->block_vector_size, l );
    else
      vector_PRECISION_copy_new( phi, x, s->block[i].start*l->num_lattice_site_var, s->block[i].start*l->num_lattice_site_var+s->block_vector_size, l );
  }
  
  // calculate D * phi with help of the almost computed residual
  // via updating the residual from odd to even
  if ( D_phi != NULL ) {
    START_LOCKED_MASTER(threading)
    g.num_vect_pass2 = latest_iter->num_vect;
    for ( mu=0; mu<4; mu++ ) {
      ghost_update_PRECISION_new( latest_iter, mu, +1, &(s->op.c), l );
      ghost_update_PRECISION_new( latest_iter, mu, -1, &(s->op.c), l );
    }
    END_LOCKED_MASTER(threading)
    
    for ( i=nb_thread_start; i<nb_thread_end; i++ ) {
      if ( 0 == s->block[i].color && s->block[i].no_comm ) {
        n_boundary_op( r, latest_iter, i, s, l );
        vector_PRECISION_minus_new( D_phi, eta, r, s->block[i].start*l->num_lattice_site_var,
                                s->block[i].start*l->num_lattice_site_var+s->block_vector_size, l );
        if ( l->relax_fac != 1.0 ) {
          vector_PRECISION_scale_new( D_phi, D_phi, relax_factor, 0, s->block[i].start*l->num_lattice_site_var,
                                  s->block[i].start*l->num_lattice_site_var+s->block_vector_size, l );
        }
      }
      if ( 1 == s->block[i].color ) {
        vector_PRECISION_minus_new( D_phi, eta, r, s->block[i].start*l->num_lattice_site_var,
                                s->block[i].start*l->num_lattice_site_var+s->block_vector_size, l );
        if ( l->relax_fac != 1.0 ) {
          vector_PRECISION_scale_new( D_phi, D_phi, relax_factor, 0, s->block[i].start*l->num_lattice_site_var,
                                  s->block[i].start*l->num_lattice_site_var+s->block_vector_size, l );
        }
      }
    }
    
    START_LOCKED_MASTER(threading)
    for ( mu=0; mu<4; mu++ ) {
      ghost_update_wait_PRECISION( latest_iter, mu, +1, &(s->op.c), l );
      ghost_update_wait_PRECISION( latest_iter, mu, -1, &(s->op.c), l );
    }
    END_LOCKED_MASTER(threading)
    
    for ( i=nb_thread_start; i<nb_thread_end; i++ ) {
      if ( 0 == s->block[i].color && !s->block[i].no_comm ) {
        n_boundary_op( r, latest_iter, i, s, l );
        vector_PRECISION_minus_new( D_phi, eta, r, s->block[i].start*l->num_lattice_site_var,
                                s->block[i].start*l->num_lattice_site_var+s->block_vector_size, l );
        if ( l->relax_fac != 1.0 ) {
          vector_PRECISION_scale_new( D_phi, D_phi, relax_factor, 0, s->block[i].start*l->num_lattice_site_var,
                                  s->block[i].start*l->num_lattice_site_var+s->block_vector_size, l );
        }
      }
    }    
  }
  SYNC_CORES(threading)
  
#ifdef SCHWARZ_RES
  START_LOCKED_MASTER(threading)
  if ( D_phi == NULL ) {
    for ( mu=0; mu<4; mu++ ) {
      g.num_vect_pass2 = latest_iter->num_vect;
      ghost_update_PRECISION_new( latest_iter, mu, +1, &(s->op.c), l );
      ghost_update_PRECISION_new( latest_iter, mu, -1, &(s->op.c), l );
    }
    
    for ( i=0; i<nb; i++ )
      if ( s->block[i].no_comm ) 
        n_boundary_op( r, latest_iter, i, s, l );

    
    for ( mu=0; mu<4; mu++ ) {
      ghost_update_wait_PRECISION( latest_iter, mu, +1, &(s->op.c), l );
      ghost_update_wait_PRECISION( latest_iter, mu, -1, &(s->op.c), l );
    }
    
    for ( i=0; i<nb; i++ )
      if ( !s->block[i].no_comm ) 
        n_boundary_op( r, latest_iter, i, s, l );
  }

  PRECISION r_norm[nvec];
  global_norm_PRECISION_new( r_norm, r, 0, l->inner_vector_size, l, no_threading );
  char number[3]; sprintf( number, "%2d", 31+l->depth ); printf0("\033[1;%2sm|", number );
  for( i=0; i<nvec; i++ )
    printf0(" ---- depth: %d, c: %d, schwarz iter %2d, vector %d, norm: %11.6le |", l->depth, s->num_colors, k, i, r_norm[i] );
  printf0("\033[0m\n"); fflush(0);
  END_LOCKED_MASTER(threading)
#endif

  END_NO_HYPERTHREADS(threading)
}

void additive_schwarz_PRECISION_new( vector_PRECISION *phi, vector_PRECISION *D_phi, vector_PRECISION *eta, const int cycles, int res, 
                                 schwarz_PRECISION_struct *s, level_struct *l, struct Thread *threading ) {
  
  START_NO_HYPERTHREADS(threading)
  
  int k, mu, i, nb = s->num_blocks, nvec = phi->num_vect_now, jj, jjj;
  vector_PRECISION *r = &(s->buf[0]), *Dphi = &(s->buf[3]), *latest_iter = &(s->buf[1]), *x = &(s->buf[2]), *latest_iter2 = &(s->buf[4]), *swap = NULL;
  void (*block_op)()      = (l->depth==0)?block_d_plus_clover_PRECISION_new:            coarse_block_operator_PRECISION_new;
  void (*boundary_op)()   = (l->depth==0)?block_PRECISION_boundary_op_new:              coarse_block_PRECISION_boundary_op_new,
       (*n_boundary_op)() = (l->depth==0)?n_block_PRECISION_boundary_op_new:            n_coarse_block_PRECISION_boundary_op_new,
       (*block_solve)()   = (l->depth==0&&g.odd_even)?block_solve_oddeven_PRECISION_new:local_minres_PRECISION_new;

  complex_PRECISION relax_factor[nvec];
  VECTOR_LOOP(jj, nvec, jjj, relax_factor[jj+jjj]=l->relax_fac;)

  int nb_thread_start;
  int nb_thread_end;
  compute_core_start_end_custom(0, nb, &nb_thread_start, &nb_thread_end, l, threading, 1);
  
  SYNC_CORES(threading)
  
  g.num_vect_pass1 = nvec;
  r->num_vect_now = nvec; Dphi->num_vect_now = nvec; latest_iter->num_vect_now = nvec; latest_iter->num_vect_now = nvec; x->num_vect_now = nvec;
  if ( res == _NO_RES ) {
    vector_PRECISION_copy_new( r, eta, nb_thread_start*s->block_vector_size, nb_thread_end*s->block_vector_size, l );
    vector_PRECISION_define_new( x, 0, nb_thread_start*s->block_vector_size, nb_thread_end*s->block_vector_size, l );
    START_MASTER(threading)
    vector_PRECISION_define_new( x, 0, l->inner_vector_size, l->schwarz_vector_size, l );
    END_MASTER(threading)
  } else {
    vector_PRECISION_copy_new( x, phi, nb_thread_start*s->block_vector_size, nb_thread_end*s->block_vector_size, l );
    vector_PRECISION_copy_new( latest_iter, phi, nb_thread_start*s->block_vector_size, nb_thread_end*s->block_vector_size, l );
  }
  
  
  SYNC_CORES(threading)
  
  for ( k=0; k<cycles; k++ ) {
    if ( res == _RES ) {
      START_LOCKED_MASTER(threading)
      for ( mu=0; mu<4; mu++ ) {
	g.num_vect_pass2 = latest_iter->num_vect;
        ghost_update_PRECISION_new( latest_iter, mu, +1, &(s->op.c), l );
        ghost_update_PRECISION_new( latest_iter, mu, -1, &(s->op.c), l );
      }
      END_LOCKED_MASTER(threading)
    } else {
      SYNC_CORES(threading)
    }
      
    for ( i=nb_thread_start; i<nb_thread_end; i++ ) {
      // for all blocks of current color NOT involved in communication
      if ( s->block[i].no_comm ) {
        // calculate block residual
        if ( res == _RES ) {
          if ( k==0 ) {
            block_op( Dphi, latest_iter, s->block[i].start*l->num_lattice_site_var, s, l, no_threading );
            boundary_op( Dphi, latest_iter, i, s, l, no_threading );
            vector_PRECISION_minus_new( r, eta, Dphi, s->block[i].start*l->num_lattice_site_var,
                                    s->block[i].start*l->num_lattice_site_var+s->block_vector_size, l );
        } else {
            n_boundary_op( r, latest_iter, i, s, l );
          }
        }
        // local minres updates x, r and latest iter
        block_solve( x, r, latest_iter2, s->block[i].start*l->num_lattice_site_var, s, l, no_threading );
      }
    }
    
    if ( res == _RES ) {
      START_LOCKED_MASTER(threading)
      for ( mu=0; mu<4; mu++ ) {
	g.num_vect_pass2 = latest_iter->num_vect;
        ghost_update_wait_PRECISION( latest_iter, mu, +1, &(s->op.c), l );
        ghost_update_wait_PRECISION( latest_iter, mu, -1, &(s->op.c), l );
      }
      END_LOCKED_MASTER(threading)
    } else {
      SYNC_CORES(threading)
    }
      
    for ( i=nb_thread_start; i<nb_thread_end; i++ ) {
      // for all blocks of current color involved in communication
      if ( !s->block[i].no_comm ) {
        // calculate block residual
        if ( res == _RES ) {
          if ( k==0 ) {
            block_op( Dphi, latest_iter, s->block[i].start*l->num_lattice_site_var, s, l, no_threading );
            boundary_op( Dphi, latest_iter, i, s, l, no_threading );
            vector_PRECISION_minus_new( r, eta, Dphi, s->block[i].start*l->num_lattice_site_var,
                                    s->block[i].start*l->num_lattice_site_var+s->block_vector_size, l );
          } else {
            n_boundary_op( r, latest_iter, i, s, l );
          }
        }
        // local minres updates x, r and latest iter
        block_solve( x, r, latest_iter2, s->block[i].start*l->num_lattice_site_var, s, l, no_threading );
      }
    }
    res = _RES;
    swap = latest_iter; latest_iter = latest_iter2; latest_iter2 = swap;
  }

  for ( i=nb_thread_start; i<nb_thread_end; i++ ) {
    if ( l->relax_fac != 1.0 )
      vector_PRECISION_scale_new( phi, x, relax_factor, 0, s->block[i].start*l->num_lattice_site_var, s->block[i].start*l->num_lattice_site_var+s->block_vector_size, l );
    else
      vector_PRECISION_copy_new( phi, x, s->block[i].start*l->num_lattice_site_var, s->block[i].start*l->num_lattice_site_var+s->block_vector_size, l );
  }
  
  // calculate D * phi with help of the almost computed residual
  if ( D_phi != NULL ) {
    START_LOCKED_MASTER(threading)
    for ( mu=0; mu<4; mu++ ) {
      g.num_vect_pass2 = latest_iter->num_vect;
      ghost_update_PRECISION_new( latest_iter, mu, +1, &(s->op.c), l );
      ghost_update_PRECISION_new( latest_iter, mu, -1, &(s->op.c), l );
    }
    END_LOCKED_MASTER(threading)
    
    for ( i=nb_thread_start; i<nb_thread_end; i++ ) {
      if ( s->block[i].no_comm ) {
        n_boundary_op( r, latest_iter, i, s, l );
        vector_PRECISION_minus_new( D_phi, eta, r,
            s->block[i].start*l->num_lattice_site_var,
            s->block[i].start*l->num_lattice_site_var+s->block_vector_size, l );
        if ( l->relax_fac != 1.0 )
          vector_PRECISION_scale_new( D_phi, D_phi, relax_factor, 0, s->block[i].start*l->num_lattice_site_var,
              s->block[i].start*l->num_lattice_site_var+s->block_vector_size, l );
      }
    }
    
    START_LOCKED_MASTER(threading)
    for ( mu=0; mu<4; mu++ ) {
      ghost_update_wait_PRECISION( latest_iter, mu, +1, &(s->op.c), l );
      ghost_update_wait_PRECISION( latest_iter, mu, -1, &(s->op.c), l );
    }
    END_LOCKED_MASTER(threading)
    
    for ( i=nb_thread_start; i<nb_thread_end; i++ ) {
      if ( !s->block[i].no_comm ) {
        n_boundary_op( r, latest_iter, i, s, l );
        vector_PRECISION_minus_new( D_phi, eta, r,
            s->block[i].start*l->num_lattice_site_var,
            s->block[i].start*l->num_lattice_site_var+s->block_vector_size, l );
        if ( l->relax_fac != 1.0 )
          vector_PRECISION_scale_new( D_phi, D_phi, relax_factor, 0, s->block[i].start*l->num_lattice_site_var,
              s->block[i].start*l->num_lattice_site_var+s->block_vector_size, l );
      }
    }
  }
  SYNC_CORES(threading)
  
#ifdef SCHWARZ_RES
  START_LOCKED_MASTER(threading)
  if ( D_phi == NULL ) {
    for ( mu=0; mu<4; mu++ ) {
      g.num_vect_pass2 = latest_iter->num_vect;
      ghost_update_PRECISION_new( latest_iter, mu, +1, &(s->op.c), l );
      ghost_update_PRECISION_new( latest_iter, mu, -1, &(s->op.c), l );
    }
    
    for ( i=0; i<nb; i++ ) {
      if ( s->block[i].no_comm ) {
        n_boundary_op( r, latest_iter, i, s, l );
      }
    }
    
    for ( mu=0; mu<4; mu++ ) {
      g.num_vect_pass2 = latest_iter->num_vect;
      ghost_update_wait_PRECISION( latest_iter, mu, +1, &(s->op.c), l );
      ghost_update_wait_PRECISION( latest_iter, mu, -1, &(s->op.c), l );
    }
    
    for ( i=0; i<nb; i++ ) {
      if ( !s->block[i].no_comm ) {
        n_boundary_op( r, latest_iter, i, s, l );
      }
    }
  }
  PRECISION r_norm[nvec];
  global_norm_PRECISION_new( r_norm, r, 0, l->inner_vector_size, l, no_threading );
  char number[3]; sprintf( number, "%2d", 31+l->depth ); printf0("\033[1;%2sm|", number );
  for( i=0; i<nvec; i++)
    printf0(" ---- depth: %d, c: %d, schwarz iter %2d, vector %d, norm: %11.6le |", l->depth, s->num_colors, k, i, r_norm[i] );
  printf0("\033[0m\n"); fflush(0);
  END_LOCKED_MASTER(threading)
#endif

  END_NO_HYPERTHREADS(threading)
}

void red_black_schwarz_PRECISION_new( vector_PRECISION *phi, vector_PRECISION *D_phi, vector_PRECISION *eta, const int cycles, int res,
				      schwarz_PRECISION_struct *s, level_struct *l, struct Thread *threading ) {
  START_NO_HYPERTHREADS(threading)
  //initial step: phi _ORDINARY eta: INNER; r at depth 0 is INNER
  int jj, jjj, k=0, mu, i, init_res = res, res_comm = res, step;
  int nvec = phi->num_vect_now;
  vector_PRECISION *r = &(s->buf[0]), *Dphi = &(s->buf[3]), *latest_iter = &(s->buf[1]), *x = &(s->buf[2]);
  void (*block_op)()       = (l->depth==0)?block_d_plus_clover_PRECISION_new:            coarse_block_operator_PRECISION_new;
  void (*boundary_op)()    = (l->depth==0)?block_PRECISION_boundary_op_new:              coarse_block_PRECISION_boundary_op_new;
  void (*n_boundary_op)()  = (l->depth==0)?n_block_PRECISION_boundary_op_new:            n_coarse_block_PRECISION_boundary_op_new;
  void (*block_solve)()    = (l->depth==0&&g.odd_even)?block_solve_oddeven_PRECISION_new:local_minres_PRECISION_new;
  void (*communicate[2])() = {ghost_update_wait_PRECISION, ghost_update_PRECISION_new};
  int commdir[8] = {+1,-1,-1,+1,-1,+1,+1,-1};

  complex_PRECISION relax_factor[nvec];
  VECTOR_LOOP(jj, nvec, jjj, relax_factor[jj+jjj]=l->relax_fac;)
       
  SYNC_CORES(threading)

  int block_thread_start[8], block_thread_end[8];
  for ( i=0; i<8; i++ )
     compute_core_start_end_custom(0, s->block_list_length[i], block_thread_start+i, block_thread_end+i, l, threading, 1 );
  int start, end;
  compute_core_start_end_custom(0, l->inner_vector_size, &start, &end, l, threading, l->num_lattice_site_var );

  //---- initialization
  g.num_vect_pass1 = nvec;
  r->num_vect_now = nvec; Dphi->num_vect_now = nvec; latest_iter->num_vect_now = nvec; x->num_vect_now = nvec;
  if ( res == _NO_RES ) {
    // if the initial guess is zero: r <- eta, x = 0????
    vector_PRECISION_copy_new( r, eta, start, end, l );
    vector_PRECISION_define_new( x, 0, start, end, l );
    START_MASTER(threading)
    vector_PRECISION_define_new( x, 0, l->inner_vector_size, l->schwarz_vector_size, l );
    END_MASTER(threading)
    SYNC_CORES(threading)
  } else {
    // if phi contains the initial guess, x <- phi and update the ghost shell
    vector_PRECISION_copy_new( x, phi, start, end, l );
    START_LOCKED_MASTER(threading)
      g.num_vect_pass2 = x->num_vect;
    for ( mu=0; mu<4; mu++ )
      ghost_update_PRECISION_new( x, mu, +1, &(s->op.c), l );
    for ( mu=0; mu<4; mu++ )
      ghost_update_PRECISION_new( x, mu, -1, &(s->op.c), l );
    END_LOCKED_MASTER(threading)
  }

  g.num_vect_pass2 = (k==0 && step < 6 && init_res == _RES)?x->num_vect:latest_iter->num_vect;
  //--- perform the Schwarz iteration, solve the block systems
  // outer iteration over Scwarz cycles
  for ( k=0; k<cycles; k++ ) {
    for ( step=0; step<8; step++ ) {
      for ( i=block_thread_start[step]; i<block_thread_end[step]; i++ ) {
        int index = s->block_list[step][i];
        START_MASTER(threading)
        PROF_PRECISION_START( _SM3 );
        END_MASTER(threading)
        if ( res == _RES ) {
          if ( k==0 && init_res == _RES ) {
            block_op( Dphi, x, s->block[index].start*l->num_lattice_site_var, s, l, no_threading );
            boundary_op( Dphi, x, index, s, l, no_threading );
            vector_PRECISION_minus_new( r, eta, Dphi, s->block[index].start*l->num_lattice_site_var,
					s->block[index].start*l->num_lattice_site_var+s->block_vector_size, l );
          } else {
            n_boundary_op( r, latest_iter, index, s, l );
          }
        }
        START_MASTER(threading)
        PROF_PRECISION_STOP( _SM3, 1 );
        PROF_PRECISION_START( _SM4 );
        END_MASTER(threading)
	// local minres updates x, r and latest iter
	block_solve( x, r, latest_iter, s->block[index].start*l->num_lattice_site_var, s, l, no_threading );
	START_MASTER(threading)
	PROF_PRECISION_STOP( _SM4, 1 );
        END_MASTER(threading)
	  }

      if ( res_comm == _RES && !(k==cycles-1 && (step==6||step==7) && D_phi==NULL) ) {
        START_LOCKED_MASTER(threading)
	for ( mu=0; mu<4; mu++ ) {
	  communicate[(step%4)/2]( (k==0 && step < 6 && init_res == _RES)?x:latest_iter, mu, commdir[step], &(s->op.c), l );
        }
        END_LOCKED_MASTER(threading)
	} else {
        SYNC_CORES(threading)
	}
      
      if ( k==0 && step == 5 ) res = _RES;
      if ( k==0 && step == 1 ) res_comm = _RES;
    }
  }
  // copy phi = x
  if ( l->relax_fac != 1.0 )
    vector_PRECISION_scale_new( phi, x, relax_factor, 0, start, end, l );
  else
    vector_PRECISION_copy_new( phi, x, start, end, l );
  // calculate D * phi from r
  if ( D_phi != NULL ) {
    for ( step=4; step<8; step++ ) {
      for ( i=block_thread_start[step]; i<block_thread_end[step]; i++ ) {
        int index = s->block_list[step][i];
        vector_PRECISION_minus_new( D_phi, eta, r, s->block[index].start*l->num_lattice_site_var,
                                s->block[index].start*l->num_lattice_site_var+s->block_vector_size, l );
        if ( l->relax_fac != 1.0 ) {
          vector_PRECISION_scale_new( D_phi, D_phi, relax_factor, 0, s->block[index].start*l->num_lattice_site_var,
                                  s->block[index].start*l->num_lattice_site_var+s->block_vector_size, l );
        }
      }
    }
    
    for ( step=0; step<4; step++ ) {
      for ( i=block_thread_start[step]; i<block_thread_end[step]; i++ ) {
        int index = s->block_list[step][i];
        
        START_MASTER(threading)
        PROF_PRECISION_START( _SM3 );
        END_MASTER(threading)
        n_boundary_op( r, latest_iter, index, s, l );
        vector_PRECISION_minus_new( D_phi, eta, r, s->block[index].start*l->num_lattice_site_var,
                                s->block[index].start*l->num_lattice_site_var+s->block_vector_size, l );
        if ( l->relax_fac != 1.0 ) {
          vector_PRECISION_scale_new( D_phi, D_phi, relax_factor, 0, s->block[index].start*l->num_lattice_site_var,
                                  s->block[index].start*l->num_lattice_site_var+s->block_vector_size, l );
        }
        START_MASTER(threading)
        PROF_PRECISION_STOP( _SM3, 1 );
        END_MASTER(threading)
      }
      if ( step == 0 || step == 1 ) {
        START_LOCKED_MASTER(threading)
	  g.num_vect_pass2 = latest_iter->num_vect;
        for ( mu=0; mu<4; mu++ ) {
          communicate[0]( latest_iter, mu, commdir[step], &(s->op.c), l );
        }
        END_LOCKED_MASTER(threading)
      } else {
        SYNC_CORES(threading)
      }
    }
  }
  SYNC_CORES(threading)

#ifdef SCHWARZ_RES
  int nb = s->num_blocks;
  START_LOCKED_MASTER(threading)
  if ( D_phi == NULL ) {
    g.num_vect_pass2 = latest_iter->num_vect;
    for ( mu=0; mu<4; mu++ ) {
      ghost_update_PRECISION_new( latest_iter, mu, +1, &(s->op.c), l );
      ghost_update_PRECISION_new( latest_iter, mu, -1, &(s->op.c), l );
    }
    
    for ( i=0; i<nb; i++ )
      if ( s->block[i].no_comm )
        n_boundary_op( r, latest_iter, i, s, l );
    
    for ( mu=0; mu<4; mu++ ) {
      ghost_update_wait_PRECISION( latest_iter, mu, +1, &(s->op.c), l );
      ghost_update_wait_PRECISION( latest_iter, mu, -1, &(s->op.c), l );
    }
    
    for ( i=0; i<nb; i++ )
      if ( !s->block[i].no_comm )
        n_boundary_op( r, latest_iter, i, s, l );

  }
  PRECISION r_norm[nvec];
  global_norm_PRECISION_new( r_norm, r, 0, l->inner_vector_size, l, no_threading );
  char number[3]; sprintf( number, "%2d", 31+l->depth ); printf0("\033[1;%2sm|", number );
  for( i=0; i<nvec; i++)
    printf0(" ---- depth: %d, c: %d, schwarz iter %2d, vector %d, norm: %11.6le |", l->depth, s->num_colors, k, i, r_norm[i] );
  printf0("\033[0m\n"); fflush(0);
  END_LOCKED_MASTER(threading)
#endif
  END_NO_HYPERTHREADS(threading)
}

void sixteen_color_schwarz_PRECISION_new( vector_PRECISION *phi, vector_PRECISION *D_phi, vector_PRECISION *eta, const int cycles, int res, 
                                      schwarz_PRECISION_struct *s, level_struct *l, struct Thread *threading ) {
  
  START_NO_HYPERTHREADS(threading)
  ASSERT( D_phi == NULL );
  
  if ( s->num_colors == 2 ) schwarz_PRECISION_new( phi, D_phi, eta, cycles, res, s, l, no_threading );
  else {
    int color, k, mu, i, nb = s->num_blocks, nvec = phi->num_vect_now, jj, jjj;
    vector_PRECISION *r = &(s->buf[0]), *Dphi = &(s->buf[3]), *latest_iter = &(s->buf[1]), *x = &(s->buf[2]);
    void (*block_op)() = (l->depth==0)?block_d_plus_clover_PRECISION_new:               coarse_block_operator_PRECISION_new;
    void (*boundary_op)() = (l->depth==0)?block_PRECISION_boundary_op_new:              coarse_block_PRECISION_boundary_op_new,
         (*n_boundary_op)() = (l->depth==0)?n_block_PRECISION_boundary_op_new:          n_coarse_block_PRECISION_boundary_op_new,
         (*block_solve)() = (l->depth==0&&g.odd_even)?block_solve_oddeven_PRECISION_new:local_minres_PRECISION_new;
        
    int color_to_comm[16][2] = { {T,-1}, {X,+1}, {Y,+1}, {X,-1}, {Z,+1}, {Y,-1}, {X,+1}, {Y,+1},
                                {T,+1}, {X,-1}, {Y,-1}, {X,+1}, {Z,-1}, {Y,+1}, {X,-1}, {Y,-1}  };

    complex_PRECISION relax_factor[nvec];
    VECTOR_LOOP(jj, nvec, jjj, relax_factor[jj+jjj]=l->relax_fac;)
    
    int nb_thread_start;
    int nb_thread_end;
    compute_core_start_end_custom(0, nb, &nb_thread_start, &nb_thread_end, l, threading, 1);
    
    SYNC_CORES(threading)
    
    g.num_vect_pass1 = nvec;
    r->num_vect_now = nvec; Dphi->num_vect_now = nvec; latest_iter->num_vect_now = nvec; x->num_vect_now = nvec;
    if ( res == _NO_RES ) {
      vector_PRECISION_copy_new( r, eta, nb_thread_start*s->block_vector_size, nb_thread_end*s->block_vector_size, l );
      vector_PRECISION_define_new( x, 0, nb_thread_start*s->block_vector_size, nb_thread_end*s->block_vector_size, l );
    } else {
      vector_PRECISION_copy_new( x, phi, nb_thread_start*s->block_vector_size, nb_thread_end*s->block_vector_size, l );
    }
    
    START_MASTER(threading)
    if ( res == _NO_RES ) {
      vector_PRECISION_define_new( x, 0, l->inner_vector_size, l->schwarz_vector_size, l );
    }
    END_MASTER(threading)
    
    SYNC_CORES(threading)

    g.num_vect_pass2 = k==0?x->num_vect:latest_iter->num_vect;
    for ( k=0; k<cycles; k++ ) {
      for ( color=0; color<s->num_colors; color++ ) {
        
        // comm start
        if ( res == _RES ) {
          START_LOCKED_MASTER(threading)
          ghost_update_PRECISION_new( k==0?x:latest_iter, color_to_comm[color][0], color_to_comm[color][1], &(s->op.c), l );
          if ( color == 0 && (k==0 || k==1) ) {
            for ( mu=1; mu<4; mu++ ) {
              ghost_update_PRECISION_new( k==0?x:latest_iter, mu, -1, &(s->op.c), l );
            }
          }
          END_LOCKED_MASTER(threading)
        } else {
          SYNC_CORES(threading)
        }
        
        // blocks which have their latest neighbor information available
        for ( i=nb_thread_start; i<nb_thread_end; i++ ) {          
          if (  color == s->block[i].color && s->block[i].no_comm ) {
            if ( res == _RES ) {
              if ( k==0 ) {
                block_op( Dphi, x, s->block[i].start*l->num_lattice_site_var, s, l, no_threading );
                boundary_op( Dphi, x, i, s, l, no_threading );
                vector_PRECISION_minus_new( r, eta, Dphi, s->block[i].start*l->num_lattice_site_var,
                                        s->block[i].start*l->num_lattice_site_var+s->block_vector_size, l );
              } else {
                n_boundary_op( r, latest_iter, i, s, l );
              }
            }
            // local minres updates x, r and latest iter
            block_solve( x, r, latest_iter, s->block[i].start*l->num_lattice_site_var, s, l, no_threading );
          }
        }
        
        // comm wait
        if ( res == _RES ) {
          START_LOCKED_MASTER(threading)
          ghost_update_wait_PRECISION( k==0?x:latest_iter, color_to_comm[color][0], color_to_comm[color][1], &(s->op.c), l );
          if ( color == 0 && (k==0 || k==1) ) {
            for ( mu=1; mu<4; mu++ ) {
              ghost_update_wait_PRECISION( k==0?x:latest_iter, mu, -1, &(s->op.c), l );
            }
          }
          END_LOCKED_MASTER(threading)
        } else {
          SYNC_CORES(threading)
        }
        
        // blocks which require certain ghost cell updates
        for ( i=nb_thread_start; i<nb_thread_end; i++ ) {
          if ( color == s->block[i].color && !s->block[i].no_comm ) {
            if ( res == _RES ) {
              if ( k==0 ) {
                block_op( Dphi, x, s->block[i].start*l->num_lattice_site_var, s, l, no_threading );
                boundary_op( Dphi, x, i, s, l, no_threading );
                vector_PRECISION_minus_new( r, eta, Dphi, s->block[i].start*l->num_lattice_site_var,
                                        s->block[i].start*l->num_lattice_site_var+s->block_vector_size, l );
              } else {
                n_boundary_op( r, latest_iter, i, s, l );
              }
            }
            // local minres updates x, r and latest iter
            block_solve( x, r, latest_iter, s->block[i].start*l->num_lattice_site_var, s, l, no_threading );
          }
        }
             
        res = _RES;
      }
    }

    SYNC_CORES(threading)
    
    if ( l->relax_fac != 1.0 )
      vector_PRECISION_scale_new( phi, x, relax_factor, 0, nb_thread_start*s->block_vector_size, nb_thread_end*s->block_vector_size, l );
    else
      vector_PRECISION_copy_new( phi, x, nb_thread_start*s->block_vector_size, nb_thread_end*s->block_vector_size, l );
    
    SYNC_CORES(threading)
    
#ifdef SCHWARZ_RES
    START_LOCKED_MASTER(threading)
    vector_PRECISION true_r;
    vector_PRECISION_init(&true_r);

    vector_PRECISION_alloc( &true_r, _ORDINARY, nvec, l, threading );
    vector_PRECISION_define_new( &true_r, 0, 0, l->inner_vector_size, l );

    if ( D_phi == NULL ) {
      g.num_vect_pass2 = x->num_vect;
      for ( mu=0; mu<4; mu++ ) {
        ghost_update_PRECISION_new( x, mu, +1, &(s->op.c), l );
        ghost_update_PRECISION_new( x, mu, -1, &(s->op.c), l );
      }
      for ( i=0; i<nb; i++ ) {
        block_op( &true_r, x, s->block[i].start*l->num_lattice_site_var, s, l, no_threading );
      }
      for ( mu=0; mu<4; mu++ ) {
        ghost_update_wait_PRECISION( x, mu, +1, &(s->op.c), l );
        ghost_update_wait_PRECISION( x, mu, -1, &(s->op.c), l );
      }
      for ( i=0; i<nb; i++ ) {
        boundary_op( &true_r, x, i, s, l );
      }
    }
    //vector_PRECISION_saxpy( &true_r, eta, &true_r, -1, 0, l->inner_vector_size, l );
    complex_PRECISION factor[nvec];
    VECTOR_LOOP( jj, nvec, jjj, factor[jj+jjj]=-1;)
    vector_PRECISION_saxpy_new( &true_r, eta, &true_r, factor, 0, 1, 0, l->inner_vector_size, l );
    //PRECISION r_norm = global_norm_PRECISION( &true_r, 0, l->inner_vector_size, l, no_threading ),
    //  den = global_norm_PRECISION( eta, 0, l->inner_vector_size, l, no_threading );
    // r_norm/=den;
    PRECISION r_norm[nvec], den[nvec];
    global_norm_PRECISION_new( r_norm, &true_r, 0, l->inner_vector_size, l, no_threading );
    global_norm_PRECISION_new( den, eta, 0, l->inner_vector_size, l, no_threading );
    VECTOR_LOOP( jj, nvec, jjj, r_norm[jj+jjj] /= den[jj+jjj];)
    char number[3]; sprintf( number, "%2d", 31+l->depth ); printf0("\033[1;%2sm|", number );
    for ( i=0; i<nvec; i++ )
      printf0(" ---- depth: %d, c: %d, schwarz iter %2d, vector %d, norm: %11.6le |", l->depth, s->num_colors, k, i, r_norm[i] );
    printf0("\033[0m\n"); fflush(0);
    vector_PRECISION_free( &true_r, l, threading );
    END_LOCKED_MASTER(threading)
#endif
  }
  
  END_NO_HYPERTHREADS(threading)
}

/*****************  TEST ROUTINES  ***************************************************************************/

void schwarz_PRECISION_mvm_testfun_new( schwarz_PRECISION_struct *s, level_struct *l, struct Thread *threading ) {
  
  START_UNTHREADED_FUNCTION(threading)

  int mu, i, nb = s->num_blocks;
  int ivs = l->inner_vector_size;
  int n_vect = num_loop;//g.num_rhs_vect;

  void (*block_op)() = (l->depth==0)?block_d_plus_clover_PRECISION_new:coarse_block_operator_PRECISION_new;
  void (*op)() = (l->depth==0)?d_plus_clover_PRECISION_new:apply_coarse_operator_PRECISION_new;
  void (*boundary_op)() = (l->depth==0)?block_PRECISION_boundary_op_new:coarse_block_PRECISION_boundary_op_new;
  
  vector_PRECISION v1, v2, v3;
  PRECISION diff1[n_vect], diff2[n_vect];
  //printf("schwarz_PRECISION_mvm_testfun: %d %d %ld\n",n_vect,nb, l->inner_vector_size); 
  g.num_vect_pass1 = n_vect; g.num_vect_pass2 = n_vect;

  vector_PRECISION_init( &v1 );
  vector_PRECISION_init( &v2 );
  vector_PRECISION_init( &v3 );

  vector_PRECISION_alloc( &v1, _SCHWARZ, n_vect, l, no_threading );
  vector_PRECISION_alloc( &v2, _ORDINARY, n_vect, l, no_threading );
  vector_PRECISION_alloc( &v3, _ORDINARY, n_vect, l, no_threading );

  v1.num_vect_now = n_vect;
  v2.num_vect_now = n_vect;
  v3.num_vect_now = n_vect;
  vector_PRECISION_define_new(&v1,0,0,v1.size,l);
  //  printf("%p %p\n",v2.vector_buffer,v3.vector_buffer);
  vector_PRECISION_define_random_new( &v1, 0, ivs, l );
  //vector_PRECISION_define_new(&v2,0,0,l->vector_size,l);
  op( &v3, &v1, &(s->op), l, no_threading );
  //    printfv_PRECISION(&v3);
  //global_norm_PRECISION_new( diff1, &v3, 0, l->inner_vector_size, l, no_threading );for( i=0; i<n_vect; i++ ) printf("gnorm: %g\n",diff1[i]);
  for ( mu=0; mu<4; mu++ ) {
    ghost_update_PRECISION_new( &v1, mu, +1, &(s->op.c), l );
    ghost_update_PRECISION_new( &v1, mu, -1, &(s->op.c), l );
  }
      
  for ( mu=0; mu<4; mu++ ) {
    ghost_update_wait_PRECISION( &v1, mu, +1, &(s->op.c), l );
    ghost_update_wait_PRECISION( &v1, mu, -1, &(s->op.c), l );
  }
  //global_norm_PRECISION_new( diff2, &v2, 0, l->inner_vector_size, l, no_threading );for( i=0; i<n_vect; i++ ) printf("gnorm v2pre: %g\n",diff2[i]);
  for ( i=0; i<nb; i++ ) {
    block_op( &v2, &v1, s->block[i].start*l->num_lattice_site_var, s, l, no_threading );//  global_norm_PRECISION_new( diff2, &v2, 0, l->inner_vector_size, l, no_threading );for( i=0; i<n_vect; i++ ) printf("gnorm: %g\n",diff2[i]);
    boundary_op( &v2, &v1, i, s, l, no_threading );//  global_norm_PRECISION_new( diff2, &v2, 0, l->inner_vector_size, l, no_threading );for( i=0; i<n_vect; i++ ) printf("gnorm: %g\n",diff2[i]);
  }
  //printf("%p %p\n",v2.vector_buffer,v3.vector_buffer);

  //global_norm_PRECISION_new( diff2, &v2, 0, l->inner_vector_size, l, no_threading );for( i=0; i<n_vect; i++ ) printf("gnorm v2post: %g\n",diff2[i]);
  //    printfv_PRECISION(&v2);
  //if(l->depth!=0)for(i=0;i<l->inner_vector_size*n_vect;i++)printf0("%d %d %13g %g %g\n",i/n_vect,i%n_vect,creal_PRECISION(v2.vector_buffer[i]-v3.vector_buffer[i]),creal_PRECISION(v2.vector_buffer[i]),creal_PRECISION(v3.vector_buffer[i]));
  vector_PRECISION_minus_new( &v3, &v3, &v2, 0, l->inner_vector_size, l );//  printfv_PRECISION(&v3);  printfv_PRECISION(&v2); 
  global_norm_PRECISION_new( diff1, &v3, 0, l->inner_vector_size, l, no_threading );//d_plus_clover_PRECISION-block_PRECISION_boundary_op*(block_d_plus_clover_PRECISION)
  global_norm_PRECISION_new( diff2, &v2, 0, l->inner_vector_size, l, no_threading );//block_PRECISION_boundary_op*(block_d_plus_clover_PRECISION)
  /*for( i=0; i<n_vect; i++ ) printf("gnormv diff: %g\n",diff1[i]);
    for( i=0; i<n_vect; i++ ) printf("gnormv dinom v2: %g\n",diff2[i]);*/
  for( i=0; i<n_vect; i++ )
    test0_PRECISION("depth: %d, correctness of local residual vector: %le\n", l->depth, diff1[i]/diff2[i] );
  
  vector_PRECISION_free( &v1, l, no_threading );
  vector_PRECISION_free( &v2, l, no_threading );
  vector_PRECISION_free( &v3, l, no_threading );    

  END_UNTHREADED_FUNCTION(threading)
}
