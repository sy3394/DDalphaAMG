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
 */

#include "main.h"

void interpolation_PRECISION_struct_init( interpolation_PRECISION_struct *is ) {

  is->agg_index[T] = NULL;
  is->agg_boundary_index[T] = NULL;
  is->agg_boundary_neighbor[T] = NULL;
  is->operator = NULL;
  vector_PRECISION_init(&(is->test_vector_vec));
  vector_PRECISION_init(&(is->interpolation_vec));
  is->eigenvalues = NULL;
  //vector_PRECISION_init(&(is->tmp));
}


void coarsening_index_table_PRECISION_alloc( interpolation_PRECISION_struct *is, level_struct *l ) {
  
  int vol, sur, mu, agg_split[4], agg_size[4];
  
  is->num_agg = 1;
  for ( mu=0; mu<4; mu++ ) {
    agg_split[mu] = l->local_lattice[mu]/l->coarsening[mu]; // #aggregates in a local lattice
    agg_size[mu]  = l->coarsening[mu];                      // dims of each aggregate
    is->num_agg  *= agg_split[mu];                          // total #aggregates in a local lattice
  }
  vol = l->local_lattice[T]*l->local_lattice[Z]*l->local_lattice[Y]*l->local_lattice[X];

  for ( mu=0; mu<4; mu++ ) {
    sur = vol/l->local_lattice[mu];
    // number of lattice sites in local volume
    // that have a neighbor site in mu-direction which belongs to the same aggregate
    is->agg_length[mu] = (agg_size[mu]-1)*sur*agg_split[mu];
    // number of lattice sites in local volume
    // that have a neighbor site in mu-direction which belongs to a different aggregate
    is->agg_boundary_length[mu] = sur*agg_split[mu];
  }
  
  // index table for contributions to the self couplings of the coarse operator
  MALLOC( is->agg_index[T], int,
          is->agg_length[T]+is->agg_length[Z]+is->agg_length[Y]+is->agg_length[X] );
  // index table for contributions to neighbor couplings of the coarse operator
  MALLOC( is->agg_boundary_index[T], int,
          is->agg_boundary_length[T]+is->agg_boundary_length[Z]+is->agg_boundary_length[Y]+is->agg_boundary_length[X] );
  // corresponging neighbors of the sites in agg_boundary_index
  MALLOC( is->agg_boundary_neighbor[T], int,
          is->agg_boundary_length[T]+is->agg_boundary_length[Z]+is->agg_boundary_length[Y]+is->agg_boundary_length[X] );
  
  // offsets for the directions
  for ( mu=1; mu<4; mu++ ) {
    is->agg_index[mu] = is->agg_index[mu-1] + is->agg_length[mu-1];
    is->agg_boundary_index[mu] = is->agg_boundary_index[mu-1] + is->agg_boundary_length[mu-1];
    is->agg_boundary_neighbor[mu] = is->agg_boundary_neighbor[mu-1] + is->agg_boundary_length[mu-1];
  }
}


void coarsening_index_table_PRECISION_free( interpolation_PRECISION_struct *is, level_struct *l ) {
  
  int mu;
  
  FREE( is->agg_index[T], int,
        is->agg_length[T]+is->agg_length[Z]+is->agg_length[Y]+is->agg_length[X] );
  FREE( is->agg_boundary_index[T], int,
        is->agg_boundary_length[T]+is->agg_boundary_length[Z]+is->agg_boundary_length[Y]+is->agg_boundary_length[X] );
  FREE( is->agg_boundary_neighbor[T], int,
        is->agg_boundary_length[T]+is->agg_boundary_length[Z]+is->agg_boundary_length[Y]+is->agg_boundary_length[X] );
  
  for ( mu=1; mu<4; mu++ ) {
    is->agg_index[mu] = NULL;
    is->agg_boundary_index[mu] = NULL;
    is->agg_boundary_neighbor[mu] = NULL;
  }
}


void coarsening_index_table_PRECISION_define( interpolation_PRECISION_struct *is, schwarz_PRECISION_struct *s,
                                              level_struct *l ) {
  
  int i, j, k, mu, t, z, y, x, a0, b0, c0, d0, a1, b1, c1, d1, stride, offset,
      agg_split[4], block_split[4], block_size[4], agg_size[4],
      *index_table = s->op.index_table, *table_dim = s->op.table_dim,
      *count[4], *index_dir, *boundary_index_dir, *boundary_neighbor_index_dir,
      *neighbor = s->op.neighbor_table;
    
  for ( mu=0; mu<4; mu++ ) {
    agg_split[mu] = l->local_lattice[mu]/l->coarsening[mu];   // #aggregates
    agg_size[mu] = l->coarsening[mu];                         // dims of an aggregate
    block_split[mu] = l->coarsening[mu]/l->block_lattice[mu]; // #blocks in an aggregate
    block_size[mu] = l->block_lattice[mu];                    // dims of a block
  }

  stride = (l->depth==0)?4:5; offset = (l->depth==0)?0:1; // See define_nt_bt_tt in data_layout.c
  count[T]=&t; count[Z]=&z; count[Y]=&y; count[X]=&x;
  // filling index tables according to the schwarz operator layout
  for ( mu=0; mu<4; mu++ ) {
    
    i = 0; j = 0;
    index_dir = is->agg_index[mu];
    boundary_index_dir = is->agg_boundary_index[mu];
    boundary_neighbor_index_dir = is->agg_boundary_neighbor[mu];

    // visit each site in local lattice in the Schwarz order (index: iter_index -> Scwarz_index)
    //  i: counts the inner pos mu aggregate boundary sites in the order encountered
    //  j: counts the other sites
    // for each aggregate
    for ( d0=0; d0<agg_split[T]; d0++ )
      for ( c0=0; c0<agg_split[Z]; c0++ )
        for ( b0=0; b0<agg_split[Y]; b0++ )
          for ( a0=0; a0<agg_split[X]; a0++ )
            // for each block in the aggregate
            for ( d1=d0*block_split[T]; d1<(d0+1)*block_split[T]; d1++ )
              for ( c1=c0*block_split[Z]; c1<(c0+1)*block_split[Z]; c1++ )
                for ( b1=b0*block_split[Y]; b1<(b0+1)*block_split[Y]; b1++ )
                  for ( a1=a0*block_split[X]; a1<(a0+1)*block_split[X]; a1++ )
                    // for each site in the block
                    for ( t=d1*block_size[T]; t<(d1+1)*block_size[T]; t++ )
                      for ( z=c1*block_size[Z]; z<(c1+1)*block_size[Z]; z++ )
                        for ( y=b1*block_size[Y]; y<(b1+1)*block_size[Y]; y++ )
                          for ( x=a1*block_size[X]; x<(a1+1)*block_size[X]; x++ )
                            // if the neighbor in the pos mu dir in the adjacent aggregate
                            if ( (*(count[mu])+1) % agg_size[mu] != 0 ) {
			      // set is->agg_index[mu]: i counts these sites in Schwarz order (table_dim = l->local_lattice[mu]+2)
                              index_dir[i] = site_index( t, z, y, x, table_dim, index_table );
                              i++;
                            } else {
                              k = site_index( t, z, y, x, table_dim, index_table );
                              boundary_index_dir[j] = k;
                              boundary_neighbor_index_dir[j] = neighbor[ stride*k + mu + offset ];
                              j++;
                            }
  }
}
