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
 * matched:11/29/2019: need to update handling of different #vectors
 * changed from sbacchio
 * checked:12/05/2019:def of gs->buffer needs to be done properly; need clearer understanding of dist and gath (parent and offset too)
 */

#include "main.h"

void gathering_PRECISION_next_level_init( gathering_PRECISION_struct *gs, level_struct *l ) {
  /**************************************************************************************************************************
   * Assume: gs is one-level below than l, e.g., gs is at the next level if l is current level
   *         So fields in gs refers to the lattice at the level whose level struct containing this gs
   * gs->dist_local_lattice[mu]: dims of local lattice if none of the processes are turned off 
   * gs->dist_inner_lattice_sites: #sites of local lattice if none of the processes are turned off 
   * gs->gather_list_length: #processes in a group (parent and child at the next level)   
   **************************************************************************************************************************/

  int mu;
  
  gs->permutation = NULL;
  gs->gather_list = NULL;
  gs->reqs = NULL;
  vector_PRECISION_init(&(gs->buffer));
  vector_PRECISION_init(&(gs->transfer_buffer));
  
  /* 
   * Some processes can be turned off when going one-level down.
   * comm_offset[mu] specifies which process is idle in the mu dir in the Cartesian grid of processes.
   * In particular, it says: every comm_offset[mu]^th process in the mu dir is idle.  So
       (#processes at the given level) = (#process at the top level)/(product of comm_offset[i]) = (product of num_process_dir[i])
   * That is, if comm_offset[i] != 1, #active processes is reduced.
   * When this happens, when going one step down, an active process needs to handle more coarse sites than the ones
     obtained from the given process.
   * l->next_level->local_lattice[mu] is #sites in the mu dir that the active process needs to take care of at the next level.
   * So, one process (parent process) take care of calculation on sites of gs->gather_list_length many processes,
     all of which except the parent process will be idle at the next_level.
     So l->next_level->local_lattice[mu] > gs->dist_local_lattice[mu] b/c a process now take care of more sites */
  gs->dist_inner_lattice_sites = 1;
  gs->gather_list_length = 1;
  for ( mu=0; mu<4; mu++ ) {
    gs->dist_local_lattice[mu] = l->local_lattice[mu] / l->coarsening[mu]; // dims of local lattice at the next level if none of the processes are turned off
    gs->dist_inner_lattice_sites *= gs->dist_local_lattice[mu];            // #sites of local lattice at the next level if none of the processes are turned off 
    gs->gather_list_length *= l->next_level->local_lattice[mu] / gs->dist_local_lattice[mu]; // #processes in a group (parent and child at the next level)
  }
}

void gathering_PRECISION_setup( gathering_PRECISION_struct *gs, level_struct *l ) {
  
  // define data merging
  // define data gathering permutation
  int i, mu, current_rank, offset, offset_sum, nvect = num_loop;//(g.num_rhs_vect < l->num_eig_vect)? l->num_eig_vect:g.num_rhs_vect;//!!!!!!!!!
  int process_coords[4] = {0,0,0,0}, parent_coords[4] = {0,0,0,0}, *process_list = NULL;
  // for sites in the child processes
  MALLOC( process_list, int, l->num_processes );
#ifdef HAVE_TM1p1
  MALLOC( gs->transfer_buffer.vector_buffer, complex_PRECISION, 2 * gs->dist_inner_lattice_sites * l->num_lattice_site_var * nvect );
#else
  MALLOC( gs->transfer_buffer.vector_buffer, complex_PRECISION, gs->dist_inner_lattice_sites * l->num_lattice_site_var * nvect );
#endif  
  gs->transfer_buffer.num_vect = nvect;
  gs->transfer_buffer.size = gs->dist_inner_lattice_sites * l->num_lattice_site_var;

  l->idle = 0;
  i = 0;
  // go over each process both active and idle
  for ( process_coords[T]=0; process_coords[T]<g.process_grid[T]; process_coords[T]++ )
    for ( process_coords[Z]=0; process_coords[Z]<g.process_grid[Z]; process_coords[Z]++ )
      for ( process_coords[Y]=0; process_coords[Y]<g.process_grid[Y]; process_coords[Y]++ )
        for ( process_coords[X]=0; process_coords[X]<g.process_grid[X]; process_coords[X]++ ) {

          g.Cart_rank( g.comm_cart, process_coords, &current_rank );

	  // find out if the current process is a parent rank or its child
          offset_sum = 0;
          for ( mu=0; mu<4; mu++ ) {
            offset = process_coords[mu] % l->comm_offset[mu]; // processes grouped into l->comm_offset[mu] processes in the mu dir
            parent_coords[mu] = process_coords[mu] - offset;  // and the first one in the group becomes the parent
            offset_sum+=offset;
          }
          // if the current rank is a child process, put the current rank idle
	  //   if #process in dir has not changed when going deeper, i.e., l->comm_offset==1(default), each rank is its own parent rank
          if ( current_rank == g.my_rank ) {
            g.Cart_rank( g.comm_cart, parent_coords, &(l->parent_rank) );
            // find out if current process is supposed to idle
            if ( offset_sum > 0 )
              l->idle = 1;
          }
          // store current rank in case the respective process is not supposed to idle
          if ( offset_sum == 0 ) {
            process_list[i] = current_rank;
            i++;
          }
        }
        
  // creating sub communicator for inner products
  MPI_Comm_group( g.comm_cart, &(gs->level_comm_group) );
  MPI_Group_incl( g.global_comm_group, l->num_processes, process_list, &(gs->level_comm_group) );
  MPI_Comm_create( g.comm_cart, gs->level_comm_group, &(gs->level_comm) );
  FREE( process_list, int, l->num_processes );
  
  if ( !l->idle ) { // only parent performs the folowing setup
    int d0, c0, b0, a0, d1, c1, b1, a1, t, z, y, x, k, j, *field1=NULL, *field2=NULL, *count[4],
    merge[4], block_size[4], block_split[4], agg_split[4];
    MALLOC( gs->gather_list, int, gs->gather_list_length );
    MALLOC( gs->permutation, int, l->num_inner_lattice_sites );
    MALLOC( gs->reqs, MPI_Request, gs->gather_list_length );
#ifdef HAVE_TM1p1
    vector_PRECISION_alloc( &(gs->buffer), _INNER, 2*nvect, l, no_threading );
#else
    vector_PRECISION_alloc( &(gs->buffer), _INNER, 1*nvect, l, no_threading );
#endif
    gs->buffer.num_vect = nvect;

    //---- define gs->permutation, gs->gather_list
    MALLOC( field1, int, l->num_inner_lattice_sites );
    MALLOC( field2, int, l->num_inner_lattice_sites );
    
    for ( mu=0; mu<4; mu++ ) { // processes are grouped into a family of the parent rank and its children
      block_size[mu] = gs->dist_local_lattice[mu];       // dims of local lattice on each process in the family
      merge[mu] = l->local_lattice[mu] / block_size[mu]; // #processes of this family in the mu dir
    }

    // visit each site on a rank (a block in the new local lattice) in the family in the order x->y->z->t
    // field1: local_lex_index of each site on the new local lattice for the parent -> local_lex_index of each site on each rank in the family
    count[0]=&d0; count[1]=&c0; count[2]=&b0; count[3]=&a0;
    i=0; j=0;
    // for each process in the family (parent and its children)
    for ( d0=0; d0<merge[T]; d0++)
      for ( c0=0; c0<merge[Z]; c0++)
        for ( b0=0; b0<merge[Y]; b0++ )
          for ( a0=0; a0<merge[X]; a0++ ) {
            // find the rank ID for each child rank in the family
	    //   Active ranks on the level one up are grouped into families of a parent and children.
	    //   So to find a global coordinate of child ranks, we need comm_offset for the level one up.
	    //   This can be found by l->comm_offset[mu]/merge[mu].
            for ( mu=0; mu<4; mu++ )
              process_coords[mu] = g.my_coords[mu] + *(count[mu]) * (l->comm_offset[mu]/merge[mu]);
            g.Cart_rank( g.comm_cart, process_coords, gs->gather_list + j );
            j++;
	    
            // for each site on the process
            for ( t=d0*block_size[T]; t<(d0+1)*block_size[T]; t++ )
              for ( z=c0*block_size[Z]; z<(c0+1)*block_size[Z]; z++ )
                for ( y=b0*block_size[Y]; y<(b0+1)*block_size[Y]; y++ )
                  for ( x=a0*block_size[X]; x<(a0+1)*block_size[X]; x++ ) {
                    k = lex_index( t, z, y, x, l->local_lattice );
                    field1[k] = i;
                    i++;
                  }
          }
    

    // visit each site in the Schwarz order (aggregate->block->lex_site) if l->level > 0
    // otherwise, visit each site on the new local lattie lexicographically (even->odd if g.odd_even)
    // field2: local_lex_index -> Scwarz_index (if l->level != 0) and local_lex_index (if l->level == 0)
    if ( l->level > 0 ) {
      
      for ( mu=0; mu<4; mu++ ) {
        agg_split[mu] = l->local_lattice[mu]/l->coarsening[mu];   // #aggregates in the local lattice in the mu dir
        block_split[mu] = l->coarsening[mu]/l->block_lattice[mu]; // #blocks in the block in the mu dir
        block_size[mu] = l->block_lattice[mu];                    // dims of blocks
      }
      
      i = 0;
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
                            for ( x=a1*block_size[X]; x<(a1+1)*block_size[X]; x++ ) {
                              k = lex_index( t, z, y, x, l->local_lattice );
                              field2[k] = i;
                              i++;
                            }
    } else {
      if ( g.odd_even ) {
        int oe_offset=0, *le = l->local_lattice;
        
        for ( mu=0; mu<4; mu++ )
          oe_offset += (le[mu]*(g.my_coords[mu]/l->comm_offset[mu]))%2;
        oe_offset = oe_offset%2;
        
        i=0;
	// for each even sites in the local lattice
        for ( t=0; t<le[T]; t++ )
          for ( z=0; z<le[Z]; z++ )
            for ( y=0; y<le[Y]; y++ )
              for ( x=0; x<le[X]; x++ )
                if ( (t+z+y+x+oe_offset)%2 == 0 ) {
                  field2[ lex_index( t, z, y, x, le ) ] = i;
                  i++;
                }
	// for each odd sites in the local lattice  
        for ( t=0; t<le[T]; t++ )
          for ( z=0; z<le[Z]; z++ )
            for ( y=0; y<le[Y]; y++ )
              for ( x=0; x<le[X]; x++ )
                if ( (t+z+y+x+oe_offset)%2 == 1 ) {
                  field2[ lex_index( t, z, y, x, le ) ] = i;
                  i++;
                }
                        
      } else {
	// for each sites in the local lattice  
        for ( i=0; i<l->num_inner_lattice_sites; i++ )
          field2[i] = i;
      }
    }      
    
    // build the permutation table process-wise_lex_index to Scwarz_index/lex_index (level>0/level==0)
    for ( i=0; i<l->num_inner_lattice_sites; i++ )
      gs->permutation[ field1[i] ] = field2[i];
    
    FREE( field1, int, l->num_inner_lattice_sites );
    FREE( field2, int, l->num_inner_lattice_sites );
  }  
}


void gathering_PRECISION_free( gathering_PRECISION_struct *gs, level_struct *l ) {
  
  int nvec = gs->transfer_buffer.num_vect;
#ifdef HAVE_TM1p1
  nvec *= 2;
#endif

  if ( !l->idle ) {
    FREE( gs->gather_list, int, gs->gather_list_length );
    FREE( gs->permutation, int, l->num_inner_lattice_sites );
    FREE( gs->reqs, MPI_Request, gs->gather_list_length );
    vector_PRECISION_free( &(gs->buffer), l, no_threading );
  }
  
  MPI_Comm_free( &(gs->level_comm) );
  MPI_Group_free( &(gs->level_comm_group) );

  FREE( gs->transfer_buffer.vector_buffer, complex_PRECISION, gs->dist_inner_lattice_sites * l->num_lattice_site_var * nvec );

}

void conf_PRECISION_gather( operator_PRECISION_struct *out, operator_PRECISION_struct *in, level_struct *l ) {
  /*
   * Descrption: copies in into out while taking care of possible redunction of #processes
   */
  
  int send_size_hopp  = l->gs_PRECISION.dist_inner_lattice_sites * 4 * SQUARE( l->num_lattice_site_var );
  int send_size_clov  = l->gs_PRECISION.dist_inner_lattice_sites * ( (l->num_lattice_site_var*(l->num_lattice_site_var+1))/2 );
  int send_size_block = l->gs_PRECISION.dist_inner_lattice_sites * ( (l->num_lattice_site_var/2*(l->num_lattice_site_var/2+1)) );
  int jj, jjj;;
#ifdef HAVE_TM
  out->mu            = in->mu;
  for ( int i=0; i<g.num_rhs_vect; i++ ) out->mu_even_shift[i] = in->mu_even_shift[i];
  out->mu_odd_shift  = in->mu_odd_shift;
  out->odd_shifted_mu = in->odd_shifted_mu;
  for ( int i=0; i<g.num_rhs_vect; i++ ) out->diff_mu_eo[i] = in->diff_mu_eo[i];
  out->even_shift_avg = in->even_shift_avg;
  out->is_even_shifted_mu_nonzero = in->is_even_shifted_mu_nonzero;
#endif
  out->m0 = in->m0;
#ifdef HAVE_TM1p1
  out->epsbar                = in->epsbar;
  out->epsbar_ig5_even_shift = in->epsbar_ig5_even_shift;
  out->epsbar_ig5_odd_shift  = in->epsbar_ig5_odd_shift;
#endif

  if ( g.my_rank != l->parent_rank ) {
    MPI_Request req, odd_req;
#ifdef HAVE_TM1p1
    MPI_Request eps_req;
    MPI_Isend( in->epsbar_term, send_size_block, MPI_COMPLEX_PRECISION, l->parent_rank, 4, g.comm_cart, &eps_req );
#endif
#ifdef HAVE_MULT_TM
    int nrt = num_loop; //g.num_rhs_vect;
    MPI_Request tm_req;
    buffer_PRECISION tm_send_buffer;
    MALLOC( tm_send_buffer, complex_PRECISION, nrt*send_size_block );
    for ( int i=0; i< send_size_block; i++ )
      VECTOR_LOOP(jj, nrt, jjj, tm_send_buffer[nrt*i+jj+jjj] = in->tm_term[i*num_loop+jj*send_size_block+jjj]) 
    MPI_Isend( tm_send_buffer, nrt*send_size_block, MPI_COMPLEX_PRECISION, l->parent_rank, 3, g.comm_cart, &tm_req );
#endif
    MPI_Isend( in->odd_proj, send_size_block, MPI_COMPLEX_PRECISION, l->parent_rank, 2, g.comm_cart, &odd_req );
    MPI_Isend( in->D, send_size_hopp, MPI_COMPLEX_PRECISION, l->parent_rank, 0, g.comm_cart, &req );
    MPI_Send( in->clover, send_size_clov, MPI_COMPLEX_PRECISION, l->parent_rank, 1, g.comm_cart );
#ifdef HAVE_TM1p1
    MPI_Wait( &eps_req, MPI_STATUS_IGNORE );
#endif
#ifdef HAVE_MULT_TM
    MPI_Wait( &tm_req, MPI_STATUS_IGNORE );
    FREE( tm_send_buffer, complex_PRECISION, nrt*send_size_block );
#endif    
    MPI_Wait( &odd_req, MPI_STATUS_IGNORE );
    MPI_Wait( &req, MPI_STATUS_IGNORE );
  } else {
    int i, j, n=l->gs_PRECISION.gather_list_length, s=l->num_inner_lattice_sites,
      t, *pi = l->gs_PRECISION.permutation;
    buffer_PRECISION buffer_hopp = NULL, buffer_clov = NULL, buffer_odd_proj = NULL;
    MPI_Request *hopp_reqs = NULL, *clov_reqs = NULL, *odd_proj_reqs = NULL;

#ifdef HAVE_TM1p1
    buffer_PRECISION buffer_eps_term = NULL;
    MPI_Request *eps_term_reqs = NULL;
    MALLOC( buffer_eps_term, complex_PRECISION, n*send_size_block );
    MALLOC( eps_term_reqs, MPI_Request, n );
#endif
#ifdef HAVE_MULT_TM
    int nrt = num_loop;
    buffer_PRECISION buffer_tm_term = NULL;
    MPI_Request *tm_term_reqs = NULL;
    MALLOC( buffer_tm_term, complex_PRECISION, n*nrt*send_size_block );
    MALLOC( tm_term_reqs, MPI_Request, n );
#endif
    MALLOC( buffer_hopp, complex_PRECISION, n*send_size_hopp );
    MALLOC( buffer_clov, complex_PRECISION, n*send_size_clov );
    MALLOC( buffer_odd_proj, complex_PRECISION, n*send_size_block );
    MALLOC( hopp_reqs, MPI_Request, n );
    MALLOC( clov_reqs, MPI_Request, n );
    MALLOC( odd_proj_reqs, MPI_Request, n );
    
    PROF_PRECISION_START( _GD_COMM );
    for ( i=1; i<n; i++ ) {
#ifdef HAVE_TM1p1
      MPI_Irecv( buffer_eps_term+i*send_size_block, send_size_block, MPI_COMPLEX_PRECISION,
                 l->gs_PRECISION.gather_list[i], 4, g.comm_cart, &(eps_term_reqs[i]) );
#endif
#ifdef HAVE_MULT_TM
      MPI_Irecv( buffer_tm_term+i*nrt*send_size_block, nrt*send_size_block, MPI_COMPLEX_PRECISION,
		 l->gs_PRECISION.gather_list[i], 3, g.comm_cart, &(tm_term_reqs[i]) );
#endif
      MPI_Irecv( buffer_hopp+i*send_size_hopp, send_size_hopp, MPI_COMPLEX_PRECISION,
                 l->gs_PRECISION.gather_list[i], 0, g.comm_cart, &(hopp_reqs[i]) );
      MPI_Irecv( buffer_clov+i*send_size_clov, send_size_clov, MPI_COMPLEX_PRECISION,
                 l->gs_PRECISION.gather_list[i], 1, g.comm_cart, &(clov_reqs[i]) );
      MPI_Irecv( buffer_odd_proj+i*send_size_block, send_size_block, MPI_COMPLEX_PRECISION,
                 l->gs_PRECISION.gather_list[i], 2, g.comm_cart, &(odd_proj_reqs[i]) );      
    }
    PROF_PRECISION_STOP( _GD_COMM, 2*n-2 );

#ifdef HAVE_TM1p1
    for ( i=0; i<send_size_block; i++ )
      buffer_eps_term[i] = in->epsbar_term[i];
#endif
#ifdef HAVE_MULT_TM
    for ( i=0; i<send_size_block; i++ )
      VECTOR_LOOP(jj, nrt, jjj, buffer_tm_term[i*nrt+jj+jjj] = in->tm_term[i*num_loop+send_size_block*jj+jjj];)
#endif
    
    for ( i=0; i<send_size_hopp; i++ )
      buffer_hopp[i] = in->D[i];
    
    for ( i=0; i<send_size_clov; i++ )
      buffer_clov[i] = in->clover[i];
    
    for ( i=0; i<send_size_block; i++ )
      buffer_odd_proj[i] = in->odd_proj[i];

#ifdef HAVE_TM1p1
    PROF_PRECISION_START( _GD_IDLE );
    for ( i=1; i<n; i++ )
      MPI_Wait( &(eps_term_reqs[i]), MPI_STATUS_IGNORE );
    PROF_PRECISION_STOP( _GD_IDLE, n-1 );
    
    t = (send_size_block*n)/s;
    for ( i=0; i<s; i++ )
      for ( j=0; j<t; j++ )
        out->epsbar_term[ t*pi[i] + j ] = buffer_eps_term[ t*i + j ];
#endif
#ifdef HAVE_MULT_TM
    PROF_PRECISION_START( _GD_IDLE );
    for ( i=1; i<n; i++ )
      MPI_Wait( &(tm_term_reqs[i]), MPI_STATUS_IGNORE );
    PROF_PRECISION_STOP( _GD_IDLE, n-1 );

    t = (send_size_block*n)/s;//==internal d.o.f. = l->num_lattice_site_var/2*(l->num_lattice_site_var/2+1)
    for ( i=0; i<s; i++ )
      for ( j=0; j<t; j++ )
	VECTOR_LOOP(jj, nrt, jjj, out->tm_term[ (t*pi[i] + j)*num_loop+send_size_block*jj+jjj ] = buffer_tm_term[ (t*i + j)*nrt+jj+jjj ];)
#endif

    PROF_PRECISION_START( _GD_IDLE );
    for ( i=1; i<n; i++ )
      MPI_Wait( &(hopp_reqs[i]), MPI_STATUS_IGNORE );
    PROF_PRECISION_STOP( _GD_IDLE, n-1 );
    
    t = (send_size_hopp*n)/s;
    for ( i=0; i<s; i++ )
      for ( j=0; j<t; j++ )
        out->D[ t*pi[i] + j ] = buffer_hopp[ t*i + j ];
    
    PROF_PRECISION_START( _GD_IDLE );
    for ( i=1; i<n; i++ )
      MPI_Wait( &(clov_reqs[i]), MPI_STATUS_IGNORE );
    PROF_PRECISION_STOP( _GD_IDLE, n-1 );
    
    t = (send_size_clov*n)/s;
    for ( i=0; i<s; i++ )
      for ( j=0; j<t; j++ )
        out->clover[ t*pi[i] + j ] = buffer_clov[ t*i + j ];

    PROF_PRECISION_START( _GD_IDLE );
    for ( i=1; i<n; i++ )
      MPI_Wait( &(odd_proj_reqs[i]), MPI_STATUS_IGNORE );
    PROF_PRECISION_STOP( _GD_IDLE, n-1 );
    
    t = (send_size_block*n)/s;
    for ( i=0; i<s; i++ )
      for ( j=0; j<t; j++ )
        out->odd_proj[ t*pi[i] + j ] = buffer_odd_proj[ t*i + j ];
      
    FREE( buffer_hopp, complex_PRECISION, n*send_size_hopp );
    FREE( buffer_clov, complex_PRECISION, n*send_size_clov );
    FREE( buffer_odd_proj, complex_PRECISION, n*send_size_block );
    FREE( hopp_reqs, MPI_Request, n );
    FREE( clov_reqs, MPI_Request, n );
    FREE( odd_proj_reqs, MPI_Request, n );
#ifdef HAVE_MULT_TM
    FREE( buffer_tm_term, complex_PRECISION, n*nrt*send_size_block );
    FREE( tm_term_reqs, MPI_Request, n );
#endif
#ifdef HAVE_TM1p1
    FREE( buffer_eps_term, complex_PRECISION, n*send_size_block );
    FREE( eps_term_reqs, MPI_Request, n );
#endif

  }
  l->dummy_p_PRECISION.op = out;
  l->dummy_p_PRECISION.v_start = 0;
  l->dummy_p_PRECISION.v_end = l->inner_vector_size;
  l->dummy_p_PRECISION.eval_operator = apply_coarse_operator_PRECISION;
}

// when #process decreases in going deeper, gather entries needed to do restriction from other ranks to the parent, and parent hold phi_c entries
// otherwise, simply gath <- dist after reordering
void vector_PRECISION_gather( vector_PRECISION *gath, vector_PRECISION *dist, level_struct *l ) {
  /*************************
   * Description: gather values on the sites also belongin to the child processes to the parent process
   * Input:
   *   vector_PRECISION *dist: local vector to be sent to the parent rank
   * Output:
   *   vector_PRECISION *gath: a new local vector consisting of entries including the ones gathered from child/idle ranks
   * note: need to send all vectors in dist irrespective of how many are currently used as the data are consecutive
   ************************/

  int nvect = dist->num_vect_now, nvect_dist = dist->num_vect, nvect_gath = gath->num_vect;
  int send_size = l->gs_PRECISION.dist_inner_lattice_sites * l->num_lattice_site_var;

  if ( l->gs_PRECISION.buffer.num_vect < nvect_dist )
    error0("vector_PRECISION_gather: potential memory overflow\n");

  if ( g.my_rank != l->parent_rank ) {
    MPI_Send( dist->vector_buffer, send_size*nvect_dist, MPI_COMPLEX_PRECISION, l->parent_rank, g.my_rank, g.comm_cart );
  } else {
    int i, j, jj, jjj;
    int n = l->gs_PRECISION.gather_list_length;
    int s = l->num_inner_lattice_sites;
    int t = l->num_lattice_site_var;
    int *pi = l->gs_PRECISION.permutation;
    vector_PRECISION buffer = l->gs_PRECISION.buffer;

    PROF_PRECISION_START( _GD_COMM );
    for ( i=1; i<n; i++ )//store recieved data from the second chank gathered from child processes
      MPI_Irecv( buffer.vector_buffer+i*send_size*nvect_dist, send_size*nvect_dist, MPI_COMPLEX_PRECISION, l->gs_PRECISION.gather_list[i],
                 l->gs_PRECISION.gather_list[i], g.comm_cart, &(l->gs_PRECISION.reqs[i]) );
    PROF_PRECISION_STOP( _GD_COMM, n-1 );

    for ( i=0; i<send_size; i++ )//copy the data at the parent rank into the first chank of buffer
      VECTOR_LOOP(jj, nvect_dist, jjj, buffer.vector_buffer[i*nvect_dist+jj+jjj] = dist->vector_buffer[i*nvect_dist+jj+jjj];)

    PROF_PRECISION_START( _GD_IDLE );
    for ( i=1; i<n; i++ )
      MPI_Wait( &(l->gs_PRECISION.reqs[i]), MPI_STATUS_IGNORE );
    PROF_PRECISION_STOP( _GD_IDLE, n-1 );
    // permute data according to desired data layout for parent process
    for ( i=0; i<s; i++ )
      for ( j=0; j<t; j++ )
        VECTOR_LOOP(jj, nvect, jjj, gath->vector_buffer[ (t*pi[i] + j)*nvect_gath+jj+jjj ] = buffer.vector_buffer[ (t*i + j)*nvect_dist+jj+jjj ];)
  }
}

void vector_PRECISION_distribute( vector_PRECISION *dist, vector_PRECISION *gath, level_struct *l ) {
    /*************************
     * Description: send values on the sites also belongin to the child processes from the parent process
     * Input:
     *   vector_PRECISION *gath: vector on the parent rank to be send around child/idle ranks
     * Output:
     *   vector_PRECISION *dist: updated local vector distributed from parent rank
     * note: need to send all vectors in dist irrespective of how many are currently used as the data are consecutive
     ************************/

  int nvec = dist->num_vect_now, nvec_dist = dist->num_vect, nvec_gath = gath->num_vect;//nvec:tmp fix as dist contains MIN!!!!!!!
  int send_size = l->gs_PRECISION.dist_inner_lattice_sites * l->num_lattice_site_var;
  
  if ( l->gs_PRECISION.buffer.num_vect < nvec_dist || nvec_dist < nvec )//???????
    error0("vector_PRECISION_gather: potential memory overflow\n");

  if ( g.my_rank != l->parent_rank ) {// if I am a child 
    MPI_Recv( dist->vector_buffer, send_size*nvec_dist, MPI_COMPLEX_PRECISION, l->parent_rank, g.my_rank, g.comm_cart, MPI_STATUS_IGNORE );//!!!!
  } else { //if I am a parent
    int i, j, jj, jjj;
    int n_p =l->gs_PRECISION.gather_list_length;
    int n_s =l->num_inner_lattice_sites;
    int n_i =l->num_lattice_site_var;
    int *pi = l->gs_PRECISION.permutation;
    vector_PRECISION buffer = l->gs_PRECISION.buffer;
    // permute data according to desired distributed data layout and store them in buffer
    for ( i=0; i<n_s; i++ )
      for ( j=0; j<n_i; j++ )
        VECTOR_LOOP(jj, nvec, jjj, buffer.vector_buffer[ (n_i*i+j)*nvec_dist+jj+jjj ] = gath->vector_buffer[ (n_i*pi[i]+j)*nvec_gath+jj+jjj ];)//!!!!
    
    // send buffer from the second chank
    PROF_PRECISION_START( _GD_COMM );
    for ( i=1; i<n_p; i++ ) // n = Prod_mu l->next_level->local_lattice[mu]???
      MPI_Isend( buffer.vector_buffer+i*send_size*nvec_dist, send_size*nvec_dist, MPI_COMPLEX_PRECISION, l->gs_PRECISION.gather_list[i], //!!!!
                 l->gs_PRECISION.gather_list[i], g.comm_cart, &(l->gs_PRECISION.reqs[i]) );
    PROF_PRECISION_STOP( _GD_COMM, n_p-1 );
    // reflect reordering in the local chank
    for ( i=0; i<send_size; i++ )
      VECTOR_LOOP(jj, nvec, jjj, dist->vector_buffer[i*nvec_dist+jj+jjj] = buffer.vector_buffer[i*nvec_dist+jj+jjj];)//!!!
    
    PROF_PRECISION_START( _GD_IDLE );
    for ( i=1; i<n_p; i++ )
      MPI_Wait( &(l->gs_PRECISION.reqs[i]), MPI_STATUS_IGNORE );
    PROF_PRECISION_STOP( _GD_IDLE, n_p-1 );
  }  
}
