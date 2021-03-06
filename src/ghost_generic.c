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
 */

#include "main.h"


void ghost_alloc_PRECISION( int buffer_size, comm_PRECISION_struct *c, level_struct *l ) {
  
  int mu, nu, factor=1, n_vect = num_loop;

  // offset is internal d.o.f. communicated.
  if ( l->depth > 0 ) {
    // full internal d.o.f. is transferred
    c->offset = l->num_lattice_site_var;
  } else {
    // Only, half spinor (01 components) is transferred at the top, as the other half is projected out
    c->offset = l->num_lattice_site_var/2;
    if ( g.method < 5 )
      factor = 2; // In ghost_update_PRECISION for smoother, we need full internal d.o.f.
  }

#ifdef HAVE_TM1p1
  n_vect *= 2;
#endif

  //--- set c->comm_start, c->length, and c->max_length and allocate c->buffer
  if ( buffer_size <= 0 ) { // default
    c->comm_start[0] = c->offset*l->num_inner_lattice_sites;
    c->comm_start[1] = c->offset*l->num_inner_lattice_sites;
    for ( mu=0; mu<4; mu++ ) {
      if ( mu > 0 ) {
        c->comm_start[2*mu]   = c->comm_start[2*(mu-1)]   + buffer_size;
        c->comm_start[2*mu+1] = c->comm_start[2*(mu-1)+1] + buffer_size;
      }
      buffer_size = c->offset;
      for ( nu=0; nu<4; nu++ ) {
        if ( nu != mu ) {
          buffer_size *= l->local_lattice[nu];
        }
      }
      c->length[2*mu]   = buffer_size;
      c->length[2*mu+1] = buffer_size;
      c->max_length[mu] = factor*buffer_size;
      MALLOC( c->buffer[2*mu],   complex_PRECISION, factor*buffer_size*n_vect );//printf("ghot aloc(%d) %d: %d %d\n",l->depth,mu,factor*buffer_size,n_vect);
      MALLOC( c->buffer[2*mu+1], complex_PRECISION, factor*buffer_size*n_vect );
      c->in_use[2*mu]   = 0;
      c->in_use[2*mu+1] = 0;
    }
  } else { // used when reallocating the below fields in case of memory shortage in ghost_sendrecv_PRECISION
    for ( mu=0; mu<4; mu++ ) {
      c->max_length[mu] = buffer_size;
      MALLOC( c->buffer[2*mu],   complex_PRECISION, buffer_size*n_vect );
      MALLOC( c->buffer[2*mu+1], complex_PRECISION, buffer_size*n_vect );
    }
  }

  //--- allocate l->vbuf_PRECISION[4]
  if ( l->vbuf_PRECISION[4].vector_buffer == NULL ) {
    // no_threading is not allocated yet
    //vector_PRECISION_alloc( &(l->vbuf_PRECISION[4]), _ORDINARY, n_vect, l, no_threading);
    MALLOC( l->vbuf_PRECISION[4].vector_buffer, complex_PRECISION, l->vector_size*n_vect );
    l->vbuf_PRECISION[4].type = _ORDINARY;
    l->vbuf_PRECISION[4].size = l->vector_size;
    l->vbuf_PRECISION[4].num_vect = n_vect;
    l->vbuf_PRECISION[4].num_vect_now = n_vect;
    l->vbuf_PRECISION[4].l = l;
  }
  c->num_vect = n_vect;
}

void ghost_free_PRECISION( comm_PRECISION_struct *c, level_struct *l ) {
  
  int mu, n_vect = c->num_vect;
  
  for ( mu=0; mu<4; mu++ ) {    
    FREE( c->buffer[2*mu],   complex_PRECISION, c->max_length[mu]*n_vect );
    FREE( c->buffer[2*mu+1], complex_PRECISION, c->max_length[mu]*n_vect );
  }
  if ( l->vbuf_PRECISION[4].vector_buffer != NULL ){
    //vector_PRECISION_free( &(l->vbuf_PRECISION[4]), l, no_threading);
    FREE(l->vbuf_PRECISION[4].vector_buffer, complex_PRECISION, l->vector_size*n_vect);
   }
}

void ghost_sendrecv_init_PRECISION( const int type, comm_PRECISION_struct *c, level_struct *l ) {
  
  int mu; 
  
  if ( type == _COARSE_GLOBAL ) {
    c->comm = 1;
    for ( mu=0; mu<4; mu++ ) {
      ASSERT( c->in_use[2*mu] == 0 );
      ASSERT( c->in_use[2*mu+1] == 0 );
    }
  }
}

void ghost_sendrecv_PRECISION( buffer_PRECISION phi, const int mu, const int dir,
			       comm_PRECISION_struct *c, const int amount, level_struct *l ) {
  /**************************
   * buffer_PRECISION phi: vectors to be sent
   * mu: specify boundary
   * dir: eiter postive (update ghose shell) or negative (update inner boundary sites)
   *   if +1: send sites of phi on pos mu boundary ghost shell to the rank in the pos mu dir and store data from the rank in the neg mu dir in buffer
   *   if -1: send neg mu inner boudnary sites to the rank in the neg mu dir and update pos mu ghost shell
   * Note: does not allow sending in both directions at the same time
   *************************/

  if( l->global_splitting[mu] > 1 ) {

    int i, j, jj, jjj, *table=NULL, mu_dir = 2*mu-MIN(dir,0), offset = c->offset, length[2] = {0,0}, comm_start = 0, table_start = 0;
    int nvec = g.num_vect_pass1, nvec_phi = g.num_vect_pass2, nvec_com=c->num_vect;//temp fix!!!!!!!!
    buffer_PRECISION buffer, phi_pt;

    // offset = (internal d.o.f.)/2 if at the top; offset = (internal d.o.f.) if not.
    if ( amount == _FULL_SYSTEM ) {
      length[0]   = (c->num_boundary_sites[2*mu])*offset;   // #inner boundary sites on the positive mu dir x offset
      length[1]   = (c->num_boundary_sites[2*mu+1])*offset; // #inner boundary sites on the negative mu dir x offset
      comm_start  = c->comm_start[mu_dir];
      table_start = 0;
    } else if ( amount == _EVEN_SITES ) {
      length[0]   = c->num_even_boundary_sites[2*mu]*offset;   // #inner boundary even sites on the positive mu dir x offset
      length[1]   = c->num_even_boundary_sites[2*mu+1]*offset; // #inner boundary even sites on the negative mu dir x offset 
      comm_start  = c->comm_start[mu_dir];
      table_start = 0;
    } else if ( amount == _ODD_SITES ) {
      length[0]   = c->num_odd_boundary_sites[2*mu]*offset;    // #inner boundary odd sites on the positive mu dir x offset
      length[1]   = c->num_odd_boundary_sites[2*mu+1]*offset;  // #inner boundary odd sites on the negative mu dir x offset
      comm_start  = c->comm_start[mu_dir]+c->num_even_boundary_sites[mu_dir]*offset;
      table_start = c->num_even_boundary_sites[mu_dir];
    }

    ASSERT( c->in_use[mu_dir] == 0 );
    c->in_use[mu_dir] = 1;

    if ( MAX(length[0],length[1]) > c->max_length[mu] ) {
      printf("CAUTION: my_rank: %d, not enough comm buffer\n", g.my_rank ); fflush(0);
      ghost_free_PRECISION( c, l );
      ghost_alloc_PRECISION( MAX(length[0],length[1]), c, l );
    }
    
    buffer = (buffer_PRECISION)c->buffer[mu_dir];
    // dir = senddir
    if ( dir == 1 ) {//printf0("send 1: %d, %d %d\n",c->max_length[mu],length[1],nvec_phi);fflush(stdout);
      // data to be communicated is stored serially in the vector phi
      // recv target is a buffer
      // afterwards (in ghost_wait) the data has to be distributed onto the correct sites
      // touching the respective boundary in -mu direction
      
      phi_pt = phi + comm_start*nvec_phi;
      if ( length[1] > 0 ) {
	// get ghost cells on the negative mu boundary from the rank in the negative mu dir
        PROF_PRECISION_START( _OP_COMM );
        MPI_Irecv( buffer, length[1]*nvec_phi, MPI_COMPLEX_PRECISION,//why use buffer in this case only????
                   l->neighbor_rank[2*mu+1], 2*mu, g.comm_cart, &(c->rreqs[2*mu]) );
        PROF_PRECISION_STOP( _OP_COMM, 1 );
      }
      if ( length[0] > 0 ) {
	// send inner boundary sites on the positive mu boundary to the rank in the positive mu dir
        PROF_PRECISION_START( _OP_COMM );
        MPI_Isend( phi_pt, length[0]*nvec_phi, MPI_COMPLEX_PRECISION,
                   l->neighbor_rank[2*mu], 2*mu, g.comm_cart, &(c->sreqs[2*mu]) );
        PROF_PRECISION_STOP( _OP_COMM, 0 );
      }
      
      
    } else if ( dir == -1 ) {
      // data to be communicated is stored on the sites touching the boundary in -mu direction
      // this data is gathered in a buffer in the correct ordering
      // which is required on the boundary of the vector phi
      int num_boundary_sites = length[1]/offset;
      
      table = c->boundary_table[2*mu+1]+table_start;
      for ( j=0; j<num_boundary_sites; j++ ) {// copy each neg mu inner boundary sites onto buffer
        phi_pt = phi + table[j]*offset*nvec_phi;
        for ( i=0; i<offset; i++ ) {//offset = half of internal d.o.f. if depth==0
          VECTOR_LOOP( jj, nvec, jjj, buffer[i*nvec_phi+jj+jjj] = phi_pt[i*nvec_phi+jj+jjj];)
        }
        buffer += offset*nvec_phi;
      }
      buffer = (buffer_PRECISION)c->buffer[mu_dir];//move back to the head
      phi_pt = phi + comm_start*nvec_phi;// move to the end of inner vector

      if ( length[0] > 0 ) {
        PROF_PRECISION_START( _OP_COMM );
        MPI_Irecv( phi_pt, length[0]*nvec_phi, MPI_COMPLEX_PRECISION,
                   l->neighbor_rank[2*mu], 2*mu+1, g.comm_cart, &(c->rreqs[2*mu+1]) );
        PROF_PRECISION_STOP( _OP_COMM, 1 );
      }
      if ( length[1] > 0 ) {
        PROF_PRECISION_START( _OP_COMM );
        MPI_Isend( buffer, length[1]*nvec_phi, MPI_COMPLEX_PRECISION,
                   l->neighbor_rank[2*mu+1], 2*mu+1, g.comm_cart, &(c->sreqs[2*mu+1]) );
        PROF_PRECISION_STOP( _OP_COMM, 0 );
      }
      
    } else ASSERT( dir == 1 || dir == -1 );
  }
}

void ghost_wait_PRECISION( buffer_PRECISION phi, const int mu, const int dir,
                           comm_PRECISION_struct *c, const int amount, level_struct *l ) {
  
  /*
   * dir
   *   if +1: wait for incoming data arriving into buffer and set neg mu inner boundary sites from buffer
   *   if -1: wait
   */
  
  if( l->global_splitting[mu] > 1 ) {    
    int mu_dir = 2*mu-MIN(dir,0);
    int i, j, jj, jjj, *table, offset = c->offset, length[2]={0,0}, table_start = 0;
    int nvec = g.num_vect_pass1 ,nvec_phi=g.num_vect_pass2, nvec_com = c->num_vect;
    buffer_PRECISION buffer, phi_pt;

#ifdef DEBUG
    if ( nvec_com < nvec )
      error0("ghost_wait_PRECISION_: potential memory overflow (%d %d)\n", nvec_com,nvec);
#endif
    
    if ( amount == _FULL_SYSTEM ) {
      length[0] = (c->num_boundary_sites[2*mu])*offset;
      length[1] = (c->num_boundary_sites[2*mu+1])*offset;
      table_start = 0;
    } else if ( amount == _EVEN_SITES ) {
      length[0] = c->num_even_boundary_sites[2*mu]*offset;
      length[1] = c->num_even_boundary_sites[2*mu+1]*offset;
      table_start = 0;
    } else if ( amount == _ODD_SITES ) {
      length[0] = c->num_odd_boundary_sites[2*mu]*offset;
      length[1] = c->num_odd_boundary_sites[2*mu+1]*offset;
      table_start = c->num_even_boundary_sites[mu_dir];
    }

    ASSERT( c->in_use[mu_dir] == 1 );
    
    if ( dir == 1 ) {
      
      int num_boundary_sites = length[0]/offset;
      
      buffer = (buffer_PRECISION)c->buffer[mu_dir];      
      table = c->boundary_table[2*mu+1] + table_start;
      
      if ( length[0] > 0 ) {
        PROF_PRECISION_START( _OP_IDLE );
        MPI_Wait( &(c->sreqs[2*mu]), MPI_STATUS_IGNORE );
        PROF_PRECISION_STOP( _OP_IDLE, 0 );
      }
      if ( length[1] > 0 ) {
        PROF_PRECISION_START( _OP_IDLE );
        MPI_Wait( &(c->rreqs[2*mu]), MPI_STATUS_IGNORE );
        PROF_PRECISION_STOP( _OP_IDLE, 1 );
      }
      
      if ( l->depth == 0 ) {
        for ( j=0; j<num_boundary_sites; j++ ) {
          phi_pt = phi + table[j]*offset*nvec_phi;
	  for ( i=0; i<offset; i++ )
            VECTOR_LOOP( jj, nvec, jjj, phi_pt[i*nvec_phi+jj+jjj] = buffer[i*nvec_phi+jj+jjj];)
          buffer += offset*nvec_phi;
        }
      } else {
        for ( j=0; j<num_boundary_sites; j++ ) {
          phi_pt = phi + table[j]*offset*nvec_phi;
          for ( i=0; i<offset; i++ )
            VECTOR_LOOP( jj, nvec, jjj, phi_pt[i*nvec_phi+jj+jjj] += buffer[i*nvec_phi+jj+jjj];)
          buffer += offset*nvec_phi;
        }
      }
    } else if ( dir == -1 ) {
      
      if ( length[1] > 0 ) {
        PROF_PRECISION_START( _OP_IDLE );
        MPI_Wait( &(c->sreqs[2*mu+1]), MPI_STATUS_IGNORE );
        PROF_PRECISION_STOP( _OP_IDLE, 0 );
      }
      if ( length[0] > 0 ) {
        PROF_PRECISION_START( _OP_IDLE );
        MPI_Wait( &(c->rreqs[2*mu+1]), MPI_STATUS_IGNORE );
        PROF_PRECISION_STOP( _OP_IDLE, 1 );
      }
    } else ASSERT( dir == 1 || dir == -1 );
    
    c->in_use[mu_dir] = 0;
  }
}

void ghost_update_PRECISION( vector_PRECISION *phi, const int mu, const int dir, comm_PRECISION_struct *c, level_struct *l ) {
  /*
   * updates the full ghost shell
   *
   * Assume: phi is of size _SCHWARZ
   * dir == -1: updates negative ghost shell
   * dir == +1: updates positive ghost shell
   * Comment: used only in Schwarz method
   */
  if( l->global_splitting[mu] > 1 ) {
    int i, j, jj, jjj, mu_dir = 2*mu-MIN(dir,0), nu, inv_mu_dir = 2*mu+1+MIN(dir,0), length, *table=NULL,
      comm_start, num_boundary_sites, site_var;
    int nvec_phi = phi->num_vect, nvec_com=c->num_vect;
    buffer_PRECISION buffer, recv_pt, phi_pt;

#ifdef DEBUG
    if ( nvec_com < nvec_phi )
      error0("ghost_update_PRECISION: potential memory overflow\n");
#endif
    
    site_var = l->num_lattice_site_var;
    length = c->num_boundary_sites[mu_dir]*l->num_lattice_site_var*nvec_phi;
    num_boundary_sites = c->num_boundary_sites[mu_dir];
    buffer = c->buffer[mu_dir];
    
    if ( dir == -1 ) // for updating negative ghost shell
      comm_start = l->vector_size;
    else             // for updating positive ghost shell
      comm_start = l->inner_vector_size;
    for ( nu=0; nu<mu; nu++ ) 
      comm_start += c->num_boundary_sites[2*nu]*l->num_lattice_site_var;
    
    ASSERT( c->in_use[mu_dir] == 0 );
    c->in_use[mu_dir] = 1;
    
    recv_pt = phi->vector_buffer + comm_start*nvec_phi;
    if ( length > 0 ) {
      PROF_PRECISION_START( _OP_COMM );
      MPI_Irecv( recv_pt, length, MPI_COMPLEX_PRECISION,
                 l->neighbor_rank[mu_dir], mu_dir, g.comm_cart, &(c->rreqs[mu_dir]) );
      PROF_PRECISION_STOP( _OP_COMM, 1 );
    }
    table = c->boundary_table[inv_mu_dir];
    for ( j=0; j<num_boundary_sites; j++ ) {
      phi_pt = phi->vector_buffer + table[j]*site_var*nvec_phi;
      
      for ( i=0; i<site_var; i++ ) {
        VECTOR_LOOP(jj, nvec_phi, jjj, buffer[i*nvec_phi+jj+jjj] = phi_pt[i*nvec_phi+jj+jjj];)
      }
      buffer += site_var*nvec_phi;
    }
    buffer = c->buffer[mu_dir];// move back the poiner to the head
    
    if ( length > 0 ) {
      PROF_PRECISION_START( _OP_COMM );
      MPI_Isend( buffer, length, MPI_COMPLEX_PRECISION,
                 l->neighbor_rank[inv_mu_dir], mu_dir, g.comm_cart, &(c->sreqs[mu_dir]) );
      PROF_PRECISION_STOP( _OP_COMM, 0 );
    }
  }
}

// used only in Schwarz method
void ghost_update_wait_PRECISION( vector_PRECISION *phi, const int mu, const int dir, comm_PRECISION_struct *c, level_struct *l ) {
  
  if( l->global_splitting[mu] > 1 ) {
    int mu_dir = 2*mu-MIN(dir,0), length = c->num_boundary_sites[mu_dir]*l->num_lattice_site_var;
    
    ASSERT( c->in_use[mu_dir] == 1 );
      
    if ( length > 0 ) {
      PROF_PRECISION_START( _OP_IDLE );
      MPI_Wait( &(c->sreqs[mu_dir]), MPI_STATUS_IGNORE );
      MPI_Wait( &(c->rreqs[mu_dir]), MPI_STATUS_IGNORE );
      PROF_PRECISION_STOP( _OP_IDLE, 1 );
    }
    c->in_use[mu_dir] = 0;
  }
}

void negative_sendrecv_PRECISION( vector_PRECISION *phi, const int mu, comm_PRECISION_struct *c, level_struct *l ) {
  /**************************
   * dir = -1
   *   send neg mu inner boundary sites of this rank to the rank in the neg mu dir
   *   recieve data from the rank in the pos mu dir and update pos mu ghost shell
   **************************/
  if( l->global_splitting[mu] > 1 ) {    
    
    int i, j, jj, boundary_start;
    int num_boundary_sites = c->num_boundary_sites[2*mu+1];
    int *boundary_table    = c->boundary_table[2*mu+1];
    int n                  = l->num_lattice_site_var;
    int nvec = phi->num_vect;

#ifdef DEBUG
    if ( l->vbuf_PRECISION[4].num_vect < nvec )
      error0("negative_sendrecv_PRECISION: potential memory overflow\n");
#endif
    
    buffer_PRECISION buffer, tmp_pt, buffer_pt;

    boundary_start = l->num_inner_lattice_sites;
    for ( i=0; i<mu; i++ )
      boundary_start += c->num_boundary_sites[2*i]; // starting index of ghost cells of the negative mu boundary

    buffer    = l->vbuf_PRECISION[4].vector_buffer+n*(boundary_start-l->num_inner_lattice_sites)*nvec;
    buffer_pt = buffer;
    // for each site on the negative mu inner boundary, set buffer_pt, i.e., buffer, to the value of phi at the site 
    for ( i=0; i<num_boundary_sites; i++ ) {
      tmp_pt = phi->vector_buffer + n*boundary_table[i]*nvec;
      VECTOR_LOOP( j, n*nvec, jj, *buffer_pt = *tmp_pt; buffer_pt++; tmp_pt++;)
    }
    
    // recieve ghost cells of positive mu boundary from the process in the pos mu dir
    MPI_Irecv( phi->vector_buffer+n*boundary_start*nvec, n*num_boundary_sites*nvec, MPI_COMPLEX_PRECISION,
               l->neighbor_rank[2*mu], 2*mu+1, g.comm_cart, &(c->rreqs[2*mu+1]) );
    // send inner boundary sites on the negative mu boundary to the rank in the neg mu dir
    MPI_Isend( buffer, n*num_boundary_sites*nvec, MPI_COMPLEX_PRECISION,
               l->neighbor_rank[2*mu+1], 2*mu+1, g.comm_cart, &(c->sreqs[2*mu+1]) );
  }
}
 
void negative_wait_PRECISION( const int mu, comm_PRECISION_struct *c, level_struct *l ) {
 
  if( l->global_splitting[mu] > 1 ) {
    MPI_Wait( &(c->sreqs[2*mu+1]), MPI_STATUS_IGNORE );
    MPI_Wait( &(c->rreqs[2*mu+1]), MPI_STATUS_IGNORE );
  }
}

