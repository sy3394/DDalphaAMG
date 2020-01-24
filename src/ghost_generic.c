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
  
  int mu, nu, factor=1, n_vect = (g.num_rhs_vect < l->num_eig_vect)? l->num_eig_vect:g.num_rhs_vect;//!!!!!!!
  
  c->num_vect = n_vect;
  if ( l->depth > 0 ) {
    c->offset = l->num_lattice_site_var;
  } else {
    c->offset = l->num_lattice_site_var/2;//why divide by 2???? probably b/c only half of d.o.f are used; other half projected out
    if ( g.method < 5 )
      factor = 2;
  }

#ifdef HAVE_TM1p1
  factor *= 2;
#endif

  //--- set c->comm_start,c->length, and c->max_length and allocate c->buffer
  if ( buffer_size <= 0 ) {
    c->comm_start[0] = c->offset*l->num_inner_lattice_sites;//*n_vect;
    c->comm_start[1] = c->offset*l->num_inner_lattice_sites;//*n_vect;
    for ( mu=0; mu<4; mu++ ) {
      if ( mu > 0 ) {
        c->comm_start[2*mu] = c->comm_start[2*(mu-1)] + buffer_size;//*n_vect;
        c->comm_start[2*mu+1] = c->comm_start[2*(mu-1)+1] + buffer_size;//*n_vect;!!!!!!!!!
      }
      buffer_size = c->offset;
      for ( nu=0; nu<4; nu++ ) {
        if ( nu != mu ) {
          buffer_size *= l->local_lattice[nu];
        }
      }//printf("ghost aloc mu:%d %d %d %d\n", mu,n_vect, factor,buffer_size);
      c->length[2*mu]  = buffer_size;
      c->length[2*mu+1] = buffer_size;
      c->max_length[mu] = factor*buffer_size;
      MALLOC( c->buffer[2*mu],   complex_PRECISION, factor*buffer_size*n_vect );
      MALLOC( c->buffer[2*mu+1], complex_PRECISION, factor*buffer_size*n_vect );
      c->in_use[2*mu] = 0;
      c->in_use[2*mu+1] = 0;
    }
  } else {
    for ( mu=0; mu<4; mu++ ) {//printf("ghost aloc mu: else\n");
      c->max_length[mu] = buffer_size;
#ifdef HAVE_TM1p1
      MALLOC( c->buffer[2*mu],   complex_PRECISION, 2*buffer_size*n_vect );
      MALLOC( c->buffer[2*mu+1], complex_PRECISION, 2*buffer_size*n_vect );
#else
      MALLOC( c->buffer[2*mu],   complex_PRECISION, buffer_size*n_vect );
      MALLOC( c->buffer[2*mu+1], complex_PRECISION, buffer_size*n_vect );
#endif
    }
  }//printf("ghost aloc %d %d\n", n_vect, buffer_size);

  //--- allocate l->vbuf_PRECISION[8]
 if ( l->vbuf_PRECISION[8].vector_buffer == NULL ) {
   // no_threading is not allocated yet
#ifdef HAVE_TM1p1
   //vector_PRECISION_alloc( &(l->vbuf_PRECISION[8]), _ORDINARY, 2*n_vect, l, no_threading);
    MALLOC( l->vbuf_PRECISION[8].vector_buffer, complex_PRECISION, 2*l->vector_size*n_vect );
    l->vbuf_PRECISION[8].type = _ORDINARY;
    l->vbuf_PRECISION[8].size = l->vector_size;
    l->vbuf_PRECISION[8].num_vect = n_vect;
    l->vbuf_PRECISION[8].l = l;
#else
    //vector_PRECISION_alloc( &(l->vbuf_PRECISION[8]), _ORDINARY, n_vect, l, no_threading);
    MALLOC( l->vbuf_PRECISION[8].vector_buffer, complex_PRECISION, l->vector_size*n_vect );
    l->vbuf_PRECISION[8].type = _ORDINARY;
    l->vbuf_PRECISION[8].size = l->vector_size;
    l->vbuf_PRECISION[8].num_vect = n_vect;
    l->vbuf_PRECISION[8].l = l;
#endif
  }
}

void ghost_free_PRECISION( comm_PRECISION_struct *c, level_struct *l ) {
  
  int mu, n_vect = (g.num_rhs_vect < l->num_eig_vect)? l->num_eig_vect:g.num_rhs_vect;//!!!!!!!!!!
  
  for ( mu=0; mu<4; mu++ ) {    
    FREE( c->buffer[2*mu],   complex_PRECISION, c->max_length[mu]*n_vect );
    FREE( c->buffer[2*mu+1], complex_PRECISION, c->max_length[mu]*n_vect );
  }
  if ( l->vbuf_PRECISION[8].vector_buffer != NULL ){
    vector_PRECISION_free( &(l->vbuf_PRECISION[8]), l, no_threading);
    /*
#ifdef HAVE_TM1p1		
     FREE( l->vbuf_PRECISION[8].vector_buffer, complex_PRECISION, 2*l->vector_size );		
 #else		
     FREE( l->vbuf_PRECISION[8].vector_buffer, complex_PRECISION, l->vector_size );		
 #endif
    */
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

/*
void ghost_sendrecv_PRECISION( buffer_PRECISION phi, const int mu, const int dir,
                               comm_PRECISION_struct *c, const int amount, level_struct *l ) {
  // does not allow sending in both directions at the same time
  if( l->global_splitting[mu] > 1 ) {
    
    int i, j, *table=NULL, mu_dir = 2*mu-MIN(dir,0), offset = c->offset,
        length[2] = {0,0}, comm_start = 0, table_start = 0;
    buffer_PRECISION buffer, phi_pt;
    
    if ( amount == _FULL_SYSTEM ) {
      length[0] = (c->num_boundary_sites[2*mu])*offset;
      length[1] = (c->num_boundary_sites[2*mu+1])*offset;
      comm_start = c->comm_start[mu_dir];
      table_start = 0;
    } else if ( amount == _EVEN_SITES ) {
      length[0] = c->num_even_boundary_sites[2*mu]*offset;
      length[1] = c->num_even_boundary_sites[2*mu+1]*offset;
      comm_start = c->comm_start[mu_dir];
      table_start = 0;
    } else if ( amount == _ODD_SITES ) {
      length[0] = c->num_odd_boundary_sites[2*mu]*offset;
      length[1] = c->num_odd_boundary_sites[2*mu+1]*offset;
      comm_start = c->comm_start[mu_dir]+c->num_even_boundary_sites[mu_dir]*offset;
      table_start = c->num_even_boundary_sites[mu_dir];
    }
    
#ifdef HAVE_TM1p1
    if ( g.n_flavours == 2 ) {
      length[0] *= 2;
      length[1] *= 2;
      comm_start *= 2;
      offset *= 2;
    }
#endif

    ASSERT( c->in_use[mu_dir] == 0 );
    c->in_use[mu_dir] = 1;
    
    if ( MAX(length[0],length[1]) > c->max_length[mu] ) {
      printf("CAUTION: my_rank: %d, not enough comm buffer\n", g.my_rank ); fflush(0);
      ghost_free_PRECISION( c, l );
      ghost_alloc_PRECISION( MAX(length[0],length[1]), c, l );
    }
    
    buffer = c->buffer[mu_dir];
    
    // dir = senddir
    if ( dir == 1 ) {
      // data to be communicated is stored serially in the vector phi
      // recv target is a buffer
      // afterwards (in ghost_wait) the data has to be distributed onto the correct sites
      // touching the respective boundary in -mu direction
      
      phi_pt = phi + comm_start;
      if ( length[1] > 0 ) {
        PROF_PRECISION_START( _OP_COMM );
        MPI_Irecv( buffer, length[1], MPI_COMPLEX_PRECISION,
                   l->neighbor_rank[2*mu+1], 2*mu, g.comm_cart, &(c->rreqs[2*mu]) );
        PROF_PRECISION_STOP( _OP_COMM, 1 );
      }
      if ( length[0] > 0 ) {
        PROF_PRECISION_START( _OP_COMM );
        MPI_Isend( phi_pt, length[0], MPI_COMPLEX_PRECISION,
                   l->neighbor_rank[2*mu], 2*mu, g.comm_cart, &(c->sreqs[2*mu]) );
        PROF_PRECISION_STOP( _OP_COMM, 0 );
      }
      
      
    } else if ( dir == -1 ) {
      // data to be communicated is stored on the sites touching the boundary in -mu direction
      // this data is gathered in a buffer in the correct ordering
      // which is required on the boundary of the vector phi
      int num_boundary_sites = length[1]/offset;
      
      table = c->boundary_table[2*mu+1]+table_start;
      for ( j=0; j<num_boundary_sites; j++ ) {
        phi_pt = phi + table[j]*offset;
        for ( i=0; i<offset; i++ ) {
          buffer[i] = phi_pt[i];
        }
        buffer += offset;
      }
      
      buffer = c->buffer[mu_dir];      
      phi_pt = phi + comm_start;
      
      if ( length[0] > 0 ) {
        PROF_PRECISION_START( _OP_COMM );
        MPI_Irecv( phi_pt, length[0], MPI_COMPLEX_PRECISION,
                   l->neighbor_rank[2*mu], 2*mu+1, g.comm_cart, &(c->rreqs[2*mu+1]) );
        PROF_PRECISION_STOP( _OP_COMM, 1 );
      }
      if ( length[1] > 0 ) {
        PROF_PRECISION_START( _OP_COMM );
        MPI_Isend( buffer, length[1], MPI_COMPLEX_PRECISION,
                   l->neighbor_rank[2*mu+1], 2*mu+1, g.comm_cart, &(c->sreqs[2*mu+1]) );
        PROF_PRECISION_STOP( _OP_COMM, 0 );
      }
      
    } else ASSERT( dir == 1 || dir == -1 );
  }
}
*/

//if dir == +1
// send sites in the dir mu direction and 
// get sites from the negative dir mu direction and store them in c->buffer[mu_dir] (Then, c->in_use[mu_dir] flag is set on)
// After the call of this function, call ghost_wait to make sure that the data is delivered and store the data in phi
//if dir == -1
// copy the values to be send onto buffer and then send it
// recieve and store the values of ghost cells of the positive mu dir sonto phi
// After the call of this function, call ghost_update to make sure that the data are transferred onto phi
void ghost_sendrecv_PRECISION_new( buffer_PRECISION phi, const int mu, const int dir,
				   comm_PRECISION_struct *c, const int amount, level_struct *l ) {
  /**************************
   * buffer_PRECISION phi: vectors to be sent
   * mu: specify boundary
   * dir: eiter postive (update ghose shell) or negative (update inner boundary sites)
   *   if +1: send sites of phi in pos mu boundary ghost shell to the rank in the pos mu dir and store data from the rank in the neg mu dir in buffer
   *   if -1: send neg mu inner boudnary sites to the rank in the neg mu dir and update pos mu ghost shell
   *************************/
  // does not allow sending in both directions at the same time
  if( l->global_splitting[mu] > 1 ) {

    int i, j, jj, jjj, *table=NULL, mu_dir = 2*mu-MIN(dir,0), offset = c->offset, length[2] = {0,0}, comm_start = 0, table_start = 0;
    int nvec = g.num_vect_pass1, nvec_phi = g.num_vect_pass2, nvec_com=c->num_vect;//temp fix!!!!!!!!
    buffer_PRECISION buffer, phi_pt;
    
    // offset is internal d.o.f. if at the top and half of that if not
    if ( amount == _FULL_SYSTEM ) {
      length[0] = (c->num_boundary_sites[2*mu])*offset;   // #inner boundary sites on the positive mu dir x offset
      length[1] = (c->num_boundary_sites[2*mu+1])*offset; // #inner boundary sites on the negative mu dir x offset
      comm_start = c->comm_start[mu_dir];
      table_start = 0;
    } else if ( amount == _EVEN_SITES ) {
      length[0] = c->num_even_boundary_sites[2*mu]*offset;
      length[1] = c->num_even_boundary_sites[2*mu+1]*offset;
      comm_start = c->comm_start[mu_dir];
      table_start = 0;
    } else if ( amount == _ODD_SITES ) {
      length[0] = c->num_odd_boundary_sites[2*mu]*offset;
      length[1] = c->num_odd_boundary_sites[2*mu+1]*offset;
      comm_start = c->comm_start[mu_dir]+c->num_even_boundary_sites[mu_dir]*offset;
      table_start = c->num_even_boundary_sites[mu_dir];
    }
    //printf0("ghost_sendrecv_PRECISION: %d %d %d %d\n",length[0],length[1],comm_start,table_start);
/*#ifdef HAVE_TM1p1
    if ( g.n_flavours == 2 ) {
      length[0] *= 2;
      length[1] *= 2;
      comm_start *= 2;
      offset *= 2;
    }
#endif*/

    ASSERT( c->in_use[mu_dir] == 0 );
    c->in_use[mu_dir] = 1;
    //    printf(" ghost_sendrecv_PRECISION_new %d %d %d %d %d %d\n",mu,length[0],length[1],c->max_length[mu], nvec_phi,c->num_vect);    
    if ( MAX(length[0],length[1]) > c->max_length[mu] ) {// I took out nvec_com
      printf("CAUTION: my_rank: %d, not enough comm buffer\n", g.my_rank ); fflush(0);
      ghost_free_PRECISION( c, l );
      ghost_alloc_PRECISION( MAX(length[0],length[1]), c, l );
    }
    
    buffer = (buffer_PRECISION)c->buffer[mu_dir];
    //    printf0("sendrev star %d %d %d %d %d %d %d\n", mu_dir,nvec,nvec_phi,length[0],length[1],c->max_length[mu],offset);
    // dir = senddir
    if ( dir == 1 ) {
      // data to be communicated is stored serially in the vector phi
      // recv target is a buffer
      // afterwards (in ghost_wait) the data has to be distributed onto the correct sites
      // touching the respective boundary in -mu direction
      
      phi_pt = phi + comm_start*nvec_phi;//!!!!!!!
      if ( length[1] > 0 ) {//printf0("gho +1 rev %d\n",length[1]*nvec_phi);
	// get ghost cells on the negative mu boundary from the rank in the negative mu dir
        PROF_PRECISION_START( _OP_COMM );
        MPI_Irecv( buffer, length[1]*nvec_phi, MPI_COMPLEX_PRECISION,//why use buffer in this case only???
                   l->neighbor_rank[2*mu+1], 2*mu, g.comm_cart, &(c->rreqs[2*mu]) );
        PROF_PRECISION_STOP( _OP_COMM, 1 );
      }
      if ( length[0] > 0 ) {//printf0("gho +1 send %d %d %d %d %d\n", l->neighbor_rank[2*mu],length[0],length[0]*nvec_phi,c->max_length[mu]*c->num_vect,c->num_vect);
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
        phi_pt = phi + table[j]*offset*nvec_phi;//printf0("send ghost table[%d] = %d\n",j,table[j]);
        for ( i=0; i<offset; i++ ) {//offset = half of internal d.o.f. if depth==0
          VECTOR_LOOP( jj, nvec, jjj, buffer[i*nvec_phi+jj+jjj] = phi_pt[i*nvec_phi+jj+jjj];)
        }
        buffer += offset*nvec_phi;
      }
      buffer = (buffer_PRECISION)c->buffer[mu_dir];//move back to the head
      phi_pt = phi + comm_start*nvec_phi;// move to the end of inner vector // nvec!!!!!!!!!
      
      if ( length[0] > 0 ) {//printf("gho -1 rev %d:%d %d %d %d\n",g.my_rank,2*mu+1,length[0],length[0]*nvec_phi,c->max_length[mu]*c->num_vect);
        PROF_PRECISION_START( _OP_COMM );
        MPI_Irecv( phi_pt, length[0]*nvec_phi, MPI_COMPLEX_PRECISION,
                   l->neighbor_rank[2*mu], 2*mu+1, g.comm_cart, &(c->rreqs[2*mu+1]) );
        PROF_PRECISION_STOP( _OP_COMM, 1 );
      }
      if ( length[1] > 0 ) {//printf("gho -1 send %d: %d %d %d %d %d\n",g.my_rank,2*mu+1,comm_start,offset,length[1]*nvec_phi,c->max_length[mu]*c->num_vect);
        PROF_PRECISION_START( _OP_COMM );
	//MPI_Request sreq;
	//for(int kk=0;kk<length[1]*nvec_phi;kk++)printf(":%d %g: ",g.my_rank,creal_PRECISION(buffer[kk]));
        MPI_Isend( buffer, length[1]*nvec_phi, MPI_COMPLEX_PRECISION,
                   l->neighbor_rank[2*mu+1], 2*mu+1, g.comm_cart, &(c->sreqs[2*mu+1]) );
	//MPI_Wait( &sreq, MPI_STATUS_IGNORE );
        PROF_PRECISION_STOP( _OP_COMM, 0 );
      }
      
    } else ASSERT( dir == 1 || dir == -1 );
  }
  //  printf0("sendrev end\n");
}

/*
void ghost_wait_PRECISION( buffer_PRECISION phi, const int mu, const int dir,
                           comm_PRECISION_struct *c, const int amount, level_struct *l ) {
  
  if( l->global_splitting[mu] > 1 ) {    
    int mu_dir = 2*mu-MIN(dir,0);
    int i, j, *table, offset = c->offset, length[2]={0,0}, table_start = 0;
    buffer_PRECISION buffer, phi_pt;

#ifdef HAVE_TM1p1
    if ( g.n_flavours == 2 )
      offset *= 2;
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
      
      buffer = c->buffer[mu_dir];      
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
          phi_pt = phi + table[j]*offset;
          
          for ( i=0; i<offset; i++ )
            phi_pt[i] = buffer[i];
          
          buffer += offset;
        }
      } else {
        for ( j=0; j<num_boundary_sites; j++ ) {
          phi_pt = phi + table[j]*offset;
          
          for ( i=0; i<offset; i++ )
            phi_pt[i] += buffer[i];
          
          buffer += offset;
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
*/

// if dir is positive, wait (make sure to recieve the data) and update the ghost cells
// if dir is negative, just wait
void ghost_wait_PRECISION_new( buffer_PRECISION phi, const int mu, const int dir,
                           comm_PRECISION_struct *c, const int amount, level_struct *l ) {
  
  /*
   * dir
   *   if +1: wait and set neg mu inner boundary sites from buffer
   *   if -1: wait
   */
  if( l->global_splitting[mu] > 1 ) {    
    int mu_dir = 2*mu-MIN(dir,0);
    int i, j, jj, jjj, *table, offset = c->offset, length[2]={0,0}, table_start = 0;
    int nvec = g.num_vect_pass1 ,nvec_phi=g.num_vect_pass2, nvec_com = c->num_vect;
    buffer_PRECISION buffer, phi_pt;

/*#ifdef HAVE_TM1p1
    if ( g.n_flavours == 2 )
      offset *= 2;
#endif*/
      
    if ( nvec_com < nvec )
      error0("ghost_wait_PRECISION_: potential memory overflow\n");

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
    //    printf0("gho wait0\n");
    ASSERT( c->in_use[mu_dir] == 1 );
    
    if ( dir == 1 ) {//printf0("gho wait+1\n");
      
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
    } else if ( dir == -1 ) {//printf0("gho wait-1\n");
      
      if ( length[1] > 0 ) {//printf0("gho wait-1201 %d\n", 2*mu+1);
        PROF_PRECISION_START( _OP_IDLE );
        MPI_Wait( &(c->sreqs[2*mu+1]), MPI_STATUS_IGNORE );//printf0("gho wait-121\n");
        PROF_PRECISION_STOP( _OP_IDLE, 0 );//printf0("gho wait-11\n");
      }
      if ( length[0] > 0 ) {
        PROF_PRECISION_START( _OP_IDLE );//printf0("gho wait-120\n");
        MPI_Wait( &(c->rreqs[2*mu+1]), MPI_STATUS_IGNORE );
        PROF_PRECISION_STOP( _OP_IDLE, 1 );//printf0("gho wait-12\n");
      }
      
    } else ASSERT( dir == 1 || dir == -1 );
    
    c->in_use[mu_dir] = 0;
  }
}

/*
void ghost_update_PRECISION( vector_PRECISION *phi, const int mu, const int dir, comm_PRECISION_struct *c, level_struct *l ) {
  
  if( l->global_splitting[mu] > 1 ) {
    int i, j, mu_dir = 2*mu-MIN(dir,0), nu, inv_mu_dir = 2*mu+1+MIN(dir,0), length, *table=NULL,
        comm_start, num_boundary_sites, site_var;
    buffer_PRECISION buffer, recv_pt, phi_pt;
    
    site_var = l->num_lattice_site_var;
    length = c->num_boundary_sites[mu_dir]*l->num_lattice_site_var;
    num_boundary_sites = c->num_boundary_sites[mu_dir];
    buffer = c->buffer[mu_dir];
    
    if ( dir == -1 )
      comm_start = l->vector_size;
    else
      comm_start = l->inner_vector_size;
    for ( nu=0; nu<mu; nu++ ) {
      comm_start += c->num_boundary_sites[2*nu]*l->num_lattice_site_var;
    }
    
    ASSERT( c->in_use[mu_dir] == 0 );
    c->in_use[mu_dir] = 1;
    
    recv_pt = phi->vector_buffer + comm_start;
    if ( length > 0 ) {
      PROF_PRECISION_START( _OP_COMM );
      MPI_Irecv( recv_pt, length, MPI_COMPLEX_PRECISION,
                 l->neighbor_rank[mu_dir], mu_dir, g.comm_cart, &(c->rreqs[mu_dir]) );
      PROF_PRECISION_STOP( _OP_COMM, 1 );
    }
    
    table = c->boundary_table[inv_mu_dir];
    for ( j=0; j<num_boundary_sites; j++ ) {
      phi_pt = phi->vector_buffer + table[j]*site_var;
      
      for ( i=0; i<site_var; i++ ) {
        buffer[i] = phi_pt[i];
      }
      buffer += site_var;
    }
    buffer = c->buffer[mu_dir];
    
    if ( length > 0 ) {
      PROF_PRECISION_START( _OP_COMM );
      MPI_Isend( buffer, length, MPI_COMPLEX_PRECISION,
                 l->neighbor_rank[inv_mu_dir], mu_dir, g.comm_cart, &(c->sreqs[mu_dir]) );
      PROF_PRECISION_STOP( _OP_COMM, 0 );
    }
  }
}
*/

//!!!!!!!!!
// used only in Schwarz method: update the full ghost shell
// assume: phi is of size _SCHWARZ
void ghost_update_PRECISION_new( vector_PRECISION *phi, const int mu, const int dir, comm_PRECISION_struct *c, level_struct *l ) {
  //error0("ghost_update_PRECISION_new: not corrected\n");
  if( l->global_splitting[mu] > 1 ) {
    //printf("ghost_update_PRECISION_new");
    int i, j, jj, jjj, mu_dir = 2*mu-MIN(dir,0), nu, inv_mu_dir = 2*mu+1+MIN(dir,0), length, *table=NULL,
      comm_start, num_boundary_sites, site_var;
    int nvec_phi = phi->num_vect, nvec_com=c->num_vect;
    buffer_PRECISION buffer, recv_pt, phi_pt;
    
    if ( nvec_com < nvec_phi )
      error0("ghost_update_PRECISION: potential memory overflow\n");

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
    //printf("ghost_update_PRECISION_new2");
    table = c->boundary_table[inv_mu_dir];
    for ( j=0; j<num_boundary_sites; j++ ) {
      phi_pt = phi->vector_buffer + table[j]*site_var*nvec_phi;
      
      for ( i=0; i<site_var; i++ ) {
        VECTOR_LOOP(jj, nvec_phi, jjj, buffer[i*nvec_phi+jj+jjj] = phi_pt[i*nvec_phi+jj+jjj];)
      }
      buffer += site_var*nvec_phi;
    }//printf("ghost_update_PRECISION_new3");
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

/*
void negative_sendrecv_PRECISION( vector_PRECISION *phi, const int mu, comm_PRECISION_struct *c, level_struct *l ) {
  // send dir = -1
  if( l->global_splitting[mu] > 1 ) {    
    
    int i, j, num_boundary_sites = c->num_boundary_sites[2*mu+1], boundary_start,
      *boundary_table = c->boundary_table[2*mu+1], n = l->num_lattice_site_var;

    vector_PRECISION buffer, tmp_pt, buffer_pt;
    
    boundary_start = l->num_inner_lattice_sites;
    for ( i=0; i<mu; i++ )
      boundary_start += c->num_boundary_sites[2*i];

    buffer.vector_buffer  = l->vbuf_PRECISION[8].vector_buffer+n*(boundary_start-l->num_inner_lattice_sites);
    buffer_pt.vector_buffer = buffer.vector_buffer;
    
    for ( i=0; i<num_boundary_sites; i++ ) {
      tmp_pt.vector_buffer = phi->vector_buffer + n*boundary_table[i];
      for ( j=0; j<n; j++, buffer_pt.vector_buffer++, tmp_pt.vector_buffer++ )
        *buffer_pt.vector_buffer = *tmp_pt.vector_buffer;
    }
    
    MPI_Irecv( phi->vector_buffer+n*boundary_start, n*num_boundary_sites, MPI_COMPLEX_PRECISION,
               l->neighbor_rank[2*mu], 2*mu+1, g.comm_cart, &(c->rreqs[2*mu+1]) );
    MPI_Isend( buffer.vector_buffer, n*num_boundary_sites, MPI_COMPLEX_PRECISION,
               l->neighbor_rank[2*mu+1], 2*mu+1, g.comm_cart, &(c->sreqs[2*mu+1]) );
  }
}
*/

// send inner boundary sites on the negative mu boundary in the negative mu dir and
// recieve inner boundary sites from rank in the positive mu dir and store the associated ghost cells
void negative_sendrecv_PRECISION_new( vector_PRECISION *phi, const int mu, comm_PRECISION_struct *c, level_struct *l ) {
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

    if ( l->vbuf_PRECISION[8].num_vect < nvec )
      error0("negative_sendrecv_PRECISION: potential memory overflow\n");
    
    buffer_PRECISION buffer, tmp_pt, buffer_pt;

    boundary_start = l->num_inner_lattice_sites;
    for ( i=0; i<mu; i++ )
      boundary_start += c->num_boundary_sites[2*i]; // starting index of ghost cells of the negative mu boundary

    //printf("negative_sendrecv_PRECISION_new %d %d %d %d %d\n",nvec, l->vbuf_PRECISION[8].num_vect, l->vbuf_PRECISION[8].size,n*(boundary_start-l->num_inner_lattice_sites),n*num_boundary_sites);fflush(stdout);
    buffer    = l->vbuf_PRECISION[8].vector_buffer+n*(boundary_start-l->num_inner_lattice_sites)*nvec;
    buffer_pt = buffer;
    
    // for each site on the negative mu inner boundary, set buffer_pt, i.e., buffer, to the value of phi at the site 
    for ( i=0; i<num_boundary_sites; i++ ) {
      tmp_pt = phi->vector_buffer + n*boundary_table[i]*nvec;
      //VECTOR_LOOP( j, n*nvec, jj, *buffer_pt = *tmp_pt; buffer_pt++; tmp_pt++;)
      
      for ( j=0; j<n*nvec; j++, buffer_pt++, tmp_pt++ )//why not VECTOR_LOOP????
	*buffer_pt = *tmp_pt;
    }
    
    // recieve ghost cells of positive mu boundary from the process in the pos mu dir
    MPI_Irecv( phi->vector_buffer+n*boundary_start*nvec, n*num_boundary_sites*nvec, MPI_COMPLEX_PRECISION,
               l->neighbor_rank[2*mu], 2*mu+1, g.comm_cart, &(c->rreqs[2*mu+1]) );
    // send inner boundary sites on the negative mu boundary to the rank in the neg mu dir
    MPI_Isend( buffer, n*num_boundary_sites*nvec, MPI_COMPLEX_PRECISION,
               l->neighbor_rank[2*mu+1], 2*mu+1, g.comm_cart, &(c->sreqs[2*mu+1]) );
  }
}/*
void negative_sendrecv_PRECISION_new( vector_PRECISION *phi, const int mu, comm_PRECISION_struct *c, level_struct *l ) {
  // send dir = -1  
  if( l->global_splitting[mu] > 1 ) {

    int i, j, num_boundary_sites = c->num_boundary_sites[2*mu+1], boundary_start,
      *boundary_table = c->boundary_table[2*mu+1], n = l->num_lattice_site_var;
    int n_vect = phi->num_vect;
    n = n *n_vect;

    vector_PRECISION buffer, tmp_pt, buffer_pt;

    boundary_start = l->num_inner_lattice_sites;
    for ( i=0; i<mu; i++ )
      boundary_start += c->num_boundary_sites[2*i];

    buffer.vector_buffer  = l->vbuf_PRECISION[8].vector_buffer+n*(boundary_start-l->num_inner_lattice_sites);
    buffer_pt.vector_buffer = buffer.vector_buffer;

    for ( i=0; i<num_boundary_sites; i++ ) {
      tmp_pt.vector_buffer = phi->vector_buffer + n*boundary_table[i];
      for ( j=0; j<n; j++, buffer_pt.vector_buffer++, tmp_pt.vector_buffer++ )
        *buffer_pt.vector_buffer = *tmp_pt.vector_buffer;
    }

    MPI_Irecv( phi->vector_buffer+n*boundary_start, n*num_boundary_sites, MPI_COMPLEX_PRECISION,
               l->neighbor_rank[2*mu], 2*mu+1, g.comm_cart, &(c->rreqs[2*mu+1]) );
    MPI_Isend( buffer.vector_buffer, n*num_boundary_sites, MPI_COMPLEX_PRECISION,
               l->neighbor_rank[2*mu+1], 2*mu+1, g.comm_cart, &(c->sreqs[2*mu+1]) );
  }
  }*/
 
void negative_wait_PRECISION( const int mu, comm_PRECISION_struct *c, level_struct *l ) {
 
  if( l->global_splitting[mu] > 1 ) {//printf("negative_wait_PRECISION %d",mu);
    MPI_Wait( &(c->sreqs[2*mu+1]), MPI_STATUS_IGNORE );//printf("negative_wait_PRECISION send %d",mu);
    MPI_Wait( &(c->rreqs[2*mu+1]), MPI_STATUS_IGNORE );//printf("negative_wait_PRECISION recv  %d",mu);
  }
}

