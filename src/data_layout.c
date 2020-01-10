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
 * not changed from sbacchio
 * glanced over: 12/08/2019
 * glanced over:12/19/2019
 * confirmed: not changed from milla:01/02/2020
 */

#include "main.h"

// Compute volume of local lattice and size of local vectors
void data_layout_init( level_struct *l ) {
  
  int i, j;
  
  l->num_inner_lattice_sites = 1;
  for ( i=0; i<4; i++ )
    l->num_inner_lattice_sites *= l->local_lattice[i];  
  l->num_lattice_sites = l->num_inner_lattice_sites;
  l->inner_vector_size = l->num_inner_lattice_sites * l->num_lattice_site_var;
  
  j = l->num_lattice_sites;
  for ( i=0; i<4; i++ ) 
    l->num_lattice_sites += j / l->local_lattice[i];
  
  l->vector_size = l->num_lattice_sites * l->num_lattice_site_var;
  l->schwarz_vector_size = 2*l->vector_size - l->inner_vector_size;
}

void data_layout_n_flavours( int nf, level_struct *l, struct Thread *threading ) {

  ASSERT(nf>0);
  ASSERT(l->depth == 0);

#ifdef HAVE_TM1p1
  ASSERT(nf<=2);
  
  if( g.n_flavours == nf )
    return;
  else {
    START_LOCKED_MASTER(threading)
    g.n_flavours = nf;

    struct level_struct *l_tmp = l;
    
    while(1) {
      if(l_tmp->depth == 0)
        l_tmp->num_lattice_site_var = nf * 12;
      else
        l_tmp->num_lattice_site_var = nf * 2 * l_tmp->num_parent_eig_vect;
      
      l_tmp->inner_vector_size = l_tmp->num_inner_lattice_sites * l_tmp->num_lattice_site_var;
      
      l_tmp->vector_size = l_tmp->num_lattice_sites * l_tmp->num_lattice_site_var;
      l_tmp->schwarz_vector_size = 2*l_tmp->vector_size - l_tmp->inner_vector_size;
      
      if(l_tmp->depth == 0) {
        g.p.v_end = l_tmp->inner_vector_size;
        g.p_MP.sp.v_end = l_tmp->inner_vector_size;
        g.p_MP.dp.v_end = l_tmp->inner_vector_size;
      }
      
      if ( g.mixed_precision ) {
        l_tmp->s_float.block_vector_size = l_tmp->s_float.num_block_sites*l_tmp->num_lattice_site_var;
        l_tmp->p_float.v_end = l_tmp->inner_vector_size;
        l_tmp->sp_float.v_end = l_tmp->inner_vector_size;
        l_tmp->dummy_p_float.v_end = l_tmp->inner_vector_size;
        if ( (g.method >= 4 && g.odd_even) || (!l_tmp->idle && l_tmp->level == 0 && g.odd_even) ) {
          if ( l_tmp->level == 0 )
            l_tmp->p_float.v_end = l_tmp->oe_op_float.num_even_sites*l_tmp->num_lattice_site_var;
          else
            l_tmp->sp_float.v_end = l_tmp->oe_op_float.num_even_sites*l_tmp->num_lattice_site_var;
        }
        
      } else {
        l_tmp->s_double.block_vector_size = l_tmp->s_double.num_block_sites*l_tmp->num_lattice_site_var;
        l_tmp->p_double.v_end = l_tmp->inner_vector_size;
        l_tmp->sp_double.v_end = l_tmp->inner_vector_size;
        l_tmp->dummy_p_double.v_end = l_tmp->inner_vector_size;
        if ( (g.method >= 4 && g.odd_even) || (!l_tmp->idle && l_tmp->level == 0 && g.odd_even) ) {
          if ( l_tmp->level == 0 )
            l_tmp->p_double.v_end = l_tmp->oe_op_double.num_even_sites*l_tmp->num_lattice_site_var;
          else
            l_tmp->sp_double.v_end = l_tmp->oe_op_double.num_even_sites*l_tmp->num_lattice_site_var;
        } 
      }
      
      if ( l->level == 0 || l_tmp->next_level == NULL )
        break;
      
      l_tmp = l_tmp->next_level;
    }
    
    update_threading( no_threading, l);
    END_LOCKED_MASTER(threading)
  }

  update_threading( threading, l);
#else
  ASSERT(nf==1);
#endif

}

void define_eot( int *eot, int *N, level_struct *l ) {
  
  int i, mu, t, z, y, x, ls[4], le[4], oe_offset=0;
      
  for ( mu=0; mu<4; mu++ )
    oe_offset += (l->local_lattice[mu]*(g.my_coords[mu]/l->comm_offset[mu]))%2;
  oe_offset = oe_offset%2;
  
  for ( mu=0; mu<4; mu++ ) {
    ls[mu] = 0;
    le[mu] = ls[mu] + l->local_lattice[mu];
  }
  
  // eot: local_lex_index -> iter_index
  i = 0; 
  // visit even sites in the local lexicographical order
  for ( t=0; t<le[T]; t++ )
    for ( z=0; z<le[Z]; z++ )
      for ( y=0; y<le[Y]; y++ )
        for ( x=0; x<le[X]; x++ )
          if ( (t+z+y+x+oe_offset)%2 == 0 ) {
            eot[ lex_index( t, z, y, x, N ) ] = i;
            i++;
          }
  // visit odd sites in the local lexicographical order
  for ( t=0; t<le[T]; t++ )
    for ( z=0; z<le[Z]; z++ )
      for ( y=0; y<le[Y]; y++ )
        for ( x=0; x<le[X]; x++ )
          if ( (t+z+y+x+oe_offset)%2 == 1 ) {
            eot[ lex_index( t, z, y, x, N ) ] = i;
            i++;
          }
  // visit outer boundary sites in pos dir
  for ( mu=0; mu<4; mu++ ) {
    ls[mu] = le[mu];
    le[mu]++;
    // first go over even sites
    for ( t=ls[T]; t<le[T]; t++ )
      for ( z=ls[Z]; z<le[Z]; z++ )
        for ( y=ls[Y]; y<le[Y]; y++ )
          for ( x=ls[X]; x<le[X]; x++ )
            if ( (t+z+y+x+oe_offset)%2 == 0 ) {
              eot[ lex_index( t, z, y, x, N ) ] = i;
              i++;
            }
    // then go over odd sites
    for ( t=ls[T]; t<le[T]; t++ )
      for ( z=ls[Z]; z<le[Z]; z++ )
        for ( y=ls[Y]; y<le[Y]; y++ )
          for ( x=ls[X]; x<le[X]; x++ )
            if ( (t+z+y+x+oe_offset)%2 == 1 ) {
              eot[ lex_index( t, z, y, x, N ) ] = i;
              i++;
            }
    
    ls[mu] = 0;
    le[mu]--;    
  }
}


void define_eo_bt( int **bt, int *eot, int *n_ebs, int *n_obs, int *n_bs, int *N, level_struct *l ) {
  
  int i, t, z, y, x, mu, nu, le[4], bs, oe_offset=0, *bt_mu;
  
  for ( mu=0; mu<4; mu++ ) {
    le[mu] = l->local_lattice[mu];
  }
  
  for ( mu=0; mu<4; mu++ )
    oe_offset += (l->local_lattice[mu]*(g.my_coords[mu]/l->comm_offset[mu]))%2;
  oe_offset = oe_offset%2;//???????
  
  for ( mu=0; mu<4; mu++ ) {
    bt_mu = bt[2*mu];
    bs = 1; 
    le[mu] = 1;
    for ( nu=0; nu<4; nu++ )
      bs *= le[nu]; // size (#sites) in the boundary in the mu dir
     
    // define negative boundary table (bt[2*mu]: iter_index -> eot_index)                                                                                       
    //   go over each negative boundaries (mu:0->3)                                                                                                                                
    //   i: counts #even boundary sites in the mu dir
    i = 0; 
    for ( t=0; t<le[T]; t++ )
      for ( z=0; z<le[Z]; z++ )
        for ( y=0; y<le[Y]; y++ )
          for ( x=0; x<le[X]; x++ )
            if ( (t+z+y+x+oe_offset)%2 == 0 ) {
              bt_mu[i] = site_index( t, z, y, x, N, eot );
              i++;
            }
    n_ebs[2*mu] = i;   // #even boundary sites in the pos mu dir
    n_ebs[2*mu+1] = i; // #even boundary sites in the neg mu dir
    
    // define negative boundary table for communication (bt[2*mu+1]: iter_index -> eot_index)                                                                                       
    //   go over each negative boundaries (mu:0->3)                                                                                                                                
    //   i: counts neg inner boundary sites   
    for ( t=0; t<le[T]; t++ )
      for ( z=0; z<le[Z]; z++ )
        for ( y=0; y<le[Y]; y++ )
          for ( x=0; x<le[X]; x++ )
            if ( (t+z+y+x+oe_offset)%2 == 1 ) {
              bt_mu[i] = site_index( t, z, y, x, N, eot );
              i++;
            }
            
    n_obs[2*mu] = i - n_ebs[2*mu];     // #odd boundary sites in the pos mu dir 
    n_obs[2*mu+1] = i - n_ebs[2*mu+1]; // #odd boundary sites in the neg mu dir
    n_bs[2*mu] = i;                    // #total boundary sites in the pos mu dir
    n_bs[2*mu+1] = i;                  // #total boundary sites in the neg mu dir
    le[mu] = l->local_lattice[mu];
  }
}


void define_nt_bt_tt( int *nt, int *backward_nt, int **bt, int *tt, int *it, int *dt, level_struct *l ) {

/*********************************************************************************
* Defines neighbor table (for the application of the entire operator), negative 
* inner boundary table (for communication) and translation table (for translation 
* to lexicographical site ordnering).
* - int *nt: neighbor table 
* - int **bt: boundary table
* - int *tt: translation table
* - int *it: index table
* - int *dt: dimension table
*********************************************************************************/
  
  ASSERT( dt != NULL && it != NULL );
  
  int i, mu, pos, t, z, y, x, ls[4], le[4], l_st[4], l_en[4],
      offset, stride, *bt_mu, *gs = l->global_splitting;
  
  for ( mu=0; mu<4; mu++ ) {
    ls[mu] = 0;
    le[mu] = ls[mu] + l->local_lattice[mu];
    l_st[mu] = ls[mu];
    l_en[mu] = le[mu];
  }
  
  // define neighbor table for each site in the local lattice (nt: it_index+mu -> it_index of neighbor in the pos mu dir; ignoring offset)
  //   stride: for each site, we have four neighbors in pos dirs (mu)
  //   offset: but if depth>0, we also store iter_index of the site as well at the beggining
  stride = (l->depth==0)?4:5; offset = (l->depth==0)?0:1;
  for ( t=ls[T]; t<le[T]; t++ )
    for ( z=ls[Z]; z<le[Z]; z++ )
      for ( y=ls[Y]; y<le[Y]; y++ )
        for ( x=ls[X]; x<le[X]; x++ ) {
          pos = site_index( t, z, y, x, dt, it );
          if ( offset )
            nt[stride*pos  ] = pos;
          nt[stride*pos+offset+T] = site_index( (gs[T]>1)?t+1:(t+1)%le[T], z, y, x, dt, it ); // T dir
          nt[stride*pos+offset+Z] = site_index( t, (gs[Z]>1)?z+1:(z+1)%le[Z], y, x, dt, it ); // Z dir
          nt[stride*pos+offset+Y] = site_index( t, z, (gs[Y]>1)?y+1:(y+1)%le[Y], x, dt, it ); // Y dir
          nt[stride*pos+offset+X] = site_index( t, z, y, (gs[X]>1)?x+1:(x+1)%le[X], dt, it ); // X dir
        }

  // define backward neighbor table for each site in the local lattice (nt: it_index+mu -> it_index of neighbor in the neg mu dir; ignoring offset)
  //   stride: for each site, we have four neighbors in pos dirs (mu)  
  //   offset: but if depth>0, we also store iter_index of the site as well at the beggining 
  stride = (l->depth==0)?4:5; offset = (l->depth==0)?0:1;
  for ( t=ls[T]; t<le[T]; t++ )
    for ( z=ls[Z]; z<le[Z]; z++ )
      for ( y=ls[Y]; y<le[Y]; y++ )
        for ( x=ls[X]; x<le[X]; x++ ) {
          pos = site_index( t, z, y, x, dt, it );
          if ( offset )
            backward_nt[stride*pos  ] = pos;
          backward_nt[stride*pos+offset+T] = site_index( (gs[T]>1)?(t-1+dt[T])%dt[T]:(t-1+le[T])%le[T], z, y, x, dt, it ); // T dir
          backward_nt[stride*pos+offset+Z] = site_index( t, (gs[Z]>1)?(z-1+dt[Z])%dt[Z]:(z-1+le[Z])%le[Z], y, x, dt, it ); // Z dir
          backward_nt[stride*pos+offset+Y] = site_index( t, z, (gs[Y]>1)?(y-1+dt[Y])%dt[Y]:(y-1+le[Y])%le[Y], x, dt, it ); // Y dir
          backward_nt[stride*pos+offset+X] = site_index( t, z, y, (gs[X]>1)?(x-1+dt[X])%dt[X]:(x-1+le[X])%le[X], dt, it ); // X dir
        }
        
  if ( bt != NULL ) {
    for ( mu=0; mu<4; mu++ ) {
      // define negative boundary table for communication (bt[2*mu+1]: iter_index -> it_index)
      //   go over each negative boundaries (mu:0->3)
      //   i: counts neg inner boundary sites
      l_en[mu] = l_st[mu]+1;
      bt_mu = bt[2*mu+1];
      i = 0;
      for ( t=l_st[T]; t<l_en[T]; t++ )
        for ( z=l_st[Z]; z<l_en[Z]; z++ )
          for ( y=l_st[Y]; y<l_en[Y]; y++ )
            for ( x=l_st[X]; x<l_en[X]; x++ ) {
              bt_mu[i] = site_index( t, z, y, x, dt, it );
              i++;
            }
      l_en[mu] = le[mu];
      
      // define positive boundary table for communication (if desired) (bt[2*mu]: iter_index -> it_index)
      //   go over each positive boundaries (mu:0->3)
      //   i: counts pos inner boundary sites 
      if ( bt[2*mu] != bt[2*mu+1] ) {
        l_st[mu] = le[mu]-1;
        bt_mu = bt[2*mu];
        i = 0;
        for ( t=l_st[T]; t<l_en[T]; t++ )
          for ( z=l_st[Z]; z<l_en[Z]; z++ )
            for ( y=l_st[Y]; y<l_en[Y]; y++ )
              for ( x=l_st[X]; x<l_en[X]; x++ ) {
                bt_mu[i] = site_index( t, z, y, x, dt, it );
                i++;
              }
              l_st[mu] = ls[mu];
      }
    }
  }
  
  // define layout translation table (tt: local_lex_index -> it_index)
  if ( tt != NULL ) {
    i = 0;
    for ( t=0; t<l->local_lattice[T]; t++ )
      for ( z=0; z<l->local_lattice[Z]; z++ )
        for ( y=0; y<l->local_lattice[Y]; y++ )
          for ( x=0; x<l->local_lattice[X]; x++ ) {
            tt[i] = site_index( t, z, y, x, dt, it );
            i++;
          }
  }
}
