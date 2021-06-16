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

/********************
 * index table: Given the lexicographic index of a site, it returns the location of the site in an array
 * There are several layouts
 * local lexicogaphic layout (normal): lists sites in the order X->Y->Z->T
 ********************/

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

static void change_to_n_flavours( gmres_MP_struct *p, int nf ) {

  int i, k;
  if ( p->sp.preconditioner != NULL ) {
    if ( p->sp.kind == _RIGHT ) {
      k = p->dp.restart_length+1;
    } else {
      k = 0;
    }
  }

  // single precision
  p->sp.w.num_vect_now = nf*num_loop;
  for ( i=0; i<p->dp.restart_length+1; i++ ) p->sp.V[i].num_vect_now = nf*num_loop;
  if ( p->sp.Z != NULL ) 
    for ( i=0; i<k; i++ ) p->sp.Z[i].num_vect_now = nf*num_loop;

  // double precision
  p->dp.x.num_vect_now = nf*num_loop;
  p->dp.r.num_vect_now = nf*num_loop;
  p->dp.b.num_vect_now = nf*num_loop;

}

/****
 * To reduce memory requirement, except for interpolation vectors, solutions, and rhs,
 * all vector fields have max duplicaiton of 2*num_loop if HAVE_TM1p1 and num_loop otherwise
 * When HAVE_TM1p1, num_vect_now is set dynamically to g.n_flavours*num_loop
 ****/
void data_layout_n_flavours( int nf, level_struct *l, struct Thread *threading ) {

  /*
   * Description: Update inner d.o.f. and vector sizes to account for multiple flavours
   * Note: flavour index is now the next fastest running index (fastest is vector index)
   * Note: g.n_flavours's default value is 1 and changed on demand
   *
   * Case: HAVE_TM1p1
   *   ASSERT(nf==1 or nf==2)
   *   if g.n_flavours == 1: 
   *     nothing needs to be done
   *   if g.n_flavours == 2:
   *     update inner d.o.f. in various structures
   * Case: not HAVE_TM1p1 
   *   n_flavours should be 1
   *
   */
  ASSERT(nf>0);
  ASSERT(l->depth == 0);

#ifdef HAVE_TM1p1
  ASSERT(nf<=2);
  
  if( g.n_flavours == nf )
    return;
  else {
    START_LOCKED_MASTER(threading)
    g.n_flavours = nf;
    if ( g.epsbar != 0 || g.epsbar_ig5_odd_shift != 0 || g.epsbar_ig5_odd_shift != 0 )
      g.num_indep_flav = 1;
    else
      g.num_indep_flav = nf;
    
    // l->gs_PRECISION.buffer.num_vect_now is not used
    // l->gs_PRECISION.transfer_buffer.num_vect_now is tmp field and so set on demand
    // comm_PRECISION_struct does not use num_vect_now
    // op->buffer[*].num_vect now is tmp field and so set on demand in (coarse_)apply_schur_complement_PRECISION
    // only vector_buffer in l->x is used. 
    
    int i, k = 0;
    struct level_struct *l_tmp = l;
    while(1) {
      if (l_tmp->depth == 0 ) {
#ifdef INIT_ONE_PREC
	if ( g.mixed_precision == 2 && g.method >= 0 ) {
#else
	if ( g.method >= 0 ) {
#endif
	  change_to_n_flavours( &(g.p_MP), nf );
#if defined(INIT_ONE_PREC) && (defined (DEBUG) || defined (TEST_VECTOR_ANALYSIS))
	  g.p.b.num_vect_now = nf*num_loop;
	  g.p.x.num_vect_now = nf*num_loop;
#endif
#ifdef INIT_ONE_PREC
	} else
#else
        }
#endif
	  change_to_n_flavours_double( &(g.p), nf, NULL );
      }

      if ( g.mixed_precision ) {
	for ( i = 0; i<5; i++ ) l_tmp->s_float.buf[i].num_vect_now                 = nf*num_loop;
	for ( i = 0; i<4; i++ ) l_tmp->s_float.oe_buf[i].num_vect_now              = nf*num_loop;
	for ( i = 0; i<3; i++ ) l_tmp->s_float.local_minres_buffer[i].num_vect_now = nf*num_loop;
      } else {
	for ( i = 0; i<5; i++ ) l_tmp->s_double.buf[i].num_vect_now                 = nf*num_loop;
	for ( i = 0; i<4; i++ ) l_tmp->s_double.oe_buf[i].num_vect_now              = nf*num_loop;
	for ( i = 0; i<3; i++ ) l_tmp->s_double.local_minres_buffer[i].num_vect_now = nf*num_loop;
      }

      if ( g.mixed_precision ) {
	change_to_n_flavours_float( &(l_tmp->p_float), nf, l_tmp );
	if ( g.method > 3 ) 
	  change_to_n_flavours_float( &(l_tmp->sp_float), nf, l_tmp );
      } else {
	change_to_n_flavours_double( &(l_tmp->p_double), nf, l_tmp );
        if ( g.method > 3 )
          change_to_n_flavours_double( &(l_tmp->sp_double), nf, l_tmp );
      }

      //int n = (g.method != -1)?2:4;//does not hurt to change num_vect_now of unused buffers for vbuf_PRECISION
      if ( g.mixed_precision ) {
	for ( i = 0; i<5; i++ ) l_tmp->vbuf_float[i].num_vect_now = nf*num_loop;
	for ( i = 0; i<2; i++ ) l_tmp->sbuf_float[i].num_vect_now = nf*num_loop;
      } else {
	for ( i = 0; i<5; i++ ) l_tmp->vbuf_double[i].num_vect_now = nf*num_loop;
	for ( i = 0; i<2; i++ ) l_tmp->sbuf_double[i].num_vect_now = nf*num_loop;
      }
      
      if ( l->level == 0 || l_tmp->next_level == NULL )
        break;
      
      l_tmp = l_tmp->next_level;
    }
    END_LOCKED_MASTER(threading)
  }
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
