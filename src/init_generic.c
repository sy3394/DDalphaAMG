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
 * checked:11/29/2019
 * changed from sbacchio
 * checked: 12/09/2019
 * 1st cleanup:12/18/2019
 */

#include "main.h"

/*****************  PROFILING  **************************************************/
void prof_PRECISION_init( level_struct *l ) {

/*********************************************************************************
* Initializes the profiling struct by specifying the name of every entry.
*********************************************************************************/
  
  if ( l != NULL ) {
    for ( int i=0; i<_NUM_PROF; i++ ) {
      l->prof_PRECISION.time[i] = 0.0;
      l->prof_PRECISION.count[i] = 0.0;
      l->prof_PRECISION.flop[i] = 0.0;
    }
    
    double level_ratio = 1;
    for ( int mu=0; mu<4; mu++ )
      level_ratio *= (double)g.global_lattice[l->depth][mu]/(double)g.global_lattice[0][mu];
    
    sprintf( l->prof_PRECISION.name[_GIP], "global inner product, PRECISION" );
    l->prof_PRECISION.flop[_GIP] = level_ratio*l->num_lattice_site_var*8.0;
    sprintf( l->prof_PRECISION.name[_PIP], "process inner product, PRECISION" );
    l->prof_PRECISION.flop[_PIP] = level_ratio*l->num_lattice_site_var*8.0;
    sprintf( l->prof_PRECISION.name[_LA2], "2 flop vector operations, PRECISION" );
    l->prof_PRECISION.flop[_LA2] = level_ratio*l->num_lattice_site_var*2.0;
    sprintf( l->prof_PRECISION.name[_LA6], "6 flop vector operations, PRECISION" );
    l->prof_PRECISION.flop[_LA6] = level_ratio*l->num_lattice_site_var*6.0;
    sprintf( l->prof_PRECISION.name[_LA8], "8 flop vector operations, PRECISION" );
    l->prof_PRECISION.flop[_LA8] = level_ratio*l->num_lattice_site_var*8.0;
    sprintf( l->prof_PRECISION.name[_LA], "other vector operations, PRECISION" );
    sprintf( l->prof_PRECISION.name[_GRAM_SCHMIDT], "Gram-Schmidt, PRECISION" );
    sprintf( l->prof_PRECISION.name[_GRAM_SCHMIDT_ON_AGGREGATES], "Gram-Schmidt on aggregates, PRECISION" );
    sprintf( l->prof_PRECISION.name[_CPY], "copy operations, PRECISION" );
    sprintf( l->prof_PRECISION.name[_RS], "real scale operations, PRECISION" );
    sprintf( l->prof_PRECISION.name[_SET], "set value operations, PRECISION" );
    sprintf( l->prof_PRECISION.name[_PR], "interpolation and restriction, PRECISION" );
    l->prof_PRECISION.flop[_PR] = level_ratio*l->num_lattice_site_var*8.0*(l->num_lattice_site_var/2);
    sprintf( l->prof_PRECISION.name[_SC], "self coupling, PRECISION" );
    l->prof_PRECISION.flop[_SC] = (l->depth==0)?552.0:level_ratio*SQUARE(l->num_lattice_site_var)*8.0;
    sprintf( l->prof_PRECISION.name[_NC], "neighbor coupling, PRECISION" );
    l->prof_PRECISION.flop[_NC] = (l->depth==0)?1368.0:level_ratio*8.0*SQUARE(l->num_lattice_site_var)*8.0;
    sprintf( l->prof_PRECISION.name[_SM], "smoother, PRECISION" );
    double ncflops = l->prof_PRECISION.flop[_SC];
    for ( int mu=0; mu<4; mu++ )
      ncflops += (l->prof_PRECISION.flop[_NC]/4.0)*((double)(l->block_lattice[mu]-1)/(double)l->block_lattice[mu]);
    l->prof_PRECISION.flop[_SM] = ncflops * (double)(g.odd_even?l->block_iter+1:l->block_iter);
    l->prof_PRECISION.flop[_SM] += (l->prof_PRECISION.flop[_NC] + l->prof_PRECISION.flop[_SC]);
    sprintf( l->prof_PRECISION.name[_OP_COMM], "operator comm init, PRECISION" );
    sprintf( l->prof_PRECISION.name[_OP_IDLE], "operator comm wait, PRECISION" );
    sprintf( l->prof_PRECISION.name[_ALLR], "allreduces, PRECISION" );
    sprintf( l->prof_PRECISION.name[_GD_COMM], "data re-distribution comm init, PRECISION" );
    sprintf( l->prof_PRECISION.name[_GD_IDLE], "data re-distribution comm wait, PRECISION" );
    sprintf( l->prof_PRECISION.name[_SM1], "smoother - pt 1, res no comm, PRECISION" );
    sprintf( l->prof_PRECISION.name[_SM2], "smoother - pt 2, solve no comm, PRECISION" );
    sprintf( l->prof_PRECISION.name[_SM3], "smoother - pt 3, res comm, PRECISION" );
    sprintf( l->prof_PRECISION.name[_SM4], "smoother - pt 4, solve comm, PRECISION" );
    
    sprintf( l->prof_PRECISION.name[_SMALL1], "Hessenberg: qr update PRECISION" );
    sprintf( l->prof_PRECISION.name[_SMALL2], "Hessenberg: bkwd subst PRECISION" );
  }
}


double prof_PRECISION_print( level_struct *l ) {
  double flop = 0;
  for ( int i=0; i<_NUM_PROF; i++ )
    if ( l->prof_PRECISION.count[i] > 0 ) {
      if ( l->prof_PRECISION.count[i] > 9999999 )
        printf0("| %37s: %8.2le(%7.1le) |\n", l->prof_PRECISION.name[i], l->prof_PRECISION.time[i], l->prof_PRECISION.count[i] );
      else
        printf0("| %37s: %8.2le(%7d) |\n", l->prof_PRECISION.name[i], l->prof_PRECISION.time[i], (int)l->prof_PRECISION.count[i] );
      flop += (double)l->prof_PRECISION.count[i] * l->prof_PRECISION.flop[i];
    }
  return flop;
}

/*********** LEVEL STRUCTURE  ******************************************************************/

void level_PRECISION_init( level_struct *l ) {

  for ( int i=0; i<9; i++ )
    vector_PRECISION_init( &(l->vbuf_PRECISION[i]) );
  for ( int i=0; i<2; i++ )
    vector_PRECISION_init( &(l->vtmp_PRECISION[i]) );
  
  operator_PRECISION_init( &(l->op_PRECISION) );
  operator_PRECISION_init( &(l->oe_op_PRECISION) );
  schwarz_PRECISION_init( &(l->s_PRECISION), l );
  interpolation_PRECISION_struct_init( &(l->is_PRECISION) );
  fgmres_PRECISION_struct_init( &(l->p_PRECISION) );
  fgmres_PRECISION_struct_init( &(l->sp_PRECISION) );
}

void fine_level_PRECISION_alloc( level_struct *l ) {//only p_PRECISION.b/x are allocated here??????
  
  int n = 8;
  int nvecs = g.num_vect_now; //!!!!!!!!!!!
#ifdef HAVE_TM1p1
  for ( int i=0; i<n; i++ ){
    vector_PRECISION_alloc( &(l->vbuf_PRECISION[i]), _ORDINARY, 2*nvecs, l, no_threading );
  }
  vector_PRECISION_alloc( &(l->p_PRECISION.b), _INNER, 2*nvecs, l, no_threading );
  vector_PRECISION_alloc( &(l->p_PRECISION.x), _INNER, 2*nvecs, l, no_threading ); 
#else
  for ( int i=0; i<n; i++ )
    vector_PRECISION_alloc( &(l->vbuf_PRECISION[i]), _ORDINARY, nvecs, l, no_threading );
  for ( int i=0; i<2; i++ )
    vector_PRECISION_alloc( &(l->vtmp_PRECISION[i]), _ORDINARY, nvecs, l, no_threading );
  vector_PRECISION_alloc( &(l->p_PRECISION.b), _INNER, nvecs, l, no_threading );
  vector_PRECISION_alloc( &(l->p_PRECISION.x), _INNER, nvecs, l, no_threading ); 
#endif
}

void fine_level_PRECISION_free( level_struct *l ) {
  
  int n = 8;

  for ( int i=0; i<n; i++ )
    vector_PRECISION_free( &(l->vbuf_PRECISION[i]), l, no_threading );
  for ( int i=0; i<2; i++ )
    vector_PRECISION_free( &(l->vtmp_PRECISION[i]), l, no_threading );
  vector_PRECISION_free( &(l->p_PRECISION.b), l, no_threading );
  vector_PRECISION_free( &(l->p_PRECISION.x), l, no_threading );
}

// allocate memory for gmres_PRECISION_struct p_PRECISION, p_PRECISION.b, p_PRECISION.x
void next_level_PRECISION_setup( level_struct *l ) {

  prof_float_init( l->next_level );
  prof_double_init( l->next_level );
  gathering_PRECISION_next_level_init( &(l->next_level->gs_PRECISION), l );  
  gathering_PRECISION_setup( &(l->next_level->gs_PRECISION), l->next_level );
  
  if ( !l->idle ) {
    g.num_vect_now = ( g.num_rhs_vect < l->next_level->num_eig_vect )? l->next_level->num_eig_vect:g.num_rhs_vect;//??????

    coarsening_index_table_PRECISION_alloc( &(l->is_PRECISION), l );//not next_level??????? although it's about coarsened lattice structure on the next level
    coarsening_index_table_PRECISION_define( &(l->is_PRECISION), &(l->s_PRECISION), l );//not next_level??????? although it's about coarsened lattice structure on the next level 

    if ( l->level == 1 && !l->next_level->idle ) {// the-next-to coarsest level and next level is not idle???
      printf("next_level_PRECISION_setup0: level=%d depth=%d\n",l->level,l->depth);
      fgmres_PRECISION_struct_alloc( g.coarse_iter, g.coarse_restart, _ORDINARY, g.coarse_tol, 
                                     _COARSE_GMRES, _NOTHING, NULL,
				     g.odd_even?coarse_apply_schur_complement_PRECISION_new:apply_coarse_operator_PRECISION_new,//!!!!!!!!
                                     &(l->next_level->p_PRECISION), l->next_level );
    } else {
      if ( g.kcycle ) { //if the current level is not the-next-to coarsest level and uses K-cycle
	printf("next_level_PRECISION_setup kcycle: level=%d depth=%d\n",l->level,l->depth);
        fgmres_PRECISION_struct_alloc( g.kcycle_restart, g.kcycle_max_restart, _ORDINARY, g.kcycle_tol, 
                                       _K_CYCLE, _RIGHT, vcycle_PRECISION_new,//!!!!!!!!!!!!
				       apply_coarse_operator_PRECISION_new,//!!!!!!
                                       &(l->next_level->p_PRECISION), l->next_level );
      } else {//otherwise only p_PRECISION.b/x are used
        vector_PRECISION_init(&(l->next_level->p_PRECISION.b));
        vector_PRECISION_init(&(l->next_level->p_PRECISION.x));
#ifdef HAVE_TM1p1
        vector_PRECISION_alloc( &(l->next_level->p_PRECISION.b), _ORDINARY, 2*g.num_vect_now, l->next_level, no_threading );//!!!!!
        vector_PRECISION_alloc( &(l->next_level->p_PRECISION.x), _ORDINARY, 2*g.num_vect_now, l->next_level, no_threading );//!!!!!!
#else
        vector_PRECISION_alloc( &(l->next_level->p_PRECISION.b), _ORDINARY, g.num_vect_now, l->next_level, no_threading );//!!!!!!!
        vector_PRECISION_alloc( &(l->next_level->p_PRECISION.x), _ORDINARY, g.num_vect_now, l->next_level, no_threading );//!!!!!
#endif
        l->next_level->p_PRECISION.v_start = 0;
        l->next_level->p_PRECISION.v_end = l->next_level->inner_vector_size;
      }
    }

    int i, n = (l->next_level->level>0)?6:4;
    for ( i=0; i<n; i++ ){
#ifdef HAVE_TM1p1
      vector_PRECISION_alloc( &(l->next_level->vbuf_PRECISION[i]), _ORDINARY, 2*l->next_level->num_eig_vect, l->next_level, no_threading );//!!!!!!!!!!!!!
#else
      vector_PRECISION_alloc( &(l->next_level->vbuf_PRECISION[i]), _ORDINARY, l->next_level->num_eig_vect, l->next_level, no_threading );//!!!!!!1
#endif
    }
    for ( int i=0; i<2; i++ )
      vector_PRECISION_alloc( &(l->next_level->vtmp_PRECISION[i]), _ORDINARY, g.num_rhs_vect, l->next_level, no_threading );// don't know where it is used????
  }
}

void next_level_PRECISION_free( level_struct *l ) {
  
  coarse_grid_correction_PRECISION_free( l );

  if ( !l->idle ) {
    if ( ( l->level == 1 && !l->next_level->idle ) || g.kcycle ) {
      fgmres_PRECISION_struct_free( &(l->next_level->p_PRECISION), l->next_level );
    } else {
      vector_PRECISION_free( &(l->next_level->p_PRECISION.b), l->next_level, no_threading );
      vector_PRECISION_free( &(l->next_level->p_PRECISION.x), l->next_level, no_threading );
    }
  
    int i, n = (l->next_level->level>0)?6:4;  
    for ( i=0; i<n; i++)
      vector_PRECISION_free( &(l->next_level->vbuf_PRECISION[i]), l->next_level, no_threading );
    coarsening_index_table_PRECISION_free( &(l->is_PRECISION), l );
  }

  gathering_PRECISION_free( &(l->next_level->gs_PRECISION), l->next_level );
}
