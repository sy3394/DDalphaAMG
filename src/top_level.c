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
 * checked: 12/05/2019
 * 1st cleanup:12/19/2019
 */

#include "main.h"

static int wilson_driver( vector_double *solution, vector_double *source, level_struct *l, struct Thread *threading );
static void solve( vector_double *solution, vector_double *source, level_struct *l, struct Thread *threading );


void solve_driver( level_struct *l, struct Thread *threading ) {
  
  vector_double solution, source;
  double minus_twisted_bc[4], norm[g.num_rhs_vect];

  vector_double_init( &solution ); solution.num_vect_now = g.num_rhs_vect;//!!!!!!!
  vector_double_init( &source ); source.num_vect_now = g.num_rhs_vect;//!!!!!!!
  vector_double_alloc( &solution, _INNER, g.num_rhs_vect, l, threading );
  vector_double_alloc( &source, _INNER, g.num_rhs_vect, l, threading );
  rhs_define( &source, l, threading ); 
  //solution.num_vect_now = g.num_rhs_vect;
 
  if(g.bc==2)
    for ( int i=0; i<4; i++ )
      minus_twisted_bc[i] = -1*g.twisted_bc[i];
  
#ifdef HAVE_TM1p1
  if( g.epsbar != 0 || g.epsbar_ig5_odd_shift != 0 || g.epsbar_ig5_odd_shift != 0 ) { 
    data_layout_n_flavours( 2, l, threading );
    printf0("inverting doublet operator\n");
  }
#endif

  START_LOCKED_MASTER(threading)//my addition
  if(g.bc==2)
      apply_twisted_bc_to_vector_double_new( &source, &source, g.twisted_bc, l);
  END_LOCKED_MASTER(threading)
  global_norm_double_new( norm, &source, 0, l->inner_vector_size, l, threading );

  START_MASTER(threading)
  for( int i=0; i<g.num_rhs_vect; i++ ){
    //norm = global_norm_double( &source, 0, l->inner_vector_size, l, threading );
    printf0("source vector %d norm: %le\n",i,norm[i]);
  }
  END_MASTER(threading)
  SYNC_MASTER_TO_ALL(threading)

#ifdef HAVE_TM1p1
  if( g.n_flavours == 1 )
#endif
#ifdef HAVE_TM
    if ( g.mu + g.mu_odd_shift != 0.0 || g.mu + g.mu_even_shift != 0.0 )
      if(g.downprop) {
	SYNC_MASTER_TO_ALL(threading)      
	START_MASTER(threading)  
	  printf0("\n\n+--------------------------- up ---------------------------+\n\n");
	END_MASTER(threading)

	solve( &solution, &source, l, threading );    
      
	if(g.bc==2)
	  apply_twisted_bc_to_vector_double_new( &solution, &solution, minus_twisted_bc, l);
      
	START_LOCKED_MASTER(threading)  
	  printf0("\n\n+-------------------------- down --------------------------+\n\n");
	  g.mu*=-1;
	  g.mu_odd_shift*=-1;
	  g.mu_even_shift*=-1;
	END_LOCKED_MASTER(threading)
  
	tm_term_update( g.mu, l, threading );
        finalize_operator_update( l, threading );
      } 
#endif

  solve( &solution, &source, l, threading );

  if(g.bc==2)
    apply_twisted_bc_to_vector_double_new( &solution, &solution, minus_twisted_bc, l);
 
  global_norm_double_new( norm, &solution, 0, l->inner_vector_size, l, threading );
  for( int i=0; i<g.num_rhs_vect; i++ ){
    //norm = global_norm_double( &solution, 0, l->inner_vector_size, l, threading );
    printf0("solution vector %d norm: %le\n",i,norm[i]);
  }

  if (g.my_rank == 0){//need to update it
  START_LOCKED_MASTER(threading)//my addition     
    //    printf0("sol %d %d",l->inner_vector_size,solution.size);
  vector_double_change_layout( &solution, &solution, _NVEC_OUTER, no_threading );
  vector_double_change_layout( &source, &source, _NVEC_OUTER, no_threading );
  FILE *f;
  f = fopen("DDsol.out","w");
  for ( int jj=0;jj<g.num_rhs_vect; jj++)
    for ( int jjj=0; jjj<solution.size; jjj++)
      fprintf(f,"vec_%d[%d] = %g %g\n", jj,jjj,creal_double(solution.vector_buffer[jj*solution.size+jjj]),cimag_double(solution.vector_buffer[jj*solution.size+jjj]));
  fclose(f);
  END_LOCKED_MASTER(threading)
    }
  vector_double_free( &solution, l, threading );
  vector_double_free( &source, l, threading );

#ifdef HAVE_TM1p1
  if( g.n_flavours == 2 ) 
    data_layout_n_flavours( 1, l, threading );
#endif
}

static void solve( vector_double *solution, vector_double *source, level_struct *l, struct Thread *threading ) {
  
  if ( g.vt.evaluation ) {//?????
    vector_double rhs = g.mixed_precision==2?g.p_MP.dp.b:g.p.b;
    // this would yield different results if we threaded it, so we don't
    START_LOCKED_MASTER(threading)
    //vector_double_define_random( &rhs, 0, l->inner_vector_size, l );
    vector_double_define_random_new( &rhs, 0, l->inner_vector_size, l ); rhs.num_vect_now = g.num_rhs_vect;
    scan_var( &(g.vt), l );
    END_LOCKED_MASTER(threading)
  } else {
    wilson_driver( solution, source, l, threading );
  }
}

void rhs_define( vector_double *rhs, level_struct *l, struct Thread *threading ) {
  
  // no hyperthreading here
  if(threading->thread != 0)
    return;

  int j, jj;
  int start = threading->start_index[l->depth];
  int end = threading->end_index[l->depth];

  rhs->num_vect_now = g.num_rhs_vect;
  if ( g.rhs == 0 ) {
    //vector_double_define( rhs, 1, start, end, l );
    vector_double_define_new( rhs, 1, start, end, l );
    START_MASTER(threading)
    if ( g.print > 0 ) printf0("rhs = ones\n");
    END_MASTER(threading)
  } else if ( g.rhs == 1 )  {
    //vector_double_define( rhs, 0, start, end, l );
    vector_double_define_new( rhs, 0, start, end, l );
    if ( g.my_rank == 0 ) {
      START_LOCKED_MASTER(threading)
      //rhs->vector_buffer[0] = 1.0;
      //for ( int i=0; i<rhs->num_vect; i++ )
        //rhs->vector_buffer[i*(rhs->size)] = 1.0;
      VECTOR_LOOP( j, rhs->num_vect_now, jj, rhs->vector_buffer[j+jj] = 1.0; )
      END_LOCKED_MASTER(threading)
    }
    START_MASTER(threading)
    if ( g.print > 0 ) printf0("rhs = first unit vector\n");
    END_MASTER(threading)
  } else if ( g.rhs == 2 ) {
    // this would yield different results if we threaded it, so we don't
    START_LOCKED_MASTER(threading)
    //vector_double_define_random( rhs, 0, l->inner_vector_size, l );
    vector_double_define_random_new( rhs, 0, l->inner_vector_size, l );
    END_LOCKED_MASTER(threading)
    START_MASTER(threading)
    if ( g.print > 0 ) printf0("rhs = random\n");
    END_MASTER(threading)
  } else if ( g.rhs == 3 ) {
    //vector_double_define( rhs, 0, start, end, l );
    vector_double_define_new( rhs, 0, start, end, l );
    if ( g.print > 0 ) printf0("rhs = 0's\n");
  } else {
    ASSERT( g.rhs >= 0 && g.rhs <= 4 );
  }
    
}

static int wilson_driver( vector_double *solution, vector_double *source, level_struct *l, struct Thread *threading ) {
  
  int iter = 0, start = threading->start_index[l->depth], end = threading->end_index[l->depth];
  vector_double rhs = (g.mixed_precision==2 && g.method >= 0)?g.p_MP.dp.b:g.p.b; rhs.num_vect_now = g.num_rhs_vect;//!!!!!!! 
  vector_double sol = (g.mixed_precision==2 && g.method >= 0)?g.p_MP.dp.x:g.p.x; sol.num_vect_now = g.num_rhs_vect;//!!!!!!! 

#ifdef WILSON_BENCHMARK
  START_MASTER(threading)
  prof_init( l );
  END_MASTER(threading)
  double t = -MPI_Wtime();
  double t_min = 1000;
  for ( int i=0; i<100; i++ ) {
    double tmp_t = -MPI_Wtime();
#endif
    //vector_double_change_layout( &sol, &sol, _NVEC_INNER, no_threading );
    //vector_double_change_layout( &rhs, &rhs, _NVEC_INNER, no_threading );

    //printf0("wilson_driver\n");//debug!!!!!
    //vector_double_copy( &rhs, source, start, end, l );
  vector_double_copy_new( &rhs, source, start, end, l );
  if ( g.method == -1 ) {
    cgn_double( &(g.p), l, threading );
  } else if ( g.mixed_precision == 2 ) {
    iter = fgmres_MP( &(g.p_MP), l, threading );
  } else {
    //printf("wil: %d %d %d\n",g.p.b.num_vect_now, g.num_vect_now,g.mixed_precision==2 && g.method >= 0);
    iter = fgmres_double( &(g.p), l, threading );
    //printf("wil2: %d %d\n",g.p.b.num_vect_now, g.num_vect_now);
  }
  //vector_double_copy( solution, &sol, start, end, l );
  vector_double_copy_new( solution, &sol, start, end, l );
  
#ifdef WILSON_BENCHMARK
    tmp_t += MPI_Wtime();
    if ( tmp_t < t_min )
      t_min = tmp_t;
  }
  t +=MPI_Wtime();
  START_MASTER(threading)
  printf0("average over 100 solves: %lf seconds\n", t/100 );
  printf0("minimum out of 100 solves: %lf seconds\n", t_min );
  prof_print( l );
  END_MASTER(threading)
#endif
  
    //vector_double_change_layout( &sol, &sol, _NVEC_OUTER, no_threading );
    //vector_double_change_layout( &rhs, &rhs, _NVEC_OUTER, no_threading );

  return iter;
}
