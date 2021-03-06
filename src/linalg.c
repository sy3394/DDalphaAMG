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
 *
 */

#include "main.h"


// results <- (phi's, psi) processwise 
// assume resutls have phi->num_vect_now fields 
void process_multi_inner_product_MP( int count, complex_double *results, vector_float *phi,
					 vector_float *psi, int start, int end, level_struct *l, struct Thread *threading ) {
  /*****************************************
   * Assume: each vector set in phi contains the same #used vectors as in psi
   * Input:
   *  int count: #basis vectors to be multiplied
   *  vector_float *phi: array of basis vectors
   *  vector_float *psi: a vector multiplying each basis vector in the array
   * Output:
   *  complex_double *results: stores the inner products (count x psi->num_vect_now) of each vector in the basis phi and the given vector psi on each process
   *****************************************/

#ifdef DEBUG
#ifdef HAVE_TM1p1
  if ( g.n_flavours == 2 && psi->num_vect_now != 2*num_loop )
    error0("process_multi_inner_product_MP: doublet error\n");
#endif
#endif
  
  int c, i, j, jj, nvec = psi->num_vect_now, nvecmf = psi->num_vect_now;
#ifdef HAVE_TM1p1
  if ( g.n_flavours == 2 && ( g.epsbar != 0 || g.epsbar_ig5_odd_shift != 0 || g.epsbar_ig5_odd_shift != 0 ) )
    nvec /= g.n_flavours;
#endif
  
  int core_start, core_end;
  compute_core_start_end_custom(start, end, &core_start, &core_end, l, threading, 1 );

  int thread = omp_get_thread_num();
  if ( thread == 0 && core_start != core_end)
    PROF_float_START( _PIP, threading );

  VECTOR_LOOP(j, count*nvec, jj, results[j+jj] = 0.0;)

  for ( c=0; c<count; c++ ) {
#ifdef DEBUG
    if ( phi[c].num_vect_now != psi->num_vect_now )
      error0("process_multi_inner_product_MP: phi[%d]->num_vect_now != psi->num_vect_now \n",c);
#endif
    for ( i=core_start; i<core_end; i++ )
      VECTOR_LOOP(j, nvecmf, jj,
		  results[c*nvec+(j+jj)%nvec] += (complex_double) conj_float(phi[c].vector_buffer[i*phi[c].num_vect+j+jj])*psi->vector_buffer[i*psi->num_vect+j+jj])
  }

  START_NO_HYPERTHREADS(threading)
  VECTOR_LOOP( j, count*nvec, jj,((complex_double *)threading->workspace)[threading->core*count*nvec+j+jj] = results[j+jj];)
  END_NO_HYPERTHREADS(threading)
  // master sums up all results
  SYNC_CORES(threading)
  START_MASTER(threading)
  for( i=1; i<threading->n_core; i++)
    VECTOR_LOOP( j, count*nvec, jj,((complex_double *)threading->workspace)[j+jj] += ((complex_double *)threading->workspace)[i*count*nvec+j+jj];)
  END_MASTER(threading)
  // all threads need the result of the norm
  SYNC_MASTER_TO_ALL(threading)
  VECTOR_LOOP( j, count*nvec, jj, results[j+jj] = ((complex_double *)threading->workspace)[j+jj];)

  if(thread == 0 && core_start != core_end)
    PROF_float_STOP( _PIP, (double)(core_end-core_start)/(double)l->inner_vector_size, threading );
}

void global_norm_MP( double *res, vector_float *x, int start, int end, level_struct *l, struct Thread *threading ) {
  /*********************
   * res <- norms of vectors in x: contains x->num_vect_now elements
   ********************/

#ifdef DEBUG
#ifdef HAVE_TM1p1
  if ( g.n_flavours == 2 && x->num_vect_now != 2*num_loop )
    error0("global_norm_MP: doublet error\n");
#endif
#endif
  
  int core_start, core_end;
  compute_core_start_end_custom(start, end, &core_start, &core_end, l, threading, 1);

  int thread = omp_get_thread_num();
  if(thread == 0 && core_start != core_end)
    PROF_float_START( _GIP, threading );

  int i, j, jj, nvec = x->num_vect_now, nvecmf = x->num_vect_now;
#ifdef HAVE_TM1p1
  if ( g.n_flavours == 2 && ( g.epsbar != 0 || g.epsbar_ig5_odd_shift != 0 || g.epsbar_ig5_odd_shift != 0 ) )
    nvec /= g.n_flavours;
#endif
  
  double global_alpha[nvec];
  VECTOR_LOOP(j, nvec, jj, res[j+jj]=0;)

  for( i=core_start; i<core_end; i++ )
    VECTOR_LOOP(j, nvecmf, jj, res[(j+jj)%nvec] += NORM_SQUARE_float(x->vector_buffer[i*x->num_vect+j+jj]);)

  // communication  -------------------------------------
  // sum over cores: be careful about overflow of workspace!!!!!!!! also potential problem using res directly
  START_NO_HYPERTHREADS(threading)
  VECTOR_LOOP(j, nvec, jj, ((double *)threading->workspace)[threading->core*nvec+j+jj] = res[j+jj];)
  END_NO_HYPERTHREADS(threading)
  // master sums up all results
  SYNC_CORES(threading)
  START_MASTER(threading)
  for ( i=1; i<threading->n_core; i++)
    VECTOR_LOOP(j, nvec, jj, ((double *)threading->workspace)[0*nvec+j+jj] += ((double *)threading->workspace)[i*nvec+j+jj];)
  END_MASTER(threading)
  SYNC_MASTER_TO_ALL(threading)
  VECTOR_LOOP(j, nvec, jj, res[j+jj] = ((double *)threading->workspace)[0*nvec+j+jj];)
  if ( g.num_processes > 1 ) {
    START_MASTER(threading)
    PROF_double_START( _ALLR );
    MPI_Allreduce( res, global_alpha, nvec, MPI_double, MPI_SUM, (l->depth==0)?g.comm_cart:l->gs_float.level_comm );
    PROF_double_STOP( _ALLR, 1 );
    VECTOR_LOOP(j, nvec, jj, ((double *)threading->workspace)[0*nvec+j+jj] = global_alpha[j+jj];)
    END_MASTER(threading)
    // all threads need the result of the norm
    SYNC_MASTER_TO_ALL(threading)
    VECTOR_LOOP(j, nvec, jj, res[j+jj] = ((double *)threading->workspace)[0*nvec+j+jj];)
  }
  
  VECTOR_LOOP(j, nvec, jj, res[j+jj] = (double)sqrt((double)res[j+jj]);)
  
  if(thread == 0 && core_start != core_end)
    PROF_float_STOP( _GIP, (double)(core_end-core_start)/(double)l->inner_vector_size, threading );
}
