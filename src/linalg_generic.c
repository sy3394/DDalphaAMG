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
 * first check:12/03/2019
 * checked: 12/07/2019
 * 1st cleanip: 12/18/2019
 */

#include "main.h"

void global_norm_PRECISION_new( PRECISION *res, vector_PRECISION *x, int start, int end, level_struct *l, struct Thread *threading ) {

  int i, j, jj, nvec = x->num_vect_now;//!!!!!!!!!!!  
  int thread = omp_get_thread_num();
  PRECISION global_alpha[nvec];

  if ( thread == 0 && start != end )
    PROF_PRECISION_START( _GIP, threading );

  int core_start, core_end;
  compute_core_start_end(start, end, &core_start, &core_end, l, threading); //thread == core here: distribute entries from start to end among cores
  printf("global_norm_PRECISION_ %d: %d %d %d %d %d %d %d %d\n", nvec,start, core_start,end, core_end,thread, threading->n_core,threading->core,g.num_processes);
  VECTOR_LOOP(j, nvec, jj, res[j+jj]=0;)//printf("norm res %g\n",res[j+jj]);)
 
  for( i=core_start; i<core_end; i++)
    VECTOR_LOOP(j, nvec, jj, res[j+jj] += NORM_SQUARE_PRECISION(x->vector_buffer[i*x->num_vect+j+jj]);)
      //      for( i=0; i<nvec; i++ ) printf("gnorm_in: %g\n",res[i]);
      //      printf("global_norm_PRECISION: %d %d %d\n",threading->n_core,nvec,end);
  /////////////// communication -------------------------------------
  // sum over cores: be careful about overflow of workspace!!!!!!!! also potential problem using res directly
  START_NO_HYPERTHREADS(threading)
    //    printf("global_norm_PRECISION: %d\n",threading->core);
      VECTOR_LOOP(j, nvec, jj, ((PRECISION *)threading->workspace)[threading->core*nvec+j+jj] = res[j+jj];) //for( i=0; i<nvec; i++ ) printf("gnorm_in: %g\n",((PRECISION *)threading->workspace)[i]);
  END_NO_HYPERTHREADS(threading)
  // master sums up all results
  SYNC_CORES(threading)
  START_MASTER(threading)
  for ( i=1; i<threading->n_core; i++)
    VECTOR_LOOP(j, nvec, jj, ((PRECISION *)threading->workspace)[0*nvec+j+jj] += ((PRECISION *)threading->workspace)[i*nvec+j+jj];)
  END_MASTER(threading)
  SYNC_MASTER_TO_ALL(threading)
  VECTOR_LOOP(j, nvec, jj, res[j+jj] = ((PRECISION *)threading->workspace)[0*nvec+j+jj];)
  if ( g.num_processes > 1 ) {
    START_MASTER(threading)
    PROF_PRECISION_START( _ALLR );
    MPI_Allreduce( res, global_alpha, nvec, MPI_PRECISION, MPI_SUM, (l->depth==0)?g.comm_cart:l->gs_PRECISION.level_comm );
    PROF_PRECISION_STOP( _ALLR, 1 );
    VECTOR_LOOP(j, nvec, jj, ((PRECISION *)threading->workspace)[0*nvec+j+jj] = global_alpha[j+jj];)
    END_MASTER(threading)
    // all threads need the result of the norm
    SYNC_MASTER_TO_ALL(threading)
    VECTOR_LOOP(j, nvec, jj, res[j+jj] = ((PRECISION *)threading->workspace)[0*nvec+j+jj];)
  }

  // take square root
  VECTOR_LOOP(j, nvec, jj, res[j+jj] = (PRECISION)sqrt((double)res[j+jj]);)

  if(thread == 0 && start != end)
    PROF_PRECISION_STOP( _GIP, (double)(end-start)/(double)l->inner_vector_size, threading );
}

void global_inner_product_PRECISION_new( complex_PRECISION *results, vector_PRECISION *phi, vector_PRECISION *psi, int start, int end, level_struct *l, struct Thread *threading ) {
  
  if ( phi->num_vect_now != psi->num_vect_now )
    error0("global_inner_product_PRECISION: phi->num_vect_now != psi->num_vect_now \n");

  PROF_PRECISION_START( _GIP, threading );
  int i, j, jj, nvec = phi->num_vect_now;
  complex_PRECISION local_alpha[nvec], global_alpha[nvec];
 
  VECTOR_LOOP( j, nvec, jj, local_alpha[j+jj]=0; global_alpha[j+jj]=0; )

  int thread_start;
  int thread_end;
  compute_core_start_end(start, end, &thread_start, &thread_end, l, threading);
  
  SYNC_CORES(threading)
  
  for( i=thread_start; i<thread_end; i++){
    VECTOR_LOOP( j, nvec, jj, local_alpha[j+jj] += conj_PRECISION(phi->vector_buffer[i*phi->num_vect+j+jj])*psi->vector_buffer[i*psi->num_vect+j+jj]; )
  }
  
  /////////////// communication -----------------------------------
  // sum over cores
  START_NO_HYPERTHREADS(threading)
  VECTOR_LOOP( j, nvec, jj, ((complex_PRECISION *)threading->workspace)[threading->core*nvec+j+jj] = local_alpha[j+jj];)
  END_NO_HYPERTHREADS(threading)
  // master sums up all results
  SYNC_CORES(threading)
  START_MASTER(threading)
  for( i=1; i<threading->n_core; i++)
    VECTOR_LOOP( j, nvec, jj, ((complex_PRECISION *)threading->workspace)[0*nvec+j+jj] += ((complex_PRECISION *)threading->workspace)[i*nvec+j+jj];)
  END_MASTER(threading)
  SYNC_MASTER_TO_ALL(threading)
  VECTOR_LOOP( j, nvec, jj, local_alpha[j+jj] = ((complex_PRECISION *)threading->workspace)[0*nvec+j+jj];)
  
  if ( g.num_processes > 1 ) {
    START_MASTER(threading)
    PROF_PRECISION_START( _ALLR );
    MPI_Allreduce( local_alpha, global_alpha, nvec, MPI_COMPLEX_PRECISION, MPI_SUM, (l->depth==0)?g.comm_cart:l->gs_PRECISION.level_comm );
    PROF_PRECISION_STOP( _ALLR, 1 );
    VECTOR_LOOP( j, nvec, jj, ((complex_PRECISION *)threading->workspace)[0*nvec+j+jj] = global_alpha[j+jj];)
    END_MASTER(threading)
    // all threads need the result of the norm
    SYNC_MASTER_TO_ALL(threading)
    VECTOR_LOOP( j, nvec, jj, global_alpha[j+jj] = ((complex_PRECISION *)threading->workspace)[0*nvec+j+jj];)
  }
  PROF_PRECISION_STOP( _GIP, (double)(end-start)/(double)l->inner_vector_size, threading );  

  if ( g.num_processes > 1 ) {//we are able to shorten using results directly in local_alpha
    VECTOR_LOOP( j, nvec, jj, results[j+jj] = global_alpha[j+jj]; )
  } else {
    VECTOR_LOOP( j, nvec, jj, results[j+jj] = local_alpha[j+jj]; )
  }
}

// results <- (phi's, psi) processwise
// assume resutls have phi->num_vect_now fields
void process_multi_inner_product_PRECISION_new( int count, complex_PRECISION *results, vector_PRECISION *phi, vector_PRECISION *psi,
						level_struct *l, struct Thread *threading ) {
  /******************************************
   * Input:
   *  int count: #basis vectors to be multiplied
   *  vector_PRECISION *phi: array of basis vectors
   *  vector_PRECISION *psi: a vector multiplying each basis vector in the array
   * Output:
   *  complex_PRECISION *results: stores the inner products of each vector in the basis phi and the given vector psi on each process
   *****************************************/

  int start, end;
  compute_core_start_end(0, psi->size, &start, &end, l, threading);
  int thread = omp_get_thread_num();
  if(thread == 0 && start != end)
    PROF_PRECISION_START( _PIP, threading );
  //  printf("process_multi_inner_product_PRECISION: num-now=%d\n",psi->num_vect_now);
  int i, j, jj, nvec = psi->num_vect_now;//change!!!!!!!!
  VECTOR_LOOP(j, count*nvec, jj, results[j+jj] = 0.0;)

  for(int c=0; c<count; c++){
    if ( phi[c].num_vect_now != psi->num_vect_now )
      error0("process_multi_inner_product_PRECISION: phi->num_vect_now != psi->num_vect_now \n");
    for ( i=start; i<end; i++ )
      VECTOR_LOOP(j, nvec, jj, results[c*nvec+j+jj] += conj_PRECISION(phi[c].vector_buffer[i*phi[c].num_vect+j+jj])*psi->vector_buffer[i*psi->num_vect+j+jj];)//printf("%g ",creal(psi->vector_buffer[i*psi->num_vect+j+jj]));)
  }

  if(thread == 0 && start != end)
    PROF_PRECISION_STOP( _PIP, (double)(end-start)/(double)l->inner_vector_size, threading );
}

// res <- <phi,psi>/<phi,phi> with each vector sliced from start to end
void local_xy_over_xx_PRECISION_new( complex_PRECISION *res, vector_PRECISION *phi, vector_PRECISION *psi, int start, int end, level_struct *l  ) {
  
  if ( phi->num_vect_now != psi->num_vect_now )
    error0("local_xy_over_xx_PRECISION: phi->num_vect_now != psi->num_vect_now \n");

  int nvec = phi->num_vect_now, i, j, jj;
  complex_PRECISION numerator[nvec]; PRECISION denominator[nvec], total_den=0.0;
  VECTOR_LOOP(j, nvec, jj, numerator[j+jj]=0.0; denominator[j+jj]=0.0;)
  
  for( i=start; i<end; i++ )
    VECTOR_LOOP(j, nvec, jj, numerator[j+jj] += conj_PRECISION(phi->vector_buffer[i*phi->num_vect+j+jj])*psi->vector_buffer[i*psi->num_vect+j+jj]; denominator[j+jj] += NORM_SQUARE_PRECISION(phi->vector_buffer[i*phi->num_vect+j+jj]);)

  VECTOR_LOOP(j, nvec, jj, total_den += denominator[j+jj];)
  if ( abs_PRECISION(total_den) < nvec*EPS_PRECISION ) {
    VECTOR_LOOP(j, nvec, jj, res[j+jj] = 0.0;)
    return;
  }
  
  VECTOR_LOOP(j, nvec, jj, res[j+jj] = numerator[j+jj]/denominator[j+jj];)
}

void gram_schmidt_PRECISION_new( vector_PRECISION *V, const int nvec, level_struct *l, struct Thread *threading ) {
  /*****************************
   * a set of vectors in V are orthonormalized
   * nvec: #first vectors in the set to be orthonormalized
   *****************************/

  // NOTE: only thread safe, if "buffer" is the same buffer for all threads belonging to a common MPI process
  START_MASTER(threading)
  PROF_PRECISION_START( _LA );
  END_MASTER(threading)
  SYNC_CORES(threading)

  int i, j, k, start, end;  
  PRECISION norm;
  complex_PRECISION local_alpha, global_alpha; 
  
  compute_core_start_end_custom( 0, l->inner_vector_size, &start, &end, l, threading, l->num_lattice_site_var );//???????

  vector_PRECISION_change_layout( V, V, _NVEC_INNER, no_threading );//just to be sure
  for ( i=0; i<nvec; i++ ) {
    for ( j=0; j<i; j++ ) {
      //take inner product <i,j>
      SYNC_CORES(threading)
      local_alpha = 0;
      for( k=start; k<end; k++) // for each core: do we need to divide into two parts like in aggregate version????????
	local_alpha += conj_PRECISION(V->vector_buffer[k*V->num_vect+j])*V->vector_buffer[k*V->num_vect+i];//is this inefficient?????
      // communication  -----------------------------------
      // sum over cores
      START_NO_HYPERTHREADS(threading)
      ((complex_PRECISION *)threading->workspace)[threading->core] = local_alpha;
      END_NO_HYPERTHREADS(threading)
      // master sums up all results
      SYNC_CORES(threading)
      START_MASTER(threading)
      for( k=1; k<threading->n_core; k++)
	((complex_PRECISION *)threading->workspace)[0] += ((complex_PRECISION *)threading->workspace)[k];
      END_MASTER(threading)
      SYNC_MASTER_TO_ALL(threading)
      local_alpha = ((complex_PRECISION *)threading->workspace)[0];
      if ( g.num_processes > 1 ) {
	START_MASTER(threading)
	PROF_PRECISION_START( _ALLR );
	MPI_Allreduce( &local_alpha, &global_alpha, 1, MPI_COMPLEX_PRECISION, MPI_SUM, (l->depth==0)?g.comm_cart:l->gs_PRECISION.level_comm );
	PROF_PRECISION_STOP( _ALLR, 1 );
	((complex_PRECISION *)threading->workspace)[0] = global_alpha;
	END_MASTER(threading)
	// all threads need the result of the norm
	SYNC_MASTER_TO_ALL(threading)
        global_alpha = ((complex_PRECISION *)threading->workspace)[0];
      }
      else {
	global_alpha = local_alpha;
      }
      
      // v_i -= alpha*v_j
      for( k=start; k<end; k++) // for each core: do we need to divide into two parts like in aggregate version????????
        V->vector_buffer[k*V->num_vect+i] -= global_alpha*V->vector_buffer[k*V->num_vect+j];
    }
    // take squared norm
    local_alpha = 0;
    for( k=start; k<end; k++) // for each core: do we need to divide into two parts like in aggregate version????????
      local_alpha += NORM_SQUARE_PRECISION(V->vector_buffer[k*V->num_vect+i]);
    
    // communication  -----------------------------------
    // sum over cores
    START_NO_HYPERTHREADS(threading)
    ((complex_PRECISION *)threading->workspace)[threading->core] = local_alpha;
    END_NO_HYPERTHREADS(threading)
    // master sums up all results
    SYNC_CORES(threading)
    START_MASTER(threading)
    for( k=1; k<threading->n_core; k++)
      ((complex_PRECISION *)threading->workspace)[0] += ((complex_PRECISION *)threading->workspace)[k];
    END_MASTER(threading)
    SYNC_MASTER_TO_ALL(threading)
    local_alpha = ((complex_PRECISION *)threading->workspace)[0];
    if ( g.num_processes > 1 ) {
      START_MASTER(threading)
      PROF_PRECISION_START( _ALLR );
      MPI_Allreduce( &local_alpha, &global_alpha, 1, MPI_COMPLEX_PRECISION, MPI_SUM, (l->depth==0)?g.comm_cart:l->gs_PRECISION.level_comm );
      PROF_PRECISION_STOP( _ALLR, 1 );
      ((complex_PRECISION *)threading->workspace)[0] = global_alpha;
      END_MASTER(threading)
      // all threads need the result of the norm
      SYNC_MASTER_TO_ALL(threading)
      global_alpha = ((complex_PRECISION *)threading->workspace)[0];
    }
    else {
      global_alpha = local_alpha;
    }
    
    // divide v_i by its norm
    //    global_alpha = 1/sqrt(global_alpha);
    for( k=start; k<end; k++) // for each core: do we need to divide into two parts like in aggregate version????????
      V->vector_buffer[k*V->num_vect+i] /= sqrt(global_alpha);
  }
  //vector_PRECISION_change_layout( V, V, _LV_SV_NV, no_threading );

  START_MASTER(threading)
  PROF_PRECISION_STOP( _LA, 1 );
  END_MASTER(threading)
  SYNC_CORES(threading)
}

void gram_schmidt_on_aggregates_PRECISION_new( vector_PRECISION *phi, const int num_vect, level_struct *l, struct Thread *threading ) {
  /************************************************
   * vector_PRECISION *phi: a set of vectors to be orthogonalized 
   * num_vect: should be equal to phi->num_vect_now
   * note: aggregate lives in the process so that no communication is necessary??????
   ***********************************************/

  PROF_PRECISION_START( _GRAM_SCHMIDT_ON_AGGREGATES, threading );
  SYNC_CORES(threading)
  SYNC_HYPERTHREADS(threading)
    int i, j, k, k1, k2, nvec_phi = phi->num_vect;
  int num_aggregates = l->s_PRECISION.num_aggregates;
  int aggregate_size = l->inner_vector_size / num_aggregates;
  int offset         = l->num_lattice_site_var/2;
      
  complex_PRECISION alpha1, alpha2;
  vector_PRECISION v_pt1, v_pt2;
  PRECISION norm1, norm2;
  int start = threading->start_index[l->depth];
  int end   = threading->end_index[l->depth];

  if ( num_vect != phi->num_vect_now )
    error0("gram_schmidt_on_aggregates_PRECISION: assumptions are not met\n");

  vector_PRECISION_change_layout( phi, phi, _NVEC_OUTER, no_threading );
  for ( j=threading->n_thread*threading->core+threading->thread; j<num_aggregates; j+=threading->n_thread*threading->n_core ) {
    for ( k1=0; k1<num_vect; k1++ ) {
      v_pt1.vector_buffer = phi->vector_buffer + (k1*(phi->size) + j*aggregate_size);
      
      for ( k2=0; k2<k1; k2++ ) {
        v_pt2.vector_buffer = phi->vector_buffer + (k2*(phi->size) + j*aggregate_size);
        alpha1 = 0; alpha2 = 0;
        // V[k1] -= <V[k2],V[k1]> V[k2] | 2*j-th and 2*j+1-st aggregate
        for ( i=0; i<aggregate_size; ) {
          for ( k=0; k<offset; k++, i++ )//?????
            alpha1 += conj_PRECISION(v_pt2.vector_buffer[i]) * v_pt1.vector_buffer[i];
          for ( k=0; k<offset; k++, i++ )
            alpha2 += conj_PRECISION(v_pt2.vector_buffer[i]) * v_pt1.vector_buffer[i];
        }
        for ( i=0; i<aggregate_size; ) {
          for ( k=0; k<offset; k++, i++ )
            v_pt1.vector_buffer[i] -=  alpha1 * v_pt2.vector_buffer[i];
          for ( k=0; k<offset; k++, i++ )
            v_pt1.vector_buffer[i] -=  alpha2 * v_pt2.vector_buffer[i];
        }
      }
      
      norm1 = 0; norm2 = 0;
      // V[k1] = V[k1]/norm(V[k1]) | 2*j-th and 2*j+1-st aggregate    
      for ( i=0; i<aggregate_size; ) {
        for ( k=0; k<offset; k++, i++ )
          norm1 += NORM_SQUARE_PRECISION(v_pt1.vector_buffer[i]);
        for ( k=0; k<offset; k++, i++ )
          norm2 += NORM_SQUARE_PRECISION(v_pt1.vector_buffer[i]);
	//printf("%d %d: norm1 norm2 %g %g\n",j, k1,norm1,norm2);
      }

      norm1 = 1/sqrt(norm1); norm2 = 1/sqrt(norm2);
      for ( i=0; i<aggregate_size; ) {
        for ( k=0; k<offset; k++, i++ )
          v_pt1.vector_buffer[i] =  norm1 * creal_PRECISION(v_pt1.vector_buffer[i]) + I*norm1* cimag_PRECISION(v_pt1.vector_buffer[i]);
        for ( k=0; k<offset; k++, i++ )
          v_pt1.vector_buffer[i] =  norm2 * creal_PRECISION(v_pt1.vector_buffer[i]) + I*norm2* cimag_PRECISION(v_pt1.vector_buffer[i]);
      }
    }
  }

  vector_PRECISION_change_layout( phi, phi, _NVEC_INNER, no_threading );
  SYNC_HYPERTHREADS(threading)
  SYNC_CORES(threading)
  PROF_PRECISION_STOP( _GRAM_SCHMIDT_ON_AGGREGATES, 1, threading );
}

/******* ????????? ***************************/

/*
void spinwise_PRECISION_skalarmultiply( vector_PRECISION *eta1, vector_PRECISION *eta2, vector_PRECISION *phi, complex_PRECISION alpha,
                                        int start, int end, level_struct *l ) {
  
  PROF_PRECISION_START( _LA6 );  
  for ( int i=start; i<end; ) {
    FOR6( eta1->vector_buffer[i] = alpha*phi->vector_buffer[i]; eta2->vector_buffer[i] = _COMPLEX_PRECISION_ZERO; i++; )
    FOR6( eta2->vector_buffer[i] = alpha*phi->vector_buffer[i]; eta1->vector_buffer[i] = _COMPLEX_PRECISION_ZERO; i++; )
  }
  PROF_PRECISION_STOP( _LA6, 1 );
}
*/

/*
void gram_schmidt_PRECISION( vector_PRECISION *V, complex_PRECISION *buffer, const int begin, const int n, level_struct *l, struct Thread *threading ) {
  
  // NOTE: only thread safe, if "buffer" is the same buffer for all threads belonging to a common MPI process
  START_MASTER(threading)
  PROF_PRECISION_START( _LA );
  END_MASTER(threading)
  SYNC_CORES(threading)
  
  PRECISION beta;
  int i, j, start, end;
  
  compute_core_start_end_custom( 0, l->inner_vector_size, &start, &end, l, threading, l->num_lattice_site_var );
  
  for ( i=begin; i<n; i++ ) {
    
    complex_PRECISION tmp[i];
    process_multi_inner_product_PRECISION( i, tmp, V, &V[i], 0, l->inner_vector_size, l, threading );
    SYNC_CORES(threading)
    START_MASTER(threading)
    for ( j=0; j<i; j++ ) {
      buffer[j] = tmp[j];
    }
    END_MASTER(threading)
    SYNC_MASTER_TO_ALL(threading)
    
    if ( i>0 ) {
      START_MASTER(threading)
      PROF_PRECISION_START( _ALLR );
      MPI_Allreduce( buffer, buffer+n, i, MPI_COMPLEX_PRECISION, MPI_SUM, (l->depth==0)?g.comm_cart:l->gs_PRECISION.level_comm );
      PROF_PRECISION_STOP( _ALLR, 1 );
      END_MASTER(threading)
      SYNC_MASTER_TO_ALL(threading)
    }
    
    for( j=0; j<i; j++ ) {
      vector_PRECISION_saxpy( &V[i], &V[i], &V[j], -(buffer+n)[j], start, end, l );
      SYNC_CORES(threading)
    }
    
    SYNC_CORES(threading)
      
    beta = global_norm_PRECISION( &V[i], 0, l->inner_vector_size, l, threading );
    SYNC_MASTER_TO_ALL(threading)
    vector_PRECISION_real_scale( &V[i], &V[i], creal(1.0/beta), start, end, l );
    SYNC_CORES(threading)
  }
  
  START_MASTER(threading)
  PROF_PRECISION_STOP( _LA, 1 );
  END_MASTER(threading)
  SYNC_CORES(threading)
}
*/

/****************** NOT USED ************************************/
#if 0
// used only in gram_schmidt_PRECISION_new 
PRECISION global_norm_PRECISION_new_2( vector_PRECISION *x, int ind1, int start, int end, level_struct *l, struct Thread *threading ) {
  
  PROF_PRECISION_START( _GIP, threading );
  
  PRECISION local_alpha = 0, global_alpha = 0;

  int thread_start;
  int thread_end;
  compute_core_start_end(start, end, &thread_start, &thread_end, l, threading);
  
  SYNC_CORES(threading)
  
  VECTOR_FOR( int i=thread_start, i<thread_end, local_alpha += NORM_SQUARE_PRECISION(x->vector_buffer[i+ind1*end]), i++, l );

  // sum over cores
  START_NO_HYPERTHREADS(threading)
  ((PRECISION *)threading->workspace)[threading->core] = local_alpha;
  END_NO_HYPERTHREADS(threading)
  // master sums up all results
  SYNC_CORES(threading)
  START_MASTER(threading)
  for(int i=1; i<threading->n_core; i++)
    ((PRECISION *)threading->workspace)[0] += ((PRECISION *)threading->workspace)[i];
  local_alpha = ((PRECISION *)threading->workspace)[0];
  END_MASTER(threading)

  if ( g.num_processes > 1 ) {
    START_MASTER(threading)
    PROF_PRECISION_START( _ALLR );
    MPI_Allreduce( &local_alpha, &global_alpha, 1, MPI_PRECISION, MPI_SUM, (l->depth==0)?g.comm_cart:l->gs_PRECISION.level_comm );
    PROF_PRECISION_STOP( _ALLR, 1 );
    ((PRECISION *)threading->workspace)[0] = global_alpha;
    END_MASTER(threading)
    // all threads need the result of the norm
    SYNC_MASTER_TO_ALL(threading)
    global_alpha = ((PRECISION *)threading->workspace)[0];
    PROF_PRECISION_STOP( _GIP, (double)(end-start)/(double)l->inner_vector_size, threading );
    return (PRECISION)sqrt((double)global_alpha);
  } else {
    // all threads need the result of the norm
    SYNC_MASTER_TO_ALL(threading)
    local_alpha = ((PRECISION *)threading->workspace)[0];
    PROF_PRECISION_STOP( _GIP, (double)(end-start)/(double)l->inner_vector_size, threading );
    return (PRECISION)sqrt((double)local_alpha);
  }
}

// used only in gram_schmidt_PRECISION_new_w
void process_multi_inner_product_PRECISION_new_2( int count, complex_PRECISION *results, vector_PRECISION *phi, vector_PRECISION *psi,
    int start, int end, level_struct *l, struct Thread *threading ) {

  PROF_PRECISION_START( _PIP, threading );
  int i;
  for(int c=0; c<count; c++)
    results[c] = 0.0;

  int thread_start;
  int thread_end;

  SYNC_CORES(threading)
  if ( l->depth == 0 ) {
    compute_core_start_end_custom(start, end, &thread_start, &thread_end, l, threading, 12);
    for(int c=0; c<count; c++)
      for ( i=thread_start; i<thread_end; )
        FOR12( results[c] += conj_PRECISION(phi->vector_buffer[i+c*end])*psi->vector_buffer[i+count*end]; i++; )
  } else {
#ifdef _M10TV
    compute_core_start_end_custom(start, end, &thread_start, &thread_end, l, threading, 20);
    for(int c=0; c<count; c++)
      for ( i=thread_start; i<thread_end; )
        FOR20( results[c] += conj_PRECISION(phi->vector_buffer[i+c*end])*psi->vector_buffer[i+count*end]; i++; )
#else
    compute_core_start_end_custom(start, end, &thread_start, &thread_end, l, threading, 2);
    for(int c=0; c<count; c++)
      for ( i=thread_start; i<thread_end; )
        FOR2( results[c] += conj_PRECISION(phi->vector_buffer[i+c*end])*psi->vector_buffer[i+count*end]; i++; )
#endif
  }

  START_NO_HYPERTHREADS(threading)
  ((complex_PRECISION **)threading->workspace)[threading->core] = results;
  END_NO_HYPERTHREADS(threading)
  // master sums up all results
  SYNC_CORES(threading)
  START_MASTER(threading)
  for(int c=0; c<count; c++)
    for(int i=1; i<threading->n_core; i++)
      ((complex_PRECISION **)threading->workspace)[0][c] += ((complex_PRECISION **)threading->workspace)[i][c];
  END_MASTER(threading)
  // all threads need the result of the norm
  SYNC_MASTER_TO_ALL(threading)
  for(int c=0; c<count; c++)
    results[c] = ((complex_PRECISION **)threading->workspace)[0][c];

  PROF_PRECISION_STOP( _PIP, (double)(end-start)/(double)l->inner_vector_size, threading );
}

// used only in gram_schmidt_PRECISION_new
void vector_PRECISION_saxpy_new_2( vector_PRECISION *z, int ind1, vector_PRECISION *x, int ind2, vector_PRECISION *y, int ind3, complex_PRECISION alpha, int start, int end, int size, level_struct *l ) {

  int thread = omp_get_thread_num();
  if (thread == 0 && start != end )
  PROF_PRECISION_START( _LA8 );
  
  VECTOR_FOR( int i=start, i<end, z->vector_buffer[i+ind1*size] = x->vector_buffer[i+ind2*size] + alpha*y->vector_buffer[i+ind3*size], i++, l );
  
  if( thread == 0 && start != end )
  PROF_PRECISION_STOP( _LA8, (double)(end-start)/(double)l->inner_vector_size );
}


void gram_schmidt_PRECISION_new_2( vector_PRECISION *V, complex_PRECISION *buffer, const int begin, const int n, level_struct *l, struct Thread *threading ) {
  /*****************************
   * a set of vectors in V are orthonormalized
   * begin: where to start
   * n: #vectors from begin to be orthonormalized
   *****************************/

  // NOTE: only thread safe, if "buffer" is the same buffer for all threads belonging to a common MPI process
  START_MASTER(threading)
  PROF_PRECISION_START( _LA );
  END_MASTER(threading)
  SYNC_CORES(threading)
  
  PRECISION beta;
  int i, j, start, end;
  
  compute_core_start_end_custom( 0, l->inner_vector_size, &start, &end, l, threading, l->num_lattice_site_var );
  vector_PRECISION_change_layout( V, V, _NV_LV_SV, no_threading );
  for ( i=begin; i<n; i++ ) {
    
    complex_PRECISION tmp[i];
    process_multi_inner_product_PRECISION_new_2( i, tmp, V, V, 0, l->inner_vector_size, l, threading );
    SYNC_CORES(threading)
    START_MASTER(threading)
    for ( j=0; j<i; j++ ) {
      buffer[j] = tmp[j];
    }
    END_MASTER(threading)
    SYNC_MASTER_TO_ALL(threading)
    
    if ( i>0 ) {
      START_MASTER(threading)
      PROF_PRECISION_START( _ALLR );
      MPI_Allreduce( buffer, buffer+n, i, MPI_COMPLEX_PRECISION, MPI_SUM, (l->depth==0)?g.comm_cart:l->gs_PRECISION.level_comm );
      PROF_PRECISION_STOP( _ALLR, 1 );
      END_MASTER(threading)
      SYNC_MASTER_TO_ALL(threading)
    }
    
    for( j=0; j<i; j++ ) {
      vector_PRECISION_saxpy_new_2( V, i, V, i, V, j, -(buffer+n)[j], start, end, l->inner_vector_size, l );
      SYNC_CORES(threading)
    }
    
    SYNC_CORES(threading)
      
    beta = global_norm_PRECISION_new_2( V, i, 0, l->inner_vector_size, l, threading );
    SYNC_MASTER_TO_ALL(threading)
    vector_PRECISION_real_scale_new_2( V, i, V, i, creal(1.0/beta), start, end, l->inner_vector_size, l );
    SYNC_CORES(threading)
  }
  vector_PRECISION_change_layout( V, V, _LV_SV_NV, no_threading );
  START_MASTER(threading)
  PROF_PRECISION_STOP( _LA, 1 );
  END_MASTER(threading)
  SYNC_CORES(threading)
}
#endif
/****************** NOT USED ************************************
// used only in setup_gram_schmidt_PRECISION
void setup_gram_schmidt_PRECISION_compute_dots(
    complex_PRECISION *thread_buffer, vector_PRECISION *V, int count, int offset,
    int start, int end, level_struct *l, struct Thread *threading) {

  int thread_start;
  int thread_end;
  int cache_block_size = 12*64;
  complex_PRECISION tmp[cache_block_size];
  vector_PRECISION tmp_vect;
  tmp_vect.vector_buffer = tmp;

  for(int i=0; i<2*offset; i++)
    thread_buffer[i] = 0.0;

  SYNC_CORES(threading)
  compute_core_start_end_custom(start, end, &thread_start, &thread_end, l, threading, cache_block_size);
  
  for ( int i=thread_start; i<thread_end; i+=cache_block_size) {
    coarse_gamma5_PRECISION( &tmp_vect, &V[count]+i, 0, cache_block_size, l );
    for ( int j=0; j<count; j++ ) {
      for ( int k=0; k<cache_block_size; k++) {
        thread_buffer[j]   += conj_PRECISION(V[j].vector_buffer[i+k])*V[count].vector_buffer[i+k];
        thread_buffer[j+offset] += conj_PRECISION(V[j].vector_buffer[i+k])*tmp[k];
      }
    }
  }

  START_NO_HYPERTHREADS(threading)
  ((complex_PRECISION **)threading->workspace)[threading->core] = thread_buffer;
  END_NO_HYPERTHREADS(threading)
  // master sums up all results
  SYNC_CORES(threading)
  START_MASTER(threading)
  for(int i=1; i<threading->n_core; i++) {
    for(int j=0; j<count; j++) {
      ((complex_PRECISION **)threading->workspace)[0][j]        += ((complex_PRECISION **)threading->workspace)[i][j];
      ((complex_PRECISION **)threading->workspace)[0][j+offset] += ((complex_PRECISION **)threading->workspace)[i][j+offset];
    }
  }
  END_MASTER(threading)
  // only master needs the result in this case (it will be distributed later)
}

// used only in setup_gram_schmidt_PRECISION  
void setup_gram_schmidt_PRECISION_axpys(
    complex_PRECISION *thread_buffer, vector_PRECISION *V, int count, int offset,
    int start, int end, level_struct *l, struct Thread *threading) {
  
  int thread_start;
  int thread_end;
  int cache_block_size = 12*64;
  complex_PRECISION tmp[cache_block_size];
  vector_PRECISION tmp_vect;
  tmp_vect.vector_buffer = tmp;

  compute_core_start_end_custom(start, end, &thread_start, &thread_end, l, threading, cache_block_size);

  for ( int i=thread_start; i<thread_end; i+=cache_block_size) {
    for ( int j=0; j<count; j++ ) {
      coarse_gamma5_PRECISION( &tmp_vect, &V[j]+i, 0, cache_block_size, l );
      for ( int k=0; k<cache_block_size; k++) {
        V[count].vector_buffer[i+k] -= thread_buffer[2*offset+j]*V[j].vector_buffer[i+k];
        V[count].vector_buffer[i+k] -= thread_buffer[3*offset+j]*tmp[k];
      }
    }
  }
}

// used only in setup_gram_schmidt_PRECISION  
complex_PRECISION process_inner_product_PRECISION( vector_PRECISION *phi, vector_PRECISION *psi, int start, int end, level_struct *l, struct Thread *threading ) {
  
  PROF_PRECISION_START( _PIP, threading );
  int i;
  complex_PRECISION local_alpha = 0;
  
  SYNC_CORES(threading)
  
  THREADED_VECTOR_FOR( i, start, end, local_alpha += conj_PRECISION(phi->vector_buffer[i])*psi->vector_buffer[i], i++, l, threading );

  START_NO_HYPERTHREADS(threading)
  ((complex_PRECISION *)threading->workspace)[threading->core] = local_alpha;
  END_NO_HYPERTHREADS(threading)
  // master sums up all results
  SYNC_CORES(threading)
  START_MASTER(threading)
  for(int i=1; i<threading->n_core; i++)
    ((complex_PRECISION *)threading->workspace)[0] += ((complex_PRECISION *)threading->workspace)[i];
  END_MASTER(threading)
  // all threads need the result of the norm
  SYNC_MASTER_TO_ALL(threading)
  local_alpha = ((complex_PRECISION *)threading->workspace)[0];

  PROF_PRECISION_STOP( _PIP, (double)(end-start)/(double)l->inner_vector_size, threading );

  return local_alpha;
}


void setup_gram_schmidt_PRECISION( vector_PRECISION *V, vector_PRECISION *g5v,
                                   complex_PRECISION *buffer, const int n, level_struct *l,
                                   struct Thread *threading ) {
  
  PROF_PRECISION_START( _GRAM_SCHMIDT, threading );
  PRECISION beta;
  int i, j;
  int start = 0;
  int end = l->inner_vector_size;
  int thread_start = threading->start_index[l->depth];
  int thread_end   = threading->end_index[l->depth];

  complex_PRECISION thread_buffer[4*n];
  
  for ( i=0; i<4*n; i++ )
    thread_buffer[i] = 0;
  
  for ( i=0; i<n; i++ ) {
    
    if ( l->depth > 0 ) {
      coarse_gamma5_PRECISION( g5v, &V[i], thread_start, thread_end, l );
      for ( j=0; j<i; j++ ) {
        thread_buffer[j] = process_inner_product_PRECISION( &V[j], &V[i], start, end, l, threading );
        thread_buffer[j+n] = process_inner_product_PRECISION( &V[j], g5v, start, end, l, threading );
      }
    }
    else
      setup_gram_schmidt_PRECISION_compute_dots( thread_buffer, V, i, n, start, end, l, threading);
    
    
    START_LOCKED_MASTER(threading)
    if ( i>0 ) {
      PROF_PRECISION_START( _ALLR );
      MPI_Allreduce( thread_buffer, thread_buffer+2*n, 2*n, MPI_COMPLEX_PRECISION, MPI_SUM, (l->depth==0)?g.comm_cart:l->gs_PRECISION.level_comm );
      PROF_PRECISION_STOP( _ALLR, 1 );
    }
    for ( j=2*n; j<4*n; j++ )
      ((complex_PRECISION *)(threading->workspace))[j] = thread_buffer[j];
    END_LOCKED_MASTER(threading)
    for ( j=2*n; j<4*n; j++ )
      thread_buffer[j] = ((complex_PRECISION *)(threading->workspace))[j];

    
    if ( l->depth > 0 ) {
      for( j=0; j<i; j++ ) {
        vector_PRECISION_saxpy( &V[i], &V[i], &V[j], -(thread_buffer+2*n)[j], thread_start, thread_end, l );
        coarse_gamma5_PRECISION( g5v, &V[j], thread_start, thread_end, l );
        vector_PRECISION_saxpy( &V[i], &V[i], g5v, -(thread_buffer+3*n)[j], thread_start, thread_end, l );
      }
    } else {
      setup_gram_schmidt_PRECISION_axpys( thread_buffer, V, i, n, start, end, l, threading);
    }
    
    beta = global_norm_PRECISION( &V[i], start, end, l, threading );
    vector_PRECISION_real_scale( &V[i], &V[i], 1.0/beta, thread_start, thread_end, l );
  }
  PROF_PRECISION_STOP( _GRAM_SCHMIDT, (double)(end-start)/(double)l->inner_vector_size, threading );
}
*/
