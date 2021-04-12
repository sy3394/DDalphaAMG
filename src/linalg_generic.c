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

void global_norm_PRECISION( PRECISION *res, vector_PRECISION *x, int start, int end, level_struct *l, struct Thread *threading ) {
  /*
   * Note: end-start needs to match with the local size of x; otherwise it computes partial norm
   */
  
  int i, j, jj, nvec = x->num_vect_now, nvecmf = x->num_vect_now;
#ifdef HAVE_TM1p1
#ifdef DEBUG
  if ( g.n_flavours == 2 && nvec != 2*num_loop )
    error0("global_norm_PRECISION: doublet error\n");
#endif
  if ( g.n_flavours == 2 && ( g.epsbar != 0 || g.epsbar_ig5_odd_shift != 0 || g.epsbar_ig5_odd_shift != 0 ) )
    nvec /= g.n_flavours;
#endif
  int thread = omp_get_thread_num();
  PRECISION global_alpha[nvec];

  if ( thread == 0 && start != end )
    PROF_PRECISION_START( _GIP, threading );

  //thread == core here: distribute entries from start to end among cores
  int core_start, core_end;
  compute_core_start_end_custom(start, end, &core_start, &core_end, l, threading, l->num_lattice_site_var );//num_loop);

  VECTOR_LOOP(j, nvec, jj, res[j+jj]=0;)

    //SYNC_CORES(threading)//????necessary? in the original   
  for( i=core_start; i<core_end; i++) {
    VECTOR_LOOP(j, nvecmf, jj, res[(j+jj)%nvec] += NORM_SQUARE_PRECISION(x->vector_buffer[i*x->num_vect+j+jj]);)
    //for ( int f = 0; f<g.n_flavours-1; f++)VECTOR_LOOP(j, nvec, jj, res[j+jj] += NORM_SQUARE_PRECISION(x->vector_buffer[i*x->num_vect+j+jj]+f*nvec);)//expensive?
  }
  
  /////////////// communication -------------------------------------
  // sum over cores: be careful about overflow of workspace!!!!!!!! also potential problem using res directly
  START_NO_HYPERTHREADS(threading)
  VECTOR_LOOP(j, nvec, jj, ((PRECISION *)threading->workspace)[threading->core*nvec+j+jj] = res[j+jj];) 
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

void global_inner_product_PRECISION( complex_PRECISION *results, vector_PRECISION *phi, vector_PRECISION *psi, int start, int end, level_struct *l, struct Thread *threading ) {

#ifdef DEBUG
  if ( phi->num_vect_now != psi->num_vect_now )
    error0("global_inner_product_PRECISION: phi->num_vect_now != psi->num_vect_now \n");
#endif
  
  PROF_PRECISION_START( _GIP, threading );
  int i, j, jj, nvec = phi->num_vect_now, nvecmf = phi->num_vect_now;
#ifdef HAVE_TM1p1
#ifdef DEBUG
  if ( g.n_flavours == 2 && nvec != 2*num_loop )
    error0("global_inner_product_PRECISION: doublet error\n");
#endif
  if ( g.n_flavours == 2 && ( g.epsbar != 0 || g.epsbar_ig5_odd_shift != 0 || g.epsbar_ig5_odd_shift != 0 ) )
    nvec /= g.n_flavours;
#endif
  complex_PRECISION local_alpha[nvec], global_alpha[nvec];//local_alpha not necessary!!!! kept for readability
 
  VECTOR_LOOP( j, nvec, jj, local_alpha[j+jj]=0; global_alpha[j+jj]=0; )

  int thread_start;
  int thread_end;
  compute_core_start_end_custom(start, end, &thread_start, &thread_end, l, threading, l->num_lattice_site_var );//num_loop);
  
  //SYNC_CORES(threading)
  for( i=thread_start; i<thread_end; i++){
    VECTOR_LOOP( j, nvecmf, jj, local_alpha[(j+jj)%nvec] += conj_PRECISION(phi->vector_buffer[i*phi->num_vect+j+jj])*psi->vector_buffer[i*psi->num_vect+j+jj]; )
  }
  
  /////////////// communication -----------------------------------
  // sum over cores
  //// copy the results into threading->workspace 
  START_NO_HYPERTHREADS(threading)
  VECTOR_LOOP( j, nvec, jj, ((complex_PRECISION *)threading->workspace)[threading->core*nvec+j+jj] = local_alpha[j+jj];)
  END_NO_HYPERTHREADS(threading)
  //// master sums up all results
  SYNC_CORES(threading)
  START_MASTER(threading)
  for( i=1; i<threading->n_core; i++)
    VECTOR_LOOP( j, nvec, jj, ((complex_PRECISION *)threading->workspace)[0*nvec+j+jj] += ((complex_PRECISION *)threading->workspace)[i*nvec+j+jj];)
  END_MASTER(threading)
  SYNC_MASTER_TO_ALL(threading)
  VECTOR_LOOP( j, nvec, jj, results[j+jj] = ((complex_PRECISION *)threading->workspace)[0*nvec+j+jj];)
  // sum over ranks
  if ( g.num_processes > 1 ) {
    START_MASTER(threading)
    PROF_PRECISION_START( _ALLR );
    MPI_Allreduce( results, global_alpha, nvec, MPI_COMPLEX_PRECISION, MPI_SUM, (l->depth==0)?g.comm_cart:l->gs_PRECISION.level_comm );
    PROF_PRECISION_STOP( _ALLR, 1 );
    VECTOR_LOOP( j, nvec, jj, ((complex_PRECISION *)threading->workspace)[0*nvec+j+jj] = global_alpha[j+jj];)
    END_MASTER(threading)
    // all threads need the result of the norm
    SYNC_MASTER_TO_ALL(threading)
    VECTOR_LOOP( j, nvec, jj, results[j+jj] = ((complex_PRECISION *)threading->workspace)[0*nvec+j+jj];)
  }

  PROF_PRECISION_STOP( _GIP, (double)(end-start)/(double)l->inner_vector_size, threading );  

}

// results <- (phi's, psi) processwise
// assume resutls have phi->num_vect_now fields
void process_multi_inner_product_PRECISION( int count, complex_PRECISION *results, vector_PRECISION *phi, vector_PRECISION *psi,
						int start, int end, level_struct *l, struct Thread *threading ) {
  /******************************************
   * Input:
   *  int count: #basis vectors to be multiplied
   *  vector_PRECISION *phi: array of basis vectors
   *  vector_PRECISION *psi: a vector multiplying each basis vector in the array
   * Output:
   *  complex_PRECISION *results: stores the inner products of each vector in the basis phi and the given vector psi on each process
   *****************************************/

  int core_start, core_end;
  compute_core_start_end_custom(start, end, &core_start, &core_end, l, threading, l->num_lattice_site_var);

  int thread = omp_get_thread_num();
  if(thread == 0 && core_start != core_end)
    PROF_PRECISION_START( _PIP, threading );

  //SYNC_CORES(threading)//????necessary? in the original  
  int i, c, j, jj, nvec = psi->num_vect_now, nvecmf = psi->num_vect_now;
#ifdef HAVE_TM1p1
#ifdef DEBUG
  if ( g.n_flavours == 2 && nvec != 2*num_loop )
    error0("process_multi_inner_product_PRECISION: doublet error\n");
#endif
  if ( g.n_flavours == 2 && ( g.epsbar != 0 || g.epsbar_ig5_odd_shift != 0 || g.epsbar_ig5_odd_shift != 0 ) )
    nvec /= g.n_flavours;
#endif
  VECTOR_LOOP(j, count*nvec, jj, results[j+jj] = 0.0;)
    //printf0("mul inn: %d %d %d\n",count,nvec,nvecmf);
  for( c = 0; c<count; c++ ){
#ifdef DEBUG
    if ( phi[c].num_vect_now != psi->num_vect_now )
      error0("process_multi_inner_product_PRECISION: phi->num_vect_now(%d) != psi->num_vect_now(%d) \n",phi[c].num_vect_now,psi->num_vect_now);
#endif
    for ( i=core_start; i<core_end; i++ )
      VECTOR_LOOP(j, nvecmf, jj, results[c*nvec+(j+jj)%nvec] += conj_PRECISION(phi[c].vector_buffer[i*phi[c].num_vect+j+jj])*psi->vector_buffer[i*psi->num_vect+j+jj];)//if(!g.in_setup)printf0("mul in: %g ",conj_PRECISION(phi[c].vector_buffer[i*phi[c].num_vect+j+jj])*psi->vector_buffer[i*psi->num_vect+j+jj]);)
  }

  // copy results into threading->workspace
  START_NO_HYPERTHREADS(threading)
    VECTOR_LOOP( j, count*nvec, jj,((complex_PRECISION *)threading->workspace)[threading->core*count*nvec+j+jj] = results[j+jj];)//printf0("mul inn: %g\n",creal_PRECISION(results[j+jj]));)
  END_NO_HYPERTHREADS(threading)
  // master sums up all results
  SYNC_CORES(threading)
  START_MASTER(threading)
  for( i=1; i<threading->n_core; i++)
    VECTOR_LOOP( j, count*nvec, jj,((complex_PRECISION *)threading->workspace)[j+jj] += ((complex_PRECISION *)threading->workspace)[i*count*nvec+j+jj];)
  END_MASTER(threading)
  SYNC_MASTER_TO_ALL(threading)
  // all threads need the result of the norm
  VECTOR_LOOP( j, count*nvec, jj, results[j+jj] = ((complex_PRECISION *)threading->workspace)[j+jj];)

  if(thread == 0 && core_start != core_end)
    PROF_PRECISION_STOP( _PIP, (double)(core_end-core_start)/(double)l->inner_vector_size, threading );
}

// res <- <phi,psi>/<phi,phi> with each vector sliced from start to end
// used only in local_minres_PRECISION, which assumes, for the moment, res to be of size, num_vect_now
void local_xy_over_xx_PRECISION( complex_PRECISION *res, vector_PRECISION *phi, vector_PRECISION *psi, int start, int end, level_struct *l  ) {

#ifdef DEBUG
  if ( phi->num_vect_now != psi->num_vect_now )
    error0("local_xy_over_xx_PRECISION: phi->num_vect_now != psi->num_vect_now \n");
#endif
  
  int nvec = phi->num_vect_now, i, j, jj, nvecmf = phi->num_vect_now;
#ifdef HAVE_TM1p1
#ifdef DEBUG
  if ( g.n_flavours == 2 && nvec != 2*num_loop )
    error0("process_multi_inner_product_PRECISION: doublet error\n");
#endif
  if ( g.n_flavours == 2 && ( g.epsbar != 0 || g.epsbar_ig5_odd_shift != 0 || g.epsbar_ig5_odd_shift != 0 ) )
    nvec /= g.n_flavours;
#endif
  complex_PRECISION numerator[nvec]; PRECISION denominator[nvec], total_den=0.0;
  VECTOR_LOOP(j, nvec, jj, numerator[j+jj]=0.0; denominator[j+jj]=0.0;)
  
  for( i=start; i<end; i++ )
    VECTOR_LOOP(j, nvecmf, jj,
		numerator[(j+jj)%nvec] += conj_PRECISION(phi->vector_buffer[i*phi->num_vect+j+jj])*psi->vector_buffer[i*psi->num_vect+j+jj]; 
		denominator[(j+jj)%nvec] += NORM_SQUARE_PRECISION(phi->vector_buffer[i*phi->num_vect+j+jj]);)

  VECTOR_LOOP(j, nvec, jj, total_den += denominator[j+jj];)
  if ( abs_PRECISION(total_den) < nvecmf*EPS_PRECISION ) {
    VECTOR_LOOP(j, nvecmf, jj, res[j+jj] = 0.0;)
    return;
  }
  
  VECTOR_LOOP(j, nvecmf, jj, res[j+jj] = numerator[(j+jj)%nvec]/denominator[(j+jj)%nvec];)
}

void gram_schmidt_PRECISION( vector_PRECISION *V, const int nvec, level_struct *l, struct Thread *threading ) {
  /*****************************
   * a set of vectors in V are orthonormalized
   * nvec: #first vectors in the set to be orthonormalized
   * NOTE: only thread safe, if "buffer" is the same buffer for all threads belonging to a common MPI process
   * NOTE: used only in setting up interpolation op
   *****************************/

  START_MASTER(threading)
  PROF_PRECISION_START( _LA );
  END_MASTER(threading)
  SYNC_CORES(threading)

  int i, j, k, start, end;  
  PRECISION norm;
  complex_PRECISION local_alpha, global_alpha; 
  
  compute_core_start_end_custom(0, l->inner_vector_size, &start, &end, l, threading, l->num_lattice_site_var );

#ifdef DEBUG
  if ( V->layout != _NVEC_INNER )
    error0("gram_schmidt_PRECISION: assumptions are not met\n");
#endif
  
  for ( i=0; i<nvec; i++ ) {
    for ( j=0; j<i; j++ ) {
      //take inner product <i,j>
      //SYNC_CORES(threading)//changed!!!!!
      local_alpha = 0;
      for( k=start; k<end; k++) // for each core: do we need to divide into two parts like in aggregate version????????
	local_alpha += conj_PRECISION(V->vector_buffer[k*V->num_vect+j])*V->vector_buffer[k*V->num_vect+i];//is this inefficient?????
      // communication  -----------------------------------
      // sum over cores
      //// copy results into threading->workspace
      START_NO_HYPERTHREADS(threading)
      ((complex_PRECISION *)threading->workspace)[threading->core] = local_alpha;
      END_NO_HYPERTHREADS(threading)
      //// master sums up all results
      SYNC_CORES(threading)
      START_MASTER(threading)
      for( k=1; k<threading->n_core; k++)
	((complex_PRECISION *)threading->workspace)[0] += ((complex_PRECISION *)threading->workspace)[k];
      END_MASTER(threading)
      SYNC_MASTER_TO_ALL(threading)
      local_alpha = ((complex_PRECISION *)threading->workspace)[0];
      // sum over ranks
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
      
      // v_i -= alpha*v_j//may use function????
      for( k=start; k<end; k++) // for each core: do we need to divide into two parts like in aggregate version????????
        V->vector_buffer[k*V->num_vect+i] -= global_alpha*V->vector_buffer[k*V->num_vect+j];
    }
    // take squared norm
    local_alpha = 0;
    for( k=start; k<end; k++) // for each core: do we need to divide into two parts like in aggregate version????????
      local_alpha += NORM_SQUARE_PRECISION(V->vector_buffer[k*V->num_vect+i]);
    
    // communication  -----------------------------------
    // sum over cores
    //// copy results into threading->workspace
    START_NO_HYPERTHREADS(threading)
    ((complex_PRECISION *)threading->workspace)[threading->core] = local_alpha;
    END_NO_HYPERTHREADS(threading)
    //// master sums up all results
    SYNC_CORES(threading)
    START_MASTER(threading)
    for( k=1; k<threading->n_core; k++)
      ((complex_PRECISION *)threading->workspace)[0] += ((complex_PRECISION *)threading->workspace)[k];
    END_MASTER(threading)
    SYNC_MASTER_TO_ALL(threading)
    local_alpha = ((complex_PRECISION *)threading->workspace)[0];
    // sum over ranks
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
    for( k=start; k<end; k++) // for each core: do we need to divide into two parts like in aggregate version????????
      V->vector_buffer[k*V->num_vect+i] /= sqrt(global_alpha);
  }

  START_MASTER(threading)
  PROF_PRECISION_STOP( _LA, 1 );
  END_MASTER(threading)
  SYNC_CORES(threading)
}

void gram_schmidt_on_aggregates_PRECISION( vector_PRECISION *phi, const int num_vect, level_struct *l, struct Thread *threading ) {
  /************************************************
   * vector_PRECISION *phi: a set of vectors to be orthogonalized 
   * num_vect: should be equal to phi->num_vect_now
   * note: aggregate lives in the process so that no communication is necessary
   * note: used only in setting up interpolation op  
   ***********************************************/

  PROF_PRECISION_START( _GRAM_SCHMIDT_ON_AGGREGATES, threading );
  SYNC_CORES(threading)
  SYNC_HYPERTHREADS(threading)
  int i, j, k, k1, k2, nvec_phi = phi->num_vect;
  int num_aggregates = l->s_PRECISION.num_aggregates;
  int aggregate_size = l->inner_vector_size / num_aggregates;
  int offset         = l->num_lattice_site_var/2;
      
  complex_PRECISION alpha1, alpha2;
  PRECISION norm1, norm2;

#ifdef DEBUG
  if ( num_vect != phi->num_vect_now )
    error0("gram_schmidt_on_aggregates_PRECISION: assumptions are not met\n");
#endif
  
  buffer_PRECISION phi_pt = phi->vector_buffer;
  for ( j=threading->n_thread*threading->core+threading->thread; j<num_aggregates; j+=threading->n_thread*threading->n_core ) {
    for ( k1=0; k1<num_vect; k1++ ) {
      for ( k2=0; k2<k1; k2++ ) {
	alpha1 = 0; alpha2 = 0;
	  // V[k1] -= <V[k2],V[k1]> V[k2] | 2*j-th and 2*j+1-st aggregate
        for ( i=0; i<aggregate_size; ) {
          for ( k=0; k<offset; k++, i++ )
	    alpha1 += conj_PRECISION(phi->vector_buffer[(j*aggregate_size+i)*nvec_phi+k2]) * phi->vector_buffer[(j*aggregate_size+i)*nvec_phi+k1];
          for ( k=0; k<offset; k++, i++ )
            alpha2 += conj_PRECISION(phi->vector_buffer[(j*aggregate_size+i)*nvec_phi+k2]) * phi->vector_buffer[(j*aggregate_size+i)*nvec_phi+k1];
        }
        for ( i=0; i<aggregate_size; ) {
          for ( k=0; k<offset; k++, i++ )
	    phi->vector_buffer[(j*aggregate_size+i)*nvec_phi+k1] -= alpha1 * phi->vector_buffer[(j*aggregate_size+i)*nvec_phi+k2];
          for ( k=0; k<offset; k++, i++ )
	    phi->vector_buffer[(j*aggregate_size+i)*nvec_phi+k1] -= alpha2 * phi->vector_buffer[(j*aggregate_size+i)*nvec_phi+k2];
        }
      }

      norm1 = 0; norm2 = 0;
      // V[k1] = V[k1]/norm(V[k1]) | 2*j-th and 2*j+1-st aggregate    
      for ( i=0; i<aggregate_size; ) {
        for ( k=0; k<offset; k++, i++ )
          norm1 += NORM_SQUARE_PRECISION(phi->vector_buffer[(j*aggregate_size+i)*nvec_phi+k1]);
        for ( k=0; k<offset; k++, i++ )
          norm2 += NORM_SQUARE_PRECISION(phi->vector_buffer[(j*aggregate_size+i)*nvec_phi+k1]);
      }

      norm1 = 1/sqrt(norm1); norm2 = 1/sqrt(norm2);
      for ( i=0; i<aggregate_size; ) {
        for ( k=0; k<offset; k++, i++ )
          phi_pt[(j*aggregate_size+i)*nvec_phi+k1] =  norm1 * creal_PRECISION(phi_pt[(j*aggregate_size+i)*nvec_phi+k1]) + I*norm1* cimag_PRECISION(phi_pt[(j*aggregate_size+i)*nvec_phi+k1]);
        for ( k=0; k<offset; k++, i++ )
          phi_pt[(j*aggregate_size+i)*nvec_phi+k1] =  norm2 * creal_PRECISION(phi_pt[(j*aggregate_size+i)*nvec_phi+k1]) + I*norm2* cimag_PRECISION(phi_pt[(j*aggregate_size+i)*nvec_phi+k1]);
      }
    }
  }

  SYNC_HYPERTHREADS(threading)
  SYNC_CORES(threading)
  PROF_PRECISION_STOP( _GRAM_SCHMIDT_ON_AGGREGATES, 1, threading );
}
