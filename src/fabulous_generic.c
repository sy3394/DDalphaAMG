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
 * checked:11/30/2019
 * changed from sbacchio
 */

#include "main.h"

static void vector_PRECISION_copy_fab( vector_PRECISION *z, vector_PRECISION *x, int vec_ind, int length, int dir, int start, int end, level_struct *l );

void fabulous_PRECISION_init( fabulous_PRECISION_struct *fab ) {
  //when dynamically alloc fab, then do here and free in the below
  vector_PRECISION *X = &(fab->X) , *B = &(fab->B), *B0 = &(fab->B0), *X0 = &(fab->X0), *C0 = &(fab->C0);
  
  vector_PRECISION_init(X);
  vector_PRECISION_init(B);
  
  vector_PRECISION_init(X0);
  vector_PRECISION_init(B0);
  vector_PRECISION_init(C0);

  fab->eigvals = NULL;
  fab->k = 0;
  fab->handle = NULL;
  
}

void setup_fabulous_PRECISION( gmres_PRECISION_struct *p, int v_type, level_struct *l, struct Thread *threading ) {

  //----- Set fabulous structure
  fabulous_PRECISION_struct *fab = &(p->fab);
  //some fields could be eliminated!!!!
  fab->nrhs =p->num_vect; 
  int dim = (g.odd_even&&l->depth!=0)?(l->num_inner_lattice_sites/2)*l->num_lattice_site_var:p->v_end, nrhs = fab->nrhs;//!!!!!!
  fab->dim = dim; fab->ldb = dim; fab->ldx = dim;
  fab->l = l;
  fab->threading = threading;//needed?????
  //    printf0("fab PRECISION set: %d %d %d, %g %d %d\n",nrhs,dim,l->depth, p->tol, FABULOUS_COMPLEX_PRECISION, FABULOUS_COMPLEX_double);
  vector_PRECISION *X = &(fab->X) , *B = &(fab->B), *B0 = &(fab->B0), *X0 = &(fab->X0), *C0 = &(fab->C0);
  // The following fields are fed to fabulous solver; perhaps will be unnecessary if fabulous can handle _ORINARY properly!!!!
  vector_PRECISION_alloc( X, (g.odd_even&&l->depth!=0)?_EVEN_INNER:_INNER, nrhs, l, threading ); X->num_vect_now = nrhs;
  vector_PRECISION_alloc( B, (g.odd_even&&l->depth!=0)?_EVEN_INNER:_INNER, nrhs, l, threading ); B->num_vect_now = nrhs;
  // The following fields are used in mvp_PRECISION as temporay storage
  vector_PRECISION_alloc( X0, v_type, nrhs, l, threading ); X0->num_vect_now = nrhs; X0->layout = _NVEC_OUTER;
  vector_PRECISION_alloc( B0, v_type, nrhs, l, threading ); B0->num_vect_now = nrhs; B0->layout = _NVEC_OUTER;
  vector_PRECISION_alloc( C0, v_type, nrhs, l, threading ); C0->num_vect_now = nrhs; C0->layout = _NVEC_OUTER;
  if ( X->size != dim )
    error0("set_fabulous_struct_PRECISION: assumptions are not met\n");

  fab->k = g.p.fab.k;
  if ( fab->k > 0 ) {
    MALLOC( fab->eigvals, complex_PRECISION, fab->k );
    for ( int i=0; i<fab->k; i++ ) ((complex_PRECISION *) fab->eigvals)[i] = 0;
  }

  //----- Set fabulous handle
  fab->handle = fabulous_create(FABULOUS_COMPLEX_PRECISION, fab->dim, (void *) p);
  if ( fab->handle == NULL )
    error0("fabulous handle has not been created!\n");

  // Set callback:
  fabulous_set_mvp(&mvp_PRECISION, fab->handle);
  fabulous_set_dot_product(&dot_product_PRECISION, fab->handle);
  if ( g.use_fab_as_outer && l->level > 0 )
    fabulous_set_rightprecond(&fabulous_rightprecond_PRECISION, fab->handle);
  
  // Setup parameters:
  fabulous_set_ortho_process(g.f_orthoscheme, g.f_orthotype, g.ortho_iter, fab->handle);
  PRECISION tolerance[1] = { p->tol };
  fabulous_set_parameters( fab->nrhs*p->restart_length*p->num_restart, p->restart_length, tolerance, 1, fab->handle );//g.max_mvp
  fabulous_set_advanced_parameters( g.max_kept_direction, g.real_residual, g.logger_user_data_size, g.quiet, fab->handle );
  //printf0("fab param: %d %d %d %g %d %d %d %d %d\n",g.f_orthoscheme, g.f_orthotype, g.ortho_iter,tolerance[0],g.max_mvp, p->restart_length,g.max_kept_direction, g.real_residual, g.logger_user_data_size);

}

void fabulous_PRECISION_free( fabulous_PRECISION_struct *fab, level_struct *l, struct Thread *threading ) {

  vector_PRECISION_free(&(fab->X), l, threading);
  vector_PRECISION_free(&(fab->B), l, threading);
  vector_PRECISION_free(&(fab->X0), l, threading);
  vector_PRECISION_free(&(fab->B0), l, threading);
  vector_PRECISION_free(&(fab->C0), l, threading);
  if ( fab->k > 0 )
     PUBLIC_FREE( fab->eigvals, PRECISION, fab->k );
  
  fabulous_destroy(fab->handle);
  fab->handle = NULL;
}

// B <- beta*B+ alpha*A*X: Matrix x BlockOfVector product CALLBACK: vectors need to be in _NVEC_INNER
int64_t mvp_PRECISION(  void *user_env, int N,
			const void *p_alpha, const void *XX, int ldx,
			const void *p_beta, void *BB, int ldb) {

  gmres_PRECISION_struct *p = user_env;
  fabulous_PRECISION_struct *fab = &(p->fab);
  level_struct *l = fab->l;
  PROF_PRECISION_START( _FMVP );
  //printf0("mvp %d: %d vs %d; %d %d %d\n", g.my_rank, N, fab->B0.num_vect, ldx, ldb, fab->B0.size );fflush(stdout);
  //  START_UNTHREADED_FUNCTION(fab.threading)
  int n = fab->X0.num_vect;//fab->nrhs;!!!!!!
  complex_PRECISION alphas[n], betas[n];
  /*
  const complex_PRECISION alpha = *((const complex_PRECISION*) p_alpha), beta = *((const complex_PRECISION*) p_beta);
  for ( int i=0; i<n; i++ ) { alphas[i] = alpha;betas[i] = beta;}*/
  const complex_PRECISION *alpha = (const complex_PRECISION*) p_alpha, *beta = (const complex_PRECISION*) p_beta;
  for ( int i=0; i<n; i++ ) {
    alphas[i] = *alpha;
    betas[i]  = *beta;
  }

  vector_PRECISION B, X, B0 = fab->B0, X0 = fab->X0, C0 = fab->C0;
  vector_PRECISION_duplicate( &B, &B0, 0, l ); B.num_vect = N; B.num_vect_now = N; B.size = ldb;
  vector_PRECISION_duplicate( &X, &X0, 0, l ); X.num_vect = N; X.num_vect_now = N; X.size = ldx;
  B.vector_buffer = BB;
  X.vector_buffer = (const buffer_PRECISION) XX;
  if ( fab->dim != X.size || N > X.num_vect )//debugging!!!!
    error0("Fabulous matrix vector product callback: assumptions are not met\n");
  vector_PRECISION_copy_fab( &B0, &B, 0, N, 1, 0, fab->dim, l );
  vector_PRECISION_copy_fab( &X0, &X, 0, N, 1, 0, fab->dim, l );

  // may need to write a new dirac op function for column major order!!!!
  vector_PRECISION_change_layout( &X0, &X0, _NVEC_INNER, no_threading );
  vector_PRECISION_change_layout( &B0, &B0, _NVEC_INNER, no_threading );
  vector_PRECISION_scale_new( &B0, &B0, betas, 0, 0, fab->dim, fab->l );
  p->eval_operator( &C0, &X0, p->op, l, no_threading );
  vector_PRECISION_saxpy_new( &B0, &B0, &C0, alphas, 0, 1, 0, fab->dim, l );
  X0.layout = _NVEC_OUTER;
  vector_PRECISION_change_layout( &B0, &B0, _NVEC_OUTER, no_threading );

  vector_PRECISION_copy_fab( &B, &B0, 0, N, -1, 0, fab->dim, l );
  //  END_UNTHREADED_FUNCTION(fab.threading)
  //    printf0("mvp2 %d: %d vs %d; %d %d %d\n", g.my_rank, N, fab->B0.num_vect, ldx, ldb, fab->B0.size );fflush(stdout);
  PROF_PRECISION_STOP( _FMVP, 1 );
    
  return sizeof(PRECISION)/4*N*ldx*ldx;//this is actually wrong!!!!!;
}
// C = A^{H} * B: Dot product CALLBACK: vectors need to be in _NVEC_INNER 
int64_t dot_product_PRECISION(void *user_env,
			      int M, int N,
			      const void *A_, int lda,
			      const void *B_, int ldb,
			      void *C_, int ldc) {

  fabulous_PRECISION_struct *fab = &(((gmres_PRECISION_struct *) user_env)->fab);
  level_struct *l = fab->l;
  struct Thread *threading = fab->threading;
  buffer_PRECISION A, B, C;
  A = (const buffer_PRECISION) A_;
  B = (const buffer_PRECISION) B_;
  C = C_;

  PROF_PRECISION_START( _FIP );
  int i, j, k, s;
  complex_PRECISION global_alpha[ldc*N];
  //  printf0("dot %d %d, %d %d %d\n", M,N,lda,ldb,ldc);fflush(stdout);
  for ( s=0; s<ldc*N; s++ ) C[s] = 0; 
  for ( s=0; s<ldc*N; s++ ) global_alpha[s] = 0;//printf0("dot0 %d %d %d %d %d %d\n", M, N, lda, ldb,ldc, fab->dim);
#if 0
  int thread_start;
  int thread_end;
  compute_core_start_end_custom(0, fab->dim, &thread_start, &thread_end, l, threading, l->num_lattice_site_var );

  for( i=thread_start; i<thread_end; i++)
#endif
    for ( j=0; j<M; j++ )
      for (k=0; k<N; k++ )
	for ( s=0; s<fab->dim; s++ ) C[k*ldc+j] += conj_PRECISION(A[j*lda+s])*B[k*ldb+s];
  //for ( s=thread_start; s<thread_end; s++ ) C[k*ldc+j] += conj_PRECISION(A[j*lda+s])*B[k*ldb+s];

  /////////////// communication -----------------------------------
#if 0 // Assume: the fabulous solver is called out of OpenMP region
  // sum over cores
  START_NO_HYPERTHREADS(threading)
    for ( s=0; s<M*N; s++ ) ((complex_PRECISION *)threading->workspace)[threading->core*ldc*N+s] = C[ldc*(s/M)+s%M];
  END_NO_HYPERTHREADS(threading)
  // master sums up all results
  SYNC_CORES(threading)
  START_MASTER(threading)
  for( i=1; i<threading->n_core; i++)
    for ( s=0; s<M*N; s++ ) ((complex_PRECISION *)threading->workspace)[s] += ((complex_PRECISION *)threading->workspace)[i*M*N+s];
  END_MASTER(threading)
  SYNC_MASTER_TO_ALL(threading)
    for ( s=0; s<M*N; s++ )  C[ldc*(s/M)+s%M] = ((complex_PRECISION *)threading->workspace)[s];
#endif  
  if ( g.num_processes > 1 ) {
    START_MASTER(threading)
    PROF_PRECISION_START( _ALLR );
    MPI_Allreduce( C, global_alpha, ldc*N, MPI_COMPLEX_PRECISION, MPI_SUM, (l->depth==0)?g.comm_cart:l->gs_PRECISION.level_comm );
    PROF_PRECISION_STOP( _ALLR, 1 );
    for ( s=0; s<M*N; s++ ) ((complex_PRECISION *)threading->workspace)[s] = global_alpha[ldc*(s/M)+s%M];
    END_MASTER(threading)
    // all threads need the result of the norm
    SYNC_MASTER_TO_ALL(threading)
    for ( s=0; s<M*N; s++ ) C[ldc*(s/M)+s%M] = ((complex_PRECISION *)threading->workspace)[s];
    //for ( s=0; s<ldc*N; s++ ) C[s] = global_alpha[s];
  }
  //  printf0("dot end\n");fflush(stdout);
  PROF_PRECISION_STOP( _FIP, 1 );
    
  return 2L*N*M*lda;
}

int64_t fabulous_rightprecond_PRECISION(void *user_env, int N,
					const void *XX, int ldx,
					void *BB, int ldb) {
  
  gmres_PRECISION_struct *p = user_env;
  fabulous_PRECISION_struct *fab = &(p->fab);
  level_struct *l = fab->l;
  struct Thread *threading = fab->threading;

  int n = fab->B0.num_vect;//!!!!!!
  if (l->level > 0 ) {
    if ( g.mixed_precision ) l->next_level->p_float.fab.nrhs = N;//!!!!!!
    else l->next_level->p_double.fab.nrhs = N;//!!!!!!
  }
  vector_PRECISION B, X, B0 = fab->B0, X0 = fab->X0;
  /*
  vector_PRECISION_init( &X0);
  vector_PRECISION_init( &B0);
  vector_PRECISION_alloc( &X0, _ORDINARY, n, l, threading ); X0.num_vect_now = n; X0.layout = _NVEC_OUTER;
  vector_PRECISION_alloc( &B0, _ORDINARY, n, l, threading ); B0.num_vect_now = n; B0.layout = _NVEC_OUTER;
  */
  vector_PRECISION_duplicate( &B, &B0, 0, l ); B.num_vect = N; B.num_vect_now = N; B.size = ldb;
  vector_PRECISION_duplicate( &X, &X0, 0, l ); X.num_vect = N; X.num_vect_now = N; X.size = ldx;
  B.vector_buffer = BB;
  X.vector_buffer = (const buffer_PRECISION) XX;
  
  if ( l->level > 0 ) {  printf0("fab prec PRECISION: %d->%d %d; %d vs %d\n",l->depth,l->next_level->depth,l->is_PRECISION.num_agg,N,n);fflush(stdout);
    vector_PRECISION_copy_fab( &B0, &B, 0, N, 1, 0, fab->dim, l );
    vector_PRECISION_copy_fab( &X0, &X, 0, N, 1, 0, fab->dim, l );
    for ( int i=N; i<n; i++ )  vector_PRECISION_copy_fab( &X0, &X, i, 1, -1, 0, fab->dim, l );
    vector_PRECISION_change_layout( &X0, &X0, _NVEC_INNER, no_threading );
    vector_PRECISION_change_layout( &B0, &B0, _NVEC_INNER, no_threading );
    if ( l->depth == 0 ) preconditioner_new( (vector_double *)(&B0), NULL, (vector_double *) (&X0), _RES, l, threading ); // can also use !g.in_setup as this is used only at top
    else vcycle_PRECISION_new( &B0, NULL, &X0, _RES, l, threading );//printf0("fab prec2: %d\n",l->depth);fflush(stdout);
    //vcycle_PRECISION_new( &B0, NULL, &X0, _NO_RES, l, threading );
    X0.layout = _NVEC_OUTER;
    vector_PRECISION_change_layout( &B0, &B0, _NVEC_OUTER, no_threading );//printfv_PRECISION(&B0);
    vector_PRECISION_copy_fab( &B, &B0, 0, N, -1, 0, fab->dim, l );
  } else {
    vector_PRECISION_copy_fab( &B, &X, 0, N, -1, 0, fab->dim, l );
  }

  //vector_PRECISION_free(&(X0), l, threading);
  //vector_PRECISION_free(&(B0), l, threading);
  printf0("end fab prec PRECISION: %d->%d %d; %d\n",l->depth,l->next_level->depth,l->is_PRECISION.num_agg,N);fflush(stdout);
  return 4; //This is a random number
    
}

// dir == 1: z[0:length] <-x[vec_ind:vec_ind+length]; dir==-1: z[vec_ind:vec_ind+length] <- x[0:length]
static void vector_PRECISION_copy_fab( vector_PRECISION *z, vector_PRECISION *x, int vec_ind, int length, int dir, int start, int end, level_struct *l ) {

  PROF_PRECISION_START( _FAB_COPY );
  int layout = x->layout;
  int i, j;
  if ( layout == _NVEC_INNER ) {
    if ( dir == 1 )
      for( i=start; i<end; i++)
	for( j=0; j<length; j++ ) z->vector_buffer[i*z->num_vect+j] = x->vector_buffer[i*x->num_vect+vec_ind+j];
    else
      for( i=start; i<end; i++)
	for( j=0; j<length; j++ ) z->vector_buffer[i*z->num_vect+vec_ind+j] = x->vector_buffer[i*x->num_vect+j];
  }
  else if ( layout == _NVEC_OUTER ) {
    if ( dir == 1 )
      for( j=0; j<length; j++ )
	for( i=start; i<end; i++) z->vector_buffer[j*z->size+i] = x->vector_buffer[(vec_ind+j)*x->size+i];
    else
      for( j=0; j<length; j++ )
	for( i=start; i<end; i++) z->vector_buffer[(vec_ind+j)*z->size+i] = x->vector_buffer[j*x->size+i];
  }
  PROF_PRECISION_STOP( _FAB_COPY, (double)(end-start)/(double)l->inner_vector_size );
}
