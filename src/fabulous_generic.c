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
 */

#ifdef HAVE_FABULOUS
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
  fab->nrhs = 0;
  fab->dim = 0;
  fab->ldb = 0;
  fab->ldx = 0;
  fab->ldu = 0;
  fab->max_iter = 0;
  fab->mvp = 0;
  fab->dsize= 0;
  fab->U = NULL;
  fab->l = NULL;
  fab->threading = NULL;
  fab->handle = NULL;
  
}

void setup_fabulous_PRECISION( gmres_PRECISION_struct *p, int v_type, level_struct *l, struct Thread *threading ) {
  /*
   * Set fabulous_PRECISION_struct
   */

  fabulous_PRECISION_struct *fab = &(p->fab);
  
  int dim = (g.odd_even&&l->level==0)?(l->num_inner_lattice_sites/2)*l->num_lattice_site_var:p->v_end, nrhs = p->num_vect;
  fab->nrhs     = num_loop;
  fab->dim      = dim;
  fab->ldb      = dim;
  fab->ldx      = dim;
  fab->mvp      = (nrhs*p->restart_length*p->num_restart < g.max_mvp[l->depth])? g.max_mvp[l->depth]:nrhs*p->restart_length*p->num_restart;
  fab->max_iter = nrhs*p->restart_length;
  fab->l        = l;
  fab->threading = threading;
  
#ifdef HAVE_TM1p1
  // In this case, fabulous is used only in the solver phase
  if ( g.epsbar != 0 || g.epsbar_ig5_odd_shift != 0 || g.epsbar_ig5_odd_shift != 0 )
    dim *= 2;
#endif
  
  vector_PRECISION *X = &(fab->X) , *B = &(fab->B), *B0 = &(fab->B0), *X0 = &(fab->X0), *C0 = &(fab->C0);
  // The following fields are fed to fabulous solver; visible outsdie of fabulous; perhaps will be unnecessary if fabulous can handle _ORINARY properly!!!!
  vector_PRECISION_alloc( X, (g.odd_even&&l->level==0)?_EVEN_INNER:_INNER, nrhs, l, threading ); X->num_vect_now = nrhs;
  vector_PRECISION_alloc( B, (g.odd_even&&l->level==0)?_EVEN_INNER:_INNER, nrhs, l, threading ); B->num_vect_now = nrhs;
  // The following fields are used in mvp_PRECISION as temporay storage; only visible inside of fabulous
  vector_PRECISION_alloc( X0, v_type, nrhs, l, threading ); X0->layout = _NVEC_OUTER; X0->num_vect_now = nrhs;
  vector_PRECISION_alloc( B0, v_type, nrhs, l, threading ); B0->layout = _NVEC_OUTER; B0->num_vect_now = nrhs;
  vector_PRECISION_alloc( C0, v_type, nrhs, l, threading ); C0->layout = _NVEC_OUTER; C0->num_vect_now = nrhs;
  if ( X->size != fab->dim )
    error0("set_fabulous_struct_PRECISION: assumptions are not met\n");

  int sol = g.solver[l->depth];
  fab->k = g.k[l->depth];
  if ( fab->k > 0 && ( sol==_GCRO || sol==_DR || sol==_IBDR || sol==_QRDR || sol==_QRIBDR ) ) {
    MALLOC( fab->eigvals, complex_PRECISION, fab->k );
    for ( int i=0; i<fab->k; i++ ) ((complex_PRECISION *) fab->eigvals)[i] = 0;
  }

  //----- Set fabulous handle
  fab->handle = fabulous_create(FABULOUS_COMPLEX_PRECISION, dim, (void *) p);
  if ( fab->handle == NULL )
    error0("fabulous handle has not been created!\n");

  // Set callbacks:
  fabulous_set_mvp(&mvpf_PRECISION, fab->handle);
  fabulous_set_dot_product(&dot_product_PRECISION, fab->handle);
  if ( g.solver[l->depth] && l->level > 0 )
    fabulous_set_rightprecond(&fabulous_rightprecond_PRECISION, fab->handle);
  if ( p->print && g.print > 0 ) {
    if ( g.solver[l->depth] != _GCR ) g.comp_residual[l->depth] = 1;
    fabulous_set_callback(&fabulous_print_PRECISION, fab->handle);
  }
  
  // Setup parameters:
  PRECISION tolerance[1] = { p->tol };
  fabulous_set_ortho_process(g.f_orthoscheme[l->depth], g.f_orthotype[l->depth], g.ortho_iter[l->depth], fab->handle);
  fabulous_set_parameters( fab->mvp, fab->max_iter, tolerance, 1, fab->handle );
  fabulous_set_advanced_parameters( g.max_kept_direction[l->depth], g.comp_residual[l->depth], g.logger_user_data_size, g.quiet, fab->handle );

}

void fabulous_PRECISION_free( fabulous_PRECISION_struct *fab, level_struct *l, struct Thread *threading ) {

  vector_PRECISION_free(&(fab->X), l, threading);
  vector_PRECISION_free(&(fab->B), l, threading);
  vector_PRECISION_free(&(fab->X0), l, threading);
  vector_PRECISION_free(&(fab->B0), l, threading);
  vector_PRECISION_free(&(fab->C0), l, threading);

  int sol = g.solver[l->depth];
  if ( fab->k > 0 && ( sol==_GCRO || sol==_DR || sol==_IBDR || sol==_QRDR || sol==_QRIBDR ))
     PUBLIC_FREE( fab->eigvals, complex_PRECISION, fab->k );
  if(fab->U != NULL) free(fab->U);
  
  fabulous_destroy(fab->handle);
  fab->handle = NULL;
}

void reset_fab_nrhs_PRECISION( level_struct *l ) {
  // This will have no effect as the only fabulous flexible solver does not use inexact breakdown.  just for future update
  int nrhs = num_loop;
  level_struct *l_tmp = l;
#ifdef HAVE_TM1p1
  if ( g.n_flavours > 1 && g.epsbar == 0 && g.epsbar_ig5_odd_shift == 0 && g.epsbar_ig5_odd_shift == 0 )
    nrhs *= g.n_flavours;
#endif
  while ( l_tmp != NULL ) {
    if ( l_tmp->p_PRECISION.fab.nrhs != 0 )
      l_tmp->p_PRECISION.fab.nrhs = nrhs;
    l_tmp = l_tmp->next_level;
  }
}

/*********************  INSIDE OF FABULOUS CALL  **************************************************/
// B <- beta*B+ alpha*A*X: Matrix x BlockOfVector product CALLBACK: vectors need to be in _NVEC_INNER
int64_t mvpf_PRECISION(  void *user_env, int N,
			 const void *p_alpha, const void *XX, int ldx,
			 const void *p_beta, void *BB, int ldb) {
  
  
  gmres_PRECISION_struct *p      = (gmres_PRECISION_struct *) user_env;
  fabulous_PRECISION_struct *fab = &(p->fab);
  level_struct *l                = fab->l;
  struct Thread *threading       = fab->threading;
  int n = fab->X0.num_vect;
  complex_PRECISION alphas[n], betas[n];
  const complex_PRECISION *alpha = (const complex_PRECISION*) p_alpha, *beta = (const complex_PRECISION*) p_beta;
  
  PROF_PRECISION_START( _FMVP, threading );

#ifdef DEBUG
  //if ( N > fab->X0.num_vect )
    //error0("mvpf_PRECISION: assumptions are not met (%d %d)\n",N,fab->X0.num_vect);
#ifdef HAVETM1p1
  if ( g.n_flavours/g.num_indep_flav*fab->dim != ldx || g.n_flavours/g.num_indep_flav*fab->dim != ldb || g.n_flavours/g.num_indep_flav*N > fab->X0.num_vect  )
    error0("mvpf_PRECISION: assumptions are not met\n");
#endif
#endif
  
  for ( int i=0; i<n; i++ ) {
    alphas[i] = *alpha;
    betas[i]  = *beta;
  }

  int start, end;
  compute_core_start_end_custom(0, fab->dim, &start, &end, l, threading, l->num_lattice_site_var);

  vector_PRECISION B, X, B0 = fab->B0, X0 = fab->X0, C0 = fab->C0;
  vector_PRECISION_duplicate( &B, &B0, 0, l ); B.size = fab->dim; B.num_vect = N;
  vector_PRECISION_duplicate( &X, &X0, 0, l ); X.size = fab->dim; X.num_vect = N;
  START_MASTER(threading)
  B.vector_buffer = (buffer_PRECISION) BB;
  X.vector_buffer = (const buffer_PRECISION) XX;
  END_MASTER(threading)

    int ns = n;
#ifdef HAVE_TM1p1
  ns /= g.n_flavours/g.num_indep_flav;
#endif
  for ( int vi = 0; vi<N; vi += ns ) {
    int s = (((((vi+ns)/N)*2-1)%2+1)/2); // <=> vi+ns > N?1:0;
    int d = (vi+ns-N)*s;//((vi+ns)/N)
    //if(N>ns)printf0("fab mat: vi=%d d=%d ns=%d N=%d\n",vi,d,ns,N);
    //    else warning0("fab small!!!%d %d\n",N,d);
    vector_PRECISION_copy_fab( &B0, &B, vi, ns-d, 1, start, end, l );
  vector_PRECISION_copy_fab( &X0, &X, vi, ns-d, 1, start, end, l );
  //for ( int i=0; i<d; i++ ) vector_PRECISION_copy_fab( &X0, &X, N%ns+i, 1, -1, start, end, l );
  //printf0("matvf: dim=%d ffac=%d; N=%d vs  (B0,X0,C0)_now=( %d,%d,%d); size ldb,ldb=%d %d vs %d %d \n",dim, g.n_flavours/g.num_indep_flav,N, B0.num_vect_now,X0.num_vect_now,C0.num_vect_now,ldx,ldb,B0.size,X0.size);

  // may need to write a new dirac op function for column major order!!!!
  START_LOCKED_MASTER(threading)
  vector_PRECISION_change_layout( &X0, &X0, _NVEC_INNER, no_threading );
  vector_PRECISION_change_layout( &B0, &B0, _NVEC_INNER, no_threading );
  C0.layout = _NVEC_INNER;
  END_LOCKED_MASTER(threading)
  vector_PRECISION_scale( &B0, &B0, betas, 0, start, end, fab->l );
  p->eval_operator( &C0, &X0, p->op, l, no_threading ); // if(g.odd_even&&l->level==0), updates only even part
  vector_PRECISION_saxpy( &B0, &B0, &C0, alphas, 0, 1, start, end, l );
  START_LOCKED_MASTER(threading)
  X0.layout = _NVEC_OUTER;
  vector_PRECISION_change_layout( &B0, &B0, _NVEC_OUTER, no_threading );
  C0.layout = _NVEC_OUTER;
  END_LOCKED_MASTER(threading)
    vector_PRECISION_copy_fab( &B, &B0, vi, ns-d, -1, start, end, l );
  }
  //printf0("mvp2 %d: %d vs %d; %d %d %d\n", g.my_rank, N, fab->B0.num_vect, ldx, ldb, fab->B0.size );fflush(stdout);
  PROF_PRECISION_STOP( _FMVP, 1, threading );
    
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
  C = C_;//Is C column major?????

  PROF_PRECISION_START( _FIP, threading );
  int i, j, k, s;
  complex_PRECISION global_alpha[ldc*N];

  for ( s=0; s<ldc*N; s++ ) C[s] = 0; 
  for ( s=0; s<ldc*N; s++ ) global_alpha[s] = 0;
  
  int dim = fab->dim;
#ifdef HAVE_TM1p1
  dim *= g.n_flavours/g.num_indep_flav;
#ifdef DEBUG
  if(dim > lda )
    error0("fabulous dot_product_PRECISION: assumptions are not met %d %d %d\n",dim,lda,ldb);
#endif
#endif
  int start;
  int end;
  compute_core_start_end_custom(0, dim, &start, &end, l, threading, l->num_lattice_site_var );
  //printf0("dot fab: %d %d %d, dim%d lda%d lfb%d %d: %d<%d?\n",M,N,fab->nrhs,dim,lda,ldb,ldc, M,ldc);
  for (k=0; k<N; k++ ) // assume column major order of C
    for ( j=0; j<M; j++ )
      for ( s=start; s<end; s++ )
	C[k*ldc+j] += conj_PRECISION(A[j*lda+s])*B[k*ldb+s];

  /////////////// communication -----------------------------------
  // sum over cores
  START_NO_HYPERTHREADS(threading)
  for ( s=0; s<M*N; s++ ) ((complex_PRECISION *)threading->workspace)[threading->core*M*N+s] = C[ldc*(s/M)+s%M];
  END_NO_HYPERTHREADS(threading)
  // master sums up all results
  SYNC_CORES(threading)
  START_MASTER(threading)
  for( i=1; i<threading->n_core; i++)
    for ( s=0; s<M*N; s++ ) ((complex_PRECISION *)threading->workspace)[s] += ((complex_PRECISION *)threading->workspace)[i*M*N+s];
  END_MASTER(threading)
  SYNC_MASTER_TO_ALL(threading)
  for ( s=0; s<M*N; s++ )  C[ldc*(s/M)+s%M] = ((complex_PRECISION *)threading->workspace)[s];
  if ( g.num_processes > 1 ) {
    START_MASTER(threading)
    PROF_PRECISION_START( _FALLR );
    MPI_Allreduce( C, global_alpha, ldc*N, MPI_COMPLEX_PRECISION, MPI_SUM, (l->depth==0)?g.comm_cart:l->gs_PRECISION.level_comm );
    PROF_PRECISION_STOP( _FALLR, 1 );
    for ( s=0; s<M*N; s++ ) ((complex_PRECISION *)threading->workspace)[s] = global_alpha[ldc*(s/M)+s%M];
    END_MASTER(threading)
    // all threads need the result of the norm
    SYNC_MASTER_TO_ALL(threading)
    for ( s=0; s<M*N; s++ ) C[ldc*(s/M)+s%M] = ((complex_PRECISION *)threading->workspace)[s];
  }
  PROF_PRECISION_STOP( _FIP, 1, threading );
  
  return 2L*N*M*lda;
}

// B = M*X = P D_c^{-1} R X
int64_t fabulous_rightprecond_PRECISION(void *user_env, int N,
					const void *XX, int ldx,
					void *BB, int ldb) {

  gmres_PRECISION_struct *p      = user_env;
  fabulous_PRECISION_struct *fab = &(p->fab);
  level_struct *l                = fab->l;
  struct Thread *threading       = fab->threading;

#ifdef DEBUG
  if ( N > fab->X0.num_vect )
    error0("fabulous_rightprecond_PRECISION: assumptions are not met (%d %d)\n",N,fab->X0.num_vect);
#ifdef HAVE_TM1p1
  if ( g.n_flavours/g.num_indep_flav*fab->dim != ldx || g.n_flavours/g.num_indep_flav*fab->dim != ldb )
    error0("mvpf_PRECISION: assumptions are not met\n");
#endif
#endif
  
  vector_PRECISION B, X, B0 = fab->B0, X0 = fab->X0;
  vector_PRECISION_duplicate( &B, &B0, 0, l ); B.num_vect_now = N; B.size = fab->dim; 
  vector_PRECISION_duplicate( &X, &X0, 0, l ); X.num_vect_now = N; X.size = fab->dim; 
  
  B.vector_buffer = BB;
  X.vector_buffer = (const buffer_PRECISION) XX;

  int res = _NO_RES, n = fab->B0.num_vect, start, end;
#ifdef HAVE_TM1p1
  n /= g.n_flavours/g.num_indep_flav;
#endif
  //printf0("fab PRECISION prec: at %d nrhs=%d vs %d, nvec(%d,%d): layout(%d %d): ld (B,X)=%d,%d vs %d, size %d, %d\n",l->depth,N, n,X0.num_vect_now,B0.num_vect_now,B.layout,X.layout,ldx,ldb,dim, B0.size,X0.size);fflush(stdout);

  if ( l->level > 0 ) {
    
    compute_core_start_end_custom(0, fab->dim, &start, &end, l, threading, l->num_lattice_site_var );    
    vector_PRECISION_copy_fab( &X0, &X, 0, N, 1, start, end, l );
    for ( int i=N; i<n; i++ )  vector_PRECISION_copy_fab( &X0, &X, i, 1, -1, start, end, l );// just to avoid overall divergence due to unused vector fields in FGMRES
    START_LOCKED_MASTER(threading)
    vector_PRECISION_change_layout( &X0, &X0, _NVEC_INNER, no_threading );
    B0.layout = _NVEC_INNER;
    END_LOCKED_MASTER(threading)
      
    if ( l->depth == 0 ) preconditioner( (vector_double *)(&B0), NULL, (vector_double *) (&X0), res, l, threading ); // can also use !g.in_setup as this is used only at top
    else vcycle_PRECISION( &B0, NULL, &X0, res, l, threading );//printf0("fab prec2: %d\n",l->depth);fflush(stdout);

#if 0
    PRECISION norm[n];
    vector_PRECISION C0 = fab->C0;
    vector_PRECISION_minus(&C0,&X0,&B0,0,fab->dim,l);
    global_norm_PRECISION(norm,&C0,0,fab->dim,l,threading);
    for(int i=0;i<n;i++)printf0("fab norm %g\n",norm[i]);
#endif
    
    START_LOCKED_MASTER(threading)
    X0.layout = _NVEC_OUTER;
    vector_PRECISION_change_layout( &B0, &B0, _NVEC_OUTER, no_threading );//printfv_PRECISION(&B0);
    END_LOCKED_MASTER(threading)
    vector_PRECISION_copy_fab( &B, &B0, 0, N, -1, start, end, l );
    
  } else {
    error0("fabulous_rightprecond_PRECISION: should not be called at the bottom\n");
  }
  //  printf0("fab prec: at %d nrhs=%d vs %d, nvec(%d,%d): layout(%d %d): ld (B,X)=%d,%d vs %d\n",l->depth,N, n,X0.num_vect_now,B0.num_vect_now,B.layout,X.layout,ldx,ldb,dim);fflush(stdout);

  return n; //This is a random number
    
}

void fabulous_print_PRECISION(void *user_env,
                              int iter, int NRHS,
                              const void *XX, int ldx,
                              const void *RR, int ldr,
                              fabulous_handle handle) {

  gmres_PRECISION_struct *p      = user_env;
  fabulous_PRECISION_struct *fab = &(p->fab);
  level_struct *l                = fab->l;
  struct Thread *threading       = fab->threading;
  //printf0("print fab: %d NRHS=%d\n",l->depth,NRHS);
  int i, n = fab->X0.num_vect;
#ifdef HAVE_TM1p1
  n /= g.n_flavours/g.num_indep_flav;
#endif
  int start, end;
  compute_core_start_end_custom(0, fab->dim, &start, &end, l, threading, l->num_lattice_site_var);

  PRECISION res[n], norm[n];
  //printf0("fab print: %d %d, %d\n",ldx,ldr,NRHS);
  vector_PRECISION R, X, C0 = fab->C0;
  vector_PRECISION_duplicate( &R, &C0, 0, l ); R.num_vect_now = NRHS; R.size = fab->dim;
  vector_PRECISION_duplicate( &X, &C0, 0, l ); X.num_vect_now = NRHS; X.size = fab->dim;
  START_MASTER(threading)
  R.vector_buffer = (const buffer_PRECISION) RR;
  X.vector_buffer = (const buffer_PRECISION) XX;
  END_MASTER(threading)

#ifdef DEBUG
    ASSERT( RR != NULL && XX != NULL );
    if ( R.layout != _NVEC_OUTER || X.layout != _NVEC_OUTER )
      error0("fabulous_print_PRECISION: assumptions are not met\n");
#endif
  
  if ( p->print && g.print > 0 ) {
    vector_PRECISION_copy_fab( &C0, &R, 0, NRHS, 1, start, end, l );
    START_LOCKED_MASTER(threading)
    vector_PRECISION_change_layout( &C0, &C0, _NVEC_INNER, no_threading );
    END_LOCKED_MASTER(threading)
    global_norm_PRECISION(res, &C0, start, end, l, threading );
    
    vector_PRECISION_copy_fab( &C0, &X, 0, NRHS, 1, start, end, l );
    START_LOCKED_MASTER(threading)
    C0.layout = _NVEC_OUTER;
    vector_PRECISION_change_layout( &C0, &C0, _NVEC_INNER, no_threading );
    C0.layout = _NVEC_OUTER;
    END_LOCKED_MASTER(threading)
    global_norm_PRECISION(norm, &C0, start, end, l, threading );
    
    for( i=0; i<NRHS; i++ )
      printf0("| vector %d, depth: %d, approx. rel. res. after  %-6d iterations: %e |\n", i, l->depth, iter, res[i]/norm[i]);
  }
}

// dir == 1: z[0:length] <-x[vec_ind:vec_ind+length]; dir==-1: z[vec_ind:vec_ind+length] <- x[0:length]
static void vector_PRECISION_copy_fab( vector_PRECISION *z, vector_PRECISION *x, int vec_ind, int length, int dir, int start, int end, level_struct *l ) {
  // In case of HAVE_TM1p1, lenght & start,end should not include flavour d.o.f
  PROF_PRECISION_START( _FAB_COPY );
  int layout = x->layout, fdof = 1;
  int i, j;
#ifdef HAVE_TM1p1
  fdof *= g.n_flavours/g.num_indep_flav; // flavour d.o.f.
#endif
  if ( layout == _NVEC_INNER ) {
    if ( dir == 1 )
      for( i=start; i<end; i++)
	for( j=0; j<length*fdof; j++ ) z->vector_buffer[i*z->num_vect+j%length+z->num_vect/2*(j/length)]
					 = x->vector_buffer[i*x->num_vect+vec_ind+j%length+x->num_vect/2*(j/length)];
    else
      for( i=start; i<end; i++)
	for( j=0; j<length*fdof; j++ ) z->vector_buffer[i*z->num_vect+vec_ind+j%length+z->num_vect/2*(j/length)]
					 = x->vector_buffer[i*x->num_vect+j%length+x->num_vect/2*(j/length)];
  }
  else if ( layout == _NVEC_OUTER ) {
    if ( dir == 1 )
      for( j=0; j<length; j++ )
	for( i=start*fdof; i<end*fdof; i++) z->vector_buffer[j*z->size*fdof+i] = x->vector_buffer[(vec_ind+j)*x->size*fdof+i];
    else
      for( j=0; j<length; j++ )
	for( i=start*fdof; i<end*fdof; i++) z->vector_buffer[(vec_ind+j)*z->size*fdof+i] = x->vector_buffer[j*x->size*fdof+i];
  }
  PROF_PRECISION_STOP( _FAB_COPY, (double)(end-start)/(double)l->inner_vector_size );
}

#endif
