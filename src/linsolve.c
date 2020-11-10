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
 * copied: 11/29/2019 
 * changed from sbacchio
 * checked:12/06/2019
 */

#include "main.h"

void fgmres_MP_struct_init( gmres_MP_struct *p ) {
  fgmres_float_struct_init( &(p->sp) );
  fgmres_double_struct_init( &(p->dp) );
}


void fgmres_MP_struct_alloc( int m, int n, const int vl_type, double tol, const int prec_kind, 
                             void (*precond)(), gmres_MP_struct *p, level_struct *l ) {

/*********************************************************************************
 * Input:
 * Allocates memory for the fgmres struct and sets its values for both precisions in *p.
 * int m: Restart length                                                            
 * int n: Number of restarts                                                        
 * int vl: System size                                                              
 * double tol: Tolerance for relative residual                                         
 * const int type: Specifies the problem for which fgmres should be applied               
 *                 (_GLOBAL_FGMRES, _K_CYCLE, _COARSE_GMRES)                              
 * const int prec_kind: type of preconditioning: _RIGHT (flexible preconditioner),
 *                                               _LEFT (stationary preconditioner)
 *                                               or _NOTHING                                     
 * void (*precond): Function pointer to the preconditioner                           
 * level_struct *l: level structure of the current level
 *
 * Output:
 * gmres_PRECISION_struct *p
*********************************************************************************/  

  long int total=0; 
  int i, k=0, nvec=num_loop;

  // double                                       // single
  p->dp.restart_length = m;                       p->sp.restart_length = m;           
  p->dp.num_restart = n;                          p->sp.num_restart = n;
  p->dp.preconditioner = NULL;                    p->sp.preconditioner = precond;
  p->dp.eval_operator = d_plus_clover_double; p->sp.eval_operator = d_plus_clover_float;

  p->dp.tol = tol;                                p->sp.tol = MAX(tol,1E-5);
  p->dp.kind = _NOTHING;                          p->sp.kind = prec_kind;
  p->dp.timing = 1;                               p->sp.timing = 1;
                                 
  p->dp.print = g.vt.evaluation?0:1;              p->sp.print = g.vt.evaluation?0:1;
  p->dp.initial_guess_zero = 1;                   p->sp.initial_guess_zero = 1;
  p->dp.v_start = 0;                              p->sp.v_start = 0;
  p->dp.v_end = l->inner_vector_size;             p->sp.v_end = l->inner_vector_size;
  
  p->dp.op = &(g.op_double);                      p->sp.op = &(l->s_float.op);
  
  //g.p.op = &(g.op_double);//this showed up in method_init too??????? ony when INIT_PREC is not enabled
  //g.p.eval_operator = d_plus_clover_double;//why here????? ony when INIT_PREC is not enabled
  
#ifdef HAVE_TM1p1
  nvec*=2;
#endif

  //--------------------------- double precision part
  p->dp.num_vect = nvec;
  total = 0;
  total += (m+1)*m*nvec;       // for Hessenberg matrix
  total += 4*(m+1)*nvec;       // for y, gamma, c, s
  p->dp.total_storage = total; // precomputed total storage amount for y, gamma, c, s

  // allocate connected memory for Hessenberg matrix, y, gamma, c, s 
  MALLOC( p->dp.H, complex_double*, m );  
  p->dp.H[0] = NULL; 
  MALLOC( p->dp.H[0], complex_double, total );

  // reserve storage
  total = 0;
  // H
  for ( i=1; i<m; i++ )
    p->dp.H[i] = p->dp.H[0] + i*(m+1)*nvec;
  total += m*(m+1)*nvec;
  // y
  p->dp.y = p->dp.H[0] + total; total += (m+1)*nvec;
  // gamma
  p->dp.gamma = p->dp.H[0] + total; total += (m+1)*nvec;
  // c
  p->dp.c = p->dp.H[0] + total; total += (m+1)*nvec;
  // s
  p->dp.s = p->dp.H[0] + total; total += (m+1)*nvec;

  // x
  vector_double_alloc( &(p->dp.x), vl_type, nvec, l, no_threading );
  // r
  vector_double_alloc( &(p->dp.r), vl_type, nvec, l, no_threading );
  // b
  vector_double_alloc( &(p->dp.b), vl_type, nvec, l, no_threading );
  
  ASSERT( p->dp.total_storage == total );
  
  
  //--------------------------- single precision part
  p->sp.num_vect = nvec;

  // allocate connected memory for V and Z
  MALLOC( p->sp.V, vector_float, m+1 );
  if ( precond != NULL ) {
    if ( prec_kind == _RIGHT ) {
      k = m+1;
    } else {
      k = m;//m debug!!!originally 0
    }
    MALLOC( p->sp.Z, vector_float, k );
  }
  
  // reserve storage for w, V, Z
  total = 0;
  // w
  vector_float_alloc( &(p->sp.w), vl_type, nvec, l, no_threading );
  // V 
  for ( i=0; i<m+1; i++ ) {
    vector_float_init( &(p->sp.V[i]) );
    vector_float_alloc( &(p->sp.V[i]), vl_type, nvec, l, no_threading );
  }
  // Z
  if ( precond != NULL ) {
    for ( i=0; i<k; i++ ) {
      vector_float_init( &(p->sp.Z[i]) );
      vector_float_alloc( &(p->sp.Z[i]), vl_type, nvec, l, no_threading );
    }
  }
  
}  
   
   
void fgmres_MP_struct_free( gmres_MP_struct *p, level_struct *l ) {

  int i, k;  
  if ( p->sp.preconditioner != NULL ) {
    if ( p->sp.kind == _RIGHT ) {
      k = p->dp.restart_length+1;
    } else {
      k = p->dp.restart_length;//m debug!!!originally 0
    }
  }

  // single precision
  vector_float_free( &(p->sp.w), l, no_threading );
  for ( i=0; i<p->dp.restart_length+1; i++ ) vector_float_free( &(p->sp.V[i]), l, no_threading );
  FREE( p->sp.V, vector_float, p->sp.restart_length+1 );
  if ( p->sp.Z != NULL ) {
    for ( i=0; i<k; i++ ) vector_float_free( &(p->sp.Z[i]), l, no_threading );
    FREE( p->sp.Z, vector_float, k);//p->sp.kind==_RIGHT?p->sp.restart_length+1:p->sp.restart_length );// debug: 1->p->sp.restart_length
  }
  
  // double precision
  FREE( p->dp.H[0], complex_double, p->dp.total_storage );
  FREE( p->dp.H, complex_double*, p->dp.restart_length );
  vector_double_free( &(p->dp.x), l, no_threading );
  vector_double_free( &(p->dp.r), l, no_threading );
  vector_double_free( &(p->dp.b), l, no_threading );
} 
    
int fgmres_MP( gmres_MP_struct *p, level_struct *l, struct Thread *threading ) {
  
  /*********************************************************************************
   * Uses FGMRES to solve the system D x = b, where b is taken from p->b and x is 
   * stored in p->x.                                                              
   *********************************************************************************/  
  
  ASSERT( g.mixed_precision );

  // start and end indices for vector functions depending on thread
  int start;
  int end;
  
  int i,jj, j=-1, finish=0, iter=0, il, ol, n_vect=num_loop;//g.num_vect_now;//!!!!!!!
  double beta[n_vect];

  double t0=0, t1=0;
  double norm_r0[n_vect], gamma_jp1[n_vect], gamma0_real[n_vect], gamma_max, gamma_max2, H_tot;
  complex_float gamma_float[n_vect];

  if ( p->sp.num_vect < num_loop )
    error0("fgmres: memory corruption\n");

  p->dp.w.num_vect_now = n_vect; p->dp.x.num_vect_now = n_vect; p->dp.r.num_vect_now = n_vect; p->dp.b.num_vect_now = n_vect;
  p->sp.w.num_vect_now = n_vect; p->sp.x.num_vect_now = n_vect; p->sp.r.num_vect_now = n_vect; p->sp.b.num_vect_now = n_vect;
  //p->sp.Z[0].num_vect_now = n_vect;p->sp.V[0].num_vect_now = n_vect;//!!!!!!!
  p->sp.V[0].num_vect_now = n_vect;

  VECTOR_LOOP(i, n_vect, jj, norm_r0[i+jj]=1;
                             gamma_jp1[i+jj]=1;)
  
  START_LOCKED_MASTER(threading)
#ifndef WILSON_BENCHMARK
  if ( l->depth==0 && ( p->dp.timing || p->dp.print ) ) prof_init( l );
#endif
  if ( l->level==0 && g.num_levels > 1 && g.interpolation ) p->dp.tol = g.coarse_tol;
  if ( l->depth > 0 ) p->dp.timing = 1;
  if ( l->depth == 0 ) t0 = MPI_Wtime();
#if defined(TRACK_RES) && !defined(WILSON_BENCHMARK)
  if ( p->dp.print && g.print > 0 ) printf0("+----------------------------------------------------------+\n");
#endif
  END_LOCKED_MASTER(threading)

  // compute start and end indices for core
  // this puts zero for all other hyperthreads, so we can call functions below with all hyperthreads
  //compute_core_start_end(p->dp.v_start, p->dp.v_end, &start, &end, l, threading);
  compute_core_start_end_custom(p->dp.v_start, p->dp.v_end, &start, &end, l, threading, l->num_lattice_site_var);
  
  // Outer loop in double precision
  for( ol=0; ol<p->dp.num_restart && finish==0; ol++ )  {// if go beyond max #restarts or has converged, exit
    
    //------- iniital setup: Compute r_0 = b − Dx_0 , β := ||r_0||_2 , and v_1 := r_0 /β
    if( ol == 0 && p->dp.initial_guess_zero ) {
      vector_double_copy( &(p->dp.r), &(p->dp.b), start, end, l );
    } else { 
      // if the initial guess is not zero, need to compute b_l-D_l phi_i 
      apply_operator_double( &(p->dp.r), &(p->dp.x), &(p->dp), l, threading );     // compute r <- D*x
      vector_double_minus( &(p->dp.r), &(p->dp.b), &(p->dp.r), start, end, l );// compute r <- b - r = b - D*x
    }
    global_norm_double( gamma0_real, &(p->dp.r), p->dp.v_start, p->dp.v_end, l, threading );// gamma_0 = norm(r)
    START_MASTER(threading)
    VECTOR_LOOP(i, n_vect, jj, p->dp.gamma[i+jj] = (complex_double) gamma0_real[i+jj];)
    END_MASTER(threading)
      
    // report
    if( ol == 0) {
      if (l->depth == 0 && !p->dp.initial_guess_zero) {
	global_norm_double( norm_r0, &(p->dp.b), start, end, l, threading );
	START_MASTER(threading)
	for( i=0; i<n_vect; i++ ) 
	  printf0("| initial guess relative residual (%d):            %le |\n", i, gamma0_real[i]/norm_r0[i]);
	END_MASTER(threading)
	SYNC_MASTER_TO_ALL(threading)
      } else {
	VECTOR_LOOP(i, n_vect, jj, norm_r0[i+jj]= gamma0_real[i+jj];)
      }
    }
#if defined(TRACK_RES) && !defined(WILSON_BENCHMARK)
    else {
      if ( p->dp.print && g.print > 0 ) {
        START_MASTER(threading)
        printf0("+----------------------------------------------------------+\n");
        for( i=0; i<n_vect; i++ ) 
          printf0("| restarting ...          true residual norm (%d): %6e |\n", i, gamma0_real[i]/norm_r0[i] );
        printf0("+----------------------------------------------------------+\n");
        END_MASTER(threading)
      }
    }
#endif
    // turn into float
    trans_float( &(p->sp.V[0]), &(p->dp.r), l->s_float.op.translation_table, l, threading );
    VECTOR_LOOP(i, n_vect, jj, gamma_float[i+jj]= (complex_float) p->dp.gamma[0*n_vect+i+jj];)
    vector_float_real_scale( &(p->sp.V[0]), &(p->sp.V[0]), gamma_float, 0, 1, start, end, l ); // V[0] <- r / gamma_0

    // inner loop in single precision
    for( il=0; il<p->dp.restart_length && finish==0; il++) { // if (dim of the Krylov subspace)>p->dp.restart_length or has converged, exit this loop  
      j = il; iter++;

      // extend the Arnoldi basis, p->sp.V, by one; Hessenberg matrix and sol vector y are of double precision
      arnoldi_step_MP( p->sp.V, p->sp.Z, &(p->sp.w), p->dp.H, p->dp.y, j, p->sp.preconditioner, &(p->sp), l, threading );

      // perhaps need to change here!!!!!
      H_tot=0;
      VECTOR_LOOP(i, n_vect, jj, H_tot += cabs( p->dp.H[j][(j+1)*n_vect+i+jj] );)//????what is this for????
      if ( H_tot > n_vect*1E-15 ) {
        qr_update_double( p->dp.H, p->dp.s, p->dp.c, p->dp.gamma, j, l, threading );
        VECTOR_LOOP(i, n_vect, jj, gamma_jp1[i+jj] = cabs( p->dp.gamma[(j+1)*n_vect+i+jj] );)

        if ( iter%10 == 0 || p->sp.preconditioner != NULL || l->depth > 0 ) {
#if defined(TRACK_RES) && !defined(WILSON_BENCHMARK)
          START_MASTER(threading)
          if ( p->sp.print && g.print > 0 )
            for( i=0; i<n_vect; i++ )
              printf0("| vector %d, approx. rel. res. after  %-6d iterations: %e |\n", i, iter, gamma_jp1[i]/norm_r0[i] );
          END_MASTER(threading)
#endif
        }

        gamma_max = gamma_jp1[0]/norm_r0[0];
	for ( i=1; i<n_vect; i++ ) if ( gamma_max < gamma_jp1[i]/norm_r0[i] ) gamma_max = gamma_jp1[i]/norm_r0[i];
        if( gamma_max < p->dp.tol || gamma_max > 1E+5 ) {
          finish = 1;
          START_MASTER(threading)
            if ( gamma_max > 1E+5 ) printf0("Divergence of fgmres_MP, iter = %d, level=%d\n", iter, l->level );
          END_MASTER(threading)
        }
        gamma_max2 = gamma_jp1[0]/gamma0_real[0];
	for ( i=1; i<n_vect; i++ ) if ( gamma_max2 < gamma_jp1[i]/gamma0_real[i] ) gamma_max2 = gamma_jp1[i]/gamma0_real[i];
	if( gamma_max2 < p->sp.tol ){  
          break;
        }
      } else {
        finish = 1;
      }
    } // end of a single restart
    compute_solution_MP( &(p->sp.w), (p->sp.preconditioner&&p->sp.kind==_RIGHT)?p->sp.Z:p->sp.V,
                         p->dp.y, p->dp.gamma, p->dp.H, j, &(p->sp), l, threading );
                           
    trans_back_float( &(p->dp.r), &(p->sp.w), l->s_float.op.translation_table, l, threading );
    if ( ol == 0 ) {
      vector_double_copy( &(p->dp.x), &(p->dp.r), start, end, l );
    } else {
      vector_double_plus( &(p->dp.x), &(p->dp.x), &(p->dp.r), start, end, l ); // p->dp.x += p->dp.r
    }
  } // end of fgmres
  
  START_LOCKED_MASTER(threading)
  if ( l->depth == 0 ) { t1 = MPI_Wtime(); g.iter_times[0] = t1-t0; g.iter_counts[0] = iter; g.norm_res = gamma_max; }
  END_LOCKED_MASTER(threading)
  
  if ( p->dp.print ) {
#ifdef FGMRES_RESTEST
    apply_operator_double( &(p->dp.r), &(p->dp.x), &(p->dp), l, threading );
    vector_double_minus( &(p->dp.r), &(p->dp.b), &(p->dp.r), start, end, l ); 
    global_norm_double( beta, &(p->dp.r), p->dp.v_start, p->dp.v_end, l, threading );
#else
    VECTOR_LOOP(i, n_vect, jj, beta[i+jj] = creal(gamma_jp1[i+jj]);)
#endif
    START_MASTER(threading)
#if defined(TRACK_RES) && !defined(WILSON_BENCHMARK)
    if ( g.print > 0 ) printf0("+----------------------------------------------------------+\n\n");
#endif
    printf0("+----------------------------------------------------------+\n");
    printf0("|              Final Relative Residuals                    |\n");
    for( i=0; i<n_vect; i++ )
      printf0("| exact relative residual, %d: ||r||/||b|| = %e   |\n",i, beta[i]/norm_r0[i] );
    printf0("+----------------------------------------------------------+\n");
    printf0("|                   Solver Statistics                      |\n");
    printf0("|    FGMRES MP iterations: %-6d coarse average: %-6.2lf   |\n", iter,
            ((double)g.iter_counts[g.num_levels-1])/((double)iter) );
    //printf0("|    elapsed wall clock time: %-8.4lf seconds                |\n", t1-t0 );
    if ( l->depth == 0) for( i=0; i<g.num_levels; i++ ) {
        printf0("|  (depth %d) solver: %d, iterations: %-6d, time below this level: %-8.4lf sec (%04.1lf%%) |\n",
		i, g.solver[i], g.iter_counts[i], g.iter_times[i], 100*(g.iter_times[i]/(t1-t0)) );
    }
    printf0("|  consumed core minutes*: %-8.2le (solve only)           |\n", ((t1-t0)*g.num_processes*MAX(1,threading->n_core))/60.0 );
    printf0("|    max used mem/MPIproc: %-8.2le GB                     |\n", g.max_storage/1024.0 );
    printf0("+----------------------------------------------------------+\n");
    printf0("*: only correct if #MPIprocs*#threads == #CPUs\n\n");
    END_MASTER(threading)
  }

  if ( l->depth == 0 && g.vt.p_end != NULL  ) {
    if ( g.vt.p_end != NULL ) {
      START_LOCKED_MASTER(threading)
      printf0("solve iter: %d\n", iter );
      printf0("solve time: %le seconds\n", t1-t0 );
      g.vt.p_end->values[_SLV_TIME] += (t1-t0)/((double)g.vt.average_over);
      g.vt.p_end->values[_SLV_ITER] += iter/((double)g.vt.average_over);
      g.vt.p_end->values[_CRS_ITER] += (((double)g.coarse_iter_count)/((double)iter))/((double)g.vt.average_over);
      g.vt.p_end->values[_CRS_TIME] += g.coarse_time/((double)g.vt.average_over);
    END_LOCKED_MASTER(threading)
    }
  }

  if ( l->depth == 0 && ( p->dp.timing || p->dp.print ) && !(g.vt.p_end != NULL )  ) {
    START_MASTER(threading)
#ifndef WILSON_BENCHMARK
    prof_print( l );
#endif
    END_MASTER(threading)
  }
  
  return iter;
}

void arnoldi_step_MP( vector_float *V, vector_float *Z, vector_float *w,
                      complex_double **H, complex_double* buffer, int j, void (*prec)(),
                      gmres_float_struct *p, level_struct *l, struct Thread *threading ) {

  int i, jj, n, n_vect = w->num_vect_now;//g.num_vect_now, n, jj;//!!!!!!!
  double H_tot;
  complex_float H_float[n_vect];
  // start and end indices for vector functions depending on thread
  int start;
  int end;
  // compute start and end indices for core
  // this puts zero for all other hyperthreads, so we can call functions below with all hyperthreads
  compute_core_start_end_custom(p->v_start, p->v_end, &start, &end, l, threading, l->num_lattice_site_var);

  V[j].num_vect_now = n_vect; V[j+1].num_vect_now = n_vect;//can move to fgmres????   

  //--- apply D (and preconditioner) to V[j]
  if ( prec != NULL ) {
    Z[j].num_vect_now = n_vect;
    if ( p->kind == _LEFT ) {
      apply_operator_float( &Z[0], &V[j], p, l, threading );
      prec( w, NULL, &Z[0], _NO_RES, l, threading );
    } else {
      if ( g.mixed_precision == 2 && (g.method >= 1 && g.method <= 2) ) {
        prec( &Z[j], w, &V[j], _NO_RES, l, threading );    // obtains w = D * Z[j] from Schwarz
      } else {
        prec( &Z[j], NULL, &V[j], _NO_RES, l, threading );
        apply_operator_float( w, &Z[j], p, l, threading ); // w = D*Z[j]
      }
    }
  } else {
    apply_operator_float( w, &V[j], p, l, threading );     // w = D*V[j]
  }

  complex_double tmp[(j+1)*n_vect];
  process_multi_inner_product_MP( j+1, tmp, V, w, p->v_start, p->v_end, l, threading );
  START_MASTER(threading)
  for( i=0; i<=j; i++ )
    VECTOR_LOOP(n, n_vect, jj, buffer[i*n_vect+n+jj] = tmp[i*n_vect+n+jj];)

  if ( g.num_processes > 1 ) {
    PROF_double_START( _ALLR );
    MPI_Allreduce( buffer, H[j], (j+1)*n_vect, MPI_COMPLEX_double, MPI_SUM, (l->depth==0)?g.comm_cart:l->gs_float.level_comm );
    PROF_double_STOP( _ALLR, 1 );
  } else {
    for( i=0; i<=j; i++ )
      VECTOR_LOOP(n, n_vect, jj, H[j][i*n_vect+n+jj] = buffer[i*n_vect+n+jj];)
  }
  END_MASTER(threading)
  SYNC_MASTER_TO_ALL(threading)

  // orthogonalization
  complex_float alpha[(j+1)*n_vect]; 
  for( i=0; i<=j; i++ )
    VECTOR_LOOP(n, n_vect, jj, alpha[i*n_vect+n+jj] = (complex_float) H[j][i*n_vect+n+jj];)
  vector_float_multi_saxpy( w, V, alpha, -1, j+1, start, end, l, threading );

  // V_j+1 = w / H_j+1,j
  double tmp2[n_vect];
  global_norm_MP( tmp2, w, p->v_start, p->v_end, l, threading );
  START_MASTER(threading)
  VECTOR_LOOP(n, n_vect, jj, H[j][(j+1)*n_vect+n+jj] = tmp2[n+jj];)
  END_MASTER(threading)
  SYNC_MASTER_TO_ALL(threading)
  VECTOR_LOOP(n, n_vect, jj,  H_float[n+jj] = ( cabs_double( H[j][(j+1)*n_vect+n+jj] ) > 1e-15 )?(float)tmp2[n+jj]:1;)
  vector_float_real_scale( &V[j+1], w, H_float, 0, 1, start, end, l );
}

void compute_solution_MP( vector_float *x, vector_float *V, complex_double *y,
                          complex_double *gamma, complex_double **H, int j,
                          gmres_float_struct *p, level_struct *l, struct Thread *threading ) {
  
  int i, k, n, jj, n_vect=num_loop;//g.num_vect_now;
  complex_float y_float[n_vect];
  // start and end indices for vector functions depending on thread
  int start;
  int end;
  // compute start and end indices for core
  // this puts zero for all other hyperthreads, so we can call functions below with all hyperthreads
  compute_core_start_end_custom(p->v_start, p->v_end, &start, &end, l, threading, l->num_lattice_site_var );

  START_MASTER(threading)
  
  PROF_double_START( _SMALL2 );
  
  // backward substitution
  for ( i=j; i>=0; i-- ) {
    VECTOR_LOOP(n, n_vect, jj, y[i*n_vect+n+jj] = gamma[i*n_vect+n+jj];)
    for ( k=i+1; k<=j; k++ )
      VECTOR_LOOP(n, n_vect, jj, y[i*n_vect+n+jj] -= H[k][i*n_vect+n+jj]*y[k*n_vect+n+jj];)
    VECTOR_LOOP(n, n_vect, jj, y[i*n_vect+n+jj] /= H[i][i*n_vect+n+jj];)
  }
  
  PROF_double_STOP( _SMALL2, ((j+1)*(j+2))/2 + j+1 );
  
  END_MASTER(threading)
  SYNC_MASTER_TO_ALL(threading)
  
  // x = V*y
  VECTOR_LOOP(n, n_vect, jj, y_float[n+jj]= (complex_float) y[0*n_vect+n+jj];)
  vector_float_scale( x, &V[0], y_float, 0, start, end, l ); // x = y_0 * V_0 where V_0 is the first basis vector 

  complex_float alpha[j*n_vect];
  for ( i=1; i<=j; i++ )
    VECTOR_LOOP(n, n_vect, jj, alpha[(i-1)*n_vect+n+jj] = (complex_float) y[i*n_vect+n+jj];)
  vector_float_multi_saxpy( x, &V[1], alpha, 1, j, start, end, l, threading );

}
