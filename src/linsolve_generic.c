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
 * copied:11/28/2019
 * changed from sbacchio
 * checked: 12/06/2019
 * checked: 12/08/2019
 */

#include "main.h"


void fgmres_PRECISION_struct_init( gmres_PRECISION_struct *p ) {

/*********************************************************************************
* Initializes all declared pointers with NULL.                              
*********************************************************************************/

  p->V = NULL;
  p->Z = NULL;
  vector_PRECISION_init(&(p->x));
  vector_PRECISION_init(&(p->b));
  vector_PRECISION_init(&(p->r));
  vector_PRECISION_init(&(p->w));
  p->H     = NULL;
  p->y     = NULL;
  p->gamma = NULL;
  p->c     = NULL;
  p->s     = NULL;
  p->preconditioner = NULL;
  p->eval_operator  = NULL;
  fabulous_PRECISION_init( &(p->fab) );
}


void fgmres_PRECISION_struct_alloc( int m, int n, const int vl_type, PRECISION tol, const int type, const int prec_kind,
                                    void (*precond)(), void (*eval_op)(), gmres_PRECISION_struct *p, level_struct *l ) {

/*********************************************************************************
* Allocates memory for the fgmres struct and sets its values.                  
* Input:
* int m: Restart length                                                            
* int n: Number of restarts                                                        
* int vl_type: Type of a vector 
*              (_ORDINARY, _SCHWARZ, _ODDEVEN, _INNER)
* PRECISION tol: Tolerance for relative residual                                         
* const int type: Specifies the problem for which fgmres should be applied               
*                 (_GLOBAL_FGMRES, _K_CYCLE, _COARSE_GMRES)                              
* const int prec_kind: type of preconditioning: _RIGHT (flexible preconditioner),
*                                               _LEFT (stationary preconditioner)
*                                               _NOTHING                                     
* void (*precond): Function pointer to the preconditioner                      
* level_struct *l: level structure of the current level
*
* Output:
* gmres_PRECISION_struct *p
*   V: a orthonormalized Krylov basis (max size = m)
*   Z:= D*V
*   x: iterative guess
*   w:= D*x
*   b: rhs
*   r: residual
*   H: Hessenberg matrix (complex_PRECISION**: the first pt points to the consecutive memory for the 1st column of H (upper triangular), y, gamma, c, s 
*   y: x expressed in the basis of V
*   gamma: rhs of projected linear eq. onto the Krylov subspace after Givens rotation
*   c: cosine of Givens rotation
*   s: sine of Givens rotation
*********************************************************************************/  

  long int total=0; // keep track of size of H allocated for H, y, gamma, c, s 
  int i, k=0, nvec=num_loop;
#ifdef HAVE_TM1p1
  nvec*=2;
#endif
  
  p->restart_length = m;
  p->num_restart    = n;
  p->preconditioner = precond;
  p->eval_operator  = eval_op; 
  p->tol            = tol;
  p->kind           = prec_kind;
  p->num_vect       = nvec;

  if ( type == _GLOBAL_FSOLVER ) {
    p->timing = 1;
    p->print = g.vt.evaluation?0:1;
    p->initial_guess_zero = 1;
    p->v_start = 0;
    p->v_end = l->inner_vector_size;
    p->op = &(g.op_PRECISION);
    if ( g.solver[l->depth] )
      setup_fabulous_PRECISION( p, vl_type, l, no_threading );//!!!!!!no_threaing->threading desireble 
  } else if ( type == _K_CYCLE ) {
    // these settings also work for GMRES as a smoother
    p->timing = 0;
    p->print = 0;//just for debugging orig0
    p->initial_guess_zero = 1;
    p->v_start = 0;
    p->v_end = l->inner_vector_size;
    p->op = &(l->s_PRECISION.op);
    if ( g.solver[l->depth] )
      setup_fabulous_PRECISION( p, vl_type, l, no_threading );//!!!!!!no_threaing->threading desireble 
  } else if ( type == _COARSE_SOLVER ) {
    p->timing = 0;
    p->print = 0;
    p->initial_guess_zero = 1;
    p->layout = -1;
    p->v_start = 0;
    p->v_end = l->inner_vector_size;
    if ( g.odd_even )
      p->op = &(l->oe_op_PRECISION);
    else  
      p->op = &(l->s_PRECISION.op);
    if ( g.solver[l->depth] )
      setup_fabulous_PRECISION( p, vl_type, l, no_threading );//!!!!!!no_threaing->threading desireble
  } else {
    ASSERT( type < 6 );
  }

  //---- allocate vector fields
  if(( !g.solver[l->depth] || g.use_only_fgrmes_at_setup) && m > 0) {
    // allocate space for H, V, Z
    total += (m+1)*m*nvec; // for Hessenberg matrix
    MALLOC( p->H, complex_PRECISION*, m );
    
    MALLOC( p->V, vector_PRECISION, m+1 );
    
    if ( precond != NULL ) {
      if ( prec_kind == _RIGHT ) {
	k = m+1;
      } else {
	k = 1;//?????
      }
      MALLOC( p->Z, vector_PRECISION, k );
    } else {
#if defined(SINGLE_ALLREDUCE_ARNOLDI) && defined(PIPELINED_ARNOLDI)
      if ( l->level == 0 && l->depth > 0 ) {
	k = m+2;
	MALLOC( p->Z, vector_PRECISION, k );
      }
#else
      k = m;//0;//temp fix!!!!!!!
      MALLOC( p->Z, vector_PRECISION, k );//tempfix!!!!!!!
#endif
    }

    // allocate connected memory
    total += 4*(m+1)*nvec; // for y, gamma, c, s
    p->H[0] = NULL; 
    MALLOC( p->H[0], complex_PRECISION, total );
    
    p->total_storage = total;
    total = 0; // for sanity check
    
    // ordering: H, y, gamma, c, s
    // H
    for ( i=1; i<m; i++ )
      p->H[i] = p->H[0] + i*(m+1)*nvec;
    total += m*(m+1)*nvec;
    // y
    p->y = p->H[0] + total; total += (m+1)*nvec;
    // gamma
    p->gamma = p->H[0] + total; total += (m+1)*nvec;
    // c
    p->c = p->H[0] + total; total += (m+1)*nvec;
    // s
    p->s = p->H[0] + total; total += (m+1)*nvec;

    // V
    for ( i=0; i<m+1; i++ ) {
      vector_PRECISION_init(&(p->V[i]));
      vector_PRECISION_alloc( &(p->V[i]), vl_type, nvec, l, no_threading );
    }
    // Z
    for ( i=0; i<k; i++ ) {
      vector_PRECISION_init(&(p->Z[i]));
      vector_PRECISION_alloc( &(p->Z[i]), vl_type, nvec, l, no_threading );
    }
    // x
    vector_PRECISION_alloc( &(p->x), vl_type, nvec, l, no_threading );
    // b
    vector_PRECISION_alloc( &(p->b), vl_type, nvec, l, no_threading );
    // w
    vector_PRECISION_alloc( &(p->w), vl_type, nvec, l, no_threading );
    // r
    vector_PRECISION_alloc( &(p->r), vl_type, nvec, l, no_threading );
    
    ASSERT( p->total_storage == total );
  } else if ( g.solver[l->depth] ) {
    // x
    vector_PRECISION_alloc( &(p->x), vl_type, nvec, l, no_threading );
    // b
    vector_PRECISION_alloc( &(p->b), vl_type, nvec, l, no_threading );
  }
  
}


void fgmres_PRECISION_struct_free( gmres_PRECISION_struct *p, level_struct *l ) {

/*********************************************************************************
* Frees the allocated space for the gmres struct p.                            
*********************************************************************************/ 
  
  int i, k=0;
  
  if ( p->preconditioner != NULL ) {
    if ( p->kind == _RIGHT )
      k = p->restart_length+1;
    else
      k = 1;
  } else {
#if defined(SINGLE_ALLREDUCE_ARNOLDI) && defined(PIPELINED_ARNOLDI)
    if ( l->level == 0 && l->depth > 0 ) {
      k = p->restart_length+2;
    }
#else
    k = p->restart_length;//0;//temp fix!!!!!!!
#endif
  }

  vector_PRECISION_free( &(p->x), l, no_threading );
  vector_PRECISION_free( &(p->b), l, no_threading );
  if((!g.solver[l->depth]||g.use_only_fgrmes_at_setup)&&p->restart_length > 0) {
    FREE( p->H[0], complex_PRECISION, p->total_storage );
    FREE( p->H, complex_PRECISION*, p->restart_length );
    for ( i=0; i<p->restart_length+1; i++ ) vector_PRECISION_free( &(p->V[i]), l, no_threading );
    FREE( p->V, vector_PRECISION, p->restart_length+1 );
    vector_PRECISION_free( &(p->w), l, no_threading );
    vector_PRECISION_free( &(p->r), l, no_threading );
    if ( p->Z != NULL ) {
      for ( i=0; i<k; i++ ) vector_PRECISION_free( &(p->Z[i]), l, no_threading );
      FREE( p->Z, vector_PRECISION, k );
    }
  }
  if ( p->fab.handle != NULL ) {
    fabulous_PRECISION_free( &(p->fab), l, no_threading );
  }
  
  p->D = NULL;
  p->clover = NULL;
}

// so far integrated interface to fgmres and fabulous
int solver_PRECISION( gmres_PRECISION_struct *p, level_struct *l, struct Thread *threading ) {

  int i, iter = 0, n_vect = num_loop;

  START_LOCKED_MASTER(threading)
  if ( l->depth==0 && ( p->timing || p->print ) ) prof_init( l );
  if ( l->level==0 && g.num_levels > 1 && g.interpolation ) p->tol = g.coarse_tol;
  if ( l->depth > 0 ) p->timing = 1;
  END_LOCKED_MASTER(threading)

  if ( p->fab.handle != NULL && !(g.use_only_fgrmes_at_setup && g.in_setup) ) {//need more generality if we are to use this as generic entry pt
    iter = fabulous_PRECISION( p, threading );
  } else {
    iter = fgmres_PRECISION( p, l, threading );
  }

  START_LOCKED_MASTER(threading)
  g.iter_counts[l->depth] += iter;
  END_LOCKED_MASTER(threading)
  if ( l->level == 0 ) {
    START_LOCKED_MASTER(threading)
    g.coarse_iter_count += iter;
    END_LOCKED_MASTER(threading)
  }

  // Report the statistics
  if ( p->print ) {
    START_LOCKED_MASTER(threading)
    printf0("+----------------------------------------------------------+\n");
    printf0("|              Final Relative Residuals                    |\n");
    for( i=0; i<n_vect; i++ )
      printf0("| exact relative residual, %d: ||r||/||b|| = %e   |\n",i, g.resids[i] );
    printf0("+----------------------------------------------------------+\n");
    printf0("|                   Solver Statistics                      |\n");
    printf0("|     %8s iterations: %-6d coarse average: %-6.2lf   |\n",
	    (g.solver[0]==0)?"FGMRES":"FABULOUS",
	    iter, ((double)g.iter_counts[g.num_levels-1])/((double)iter) ); //((double) iter)/((double) fab->nrhs)
    //printf0("|     elapsed wall clock time: %-8.4lf seconds                |\n", g.iter_times[0] );
    if ( l->depth == 0) for( i=0; i<g.num_levels; i++ ) {
        printf0("|  (depth %d) solver: %d, iterations: %-6d, time below this level: %-8.4lf sec (%04.1lf%%) |\n", i, g.solver[i], g.iter_counts[i],
		g.iter_times[i], 100*(g.iter_times[i]/g.iter_times[0]) );
    }
    printf0("|  consumed core minutes*: %-8.2le (solve only)           |\n", (g.iter_times[0]*g.num_processes*MAX(1,threading->n_core))/60.0 );
    printf0("|    max used mem/MPIproc: %-8.2le GB                     |\n", g.max_storage/1024.0 );
    printf0("+----------------------------------------------------------+\n");
    printf0("*: only correct if #MPIprocs*#threads == #CPUs\n\n");
    END_MASTER(threading)
  }
  if ( l->depth == 0 && g.vt.p_end != NULL  ) {
    if ( g.vt.p_end != NULL ) {
      START_LOCKED_MASTER(threading)
      printf0("solve iter: %d\n", iter );
      printf0("solve time: %le seconds\n", g.total_time );
      g.vt.p_end->values[_SLV_TIME] += g.total_time/((double)g.vt.average_over);
      g.vt.p_end->values[_SLV_ITER] += iter/((double)g.vt.average_over);
      g.vt.p_end->values[_CRS_ITER] += (((double)g.coarse_iter_count)/((double)iter))/((double)g.vt.average_over);
      g.vt.p_end->values[_CRS_TIME] += g.coarse_time/((double)g.vt.average_over);
      END_LOCKED_MASTER(threading)
    }
  }
  // Report profiling
  if ( l->depth == 0 && ( p->timing || p->print ) && !(g.vt.p_end != NULL )  ) {
    START_MASTER(threading)
    prof_print( l );
    END_MASTER(threading)
  }
  
  return iter;
}

int fgmres_PRECISION( gmres_PRECISION_struct *p, level_struct *l, struct Thread *threading ) {

/*********************************************************************************
* Uses FGMRES to solve the system D x = b, where b is taken from p->b and x is 
* stored in p->x.                                                              
*********************************************************************************/  


  // start and end indices for vector functions depending on thread
  int start;
  int end;

  SYNC_MASTER_TO_ALL(threading)
  SYNC_CORES(threading)

  int j=-1, finish=0, iter=0, il, ol, res;// related to the above change
  int n_vect=num_loop, i, jj;//!!!!!!!!!!!g.num_vect_now
  PRECISION beta[n_vect];

  double H_tot;
  PRECISION norm_r0[n_vect], gamma_jp1[n_vect], gamma_max, gamma0_real[n_vect], t0=0, t1=0;

  if ( p->num_vect < num_loop )//g.num_vect_now ) 
    error0("fgmres_PRECISION: memory corruption\n");

  p->w.num_vect_now = n_vect; p->x.num_vect_now = n_vect; p->r.num_vect_now = n_vect; p->b.num_vect_now = n_vect;
  p->Z[0].num_vect_now = n_vect; p->V[0].num_vect_now = n_vect;

  VECTOR_LOOP(i, n_vect, jj, norm_r0[i+jj]=1; 
                             gamma_jp1[i+jj]=1;)
  
  START_LOCKED_MASTER(threading)
  if ( l->depth == 0 ) t0 = MPI_Wtime();
#if defined(TRACK_RES) && !defined(WILSON_BENCHMARK)
  if ( p->print && g.print > 0 ) printf0("+----------------------------------------------------------+\n");
#endif
  END_LOCKED_MASTER(threading)

  // compute start and end indices for core
  // this puts zero for all other hyperthreads, so we can call functions below with all hyperthreads
  compute_core_start_end_custom(p->v_start, p->v_end, &start, &end, l, threading,l->num_lattice_site_var);

  // outer loop over restartings: if go beyond max #restarts or has converged, exit
  for( ol=0; ol<p->num_restart && finish==0; ol++ ) {
    
    //------- iniital setup: Compute r_0 = b − Dx_0 , β := ||r_0||_2 , and v_1 := r_0 /β
    if( ol == 0 && p->initial_guess_zero ) {
      res = _NO_RES;
      vector_PRECISION_copy( &(p->r), &(p->b), start, end, l ); // p->r <- p->b
    } else { // if the initial guess is not zero, need to compute b_l-D_l phi_i 
      res = _RES;
      if ( p->kind == _LEFT && p->preconditioner ) {// what is LEFT?????
        apply_operator_PRECISION( &(p->Z[0]), &(p->x), p, l, threading ); // Z[0] <- D*x
        if ( g.method == 5 ) {
          START_LOCKED_MASTER(threading)
	  g.bicgstab_tol = (!g.mixed_precision)?p->tol:MAX( 1E-3, (p->tol/find_max_PRECISION(n_vect,gamma_jp1,norm_r0))*5E-1 );
	  END_LOCKED_MASTER(threading)
	}
        p->preconditioner( &(p->w), NULL, &(p->Z[0]), _NO_RES, l, threading );
      } else {
        apply_operator_PRECISION( &(p->w), &(p->x), p, l, threading ); // compute w = D*x
      }
      vector_PRECISION_minus( &(p->r), &(p->b), &(p->w), start, end, l ); // compute r = b - w = b - D*x
    }
    global_norm_PRECISION( gamma0_real, &(p->r), p->v_start, p->v_end, l, threading ); // gamma_0 = norm(r)  
    START_MASTER(threading)
    VECTOR_LOOP(i, n_vect, jj, p->gamma[i+jj] = (complex_PRECISION) gamma0_real[i+jj];)
    END_MASTER(threading);
    SYNC_MASTER_TO_ALL(threading);
    
    // report
    if ( ol == 0 ) {
      if (l->depth == 0 && !p->initial_guess_zero) {
	global_norm_PRECISION( norm_r0, &(p->b), p->v_start, p->v_end, l, threading );
	for( i=0; i<n_vect; i++ )
	  printf0("| initial guess relative residual (%d):        %le |\n", i, gamma0_real[i]/norm_r0[i]);
      } else {
	VECTOR_LOOP(i, n_vect, jj, norm_r0[i+jj] = gamma0_real[i+jj];)
	}
    }
    vector_PRECISION_real_scale( &(p->V[0]), &(p->r), p->gamma, 0, 1, start, end, l ); // v_0 = r / gamma_0

    //----- main body: Expand the Krylov space possibly with preconditioner until one of the stop conditions is met
#if defined(SINGLE_ALLREDUCE_ARNOLDI) && defined(PIPELINED_ARNOLDI)
    if ( l->level == 0 && l->depth > 0 ) {
      arnoldi_step_PRECISION( p->V, p->Z, &(p->w), p->H, p->y, 0, p->preconditioner, p, l, threading );
    }
#endif   
    
    // inner loop over Arnoldi steps: if (dim of the Krylov subspace)>p->restart_length or has converged, exit this loop
    for( il=0; il<p->restart_length && finish==0; il++) { 
      j = il; iter++;
      if ( g.method == 6 ) {
        START_LOCKED_MASTER(threading)
	g.bicgstab_tol = (!g.mixed_precision)?p->tol:MAX( 1E-3, (p->tol/find_max_PRECISION(n_vect,gamma_jp1,norm_r0))*5E-1 );
        END_LOCKED_MASTER(threading)
      }
      
      // one step of Arnoldi: extend the Krylov basis and define the (m + 1) x m Hessenberg matrix
#if defined(SINGLE_ALLREDUCE_ARNOLDI) && defined(PIPELINED_ARNOLDI)
      if ( l->level == 0 && l->depth > 0 ) {
        if ( !arnoldi_step_PRECISION( p->V, p->Z, &(p->w), p->H, p->y, j+1, p->preconditioner, p, l, threading ) ) {
          printf0("| -------------- iteration %d, restart due to H(%d,%d) < 0 |\n", iter, j+2, j+1 );
          break;
        }
      } else {
        if ( !arnoldi_step_PRECISION( p->V, p->Z, &(p->w), p->H, p->y, j, p->preconditioner, p, l, threading ) ) {
          printf0("| -------------- iteration %d, restart due to H(%d,%d) < 0 |\n", iter, j+1, j );
          break;
        }
      }
#else
      if ( !arnoldi_step_PRECISION( p->V, p->Z, &(p->w), p->H, p->y, j, p->preconditioner, p, l, threading ) ) {
        printf0("| -------------- iteration %d, restart due to H(%d,%d) < 0 |\n", iter, j+1, j );
        break;//exit: for( il=0; il<p->restart_length && finish==0; il++)
      }
#endif

      H_tot=0;
      VECTOR_LOOP(i, n_vect, jj, H_tot += (double) cabs_PRECISION( p->H[j][(j+1)*n_vect+i+jj] );)//what is this for????
      if ( H_tot > n_vect*p->tol/10 ) {
        qr_update_PRECISION( p->H, p->s, p->c, p->gamma, j, l, threading );
        VECTOR_LOOP(i, n_vect, jj, gamma_jp1[i+jj] = cabs_PRECISION( p->gamma[(j+1)*n_vect+i+jj] );)

#if defined(TRACK_RES) && !defined(WILSON_BENCHMARK)
        if ( iter%10 == 0 || p->preconditioner != NULL || l->depth > 0 ) {
          START_MASTER(threading)
	  if ( p->print && g.print > 0 )//printed only at the top
            for( i=0; i<n_vect; i++ )
	      printf0("| vector %d, depth: %d, approx. rel. res. after  %-6d iterations: %e |\n", i, l->depth, iter, gamma_jp1[i]/norm_r0[i]);
          END_MASTER(threading)
        }
#endif
	// The last entry of gamma at (j+1)^th step is the minimized residual of the solution to the projected lin. eq. 
        gamma_max = gamma_jp1[0]/norm_r0[0];
	for ( i=1; i<n_vect; i++ ) if ( gamma_max < gamma_jp1[i]/norm_r0[i] ) gamma_max = gamma_jp1[i]/norm_r0[i];
	if( gamma_max < p->tol || gamma_max > 1E+5 ) { // if satisfied ... stop
	  finish = 1;
	  START_MASTER(threading)
	  if ( gamma_max > 1E+5 ) printf0("Divergence of fgmres_PRECISION, iter = %d, level=%d\n", iter, l->level );
          END_MASTER(threading)
	}
      } else {
        for( i=0; i<n_vect; i++ )
          printf0("vector: %d, depth: %d, iter: %d, p->H(%d,%d) = %+lf+%lfi\n", i, l->depth, iter, j+1, j, CSPLIT( p->H[j][(j+1)*n_vect+i] ) );
        finish = 1;
        break;
      }
    } // end of a single restart
    // Compute the minimizer y of the residual on the Krylov subspace and update the iterate: reach this pt. if converged or go beyond restart length
    compute_solution_PRECISION( &(p->x), (p->preconditioner&&p->kind==_RIGHT)?(p->Z):(p->V),
				p->y, p->gamma, p->H, j, (res==_NO_RES)?ol:1, p, l, threading );
  } // end of fgmres:for( ol=0; ol<p->num_restart && finish==0; ol++ ) 
  
  //-----------  Compute the statistics
  START_LOCKED_MASTER(threading)
  if ( l->depth == 0 ) { t1 = MPI_Wtime(); g.iter_times[0] = t1-t0; g.norm_res = gamma_max ; }//g.iter_counts[0] = iter; g.norm_res = gamma_max ; }
  END_LOCKED_MASTER(threading)
  if ( p->print ) {
#ifdef FGMRES_RESTEST
    apply_operator_PRECISION( &(p->w), &(p->x), p, l, threading );
    vector_PRECISION_minus( &(p->r), &(p->b), &(p->w), start, end, l );
    global_norm_PRECISION( beta, &(p->r), p->v_start, p->v_end, l, threading );
#else
    VECTOR_LOOP(i, n_vect, jj, beta[i+jj] = creal_PRECISION(gamma_jp1[i+jj]);)
#endif
    START_MASTER(threading)
    g.norm_res = 0;
    VECTOR_LOOP(i, n_vect, jj, g.norm_res += beta[i+jj]/norm_r0[i+jj];)
    VECTOR_LOOP(i, n_vect, jj, g.resids[i+jj] = beta[i]/norm_r0[i] );
#if defined(TRACK_RES) && !defined(WILSON_BENCHMARK)
    if ( g.print > 0 ) printf0("+----------------------------------------------------------+\n\n");
#endif
    END_MASTER(threading)
  }
  // This is FGMRES specific
#ifdef COARSE_RES
  if ( l->depth > 0 ) {
    START_MASTER(threading)
    char number[3]; sprintf( number, "%2d", 31+l->depth ); printf0("\033[1;%2sm|", number );
    for (i=0; i<num_loop; i++ ) printf0(" - depth: %d, gmres iter: %2d, approx rel res: %le |", l->depth, iter, gamma_jp1[i]/norm_r0[i] );
    printf0("\033[0m\n"); fflush(0);
    END_MASTER(threading)
  }
#endif

  return iter;
}

//#define FAB_OPENMP
int fabulous_PRECISION( gmres_PRECISION_struct *p, struct Thread *threading ) {

  fabulous_PRECISION_struct *fab = &(p->fab);
  level_struct * l = fab->l;
  int iter = 0, even_size = p->op->num_even_sites*l->num_lattice_site_var;
  PRECISION t0, t1;
  //printf0("begin PRECISION fab: solver %d depth %d: %d=%d? %d, sizes %d %d %d even sts %d\n", g.solver, l->depth, fab->nrhs,fab->B.num_vect,p->b.num_vect_now, p->b.size, fab->B.size, fab->dim, even_size);fflush(stdout);
    
  p->b.num_vect_now = num_loop; p->x.num_vect_now = num_loop;

  START_LOCKED_MASTER(threading)
  if ( l->depth == 0 ) t0 = MPI_Wtime();
  END_LOCKED_MASTER(threading)

#ifdef FAB_OPENMP
  fab->threading = threading; // if fabulous is OpenMP safe.  This can be uncommented.
#endif
  vector_PRECISION_copy( &(fab->B), &(p->b), 0, fab->dim, l );
  if( p->initial_guess_zero ) {
    vector_PRECISION_define( &(fab->X), 0, 0, fab->dim, l );
  } else {
    vector_PRECISION_copy( &(fab->X), &(p->x), 0, fab->dim, l );
  }
  START_LOCKED_MASTER(threading)
  vector_PRECISION_change_layout( &(fab->B), &(fab->B), _NVEC_OUTER, no_threading );
  vector_PRECISION_change_layout( &(fab->X), &(fab->X), _NVEC_OUTER, no_threading );
#ifdef FAB_OPENMP
  END_LOCKED_MASTER(threading)
#endif
  void *X = fab->X.vector_buffer, *B = fab->B.vector_buffer;
  switch ( g.solver[l->depth] ) {
  case _GCR:    iter = fabulous_solve_GCR( fab->nrhs, B, fab->ldb, X, fab->ldx, fab->handle ); break;
  case _IB:     iter = fabulous_solve_IB( fab->nrhs, B, fab->ldb, X, fab->ldx, fab->handle ); break;
  case _DR:     iter = fabulous_solve_DR( fab->nrhs, B, fab->ldb, X, fab->ldx, fab->k, fab->eigvals, fab->handle ); break;
  case _IBDR:   iter = fabulous_solve_IBDR( fab->nrhs, B, fab->ldb, X, fab->ldx, fab->k, fab->eigvals, fab->handle ); break;
  case _QR:     iter = fabulous_solve_QR( fab->nrhs, B, fab->ldb, X, fab->ldx, fab->handle ); break;
  case _QRIB:   iter = fabulous_solve_QRIB( fab->nrhs, B, fab->ldb, X, fab->ldx, fab->handle ); break;
  case _QRDR:   iter = fabulous_solve_QRDR( fab->nrhs, B, fab->ldb, X, fab->ldx, fab->k, fab->eigvals, fab->handle ); break;
  case _QRIBDR: iter = fabulous_solve_QRIBDR( fab->nrhs, B, fab->ldb, X, fab->ldx, fab->k, fab->eigvals, fab->handle ); break;
  default:      iter = fabulous_solve( fab->nrhs, B, fab->ldb, X, fab->ldx, fab->handle ); break;
  }
#ifdef FAB_OPENMP
  START_LOCKED_MASTER(threading)
#endif
  fab->B.layout = _NVEC_INNER;
  vector_PRECISION_change_layout( &(fab->X), &(fab->X), _NVEC_INNER, no_threading );
  END_LOCKED_MASTER(threading)
  vector_PRECISION_copy( &(p->x), &(fab->X), 0, fab->dim, l );

  // Compute the statitics
  START_LOCKED_MASTER(threading)
    if ( l->depth == 0 ) { t1 = MPI_Wtime(); g.iter_times[0] = t1-t0;}// g.iter_counts[0] = iter; }
  END_LOCKED_MASTER(threading)
  if ( p->print > 0 ) {
    int start, end, j, jj;
    PRECISION norm[num_loop], norm2[num_loop];
    compute_core_start_end_custom(p->v_start, p->v_end, &start, &end, l, threading,l->num_lattice_site_var);
    apply_operator_PRECISION( &(fab->B0), &(p->x), p, l, threading ); // compute w = D*x
    vector_PRECISION_minus( &(fab->B0), &(p->b), &(fab->B0), 0, p->v_end, l ); // compute r = b - w = b - D*x
    global_norm_PRECISION( norm, &(fab->B0), p->v_start, p->v_end, l, threading ); // gamma_0 = norm(r)
    global_norm_PRECISION( norm2, &(p->b), p->v_start, p->v_end, l, threading );
    VECTOR_LOOP(j, num_loop, jj, g.resids[j+jj] = norm[j+jj]/norm2[j+jj];)
  }
      
  return iter;
}

void bicgstab_PRECISION( gmres_PRECISION_struct *ps, level_struct *l, struct Thread *threading ) {
  /*********************************************************************************
   * Uses BiCGstab to solve the system D x = b, where b is taken from ps->b and x is 
   * stored in ps->x.                                                              
   *********************************************************************************/
  warning0("bicgstab_PRECISION: this function has not been tested\n");

  int iter=0, maxiter, j, jj, nvec = num_loop;
  vector_PRECISION x, b, r, r_tilde, p, pp, v, s, t; // Krylov subspace size: 5
  complex_PRECISION alpha[nvec], beta[nvec], rho[nvec], rho_old[nvec], omega[nvec], norm1[nvec], norm2[nvec];
  PRECISION tol, b_norm[nvec], r_norm[nvec], s_norm[nvec];
  // start and end indices for vector functions depending on thread
  int start;
  int end;
  // compute start and end indices for core
  // this puts zero for all other hyperthreads, so we can call functions below with all hyperthreads
  compute_core_start_end_custom(ps->v_start, ps->v_end, &start, &end, l, threading, l->num_lattice_site_var );
  
  tol = (l->level==0 && g.num_levels > 1 && g.interpolation )?g.coarse_tol:g.bicgstab_tol;
  maxiter = 1000000; 
  r = ps->r; b = ps->b; x = ps->x; p = ps->w;
  pp = ps->V[0]; r_tilde = ps->V[1]; v = ps->V[2]; s = ps->V[3]; t = ps->V[4];

  VECTOR_LOOP(j, nvec, jj, alpha[j+jj] = 1; beta[j+jj] = 1; );
  vector_PRECISION_copy( &r, &b, start, end, l );
  vector_PRECISION_copy( &r_tilde, &b, start, end, l );
  vector_PRECISION_define( &x, 0, start, end, l );
  vector_PRECISION_define( &v, 0, start, end, l );
  vector_PRECISION_define( &s, 0, start, end, l );
  vector_PRECISION_define( &t, 0, start, end, l );
  global_norm_PRECISION( b_norm, &b, ps->v_start, ps->v_end, l, threading );

  VECTOR_LOOP(j, nvec, jj, r_norm[j+jj] = b_norm[j+jj];)
#if defined(TRACK_RES) && !defined(WILSON_BENCHMARK)  
  START_MASTER(threading)
  printf0("+----------------------------------------------------------+\n");
  END_MASTER(threading)
#endif
  while ( find_max_PRECISION( nvec, r_norm, b_norm ) > tol && iter < maxiter ) {
    iter++;
    
    VECTOR_LOOP(j, nvec, jj, rho_old[j+jj] = rho[j+jj]; )
    global_inner_product_PRECISION( rho, &r_tilde, &r, ps->v_start, ps->v_end, l, threading );

    if ( rho == 0 ) {
      START_MASTER(threading)
      printf0("rho = 0: BiCGstab did not converge.\n");
      END_MASTER(threading)
      break;
    }
    
    if ( iter == 1 ) {
      vector_PRECISION_copy( &p, &r, start, end, l );
    } else {
      VECTOR_LOOP(j, nvec, jj, beta[j+jj] = (rho[j+jj]/rho_old[j+jj])*(alpha[j+jj]/omega[j+jj]);)
	vector_PRECISION_saxpy( &pp, &p,  &v, omega, 0, -1, start, end, l );
      vector_PRECISION_saxpy( &p,  &r, &pp  , beta, 0, 1, start, end, l );
    }    
    apply_operator_PRECISION( &v, &p, ps, l, threading );
    global_inner_product_PRECISION( norm1, &r_tilde, &v, ps->v_start, ps->v_end, l, threading );
    VECTOR_LOOP(j, nvec, jj, alpha[j+jj] = rho[j+jj] / norm1[j+jj];)
    vector_PRECISION_saxpy( &s, &r, &v, alpha, 0, -1, start, end, l );
    global_norm_PRECISION( s_norm, &s, ps->v_start, ps->v_end, l, threading );

    if ( find_max_PRECISION( nvec, s_norm,b_norm) < tol ) {
      vector_PRECISION_saxpy( &x, &x, &p, alpha, 0, 1, start, end, l );
      break;
    }
    
    apply_operator_PRECISION( &t, &s, ps, l, threading );
    global_inner_product_PRECISION( norm1, &t, &s, ps->v_start, ps->v_end, l, threading );
    global_inner_product_PRECISION( norm2, &t, &t, ps->v_start, ps->v_end, l, threading );
    VECTOR_LOOP(j, nvec, jj, omega[j+jj] = norm1[j+jj]/norm2[j+jj];)
    
    vector_PRECISION_saxpy( &x, &x, &p,  alpha, 0, 1, start, end, l );
    vector_PRECISION_saxpy( &x, &x, &s,  omega, 0, 1, start, end, l );
    vector_PRECISION_saxpy( &r, &s, &t,  omega, 0, -1, start, end, l );

    global_norm_PRECISION( r_norm, &r, ps->v_start, ps->v_end, l, threading );

#if defined(TRACK_RES) && !defined(WILSON_BENCHMARK)
    START_MASTER(threading)
    if ( iter % 100 == 0 ) for ( int i=0; i<nvec; i++ ) printf0("| vector %d, biCGstab relres: %12.6le,  iterations: %-8d     |\n", i, r_norm[i]/b_norm[i], iter );
    END_MASTER(threading)
#endif
  }
#if defined(TRACK_RES) && !defined(WILSON_BENCHMARK)
  START_MASTER(threading)
  for ( int i=0; i<nvec; i++ ) printf0("| vector %d, biCGstab relres: %12.6le,  iterations: %-8d     |\n", i, r_norm[i]/b_norm[i], iter );
  printf0("+----------------------------------------------------------+\n");
  END_MASTER(threading)
#endif
}


void cgn_PRECISION( gmres_PRECISION_struct *ps, level_struct *l, struct Thread *threading ) {
  /*********************************************************************************
   * Uses CGN to solve the system D x = b, where b is taken from ps->b and x is 
   * stored in ps->x.                                                              
   *********************************************************************************/
  warning0("cgn_PRECISION: this function has not been tested\n");

  int maxiter, iter=0, nvec = num_loop, j, jj;
  vector_PRECISION r_old, r_new, r_true, p, pp, Dp, x, b;
  complex_PRECISION alpha[nvec], beta[nvec], gamma[nvec];
  PRECISION tol, norm[nvec], r0_norm[nvec], r_norm[nvec], prod_rr_old[nvec], t0=0, t1=0;
  // start and end indices for vector functions depending on thread
  int start;
  int end;
  
  b = ps->b; x = ps->x;
  r_old = ps->V[2]; r_new = ps->V[3]; r_true = ps->r;
  p = ps->w; pp = ps->V[0]; Dp = ps->V[1];
  tol = (l->level==0 && g.num_levels > 1 && g.interpolation )?g.coarse_tol:ps->tol;
  maxiter = ps->num_restart;
  
  START_MASTER(threading)
  if ( ps->timing || ps->print ) t0 = MPI_Wtime();
  END_MASTER(threading)

  // compute start and end indices for core
  // this puts zero for all other hyperthreads, so we can call functions below with all hyperthreads
  compute_core_start_end_custom(ps->v_start, ps->v_end, &start, &end, l, threading, l->num_lattice_site_var);

  vector_PRECISION_define( &x, 0, start, end, l );
  apply_operator_PRECISION( &Dp, &x, ps, l, threading );
  vector_PRECISION_minus( &pp, &b, &Dp, start, end, l );
  apply_operator_dagger_PRECISION( &r_old, &pp, ps, l, threading );
  
  vector_PRECISION_copy( &p, &r_old, start, end, l );
  global_norm_PRECISION( r0_norm, &r_old, ps->v_start, ps->v_end, l, threading );
  VECTOR_LOOP(j, nvec, jj, prod_rr_old[j+jj] = r0_norm[j+jj]*r0_norm[j+jj];)
#if defined(TRACK_RES) && !defined(WILSON_BENCHMARK)
  if ( ps->print ) {
    START_MASTER(threading)
    printf0("\n+----------------------------------------------------------+\n");
    END_MASTER(threading)
  }  
#endif
  while ( find_max_PRECISION( nvec, norm, r0_norm ) > tol && iter < maxiter ) {
    iter++;
    
    apply_operator_PRECISION( &pp, &p, ps, l, threading );
    apply_operator_dagger_PRECISION( &Dp, &pp, ps, l, threading );
    
    global_inner_product_PRECISION( gamma, &p, &Dp, ps->v_start, ps->v_end, l, threading );
    VECTOR_LOOP(j, nvec, jj, alpha[j+jj] = prod_rr_old[j+jj] / gamma[j+jj];)
    vector_PRECISION_saxpy( &x, &x, &p, alpha, 0, 1, start, end, l );
    vector_PRECISION_saxpy( &r_new, &r_old, &Dp, alpha, 0, -1, start, end, l );
    
    global_inner_product_PRECISION( gamma, &r_new, &r_new, ps->v_start, ps->v_end, l, threading );
    VECTOR_LOOP(j, nvec, jj, beta[j+jj] = gamma[j+jj] / prod_rr_old[j+jj];)
    
    vector_PRECISION_saxpy( &p, &r_new, &p, beta, 0, 1, start, end, l );
    vector_PRECISION_copy( &r_old, &r_new, start, end, l );
    VECTOR_LOOP(j, nvec, jj, prod_rr_old[j+jj] = gamma[j+jj];)
    VECTOR_LOOP(j, nvec, jj, norm[j+jj] = sqrt(prod_rr_old[j+jj]);)
#if defined(TRACK_RES) && !defined(WILSON_BENCHMARK)    
    if ( iter%100 == 0 && ps->print >=1 ) {
      START_MASTER(threading)
      for( int i=0; i<nvec; i++ ) printf0("|      vector %d, NE rel. res. after  %-6d iterations: %e |\n", i, iter, sqrt(prod_rr_old[i])/r0_norm[i] );
      END_MASTER(threading)
    }
#endif
  }
  
  global_norm_PRECISION( r0_norm, &b, ps->v_start, ps->v_end, l, threading );
  apply_operator_PRECISION( &Dp, &x, ps, l, threading );
  vector_PRECISION_minus( &r_true, &b, &Dp, start, end, l );
  global_norm_PRECISION( r_norm, &r_true, ps->v_start, ps->v_end, l, threading );

#if defined(TRACK_RES) && !defined(WILSON_BENCHMARK)  
  if ( ps->print ) {
    START_MASTER(threading)
    printf0("+----------------------------------------------------------+\n");
    for( int i=0; i<nvec; i++ ) printf0("| vector %d, switching to CGNR, iter  %-6d true r res: %e |\n", i, iter, r_norm[i]/r0_norm[i] );
    printf0("+----------------------------------------------------------+\n");
    END_MASTER(threading)
  }
#endif
  
  while ( find_max_PRECISION(nvec, r_norm, r0_norm) > tol && iter < maxiter ) {
    iter++;
    
    apply_operator_PRECISION( &pp, &p, ps, l, threading );
    apply_operator_dagger_PRECISION( &Dp, &pp, ps, l, threading );
    
    global_inner_product_PRECISION( gamma, &p, &Dp, ps->v_start, ps->v_end, l, threading );
    VECTOR_LOOP(j, nvec, jj, alpha[j+jj] = prod_rr_old[j+jj] / gamma[j+jj];)
    vector_PRECISION_saxpy( &x, &x, &p, alpha, 0, 1, start, end, l );
    vector_PRECISION_saxpy( &r_new, &r_old, &Dp, alpha, 0, -1, start, end, l );
    
    // residual update
    vector_PRECISION_saxpy( &r_true, &r_true, &pp, alpha, 0, -1, start, end, l );
    global_norm_PRECISION( r_norm, &r_true, ps->v_start, ps->v_end, l, threading );
    global_inner_product_PRECISION( gamma, &r_new, &r_new, ps->v_start, ps->v_end, l, threading );
    VECTOR_LOOP(j, nvec, jj, beta[j+jj] = gamma[j+jj] / prod_rr_old[j+jj];)
    
    vector_PRECISION_saxpy( &p, &r_new, &p, beta, 0, 1, start, end, l );
    vector_PRECISION_copy( &r_old, &r_new, start, end, l );
    VECTOR_LOOP(j, nvec, jj, prod_rr_old[j+jj] = gamma[j+jj];)
#if defined(TRACK_RES) && !defined(WILSON_BENCHMARK)    
    if ( iter%100 ==  0 && ps->print >=1 ) {
      START_MASTER(threading)
      for ( int i=0; i<nvec; i++ ) printf0("|         vector %d, rel. res. after  %-6d iterations: %e |\n", i, iter, r_norm[i]/r0_norm[i] );
      END_MASTER(threading)
    }
#endif
  }
  
  if ( ps->timing || ps->print ) t1 = MPI_Wtime();
  if ( ps->print ) {
    START_MASTER(threading)
    printf0("+----------------------------------------------------------+\n");
    printf0("|          CGN iterations: %-6d                          |\n", iter );
    END_MASTER(threading)
    apply_operator_PRECISION( &Dp, &x, ps, l, threading );
    vector_PRECISION_minus( &pp, &b, &Dp, start, end, l );

    global_norm_PRECISION( norm, &pp, ps->v_start, ps->v_end, l, threading );
    VECTOR_LOOP(j, nvec, jj, beta[j+jj] = norm[j+jj];)
    START_MASTER(threading)
    if ( ps->timing ) for( int i=0; i<nvec; i++ ) printf0("| vector %d, exact relative residual: ||r||/||b|| = %e      |\n", i, creal(beta[i]/r0_norm[i]) );
    printf0("| elapsed wall clock time: %-12g seconds            |\n", t1-t0 );
    printf0("|  consumed core minutes*: %-8.2le (solve only)           |\n", ((t1-t0)*g.num_processes*MAX(1,threading->n_core))/60.0 );
    printf0("|    max used mem/MPIproc: %-8.2le GB                     |\n", g.max_storage/1024.0 );
    printf0("+----------------------------------------------------------+\n");
    printf0("*: only correct if #MPIprocs*#threads == #CPUs\n\n");
    END_MASTER(threading)
  }
  
  START_LOCKED_MASTER(threading)
  if ( l->level == 0 )
    g.coarse_iter_count += iter;

  if ( l->depth == 0 && g.vt.p_end != NULL  ) {
    if ( g.vt.p_end != NULL ) {
      g.vt.p_end->values[_SLV_TIME] += (t1-t0)/((double)g.vt.average_over);
      g.vt.p_end->values[_SLV_ITER] += (iter)/((double)g.vt.average_over);
    }
  }
  END_LOCKED_MASTER(threading)
}

//   * NOTE: not used anywhere!!!!
void fgcr_PRECISION( gmres_PRECISION_struct *p, level_struct *l ) { 
  /*********************************************************************************
   * Uses FGCR to solve the system D x = b, where b is taken from p->b and x is 
   * stored in p->x.                                                              
   * Assume: should be called only by master
   *********************************************************************************/
  warning0("fgcr_PRECISION: this function has not been tested\n");
  int i, j=-1, finish=0, iter=0, il, ol, nvec = num_loop, jj, jjj;
  complex_PRECISION beta[nvec], alpha[nvec];
  PRECISION norm[nvec], r0_norm[nvec], t0=0, t1=0;
  
  if ( p->timing || p->print ) t0 = MPI_Wtime();
#if defined(TRACK_RES) && !defined(WILSON_BENCHMARK)  
  if ( p->print ) printf0("+----------------------------------------------------------+\n");
#endif
  for( ol=0; ol<p->num_restart && finish==0; ol++ )  {
  
    if( ol == 0 && p->initial_guess_zero ) {
      vector_PRECISION_copy( &(p->r), &(p->b), p->v_start, p->v_end, l );

    } else {
      apply_operator_PRECISION( &(p->w), &(p->x), p, l, no_threading ); // compute w = D*x
      vector_PRECISION_minus( &(p->r), &(p->b), &(p->w), p->v_start, p->v_end, l ); // compute r = b - w
    }
    
    if( ol == 0) {
      global_norm_PRECISION( r0_norm, &(p->r), p->v_start, p->v_end, l, no_threading );
    }
    
    for( il=0; il<p->restart_length && finish==0; il++ ) {
      
      j = il; iter++;
      
      p->preconditioner( &(p->V[j]), &(p->r), _NO_RES, l, no_threading );
      apply_operator_PRECISION( &(p->Z[j]), &(p->V[j]), p, l, no_threading );
      
      for( i=0; i<j; i++ ) {
        global_inner_product_PRECISION( beta, &(p->Z[i]), &(p->Z[j]), p->v_start, p->v_end, l, no_threading );
	VECTOR_LOOP(jj, nvec, jjj, beta[jj+jjj] /= p->gamma[j*nvec+jj+jjj];)
        vector_PRECISION_saxpy( &(p->V[j]), &(p->V[j]), &(p->V[i]), beta, 0, -1, p->v_start, p->v_end, l );
        vector_PRECISION_saxpy( &(p->Z[j]), &(p->Z[j]), &(p->Z[i]), beta, 0, -1, p->v_start, p->v_end, l );
      }
      
      global_inner_product_PRECISION( p->gamma+j*nvec, &(p->Z[j]), &(p->Z[j]), p->v_start, p->v_end, l, no_threading );
      global_inner_product_PRECISION( alpha, &(p->Z[j]), &(p->r), p->v_start, p->v_end, l, no_threading );
      VECTOR_LOOP(jj, nvec, jjj, alpha[jj+jjj] /= p->gamma[j*nvec+jj+jjj];)
      vector_PRECISION_saxpy( &(p->x), &(p->x), &(p->V[j]), alpha, 0, 1, p->v_start, p->v_end, l );
      vector_PRECISION_saxpy( &(p->r), &(p->r), &(p->Z[j]), alpha, 0, -1, p->v_start, p->v_end, l );
      
      global_norm_PRECISION( norm, &(p->r), p->v_start, p->v_end, l, no_threading );
      VECTOR_LOOP(jj, nvec, jjj, alpha[jj+jjj] = norm[jj+jjj]/r0_norm[jj+jjj]; norm[jj+jjj] = creal(alpha[jj+jjj]);)
	if ( find_max_PRECISION(nvec, norm, NULL) < p->tol ) {
        finish = 1;
        break;
      } else {
#if defined(TRACK_RES) && !defined(WILSON_BENCHMARK)
        if ( iter%10 == 0 && p->print  ) for( int i=0; i<nvec; i++ ) printf0("| vector %d, approx. rel. res. after  %-6d iterations: %e |\n", i, iter, alpha[i] );
#endif
      }
    } // end of restart
  } // end of fgcr
  
  if ( p->timing || p->print ) t1 = MPI_Wtime();
  if ( p->print ) {
    apply_operator_PRECISION( &(p->w), &(p->x), p, l, no_threading );
    vector_PRECISION_minus( &(p->r), &(p->b), &(p->w), p->v_start, p->v_end, l );
    global_norm_PRECISION( norm, &(p->r), p->v_start, p->v_end, l, no_threading );
#if defined(TRACK_RES) && !defined(WILSON_BENCHMARK)
    printf0("+----------------------------------------------------------+\n");
    printf0("\n");
#endif
    printf0("+----------------------------------------------------------+\n");
    printf0("|         FGCR iterations: %-6d                          |\n", iter );
    for( int i=0; i<nvec; i++ ) printf0("| vector %d, exact relative residual: ||r||/||b|| = %e      |\n", i, norm[i]/r0_norm[i] );
    printf0("| elapsed wall clock time: %-7lf seconds                |\n", t1-t0 );
    if ( g.coarse_time > 0 ) 
      printf0("|        coarse grid time: %-7lf seconds (%04.1lf%%)        |\n",
              g.coarse_time, 100*(g.coarse_time/(t1-t0)) );
    printf0("+----------------------------------------------------------+\n\n");
  }
}

int arnoldi_step_PRECISION( vector_PRECISION *V, vector_PRECISION *Z, vector_PRECISION *w,
                            complex_PRECISION **H, complex_PRECISION* buffer, int j, void (*prec)(),
                            gmres_PRECISION_struct *p, level_struct *l, struct Thread *threading ) {

  /*********************************************************************************
   * Extends the Arnoldi basis by one vector (with application of prec if available),
   * while computing Hessenberg matrix on the extended basis
   * - vector_PRECISION **V: Contains the Arnoldi basis vectors.
   * - vector_PRECISION **Z: If a right precond. P is used, contains P*V[j] for all j.
   * - vector_PRECISION *w: Will be appended to existing Arnoldi basis at 
   *   position j+1.
   * - complex_PRECISION **H: Contains full Hessenberg matrix from the Arnoldi 
   *   decomposition (columnmajor!)
   * - complex_PRECISION* buffer: Buffer for local inner products.
   * - int j: index of the new Arnoldi vector to be orthonormalized
   *   against all previous ones.
   * - void (*prec)(): Function pointer to preconditioner (can be NULL if no 
   *   preconditioning is used).
   *********************************************************************************/
  
  int start, end, i, jj, jjj, nvec = w->num_vect_now;//!!!!!!!!!!!
#ifdef SINGLE_ALLREDUCE_ARNOLDI
  complex_PRECISION sigma[nvec], sigma_sum = 0;//originally decleared as const!!!!!!!!what does this do?????
  VECTOR_LOOP(jj, nvec, jjj, sigma[jj+jjj]=0;)
#ifdef PIPELINED_ARNOLDI
  warning0("arnoldi_step_PRECISION: with single step Arnoldi and Pipelined, this has not been tested\n");
  if ( l->level == 0 && l->depth > 0 ) {
    SYNC_MASTER_TO_ALL(threading)
    MPI_Request req;
    MPI_Status stat;
    compute_core_start_end_cutom(p->v_start, p->v_end, &start, &end, l, threading, l->num_lattice_site_var);

    V[j].num_vect_now = nvec;  V[j+1].num_vect_now = nvec; Z[j].num_vect_now = nvec; Z[j+1].num_vect_now = nvec;
    if ( j == 0 )
      vector_PRECISION_copy( &Z[0], &V[0], start, end, l );
    else
      vector_PRECISION_copy( &V[j], &Z[j], start, end, l );

    complex_PRECISION tmp[(j+1)*nvec];
    process_multi_inner_product_PRECISION( j+1, tmp, V, &V[j], p->v_start, p->v_end, l, threading );
    START_MASTER(threading)
    PROF_PRECISION_START( _ALLR );
    for( i=0; i<=j; i++ )
      VECTOR_LOOP(jj, nvec, jjj, buffer[i*nvec+jj+jjj] = tmp[i*nvec+jj+jjj];)
    if ( g.num_processes > 1 ) {
      MPI_Iallreduce( buffer, H[MAX(0,j-1)], (j+1)*nvec, MPI_COMPLEX_PRECISION, MPI_SUM,
                      (l->depth==0)?g.comm_cart:l->gs_PRECISION.level_comm, &req );
    } else {
      for( i=0; i<=j; i++ )
        VECTOR_LOOP(jj, nvec, jjj, H[MAX(0,j-1)][i*nvec+jj+jjj] = buffer[i*nvec+jj+jjj];)
    }
    PROF_PRECISION_STOP( _ALLR, 1 );
    END_MASTER(threading)
    
    apply_operator_PRECISION( &Z[j+1], &Z[j], p, l, threading );
    
    START_MASTER(threading)
    PROF_PRECISION_START( _ALLR );
    if ( g.num_processes > 1 ) {
      MPI_Wait( &req, &stat );
    }
    PROF_PRECISION_STOP( _ALLR, 0 );
    if ( j > 0 ) {
      for ( i=0; i<j; i++ )
        VECTOR_LOOP(jj, nvec, jjj, H[j-1][j*nvec+jj+jjj] -= conj( H[j-1][i*nvec+jj+jjj] )*H[j-1][i*nvec+jj+jjj];)
    }
    VECTOR_LOOP(jj, nvec, jjj, H[MAX(0,j-1)][j*nvec+jj+jjj] = sqrt( creal( H[MAX(0,j-1)][j*nvec+jj+jjj] ) );)
    END_MASTER(threading)
    SYNC_MASTER_TO_ALL(threading) 
    
    for( i=0; i<j; i++ )
      vector_PRECISION_saxpy( &V[j], &V[j], &V[i], H[j-1], i, -1, start, end, l );
    vector_PRECISION_real_scale( &V[j], &V[j], H[MAX(0,j-1)], j, 1, start, end, l );
    
    START_MASTER(threading)
    if ( j > 0 ) {
      VECTOR_LOOP(jj, nvec, jjj, H[j-1][(j-1)*nvec+jj+jjj] += sigma[jj+jjj];)
    }
    END_MASTER(threading)
    SYNC_MASTER_TO_ALL(threading)
    
    if ( j == 0 ) {
      VECTOR_LOOP(jj, nvec, jjj, sigma_sum += sigma[jj+jjj];)//or simply eliminate if-clause
      if ( sigma_sum ) vector_PRECISION_saxpy( &Z[j+1], &Z[j+1], &Z[j], sigma, 0, -1, start, end, l );
    } else {
      for( i=0; i<j; i++ )
        vector_PRECISION_saxpy( &Z[j+1], &Z[j+1], &Z[i+1], H[j-1], i, -1, start, end, l );
    }
    vector_PRECISION_real_scale( &Z[j+1], &Z[j+1], H[MAX(0,j-1)], j, 1, start, end, l );

  } else {
#endif
    warning0("arnoldi_step_PRECISION: with single step Arnoldi but not Pipelined, this has not been tested\n");
    SYNC_MASTER_TO_ALL(threading)

    PRECISION H_tot; sigma_sum = 0; 
    VECTOR_LOOP(jj, nvec, jjj, sigma_sum += sigma[jj+jjj];)//or simply eliminate if-clause 
    compute_core_start_end_custom(p->v_start, p->v_end, &start, &end, l, threading, l->num_lattice_site_var);

    V[j].num_vect_now = nvec; V[j+1].num_vect_now = nvec;  Z[j].num_vect_now = nvec;
    if ( prec != NULL ) {
      if ( p->kind == _LEFT ) {
        apply_operator_PRECISION( &Z[0], &V[j], p, l, threading );
        prec( &V[j+1], NULL, &Z[0], _NO_RES, l, threading );
	if ( sigma_sum ) vector_PRECISION_saxpy( &V[j+1], &V[j+1], &V[j], sigma, 0, -1, start, end, l );
      } else {
        if ( l->level == 0 ) {
          prec( &Z[j], NULL, &V[j], _NO_RES, l, threading );
          apply_operator_PRECISION( &V[j+1], &Z[j], p, l, threading );
        } else {
          if ( g.mixed_precision == 2 && (g.method >= 1 && g.method <= 2 ) ) {
            prec( &Z[j], &V[j+1], &V[j], _NO_RES, l, threading );
            // obtains w = D * Z[j] from Schwarz
          } else {
            prec( &Z[j], NULL, &V[j], _NO_RES, l, threading );
            apply_operator_PRECISION( &V[j+1], &Z[j], p, l, threading ); // w = D*Z[j]
          }
        }
	if ( sigma_sum ) vector_PRECISION_saxpy( &V[j+1], &V[j+1], &V[j], sigma, 0, -1, start, end, l );
      }
    } else {
      apply_operator_PRECISION( &V[j+1], &V[j], p, l, threading ); // w = D*V[j]
      if ( sigma_sum ) vector_PRECISION_saxpy( &V[j+1], &V[j+1], &V[j], sigma, 0, -1, start, end, l );
    }
    
    complex_PRECISION tmp[(j+2)*nvec];
    process_multi_inner_product_PRECISION( j+2, tmp, V, &V[j+1], p->v_start, p->v_end, l, threading );
    START_MASTER(threading)
    for( i=0; i<=j+1; i++ )
       VECTOR_LOOP(jj, nvec, jjj, buffer[i*nvec+jj+jjj] = tmp[i*nvec+jj+jjj];)
    
    if ( g.num_processes > 1 ) {
      PROF_PRECISION_START( _ALLR );
      MPI_Allreduce( buffer, H[j], (j+2)*nvec, MPI_COMPLEX_PRECISION, MPI_SUM, (l->depth==0)?g.comm_cart:l->gs_PRECISION.level_comm );
      PROF_PRECISION_STOP( _ALLR, 1 );
    } else {
      for( i=0; i<=j+1; i++ )
	VECTOR_LOOP(jj, nvec, jjj, H[j][i*nvec+jj+jjj] = buffer[i*nvec+jj+jjj];)
    }  
    for ( i=0; i<=j; i++ )
      VECTOR_LOOP(jj, nvec, jjj, H[j][(j+1)*nvec+jj+jjj] -= conj( H[j][i*nvec+jj+jjj] )*H[j][i*nvec+jj+jjj];)
    END_MASTER(threading)
    SYNC_MASTER_TO_ALL(threading)
    H_tot=0;
    VECTOR_LOOP(jj, nvec, jjj, H_tot += creal( H[j][(j+1)*nvec+jj+jjj] );)
    if ( H_tot < 0 )
      return 0;
    START_MASTER(threading)
    VECTOR_LOOP(jj, nvec, jjj, H[j][(j+1)*nvec+jj+jjj] = sqrt( creal( H[j][(j+1)*nvec+jj+jjj] ) );)
    END_MASTER(threading)
    SYNC_MASTER_TO_ALL(threading)

    for( i=0; i<=j; i++ )
      vector_PRECISION_saxpy( &V[j+1], &V[j+1], &V[i], H[j], i, -1, start, end, l );
    vector_PRECISION_real_scale( &V[j+1], &V[j+1], H[j], j+1, 1, start, end, l );
    START_LOCKED_MASTER(threading)
    VECTOR_LOOP(jj, nvec, jjj, H[j][j*nvec+jj+jjj] += sigma[jj+jjj];)
    END_LOCKED_MASTER(threading)
#ifdef PIPELINED_ARNOLDI
  }
#endif
#else
  SYNC_MASTER_TO_ALL(threading)
  // start and end indices for vector functions depending on thread
  // compute start and end indices for core
  // this puts zero for all other hyperthreads, so we can call functions below with all hyperthreads????
  compute_core_start_end_custom(p->v_start, p->v_end, &start, &end, l, threading, l->num_lattice_site_var);

  V[j].num_vect_now = nvec;  V[j+1].num_vect_now = nvec; Z[j].num_vect_now = nvec; // can move to fgmres????
  // apply preconditioner if provided
  if ( prec != NULL ) { 
    if ( p->kind == _LEFT ) {
      apply_operator_PRECISION( &Z[0], &V[j], p, l, threading ); 
      prec( w, NULL, &Z[0], _NO_RES, l, threading );
    } else {
      if ( l->level == 0 ) { 
        apply_operator_PRECISION( w, &Z[j], p, l, threading );
      } else {
        if ( g.mixed_precision == 2 && (g.method >= 1 && g.method <= 2 || g.method ==4) ) {
          prec( &Z[j], w, &V[j], _NO_RES, l, threading ); // obtains w = D * Z[j] from Schwarz
        } else {
          prec( &Z[j], NULL, &V[j], _NO_RES, l, threading );     // Z[j] = prec(V[j])
          apply_operator_PRECISION( w, &Z[j], p, l, threading ); // w    = D*Z[j]
        }
      }
    }
  // if precontioner is not set, just apply D to V[j]
  } else { 
    apply_operator_PRECISION( w, &V[j], p, l, threading ); // w = D*V[j]
  }

  // Compute the (j+1)^th column of Hessenberg matrix
  complex_PRECISION tmp[(j+1)*nvec];
  process_multi_inner_product_PRECISION( j+1, tmp, V, w, p->v_start, p->v_end, l, threading );
  START_MASTER(threading)
  for( i=0; i<=j; i++ )
    VECTOR_LOOP(jj, nvec, jjj, buffer[i*nvec+jj+jjj] = tmp[i*nvec+jj+jjj];)

  if ( g.num_processes > 1 ) {
    PROF_PRECISION_START( _ALLR );
    MPI_Allreduce( buffer, H[j], (j+1)*nvec, MPI_COMPLEX_PRECISION, MPI_SUM, (l->depth==0)?g.comm_cart:l->gs_PRECISION.level_comm );
    PROF_PRECISION_STOP( _ALLR, 1 );
  } else {
    for( i=0; i<=j; i++ )
      VECTOR_LOOP(jj, nvec, jjj, H[j][i*nvec+jj+jjj] = buffer[i*nvec+jj+jjj];)
  }
  END_MASTER(threading)
  SYNC_MASTER_TO_ALL(threading)
  for( i=0; i<=j; i++ )
    vector_PRECISION_saxpy( w, w, &V[i], H[j], i, -1, start, end, l );// w -= \sum_{i=0}^{j}H[j][i]*V[i]

#ifdef REORTH
  // re-orthogonalization
  process_multi_inner_product_PRECISION( j+1, tmp, V, w, p->v_start, p->v_end, l, threading );
  START_MASTER(threading)
  for( i=0; i<=j; i++ )
    VECTOR_LOOP(jj, nvec, jjj, buffer[i*nvec+jj+jjj] = tmp[i*nvec+jj+jjj];)
  
  if ( g.num_processes > 1 ) {
    PROF_PRECISION_START( _ALLR );
    MPI_Allreduce( buffer, tmp, (j+1)*nvec, MPI_COMPLEX_PRECISION, MPI_SUM, (l->depth==0)?g.comm_cart:l->gs_PRECISION.level_comm );
    PROF_PRECISION_STOP( _ALLR, 1 );
  }
  
  for( i=0; i<=j; i++ )
    VECTOR_LOOP(jj, nvec, jjj, H[j][i*nvec+jj+jjj] += tmp[i*nvec+jj+jjj];)

  END_MASTER(threading)
  SYNC_MASTER_TO_ALL(threading)
  for( i=0; i<=j; i++ )
    vector_PRECISION_saxpy( w, w, &V[i], tmp, i, -1, start, end, l );
#endif

  // Compute H[j+1][j]
  PRECISION tmp2[nvec]; 
  global_norm_PRECISION( tmp2, w, p->v_start, p->v_end, l, threading );

  START_MASTER(threading)
  VECTOR_LOOP(jj, nvec, jjj, H[j][(j+1)*nvec+jj+jjj] = tmp2[jj+jjj];)
  END_MASTER(threading)
  SYNC_MASTER_TO_ALL(threading)

  // compute V[j+1] = w / H[j+1][j]
  complex_PRECISION tmp3[nvec];
  VECTOR_LOOP(jj, nvec, jjj, tmp3[jj+jjj] = ( cabs_PRECISION( H[j][(j+1)*nvec+jj+jjj] ) > 1e-15 )?tmp2[jj+jjj]:1;)
  vector_PRECISION_real_scale( &V[j+1], w, tmp3, 0, 1, start, end, l );
#endif
  return 1;
}


void qr_update_PRECISION( complex_PRECISION **H, complex_PRECISION *s,
                          complex_PRECISION *c, complex_PRECISION *gamma, int j,
                          level_struct *l, struct Thread *threading ) {

/*********************************************************************************
* Applies one Givens rotation to the Hessenberg matrix H in order to solve the 
* least squares problem in (F)GMRES for computing the solution.
* - complex_PRECISION **H: Hessenberg matrix from Arnoldi decomposition
* - complex_PRECISION *s: sin values from givens rotations
* - complex_PRECISION *c: cos valies from givens rotations
* - complex_PRECISION *gamma: Approximation to residual from every step
* - int j: Denotes current iteration.
*********************************************************************************/  

  SYNC_HYPERTHREADS(threading)
  SYNC_CORES(threading)
  START_MASTER(threading)
  
  PROF_PRECISION_START( _SMALL1 );
  
  int i, n, jj, n_vect=num_loop;//g.num_vect_now;// may want to find a better way!!!!!
  complex_PRECISION beta[n_vect];
  
  // update QR factorization
  // apply previous Givens rotation
  for ( i=0; i<j; i++ ) {
    VECTOR_LOOP(n, n_vect, jj, beta[n+jj] = (-s[i*n_vect+n+jj])*H[j][i*n_vect+n+jj] + (c[i*n_vect+n+jj])*H[j][(i+1)*n_vect+n+jj];)
    VECTOR_LOOP(n, n_vect, jj, H[j][i*n_vect+n+jj] = conj_PRECISION(c[i*n_vect+n+jj])*H[j][i*n_vect+n+jj] + conj_PRECISION(s[i*n_vect+n+jj])*H[j][(i+1)*n_vect+n+jj];)
    VECTOR_LOOP(n, n_vect, jj, H[j][(i+1)*n_vect+n+jj] = beta[n+jj];)
  }

  // compute current Givens rotation angles
  VECTOR_LOOP(n, n_vect, jj, beta[n+jj] = (complex_PRECISION) sqrt( NORM_SQUARE_PRECISION(H[j][j*n_vect+n+jj]) + NORM_SQUARE_PRECISION(H[j][(j+1)*n_vect+n+jj]) );)
  VECTOR_LOOP(n, n_vect, jj, s[j*n_vect+n+jj] = H[j][(j+1)*n_vect+n+jj]/beta[n+jj];)
  VECTOR_LOOP(n, n_vect, jj, c[j*n_vect+n+jj] = H[j][j*n_vect+n+jj]/beta[n+jj];)
  // update right column
  VECTOR_LOOP(n, n_vect, jj, gamma[(j+1)*n_vect+n+jj] = (-s[j*n_vect+n+jj])*gamma[j*n_vect+n+jj];)
  VECTOR_LOOP(n, n_vect, jj, gamma[j*n_vect+n+jj] = conj_PRECISION(c[j*n_vect+n+jj])*gamma[j*n_vect+n+jj];)
  // apply current Givens rotation
  VECTOR_LOOP(n, n_vect, jj, H[j][j*n_vect+n+jj] = beta[n+jj];)
  VECTOR_LOOP(n, n_vect, jj, H[j][(j+1)*n_vect+n+jj] = 0;)
  
  PROF_PRECISION_STOP( _SMALL1, 6*j+6 );
  
  END_MASTER(threading)
  SYNC_MASTER_TO_ALL(threading)
    // what about hyperthread????
}

void compute_solution_PRECISION( vector_PRECISION *x, vector_PRECISION *V, complex_PRECISION *y,
                                 complex_PRECISION *gamma, complex_PRECISION **H, int j, int ol,
                                 gmres_PRECISION_struct *p, level_struct *l, struct Thread *threading ) {
  
  int i, k, n, jj, n_vect = num_loop;
  // start and end indices for vector functions depending on thread
  int start;
  int end;
  // compute start and end indices for core
  // this puts zero for all other hyperthreads, so we can call functions below with all hyperthreads
  compute_core_start_end_custom(p->v_start, p->v_end, &start, &end, l, threading, l->num_lattice_site_var );//!!!!!!!

  START_MASTER(threading)
  
  PROF_PRECISION_START( _SMALL2 );

  // backward substitution
  for ( i=j; i>=0; i-- ) {
    VECTOR_LOOP(n, n_vect, jj, y[i*n_vect+n+jj] = gamma[i*n_vect+n+jj];)
      for ( k=i+1; k<=j; k++ )
	VECTOR_LOOP(n, n_vect, jj, y[i*n_vect+n+jj] -= H[k][i*n_vect+n+jj]*y[k*n_vect+n+jj];)
    VECTOR_LOOP(n, n_vect, jj, y[i*n_vect+n+jj] /= H[i][i*n_vect+n+jj];)
  }
  
  PROF_PRECISION_STOP( _SMALL2, ((j+1)*(j+2))/2 + j+1 );
  
  END_MASTER(threading)
  SYNC_MASTER_TO_ALL(threading)
  
  // x = x + V*y
  if ( ol ) {
    for ( i=0; i<=j; i++ ) {
      vector_PRECISION_saxpy( x, x, &V[i], y, i, 1, start, end, l ); // x += y_i * V_i where V_i is the i^th basis vector
    }
  } else {
    vector_PRECISION_scale( x, &V[0], y, 0, start, end, l );         // x  = y_0 * V_0 where V_0 is the first basis vector
    for ( i=1; i<=j; i++ ) {
      vector_PRECISION_saxpy( x, x, &V[i], y, i, 1, start, end, l ); // x += y_i * V_i where V_i is the i^th basis vector
    }
  }
}

void local_minres_PRECISION( vector_PRECISION *phi, vector_PRECISION *eta, vector_PRECISION *latest_iter,
                             int start, schwarz_PRECISION_struct *s, level_struct *l, struct Thread *threading ) {
  
/*********************************************************************************
* Minimal Residual iteration solver used to solve the block systems
*     blockD varphi = r = eta ... (1)
* within the Schwarz method.
* In SAP, we solve D*phi = eta_0 ... (2) via phi <- phi + varphi where varphi = blockD^{-1}*r
* In SAP, eta = r = eta_0 - D*phi.  This function computes the increment varphi and update r
* Assume: initial guess for Eq. (1)is zero.
* phi: (in) current estimate of Eq. (2); (return) updated estimate after l->block_iter steps of MinRes
* eta: (in) initial residual; (return) overwritten by the block residual r.
* latest_iter: the increment to the estimate in the last step of MinRes applied to Eq. (1) = varphi
* Note: The increment in the update eq. for Eq. (2), i.e., latest_iter = "phi - phi_old"
*       is returned to calculate the missing contributions to r on the current Schwarz block
*       coming from outside of the block, i.e., block boundaries and complete the residual 
*       update in SAP at a cheaper cost.
*********************************************************************************/
  
  START_UNTHREADED_FUNCTION(threading)

  int i, nvec = eta->num_vect_now;
  int nv  = l->num_lattice_site_var;
  int n   = l->block_iter;
  int end = (g.odd_even&&l->depth==0)?(start+nv*s->num_block_even_sites):(start+s->block_vector_size);
  vector_PRECISION Dr, r, lphi;
 
  if ( nvec != num_loop ) //g.num_vect_now )//!!!!!!!!
    error0("local_minres_PRECISION: incosistent number of vectors %d %d\n",nvec,g.num_vect_now);

  complex_PRECISION alpha[nvec];
  void (*block_op)() = (l->depth==0)?(g.odd_even?apply_block_schur_complement_PRECISION:block_d_plus_clover_PRECISION)
                                    :coarse_block_operator_PRECISION;
  Dr.num_vect   = s->buf[0].num_vect;
  r.num_vect    = s->buf[0].num_vect;
  lphi.num_vect = s->buf[0].num_vect;

  Dr.num_vect_now   = nvec;
  r.num_vect_now    = nvec;
  lphi.num_vect_now = nvec;

  Dr.vector_buffer   = s->local_minres_buffer[0]; // _SCHWARZ size 
  r.vector_buffer    = s->local_minres_buffer[1]; // _SCHWARZ size 
  lphi.vector_buffer = s->local_minres_buffer[2]; // _SCHWARZ size 

  vector_PRECISION_copy( &r, eta, start, end, l );    // r = eta
  vector_PRECISION_define( &lphi, 0, start, end, l ); // phi = 0
  
  for ( i=0; i<n; i++ ) {
    // Dr = blockD*r
    block_op( &Dr, &r, start, s, l, no_threading );
    // alpha = <Dr,r>/<Dr,Dr>
    local_xy_over_xx_PRECISION( alpha, &Dr, &r, start, end, l );
    // phi += alpha * r
    vector_PRECISION_saxpy( &lphi, &lphi, &r, alpha, 0, 1, start, end, l );
    // r -= alpha * Dr
    vector_PRECISION_saxpy( &r, &r, &Dr, alpha, 0, -1, start, end, l );
  }
  
  if ( latest_iter != NULL ) vector_PRECISION_copy( latest_iter, &lphi, start, end, l );
  if ( phi != NULL ) vector_PRECISION_plus( phi, phi, &lphi, start, end, l );
  vector_PRECISION_copy( eta, &r, start, end, l );

  END_UNTHREADED_FUNCTION(threading)
}
