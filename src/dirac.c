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

void dirac_setup( config_double hopp, level_struct *l ) {

/*********************************************************************************
* Sets up the gauge matrices + clover term for the dirac operator and calculates 
* the plaquette of the configuration. (only those that need to be precomputed)
* config_double hopp: Vector containing all entries of the gauge configuration.
*********************************************************************************/

  double t0, t1;
  int i, j, t, z, y, x, mu;
  SU3_storage U = NULL;
  complex_double phase[4];
  int *ll=l->local_lattice, onb[4];
  for (i=0; i<4; i++)
    onb[i] = (g.my_coords[i]==g.process_grid[i]-1)?1:0;
  
  if ( g.print > 0 ) printf0("%s\n", CLIFFORD_BASIS );
  if ( g.bc == _ANTIPERIODIC ) printf0("antiperiodic in time");
  else if ( g.bc == _TWISTED ) printf0("twisted (%.2f, %.2f, %.2f, %.2f)", g.twisted_bc[0], 
               g.twisted_bc[1], g.twisted_bc[2], g.twisted_bc[3]);
  else printf0("periodic in time");
  printf0(" boundary conditions\n");

  t0 = MPI_Wtime();

  // read and store configuration

  SU3_storage_alloc( &U, l );

  if( g.bc == _ANTIPERIODIC && onb[T] ) {
    phase[Z] = 1; phase[Y] = 1; phase[X] = 1;
    for ( t=1, i=0; t<ll[T]+1; t++ ) {
      if (t<ll[T]) phase[T] = 1; 
      else phase[T] = -1;
      for ( z=1; z<ll[Z]+1; z++ )
        for ( y=1; y<ll[Y]+1; y++ )
          for ( x=1; x<ll[X]+1; x++ )
            for ( mu=0; mu<4; mu++ )
              for (j=0; j<9; j++, i++) {
                g.op_double.D[i] = 0.5*phase[mu]*hopp[i];
                U[t][z][y][x][mu][j] = phase[mu]*hopp[i];
              }
    }
  }
  else if( g.bc == _TWISTED && ( onb[T] || onb[Z] || onb[Y] || onb[X] )) {
    //TODO
    warning0("Twisted boundary conditions not supported outside the library.\n");
    for ( t=1, i=0; t<ll[T]+1; t++ ) {
      if ( !onb[T] || t<ll[T] || g.twisted_bc[T]==0) phase[T] = 1; 
      else phase[T] = cexp(I*g.twisted_bc[T]);
      for ( z=1; z<ll[Z]+1; z++ ) {
        if ( !onb[Z] || z<ll[Z] || g.twisted_bc[Z]==0) phase[Z] = 1; 
        else phase[Z] = cexp(I*g.twisted_bc[Z]);
        for ( y=1; y<ll[Y]+1; y++ ) {
          if ( !onb[Y] || y<ll[Y] || g.twisted_bc[Y]==0) phase[Y] = 1; 
          else phase[Y] = cexp(-I*g.twisted_bc[Y]);
          for ( x=1; x<ll[X]+1; x++ ) {
            if ( !onb[X] || x<ll[X] || g.twisted_bc[X]==0) phase[X] = 1; 
            else phase[X] = cexp(-I*g.twisted_bc[X]);
            for ( mu=0; mu<4; mu++ ) 
              for (j=0; j<9; j++, i++) {
                g.op_double.D[i] = 0.5*phase[mu]*hopp[i];
                U[t][z][y][x][mu][j] = phase[mu]*hopp[i];
              }
          }
        }
      }
    }
  }

  else
    for ( t=1, i=0; t<ll[T]+1; t++ )
      for ( z=1; z<ll[Z]+1; z++ )
        for ( y=1; y<ll[Y]+1; y++ )
          for ( x=1; x<ll[X]+1; x++ )
            for ( mu=0; mu<4; mu++ )
              for (j=0; j<9; j++, i++) {
                g.op_double.D[i] = 0.5*hopp[i];
                U[t][z][y][x][mu][j] = hopp[i];
              }
  
  SU3_ghost_update( &U, l );
  if ( g.print > 0 ) printf0("Configuration stored...\n");
  
  compute_clover_term( U, l );
  
  // calculate the plaquette
  g.plaq_clov = calc_plaq( U, l );
  if(g.print > 0) printf0("average plaquette: %.13lf\n", g.plaq_clov);
    
  SU3_storage_free( &U, l );
  
  t1 = MPI_Wtime();
  if ( g.print > 0 ) {    
    printf0("\n+----------------------------------------------------------+\n");
    printf0("| read in and set up the parallel dirac operator           |\n");
    printf0("| elapsed wall clock time: %-12g seconds            |\n", t1-t0 );
    printf0("+----------------------------------------------------------+\n");
  }
}


double *dirac_setup_get_gauge_pointer() {
  return (double *)g.op_double.D;
}
double *dirac_setup_get_clover_pointer() {
  return (double *)g.op_double.clover;
}


void mat_alloc( complex_double **A, int n )  {
  *A = NULL;
  MALLOC( (*A), complex_double, n*n );
}


void mat_free( complex_double **A, int n ) {
  FREE( (*A), complex_double, n*n );
}


void spin_alloc( int num_spin, int n ) {
  int i;
  MALLOC( g.gamma, complex_double*, num_spin );
  for ( i=0; i<num_spin; i++ )
    mat_alloc( &(g.gamma[i]), n );
}


void spin_free( int num_spin, int n ) {
  int i;
  for ( i=0; i<num_spin; i++ )
    mat_free( &(g.gamma[i]), n );
  FREE( g.gamma, complex_double*, num_spin );
}


void spin_define( void ) {
  
  // assertion: gamma5 = (+/-)diag( -1, -1, 1, 1 )
  int i, j;
  
  for ( i=0; i<4; i++ )
    for ( j=0; j<16; j++ )
      g.gamma[i][j] = _COMPLEX_double_ZERO;
  
  g.gamma[T][   GAMMA_T_SPIN0_CO] = GAMMA_T_SPIN0_VAL;
  g.gamma[T][ 4+GAMMA_T_SPIN1_CO] = GAMMA_T_SPIN1_VAL;
  g.gamma[T][ 8+GAMMA_T_SPIN2_CO] = GAMMA_T_SPIN2_VAL;
  g.gamma[T][12+GAMMA_T_SPIN3_CO] = GAMMA_T_SPIN3_VAL;
  
  g.gamma[Z][   GAMMA_Z_SPIN0_CO] = GAMMA_Z_SPIN0_VAL;
  g.gamma[Z][ 4+GAMMA_Z_SPIN1_CO] = GAMMA_Z_SPIN1_VAL;
  g.gamma[Z][ 8+GAMMA_Z_SPIN2_CO] = GAMMA_Z_SPIN2_VAL;
  g.gamma[Z][12+GAMMA_Z_SPIN3_CO] = GAMMA_Z_SPIN3_VAL;
  
  g.gamma[Y][   GAMMA_Y_SPIN0_CO] = GAMMA_Y_SPIN0_VAL;
  g.gamma[Y][ 4+GAMMA_Y_SPIN1_CO] = GAMMA_Y_SPIN1_VAL;
  g.gamma[Y][ 8+GAMMA_Y_SPIN2_CO] = GAMMA_Y_SPIN2_VAL;
  g.gamma[Y][12+GAMMA_Y_SPIN3_CO] = GAMMA_Y_SPIN3_VAL;
  
  g.gamma[X][   GAMMA_X_SPIN0_CO] = GAMMA_X_SPIN0_VAL;
  g.gamma[X][ 4+GAMMA_X_SPIN1_CO] = GAMMA_X_SPIN1_VAL;
  g.gamma[X][ 8+GAMMA_X_SPIN2_CO] = GAMMA_X_SPIN2_VAL;
  g.gamma[X][12+GAMMA_X_SPIN3_CO] = GAMMA_X_SPIN3_VAL;
}

// tensor += lambda * A (x) B --- where A n times n and B m times m
static inline void mat_tensor( complex_double *tensor, complex_double lambda, complex_double *A, int n, complex_double *B, int m ) {
  int i1, i2, j1, j2;
  for ( i1=0; i1<n; i1++ )
    for ( i2=0; i2<n; i2++ )
      for ( j1=0; j1<m; j1++ )
        for ( j2=0; j2<m; j2++ )
          tensor[n*m*(i1*m+j1)+i2*m+j2] += lambda*A[n*i1+i2]*B[m*j1+j2];
}

// C += A*B
static inline void addMatMul( complex_double *C, complex_double *A, complex_double *B, int n ) {
  int i, j, k;
  for ( i=0; i<n; i++ )
    for ( j=0; j<n; j++ )
      for ( k=0; k<n; k++ )
        C[n*i+j] += A[n*i+k]*B[n*k+j];
}


static inline void addMatMulHermit2( complex_double *C, complex_double *A, complex_double *B, int n ) {
  int i, j, k;
  for ( i=0; i<n; i++ )
    for ( j=0; j<n; j++ )
      for ( k=0; k<n; k++ )
        C[n*i+j] += A[n*i+k]*conj(B[n*j+k]);
}


static inline void addMatMulHermit1( complex_double *C, complex_double *A, complex_double *B, int n ) {
  int i, j, k;
  for ( i=0; i<n; i++ )
    for ( j=0; j<n; j++ )
      for ( k=0; k<n; k++ )
        C[n*i+j] += conj(A[n*k+i])*B[n*k+j];
}


static inline void addMatMulHermit12( complex_double *C, complex_double *A, complex_double *B, int n ) {
  int i, j, k;
  for ( i=0; i<n; i++ )
    for ( j=0; j<n; j++ )
      for ( k=0; k<n; k++ )
        C[n*i+j] += conj(A[n*k+i])*conj(B[n*j+k]);
}

// C := alpha*A
static inline void scaleMat( complex_double *C, complex_double *A, complex_double alpha, int n ) {
  int i, m=n*n;
  for ( i=0; i<m; i++ )
    C[i] = alpha*A[i];
}

// C := A-B
static inline void subtMat( complex_double *C, complex_double *A, complex_double *B, int n ) {
  int i, m=n*n;
  for ( i=0; i<m; i++ )
    C[i] = A[i]-B[i];
}

// A := 0
void zeroMat( complex_double *A, int n ) {
  int i, m=n*n;
  for ( i=0; i<m; i++ )
    A[i]=0;
}


void Q( complex_double *Qstore, int mu, int nu, int t, int z, int y, int x, SU3_storage U ) {
    
  complex_double *tmp1 = NULL, *tmp2 = NULL;
  int mH[4][4] = {{ 1,0,0,0 },{ 0,1,0,0 },{ 0,0,1,0 },{ 0,0,0,1 }};
  
  mat_alloc( &tmp1, 3 );
  mat_alloc( &tmp2, 3 );  
  
  // plaquette no1: mu, nu
  zeroMat( Qstore, 3 );
  zeroMat( tmp1, 3 );
  zeroMat( tmp2, 3 );
  addMatMul( tmp1, U[t][z][y][x][mu],
             U[t+mH[mu][T] ][z+mH[mu][Z] ][y+mH[mu][Y] ][x+mH[mu][X] ][nu], 3 );
  addMatMulHermit2( tmp2, tmp1,
                    U[t+mH[nu][T] ][z+mH[nu][Z] ][y+mH[nu][Y] ][x+mH[nu][X] ][mu], 3 );
  addMatMulHermit2( Qstore, tmp2, U[t][z][y][x][nu], 3 );
  
  // plaquette no2: nu, -mu
  zeroMat( tmp1, 3 );
  zeroMat( tmp2, 3 );
  addMatMulHermit2( tmp1, U[t][z][y][x][nu],
                    U[t+mH[nu][T]-mH[mu][T] ][z+mH[nu][Z]-mH[mu][Z] ][y+mH[nu][Y]-mH[mu][Y] ][x+mH[nu][X]-mH[mu][X] ][mu], 3 );
  addMatMulHermit2( tmp2, tmp1,
                    U[t-mH[mu][T] ][z-mH[mu][Z] ][y-mH[mu][Y] ][x-mH[mu][X] ][nu], 3 );
  addMatMul( Qstore, tmp2,
             U[t-mH[mu][T] ][z-mH[mu][Z] ][y-mH[mu][Y] ][x-mH[mu][X] ][mu], 3 );
  
  // plaquette no3: -mu, -nu
  zeroMat( tmp1, 3 );
  zeroMat( tmp2, 3 );
  addMatMulHermit12( tmp1, U[t-mH[mu][T] ][z-mH[mu][Z] ][y-mH[mu][Y] ][x-mH[mu][X] ][mu],
                     U[t-mH[mu][T]-mH[nu][T] ][z-mH[mu][Z]-mH[nu][Z] ][y-mH[mu][Y]-mH[nu][Y] ][x-mH[mu][X]-mH[nu][X] ][nu], 3 );
  addMatMul( tmp2, tmp1,
             U[t-mH[mu][T]-mH[nu][T] ][z-mH[mu][Z]-mH[nu][Z] ][y-mH[mu][Y]-mH[nu][Y] ][x-mH[mu][X]-mH[nu][X] ][mu], 3 );   
  addMatMul( Qstore, tmp2,
             U[t-mH[nu][T] ][z-mH[nu][Z] ][y-mH[nu][Y] ][x-mH[nu][X] ][nu], 3 );
  
  // plaquette no4: -nu, mu
  zeroMat( tmp1, 3 );
  zeroMat( tmp2, 3 );
  addMatMulHermit1( tmp1, U[t-mH[nu][T] ][z-mH[nu][Z] ][y-mH[nu][Y] ][x-mH[nu][X] ][nu],
                    U[t-mH[nu][T] ][z-mH[nu][Z] ][y-mH[nu][Y] ][x-mH[nu][X] ][mu], 3 );
  addMatMul( tmp2, tmp1,
             U[t-mH[nu][T]+mH[mu][T] ][z-mH[nu][Z]+mH[mu][Z] ][y-mH[nu][Y]+mH[mu][Y] ][x-mH[nu][X]+mH[mu][X] ][nu], 3 );     
  addMatMulHermit2( Qstore, tmp2,
                    U[t][z][y][x][mu], 3 );
  scaleMat( Qstore, Qstore, 1.0/16.0, 3 );
  
  mat_free( &tmp1, 3 );
  mat_free( &tmp2, 3 );
}


void Qdiff( complex_double *Qstore, int mu, int nu, int t, int z, int y, int x, SU3_storage U ) {
  
  complex_double *tmp1 = NULL, *tmp2 = NULL;
  
  mat_alloc( &tmp1, 3 );
  mat_alloc( &tmp2, 3 );
  
  Q( tmp1, mu, nu, t, z, y, x, U );
  Q( tmp2, nu, mu, t, z, y, x, U );
  subtMat( Qstore, tmp1, tmp2, 3 );
  
  mat_free( &tmp1, 3 );
  mat_free( &tmp2, 3 );
}


void set_clover( complex_double *Qstore, int mu, int nu, int index, config_double clover ) {
  
  complex_double *gamma_mu_gamma_nu = NULL, *tmp = NULL;
  int i, j, k;
  
  mat_alloc( &tmp, 12 );
  mat_alloc( &gamma_mu_gamma_nu, 4 );
  zeroMat( tmp, 12 );
  zeroMat( gamma_mu_gamma_nu, 4 );
  addMatMul( gamma_mu_gamma_nu, g.gamma[mu], g.gamma[nu], 4 );
  mat_tensor( tmp, -g.csw, gamma_mu_gamma_nu, 4, Qstore, 3 );
  
  // diagonal including the shift
  for ( k=0; k<12; k++)
    clover[42*index+k] += tmp[13*k];
  
  // spin 1, spin 2
  for ( i=0, k=12; i<6; i++)
    for ( j=i+1; j<6; j++, k++ )
      clover[42*index+k] += tmp[12*i+j];
  
  // spin 3, spin 4
  for ( i=6, k=27; i<12; i++)
    for ( j=i+1; j<12; j++, k++ )
      clover[42*index+k] += tmp[12*i+j];
      
  mat_free( &gamma_mu_gamma_nu, 4 );
  mat_free( &tmp, 12 );
}


// compute self-coupling terms
void compute_clover_term ( SU3_storage U, level_struct *l ) {
  int i, j, t, z, y, x, mu, nu;
  operator_double_struct *op = &(g.op_double);

  op->m0 = g.m0;
  
  for ( mu=0; mu<4; mu++ )  
    op->oe_offset += (l->local_lattice[mu]*(g.my_coords[mu]/l->comm_offset[mu]))%2;
  op->oe_offset = op->oe_offset%2;
  
  for ( i=0,t=0; t<l->local_lattice[T]; t++ )
    for ( z=0; z<l->local_lattice[Z]; z++ )
      for ( y=0; y<l->local_lattice[Y]; y++ )
        for ( x=0; x<l->local_lattice[X]; x++ ){
          if((t+z+y+x+op->oe_offset)%2) { //odd
            FOR12(op->odd_proj[i] = 1; i++;);
          } else {
            FOR12(op->odd_proj[i] = _COMPLEX_double_ZERO; i++;);
          }
        }
  
#ifdef HAVE_TM
#ifdef HAVE_MULT_TM
  if ( !g.is_even_shifted_mu_nonzero && g.mu + g.mu_odd_shift == 0 ) {// tm_term is zero in this case
    buffer_double_define( op->tm_term, _COMPLEX_double_ZERO, 0, l->inner_vector_size, l );
  }
  else
#endif
    tm_term_double_setup( g.mu, g.mu_even_shift, g.mu_odd_shift, g.mu_factor[l->depth], op, l, no_threading );
#endif

#ifdef HAVE_TM1p1
  if ( g.epsbar == 0 && g.epsbar_ig5_even_shift == 0 && g.epsbar_ig5_odd_shift == 0 ) 
    buffer_double_define( op->epsbar_term, _COMPLEX_double_ZERO, 0, l->inner_vector_size, l );
  else
    epsbar_term_double_setup( g.epsbar, g.epsbar_ig5_even_shift, g.epsbar_ig5_odd_shift, op, l, no_threading );  
#endif

  // generate clover term
  if ( g.csw != 0.0 ) {
    spin_alloc( 4, 4 );
    spin_define();
    complex_double *Qstore = NULL;
    mat_alloc( &Qstore, 3 );
    
    j = 42*l->num_inner_lattice_sites;
    for ( i=0; i<j; i++ )
      op->clover[i] = 0;
    i = 0;
    for ( t=1; t<l->local_lattice[T]+1; t++ )
      for ( z=1; z<l->local_lattice[Z]+1; z++ )
        for ( y=1; y<l->local_lattice[Y]+1; y++ )
          for ( x=1; x<l->local_lattice[X]+1; x++ ) {
            // diagonal including the shift
            for ( j=0; j<12; j++)
              op->clover[42*i+j] = 4+op->m0;
            
            for ( mu=0; mu<4; mu++ )
              for ( nu=mu+1; nu<4; nu++ ) {
                Qdiff( Qstore, mu, nu, t, z, y, x, U );
                set_clover( Qstore, mu, nu, i, op->clover );
              }
              i++;
          }
    
    mat_free( &Qstore, 3 );
    spin_free( 4, 4 );
  } else {
    buffer_double_define( op->clover, 4+op->m0, 0, l->inner_vector_size, l );
  }
}


void SU3_ghost_update( SU3_storage *U, level_struct *l ) {
  
/*********************************************************************************
* Updates ghost cells of a gauge matrix U.
* - SU3_storage *U: Current gauge matrix which sends and recieves ghost cells.
*********************************************************************************/  
  
  int t, z, y, x, mu, nu, *ll = l->local_lattice, ls[4], le[4];
  long int i, j, send_size, max_size;
  buffer_double buffer1 = NULL, buffer2 = NULL, buffer3 = NULL, buffer4 = NULL;

  max_size = 0;
  for ( mu=0; mu<4; mu++ ) {
    send_size=36;
    for ( nu=0; nu<4; nu++ ) {
      if (nu < mu) send_size *= ll[nu]+2;
      if (nu > mu) send_size *= ll[nu];
    }
    if (send_size > max_size) max_size = send_size;
  }
  
  MALLOC( buffer1, complex_double, max_size );
  MALLOC( buffer2, complex_double, max_size );
  MALLOC( buffer3, complex_double, max_size );
  MALLOC( buffer4, complex_double, max_size );
  
  for ( mu=0; mu<4; mu++ ) {
    ls[mu] = 1;
    le[mu] = ll[mu]+1;
  }
  
  for ( mu=0; mu<4; mu++ ) {
    // send own negative inner boundary
    le[mu] = 2;
    i = 0;
    for ( t=ls[T]; t<le[T]; t++ )
      for ( z=ls[Z]; z<le[Z]; z++ )
        for ( y=ls[Y]; y<le[Y]; y++ )
          for ( x=ls[X]; x<le[X]; x++ ) {
            j = 0;
            FOR36( buffer1[i] = *((*U)[t][z][y][x][0]+j); i++; j++; )
          }
    le[mu] = ll[mu]+1;
    send_size = i;
    ASSERT(send_size<=max_size);
    MPI_Irecv( buffer3, send_size, MPI_COMPLEX_double, l->neighbor_rank[2*mu], 2*mu, g.comm_cart, &(g.rreqs[2*mu]) );
    MPI_Isend( buffer1, send_size, MPI_COMPLEX_double, l->neighbor_rank[2*mu+1], 2*mu, g.comm_cart, &(g.sreqs[2*mu]) );
    
    // send own positive inner boundary
    ls[mu] = ll[mu];
    i = 0;
    for ( t=ls[T]; t<le[T]; t++ )
      for ( z=ls[Z]; z<le[Z]; z++ )
        for ( y=ls[Y]; y<le[Y]; y++ )
          for ( x=ls[X]; x<le[X]; x++ ) {
            j = 0;
            FOR36( buffer2[i] = *((*U)[t][z][y][x][0]+j); i++; j++; )
          }
    ls[mu] = 1;
    send_size = i;
    ASSERT(send_size<=max_size);
    MPI_Irecv( buffer4, send_size, MPI_COMPLEX_double, l->neighbor_rank[2*mu+1], 2*mu+1, g.comm_cart, &(g.rreqs[2*mu+1]) );
    MPI_Isend( buffer2, send_size, MPI_COMPLEX_double, l->neighbor_rank[2*mu], 2*mu+1, g.comm_cart, &(g.sreqs[2*mu+1]) );
    
    //recv own positive boundary
    MPI_Wait( &(g.sreqs[2*mu]), MPI_STATUS_IGNORE );
    MPI_Wait( &(g.rreqs[2*mu]), MPI_STATUS_IGNORE );
    le[mu] = ll[mu]+2;
    ls[mu] = ll[mu]+1;
    i = 0;
    for ( t=ls[T]; t<le[T]; t++ )
      for ( z=ls[Z]; z<le[Z]; z++ )
        for ( y=ls[Y]; y<le[Y]; y++ )
          for ( x=ls[X]; x<le[X]; x++ ) {
            j = 0;
            FOR36( *((*U)[t][z][y][x][0]+j) = buffer3[i]; i++; j++; )
          }
    le[mu] = ll[mu]+1;
    ls[mu] = 1;
    send_size = i;
    ASSERT(send_size<=max_size);
    
    // recv own negative boundary
    MPI_Wait( &(g.sreqs[2*mu+1]), MPI_STATUS_IGNORE );
    MPI_Wait( &(g.rreqs[2*mu+1]), MPI_STATUS_IGNORE );
    le[mu] = 1;
    ls[mu] = 0;
    i = 0;
    for ( t=ls[T]; t<le[T]; t++ )
      for ( z=ls[Z]; z<le[Z]; z++ )
        for ( y=ls[Y]; y<le[Y]; y++ )
          for ( x=ls[X]; x<le[X]; x++ ) {
            j = 0;
            FOR36( *((*U)[t][z][y][x][0]+j) = buffer4[i]; i++; j++; )
          }
    // extending, we send also the corners
    le[mu] = ll[mu]+2;
    ls[mu] = 0;
    send_size = i;
    ASSERT(send_size<=max_size);
  }
  FREE( buffer1, complex_double, max_size );		
  FREE( buffer2, complex_double, max_size );		
  FREE( buffer3, complex_double, max_size );		
  FREE( buffer4, complex_double, max_size );

}


void SU3_storage_alloc( SU3_storage *U, level_struct *l ) {

  int t, z, y, x, mu;
  long int lsize[4];
  complex_double *field = NULL;
  
  for (mu=0; mu<4; mu++)
    lsize[mu]=l->local_lattice[mu]+2;

  MALLOC( field, complex_double, lsize[T]*lsize[Z]*lsize[Y]*lsize[X]*36);

  MALLOC( (*U), complex_double*****, lsize[T] );
  for (t=0; t<lsize[T]; t++) {
    (*U)[t] = NULL;
    MALLOC( (*U)[t], complex_double****, lsize[Z] );
    for (z=0; z<lsize[Z]; z++) {
      (*U)[t][z] = NULL;
      MALLOC( (*U)[t][z], complex_double***, lsize[Y] );
      for (y=0; y<lsize[Y]; y++) {
        (*U)[t][z][y] = NULL;
        MALLOC( (*U)[t][z][y], complex_double**, lsize[X] );
        for (x=0; x<lsize[X]; x++) {
          (*U)[t][z][y][x] = NULL;
          MALLOC( (*U)[t][z][y][x], complex_double*, 4 );
          for (mu=0; mu<4; mu++) {
            (*U)[t][z][y][x][mu] = 
            &field[ t*lsize[Z]*lsize[Y]*lsize[X]*36
            + z*lsize[Y]*lsize[X]*36
            + y*lsize[X]*36
            + x*36
            + mu*9 ];
          }
        }
      }
    }
  }
}


void SU3_storage_free( SU3_storage *U, level_struct *l ) {

  int t,z,y,x,mu;
  int lsize[4];
  for (mu=0;mu<4; mu++)
    lsize[mu]=l->local_lattice[mu]+2;

  FREE( (*U)[0][0][0][0][0], complex_double, lsize[T]*lsize[Z]*lsize[Y]*lsize[X]*36);
  for (t=0; t<lsize[T]; t++) {
    for (z=0; z<lsize[Z]; z++) {
      for (y=0; y<lsize[Y]; y++) {
        for (x=0; x<lsize[X]; x++) {
          FREE( (*U)[t][z][y][x], complex_double*, 4 );
        }
        FREE( (*U)[t][z][y], complex_double**, lsize[X] );
      }
      FREE( (*U)[t][z], complex_double***, lsize[Y] );
    }
    FREE( (*U)[t], complex_double****, lsize[Z] );
  }
  FREE( (*U), complex_double*****, lsize[T] );
}

double calc_plaq( SU3_storage U, level_struct *l ) {
  
  int t, z, y, x, mu, nu, ls[4], mH[4][4]  = {{ 1,0,0,0 },{ 0,1,0,0 },{ 0,0,1,0 },{ 0,0,0,1 }};
  complex_double plaq, averagePlaq, *su3Prod, *tmp1, *tmp2;
  long int k;
  
  mat_alloc( &su3Prod, 3 );
  mat_alloc( &tmp1, 3 );
  mat_alloc( &tmp2, 3 );
  
  for (mu=0; mu<4; mu++) {
    ls[mu] = 1;
  }
  
  plaq = 0; averagePlaq=0;
  k = l->global_lattice[0]*l->global_lattice[1]*l->global_lattice[2]*l->global_lattice[3];
  
  for (t=ls[0]; t<l->local_lattice[0]+ls[0]; t++)
    for (z=ls[1]; z<l->local_lattice[1]+ls[1]; z++)
      for (y=ls[2]; y<l->local_lattice[2]+ls[2]; y++)
        for (x=ls[3]; x<l->local_lattice[3]+ls[3]; x++) {
          
          zeroMat( su3Prod, 3 );
          for (mu=0; mu<4; mu++)
            for (nu=mu+1; nu<4; nu++) {
              zeroMat( tmp1, 3 ); zeroMat( tmp2, 3 );
              addMatMul( tmp1, U[t][z][y][x][mu],
                         U[t+mH[mu][0] ][z+mH[mu][1] ][y+mH[mu][2] ][x+mH[mu][3] ][nu], 3 );
              addMatMulHermit2( tmp2, tmp1,
                                U[t+mH[nu][0] ][z+mH[nu][1] ][y+mH[nu][2] ][x+mH[nu][3] ][mu], 3 );
              addMatMulHermit2( su3Prod, tmp2,
                                U[t][z][y][x][nu], 3 );
            }
#ifndef HAVE_LIME
            plaq += (su3Prod[0] + su3Prod[4] + su3Prod[8])/(k*6.0);
#else
            /* In LIME format:
              *  plaquette = trace / (volume * N_color * number of mu-nu combinations)
              */
            plaq += (su3Prod[0] + su3Prod[4] + su3Prod[8])/(k*3.0*6.0); 
#endif
        }
  
  MPI_Allreduce( &plaq, &averagePlaq, 1, MPI_COMPLEX_double, MPI_SUM, g.comm_cart );
  g.plaq = creal(averagePlaq);
  
  mat_free( &su3Prod, 3 );
  mat_free( &tmp1, 3 );
  mat_free( &tmp2, 3 );

  return creal(averagePlaq);
}


void define_odd_even_table( level_struct *l ) {
  
  int t, z, y, x, mu, i, oe_offset = 0, *le = l->local_lattice, *odd_even_table = g.odd_even_table;
  
  // Determine if whole block is even or odd
  for ( mu=0; mu<4; mu++ )
    oe_offset += (l->local_lattice[mu]*(g.my_coords[mu]/l->comm_offset[mu]))%2;
  oe_offset = oe_offset%2;
  
  // Determine for every site of block
  i = 0;
  for ( t=0; t<le[T]; t++ )
    for ( z=0; z<le[Z]; z++ )
      for ( y=0; y<le[Y]; y++ )
        for ( x=0; x<le[X]; x++ ) {
          odd_even_table[i] = ((t+z+y+x+oe_offset)%2 == 1)?_ODD:_EVEN;
          i++;
        }
}


void m0_update( double m0, level_struct *l, struct Thread *threading ) {

  if (l->depth == 0) {
    m0_update_double( m0, &(g.op_double), l, threading );
    m0_update_float( m0, &(g.op_float), l, threading );
  } else {
    if ( g.mixed_precision )
      m0_update_float( m0, &(l->op_float), l, threading );
    else
      m0_update_double( m0, &(l->op_double), l, threading );
  }
  
  if ( g.mixed_precision ) {
      m0_update_float( m0, &(l->oe_op_float), l, threading );
      m0_update_float( m0, &(l->s_float.op), l, threading );      
  } else {
      m0_update_double( m0, &(l->oe_op_double), l, threading );
      m0_update_double( m0, &(l->s_double.op), l, threading );
  }  

  START_LOCKED_MASTER(threading)
  if(g.print>0) printf0("depth: %d, kappa updated to %f \n", (l->depth), 0.5/(m0 + 4.));
  END_LOCKED_MASTER(threading)
  
  if ( g.interpolation && l->level > 0 && l->next_level != NULL )
    m0_update(m0, l->next_level, threading);
}


void tm_term_update( double mu, level_struct *l, struct Thread *threading ) {
  /***********************
   * tm term data are contained in operator_PRECISION_struct.
   * operator_PRECISION_struct is contained in g, l, gmres_PRECISION_struct & schwarz_PRECISION_struct
   * These structure are in turn contained in g or l
   *   g: op_PRECISION (<- op to be inverted), p or p_MP
   *   l: op_PRECISION, oe_op_PRECISION, s_PRECISION,  p_PRECISION (if g.method < 4)
   * g.op_PRECISION & l->op_PRECISION are constructed at the setup
   * s_PRECISION.op are defined as rearrengement of sites of l->op_PRECISION
   * oe_op_PRECISION is defined as rearrengement of sites of s_PRECISION.op
   * op in p & p_MP in g and p_PRECISION in l points to op in g or op in op_PRECISION, oe_op_PRECISION, s_PRECISION in l
   * ? l->op_PRECISION is actually not used????
   **********************/
  
#ifdef HAVE_TM
  int i;
  double factor = g.mu_factor[l->depth];
  double *even_shift = g.mu_even_shift, odd_shift = g.mu_odd_shift;
  if (l->depth == 0) { // we don't use the multiplicative factor at the top
    tm_term_double_setup( mu, even_shift, odd_shift, 1, &(g.op_double), l, threading ); 
    tm_term_float_setup( mu, even_shift, odd_shift, 1, &(g.op_float), l, threading );
  } else {// are these used after update????
    if ( g.mixed_precision )
      tm_term_float_setup( mu, even_shift, odd_shift, factor, &(l->op_float), l, threading );
    else
      tm_term_double_setup( mu, even_shift, odd_shift, factor, &(l->op_double), l, threading );
  }
  
  if ( g.mixed_precision ) {// if oe_op_PRECISION is not set, setup will be skipped
    tm_term_float_setup( mu, even_shift, odd_shift, factor, &(l->oe_op_float), l, threading );
    tm_term_float_setup( mu, even_shift, odd_shift, factor, &(l->s_float.op), l, threading );
  } else {
    tm_term_double_setup( mu, even_shift, odd_shift, factor, &(l->oe_op_double), l, threading );   
    tm_term_double_setup( mu, even_shift, odd_shift, factor, &(l->s_double.op), l, threading );   
  }

  START_MASTER(threading)
  if(g.print>0)
    for( i=0; i<g.num_rhs_vect; i++) {
      double even = (g.in_setup && i<num_loop)?even_shift[0]:even_shift[i];
      printf0("depth: %d (%dth rhs), mu updated to %f on even sites and %f on odd sites \n",
	      l->depth, i, factor*(mu+even), factor*(mu+odd_shift));
    }
  END_MASTER(threading)

  if ( g.interpolation && l->level > 0 && l->next_level != NULL )
    tm_term_update( mu, l->next_level, threading );
#endif
}

void epsbar_term_update( level_struct *l, struct Thread *threading ) {

#ifdef HAVE_TM1p1
  double factor = g.epsbar_factor[l->depth];
  double epsbar = g.epsbar;
  double even_shift = g.epsbar_ig5_even_shift, odd_shift = g.epsbar_ig5_odd_shift;
    
  if (l->depth == 0) {
    epsbar_term_double_setup( epsbar, even_shift, odd_shift, &(g.op_double), l, threading ); 
    epsbar_term_float_setup( epsbar, even_shift, odd_shift, &(g.op_float), l, threading );
  } else {
    if ( g.mixed_precision )
      epsbar_term_float_setup( factor*epsbar, factor*even_shift, factor*odd_shift, &(l->op_float), l, threading );
    else
      epsbar_term_double_setup( factor*epsbar, factor*even_shift, factor*odd_shift, &(l->op_double), l, threading );
  }
  
  if ( g.mixed_precision ) {
      epsbar_term_float_setup( factor*epsbar, factor*even_shift, factor*odd_shift, &(l->oe_op_float),l, threading );
      epsbar_term_float_setup( factor*epsbar, factor*even_shift, factor*odd_shift, &(l->s_float.op), l, threading );
  } else {
      epsbar_term_double_setup( factor*epsbar, factor*even_shift, factor*odd_shift, &(l->oe_op_double),l, threading );
      epsbar_term_double_setup( factor*epsbar, factor*even_shift, factor*odd_shift, &(l->s_double.op), l, threading );
  }

  START_MASTER(threading)
  if(g.print>0) {
    if( even_shift == odd_shift )
      printf0("depth: %d, epsbar term updated to %f + ig5 %f \n", l->depth, factor*epsbar, factor*even_shift);
    else  
      printf0("depth: %d, epsbar term updated to %f + ig5 %f on even sites and + ig5 %f on odd sites \n", l->depth,
              factor*epsbar, factor*even_shift, factor*odd_shift);
  }
  END_MASTER(threading)

  if ( g.interpolation && l->level > 0 && l->next_level != NULL )
    epsbar_term_update( l->next_level, threading );
#endif
}

void finalize_operator_update( level_struct *l, struct Thread *threading ) {
  // update tm_term or eps_term dependent fiedls
  
  if ( g.odd_even ) {
    if (l->depth == 0) {
      START_LOCKED_MASTER(threading)  
      if(l->s_double.op.clover != NULL )
        schwarz_double_oddeven_setup( &(l->s_double), l );
      if ( l->s_float.op.clover != NULL )
        schwarz_float_oddeven_setup( &(l->s_float), l );
      END_LOCKED_MASTER(threading)
    } else {
      SYNC_CORES(threading)
      if ( g.mixed_precision && !l->idle && (g.method >= 4 || l->level == 0) )
	coarse_oddeven_float_set_self_couplings( l, threading );
      else if ( !l->idle && (g.method >= 4 || l->level == 0) )
	coarse_oddeven_double_set_self_couplings( l, threading );
    }
  }

  if ( g.interpolation && l->level > 0 && l->next_level != NULL )
    finalize_operator_update( l->next_level, threading );
         
#ifdef DEBUG
  if (l->depth == 0 && l->next_level != NULL) 
    test_routine( l, threading );
#endif

}
