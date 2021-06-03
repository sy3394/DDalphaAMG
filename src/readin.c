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
#include "DDalphaAMG.h"


static int read_parameter( void **save_at, char *search_pattern, char *read_format, int number, FILE *read_from, int set_default ) {

/********************************************************************************* 
* Search for the patter and Reads a value for an input parameter from the provided inputfile.
* - void **save_at: points to variable, where parameters are stored.               
* - char *search_pattern: Gives the name of parameter to be read.                      
* - char *read_format: Specifies datatype of parameter.                                
* - int number: Specifies how many numbers need to be read.  
* - FILE *read_from: The inputfile.                         
* - int set_default = _NO_DEFAULT_SET or DEFAULT_SET
*     _NO_DEFAULT_SET: for no default value, i.e., parameters MUST be set in the
*                      inputfile. 
*     DEFAULT_SET:     A default value is assigned to parameters.
*                                                                   
* To see how it is used, see lg_in(...) in the futher below
*********************************************************************************/ 
  
  int i=0, j, k, n=strlen( search_pattern ), match=0;
  char read_pattern[100000], *read_pattern_pt, buffer[50];
  var_table_entry e;
  
  if ( read_from == NULL ) {
    if ( !set_default )
      error0("FILE NULL, unable to find string \"%s\" --- fatal error\n", search_pattern);
    else
      return match;
  }
    
  fseek( read_from, 0L, SEEK_SET );

  // Search through the input file line by line to find search_pattern
  // Note: as it is a while loop, only the first occurance of search_pattern counts
  while ( !match && fgets( read_pattern, 100000, read_from ) ) {

    // Find out if the current line, read_pattern, mathces with earch_pattern
    k = strlen( read_pattern );
    if(k>n) { // Note: read_pattern from fgets contains a new line character
      match = 1;
      i = 0;
      while ( i<n && match ) { 
        if ( search_pattern[i] != read_pattern[i] )
          match = 0;
        i++;
      }
    }
  }

  // deal with initial white spaces
  read_pattern_pt = read_pattern+i;
  while ( *read_pattern_pt == ' ' )
    read_pattern_pt++;

  // if a matched line is found
  if ( match ) {
    if ( strcmp(read_format,"%s") != 0 ) {
      e.pt = *save_at;
      for ( i=0; i<n-1; i++ ) {
        e.name[i] = search_pattern[i];
      }
      e.name[n-1]='\0';
      if ( strcmp(read_format,"%d") == 0 ) {
        // int
        for ( j=0; j<number; j++ ) {
          sscanf( read_pattern_pt, read_format, &(((int*)*save_at)[j]) );//sscanf does not put anyhting to save_at if read pattern is empty or shorter
          sscanf( read_pattern_pt, "%s", buffer );
          read_pattern_pt += strlen( buffer );
          while ( *read_pattern_pt == ' ' )
            read_pattern_pt++;
        }
        sprintf( e.datatype, "int" );
      } else {
        // double
        for ( j=0; j<number; j++ ) {
          sscanf( read_pattern_pt, read_format, &(((double*)*save_at)[j]) );
          sscanf( read_pattern_pt, "%s", buffer );
          read_pattern_pt += strlen( buffer );
          while ( *read_pattern_pt == ' ' )
            read_pattern_pt++;
        }
        sprintf( e.datatype, "double" );
      }
      if ( number == 1 ) {
        var_table_insert( &(g.vt), e );
      }
    } else {
      // string
      sprintf( ((char*)*save_at), "%s", read_pattern_pt );
      ((char*)*save_at)[strlen(read_pattern_pt)-1] = '\0';
    }
  } else {
    if ( !set_default )
      error0("unable to find string \"%s\" --- fatal error\n", search_pattern);
  }

  return match;
}

static void read_global_info( FILE *in ) {

  void *save_pt;
    
  // Note: There is actually no default set for the three following values
  // Though, when using the code as a library, no configuration paths are required.
  save_pt = &(g.in); g.in[0] = '\0';
  read_parameter( &save_pt, "configuration:", "%s", 1, in, _NO_DEFAULT_SET );
  
  save_pt = &(g.in_format); g.in_format = _STANDARD;
  read_parameter( &save_pt, "format:", "%d", 1, in, _DEFAULT_SET );

  // right hand side
  save_pt = &(g.rhs);  g.rhs = 1;
  read_parameter( &save_pt, "right hand side:", "%d", 1, in, _DEFAULT_SET );
  // TODO: support the case below
  if ( g.rhs == 4 ) {
    save_pt = &(g.source_list);
    read_parameter( &save_pt, "source list:", "%s", 1, in, _NO_DEFAULT_SET );
  }
  save_pt = &(g.num_rhs_vect); g.num_rhs_vect=1;
  read_parameter( &save_pt, "number of rhs vectors:", "%d", 1, in, _DEFAULT_SET );

  
  save_pt = &(g.num_levels); g.num_levels = 2;
  read_parameter( &save_pt, "number of levels:", "%d", 1, in, _DEFAULT_SET );
  g.num_desired_levels = g.num_levels;
    
  save_pt = &(g.bc); g.bc = _PERIODIC;
  read_parameter( &save_pt, "boundary conditions:", "%d", 1, in, _DEFAULT_SET );

  if(g.bc==_TWISTED) {
    save_pt = g.twisted_bc; for(int i=0; i<4; i++) g.twisted_bc[i]=0;
    read_parameter( &save_pt, "twisted boundary conditions:", "%d", 4, in, _DEFAULT_SET );
    for(int i=0; i<4; i++) g.twisted_bc[i]*=M_PI;
  }
  
  save_pt = &(g.num_openmp_processes); g.num_openmp_processes = 1;
  read_parameter( &save_pt, "number of openmp threads:", "%d", 1, in, _DEFAULT_SET );

}

// Default values are not available for the following variables.
static void read_no_default_info( FILE *in ) {

  void *save_pt;

  // global lattice
  save_pt = g.global_lattice[0];
  read_parameter( &save_pt, "d0 global lattice:", "%d", 4, in, _NO_DEFAULT_SET );

  // local lattice
  save_pt = g.local_lattice[0];
  read_parameter( &save_pt, "d0 local lattice:", "%d", 4, in, _NO_DEFAULT_SET );

  // block decomposition
  save_pt = g.block_lattice[0];
  read_parameter( &save_pt, "d0 block lattice:", "%d", 4, in, _NO_DEFAULT_SET );

    // Wilson mass
  save_pt = &(g.m0); g.m0 = 0;
  read_parameter( &save_pt, "m0:", "%lf", 1, in, _DEFAULT_SET ); 
  if ( g.m0 == 0 ) {
    double kappa=0; save_pt = &(kappa);    
    read_parameter( &save_pt, "kappa:", "%lf", 1, in, _DEFAULT_SET );
    ASSERT(kappa != 0);
    g.m0 = 1./(2.*kappa)-4.; //setting m0 from kappa
  }
  save_pt = &(g.csw);
  read_parameter( &save_pt, "csw:", "%lf", 1, in, _NO_DEFAULT_SET );
  
#ifdef HAVE_TM
  save_pt = &(g.mu);g.mu=0;
  read_parameter( &save_pt, "mu:", "%lf", 1, in, _DEFAULT_SET );
  if ( g.mu == 0 ) {
    read_parameter( &save_pt, "2KappaMu:", "%lf", 1, in, _DEFAULT_SET );
    g.mu = g.mu*(4.+g.m0);
  }
#endif

#ifdef HAVE_TM1p1
  save_pt = &(g.epsbar); g.epsbar = 0;
  read_parameter( &save_pt, "epsbar:", "%lf", 1, in, _DEFAULT_SET );
  save_pt = &(g.force_2flavours); g.force_2flavours = 0;
  read_parameter( &save_pt, "force_2flavours:", "%d", 1, in, _DEFAULT_SET );
#endif
}

/*  helper functions for read_geometry_data() */
static int shortest_dir( int* data ) {
  int min=0, mu;
  for ( mu=1; mu<4; mu++ )
    if ( data[mu] < data[min] )
      min = mu;
  return min;
}
static int gcd( int a, int b ) {
  if ( b==0 )
    return a;
  return gcd( b, a%b );
}
static int lcm( int a, int b ) {
  return ( a*b / gcd( a, b ) );
}

// ls = max(number of levels,2)
static void read_geometry_data( FILE *in, int ls ) {

  void *save_pt;
  char inputstr[STRINGLENGTH];
  int i, mu, nb, nls, nlls, flag;

  for ( i=0; i<ls; i++ ) {
    
    if(i>0) {
      // global lattice
      sprintf( inputstr, "d%d global lattice:", i );
      save_pt = g.global_lattice[i];
      
      if ( ! read_parameter( &save_pt, inputstr, "%d", 4, in, _DEFAULT_SET ) ) {
        nls = 1;
        for ( mu=0; mu<4; mu++ ) {
          g.global_lattice[i][mu] = g.global_lattice[i-1][mu]/g.block_lattice[i-1][mu];
          nls *= g.global_lattice[i][mu];
        }
        if ( g.odd_even && nls < 2 ) {
          warning0("lattice dimensions not valid for a %d-level method, choosing a %d-level method\n", g.num_levels, i );
          g.num_levels = i; ls = i;
          break;
        }
      }
      
      // local lattice
      sprintf( inputstr, "d%d local lattice:", i );
      save_pt = g.local_lattice[i];
      
      if ( ! read_parameter( &save_pt, inputstr, "%d", 4, in, _DEFAULT_SET ) ) {
        nls = 1;
        nlls = 1;
        for ( mu=0; mu<4; mu++ ) {
          g.local_lattice[i][mu] = g.local_lattice[i-1][mu]/g.block_lattice[i-1][mu];
          nlls *= g.local_lattice[i][mu];
          nls *= g.global_lattice[i][mu];
        }
        if ( g.odd_even && nlls < 2 ) {
          if ( nls/nlls > 1 ) {
            mu = shortest_dir( g.local_lattice[i] );
            if ( g.global_lattice[i][mu] > g.local_lattice[i][mu] ) {
              g.local_lattice[i][mu] *= lcm( g.local_lattice[i][mu],
                                             g.global_lattice[i][mu]/g.local_lattice[i][mu] );
            }
          }
        }
      }
      
      // block lattice
      for ( mu=0; mu<4; mu++ )
        g.block_lattice[i][mu] = 1;
      if ( i<ls-1 ) {
        sprintf( inputstr, "d%d block lattice:", i );
        save_pt = g.block_lattice[i];
        if ( ! read_parameter( &save_pt, inputstr, "%d", 4, in, _DEFAULT_SET ) ) {
          nls = 1;
          nb = 1;
          flag = 1;
          for ( mu=0; mu<4; mu++ )  {
            if ( DIVIDES( 2, g.global_lattice[i][mu] ) ) {
              g.block_lattice[i][mu] = 2;
            } else if ( DIVIDES( 3, g.global_lattice[i][mu] ) ) {
              g.block_lattice[i][mu] = 3;
            } else {
              warning0("lattice dimensions not valid for a %d-level method, choosing a %d-level method\n", g.num_levels, i+1 );
              g.num_levels = i+1; ls=i+1;
              g.block_lattice[i][mu] = 1;
              flag = 0;
              break;
            }
            nb *= g.local_lattice[i][mu]/g.block_lattice[i][mu];
            
            if ( g.local_lattice[i][mu] < g.block_lattice[i][mu] ) {
              g.local_lattice[i][mu] *= g.block_lattice[i][mu];
              if ( ! DIVIDES( g.local_lattice[i][mu], g.global_lattice[i][mu] ) ) {
                g.local_lattice[i][mu] /= g.block_lattice[i][mu];
              }
              warning0("lattice dimensions not valid for a %d-level method, choosing a %d-level method\n", g.num_levels, i+1 );
              g.num_levels = i+1; ls=i+1;
              g.block_lattice[i][mu] = 1;
              flag = 0;
              break;
            }
          }
          
          if ( flag == 1 && g.method == 2 && nb == 1 ) {
            mu = shortest_dir( g.local_lattice[i] );
            if ( g.global_lattice[i][mu] > g.local_lattice[i][mu] ) {
              g.local_lattice[i][mu] *= lcm( g.local_lattice[i][mu],
                                             g.global_lattice[i][mu]/g.local_lattice[i][mu] );
            }
          }
          
        }
      }
    }
#ifdef DEBUG
    printf00("level: %d, gl: %3d %3d %3d %3d\n", i, g.global_lattice[i][0],
             g.global_lattice[i][1],g.global_lattice[i][2],g.global_lattice[i][3] );
    
    printf00("level: %d, ll: %3d %3d %3d %3d\n", i, g.local_lattice[i][0],
             g.local_lattice[i][1],g.local_lattice[i][2],g.local_lattice[i][3] );
    
    printf00("level: %d, bl: %3d %3d %3d %3d\n\n", i, g.block_lattice[i][0],
             g.block_lattice[i][1],g.block_lattice[i][2],g.block_lattice[i][3] );
#endif
    
    
    sprintf( inputstr, "d%d post smooth iter:", i );
    save_pt = &(g.post_smooth_iter[i]); g.post_smooth_iter[i] = 4;
    read_parameter( &save_pt, inputstr, "%d", 1, in, _DEFAULT_SET );
    
    sprintf( inputstr, "d%d preconditioner cycles:", i );
    save_pt = &(g.ncycle[i]); g.ncycle[i] = 1;
    read_parameter( &save_pt, inputstr, "%d", 1, in, _DEFAULT_SET );
    
    
    sprintf( inputstr, "d%d relaxation factor:", i );
    save_pt = &(g.relax_fac[i]); g.relax_fac[i] = 1.0;
    read_parameter( &save_pt, inputstr, "%lf", 1, in, _DEFAULT_SET );
    
    sprintf( inputstr, "d%d block iter:", i );
    save_pt = &(g.block_iter[i]);
    g.block_iter[i] = 4;
    read_parameter( &save_pt, inputstr, "%d", 1, in, _DEFAULT_SET );
    
    sprintf( inputstr, "d%d setup iter:", i );
    save_pt = &(g.setup_iter[i]);
    if ( i==0 ) g.setup_iter[i] = 5;
    else if ( i==1 ) g.setup_iter[i] = 3;
    else if ( i>1 ) g.setup_iter[i] = 2;
    read_parameter( &save_pt, inputstr, "%d", 1, in, _DEFAULT_SET );
    
#ifdef HAVE_TM
    sprintf( inputstr, "d%d mu factor:", i );
    save_pt = &(g.mu_factor[i]); g.mu_factor[i] = 1;
    read_parameter( &save_pt, inputstr, "%lf", 1, in, _DEFAULT_SET );
#endif

#ifdef HAVE_TM1p1
    sprintf( inputstr, "d%d epsbar factor:", i );
    save_pt = &(g.epsbar_factor[i]); g.epsbar_factor[i] = 1;
    read_parameter( &save_pt, inputstr, "%lf", 1, in, _DEFAULT_SET );
#endif
    
    sprintf( inputstr, "d%d test vectors:", i );
    save_pt = &(g.num_eig_vect[i]);
    if ( i==0 ) g.num_eig_vect[i] = 24;
    else g.num_eig_vect[i] = 28;
    read_parameter( &save_pt, inputstr, "%d", 1, in, _DEFAULT_SET );

    sprintf( inputstr, "d%d solver:", i );
    save_pt = &(g.solver[i]); g.solver[i] = _FGMRES;
    read_parameter( &save_pt, inputstr, "%d", 1, in, _DEFAULT_SET );

#ifdef HAVE_FABULOUS
    sprintf( inputstr, "d%d fabulous orthogonalization scheme:", i );
    save_pt = &(g.f_orthoscheme[i]); g.f_orthoscheme[i] = FABULOUS_MGS;
    read_parameter( &save_pt, inputstr, "%d", 1, in, _DEFAULT_SET );

    sprintf( inputstr, "d%d fabulous orthogonalization type:", i );
    save_pt = &(g.f_orthotype[i]); g.f_orthotype[i] = FABULOUS_RUHE;
    read_parameter( &save_pt, inputstr, "%d", 1, in, _DEFAULT_SET );

    sprintf( inputstr, "d%d fabulous orthogonalization iter:", i );
    save_pt = &(g.ortho_iter[i]); g.ortho_iter[i] = 2;
    read_parameter( &save_pt, inputstr, "%d", 1, in, _DEFAULT_SET );

    sprintf( inputstr, "d%d fabulous max kept dir:", i );
    save_pt = &(g.max_kept_direction[i]); g.max_kept_direction[i] = -1;
    read_parameter( &save_pt, inputstr, "%d", 1, in, _DEFAULT_SET );

    sprintf( inputstr, "d%d number of deflating eigenvectors for fabulous:", i );
    save_pt = &(g.k[i]); g.k[i] = 0;
    read_parameter( &save_pt, inputstr, "%d", 1, in, _DEFAULT_SET );

    sprintf( inputstr, "d%d max number of mat-vec product for fabulous:", i );
    save_pt = &(g.max_mvp[i]); g.max_mvp[i] = 0;
    read_parameter( &save_pt, inputstr, "%d", 1, in, _DEFAULT_SET );
    
    sprintf( inputstr, "d%d fabulous compute real residual:", i );
    save_pt = &(g.real_residual[i]); g.real_residual[i] = 0;
    read_parameter( &save_pt, inputstr, "%d", 1, in, _DEFAULT_SET );
#endif
  }
}

static void read_solver_parameters( FILE *in ) {

  void *save_pt;

  save_pt = &(g.mixed_precision); g.mixed_precision = 2;
  read_parameter( &save_pt, "mixed precision:", "%d", 1, in, _DEFAULT_SET );
  if ( g.num_levels == 1 ) g.interpolation = 0; else {
    save_pt = &(g.interpolation); g.interpolation = 2;
    read_parameter( &save_pt, "interpolation:", "%d", 1, in, _DEFAULT_SET );
  }

  int db = g.num_levels-1;
  save_pt = &(g.randomize); g.randomize = 0;
  read_parameter( &save_pt, "randomize test vectors:", "%d", 1, in, _DEFAULT_SET );
  save_pt = &(g.max_iter[db]); g.max_iter[db] = 200;
  read_parameter( &save_pt, "coarse grid iterations:", "%d", 1, in, _DEFAULT_SET );
  save_pt = &(g.max_restart[db]); g.max_restart[db] = 10;
  read_parameter( &save_pt, "coarse grid restarts:", "%d", 1, in, _DEFAULT_SET );
  save_pt = &(g.tol[db]); g.tol[db] = 1E-1;
  read_parameter( &save_pt, "coarse grid tolerance:", "%le", 1, in, _DEFAULT_SET );
  save_pt = &(g.setup_tol); g.setup_tol = g.tol[db];
  read_parameter( &save_pt, "adaptive setup tolerance:", "%le", 1, in, _DEFAULT_SET );
  save_pt = &(g.odd_even); g.odd_even = 1;
  read_parameter( &save_pt, "odd even preconditioning:", "%d", 1, in, _DEFAULT_SET );

  save_pt = &(g.setup_m0); g.setup_m0 = g.m0;
  read_parameter( &save_pt, "setup m0:", "%lf", 1, in, _DEFAULT_SET );
#ifdef HAVE_TM
  //TODO: multi-mass shifts solveron odd-sites is not implemented
  save_pt = &(g.mu_odd_shift); g.mu_odd_shift = 0.;
  read_parameter( &save_pt, "mu odd shift:", "%lf", 1, in, _DEFAULT_SET );
  save_pt = g.mu_even_shift; for ( int j=0; j<g.num_rhs_vect; j++ ) g.mu_even_shift[j] = 0.;
  read_parameter( &save_pt, "mu even shift:", "%lf", g.num_rhs_vect, in, _DEFAULT_SET );
  save_pt = &(g.setup_mu); g.setup_mu = g.mu;
  read_parameter( &save_pt, "setup mu:", "%lf", 1, in, _DEFAULT_SET );
#endif

#ifdef HAVE_TM1p1
  save_pt = &(g.epsbar_ig5_odd_shift);g.epsbar_ig5_odd_shift=0;
  read_parameter( &save_pt, "epsbar odd shift:", "%lf", 1, in, _DEFAULT_SET );
  save_pt = &(g.epsbar_ig5_even_shift);g.epsbar_ig5_even_shift=0;
  read_parameter( &save_pt, "epsbar even shift:", "%lf", 1, in, _DEFAULT_SET );
#endif
  
  save_pt = &(g.method); g.method = 2;
  read_parameter( &save_pt, "method:", "%d", 1, in, _DEFAULT_SET );
  save_pt = &(g.max_iter[0]); g.max_iter[0] = 30;
  read_parameter( &save_pt, "iterations between restarts:", "%d", 1, in, _DEFAULT_SET );
  save_pt = &(g.max_restart[0]); g.max_restart[0] = 20;
  read_parameter( &save_pt, "maximum of restarts:", "%d", 1, in, _DEFAULT_SET );
  save_pt = &(g.tol[0]); g.tol[0] = 1E-10;
  read_parameter( &save_pt, "tolerance for relative residual:", "%le", 1, in, _DEFAULT_SET );
  save_pt = &(g.print); g.print = 0;
  read_parameter( &save_pt, "print mode:", "%d", 1, in, _DEFAULT_SET );
#ifdef HAVE_TM
  save_pt = &(g.downprop); g.downprop=1;
  read_parameter( &save_pt, "addDownPropagator:", "%d", 1, in, _DEFAULT_SET );
#endif

  save_pt = &(g.use_only_fgrmes_at_setup); g.use_only_fgrmes_at_setup = 0;
  read_parameter( &save_pt, "use only FGMRES at setup:", "%d", 1, in, _DEFAULT_SET );
#ifdef HAVE_FABULOUS
  save_pt = &(g.n_defl_updates); g.n_defl_updates = 5;
  read_parameter( &save_pt, "fabulous num of deflation space updates:", "%d", 1, in, _DEFAULT_SET );
  save_pt = &(g.logger_user_data_size); g.logger_user_data_size = 0;
  read_parameter( &save_pt, "fabulous user data size for log:", "%d", 1, in, _DEFAULT_SET );
  save_pt = &(g.quiet); g.quiet = 1;
  read_parameter( &save_pt, "fabulous silent run:", "%d", 1, in, _DEFAULT_SET );
#endif
  if ( g.randomize ) {
    srand( time( 0 ) + 1000*g.my_rank );
  } else 
    srand( 1000*g.my_rank );
}

static void read_testvector_io_data_if_necessary( FILE *in ) {

  void *save_pt;
  if ( g.interpolation == 4 ) {
    save_pt = &(g.tv_io_single_file); g.tv_io_single_file = 1;
    read_parameter( &save_pt, "test vector io from single file:", "%d", 1, in, _DEFAULT_SET );
    save_pt = &(g.tv_io_file_name);
    read_parameter( &save_pt, "test vector io file name:", "%s", 1, in, _NO_DEFAULT_SET );
  }
}

static void read_evaluation_parameters_if_necessary( FILE *in ) {

  void *save_pt;
  save_pt = &(g.vt.evaluation); g.vt.evaluation = 0;
  read_parameter( &save_pt, "evaluation:", "%d", 1, in, _DEFAULT_SET );
  if ( g.vt.evaluation ) {
    save_pt = &(g.vt.scan_var);
    read_parameter( &save_pt, "scan variable:", "%s", 1, in, _NO_DEFAULT_SET );  
    save_pt = &(g.vt.start_val);
    read_parameter( &save_pt, "start value:", "%lf", 1, in, _NO_DEFAULT_SET );
    save_pt = &(g.vt.end_val);
    read_parameter( &save_pt, "end value:", "%lf", 1, in, _NO_DEFAULT_SET );
    save_pt = &(g.vt.step_size);
    read_parameter( &save_pt, "step size:", "%lf", 1, in, _NO_DEFAULT_SET );
    save_pt = &(g.vt.multiplicative);
    read_parameter( &save_pt, "multiplicative:", "%d", 1, in, _NO_DEFAULT_SET );
    save_pt = &(g.vt.shift_update);
    read_parameter( &save_pt, "shift update:", "%d", 1, in, _NO_DEFAULT_SET );
    save_pt = &(g.vt.re_setup);
    read_parameter( &save_pt, "setup update:", "%d", 1, in, _NO_DEFAULT_SET );
    save_pt = &(g.vt.track_error); g.vt.track_error = 0;
    read_parameter( &save_pt, "track error:", "%d", 1, in, _NO_DEFAULT_SET );
    save_pt = &(g.vt.track_cgn_error); g.vt.track_cgn_error = 0;
    read_parameter( &save_pt, "compare with CGN error:", "%d", 1, in, _DEFAULT_SET );
    save_pt = &(g.vt.average_over); g.vt.average_over = 1;
    read_parameter( &save_pt, "average over:", "%d", 1, in, _DEFAULT_SET );
  }
}

static void read_kcycle_data( FILE *in ) {

  void *save_pt;
  save_pt = &(g.kcycle); g.kcycle = 1;
  read_parameter( &save_pt, "kcycle:", "%d", 1, in, _DEFAULT_SET );
  save_pt = &(g.kcycle_restart); g.kcycle_restart = 5;
  read_parameter( &save_pt, "kcycle length:", "%d", 1, in, _DEFAULT_SET );
  save_pt = &(g.kcycle_max_restart); g.kcycle_max_restart = 2;
  read_parameter( &save_pt, "kcycle restarts:", "%d", 1, in, _DEFAULT_SET );
  save_pt = &(g.kcycle_tol); g.kcycle_tol = 1E-1;
  read_parameter( &save_pt, "kcycle tolerance:", "%le", 1, in, _DEFAULT_SET );

  for ( int i=1; i<g.num_levels-1; i++ ) {
    g.max_iter[i] = g.kcycle_restart;
    g.max_restart[i] = g.kcycle_max_restart;
    g.tol[i] = g.kcycle_tol;
  }
}

static void validate_parameters( int ls, level_struct *l ) {

  int i;
  int mu;

  if ( g.method > 3 ) {
    if( g.interpolation != 0 ) {
      warning0("Multigrid with GMRES/BiCGstab smoothing is not supported.\n         Switching to FGMRES preconditioned with GMRES/BiCGstab (g.interpolation=0).\n");
      g.interpolation = 0;
    }
    if ( g.mixed_precision == 2 ) {
      warning0("FGMRES+GMRES or BiCGstab with MP is not currently supported.\n         Switching to double precision\n");
      g.mixed_precision = 0;
    }
    if ( g.odd_even ) {//to be removed
    warning0("BiCGstab w/t even-odd prec. is not currently supported.\n         Switching to even odd prec..\n");
    g.odd_even = 1;
    }
  }

  ASSERT( ASCENDING( 0, g.rhs, 2 ) );
  ASSERT( ASCENDING( -1, g.method, 5 ) );
  if ( g.method < 1 ) {
    warning0("Multigrid is not supported.\n         Switching to the chosen method with no AMG (g.interpolation=0).\n");
    g.interpolation = 0;
    //ASSERT ( g.interpolation == 0 );
  }
  
  if ( g.method == 0 && g.mixed_precision == 1 ) {
    warning0("Pure GMRES uses either mixed precision solver or double precision solver.\n         Switching to doule precision\n");
    g.mixed_precision = 0;
  }

  int solver = 0;
  for ( i=0; i<g.num_levels; i++ ) solver += g.solver[i];
#ifndef HAVE_FABULOUS
  if ( solver != 0 )
    error0("Non-zero solver options indicate fabulous solver, which is switched off by the user\n");
  if ( g.use_only_fgrmes_at_setup != 0 ) {
    warning0("As fabulous solvers are not available in the current parameter choice, we need to use FGMRES at the setup\n         Switching to FGMRES\n");
    g.use_only_fgrmes_at_setup = 0;
  }
#else
  if ( g.method > 0 && g.method < 3 ) {
    if ( g.mixed_precision == 2 && g.solver[0] ) {
      warning0("Fabulous solvers do not support mixed precision.\n         Switching to single precision.\n");
      g.mixed_precision = 1;
    }
    for ( i=0; i<g.num_levels; i++ ) {
      if ( !(g.solver[i] < _NUM_SOLVER) ) {
	warning0("The solver must be smaller than %d (depth==%d).\n         Switching to fabulous BGCR method.\n", _NUM_SOLVER, i);
	g.solver[i] = _GCR;
      }
      if ( i < g.num_levels-1  && g.solver[i] != _GCR && g.solver[i] != _FGMRES ) {
	warning0("The solver at the depth %d needs to be flexible to employ AMG as its right preconditioner.\n         Switching to fabulous BGCR method.\n", i);
	g.solver[i] = _GCR;
      }
      if ( g.solver[i] == _GCR && g.f_orthotype[i] != FABULOUS_BLOCK ) {
	warning0("Only BLOCK-wise orthogonalization is currently implemented for BCGR. The BLOCK-wise version will be used at depth %d\n", i);
	g.f_orthotype[i] = FABULOUS_BLOCK;
      }
      if ( g.k[i] < 0 ) {
	warning0("Number of deflating eigenvectors for FABULOUS solvers (depth %d) must be non-negative.\n         Setting it to 0.\n", i);
	g.k[i] = 0;
      }
      int nvec = num_loop;
#ifdef HAVE_TM1p1
      if ( g.epsbar == 0 && g.epsbar_ig5_odd_shift == 0 && g.epsbar_ig5_odd_shift == 0 )
	nvec *= 2;
#endif
      if ( g.solver[i] == _GCRO && g.k[i] > 2*g.max_iter[i]-nvec ) {
	warning0("Deflation space size should be smaller than %d.\n         Setting it to %d.\n", 2*g.max_iter[i]-nvec, 2*g.max_iter[i]-nvec);
	g.k[i] = 2*g.max_iter[i]-nvec;
      }
    }
  } else if ( solver != 0 ) {
    // TODO: for g.method > 3, we can write fabulous version
    error0("Fabulous solvers are not supported when g.method != 0, 1, 2, 3\n");
  }
#endif

  ASSERT( IMPLIES( g.vt.evaluation, g.rhs <= 2 ) );

  ASSERT( DIVIDES( num_loop, g.num_rhs_vect ) );
  for ( i=0; i<g.num_levels-1; i++ )
    ASSERT( DIVIDES( num_loop, g.num_eig_vect[i] ) );

  for ( i=0; i<g.num_levels; i++ )
    for ( mu=0; mu<4; mu++)
      ASSERT( DIVIDES( g.local_lattice[i][mu], g.global_lattice[i][mu] ) );
  
  for ( i=0; i<g.num_levels-1; i++ )
    for ( mu=0; mu<4; mu++) {
      ASSERT( DIVIDES( g.global_lattice[i+1][mu], g.global_lattice[i][mu] ) );
      ASSERT( DIVIDES( g.block_lattice[i][mu], g.local_lattice[i][mu] ) );
      ASSERT( DIVIDES( g.global_lattice[i][mu]/g.global_lattice[i+1][mu], g.local_lattice[i][mu] ) ); 
      ASSERT( DIVIDES( g.block_lattice[i][mu], g.global_lattice[i][mu]/g.global_lattice[i+1][mu] ) );
    }
    
  if ( g.odd_even ) {
    int coarse_sites_per_core = 1;
    for ( mu=0; mu<4; mu++ ) {
      ASSERT( DIVIDES( 2, g.global_lattice[g.num_levels-1][mu] ) );
      coarse_sites_per_core *= g.local_lattice[g.num_levels-1][mu];
    }
    ASSERT( DIVIDES( 2, coarse_sites_per_core ) );
  }

  if ( g.method == 2 ) {
    for ( i=0; i<ls-1; i++ ) {
      int num_blocks = 1;
      for ( mu=0; mu<4; mu++) {
        num_blocks *= ( g.local_lattice[i][mu]/g.block_lattice[i][mu] );
      }
      ASSERT( num_blocks >= 2 );
    }
  }
  /* TODO: perhaps support this
  for ( mu=0; mu<4; mu++ )
    ASSERT( IMPLIES( g.rhs == 3, ASCENDING( 0, g.propagator_coords[mu], l->global_lattice[mu]-1 ) ) );
  */
  
  ASSERT( IMPLIES( g.method > 0 && g.interpolation > 0, g.max_iter[g.num_levels-1] > 0 ) );
  ASSERT( IMPLIES( g.method > 0 && g.interpolation > 0, g.max_restart[g.num_levels-1] > 0 ) );
  ASSERT( IMPLIES( g.method > 0 && g.interpolation > 0, 0 < g.tol[g.num_levels-1] && g.tol[g.num_levels-1] < 1 ) );
  ASSERT( IMPLIES( g.method > 0, l->n_cy > 0 ) );
  ASSERT( g.max_restart[0] > 0 );
  ASSERT( 0 < g.tol[0] && g.tol[0] < 1 );
  ASSERT( ASCENDING( 0, g.kcycle, 1 ) );
  ASSERT( IMPLIES( g.kcycle && g.method > 0, g.kcycle_restart > 0 ) );
  ASSERT( IMPLIES( g.kcycle && g.method > 0, g.kcycle_max_restart > 0 ) );
  ASSERT( IMPLIES( g.kcycle && g.method > 0, 0 < g.kcycle_tol && g.kcycle_tol < 1 ) );
  
#ifdef HAVE_TM1p1
  if ( g.epsbar != 0 || g.epsbar_ig5_odd_shift != 0 || g.epsbar_ig5_odd_shift != 0 ) {
    warning0("force_2flavours is set to 0 when eps term is non-zero to aovid confusion\n");
    g.force_2flavours = 0;
#ifdef HAVE_FABULOUS
    warning0("With non-zero eps term, fabulous is used only during the solver phase\n         Switching use_only_fgrmes_at_setup to true\n");
    g.use_only_fgrmes_at_setup = 1;
#endif
  }
#endif
  //LIST OF CASES WHICH SHOULD WORK, BUT DO NOT (TODO)

  //TODO: Could work without, but you need to fix the setup phase.    
  for ( i=0; i<g.num_levels-2; i++ )
    ASSERT( g.num_eig_vect[i] <= g.num_eig_vect[i+1] );

  //TODO: for some reason g.mixed_precision=0 do not work with g.num_levels>2
  if ( g.num_levels>2 && g.interpolation )
    ASSERT( g.mixed_precision );

}

static void allocate_for_global_struct_after_read_global_info( int ls ) {

  int i;
  MALLOC( g.global_lattice, int*, ls );
  MALLOC( g.local_lattice, int*, ls );
  MALLOC( g.block_lattice, int*, ls );
  g.global_lattice[0] = NULL;
  g.local_lattice[0] = NULL;
  g.block_lattice[0] = NULL;
  MALLOC( g.global_lattice[0], int, 4*ls );
  MALLOC( g.local_lattice[0], int, 4*ls );
  MALLOC( g.block_lattice[0], int, 4*ls );
  MALLOC( g.post_smooth_iter, int, ls );
  MALLOC( g.ncycle, int, ls );
  MALLOC( g.relax_fac, double, ls );
#ifdef HAVE_TM
  MALLOC( g.mu_factor, double, ls );
  MALLOC( g.mu_even_shift, double, g.num_rhs_vect );
#endif
#ifdef HAVE_TM1p1
  MALLOC( g.epsbar_factor, double, ls );
#endif
  MALLOC( g.max_iter, int, ls );
  MALLOC( g.max_restart, int, ls );
  MALLOC( g.tol, double, ls );
  MALLOC( g.iter_times, double, ls );
  MALLOC( g.iter_counts, int, ls );
  MALLOC( g.block_iter, int, ls );
  MALLOC( g.setup_iter, int, ls );
  MALLOC( g.num_eig_vect, int, ls );
  MALLOC( g.solver, int, ls );
#ifdef HAVE_FABULOUS
  MALLOC( g.f_orthoscheme, fabulous_orthoscheme, ls );
  MALLOC( g.f_orthotype, fabulous_orthotype, ls );
  MALLOC( g.ortho_iter, int, ls );
  MALLOC( g.max_kept_direction, int, ls );
  MALLOC( g.k, int, ls );
  MALLOC( g.max_mvp, int, ls );
  MALLOC( g.real_residual, int, ls );
#endif
  for ( i=1; i<ls; i++ ) {
    g.global_lattice[i] = g.global_lattice[0] + i*4;
    g.local_lattice[i] = g.local_lattice[0] + i*4;
    g.block_lattice[i] = g.block_lattice[0] + i*4;
  }
}

// This is for the level structure associated with the top level
static void set_level_and_global_structs_according_to_global_struct( level_struct *l ) {

  int mu;
  
  l->level = g.num_levels-1; l->depth = 0; l->idle = 0;
  l->global_lattice = g.global_lattice[0];
  l->local_lattice = g.local_lattice[0];
  l->block_lattice = g.block_lattice[0];
  l->post_smooth_iter = g.post_smooth_iter[0];
  l->n_cy = g.ncycle[0];
  l->relax_fac = g.relax_fac[0];
  l->block_iter = g.block_iter[0];
  l->setup_iter = g.setup_iter[0];
  l->num_eig_vect = g.num_eig_vect[0];
  l->num_parent_eig_vect = 6; //for consistency sake
    
  // compute some additional values
  l->num_lattice_site_var = 12;
  g.num_processes = 1;
  for ( mu=0; mu<4; mu++ ) {
    l->comm_offset[mu] = 1;
    l->coarsening[mu] = g.global_lattice[0][mu]/MAX(1,g.global_lattice[1][mu]);
    l->global_splitting[mu] = l->global_lattice[mu]/l->local_lattice[mu];
    g.process_grid[mu] = l->global_splitting[mu];
    l->periodic_bc[mu] = 1;
    g.num_processes *= l->global_splitting[mu];
  }
  
  g.setup_m0 = g.m0;
}

/*************        Public Functions    **************/
// This is used in method_init
void lg_in( char *inputfile, level_struct *l ) {

  FILE *in;

  ASSERT( (in = fopen( inputfile, "r" )) != NULL );

  read_global_info( in );
  
  int ls = MAX(g.num_levels,2);
  allocate_for_global_struct_after_read_global_info( ls );

  read_no_default_info( in ); 
  read_solver_parameters( in ); 
  read_geometry_data( in, ls );

  ls = MAX(g.num_levels,2); // update ls; could be changed in read_geometry_data

  set_level_and_global_structs_according_to_global_struct( l );

  read_testvector_io_data_if_necessary( in );
  read_evaluation_parameters_if_necessary( in );
  read_kcycle_data( in );

#ifdef HAVE_TM
  // sort the even shifts to ascending order in magnitude
  double mu;
  for ( int i=0; i<g.num_rhs_vect; i++ ) 
    for ( int j=i+1; j<g.num_rhs_vect; j++ )
      if (  g.mu_even_shift[i]*g.mu_even_shift[i] >  g.mu_even_shift[j]*g.mu_even_shift[j] ) {
	mu = g.mu_even_shift[i];
	g.mu_even_shift[i] = g.mu_even_shift[j];
	g.mu_even_shift[j] = mu;
      }
  // compute meta info for even shifts
  for ( int i=0; i<g.num_rhs_vect; i++ ) {
    g.even_shifted += g.mu_even_shift[i]?1:0;
    g.is_even_shifted_mu_nonzero += (g.mu + g.mu_even_shift[i]!=0.0)?1:0;
  }
#endif
  validate_parameters( ls, l );

  printf00("configuration: %s\n", g.in );
  if( g.rhs == 4 )
    printf00("source list: %s\n", g.source_list );
  fclose(in);
}

/*****************  Interface  ***************************/

// Helper function for set_DDalphaAMG_parameters
static void set_global_info( struct init *params ) {

  // global lattice
  for( int i=0; i<4; i++ ) {
    g.global_lattice[0][i] = params->global_lattice[i];
    g.local_lattice[0][i] = params->global_lattice[i]/params->procs[i];
    g.block_lattice[0][i] = params->block_lattice[i];
  }

  g.bc = params->bc;

  if(g.bc==_TWISTED) {
    for(int i=0; i<4; i++){
      g.twisted_bc[i]=params->theta[i];
      g.twisted_bc[i]*=M_PI;
    }
  }

  // rhs
  g.num_rhs_vect = params->nrhs;
  
  // Operator
  g.m0 = 1./(2.*params->kappa)-4.;
  g.csw = params->csw;
#ifdef HAVE_TM
  g.mu = params->mu;
#endif
  
  g.num_openmp_processes = params->number_openmp_threads;

}

// For interface; replaces lg_in
void set_DDalphaAMG_parameters( struct init *params, level_struct *l ) {

  FILE *in=NULL;

  if (params->init_file != NULL) 
    ASSERT( (in = fopen( params->init_file, "r" )) != NULL );
  
  g.num_levels = params->number_of_levels;
  g.num_desired_levels = g.num_levels;

  int ls = MAX(g.num_levels,2);
  allocate_for_global_struct_after_read_global_info( ls );

  set_global_info( params );
  read_solver_parameters( in );
  read_geometry_data( in, ls );
  
  ls = MAX(g.num_levels,2); // update ls; could be changed in read_geometry_data
  
  set_level_and_global_structs_according_to_global_struct( l );

  read_testvector_io_data_if_necessary( in );
  read_evaluation_parameters_if_necessary( in );
  read_kcycle_data( in );

#ifdef HAVE_TM
  // sort the even shifts to ascending order in magnitude 
  double mu;
  for ( int i=0; i<g.num_rhs_vect; i++ )
    for ( int j=i+1; j<g.num_rhs_vect; j++ )
      if (  g.mu_even_shift[i]*g.mu_even_shift[i] >  g.mu_even_shift[j]*g.mu_even_shift[j] ) {
        mu = g.mu_even_shift[i];
        g.mu_even_shift[i] = g.mu_even_shift[j];
        g.mu_even_shift[j] = mu;
      }
  // compute meta info for even shifts
  for ( int i=0; i<g.num_rhs_vect; i++ ) {
    g.even_shifted += g.mu_even_shift[i]?1:0;
    g.is_even_shifted_mu_nonzero += (g.mu + g.mu_even_shift[i]!=0.0)?1:0;
  }
#endif
#ifdef HAVE_FABULOUS
  // For interface solver, we use FABULOUS only during the inversion
  g.use_only_fgrmes_at_setup = 1;
#endif
  
  validate_parameters( ls, l );

  if (params->init_file != NULL) 
    fclose(in);

}

// this function is used in var_table.h too
void parameter_update( level_struct *l ) {
  
  if(l->depth==0) {
    int ls = MAX(g.num_levels,2);
    set_level_and_global_structs_according_to_global_struct( l );
    validate_parameters( ls, l );
  }

  l->level = g.num_levels-1-l->depth;
  l->post_smooth_iter = g.post_smooth_iter[l->depth];
  l->block_iter = g.block_iter[l->depth];
  l->setup_iter = g.setup_iter[l->depth];
  l->num_eig_vect = g.num_eig_vect[l->depth];
  if(l->depth>0)
    l->num_parent_eig_vect = g.num_eig_vect[l->depth-1];
  else
    l->num_parent_eig_vect = 6;
  
  if ( l->level > 0 && l->next_level != NULL ) 
    parameter_update( l->next_level );
}
