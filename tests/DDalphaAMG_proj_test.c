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

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include <mpi.h>
#include <omp.h>

//library to include
#include "DDalphaAMG.h"

enum { T, Z, Y, X };
#define printf0(...) do{if(rank==0)printf(__VA_ARGS__);}while(0)

int rank, max_setup=0, conf_format;
float residual=1e-9;
MPI_Comm comm_cart;
DDalphaAMG_init init;
DDalphaAMG_parameters params;
DDalphaAMG_status status;
char * conf_file = "../conf/8x8x8x8b6.0000id3n1";
char * options = "c:f:i:L:p:B:n:l:k:w:u:t:hr:C:K:V:m:s:S:v:D:b:g:I:M:e";
/*
 * Setting standard values for DDalphaAMG_init
 */
void standard_init() {
  init.global_lattice[T] = 8;
  init.global_lattice[Z] = 8;
  init.global_lattice[Y] = 8;
  init.global_lattice[X] = 8;

  init.nrhs = BASE_LOOP_COUNT;
  
  init.kappa = 0.142857143;
  init.mu = 0.1;
  init.csw = 1;

  init.bc = 1;

  init.number_of_levels = 3;

  init.procs[T] = 2;
  init.procs[Z] = 1;
  init.procs[Y] = 1;
  init.procs[X] = 1;

  init.number_openmp_threads = 2;

  init.block_lattice[T]=2;
  init.block_lattice[Z]=2;
  init.block_lattice[Y]=2;
  init.block_lattice[X]=2;

  init.comm_cart = MPI_COMM_WORLD;
  init.Cart_coords = NULL;
  init.Cart_rank = NULL;
  init.init_file = NULL;
}

/*
 * Printing implemented parameters
 */
void help( char * arg0 ) {
  static int printed = 0;
  if(!printed) {
    printf0("\n\n");
    printf0("Usage: %s -c <conf> [<option(s)>]\n", arg0);
    printf0("   -c PATH      Configuration to load\n");
    printf0("   -f #         Configuration format (0 -> DDalphaAMG, 1 -> Lime)\n");
    printf0("   -i PATH      Input file (optional)\n");
    printf0("   -L T Y X Z   Lattice size in each direction\n");
    printf0("   -p T Y X Z   Processors in each direction\n");
    printf0("   -B T Y X Z   Block size in each direction on first level.\n");
    printf0("   -k #         kappa for the configuration\n");
    printf0("   -w #         c_sw for the configuration\n");
    printf0("   -u #         mu for the configuration\n");
    printf0("   -t #         Number of OpenMp threads\n");
    printf0("   -r #         Relative residual\n");
    printf0("   -K #         K-cycle tolerance\n");
    printf0("   -C #         Tolerance on coarsest grid\n");
    printf0("   -l #         Number of levels, l (from 1 to 4)\n");
    printf0("   -V 1 [2] [3] Basis vectors between each level (l-1)\n");
    printf0("   -m 2 [3] [4] Factor for mu on coarse levels\n");
    printf0("   -s 1 [2] [3] Setup iterations on each level (l-1)\n");
    printf0("   -b 1 [2] [3] [4] Orthogonalization type for a FABULOUS solver on each level\n");
    printf0("   -g 1 [2] [3] [4] Orthogonalization scheme for a FABULOUS solver on each level\n");
    printf0("   -I 1 [2] [3] [4] Number of iteration for iterative variant of Gram-Schmidt procedure for a FABULOUS solver on each level\n");
    printf0("   -M 1 [2] [3] [4] Max kept direction for some Fabulous solvers  on each level\n");
    printf0("   -e 1 [2] [3] [4] Number of deflating eigenvectors for some Fabulous solvers on each level\n");
    printf0("   -v           Verbose\n");
  }
  printf0("\n\n");
  printed++;
} 

/*
 * Printing implemented parameters
 */
void read_init_arg(int argc, char *argv[] ) {

  int opt, p=0, fail=0, mu;
  optind = 0;
  while ((opt = getopt(argc, argv, options)) != -1) {
    switch (opt) {
    case 'c':
      conf_file = optarg;
      break;
    case 'f':
      conf_format = atoi(optarg);
      break;
    case 'i':
      init.init_file = optarg;
      break;
    case 'L':
      optind--;
      mu=0;
      for ( ; optind < argc && *argv[optind] != '-'; optind++){
	if(mu > 3) {
	  printf0("Error: too many arguments in -L.\n");
	  p++;
	  fail++;
	  break;
	}
	init.global_lattice[mu] = atoi(argv[optind]);
	mu++;        
      }
      if(mu < 4) {
	printf0("Warning: too few arguments in -L.\n");
	p++;
	fail++;
      }
      break;
    case 'p':
      optind--;
      mu=0;
      for ( ; optind < argc && *argv[optind] != '-'; optind++){
	if(mu > 3) {
	  printf0("Error: too many arguments in -p.\n");
	  p++;
	  fail++;
	  break;
	}
	init.procs[mu] = atoi(argv[optind]);
	mu++;        
      }
      if(mu < 4) {
	printf0("Warning: too few arguments in -p.\n");
	p++;
      }
      break;
    case 'B':
      optind--;
      mu=0;
      for ( ; optind < argc && *argv[optind] != '-'; optind++){
	if(mu > 3) {
	  printf0("Error: too many arguments in -B.\n");
	  p++;
	  fail++;
	  break;
	}
	init.block_lattice[mu] = atoi(argv[optind]);
	mu++;        
      }
      if(mu < 4) {
	printf0("Warning: too few arguments in -B.\n");
	p++;
      }
      break;
    case 'n':
      init.nrhs = atoi(optarg);
      break;
    case 'l':
      init.number_of_levels = atoi(optarg);
      break;
    case 'k':
      init.kappa = atof(optarg);
      break;
    case 'w':
      init.csw = atof(optarg);
      break;
    case 'u':
      init.mu = atof(optarg);
      break;
    case 't':
      init.number_openmp_threads = atoi(optarg);
      break;
    case '?':
    case 'h':
      p++;
      break;
    default: 
      break;
    }
  }

  if(conf_file == NULL) {
    printf0("Error: configuration file is missing (use -c PATH).\n");
    p++;
    fail++;
  }
  
  if(p) {
    help(argv[0]);
    MPI_Abort(MPI_COMM_WORLD,0);
  }
  if(fail)
    MPI_Abort(MPI_COMM_WORLD,0);
}

void read_params_arg(int argc, char *argv[] ) {

  int opt, p=0, fail=0, mu;
  optind = 0;
  while ((opt = getopt(argc, argv, options)) != -1) {
    switch (opt) {
    case 'r':
      residual = atof(optarg);
      break;
    case 'K':
      params.kcycle_tolerance = atof(optarg);
      break;
    case 'C':
      params.coarse_tolerance = atof(optarg);
      break;
    case 'V':
      optind--;
      mu=0;
      for ( ; optind < argc && *argv[optind] != '-'; optind++){
	if(mu > 2) {
	  printf0("Error: too many arguments in -V.\n");
	  p++;
	  fail++;
	  break;
	}
	params.mg_basis_vectors[mu] = atoi(argv[optind]);
	mu++;        
      }
      break;
    case 'm':
      optind--;
      params.mu_factor[0]=1;
      mu=1;
      for ( ; optind < argc && *argv[optind] != '-'; optind++){
	if(mu > 3) {
	  printf0("Error: too many arguments in -m.\n");
	  p++;
	  fail++;
	  break;
	}
	params.mu_factor[mu] = atof(argv[optind]);
	mu++;        
      }
      break;
    case 's':
      optind--;
      mu=0;
      for ( ; optind < argc && *argv[optind] != '-'; optind++){
	if(mu > 2) {
	  printf0("Error: too many arguments in -s.\n");
	  p++;
	  fail++;
	  break;
	}
	params.setup_iterations[mu] = atoi(argv[optind]);
	mu++;        
      }
      break;
    case 'v':
      params.print = 1;
      break;
    case 'D':
      optind--;
      mu=0;
      for ( ; optind < argc && *argv[optind] != '-'; optind++){
        if(mu > 3) {
          printf0("Error: too many arguments in -so.\n");
          p++;
          fail++;
          break;
        }
        params.solver[mu] = atoi(argv[optind]);
        mu++;
      }
      break;
    case 'g':
      optind--;
      mu=0;
      for ( ; optind < argc && *argv[optind] != '-'; optind++){
        if(mu > 3) {
          printf0("Error: too many arguments in -ot.\n");
          p++;
          fail++;
          break;
        }
        params.fab_orthoscheme[mu] = atoi(argv[optind]);
        mu++;
      }
      break;
      optind--;
      mu=0;
      for ( ; optind < argc && *argv[optind] != '-'; optind++){
        if(mu > 3) {
          printf0("Error: too many arguments in -os.\n");
          p++;
          fail++;
          break;
        }
        params.fab_orthoscheme[mu] = atoi(argv[optind]);
        mu++;
      }
      break;
    case 'I':
      optind--;
      mu=0;
      for ( ; optind < argc && *argv[optind] != '-'; optind++){
        if(mu > 3) {
          printf0("Error: too many arguments in -oi.\n");
          p++;
          fail++;
          break;
        }
        params.fab_ortho_iter[mu] = atoi(argv[optind]);
        mu++;
      }
      break;
    case 'M':
      optind--;
      mu=0;
      for ( ; optind < argc && *argv[optind] != '-'; optind++){
	if(mu > 3) {
          printf0("Error: too many arguments in -fm.\n");
          p++;
          fail++;
          break;
        }
        params.fab_max_kept_direction[mu] = atoi(argv[optind]);
        mu++;
      }
      break;
    case 'e':
      optind--;
      mu=0;
      for ( ; optind < argc && *argv[optind] != '-'; optind++){
        if(mu > 3) {
          printf0("Error: too many arguments in -fk.\n");
          p++;
          fail++;
          break;
        }
        params.fab_num_deflating_eig[mu] = atoi(argv[optind]);
        mu++;
      }
      break;
    case '?':
    default: 
      break;
    }
  }

  if(conf_file == NULL) {
    printf0("Error: configuration file is missing (use -c PATH).\n");
    p++;
    fail++;
  }
  
  if(p) {
    help(argv[0]);
    MPI_Abort(MPI_COMM_WORLD,0);
  }
  if(fail)
    MPI_Abort(MPI_COMM_WORLD,0);
}

int main( int argc, char *argv[] ) {

  
  MPI_Init( &argc, &argv );
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  standard_init();
  read_init_arg(argc, argv);
  
  printf0("Running initialization...\n");
  DDalphaAMG_initialize( &init, &params, &status );
  printf0("Initialized %d levels in %.2f sec\n", status.success, status.time);

  int nlvl = status.success;
  read_params_arg(argc, argv);

  comm_cart =  DDalphaAMG_get_communicator();
  MPI_Comm_rank( comm_cart, &rank );

  printf0("Running updating\n");
  DDalphaAMG_update_parameters( &params, &status );
  if (status.success)
    printf0("Updating time %.2f sec\n", status.time);

  /*
   * Reading the configuration. In plaq, it returns the plaquette value
   *  if provided in the configuration header.
   */
  double *gauge_field;
  int vol = init.global_lattice[T] * init.global_lattice[X] * init.global_lattice[Y] *
    init.global_lattice[Z] / init.procs[T] / init.procs[X] / init.procs[Y] / init.procs[Z];
  gauge_field = (double *) malloc(18*4*vol*sizeof(double));

  printf0("Reading config.\n");
  DDalphaAMG_read_configuration( gauge_field, conf_file, conf_format, &status );
  printf0("Reading configuration time %.2f sec\n", status.time);
  printf0("Desired plaquette %.13lf\n", status.info);
      
  printf0("Setting config.\n");
  DDalphaAMG_set_configuration( gauge_field, &status );
  printf0("Setting configuration time %.2f sec\n", status.time);
  printf0("Computed plaquette %.13lf\n", status.info);
  

  /*
   * Setup the solver
   */
  int i, nrhs = init.nrhs;
  printf0("Running setup\n");
  DDalphaAMG_setup( &status );
  for( i=0; i<nlvl; i++ ) {
    printf0("Run %d setup iterations in %.2f sec (%.1f %% on coarse grid at and below depth = %d)\n", status.success,
            status.time, 100.*(status.iter_times[i]/status.time), i);
    printf0("Total iterations at depth %d %d\n", status.iter_counts[i],i);
  }
 
  /*
   * Defining fine and coarse vector randomly.
   */
  double *vector1[nlvl], *vector2[nlvl];
  int vols[nlvl], vars[nlvl], sizes [nlvl];
  vols[0]=vol;
  vars[0]=3*4*2;
  sizes[0] = vols[0]*vars[0];

  for ( int i=1; i<nlvl; i++ ) {
    vols[i] = vols[i-1] / params.block_lattice[i-1][T] / params.block_lattice[i-1][X] / params.block_lattice[i-1][Y] / params.block_lattice[i-1][Z];
    vars[i]=params.mg_basis_vectors[i-1]*2*2; // a factor of 2 is for the spin, the other for the complex
    sizes[i] = vols[i]*vars[i];
  }
  
  for ( int i=0; i<nlvl; i++ ) {
    vector1[i] = (double *) malloc(vars[i]*vols[i]*sizeof(double)*nrhs);
    vector2[i] = (double *) malloc(vars[i]*vols[i]*sizeof(double)*nrhs);
  }

  for ( int i=0; i<nlvl; i++ )
    for ( int j=0; j<vars[i]*vols[i]*nrhs; j++ )
      vector1[i][j] = ((double)rand()/(double)RAND_MAX)-0.5;


  for ( int i=1; i<nlvl; i++ ) {
    printf0("Testing RP=1 on level %d\n",i);
    DDalphaAMG_prolongate(vector2[i-1], vector1[i], i-1, &status);
    DDalphaAMG_restrict(vector2[i], vector2[i-1], i-1, &status);

    for ( int n=0; n<nrhs; n++ ) {
      double num=0, den=0;
      for ( int j=0; j<sizes[i]; j++ ) {
	vector2[i][j+n*sizes[i]] -= vector1[i][j+n*sizes[i]];
	num += vector2[i][j+n*sizes[i]]*vector2[i][j+n*sizes[i]];
	den += vector1[i][j+n*sizes[i]]*vector1[i][j+n*sizes[i]];
      }
      printf0("Restult (1-RP)v[%d] = %e\n\n", n, num/den);
    }
  }

  for ( int i=1; i<nlvl; i++ ) {
    printf0("Testing coarse operator on level %d\n",i);
    DDalphaAMG_prolongate(vector1[i-1], vector1[i], i-1, &status);
    DDalphaAMG_apply_coarse_operator(vector2[i-1], vector1[i-1], i-1, &status);
    DDalphaAMG_restrict(vector2[i], vector2[i-1], i-1, &status);

    DDalphaAMG_apply_coarse_operator(vector1[i], vector1[i], i, &status);
    for ( int n=0; n<nrhs; n++ ) {
      double num=0, den=0;
      for ( int j=0; j<sizes[i]; j++ ) {
        vector2[i][j+n*sizes[i]] -= vector1[i][j+n*sizes[i]];
        num += vector2[i][j+n*sizes[i]]*vector2[i][j+n*sizes[i]];
        den += vector1[i][j+n*sizes[i]]*vector1[i][j+n*sizes[i]];
      }
      printf0("Restult (D_c-RDP)v[%d] = %e\n\n", n, num/den);
    }
  }
  
  //  free(vector_in);
  // free(vector_out);
  //free(gauge_field);
  //  DDalphaAMG_finalize();
  MPI_Finalize();
}
