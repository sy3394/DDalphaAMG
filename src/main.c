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
#include <unistd.h>

global_struct g;
#ifdef HAVE_HDF5
Hdf5_fileinfo h5info;
#endif
struct common_thread_data *commonthreaddata;
struct Thread *no_threading; // no_threading is global both in serial and thread pralel region

int main( int argc, char **argv ) {
#if 0 // attempt to enable debugging MPI process via gdb
  int ifl = 0;
  char hostname[256];
  gethostname(hostname, sizeof(hostname));
  printf("PID %d on %s ready for attach\n", getpid(), hostname);
  fflush(stdout);
  while (0 == ifl)
    sleep(5);
  #endif
 
#ifdef HAVE_HDF5
  h5info.filename=NULL;
  h5info.file_id=-1; 
  h5info.rootgroup_id=-1; 
  h5info.configgroup_id=-1;
  h5info.eigenmodegroup_id=-1;
  h5info.thiseigenmodegroup_id=-1;
  h5info.isOpen=0;
  h5info.mode=-1;
#endif
  level_struct l;
  config_double hopp = NULL;

  int provided;
  MPI_Init_thread( &argc, &argv, MPI_THREAD_FUNNELED, &provided);
  MPI_Comm_rank( MPI_COMM_WORLD, &(g.my_rank) );
  
  if ( g.my_rank == 0 ) {
    if ( provided < MPI_THREAD_FUNNELED )
      warning("MPI_THREAD_FUNNELED is not ensured\n");
    printf("\n\n+----------------------------------------------------------+\n");
    printf("| The DDalphaAMG solver library.                           |\n");
    printf("| Copyright (C) 2016, Matthias Rottmann, Artur Strebel,    |\n");
    printf("|       Simon Heybrock, Simone Bacchio, Bjoern Leder,      |\n");
    printf("|       Shuhei Yamamoto                             .      |\n");
    printf("|                                                          |\n");
    printf("| This program comes with ABSOLUTELY NO WARRANTY.          |\n");
    printf("+----------------------------------------------------------+\n\n");
  }
  
  //------------ initialize and setup g and l according to the inputfile
  method_init( &argc, &argv, &l );

  // setup threading structures
  no_threading = (struct Thread *)malloc(sizeof(struct Thread));
  setup_no_threading(no_threading, &l);

  commonthreaddata = (struct common_thread_data *)malloc(sizeof(struct common_thread_data));
  init_common_thread_data(commonthreaddata);

  // store gauge configuration in hopp
  MALLOC( hopp, complex_double, 3*l.inner_vector_size );
  if(g.in_format == _LIME)
    lime_read_conf( (double*)(hopp), g.in, &(g.plaq_hopp) );
  else 
    read_conf( (double*)(hopp), g.in, &(g.plaq_hopp), &l );

  // setup Dirac matrices (gauge matrices, self-coupling terms) and compute plaquettes
  dirac_setup( hopp, &l );
  FREE( hopp, complex_double, 3*l.inner_vector_size );

  THREADED(g.num_openmp_processes)
  {
    //------------------ allocate memory and set up for multigrid solver
    struct Thread threading;
    setup_threading(&threading, commonthreaddata, &l);

    // setup up initial MG hierarchy (null vector setup)
    method_setup( NULL, &l, &threading );
    
    // iteratively update P and D_c
    method_iterative_setup( l.setup_iter, &l, &threading );
    //        error0("STOP\n");
    //----------------- solve
    solve_driver( &l, &threading );
  }
  printf0("Number of rhs vectors = %d\n", g.num_rhs_vect);
  method_free( &l );
  finalize_common_thread_data(commonthreaddata); // free workspace in commonthreaddata
  finalize_no_threading(no_threading);           // free workspace in no_threading
  free(commonthreaddata);
  free(no_threading);
  method_finalize( &l );
  
  MPI_Finalize();
  
  return 0;
}
