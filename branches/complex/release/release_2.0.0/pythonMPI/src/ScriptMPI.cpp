
/*******************************************************
*
* Copyright (c) 2003-2008 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


#include <Python.h>
#ifdef PASO_MPI
#include <mpi.h>
#endif
#include <iostream>
#include <stdexcept>

extern "C"{
#include "paso/Paso_MPI.h"
}
#ifdef PASO_MPI

int main( int argc, char **argv ) {
  int status = 0;
  int provided;
  Paso_MPIInfo *mpi_info=NULL;
  try
  {
    /*
     * Initialise MPI
     */
    /* status = MPI_Init(&argc, &argv); */
    status = MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided );
    if (status != MPI_SUCCESS) {
      std::cerr << argv[0] << ": MPI_Init failed, exiting." << std::endl;
      return status;
    }
    mpi_info = Paso_MPIInfo_alloc( MPI_COMM_WORLD );

    if( mpi_info->rank )
    {
      FILE *fp_out;
      fp_out = freopen( "/dev/null", "w+", stdout );
    }
    /*
     * Start the python parser
     */
    status = Py_Main(argc, argv);
  
    /*
     * Finalise MPI for a clean exit.
     */
    MPI_Finalize();

    Paso_MPIInfo_free( mpi_info );
  }
  catch (std::runtime_error &e)
  {
    std::cerr << "EXCEPTION: " << e.what() << std::endl;
    throw;
  }
  catch (char *e)
  {
    std::cerr << "EXCEPTION: " << e << std::endl;
    throw;
  }
  catch (...)
  {
    std::cerr << "EXCEPTION: " << "UNKNOWN." << std::endl;
    throw;
  }

  return status;
}

#else
int main( int argc, char **argv ) {
	printf( "Esys must be compiled with PASO_MPI defined to make the MPI version available\n\n" );
	return 0;
}
#endif
