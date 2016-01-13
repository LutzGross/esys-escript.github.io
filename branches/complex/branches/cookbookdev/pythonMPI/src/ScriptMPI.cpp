
/*******************************************************
*
* Copyright (c) 2003-2009 by University of Queensland
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
      #ifdef _WIN32
      fp_out = freopen( "NUL", "w+", stdout );
      #else
      fp_out = freopen( "/dev/null", "w+", stdout );
      #endif
    }
    /*
     * Start the python parser
     */
    status = Py_Main(argc, argv);
 
    /*
     * Close down MPI.
     * status==1 : uncaught python exception
     * status==2 : invalid python cmd line
     * status>2 : supposed to be param of sys.exit()
     *            sys.exit doesn't return as it should.
     *
     * I have made an exception for 2 because calling MPI_Abort
     * can display pretty ugly messages for not typing params
     * properly.
     */ 
    if ((status!=0) && (status!=2))
    {
	MPI_Abort(MPI_COMM_WORLD,status);
    }
    else
    { 
    /*
     * Finalise MPI for a clean exit.
     */
        MPI_Finalize();
    }
    Paso_MPIInfo_free( mpi_info );
  }
  catch (std::runtime_error &e)
  {
    std::cerr << "EXCEPTION: " << e.what() << std::endl;
    MPI_Abort(MPI_COMM_WORLD,1);
    throw;
  }
  catch (char *e)
  {
    std::cerr << "EXCEPTION: " << e << std::endl;
    MPI_Abort(MPI_COMM_WORLD,1);
    throw;
  }
  catch (...)
  {
    std::cerr << "EXCEPTION: " << "UNKNOWN." << std::endl;
    MPI_Abort(MPI_COMM_WORLD,1);
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
