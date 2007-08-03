#include <Python.h>
#include <mpi.h>
#include <iostream>
#include <stdexcept>

extern "C"{
#include "paso/Paso_MPI.h"
}
#ifdef PASO_MPI

int main( int argc, char **argv ) {
  int status = 0;
  Paso_MPIInfo *mpi_info=NULL;
  try
  {
    /*
     * Initialise MPI
     */
    status = MPI_Init(&argc, &argv);
    if (status != MPI_SUCCESS) {
      std::cerr << argv[0] << ": MPI_Init failed, exiting." << std::endl;
      return status;
    }
    mpi_info = Paso_MPIInfo_alloc( MPI_COMM_WORLD );

    if( mpi_info->rank )
    {
      char fname[256];
      sprintf( fname, "stdout_cpu_%04d.out", mpi_info->rank );
      FILE *fp_out = freopen( fname, "w+", stdout );
      sprintf( fname, "stdout_cpu_%04d.err", mpi_info->rank );
      FILE *fp_err = freopen( fname, "w+", stderr );
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
