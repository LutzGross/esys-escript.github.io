#include <Python.h>
#include <mpi.h>
#include <iostream>
#include <stdexcept>

extern "C"{
#include "paso/Paso.h"
#include "finley/Finley.h"
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

      sprintf( fname, "log_P%d.txt", mpi_info->rank );
      FILE *fp = freopen( fname, "w+", stdout );
    }
    /*
     * Start the python parser
     */
    status = Py_Main(argc, argv);
  
    /*
     * Finalise MPI for a clean exit.
     */
    MPI_Finalize();

    Paso_MPIInfo_dealloc( mpi_info );
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
