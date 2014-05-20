
/*******************************************************
*
* Copyright (c) 2003-2011 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


#include <Python.h>
#include <iostream>
#include <stdexcept>

extern "C"{
#include "esysUtils/Esys_MPI.h"
}

#ifdef ESYS_MPI

int main( int argc, char **argv ) {
  int status = 0;
  int provided;
  Esys_MPIInfo *mpi_info=NULL;
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
    mpi_info = Esys_MPIInfo_alloc( MPI_COMM_WORLD );

    if( mpi_info->rank )
    {
      char fname[256];
      sprintf( fname, "stdout_%04d.out", mpi_info->rank );
      /*FILE *fp_out =*/ freopen( fname, "w+", stdout );
      sprintf( fname, "stdout_%04d.err", mpi_info->rank );
      /*FILE *fp_err =*/ freopen( fname, "w+", stderr );
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
     * Yes that means you probably shouldn't use 2 as an exit code
     * but we don't really recommend sys.exit anyway.
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

    Esys_MPIInfo_free( mpi_info );
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
	printf( "Esys must be compiled with ESYS_MPI defined to make the MPI version available\n\n" );
	return 0;
}
#endif // ESYS_MPI

