
/*****************************************************************************
*
* Copyright (c) 2003-2016 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/


#include <Python.h>
#include <iostream>
#include <stdexcept>

#include "esysUtils/Esys_MPI.h"

#ifdef ESYS_MPI

int main( int argc, char **argv ) {
  int status = 0;
  int provided;
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
    int rank=0;
    if (MPI_Comm_rank(MPI_COMM_WORLD, &rank)!=MPI_SUCCESS)
    {
	MPI_Abort(MPI_COMM_WORLD, 1);
    }

    if( rank )
    {
      #ifdef _WIN32
         if ( freopen( "NUL", "w+", stdout ) == NULL ) {
           exit(EXIT_FAILURE);
         }
      #else
         if ( freopen( "/dev/null", "w+", stdout ) == NULL ) {
           exit(EXIT_FAILURE);
         }
      #endif
    }
    /*
     * Start the python parser
     */

#ifdef ESPYTHON3 
    wchar_t** wargv=new wchar_t*[argc+1];
    for (int i=0;i<argc;++i)
    {
        int len=strlen(argv[i]);
	wargv[i]=new wchar_t[len+1];
	for (int j=0;j<len;++j) 
	{
            wargv[i][j]=wchar_t(argv[i][j]);
	}
	wargv[i][len]=0;
    }
    wargv[argc]=0;
    status = Py_Main(argc, wargv);
#else    


    status = Py_Main(argc, argv);
#endif
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

