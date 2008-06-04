
/* $Id$ */

/*******************************************************
 *
 *           Copyright 2003-2007 by ACceSS MNRF
 *       Copyright 2007 by University of Queensland
 *
 *                http://esscc.uq.edu.au
 *        Primary Business: Queensland, Australia
 *  Licensed under the Open Software License version 3.0
 *     http://www.opensource.org/licenses/osl-3.0.php
 *
 *******************************************************/

#include "Utils.h"
#include "DataVector.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef PASO_MPI
#include <mpi.h>
#endif

#ifdef  _WIN32
#include <WinSock2.h>
#endif

namespace escript {

int getSvnVersion() 
{
#ifdef SVN_VERSION
  return SVN_VERSION;
#else
  return 0;
#endif
}

void printParallelThreadCnt() 
{
  int mpi_iam=0, mpi_num=1;
  char hname[64];

  gethostname(hname, 64);

  #ifdef PASO_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_iam);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_num);
  #endif

  #ifdef _OPENMP
  #pragma omp parallel
  #endif
  {
    int omp_iam=0, omp_num=1;
    #ifdef _OPENMP
    omp_iam = omp_get_thread_num(); /* Call in a parallel region */
    omp_num = omp_get_num_threads();
    #endif
    printf("printParallelThreadCounts: MPI=%d/%d OpenMP=%d/%d running on %s\n", mpi_iam, mpi_num, omp_iam, omp_num, hname);
  }
}

void setNumberOfThreads(const int num_threads) 
{

   #ifdef _OPENMP
   omp_set_num_threads(num_threads);
   #endif

}

int getNumberOfThreads() 
{
   #ifdef _OPENMP
   return omp_get_max_threads();
   #else
   return 1;
   #endif

}

}  // end of namespace
