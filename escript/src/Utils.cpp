
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

#include <unistd.h>
#include <string.h>

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

/* This is probably not very robust, but it works on Savanna today and is useful for performance analysis */
int get_core_id() {
  int processor_num=-1;
#ifdef CORE_ID1
  FILE *fp;
  int i, count_spaces=0;
  char fname[100];
  char buf[1000];

  sprintf(fname, "/proc/%d/stat", getpid());
  fp = fopen(fname, "r");
  if (fp == NULL) return(-1);
  fgets(buf, 1000, fp);
  fclose(fp);

  for (i=strlen(buf)-1; i>=0; i--) {
    if (buf[i] == ' ') count_spaces++;
    if (count_spaces == 4) break;
  }
  processor_num = atoi(&buf[i+1]);
#endif
  return(processor_num);
}


void printParallelThreadCnt() 
{
  int mpi_iam=0, mpi_num=1;
  char hname[64];

#ifdef HAVE_GETHOSTNAME
  gethostname(hname, 64);
  hname[63] = '\0';
#else
  strcpy(hname, "unknown host");
#endif

  #ifdef PASO_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_iam);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_num);
  #endif

  #pragma omp parallel
  {
    int omp_iam=0, omp_num=1;
    #ifdef _OPENMP
    omp_iam = omp_get_thread_num(); /* Call in a parallel region */
    omp_num = omp_get_num_threads();
    #endif
    printf("printParallelThreadCounts: MPI=%03d/%03d OpenMP=%03d/%03d running on %s core %d\n",
      mpi_iam, mpi_num, omp_iam, omp_num, hname, get_core_id());
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

int getMPISizeWorld() {
  int mpi_num = 1;
  #ifdef PASO_MPI
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_num);
  #endif
  return mpi_num;
}

int getMPIRankWorld() {
  int mpi_iam = 0;
  #ifdef PASO_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_iam);
  #endif
  return mpi_iam;
}

}  // end of namespace
