
/*******************************************************
*
* Copyright (c) 2003-2010 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


#include <stdlib.h>
#include <stdio.h>
#include <string.h>


#include "Esys_MPI.h"
#include "index.h"
#include "mem.h"
#include "error.h"


/* allocate memory for an mpi_comm, and find the communicator details */
Esys_MPIInfo* Esys_MPIInfo_alloc( MPI_Comm comm )
{
  #ifdef ESYS_MPI
    int error;
  #endif

  Esys_MPIInfo *out=NULL;

  out = MEMALLOC( 1, Esys_MPIInfo );
  
  out->reference_counter = 0;
  out->msg_tag_counter = 0;
  #ifdef ESYS_MPI
     error = MPI_Comm_rank( comm, &out->rank )==MPI_SUCCESS && MPI_Comm_size( comm, &out->size )==MPI_SUCCESS;
     if( !error ) {
       Esys_setError( ESYS_MPI_ERROR, "Esys_MPIInfo_alloc : error finding comm rank/size" );
     }
  
     out->comm = comm;
  #else
     out->rank=0;
     out->size=1;
     out->comm=-1;
  #endif
  out->reference_counter++;

  return out;
}

/* free memory for an mpi_comm */
void Esys_MPIInfo_free( Esys_MPIInfo *in )
{
  if( in && !(--in->reference_counter) )
    MEMFREE( in );
}

Esys_MPIInfo *Esys_MPIInfo_getReference( Esys_MPIInfo* in )
{
  if (in!=NULL) 
    ++(in->reference_counter);
  
  return in;
}
/* N = #CPUs, k is a CPU number but out of range or even negative. Return a CPU number in 0...n-1. */
index_t Esys_MPIInfo_mod(index_t n, index_t k) 
{
    index_t q, out=0;
    if (n>1) {
        q=k/n;
        if (k>0) {
           out=k-n*q;
        } else if (k<0) {
           out=k-n*(q-1);
        }
    }
    return out;
}

void Esys_MPIInfo_Split( Esys_MPIInfo *mpi_info, dim_t N, dim_t* local_N,index_t* offset) 
{
   int rest=0;
   int s=mpi_info->size;
   int r=mpi_info->rank;
   *local_N=N/s;
   rest=N-(*local_N)*s;
   if (r<rest) {
       (*local_N)++;
       (*offset)=(*local_N)*r;
   } else {
       (*offset)=(*local_N)*r+rest;
   }
}


dim_t Esys_MPIInfo_setDistribution(Esys_MPIInfo* mpi_info ,index_t min_id,index_t max_id,index_t* distribution) {
   int rest=0, p;
   dim_t out;
   int s=mpi_info->size;
   dim_t N=max_id-min_id+1;
   if (N>0) {
      int local_N=N/s;
      rest=N-local_N*s;
      for (p=0; p<s; ++p) {
         if (p<rest) {
             distribution[p]=min_id+(local_N+1)*p;
             out=local_N+1;
         } else {
             distribution[p]=min_id+rest+local_N*p;
         }
      }
      distribution[s]=max_id+1;
      if (rest==0) {
         return local_N;
      } else {
         return local_N+1;
      }
  } else {
      for (p=0; p<s+1; ++p) distribution[p]=min_id;
      return 0;
  }
}

/* checks that there is no error accross all processes in a communicator */
/* NOTE : does not make guarentee consistency of error string on each process */
bool_t Esys_MPIInfo_noError( Esys_MPIInfo *mpi_info )
{
  int errorLocal = Esys_noError() ? 0 : 1;
  int errorGlobal = errorLocal;

  return (errorGlobal==0);
}

/**************************************************
                 WRAPPERS 
**************************************************/

int Esys_MPIInfo_initialized( void )
{
  #ifdef ESYS_MPI
     int error=0, initialised=0;
     error = MPI_Initialized( &initialised );
     if( error!=MPI_SUCCESS )
         Esys_setError( ESYS_MPI_ERROR, "mpi_initialised : MPI error" );
     return initialised;
  #else
     return TRUE;
  #endif
}

/* Append MPI rank to file name if multiple MPI processes */
char *Esys_MPI_appendRankToFileName(const char *fileName, int mpi_size, int mpi_rank) {
  /* Make plenty of room for the mpi_rank number and terminating '\0' */
  char *newFileName = TMPMEMALLOC(strlen(fileName)+20,char);
  strncpy(newFileName, fileName, strlen(fileName)+1);
  if (mpi_size>1) sprintf(newFileName+strlen(newFileName), ".%04d", mpi_rank);
  return(newFileName);
}

#ifndef _OPENMP 
int omp_get_max_threads(void) {
   return 1;
}
int omp_get_thread_num(void) {
   return 0;
}
#endif

