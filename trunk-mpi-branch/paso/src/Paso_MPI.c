#include <stdlib.h>
#include <stdio.h>


#include "Paso_MPI.h"


/* allocate memory for an mpi_comm, and find the communicator details */
Paso_MPIInfo* Paso_MPIInfo_alloc( MPI_Comm comm )
{
  int error;
  Paso_MPIInfo *out=NULL;

  out = MEMALLOC( 1, Paso_MPIInfo );
  
  out->reference_counter = 0;
  out->msg_tag_counter = 0;
  #ifdef PASO_MPI
     error = MPI_Comm_rank( comm, &out->rank )==MPI_SUCCESS && MPI_Comm_size( comm, &out->size )==MPI_SUCCESS;
     if( !error ) {
       Paso_setError( PASO_MPI_ERROR, "Paso_MPIInfo_alloc : error finding comm rank/size" );
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
void Paso_MPIInfo_free( Paso_MPIInfo *in )
{
  if( in && !(--in->reference_counter) )
    MEMFREE( in );
}

Paso_MPIInfo *Paso_MPIInfo_getReference( Paso_MPIInfo* in )
{
  if (in!=NULL) 
    ++(in->reference_counter);
  
  return in;
}
/* N = #CPUs, k is a CPU number but out of range or even negative. Return a CPU number in 0...n-1. */
index_t Paso_MPIInfo_mod(index_t n, index_t k) 
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

void Paso_MPIInfo_Split( Paso_MPIInfo *mpi_info, dim_t N, dim_t* local_N,index_t* offset) 
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


dim_t Paso_MPIInfo_setDistribution(Paso_MPIInfo* mpi_info ,index_t min_id,index_t max_id,index_t* distribution) {
   int rest=0, p;
   dim_t out;
   int s=mpi_info->size;
   dim_t N=max_id-min_id+1;
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
}

/* checks that there is no error accross all processes in a communicator */
/* NOTE : does not make guarentee consistency of error string on each process */
bool_t Paso_MPIInfo_noError( Paso_MPIInfo *mpi_info )
{
  int errorLocal = 0;
  int errorGlobal= 0;
  errorLocal= Paso_noError() ? 0 : 1;
  if (mpi_info->size>1) {
     #ifdef PASO_MPI
#if 1 /* ksteube disable error checking during benchmarking activities */
     MPI_Allreduce( &errorLocal, &errorGlobal, 1, MPI_INT, MPI_MAX, mpi_info->comm  );
#else
     errorGlobal=errorLocal;
#endif
     #else
     errorGlobal=errorLocal;
     #endif
     /* take care of the case where the error was on another processor */
     if( (errorLocal==0) && (errorGlobal==1) ) {
         Paso_setError( PASO_MPI_ERROR, "Paso_MPI_noError() : there was an error on another MPI process" );
     }
  }
  return (bool_t) errorGlobal;
}


/**************************************************
                 WRAPPERS 
**************************************************/

int Paso_MPIInfo_initialized( void )
{
  int error=0, initialised=0;

  #ifdef PASO_MPI
     error = MPI_Initialized( &initialised );
     if( error!=MPI_SUCCESS )
         Paso_setError( PASO_MPI_ERROR, "mpi_initialised : MPI error" );
     return initialised;
  #else
     return TRUE;
  #endif
}
