#include <stdlib.h>
#include <stdio.h>


#include "Paso.h"

#ifdef PASO_MPI

/* allocate memory for an mpi_comm, and find the communicator details */
Paso_MPIInfo* Paso_MPIInfo_alloc( MPI_Comm comm )
{
  int error;
  Paso_MPIInfo *out=NULL;

  out = MEMALLOC( 1, Paso_MPIInfo );
  
  out->reference_counter = 0;
  error = MPI_Comm_rank( comm, &out->rank )==MPI_SUCCESS && MPI_Comm_size( comm, &out->size )==MPI_SUCCESS;
  if( !error ) {
    Paso_setError( PASO_MPI_ERROR, "Paso_MPIInfo_alloc : error finding comm rank/size" );
  }
  
  out->comm = comm;
  out->reference_counter++;

  return out;
}

/* free memory for an mpi_comm */
void Paso_MPIInfo_dealloc( Paso_MPIInfo *in )
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

/**************************************************
                 WRAPPERS 
**************************************************/

int Paso_MPI_initialized( void )
{
  int error=0, initialised=0;

  error = MPI_Initialized( &initialised );
  if( error!=MPI_SUCCESS )
    Paso_setError( PASO_MPI_ERROR, "mpi_initialised : MPI error" );

  return initialised;
}

#endif
