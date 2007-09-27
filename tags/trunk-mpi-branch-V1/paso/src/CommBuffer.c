/*******************************************************************

  a buffer type for communication of "boundary" data between domains.
  send and receive buffers are allocated for each neighbouring domain,
  and can be used for sending any type of data distributed over the domains.
  Looks after communication completion, and safeguards against message collision
  etc.

*******************************************************************/

#include "Paso.h"

#ifdef PASO_MPI
#include <string.h>

#include "CommBuffer.h"

static index_t Paso_CommBuffer_checkDomain( Paso_CommBuffer *in, index_t dom )
{
  index_t position;

  if( !in || !in->numDomains )
  { 
    Paso_setError( VALUE_ERROR, "Paso_CommBuffer_checkDomain() : invalid or improperly uninitialised CommBuffer"); 
    return -1; 
  } 
  if( dom<0 || dom>(in->MPIInfo->size-1) ) 
  { 
    Paso_setError( VALUE_ERROR, "Paso_CommBuffer_checkDomain() : domain not in the communicator"); 
    return -1; 
  } 
  position = in->indexDomains[dom]; 
  if( position<0 ) 
  { 
    Paso_setError(  VALUE_ERROR,"Paso_CommBuffer_checkDomain() : domain not a neighbour of this domain" );  
    return -1;
  }

  return position;
}

/***************************************
 *  alloc/dealloc for CommBuffer type
 ***************************************/
Paso_CommBuffer *Paso_CommBuffer_alloc( Paso_MPIInfo *MPIInfo, index_t tag )
{
  Paso_CommBuffer *out=NULL;
  index_t i;

  out = MEMALLOC( 1, Paso_CommBuffer );

  if( Paso_checkPtr(out) )
    return NULL;

  out->reference_counter = 0;

  out->MPIInfo = Paso_MPIInfo_getReference( MPIInfo );
  if( !Paso_noError() )
  {
    MEMFREE(out);
    return NULL;
  }

  out->bufferForward = out->bufferBackward = NULL;
  out->numForward = out->numBackward = NULL;
  out->domains = NULL;
  out->indexDomains = MEMALLOC( MPIInfo->size, index_t );
  for( i=0; i<MPIInfo->size; i++ )
    out->indexDomains[i] = -1;
  out->statusForward = NULL;  
  out->requestForward = NULL;
  out->statusBackward = NULL;  
  out->requestBackward = NULL;
	out->requestedRecvLength = NULL;
  out->tag = tag;
  out->numDomains = 0;
  out->maxItemSize = 0;
  out->reference_counter++;

  return out;
}

void Paso_CommBuffer_dealloc( Paso_CommBuffer *in )
{
  if( in && !(--in->reference_counter) )
  {
    index_t i;

    if( in->bufferForward )
      for( i=0; i<in->numDomains; i++ ) {
        MEMFREE( in->bufferForward[i] );
      }
    MEMFREE( in->bufferForward );

    if( in->bufferBackward )
      for( i=0; i<in->numDomains; i++ ) {
        MEMFREE( in->bufferBackward[i] );
      }
    MEMFREE( in->bufferBackward );
		MEMFREE( in->requestedRecvLength);
    MEMFREE( in->numForward );
    MEMFREE( in->numBackward );
    MEMFREE( in->domains );
    MEMFREE( in->indexDomains );
    MEMFREE( in->statusForward );
    MEMFREE( in->requestForward );
    MEMFREE( in->statusBackward );
    MEMFREE( in->requestBackward );

    Paso_MPIInfo_dealloc( in->MPIInfo );

    MEMFREE( in );
  }
}

Paso_CommBuffer *Paso_CommBuffer_getReference( Paso_CommBuffer *in )
{
  if( !in )
    return NULL;

  in->reference_counter++;

  return in;
}

void Paso_CommBuffer_allocTable( Paso_CommBuffer *in, size_t itemSize, index_t *numForward, index_t *numBackward, dim_t numDomains, index_t *domains )
{
  index_t i;

  if( in && ( !in->domains || itemSize>in->maxItemSize ) )
  {
    /* Setup the domain info if it has not already been done.
       In effect, this says that the domain info is fixed for the lifespan of the
       CommBuffer. For example, one cannot reconfigure the structure for DOF 
       data to element data. */
    if( !in->domains ) 
    {
      if( numDomains>0 )
      {
        in->numForward  = MEMALLOC( numDomains, index_t );
        in->numBackward = MEMALLOC( numDomains, index_t );
        in->domains     = MEMALLOC( numDomains, index_t );
        in->requestedRecvLength= MEMALLOC( numDomains, dim_t);
        in->requestForward  = MEMALLOC( numDomains, MPI_Request );
        in->statusForward   = MEMALLOC( numDomains, MPI_Status );
        in->requestBackward = MEMALLOC( numDomains, MPI_Request );
        in->statusBackward  = MEMALLOC( numDomains, MPI_Status );

        if( ( Paso_checkPtr(in->numForward) || Paso_checkPtr(in->requestedRecvLength) || Paso_checkPtr(in->numBackward) || Paso_checkPtr(in->domains) || Paso_checkPtr(in->requestForward) || Paso_checkPtr(in->statusForward) ) )
        {
					 MEMFREE( in->requestedRecvLength );
           MEMFREE( in->numForward );
           MEMFREE( in->numBackward );
           MEMFREE( in->domains );
           MEMFREE( in->requestForward );
           MEMFREE( in->statusForward );
           MEMFREE( in->requestBackward );
           MEMFREE( in->statusBackward );
           return;
        }

        memcpy( in->numForward, numForward, numDomains*sizeof(index_t) );
        memcpy( in->numBackward, numBackward, numDomains*sizeof(index_t) );
        memcpy( in->domains, domains, numDomains*sizeof(index_t) );
        for( i=0; i<numDomains; i++ )
          in->requestForward[i] = in->requestBackward[i] = MPI_REQUEST_NULL;

        for( i=0; i<numDomains; i++ )   
          in->requestedRecvLength[i] =0;

        for( i=0; i<numDomains; i++ )   
          in->indexDomains[domains[i]] = i;
        
      }
      in->numDomains = numDomains;
    }

    /* setup the buffers, which may need to be reallocated in the case where
       the maximum itemSize has been reset */
    if( in->numDomains>0 )
    {
      if( !in->bufferForward )
      {
        in->bufferForward = MEMALLOC( in->numDomains, void* );
        for( i=0; i<in->numDomains; i++ )
          in->bufferForward[i] = NULL;
      }
      for( i=0; i<in->numDomains; i++ )
        if( in->numForward[i]>0 && itemSize>0 )
        {
           MEMREALLOC( in->bufferForward[i], in->numForward[i]*itemSize, char );
        }

      if( !in->bufferBackward )
      {
        in->bufferBackward = MEMALLOC( numDomains, void* );
        for( i=0; i<in->numDomains; i++ )
          in->bufferBackward[i] = NULL;
      }
      for( i=0; i<in->numDomains; i++ )
        if( in->numBackward[i]>0 && itemSize>0 )
        {
          MEMREALLOC( in->bufferBackward[i], in->numBackward[i]*itemSize, char );
        }
      
      in->maxItemSize = itemSize;
    }
  }
}

bool_t Paso_CommBuffer_recv( Paso_CommBuffer *in, index_t dom, size_t itemSize  )
{
  index_t position=0;
  int result;
  MPI_Status status;
  MPI_Request request;

  if( (position=Paso_CommBuffer_checkDomain( in, dom )) == -1 )
    return FALSE;
  
    /* ensure that the buffers are large enough */
  if( itemSize>in->maxItemSize )
    Paso_CommBuffer_allocTable( in, itemSize, NULL, NULL, -1, NULL );

  /* check that there isn't already a pending receive on the buffer */
  if( in->requestBackward[position]!=MPI_REQUEST_NULL )
  {
    Paso_setError( PASO_MPI_ERROR, "Paso_CommBuffer_recv() : cannot request a receive in a buffer that has an unresolved request" );    
    return FALSE;
  }

  result = MPI_Irecv( in->bufferBackward[position], in->numBackward[position]*itemSize, MPI_BYTE, dom, in->tag, in->MPIInfo->comm, in->requestBackward+position );

  if( result!=MPI_SUCCESS )
    return FALSE; 
  
	/* set the size of the requested receive for error checking when testing if the receive completed */
	in->requestedRecvLength[position] = itemSize*in->numBackward[position];

  return TRUE;
}

/* Assumes that there is/will be a message sent */
bool_t Paso_CommBuffer_recvAny( Paso_CommBuffer *in, index_t *dom, size_t itemSize )
{
  int result;
  MPI_Status status;
  
  /* ensure that the buffers are large enough */
  if( itemSize>in->maxItemSize )
    Paso_CommBuffer_allocTable( in, itemSize, NULL, NULL, -1, NULL );
    

  /* probe for a message. Note, this is a blocking call */
  MPI_Probe( MPI_ANY_SOURCE, in->tag, in->MPIInfo->comm, &status);
  *dom = status.MPI_SOURCE;

  return Paso_CommBuffer_recv( in, *dom, itemSize );
}

bool_t Paso_CommBuffer_send( Paso_CommBuffer *in, index_t dom, size_t itemSize )
{
  index_t position=0;
  int result, i;

  if( (position = Paso_CommBuffer_checkDomain( in, dom ))==-1 )
    return FALSE;

  /* check to make sure that the send buffer doesn't already have a pending send */
  if( Paso_CommBuffer_waitSend( in, dom )==FALSE )
    return FALSE;

  /* is the buffer large enough? if so, it is too late to resize it
     and we have to flag an error. */
  if( itemSize>in->maxItemSize )
  {
    Paso_setError( VALUE_ERROR, "Paso_CommBuffer_send() : the requested itemSize exceeds in->maxItemSize" );
    return FALSE;
  }

  /* send the buffer */
  result = MPI_Isend( in->bufferForward[position], itemSize*in->numForward[position], MPI_BYTE, dom, in->tag, in->MPIInfo->comm, &in->requestForward[position] );

  if( result!=MPI_SUCCESS )
    return FALSE;
  return TRUE;
}

/*
  waits for all pending sends in the buffer to complete. 
  returns immediately if an error is encountered.
*/
bool_t Paso_CommBuffer_waitSendPending( Paso_CommBuffer *in )
{
  index_t i;
  int success;

  for( i=0; i<in->numDomains; i++ )
    if( Paso_CommBuffer_waitSend( in, in->domains[i] )==FALSE )
      return FALSE;

  return TRUE;
}

/*
  waits for all a pending send in the buffer before continueing
*/
bool_t Paso_CommBuffer_waitSend( Paso_CommBuffer *in, index_t dom )
{
  index_t position;
  int success;

  if( (position = Paso_CommBuffer_checkDomain( in, dom ))==-1 )
    return FALSE;

  if( in->requestForward[position]!=MPI_REQUEST_NULL )
  {
    success = MPI_Wait( &in->requestForward[position], &in->statusForward[position] );

    if( success!=MPI_SUCCESS )
    {
      Paso_setError( PASO_MPI_ERROR, "Paso_CommBuffer_waitSend() : failed MPI_Isend" );
      return FALSE;
    }
  }

  return TRUE;
}

/*
  waits for all pending receives in the buffer to complete. 
  returns immediately if an error is encountered.
*/
bool_t Paso_CommBuffer_waitRecvPending( Paso_CommBuffer *in )
{
  index_t i;
  int success;

  for( i=0; i<in->numDomains; i++ )
    if( Paso_CommBuffer_waitRecv( in, in->domains[i] )==FALSE )
      return FALSE;

  return TRUE;
}

/*
  waits for all a pending receives in the buffer before continueing
*/
bool_t Paso_CommBuffer_waitRecv( Paso_CommBuffer *in, index_t dom )
{
  index_t position;
  int success;

  if( (position = Paso_CommBuffer_checkDomain( in, dom ))==-1 )
    return FALSE;

  if( in->requestBackward[position]!=MPI_REQUEST_NULL )
  {
    success = MPI_Wait( &in->requestBackward[position], &in->statusBackward[position] );

    if( success!=MPI_SUCCESS )
    {
      Paso_setError( PASO_MPI_ERROR, "Paso_CommBuffer_waitRecv() : failed MPI_Irecv" );
      return FALSE;
    }

		/* verify that the received buffer was the length of the requested buffer  */			
		MPI_Get_count( in->statusBackward + position, MPI_BYTE, &success );
		if( success!=in->requestedRecvLength[position] )
		{
      Paso_setError( PASO_MPI_ERROR, "Paso_CommBuffer_waitRecv() : size of received buffer and backward count are not equal" );
      return FALSE;
		}
  }

  return TRUE;
}


/*
  pack information into a send buffer
*/
void Paso_CommBuffer_pack( Paso_CommBuffer *in, index_t dom, index_t *index, void *data, size_t itemSize, dim_t offset )
{
  long i;
  index_t position;
  unsigned char *from, *to;

  /* wait for pending sends from the buffer */
  if( Paso_CommBuffer_waitSend( in, dom )==FALSE )
      return;
  
  /* ensure that the buffers are large enough */
  if( itemSize>in->maxItemSize )
    Paso_CommBuffer_allocTable( in, itemSize, NULL, NULL, -1, NULL );

  /* setup pointers to regions for copying */
  if( (position = Paso_CommBuffer_checkDomain( in, dom ))==-1 )
    return;
  if( in->numForward[position]==0 )
    return;
  from = (unsigned char*)data;
  to =   (unsigned char*)in->bufferForward[position];
  from += offset*itemSize;

  if( index==NULL && in->numForward[position]>0 )
  {
    /* if there is no index, pack the data as is */
    memcpy( to, from, itemSize*in->numForward[position] );
  }
  else
  {
    /* pack the data according to index */
    for( i=0; i<in->numForward[position]; i++, to+=itemSize )
      memcpy( to , from + (index[i]*itemSize), itemSize );
  }
} 

void Paso_CommBuffer_unpack( Paso_CommBuffer *in, index_t dom, index_t *index, void *data, size_t itemSize, dim_t offset )
{
  long i;
  index_t position;
  unsigned char *from, *to;

  /* wait for pending receives to the buffer */
  if( Paso_CommBuffer_waitRecv( in, dom )==FALSE )
      return;

  /* ensure that the buffers are large enough */
  if( itemSize>in->maxItemSize )
    Paso_CommBuffer_allocTable( in, itemSize, NULL, NULL, -1, NULL );

  /* setup pointers to regions for copying */
  if( (position = Paso_CommBuffer_checkDomain( in, dom ))==-1 )
    return;
  if( in->numBackward[position]==0 )
    return;

  from = (unsigned char*)in->bufferBackward[position];
  to =   (unsigned char*)data;
  to += offset*itemSize;

  if( index==NULL && in->numBackward[position]>0 )
  {
    /* if there is no index, unpack the data as is */
    memcpy( to, from, itemSize*in->numBackward[position] );
  }
  else
  {
    /* unpack the data according to index */
    for( i=0; i<in->numBackward[position]; i++ )
			memcpy( to + itemSize*index[i], from + (i*itemSize), itemSize );
  }
}

/*
  verify that the information stored accross the processors on the communicator
  over which the CommBuffer is allocated is valid
*/
bool_t Paso_CommBuffer_validate( Paso_CommBuffer *in )
{
  index_t *tmpForwardMap=NULL, *tmpBackwardMap=NULL, *tmpConnectionMap=NULL, *tmpForward=NULL, *tmpBackward=NULL, *tmpConnection=NULL;  
  dim_t size, rank, i;

  if( in && in->MPIInfo->size>0 )  
  {
    size = in->MPIInfo->size;
    rank = in->MPIInfo->rank;

    tmpForwardMap = MEMALLOC( size*size, index_t );
    tmpBackwardMap = MEMALLOC( size*size, index_t );
    tmpConnectionMap = MEMALLOC( size*size, index_t );
    tmpForward = MEMALLOC( size, index_t );
    tmpBackward = MEMALLOC( size, index_t );
    tmpConnection = MEMALLOC( size, index_t );

    /* collect global image of the distribution */
    for( i=0; i<size; i++ )
      tmpConnection[i] = 0;
    for( i=0; i<in->numDomains; i++ )
      tmpConnection[in->domains[i]] = 1;
    MPI_Allgather( tmpConnection, size, MPI_INT, tmpConnectionMap, size, MPI_INT, in->MPIInfo->comm );

    for( i=0; i<size; i++ )
      tmpForward[i] = 0;
    for( i=0; i<in->numDomains; i++ )
      tmpForward[in->domains[i]] = in->numForward[i];
    MPI_Allgather( tmpForward, size, MPI_INT, tmpForwardMap, size, MPI_INT, in->MPIInfo->comm );  


    for( i=0; i<size; i++ )
      tmpBackward[i] = 0;
    for( i=0; i<in->numDomains; i++ )
      tmpBackward[in->domains[i]] = in->numBackward[i];
    MPI_Allgather( tmpBackward, size, MPI_INT, tmpBackwardMap, size, MPI_INT, in->MPIInfo->comm );

    /* verify that information on different processors is consisent */
    for( i=0; i<size; i++ )
    {
      if( tmpConnection[i] != tmpConnectionMap[rank+i*size] )
      {
        Paso_setError( VALUE_ERROR, "Paso_CommBuffer_validate() : neighbour connection map is inconsistent" );
        goto clean;
      }
    }
    for( i=0; i<size; i++ )
    {
      if( tmpForward[i] != tmpBackwardMap[rank+i*size] )
      {
        Paso_setError( VALUE_ERROR, "Paso_CommBuffer_validate() : neighbour forward map is inconsistent" );
        goto clean;
      }
    }
    for( i=0; i<size; i++ )
    {
      if( tmpBackward[i] != tmpForwardMap[rank+i*size] )
      {
        Paso_setError( VALUE_ERROR, "Paso_CommBuffer_validate() : neighbour backward map is inconsistent" );
        goto clean;
      }
    }    

  }

clean :
  MEMFREE( tmpBackwardMap );
  MEMFREE( tmpForwardMap );
  MEMFREE( tmpConnectionMap );
  MEMFREE( tmpConnection );
  MEMFREE( tmpForward );
  MEMFREE( tmpBackward );

  return Paso_MPI_noError( in->MPIInfo );
}

#endif
