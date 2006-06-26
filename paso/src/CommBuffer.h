#ifndef INC_COMMBUFFER
#define INC_COMMBUFFER

#include "Paso.h"

#ifdef PASO_MPI
      
#define TRACE_COMMBUFFER

/***********************************
 * TYPE DEFINITIONS
 ***********************************/

/* buffer for communication */
typedef struct Paso_CommBuffer
{
  Paso_MPIInfo *MPIInfo;
  index_t reference_counter;

  void **bufferForward;
  void **bufferBackward;
  dim_t *numForward; 
  dim_t *numBackward;
  size_t maxItemSize;
	dim_t *requestedRecvLength;

  /* the source/target domains for the buffer */
  dim_t numDomains;
  index_t *domains;
  index_t *indexDomains;

  MPI_Status  *statusForward;  
  MPI_Request *requestForward;
  MPI_Status  *statusBackward;  
  MPI_Request *requestBackward;

  int tag; /* tag used for sends a receives for this buffer */
} Paso_CommBuffer;



/***********************************
 * PROTOTYPES
 ***********************************/

Paso_CommBuffer* Paso_CommBuffer_alloc( Paso_MPIInfo *MPIInfo, int tag );
Paso_CommBuffer* Paso_CommBuffer_getReference( Paso_CommBuffer *in );
void   Paso_CommBuffer_dealloc( Paso_CommBuffer *in );
void   Paso_CommBuffer_allocTable( Paso_CommBuffer *in, size_t itemSize, index_t *numForward, index_t *numBackward, dim_t numDomains, index_t *domains );
bool_t Paso_CommBuffer_recv( Paso_CommBuffer *in, index_t dom, size_t itemSize );
bool_t Paso_CommBuffer_recvAny( Paso_CommBuffer *in, index_t *dom, size_t itemSize );
bool_t Paso_CommBuffer_send( Paso_CommBuffer *in, index_t dom, size_t itemSize );
bool_t Paso_CommBuffer_waitSendPending( Paso_CommBuffer *in );
bool_t Paso_CommBuffer_waitSend( Paso_CommBuffer *in, index_t dom );
bool_t Paso_CommBuffer_waitRecvPending( Paso_CommBuffer *in );
bool_t Paso_CommBuffer_waitRecv( Paso_CommBuffer *in, index_t dom );
void*  Paso_CommBuffer_getBufferForward( Paso_CommBuffer *in, index_t dom );
void   Paso_CommBuffer_pack( Paso_CommBuffer *in, index_t dom, index_t *index, void *data, size_t itemSize, dim_t offset );
void   Paso_CommBuffer_unpack( Paso_CommBuffer *in, index_t dom, index_t *index, void *data, size_t itemSize, dim_t offset );
bool_t Paso_CommBuffer_validate( Paso_CommBuffer *in );

#endif

#endif
