#include "Distribution.h"

#ifdef PASO_MPI
void Finley_NodeDistribution_calculateIndexExternal( Finley_NodeDistribution *Distribution, Paso_CommBuffer *CommBuffer )
{
  index_t i, n, thisDom=Distribution->MPIInfo->rank;
  index_t *buffer=NULL;
  dim_t bufferLength=0;

  for( n=0; n<Distribution->numNeighbours; n++ ) {
    if( Distribution->edges[n]->numForward>bufferLength )
      bufferLength = Distribution->edges[n]->numForward;
    if( Distribution->edges[n]->numBackward>bufferLength )
      bufferLength = Distribution->edges[n]->numBackward;
  }

  buffer = TMPMEMALLOC( bufferLength, index_t );

  /* send updated global indices for local DOF that are to neighbours for whom they are external */
  for( n=0; n<Distribution->numNeighbours; n++ ){
    for( i=0; i<Distribution->edges[n]->numForward; i++ )
      buffer[i] = Distribution->vtxdist[thisDom] + Distribution->edges[n]->indexForward[i];
    Paso_CommBuffer_pack( CommBuffer, Distribution->neighbours[n], NULL, buffer, sizeof(index_t), 0 );
    Paso_CommBuffer_send( CommBuffer, Distribution->neighbours[n], sizeof(index_t) );
  }

  /* receive updated global indices for external DOF from neighbours */
  for( n=0; n<Distribution->numNeighbours; n++ ){
    Paso_CommBuffer_recv( CommBuffer, Distribution->neighbours[n], sizeof(index_t)  );      
    Paso_CommBuffer_unpack( CommBuffer, Distribution->neighbours[n], Distribution->edges[n]->indexBackward, Distribution->indexExternal, sizeof(index_t), -Distribution->numLocal );
  }

  TMPMEMFREE( buffer );
}
#endif

