#include "Distribution.h"

#ifdef PASO_MPI

void Finley_NodeDistribution_formCommBuffer( Finley_NodeDistribution *in, Paso_CommBuffer *CommBuffer )
{
  index_t *numForward=NULL, *numBackward=NULL;
  index_t i;  

  numForward  = MEMALLOC( in->numNeighbours, index_t );
  numBackward = MEMALLOC( in->numNeighbours, index_t );

  if( Finley_checkPtr(numForward) || Finley_checkPtr(numBackward) ){
    MEMFREE( numForward );  
    MEMFREE( numBackward );
    return; 
  }

  for( i=0; i<in->numNeighbours; i++ ){
    numForward[i] = in->edges[i]->numForward;
    numBackward[i] = in->edges[i]->numBackward;
  }

  Paso_CommBuffer_allocTable( CommBuffer, FINLEY_INIT_ITEMSIZE, numForward, numBackward, in->numNeighbours, in->neighbours );

  MEMFREE( numForward );
  MEMFREE( numBackward );
}

#endif
