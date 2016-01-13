/* created by Ben Cumming on 26/04/2006 */
#include "Distribution.h"

#ifdef PASO_MPI

Finley_NodeDistribution* Finley_NodeDistribution_alloc( Paso_MPIInfo *MPIInfo )
{
  Finley_NodeDistribution *out=NULL;

  out = MEMALLOC( 1, Finley_NodeDistribution );
  if (Finley_checkPtr(out)) return NULL;

  out->reference_counter = 0;

  out->MPIInfo = Paso_MPIInfo_getReference(MPIInfo);

  out->numLocal = 0;
  out->numGlobal = 0;
  out->numInternal = 0;
  out->numBoundary = 0;
  out->numExternal = 0;

  out->numNeighbours = 0;
  out->neighbours = NULL;
  out->edges = NULL;
  out->indexExternal = NULL;
  out->vtxdist=NULL;

  out->reference_counter++;

  return out;
}

void Finley_NodeDistribution_dealloc( Finley_NodeDistribution *in )
{
  index_t i;

  if( in && !(--in->reference_counter) )
  {
    Paso_MPIInfo_dealloc( in->MPIInfo );

    for( i=0; i<in->numNeighbours; i++ )
      Finley_NodeGhostEdge_dealloc( in->edges[i] );
    MEMFREE( in->edges );
    MEMFREE( in->vtxdist );
    MEMFREE( in->indexExternal );
    MEMFREE( in->neighbours );

    MEMFREE( in );
  } 
}

Finley_NodeDistribution* Finley_NodeDistribution_getReference( Finley_NodeDistribution *in )
{
  if( in ) 
    in->reference_counter++;
  
  return in;
}

Finley_NodeGhostEdge* Finley_NodeGhostEdge_alloc( void )
{
  Finley_NodeGhostEdge *out=NULL;

  out = MEMALLOC( 1, Finley_NodeGhostEdge );
  if (Finley_checkPtr(out)) return NULL;
  
  out->reference_counter = 0;
  out->numForward = out->numBackward = 0;
  out->indexForward = out->indexBackward = NULL;
  
  out->reference_counter++;  

  return out;
}

void Finley_NodeGhostEdge_dealloc( Finley_NodeGhostEdge *in )
{
  if( in && !(--in->reference_counter) )
  {
    MEMFREE( in->indexForward );
    MEMFREE( in->indexBackward );
    MEMFREE( in );
  }
}

Finley_NodeGhostEdge* Finley_NodeGhostEdge_getReference( Finley_NodeGhostEdge *in )
{
  if( in ) 
    in->reference_counter++;
  
  return in;
}

#endif
