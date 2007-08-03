#include "Distribution.h"

#ifdef PASO_MPI

/* 
 *  allocates memory arrays for a node distribution.
 *  assumes that numLocal and numExternal are greater than zero.
 *  numNeigbhours may be zero.
 */
void Finley_NodeDistribution_allocTable( Finley_NodeDistribution *in, dim_t numLocal, dim_t numExternal, dim_t numNeighbours )
{
  index_t i;
  index_t *neighbours=NULL, *vtxdist=NULL, *indexExternal=NULL;
  Finley_NodeGhostEdge **edges=NULL;

  if( Finley_checkPtr(in) )
    return;

	neighbours = numNeighbours ? MEMALLOC( numNeighbours, index_t ) : NULL;
	edges      = numNeighbours ? MEMALLOC( numNeighbours, Finley_NodeGhostEdge* ) : NULL;
  vtxdist = MEMALLOC( in->MPIInfo->size+1, index_t ); 
 	indexExternal = numExternal ? MEMALLOC( numExternal, index_t ) : NULL;
    
  if( (numNeighbours && ( Finley_checkPtr(neighbours) || Finley_checkPtr(edges) )) || Finley_checkPtr(vtxdist) || (numExternal &&  Finley_checkPtr(indexExternal)) )
  {
    MEMFREE( neighbours );
    MEMFREE( edges );
    MEMFREE( vtxdist );
    MEMFREE( indexExternal );
  }
  else
  {
    Finley_NodeDistribution_deallocTable( in );
    in->neighbours = neighbours;
    in->edges = edges;
    in->vtxdist = vtxdist;
    in->indexExternal = indexExternal;

    /* initialise dummy values */
    for( i=0; i<numNeighbours; i++ )
    {
      neighbours[i] = -1;
      edges[i] = NULL;
    }
    for( i=0; i<numExternal; i++ )
      indexExternal[i] = -1;

    /* MPI Communication to determine global info from local info */
    MPI_Allreduce( &numLocal, &(in->numGlobal), 1, MPI_INT, MPI_SUM, in->MPIInfo->comm );
    vtxdist[0] = 0;
    MPI_Allgather( &numLocal, 1, MPI_INT, vtxdist+1, 1, MPI_INT, in->MPIInfo->comm );
    for( i=1; i<in->MPIInfo->size+1; i++ )
      vtxdist[i] += vtxdist[i-1];

    /* table dimensions */
    in->numLocal = numLocal;
    in->numInternal = 0;
    in->numExternal = numExternal;
    in->numBoundary = 0;
    in->numNeighbours = numNeighbours;
  }

  return;
}

void Finley_NodeDistribution_deallocTable( Finley_NodeDistribution *in )
{
  index_t i;

  if( in!=NULL )
  {
    MEMFREE( in->neighbours );
    MEMFREE( in->vtxdist );
    MEMFREE( in->indexExternal );
    for( i=0; i<in->numNeighbours; i++ )
      Finley_NodeGhostEdge_dealloc( in->edges[i] );
    MEMFREE( in->edges ); 

    in->numBoundary = 0;
    in->numInternal = 0;
    in->numLocal = 0;
    in->numNeighbours = 0;
    in->numExternal = 0;
    in->numGlobal = 0;
  }
}

/* TODO */
/* maybe these should be replaced by a routine to compile the edge list automatically */
void Finley_NodeGhostEdge_allocTable( Finley_NodeGhostEdge *in, dim_t numForward, dim_t numBackward )
{
  index_t *indexForward=NULL, *indexBackward=NULL;

  if( Finley_checkPtr(in) )
    return;

  if( numForward )
  {
    indexForward  = MEMALLOC( numForward,  index_t );
  }
  if( numBackward )
  {
    indexBackward = MEMALLOC( numBackward, index_t );
  }

  if( (numForward && Finley_checkPtr(indexForward)) || (numBackward && Finley_checkPtr(indexBackward)) )
  {
    MEMFREE( indexForward );
    MEMFREE( indexBackward );
  }
  else
  {
    Finley_NodeGhostEdge_deallocTable( in );
    in->indexForward = indexForward;
    in->indexBackward = indexBackward;
    in->numForward = numForward;
    in->numBackward = numBackward;
  }
  return;
}

void Finley_NodeGhostEdge_deallocTable( Finley_NodeGhostEdge *in )
{
  if( in!=NULL )
  {
    in->numForward = in->numBackward = 0;
    MEMFREE( in->indexForward);
    MEMFREE( in->indexBackward );
  }
}


#endif
