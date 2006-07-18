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

/* print the local distribution data to fid */
void Finley_NodeDistribution_print( Finley_NodeDistribution *in, FILE *fid )
{
  int n, i;

  fprintf( fid, "=======================================================================================\nNode Distribution\n---------------------\n\nInternal %d\nBoundary %d\nLocal %d\nExternal %d\nGlobal %d\n", in->numInternal, in->numBoundary, in->numLocal, in->numExternal, in->numGlobal );
  fprintf( fid, "vtxdist [ " );
  for( i=0; i<in->MPIInfo->size+1; i++ )
    fprintf( fid, "%d ", in->vtxdist[i] );
  fprintf( fid, "]\n" );
  fprintf( fid, "indexExternal [ " );
  for( i=0; i<in->numExternal; i++ )
    fprintf( fid, "%d ", in->indexExternal[i] );
  fprintf( fid, "]\n" );
  fprintf( fid, "NumNeighbours %d\n", in->numNeighbours );
  /* list each neighbour */
  for( n=0; n<in->numNeighbours; n++ )
  {
		if( in->neighbours[n]!=-1 ){
			fprintf( fid, "\nNeighbour %d - Domain %d\n=====================\n\n", n, in->neighbours[n] );
			fprintf( fid, "numForward = %d\tnumBackward = %d\n", in->edges[n]->numForward, in->edges[n]->numBackward );
			fprintf( fid, "indexForward\n[ " );
			for( i=0; i<in->edges[n]->numForward; i++ )
				fprintf( fid, "%d ", in->edges[n]->indexForward[i] );
			fprintf( fid, "]\n" );
			fprintf( fid, "indexBackward\n[ " );
			for( i=0; i<in->edges[n]->numBackward; i++ )
				fprintf( fid, "%d ", in->edges[n]->indexBackward[i] );
			fprintf( fid, "]\n=======================================================================================\n" );
		}else
			fprintf( fid, "\nNeighbour %d - Empty\n=====================\n\n", n );
  }
}

#endif
