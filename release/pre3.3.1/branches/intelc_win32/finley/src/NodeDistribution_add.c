/* created by Ben Cumming on 2/5/2006 */

#include "Distribution.h"

#ifdef PASO_MPI

/* returns position of first occurence of val in list of lenght len 
   returns -1 if val is not in list */
static index_t findInList( index_t val, index_t* list, dim_t len )
{
  index_t i;

  for( i=0; i<len; i++ )
    if( list[i]==val )
      return i;
  return -1;
}

void Finley_NodeDistribution_addForward( Finley_NodeDistribution *in, index_t domain, dim_t numForward, index_t* indexLocal  )
{
    index_t posDomain=0;

    if( !numForward )
      return;

    if( Finley_checkPtr(in) || Finley_checkPtr( indexLocal ) )
      return;
      
    /* determine if we already share DOF with this domain */
    posDomain = findInList( domain, in->neighbours, in->numNeighbours );

    /* we do */
    if( posDomain>=0 )
    {
      MEMREALLOC( in->edges[posDomain]->indexForward, in->edges[posDomain]->numForward + numForward, index_t );
      if( Finley_checkPtr(in->edges[posDomain]->indexForward) )
      {
        Finley_setError( PASO_MPI_ERROR, "Finley_NodeDistribution_addForward() : unable to realloc() memory for aditional degrees of freedom" ); 
        return;
      }
      memcpy( in->edges[posDomain]->indexForward + in->edges[posDomain]->numForward, indexLocal, sizeof(index_t)*numForward );
      in->edges[posDomain]->numForward += numForward;
    }
    /* we don't */
    else
    {
      posDomain = findInList( -1, in->neighbours, in->numNeighbours );
      if( posDomain==-1 )
      {     
        posDomain = in->numNeighbours++;
        MEMREALLOC( in->edges, in->numNeighbours, Finley_NodeGhostEdge* );
        MEMREALLOC( in->neighbours, in->numNeighbours, index_t );
        if( Finley_checkPtr( in->edges ) || Finley_checkPtr( in->neighbours ) )
          return;
        in->edges[posDomain] = NULL;
      }
      in->neighbours[posDomain] = domain;  
      in->edges[posDomain] = Finley_NodeGhostEdge_alloc();
      if( Finley_checkPtr(in->edges[posDomain]) )
      {
        Finley_setError( MEMORY_ERROR, "Finley_NodeDistribution_addForward() : invalid pointer to ghost edge" );
        return;
      }
      Finley_NodeGhostEdge_allocTable( in->edges[posDomain], numForward, 0 );
      memcpy( in->edges[posDomain]->indexForward, indexLocal, sizeof(index_t)*numForward );
    }
}

void Finley_NodeDistribution_addBackward( Finley_NodeDistribution *in, index_t domain, dim_t numBackward, index_t* indexLocal  )
{
    index_t posDomain=0;

    if( !numBackward )
      return;

    if( Finley_checkPtr(in) || Finley_checkPtr( indexLocal ) )
      return;

    /* determine if we already share DOF with this domain */
    posDomain = findInList( domain, in->neighbours, in->numNeighbours );

    /* we do */
    if( posDomain>=0 )
    {
      MEMREALLOC( in->edges[posDomain]->indexBackward, in->edges[posDomain]->numBackward + numBackward, index_t );
      if( Finley_checkPtr(in->edges[posDomain]->indexBackward) )
        return;
      memcpy( in->edges[posDomain]->indexBackward + in->edges[posDomain]->numBackward, indexLocal, sizeof(index_t)*numBackward );
      in->edges[posDomain]->numBackward += numBackward;
    }
    /* we don't */
    else
    {
      posDomain = in->numNeighbours++;

      if( posDomain==-1 )
      {     
        posDomain = in->numNeighbours++;
        MEMREALLOC( in->edges, in->numNeighbours, Finley_NodeGhostEdge* );
        MEMREALLOC( in->neighbours, in->numNeighbours, index_t );
        if( Finley_checkPtr( in->edges ) || Finley_checkPtr( in->neighbours ) )
          return;
        in->edges[posDomain] = NULL;
      }
      in->neighbours[posDomain] = domain;
      in->edges[posDomain] = Finley_NodeGhostEdge_alloc();
      if( Finley_checkPtr(in->edges[posDomain]))
      {
        Finley_setError( MEMORY_ERROR, "Finley_NodeDistribution_addForward() : invalid pointer to ghost edge" );
        return;
      }
      Finley_NodeGhostEdge_allocTable( in->edges[posDomain], 0, numBackward );
      memcpy( in->edges[posDomain]->indexBackward, indexLocal, sizeof(index_t)*numBackward );
    }
} 

#endif
