/* created by Ben Cumming on 27/04/2006 */

#include "Distribution.h"

#ifdef PASO_MPI

/*
struct Finley_ElementDistribution
{
  dim_t reference_counter;
  Paso_MPIInfo *MPIInfo;
  index_t numLocal;     
  index_t numInternal;    
  index_t numBoundary;    
  index_t numNeighbours; 
  index_t *neighbours;    
};
*/

Finley_ElementDistribution* Finley_ElementDistribution_alloc( Paso_MPIInfo *MPIInfo )
{
  Finley_ElementDistribution *out = NULL;

  out = MEMALLOC( 1, Finley_ElementDistribution );
  if (Finley_checkPtr(out)) return NULL;

  out->reference_counter = 0;

  out->numLocal = 0;
  out->numInternal = 0;
  out->numBoundary = 0;
  //out->numNeighbours = 0;
  //out->neighbours = NULL;
  
  out->reference_counter++;

  return out;
}

void Finley_ElementDistribution_dealloc( Finley_ElementDistribution* in )
{
  if( in && !(--in->reference_counter) )
  {
    Paso_MPIInfo_dealloc( in->MPIInfo );
    //MEMFREE( in->neighbours );

    MEMFREE( in );
  }
}

Finley_ElementDistribution* Finley_ElementDistribution_getReference( Finley_ElementDistribution* in )
{
  if( in ) 
    in->reference_counter++;
  
  return in;
}

#endif
