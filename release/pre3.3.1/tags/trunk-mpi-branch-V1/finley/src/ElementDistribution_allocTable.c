
#ifdef PASO_MPI
#include "Distribution.h"

/* note that numElementsThis is the number of elements assigned to this process */
void Finley_ElementDistribution_allocTable( Finley_ElementDistribution *in, dim_t numElements, dim_t numElementsThis ){
	index_t i;

	in->numLocal = numElements;
	in->numBoundary = 0;
	in->numInternal = 0;
	in->vtxdist = MEMALLOC(in->MPIInfo->size+1,dim_t);
	in->vtxdist[0] = 0;
	if( in->MPIInfo->size>1 )
		MPI_Allgather( &numElementsThis, 1, MPI_INT, in->vtxdist+1, 1, MPI_INT, in->MPIInfo->comm );
	else
		in->vtxdist[1] = numElements;
	for( i=0; i<in->MPIInfo->size; i++ )	
		in->vtxdist[i+1] += in->vtxdist[i];
}
#else
void Finley_ElementDistribution_allocTable(void){
}
#endif
