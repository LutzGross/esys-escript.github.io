/* created on 3-7-2006 by Ben Cumming */

/* Take a distributed mesh and relable the degrees of freedom
   in the order Internal-Boundary-External. For a 2nd order mesh
	 it is assumed that degrees of freedom are ordered local-external,
	 so as to avoid ambiguity */
	 
#include "Mesh.h"
#include "Util.h"

#ifdef PASO_MPI	 
#include "Distribution.h"

void Finley_Mesh_resolveDegreeOfFreedomOrder( Finley_Mesh *in, bool_t doReduced )
{
	index_t *mask;	
	index_t i, n, iI=0, iB=0, iE=0, minDOF, maxDOF, len;
	dim_t numInternal, numBoundary, numLocal, numExternal, totalDOF;
	int status;
	Finley_NodeDistribution *distribution=NULL;

	minDOF = Finley_Util_getMinInt( 1, in->Nodes->numNodes, in->Nodes->degreeOfFreedom ); 
	maxDOF = Finley_Util_getMaxInt( 1, in->Nodes->numNodes, in->Nodes->degreeOfFreedom ); 
	len = maxDOF - 2*minDOF + 1;

	/* generate the DOF mask */
  mask = MEMALLOC( len, index_t );
  for( n=0; n<len; n++ )
    mask[n] = -1;

	for( n=0; n<in->Nodes->numNodes; n++ )
		mask[in->Nodes->degreeOfFreedom[n]-minDOF] = in->Nodes->Dom[n];

	/* determine the number of internal/boundary/external DOF */
  numInternal = numBoundary = numLocal = numExternal = 0;
	for( n=0; n<len; n++ )
		switch( mask[n] ) {
			case NODE_INTERNAL :
				numInternal++;
				break;
			case NODE_BOUNDARY :
				numBoundary++;
				break;
			case NODE_EXTERNAL :
				numExternal++;
				break;
		}
	numLocal = numInternal+numBoundary;

	/* find the new DOF labels according to the internal/boundary/external info */
	iI = iB = iE = 0;
	for( n=0; n<len; n++ )
		switch( mask[n] ){
			case NODE_INTERNAL :
				mask[n] = iI++;
				break;
			case NODE_BOUNDARY :
				mask[n] = numInternal + iB++;
				break;
			case NODE_EXTERNAL :
				mask[n] = numLocal + iE++;
				break;
		}

	/* update the NodeFile to the new DOF labels */
	for( n=0; n<in->Nodes->numNodes; n++ )
		in->Nodes->degreeOfFreedom[n] = mask[ in->Nodes->degreeOfFreedom[n]-minDOF ];

	/* update the distribution information */
	in->Nodes->degreeOfFreedomDistribution->numLocal = numLocal;
	in->Nodes->degreeOfFreedomDistribution->numInternal = numInternal;
	in->Nodes->degreeOfFreedomDistribution->numBoundary = numBoundary;
	in->Nodes->degreeOfFreedomDistribution->numExternal = numExternal;
	
	/* update the forward and backward indices for each neighbour */
	for( n=0; n<in->Nodes->degreeOfFreedomDistribution->numNeighbours; n++ ){
		for( i=0, iI=0; i<in->Nodes->degreeOfFreedomDistribution->edges[n]->numForward; i++ ){
				in->Nodes->degreeOfFreedomDistribution->edges[n]->indexForward[iI++] = mask[in->Nodes->degreeOfFreedomDistribution->edges[n]->indexForward[i]];
		}
		for( i=0, iI=0; i<in->Nodes->degreeOfFreedomDistribution->edges[n]->numBackward; i++ ){
				in->Nodes->degreeOfFreedomDistribution->edges[n]->indexBackward[iI++] = mask[in->Nodes->degreeOfFreedomDistribution->edges[n]->indexBackward[i]];
		}
	}

	/* recalculate vtxdist */
	distribution = in->Nodes->degreeOfFreedomDistribution;
	
	status = MPI_Allgather( &distribution->numLocal, 1, MPI_INT, distribution->vtxdist+1, 1, MPI_INT, in->MPIInfo->comm );
	if( status!=MPI_SUCCESS ){
		Finley_setError( PASO_MPI_ERROR, "Mesh_resolveDegreeOfFreedomOrder() : error performing global MPI operation : MPI_Allgather()" );
		goto clean;
	}

	distribution->vtxdist[0] = 0;
	for( i=2; i<in->MPIInfo->size+1; i++ )
		distribution->vtxdist[i] += distribution->vtxdist[i-1];
		
	status = MPI_Allreduce( &distribution->numLocal, &distribution->numGlobal, 1, MPI_INT, MPI_SUM, in->MPIInfo->comm );
	if( status!=MPI_SUCCESS ){
		Finley_setError( PASO_MPI_ERROR, "Mesh_resolveDegreeOfFreedomOrder() : error performing global MPI operation : MPI_Allreduce()" );
		goto clean;
	}
clean:
	/* make the global error status local on each MPI process */
	Finley_MPI_noError(in->MPIInfo);
  MEMFREE( mask);
}
#endif

