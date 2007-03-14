#include "Mesh.h"

/*
  takes the ids for the elements in a mesh and updates to global values. After
	sorting the elements are ordered [elements,faceElements,contactElements,points] 
	*/

void Finley_Mesh_prepareElementDistribution( Finley_Mesh *in ){
#ifdef PASO_MPI
	int i;
	
	if(!Finley_checkPtr(in) || !Finley_checkPtr(in->MPIInfo) ){
		Finley_setError( MEMORY_ERROR, "Finley_Mesh_prepareElementDistribution() : mesh not initialised" );
		return;
	}
	
	for( i=0; i<in->MPIInfo->size+1; i++ )
		in->FaceElements->elementDistribution->vtxdist[i]+=in->Elements->elementDistribution->vtxdist[in->MPIInfo->size];
	for( i=0; i<in->FaceElements->numElements; i++ )
		in->FaceElements->Id[i] += in->Elements->elementDistribution->vtxdist[in->MPIInfo->size];

	for( i=0; i<in->MPIInfo->size+1; i++ )
		in->ContactElements->elementDistribution->vtxdist[i]+=in->FaceElements->elementDistribution->vtxdist[in->MPIInfo->size];
	for( i=0; i<in->ContactElements->numElements; i++ )
		in->ContactElements->Id[i] += in->FaceElements->elementDistribution->vtxdist[in->MPIInfo->size];

	for( i=0; i<in->MPIInfo->size+1; i++ )
		in->Points->elementDistribution->vtxdist[i]+=in->ContactElements->elementDistribution->vtxdist[in->MPIInfo->size];
	for( i=0; i<in->Points->numElements; i++ )
		in->Points->Id[i] += in->ContactElements->elementDistribution->vtxdist[in->MPIInfo->size];
#endif
}
