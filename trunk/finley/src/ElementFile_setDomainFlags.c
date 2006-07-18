#include "ElementFile.h"

#ifdef PASO_MPI
void Finley_ElementFile_setDomainFlags( Finley_ElementFile *in  )
{
	int i, internal, boundary;

	internal = boundary = 0;

	for( i=0; i<in->numElements; i++ )
		switch( in->Dom[i] ){
			case ELEMENT_INTERNAL :
				internal++;
				break;
			case ELEMENT_BOUNDARY :
			  boundary++;
				break;
		}
	in->elementDistribution->numInternal = internal;
	in->elementDistribution->numBoundary = boundary;
	in->elementDistribution->numLocal = internal+boundary;
}
#endif

