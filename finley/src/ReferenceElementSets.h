
/*******************************************************
*
* Copyright (c) 2003-2012 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


/***************************************************************************************************************

    Finley: Reference element sets manage the reference elements for the full and reduced integration order

**************************************************************************************************************/

#ifndef INC_FINLEY_REFERENCEELEMENTSETS
#define INC_FINLEY_REFERENCEELEMENTSETS


/**************************************************************/

#include "ReferenceElements.h"

/**************************************************************/

   
typedef struct Finley_ReferenceElementSet {
	Finley_ReferenceElement* referenceElementReducedQuadrature;
	Finley_ReferenceElement* referenceElement;
	dim_t numNodes;
	index_t reference_counter;
} Finley_ReferenceElementSet;



Finley_ReferenceElementSet* Finley_ReferenceElementSet_alloc(Finley_ElementTypeId id, index_t order, index_t reduced_order);
void Finley_ReferenceElementSet_dealloc(Finley_ReferenceElementSet* in);
Finley_ReferenceElementSet* Finley_ReferenceElementSet_reference(Finley_ReferenceElementSet* in);
Finley_ShapeFunction* Finley_ReferenceElementSet_borrowBasisFunctions(Finley_ReferenceElementSet* in, bool_t reducedShapefunction, bool_t reducedIntegrationOrder);
Finley_ShapeFunction* Finley_ReferenceElementSet_borrowParametrization(Finley_ReferenceElementSet* in, bool_t reducedIntegrationOrder);
Finley_ReferenceElement* Finley_ReferenceElementSet_borrowReferenceElement(Finley_ReferenceElementSet* in, bool_t reducedIntegrationOrder);
#define Finley_ReferenceElementSet_getNumNodes(__IN__) ((__IN__)->numNodes)

#endif /* #ifndef INC_FINLEY_REFERENCEELEMENTSETS */
