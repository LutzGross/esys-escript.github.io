
/*******************************************************
*
* Copyright (c) 2003-2010 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


/***************************************************************************************************************

    Dudley: Reference elements set managing the reference elements for the full and reduced intergation order

**************************************************************************************************************/

#ifndef INC_DUDLEY_REFERENCEELEMENTSETS
#define INC_DUDLEY_REFERENCEELEMENTSETS


/**************************************************************/

#include "ReferenceElements.h"

/**************************************************************/

   
typedef struct Dudley_ReferenceElementSet {
	Dudley_ReferenceElement* referenceElementReducedQuadrature;
	Dudley_ReferenceElement* referenceElement;
	dim_t numNodes;
	index_t reference_counter;
} Dudley_ReferenceElementSet;



Dudley_ReferenceElementSet* Dudley_ReferenceElementSet_alloc(ElementTypeId id, index_t order, index_t reduced_order);
void Dudley_ReferenceElementSet_dealloc(Dudley_ReferenceElementSet* in);
Dudley_ReferenceElementSet* Dudley_ReferenceElementSet_reference(Dudley_ReferenceElementSet* in);
Dudley_ShapeFunction* Dudley_ReferenceElementSet_borrowBasisFunctions(Dudley_ReferenceElementSet* in, bool_t reducedShapefunction, bool_t reducedIntegrationOrder);
Dudley_ShapeFunction* Dudley_ReferenceElementSet_borrowParametrization(Dudley_ReferenceElementSet* in, bool_t reducedIntegrationOrder);
Dudley_ReferenceElement* Dudley_ReferenceElementSet_borrowReferenceElement(Dudley_ReferenceElementSet* in, bool_t reducedIntegrationOrder);
#define Dudley_ReferenceElementSet_getNumNodes(__IN__) ((__IN__)->numNodes)

#endif /* #ifndef INC_DUDLEY_REFERENCEELEMENTSETS */
