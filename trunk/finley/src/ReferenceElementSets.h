
/*****************************************************************************
*
* Copyright (c) 2003-2013 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development since 2012 by School of Earth Sciences
*
*****************************************************************************/


/****************************************************************************

  Finley: Reference element sets manage the reference elements for the full
  and reduced integration order.

*****************************************************************************/

#ifndef __FINLEY_REFERENCEELEMENTSETS_H__
#define __FINLEY_REFERENCEELEMENTSETS_H__

#include "ReferenceElements.h"

namespace finley {

struct ReferenceElementSet {
    ReferenceElement* referenceElementReducedQuadrature;
    ReferenceElement* referenceElement;
    int numNodes;
    int reference_counter;
};


ReferenceElementSet* ReferenceElementSet_alloc(ElementTypeId id, int order,
                                               int reduced_order);

void ReferenceElementSet_dealloc(ReferenceElementSet* in);

ReferenceElementSet* ReferenceElementSet_reference(ReferenceElementSet* in);
ShapeFunction* ReferenceElementSet_borrowBasisFunctions(
        ReferenceElementSet* in, bool reducedShapefunction,
        bool reducedIntegrationOrder);

ShapeFunction* ReferenceElementSet_borrowParametrization(
        ReferenceElementSet* in, bool reducedIntegrationOrder);

ReferenceElement* ReferenceElementSet_borrowReferenceElement(
        ReferenceElementSet* in, bool reducedIntegrationOrder);

inline int ReferenceElementSet_getNumNodes(const ReferenceElementSet* in)
{
    return in->numNodes;
}


} // namespace finley

#endif // __FINLEY_REFERENCEELEMENTSETS_H__

