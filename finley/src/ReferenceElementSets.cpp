
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

#include "ReferenceElementSets.h"
#include "esysUtils/mem.h"
#include <algorithm> // std::max

namespace finley {

ReferenceElementSet* ReferenceElementSet_alloc(ElementTypeId id, int order, int reduced_order)
{
    ReferenceElementInfo* id_info=NULL;
    ShapeFunctionInfo* bf_info=NULL;
    ReferenceElementSet *out=NULL;
    id_info=ReferenceElement_getInfo(id);
    if (! Finley_noError()) return NULL;
    bf_info=ShapeFunction_getInfo(id_info->BasisFunctions);
    if (! Finley_noError()) return NULL;

    out=new ReferenceElementSet;
    if (Finley_checkPtr(out)) return NULL;
    out->reference_counter=0;
    out->referenceElement =NULL;
    out->referenceElementReducedQuadrature =NULL;

    if (Finley_noError()) {
        if (order<0) order=std::max(2*(bf_info->numOrder),0);
        out->referenceElement=ReferenceElement_alloc(id,  order);
    }
    if (Finley_noError())  {
        if (reduced_order<0) reduced_order=std::max(2*(bf_info->numOrder-1),0);
        out->referenceElementReducedQuadrature=ReferenceElement_alloc(id,  reduced_order);
    }

    if (Finley_noError()) {
        if (! (ReferenceElement_getNumNodes(out->referenceElement) == ReferenceElement_getNumNodes(out->referenceElementReducedQuadrature) ) ) {
            Finley_setError(VALUE_ERROR,"ReferenceElementSet_alloc: numNodes in referenceElement and referenceElementReducedQuadrature don't match.");
        }
    }

    if (! Finley_noError()) {
        ReferenceElementSet_dealloc(out);
        return NULL;
    } else {
        out->numNodes=ReferenceElement_getNumNodes(out->referenceElement);
        return ReferenceElementSet_reference(out);
    }
}


void ReferenceElementSet_dealloc(ReferenceElementSet* in)
{
    if (in!=NULL) {
        in->reference_counter--;
        if (in->reference_counter<1) {
            ReferenceElement_dealloc(in->referenceElement);
            ReferenceElement_dealloc(in->referenceElementReducedQuadrature);
            delete in;
        }
    }
}
ReferenceElementSet* ReferenceElementSet_reference(ReferenceElementSet* in)
{
     if (in!=NULL) ++(in->reference_counter);
     return in;
}

ReferenceElement* ReferenceElementSet_borrowReferenceElement(
        ReferenceElementSet* in, bool reducedIntegrationOrder)
{
    ReferenceElement* out=NULL;
    if (in !=NULL) {              
        if (reducedIntegrationOrder) {
            out=in->referenceElementReducedQuadrature;
        } else {
            out=in->referenceElement;
        }
    }
    return out;
}
    
ShapeFunction* ReferenceElementSet_borrowBasisFunctions(
        ReferenceElementSet* in, bool reducedShapefunction,
        bool reducedIntegrationOrder)
{
    ShapeFunction* basis=NULL;
    if (in != NULL) {    
        if (reducedShapefunction) {
            if (reducedIntegrationOrder) {
                basis=in->referenceElementReducedQuadrature->LinearBasisFunctions;
            } else {
                basis=in->referenceElement->LinearBasisFunctions;
            }
        } else {
            if (reducedIntegrationOrder) {
                basis=in->referenceElementReducedQuadrature->BasisFunctions;
            } else {
                basis=in->referenceElement->BasisFunctions;
            }
        }
    }
    return basis;
}

ShapeFunction* ReferenceElementSet_borrowParametrization(
        ReferenceElementSet* in, bool reducedIntegrationOrder)
{
    ShapeFunction* shape=NULL;
    if (in !=NULL) {        
        if (reducedIntegrationOrder) {
            shape=in->referenceElementReducedQuadrature->Parametrization;
        } else {
            shape=in->referenceElement->Parametrization;
        }
    }
    return shape;
}

} // namespace finley

