
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

#ifndef __FINLEY_REFERENCEELEMENTSETS_H__
#define __FINLEY_REFERENCEELEMENTSETS_H__

#include "ReferenceElements.h"

namespace finley {

///  A reference element set manages the reference elements for the full
///  and reduced integration order
struct ReferenceElementSet {
    ReferenceElementSet(ElementTypeId id, int order, int reduced_order)
    {
        const ReferenceElementInfo* id_info = ReferenceElement::getInfo(id);
        const ShapeFunctionInfo* bf_info = ShapeFunction::getInfo(
                                                    id_info->BasisFunctions);
        if (!noError()) return;

        if (order<0)
            order=std::max(2*bf_info->numOrder, 0);

        referenceElement.reset(new ReferenceElement(id, order));
        if (noError()) {
            if (reduced_order<0)
                reduced_order=std::max(2*(bf_info->numOrder-1), 0);
            referenceElementReducedQuadrature.reset(
                    new ReferenceElement(id, reduced_order));
        }

        if (noError()) {
            if (referenceElement->getNumNodes() != referenceElementReducedQuadrature->getNumNodes()) {
                setError(VALUE_ERROR, "ReferenceElementSet: numNodes in referenceElement and referenceElementReducedQuadrature don't match.");
            }
        }
    }

    const_ShapeFunction_ptr borrowBasisFunctions(bool reducedShapefunction,
                                                 bool reducedIntegrationOrder) const
    {
        if (reducedShapefunction) {
            return (reducedIntegrationOrder ?
                    referenceElementReducedQuadrature->LinearBasisFunctions :
                    referenceElement->LinearBasisFunctions);
        }
        return (reducedIntegrationOrder ?
                referenceElementReducedQuadrature->BasisFunctions :
                referenceElement->BasisFunctions);
    }

    const_ShapeFunction_ptr borrowParametrization(bool reducedIntegrationOrder) const
    {
        return (reducedIntegrationOrder ?
                referenceElementReducedQuadrature->Parametrization :
                referenceElement->Parametrization);
    }

    const_ReferenceElement_ptr borrowReferenceElement(bool reducedIntOrder) const
    {
        return (reducedIntOrder ? referenceElementReducedQuadrature :
                                  referenceElement);
    }

    inline int getNumNodes() const { return referenceElement->getNumNodes(); }

    ReferenceElement_ptr referenceElementReducedQuadrature;
    ReferenceElement_ptr referenceElement;
};


typedef boost::shared_ptr<const ReferenceElementSet> const_ReferenceElementSet_ptr;


} // namespace finley

#endif // __FINLEY_REFERENCEELEMENTSETS_H__

