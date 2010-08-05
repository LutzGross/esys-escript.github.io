
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

#include "ReferenceElementSets.h"

/**************************************************************/


Dudley_ReferenceElementSet* Dudley_ReferenceElementSet_alloc(ElementTypeId id, index_t order, index_t reduced_order) {
        Dudley_ReferenceElementInfo* id_info=NULL;
        Dudley_ShapeFunctionInfo* bf_info=NULL;
	Dudley_ReferenceElementSet *out=NULL;
        id_info=Dudley_ReferenceElement_getInfo(id);
	if (! Dudley_noError()) return NULL;
        bf_info=Dudley_ShapeFunction_getInfo(id_info->BasisFunctions);
	if (! Dudley_noError()) return NULL;

	out=MEMALLOC(1,Dudley_ReferenceElementSet);
	if (Dudley_checkPtr(out)) return NULL;
	out->reference_counter=0;
        out->referenceElement =NULL;
        out->referenceElementReducedQuadrature =NULL;

	if (Dudley_noError()) {
            if (order<0) order=MAX(2*(bf_info->numOrder),0);
	    out->referenceElement                 =Dudley_ReferenceElement_alloc(id,  order);
        }
        if (Dudley_noError())  {
               if (reduced_order<0) reduced_order=MAX(2*(bf_info->numOrder-1),0);
	       out->referenceElementReducedQuadrature=Dudley_ReferenceElement_alloc(id,  reduced_order);
        }

	if (Dudley_noError()) {
	     if (! (Dudley_ReferenceElement_getNumNodes(out->referenceElement) == Dudley_ReferenceElement_getNumNodes(out->referenceElementReducedQuadrature) ) ) {
		Dudley_setError(VALUE_ERROR,"Dudley_ReferenceElementSet_alloc: numNodes in referenceElement  and referenceElementReducedQuadrature don't match.");
             }
        }

	if (! Dudley_noError()) {
		Dudley_ReferenceElementSet_dealloc(out);
		return NULL;
	} else {
		out->numNodes=Dudley_ReferenceElement_getNumNodes(out->referenceElement);
		return Dudley_ReferenceElementSet_reference(out);
	}
}

/**************************************************************/

void Dudley_ReferenceElementSet_dealloc(Dudley_ReferenceElementSet* in) {
	if (in!=NULL) {
		in->reference_counter--;
		if (in->reference_counter<1) {
			Dudley_ReferenceElement_dealloc(in->referenceElement);
			Dudley_ReferenceElement_dealloc(in->referenceElementReducedQuadrature);
		}
	}
}
Dudley_ReferenceElementSet* Dudley_ReferenceElementSet_reference(Dudley_ReferenceElementSet* in) {
     if (in!=NULL) ++(in->reference_counter);
     return in;
}

Dudley_ReferenceElement* Dudley_ReferenceElementSet_borrowReferenceElement(Dudley_ReferenceElementSet* in, bool_t reducedIntegrationOrder) {
		Dudley_ReferenceElement* out=NULL;
    if (in !=NULL) {			  
		 if (reducedIntegrationOrder) {
             out=in->referenceElementReducedQuadrature;
		  } else {
			 out=in->referenceElement;
		  }
	}
	return out;
}
	
Dudley_ShapeFunction* Dudley_ReferenceElementSet_borrowBasisFunctions(Dudley_ReferenceElementSet* in, bool_t reducedShapefunction, bool_t reducedIntegrationOrder) {
	Dudley_ShapeFunction* basis=NULL;
    if (in !=NULL) {	
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

Dudley_ShapeFunction* Dudley_ReferenceElementSet_borrowParametrization(Dudley_ReferenceElementSet* in, bool_t reducedIntegrationOrder) {
	Dudley_ShapeFunction* shape=NULL;
    if (in !=NULL) {		
		  if (reducedIntegrationOrder) {
			 shape=in->referenceElementReducedQuadrature->Parametrization;
		  } else {
			  shape=in->referenceElement->Parametrization;
		  }
	}
	return shape;
}

