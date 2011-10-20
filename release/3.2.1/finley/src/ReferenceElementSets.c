
/*******************************************************
*
* Copyright (c) 2003-2011 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


/***************************************************************************************************************

    Finley: Reference elements set managing the reference elements for the full and reduced intergation order

**************************************************************************************************************/

#include "ReferenceElementSets.h"
#include "esysUtils/mem.h"

#define MAX(X,Y) ((X)>(Y)?(X):(Y))

/**************************************************************/


Finley_ReferenceElementSet* Finley_ReferenceElementSet_alloc(Finley_ElementTypeId id, index_t order, index_t reduced_order) {
        Finley_ReferenceElementInfo* id_info=NULL;
        Finley_ShapeFunctionInfo* bf_info=NULL;
	Finley_ReferenceElementSet *out=NULL;
        id_info=Finley_ReferenceElement_getInfo(id);
	if (! Finley_noError()) return NULL;
        bf_info=Finley_ShapeFunction_getInfo(id_info->BasisFunctions);
	if (! Finley_noError()) return NULL;

	out=MEMALLOC(1,Finley_ReferenceElementSet);
	if (Finley_checkPtr(out)) return NULL;
	out->reference_counter=0;
        out->referenceElement =NULL;
        out->referenceElementReducedQuadrature =NULL;

	if (Finley_noError()) {
            if (order<0) order=MAX(2*(bf_info->numOrder),0);
	    out->referenceElement                 =Finley_ReferenceElement_alloc(id,  order);
        }
        if (Finley_noError())  {
               if (reduced_order<0) reduced_order=MAX(2*(bf_info->numOrder-1),0);
	       out->referenceElementReducedQuadrature=Finley_ReferenceElement_alloc(id,  reduced_order);
        }

	if (Finley_noError()) {
	     if (! (Finley_ReferenceElement_getNumNodes(out->referenceElement) == Finley_ReferenceElement_getNumNodes(out->referenceElementReducedQuadrature) ) ) {
		Finley_setError(VALUE_ERROR,"Finley_ReferenceElementSet_alloc: numNodes in referenceElement  and referenceElementReducedQuadrature don't match.");
             }
        }

	if (! Finley_noError()) {
		Finley_ReferenceElementSet_dealloc(out);
		return NULL;
	} else {
		out->numNodes=Finley_ReferenceElement_getNumNodes(out->referenceElement);
		return Finley_ReferenceElementSet_reference(out);
	}
}

/**************************************************************/

void Finley_ReferenceElementSet_dealloc(Finley_ReferenceElementSet* in) {
	if (in!=NULL) {
		in->reference_counter--;
		if (in->reference_counter<1) {
			Finley_ReferenceElement_dealloc(in->referenceElement);
			Finley_ReferenceElement_dealloc(in->referenceElementReducedQuadrature);
			MEMFREE(in);
		}
	}
}
Finley_ReferenceElementSet* Finley_ReferenceElementSet_reference(Finley_ReferenceElementSet* in) {
     if (in!=NULL) ++(in->reference_counter);
     return in;
}

Finley_ReferenceElement* Finley_ReferenceElementSet_borrowReferenceElement(Finley_ReferenceElementSet* in, bool_t reducedIntegrationOrder) {
		Finley_ReferenceElement* out=NULL;
    if (in !=NULL) {			  
		 if (reducedIntegrationOrder) {
             out=in->referenceElementReducedQuadrature;
		  } else {
			 out=in->referenceElement;
		  }
	}
	return out;
}
	
Finley_ShapeFunction* Finley_ReferenceElementSet_borrowBasisFunctions(Finley_ReferenceElementSet* in, bool_t reducedShapefunction, bool_t reducedIntegrationOrder) {
	Finley_ShapeFunction* basis=NULL;
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

Finley_ShapeFunction* Finley_ReferenceElementSet_borrowParametrization(Finley_ReferenceElementSet* in, bool_t reducedIntegrationOrder) {
	Finley_ShapeFunction* shape=NULL;
    if (in !=NULL) {		
		  if (reducedIntegrationOrder) {
			 shape=in->referenceElementReducedQuadrature->Parametrization;
		  } else {
			  shape=in->referenceElement->Parametrization;
		  }
	}
	return shape;
}

