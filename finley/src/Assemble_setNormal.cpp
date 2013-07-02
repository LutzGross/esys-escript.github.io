
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

  Assemblage routines: calculates the normal vector at quadrature points on
  face elements

*****************************************************************************/

#include "Assemble.h"
#include "Util.h"

namespace finley {

void Assemble_setNormal(NodeFile* nodes, ElementFile* elements, escript::Data& normal)
{
  Finley_resetError();
  if (!nodes || !elements)
      return;

    ReferenceElement* reference_element =
        ReferenceElementSet_borrowReferenceElement(
                elements->referenceElementSet,
                util::hasReducedIntegrationOrder(normal));
    const int NN=elements->numNodes;
    const int numDim=nodes->numDim;
    const int numQuad=reference_element->Parametrization->numQuadNodes;
    const int numDim_local=reference_element->Parametrization->Type->numDim;
    const int NS=reference_element->Parametrization->Type->numShapes;
  
    int sign, node_offset;
    if (normal.getFunctionSpace().getTypeCode()==FINLEY_CONTACT_ELEMENTS_2) {
        node_offset=reference_element->Type->offsets[1];
        sign=-1;
    } else {
        node_offset=reference_element->Type->offsets[0];
        sign=1;
    }

    // check the dimensions of normal
    if (!(numDim==numDim_local || numDim-1==numDim_local)) {
        Finley_setError(TYPE_ERROR, "Assemble_setNormal: Cannot calculate normal vector");
    } else if (!normal.numSamplesEqual(numQuad, elements->numElements)) {
        Finley_setError(TYPE_ERROR, "Assemble_setNormal: illegal number of samples of normal Data object");
    } else if (!normal.isDataPointShapeEqual(1, &numDim)) {
        Finley_setError(TYPE_ERROR, "Assemble_setNormal: illegal point data shape of normal Data object");
    } else if (!normal.actsExpanded()) {
        Finley_setError(TYPE_ERROR, "Assemble_setNormal: expanded Data object is expected for normal.");
    }
   
    if (Finley_noError()) {
        normal.requireWrite();
#pragma omp parallel
        {
            std::vector<double> local_X(NS*numDim); 
            std::vector<double> dVdv(numQuad*numDim*numDim_local);
            // open the element loop
#pragma omp for
            for (int e=0; e<elements->numElements; e++) {
                // gather local coordinates of nodes into local_X:
                util::gather(NS, &(elements->Nodes[INDEX2(node_offset,e,NN)]),
                             numDim, nodes->Coordinates, &local_X[0]);
                // calculate dVdv(i,j,q)=local_X(i,n)*DSDv(n,j,q)
                util::smallMatMult(numDim, numDim_local*numQuad, &dVdv[0], NS,
                        &local_X[0], reference_element->Parametrization->dSdv);
                double *normal_array=normal.getSampleDataRW(e);
                util::normalVector(numQuad, numDim, numDim_local, &dVdv[0],
                                   normal_array);
                for (int q=0; q<numQuad*numDim; q++)
                    normal_array[q]*=sign;
            }
        }
    }
}

} // namespace finley
