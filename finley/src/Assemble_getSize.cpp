
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

  Assemblage routines: calculates the minimum distance between two vertices
  of elements and assigns the value to each quadrature point in element_size.

*****************************************************************************/

#include "Assemble.h"
#include "Util.h"

#include <vector>

void Finley_Assemble_getSize(Finley_NodeFile* nodes,
                             Finley_ElementFile* elements,
                             escriptDataC* element_size)
{
    Finley_resetError();

    if (!nodes || !elements)
        return;

    Finley_ReferenceElement *refElement =
        Finley_ReferenceElementSet_borrowReferenceElement(
                elements->referenceElementSet,
                Finley_Assemble_reducedIntegrationOrder(element_size));

    const dim_t numDim=nodes->numDim;
    const dim_t numQuad=refElement->Parametrization->numQuadNodes;
    const dim_t NN=elements->numNodes;
    const dim_t NS=refElement->Parametrization->Type->numShapes;
    const dim_t NVertices=refElement->Parametrization->Type->numVertices;

    // check the dimensions of element_size
    if (!numSamplesEqual(element_size, numQuad, elements->numElements)) {
        Finley_setError(TYPE_ERROR, "Finley_Assemble_getSize: illegal number of samples of element_size Data object");
    } else if (!isDataPointShapeEqual(element_size, 0, &numDim)) {
        Finley_setError(TYPE_ERROR, "Finley_Assemble_getSize: illegal data point shape of element_size Data object");
    } else if (!isExpanded(element_size)) {
        Finley_setError(TYPE_ERROR, "Finley_Assemble_getSize: expanded Data object is expected for element size.");
    }

    if (!Finley_noError())
        return;

    // now we can start
    int node_offset;
    if (getFunctionSpaceType(element_size)==FINLEY_CONTACT_ELEMENTS_2) {
        node_offset=refElement->Type->offsets[1];
    } else {
        node_offset=refElement->Type->offsets[0];
    }
    const double f=pow(0.5, pow((double)(refElement->Type->numSubElements),
                1./(double)(numDim))-1);

    requireWrite(element_size);
#pragma omp parallel
    {
        std::vector<double> local_X(NN*numDim);
#pragma omp parallel for
        for (int e=0; e<elements->numElements; e++) {
            // gather local coordinates of nodes into
            // local_X(numDim,NN):
            Finley_Util_Gather_double(NS,
                    &(elements->Nodes[INDEX2(node_offset,e,NN)]),
                    numDim, nodes->Coordinates, &local_X[0]);
            // calculate minimal differences:
            double max_diff=0.;
            for (dim_t n0=0; n0<NVertices; n0++) {
                for (dim_t n1=n0+1; n1<NVertices; n1++) {
                    double diff=0;
                    for (dim_t i=0; i<numDim; i++) {
                        const double d=local_X[INDEX2(i,n0,numDim)]-local_X[INDEX2(i,n1,numDim)];
                        diff += d*d;
                    }
                    max_diff=MAX(max_diff,diff);
                }
            }
            max_diff=sqrt(max_diff)*f;
            // set all values to max_diff
            double *element_size_array=getSampleDataRW(element_size,e);
            for (dim_t q=0; q<numQuad; q++)
                element_size_array[q]=max_diff;
        }
    } // end of parallel region
}

