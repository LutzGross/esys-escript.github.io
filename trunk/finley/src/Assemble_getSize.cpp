
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
  of elements and assigns the value to each quadrature point in out.

*****************************************************************************/

#include "Assemble.h"
#include "Util.h"

namespace finley {

void Assemble_getSize(NodeFile* nodes, ElementFile* elements, escript::Data& out)
{
    Finley_resetError();

    if (!nodes || !elements)
        return;

    ReferenceElement *refElement = ReferenceElementSet_borrowReferenceElement(
                elements->referenceElementSet,
                util::hasReducedIntegrationOrder(out));

    const int numDim=nodes->numDim;
    const int numQuad=refElement->Parametrization->numQuadNodes;
    const int NN=elements->numNodes;
    const int NS=refElement->Parametrization->Type->numShapes;
    const int NVertices=refElement->Parametrization->Type->numVertices;

    // check the dimensions of out
    if (!out.numSamplesEqual(numQuad, elements->numElements)) {
        Finley_setError(TYPE_ERROR, "Assemble_getSize: illegal number of samples of out Data object");
    } else if (!out.isDataPointShapeEqual(0, &numDim)) {
        Finley_setError(TYPE_ERROR, "Assemble_getSize: illegal data point shape of out Data object");
    } else if (!out.actsExpanded()) {
        Finley_setError(TYPE_ERROR, "Assemble_getSize: expanded Data object is expected for element size.");
    }

    if (!Finley_noError())
        return;

    // now we can start
    int node_offset;
    if (out.getFunctionSpace().getTypeCode()==FINLEY_CONTACT_ELEMENTS_2) {
        node_offset=refElement->Type->offsets[1];
    } else {
        node_offset=refElement->Type->offsets[0];
    }
    const double f=pow(0.5, pow((double)(refElement->Type->numSubElements),
                1./(double)(numDim))-1);

    out.requireWrite();
#pragma omp parallel
    {
        std::vector<double> local_X(NN*numDim);
#pragma omp parallel for
        for (int e=0; e<elements->numElements; e++) {
            // gather local coordinates of nodes into
            // local_X(numDim,NN):
            util::gather(NS, &(elements->Nodes[INDEX2(node_offset,e,NN)]),
                         numDim, nodes->Coordinates, &local_X[0]);
            // calculate minimal differences:
            double max_diff=0.;
            for (int n0=0; n0<NVertices; n0++) {
                for (int n1=n0+1; n1<NVertices; n1++) {
                    double diff=0;
                    for (int i=0; i<numDim; i++) {
                        const double d=local_X[INDEX2(i,n0,numDim)]-local_X[INDEX2(i,n1,numDim)];
                        diff += d*d;
                    }
                    max_diff=std::max(max_diff,diff);
                }
            }
            max_diff=sqrt(max_diff)*f;
            // set all values to max_diff
            double *out_array=out.getSampleDataRW(e);
            for (int q=0; q<numQuad; q++)
                out_array[q]=max_diff;
        }
    } // end of parallel region
}

} // namespace finley

