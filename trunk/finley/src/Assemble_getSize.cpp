
/*****************************************************************************
*
* Copyright (c) 2003-2015 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/


/****************************************************************************

  Assemblage routines: calculates the minimum distance between two vertices
  of elements and assigns the value to each quadrature point in out.

*****************************************************************************/
#define ESNEEDPYTHON
#include "esysUtils/first.h"

#include "Assemble.h"
#include "Util.h"

namespace finley {

void Assemble_getSize(const NodeFile* nodes, const ElementFile* elements,
                      escript::Data& out)
{
    resetError();

    if (!nodes || !elements)
        return;

    const_ReferenceElement_ptr refElement(elements->referenceElementSet->
            borrowReferenceElement(
                util::hasReducedIntegrationOrder(out)));

    const int numDim=nodes->numDim;
    const int numQuad=refElement->Parametrization->numQuadNodes;
    const int NN=elements->numNodes;
    const int NS=refElement->Parametrization->Type->numShapes;
    const int NVertices=refElement->Parametrization->Type->numVertices;

    // check the dimensions of out
    if (!out.numSamplesEqual(numQuad, elements->numElements)) {
        setError(TYPE_ERROR, "Assemble_getSize: illegal number of samples of out Data object");
    } else if (!out.isDataPointShapeEqual(0, &numDim)) {
        setError(TYPE_ERROR, "Assemble_getSize: illegal data point shape of out Data object");
    } else if (!out.actsExpanded()) {
        setError(TYPE_ERROR, "Assemble_getSize: expanded Data object is expected for element size.");
    }

    if (!noError())
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

