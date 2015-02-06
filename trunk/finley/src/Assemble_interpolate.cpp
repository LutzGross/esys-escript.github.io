
/*****************************************************************************
*
* Copyright (c) 2003-2015 by University of Queensland
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

  Assemblage routines: interpolates nodal data in a data array onto elements
  (=integration points)

*****************************************************************************/

#include "Assemble.h"
#include "Util.h"

namespace finley {

void Assemble_interpolate(const NodeFile* nodes, const ElementFile* elements,
                          const escript::Data& data,
                          escript::Data& interpolated_data)
{
    resetError();
    if (!nodes || !elements)
        return;

    const int data_type=data.getFunctionSpace().getTypeCode();
    const bool reducedOrder = util::hasReducedIntegrationOrder(interpolated_data);
    const_ReferenceElement_ptr refElement(elements->referenceElementSet->
                                        borrowReferenceElement(reducedOrder));

    const int *resort_nodes = NULL, *map = NULL;
    int numSub = 0, numNodes = 0;
    const_ShapeFunction_ptr basis;
    int dof_offset = 0;

    if (data_type==FINLEY_NODES) {
        numSub=refElement->Type->numSubElements;
        resort_nodes=refElement->Type->subElementNodes;
        basis=refElement->BasisFunctions;
        numNodes=nodes->getNumNodes();
        map=nodes->borrowTargetNodes();
        if (interpolated_data.getFunctionSpace().getTypeCode()==FINLEY_CONTACT_ELEMENTS_2) {
            dof_offset=refElement->Type->offsets[1];
        } else {
            dof_offset=refElement->Type->offsets[0];
        }
    } else if (data_type==FINLEY_REDUCED_NODES) {
        numSub=1;
        resort_nodes=refElement->Type->linearNodes;
        basis=refElement->LinearBasisFunctions;
        numNodes=nodes->getNumReducedNodes();
        map=nodes->borrowTargetReducedNodes();
        if (interpolated_data.getFunctionSpace().getTypeCode()==FINLEY_CONTACT_ELEMENTS_2) {
            dof_offset=refElement->LinearType->offsets[1];
        } else {
            dof_offset=refElement->LinearType->offsets[0];
        }
    } else if (data_type==FINLEY_DEGREES_OF_FREEDOM) {
        if (elements->MPIInfo->size > 1) {
            setError(TYPE_ERROR,"Assemble_interpolate: for more than one processor DEGREES_OF_FREEDOM data are not accepted as input.");
            return;
        }
        numSub=refElement->Type->numSubElements;
        resort_nodes=refElement->Type->subElementNodes;
        basis=refElement->BasisFunctions;
        numNodes=nodes->getNumDegreesOfFreedom();
        map=nodes->borrowTargetDegreesOfFreedom();
        if (interpolated_data.getFunctionSpace().getTypeCode()==FINLEY_CONTACT_ELEMENTS_2) {
            dof_offset=refElement->Type->offsets[1];
        } else {
            dof_offset=refElement->Type->offsets[0];
        }
    } else if (data_type==FINLEY_REDUCED_DEGREES_OF_FREEDOM) {
        if (elements->MPIInfo->size > 1) {
            setError(TYPE_ERROR, "Assemble_interpolate: for more than one processor REDUCED_DEGREES_OF_FREEDOM data are not accepted as input.");
            return;
        }
        numSub=1;
        resort_nodes=refElement->Type->linearNodes;
        basis=refElement->LinearBasisFunctions;
        numNodes=nodes->getNumReducedDegreesOfFreedom();
        map=nodes->borrowTargetReducedDegreesOfFreedom();
        if (interpolated_data.getFunctionSpace().getTypeCode()==FINLEY_CONTACT_ELEMENTS_2) {
            dof_offset=refElement->LinearType->offsets[1];
        } else {
            dof_offset=refElement->LinearType->offsets[0];
        }
    } else {
        setError(TYPE_ERROR,"Assemble_interpolate: Cannot interpolate data");
        return;
    }

    const int numComps=data.getDataPointSize();
    const int numQuad=basis->numQuadNodes;
    const int numShapesTotal=basis->Type->numShapes*refElement->Type->numSides;
    const int NN=elements->numNodes;
    const int NS_DOF=basis->Type->numShapes;

    // check the dimensions of interpolated_data and data
    if (!interpolated_data.numSamplesEqual(numQuad*numSub, elements->numElements)) {
        setError(TYPE_ERROR, "Assemble_interpolate: illegal number of samples of output Data object");
    } else if (!data.numSamplesEqual(1,numNodes)) {
        setError(TYPE_ERROR, "Assemble_interpolate: illegal number of samples of input Data object");
    } else if (numComps != interpolated_data.getDataPointSize()) {
        setError(TYPE_ERROR, "Assemble_interpolate: number of components of input and interpolated Data do not match.");
    }  else if (!interpolated_data.actsExpanded()) {
        setError(TYPE_ERROR, "Assemble_interpolate: expanded Data object is expected for output data.");
    }

    if (noError()) {
        interpolated_data.requireWrite();
#pragma omp parallel
        {
            // allocation of work array
            std::vector<double> local_data(NS_DOF*numComps*numSub);
            const size_t numComps_size=numComps*sizeof(double);
            // open the element loop
#pragma omp for
            for (int e=0; e<elements->numElements; e++) {
                for (int isub=0; isub<numSub; isub++) {
                    for (int q=0; q<NS_DOF; q++) {
                        const int i=elements->Nodes[INDEX2(resort_nodes[INDEX2(dof_offset+q,isub,numShapesTotal)],e,NN)];
                        const double *data_array=data.getSampleDataRO(map[i]);
                        memcpy(&local_data[INDEX3(0,q,isub, numComps,NS_DOF)], data_array, numComps_size);
                    }
                }
                // calculate interpolated_data=local_data*S
                util::smallMatSetMult1(numSub, numComps, numQuad,
                      interpolated_data.getSampleDataRW(e), NS_DOF,
                      local_data, basis->S);
            } // end of element loop
        } // end of parallel region
    } // no error
}

} // namespace finley

