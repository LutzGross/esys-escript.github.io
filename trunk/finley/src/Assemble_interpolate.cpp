
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

  Assemblage routines: interpolates nodal data in a data array onto elements
  (=integration points)

*****************************************************************************/

#include "Assemble.h"
#include "Util.h"

#include <vector>

void Finley_Assemble_interpolate(finley::NodeFile *nodes,
                                 Finley_ElementFile *elements,
                                 escriptDataC *data,
                                 escriptDataC *interpolated_data)
{
    Finley_resetError();
    if (!nodes || !elements)
        return;

    const type_t data_type=getFunctionSpaceType(data);
    const bool_t reduced_integration = Finley_Assemble_reducedIntegrationOrder(interpolated_data);
    Finley_ReferenceElement *reference_element =
        Finley_ReferenceElementSet_borrowReferenceElement(
                elements->referenceElementSet, reduced_integration);

    int *resort_nodes = NULL, *map = NULL;
    int numSub = 0, numNodes = 0;
    Finley_ShapeFunction *basis = NULL;
    int dof_offset = 0;

    if (data_type==FINLEY_NODES) {
        numSub=reference_element->Type->numSubElements;
        resort_nodes=reference_element->Type->subElementNodes;
        basis=reference_element->BasisFunctions;
        numNodes=nodes->getNumNodes();
        map=nodes->borrowTargetNodes();
        if (getFunctionSpaceType(interpolated_data)==FINLEY_CONTACT_ELEMENTS_2) {
            dof_offset=reference_element->Type->offsets[1];
        } else {
            dof_offset=reference_element->Type->offsets[0];
        }
    } else if (data_type==FINLEY_REDUCED_NODES) {
        numSub=1;
        resort_nodes=reference_element->Type->linearNodes;
        basis=reference_element->LinearBasisFunctions;
        numNodes=nodes->getNumReducedNodes();
        map=nodes->borrowTargetReducedNodes();
        if (getFunctionSpaceType(interpolated_data)==FINLEY_CONTACT_ELEMENTS_2) {
            dof_offset=reference_element->LinearType->offsets[1];
        } else {
            dof_offset=reference_element->LinearType->offsets[0];
        }
    } else if (data_type==FINLEY_DEGREES_OF_FREEDOM) {
        if (elements->MPIInfo->size > 1) {
            Finley_setError(TYPE_ERROR,"Finley_Assemble_interpolate: for more than one processor DEGREES_OF_FREEDOM data are not accepted as input.");
            return;
        }
        numSub=reference_element->Type->numSubElements;
        resort_nodes=reference_element->Type->subElementNodes;
        basis=reference_element->BasisFunctions;
        numNodes=nodes->getNumDegreesOfFreedom();
        map=nodes->borrowTargetDegreesOfFreedom();
        if (getFunctionSpaceType(interpolated_data)==FINLEY_CONTACT_ELEMENTS_2) {
            dof_offset=reference_element->Type->offsets[1];
        } else {
            dof_offset=reference_element->Type->offsets[0];
        }
    } else if (data_type==FINLEY_REDUCED_DEGREES_OF_FREEDOM) {
        if (elements->MPIInfo->size > 1) {
            Finley_setError(TYPE_ERROR, "Finley_Assemble_interpolate: for more than one processor REDUCED_DEGREES_OF_FREEDOM data are not accepted as input.");
            return;
        }
        numSub=1;
        resort_nodes=reference_element->Type->linearNodes;
        basis=reference_element->LinearBasisFunctions;
        numNodes=nodes->getNumReducedDegreesOfFreedom();
        map=nodes->borrowTargetReducedDegreesOfFreedom();
        if (getFunctionSpaceType(interpolated_data)==FINLEY_CONTACT_ELEMENTS_2) {
            dof_offset=reference_element->LinearType->offsets[1];
        } else {
            dof_offset=reference_element->LinearType->offsets[0];
        }
    } else {
        Finley_setError(TYPE_ERROR,"Finley_Assemble_interpolate: Cannot interpolate data");
        return;
    }

    const int numComps=getDataPointSize(data);
    const int numQuad=basis->numQuadNodes;
    const int numShapesTotal=basis->Type->numShapes*reference_element->Type->numSides;
    const int NN=elements->numNodes;
    const int NS_DOF=basis->Type->numShapes;

    // check the dimensions of interpolated_data and data
    if (!numSamplesEqual(interpolated_data, numQuad*numSub, elements->numElements)) {
        Finley_setError(TYPE_ERROR, "Finley_Assemble_interpolate: illegal number of samples of output Data object");
    } else if (! numSamplesEqual(data,1,numNodes)) {
        Finley_setError(TYPE_ERROR, "Finley_Assemble_interpolate: illegal number of samples of input Data object");
    } else if (numComps!=getDataPointSize(interpolated_data)) {
        Finley_setError(TYPE_ERROR, "Finley_Assemble_interpolate: number of components of input and interpolated Data do not match.");
    }  else if (!isExpanded(interpolated_data)) {
        Finley_setError(TYPE_ERROR, "Finley_Assemble_interpolate: expanded Data object is expected for output data.");
    }

    // now we can start
    if (Finley_noError()) {
        requireWrite(interpolated_data);
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
                        const double *data_array=getSampleDataRO(data, map[i]);
                        memcpy(&(local_data[INDEX3(0,q,isub, numComps,NS_DOF)]), data_array, numComps_size);
                    }
                }
                // calculate interpolated_data=local_data*S
                Finley_Util_SmallMatSetMult1(numSub, numComps, numQuad,
                      getSampleDataRW(interpolated_data,e),
                      NS_DOF, &local_data[0], basis->S);
            } // end of element loop
        } // end of parallel region
    } // no error
}

