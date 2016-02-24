
/*****************************************************************************
*
* Copyright (c) 2003-2016 by The University of Queensland
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

#define ESNEEDPYTHON
#include "esysUtils/first.h"

#include "Assemble.h"
#include "ElementFile.h"
#include "ShapeTable.h"

Dudley_ElementFile_Jacobeans *Dudley_ElementFile_Jacobeans_alloc(void)
{
    Dudley_ElementFile_Jacobeans *out = new Dudley_ElementFile_Jacobeans;
    out->status = DUDLEY_INITIAL_STATUS - 1;
    out->numDim = 0;
    out->numQuad = 0;
    out->numElements = 0;
    out->absD = NULL;
    out->quadweight = 0;
    out->DSDX = NULL;
    return out;
}

/************************************************************************************/

void Dudley_ElementFile_Jacobeans_dealloc(Dudley_ElementFile_Jacobeans * in)
{
    if (in != NULL)
    {
        delete[] in->DSDX;
        delete[] in->absD;
        delete in;
    }
}

/************************************************************************************/

Dudley_ElementFile_Jacobeans *Dudley_ElementFile_borrowJacobeans(Dudley_ElementFile * self, Dudley_NodeFile * nodes,
                                                                 bool reducedIntegrationOrder)
{
    Dudley_ElementFile_Jacobeans *out = NULL;

    dim_t numNodes = self->numNodes;

    if (reducedIntegrationOrder)
    {
        out = self->jacobeans_reducedQ;
    }
    else
    {
        out = self->jacobeans;
    }
    if (out->status < nodes->status)
    {
        out->numDim = nodes->numDim;
        out->numQuad = QuadNums[self->numDim][!reducedIntegrationOrder];
        out->numShapes = self->numDim + 1;
        out->numElements = self->numElements;
        if (out->DSDX == NULL)
            out->DSDX = new  double[(out->numElements) * (out->numShapes) * (out->numDim) * (out->numQuad)];
        if (out->absD == NULL)
            out->absD = new  double[out->numElements];

        /*========================== dim = 1 ============================= */
        if (out->numDim == 1)
        {
            Dudley_setError(SYSTEM_ERROR, "Dudley does not support 1D domains.");
        /*========================== dim = 2 ============================= */
        }
        else if (out->numDim == 2)
        {
            if (self->numLocalDim == 0)
            {
                Dudley_setError(SYSTEM_ERROR,
                                "Dudley_ElementFile_borrowJacobeans: 2D does not support local dimension 0.");
            }
            else if (self->numLocalDim == 1)
            {
                Dudley_Assemble_jacobeans_2D_M1D_E1D(nodes->Coordinates, out->numQuad, self->numElements, numNodes,
                                              self->Nodes, out->DSDX, out->absD, &(out->quadweight), self->Id);
            }
            else if (self->numLocalDim == 2)
            {
                Dudley_Assemble_jacobeans_2D(nodes->Coordinates, out->numQuad, self->numElements, numNodes, self->Nodes,
                                      out->DSDX, out->absD, &(out->quadweight), self->Id);
            }
            else
            {
                Dudley_setError(SYSTEM_ERROR,
                                "Dudley_ElementFile_borrowJacobeans: local dimension in a 2D domain has to be 1 or 2.");
            }
        /*========================== dim = 3 ============================= */
        }
        else if (out->numDim == 3)
        {
            if (self->numLocalDim == 0)
            {
                Dudley_setError(SYSTEM_ERROR,
                                "Dudley_ElementFile_borrowJacobeans: 3D does not support local dimension 0.");
            }
            else if (self->numLocalDim == 2)
            {
                Dudley_Assemble_jacobeans_3D_M2D_E2D(nodes->Coordinates, out->numQuad, self->numElements, numNodes,
                                              self->Nodes, out->DSDX, out->absD, &(out->quadweight), self->Id);
            }
            else if (self->numLocalDim == 3)
            {
                Dudley_Assemble_jacobeans_3D(nodes->Coordinates, out->numQuad, self->numElements, numNodes, self->Nodes,
                                      out->DSDX, out->absD, &(out->quadweight), self->Id);
            }
            else
            {
                Dudley_setError(SYSTEM_ERROR,
                                "Dudley_ElementFile_borrowJacobeans: local dimension in a 3D domain has to be 2 or 3.");
            }
        }
        else
        {
            Dudley_setError(SYSTEM_ERROR,
                            "Dudley_ElementFile_borrowJacobeans: spatial dimension has to be 1, 2 or 3.");
        }
        if (Dudley_noError())
        {
            out->status = nodes->status;
        }
        else
        {
            out = NULL;
        }
    }
    return out;
}

