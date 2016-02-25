
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

/************************************************************************************/

/*   Dudley: ElementFile */

/*   allocates an element file to hold elements of type id and with integration order order. */
/*   use Dudley_Mesh_allocElementTable to allocate the element table (Id,Nodes,Tag,Owner). */

/************************************************************************************/

#define ESNEEDPYTHON
#include "esysUtils/first.h"

#include "ElementFile.h"
#include "ShapeTable.h"

Dudley_ElementFile *Dudley_ElementFile_alloc(Dudley_ElementTypeId etype, esysUtils::JMPI& MPIInfo)
{
    Dudley_ElementFile *out;

    if (!Dudley_noError())
        return NULL;

    /*  allocate the return value */

    out = new Dudley_ElementFile;
    out->numElements = 0;
    out->Id = NULL;
    out->Nodes = NULL;
    out->Tag = NULL;
    out->Color = NULL;
    out->minColor = 0;
    out->maxColor = -1;
    out->jacobians = NULL;
    out->jacobians_reducedQ = NULL;

    out->Owner = NULL;
    out->numTagsInUse = 0;
    out->tagsInUse = NULL;

    out->MPIInfo = MPIInfo;

    out->jacobians = Dudley_ElementFile_Jacobians_alloc();
    out->jacobians_reducedQ = Dudley_ElementFile_Jacobians_alloc();

    if (!Dudley_noError())
    {
        Dudley_ElementFile_free(out);
        return NULL;
    }
    out->etype = etype;
    out->numDim = Dims[out->etype];
    out->numNodes = out->numDim + 1;
    out->numLocalDim = localDims[out->etype];
    out->numShapes = out->numLocalDim + 1;
    out->ename = getElementName(out->etype);
    return out;
}

/*  deallocates an element file: */

void Dudley_ElementFile_free(Dudley_ElementFile * in)
{
    if (in != NULL)
    {
        Dudley_ElementFile_freeTable(in);
        Dudley_ElementFile_Jacobians_dealloc(in->jacobians);
        Dudley_ElementFile_Jacobians_dealloc(in->jacobians_reducedQ);
        delete in;
    }
}

void Dudley_ElementFile_setElementDistribution(Dudley_ElementFile * in, dim_t * distribution)
{
    dim_t local_num_elements, e, num_elements = 0;
    int myRank;
    if (in == NULL)
    {
        distribution[0] = num_elements;
    }
    else
    {
        if (in->MPIInfo->size > 1)
        {
            num_elements = 0;
            myRank = in->MPIInfo->rank;
#pragma omp parallel private(local_num_elements)
            {
                local_num_elements = 0;
#pragma omp for private(e)
                for (e = 0; e < in->numElements; e++)
                {
                    if (in->Owner[e] == myRank)
                        local_num_elements++;
                }
#pragma omp critical
                num_elements += local_num_elements;
            }
#ifdef ESYS_MPI
            MPI_Allgather(&num_elements, 1, MPI_INT, distribution, 1, MPI_INT, in->MPIInfo->comm);
#else
            distribution[0] = num_elements;
#endif
        }
        else
        {
            distribution[0] = in->numElements;
        }
    }
}

dim_t Dudley_ElementFile_getGlobalNumElements(Dudley_ElementFile * in)
{
    dim_t size, *distribution = NULL, out, p;
    if (in == NULL)
    {
        return 0;
    }
    else
    {
        size = in->MPIInfo->size;
        distribution = new  dim_t[size];
        Dudley_ElementFile_setElementDistribution(in, distribution);
        out = 0;
        for (p = 0; p < size; ++p)
            out += distribution[p];
        delete[] distribution;
        return out;
    }
}

dim_t Dudley_ElementFile_getMyNumElements(Dudley_ElementFile * in)
{
    dim_t size, *distribution = NULL, out;
    if (in == NULL)
    {
        return 0;
    }
    else
    {
        size = in->MPIInfo->size;
        distribution = new  dim_t[size];
        Dudley_ElementFile_setElementDistribution(in, distribution);
        out = distribution[in->MPIInfo->rank];
        delete[] distribution;
        return out;
    }

}

index_t Dudley_ElementFile_getFirstElement(Dudley_ElementFile * in)
{
    dim_t size, *distribution = NULL, out, p;
    if (in == NULL)
    {
        return 0;
    }
    else
    {
        size = in->MPIInfo->size;
        distribution = new  dim_t[size];
        Dudley_ElementFile_setElementDistribution(in, distribution);
        out = 0;
        for (p = 0; p < in->MPIInfo->rank; ++p)
            out += distribution[p];
        delete[] distribution;
        return out;
    }
}

