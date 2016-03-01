
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

/****************************************************************************/

/*   Dudley: ElementFile */

/*   allocates and deallocates element table                  */

/****************************************************************************/

#include "ElementFile.h"
#include "Util.h"

namespace dudley {

/*  allocates the element table within an element file to hold numElements: */

void Dudley_ElementFile_allocTable(Dudley_ElementFile* in, dim_t numElements)
{
    index_t *Id2 = NULL, *Nodes2 = NULL, *Tag2 = NULL, *Color2 = NULL;
    int *Owner2 = NULL;
    dim_t numNodes, e, i;

    /*  allocate memory: */
    numNodes = in->numNodes;
    Owner2 = new int[numElements];
    Id2 = new index_t[numElements];
    Nodes2 = new index_t[numElements * in->numNodes];
    Tag2 = new index_t[numElements];
    Color2 = new index_t[numElements];

    /*  if fine, deallocate the old table and replace by new: */
    Dudley_ElementFile_freeTable(in);
    in->Owner = Owner2;
    in->numElements = numElements;
    in->Id = Id2;
    in->Nodes = Nodes2;
    in->Tag = Tag2;
    in->Color = Color2;

    /* this initialization makes sure that data are located on the right processor */

#pragma omp parallel for private(e,i) schedule(static)
    for (e = 0; e < numElements; e++)
    {
        for (i = 0; i < numNodes; i++)
            in->Nodes[INDEX2(i, e, numNodes)] = -1;
        in->Owner[e] = -1;
        in->Id[e] = -1;
        in->Tag[e] = -1;
        in->Color[e] = -1;
    }
    in->maxColor = -1;
    in->minColor = 0;
}

void Dudley_ElementFile_setTagsInUse(Dudley_ElementFile * in)
{
    index_t *tagsInUse = NULL;
    dim_t numTagsInUse;
    if (in != NULL)
    {
        Dudley_Util_setValuesInUse(in->Tag, in->numElements, &numTagsInUse, &tagsInUse, in->MPIInfo);
        delete[] in->tagsInUse;
        in->tagsInUse = tagsInUse;
        in->numTagsInUse = numTagsInUse;
    }
}

/*  deallocates the element table within an element file: */

void Dudley_ElementFile_freeTable(Dudley_ElementFile * in)
{
    delete[] in->Owner;
    delete[] in->Id;
    delete[] in->Nodes;
    delete[] in->Tag;
    delete[] in->Color;
    delete[] in->tagsInUse;
    in->numTagsInUse = 0;
    in->numElements = 0;
    in->maxColor = -1;
    in->minColor = 0;
}

} // namespace dudley

