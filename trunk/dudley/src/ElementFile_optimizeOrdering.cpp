
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
/*                                                                                                         */
/*   Dudley: ElementFile                                                                                   */
/*                                                                                                         */
/*  reorders the elements in the element file such that the elements are stored close to the nodes         */
/*                                                                                                         */
/************************************************************************************/

#define ESNEEDPYTHON
#include "esysUtils/first.h"

#include "Util.h"
#include "ElementFile.h"

/************************************************************************************/

void Dudley_ElementFile_optimizeOrdering(Dudley_ElementFile** in)
{
    if (*in != NULL)
    {
        if ((*in)->numElements < 1)
            return;

        Dudley_Util_ValueAndIndex *item_list = NULL;
        Dudley_ElementFile *out = NULL;
        dim_t e, i, NN;
        index_t *index = NULL;
        NN = (*in)->numNodes;
        item_list = new Dudley_Util_ValueAndIndex[(*in)->numElements];
        index = new index_t[(*in)->numElements];
        out = Dudley_ElementFile_alloc((*in)->etype, (*in)->MPIInfo);
        if (Dudley_noError())
        {
            Dudley_ElementFile_allocTable(out, (*in)->numElements);
            if (Dudley_noError())
            {
#pragma omp parallel for private(e,i) schedule(static)
                for (e = 0; e < (*in)->numElements; e++)
                {
                    item_list[e].index = e;
                    item_list[e].value = (*in)->Nodes[INDEX2(0, e, NN)];
                    for (i = 1; i < NN; i++)
                        item_list[e].value = MIN(item_list[e].value, (*in)->Nodes[INDEX2(i, e, NN)]);
                }
                Dudley_Util_sortValueAndIndex((*in)->numElements, item_list);
#pragma omp parallel for private(e) schedule(static)
                for (e = 0; e < (*in)->numElements; e++)
                    index[e] = item_list[e].index;
                Dudley_ElementFile_gather(index, *in, out);
                Dudley_ElementFile_free(*in);
                *in = out;
            }
            else
            {
                Dudley_ElementFile_free(out);
            }
        }
        delete[] item_list;
        delete[] index;
    }
}

