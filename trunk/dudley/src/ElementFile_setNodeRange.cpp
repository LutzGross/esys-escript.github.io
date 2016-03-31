
/*****************************************************************************
*
* Copyright (c) 2003-2016 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/

/****************************************************************************/
/*                                                                                            */
/*   Dudley: ElementFile                                                                      */
/*                                                                                            */
/*   returns the maximum and minimum node reference number of nodes describing the elements:; */
/*                                                                                            */
/*                                                                                            */
/****************************************************************************/

#include "ElementFile.h"
#include "Util.h"

namespace dudley {

void Dudley_ElementFile_setNodeRange(index_t* min_id, index_t* max_id, Dudley_ElementFile* in)
{
    if (in != NULL)
    {
        *min_id = Dudley_Util_getMinInt(in->numNodes, in->numElements, in->Nodes);
        *max_id = Dudley_Util_getMaxInt(in->numNodes, in->numElements, in->Nodes);
    }
    else
    {
        *min_id = escript::DataTypes::index_t_max();
        *max_id = -escript::DataTypes::index_t_max();
    }
}

} // namespace dudley

