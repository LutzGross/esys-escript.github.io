
/*****************************************************************************
*
* Copyright (c) 2003-2012 by University of Queensland
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

/************************************************************************************/

/*   Dudley: Mesh: NodeFile */

/*   returns the maximum and minimum node id number of nodes: */

/************************************************************************************/

#include "NodeFile.h"
#include "Util.h"

/************************************************************************************/

void Dudley_NodeFile_setGlobalIdRange(index_t * min_id, index_t * max_id, Dudley_NodeFile * in)
{
    index_t min_id_local, max_id_local;
#ifdef ESYS_MPI
    index_t global_id_range[2], id_range[2];
#endif

    min_id_local = Dudley_Util_getMinInt(1, in->numNodes, in->Id);
    max_id_local = Dudley_Util_getMaxInt(1, in->numNodes, in->Id);

#ifdef ESYS_MPI
    id_range[0] = -min_id_local;
    id_range[1] = max_id_local;
    MPI_Allreduce(id_range, global_id_range, 2, MPI_INT, MPI_MAX, in->MPIInfo->comm);
    *min_id = -global_id_range[0];
    *max_id = global_id_range[1];
#else
    *min_id = min_id_local;
    *max_id = max_id_local;
#endif
    if (*max_id < *min_id)
    {
	*max_id = 0;
	*min_id = -1;
    }
}

void Dudley_NodeFile_setIdRange(index_t * min_id, index_t * max_id, Dudley_NodeFile * in)
{
    *min_id = Dudley_Util_getMinInt(1, in->numNodes, in->Id);
    *max_id = Dudley_Util_getMaxInt(1, in->numNodes, in->Id);
    if (*max_id < *min_id)
    {
	*max_id = 0;
	*min_id = -1;
    }
}

void Dudley_NodeFile_setGlobalDOFRange(index_t * min_id, index_t * max_id, Dudley_NodeFile * in)
{
    index_t min_id_local, max_id_local;
#ifdef ESYS_MPI
    index_t global_id_range[2], id_range[2];
#endif

    min_id_local = Dudley_Util_getMinInt(1, in->numNodes, in->globalDegreesOfFreedom);
    max_id_local = Dudley_Util_getMaxInt(1, in->numNodes, in->globalDegreesOfFreedom);

#ifdef ESYS_MPI
    id_range[0] = -min_id_local;
    id_range[1] = max_id_local;
    MPI_Allreduce(id_range, global_id_range, 2, MPI_INT, MPI_MAX, in->MPIInfo->comm);
    *min_id = -global_id_range[0];
    *max_id = global_id_range[1];
#else
    *min_id = min_id_local;
    *max_id = max_id_local;
#endif
    if (*max_id < *min_id)
    {
	*max_id = 0;
	*min_id = -1;
    }
}

void Dudley_NodeFile_setDOFRange(index_t * min_id, index_t * max_id, Dudley_NodeFile * in)
{
    *min_id = Dudley_Util_getMinInt(1, in->numNodes, in->globalDegreesOfFreedom);
    *max_id = Dudley_Util_getMaxInt(1, in->numNodes, in->globalDegreesOfFreedom);
    if (*max_id < *min_id)
    {
	*max_id = 0;
	*min_id = -1;
    }
}

void Dudley_NodeFile_setReducedDOFRange(index_t * min_id, index_t * max_id, Dudley_NodeFile * in)
{
    *min_id = Dudley_Util_getFlaggedMinInt(1, in->numNodes, in->globalReducedDOFIndex, -1);
    *max_id = Dudley_Util_getFlaggedMaxInt(1, in->numNodes, in->globalReducedDOFIndex, -1);
    if (*max_id < *min_id)
    {
	*max_id = 0;
	*min_id = -1;
    }
}

index_t Dudley_NodeFile_maxGlobalDegreeOfFreedomIndex(Dudley_NodeFile * in)
{
    index_t min_id, max_id;
    Dudley_NodeFile_setGlobalDOFRange(&min_id, &max_id, in);
    return max_id;
}

index_t Dudley_NodeFile_maxGlobalReducedDegreeOfFreedomIndex(Dudley_NodeFile * in)
{
    index_t min_id, max_id;
    Dudley_NodeFile_setGlobalReducedDegreeOfFreedomRange(&min_id, &max_id, in);
    return max_id;
}

void Dudley_NodeFile_setGlobalReducedDegreeOfFreedomRange(index_t * min_id, index_t * max_id, Dudley_NodeFile * in)
{
    index_t min_id_local, max_id_local;
#ifdef ESYS_MPI
    index_t global_id_range[2], id_range[2];
#endif

    min_id_local = Dudley_Util_getFlaggedMaxInt(1, in->numNodes, in->globalReducedDOFIndex, -1);
    max_id_local = Dudley_Util_getFlaggedMinInt(1, in->numNodes, in->globalReducedDOFIndex, -1);

#ifdef ESYS_MPI
    id_range[0] = -min_id_local;
    id_range[1] = max_id_local;
    MPI_Allreduce(id_range, global_id_range, 2, MPI_INT, MPI_MAX, in->MPIInfo->comm);
    *min_id = -global_id_range[0];
    *max_id = global_id_range[1];
#else
    *min_id = min_id_local;
    *max_id = max_id_local;
#endif
    if (*max_id < *min_id)
    {
	*max_id = 0;
	*min_id = -1;
    }
}

index_t Dudley_NodeFile_maxGlobalNodeIDIndex(Dudley_NodeFile * in)
{
    index_t min_id, max_id;
    Dudley_NodeFile_setGlobalNodeIDIndexRange(&min_id, &max_id, in);
    return max_id;
}

void Dudley_NodeFile_setGlobalNodeIDIndexRange(index_t * min_id, index_t * max_id, Dudley_NodeFile * in)
{
    index_t min_id_local, max_id_local;
#ifdef ESYS_MPI
    index_t global_id_range[2], id_range[2];
#endif

    max_id_local = Dudley_Util_getMaxInt(1, in->numNodes, in->globalNodesIndex);
    min_id_local = Dudley_Util_getMinInt(1, in->numNodes, in->globalNodesIndex);

#ifdef ESYS_MPI
    id_range[0] = -min_id_local;
    id_range[1] = max_id_local;
    MPI_Allreduce(id_range, global_id_range, 2, MPI_INT, MPI_MAX, in->MPIInfo->comm);
    *min_id = -global_id_range[0];
    *max_id = global_id_range[1];
#else
    *min_id = min_id_local;
    *max_id = max_id_local;
#endif
    if (*max_id < *min_id)
    {
	*max_id = 0;
	*min_id = -1;
    }
}

index_t Dudley_NodeFile_maxGlobalReducedNodeIDIndex(Dudley_NodeFile * in)
{
    index_t min_id, max_id;
    Dudley_NodeFile_setGlobalReducedNodeIDIndexRange(&min_id, &max_id, in);
    return max_id;
}

void Dudley_NodeFile_setGlobalReducedNodeIDIndexRange(index_t * min_id, index_t * max_id, Dudley_NodeFile * in)
{
    index_t min_id_local, max_id_local;
#ifdef ESYS_MPI
    index_t global_id_range[2], id_range[2];
#endif

    max_id_local = Dudley_Util_getFlaggedMaxInt(1, in->numNodes, in->globalReducedNodesIndex, -1);
    min_id_local = Dudley_Util_getFlaggedMinInt(1, in->numNodes, in->globalReducedNodesIndex, -1);

#ifdef ESYS_MPI
    id_range[0] = -min_id_local;
    id_range[1] = max_id_local;
    MPI_Allreduce(id_range, global_id_range, 2, MPI_INT, MPI_MAX, in->MPIInfo->comm);
    *min_id = -global_id_range[0];
    *max_id = global_id_range[1];
#else
    *min_id = min_id_local;
    *max_id = max_id_local;
#endif
    if (*max_id < *min_id)
    {
	*max_id = 0;
	*min_id = -1;
    }
}
