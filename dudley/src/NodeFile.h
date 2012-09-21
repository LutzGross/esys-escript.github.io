
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

#ifndef INC_DUDLEY_NODEFILE
#define INC_DUDLEY_NODEFILE

#define MAX_numDim 3

#include "Dudley.h"
#include "NodeMapping.h"
#include "escript/DataC.h"
#include "paso/Distribution.h"
#include "paso/Coupler.h"

struct Dudley_NodeFile {
    Esys_MPIInfo *MPIInfo;	/* MPI information */

    dim_t numNodes;		/* number of nodes */
    dim_t numDim;		/* spatial dimension */
    index_t *Id;		/* Id[i] is the id number of node i. It need to be unique. */
    index_t *Tag;		/* Tag[i] is the tag of node i. */
    index_t *tagsInUse;		/* array of tags which are actually used */
    dim_t numTagsInUse;		/* number of tags used */

    index_t *globalDegreesOfFreedom;	/* globalDegreesOfFreedom[i] is the global degree of freedom assigned to node i */
    /* this index is used to consider periodic boundary conditions by assigning */
    /* the same degreesOfFreedom to the same node */
    double *Coordinates;	/* Coordinates[INDEX2(k,i,numDim)] is the k-th coordinate of the */
    /* node i. */
    index_t *globalReducedDOFIndex;	/* assigns each local node a global unique Id in a dens labeling of reduced DOF */
    /* value <0 indicates that the DOF is not used */
    index_t *globalReducedNodesIndex;	/* assigns each local node a global unique Id in a dens labeling */
    /* value <0 indicates that the DOF is not used */
    index_t *globalNodesIndex;	/* assigns each local reduced node a global unique Id in a dens labeling */

    Dudley_NodeMapping *nodesMapping;
    Dudley_NodeMapping *reducedNodesMapping;
    Dudley_NodeMapping *degreesOfFreedomMapping;
    Dudley_NodeMapping *reducedDegreesOfFreedomMapping;

    Paso_Distribution *nodesDistribution;
    Paso_Distribution *reducedNodesDistribution;
    Paso_Distribution *degreesOfFreedomDistribution;
    Paso_Distribution *reducedDegreesOfFreedomDistribution;

    Paso_Connector *degreesOfFreedomConnector;
    Paso_Connector *reducedDegreesOfFreedomConnector;

    /* these a the packed versions of Id */
    index_t *reducedNodesId;
    index_t *degreesOfFreedomId;
    index_t *reducedDegreesOfFreedomId;

    int status;			/* the status counts the updates done on the node coordinates */
    /* the value of status is increased by when the node coordinates are updated. */

};

typedef struct Dudley_NodeFile Dudley_NodeFile;


Dudley_NodeFile *Dudley_NodeFile_alloc(dim_t, Esys_MPIInfo * MPIInfo);
index_t Dudley_NodeFile_getFirstReducedNode(Dudley_NodeFile * in);
index_t Dudley_NodeFile_getLastReducedNode(Dudley_NodeFile * in);
dim_t Dudley_NodeFile_getGlobalNumReducedNodes(Dudley_NodeFile * in);
index_t *Dudley_NodeFile_borrowGlobalReducedNodesIndex(Dudley_NodeFile * in);
index_t Dudley_NodeFile_maxGlobalNodeIDIndex(Dudley_NodeFile * in);
index_t Dudley_NodeFile_maxGlobalReducedNodeIDIndex(Dudley_NodeFile * in);
index_t Dudley_NodeFile_GlobalDegreeOfFreedomIndex(Dudley_NodeFile * in);
index_t Dudley_NodeFile_GlobalReducedDegreeOfFreedomIndex(Dudley_NodeFile * in);

index_t Dudley_NodeFile_getFirstNode(Dudley_NodeFile * in);
index_t Dudley_NodeFile_getLastNode(Dudley_NodeFile * in);
dim_t Dudley_NodeFile_getGlobalNumNodes(Dudley_NodeFile * in);
index_t *Dudley_NodeFile_borrowGlobalNodesIndex(Dudley_NodeFile * in);

/* returns the number of target */
dim_t Dudley_NodeFile_getNumReducedNodes(Dudley_NodeFile * in);
dim_t Dudley_NodeFile_getNumDegreesOfFreedom(Dudley_NodeFile * in);
dim_t Dudley_NodeFile_getNumNodes(Dudley_NodeFile * in);
dim_t Dudley_NodeFile_getNumReducedDegreesOfFreedom(Dudley_NodeFile * in);

/* returns the mapping from local nodes to a target */
index_t *Dudley_NodeFile_borrowTargetReducedNodes(Dudley_NodeFile * in);
index_t *Dudley_NodeFile_borrowTargetDegreesOfFreedom(Dudley_NodeFile * in);
index_t *Dudley_NodeFile_borrowTargetNodes(Dudley_NodeFile * in);
index_t *Dudley_NodeFile_borrowTargetReducedDegreesOfFreedom(Dudley_NodeFile * in);
/* returns the mapping from target to the local nodes */
index_t *Dudley_NodeFile_borrowReducedNodesTarget(Dudley_NodeFile * in);
index_t *Dudley_NodeFile_borrowDegreesOfFreedomTarget(Dudley_NodeFile * in);
index_t *Dudley_NodeFile_borrowNodesTarget(Dudley_NodeFile * in);
index_t *Dudley_NodeFile_borrowReducedDegreesOfFreedomTarget(Dudley_NodeFile * in);

void Dudley_NodeFile_allocTable(Dudley_NodeFile *, dim_t);
void Dudley_NodeFile_free(Dudley_NodeFile *);
void Dudley_NodeFile_freeTable(Dudley_NodeFile *);
void Dudley_NodeFile_setIdGlobalRange(index_t *, index_t *, Dudley_NodeFile *);
void Dudley_NodeFile_setIdRange(index_t *, index_t *, Dudley_NodeFile *);
void Dudley_NodeFile_setDOFGlobalRange(index_t *, index_t *, Dudley_NodeFile *);
void Dudley_NodeFile_setDOFRange(index_t *, index_t *, Dudley_NodeFile *);

void Dudley_NodeFile_setGlobalDOFRange(index_t *, index_t *, Dudley_NodeFile *);
void Dudley_NodeFile_setGlobalIdRange(index_t *, index_t *, Dudley_NodeFile *);
index_t Dudley_NodeFile_maxGlobalDegreeOfFreedomIndex(Dudley_NodeFile *);
index_t Dudley_NodeFile_maxGlobalReducedDegreeOfFreedomIndex(Dudley_NodeFile *);

void Dudley_NodeFile_setReducedDOFRange(index_t *, index_t *, Dudley_NodeFile *);
dim_t Dudley_NodeFile_createDenseDOFLabeling(Dudley_NodeFile *);
dim_t Dudley_NodeFile_createDenseNodeLabeling(Dudley_NodeFile * in, index_t * node_distribution,
					      const index_t * dof_distribution);
dim_t Dudley_NodeFile_createDenseReducedNodeLabeling(Dudley_NodeFile * in, index_t * reducedNodeMask);
dim_t Dudley_NodeFile_createDenseReducedDOFLabeling(Dudley_NodeFile * in, index_t * reducedNodeMask);
void Dudley_NodeFile_assignMPIRankToDOFs(Dudley_NodeFile * in, Esys_MPI_rank * mpiRankOfDOF, index_t * distribution);
void Dudley_NodeFile_gather(index_t *, Dudley_NodeFile *, Dudley_NodeFile *);
void Dudley_NodeFile_gather_global(index_t *, Dudley_NodeFile *, Dudley_NodeFile *);
void Dudley_NodeFile_gatherEntries(dim_t, index_t *, index_t, index_t, index_t *, index_t *, index_t *, index_t *,
				   index_t *, index_t *, dim_t numDim, double *, double *);
void Dudley_NodeFile_copyTable(dim_t, Dudley_NodeFile *, dim_t, dim_t, Dudley_NodeFile *);
void Dudley_NodeFile_scatter(index_t *, Dudley_NodeFile *, Dudley_NodeFile *);
void Dudley_NodeFile_scatterEntries(dim_t, index_t *, index_t, index_t, index_t *, index_t *, index_t *, index_t *,
				    index_t *, index_t *, dim_t numDim, double *, double *);
void Dudley_NodeFile_copyTable(dim_t, Dudley_NodeFile *, dim_t, dim_t, Dudley_NodeFile *);
void Dudley_NodeFile_setGlobalReducedDegreeOfFreedomRange(index_t * min_id, index_t * max_id, Dudley_NodeFile * in);
void Dudley_NodeFile_setGlobalNodeIDIndexRange(index_t * min_id, index_t * max_id, Dudley_NodeFile * in);
void Dudley_NodeFile_setGlobalReducedNodeIDIndexRange(index_t * min_id, index_t * max_id, Dudley_NodeFile * in);

/* ===================== */
void Dudley_NodeFile_setCoordinates(Dudley_NodeFile *, escriptDataC *);
void Dudley_NodeFile_setTags(Dudley_NodeFile *, const int, escriptDataC *);
void Dudley_NodeFile_setTagsInUse(Dudley_NodeFile * in);

#endif
