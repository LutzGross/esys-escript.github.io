
/* $Id$ */

/*******************************************************
 *
 *           Copyright 2003-2007 by ACceSS MNRF
 *       Copyright 2007 by University of Queensland
 *
 *                http://esscc.uq.edu.au
 *        Primary Business: Queensland, Australia
 *  Licensed under the Open Software License version 3.0
 *     http://www.opensource.org/licenses/osl-3.0.php
 *
 *******************************************************/

#ifndef INC_FINLEY_NODEFILE
#define INC_FINLEY_NODEFILE

#define MAX_numDim 3

#include "Finley.h"
#include "NodeMapping.h"
#include "escript/DataC.h"
#include "paso/Distribution.h"
#include "paso/Coupler.h"

struct Finley_NodeFile {
  Paso_MPIInfo *MPIInfo;              /* MPI information */

  dim_t numNodes;                      /* number of nodes */
  dim_t numDim;                        /* spatial dimension */
  index_t *Id;                         /* Id[i] is the id number of node i. It need to be unique. */
  index_t *Tag;                        /* Tag[i] is the tag of node i. */

  index_t* globalDegreesOfFreedom;      /* globalDegreesOfFreedom[i] is the global degree of freedom assigned to node i */
                                       /* this index is used to consider periodic boundary conditions by assigning */
                                       /* the same degreesOfFreedom to the same node */
  double *Coordinates;                 /* Coordinates[INDEX2(k,i,numDim)] is the k-th coordinate of the */
                                       /* node i. */
  index_t *globalReducedDOFIndex;    /* assigns each local node a global unique Id in a dens labeling of reduced DOF*/
                                     /* value <0 indicates that the DOF is not used */
  index_t *globalReducedNodesIndex;    /* assigns each local node a global unique Id in a dens labeling */
                                     /* value <0 indicates that the DOF is not used */
  index_t *globalNodesIndex;           /* assigns each local reduced node a global unique Id in a dens labeling */


 Finley_NodeMapping *nodesMapping;
 Finley_NodeMapping *reducedNodesMapping;
 Finley_NodeMapping *degreesOfFreedomMapping;
 Finley_NodeMapping *reducedDegreesOfFreedomMapping;
 
 Paso_Distribution *nodesDistribution;
 Paso_Distribution *reducedNodesDistribution;
 Paso_Distribution *degreesOfFreedomDistribution;
 Paso_Distribution *reducedDegreesOfFreedomDistribution;

 Paso_Coupler* degreesOfFreedomCoupler;
 Paso_Coupler *reducedDegreesOfFreedomCoupler;
  
                     /* these a the packed versions of Id */
 index_t *reducedNodesId;        
 index_t *degreesOfFreedomId;
 index_t *reducedDegreesOfFreedomId;


 int status; /* the status counts the updates done on the node coordinates */
              /* the value of status is increased by when the node coordinates are updated.*/
                                                                                                                                                                                                 
                                                                                                                                                                                                 
};

typedef struct Finley_NodeFile Finley_NodeFile;



Finley_NodeFile* Finley_NodeFile_alloc(dim_t, Paso_MPIInfo *MPIInfo);
index_t Finley_NodeFile_getFirstReducedNode(Finley_NodeFile* in);
index_t Finley_NodeFile_getLastReducedNode(Finley_NodeFile* in);
dim_t Finley_NodeFile_getGlobalNumReducedNodes(Finley_NodeFile* in);
index_t* Finley_NodeFile_borrowGlobalReducedNodesIndex(Finley_NodeFile* in);
index_t Finley_NodeFile_maxGlobalNodeIDIndex(Finley_NodeFile* in);
index_t Finley_NodeFile_maxGlobalReducedNodeIDIndex(Finley_NodeFile* in);
index_t Finley_NodeFile_GlobalDegreeOfFreedomIndex(Finley_NodeFile* in);
index_t Finley_NodeFile_GlobalReducedDegreeOfFreedomIndex(Finley_NodeFile* in);

index_t Finley_NodeFile_getFirstNode(Finley_NodeFile* in);
index_t Finley_NodeFile_getLastNode(Finley_NodeFile* in);
dim_t Finley_NodeFile_getGlobalNumNodes(Finley_NodeFile* in);
index_t* Finley_NodeFile_borrowGlobalNodesIndex(Finley_NodeFile* in);

/* returns the number of target */
dim_t Finley_NodeFile_getNumReducedNodes(Finley_NodeFile* in);
dim_t Finley_NodeFile_getNumDegreesOfFreedom(Finley_NodeFile* in);
dim_t Finley_NodeFile_getNumNodes(Finley_NodeFile* in);
dim_t Finley_NodeFile_getNumReducedDegreesOfFreedom(Finley_NodeFile* in);

/* returns the mapping from local nodes to a target */
index_t* Finley_NodeFile_borrowTargetReducedNodes(Finley_NodeFile* in);
index_t* Finley_NodeFile_borrowTargetDegreesOfFreedom(Finley_NodeFile* in);
index_t* Finley_NodeFile_borrowTargetNodes(Finley_NodeFile* in);
index_t* Finley_NodeFile_borrowTargetReducedDegreesOfFreedom(Finley_NodeFile* in);
/* returns the mapping from target to the local nodes */
index_t* Finley_NodeFile_borrowReducedNodesTarget(Finley_NodeFile* in);
index_t* Finley_NodeFile_borrowDegreesOfFreedomTarget(Finley_NodeFile* in);
index_t* Finley_NodeFile_borrowNodesTarget(Finley_NodeFile* in);
index_t* Finley_NodeFile_borrowReducedDegreesOfFreedomTarget(Finley_NodeFile* in);

void Finley_NodeFile_allocTable(Finley_NodeFile*,dim_t);
void Finley_NodeFile_free(Finley_NodeFile*);
void Finley_NodeFile_freeTable(Finley_NodeFile*);
void Finley_NodeFile_setIdGlobalRange(index_t*,index_t*,Finley_NodeFile*);
void Finley_NodeFile_setIdRange(index_t*,index_t*,Finley_NodeFile*);
void Finley_NodeFile_setDOFGlobalRange(index_t*,index_t*,Finley_NodeFile*);
void Finley_NodeFile_setDOFRange(index_t*,index_t*,Finley_NodeFile*);

void Finley_NodeFile_setGlobalDOFRange(index_t*,index_t*,Finley_NodeFile*);
void Finley_NodeFile_setGlobalIdRange(index_t*,index_t*,Finley_NodeFile*);
index_t Finley_NodeFile_maxGlobalDegreeOfFreedomIndex(Finley_NodeFile*);
index_t Finley_NodeFile_maxGlobalReducedDegreeOfFreedomIndex(Finley_NodeFile*);

void Finley_NodeFile_setReducedDOFRange(index_t*,index_t*,Finley_NodeFile*);
dim_t Finley_NodeFile_createDenseDOFLabeling(Finley_NodeFile*);
dim_t Finley_NodeFile_createDenseNodeLabeling(Finley_NodeFile* in);
dim_t Finley_NodeFile_createDenseReducedNodeLabeling(Finley_NodeFile* in, index_t* reducedNodeMask);
dim_t Finley_NodeFile_createDenseReducedDOFLabeling(Finley_NodeFile* in, index_t* reducedNodeMask);
void Finley_NodeFile_assignMPIRankToDOFs(Finley_NodeFile* in,Paso_MPI_rank* mpiRankOfDOF, index_t *distribution);
void Finley_NodeFile_gather(index_t*,Finley_NodeFile*,Finley_NodeFile*);
void Finley_NodeFile_gather_global(index_t*,Finley_NodeFile*,Finley_NodeFile*);
void Finley_NodeFile_gatherEntries(dim_t, index_t*, index_t, index_t, index_t*, index_t*, index_t*, index_t*, index_t*, index_t*, dim_t numDim, double*, double*);
void Finley_NodeFile_copyTable(dim_t,Finley_NodeFile*,dim_t,dim_t,Finley_NodeFile*);
void Finley_NodeFile_scatter(index_t*,Finley_NodeFile*,Finley_NodeFile*);
void Finley_NodeFile_scatterEntries(dim_t, index_t*, index_t, index_t, index_t*, index_t*, index_t*, index_t*, index_t*, index_t*, dim_t numDim, double*, double*);
void Finley_NodeFile_copyTable(dim_t,Finley_NodeFile*,dim_t,dim_t,Finley_NodeFile*);

/* ===================== */
void Finley_NodeFile_setCoordinates(Finley_NodeFile*,escriptDataC*);
void Finley_NodeFile_setTags(Finley_NodeFile*,const int,escriptDataC*);

#endif

