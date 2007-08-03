/*
 ************************************************************
 *          Copyright 2007 by ACcESS MNRF                   *
 *                                                          *
 *              http://www.access.edu.au                    *
 *       Primary Business: Queensland, Australia            *
 *  Licensed under the Open Software License version 3.0    *
 *     http://www.opensource.org/licenses/osl-3.0.php       *
 *                                                          *
 ************************************************************
*/
/*                                                                                                                     */
/* DOFMap provides a mapping from the local FEM nodes onto the local as well as global, distributed degrees of freedom */
/*                                                                                                                     */

/* Version: $Id$ */

#ifndef INC_FINLEY_DOFMAP
#define INC_FINLEY_DOFMAP

#include "paso/Paso_MPI.h"
#include "paso/Distribution.h"


struct Finley_DOFMap {
  Paso_MPIInfo *MPIInfo; 
  Paso_Distribution* distribution;  /* describes the global distribution of the DOFs */
  dim_t numNodes; /* number of FEM nodes */
  index_t* ID; /* ID[i] is the local DOF assigned to FEM node i in 0,....,numNodes-1 */
  dim_t myNumDOFs; /* number of DOFs controled by this process */
  dim_t numDOFs; /* number pf DOFs used by this process: numDOFs-myNumDOFs=number of nodes controled by other processes */
  dim_t numRemotes;  /* number of DOFs controled by other processes = numDOFs-myNumDOFs */
  index_t* remoteID; /* remoteID[i] is the id of DOF i on the processor controling it (i=0,...,numRemotes-1) */
  Paso_MPI_rank* remoteProcessor; /* remoteProcessor[i] is the processor contoling DOF i (i=0,...,numRemotes-1) */
  index_t* offsetInRemoteID; /* offsetInRemoteID[i]  points to the first input value in remoteDOFID
                                 used by processor neighbours[i]. Has length numNeighbours+1 */
  dim_t numNeighbours;        /* number of processor controling DOFs used by this process */
  Paso_MPI_rank* neighbours;  /* array of processor controlling DOFs on this process */
  dim_t reference_counter;
};
typedef struct Finley_DOFMap Finley_DOFMap;

/* some short cuts to maps */
#define Finley_DOFMap_mapNodeToLocalDOF(_in_,_i_)  (_in_)->ID[(_i_)]
#define Finley_DOFMap_mapLocalDOFToGlobal(_in_,_i_)  ( (_i_) <  (_in_)->myNumDOFs \
                                                    ? (_i_)+(_in_)->distribution->myFirstComponent \
                                                    : (_in_)->remoteID[(_i_)-(_in_)->myNumDOFs]+(_in_)->distribution->first_component[(_in_)->remoteProcessor[(_i_)-(_in_)->myNumDOFs]])
#define Finley_DOFMap_mapNodeToToGlobalDOF(_in_,_i_) Finley_DOFMap_mapLocalDOFToGlobal(_in_,Finley_DOFMap_mapNodeToLocalDOF(_in_,_i_))


Finley_DOFMap* Finley_DOFMap_alloc(dim_t numNodes, index_t* globalID, Paso_Distribution* distribution);
void Finley_DOFMap_free(Finley_DOFMap*);
Finley_DOFMap*  DOFMap_getReference(Finley_DOFMap *in );

#endif
