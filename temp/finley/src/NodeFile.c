
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

/**************************************************************/
/*                                                             */
/*   Finley: Mesh : NodeFile */
/*                                                             */
/*   allocates and freeates node files                      */
/*                                                             */
/**************************************************************/

#include "NodeFile.h"

/**************************************************************/

/*   allocates a node file to hold nodes */
/*   use Finley_NodeFile_allocTable to allocate the node table (Id,Coordinatess). */

Finley_NodeFile* Finley_NodeFile_alloc(dim_t numDim, Paso_MPIInfo *MPIInfo)
{
  Finley_NodeFile *out;
  
  /*  allocate the return value */
  
  out=MEMALLOC(1,Finley_NodeFile);
  if (Finley_checkPtr(out)) return NULL;
  out->numNodes=0;
  out->numDim=numDim;
  out->Id=NULL;
  out->globalDegreesOfFreedom=NULL;
  out->Tag=NULL;
  out->Coordinates=NULL;
  out->status=FINLEY_INITIAL_STATUS;

  out->nodesMapping=NULL;
  out->reducedNodesMapping=NULL;
  out->degreesOfFreedomMapping=NULL;
  out->reducedDegreesOfFreedomMapping=NULL;

  out->globalReducedDOFIndex=NULL;
  out->globalReducedNodesIndex=NULL;
  out->globalNodesIndex=NULL;
  out->reducedNodesId=NULL;
  out->degreesOfFreedomId=NULL;
  out->reducedDegreesOfFreedomId=NULL;
  out->nodesDistribution=NULL;
  out->reducedNodesDistribution=NULL;
  out->degreesOfFreedomDistribution=NULL;
  out->reducedDegreesOfFreedomDistribution=NULL;
  out->degreesOfFreedomCoupler=NULL;
  out->reducedDegreesOfFreedomCoupler=NULL;

  out->MPIInfo = Paso_MPIInfo_getReference( MPIInfo );
  return out;
}

/*  frees a node file: */

void Finley_NodeFile_free(Finley_NodeFile* in) {
  if (in!=NULL) {
     Finley_NodeFile_freeTable(in);
     Paso_MPIInfo_free( in->MPIInfo );
     MEMFREE(in);      
  }
}

index_t Finley_NodeFile_getFirstReducedNode(Finley_NodeFile* in) {
  if (in!=NULL) {
    return Paso_Distribution_getFirstComponent(in->reducedNodesDistribution);
  } else {
    return 0;
  }
}
index_t Finley_NodeFile_getLastReducedNode(Finley_NodeFile* in){
  if (in!=NULL) {
    return Paso_Distribution_getLastComponent(in->reducedNodesDistribution);
  } else {
    return 0;
  }

}

dim_t Finley_NodeFile_getGlobalNumReducedNodes(Finley_NodeFile* in){
  if (in!=NULL) {
    return Paso_Distribution_getGlobalNumComponents(in->reducedNodesDistribution);
  } else {
    return 0;
  }

}
index_t* Finley_NodeFile_borrowGlobalReducedNodesIndex(Finley_NodeFile* in){
  if (in!=NULL) {
    return in->globalReducedNodesIndex;
  } else {
    return NULL;
  }
}
index_t Finley_NodeFile_getFirstNode(Finley_NodeFile* in) {
  if (in!=NULL) {
    return Paso_Distribution_getFirstComponent(in->nodesDistribution);
  } else {
    return 0;
  }
}
index_t Finley_NodeFile_getLastNode(Finley_NodeFile* in){
  if (in!=NULL) {
    return Paso_Distribution_getLastComponent(in->nodesDistribution);
  } else {
    return 0;
  }

}
dim_t Finley_NodeFile_getGlobalNumNodes(Finley_NodeFile* in){
  if (in!=NULL) {
    return Paso_Distribution_getGlobalNumComponents(in->nodesDistribution);
  } else {
    return 0;
  }

}
index_t* Finley_NodeFile_borrowGlobalNodesIndex(Finley_NodeFile* in){
  if (in!=NULL) {
    return in->globalNodesIndex;
  } else {
    return NULL;
  }
}

dim_t Finley_NodeFile_getNumReducedNodes(Finley_NodeFile* in) {
  if (in!=NULL) {
       return in->reducedNodesMapping->numTargets;
  } else {
    return 0;
  }

}
dim_t Finley_NodeFile_getNumDegreesOfFreedom(Finley_NodeFile* in) {
  if (in!=NULL) {
      return Paso_Distribution_getMyNumComponents(in->degreesOfFreedomDistribution);
  } else {
    return 0;
  }
}
dim_t Finley_NodeFile_getNumNodes(Finley_NodeFile* in) {
  if (in!=NULL) {
        return in->nodesMapping->numNodes;
  } else {
    return 0;
  }
}
dim_t Finley_NodeFile_getNumReducedDegreesOfFreedom(Finley_NodeFile* in) {
  if (in!=NULL) {
      return Paso_Distribution_getMyNumComponents(in->reducedDegreesOfFreedomDistribution);
  } else {
    return 0;
  }
}


index_t* Finley_NodeFile_borrowTargetReducedNodes(Finley_NodeFile* in){
  if (in!=NULL) {
    return in->reducedNodesMapping->target;
  } else {
    return NULL;
  }
}

index_t* Finley_NodeFile_borrowTargetDegreesOfFreedom(Finley_NodeFile* in){
  if (in!=NULL) {
    return in->degreesOfFreedomMapping->target;
  } else {
    return NULL;
  }
}

index_t* Finley_NodeFile_borrowTargetNodes(Finley_NodeFile* in){
  if (in!=NULL) {
    return in->nodesMapping->target;
  } else {
    return NULL;
  }
}

index_t* Finley_NodeFile_borrowTargetReducedDegreesOfFreedom(Finley_NodeFile* in){
  if (in!=NULL) {
    return in->reducedDegreesOfFreedomMapping->target;
  } else {
    return NULL;
  }
}

index_t* Finley_NodeFile_borrowReducedNodesTarget(Finley_NodeFile* in){
  if (in!=NULL) {
    return in->reducedNodesMapping->map;
  } else {
    return NULL;
  }
}

index_t* Finley_NodeFile_borrowDegreesOfFreedomTarget(Finley_NodeFile* in){
  if (in!=NULL) {
    return in->degreesOfFreedomMapping->map;
  } else {
    return NULL;
  }
}

index_t* Finley_NodeFile_borrowNodesTarget(Finley_NodeFile* in){
  if (in!=NULL) {
    return in->nodesMapping->map;
  } else {
    return NULL;
  }
}

index_t* Finley_NodeFile_borrowReducedDegreesOfFreedomTarget(Finley_NodeFile* in){
  if (in!=NULL) {
    return in->reducedDegreesOfFreedomMapping->map;
  } else {
    return NULL;
  }
}


