/*
 ************************************************************
 *          Copyright 2006 by ACcESS MNRF                   *
 *                                                          *
 *              http://www.access.edu.au                    *
 *       Primary Business: Queensland, Australia            *
 *  Licensed under the Open Software License version 3.0    *
 *     http://www.opensource.org/licenses/osl-3.0.php       *
 *                                                          *
 ************************************************************
*/

/**************************************************************/

/*   Finley: Mesh: NodeFile */

/**************************************************************/

/*   Author: gross@access.edu.au */
/*   Version: $Id$ */

/**************************************************************/

#include "NodeFile.h"

/**************************************************************/

/*  allocates the node table within an node file to hold numNodes of nodes. The LinearTo mapping, if it exists, */
/*  is frees. use Finley_Mesh_setLinearMesh to create a new one. */

void Finley_NodeFile_allocTable(Finley_NodeFile* in ,dim_t numNodes) 
{
  index_t *Id2=NULL, *Tag2=NULL, *globalDegreesOfFreedom2=NULL;
  double *Coordinates2=NULL;
  dim_t n,i;
  
  /*  allocate memory: */
  Id2=MEMALLOC(numNodes,index_t);
  Coordinates2=MEMALLOC(numNodes*in->numDim,double);
  Tag2=MEMALLOC(numNodes,index_t);
  globalDegreesOfFreedom2=MEMALLOC(numNodes,index_t);
  
  /*  if fine, freeate the old table and replace by new: */
  if (Finley_checkPtr(Id2) || Finley_checkPtr(Coordinates2) || Finley_checkPtr(Tag2) || Finley_checkPtr(globalDegreesOfFreedom2) ) {
    MEMFREE(Id2);
    MEMFREE(Coordinates2);
    MEMFREE(Tag2);
    MEMFREE(globalDegreesOfFreedom2);
  } else { 
    Finley_NodeFile_freeTable(in);
    in->Id=Id2;
    in->Coordinates=Coordinates2;
    in->globalDegreesOfFreedom=globalDegreesOfFreedom2;
    in->Tag=Tag2;
    in->numNodes=numNodes;
    /* this initialization makes sure that data are located on the right processor */
    #pragma omp parallel for private(n,i) schedule(static)
    for (n=0;n<numNodes;n++) {
       for (i=0;i<in->numDim;i++) in->Coordinates[INDEX2(i,n,in->numDim)]=0.;
       in->Id[n]=-1;
       in->Tag[n]=-1;
       in->globalDegreesOfFreedom[n]=-1;
    }
  }
  return;
}

/*  frees the node table within an node file: */

void Finley_NodeFile_freeTable(Finley_NodeFile* in) {
  if (in!=NULL) {
    MEMFREE(in->Id);
    MEMFREE(in->Coordinates);
    MEMFREE(in->globalDegreesOfFreedom);
    MEMFREE(in->Tag);
    Finley_NodeMapping_free(in->nodesMapping);
    Finley_NodeMapping_free(in->reducedNodesMapping);
    Finley_NodeMapping_free(in->degreesOfFreedomMapping);
    Finley_NodeMapping_free(in->reducedDegreesOfFreedomMapping);
    in->numNodes=0;
  }
}
