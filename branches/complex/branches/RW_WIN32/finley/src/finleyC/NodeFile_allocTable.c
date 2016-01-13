/*
 ******************************************************************************
 *                                                                            *
 *       COPYRIGHT  ACcESS 2003,2004,2005 -  All Rights Reserved              *
 *                                                                            *
 * This software is the property of ACcESS. No part of this code              *
 * may be copied in any form or by any means without the expressed written    *
 * consent of ACcESS.  Copying, use or modification of this software          *
 * by any unauthorised person is illegal unless that person has a software    *
 * license agreement with ACcESS.                                             *
 *                                                                            *
 ******************************************************************************
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
/*  is deallocated. use Finley_Mesh_setLinearMesh to create a new one. */

void Finley_NodeFile_allocTable(Finley_NodeFile* in ,int numNodes) {
  index_t *Id2=NULL, *Tag2=NULL, *degreeOfFreedom2=NULL, *reducedDegreeOfFreedom2=NULL, *toReduced2=NULL;
  double *Coordinates2=NULL;
  dim_t n,i;
  
  /*  allocate memory: */
  
  Id2=MEMALLOC(numNodes,index_t);
  Coordinates2=MEMALLOC(numNodes*in->numDim,double);
  Tag2=MEMALLOC(numNodes,index_t);
  degreeOfFreedom2=MEMALLOC(numNodes,index_t);
  reducedDegreeOfFreedom2=MEMALLOC(numNodes,index_t);
  toReduced2=MEMALLOC(numNodes,index_t);
  
  /*  if fine, deallocate the old table and replace by new: */
  
  if (Finley_checkPtr(Id2) || Finley_checkPtr(Coordinates2) || Finley_checkPtr(Tag2) ||
      Finley_checkPtr(degreeOfFreedom2) || Finley_checkPtr(reducedDegreeOfFreedom2) || Finley_checkPtr(toReduced2) ) {
    MEMFREE(Id2);
    MEMFREE(Coordinates2);
    MEMFREE(Tag2);
    MEMFREE(degreeOfFreedom2);
    MEMFREE(reducedDegreeOfFreedom2);
    MEMFREE(toReduced2);
  } else { 
    Finley_NodeFile_deallocTable(in);
    in->Id=Id2;
    in->Coordinates=Coordinates2;
    in->Tag=Tag2;
    in->degreeOfFreedom=degreeOfFreedom2;
    in->reducedDegreeOfFreedom=reducedDegreeOfFreedom2;
    in->toReduced=toReduced2;
    in->numNodes=numNodes;
    in->numDegreesOfFreedom=numNodes;
    in->reducedNumDegreesOfFreedom=numNodes;
    in->reducedNumNodes=numNodes;
    /* this initialization makes sure that data are located on the right processor */
    #pragma omp parallel for private(n,i) schedule(static)
    for (n=0;n<numNodes;n++) {
       for (i=0;i<in->numDim;i++) in->Coordinates[INDEX2(i,n,in->numDim)]=0.;
       in->Id[n]=-1;
       in->Tag[n]=-1;
       in->degreeOfFreedom[n]=n;
       in->reducedDegreeOfFreedom[n]=n;
       in->toReduced[n]=n;
    }
  }
  return;
}

/*  deallocates the node table within an node file: */

void Finley_NodeFile_deallocTable(Finley_NodeFile* in) {
  if (in!=NULL) {
    MEMFREE(in->Id);
    MEMFREE(in->Coordinates);
    MEMFREE(in->Tag);
    MEMFREE(in->degreeOfFreedom);
    MEMFREE(in->reducedDegreeOfFreedom);
    MEMFREE(in->toReduced);
    in->numNodes=0;
    in->numDegreesOfFreedom=0;
    in->reducedNumDegreesOfFreedom=0;
    in->reducedNumNodes=0;
  }
}
/* 
* $Log$
* Revision 1.6  2005/09/15 03:44:23  jgs
* Merge of development branch dev-02 back to main trunk on 2005-09-15
*
* Revision 1.5.2.1  2005/09/07 06:26:20  gross
* the solver from finley are put into the standalone package paso now
*
* Revision 1.5  2005/07/08 04:07:55  jgs
* Merge of development branch back to main trunk on 2005-07-08
*
* Revision 1.4  2004/12/15 07:08:33  jgs
* *** empty log message ***
* Revision 1.1.1.1.2.2  2005/06/29 02:34:54  gross
* some changes towards 64 integers in finley
*
* Revision 1.1.1.1.2.1  2004/11/24 01:37:15  gross
* some changes dealing with the integer overflow in memory allocation. Finley solves 4M unknowns now
*
*
*
*/
