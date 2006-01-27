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

/*   Finley: Mesh */
/*   prepare nodes does */

/*   - creates a dense labeling of  degressOfFreedom */
/*   - creates/overwrites in->Nodes->reducedDegressOfFreedom */
/*   - creates/overwrites in->Nodes->reducedTo */

/**************************************************************/

/*   Author: gross@access.edu.au */
/*   Version: $Id$ */

/**************************************************************/

#include "Mesh.h"
#include "Util.h"

/**************************************************************/

void Finley_Mesh_prepareNodes(Finley_Mesh* in) {
  dim_t n,len;
  index_t id,max_id,min_id,*maskReducedDOF=NULL,*maskDOF=NULL,*reducedNodesMask=NULL,*index=NULL;

  max_id=Finley_Util_getMaxInt(1,in->Nodes->numNodes,in->Nodes->degreeOfFreedom);
  min_id=Finley_Util_getMinInt(1,in->Nodes->numNodes,in->Nodes->degreeOfFreedom);
  len=max_id-min_id+1;

  reducedNodesMask=TMPMEMALLOC(in->Nodes->numNodes,index_t);
  maskDOF=TMPMEMALLOC(len,index_t);
  maskReducedDOF=TMPMEMALLOC(len,index_t);
  index=TMPMEMALLOC(MAX(in->Nodes->numNodes,len),index_t);

  if  (! (Finley_checkPtr(maskDOF) || Finley_checkPtr(maskReducedDOF) 
                                        || Finley_checkPtr(reducedNodesMask) || Finley_checkPtr(index) ) ) {
      /* initialize everything */
      #pragma omp parallel
      {
         #pragma omp for private(n) schedule(static)
         for (n=0;n<in->Nodes->numNodes;n++) {
              in->Nodes->toReduced[n]=-1;
              in->Nodes->reducedDegreeOfFreedom[n]=-1;
              reducedNodesMask[n]=-1;
         }
         #pragma omp for private(n) schedule(static)
         for (n=0;n<len;n++) {
              maskDOF[n]=-1;
              maskReducedDOF[n]=-1;
         }
      }
      /* mark all nodes used by reduced elements */
      Finley_Mesh_markNodes(reducedNodesMask,0,in,TRUE);
      
      /* mark used degrees of freedom */
      /* OMP */
      for (n=0;n<in->Nodes->numNodes;n++) {
              id=in->Nodes->degreeOfFreedom[n]-min_id;
              maskDOF[id]=1;
              if (reducedNodesMask[n]>=0) maskReducedDOF[id]=1;
      }
      /* get a list of all nodes used in the reduced mesh and convert into in->Nodes->toReduced: */
      in->Nodes->reducedNumNodes=Finley_Util_packMask(in->Nodes->numNodes,reducedNodesMask,index);
      #pragma omp parallel for private(n) schedule(static)
      for (n=0;n<in->Nodes->reducedNumNodes;n++) in->Nodes->toReduced[index[n]]=n;
      /* get a list of the DOFs in the reduced mesh and convert it into reducedDegreeOfFreedom */
      in->Nodes->reducedNumDegreesOfFreedom=Finley_Util_packMask(len,maskReducedDOF,index);
      #pragma omp parallel for private(n) schedule(static)
      for (n=0;n<in->Nodes->reducedNumDegreesOfFreedom;n++) maskReducedDOF[index[n]]=n;
      /* get a list of the DOFs and convert it into degreeOfFreedom */
      in->Nodes->numDegreesOfFreedom=Finley_Util_packMask(len,maskDOF,index);
      #pragma omp parallel 
      {
          #pragma omp for private(n) schedule(static)
          for (n=0;n<in->Nodes->numDegreesOfFreedom;n++) maskDOF[index[n]]=n;
          #pragma omp for private(n,id) schedule(static)
          for (n=0;n<in->Nodes->numNodes;n++) {
                id=in->Nodes->degreeOfFreedom[n]-min_id;
                in->Nodes->degreeOfFreedom[n]=maskDOF[id];
                in->Nodes->reducedDegreeOfFreedom[n]=maskReducedDOF[id];
          }
      }
   }
   TMPMEMFREE(reducedNodesMask);
   TMPMEMFREE(maskDOF);
   TMPMEMFREE(maskReducedDOF);
   TMPMEMFREE(index);
}

/*
* $Log$
* Revision 1.6  2005/09/15 03:44:22  jgs
* Merge of development branch dev-02 back to main trunk on 2005-09-15
*
* Revision 1.5.2.1  2005/09/07 06:26:19  gross
* the solver from finley are put into the standalone package paso now
*
* Revision 1.5  2005/07/08 04:07:53  jgs
* Merge of development branch back to main trunk on 2005-07-08
*
* Revision 1.4  2004/12/15 07:08:33  jgs
* *** empty log message ***
* Revision 1.1.1.1.2.3  2005/06/29 02:34:52  gross
* some changes towards 64 integers in finley
*
* Revision 1.1.1.1.2.2  2004/11/24 01:37:14  gross
* some changes dealing with the integer overflow in memory allocation. Finley solves 4M unknowns now
*
*
*
*/

