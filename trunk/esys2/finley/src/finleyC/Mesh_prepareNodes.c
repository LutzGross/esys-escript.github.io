/**************************************************************/

/*   Finley: Mesh */
/*   prepare nodes does */

/*   - creates a dense labeling of  degressOfFreedom */
/*   - creates/overwrites in->Nodes->reducedDegressOfFreedom */
/*   - creates/overwrites in->Nodes->reducedTo */

/**************************************************************/

/*   Copyrights by ACcESS Australia 2003 */
/*   Author: gross@access.edu.au */
/*   Version: $Id$ */

/**************************************************************/

#include "Finley.h"
#include "Mesh.h"
#include "Util.h"

/**************************************************************/

void Finley_Mesh_prepareNodes(Finley_Mesh* in) {
  int n,id,max_id,min_id,len;
  maybelong *maskReducedDOF=NULL,*maskDOF=NULL,*reducedNodesMask=NULL,*index=NULL;

  max_id=Finley_Util_getMaxInt(1,in->Nodes->numNodes,in->Nodes->degreeOfFreedom);
  min_id=Finley_Util_getMinInt(1,in->Nodes->numNodes,in->Nodes->degreeOfFreedom);
  len=max_id-min_id+1;

  reducedNodesMask=TMPMEMALLOC(in->Nodes->numNodes,maybelong);
  maskDOF=TMPMEMALLOC(len,maybelong);
  maskReducedDOF=TMPMEMALLOC(len,maybelong);
  index=TMPMEMALLOC(MAX(in->Nodes->numNodes,len),maybelong);

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
* Revision 1.2  2004/12/14 05:39:30  jgs
* *** empty log message ***
*
* Revision 1.1.1.1.2.2  2004/11/24 01:37:14  gross
* some changes dealing with the integer overflow in memory allocation. Finley solves 4M unknowns now
*
* Revision 1.1.1.1.2.1  2004/11/12 06:58:18  gross
* a lot of changes to get the linearPDE class running: most important change is that there is no matrix format exposed to the user anymore. the format is chosen by the Domain according to the solver and symmetry
*
* Revision 1.1.1.1  2004/10/26 06:53:57  jgs
* initial import of project esys2
*
* Revision 1.2  2004/07/30 04:37:06  gross
* escript and finley are linking now and RecMeshTest.py has been passed
*
* Revision 1.1.1.1  2004/06/24 04:00:40  johng
* Initial version of eys using boost-python.
*
*
*/

