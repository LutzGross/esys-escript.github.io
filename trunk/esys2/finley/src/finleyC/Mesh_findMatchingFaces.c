/**************************************************************/

/*   Finley: Mesh */

/* searches for faces in the mesh which are matching */

/**************************************************************/

/*   Copyrights by ACcESS Australia 2003/04 */
/*   Author: gross@access.edu.au */
/*   Version: $Id$ */

/**************************************************************/

#include "Common.h"
#include "Finley.h"
#include "Util.h"
#include "Mesh.h"

/**************************************************************/

static double  Finley_Mesh_lockingGridSize=0;

int Finley_Mesh_findMatchingFaces_compar(const void *arg1 , const void *arg2 ) {
   Finley_Mesh_findMatchingFaces_center *e1,*e2;
   int i,l,g;
   e1=(Finley_Mesh_findMatchingFaces_center*) arg1;
   e2=(Finley_Mesh_findMatchingFaces_center*) arg2;
   for (i=0;i<MAX_numDim;i++) {
     l= (e1->x[i]<e2->x[i]+Finley_Mesh_lockingGridSize) ? TRUE : FALSE;
     g= (e2->x[i]<e1->x[i]+Finley_Mesh_lockingGridSize) ? TRUE : FALSE;
     if (! (l && g)) {
        if (l) return -1;
        if (g) return 1;
     }
   }
   return 0;
}

void Finley_Mesh_findMatchingFaces(Finley_NodeFile *nodes, Finley_ElementFile *faces, double safety_factor,double tolerance,
                                   int* numPairs, int* elem0,int* elem1,int* matching_nodes_in_elem1) {
#define getDist(_dist_,_e0_,_i0_,_e1_,_i1_) \
      {int i;   \
      _dist_=0; \
      for (i=0;i<numDim;i++) _dist_=MAX(_dist_,ABS(X[INDEX3(i,_i0_,_e0_,numDim,NN)]-X[INDEX3(i,_i1_,_e1_,numDim,NN)])); \
      } 

#define SWAP(_i1_,_i2_) \
            {int* i;  \
              i=(_i2_); \
              (_i2_)=(_i1_); \
              (_i1_)=i; \
             }

    int e,i,i0,i1,e_0,e_1;
    double h=DBLE(HUGE_VAL),h_local,dist,*X=NULL;
    int NN=faces->ReferenceElement->Type->numNodes;
    int numDim=nodes->numDim;
    int a1[NN],a2[NN],*perm,*perm_tmp,n;
    Finley_Mesh_findMatchingFaces_center *center;

    X=(double*) TMPMEMALLOC(sizeof(double)*NN*numDim*faces->numElements);
    center=(Finley_Mesh_findMatchingFaces_center*) TMPMEMALLOC(sizeof(Finley_Mesh_findMatchingFaces_center)*faces->numElements);
    if (!(Finley_checkPtr(X) || Finley_checkPtr(center)) ) {
       /* OMP */
       for (e=0;e<faces->numElements;e++) {
            /* get the coordinates of the nodes */
            Finley_Util_Gather_double(NN,&(faces->Nodes[INDEX2(0,e,NN)]),numDim,nodes->Coordinates,&(X[INDEX3(0,0,e,numDim,NN)]));
            /* get the element center */
            center[e].refId=e;
            for (i=0;i<MAX_numDim;i++) center[e].x[i]=0;
            for (i0=0;i0<faces->ReferenceElement->Type->numNodesOnFace;i0++) {
               for (i=0;i<numDim;i++) center[e].x[i]+=X[INDEX3(i,faces->ReferenceElement->Type->faceNode[i0],e,numDim,NN)];
            }
            for (i=0;i<numDim;i++) center[e].x[i]/=faces->ReferenceElement->Type->numNodesOnFace;
            /* get the minimum distance between nodes in the element */
            for (i0=0;i0<faces->ReferenceElement->Type->numNodesOnFace;i0++) {
               for (i1=i0+1;i1<faces->ReferenceElement->Type->numNodesOnFace;i1++) {
                  getDist(h_local,e,faces->ReferenceElement->Type->faceNode[i0],e,faces->ReferenceElement->Type->faceNode[i1]);
                  h=MIN(h,h_local);
               }
            }
       }
       /* set the */
       Finley_Mesh_lockingGridSize=h*MAX(safety_factor,0);
       #ifdef Finley_TRACE
       printf("locking grid size is %e\n",Finley_Mesh_lockingGridSize);
       #endif
       /* sort the elements by center center coordinates (lexigraphical)*/
       qsort(center,faces->numElements,sizeof(Finley_Mesh_findMatchingFaces_center),Finley_Mesh_findMatchingFaces_compar);
       /* find elements with matching center */
       *numPairs=0;
       /* OMP */
       for (e=0;e<faces->numElements-1 && Finley_ErrorCode==NO_ERROR;e++) {
          if (Finley_Mesh_findMatchingFaces_compar((void*) &(center[e]),(void*) &(center[e+1]))==0) {
              e_0=center[e].refId;
              e_1=center[e+1].refId;
              elem0[*numPairs]=e_0;
              elem1[*numPairs]=e_1;
              /* now the element e_1 is rotated such that the first node in element e_0 and e_1 have the same coordinates */
              perm=a1;
              perm_tmp=a2;
              for (i=0;i<NN;i++) perm[i]=i;
              while (Finley_ErrorCode==NO_ERROR) {
                 /* if node 0 and perm[0] are the same we are ready */
                 getDist(dist,e_0,0,e_1,perm[0]);
                 if (dist<=h*tolerance) break;
                 if (faces->ReferenceElement->Type->shiftNodes[0]>=0) {
                    /* rotate the nodes */
                    for (i=0;i<NN;i++) perm_tmp[i]=perm[faces->ReferenceElement->Type->shiftNodes[i]];
                    SWAP(perm,perm_tmp);
                 }
                 /* if the permutation is back at the identity, ie. perm[0]=0, the faces don't match: */
                 if (perm[0]==0) {
                       Finley_ErrorCode=VALUE_ERROR;
                       sprintf(Finley_ErrorMsg,"couldn't match first node of element %d to touching element %d",e_0,e_1);
                 }
              }
              /* now we check if the second nodes match */
              if (Finley_ErrorCode==NO_ERROR) {
                 if (faces->ReferenceElement->Type->numNodesOnFace>1) {
                    getDist(dist,e_0,1,e_1,perm[faces->ReferenceElement->Type->faceNode[1]]);
                    /* if the second node does not match we reverse the direction of the nodes */
                    if (dist>h*tolerance) {
                          /* rotate the nodes */
                          if (faces->ReferenceElement->Type->reverseNodes[0]<0) {
                             Finley_ErrorCode=VALUE_ERROR;
                             sprintf(Finley_ErrorMsg,"couldn't match the second node of element %d to touching element %d",e_0,e_1);
                          } else {
                             for (i=0;i<NN;i++) perm_tmp[i]=perm[faces->ReferenceElement->Type->reverseNodes[i]];
                             SWAP(perm,perm_tmp);
                             getDist(dist,e_0,1,e_1,perm[faces->ReferenceElement->Type->faceNode[1]]);
                             if (dist>h*tolerance) {
                                 Finley_ErrorCode=VALUE_ERROR;
                                 sprintf(Finley_ErrorMsg,"couldn't match the second node of element %d to touching element %d",e_0,e_1);
                             }
                          }
                    }
                 }
              }
              /* we check if the rest of the face nodes match: */
              if (Finley_ErrorCode==NO_ERROR) {
                 for (i=2;i<faces->ReferenceElement->Type->numNodesOnFace;i++) {
                    n=faces->ReferenceElement->Type->faceNode[i];
                    getDist(dist,e_0,n,e_1,perm[n]);
                    if (dist>h*tolerance) {
                       Finley_ErrorCode=VALUE_ERROR;
                       sprintf(Finley_ErrorMsg,"couldn't match the %d-th node of element %d to touching element %d",i,e_0,e_1);
                       break;
                    }
                 }
              }
              /* copy over the permuted nodes of e_1 into matching_nodes_in_elem1 */
              if (Finley_ErrorCode==NO_ERROR) {
                 for (i=0;i<NN;i++)  matching_nodes_in_elem1[INDEX2(i,*numPairs,NN)]=faces->Nodes[INDEX2(perm[i],e_1,NN)];
              }
              (*numPairs)++;
          }
       }
       #ifdef Finley_TRACE
       printf("number of pairs of matching faces %d\n",*numPairs);
       #endif
    }
    /* clean up */
    TMPMEMFREE(X);
    TMPMEMFREE(center);

#undef getDist
#undef SWAP
}

/*
* $Log$
* Revision 1.3  2004/12/15 03:48:45  jgs
* *** empty log message ***
*
* Revision 1.1.1.1  2004/10/26 06:53:57  jgs
* initial import of project esys2
*
* Revision 1.2  2004/07/02 04:21:13  gross
* Finley C code has been included
*
* Revision 1.1.1.1  2004/06/24 04:00:40  johng
* Initial version of eys using boost-python.
*
*
*/

