
/*****************************************************************************
*
* Copyright (c) 2003-2013 by University of Queensland
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

/*   Finley: Mesh */

/* searches for faces in the mesh which are matching */

/************************************************************************************/

#include "Util.h"
#include "Mesh.h"

/************************************************************************************/

using namespace finley;

static double Finley_Mesh_lockingGridSize=0;

int Finley_Mesh_findMatchingFaces_compar(const void *arg1 , const void *arg2 ) {
   Finley_Mesh_findMatchingFaces_center *e1,*e2;
   bool_t l,g;
   dim_t i;
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
   if ( e1->refId < e2->refId ) {
       return -1;
   } else if ( e1->refId > e2->refId ) {
       return 1;
   } else {
       return 0;
   }
}

void Finley_Mesh_findMatchingFaces(NodeFile *nodes, ElementFile *faces,
                                   double safety_factor,double tolerance,
                                   int* numPairs, int* elem0, int* elem1,
                                   int* matching_nodes_in_elem1)
{
#define getDist(_dist_,_e0_,_i0_,_e1_,_i1_) \
  {\
   _dist_=0;\
   for (int i=0;i<numDim;i++) _dist_=MAX(_dist_,ABS(X[INDEX3(i,_i0_,_e0_,numDim,NN)]-X[INDEX3(i,_i1_,_e1_,numDim,NN)]));\
  }
    Finley_ReferenceElement* refElement=NULL;
    char error_msg[LenErrorMsg_MAX];
    double h=DBLE(HUGE_VAL),h_local,dist,*X=NULL;
    Finley_Mesh_findMatchingFaces_center *center;
    index_t e_0,e_1,*a1=NULL,*a2=NULL,*perm=NULL,*perm_tmp=NULL,*itmp_ptr=NULL, *faceNodes=NULL, *shiftNodes=NULL, *reverseNodes=NULL ;
    dim_t e,i,i0,i1,n,NN,numNodesOnFace;

    dim_t numDim=nodes->numDim;

    refElement= Finley_ReferenceElementSet_borrowReferenceElement(faces->referenceElementSet, FALSE);
    NN=faces->numNodes;

    numNodesOnFace=refElement->Type->numNodesOnFace;
    faceNodes=refElement->Type->faceNodes;
    shiftNodes=refElement->Type->shiftNodes;
    reverseNodes=refElement->Type->reverseNodes;

    if (numNodesOnFace<=0) {
        sprintf(error_msg,"Finley_Mesh_findMatchingFaces: matching faces cannot be applied to face elements of type %s",refElement->Type->Name);
        Finley_setError(TYPE_ERROR,error_msg);
        return;
    }
    X=new double[NN*numDim*faces->numElements];
    center=new Finley_Mesh_findMatchingFaces_center[faces->numElements];
    a1=new int[NN];
    a2=new int[NN];
    /* TODO: OMP */
    for (e=0;e<faces->numElements;e++) {
        /* get the coordinates of the nodes */
        util::gather(NN, &(faces->Nodes[INDEX2(0,e,NN)]), numDim,
                     nodes->Coordinates, &(X[INDEX3(0,0,e,numDim,NN)]));
        /* get the element center */
        center[e].refId=e;
        for (i=0;i<MAX_numDim;i++) center[e].x[i]=0;
        for (i0=0;i0<numNodesOnFace;i0++) {
            for (i=0;i<numDim;i++) center[e].x[i]+=X[INDEX3(i,faceNodes[i0],e,numDim,NN)];
        }
        for (i=0;i<numDim;i++) center[e].x[i]/=numNodesOnFace;
        /* get the minimum distance between nodes in the element */
        for (i0=0;i0<numNodesOnFace;i0++) {
            for (i1=i0+1;i1<numNodesOnFace;i1++) {
                getDist(h_local,e,faceNodes[i0],e,faceNodes[i1]);
                h=MIN(h,h_local);
            }
        }
    }
    Finley_Mesh_lockingGridSize=h*MAX(safety_factor,0);
#ifdef Finley_TRACE
    printf("locking grid size is %e\n",Finley_Mesh_lockingGridSize);
    printf("absolute tolerance is %e.\n",h * tolerance);
#endif
    /* sort the elements by center coordinates (lexicographical)*/
    qsort(center,faces->numElements,sizeof(Finley_Mesh_findMatchingFaces_center),Finley_Mesh_findMatchingFaces_compar);
    /* find elements with matching center */
    *numPairs=0;
    /* TODO: OMP */
    for (e=0;e<faces->numElements-1 && Finley_noError();e++) {
        dist=0;
        for (i=0;i<numDim;i++) dist=MAX(dist,ABS(center[e].x[i]-center[e+1].x[i]));
        if (dist < h * tolerance) {
            e_0=center[e].refId;
            e_1=center[e+1].refId;
            elem0[*numPairs]=e_0;
            elem1[*numPairs]=e_1;
            /* now the element e_1 is rotated such that the first node in element e_0 and e_1 have the same coordinates */
            perm=a1;
            perm_tmp=a2;
            for (i=0;i<NN;i++) perm[i]=i;
            while (Finley_noError()) {
                /* if node 0 and perm[0] are the same we are ready */
                getDist(dist,e_0,0,e_1,perm[0]);
                if (dist<=h*tolerance) break;
                if (shiftNodes[0]>=0) {
                    /* rotate the nodes */
                    itmp_ptr=perm;
                    perm=perm_tmp;
                    perm_tmp=itmp_ptr;
                    #pragma ivdep
                    for (i=0;i<NN;i++) perm[i]=perm_tmp[shiftNodes[i]];
                }
                /* if the permutation is back at the identity, i.e. perm[0]=0, the faces don't match: */
                if (perm[0]==0) {
                    sprintf(error_msg,"Mesh_findMatchingFaces: couldn't match first node of element %d to touching element %d",e_0,e_1);
                    Finley_setError(VALUE_ERROR,error_msg);
                }
            }
            /* now we check if the second nodes match */
            if (Finley_noError()) {
                if (numNodesOnFace>1) {
                    getDist(dist,e_0,1,e_1,perm[faceNodes[1]]);
                    /* if the second node does not match we reverse the direction of the nodes */
                    if (dist>h*tolerance) {
                        /* rotate the nodes */
                        if (reverseNodes[0]<0) {
                            sprintf(error_msg,"Mesh_findMatchingFaces: couldn't match the second node of element %d to touching element %d",e_0,e_1);
                            Finley_setError(VALUE_ERROR,error_msg);
                        } else {
                            itmp_ptr=perm;
                            perm=perm_tmp;
                            perm_tmp=itmp_ptr;
                            #pragma ivdep
                            for (i=0;i<NN;i++) perm[i]=perm_tmp[reverseNodes[i]];
                            getDist(dist,e_0,1,e_1,perm[faceNodes[1]]);
                            if (dist>h*tolerance) {
                                sprintf(error_msg,"Mesh_findMatchingFaces: couldn't match the second node of element %d to touching element %d",e_0,e_1);
                                Finley_setError(VALUE_ERROR,error_msg);
                            }
                        }
                    }
                }
            }
            /* we check if the rest of the face nodes match: */
            if (Finley_noError()) {
                for (i=2;i<numNodesOnFace;i++) {
                    n=faceNodes[i];
                    getDist(dist,e_0,n,e_1,perm[n]);
                    if (dist>h*tolerance) {
                        sprintf(error_msg,"Mesh_findMatchingFaces: couldn't match the %d-th node of element %d to touching element %d",i,e_0,e_1);
                        Finley_setError(VALUE_ERROR,error_msg);
                        break;
                    }
                }
            }
            // copy over the permuted nodes of e_1 into matching_nodes_in_elem1
            if (Finley_noError()) {
                for (i=0;i<NN;i++)
                    matching_nodes_in_elem1[INDEX2(i,*numPairs,NN)]=faces->Nodes[INDEX2(perm[i],e_1,NN)];
            }
            (*numPairs)++;
        }
    }
#ifdef Finley_TRACE
    printf("number of pairs of matching faces %d\n",*numPairs);
#endif

    /* clean up */
    delete[] X;
    delete[] center;
    delete[] a1;
#undef getDist
}

