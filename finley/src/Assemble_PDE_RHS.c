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

/*    assembles the system of numEq PDEs into side F */

/*     -div X + Y */

/*    Shape of the coefficients: */

/*      X = numEqu x numDim   */
/*      Y = numEqu */

/*    The coefficients X and Y have to be defined on the integartion points or not present (=NULL). */

/*    F has to be initialized before the routine is called. F can be NULL.  */

/*    The routine does not consider any boundary conditions. */

/**************************************************************/

/*  Author: gross@access.edu.au */
/*  Version: $Id$ */

/**************************************************************/

#include "Assemble.h"
#include "Util.h"
#ifdef _OPENMP
#include <omp.h>
#endif


/**************************************************************/

void Finley_Assemble_PDE_RHS(Finley_NodeFile* nodes,Finley_ElementFile* elements,escriptDataC* F,escriptDataC* X, escriptDataC* Y ) {

  char error_msg[LenErrorMsg_MAX];
  double *EM_F=NULL,*V=NULL,*dVdv=NULL,*dSdV=NULL,*Vol=NULL,*dvdV=NULL;
  double time0;
  dim_t dimensions[ESCRIPT_MAX_DATA_RANK],e,q;
  Assemble_Parameters p;
  index_t *index_row=NULL,color;
  bool_t assemble;
  Finley_resetError();

  if (nodes==NULL || elements==NULL) return;
  if (isEmpty(F) || (isEmpty(Y) && isEmpty(X)) ) return;

  /* set all parameters in p*/
  Assemble_getAssembleParameters(nodes,elements,NULL,F,&p);
  if (! Finley_noError()) return;

  /*  this function assumes NS=NN */
  if (p.NN!=p.NS) {
    Finley_setError(SYSTEM_ERROR,"__FILE__: for Finley_Assemble_PDE numNodes and numShapes have to be identical.");
    return;
  } 
  if (p.numDim!=p.numElementDim) {
    Finley_setError(SYSTEM_ERROR,"__FILE__: Finley_Assemble_PDE accepts volume elements only.");
    return;
  }
  /*  get a functionspace */
  type_t funcspace=UNKNOWN;
  updateFunctionSpaceType(funcspace,X);
  updateFunctionSpaceType(funcspace,Y);
  if (funcspace==UNKNOWN) return; /* all  data are empty */

  /* check if all function spaces are the same */

  if (! functionSpaceTypeEqual(funcspace,X) ) {
        Finley_setError(TYPE_ERROR,"__FILE__: unexpected function space type for coefficient X");
  }
  if (! functionSpaceTypeEqual(funcspace,Y) ) {
        Finley_setError(TYPE_ERROR,"__FILE__: unexpected function space type for coefficient Y");
  }

  /* check if all function spaces are the same */

  if (! numSamplesEqual(X,p.numQuad,elements->numElements) ) {
        sprintf(error_msg,"__FILE__: sample points of coefficient X don't match (%d,%d)",p.numQuad,elements->numElements);
        Finley_setError(TYPE_ERROR,error_msg);
  }

  if (! numSamplesEqual(Y,p.numQuad,elements->numElements) ) {
        sprintf(error_msg,"__FILE__: sample points of coefficient Y don't match (%d,%d)",p.numQuad,elements->numElements);
        Finley_setError(TYPE_ERROR,error_msg);
  }

  /*  check the dimensions: */
  
  if (p.numEqu==1 && p.numComp==1) {
    if (!isEmpty(X)) {
       dimensions[0]=p.numDim;
       if (!isDataPointShapeEqual(X,1,dimensions)) {
          sprintf(error_msg,"__FILE__: coefficient X, expected shape (%d,",dimensions[0]);
          Finley_setError(TYPE_ERROR,error_msg);
       }
    }
    if (!isEmpty(Y)) {
       if (!isDataPointShapeEqual(Y,0,dimensions)) {
          Finley_setError(TYPE_ERROR,"__FILE__: coefficient Y, rank 0 expected.");
       }
    } 
  } else {
    if (!isEmpty(X)) {
      dimensions[0]=p.numEqu;
      dimensions[1]=p.numDim;
      if (!isDataPointShapeEqual(X,2,dimensions)) {
          sprintf(error_msg,"__FILE__: coefficient X, expected shape (%d,%d)",dimensions[0],dimensions[1]);
          Finley_setError(TYPE_ERROR,error_msg);
      }
    }
    if (!isEmpty(Y)) {
      dimensions[0]=p.numEqu;
      if (!isDataPointShapeEqual(Y,1,dimensions)) {
          sprintf(error_msg,"__FILE__: coefficient Y, expected shape (%d,)",dimensions[0]);
          Finley_setError(TYPE_ERROR,error_msg);
      }
    }
  }

  if (Finley_noError()) {
     time0=Finley_timer();
     #pragma omp parallel private(index_row,EM_F,V,dVdv,dSdV,Vol,dvdV,color,q,assemble) 
     {
         EM_F=V=dVdv=dSdV=Vol=dvdV=NULL;
         index_row=NULL;

         /* allocate work arrays: */

         EM_F=(double*) THREAD_MEMALLOC(p.NN_row*p.numEqu,double);
         V=(double*) THREAD_MEMALLOC(p.NN*p.numDim,double);
         dVdv=(double*) THREAD_MEMALLOC(p.numDim*p.numDim*p.numQuad,double);
         dvdV=(double*) THREAD_MEMALLOC(p.numDim*p.numDim*p.numQuad,double);
         dSdV=(double*) THREAD_MEMALLOC(p.NS_row*p.numQuad*p.numDim,double);
         Vol=(double*) THREAD_MEMALLOC(p.numQuad,double);
         index_row=(index_t*) THREAD_MEMALLOC(p.NN_row,index_t);

         if (! (Finley_checkPtr(EM_F) || Finley_checkPtr(V) || 
                Finley_checkPtr(index_row) || Finley_checkPtr(dVdv) || Finley_checkPtr(dSdV) || Finley_checkPtr(Vol) ))  {

           /*  open loop over all colors: */
           for (color=elements->minColor;color<=elements->maxColor;color++) {
              /*  open loop over all elements: */
              #pragma omp for private(e) schedule(dynamic) 
              for(e=0;e<elements->numElements;e++){
                if (elements->Color[e]==color) {
                  /* check if element has to be processed */
                  if (isEmpty(X)) {
                      assemble=FALSE;
                  } else {
                      if (isExpanded(X)) {
                         assemble=Finley_Util_anyNonZeroDouble(p.numEqu*p.numDim*p.numQuad,getSampleData(X,e));
                     } else {
                         assemble=Finley_Util_anyNonZeroDouble(p.numEqu*p.numDim,getSampleData(X,e));
                     } 
                  }
                  if (!assemble && !isEmpty(Y)) {
                      if (isExpanded(Y)) {
                         assemble=Finley_Util_anyNonZeroDouble(p.numEqu*p.numQuad,getSampleData(Y,e));
                      } else {
                         assemble=Finley_Util_anyNonZeroDouble(p.numEqu,getSampleData(Y,e));
                      }
                  }
                  if (assemble) {
                       for (q=0;q<p.NN_row;q++) index_row[q]=p.label_row[elements->Nodes[INDEX2(p.row_node[q],e,p.NN)]];
                       /* gather V-coordinates of nodes into V: */
		       Finley_Util_Gather_double(p.NN,&(elements->Nodes[INDEX2(0,e,p.NN)]),p.numDim,nodes->Coordinates,V);
                       /*  calculate dVdv(i,j,q)=V(i,k)*DSDv(k,j,q) */
		       Finley_Util_SmallMatMult(p.numDim,p.numDim*p.numQuad,dVdv,p.NS,V,p.referenceElement->dSdv);
                       /*  dvdV=invert(dVdv) inplace */
		       Finley_Util_InvertSmallMat(p.numQuad,p.numDim,dVdv,dvdV,Vol);
                       /*  calculate dSdV=DSDv*DvDV */
		       Finley_Util_SmallMatSetMult(p.numQuad,p.NS_row,p.numDim,dSdV,p.numDim,p.referenceElement_row->dSdv,dvdV);
                       /*  scale volume: */
		       for (q=0;q<p.numQuad;q++) Vol[q]=ABS(Vol[q]*p.referenceElement->QuadWeights[q]);
                       for (q=0;q<p.NN_row*p.numEqu;q++) EM_F[q]=0;
                       if (p.numEqu==1) {
  	                        Finley_Assemble_RHSMatrix_Single(p.NS_row,p.numDim,p.numQuad,
                                                                      p.referenceElement_row->S,dSdV,Vol,p.NN_row,EM_F,
                                                                      getSampleData(X,e),isExpanded(X),
                                                                      getSampleData(Y,e),isExpanded(Y));
                            } else {
  	                        Finley_Assemble_RHSMatrix_System(p.NS_row,p.numDim,p.numQuad,
                                                                      p.numEqu,p.referenceElement_row->S,dSdV,Vol,p.NN_row,EM_F,
                                                                      getSampleData(X,e),isExpanded(X),
                                                                      getSampleData(Y,e),isExpanded(Y));
                       }
                       /* add  */
#ifndef PASO_MPI
                       Finley_Util_AddScatter(p.NN_row,index_row,p.numEqu,EM_F,getSampleData(F,0));
#else
                       /* add influence to local degrees of freedom */
                       /* it might be more efficient to use a local buffer for copying into using Finley_Util_AddScatter()
                          then doing a memcpy at the end */
                       Finley_Util_AddScatter_upperBound(p.NN_row,index_row,p.numEqu,EM_F,getSampleData(F,0),p.degreeOfFreedomUpperBound);
#endif
                    }
                  }
              }
            }
         }
         /* clean up */
         THREAD_MEMFREE(EM_F);
         THREAD_MEMFREE(V);
         THREAD_MEMFREE(dVdv);
         THREAD_MEMFREE(dvdV);
         THREAD_MEMFREE(dSdV);
         THREAD_MEMFREE(Vol); 
         THREAD_MEMFREE(index_row); 
     }
     #ifdef Finley_TRACE
     printf("timing: assemblage PDE Right hand Side: %.4e sec\n",Finley_timer()-time0);
     #endif
  }
}
/*
 * $Log$
 * Revision 1.4  2005/09/15 03:44:21  jgs
 * Merge of development branch dev-02 back to main trunk on 2005-09-15
 *
 * Revision 1.3  2005/09/01 03:31:35  jgs
 * Merge of development branch dev-02 back to main trunk on 2005-09-01
 *
 * Revision 1.2  2005/08/12 01:45:42  jgs
 * erge of development branch dev-02 back to main trunk on 2005-08-12
 *
 * Revision 1.1.2.3  2005/09/07 06:26:17  gross
 * the solver from finley are put into the standalone package paso now
 *
 * Revision 1.1.2.2  2005/08/24 02:02:18  gross
 * timing output switched off. solver output can be swiched through getSolution(verbose=True) now.
 *
 * Revision 1.1.2.1  2005/08/04 22:41:11  gross
 * some extra routines for finley that might speed-up RHS assembling in some cases (not actived right now)
 *
 *
 */
