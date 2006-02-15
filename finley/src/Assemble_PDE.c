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

/*    assembles the system of numEq PDEs into the stiffness matrix S and right hand side F */

/*     -div(A*grad u)-div(B*u)+C*grad u + D*u= -div X + Y */

/*      -(A_{k,i,m,j} u_m,j)_i-(B_{k,i,m} u_m)_i+C_{k,m,j} u_m,j-D_{k,m} u_m = -(X_{k,i})_i + Y_k */

/*    u has numComp components. */

/*    Shape of the coefficients: */

/*      A = numEqu x numDim x numComp x numDim */
/*      B = numDim x numEqu x numComp  */
/*      C = numEqu x numDim x numComp  */
/*      D = numEqu x numComp  */
/*      X = numEqu x numDim   */
/*      Y = numEqu */

/*    The coefficients A,B,C,D,X and Y have to be defined on the integartion points or not present (=NULL). */

/*    S and F have to be initialized before the routine is called. S or F can be NULL. In this case the left or */
/*    the right hand side of the PDE is not processed.  */

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

void Finley_Assemble_PDE(Finley_NodeFile* nodes,Finley_ElementFile* elements,Paso_SystemMatrix* S, escriptDataC* F,
			 escriptDataC* A, escriptDataC* B, escriptDataC* C, escriptDataC* D, escriptDataC* X, escriptDataC* Y ) {

  char error_msg[LenErrorMsg_MAX];
  double *EM_S=NULL,*EM_F=NULL,*V=NULL,*dVdv=NULL,*dSdV=NULL,*Vol=NULL,*dvdV=NULL;
  double time0;
  dim_t dimensions[ESCRIPT_MAX_DATA_RANK],e,q;
  Assemble_Parameters p;
  index_t *index_row=NULL,*index_col=NULL,color;
  Finley_resetError();

  if (nodes==NULL || elements==NULL) return;
  if (S==NULL && isEmpty(F)) return;

  /* set all parameters in p*/
  Assemble_getAssembleParameters(nodes,elements,S,F,&p);
  if (! Finley_noError()) return;

  /*  this function assumes NS=NN */
  if (p.NN!=p.NS) {
    Finley_setError(SYSTEM_ERROR,"Finley_Assemble_PDE: for Finley_Assemble_PDE numNodes and numShapes have to be identical.");
    return;
  } 
  if (p.numDim!=p.numElementDim) {
    Finley_setError(SYSTEM_ERROR,"Finley_Assemble_PDE: Finley_Assemble_PDE accepts volume elements only.");
    return;
  }
  /*  get a functionspace */
  type_t funcspace=UNKNOWN;
  updateFunctionSpaceType(funcspace,A);
  updateFunctionSpaceType(funcspace,B);
  updateFunctionSpaceType(funcspace,C);
  updateFunctionSpaceType(funcspace,D);
  updateFunctionSpaceType(funcspace,X);
  updateFunctionSpaceType(funcspace,Y);
  if (funcspace==UNKNOWN) return; /* all  data are empty */

  /* check if all function spaces are the same */

  if (! functionSpaceTypeEqual(funcspace,A) ) {
        Finley_setError(TYPE_ERROR,"Finley_Assemble_PDE: unexpected function space type for coefficient A");
  }
  if (! functionSpaceTypeEqual(funcspace,B) ) {
        Finley_setError(TYPE_ERROR,"Finley_Assemble_PDE: unexpected function space type for coefficient B");
  }
  if (! functionSpaceTypeEqual(funcspace,C) ) {
        Finley_setError(TYPE_ERROR,"Finley_Assemble_PDE: unexpected function space type for coefficient C");
  }
  if (! functionSpaceTypeEqual(funcspace,D) ) {
        Finley_setError(TYPE_ERROR,"Finley_Assemble_PDE: unexpected function space type for coefficient D");
  }
  if (! functionSpaceTypeEqual(funcspace,X) ) {
        Finley_setError(TYPE_ERROR,"Finley_Assemble_PDE: unexpected function space type for coefficient X");
  }
  if (! functionSpaceTypeEqual(funcspace,Y) ) {
        Finley_setError(TYPE_ERROR,"Finley_Assemble_PDE: unexpected function space type for coefficient Y");
  }

  /* check if all function spaces are the same */

  if (! numSamplesEqual(A,p.numQuad,elements->numElements) ) {
        sprintf(error_msg,"Finley_Assemble_PDE: sample points of coefficient A don't match (%d,%d)",p.numQuad,elements->numElements);
        Finley_setError(TYPE_ERROR,error_msg);
  }

  if (! numSamplesEqual(B,p.numQuad,elements->numElements) ) {
        sprintf(error_msg,"Finley_Assemble_PDE: sample points of coefficient B don't match (%d,%d)",p.numQuad,elements->numElements);
        Finley_setError(TYPE_ERROR,error_msg);
  }

  if (! numSamplesEqual(C,p.numQuad,elements->numElements) ) {
        sprintf(error_msg,"Finley_Assemble_PDE: sample points of coefficient C don't match (%d,%d)",p.numQuad,elements->numElements);
        Finley_setError(TYPE_ERROR,error_msg);
  }

  if (! numSamplesEqual(D,p.numQuad,elements->numElements) ) {
        sprintf(error_msg,"Finley_Assemble_PDE: sample points of coefficient D don't match (%d,%d)",p.numQuad,elements->numElements);
        Finley_setError(TYPE_ERROR,error_msg);
  }

  if (! numSamplesEqual(X,p.numQuad,elements->numElements) ) {
        sprintf(error_msg,"Finley_Assemble_PDE: sample points of coefficient X don't match (%d,%d)",p.numQuad,elements->numElements);
        Finley_setError(TYPE_ERROR,error_msg);
  }

  if (! numSamplesEqual(Y,p.numQuad,elements->numElements) ) {
        sprintf(error_msg,"Finley_Assemble_PDE: sample points of coefficient Y don't match (%d,%d)",p.numQuad,elements->numElements);
        Finley_setError(TYPE_ERROR,error_msg);
  }

  /*  check the dimensions: */
  
  if (p.numEqu==1 && p.numComp==1) {
    if (!isEmpty(A)) {
      dimensions[0]=p.numDim;
      dimensions[1]=p.numDim;
      if (!isDataPointShapeEqual(A,2,dimensions)) {
          sprintf(error_msg,"Finley_Assemble_PDE: coefficient A: illegal shape, expected shape (%d,%d)",dimensions[0],dimensions[1]);
          Finley_setError(TYPE_ERROR,error_msg);
      }
    }
    if (!isEmpty(B)) {
       dimensions[0]=p.numDim;
       if (!isDataPointShapeEqual(B,1,dimensions)) {
          sprintf(error_msg,"Finley_Assemble_PDE: coefficient B: illegal shape (%d,)",dimensions[0]);
          Finley_setError(TYPE_ERROR,error_msg);
       }
    }
    if (!isEmpty(C)) {
       dimensions[0]=p.numDim;
       if (!isDataPointShapeEqual(C,1,dimensions)) {
          sprintf(error_msg,"Finley_Assemble_PDE: coefficient C, expected shape (%d,)",dimensions[0]);
          Finley_setError(TYPE_ERROR,error_msg);
       }
    }
    if (!isEmpty(D)) {
       if (!isDataPointShapeEqual(D,0,dimensions)) {
          Finley_setError(TYPE_ERROR,"Finley_Assemble_PDE: coefficient D, rank 0 expected.");
       }
    }
    if (!isEmpty(X)) {
       dimensions[0]=p.numDim;
       if (!isDataPointShapeEqual(X,1,dimensions)) {
          sprintf(error_msg,"Finley_Assemble_PDE: coefficient X, expected shape (%d,",dimensions[0]);
          Finley_setError(TYPE_ERROR,error_msg);
       }
    }
    if (!isEmpty(Y)) {
       if (!isDataPointShapeEqual(Y,0,dimensions)) {
          Finley_setError(TYPE_ERROR,"Finley_Assemble_PDE: coefficient Y, rank 0 expected.");
       }
    } 
  } else {
    if (!isEmpty(A)) {
      dimensions[0]=p.numEqu;
      dimensions[1]=p.numDim;
      dimensions[2]=p.numComp;
      dimensions[3]=p.numDim;
      if (!isDataPointShapeEqual(A,4,dimensions)) {
          sprintf(error_msg,"Finley_Assemble_PDE: coefficient A, expected shape (%d,%d,%d,%d)",dimensions[0],dimensions[1],dimensions[2],dimensions[3]);
          Finley_setError(TYPE_ERROR,error_msg);
      }
    }
    if (!isEmpty(B)) {
      dimensions[0]=p.numEqu;
      dimensions[1]=p.numDim;
      dimensions[2]=p.numComp;
      if (!isDataPointShapeEqual(B,3,dimensions)) {
          sprintf(error_msg,"Finley_Assemble_PDE: coefficient B, expected shape (%d,%d,%d)",dimensions[0],dimensions[1],dimensions[2]);
          Finley_setError(TYPE_ERROR,error_msg);
      }
    }
    if (!isEmpty(C)) {
      dimensions[0]=p.numEqu;
      dimensions[1]=p.numComp;
      dimensions[2]=p.numDim;
      if (!isDataPointShapeEqual(C,3,dimensions)) {
          sprintf(error_msg,"Finley_Assemble_PDE: coefficient C, expected shape (%d,%d,%d)",dimensions[0],dimensions[1],dimensions[2]);
          Finley_setError(TYPE_ERROR,error_msg);
      }
    }
    if (!isEmpty(D)) {
      dimensions[0]=p.numEqu;
      dimensions[1]=p.numComp;
      if (!isDataPointShapeEqual(D,2,dimensions)) {
          sprintf(error_msg,"Finley_Assemble_PDE: coefficient D, expected shape (%d,%d)",dimensions[0],dimensions[1]);
          Finley_setError(TYPE_ERROR,error_msg);
      }
    }
    if (!isEmpty(X)) {
      dimensions[0]=p.numEqu;
      dimensions[1]=p.numDim;
      if (!isDataPointShapeEqual(X,2,dimensions)) {
          sprintf(error_msg,"Finley_Assemble_PDE: coefficient X, expected shape (%d,%d)",dimensions[0],dimensions[1]);
          Finley_setError(TYPE_ERROR,error_msg);
      }
    }
    if (!isEmpty(Y)) {
      dimensions[0]=p.numEqu;
      if (!isDataPointShapeEqual(Y,1,dimensions)) {
          sprintf(error_msg,"Finley_Assemble_PDE: coefficient Y, expected shape (%d,)",dimensions[0]);
          Finley_setError(TYPE_ERROR,error_msg);
      }
    }
  }

  if (Finley_noError()) {
     time0=Finley_timer();
     #pragma omp parallel private(index_col,index_row,EM_S,EM_F,V,dVdv,dSdV,Vol,dvdV,color,q) \
            firstprivate(elements,nodes,S,F,A,B,C,D,X,Y)
     {
         EM_S=EM_F=V=dVdv=dSdV=Vol=dvdV=NULL;
         index_row=index_col=NULL;

         /* allocate work arrays: */

         EM_S=(double*) THREAD_MEMALLOC(p.NN_row*p.NN_col*p.numEqu*p.numComp,double);
         EM_F=(double*) THREAD_MEMALLOC(p.NN_row*p.numEqu,double);
         V=(double*) THREAD_MEMALLOC(p.NN*p.numDim,double);
         dVdv=(double*) THREAD_MEMALLOC(p.numDim*p.numDim*p.numQuad,double);
         dvdV=(double*) THREAD_MEMALLOC(p.numDim*p.numDim*p.numQuad,double);
         dSdV=(double*) THREAD_MEMALLOC(p.NS_row*p.numQuad*p.numDim,double);
         Vol=(double*) THREAD_MEMALLOC(p.numQuad,double);
         index_col=(index_t*) THREAD_MEMALLOC(p.NN_col,index_t);
         index_row=(index_t*) THREAD_MEMALLOC(p.NN_row,index_t);

         if (! (Finley_checkPtr(EM_S) || Finley_checkPtr(EM_F) || Finley_checkPtr(V) || Finley_checkPtr(index_col) ||
                Finley_checkPtr(index_row) || Finley_checkPtr(dVdv) || Finley_checkPtr(dSdV) || Finley_checkPtr(Vol) ))  {

           /*  open loop over all colors: */
           for (color=elements->minColor;color<=elements->maxColor;color++) {
              /*  open loop over all elements: */
              #pragma omp for private(e) schedule(static) 
              for(e=0;e<elements->numElements;e++){
                if (elements->Color[e]==color) {
//============================
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
//============================
    
                   /*   integration for the stiffness matrix: */
                   /*   in order to optimze the number of operations the case of constants coefficience needs a bit more work */
                   /*   to make use of some symmetry. */

                     if (S!=NULL) {
                       for (q=0;q<p.NN_row*p.NN_col*p.numEqu*p.numComp;q++) EM_S[q]=0;
                       if (p.numEqu==1 && p.numComp==1) {
  	                   Finley_Assemble_PDEMatrix_Single2(p.NS_row,p.numDim,p.numQuad,
                                                                     p.referenceElement_row->S,dSdV,Vol,p.NN_row,EM_S,
                                                                     getSampleData(A,e),isExpanded(A),
                                                                     getSampleData(B,e),isExpanded(B),
                                                                     getSampleData(C,e),isExpanded(C),
                                                                     getSampleData(D,e),isExpanded(D));
                       } else {
  	                   Finley_Assemble_PDEMatrix_System2(p.NS_row,p.numDim,p.numQuad,p.numEqu,p.numComp,
                                                                     p.referenceElement_row->S,dSdV,Vol,p.NN_row,EM_S,
                                                                     getSampleData(A,e),isExpanded(A),
                                                                     getSampleData(B,e),isExpanded(B),
                                                                     getSampleData(C,e),isExpanded(C),
                                                                     getSampleData(D,e),isExpanded(D));
                       }
                       for (q=0;q<p.NN_col;q++) index_col[q]=p.label_col[elements->Nodes[INDEX2(p.col_node[q],e,p.NN)]];
                       Finley_Assemble_addToSystemMatrix(S,p.NN_row,index_row,p.numEqu,p.NN_col,index_col,p.numComp,EM_S);
                     }
                     if (!isEmpty(F)) {
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
                       Finley_Util_AddScatter(p.NN_row,index_row,p.numEqu,EM_F,getSampleData(F,0));
                    }
                }
              }
            }
         }
         /* clean up */
         THREAD_MEMFREE(EM_S);
         THREAD_MEMFREE(EM_F);
         THREAD_MEMFREE(V);
         THREAD_MEMFREE(dVdv);
         THREAD_MEMFREE(dvdV);
         THREAD_MEMFREE(dSdV);
         THREAD_MEMFREE(Vol); 
         THREAD_MEMFREE(index_col); 
         THREAD_MEMFREE(index_row); 
     }
     #ifdef Finley_TRACE
     printf("timing: assemblage PDE: %.4e sec\n",Finley_timer()-time0);
     #endif
  }
}
/*
 * $Log$
 * Revision 1.8  2005/09/15 03:44:21  jgs
 * Merge of development branch dev-02 back to main trunk on 2005-09-15
 *
 * Revision 1.7  2005/09/01 03:31:35  jgs
 * Merge of development branch dev-02 back to main trunk on 2005-09-01
 *
 * Revision 1.6  2005/08/12 01:45:42  jgs
 * erge of development branch dev-02 back to main trunk on 2005-08-12
 *
 * Revision 1.5.2.3  2005/09/07 06:26:17  gross
 * the solver from finley are put into the standalone package paso now
 *
 * Revision 1.5.2.2  2005/08/24 02:02:18  gross
 * timing output switched off. solver output can be swiched through getSolution(verbose=True) now.
 *
 * Revision 1.5.2.1  2005/08/03 08:54:27  gross
 * contact element assemblage was called with wrong element table pointer
 *
 * Revision 1.5  2005/07/08 04:07:46  jgs
 * Merge of development branch back to main trunk on 2005-07-08
 *
 * Revision 1.4  2004/12/15 07:08:32  jgs
 * *** empty log message ***
 * Revision 1.1.1.1.2.2  2005/06/29 02:34:47  gross
 * some changes towards 64 integers in finley
 *
 * Revision 1.1.1.1.2.1  2004/11/24 01:37:12  gross
 * some changes dealing with the integer overflow in memory allocation. Finley solves 4M unknowns now
 *
 *
 *
 */
