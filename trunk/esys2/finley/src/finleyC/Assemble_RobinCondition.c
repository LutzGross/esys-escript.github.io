/* $Id$ */

/**************************************************************/

/*    assembles a system of numEq natural boundary condition into the stiffness matrix S and right hand side F: */

/*     d*u = y */

/*    u has numComp components. */

/*    Shape of the coefficients: */

/*      d = numEqu x numComp  */
/*      y = numEqu */

/*    The coefficients d and y have to be defined on the integartion points on face elements or not present (=NULL). */

/*    S and F have to be initialized before the routine is called. S or F can be NULL. In this case the left or */
/*    the right hand side of the natural boundary conditions is not processed. */


/**************************************************************/

/*   Copyrights by ACcESS Australia, 2003,2004 */
/*   author: gross@access.edu.au */
/*   Version: $Id$ */

/**************************************************************/

#include "escript/Data/DataC.h"
#include "Assemble.h"
#include "NodeFile.h"
#include "ElementFile.h"
#include "Finley.h"
#include "Util.h"
#ifdef _OPENMP
#include <omp.h>
#endif

/**************************************************************/

void Finley_Assemble_RobinCondition(Finley_NodeFile* nodes,Finley_ElementFile* elements,Finley_SystemMatrix* S, escriptDataC* F, escriptDataC* d, escriptDataC* y, Finley_Assemble_handelShapeMissMatch handelShapeMissMatchForEM) {
  double *EM_S=NULL,*EM_F=NULL,*V=NULL,*dVdv=NULL,*dSdV=NULL,*Area=NULL;
  double time0;
  Assemble_Parameters p;
  index_t *index_row=NULL,*index_col=NULL,color;
  dim_t dimensions[ESCRIPT_MAX_DATA_RANK],e,q;

  if (nodes==NULL || elements==NULL) return;
  if (S==NULL && isEmpty(F)) return;

  /* set all parameters in p*/
  Assemble_getAssembleParameters(nodes,elements,S,F,&p);
  if (Finley_ErrorCode!=NO_ERROR) return;

  /*  get a functionspace */
  type_t funcspace=UNKNOWN;
  updateFunctionSpaceType(funcspace,d);
  updateFunctionSpaceType(funcspace,y);
  if (funcspace==UNKNOWN) return; /* all  data are empty */

  /* check if all function spaces are the same */

  if (! functionSpaceTypeEqual(funcspace,d) ) {
        Finley_ErrorCode=TYPE_ERROR; 
        sprintf(Finley_ErrorMsg,"unexpected function space typ for coefficient d");
  }
  if (! functionSpaceTypeEqual(funcspace,y) ) {
        Finley_ErrorCode=TYPE_ERROR; 
        sprintf(Finley_ErrorMsg,"unexpected function space typ for coefficient y");
  }

  /* check if all function spaces are the same */

  if (! numSamplesEqual(d,p.numQuad,elements->numElements) ) {
        Finley_ErrorCode=TYPE_ERROR; 
        sprintf(Finley_ErrorMsg,"sample points of coefficient d don't match (%d,%d)",p.numQuad,elements->numElements);
  }

  if (! numSamplesEqual(y,p.numQuad,elements->numElements) ) {
        Finley_ErrorCode=TYPE_ERROR; 
        sprintf(Finley_ErrorMsg,"sample points of coefficient y don't match (%d,%d)",p.numQuad,elements->numElements);
  }

  
  /*  check coefficient */
  
  if (p.numEqu==1 && p.numComp==1) {
    if (!isEmpty(d)) {
      if (! isDataPointShapeEqual(d,0,dimensions)) {
          Finley_ErrorCode=TYPE_ERROR;
          sprintf(Finley_ErrorMsg,"coefficient d, rank 0 expected.");
      }
    }
    if (!isEmpty(y)) {
      if (! isDataPointShapeEqual(y,0,dimensions)) {
          Finley_ErrorCode=TYPE_ERROR;
          sprintf(Finley_ErrorMsg,"coefficient y, rank 0 expected.");
      }
    }
  } else {
    if (!isEmpty(d)) {
      dimensions[0]=p.numEqu;
      dimensions[1]=p.numComp;
      if (! isDataPointShapeEqual(d,2,dimensions)) {
          Finley_ErrorCode=TYPE_ERROR;
          sprintf(Finley_ErrorMsg,"coefficient D, expected shape (%d,%d)",dimensions[0],dimensions[1]);
      }
    }
    if (!isEmpty(y)) {
      dimensions[0]=p.numEqu;
      if (! isDataPointShapeEqual(y,1,dimensions)) {
          Finley_ErrorCode=TYPE_ERROR;
          sprintf(Finley_ErrorMsg,"coefficient D, expected shape (%d,)",dimensions[0]);
      }
    }
  }
  
  if (Finley_ErrorCode==NO_ERROR) {
     time0=Finley_timer();
     #pragma omp parallel private(index_row,index_col,EM_S,EM_F,V,dVdv,dSdV,Area,color)
     {
        EM_S=EM_F=V=dVdv=dSdV=Area=NULL;
        index_row=index_col=NULL;
        /* allocate work arrays: */
        EM_S=THREAD_MEMALLOC(p.NN_col*p.NN_row*p.numEqu*p.numComp,double);
        EM_F=THREAD_MEMALLOC(p.NN_row*p.numEqu,double);
        V=THREAD_MEMALLOC(p.NN*p.numDim,double);
        dVdv=THREAD_MEMALLOC(p.numDim*p.numElementDim*p.numQuad,double);
        Area=THREAD_MEMALLOC(p.numQuad,double);
        index_row=THREAD_MEMALLOC(p.NN_row,index_t);
        index_col=THREAD_MEMALLOC(p.NN_col,index_t);
        if (! ( Finley_checkPtr(EM_S) || Finley_checkPtr(EM_F) || Finley_checkPtr(index_row) || Finley_checkPtr(index_col) || 
                Finley_checkPtr(V) || Finley_checkPtr(dVdv) || Finley_checkPtr(Area))) {
           /*  open loop over all colors: */
           for (color=elements->minColor;color<=elements->maxColor;color++) {
              /*  open loop over all elements: */
              #pragma omp for private(e,q) schedule(static)
              for(e=0;e<elements->numElements;e++){
                if (elements->Color[e]==color) {
                   for (q=0;q<p.NN_row;q++) index_row[q]=p.label_row[elements->Nodes[INDEX2(p.row_node[q],e,p.NN)]];
                   /* gather V-coordinates of nodes into V: */
                   Finley_Util_Gather_double(p.NN,&(elements->Nodes[INDEX2(0,e,p.NN)]),p.numDim,nodes->Coordinates,V);
                   /* handel the case where NS!=p.NN (contact elements) */
                   Finley_Assemble_handelShapeMissMatch_Mean_in(p.numDim,p.NN,1,V,p.NS,1);
                   /*  calculate dVdv(i,j,q)=V(i,k)*DSDv(k,j,q) */
                   Finley_Util_SmallMatMult(p.numDim,p.numElementDim*p.numQuad,dVdv,p.NS,V,p.referenceElement->dSdv);
                   /*  */
                   Finley_LengthOfNormalVector(p.numQuad,p.numDim,p.numElementDim,dVdv,Area);
                   /*  scale area: */
                   for (q=0;q<p.numQuad;q++) Area[q]=ABS(Area[q]*p.referenceElement->QuadWeights[q]);
                   /*   integration for the stiffness matrix: */
                   /*   in order to optimze the number of operations the case of constants coefficience needs a bit more work */
                   /*   to make use of some symmetry. */
                   if (S!=NULL) {
                       for (q=0;q<p.NN_col*p.NN_row*p.numEqu*p.numComp;q++) EM_S[q]=0;
                       if (p.numEqu==1 && p.numComp==1) {
                           Finley_Assemble_PDEMatrix_Single2(p.NS_row,p.numElementDim,p.numQuad,
                                                                     p.referenceElement_row->S,dSdV,Area,p.NN_row,EM_S,
                                                                     NULL,TRUE,
                                                                     NULL,TRUE,
                                                                     NULL,TRUE,
                                                                     getSampleData(d,e),isExpanded(d));
                       } else {
                           Finley_Assemble_PDEMatrix_System2(p.NS_row,p.numElementDim,p.numQuad,
                                                                     p.numEqu,p.numComp,p.referenceElement_row->S,dSdV,Area,p.NN_row,EM_S,
                                                                     NULL,TRUE,
                                                                     NULL,TRUE,
                                                                     NULL,TRUE,
                                                                     getSampleData(d,e),isExpanded(d));
                       }
                       /* handel the case of p.NS<NN */
                       handelShapeMissMatchForEM(p.numEqu*p.numComp,p.NN_row,p.NN_col,EM_S,p.NS_row,p.NS_col);
                  /*
                       {int g;
	               for(g=0;g<p.NN*p.NN*p.numEqu*p.numComp;g++) printf("%f ",EM_S[g]);
	               printf("\n");
                       }
                     */
                       /* add  */
                       for (q=0;q<p.NN_col;q++) index_col[q]=p.label_col[elements->Nodes[INDEX2(p.col_node[q],e,p.NN)]];
                       Finley_SystemMatrix_add(S,p.NN_row,index_row,p.numEqu,p.NN_col,index_col,p.numComp,EM_S);
                   }
                   if (!isEmpty(F)) {
                     for (q=0;q<p.NN_row*p.numEqu;q++) EM_F[q]=0;
                     if (p.numEqu==1) {
	                Finley_Assemble_RHSMatrix_Single(p.NS_row,p.numElementDim,p.numQuad,
                                                              p.referenceElement_row->S,dSdV,Area,p.NN_row,EM_F,
                                                              NULL,TRUE,
                                                              getSampleData(y,e),isExpanded(y));
                     } else {
	                Finley_Assemble_RHSMatrix_System(p.NS_row,p.numElementDim,p.numQuad,
                                                              p.numEqu,p.referenceElement_row->S,dSdV,Area,p.NN_row,EM_F,
                                                              NULL,TRUE,
                                                              getSampleData(y,e),isExpanded(y));
                     }
                     /* handel the case of NS<NN */
                     handelShapeMissMatchForEM(p.numEqu,p.NN_row,1,EM_F,p.NS_row,1);
                     /*
                       {int g;
	               for(g=0;g<p.NN*p.numEqu;g++) printf("%f ",EM_F[g]);
	               printf("\n");
                       }
                     */
                     /* add  */
                     Finley_Util_AddScatter(p.NN_row,index_row,p.numEqu,EM_F,getSampleData(F,0));
                   }
                 }
              } /*  end of element loop: */
           } /*  end of color loop: */
        }
        THREAD_MEMFREE(EM_S);
        THREAD_MEMFREE(EM_F);
        THREAD_MEMFREE(V);
        THREAD_MEMFREE(dVdv);
        THREAD_MEMFREE(dSdV);
        THREAD_MEMFREE(Area);
        THREAD_MEMFREE(index_row);
        THREAD_MEMFREE(index_col);
     }
     printf("timing: assemblage natural BC: %.4e sec\n",Finley_timer()-time0);
  }
}
/*
 * $Log$
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
