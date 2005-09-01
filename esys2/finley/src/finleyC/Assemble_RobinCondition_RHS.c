/* $Id$ */

/**************************************************************/

/*    assembles a system of numEq natural boundary condition into the stiffness matrix S and right hand side F: */

/*      y */

/*    Shape of the coefficients: */

/*      y = numEqu */

/*    The coefficient y have to be defined on the integartion points on face elements or not present (=NULL). */

/*    F have to be initialized before the routine is called.F can be NULL. */


/**************************************************************/

/*   Copyrights by ACcESS Australia, 2003,2004,2005 */
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

void Finley_Assemble_RobinCondition_RHS(Finley_NodeFile* nodes,Finley_ElementFile* elements,escriptDataC* F, escriptDataC* y, Finley_Assemble_handelShapeMissMatch handelShapeMissMatchForEM) {
  double *EM_F=NULL,*V=NULL,*dVdv=NULL,*dSdV=NULL,*Area=NULL;
  double time0;
  Assemble_Parameters p;
  index_t *index_row=NULL,color;
  dim_t dimensions[ESCRIPT_MAX_DATA_RANK],e,q;
  bool_t assemble;

  if (nodes==NULL || elements==NULL) return;
  if (isEmpty(F) || isEmpty(y) ) return;

  /* set all parameters in p*/
  Assemble_getAssembleParameters(nodes,elements,NULL,F,&p);
  if (Finley_ErrorCode!=NO_ERROR) return;

  /*  get a functionspace */
  type_t funcspace=UNKNOWN;
  updateFunctionSpaceType(funcspace,y);
  if (funcspace==UNKNOWN) return; /* all  data are empty */

  /* check if all function spaces are the same */

  if (! functionSpaceTypeEqual(funcspace,y) ) {
        Finley_ErrorCode=TYPE_ERROR; 
        sprintf(Finley_ErrorMsg,"unexpected function space type for coefficient y");
  }

  /* check if all function spaces are the same */

  if (! numSamplesEqual(y,p.numQuad,elements->numElements) ) {
        Finley_ErrorCode=TYPE_ERROR; 
        sprintf(Finley_ErrorMsg,"sample points of coefficient y don't match (%d,%d)",p.numQuad,elements->numElements);
  }

  
  /*  check coefficient */
  
  if (p.numEqu==1 && p.numComp==1) {
    if (!isEmpty(y)) {
      if (! isDataPointShapeEqual(y,0,dimensions)) {
          Finley_ErrorCode=TYPE_ERROR;
          sprintf(Finley_ErrorMsg,"coefficient y, rank 0 expected.");
      }
    }
  } else {
    if (!isEmpty(y)) {
      dimensions[0]=p.numEqu;
      if (! isDataPointShapeEqual(y,1,dimensions)) {
          Finley_ErrorCode=TYPE_ERROR;
          sprintf(Finley_ErrorMsg,"coefficient y, expected shape (%d,)",dimensions[0]);
      }
    }
  }
  
  if (Finley_ErrorCode==NO_ERROR) {
     time0=Finley_timer();
     #pragma omp parallel private(assemble,index_row,EM_F,V,dVdv,dSdV,Area,color)
     {
        EM_F=V=dVdv=dSdV=Area=NULL;
        index_row=NULL;
        /* allocate work arrays: */
        EM_F=THREAD_MEMALLOC(p.NN_row*p.numEqu,double);
        V=THREAD_MEMALLOC(p.NN*p.numDim,double);
        dVdv=THREAD_MEMALLOC(p.numDim*p.numElementDim*p.numQuad,double);
        Area=THREAD_MEMALLOC(p.numQuad,double);
        index_row=THREAD_MEMALLOC(p.NN_row,index_t);
        if (! ( Finley_checkPtr(EM_F) || Finley_checkPtr(index_row) || Finley_checkPtr(V) || Finley_checkPtr(dVdv) || Finley_checkPtr(Area))) {
           /*  open loop over all colors: */
           for (color=elements->minColor;color<=elements->maxColor;color++) {
              /*  open loop over all elements: */
              #pragma omp for private(e,q) schedule(dynamic)
              for(e=0;e<elements->numElements;e++){
                if (elements->Color[e]==color) { 
                  if (isExpanded(y)) {
                      assemble=Finley_Util_anyNonZeroDouble(p.numEqu*p.numQuad,getSampleData(y,e));
                  } else {
                     assemble=Finley_Util_anyNonZeroDouble(p.numEqu,getSampleData(y,e));
                  }
                  if (assemble) {
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
                     /* add  */
                     Finley_Util_AddScatter(p.NN_row,index_row,p.numEqu,EM_F,getSampleData(F,0));
                   }
                 }
              } /*  end of element loop: */
           } /*  end of color loop: */
        }
        THREAD_MEMFREE(EM_F);
        THREAD_MEMFREE(V);
        THREAD_MEMFREE(dVdv);
        THREAD_MEMFREE(dSdV);
        THREAD_MEMFREE(Area);
        THREAD_MEMFREE(index_row);
     }
     #ifdef Finley_TRACE
     printf("timing: assemblage natural BC right hand side: %.4e sec\n",Finley_timer()-time0);
     #endif
  }
}
/*
 * $Log$
 * Revision 1.3  2005/09/01 03:31:35  jgs
 * Merge of development branch dev-02 back to main trunk on 2005-09-01
 *
 * Revision 1.2  2005/08/12 01:45:43  jgs
 * erge of development branch dev-02 back to main trunk on 2005-08-12
 *
 * Revision 1.1.2.2  2005/08/24 02:02:18  gross
 * timing output switched off. solver output can be swiched through getSolution(verbose=True) now.
 *
 * Revision 1.1.2.1  2005/08/04 22:41:11  gross
 * some extra routines for finley that might speed-up RHS assembling in some cases (not actived right now)
 *
 *
 */
