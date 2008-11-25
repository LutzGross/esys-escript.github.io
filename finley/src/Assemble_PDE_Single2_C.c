
/*******************************************************
*
* Copyright (c) 2003-2008 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


/**************************************************************/

/*    assembles a single PDEs into the stiffness matrix S right hand side F  */

/*      D u_m  and Y */

/*    in a 3D domain. The shape functions for test and solution must be identical  */
/*    and 2* row_NS == row_NN (contact elements)                                   */

/*    Shape of the coefficients: */

/*      D = scalar    */
/*      Y = scalar    */


/**************************************************************/


#include "Assemble.h"
#include "Util.h"
#ifdef _OPENMP
#include <omp.h>
#endif


/**************************************************************/

void  Finley_Assemble_PDE_Single2_C(Assemble_Parameters p, Finley_ElementFile* elements,
                                    Paso_SystemMatrix* Mat, escriptDataC* F, escriptDataC* D, escriptDataC* Y) {

    index_t color;
    dim_t e;
    __const double *D_p, *Y_p;
    double *EM_S, *EM_F, *Vol;
    index_t *row_index;
    register dim_t q, s,r;
    register double rtmp, rtmp_D;
    bool_t add_EM_F, add_EM_S;

    bool_t extendedD=isExpanded(D);
    bool_t extendedY=isExpanded(Y);
    double *F_p=getSampleDataRW(F,0);
    double *S=p.row_jac->ReferenceElement->S;


    #pragma omp parallel private(color,EM_S, EM_F, Vol, D_p, Y_p,row_index,q, s,r,rtmp, rtmp_D,add_EM_F, add_EM_S)
    {
       EM_S=THREAD_MEMALLOC(p.row_NN*p.col_NN,double);
       EM_F=THREAD_MEMALLOC(p.row_NN,double);
       row_index=THREAD_MEMALLOC(p.row_NN,index_t);
                                                                                                                                                                                                     
       if (!Finley_checkPtr(EM_S) && !Finley_checkPtr(EM_F) && !Finley_checkPtr(row_index) ) {

          for (color=elements->minColor;color<=elements->maxColor;color++) {
             /*  open loop over all elements: */
             #pragma omp for private(e) schedule(static)
             for(e=0;e<elements->numElements;e++){
                if (elements->Color[e]==color) {
                   Vol=&(p.row_jac->volume[INDEX2(0,e,p.numQuad)]);
                   add_EM_F=FALSE;
                   add_EM_S=FALSE;
                   /************************************************************* */
                   /* process D */
                   /**************************************************************/
                   D_p=getSampleDataRO(D,e);
                   if (NULL!=D_p) {
                     add_EM_S=TRUE;
                     if (extendedD) {
                         for (s=0;s<p.row_NS;s++) {
                           for (r=0;r<p.col_NS;r++) {
                              rtmp=0;
                              for (q=0;q<p.numQuad;q++) rtmp+=Vol[q]*S[INDEX2(s,q,p.row_NS)]*D_p[q]*S[INDEX2(r,q,p.row_NS)];
                              EM_S[INDEX4(0,0,s         ,r         ,p.numEqu,p.numComp,p.row_NN)]= rtmp;
                              EM_S[INDEX4(0,0,s         ,r+p.col_NS,p.numEqu,p.numComp,p.row_NN)]=-rtmp;
                              EM_S[INDEX4(0,0,s+p.row_NS,r         ,p.numEqu,p.numComp,p.row_NN)]=-rtmp;
                              EM_S[INDEX4(0,0,s+p.row_NS,r+p.col_NS,p.numEqu,p.numComp,p.row_NN)]= rtmp;
                           }
                         }
                     } else {
                         for (s=0;s<p.row_NS;s++) {
                           for (r=0;r<p.col_NS;r++) {
                               rtmp=0;
                               for (q=0;q<p.numQuad;q++) rtmp+=Vol[q]*S[INDEX2(s,q,p.row_NS)]*S[INDEX2(r,q,p.row_NS)];
                               rtmp_D=rtmp*D_p[0];
                               EM_S[INDEX4(0,0,s         ,r         ,p.numEqu,p.numComp,p.row_NN)]= rtmp_D;
                               EM_S[INDEX4(0,0,s         ,r+p.col_NS,p.numEqu,p.numComp,p.row_NN)]=-rtmp_D;
                               EM_S[INDEX4(0,0,s+p.row_NS,r         ,p.numEqu,p.numComp,p.row_NN)]=-rtmp_D;
                               EM_S[INDEX4(0,0,s+p.row_NS,r+p.col_NS,p.numEqu,p.numComp,p.row_NN)]= rtmp_D;
                           }
                         }
                     }
                   }
                  /**************************************************************/
                  /*   process Y: */
                  /**************************************************************/
                   Y_p=getSampleDataRO(Y,e);
                   if (NULL!=Y_p) {
                     add_EM_F=TRUE;
                     if (extendedY) {
                        for (s=0;s<p.row_NS;s++) {
                          rtmp=0;
                          for (q=0;q<p.numQuad;q++) rtmp+=Vol[q]*S[INDEX2(s,q,p.row_NS)]*Y_p[q];
                          EM_F[INDEX2(0,s         ,p.numEqu)]=-rtmp;
                          EM_F[INDEX2(0,s+p.row_NS,p.numEqu)]= rtmp;
                        }
                      } else {
                        for (s=0;s<p.row_NS;s++) {
                            rtmp=0;
                            for (q=0;q<p.numQuad;q++) rtmp+=Vol[q]*S[INDEX2(s,q,p.row_NS)];
                            rtmp_D=rtmp*Y_p[0];
                            EM_F[INDEX2(0,s         ,p.numEqu)]=-rtmp_D;
                            EM_F[INDEX2(0,s+p.row_NS,p.numEqu)]= rtmp_D;
                        }
                      }
                    }
                    /***********************************************************************************************/
                    /* add the element matrices onto the matrix and right hand side                                */
                    /***********************************************************************************************/
                    for (q=0;q<p.row_NN;q++) row_index[q]=p.row_DOF[elements->Nodes[INDEX2(p.row_node[q],e,p.row_NN)]];
   
                    if (add_EM_F) Finley_Util_AddScatter(p.row_NN,row_index,p.numEqu,EM_F,F_p, p.row_DOF_UpperBound);
                    if (add_EM_S) Finley_Assemble_addToSystemMatrix(Mat,p.row_NN,row_index,p.numEqu,p.col_NN,row_index,p.numComp,EM_S);
   
                } /* end color check */
             } /* end element loop */
         } /* end color loop */
           
         THREAD_MEMFREE(EM_S);
         THREAD_MEMFREE(EM_F);
         THREAD_MEMFREE(row_index);

      } /* end of pointer check */
   } /* end parallel region */
}
/*
 * $Log$
 */
