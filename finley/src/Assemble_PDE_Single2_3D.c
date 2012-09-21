
/*****************************************************************************
*
* Copyright (c) 2003-2012 by University of Queensland
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

/*    assembles the system of numEq PDEs into the stiffness matrix S and right hand side F  */

/*      -(A_{i,j} u_,j)_i-(B_{i} u)_i+C_{j} u_,j-D u_m  and -(X_,i)_i + Y */

/*    in a 3D domain. The shape functions for test and solution must be identical  */
/*    and row_NS == row_NN                                                         */

/*    Shape of the coefficients: */

/*      A = 3 x 3 */
/*      B = 3   */
/*      C = 3   */
/*      D = scalar  */
/*      X = 3  */
/*      Y = scalar   */


/************************************************************************************/


#include "Assemble.h"
#include "Util.h"
#ifdef _OPENMP
#include <omp.h>
#endif


/************************************************************************************/

void  Finley_Assemble_PDE_Single2_3D(Finley_Assemble_Parameters p,
                                     Finley_ElementFile* elements,
                                     Paso_SystemMatrix* Mat, escriptDataC* F,
                                     escriptDataC* A, escriptDataC* B, escriptDataC* C, escriptDataC* D, escriptDataC* X, escriptDataC* Y) {

    #define DIM 3
    index_t color;
    dim_t e, isub;
    __const double *A_p, *B_p, *C_p, *D_p, *X_p, *Y_p, *A_q, *B_q, *C_q, *D_q, *X_q, *Y_q;
    double *EM_S, *EM_F, *Vol, *DSDX;
    index_t *row_index;
    register dim_t q, s,r;
    register double rtmp, rtmp00, rtmp01, rtmp02, rtmp10, rtmp11, rtmp12, rtmp20, rtmp21, rtmp22, rtmp0, rtmp1, rtmp2;
    bool_t add_EM_F, add_EM_S;

    bool_t extendedA=isExpanded(A);
    bool_t extendedB=isExpanded(B);
    bool_t extendedC=isExpanded(C);
    bool_t extendedD=isExpanded(D);
    bool_t extendedX=isExpanded(X);
    bool_t extendedY=isExpanded(Y);
    double *F_p=(requireWrite(F), getSampleDataRW(F,0));	/* use comma, to get around the mixed code and declarations thing */
    double *S=p.row_jac->BasisFunctions->S;
    dim_t len_EM_S=p.row_numShapesTotal*p.col_numShapesTotal;
    dim_t len_EM_F=p.row_numShapesTotal;

    #pragma omp parallel private(color,EM_S, EM_F, Vol, DSDX, A_p, B_p, C_p, D_p, X_p, Y_p, A_q, B_q, C_q, D_q, X_q, Y_q,row_index,q, s,r,rtmp, rtmp00, rtmp01, rtmp02, rtmp10, rtmp11, rtmp12, rtmp20, rtmp21, rtmp22, rtmp0, rtmp1, rtmp2,add_EM_F, add_EM_S, isub)
    {
       EM_S=THREAD_MEMALLOC(len_EM_S,double);
       EM_F=THREAD_MEMALLOC(len_EM_F,double);
       row_index=THREAD_MEMALLOC(p.row_numShapesTotal,index_t);
                                                                                                                                                                                                     
       if (!Finley_checkPtr(EM_S) && !Finley_checkPtr(EM_F) && !Finley_checkPtr(row_index) ) {

          for (color=elements->minColor;color<=elements->maxColor;color++) {
             /*  open loop over all elements: */
             #pragma omp for private(e) schedule(static)
             for(e=0;e<elements->numElements;e++){
                if (elements->Color[e]==color) {

                   A_p=getSampleDataRO(A,e);
                   B_p=getSampleDataRO(B,e);
                   C_p=getSampleDataRO(C,e);
                   D_p=getSampleDataRO(D,e);
                   X_p=getSampleDataRO(X,e);
                   Y_p=getSampleDataRO(Y,e);

                   for (isub=0; isub<p.numSub; isub++) {

                      Vol=&(p.row_jac->volume[INDEX3(0,isub,e, p.numQuadSub,p.numSub)]);
                      DSDX=&(p.row_jac->DSDX[INDEX5(0,0,0,isub,e,p.row_numShapesTotal,DIM,p.numQuadSub,p.numSub)]);
                      for (q=0;q<len_EM_S;++q) EM_S[q]=0;
                      for (q=0;q<len_EM_F;++q) EM_F[q]=0;
                      add_EM_F=FALSE;
                      add_EM_S=FALSE;

                      /************************************************************************************/
                      /*   process A: */
                      /************************************************************************************/
                      if (NULL!=A_p) {
                         add_EM_S=TRUE;
                         if (extendedA) {
			    A_q=&(A_p[INDEX4(0,0,0,isub, DIM,DIM,p.numQuadSub)]);
                            for (s=0;s<p.row_numShapes;s++) {
                              for (r=0;r<p.col_numShapes;r++) {
                                rtmp=0;
                                for (q=0;q<p.numQuadSub;q++) {
                                      rtmp+=Vol[q]*( DSDX[INDEX3(s,0,q,p.row_numShapesTotal,DIM)]*A_q[INDEX3(0,0,q,DIM,DIM)]*DSDX[INDEX3(r,0,q,p.row_numShapesTotal,DIM)]
                                                   + DSDX[INDEX3(s,0,q,p.row_numShapesTotal,DIM)]*A_q[INDEX3(0,1,q,DIM,DIM)]*DSDX[INDEX3(r,1,q,p.row_numShapesTotal,DIM)]
                                                   + DSDX[INDEX3(s,0,q,p.row_numShapesTotal,DIM)]*A_q[INDEX3(0,2,q,DIM,DIM)]*DSDX[INDEX3(r,2,q,p.row_numShapesTotal,DIM)]
                                                   + DSDX[INDEX3(s,1,q,p.row_numShapesTotal,DIM)]*A_q[INDEX3(1,0,q,DIM,DIM)]*DSDX[INDEX3(r,0,q,p.row_numShapesTotal,DIM)]
                                                   + DSDX[INDEX3(s,1,q,p.row_numShapesTotal,DIM)]*A_q[INDEX3(1,1,q,DIM,DIM)]*DSDX[INDEX3(r,1,q,p.row_numShapesTotal,DIM)]
                                                   + DSDX[INDEX3(s,1,q,p.row_numShapesTotal,DIM)]*A_q[INDEX3(1,2,q,DIM,DIM)]*DSDX[INDEX3(r,2,q,p.row_numShapesTotal,DIM)]
                                                   + DSDX[INDEX3(s,2,q,p.row_numShapesTotal,DIM)]*A_q[INDEX3(2,0,q,DIM,DIM)]*DSDX[INDEX3(r,0,q,p.row_numShapesTotal,DIM)]
                                                   + DSDX[INDEX3(s,2,q,p.row_numShapesTotal,DIM)]*A_q[INDEX3(2,1,q,DIM,DIM)]*DSDX[INDEX3(r,1,q,p.row_numShapesTotal,DIM)]
                                                   + DSDX[INDEX3(s,2,q,p.row_numShapesTotal,DIM)]*A_q[INDEX3(2,2,q,DIM,DIM)]*DSDX[INDEX3(r,2,q,p.row_numShapesTotal,DIM)]);
                                }
                                EM_S[INDEX4(0,0,s,r,p.numEqu,p.numComp,p.row_numShapesTotal)]+=rtmp;
                              }
                            }
                         } else {
                            for (s=0;s<p.row_numShapes;s++) {
                              for (r=0;r<p.col_numShapes;r++) {
                                       rtmp00=0;
                                       rtmp01=0;
                                       rtmp02=0;
                                       rtmp10=0;
                                       rtmp11=0;
                                       rtmp12=0;
                                       rtmp20=0;
                                       rtmp21=0;
                                       rtmp22=0;
                                       for (q=0;q<p.numQuadSub;q++) {
      
                                           rtmp0=Vol[q]*DSDX[INDEX3(s,0,q,p.row_numShapesTotal,DIM)];
                                           rtmp00+=rtmp0*DSDX[INDEX3(r,0,q,p.row_numShapesTotal,DIM)];
                                           rtmp01+=rtmp0*DSDX[INDEX3(r,1,q,p.row_numShapesTotal,DIM)];
                                           rtmp02+=rtmp0*DSDX[INDEX3(r,2,q,p.row_numShapesTotal,DIM)];
      
                                           rtmp1=Vol[q]*DSDX[INDEX3(s,1,q,p.row_numShapesTotal,DIM)];
                                           rtmp10+=rtmp1*DSDX[INDEX3(r,0,q,p.row_numShapesTotal,DIM)];
                                           rtmp11+=rtmp1*DSDX[INDEX3(r,1,q,p.row_numShapesTotal,DIM)];
                                           rtmp12+=rtmp1*DSDX[INDEX3(r,2,q,p.row_numShapesTotal,DIM)];
      
                                           rtmp2=Vol[q]*DSDX[INDEX3(s,2,q,p.row_numShapesTotal,DIM)];
                                           rtmp20+=rtmp2*DSDX[INDEX3(r,0,q,p.row_numShapesTotal,DIM)];
                                           rtmp21+=rtmp2*DSDX[INDEX3(r,1,q,p.row_numShapesTotal,DIM)];
                                           rtmp22+=rtmp2*DSDX[INDEX3(r,2,q,p.row_numShapesTotal,DIM)];
                                       }
                                       EM_S[INDEX4(0,0,s,r,p.numEqu,p.numComp,p.row_numShapesTotal)]+=rtmp00*A_p[INDEX2(0,0,DIM)]
                                                                                         +rtmp01*A_p[INDEX2(0,1,DIM)]
                                                                                         +rtmp02*A_p[INDEX2(0,2,DIM)]
                                                                                         +rtmp10*A_p[INDEX2(1,0,DIM)]
                                                                                         +rtmp11*A_p[INDEX2(1,1,DIM)]
                                                                                         +rtmp12*A_p[INDEX2(1,2,DIM)]
                                                                                         +rtmp20*A_p[INDEX2(2,0,DIM)]
                                                                                         +rtmp21*A_p[INDEX2(2,1,DIM)]
                                                                                         +rtmp22*A_p[INDEX2(2,2,DIM)];
                              }
                            }
                          }
                      }
                      /************************************************************************************/
                      /*   process B: */
                      /************************************************************************************/
                      if (NULL!=B_p) {
                         add_EM_S=TRUE;
                         if (extendedB) {
		            B_q=&(B_p[INDEX3(0,0,isub, DIM,p.numQuadSub)]);
                            for (s=0;s<p.row_numShapes;s++) {
                              for (r=0;r<p.col_numShapes;r++) {
                                rtmp=0;
                                for (q=0;q<p.numQuadSub;q++) {
                                    rtmp+=Vol[q]*S[INDEX2(r,q,p.row_numShapes)]*
                                                 ( DSDX[INDEX3(s,0,q,p.row_numShapesTotal,DIM)]*B_q[INDEX2(0,q,DIM)]
                                                 + DSDX[INDEX3(s,1,q,p.row_numShapesTotal,DIM)]*B_q[INDEX2(1,q,DIM)]
                                                 + DSDX[INDEX3(s,2,q,p.row_numShapesTotal,DIM)]*B_q[INDEX2(2,q,DIM)]);
                                }
                                EM_S[INDEX4(0,0,s,r,p.numEqu,p.numComp,p.row_numShapesTotal)]+=rtmp;
                              }
                            }
                         } else {
                            for (s=0;s<p.row_numShapes;s++) {
                              for (r=0;r<p.col_numShapes;r++) {
                                  rtmp0=0;
                                  rtmp1=0;
                                  rtmp2=0;
                                  for (q=0;q<p.numQuadSub;q++) {
                                      rtmp=Vol[q]*S[INDEX2(r,q,p.row_numShapes)];
                                      rtmp0+=rtmp*DSDX[INDEX3(s,0,q,p.row_numShapesTotal,DIM)];
                                      rtmp1+=rtmp*DSDX[INDEX3(s,1,q,p.row_numShapesTotal,DIM)];
                                      rtmp2+=rtmp*DSDX[INDEX3(s,2,q,p.row_numShapesTotal,DIM)];
                                  }
                                  EM_S[INDEX4(0,0,s,r,p.numEqu,p.numComp,p.row_numShapesTotal)]+=rtmp0*B_p[0]+rtmp1*B_p[1]+rtmp2*B_p[2];
                              }
                            }
                         }
                      }
                      /************************************************************************************/
                      /*   process C: */
                      /************************************************************************************/
                      if (NULL!=C_p) {
                        add_EM_S=TRUE;
                        if (extendedC) {
							C_q=&(C_p[INDEX3(0,0,isub, DIM, p.numQuadSub)]);
                            for (s=0;s<p.row_numShapes;s++) {
                              for (r=0;r<p.col_numShapes;r++) {
                                rtmp=0;
                                for (q=0;q<p.numQuadSub;q++) {
                                   rtmp+=Vol[q]*S[INDEX2(s,q,p.row_numShapes)]* 
                                         ( C_q[INDEX2(0,q,DIM)]*DSDX[INDEX3(r,0,q,p.row_numShapesTotal,DIM)]
                                         + C_q[INDEX2(1,q,DIM)]*DSDX[INDEX3(r,1,q,p.row_numShapesTotal,DIM)]
                                         + C_q[INDEX2(2,q,DIM)]*DSDX[INDEX3(r,2,q,p.row_numShapesTotal,DIM)]);
                                }
                                EM_S[INDEX4(0,0,s,r,p.numEqu,p.numComp,p.row_numShapesTotal)]+=rtmp;
                              }
                            }
                        } else {
                            for (s=0;s<p.row_numShapes;s++) {
                              for (r=0;r<p.col_numShapes;r++) {
                                  rtmp0=0;
                                  rtmp1=0;
                                  rtmp2=0;
                                  for (q=0;q<p.numQuadSub;q++) {
                                      rtmp=Vol[q]*S[INDEX2(s,q,p.row_numShapes)];
                                      rtmp0+=rtmp*DSDX[INDEX3(r,0,q,p.row_numShapesTotal,DIM)];
                                      rtmp1+=rtmp*DSDX[INDEX3(r,1,q,p.row_numShapesTotal,DIM)];
                                      rtmp2+=rtmp*DSDX[INDEX3(r,2,q,p.row_numShapesTotal,DIM)];
                                  }
                                  EM_S[INDEX4(0,0,s,r,p.numEqu,p.numComp,p.row_numShapesTotal)]+=rtmp0*C_p[0]+rtmp1*C_p[1]+rtmp2*C_p[2];
                              }
                            }
                        }
                      }
                      /*********************************************************************************** */
                      /* process D */
                      /************************************************************************************/
                      if (NULL!=D_p) {
                        add_EM_S=TRUE;
                        if (extendedD) {
			    D_q=&(D_p[INDEX2(0,isub, p.numQuadSub)]);
                            for (s=0;s<p.row_numShapes;s++) {
                              for (r=0;r<p.col_numShapes;r++) {
                                rtmp=0;
                                for (q=0;q<p.numQuadSub;q++) rtmp+=Vol[q]*S[INDEX2(s,q,p.row_numShapes)]*D_q[q]*S[INDEX2(r,q,p.row_numShapes)];
                                EM_S[INDEX4(0,0,s,r,p.numEqu,p.numComp,p.row_numShapesTotal)]+=rtmp;
                              }
                            }
                        } else {
                            for (s=0;s<p.row_numShapes;s++) {
                              for (r=0;r<p.col_numShapes;r++) {
                                  rtmp=0;
                                  for (q=0;q<p.numQuadSub;q++) rtmp+=Vol[q]*S[INDEX2(s,q,p.row_numShapes)]*S[INDEX2(r,q,p.row_numShapes)];
                                  EM_S[INDEX4(0,0,s,r,p.numEqu,p.numComp,p.row_numShapesTotal)]+=rtmp*D_p[0];
                              }
                            }
                        }
                      }
                      /************************************************************************************/
                      /*   process X: */
                      /************************************************************************************/
                      if (NULL!=X_p) {
                        add_EM_F=TRUE;
                        if (extendedX) {
			   X_q=&(X_p[INDEX3(0,0,isub, DIM,p.numQuadSub)]);
                           for (s=0;s<p.row_numShapes;s++) {
                              rtmp=0;
                              for (q=0;q<p.numQuadSub;q++) {
                                   rtmp+=Vol[q]*( DSDX[INDEX3(s,0,q,p.row_numShapesTotal,DIM)]*X_q[INDEX2(0,q,DIM)]
                                                + DSDX[INDEX3(s,1,q,p.row_numShapesTotal,DIM)]*X_q[INDEX2(1,q,DIM)]
                                                + DSDX[INDEX3(s,2,q,p.row_numShapesTotal,DIM)]*X_q[INDEX2(2,q,DIM)]);
                              }
                              EM_F[INDEX2(0,s,p.numEqu)]+=rtmp;
                           }
                         } else {
                           for (s=0;s<p.row_numShapes;s++) {
                              rtmp0=0;
                              rtmp1=0;
                              rtmp2=0;
                              for (q=0;q<p.numQuadSub;q++) {
                                   rtmp0+=Vol[q]*DSDX[INDEX3(s,0,q,p.row_numShapesTotal,DIM)];
                                   rtmp1+=Vol[q]*DSDX[INDEX3(s,1,q,p.row_numShapesTotal,DIM)];
                                   rtmp2+=Vol[q]*DSDX[INDEX3(s,2,q,p.row_numShapesTotal,DIM)];
                              }
                              EM_F[INDEX2(0,s,p.numEqu)]+=rtmp0*X_p[0]+rtmp1*X_p[1]+rtmp2*X_p[2];
                           }
                        }
                     }
                     /************************************************************************************/
                     /*   process Y: */
                     /************************************************************************************/
                      if (NULL!=Y_p) {
                        add_EM_F=TRUE;
                        if (extendedY) {
			   Y_q=&(Y_p[INDEX2(0,isub, p.numQuadSub)]);
                           for (s=0;s<p.row_numShapes;s++) {
                              rtmp=0;
                              for (q=0;q<p.numQuadSub;q++) rtmp+=Vol[q]*S[INDEX2(s,q,p.row_numShapes)]*Y_q[q];
                              EM_F[INDEX2(0,s,p.numEqu)]+=rtmp;
                           }
                         } else {
                           for (s=0;s<p.row_numShapes;s++) {
                               rtmp=0;
                               for (q=0;q<p.numQuadSub;q++) rtmp+=Vol[q]*S[INDEX2(s,q,p.row_numShapes)];
                               EM_F[INDEX2(0,s,p.numEqu)]+=rtmp*Y_p[0];
                           }
                         }
                       }
                       /*********************************************************************************************************************/
                       /* add the element matrices onto the matrix and right hand side                                */
                       /*********************************************************************************************************************/

                       for (q=0;q<p.row_numShapesTotal;q++) row_index[q]=p.row_DOF[elements->Nodes[INDEX2(p.row_node[INDEX2(q,isub,p.row_numShapesTotal)],e,p.NN)]];

                       if (add_EM_F) Finley_Util_AddScatter(p.row_numShapesTotal,row_index,p.numEqu,EM_F,F_p, p.row_DOF_UpperBound);
                       if (add_EM_S) Finley_Assemble_addToSystemMatrix(Mat,p.row_numShapesTotal,row_index,p.numEqu,p.col_numShapesTotal,row_index,p.numComp,EM_S);

                  } /* end of isub */
                } /* end color check */
             } /* end element loop */
         } /* end color loop */
           
         THREAD_MEMFREE(EM_S);
         THREAD_MEMFREE(EM_F);
         THREAD_MEMFREE(row_index);

      } /* end of pointer check */
   } /* end parallel region */
}
