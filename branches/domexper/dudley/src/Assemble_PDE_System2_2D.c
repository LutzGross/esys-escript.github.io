
/*******************************************************
*
* Copyright (c) 2003-2010 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


/**************************************************************/

/*    assembles the system of numEq PDEs into the stiffness matrix S right hand side F  */
/*    the shape functions for test and solution must be identical */


/*      -(A_{k,i,m,j} u_m,j)_i-(B_{k,i,m} u_m)_i+C_{k,m,j} u_m,j-D_{k,m} u_m  and -(X_{k,i})_i + Y_k */

/*    u has p.numComp components in a 2D domain. The shape functions for test and solution must be identical  */
/*    and row_NS == row_NN                                                                                  */

/*    Shape of the coefficients: */

/*      A = p.numEqu x 2 x p.numComp x 2 */
/*      B = 2 x p.numEqu x p.numComp  */
/*      C = p.numEqu x 2 x p.numComp  */
/*      D = p.numEqu x p.numComp  */
/*      X = p.numEqu x 2  */
/*      Y = p.numEqu   */


/**************************************************************/


#include "Assemble.h"
#include "Util.h"
#ifdef _OPENMP
#include <omp.h>
#endif


/**************************************************************/

void  Dudley_Assemble_PDE_System2_2D(Assemble_Parameters p, Dudley_ElementFile* elements,
                                     Paso_SystemMatrix* Mat, escriptDataC* F,
                                     escriptDataC* A, escriptDataC* B, escriptDataC* C, escriptDataC* D, escriptDataC* X, escriptDataC* Y) {

    #define DIM 2
    index_t color;
    dim_t e, isub;
    __const double *A_p, *B_p, *C_p, *D_p, *X_p, *Y_p, *A_q, *B_q, *C_q, *D_q, *X_q, *Y_q;
    double *EM_S, *EM_F, *Vol, *DSDX; 
    index_t *row_index;
    register dim_t q, s,r,k,m;
    register double rtmp, rtmp0, rtmp1, rtmp00, rtmp10, rtmp01, rtmp11;
    bool_t add_EM_F, add_EM_S;

    bool_t extendedA=isExpanded(A);
    bool_t extendedB=isExpanded(B);
    bool_t extendedC=isExpanded(C);
    bool_t extendedD=isExpanded(D);
    bool_t extendedX=isExpanded(X);
    bool_t extendedY=isExpanded(Y);
    double *F_p=(requireWrite(F), getSampleDataRW(F,0));	/* use comma, to get around the mixed code and declarations thing */
    double *S=p.row_jac->BasisFunctions->S;
    dim_t len_EM_S=p.row_numShapesTotal*p.col_numShapesTotal*p.numEqu*p.numComp;
    dim_t len_EM_F=p.row_numShapesTotal*p.numEqu;

    #pragma omp parallel private(color,EM_S, EM_F, Vol, DSDX, A_p, B_p, C_p, D_p, X_p, Y_p, A_q, B_q, C_q, D_q, X_q, Y_q,row_index,q, s,r,k,m,rtmp, rtmp0, rtmp1, rtmp00, rtmp10, rtmp01, rtmp11,add_EM_F, add_EM_S, isub)
    {

       EM_S=THREAD_MEMALLOC(len_EM_S,double);
       EM_F=THREAD_MEMALLOC(len_EM_F,double);
       row_index=THREAD_MEMALLOC(p.row_numShapesTotal,index_t);
                                                                                                                                                                                                     
       if (!Dudley_checkPtr(EM_S) && !Dudley_checkPtr(EM_F) && !Dudley_checkPtr(row_index) ) {

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

                      /**************************************************************/
                      /*   process A: */
                      /**************************************************************/
                      A_p=getSampleDataRO(A,e);
                      if (NULL!=A_p) {
                         add_EM_S=TRUE;
                         if (extendedA) {
			    A_q=&(A_p[INDEX6(0,0,0,0,0,isub, p.numEqu,DIM,p.numComp,DIM,p.numQuadSub)]);
                            for (s=0;s<p.row_numShapes;s++) {
                              for (r=0;r<p.col_numShapes;r++) {
                                for (k=0;k<p.numEqu;k++) {
                                  for (m=0;m<p.numComp;m++) {
                                    rtmp=0;
                                    for (q=0;q<p.numQuadSub;q++) {
                                       rtmp+=Vol[q]* (
                                           DSDX[INDEX3(s,0,q,p.row_numShapesTotal,DIM)]*A_q[INDEX5(k,0,m,0,q,p.numEqu,DIM,p.numComp,DIM)]*DSDX[INDEX3(r,0,q,p.row_numShapesTotal,DIM)]
                                          +DSDX[INDEX3(s,0,q,p.row_numShapesTotal,DIM)]*A_q[INDEX5(k,0,m,1,q,p.numEqu,DIM,p.numComp,DIM)]*DSDX[INDEX3(r,1,q,p.row_numShapesTotal,DIM)]
                                          +DSDX[INDEX3(s,1,q,p.row_numShapesTotal,DIM)]*A_q[INDEX5(k,1,m,0,q,p.numEqu,DIM,p.numComp,DIM)]*DSDX[INDEX3(r,0,q,p.row_numShapesTotal,DIM)]
                                          +DSDX[INDEX3(s,1,q,p.row_numShapesTotal,DIM)]*A_q[INDEX5(k,1,m,1,q,p.numEqu,DIM,p.numComp,DIM)]*DSDX[INDEX3(r,1,q,p.row_numShapesTotal,DIM)]);
                                    }
                                    EM_S[INDEX4(k,m,s,r,p.numEqu,p.numComp,p.row_numShapesTotal)]+=rtmp;
                                  }
                                }
                              }
                            }
                         } else {
                            for (s=0;s<p.row_numShapes;s++) {
                              for (r=0;r<p.col_numShapes;r++) {
                                 rtmp00=0;
                                 rtmp01=0;
                                 rtmp10=0;
                                 rtmp11=0;
                                 for (q=0;q<p.numQuadSub;q++) {
                                       rtmp0=Vol[q]*DSDX[INDEX3(s,0,q,p.row_numShapesTotal,DIM)];
                                       rtmp1=Vol[q]*DSDX[INDEX3(s,1,q,p.row_numShapesTotal,DIM)];
                                       rtmp00+=rtmp0*DSDX[INDEX3(r,0,q,p.row_numShapesTotal,DIM)];
                                       rtmp01+=rtmp0*DSDX[INDEX3(r,1,q,p.row_numShapesTotal,DIM)];
                                       rtmp10+=rtmp1*DSDX[INDEX3(r,0,q,p.row_numShapesTotal,DIM)];
                                       rtmp11+=rtmp1*DSDX[INDEX3(r,1,q,p.row_numShapesTotal,DIM)];
                                 }
                                 for (k=0;k<p.numEqu;k++) {
                                     for (m=0;m<p.numComp;m++) {
                                        EM_S[INDEX4(k,m,s,r,p.numEqu,p.numComp,p.row_numShapesTotal)]+=
                                                   rtmp00*A_p[INDEX4(k,0,m,0,p.numEqu,DIM,p.numComp)]
                                                  +rtmp01*A_p[INDEX4(k,0,m,1,p.numEqu,DIM,p.numComp)]
                                                  +rtmp10*A_p[INDEX4(k,1,m,0,p.numEqu,DIM,p.numComp)]
                                                  +rtmp11*A_p[INDEX4(k,1,m,1,p.numEqu,DIM,p.numComp)];
                                     }
                                 }
                              }
                            }
                        }
                      }
                      /**************************************************************/
                      /*   process B: */
                      /**************************************************************/
                      B_p=getSampleDataRO(B,e);
                      if (NULL!=B_p) {
                         add_EM_S=TRUE;
                         if (extendedB) {
		            B_q=&(B_p[INDEX5(0,0,0,0,isub ,p.numEqu,DIM,p.numComp,p.numQuadSub)]);
                            for (s=0;s<p.row_numShapes;s++) {
                              for (r=0;r<p.col_numShapes;r++) {
                                for (k=0;k<p.numEqu;k++) {
                                  for (m=0;m<p.numComp;m++) {
                                     rtmp=0;
                                     for (q=0;q<p.numQuadSub;q++) {
                                        rtmp+=Vol[q]*S[INDEX2(r,q,p.row_numShapes)]* 
                                             ( DSDX[INDEX3(s,0,q,p.row_numShapesTotal,DIM)]*B_q[INDEX4(k,0,m,q,p.numEqu,DIM,p.numComp)]
                                             + DSDX[INDEX3(s,1,q,p.row_numShapesTotal,DIM)]*B_q[INDEX4(k,1,m,q,p.numEqu,DIM,p.numComp)]);
                                     }
                                     EM_S[INDEX4(k,m,s,r,p.numEqu,p.numComp,p.row_numShapesTotal)]+=rtmp;
                                  }
                                }
                              }
                            }
                         } else {
                            for (s=0;s<p.row_numShapes;s++) {
                              for (r=0;r<p.col_numShapes;r++) {
                                  rtmp0=0;
                                  rtmp1=0;
                                  for (q=0;q<p.numQuadSub;q++) {
                                      rtmp=Vol[q]*S[INDEX2(r,q,p.row_numShapes)];
                                      rtmp0+=rtmp*DSDX[INDEX3(s,0,q,p.row_numShapesTotal,DIM)];
                                      rtmp1+=rtmp*DSDX[INDEX3(s,1,q,p.row_numShapesTotal,DIM)];
                                  }
                                  for (k=0;k<p.numEqu;k++) {
                                     for (m=0;m<p.numComp;m++) {
                                        EM_S[INDEX4(k,m,s,r,p.numEqu,p.numComp,p.row_numShapesTotal)]+= rtmp0*B_p[INDEX3(k,0,m,p.numEqu,DIM)]
                                                                                          + rtmp1*B_p[INDEX3(k,1,m,p.numEqu,DIM)];
                                     }
                                  }
                              }
                            }
                         }
                      }
                      /**************************************************************/
                      /*   process C: */
                      /**************************************************************/
                      C_p=getSampleDataRO(C,e);
                      if (NULL!=C_p) {
                        add_EM_S=TRUE;
                        if (extendedC) {
			    C_q=&(C_p[INDEX5(0,0,0,0,isub, p.numEqu,p.numComp,DIM,p.numQuadSub)]);
                            for (s=0;s<p.row_numShapes;s++) {
                              for (r=0;r<p.col_numShapes;r++) {
                                for (k=0;k<p.numEqu;k++) {
                                  for (m=0;m<p.numComp;m++) {
                                     rtmp=0;
                                     for (q=0;q<p.numQuadSub;q++) {
                                         rtmp+=Vol[q]*S[INDEX2(s,q,p.row_numShapes)]*
                                              ( C_q[INDEX4(k,m,0,q,p.numEqu,p.numComp,DIM)]*DSDX[INDEX3(r,0,q,p.row_numShapesTotal,DIM)]
                                              + C_q[INDEX4(k,m,1,q,p.numEqu,p.numComp,DIM)]*DSDX[INDEX3(r,1,q,p.row_numShapesTotal,DIM)]);
                                     }
                                     EM_S[INDEX4(k,m,s,r,p.numEqu,p.numComp,p.row_numShapesTotal)]+=rtmp;
                                  }
                                }
                              }
                            }
                        } else {
                            for (s=0;s<p.row_numShapes;s++) {
                              for (r=0;r<p.col_numShapes;r++) {
                                 rtmp0=0;
                                 rtmp1=0;
                                 for (q=0;q<p.numQuadSub;q++) {
                                         rtmp=Vol[q]*S[INDEX2(s,q,p.row_numShapes)];
                                         rtmp0+=rtmp*DSDX[INDEX3(r,0,q,p.row_numShapesTotal,DIM)];
                                         rtmp1+=rtmp*DSDX[INDEX3(r,1,q,p.row_numShapesTotal,DIM)];
                                 }
                                 for (k=0;k<p.numEqu;k++) {
                                    for (m=0;m<p.numComp;m++) {
                                           EM_S[INDEX4(k,m,s,r,p.numEqu,p.numComp,p.row_numShapesTotal)]+=rtmp0*C_p[INDEX3(k,m,0,p.numEqu,p.numComp)]
                                                                                             +rtmp1*C_p[INDEX3(k,m,1,p.numEqu,p.numComp)];
                                    }
                                 }
                              }
                            }
                        }
                      }
                      /************************************************************* */
                      /* process D */
                      /**************************************************************/
                      D_p=getSampleDataRO(D,e);
                      if (NULL!=D_p) {
                        add_EM_S=TRUE;
                        if (extendedD) {
			    D_q=&(D_p[INDEX4(0,0,0,isub, p.numEqu,p.numComp,p.numQuadSub)]);
                            for (s=0;s<p.row_numShapes;s++) {
                              for (r=0;r<p.col_numShapes;r++) {
                                for (k=0;k<p.numEqu;k++) {
                                  for (m=0;m<p.numComp;m++) {
                                    rtmp=0;
                                    for (q=0;q<p.numQuadSub;q++) {
                                        rtmp+=Vol[q]*S[INDEX2(s,q,p.row_numShapes)]*D_q[INDEX3(k,m,q,p.numEqu,p.numComp)]*S[INDEX2(r,q,p.row_numShapes)];
                                    }
                                    EM_S[INDEX4(k,m,s,r,p.numEqu,p.numComp,p.row_numShapesTotal)]+=rtmp;
                                  }
                                }
                              }
                            }
                        } else {
                            for (s=0;s<p.row_numShapes;s++) {
                              for (r=0;r<p.col_numShapes;r++) {
                                  rtmp=0;
                                  for (q=0;q<p.numQuadSub;q++) rtmp+=Vol[q]*S[INDEX2(s,q,p.row_numShapes)]*S[INDEX2(r,q,p.row_numShapes)];
                                  for (k=0;k<p.numEqu;k++) {
                                      for (m=0;m<p.numComp;m++) {
                                        EM_S[INDEX4(k,m,s,r,p.numEqu,p.numComp,p.row_numShapesTotal)]+=rtmp*D_p[INDEX2(k,m,p.numEqu)];
                                     }
                                  }
                              }
                            }
                        }
                      }
                      /**************************************************************/
                      /*   process X: */
                      /**************************************************************/
                      X_p=getSampleDataRO(X,e);
                      if (NULL!=X_p) {
                        add_EM_F=TRUE;
                        if (extendedX) {
			   X_q=&(X_p[INDEX4(0,0,0,isub, p.numEqu,DIM,p.numQuadSub)]);
                           for (s=0;s<p.row_numShapes;s++) {
                              for (k=0;k<p.numEqu;k++) {
                                rtmp=0;
                                for (q=0;q<p.numQuadSub;q++) {
                                     rtmp+=Vol[q]*(DSDX[INDEX3(s,0,q,p.row_numShapesTotal,DIM)]*X_q[INDEX3(k,0,q,p.numEqu,DIM)]
                                                  +DSDX[INDEX3(s,1,q,p.row_numShapesTotal,DIM)]*X_q[INDEX3(k,1,q,p.numEqu,DIM)]);
                                }
                                EM_F[INDEX2(k,s,p.numEqu)]+=rtmp;
                              }
                           }
                        } else {
                           for (s=0;s<p.row_numShapes;s++) {
                                 rtmp0=0;
                                 rtmp1=0;
                                 for (q=0;q<p.numQuadSub;q++) {
                                     rtmp0+=Vol[q]*DSDX[INDEX3(s,0,q,p.row_numShapesTotal,DIM)];
                                     rtmp1+=Vol[q]*DSDX[INDEX3(s,1,q,p.row_numShapesTotal,DIM)];
                                 }
                                 for (k=0;k<p.numEqu;k++) EM_F[INDEX2(k,s,p.numEqu)]+=rtmp0*X_p[INDEX2(k,0,p.numEqu)]+rtmp1*X_p[INDEX2(k,1,p.numEqu)];
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
			   Y_q=&(Y_p[INDEX3(0,0,isub, p.numEqu,p.numQuadSub)]);
                           for (s=0;s<p.row_numShapes;s++) {
                              for (k=0;k<p.numEqu;k++) {
                                 rtmp=0;
                                 for (q=0;q<p.numQuadSub;q++) rtmp+=Vol[q]*S[INDEX2(s,q,p.row_numShapes)]*Y_q[INDEX2(k,q,p.numEqu)];
                                 EM_F[INDEX2(k,s,p.numEqu)]+=rtmp;
                              }
                           }
                         } else {
                           for (s=0;s<p.row_numShapes;s++) {
                               rtmp=0;
                               for (q=0;q<p.numQuadSub;q++) rtmp+=Vol[q]*S[INDEX2(s,q,p.row_numShapes)];
                               for (k=0;k<p.numEqu;k++) EM_F[INDEX2(k,s,p.numEqu)]+=rtmp*Y_p[k];
                           }
                         }
                       }
                       /***********************************************************************************************/
                       /* add the element matrices onto the matrix and right hand side                                */
                       /***********************************************************************************************/
                       for (q=0;q<p.row_numShapesTotal;q++) row_index[q]=p.row_DOF[elements->Nodes[INDEX2(p.row_node[INDEX2(q,isub,p.row_numShapesTotal)],e,p.NN)]];

                       if (add_EM_F) Dudley_Util_AddScatter(p.row_numShapesTotal,row_index,p.numEqu,EM_F,F_p, p.row_DOF_UpperBound);
                       if (add_EM_S) Dudley_Assemble_addToSystemMatrix(Mat,p.row_numShapesTotal,row_index,p.numEqu,p.col_numShapesTotal,row_index,p.numComp,EM_S);

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
/*
 * $Log$
 */
