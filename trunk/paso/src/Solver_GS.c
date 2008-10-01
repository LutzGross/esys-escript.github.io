
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

/* Paso: GS preconditioner with reordering                 */

/**************************************************************/

/* Copyrights by ACcESS Australia 2003,2004,2005,2006,2007,2008  */
/* Author: artak@uq.edu.au                                   */

/**************************************************************/

#include "Paso.h"
#include "Solver.h"
#include "PasoUtil.h"

#include <stdio.h>

/**************************************************************/

/* free all memory used by GS                                */

void Paso_Solver_GS_free(Paso_Solver_GS * in) {
     if (in!=NULL) {
        MEMFREE(in->colorOf);
        Paso_SparseMatrix_free(in->factors);
        MEMFREE(in->diag);
        MEMFREE(in->main_iptr);
        Paso_Pattern_free(in->pattern);   
        MEMFREE(in);
     }
}

/**************************************************************/

/*   constructs the incomplete block factorization of 

*/
Paso_Solver_GS* Paso_Solver_getGS(Paso_SparseMatrix * A,bool_t verbose) {
  dim_t n=A->numRows;
  dim_t n_block=A->row_block_size;
  dim_t block_size=A->block_size;
  index_t num_colors=0, *mis_marker=NULL;
  register index_t i,iptr_main,iPtr,iptr_ik,k,iptr_kj,j,iptr_ij,color,color2;
  double time0,time_color,time_fac;
  /* allocations: */  
  Paso_Solver_GS* out=MEMALLOC(1,Paso_Solver_GS);
  if (Paso_checkPtr(out)) return NULL;
  out->colorOf=MEMALLOC(n,index_t);
  out->diag=MEMALLOC( ((size_t) n) * ((size_t) block_size),double);
  /*out->diag=MEMALLOC(A->len,double);*/
  out->main_iptr=MEMALLOC(n,index_t);
  out->pattern=Paso_Pattern_getReference(A->pattern);
  out->factors=Paso_SparseMatrix_getReference(A);
  out->n_block=n_block;
  out->n=n;
 
  if ( !(Paso_checkPtr(out->colorOf) || Paso_checkPtr(out->main_iptr) || Paso_checkPtr(out->factors)) ) {
    time0=Paso_timer();
    Paso_Pattern_color(A->pattern,&out->num_colors,out->colorOf);
    time_color=Paso_timer()-time0;

    if (Paso_noError()) {
       time0=Paso_timer();

          if (! (Paso_checkPtr(out->diag))) {
             if (n_block==1) {
                #pragma omp parallel for private(i, iPtr) schedule(static)
                for (i = 0; i < A->pattern->numOutput; i++) {
                   out->diag[i]=1.;
                   /* find main diagonal */
                   for (iPtr = A->pattern->ptr[i]; iPtr < A->pattern->ptr[i + 1]; iPtr++) {
                       if (A->pattern->index[iPtr]==i) {
                           iptr_main=iPtr;
                           out->diag[i]=A->val[iPtr];
                           break;
                       }
                   }
                   out->main_iptr[i]=iptr_main;
                }
             } else if (n_block==2) {
                #pragma omp parallel for private(i, iPtr) schedule(static)
                for (i = 0; i < A->pattern->numOutput; i++) {
                   out->diag[i*4+0]= 1.;
                   out->diag[i*4+1]= 0.;
                   out->diag[i*4+2]= 0.;
                   out->diag[i*4+3]= 1.;
                   /* find main diagonal */
                   for (iPtr = A->pattern->ptr[i]; iPtr < A->pattern->ptr[i + 1]; iPtr++) {
                       if (A->pattern->index[iPtr]==i) {
                             iptr_main=iPtr;
                             out->diag[i*4]= A->val[iPtr*4];
                             out->diag[i*4+1]=A->val[iPtr*4+1];
                             out->diag[i*4+2]=A->val[iPtr*4+2];
                             out->diag[i*4+3]= A->val[iPtr*4+3];
                          break;
                       }
                   }
                   out->main_iptr[i]=iptr_main;
                }  
             } else if (n_block==3) {
                #pragma omp parallel for private(i, iPtr) schedule(static)
                for (i = 0; i < A->pattern->numOutput; i++) {
                   out->diag[i*9  ]=1.;
                   out->diag[i*9+1]=0.;
                   out->diag[i*9+2]=0.;
                   out->diag[i*9+3]=0.;
                   out->diag[i*9+4]=1.;
                   out->diag[i*9+5]=0.;
                   out->diag[i*9+6]=0.;
                   out->diag[i*9+7]=0.;
                   out->diag[i*9+8]=1.;
                   /* find main diagonal */
                   for (iPtr = A->pattern->ptr[i]; iPtr < A->pattern->ptr[i + 1]; iPtr++) {
                       if (A->pattern->index[iPtr]==i) {
                           iptr_main=iPtr;
                           out->diag[i*9  ]=A->val[iPtr*9  ];
                           out->diag[i*9+1]=A->val[iPtr*9+1];
                           out->diag[i*9+2]=A->val[iPtr*9+2];
                           out->diag[i*9+3]=A->val[iPtr*9+3];
                           out->diag[i*9+4]=A->val[iPtr*9+4];
                           out->diag[i*9+5]=A->val[iPtr*9+5];
                           out->diag[i*9+6]=A->val[iPtr*9+6];
                           out->diag[i*9+7]=A->val[iPtr*9+7];
                           out->diag[i*9+8]=A->val[iPtr*9+8];
                           break;
                       }
                   }
                   out->main_iptr[i]=iptr_main;
                }
             }
           }

     time_fac=Paso_timer()-time0;
     }
  }
  if (Paso_noError()) {
      if (verbose) {
         printf("GS: %d color used \n",out->num_colors);
         printf("timing: GS: coloring/elemination : %e/%e\n",time_color,time_fac);
     }
     return out;
  } else  {
     Paso_Solver_GS_free(out);
     return NULL;
  }
}

/**************************************************************/

/* apply GS precondition b-> x                               

     in fact it solves LD^{-1}Ux=b in the form x= U^{-1} D L^{-1}b 

 should be called within a parallel region                                              
 barrier synconization should be performed to make sure that the input vector available 

*/

void Paso_Solver_solveGS(Paso_Solver_GS * gs, double * x, double * b) {
     register dim_t i,k;
     register index_t color,iptr_ik,iptr_main;
     register double A11,A12,A21,A22,A13,A23,A33,A32,A31,S1,S2,S3,R1,R2,R3,D,S21,S22,S12,S11,S13,S23,S33,S32,S31;
     dim_t n_block=gs->n_block;
     dim_t n=gs->n;
     index_t* pivot=NULL;
     
     /* copy x into b*/
     #pragma omp parallel for private(i) schedule(static)
     for (i=0;i<n*n_block;++i) x[i]=b[i];
     /* forward substitution */
     for (color=0;color<gs->num_colors;++color) {
           if (n_block==1) {
              #pragma omp parallel for schedule(static) private(i,iptr_ik,k,S1,R1,iptr_main)
              for (i = 0; i < n; ++i) {
                   if (gs->colorOf[i]==color) {
                     /* x_i=x_i-a_ik*x_k */                     
                     S1=x[i];
                     for (iptr_ik=gs->pattern->ptr[i];iptr_ik<gs->pattern->ptr[i+1]; ++iptr_ik) {
                          k=gs->pattern->index[iptr_ik];                          
                          if (gs->colorOf[k]<color) { 
                             R1=x[k];                              
                             S1-=gs->factors->val[iptr_ik]*R1;
                           }
                     }
                     iptr_main=gs->main_iptr[i];
                     x[i]=(1/gs->factors->val[iptr_main])*S1;
                   }
              }
           } else if (n_block==2) {
              #pragma omp parallel for schedule(static) private(i,iptr_ik,k,iptr_main,S1,S2,R1,R2)
              for (i = 0; i < n; ++i) {
                   if (gs->colorOf[i]==color) {
                     /* x_i=x_i-a_ik*x_k */
                     S1=x[2*i];
                     S2=x[2*i+1];
                     for (iptr_ik=gs->pattern->ptr[i];iptr_ik<gs->pattern->ptr[i+1]; ++iptr_ik) {
                          k=gs->pattern->index[iptr_ik];                          
                          if (gs->colorOf[k]<color) {
                             R1=x[2*k];
                             R2=x[2*k+1];
                             S1-=gs->factors->val[4*iptr_ik  ]*R1+gs->factors->val[4*iptr_ik+2]*R2;
                             S2-=gs->factors->val[4*iptr_ik+1]*R1+gs->factors->val[4*iptr_ik+3]*R2;
                          }
                     }
                     iptr_main=gs->main_iptr[i];
                     A11=gs->factors->val[iptr_main*4];
                     A21=gs->factors->val[iptr_main*4+1];
                     A12=gs->factors->val[iptr_main*4+2];
                     A22=gs->factors->val[iptr_main*4+3];
                     D = A11*A22-A12*A21;
                     if (ABS(D)>0.) {
                          D=1./D;
                          S11= A22*D;
                          S21=-A21*D;
                          S12=-A12*D;
                          S22= A11*D;
                          x[2*i  ]=S11*S1+S12*S2;
                          x[2*i+1]=S21*S1+S22*S2;
                     } else {
                            Paso_setError(ZERO_DIVISION_ERROR, "Paso_Solver_getGS: non-regular main diagonal block.");
                       }
                   }

              }
           } else if (n_block==3) {
              #pragma omp parallel for schedule(static) private(i,iptr_ik,iptr_main,k,S1,S2,S3,R1,R2,R3)
              for (i = 0; i < n; ++i) {
                   if (gs->colorOf[i]==color) {
                     /* x_i=x_i-a_ik*x_k */
                     S1=x[3*i];
                     S2=x[3*i+1];
                     S3=x[3*i+2];
                     for (iptr_ik=gs->pattern->ptr[i];iptr_ik<gs->pattern->ptr[i+1]; ++iptr_ik) {
                          k=gs->pattern->index[iptr_ik];                          
                          if (gs->colorOf[k]<color) {
                             R1=x[3*k];
                             R2=x[3*k+1];
                             R3=x[3*k+2];
                             S1-=gs->factors->val[9*iptr_ik  ]*R1+gs->factors->val[9*iptr_ik+3]*R2+gs->factors->val[9*iptr_ik+6]*R3;
                             S2-=gs->factors->val[9*iptr_ik+1]*R1+gs->factors->val[9*iptr_ik+4]*R2+gs->factors->val[9*iptr_ik+7]*R3;
                             S3-=gs->factors->val[9*iptr_ik+2]*R1+gs->factors->val[9*iptr_ik+5]*R2+gs->factors->val[9*iptr_ik+8]*R3;
                          }
                     }
                     iptr_main=gs->main_iptr[i];
                     A11=gs->factors->val[iptr_main*9  ];
                     A21=gs->factors->val[iptr_main*9+1];
                     A31=gs->factors->val[iptr_main*9+2];
                     A12=gs->factors->val[iptr_main*9+3];
                     A22=gs->factors->val[iptr_main*9+4];
                     A32=gs->factors->val[iptr_main*9+5];
                     A13=gs->factors->val[iptr_main*9+6];
                     A23=gs->factors->val[iptr_main*9+7];
                     A33=gs->factors->val[iptr_main*9+8];
                     D  =  A11*(A22*A33-A23*A32)+ A12*(A31*A23-A21*A33)+A13*(A21*A32-A31*A22);
                     if (ABS(D)>0.) {
                          D=1./D;
                          S11=(A22*A33-A23*A32)*D;
                          S21=(A31*A23-A21*A33)*D;
                          S31=(A21*A32-A31*A22)*D;
                          S12=(A13*A32-A12*A33)*D;
                          S22=(A11*A33-A31*A13)*D;
                          S32=(A12*A31-A11*A32)*D;
                          S13=(A12*A23-A13*A22)*D;
                          S23=(A13*A21-A11*A23)*D;
                          S33=(A11*A22-A12*A21)*D;
                             /* a_ik=a_ii^{-1}*a_ik */
                          x[3*i  ]=S11*S1+S12*S2+S13*S3;
                          x[3*i+1]=S21*S1+S22*S2+S23*S3;
                          x[3*i+2]=S31*S1+S32*S2+S33*S3;
                       } else {
                            Paso_setError(ZERO_DIVISION_ERROR, "Paso_Solver_getGS: non-regular main diagonal block.");
                       }
                }
              }
           }
     }
     
     /* Multipling with diag(A) */
     Paso_Solver_applyBlockDiagonalMatrix(gs->n_block,gs->n,gs->diag,pivot,x,x);

     /* backward substitution */
     for (color=(gs->num_colors)-1;color>-1;--color) {
           if (n_block==1) {
              #pragma omp parallel for schedule(static) private(i,iptr_ik,k,S1,R1)
              for (i = 0; i < n; ++i) {
                   if (gs->colorOf[i]==color) {
                     /* x_i=x_i-a_ik*x_k */
                     S1=x[i];
                     for (iptr_ik=gs->pattern->ptr[i];iptr_ik<gs->pattern->ptr[i+1]; ++iptr_ik) {
                          k=gs->pattern->index[iptr_ik];                          
                          if (gs->colorOf[k]>color) {
                             R1=x[k]; 
                             S1-=gs->factors->val[iptr_ik]*R1;
                          }
                     }
                     /*x[i]=S1;*/
                     iptr_main=gs->main_iptr[i];
                     x[i]=(1/gs->factors->val[iptr_main])*S1;
                   }
              }
           } else if (n_block==2) {
              #pragma omp parallel for schedule(static) private(i,iptr_ik,k,S1,S2,R1,R2)
              for (i = 0; i < n; ++i) {
                   if (gs->colorOf[i]==color) {
                     /* x_i=x_i-a_ik*x_k */
                     S1=x[2*i];
                     S2=x[2*i+1];
                     for (iptr_ik=gs->pattern->ptr[i];iptr_ik<gs->pattern->ptr[i+1]; ++iptr_ik) {
                          k=gs->pattern->index[iptr_ik];                          
                          if (gs->colorOf[k]>color) {
                             R1=x[2*k];
                             R2=x[2*k+1];
                             S1-=gs->factors->val[4*iptr_ik  ]*R1+gs->factors->val[4*iptr_ik+2]*R2;
                             S2-=gs->factors->val[4*iptr_ik+1]*R1+gs->factors->val[4*iptr_ik+3]*R2;
                          }
                     }
                     /*x[2*i]=S1;
                     x[2*i+1]=S2;*/
                     iptr_main=gs->main_iptr[i];
                     A11=gs->factors->val[iptr_main*4];
                     A21=gs->factors->val[iptr_main*4+1];
                     A12=gs->factors->val[iptr_main*4+2];
                     A22=gs->factors->val[iptr_main*4+3];
                     D = A11*A22-A12*A21;
                     if (ABS(D)>0.) {
                          D=1./D;
                          S11= A22*D;
                          S21=-A21*D;
                          S12=-A12*D;
                          S22= A11*D;
                          x[2*i  ]=S11*S1+S12*S2;
                          x[2*i+1]=S21*S1+S22*S2;
                     } else {
                            Paso_setError(ZERO_DIVISION_ERROR, "Paso_Solver_getGS: non-regular main diagonal block.");
                       }
 
                    }
              }
           } else if (n_block==3) {
              #pragma omp parallel for schedule(static) private(i,iptr_ik,k,S1,S2,S3,R1,R2,R3)
              for (i = 0; i < n; ++i) {
                   if (gs->colorOf[i]==color) {
                     /* x_i=x_i-a_ik*x_k */
                     S1=x[3*i  ];
                     S2=x[3*i+1];
                     S3=x[3*i+2];
                     for (iptr_ik=gs->pattern->ptr[i];iptr_ik<gs->pattern->ptr[i+1]; ++iptr_ik) {
                          k=gs->pattern->index[iptr_ik];                          
                          if (gs->colorOf[k]>color) {
                             R1=x[3*k];
                             R2=x[3*k+1];
                             R3=x[3*k+2];
                             S1-=gs->factors->val[9*iptr_ik  ]*R1+gs->factors->val[9*iptr_ik+3]*R2+gs->factors->val[9*iptr_ik+6]*R3;
                             S2-=gs->factors->val[9*iptr_ik+1]*R1+gs->factors->val[9*iptr_ik+4]*R2+gs->factors->val[9*iptr_ik+7]*R3;
                             S3-=gs->factors->val[9*iptr_ik+2]*R1+gs->factors->val[9*iptr_ik+5]*R2+gs->factors->val[9*iptr_ik+8]*R3;
                          }
                     }
/*                     x[3*i]=S1;
                     x[3*i+1]=S2;
                     x[3*i+2]=S3;
*/                   iptr_main=gs->main_iptr[i];
                     A11=gs->factors->val[iptr_main*9  ];
                     A21=gs->factors->val[iptr_main*9+1];
                     A31=gs->factors->val[iptr_main*9+2];
                     A12=gs->factors->val[iptr_main*9+3];
                     A22=gs->factors->val[iptr_main*9+4];
                     A32=gs->factors->val[iptr_main*9+5];
                     A13=gs->factors->val[iptr_main*9+6];
                     A23=gs->factors->val[iptr_main*9+7];
                     A33=gs->factors->val[iptr_main*9+8];
                     D  =  A11*(A22*A33-A23*A32)+ A12*(A31*A23-A21*A33)+A13*(A21*A32-A31*A22);
                     if (ABS(D)>0.) {
                          D=1./D;
                          S11=(A22*A33-A23*A32)*D;
                          S21=(A31*A23-A21*A33)*D;
                          S31=(A21*A32-A31*A22)*D;
                          S12=(A13*A32-A12*A33)*D;
                          S22=(A11*A33-A31*A13)*D;
                          S32=(A12*A31-A11*A32)*D;
                          S13=(A12*A23-A13*A22)*D;
                          S23=(A13*A21-A11*A23)*D;
                          S33=(A11*A22-A12*A21)*D;
                          x[3*i  ]=S11*S1+S12*S2+S13*S3;
                          x[3*i+1]=S21*S1+S22*S2+S23*S3;
                          x[3*i+2]=S31*S1+S32*S2+S33*S3;
                       } else {
                            Paso_setError(ZERO_DIVISION_ERROR, "Paso_Solver_getGS: non-regular main diagonal block.");
                       }
                   }
              }
         }
     }
     
     if (gs->sweeps>1) {
     /* Compute the residual b=b-Ax*/
     Paso_SparseMatrix_MatrixVector_CSR_OFFSET0(DBLE(-1), gs->factors, x, DBLE(2), b);
     /* Go round again*/
     gs->sweeps=gs->sweeps-1;
     Paso_Solver_solveGS(gs,x,b);
     }
    
     return;
}

