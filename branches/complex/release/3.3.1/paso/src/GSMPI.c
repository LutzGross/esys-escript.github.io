
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

/* Paso: GS preconditioner with reordering                    */

/************************************************************************************/

/* Author: l.gao@uq.edu.au                                    */

/************************************************************************************/

#include "Paso.h"
#include "SystemMatrix.h"
#include "Solver.h"
#include "PasoUtil.h"

#include <stdio.h>

/************************************************************************************/

/* free all memory used by GS                                */

void Paso_Solver_GSMPI_free(Paso_Solver_GS * in) {
     if (in!=NULL) {
        MEMFREE(in->colorOf);
        Paso_SparseMatrix_free(in->factors);
        MEMFREE(in->diag);
        MEMFREE(in->main_iptr);
        Paso_Pattern_free(in->pattern);   
        MEMFREE(in);
     }
}

==========================================================================
/************************************************************************************/

/*   gs->diag saves the matrix of D{-1}
     This is different from Paso_Solver_getGS(), in which, gs->diag
     is the matrix D. 
*/
Paso_Solver_GS* Paso_Solver_getGSMPI(Paso_SparseMatrix * A,bool_t verbose) {
  dim_t n=A->numRows;
  dim_t n_block=A->row_block_size;
  dim_t block_size=A->block_size;
  register index_t i,iptr_main=0,iPtr;
  double time0=0,time_color=0,time_fac=0;
  double D, A11, A21, A31, A12, A22, A32, A13, A23, A33;

  /* allocations: */  
/*  printf("n_block= %d, n= %d\n", n_block, n); */
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
                #pragma omp parallel for private(i,iPtr,iptr_main) schedule(static)
                for (i = 0; i < A->pattern->numOutput; i++) {
                   iptr_main=0; 
                   out->diag[i]=1.;
                   /* find main diagonal */
                   for (iPtr = A->pattern->ptr[i]; iPtr < A->pattern->ptr[i + 1]; iPtr++) {
                       if (A->pattern->index[iPtr]==i) {
                           iptr_main=iPtr;
                           if (ABS(A->val[iPtr]) > 0.) {
                                out->diag[i]=1./(A->val[iPtr]);
                           } else {
                                Paso_setError(ZERO_DIVISION_ERROR, "Paso_Solver_getGSMPI: non-regular main diagonal block.");
                           }
                           break;
                       }
                   }
                   out->main_iptr[i]=iptr_main;
                }
             } else if (n_block==2) {
                #pragma omp parallel for private(i,iPtr,iptr_main) schedule(static)
                for (i = 0; i < A->pattern->numOutput; i++) {
                   out->diag[i*4+0]= 1.;
                   out->diag[i*4+1]= 0.;
                   out->diag[i*4+2]= 0.;
                   out->diag[i*4+3]= 1.;
                   iptr_main=0;
                   /* find main diagonal */
                   for (iPtr = A->pattern->ptr[i]; iPtr < A->pattern->ptr[i + 1]; iPtr++) {
                       if (A->pattern->index[iPtr]==i) {
                           iptr_main=iPtr;
                           A11=A->val[iPtr*4];
                           A21=A->val[iPtr*4+1];
                           A12=A->val[iPtr*4+2];
                           A22=A->val[iPtr*4+3];
                           D = A11*A22-A12*A21;
                           if (ABS(D)>0.) {
                                D=1./D;
                                out->diag[i*4  ]=  A22*D;
                                out->diag[i*4+1]= -A21*D;
                                out->diag[i*4+2]= -A12*D;
                                out->diag[i*4+3]=  A11*D;
                           } else {
                                Paso_setError(ZERO_DIVISION_ERROR, "Paso_Solver_getGSMPI: non-regular main diagonal block.");
                           }
                           break;
                       }
                   }
                   out->main_iptr[i]=iptr_main;
                }  
             } else if (n_block==3) {
                #pragma omp parallel for private(i, iPtr,iptr_main) schedule(static)
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
                   iptr_main=0;
                   /* find main diagonal */
                   for (iPtr = A->pattern->ptr[i]; iPtr < A->pattern->ptr[i + 1]; iPtr++) {
                       if (A->pattern->index[iPtr]==i) {
                           iptr_main=iPtr;
                           A11=A->val[iPtr*9  ];
                           A21=A->val[iPtr*9+1];
                           A31=A->val[iPtr*9+2];
                           A12=A->val[iPtr*9+3];
                           A22=A->val[iPtr*9+4];
                           A32=A->val[iPtr*9+5];
                           A13=A->val[iPtr*9+6];
                           A23=A->val[iPtr*9+7];
                           A33=A->val[iPtr*9+8];
                           D = A11*(A22*A33-A23*A32) + A12*(A31*A23-A21*A33) + A13*(A21*A32-A31*A22);
                           if (ABS(D)>0.) {
                                  D=1./D;
                                  out->diag[i*9  ]= (A22*A33-A23*A32)*D;
                                  out->diag[i*9+1]= (A31*A23-A21*A33)*D;
                                  out->diag[i*9+2]= (A21*A32-A31*A22)*D;
                                  out->diag[i*9+3]= (A13*A32-A12*A33)*D;
                                  out->diag[i*9+4]= (A11*A33-A31*A13)*D;
                                  out->diag[i*9+5]= (A12*A31-A11*A32)*D;
                                  out->diag[i*9+6]= (A12*A23-A13*A22)*D;
                                  out->diag[i*9+7]= (A13*A21-A11*A23)*D;
                                  out->diag[i*9+8]= (A11*A22-A12*A21)*D;
                           } else {
                                Paso_setError(ZERO_DIVISION_ERROR, "Paso_Solver_getGSMPI: non-regular main diagonal block.");
                           }
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
         printf("GS_MPI: %d color used \n",out->num_colors);
         printf("timing: GS_MPI: coloring/elimination : %e/%e\n",time_color,time_fac);
     }
     return out;
  } else  {
     Paso_Solver_GSMPI_free(out);
     return NULL;
  }
}

void Paso_Solver_GS_local(Paso_SystemMatrix* A, Paso_Solver_GS * gs, double * x, double * b);

/************************************************************************************/

/* Applies MPI versioned GS

     In fact it solves Ax=b in two steps:
     step 1: among different nodes (MPI ranks), we use block Jacobi
	   x{k} = x{k-1} + D{-1}(b-A*x{k-1})
        => D*x{k} = b - (E+F)x{k-1} 
	      where matrix D is (let p be the number of nodes): 
	       --------------------
	       |A1|  |  | ...  |  |
	       --------------------
	       |  |A2|  | ...  |  |
	       --------------------
	       |  |  |A3| ...  |  |
	       --------------------
               |          ...     |
	       --------------------
               |  |  |  | ...  |Ap|
	       --------------------
	      and Ai (i \in [1,p]) represents the mainBlock of matrix 
	      A on node i. Matrix (E+F) is represented as the coupleBlock
              of matrix A on each node (annotated as ACi). 
           Therefore, step 1 can be turned into the following for node i:
       => Ai * x{k} = b - ACi * x{k-1} 
           where both x{k} and b are the segment of x and b on node i, 
	   and x{k-1} is the old segment values of x on all other nodes. 

     step 2: inside node i, we use Gauss-Seidel
         let b'= b - ACi * x{k-1} we have Ai * x{k} = b' for node i
         by using symmetric Gauss-Seidel, this can be solved in a forward
         phase and a backward phase:
	   forward phase:  x{m} = diag(Ai){-1} (b' - E*x{m} - F*x{m-1})
           backward phase: x{m+1} = diag(Ai){-1} (b' - F*{m+1} - E*x{m})
*/

void Paso_Solver_solveGSMPI(Paso_SystemMatrix* A, Paso_Solver_GS * gs, double * x, double * b) {
     register dim_t i;
     dim_t n_block=gs->n_block;
     dim_t n=gs->n;
     dim_t sweeps=gs->sweeps;

     /*xi{0} = 0
       xi{1} = Ai{-1} * bi
       xi{2} = Ai{-1} * (bi - ACi * xj{1})
       ...
       xi{k} = Ai{-1} * (bi - ACi * xj{k-1}) */
     #pragma omp parallel for private(i) schedule(static)
     for (i=0;i<n*n_block;++i) x[i]=0;

     Paso_Solver_GS_local(A,gs,x,b);

     if (sweeps > 1) {
          double *new_b=MEMALLOC(n*n_block,double);
          double *remote_x=NULL;

          while (sweeps > 1) {
               /* calculate new_b = b - ACi * x{k-1}, where x{k-1} are remote
                  values of x, which requires MPI communication */
               #pragma omp parallel for private(i) schedule(static)
               for (i=0;i<n*n_block;++i) new_b[i]=b[i];

               if (A->col_coupleBlock->pattern->ptr!=NULL){
                    Paso_SystemMatrix_startCollect(A,x);
                    remote_x=Paso_SystemMatrix_finishCollect(A);
                    Paso_SparseMatrix_MatrixVector_CSR_OFFSET0(DBLE(-1),A->col_coupleBlock,remote_x,DBLE(1), new_b);
               }

               Paso_Solver_GS_local(A,gs,x,new_b);
               sweeps --;
          }
          MEMFREE(new_b);
     }

     return;
}

/* Locally solve A'x=b, where A' is the mainBlock of global system matrix A */
void Paso_Solver_GS_local(Paso_SystemMatrix* A, Paso_Solver_GS * gs, double * x, double * b) {
     dim_t n_block=gs->n_block;
     dim_t n=gs->n;
     double sum0, sum1, sum2, X0, X1, X2;
     double *val=A->mainBlock->val;
     double *diag=gs->diag;
     index_t *ptr=gs->pattern->ptr;
     index_t *index=gs->pattern->index;
     dim_t i, j, iptr, xi, ai, xj, aj;
#ifdef _OPENMP
     dim_t nt, len, rest, t, istart, iend;

     nt=omp_get_max_threads();
     len=n/nt;
     rest=n-len*nt;
#endif
     /* TO BE DONE: add handler to deal with the case "n is too small"
                    to be worth run in threads. */

#ifdef _OPENMP
     /* calculate new_b = b - ACi * x{k-1}, where x{k-1} are x values
        computed by other threads in previous sweep */
     if (nt > 1) {
     if (n_block == 1){
         #pragma omp parallel for private(t,istart,iend,i,sum0,iptr,j) schedule(static)
         for (t=0; t<nt; t++) {
              istart=len*t+MIN(t,rest);
              iend=istart+len+(t<rest ? 1:0);
              for (i=istart; i<iend; i++){
                   sum0=b[i];
                   for (iptr=ptr[i]; iptr<ptr[i+1]; iptr++){
                        j=index[iptr];
                        if (j<istart || j>=iend){
                            sum0 = sum0 - val[iptr] * x[j];
                        }
                   }
                   b[i]=sum0;
              }
         }
     } else if (n_block == 2) {
         #pragma omp parallel for private(t,istart,iend,i,xi,sum0,sum1,iptr,j,xj,aj,X0,X1) schedule(static)
         for (t=0; t<nt; t++) {
              istart=len*t+MIN(t,rest);
              iend=istart+len+(t<rest ? 1:0);
              for (i=istart; i<iend; i++){
                   xi=2*i;
                   sum0=b[xi];
                   sum1=b[xi+1];
                   for (iptr=ptr[i]; iptr<ptr[i+1]; iptr++){
                        j=index[iptr];
                        if (j<istart || j>=iend){
                            xj=2*j;
                            aj=4*iptr;
                            X0=x[xj];
                            X1=x[xj+1];
                            sum0 = sum0 - val[aj  ]*X0 - val[aj+2]*X1;
                            sum1 = sum1 - val[aj+1]*X0 - val[aj+3]*X1;
                        }
                   }
                   b[xi]=sum0;
                   b[xi+1]=sum1;
              }
         }
     } else if (n_block == 3) {
         #pragma omp parallel for private(t,istart,iend,i,xi,sum0,sum1,sum2,iptr,j,xj,aj,X0,X1,X2) schedule(static)
         for (t=0; t<nt; t++) {
              istart=len*t+MIN(t,rest);
              iend=istart+len+(t<rest ? 1:0);
              for (i=istart; i<iend; i++){
                   xi=3*i;
                   sum0=b[xi];
                   sum1=b[xi+1];
                   sum2=b[xi+2];
                   for (iptr=ptr[i]; iptr<ptr[i+1]; iptr++){
                        j=index[iptr];
                        if (j<istart || j>=iend){
                            xj=3*j;
                            aj=9*iptr;
                            X0=x[xj];
                            X1=x[xj+1];
                            X2=x[xj+2];
                            sum0 = sum0 - val[aj  ]*X0 - val[aj+3]*X1 - val[aj+6]*X2;
                            sum1 = sum1 - val[aj+1]*X0 - val[aj+4]*X1 - val[aj+7]*X2;
                            sum2 = sum2 - val[aj+2]*X0 - val[aj+5]*X1 - val[aj+8]*X2;
                        }
                   }
                   b[xi]=sum0;
                   b[xi+1]=sum1;
                   b[xi+2]=sum2;
              }
         }
     }
     }
#endif

     /* step 1: forward iteration
               x{k} = D{-1}(b - E*x{k} - F*x{k-1}) */
     /* One Gauss-Seidel iteration
        In case of forward iteration x{k} = D{-1}(b - E*x{k} - F*x{k-1})
           => into a loop (without coloring):
            for i in [0,n-1] do     
               x_i = (1/a_ii) *
                 (b_i - \sum{j=0}{i-1}(a_ij*x_j) - \sum{j=i+1}{n-1}(a_ij*x_j))
        where the first "\sum" sums up newly updated values of x elements 
        while the second "\sum" sums up previous (old) values of x elements.
        In case of backward iteration x{k} = D{-1}(b - F*x{k} - E*x{k-1})
     */
     if (n_block == 1){
#ifdef _OPENMP
         #pragma omp parallel for private(t,istart,iend,i,sum0,iptr,j) schedule(static)
         for (t=0; t<nt; t++) {
           istart=len*t+MIN(t,rest);
           iend=istart+len+(t<rest ? 1:0);
           for (i=istart; i<iend; i++){
#else
         for (i=0; i<n; i++) {
#endif
              sum0 = b[i];
              for (iptr=ptr[i]; iptr<ptr[i+1]; ++iptr) {
                   j=index[iptr];
#ifdef _OPENMP
                   if (j >= istart && j < iend && i != j){
#else
                   if (i != j) {
#endif
                       sum0 = sum0 - val[iptr] * x[j];
                   }
              }
              x[i] = sum0*diag[i];
#ifdef _OPENMP
           }
         }
#else
         }
#endif
     } else if (n_block == 2) {
#ifdef _OPENMP
         #pragma omp parallel for private(t,istart,iend,i,xi,ai,sum0,sum1,iptr,j,xj,aj,X0,X1) schedule(static)
         for (t=0; t<nt; t++) {
           istart=len*t+MIN(t,rest);
           iend=istart+len+(t<rest ? 1:0);
           for (i=istart; i<iend; i++){
#else
         for (i=0; i<n; i++) {
#endif
              xi=2*i;
              ai=4*i;
              sum0 = b[xi];
              sum1 = b[xi+1];
              for (iptr=ptr[i]; iptr<ptr[i+1]; ++iptr) {
                   j=index[iptr];
#ifdef _OPENMP
                   if (j >= istart && j < iend && i != j){
#else
                   if (i != j) {
#endif
                       xj=2*j;
                       aj=4*iptr;
                       X0=x[xj];
                       X1=x[xj+1];
                       sum0 = sum0 - val[aj  ]*X0 - val[aj+2]*X1;
                       sum1 = sum1 - val[aj+1]*X0 - val[aj+3]*X1;
                   }
              }
              x[xi  ]=diag[ai  ]*sum0 + diag[ai+2]*sum1;
              x[xi+1]=diag[ai+1]*sum0 + diag[ai+3]*sum1;
#ifdef _OPENMP
           }
         }
#else
         }
#endif 
     } else if (n_block == 3) {
#ifdef _OPENMP
         #pragma omp parallel for private(t,istart,iend,i,xi,ai,sum0,sum1,sum2,iptr,j,xj,aj,X0,X1,X2) schedule(static)
         for (t=0; t<nt; t++) {
           istart=len*t+MIN(t,rest);
           iend=istart+len+(t<rest ? 1:0);
           for (i=istart; i<iend; i++){
#else
         for (i=0; i<n; i++) {
#endif
              xi=3*i;
              ai=9*i;
              sum0 = b[xi];
              sum1 = b[xi+1];
              sum2 = b[xi+2];
              for (iptr=ptr[i]; iptr<ptr[i+1]; ++iptr) {
                   j=index[iptr];
#ifdef _OPENMP
                   if (j >= istart && j < iend && i != j){
#else
                   if (i != j) {
#endif
                       xj=3*j;
                       aj=9*iptr;
                       X0=x[xj];
                       X1=x[xj+1];
                       X2=x[xj+2];
                       sum0 = sum0 - val[aj  ]*X0 - val[aj+3]*X1 - val[aj+6]*X2;
                       sum1 = sum1 - val[aj+1]*X0 - val[aj+4]*X1 - val[aj+7]*X2;
                       sum2 = sum2 - val[aj+2]*X0 - val[aj+5]*X1 - val[aj+8]*X2;
                   }
              }
              x[xi  ] = diag[ai  ]*sum0 + diag[ai+3]*sum1 + diag[ai+6]*sum2;
              x[xi+1] = diag[ai+1]*sum0 + diag[ai+4]*sum1 + diag[ai+7]*sum2;
              x[xi+2] = diag[ai+2]*sum0 + diag[ai+5]*sum1 + diag[ai+8]*sum2;
#ifdef _OPENMP
           }
         }
#else
         }
#endif 
     }

     /* step 2: backward iteration 
               x{k} = D{-1}(b - F*x{k} - E*x{k-1}) */
     if (n_block == 1){
#ifdef _OPENMP
         #pragma omp parallel for private(t,istart,iend,i,sum0,iptr,j) schedule(static)
         for (t=nt-1; t>=0; t--) {
           istart=len*t+MIN(t,rest);
           iend=istart+len+(t<rest ? 1:0);
           for (i=iend-1; i>=istart; i--){
#else
         for (i=n-1; i>=0; i--) {
#endif
              sum0 = b[i];
              for (iptr=ptr[i]; iptr<ptr[i+1]; ++iptr) {
                   j=index[iptr];
#ifdef _OPENMP
                   if (j >= istart && j < iend && i != j){
#else
                   if (i != j) {
#endif
                       sum0 = sum0 - val[iptr] * x[j];
                   }
              }
              x[i] = sum0*diag[i];
#ifdef _OPENMP
           }
         }
#else
         }
#endif
     } else if (n_block == 2) {
#ifdef _OPENMP
         #pragma omp parallel for private(t,istart,iend,i,xi,ai,sum0,sum1,iptr,j,xj,aj,X0,X1) schedule(static)
         for (t=nt-1; t>=0; t--) {
           istart=len*t+MIN(t,rest);
           iend=istart+len+(t<rest ? 1:0);
           for (i=iend-1; i>=istart; i--){
#else
         for (i=n-1; i>=0; i--) {
#endif
              xi=2*i;
              ai=4*i;
              sum0 = b[xi];
              sum1 = b[xi+1];
              for (iptr=ptr[i]; iptr<ptr[i+1]; ++iptr) {
                   j=index[iptr];
#ifdef _OPENMP
                   if (j >= istart && j < iend && i != j){
#else
                   if (i != j) {
#endif
                       xj=2*j;
                       aj=4*iptr;
                       X0=x[xj];
                       X1=x[xj+1];
                       sum0 = sum0 - val[aj  ]*X0 - val[aj+2]*X1;
                       sum1 = sum1 - val[aj+1]*X0 - val[aj+3]*X1;
                   }
              }
              x[xi  ]=diag[ai  ]*sum0 + diag[ai+2]*sum1;
              x[xi+1]=diag[ai+1]*sum0 + diag[ai+3]*sum1;
#ifdef _OPENMP
           }
         }
#else
         }
#endif
     } else if (n_block == 3) {
#ifdef _OPENMP
         #pragma omp parallel for private(t,istart,iend,i,xi,ai,sum0,sum1,sum2,iptr,j,xj,aj,X0,X1,X2) schedule(static)
         for (t=nt-1; t>=0; t--) {
           istart=len*t+MIN(t,rest);
           iend=istart+len+(t<rest ? 1:0);
           for (i=iend-1; i>=istart; i--){
#else
         for (i=n-1; i>=0; i--) {
#endif
              xi=3*i;
              ai=9*i;
              sum0 = b[xi];
              sum1 = b[xi+1];
              sum2 = b[xi+2];
              for (iptr=ptr[i]; iptr<ptr[i+1]; ++iptr) {
                   j=index[iptr];
#ifdef _OPENMP
                   if (j >= istart && j < iend && i != j){
#else
                   if (i != j) {
#endif
                       xj=3*j;
                       aj=9*iptr;
                       X0=x[xj];
                       X1=x[xj+1];
                       X2=x[xj+2];
                       sum0 = sum0 - val[aj  ]*X0 - val[aj+3]*X1 - val[aj+6]*X2;
                       sum1 = sum1 - val[aj+1]*X0 - val[aj+4]*X1 - val[aj+7]*X2;
                       sum2 = sum2 - val[aj+2]*X0 - val[aj+5]*X1 - val[aj+8]*X2;
                   }
              }
              x[xi  ] = diag[ai  ]*sum0 + diag[ai+3]*sum1 + diag[ai+6]*sum2;
              x[xi+1] = diag[ai+1]*sum0 + diag[ai+4]*sum1 + diag[ai+7]*sum2;
              x[xi+2] = diag[ai+2]*sum0 + diag[ai+5]*sum1 + diag[ai+8]*sum2;
#ifdef _OPENMP
           }
         }
#else
         }
#endif
     }
}

