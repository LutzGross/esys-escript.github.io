
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
	Paso_Solver_LocalGS_free(in->localGS);
        MEMFREE(in);
     }
}
void Paso_Solver_LocalGS_free(Paso_Solver_LocalGS * in) {
   if (in!=NULL) {
      MEMFREE(in->diag);
      MEMFREE(in->pivot); 
      MEMFREE(in);
   }
}
/**************************************************************/

/*   constructs the symmetric Gauss-Seidel preconditioner     

*/
Paso_Solver_GS* Paso_Solver_getGS(Paso_SystemMatrix * A_p, dim_t sweeps, bool_t is_local, bool_t verbose) 
{
  
  /* allocations: */  
  Paso_Solver_GS* out=MEMALLOC(1,Paso_Solver_GS);
  if (! Paso_checkPtr(out)) {
     out->localGS=Paso_Solver_getLocalGS(A_p->mainBlock,sweeps,verbose);
     out->is_local=is_local;
  }
  if (Paso_MPIInfo_noError(A_p->mpi_info)) {
     return out;
  } else {
     Paso_Solver_GS_free(out);
     return NULL;
  }
}
Paso_Solver_LocalGS* Paso_Solver_getLocalGS(Paso_SparseMatrix * A_p, dim_t sweeps, bool_t verbose) {
   
   dim_t n=A_p->numRows;
   dim_t n_block=A_p->row_block_size;
   dim_t block_size=A_p->block_size;
   
   double time0=Paso_timer();
   /* allocations: */  
   Paso_Solver_LocalGS* out=MEMALLOC(1,Paso_Solver_LocalGS);
   if (! Paso_checkPtr(out)) {
      
      out->diag=MEMALLOC( ((size_t) n) * ((size_t) block_size),double);
      out->pivot=MEMALLOC( ((size_t) n) * ((size_t)  n_block), index_t);
      out->sweeps=sweeps;
      
      if ( ! ( Paso_checkPtr(out->diag) || Paso_checkPtr(out->pivot) ) ) {
	 Paso_SparseMatrix_invMain(A_p, out->diag, out->pivot );
      }
      
   }
   time0=Paso_timer()-time0;
   
   if (Paso_noError()) {
      if (verbose) printf("timing: Gauss-Seidel preparation: elemination : %e\n",time0);
      return out;
   } else {
      Paso_Solver_LocalGS_free(out);
      return NULL;
   }
}

/**************************************************************/

/* Apply Gauss-Seidel                                

   in fact it solves Ax=b in two steps:
   
   step1: among different MPI ranks, we use block Jacobi

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
A on rank i. Matrix (E+F) is represented as the coupleBlock
of matrix A on each rank (annotated as ACi). 

Therefore, step1 can be turned into the following for rank i:

=> Ai * x{k} = b - ACi * x{k-1} 

where both x{k} and b are the segment of x and b on node i, 
and x{k-1} is the old segment values of x on all other nodes. 

step2: inside rank i, we use Gauss-Seidel

let b'= b - ACi * x{k-1} we have Ai * x{k} = b' for rank i
by using symetrix Gauss-Seidel, 

this can be solved in a forward phase and a backward phase:

   forward phase:  x{m} = diag(Ai){-1} (b' - E*x{m} - F*x{m-1})
   backward phase: x{m+1} = diag(Ai){-1} (b' - F*{m+1} - E*x{m})
   
*/

void Paso_Solver_solveGS(Paso_SystemMatrix* A_p, Paso_Solver_GS * gs, double * x, const double * b) 
{
   register dim_t i;
   const dim_t n=A_p->mainBlock->numRows;
   const dim_t n_block=A_p->mainBlock->row_block_size;
   double *remote_x=NULL;
   const dim_t sweeps_ref=gs->localGS->sweeps;
   dim_t sweeps=gs->localGS->sweeps;
   const bool_t remote = (!gs->is_local) && (A_p->mpi_info->size > 1);
   double *new_b = NULL;
   
   if (remote) {
      new_b=TMPMEMALLOC(n*n_block,double);
      gs->localGS->sweeps=(dim_t) ceil( sqrt(DBLE(gs->localGS->sweeps)) );
   }

   
   Paso_Solver_solveLocalGS(A_p->mainBlock,gs->localGS,x,b);
   sweeps-=gs->localGS->sweeps;
   
   while (sweeps > 0 && remote ) {
         Paso_SystemMatrix_startCollect(A_p,x);
	 /* calculate new_b = b - ACi * x{k-1}, where x{k-1} are remote
	    value of x, which requires MPI communications */
	 #pragma omp parallel for private(i) schedule(static)
	 for (i=0;i<n*n_block;++i) new_b[i]=b[i];
	    
	 remote_x=Paso_SystemMatrix_finishCollect(A_p);
	 /*new_b = (-1) * AC * x + 1. * new_b */
	 Paso_SparseMatrix_MatrixVector_CSR_OFFSET0(DBLE(-1),A_p->col_coupleBlock,remote_x,DBLE(1), new_b);
	 
	 Paso_Solver_solveLocalGS(A_p->mainBlock,gs->localGS,x,new_b);
	 sweeps-=gs->localGS->sweeps;
   }
   TMPMEMFREE(new_b);
   gs->localGS->sweeps=sweeps_ref;
   return;
}

void Paso_Solver_solveLocalGS(Paso_SparseMatrix* A, Paso_Solver_LocalGS * gs, double * x, const double * b) 
{
   dim_t i;
   #ifdef _OPENMP
   const dim_t nt=omp_get_max_threads();
   #else
   const dim_t nt = 1;
   #endif
   gs->sweeps=MAX(gs->sweeps,1);
   
   for (i =0 ; i<gs->sweeps; i++) {
      if (nt > 1) {
	 Paso_Solver_solveLocalGS_sequential(A,gs,x,b);
      } else {
	 Paso_Solver_solveLocalGS_colored(A,gs,x,b);
	    /* Paso_Solver_solveLocalGS_tiled(A,gs,x,b); LIN: ADD YOUR STUFF */
      }
   }
}

/*

   applies symmetric Gauss Seidel with coloring = (U+D)^{-1}*D* (L+D)^{-1} 

*/
   
       
void Paso_Solver_solveLocalGS_colored(Paso_SparseMatrix* A_p, Paso_Solver_LocalGS * gs, double * x, const double * b) 
{
   const dim_t n=A_p->numRows;
   const dim_t n_block=A_p->row_block_size;
   const double *diag = gs->diag;
   /* const index_t* pivot = gs->pivot;
      const dim_t block_size=A_p->block_size;  use for block size >3*/
   
   register dim_t i,k;
   register index_t color,iptr_ik, mm;
   register double A11,A12,A21,A22,A13,A23,A33,A32,A31,S1,S2,S3,R1,R2,R3;
   
   const index_t* coloring = Paso_Pattern_borrowColoringPointer(A_p->pattern);
   const dim_t num_colors = Paso_Pattern_getNumColors(A_p->pattern);
   const index_t* ptr_main = Paso_SparseMatrix_borrowMainDiagonalPointer(A_p);
   
   /* forward substitution */
   
   /* color = 0 */
   if (n_block==1) { 
      #pragma omp parallel for schedule(static) private(i,S1, A11)
      for (i = 0; i < n; ++i) {
	   if (coloring[i]==0) {
	       /* x_i=x_i-a_ik*x_k */                     
	       S1=b[i];
	       A11=diag[i];
	       x[i]=A11*S1;
	    }
	 }
   } else if (n_block==2) {
	 #pragma omp parallel for schedule(static) private(i,S1,S2,A11,A21,A12,A22)
	 for (i = 0; i < n; ++i) {
	    if (coloring[i]== 0 ) {
	       /* x_i=x_i-a_ik*x_k */
	       S1=b[2*i];
	       S2=b[2*i+1];

	       A11=diag[i*4];
	       A12=diag[i*4+2];
	       A21=diag[i*4+1];
	       A22=diag[i*4+3];
	       
	       x[2*i  ]=A11 * S1 + A12 * S2;
	       x[2*i+1]=A21 * S1 + A22 * S2;
	       
	    }
	 }
      } else if (n_block==3) {
	 #pragma omp parallel for schedule(static) private(i,S1,S2,S3,A11,A21,A31,A12,A22,A32,A13,A23,A33)
	 for (i = 0; i < n; ++i) {
	    if (coloring[i]==0) {
	       /* x_i=x_i-a_ik*x_k */
	       S1=b[3*i];
	       S2=b[3*i+1];
	       S3=b[3*i+2];
	       A11=diag[i*9  ];
	       A21=diag[i*9+1];
	       A31=diag[i*9+2];
	       A12=diag[i*9+3];
	       A22=diag[i*9+4];
	       A32=diag[i*9+5];
	       A13=diag[i*9+6];
	       A23=diag[i*9+7];
	       A33=diag[i*9+8];
	       x[3*i  ]=A11 * S1 + A12 * S2 + A13 * S3;
	       x[3*i+1]=A21 * S1 + A22 * S2 + A23 * S3;
	       x[3*i+2]=A31 * S1 + A32 * S2 + A33 * S3;
	    }
	 }
   } /* add block size >3 */
   
   for (color=1;color<num_colors;++color) {
      if (n_block==1) {
	 #pragma omp parallel for schedule(static) private(i,iptr_ik,k,S1, A11,R1)
	 for (i = 0; i < n; ++i) {
	    if (coloring[i]==color) {
	       /* x_i=x_i-a_ik*x_k */                     
	       S1=b[i];
	       for (iptr_ik=A_p->pattern->ptr[i];iptr_ik<A_p->pattern->ptr[i+1]; ++iptr_ik) {
		  k=A_p->pattern->index[iptr_ik];                          
		  if (coloring[k]<color) { 
		     R1=x[k];                              
		     S1-=A_p->val[iptr_ik]*R1;
		  }
	       }
	       A11=diag[i];
	       x[i]=A11*S1;
	    }
	 }
      } else if (n_block==2) {
	 #pragma omp parallel for schedule(static) private(i,iptr_ik,k,S1,S2,R1,R2,A11,A21,A12,A22)
	 for (i = 0; i < n; ++i) {
	    if (coloring[i]==color) {
	       /* x_i=x_i-a_ik*x_k */
	       S1=b[2*i];
	       S2=b[2*i+1];
	       for (iptr_ik=A_p->pattern->ptr[i];iptr_ik<A_p->pattern->ptr[i+1]; ++iptr_ik) {
		  k=A_p->pattern->index[iptr_ik];                          
		  if (coloring[k]<color) {
		     R1=x[2*k];
		     R2=x[2*k+1];
		     S1-=A_p->val[4*iptr_ik  ]*R1+A_p->val[4*iptr_ik+2]*R2;
		     S2-=A_p->val[4*iptr_ik+1]*R1+A_p->val[4*iptr_ik+3]*R2;
		  }
	       }
	       A11=diag[i*4];
	       A12=diag[i*4+2];
	       A21=diag[i*4+1];
	       A22=diag[i*4+3];

               x[2*i  ]=A11 * S1 + A12 * S2;
     	       x[2*i+1]=A21 * S1 + A22 * S2;

	    }
	    
	 }
      } else if (n_block==3) {
	 #pragma omp parallel for schedule(static) private(i,iptr_ik,k,S1,S2,S3,R1,R2,R3,A11,A21,A31,A12,A22,A32,A13,A23,A33)
	 for (i = 0; i < n; ++i) {
	    if (coloring[i]==color) {
	       /* x_i=x_i-a_ik*x_k */
	       S1=b[3*i];
	       S2=b[3*i+1];
	       S3=b[3*i+2];
	       for (iptr_ik=A_p->pattern->ptr[i];iptr_ik<A_p->pattern->ptr[i+1]; ++iptr_ik) {
		  k=A_p->pattern->index[iptr_ik];                          
		  if (coloring[k]<color) {
		     R1=x[3*k];
		     R2=x[3*k+1];
		     R3=x[3*k+2];
		     S1-=A_p->val[9*iptr_ik  ]*R1+A_p->val[9*iptr_ik+3]*R2+A_p->val[9*iptr_ik+6]*R3;
		     S2-=A_p->val[9*iptr_ik+1]*R1+A_p->val[9*iptr_ik+4]*R2+A_p->val[9*iptr_ik+7]*R3;
		     S3-=A_p->val[9*iptr_ik+2]*R1+A_p->val[9*iptr_ik+5]*R2+A_p->val[9*iptr_ik+8]*R3;
		  }
	       }
	       A11=diag[i*9  ];
	       A21=diag[i*9+1];
	       A31=diag[i*9+2];
	       A12=diag[i*9+3];
	       A22=diag[i*9+4];
	       A32=diag[i*9+5];
	       A13=diag[i*9+6];
	       A23=diag[i*9+7];
	       A33=diag[i*9+8];
	       x[3*i  ]=A11 * S1 + A12 * S2 + A13 * S3;
	       x[3*i+1]=A21 * S1 + A22 * S2 + A23 * S3;
	       x[3*i+2]=A31 * S1 + A32 * S2 + A33 * S3;
	    }
	 }
      } /* add block size >3 */
   } /* end of coloring loop */
   

   /* backward substitution */
   for (color=(num_colors)-2 ;color>-1;--color) { /* Note: color=(num_colors)-1 is not required */
      if (n_block==1) {
	 #pragma omp parallel for schedule(static) private(mm, i,iptr_ik,k,S1,R1)
	 for (i = 0; i < n; ++i) {
	    if (coloring[i]==color) {
	       
	       mm=ptr_main[i];
	       R1=x[i];
	       A11=A_p->val[mm];
	       S1 = A11 * R1; 
	       
	       for (iptr_ik=A_p->pattern->ptr[i];iptr_ik<A_p->pattern->ptr[i+1]; ++iptr_ik) {
		  k=A_p->pattern->index[iptr_ik];                          
		  if (coloring[k]>color) {
		     R1=x[k]; 
		     S1-=A_p->val[iptr_ik]*R1;
		  }
	       }
	       
	       A11=diag[i];
	       x[i]=A11*S1;
	       
	    }
	 }
      } else if (n_block==2) {
	 #pragma omp parallel for schedule(static) private(mm, i,iptr_ik,k,S1,S2,R1,R2,A11,A21,A12,A22)
	 for (i = 0; i < n; ++i) {
	    if (coloring[i]==color) {
	       
	       mm=ptr_main[i];
	       
	       R1=x[2*i];
	       R2=x[2*i+1];
	       
	       A11=A_p->val[mm*4  ];
	       A21=A_p->val[mm*4+1];
	       A12=A_p->val[mm*4+2];
	       A22=A_p->val[mm*4+3];
	       
	       S1 = A11 * R1 + A12 * R2;
	       S2 = A21 * R1 + A22 * R2;
	       
	       for (iptr_ik=A_p->pattern->ptr[i];iptr_ik<A_p->pattern->ptr[i+1]; ++iptr_ik) {
		  k=A_p->pattern->index[iptr_ik];                          
		  if (coloring[k]>color) {
		     R1=x[2*k];
		     R2=x[2*k+1];
		     S1-=A_p->val[4*iptr_ik  ]*R1+A_p->val[4*iptr_ik+2]*R2;
		     S2-=A_p->val[4*iptr_ik+1]*R1+A_p->val[4*iptr_ik+3]*R2;
		  }
	       }
	       
	       A11=diag[i*4];
	       A12=diag[i*4+2];
	       A21=diag[i*4+1];
	       A22=diag[i*4+3];
	       
	       x[2*i  ]=A11 * S1 + A12 * S2;
	       x[2*i+1]=A21 * S1 + A22 * S2;
	       
	    }
	 }
      } else if (n_block==3) {
	 #pragma omp parallel for schedule(static) private(mm, i,iptr_ik,k,S1,S2,S3,R1,R2,R3,A11,A21,A31,A12,A22,A32,A13,A23,A33)
	 for (i = 0; i < n; ++i) {
	    if (coloring[i]==color) {

	       mm=ptr_main[i];
	       R1=x[3*i];
	       R2=x[3*i+1];
	       R3=x[3*i+2];
	       
	       A11=A_p->val[mm*9  ];
	       A21=A_p->val[mm*9+1];
	       A31=A_p->val[mm*9+2];
	       A12=A_p->val[mm*9+3];
	       A22=A_p->val[mm*9+4];
	       A32=A_p->val[mm*9+5];
	       A13=A_p->val[mm*9+6];
	       A23=A_p->val[mm*9+7];
	       A33=A_p->val[mm*9+8];
	       
	       S1 =A11 * R1 + A12 * R2 + A13 * R3;
	       S2 =A21 * R1 + A22 * R2 + A23 * R3;
	       S3 =A31 * R1 + A32 * R2 + A33 * R3;

	       for (iptr_ik=A_p->pattern->ptr[i];iptr_ik<A_p->pattern->ptr[i+1]; ++iptr_ik) {
		  k=A_p->pattern->index[iptr_ik];                          
		  if (coloring[k]>color) {
		     R1=x[3*k];
		     R2=x[3*k+1];
		     R3=x[3*k+2];
		     S1-=A_p->val[9*iptr_ik  ]*R1+A_p->val[9*iptr_ik+3]*R2+A_p->val[9*iptr_ik+6]*R3;
		     S2-=A_p->val[9*iptr_ik+1]*R1+A_p->val[9*iptr_ik+4]*R2+A_p->val[9*iptr_ik+7]*R3;
		     S3-=A_p->val[9*iptr_ik+2]*R1+A_p->val[9*iptr_ik+5]*R2+A_p->val[9*iptr_ik+8]*R3;
		  }
	       }
	       
	       A11=diag[i*9  ];
	       A21=diag[i*9+1];
	       A31=diag[i*9+2];
	       A12=diag[i*9+3];
	       A22=diag[i*9+4];
	       A32=diag[i*9+5];
	       A13=diag[i*9+6];
	       A23=diag[i*9+7];
	       A33=diag[i*9+8];
	       
	       x[3*i  ]=A11 * S1 + A12 * S2 + A13 * S3;
	       x[3*i+1]=A21 * S1 + A22 * S2 + A23 * S3;
	       x[3*i+2]=A31 * S1 + A32 * S2 + A33 * S3;
	       
	    }
	 }
      } /* add block size >3 */      
   }
   return;
}

void Paso_Solver_solveLocalGS_sequential(Paso_SparseMatrix* A_p, Paso_Solver_LocalGS * gs, double * x, const double * b)
{
      const dim_t n=A_p->numRows;
      const dim_t n_block=A_p->row_block_size;
      const double *diag = gs->diag;
      /* const index_t* pivot = gs->pivot;
      const dim_t block_size=A_p->block_size;  use for block size >3*/
      
      register dim_t i,k;
      register index_t iptr_ik, mm;
      register double A11,A12,A21,A22,A13,A23,A33,A32,A31,S1,S2,S3,R1,R2,R3;
      
      const index_t* ptr_main = Paso_SparseMatrix_borrowMainDiagonalPointer(A_p);
      
      /* forward substitution */
      
	 if (n_block==1) {
	    for (i = 0; i < n; ++i) {
		  /* x_i=x_i-a_ik*x_k */                     
		  S1=b[i];
		  for (iptr_ik=A_p->pattern->ptr[i];iptr_ik<A_p->pattern->ptr[i+1]; ++iptr_ik) {
		     k=A_p->pattern->index[iptr_ik];                          
		     if (k<i) { 
			R1=x[k];                              
			S1-=A_p->val[iptr_ik]*R1;
		     } else {
			break; /* index is ordered */
		     }
		  }
		  A11=diag[i];
		  x[i]=A11*S1;
	       }
	 } else if (n_block==2) {
	    for (i = 0; i < n; ++i) {
		  S1=b[2*i];
		  S2=b[2*i+1];
		  for (iptr_ik=A_p->pattern->ptr[i];iptr_ik<A_p->pattern->ptr[i+1]; ++iptr_ik) {
		     k=A_p->pattern->index[iptr_ik];                          
		     if (k<i) {
			R1=x[2*k];
			R2=x[2*k+1];
			S1-=A_p->val[4*iptr_ik  ]*R1+A_p->val[4*iptr_ik+2]*R2;
			S2-=A_p->val[4*iptr_ik+1]*R1+A_p->val[4*iptr_ik+3]*R2;
		     } else {
			break; /* index is ordered */
		     }
		  }
		  A11=diag[i*4];
		  A12=diag[i*4+2];
		  A21=diag[i*4+1];
		  A22=diag[i*4+3];
		  
		  x[2*i  ]=A11 * S1 + A12 * S2;
		  x[2*i+1]=A21 * S1 + A22 * S2;
		  
	    }
	 } else if (n_block==3) {
	    for (i = 0; i < n; ++i) {
		  /* x_i=x_i-a_ik*x_k */
		  S1=b[3*i];
		  S2=b[3*i+1];
		  S3=b[3*i+2];
		  for (iptr_ik=A_p->pattern->ptr[i];iptr_ik<A_p->pattern->ptr[i+1]; ++iptr_ik) {
		     k=A_p->pattern->index[iptr_ik];                          
		     if ( k<i ) {
			R1=x[3*k];
			R2=x[3*k+1];
			R3=x[3*k+2];
			S1-=A_p->val[9*iptr_ik  ]*R1+A_p->val[9*iptr_ik+3]*R2+A_p->val[9*iptr_ik+6]*R3;
			S2-=A_p->val[9*iptr_ik+1]*R1+A_p->val[9*iptr_ik+4]*R2+A_p->val[9*iptr_ik+7]*R3;
			S3-=A_p->val[9*iptr_ik+2]*R1+A_p->val[9*iptr_ik+5]*R2+A_p->val[9*iptr_ik+8]*R3;
		     } else {
			break; /* index is ordered */
		     }
		  }
		  A11=diag[i*9  ];
		  A21=diag[i*9+1];
		  A31=diag[i*9+2];
		  A12=diag[i*9+3];
		  A22=diag[i*9+4];
		  A32=diag[i*9+5];
		  A13=diag[i*9+6];
		  A23=diag[i*9+7];
		  A33=diag[i*9+8];
		  x[3*i  ]=A11 * S1 + A12 * S2 + A13 * S3;
		  x[3*i+1]=A21 * S1 + A22 * S2 + A23 * S3;
		  x[3*i+2]=A31 * S1 + A32 * S2 + A33 * S3;
	       }
	   
	 } /* add block size >3 */

      
      
      /* backward substitution */

	 if (n_block==1) {
	    for (i = n-2; i > -1; ++i) {
		  
		  mm=ptr_main[i];
		  R1=x[i];
		  A11=A_p->val[mm];
		  S1 = A11 * R1; 
		  
		  for (iptr_ik=A_p->pattern->ptr[i+1]-1; iptr_ik > A_p->pattern->ptr[i]-1; ++iptr_ik) {
		     k=A_p->pattern->index[iptr_ik];                          
		     if (k > i) {
			R1=x[k]; 
			S1-=A_p->val[iptr_ik]*R1;
		     } else {
			break ;
		     }
		  }
		  
		  A11=diag[i];
		  x[i]=A11*S1;
		  
	    }
	 } else if (n_block==2) {
	    for (i = n-2; i > -1; ++i) {
		  
		  mm=ptr_main[i];
		  
		  R1=x[2*i];
		  R2=x[2*i+1];
		  
		  A11=A_p->val[mm*4  ];
		  A21=A_p->val[mm*4+1];
		  A12=A_p->val[mm*4+2];
		  A22=A_p->val[mm*4+3];
		  
		  S1 = A11 * R1 + A12 * R2;
		  S2 = A21 * R1 + A22 * R2;
		  
		  for (iptr_ik=A_p->pattern->ptr[i+1]-1; iptr_ik > A_p->pattern->ptr[i]-1; ++iptr_ik) {
		     k=A_p->pattern->index[iptr_ik];                          
		     if (k > i) {
			R1=x[2*k];
			R2=x[2*k+1];
			S1-=A_p->val[4*iptr_ik  ]*R1+A_p->val[4*iptr_ik+2]*R2;
			S2-=A_p->val[4*iptr_ik+1]*R1+A_p->val[4*iptr_ik+3]*R2;
		     } else {
			break ;
		     }
		  }
		  
		  A11=diag[i*4];
		  A12=diag[i*4+2];
		  A21=diag[i*4+1];
		  A22=diag[i*4+3];
		  
		  x[2*i  ]=A11 * S1 + A12 * S2;
		  x[2*i+1]=A21 * S1 + A22 * S2;
		  
	     
	    }
	 } else if (n_block==3) {
	    for (i = n-2; i > -1; ++i) {
		  
		  mm=ptr_main[i];
		  R1=x[3*i];
		  R2=x[3*i+1];
		  R3=x[3*i+2];
		  
		  A11=A_p->val[mm*9  ];
		  A21=A_p->val[mm*9+1];
		  A31=A_p->val[mm*9+2];
		  A12=A_p->val[mm*9+3];
		  A22=A_p->val[mm*9+4];
		  A32=A_p->val[mm*9+5];
		  A13=A_p->val[mm*9+6];
		  A23=A_p->val[mm*9+7];
		  A33=A_p->val[mm*9+8];
		  
		  S1 =A11 * R1 + A12 * R2 + A13 * R3;
		  S2 =A21 * R1 + A22 * R2 + A23 * R3;
		  S3 =A31 * R1 + A32 * R2 + A33 * R3;
		  
		  for (iptr_ik=A_p->pattern->ptr[i+1]-1; iptr_ik > A_p->pattern->ptr[i]-1; ++iptr_ik) {
		     k=A_p->pattern->index[iptr_ik];                          
		     if (k > i) {
			R1=x[3*k];
			R2=x[3*k+1];
			R3=x[3*k+2];
			S1-=A_p->val[9*iptr_ik  ]*R1+A_p->val[9*iptr_ik+3]*R2+A_p->val[9*iptr_ik+6]*R3;
			S2-=A_p->val[9*iptr_ik+1]*R1+A_p->val[9*iptr_ik+4]*R2+A_p->val[9*iptr_ik+7]*R3;
			S3-=A_p->val[9*iptr_ik+2]*R1+A_p->val[9*iptr_ik+5]*R2+A_p->val[9*iptr_ik+8]*R3;
		     } else {
			break ;
		     }
		  }
		  
		  A11=diag[i*9  ];
		  A21=diag[i*9+1];
		  A31=diag[i*9+2];
		  A12=diag[i*9+3];
		  A22=diag[i*9+4];
		  A32=diag[i*9+5];
		  A13=diag[i*9+6];
		  A23=diag[i*9+7];
		  A33=diag[i*9+8];
		  
		  x[3*i  ]=A11 * S1 + A12 * S2 + A13 * S3;
		  x[3*i+1]=A21 * S1 + A22 * S2 + A23 * S3;
		  x[3*i+2]=A31 * S1 + A32 * S2 + A33 * S3;
		  
	       }
	    
	 } /* add block size >3 */      

      return;
}
 