
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

/* Paso: Gauss-Seidel                */

/**************************************************************/

/* Author: artak@uq.edu.au                                   */

/**************************************************************/

#include "Paso.h"
#include "Solver.h"
#include "PasoUtil.h"
#include "BlockOps.h"

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
      MEMFREE(in->buffer);
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
      out->buffer=MEMALLOC( ((size_t) n) * ((size_t)  n_block), double);
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

/*

performs a few sweeps of the  from

S (x_{k} -  x_{k-1}) = b - A x_{k-1}

where x_{0}=0 and S provides some approximatioon of A.

Under MPI S is build on using A_p->mainBlock only.
if GS is local the defect b - A x_{k-1} is calculated using A_p->mainBlock only.

*/

void Paso_Solver_solveGS(Paso_SystemMatrix* A_p, Paso_Solver_GS * gs, double * x, const double * b) 
{
   register dim_t i;
   const dim_t n= (A_p->mainBlock->numRows) * (A_p->mainBlock->row_block_size);
   
   double *b_new = gs->localGS->buffer;
   dim_t sweeps=gs->localGS->sweeps;
   
   if (gs->is_local) {
      Paso_Solver_solveLocalGS(A_p->mainBlock,gs->localGS,x,b);
   } else {
      #pragma omp parallel for private(i) schedule(static)
      for (i=0;i<n;++i) x[i]=b[i];
      
      Paso_Solver_localGSSweep(A_p->mainBlock,gs->localGS,x);
      
      while (sweeps > 1 ) {
	 #pragma omp parallel for private(i) schedule(static)
	 for (i=0;i<n;++i) b_new[i]=b[i];

         Paso_SystemMatrix_MatrixVector((-1.), A_p, x, 1., b_new); /* b_new = b - A*x */
	 
	 Paso_Solver_localGSSweep(A_p->mainBlock,gs->localGS,b_new);
	 
	 #pragma omp parallel for private(i) schedule(static)
	 for (i=0;i<n;++i) x[i]+=b_new[i]; 
	 sweeps--;
      }
      
   }
}
void Paso_Solver_solveLocalGS(Paso_SparseMatrix* A_p, Paso_Solver_LocalGS * gs, double * x, const double * b) 
{
   register dim_t i;
   const dim_t n= (A_p->numRows) * (A_p->row_block_size);
   double *b_new = gs->buffer;
   dim_t sweeps=gs->sweeps;
   
   #pragma omp parallel for private(i) schedule(static)
   for (i=0;i<n;++i) x[i]=b[i];
   
   Paso_Solver_localGSSweep(A_p,gs,x);
   
   while (sweeps > 1 ) {
	 #pragma omp parallel for private(i) schedule(static)
	 for (i=0;i<n;++i) b_new[i]=b[i];
	 
	 Paso_SparseMatrix_MatrixVector_CSC_OFFSET0((-1.), A_p, x, 1., b_new); /* b_new = b - A*x */
	 
	 Paso_Solver_localGSSweep(A_p,gs,b_new);
	 
	 #pragma omp parallel for private(i) schedule(static)
	 for (i=0;i<n;++i) x[i]+=b_new[i];
	 
	 sweeps--;
   }
}

void Paso_Solver_localGSSweep(Paso_SparseMatrix* A, Paso_Solver_LocalGS * gs, double * x) 
{
   #ifdef _OPENMP
   const dim_t nt=omp_get_max_threads();
   #else
   const dim_t nt = 1;
   #endif
   if (nt < 2) {
      Paso_Solver_localGSSweep_sequential(A,gs,x);
   } else {
      Paso_Solver_localGSSweep_colored(A,gs,x);
   }
}

/* inplace Gauss-Seidel sweep in seqential mode: */

void Paso_Solver_localGSSweep_sequential(Paso_SparseMatrix* A_p, Paso_Solver_LocalGS * gs, double * x)
{
   const dim_t n=A_p->numRows;
   const dim_t n_block=A_p->row_block_size;
   const double *diag = gs->diag;
   /* const index_t* pivot = gs->pivot;
   const dim_t block_size=A_p->block_size;  use for block size >3*/
   
   register dim_t i,k;
   register index_t iptr_ik, mm;
   
   const index_t* ptr_main = Paso_SparseMatrix_borrowMainDiagonalPointer(A_p);
   /* forward substitution */
   
   if (n_block==1) {
      Paso_BlockOps_MV_1(&x[0], &diag[0], &x[0]);
      for (i = 1; i < n; ++i) {
	 mm=ptr_main[i];
	 /* x_i=x_i-a_ik*x_k  (with k<i) */
	 for (iptr_ik=A_p->pattern->ptr[i];iptr_ik<mm; ++iptr_ik) {
	    k=A_p->pattern->index[iptr_ik];  
	    Paso_BlockOps_SMV_1(&x[i], &A_p->val[iptr_ik], &x[k]); 
	 }
	 Paso_BlockOps_MV_1(&x[i], &diag[i], &x[i]);
      }
   } else if (n_block==2) {
      Paso_BlockOps_MV_2(&x[0], &diag[0], &x[0]);
      for (i = 1; i < n; ++i) {
	 mm=ptr_main[i];
	 for (iptr_ik=A_p->pattern->ptr[i];iptr_ik<mm; ++iptr_ik) {
	    k=A_p->pattern->index[iptr_ik];                          
	    Paso_BlockOps_SMV_2(&x[2*i], &A_p->val[4*iptr_ik], &x[2*k]);
	 }
	 Paso_BlockOps_MV_2(&x[2*i], &diag[4*i], &x[2*i]);
      }
   } else if (n_block==3) {
      Paso_BlockOps_MV_3(&x[0], &diag[0], &x[0]);
      for (i = 1; i < n; ++i) {
	 mm=ptr_main[i];
	 /* x_i=x_i-a_ik*x_k */
	 for (iptr_ik=A_p->pattern->ptr[i];iptr_ik<mm; ++iptr_ik) {
	    k=A_p->pattern->index[iptr_ik];
	    Paso_BlockOps_SMV_3(&x[3*i], &A_p->val[9*iptr_ik], &x[3*k]);
	 }
	 Paso_BlockOps_MV_3(&x[3*i], &diag[9*i], &x[3*i]); 
      }
   } /* add block size >3 */

   /* backward substitution */
   
   if (n_block==1) {
      for (i = n-2; i > -1; --i) {	       
	    mm=ptr_main[i];
	    Paso_BlockOps_MV_1(&x[i], &A_p->val[mm], &x[i]);
	    for (iptr_ik=mm+1; iptr_ik < A_p->pattern->ptr[i+1]; ++iptr_ik) {
	       k=A_p->pattern->index[iptr_ik];  
	       Paso_BlockOps_SMV_1(&x[i], &A_p->val[iptr_ik], &x[k]);
	    }
	    Paso_BlockOps_MV_1(&x[i], &diag[i], &x[i]);
      }
      
   } else if (n_block==2) {
      for (i = n-2; i > -1; --i) {
	    mm=ptr_main[i];
	    Paso_BlockOps_MV_2(&x[2*i], &A_p->val[4*mm], &x[2*i]);
	    for (iptr_ik=mm+1; iptr_ik < A_p->pattern->ptr[i+1]; ++iptr_ik) {
	       k=A_p->pattern->index[iptr_ik]; 
	       Paso_BlockOps_SMV_2(&x[2*i], &A_p->val[4*iptr_ik], &x[2*k]);
	    }
	    Paso_BlockOps_MV_2(&x[2*i], &diag[i*4], &x[2*i]);
      }
   } else if (n_block==3) {
      for (i = n-2; i > -1; --i) {
	    mm=ptr_main[i];
	    Paso_BlockOps_MV_3(&x[3*i], &A_p->val[9*mm], &x[3*i]);
	 
	    for (iptr_ik=mm+1; iptr_ik < A_p->pattern->ptr[i+1]; ++iptr_ik) {
	       k=A_p->pattern->index[iptr_ik];    
	       Paso_BlockOps_SMV_3(&x[3*i], &A_p->val[9*iptr_ik], &x[3*k]);
	    }
	    Paso_BlockOps_MV_3(&x[3*i], &diag[i*9], &x[3*i]);
      }
      
   } /* add block size >3 */      
   
   return;
}
       
void Paso_Solver_localGSSweep_colored(Paso_SparseMatrix* A_p, Paso_Solver_LocalGS * gs, double * x) 
{
   const dim_t n=A_p->numRows;
   const dim_t n_block=A_p->row_block_size;
   const double *diag = gs->diag;
   index_t* pivot = gs->pivot; 
   const dim_t block_size=A_p->block_size;
   
   register dim_t i,k;
   register index_t color,iptr_ik, mm;
   
   const index_t* coloring = Paso_Pattern_borrowColoringPointer(A_p->pattern);
   const dim_t num_colors = Paso_Pattern_getNumColors(A_p->pattern);
   const index_t* ptr_main = Paso_SparseMatrix_borrowMainDiagonalPointer(A_p);
   
   /* forward substitution */

   
   /* color = 0 */
   if (n_block==1) { 
      #pragma omp parallel for schedule(static) private(i)
      for (i = 0; i < n; ++i) {
	 if (coloring[i]== 0 ) Paso_BlockOps_MV_1(&x[i], &diag[i], &x[i]);
      }
   } else if (n_block==2) {
         #pragma omp parallel for schedule(static) private(i)
	 for (i = 0; i < n; ++i) {
	    if (coloring[i]== 0 ) Paso_BlockOps_MV_2(&x[2*i], &diag[i*4], &x[2*i]);
	 }
    } else if (n_block==3) {
	 #pragma omp parallel for schedule(static) private(i)
	 for (i = 0; i < n; ++i) {
	    if (coloring[i]==0) Paso_BlockOps_MV_3(&x[3*i], &diag[i*9], &x[3*i]);
	 }
   } else {
      #pragma omp parallel for schedule(static) private(i)
      for (i = 0; i < n; ++i) {
	 if (coloring[i]==0) Paso_BlockOps_Solve_N(n_block, &x[n_block*i], &diag[block_size*i], &pivot[n_block*i], &x[n_block*i]);
      }
   }
   
   for (color=1;color<num_colors;++color) {
 
      if (n_block==1) {
	 #pragma omp parallel for schedule(static) private(i,iptr_ik,k)
	 for (i = 0; i < n; ++i) {
	    if (coloring[i]==color) {
	       /* x_i=x_i-a_ik*x_k */                     
	       for (iptr_ik=A_p->pattern->ptr[i];iptr_ik<A_p->pattern->ptr[i+1]; ++iptr_ik) {
		  k=A_p->pattern->index[iptr_ik];
		  if (coloring[k]<color) Paso_BlockOps_SMV_1(&x[i], &A_p->val[iptr_ik], &x[k]); 
	       }
	       Paso_BlockOps_MV_1(&x[i], &diag[i], &x[i]);
	    }
	 }
      } else if (n_block==2) {
	 #pragma omp parallel for schedule(static) private(i,iptr_ik,k)
	 for (i = 0; i < n; ++i) {
	    if (coloring[i]==color) {
	       for (iptr_ik=A_p->pattern->ptr[i];iptr_ik<A_p->pattern->ptr[i+1]; ++iptr_ik) {
		  k=A_p->pattern->index[iptr_ik];                          
		  if (coloring[k]<color) Paso_BlockOps_SMV_2(&x[2*i], &A_p->val[4*iptr_ik], &x[2*k]); 
	       }
	       Paso_BlockOps_MV_2(&x[2*i], &diag[4*i], &x[2*i]);
	    }
	 }
      } else if (n_block==3) {
	 #pragma omp parallel for schedule(static) private(i,iptr_ik,k)
	 for (i = 0; i < n; ++i) {
	    if (coloring[i]==color) {
	       for (iptr_ik=A_p->pattern->ptr[i];iptr_ik<A_p->pattern->ptr[i+1]; ++iptr_ik) {
		  k=A_p->pattern->index[iptr_ik];                          
		  if (coloring[k]<color) Paso_BlockOps_SMV_3(&x[3*i], &A_p->val[9*iptr_ik], &x[3*k]);
	       }
	       Paso_BlockOps_MV_3(&x[3*i], &diag[9*i], &x[3*i]);
	    }
	 }
      } else {
	 #pragma omp parallel for schedule(static) private(i,iptr_ik,k)
	 for (i = 0; i < n; ++i) {
	    if (coloring[i] == color) {
	       for (iptr_ik=A_p->pattern->ptr[i];iptr_ik<A_p->pattern->ptr[i+1]; ++iptr_ik) {
		  k=A_p->pattern->index[iptr_ik];                          
		  if (coloring[k]<color) Paso_BlockOps_SMV_N(n_block, &x[n_block*i], &A_p->val[block_size*iptr_ik], &x[n_block*k]);
	       }
	       Paso_BlockOps_Solve_N(n_block, &x[n_block*i], &diag[block_size*i], &pivot[n_block*i], &x[n_block*i]);
	    }
	 }
      }
   } /* end of coloring loop */
   
   /* backward substitution */
   for (color=(num_colors)-2 ;color>-1;--color) { /* Note: color=(num_colors)-1 is not required */
      if (n_block==1) {
	 #pragma omp parallel for schedule(static) private(mm, i,iptr_ik,k)
	 for (i = 0; i < n; ++i) {
	    if (coloring[i]==color) {
	       mm=ptr_main[i];
	       Paso_BlockOps_MV_1(&x[i], &A_p->val[mm], &x[i]);
	       for (iptr_ik=A_p->pattern->ptr[i];iptr_ik<A_p->pattern->ptr[i+1]; ++iptr_ik) {
		  k=A_p->pattern->index[iptr_ik];                          
		  if (coloring[k]>color) Paso_BlockOps_SMV_1(&x[i], &A_p->val[iptr_ik], &x[k]);
	       }
	       Paso_BlockOps_MV_1(&x[i], &diag[i], &x[i]);
	    }
	 }
      } else if (n_block==2) {
	 #pragma omp parallel for schedule(static) private(mm, i,iptr_ik,k)
	 for (i = 0; i < n; ++i) {
	    if (coloring[i]==color) {
	       mm=ptr_main[i];
	       Paso_BlockOps_MV_2(&x[2*i], &A_p->val[4*mm], &x[2*i]);
	       for (iptr_ik=A_p->pattern->ptr[i];iptr_ik<A_p->pattern->ptr[i+1]; ++iptr_ik) {
		  k=A_p->pattern->index[iptr_ik];                          
		  if (coloring[k]>color) Paso_BlockOps_SMV_2(&x[2*i], &A_p->val[4*iptr_ik], &x[2*k]);
	       }
	       Paso_BlockOps_MV_2(&x[2*i], &diag[4*i], &x[2*i]);
	    }
	 }
      } else if (n_block==3) {
	 #pragma omp parallel for schedule(static) private(mm, i,iptr_ik,k)
	 for (i = 0; i < n; ++i) {
	    if (coloring[i]==color) {
	       mm=ptr_main[i];
	       Paso_BlockOps_MV_3(&x[3*i], &A_p->val[9*mm], &x[3*i]);
	       for (iptr_ik=A_p->pattern->ptr[i];iptr_ik<A_p->pattern->ptr[i+1]; ++iptr_ik) {
		  k=A_p->pattern->index[iptr_ik];                          
		  if (coloring[k]>color) Paso_BlockOps_SMV_3(&x[3*i], &A_p->val[9*iptr_ik], &x[3*k]);
	       }
	       Paso_BlockOps_MV_3(&x[3*i], &diag[9*i], &x[3*i]);
	    }
	 }
      } else {
	 #pragma omp parallel for schedule(static) private(mm, i,iptr_ik,k)
	 for (i = 0; i < n; ++i) {
	    if (coloring[i]==color) {
	       mm=ptr_main[i];
	       Paso_BlockOps_MV_N( n_block, &x[n_block*i], &A_p->val[block_size*mm], &x[n_block*i] );
	       for (iptr_ik=A_p->pattern->ptr[i];iptr_ik<A_p->pattern->ptr[i+1]; ++iptr_ik) {
		  k=A_p->pattern->index[iptr_ik];                          
		  if (coloring[k]>color) Paso_BlockOps_SMV_N(n_block, &x[n_block*i], &A_p->val[block_size*iptr_ik], &x[n_block*k]);
	       }
	       Paso_BlockOps_Solve_N(n_block, &x[n_block*i], &diag[block_size*i], &pivot[n_block*i], &x[n_block*i]);
	    }
	 }
      }
   }
   return;
}


 