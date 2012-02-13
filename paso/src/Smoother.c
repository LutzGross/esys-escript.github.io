
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

/* Paso: Gauss-Seidel                                         */

/**************************************************************/

/* Author: artak@uq.edu.au                                    */

/**************************************************************/

#include "Paso.h"
#include "Preconditioner.h"
#include "PasoUtil.h"
#include "BlockOps.h"

#include <stdio.h>


/**************************************************************/

/* free all memory used by Smoother                           */

void Paso_Preconditioner_Smoother_free(Paso_Preconditioner_Smoother * in) {
     if (in!=NULL) {
	Paso_Preconditioner_LocalSmoother_free(in->localSmoother);
        MEMFREE(in);
     }
}
void Paso_Preconditioner_LocalSmoother_free(Paso_Preconditioner_LocalSmoother * in) {
   if (in!=NULL) {
      MEMFREE(in->diag);
      MEMFREE(in->pivot); 
      MEMFREE(in->buffer);
      MEMFREE(in);
   }
}
/**************************************************************/

/*   constructs the symmetric Gauss-Seidel preconditioner     */

Paso_Preconditioner_Smoother* Paso_Preconditioner_Smoother_alloc(Paso_SystemMatrix * A_p, const bool_t jacobi, const bool_t is_local, const bool_t verbose) 
{
  
  /* allocations: */  
  Paso_Preconditioner_Smoother* out=MEMALLOC(1,Paso_Preconditioner_Smoother);
  if (! Esys_checkPtr(out)) {
     out->localSmoother=Paso_Preconditioner_LocalSmoother_alloc(A_p->mainBlock,jacobi,verbose);
     out->is_local=is_local;
  }
  if (Esys_MPIInfo_noError(A_p->mpi_info)) {
     return out;
  } else {
     Paso_Preconditioner_Smoother_free(out);
     return NULL;
  }
}
Paso_Preconditioner_LocalSmoother* Paso_Preconditioner_LocalSmoother_alloc(Paso_SparseMatrix * A_p, const bool_t jacobi, bool_t verbose)
{
   
   const dim_t n=A_p->numRows;
   const dim_t n_block=A_p->row_block_size;
   const dim_t block_size=A_p->block_size;
   
   double time0=Esys_timer();
   /* allocations: */  
   Paso_Preconditioner_LocalSmoother* out=MEMALLOC(1,Paso_Preconditioner_LocalSmoother);
   if (! Esys_checkPtr(out)) {
      
      out->diag=MEMALLOC( ((size_t) n) * ((size_t) block_size),double);
      out->pivot=MEMALLOC( ((size_t) n) * ((size_t)  n_block), index_t);
      out->buffer=MEMALLOC( ((size_t) n) * ((size_t)  n_block), double);
      out->Jacobi=jacobi;
      
      if ( ! ( Esys_checkPtr(out->diag) || Esys_checkPtr(out->pivot) ) ) {
	 Paso_SparseMatrix_invMain(A_p, out->diag, out->pivot );
      }
      
   }
   time0=Esys_timer()-time0;
   
   if (Esys_noError()) {
      return out;
   } else {
      Paso_Preconditioner_LocalSmoother_free(out);
      return NULL;
   }
}

/*

performs a few sweeps of the form

S (x_{k} -  x_{k-1}) = b - A x_{k-1}

where x_{0}=0 and S provides some approximation of A.

Under MPI S is built using A_p->mainBlock only.
If Smoother is local the defect b - A x_{k-1} is calculated using A_p->mainBlock only.
*/

void Paso_Preconditioner_Smoother_solve(Paso_SystemMatrix* A_p, Paso_Preconditioner_Smoother * smoother, double * x, const double * b, 
					const dim_t sweeps, const bool_t x_is_initial) 
{
   const dim_t n= (A_p->mainBlock->numRows) * (A_p->mainBlock->row_block_size);
   
   double *b_new = smoother->localSmoother->buffer;
   dim_t nsweeps=sweeps;
   if (smoother->is_local) {
      Paso_Preconditioner_LocalSmoother_solve(A_p->mainBlock,smoother->localSmoother,x,b,sweeps,x_is_initial);
   } else {
      if (! x_is_initial) {
	 Paso_Copy(n, x, b);
	 Paso_Preconditioner_LocalSmoother_Sweep(A_p->mainBlock,smoother->localSmoother,x);
	 nsweeps--;
      }
      while (nsweeps > 0 ) {
	 Paso_Copy(n, b_new, b);
         Paso_SystemMatrix_MatrixVector((-1.), A_p, x, 1., b_new); /* b_new = b - A*x */
	 Paso_Preconditioner_LocalSmoother_Sweep(A_p->mainBlock,smoother->localSmoother,b_new);	 
	 Paso_AXPY(n, x, 1., b_new); 
	 nsweeps--;
      }
      
   }
}

err_t Paso_Preconditioner_Smoother_solve_byTolerance(Paso_SystemMatrix* A_p, Paso_Preconditioner_Smoother * smoother, 
						    double * x, const double * b, 
						    const double atol, dim_t *sweeps, const bool_t x_is_initial)
{
   const dim_t n= (A_p->mainBlock->numRows) * (A_p->mainBlock->row_block_size);
   double *b_new = smoother->localSmoother->buffer;
   const dim_t max_sweeps=*sweeps;
   dim_t s=0;
   double norm_dx = 2. * atol;
   err_t errorCode = PRECONDITIONER_NO_ERROR;
   
   if (! x_is_initial) {
	 Paso_Copy(n, x, b);
	 Paso_Preconditioner_LocalSmoother_Sweep(A_p->mainBlock,smoother->localSmoother,x);
	 norm_dx=Paso_lsup(n,x,A_p->mpi_info);
	 s++;
   } 
   while (norm_dx > atol) {
	 Paso_Copy(n, b_new, b);
         Paso_SystemMatrix_MatrixVector((-1.), A_p, x, 1., b_new); /* b_new = b - A*x */

	 Paso_Preconditioner_LocalSmoother_Sweep(A_p->mainBlock,smoother->localSmoother,b_new);	
         norm_dx=Paso_lsup(n,b_new,A_p->mpi_info);
	 Paso_AXPY(n, x, 1., b_new); 
	 if (s >= max_sweeps) {
	      errorCode = PRECONDITIONER_MAXITER_REACHED;
	      break;
	 }
	 s++;
   }
   *sweeps=s;
   return errorCode;
}



void Paso_Preconditioner_LocalSmoother_solve(Paso_SparseMatrix* A_p, Paso_Preconditioner_LocalSmoother * smoother, double * x, const double * b, 
					     const dim_t sweeps, const bool_t x_is_initial) 
{
   const dim_t n= (A_p->numRows) * (A_p->row_block_size);
   double *b_new = smoother->buffer;
   dim_t nsweeps=sweeps;
   
   if (! x_is_initial) {
      Paso_Copy(n, x, b);
      Paso_Preconditioner_LocalSmoother_Sweep(A_p,smoother,x);
      nsweeps--;
   }
   
   while (nsweeps > 0 ) {
	 Paso_Copy(n, b_new, b);

	 Paso_SparseMatrix_MatrixVector_CSR_OFFSET0((-1.), A_p, x, 1., b_new); /* b_new = b - A*x */
	 Paso_Preconditioner_LocalSmoother_Sweep(A_p,smoother,b_new);
	 Paso_AXPY(n, x, 1., b_new);
	 nsweeps--;
   }
}

void Paso_Preconditioner_LocalSmoother_Sweep(Paso_SparseMatrix* A, Paso_Preconditioner_LocalSmoother * smoother, double * x) 
{
   const dim_t nt=omp_get_max_threads();
   if (smoother->Jacobi) {
      Paso_BlockOps_solveAll(A->row_block_size,A->numRows,smoother->diag,smoother->pivot,x);
   } else {
      if (nt < 2) {
	 Paso_Preconditioner_LocalSmoother_Sweep_sequential(A,smoother,x);
      } else {
	 Paso_Preconditioner_LocalSmoother_Sweep_colored(A,smoother,x);
      }
   }
}

/* inplace Gauss-Seidel sweep in sequential mode: */

void Paso_Preconditioner_LocalSmoother_Sweep_sequential(Paso_SparseMatrix* A_p, Paso_Preconditioner_LocalSmoother * smoother, double * x)
{
   const dim_t n=A_p->numRows;
   const dim_t n_block=A_p->row_block_size;
   double *diag = smoother->diag;
   index_t* pivot = smoother->pivot; 
   const dim_t block_len=A_p->block_size; 

   
   register dim_t i,k;
   register index_t iptr_ik, mm;
   register double rtmp;
   int failed = 0;
   
   const index_t* ptr_main = Paso_SparseMatrix_borrowMainDiagonalPointer(A_p);
   
   (void)pivot;		/* silence warning from var being unused by macros */
   (void)block_len;    
   /* forward substitution */
   
   if (n_block==1) {
      x[0]*=diag[0];
      for (i = 1; i < n; ++i) {
	 mm=ptr_main[i];
	 /* x_i=x_i-a_ik*x_k  (with k<i) */
	 rtmp=x[i];
	 for (iptr_ik=A_p->pattern->ptr[i];iptr_ik<mm; ++iptr_ik) {
	    k=A_p->pattern->index[iptr_ik]; 
	    rtmp-=A_p->val[iptr_ik]*x[k];
	 }
	 x[i]=rtmp*diag[i];
      }
   } else if (n_block==2) {
      Paso_BlockOps_MViP_2(&diag[0], &x[0]);
      for (i = 1; i < n; ++i) {
	 mm=ptr_main[i];
	 for (iptr_ik=A_p->pattern->ptr[i];iptr_ik<mm; ++iptr_ik) {
	    k=A_p->pattern->index[iptr_ik];                          
	    Paso_BlockOps_SMV_2(&x[2*i], &A_p->val[4*iptr_ik], &x[2*k]);
	 }
	 Paso_BlockOps_MViP_2(&diag[4*i], &x[2*i]);
      }
   } else if (n_block==3) {
      Paso_BlockOps_MViP_3(&diag[0], &x[0]);
      for (i = 1; i < n; ++i) {
	 mm=ptr_main[i];
	 /* x_i=x_i-a_ik*x_k */
	 for (iptr_ik=A_p->pattern->ptr[i];iptr_ik<mm; ++iptr_ik) {
	    k=A_p->pattern->index[iptr_ik];
	    Paso_BlockOps_SMV_3(&x[3*i], &A_p->val[9*iptr_ik], &x[3*k]);
	 }
	 Paso_BlockOps_MViP_3(&diag[9*i], &x[3*i]); 
      }
   } else {
      Paso_BlockOps_solve_N(n_block, &x[0], &diag[0], &pivot[0], &failed);
      for (i = 1; i < n; ++i) {
	 mm=ptr_main[i];
	 /* x_i=x_i-a_ik*x_k */
	 for (iptr_ik=A_p->pattern->ptr[i];iptr_ik<mm; ++iptr_ik) {
	    k=A_p->pattern->index[iptr_ik];
	    Paso_BlockOps_SMV_N(n_block, &x[n_block*i], &A_p->val[block_len*iptr_ik], &x[n_block*k]);
	 }
	 Paso_BlockOps_solve_N(n_block, &x[n_block*i], &diag[block_len*i], &pivot[n_block*i], &failed);
      }
   }
   /* backward substitution */
   
   if (n_block==1) {
      for (i = n-2; i > -1; --i) {	       
	    mm=ptr_main[i];
	    rtmp=x[i]*A_p->val[mm];
	    for (iptr_ik=mm+1; iptr_ik < A_p->pattern->ptr[i+1]; ++iptr_ik) {
	       k=A_p->pattern->index[iptr_ik];
	       rtmp-=A_p->val[iptr_ik]*x[k];
	    }
	    x[i]=diag[i]*rtmp;
      }
      
   } else if (n_block==2) {
      for (i = n-2; i > -1; --i) {
	    mm=ptr_main[i];
	    Paso_BlockOps_MViP_2(&A_p->val[4*mm], &x[2*i]);
	    for (iptr_ik=mm+1; iptr_ik < A_p->pattern->ptr[i+1]; ++iptr_ik) {
	       k=A_p->pattern->index[iptr_ik]; 
	       Paso_BlockOps_SMV_2(&x[2*i], &A_p->val[4*iptr_ik], &x[2*k]);
	    }
	    Paso_BlockOps_MViP_2(&diag[i*4], &x[2*i]);
	    
      }
   } else if (n_block==3) {
      for (i = n-2; i > -1; --i) {
	    mm=ptr_main[i];
	    Paso_BlockOps_MViP_3(&A_p->val[9*mm], &x[3*i]);
	    for (iptr_ik=mm+1; iptr_ik < A_p->pattern->ptr[i+1]; ++iptr_ik) {
	       k=A_p->pattern->index[iptr_ik];    
	       Paso_BlockOps_SMV_3(&x[3*i], &A_p->val[9*iptr_ik], &x[3*k]);
	    }
	    Paso_BlockOps_MViP_3(&diag[i*9], &x[3*i]);
      }
      
   } else {
      double *y=TMPMEMALLOC(n_block, double);
      for (i = n-2; i > -1; --i) {
	 mm=ptr_main[i];
	 Paso_BlockOps_MV_N(n_block, &y[0], &A_p->val[block_len*mm], &x[n_block*i]);
	 for (iptr_ik=mm+1; iptr_ik < A_p->pattern->ptr[i+1]; ++iptr_ik) {
	    k=A_p->pattern->index[iptr_ik];    
	    Paso_BlockOps_SMV_N(n_block, &y[0], &A_p->val[block_len*iptr_ik], &x[n_block*k]);
	 }
	 Paso_BlockOps_Cpy_N(n_block ,&x[n_block*i], &y[0]);
	 Paso_BlockOps_solve_N(n_block, &x[n_block*i], &diag[i*block_len], &pivot[i*n_block], &failed);
      }   
      TMPMEMFREE(y);
   }
   
   if (failed > 0) {
      Esys_setError(ZERO_DIVISION_ERROR, "Paso_Preconditioner_LocalSmoother_Sweep_sequential: non-regular main diagonal block.");
   }
   
   return;
}
       
void Paso_Preconditioner_LocalSmoother_Sweep_colored(Paso_SparseMatrix* A_p, Paso_Preconditioner_LocalSmoother * smoother, double * x) 
{
   const dim_t n=A_p->numRows;
   const dim_t n_block=A_p->row_block_size;
   double *diag = smoother->diag;
   index_t* pivot = smoother->pivot;   
   const dim_t block_len=A_p->block_size;   
   double *y;
   
   register dim_t i,k;
   register index_t color,iptr_ik, mm;
   register double rtmp;
   int failed = 0;
   
   const index_t* coloring = Paso_Pattern_borrowColoringPointer(A_p->pattern);
   const dim_t num_colors = Paso_Pattern_getNumColors(A_p->pattern);
   const index_t* ptr_main = Paso_SparseMatrix_borrowMainDiagonalPointer(A_p);
   
   (void)pivot;			/* These vars are dropped by some macros*/
   (void)block_len;

   #pragma omp parallel  private(mm, i,iptr_ik,k,rtmp, color, y)
   {
      if (n_block>3) {
	 y=TMPMEMALLOC(n_block, double);
      } else {
	 y=NULL;
      }
      /* forward substitution */

      /* color = 0 */
      if (n_block==1) { 
	 #pragma omp  for schedule(static) 
	 for (i = 0; i < n; ++i) {
	    if (coloring[i]== 0 ) x[i]*=diag[i];
	 }
      } else if (n_block==2) {
         #pragma omp for schedule(static)
	 for (i = 0; i < n; ++i) {
	    if (coloring[i]== 0 ) Paso_BlockOps_MViP_2(&diag[i*4], &x[2*i]);
	 }
      } else if (n_block==3) {
	 #pragma omp for schedule(static)
	 for (i = 0; i < n; ++i) {
	    if (coloring[i]==0) Paso_BlockOps_MViP_3(&diag[i*9], &x[3*i]);
	 }
      } else {
	 #pragma omp for schedule(static)
	 for (i = 0; i < n; ++i) {
	    if (coloring[i]==0) Paso_BlockOps_solve_N(n_block, &x[n_block*i], &diag[block_len*i], &pivot[n_block*i], &failed);
	 }
      }
   

      for (color=1;color<num_colors;++color) {
 
	 if (n_block==1) {
	    #pragma omp for schedule(static)
	    for (i = 0; i < n; ++i) {
	       if (coloring[i]==color) {
		  /* x_i=x_i-a_ik*x_k */  
		  rtmp=x[i];
		  for (iptr_ik=A_p->pattern->ptr[i];iptr_ik<A_p->pattern->ptr[i+1]; ++iptr_ik) {
		     k=A_p->pattern->index[iptr_ik];
		     if (coloring[k]<color) rtmp-=A_p->val[iptr_ik]*x[k]; 
		  }
		  x[i]=diag[i]*rtmp;
	       }
	    }
	 } else if (n_block==2) {
	    #pragma omp for schedule(static)
	    for (i = 0; i < n; ++i) {
	       if (coloring[i]==color) {
		  for (iptr_ik=A_p->pattern->ptr[i];iptr_ik<A_p->pattern->ptr[i+1]; ++iptr_ik) {
		     k=A_p->pattern->index[iptr_ik];                          
		     if (coloring[k]<color) Paso_BlockOps_SMV_2(&x[2*i], &A_p->val[4*iptr_ik], &x[2*k]); 
		  }
		  Paso_BlockOps_MViP_2(&diag[4*i], &x[2*i]);
	       }
	    }
	 } else if (n_block==3) {
	    #pragma omp for schedule(static)
	    for (i = 0; i < n; ++i) {
	       if (coloring[i]==color) {
		  for (iptr_ik=A_p->pattern->ptr[i];iptr_ik<A_p->pattern->ptr[i+1]; ++iptr_ik) {
		     k=A_p->pattern->index[iptr_ik];                          
		     if (coloring[k]<color) Paso_BlockOps_SMV_3(&x[3*i], &A_p->val[9*iptr_ik], &x[3*k]);
		  }
		  Paso_BlockOps_MViP_3(&diag[9*i], &x[3*i]);
	       }
	    }
	 } else {
	    #pragma omp for schedule(static)
	    for (i = 0; i < n; ++i) {
	       if (coloring[i] == color) {
		  for (iptr_ik=A_p->pattern->ptr[i];iptr_ik<A_p->pattern->ptr[i+1]; ++iptr_ik) {
		     k=A_p->pattern->index[iptr_ik];                          
		     if (coloring[k]<color) Paso_BlockOps_SMV_N(n_block, &x[n_block*i], &A_p->val[block_len*iptr_ik], &x[n_block*k]);
		  }
		  Paso_BlockOps_solve_N(n_block, &x[n_block*i], &diag[block_len*i], &pivot[n_block*i], &failed);
	       }
	    }
	 }
      } /* end of coloring loop */
      
      /* backward substitution */
      for (color=(num_colors)-2 ;color>-1;--color) { /* Note: color=(num_colors)-1 is not required */
	 if (n_block==1) {
	    #pragma omp for schedule(static)
	    for (i = 0; i < n; ++i) {
	       if (coloring[i]==color) {
		  mm=ptr_main[i];
		  rtmp=A_p->val[mm]*x[i];
		  for (iptr_ik=A_p->pattern->ptr[i];iptr_ik<A_p->pattern->ptr[i+1]; ++iptr_ik) {
		     k=A_p->pattern->index[iptr_ik];                          
		     if (coloring[k]>color) rtmp-=A_p->val[iptr_ik]*x[k];
		  }
		  x[i]= rtmp*diag[i];
	       }
	    }
	 } else if (n_block==2) {
	    #pragma omp for schedule(static)
	    for (i = 0; i < n; ++i) {
	       if (coloring[i]==color) {
		  mm=ptr_main[i];
		  Paso_BlockOps_MViP_2(&A_p->val[4*mm], &x[2*i]);
		  for (iptr_ik=A_p->pattern->ptr[i];iptr_ik<A_p->pattern->ptr[i+1]; ++iptr_ik) {
		     k=A_p->pattern->index[iptr_ik];                          
		     if (coloring[k]>color) Paso_BlockOps_SMV_2(&x[2*i], &A_p->val[4*iptr_ik], &x[2*k]);
		  }
		  Paso_BlockOps_MViP_2(&diag[4*i], &x[2*i]);
	       }
	    }
	 } else if (n_block==3) {
	    #pragma omp for schedule(static)
	    for (i = 0; i < n; ++i) {
	       if (coloring[i]==color) {
		  mm=ptr_main[i];
		  Paso_BlockOps_MViP_3(&A_p->val[9*mm], &x[3*i]);
		  for (iptr_ik=A_p->pattern->ptr[i];iptr_ik<A_p->pattern->ptr[i+1]; ++iptr_ik) {
		     k=A_p->pattern->index[iptr_ik];                          
		     if (coloring[k]>color) Paso_BlockOps_SMV_3(&x[3*i], &A_p->val[9*iptr_ik], &x[3*k]);
		  }
		  Paso_BlockOps_MViP_3(&diag[9*i], &x[3*i]);
	       }
	    }
	 } else {
	    #pragma omp for schedule(static)
	    for (i = 0; i < n; ++i) {
	       if (coloring[i]==color) {
		  mm=ptr_main[i];
		  Paso_BlockOps_MV_N(n_block, &y[0], &A_p->val[block_len*mm], &x[n_block*i]);
		  for (iptr_ik=A_p->pattern->ptr[i];iptr_ik<A_p->pattern->ptr[i+1]; ++iptr_ik) {
		     k=A_p->pattern->index[iptr_ik];                          
		     if (coloring[k]>color) Paso_BlockOps_SMV_N(n_block, &y[0], &A_p->val[block_len*iptr_ik], &x[n_block*k]);
		  }
		  Paso_BlockOps_Cpy_N(n_block ,&x[n_block*i], &y[0]);
		  Paso_BlockOps_solve_N(n_block, &x[n_block*i], &diag[i*block_len], &pivot[i*n_block], &failed);
	       }
	    }
	 }
      }
      TMPMEMFREE(y);
   }
   if (failed > 0) {
      Esys_setError(ZERO_DIVISION_ERROR, "Paso_Preconditioner_LocalSmoother_Sweep_colored: non-regular main diagonal block.");
   }
   return;
}
