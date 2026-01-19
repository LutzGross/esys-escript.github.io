
/*****************************************************************************
*
* Copyright (c) 2003-2018 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
*
*****************************************************************************/


/****************************************************************************/

/* Paso: Gauss-Seidel                                         */

/****************************************************************************/

/* Author: artak@uq.edu.au                                    */

/****************************************************************************/

#include "Preconditioner.h"
#include "BlockOps.h"
#include "PasoUtil.h"

namespace paso {

void Preconditioner_Smoother_free(Preconditioner_Smoother* in)
{
    if (in!=NULL) {
        Preconditioner_LocalSmoother_free(in->localSmoother);
        delete in;
    }
}

void Preconditioner_LocalSmoother_free(Preconditioner_LocalSmoother* in)
{
    if (in!=NULL) {
        delete[] in->diag;
        delete[] in->pivot;
        delete[] in->buffer;
        delete in;
    }
}


/// constructs the symmetric Gauss-Seidel preconditioner
Preconditioner_Smoother* Preconditioner_Smoother_alloc(SystemMatrix_ptr<double> A,
        bool jacobi, bool is_local, bool verbose)
{
    Preconditioner_Smoother* out=new Preconditioner_Smoother;
    out->localSmoother=Preconditioner_LocalSmoother_alloc(A->mainBlock,
                                                jacobi, verbose);
    out->is_local=is_local;
    return out;
}

Preconditioner_LocalSmoother* Preconditioner_LocalSmoother_alloc(
        SparseMatrix_ptr<double> A, bool jacobi, bool verbose)
{
    const dim_t n=A->numRows;
    const dim_t n_block=A->row_block_size;
    const dim_t block_size=A->block_size;
    double time0=escript::gettime();
    Preconditioner_LocalSmoother* out=new Preconditioner_LocalSmoother;

    out->diag=new double[((size_t) n) * ((size_t) block_size)];
    out->pivot=new index_t[ ((size_t) n) * ((size_t)  n_block)];
    out->buffer=new double[((size_t) n) * ((size_t)  n_block)];
    out->Jacobi=jacobi;
    A->invMain(out->diag, out->pivot);
    time0=escript::gettime()-time0;
    return out;
}

/*
performs a few sweeps of the form

S (x_{k} -  x_{k-1}) = b - A x_{k-1}

where x_{0}=0 and S provides some approximation of A.

Under MPI S is built using A->mainBlock only.
If Smoother is local the defect b - A x_{k-1} is calculated using A->mainBlock only.

The MPI-versioned smoother works in the following way:
 (1) initialize x
 (2) calculate residual r_n = b - A * x_n
 (3) store residual r_n in b_new and pass to Preconditioner_LocalSmoother_Sweep() as input
 (4) /delta x_{n+1} is returned (stored in b_new) as the output of Preconditioner_LocalSmoother_Sweep()
 (5) recover x_{n+1} by /delta x_{n+1} + x_n
 (6) repeat steps 2-5 for until sweeps number reduce to 0
*/
void Preconditioner_Smoother_solve(SystemMatrix_ptr<double> A,
        Preconditioner_Smoother* smoother, double* x, const double* b,
        dim_t sweeps, bool x_is_initial)
{
    const dim_t n = A->mainBlock->numRows * A->mainBlock->row_block_size;
    double *b_new = smoother->localSmoother->buffer;
    dim_t nsweeps=sweeps;
    if (smoother->is_local) {
        Preconditioner_LocalSmoother_solve(A->mainBlock,smoother->localSmoother,x,b,sweeps,x_is_initial);
    } else {
        if (! x_is_initial) {
            util::copy(n, x, b);

            Preconditioner_LocalSmoother_Sweep(A->mainBlock,smoother->localSmoother,x);
            nsweeps--;
        }
        while (nsweeps > 0 ) {
            util::copy(n, b_new, b);
            SparseMatrix_MatrixVector_CSR_OFFSET0(-1., A->mainBlock, x, 1., b_new); /* b_new = b - A*x */
            //A->MatrixVector_CSR_OFFSET0(-1., x, 1., b_new); /* b_new = b - A*x */
            Preconditioner_LocalSmoother_Sweep(A->mainBlock,smoother->localSmoother,b_new);
            util::AXPY(n, x, 1., b_new);
            nsweeps--;
        }
    }
}

SolverResult Preconditioner_Smoother_solve_byTolerance(SystemMatrix_ptr<double> A,
            Preconditioner_Smoother* smoother, double* x, const double* b,
            double atol, dim_t* sweeps, bool x_is_initial)
{
   const dim_t n = A->mainBlock->numRows * A->mainBlock->row_block_size;
   double *b_new = smoother->localSmoother->buffer;
   const dim_t max_sweeps=*sweeps;
   dim_t s=0;
   double norm_dx = atol * 2.;
   SolverResult errorCode = NoError;

   if (! x_is_initial) {
        util::copy(n, x, b);
        Preconditioner_LocalSmoother_Sweep(A->mainBlock,smoother->localSmoother,x);
        norm_dx=util::lsup(n,x,A->mpi_info);
        s++;
   }
   while (norm_dx > atol) {
        util::copy(n, b_new, b);
        SparseMatrix_MatrixVector_CSR_OFFSET0(-1., A->mainBlock, x, 1., b_new); /* b_new = b - A*x */
        //A->MatrixVector(-1., x, 1., b_new); /* b_new = b - A*x */
        Preconditioner_LocalSmoother_Sweep(A->mainBlock,smoother->localSmoother,b_new);
        norm_dx=util::lsup(n,b_new,A->mpi_info);
        util::AXPY(n, x, 1., b_new);
        if (s >= max_sweeps) {
            errorCode = MaxIterReached;
            break;
        }
        s++;
   }
   *sweeps=s;
   return errorCode;
}

void Preconditioner_LocalSmoother_solve(SparseMatrix_ptr<double> A,
                                        Preconditioner_LocalSmoother* smoother,
                                        double* x, const double* b,
                                        dim_t sweeps, bool x_is_initial)
{
   const dim_t n = A->numRows * A->row_block_size;
   double *b_new = smoother->buffer;
   dim_t nsweeps=sweeps;

   if (! x_is_initial) {
        util::copy(n, x, b);
        Preconditioner_LocalSmoother_Sweep(A, smoother, x);
        nsweeps--;
   }

   while (nsweeps > 0 ) {
       util::copy(n, b_new, b);

        SparseMatrix_MatrixVector_CSR_OFFSET0((-1.), A, x, 1., b_new); /* b_new = b - A*x */
        Preconditioner_LocalSmoother_Sweep(A, smoother, b_new);
        util::AXPY(n, x, 1., b_new);
        nsweeps--;
   }
}

/*
  Gauss-Seidel sweep to calculate /delta x_{n+1} from matrix A and
  residual r_n. It has two steps: forward substitution and backward
  substitution.
  Forward substitution: (D+L) * /delta x_{n+1/2} = r_n
  Backward substitution: (D+U) * /delta x_{n+1} = r_n - L * /delta x_{n+1/2}
                                                = D * /delta x_{n+1/2}
  where A = D + L + U
        /delta x_{n+1/2} = x_{n+1/2} - x_n
        /delta x_{n+1} = x_{n+1} - x_n

  Input: Matrix A
         residual r_n (in x)
  Output: /delta x_{n+1} (in x)
*/
void Preconditioner_LocalSmoother_Sweep(SparseMatrix_ptr<double> A,
        Preconditioner_LocalSmoother* smoother, double* x)
{
#ifdef _OPENMP
    const dim_t nt=omp_get_max_threads();
#else
    const dim_t nt=1;
#endif
    if (smoother->Jacobi) {
        BlockOps_solveAll(A->row_block_size,A->numRows,smoother->diag,smoother->pivot,x);
    } else {
        if (nt < 2) {
            Preconditioner_LocalSmoother_Sweep_sequential(A,smoother,x);
        } else {
            Preconditioner_LocalSmoother_Sweep_colored(A,smoother,x);
        }
    }
}

/// inplace Gauss-Seidel sweep in sequential mode
void Preconditioner_LocalSmoother_Sweep_sequential(SparseMatrix_ptr<double> A,
        Preconditioner_LocalSmoother* smoother, double* x)
{
    const dim_t n=A->numRows;
    if (n==0)
        return;

    const dim_t n_block=A->row_block_size;
    double *diag = smoother->diag;
    index_t* pivot = smoother->pivot;
    const dim_t block_len=A->block_size;
    dim_t i,k;
    index_t iptr_ik, mm;
    double rtmp;
    int failed = 0;
    const index_t* ptr_main = A->borrowMainDiagonalPointer();

    // silence warning from var being unused by macros
    (void)pivot;
    (void)block_len;

    // forward substitution
    if (n_block==1) {
        x[0]*=diag[0];
        for (i = 1; i < n; ++i) {
            mm=ptr_main[i];
            /* x_i=x_i-a_ik*x_k (with k<i) */
            rtmp=x[i];
            for (iptr_ik=A->pattern->ptr[i];iptr_ik<mm; ++iptr_ik) {
                k=A->pattern->index[iptr_ik];
                rtmp-=A->val[iptr_ik]*x[k];
            }
            x[i]=rtmp*diag[i];
        }
    } else if (n_block==2) {
        BlockOps_MViP_2(&diag[0], &x[0]);
        for (i = 1; i < n; ++i) {
            mm=ptr_main[i];
            for (iptr_ik=A->pattern->ptr[i];iptr_ik<mm; ++iptr_ik) {
                k=A->pattern->index[iptr_ik];
                BlockOps_SMV_2(&x[2*i], &A->val[4*iptr_ik], &x[2*k]);
            }
            BlockOps_MViP_2(&diag[4*i], &x[2*i]);
        }
    } else if (n_block==3) {
        BlockOps_MViP_3(&diag[0], &x[0]);
        for (i = 1; i < n; ++i) {
            mm=ptr_main[i];
            for (iptr_ik=A->pattern->ptr[i];iptr_ik<mm; ++iptr_ik) {
                k=A->pattern->index[iptr_ik];
                BlockOps_SMV_3(&x[3*i], &A->val[9*iptr_ik], &x[3*k]);
            }
            BlockOps_MViP_3(&diag[9*i], &x[3*i]);
        }
    } else {
        BlockOps_solve_N(n_block, &x[0], &diag[0], &pivot[0], &failed);
        for (i = 1; i < n; ++i) {
            mm=ptr_main[i];
            for (iptr_ik=A->pattern->ptr[i];iptr_ik<mm; ++iptr_ik) {
                k=A->pattern->index[iptr_ik];
                BlockOps_SMV_N(n_block, &x[n_block*i], &A->val[block_len*iptr_ik], &x[n_block*k]);
            }
            BlockOps_solve_N(n_block, &x[n_block*i], &diag[block_len*i], &pivot[n_block*i], &failed);
        }
    }

    // backward sweeps
    if (n_block==1) {
        for (i = n-2; i > -1; --i) {
            mm=ptr_main[i];
            rtmp=x[i]*A->val[mm];
            for (iptr_ik=mm+1; iptr_ik < A->pattern->ptr[i+1]; ++iptr_ik) {
                k=A->pattern->index[iptr_ik];
                rtmp-=A->val[iptr_ik]*x[k];
            }
            x[i]=diag[i]*rtmp;
        }
    } else if (n_block==2) {
        for (i = n-2; i > -1; --i) {
            mm=ptr_main[i];
            BlockOps_MViP_2(&A->val[4*mm], &x[2*i]);
            for (iptr_ik=mm+1; iptr_ik < A->pattern->ptr[i+1]; ++iptr_ik) {
                k=A->pattern->index[iptr_ik];
                BlockOps_SMV_2(&x[2*i], &A->val[4*iptr_ik], &x[2*k]);
            }
            BlockOps_MViP_2(&diag[i*4], &x[2*i]);
        }
    } else if (n_block==3) {
        for (i = n-2; i > -1; --i) {
            mm=ptr_main[i];
            BlockOps_MViP_3(&A->val[9*mm], &x[3*i]);
            for (iptr_ik=mm+1; iptr_ik < A->pattern->ptr[i+1]; ++iptr_ik) {
                k=A->pattern->index[iptr_ik];
                BlockOps_SMV_3(&x[3*i], &A->val[9*iptr_ik], &x[3*k]);
            }
            BlockOps_MViP_3(&diag[i*9], &x[3*i]);
        }
    } else {
        double *y=new double[n_block];
        for (i = n-2; i > -1; --i) {
            mm=ptr_main[i];
            BlockOps_MV_N(n_block, &y[0], &A->val[block_len*mm], &x[n_block*i]);
            for (iptr_ik=mm+1; iptr_ik < A->pattern->ptr[i+1]; ++iptr_ik) {
                k=A->pattern->index[iptr_ik];
                BlockOps_SMV_N(n_block, &y[0], &A->val[block_len*iptr_ik], &x[n_block*k]);
            }
            BlockOps_Cpy_N(n_block ,&x[n_block*i], &y[0]);
            BlockOps_solve_N(n_block, &x[n_block*i], &diag[i*block_len], &pivot[i*n_block], &failed);
        }
        delete[] y;
    }

    if (failed > 0) {
        throw PasoException("Preconditioner_LocalSmoother_Sweep_sequential: non-regular main diagonal block.");
    }
}

void Preconditioner_LocalSmoother_Sweep_colored(SparseMatrix_ptr<double> A,
        Preconditioner_LocalSmoother* smoother, double* x)
{
    const dim_t n=A->numRows;
    const dim_t n_block=A->row_block_size;
    double *diag = smoother->diag;
    index_t* pivot = smoother->pivot;
    const dim_t block_len=A->block_size;
    double *y;

    dim_t i,k;
    index_t color,iptr_ik, mm;
    double rtmp;
    int failed = 0;

    const index_t* coloring = A->pattern->borrowColoringPointer();
    const dim_t num_colors = A->pattern->getNumColors();
    const index_t* ptr_main = A->borrowMainDiagonalPointer();

    (void)pivot;                 /* These vars are dropped by some macros*/
    (void)block_len;

    #pragma omp parallel  private(mm, i,iptr_ik,k,rtmp, color, y)
    {
        if (n_block>3) {
            y=new double[n_block];
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
                if (coloring[i]== 0 ) BlockOps_MViP_2(&diag[i*4], &x[2*i]);
            }
        } else if (n_block==3) {
            #pragma omp for schedule(static)
            for (i = 0; i < n; ++i) {
                if (coloring[i]==0) BlockOps_MViP_3(&diag[i*9], &x[3*i]);
            }
        } else {
            #pragma omp for schedule(static)
            for (i = 0; i < n; ++i) {
                if (coloring[i]==0) BlockOps_solve_N(n_block, &x[n_block*i], &diag[block_len*i], &pivot[n_block*i], &failed);
            }
        }

        for (color=1;color<num_colors;++color) {
            if (n_block==1) {
                #pragma omp for schedule(static)
                for (i = 0; i < n; ++i) {
                    if (coloring[i]==color) {
                        /* x_i=x_i-a_ik*x_k */
                        rtmp=x[i];
                        for (iptr_ik=A->pattern->ptr[i];iptr_ik<A->pattern->ptr[i+1]; ++iptr_ik) {
                            k=A->pattern->index[iptr_ik];
                            if (coloring[k]<color) rtmp-=A->val[iptr_ik]*x[k];
                        }
                        x[i]=diag[i]*rtmp;
                    }
                }
            } else if (n_block==2) {
                #pragma omp for schedule(static)
                for (i = 0; i < n; ++i) {
                    if (coloring[i]==color) {
                        for (iptr_ik=A->pattern->ptr[i];iptr_ik<A->pattern->ptr[i+1]; ++iptr_ik) {
                            k=A->pattern->index[iptr_ik];
                            if (coloring[k]<color) BlockOps_SMV_2(&x[2*i], &A->val[4*iptr_ik], &x[2*k]);
                        }
                        BlockOps_MViP_2(&diag[4*i], &x[2*i]);
                    }
                }
            } else if (n_block==3) {
                #pragma omp for schedule(static)
                for (i = 0; i < n; ++i) {
                    if (coloring[i]==color) {
                        for (iptr_ik=A->pattern->ptr[i];iptr_ik<A->pattern->ptr[i+1]; ++iptr_ik) {
                            k=A->pattern->index[iptr_ik];
                            if (coloring[k]<color) BlockOps_SMV_3(&x[3*i], &A->val[9*iptr_ik], &x[3*k]);
                        }
                        BlockOps_MViP_3(&diag[9*i], &x[3*i]);
                    }
                }
            } else {
                #pragma omp for schedule(static)
                for (i = 0; i < n; ++i) {
                    if (coloring[i] == color) {
                        for (iptr_ik=A->pattern->ptr[i];iptr_ik<A->pattern->ptr[i+1]; ++iptr_ik) {
                            k=A->pattern->index[iptr_ik];
                            if (coloring[k]<color) BlockOps_SMV_N(n_block, &x[n_block*i], &A->val[block_len*iptr_ik], &x[n_block*k]);
                        }
                        BlockOps_solve_N(n_block, &x[n_block*i], &diag[block_len*i], &pivot[n_block*i], &failed);
                    }
                }
            }
        } // end of coloring loop

        // backward substitution

        // Note: color=(num_colors)-1 is not required
        for (color=(num_colors)-2 ;color>-1;--color) {
            if (n_block==1) {
                #pragma omp for schedule(static)
                for (i = 0; i < n; ++i) {
                    if (coloring[i]==color) {
                        mm=ptr_main[i];
                        rtmp=A->val[mm]*x[i];
                        for (iptr_ik=A->pattern->ptr[i];iptr_ik<A->pattern->ptr[i+1]; ++iptr_ik) {
                            k=A->pattern->index[iptr_ik];
                            if (coloring[k]>color) rtmp-=A->val[iptr_ik]*x[k];
                        }
                        x[i]= rtmp*diag[i];
                    }
                }
            } else if (n_block==2) {
                #pragma omp for schedule(static)
                for (i = 0; i < n; ++i) {
                    if (coloring[i]==color) {
                        mm=ptr_main[i];
                        BlockOps_MViP_2(&A->val[4*mm], &x[2*i]);
                        for (iptr_ik=A->pattern->ptr[i];iptr_ik<A->pattern->ptr[i+1]; ++iptr_ik) {
                            k=A->pattern->index[iptr_ik];
                            if (coloring[k]>color) BlockOps_SMV_2(&x[2*i], &A->val[4*iptr_ik], &x[2*k]);
                        }
                        BlockOps_MViP_2(&diag[4*i], &x[2*i]);
                    }
                }
            } else if (n_block==3) {
                #pragma omp for schedule(static)
                for (i = 0; i < n; ++i) {
                    if (coloring[i]==color) {
                        mm=ptr_main[i];
                        BlockOps_MViP_3(&A->val[9*mm], &x[3*i]);
                        for (iptr_ik=A->pattern->ptr[i];iptr_ik<A->pattern->ptr[i+1]; ++iptr_ik) {
                            k=A->pattern->index[iptr_ik];
                            if (coloring[k]>color) BlockOps_SMV_3(&x[3*i], &A->val[9*iptr_ik], &x[3*k]);
                        }
                        BlockOps_MViP_3(&diag[9*i], &x[3*i]);
                    }
                }
            } else {
                #pragma omp for schedule(static)
                for (i = 0; i < n; ++i) {
                    if (coloring[i]==color) {
                        mm=ptr_main[i];
                        BlockOps_MV_N(n_block, &y[0], &A->val[block_len*mm], &x[n_block*i]);
                        for (iptr_ik=A->pattern->ptr[i];iptr_ik<A->pattern->ptr[i+1]; ++iptr_ik) {
                            k=A->pattern->index[iptr_ik];
                            if (coloring[k]>color) BlockOps_SMV_N(n_block, &y[0], &A->val[block_len*iptr_ik], &x[n_block*k]);
                        }
                        BlockOps_Cpy_N(n_block ,&x[n_block*i], &y[0]);
                        BlockOps_solve_N(n_block, &x[n_block*i], &diag[i*block_len], &pivot[i*n_block], &failed);
                    }
                }
            }
        }
        delete[] y;
    }
    if (failed > 0) {
        throw PasoException("Preconditioner_LocalSmoother_Sweep_colored: non-regular main diagonal block.");
    }
}

} // namespace paso

