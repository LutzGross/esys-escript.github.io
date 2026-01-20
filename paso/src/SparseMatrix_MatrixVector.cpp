
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* https://github.com/LutzGross/esys-escript.github.io
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
*
*****************************************************************************/


/****************************************************************************
 *
 * Paso: raw scaled vector update operation:
 *                  out = alpha * A * in + beta * out
 ****************************************************************************/

/* Author: Lutz Gross, l.gross@uq.edu.au */

/****************************************************************************/

#include "SparseMatrix.h"

namespace paso {

// forward declaration
void SparseMatrix_MatrixVector_CSR_OFFSET0_stripe(double alpha, dim_t nRows,
                                                  dim_t row_block_size,
                                                  dim_t col_block_size,
                                                  const index_t* ptr,
                                                  const index_t* index,
                                                  const double* val,
                                                  const double* in,
                                                  double beta, double* out);

/* CSC format with offset 0 */
void SparseMatrix_MatrixVector_CSC_OFFSET0(double alpha,
                                           const_SparseMatrix_ptr<double> A,
                                           const double* in,
                                           double beta, double* out)
{
    const int totalRowSize = A->numRows * A->row_block_size;
    if (std::abs(beta) > 0) {
        if (beta != 1.) {
#pragma omp parallel for schedule(static)
            for (index_t irow=0; irow < totalRowSize; irow++) {
                out[irow] *= beta;
            }
        }
    } else {
#pragma omp parallel for schedule(static)
        for (index_t irow=0; irow < totalRowSize; irow++) {
            out[irow] = 0;
        }
    }

    if (A->pattern->isEmpty())
        return;

    if (std::abs(alpha) > 0) {
        if (A->col_block_size==1 && A->row_block_size==1) {
            /* TODO: parallelize (good luck!) */
            for (index_t icol=0; icol < A->pattern->numOutput; ++icol) {
                #pragma ivdep
                for (index_t iptr=A->pattern->ptr[icol];
                        iptr < A->pattern->ptr[icol+1]; ++iptr) {
                    out[A->pattern->index[iptr]] += alpha * A->val[iptr] * in[icol];
                }
            }
        } else if (A->col_block_size==2 && A->row_block_size==2) {
            /* TODO: parallelize */
            for (index_t ic=0; ic < A->pattern->numOutput; ic++) {
                #pragma ivdep
                for (index_t iptr=A->pattern->ptr[ic];
                        iptr < A->pattern->ptr[ic+1]; iptr++) {
                    const index_t ir = 2*A->pattern->index[iptr];
                    out[  ir] += alpha * (A->val[iptr*4  ]*in[ic] + A->val[iptr*4+2]*in[1+ic]);
                    out[1+ir] += alpha * (A->val[iptr*4+1]*in[ic] + A->val[iptr*4+3]*in[1+ic]);
                }
            }
        } else if (A->col_block_size==3 && A->row_block_size==3) {
            /* TODO: parallelize */
            for (index_t ic=0; ic < A->pattern->numOutput; ic++) {
                #pragma ivdep
                for (index_t iptr=A->pattern->ptr[ic];
                        iptr < A->pattern->ptr[ic+1]; iptr++) {
                    const index_t ir = 3*A->pattern->index[iptr];
                    out[  ir] += alpha * (A->val[iptr*9  ]*in[ic] + A->val[iptr*9+3]*in[1+ic] + A->val[iptr*9+6]*in[2+ic]);
                    out[1+ir] += alpha * (A->val[iptr*9+1]*in[ic] + A->val[iptr*9+4]*in[1+ic] + A->val[iptr*9+7]*in[2+ic]);
                    out[2+ir] += alpha * (A->val[iptr*9+2]*in[ic] + A->val[iptr*9+5]*in[1+ic] + A->val[iptr*9+8]*in[2+ic]);
                }
            }
        } else {
            /* TODO: parallelize */
            for (index_t ic=0; ic < A->pattern->numOutput; ic++) {
                for (index_t iptr=A->pattern->ptr[ic];
                        iptr < A->pattern->ptr[ic+1]; iptr++) {
                    for (index_t irb=0; irb < A->row_block_size; irb++) {
                        const index_t irow = irb+A->row_block_size*A->pattern->index[iptr];
                        #pragma ivdep
                        for (index_t icb=0; icb < A->col_block_size; icb++) {
                            const index_t icol=icb+A->col_block_size*ic;
                            out[irow] += alpha * A->val[iptr*A->block_size+irb+A->row_block_size*icb] * in[icol];
                        }
                    }
                }
            }
        } // blocksizes
    } // alpha > 0
}

/* CSC format with offset 1 */
void SparseMatrix_MatrixVector_CSC_OFFSET1(double alpha,
                                           const_SparseMatrix_ptr<double> A,
                                           const double* in,
                                           double beta, double* out)
{
    const int totalRowSize = A->numRows * A->row_block_size;
    if (std::abs(beta) > 0) {
        if (beta != 1.) {
#pragma omp parallel for schedule(static)
            for (index_t irow=0; irow < totalRowSize; irow++) {
                out[irow] *= beta;
            }
        }
    } else {
#pragma omp parallel for schedule(static)
        for (index_t irow=0; irow < totalRowSize; irow++) {
            out[irow] = 0;
        }
    }

    if (std::abs(alpha) > 0) {
        if (A->col_block_size==1 && A->row_block_size==1) {
            /* TODO: parallelize (good luck!) */
            for (index_t icol=0; icol < A->pattern->numOutput; icol++) {
                #pragma ivdep
                for (index_t iptr=A->pattern->ptr[icol]-1;
                        iptr<A->pattern->ptr[icol+1]-1; iptr++) {
                    out[A->pattern->index[iptr]-1]+= alpha * A->val[iptr] * in[icol];
                }
            }
        } else if (A->col_block_size==2 && A->row_block_size==2) {
            /* TODO: parallelize */
            for (index_t ic=0; ic < A->pattern->numOutput; ic++) {
                for (index_t iptr=A->pattern->ptr[ic]-1;
                        iptr < A->pattern->ptr[ic+1]-1; iptr++) {
                    const index_t ir=2*(A->pattern->index[iptr]-1);
                    out[  ir] += alpha * (A->val[iptr*4  ]*in[ic] + A->val[iptr*4+2]*in[1+ic]);
                    out[1+ir] += alpha * (A->val[iptr*4+1]*in[ic] + A->val[iptr*4+3]*in[1+ic]);
                }
            }
        } else if (A->col_block_size==3 && A->row_block_size==3) {
            /* TODO: parallelize */
            for (index_t ic=0; ic < A->pattern->numOutput; ic++) {
                #pragma ivdep
                for (index_t iptr=A->pattern->ptr[ic]-1;
                        iptr < A->pattern->ptr[ic+1]-1; iptr++) {
                    const index_t ir=3*(A->pattern->index[iptr]-1);
                    out[  ir] += alpha * (A->val[iptr*9  ]*in[ic] + A->val[iptr*9+3]*in[1+ic] + A->val[iptr*9+6]*in[2+ic]);
                    out[1+ir] += alpha * (A->val[iptr*9+1]*in[ic] + A->val[iptr*9+4]*in[1+ic] + A->val[iptr*9+7]*in[2+ic]);
                    out[2+ir] += alpha * (A->val[iptr*9+2]*in[ic] + A->val[iptr*9+5]*in[1+ic] + A->val[iptr*9+8]*in[2+ic]);
                }
            }
        } else {
            /* TODO: parallelize */
            for (index_t ic=0; ic < A->pattern->numOutput; ic++) {
                for (index_t iptr=A->pattern->ptr[ic]-1;
                        iptr < A->pattern->ptr[ic+1]-1; iptr++) {
                    for (index_t irb=0; irb < A->row_block_size; irb++) {
                        const index_t irow=irb+A->row_block_size*(A->pattern->index[iptr]-1);
                        #pragma ivdep
                        for (index_t icb=0; icb < A->col_block_size; icb++) {
                            const index_t icol=icb+A->col_block_size*ic;
                            out[irow] += alpha * A->val[iptr*A->block_size+irb+A->row_block_size*icb] * in[icol];
                        }
                    }
                }
            }
        } // blocksizes
    } // alpha > 0
}

/* CSR format with offset 0 */
void SparseMatrix_MatrixVector_CSR_OFFSET0(double alpha,
                                           const_SparseMatrix_ptr<double> A,
                                           const double* in,
                                           double beta, double* out)
{
//#define PASO_DYNAMIC_SCHEDULING_MVM
#if defined PASO_DYNAMIC_SCHEDULING_MVM && defined _OPENMP
#define USE_DYNAMIC_SCHEDULING
#endif

    const dim_t nrow = A->numRows;
#ifdef _OPENMP
    const dim_t np = omp_get_max_threads();
#else
    const dim_t np = 1;
#endif
    const dim_t len = nrow/np;

#ifdef USE_DYNAMIC_SCHEDULING
    dim_t chunk_size=1;
    char* chksz_chr=getenv("PASO_CHUNK_SIZE_MVM");
    if (chksz_chr!=NULL) sscanf(chksz_chr, "%d",&chunk_size);
    chunk_size=std::min(std::max(1,chunk_size),len);
    dim_t n_chunks=nrow/chunk_size;
    if (n_chunks*chunk_size<nrow) n_chunks+=1;

#pragma omp parallel for schedule(dynamic,1)
    for (dim_t p=0; p < n_chunks; p++) {
        const dim_t irow=chunk_size*p;
        const dim_t local_n=std::min(chunk_size,nrow-chunk_size*p);
        SparseMatrix_MatrixVector_CSR_OFFSET0_stripe(alpha, local_n,
            A->row_block_size, A->col_block_size, &(A->pattern->ptr[irow]),
            A->pattern->index, A->val, in, beta,
            &out[irow*A->row_block_size]);
    }

#else // static scheduling
    const dim_t rest=nrow-len*np;

#pragma omp parallel for
    for (dim_t p=0; p < np; p++) {
        const dim_t irow=len*p+std::min(p,rest);
        const dim_t local_n=len+(p<rest ? 1 :0 );
        SparseMatrix_MatrixVector_CSR_OFFSET0_stripe(alpha, local_n,
            A->row_block_size, A->col_block_size, &(A->pattern->ptr[irow]),
            A->pattern->index, A->val, in, beta,
            &out[irow*A->row_block_size]);
    }
#endif // scheduling
}

/* CSR format with offset 0 */
void SparseMatrix_MatrixVector_CSR_OFFSET0_stripe(double alpha, dim_t nRows,
        dim_t row_block_size, dim_t col_block_size, const index_t* ptr,
        const index_t* index, const double* val, const double* in,
        double beta, double* out)
{
    if (std::abs(beta) > 0) {
        if (beta != 1.) {
            for (index_t irow=0; irow < nRows*row_block_size; irow++)
                out[irow] *= beta;
        }
    } else {
        for (index_t irow=0; irow < nRows*row_block_size; irow++)
            out[irow] = 0;
    }

    if (std::abs(alpha) > 0) {
        if (col_block_size==1 && row_block_size ==1) {
            for (index_t irow=0; irow < nRows; ++irow) {
                double reg=0.;
                #pragma ivdep
                for (index_t iptr=ptr[irow]; iptr<ptr[irow+1]; ++iptr) {
                    reg += val[iptr] * in[index[iptr]];
                }
                out[irow] += alpha * reg;
            }
        } else if (col_block_size==2 && row_block_size==2) {
            for (index_t ir=0; ir < nRows; ir++) {
                double reg1=0.;
                double reg2=0.;
                #pragma ivdep
                for (index_t iptr=ptr[ir]; iptr<ptr[ir+1]; iptr++) {
                    const index_t ic=2*index[iptr];
                    const index_t Aiptr=iptr*4;
                    const double in1=in[ic];
                    const double in2=in[1+ic];
                    const double A00=val[Aiptr  ];
                    const double A10=val[Aiptr+1];
                    const double A01=val[Aiptr+2];
                    const double A11=val[Aiptr+3];
                    reg1 += A00*in1 + A01*in2;
                    reg2 += A10*in1 + A11*in2;
                }
                out[  2*ir] += alpha * reg1;
                out[1+2*ir] += alpha * reg2;
            }
        } else if (col_block_size==3 && row_block_size==3) {
            for (index_t ir=0; ir < nRows; ir++) {
                double reg1=0.;
                double reg2=0.;
                double reg3=0.;
                #pragma ivdep
                for (index_t iptr=ptr[ir]; iptr < ptr[ir+1]; iptr++) {
                    const index_t ic=3*index[iptr];
                    const index_t Aiptr=iptr*9;
                    const double in1=in[ic];
                    const double in2=in[1+ic];
                    const double in3=in[2+ic];
                    const double A00=val[Aiptr  ];
                    const double A10=val[Aiptr+1];
                    const double A20=val[Aiptr+2];
                    const double A01=val[Aiptr+3];
                    const double A11=val[Aiptr+4];
                    const double A21=val[Aiptr+5];
                    const double A02=val[Aiptr+6];
                    const double A12=val[Aiptr+7];
                    const double A22=val[Aiptr+8];
                    reg1 += A00*in1 + A01*in2 + A02*in3;
                    reg2 += A10*in1 + A11*in2 + A12*in3;
                    reg3 += A20*in1 + A21*in2 + A22*in3;
                }
                out[  3*ir] += alpha * reg1;
                out[1+3*ir] += alpha * reg2;
                out[2+3*ir] += alpha * reg3;
            }
        } else { // blocksizes > 3
            const dim_t block_size=col_block_size*row_block_size;
            for (index_t ir=0; ir < nRows; ir++) {
                for (index_t iptr=ptr[ir]; iptr < ptr[ir+1]; iptr++) {
                    for (index_t irb=0; irb < row_block_size; irb++) {
                        double reg=0.;
                        #pragma ivdep
                        for (index_t icb=0; icb < col_block_size; icb++) {
                            const index_t icol=icb+col_block_size*index[iptr];
                            reg += val[iptr*block_size+irb+row_block_size*icb] * in[icol];
                        }
                        const index_t irow=irb+row_block_size*ir;
                        out[irow] += alpha * reg;
                    }
                }
            }
        }
    }
}

/* CSR format with offset 0 (diagonal only) */
void SparseMatrix_MatrixVector_CSR_OFFSET0_DIAG(double alpha,
                                                const_SparseMatrix_ptr<double> A,
                                                const double* in,
                                                double beta, double* out)
{
    const int totalRowSize = A->numRows * A->row_block_size;
    if (std::abs(beta) > 0) {
        if (beta != 1.) {
#pragma omp parallel for schedule(static)
            for (index_t irow=0; irow < totalRowSize; irow++)
                out[irow] *= beta;
        }
    } else {
#pragma omp parallel for schedule(static)
        for (index_t irow=0; irow < totalRowSize; irow++)
            out[irow] = 0;
    }

    if (std::abs(alpha) > 0) {
        const int nRows = A->pattern->numOutput;
        if (A->block_size == 1) {
#pragma omp parallel for schedule(static)
            for (index_t irow=0; irow < nRows; irow++) {
                double reg=0.;
                #pragma ivdep
                for (index_t iptr=A->pattern->ptr[irow];
                        iptr < A->pattern->ptr[irow+1]; ++iptr) {
                    reg += A->val[iptr] * in[A->pattern->index[iptr]];
                }
                out[irow] += alpha * reg;
            }
        } else if (A->block_size == 2) {
#pragma omp parallel for schedule(static)
            for (index_t ir=0; ir < nRows; ir++) {
                double reg1=0.;
                double reg2=0.;
                #pragma ivdep
                for (index_t iptr=A->pattern->ptr[ir];
                        iptr < A->pattern->ptr[ir+1]; iptr++) {
                    const index_t ic = 2*A->pattern->index[iptr];
                    reg1 += A->val[iptr*2  ]*in[  ic];
                    reg2 += A->val[iptr*2+1]*in[1+ic];
                }
                out[  2*ir] += alpha * reg1;
                out[1+2*ir] += alpha * reg2;
            }
        } else if (A->block_size == 3) {
#pragma omp parallel for schedule(static)
            for (index_t ir=0; ir < nRows; ir++) {
                double reg1=0.;
                double reg2=0.;
                double reg3=0.;
                #pragma ivdep
                for (index_t iptr=A->pattern->ptr[ir];
                        iptr < A->pattern->ptr[ir+1]; iptr++) {
                    const index_t ic = 3*A->pattern->index[iptr];
                    reg1 += A->val[iptr*3  ]*in[  ic];
                    reg2 += A->val[iptr*3+1]*in[1+ic];
                    reg3 += A->val[iptr*3+2]*in[2+ic];
                }
                out[  3*ir] += alpha * reg1;
                out[1+3*ir] += alpha * reg2;
                out[2+3*ir] += alpha * reg3;
            }
        } else if (A->block_size == 4) {
#pragma omp parallel for schedule(static)
            for (index_t ir=0; ir < nRows; ir++) {
                double reg1=0.;
                double reg2=0.;
                double reg3=0.;
                double reg4=0.;
                #pragma ivdep
                for (index_t iptr=A->pattern->ptr[ir];
                        iptr < A->pattern->ptr[ir+1]; iptr++) {
                    const index_t ic = 4 * A->pattern->index[iptr];
                    reg1 += A->val[iptr*4  ]*in[  ic];
                    reg2 += A->val[iptr*4+1]*in[1+ic];
                    reg3 += A->val[iptr*4+2]*in[2+ic];
                    reg4 += A->val[iptr*4+3]*in[3+ic];
                }
                out[  4*ir] += alpha * reg1;
                out[1+4*ir] += alpha * reg2;
                out[2+4*ir] += alpha * reg3;
                out[3+4*ir] += alpha * reg4;
            }
        } else { // blocksize > 4
#pragma omp parallel for schedule(static)
            for (index_t ir=0; ir < nRows; ir++) {
                for (index_t iptr=A->pattern->ptr[ir];
                        iptr < A->pattern->ptr[ir+1]; iptr++) {
                    for (index_t ib=0; ib < A->block_size; ib++) {
                        const index_t irow = ib+A->row_block_size*ir;
                        const index_t icol = ib+A->col_block_size*A->pattern->index[iptr];
                        out[irow] += alpha * A->val[iptr*A->block_size+ib] * in[icol];
                    }
                }
            }
        } // blocksize
    } // alpha > 0
}

} // namespace paso

