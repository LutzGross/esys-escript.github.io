
/*****************************************************************************
*
* Copyright (c) 2003-2014 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/


/****************************************************************************/

/*   Paso: SystemMatrix */

/****************************************************************************/

/*   Copyrights by ACcESS Australia 2003,2004,2005,2006 */
/*   Author: Lutz Gross, l.gross@uq.edu.au */

/****************************************************************************/

#ifndef __PASO_SYSTEMMATRIX_H__
#define __PASO_SYSTEMMATRIX_H__

#include "Common.h"
#include "Coupler.h"
#include "SparseMatrix.h"
#include "SystemMatrixPattern.h"
#include "Options.h"
#include "esysUtils/Esys_MPI.h"

namespace paso {

struct SystemMatrix;
typedef boost::shared_ptr<SystemMatrix> SystemMatrix_ptr;
typedef boost::shared_ptr<const SystemMatrix> const_SystemMatrix_ptr;

typedef int SystemMatrixType;

//  this struct holds a (distributed) stiffness matrix
PASO_DLL_API
struct SystemMatrix : boost::enable_shared_from_this<SystemMatrix>
{
    SystemMatrix(SystemMatrixType, SystemMatrixPattern_ptr, dim_t, dim_t,
                 bool patternIsUnrolled);

    ~SystemMatrix();

    /// Nullifies rows and columns in the matrix.
    /// The rows and columns are marked by positive values in mask_row and
    /// mask_col. Values on the main diagonal which are marked to set to
    /// zero by both mask_row and mask_col are set to main_diagonal_value.
    void nullifyRowsAndCols(double* mask_row, double* mask_col,
                            double main_diagonal_value);

    /// Nullifies rows in the matrix.
    /// The rows are marked by positive values in mask_row. Values on the
    /// main diagonal which are marked to set to zero by mask_row are set
    /// to main_diagonal_value.
    void nullifyRows(double* mask_row, double main_diagonal_value);

    void add(dim_t, index_t*, dim_t, dim_t, index_t*, dim_t, double*);

    void makeZeroRowSums(double* left_over); 

    /// copies the col_coupleBlock into row_coupleBlock. 
    /// WARNING: this method uses mpi_requests of the coupler attached to the
    /// matrix. No reordering on the received columns is performed.
    /// In practice this means that components in
    /// row_coupleBlock->pattern->index  and
    /// row_coupler->connector->recv->shared
    /// are ordered by increasing value.
    /// Note that send and receive row_coupler->connectors are swapping roles.
    void copyColCoupleBlock();

    void copyRemoteCoupleBlock(bool recreatePattern);

    void fillWithGlobalCoordinates(double f1);

    void print() const;

    /// Merges the system matrix which is distributed on several MPI ranks
    /// into a complete sparse matrix on rank 0. Used by the Merged Solver.
    SparseMatrix_ptr mergeSystemMatrix() const;

    void mergeMainAndCouple(index_t** p_ptr, index_t** p_idx, double** p_val) const;

    void mergeMainAndCouple_CSR_OFFSET0(index_t** p_ptr, index_t** p_idx, double** p_val) const;
    void mergeMainAndCouple_CSR_OFFSET0_Block(index_t** p_ptr, index_t** p_idx, double** p_val) const;

    void mergeMainAndCouple_CSC_OFFSET1(index_t** p_ptr, index_t** p_idx, double** p_val) const;

    void copyMain_CSC_OFFSET1(index_t** p_ptr, index_t** p_idx, double** p_val);

    void extendedRowsForST(dim_t* degree_ST, index_t* offset_ST, index_t* ST);

    void applyBalanceInPlace(double* x, bool RHS) const;

    void applyBalance(double* x_out, const double* x, bool RHS) const;

    void balance();

    double getGlobalSize() const;

    void setPreconditioner(Options* options);

    /// Applies the preconditioner.
    /// This method needs to be called within a parallel region. 
    /// Barrier synchronization is performed before the evaluation to make
    /// sure that the input vector is available
    void solvePreconditioner(double* x, double* b);

    void freePreconditioner();

    index_t* borrowMainDiagonalPointer() const;

    inline void startCollect(const double* in)
    {
        startColCollect(in);
    }

    inline double* finishCollect()
    {
        return finishColCollect();
    }

    inline void startColCollect(const double* in)
    {
        col_coupler->startCollect(in);
    }

    inline double* finishColCollect()
    {
        return col_coupler->finishCollect();
    }

    inline void startRowCollect(const double* in)
    {
        row_coupler->startCollect(in);
    }

    inline double* finishRowCollect()
    {
        return row_coupler->finishCollect();
    }

    inline dim_t getNumRows() const
    {
        return mainBlock->numRows;
    }

    inline dim_t getNumCols() const
    {
        return mainBlock->numCols;
    }

    inline dim_t getTotalNumRows() const
    {
        return getNumRows() * row_block_size;
    }

    inline dim_t getTotalNumCols() const
    {
        return getNumCols() * col_block_size;
    }

    inline dim_t getRowOverlap()  const
    {
        return row_coupler->getNumOverlapComponents();
    }

    inline dim_t getColOverlap() const
    {
        return col_coupler->getNumOverlapComponents();
    }

    inline dim_t getGlobalNumRows() const
    {
        if (type & MATRIX_FORMAT_CSC) {
            return pattern->input_distribution->getGlobalNumComponents();
        }
        return pattern->output_distribution->getGlobalNumComponents();
    }

    inline dim_t getGlobalNumCols() const
    {
        if (type & MATRIX_FORMAT_CSC) {
            return pattern->output_distribution->getGlobalNumComponents();
        }
        return pattern->input_distribution->getGlobalNumComponents();
    }

    inline dim_t getGlobalTotalNumRows() const
    {
        return getGlobalNumRows() * row_block_size;
    }

    inline dim_t getGlobalTotalNumCols() const
    {
        return getGlobalNumCols() * col_block_size;
    }

    inline double getSparsity() const
    {
        return getGlobalSize() /
                 (DBLE(getGlobalTotalNumRows())*getGlobalTotalNumCols());
    }

    inline dim_t getNumOutput() const
    {
       return pattern->getNumOutput();
    }

    inline void copyBlockFromMainDiagonal(double* out) const
    {
        mainBlock->copyBlockFromMainDiagonal(out);
    }

    inline void copyBlockToMainDiagonal(const double* in)
    {
        mainBlock->copyBlockToMainDiagonal(in);
    }

    inline void copyFromMainDiagonal(double* out) const
    {
        mainBlock->copyFromMainDiagonal(out);
    }

    inline void copyToMainDiagonal(const double* in)
    {
        mainBlock->copyToMainDiagonal(in);
    }

    inline void setValues(double value)
    {
        mainBlock->setValues(value);
        col_coupleBlock->setValues(value);
        row_coupleBlock->setValues(value);
        is_balanced = false;
    }

    inline void saveMM(const char* filename) const
    {
        if (mpi_info->size > 1) {
            Esys_setError(IO_ERROR, "SystemMatrix::saveMM: Only single rank supported.");
        } else {
            mainBlock->saveMM(filename);
        }
    }

    inline void saveHB(const char *filename) const
    {
        if (mpi_info->size > 1) {
            Esys_setError(TYPE_ERROR, "SystemMatrix::saveHB: Only single rank supported.");
        } else if (!(type & MATRIX_FORMAT_CSC)) {
            Esys_setError(TYPE_ERROR, "SystemMatrix::saveHB: Only CSC format supported.");
        } else {
            mainBlock->saveHB_CSC(filename);
        }
    }

    inline void rowSum(double* row_sum) const
    {
        if ((type & MATRIX_FORMAT_CSC) || (type & MATRIX_FORMAT_OFFSET1)) {
            Esys_setError(TYPE_ERROR, "SystemMatrix::rowSum: No normalization "
                  "available for compressed sparse column or index offset 1.");
        } else {
            const dim_t nrow = mainBlock->numRows*row_block_size;
#pragma omp parallel for
            for (index_t irow=0; irow<nrow; ++irow) {
                row_sum[irow]=0.;
            }
            mainBlock->addRow_CSR_OFFSET0(row_sum);
            col_coupleBlock->addRow_CSR_OFFSET0(row_sum);
        }
    }

    static SystemMatrix_ptr loadMM_toCSR(const char* filename);

    static SystemMatrix_ptr loadMM_toCSC(const char* filename);

    static index_t getSystemMatrixTypeId(index_t solver,
                                         index_t preconditioner,
                                         index_t package, bool symmetry,
                                         Esys_MPIInfo* mpi_info);

    SystemMatrixType type;
    SystemMatrixPattern_ptr pattern;

    dim_t logical_row_block_size;
    dim_t logical_col_block_size;

    dim_t row_block_size;
    dim_t col_block_size;
    dim_t block_size;

    Distribution_ptr row_distribution;
    Distribution_ptr col_distribution;
    Esys_MPIInfo *mpi_info;

    Coupler_ptr col_coupler;
    Coupler_ptr row_coupler;

    /// main block
    SparseMatrix_ptr mainBlock;
    /// coupling to neighbouring processors (row - col)
    SparseMatrix_ptr col_coupleBlock;
    /// coupling to neighbouring processors (col - row)
    SparseMatrix_ptr row_coupleBlock;
    /// coupling of rows-cols on neighbouring processors (may not be valid)
    SparseMatrix_ptr remote_coupleBlock;

    bool is_balanced;

    /// matrix may be balanced by a diagonal matrix D=diagonal(balance_vector)
    /// if is_balanced is true, the matrix stored is D*A*D where A is the
    /// original matrix.
    /// When the system of linear equations is solved we solve D*A*D*y=c.
    /// So to solve A*x=b one needs to set c=D*b and x=D*y.
    double* balance_vector;

    /// stores the global ids for all cols in col_coupleBlock
    mutable index_t* global_id;

    /// package code controlling the solver pointer
    index_t solver_package;

    /// pointer to data needed by a solver
    void* solver_p;

    /// this is only used for a trilinos matrix
    void* trilinos_data; 
};


void SystemMatrix_MatrixVector(double alpha, SystemMatrix_ptr A, const double* in, double beta, double* out);

void SystemMatrix_MatrixVector_CSR_OFFSET0(double alpha, SystemMatrix_ptr A, const double* in, double beta, double* out);

void Paso_RHS_loadMM_toCSR(const char* filename, double* b, dim_t size);


} // namespace paso
  
#endif // __PASO_SYSTEMMATRIX_H__

