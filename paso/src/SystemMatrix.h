
/*****************************************************************************
*
* Copyright (c) 2003-2018 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
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

#include "SparseMatrix.h"
#include "SystemMatrixPattern.h"

#include <escript/AbstractSystemMatrix.h>

namespace paso {

struct Options;
class SystemMatrix;
typedef boost::shared_ptr<SystemMatrix> SystemMatrix_ptr;
typedef boost::shared_ptr<const SystemMatrix> const_SystemMatrix_ptr;

typedef int SystemMatrixType;

/// this class holds a (distributed) stiffness matrix
class SystemMatrix : public escript::AbstractSystemMatrix
{
public:
    /// default constructor - throws exception.
    SystemMatrix();

    SystemMatrix(SystemMatrixType type, SystemMatrixPattern_ptr pattern,
                 dim_t rowBlockSize, dim_t columnBlockSize,
                 bool patternIsUnrolled, const escript::FunctionSpace& rowFS,
                 const escript::FunctionSpace& colFS);

    ~SystemMatrix();

    /// Nullifies rows and columns in the matrix.
    /// The rows and columns are marked by positive values in mask_row and
    /// mask_col. Values on the main diagonal which are marked to set to
    /// zero by both mask_row and mask_col are set to main_diagonal_value.
    virtual void nullifyRowsAndCols(escript::Data& mask_row,
                                    escript::Data& mask_col,
                                    double main_diagonal_value);

    virtual inline void saveMM(const std::string& filename) const
    {
        if (mpi_info->size > 1) {
            //throw PasoException("SystemMatrix::saveMM: Only single rank supported.");
            SparseMatrix_ptr merged(mergeSystemMatrix());
            if (mpi_info->rank == 0)
                merged->saveMM(filename.c_str());
        } else {
            mainBlock->saveMM(filename.c_str());
        }
    }

    virtual inline void saveHB(const std::string& filename) const
    {
        if (mpi_info->size > 1) {
            throw PasoException("SystemMatrix::saveHB: Only single rank supported.");
        } else if (!(type & MATRIX_FORMAT_CSC)) {
            throw PasoException("SystemMatrix::saveHB: Only CSC format supported.");
        } else {
            mainBlock->saveHB_CSC(filename.c_str());
        }
    }

    virtual void resetValues(bool preserveSolverData = false);

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

    inline void startCollect(const double* in) const
    {
        startColCollect(in);
    }

    inline double* finishCollect() const
    {
        return finishColCollect();
    }

    inline void startColCollect(const double* in) const
    {
        col_coupler->startCollect(in);
    }

    inline double* finishColCollect() const
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
                 ((double)getGlobalTotalNumRows()*getGlobalTotalNumCols());
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

    inline void rowSum(double* row_sum) const
    {
        if ((type & MATRIX_FORMAT_CSC) || (type & MATRIX_FORMAT_OFFSET1)) {
            throw PasoException("SystemMatrix::rowSum: No normalization "
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

    void MatrixVector(double alpha, const double* in, double beta,
                      double* out) const;

    void MatrixVector_CSR_OFFSET0(double alpha, const double* in, double beta,
                                  double* out) const;

    static SystemMatrix_ptr loadMM_toCSR(const char* filename);

    static SystemMatrix_ptr loadMM_toCSC(const char* filename);

    static int getSystemMatrixTypeId(int solver, int preconditioner,
                                     int package, bool symmetry,
                                     const escript::JMPI& mpi_info);

    SystemMatrixType type;
    SystemMatrixPattern_ptr pattern;

    dim_t logical_row_block_size;
    dim_t logical_col_block_size;

    dim_t row_block_size;
    dim_t col_block_size;
    dim_t block_size;

    escript::Distribution_ptr row_distribution;
    escript::Distribution_ptr col_distribution;
    escript::JMPI mpi_info;

    Coupler_ptr<real_t> col_coupler;
    Coupler_ptr<real_t> row_coupler;

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
    mutable index_t solver_package;

    /// pointer to data needed by a solver
    void* solver_p;

private:
    virtual void setToSolution(escript::Data& out, escript::Data& in,
                               boost::python::object& options) const;

    virtual void ypAx(escript::Data& y, escript::Data& x) const;

    void solve(double* out, double* in, Options* options) const;
};


void RHS_loadMM_toCSR(const char* filename, double* b, dim_t size);


} // namespace paso

#endif // __PASO_SYSTEMMATRIX_H__

