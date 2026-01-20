
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
**
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
template <class T> class SystemMatrix;
template <typename T> using SystemMatrix_ptr = boost::shared_ptr<SystemMatrix<T> >;
template <typename T> using const_SystemMatrix_ptr = boost::shared_ptr<const SystemMatrix<T> >;

typedef int SystemMatrixType;

/// this class holds a (distributed) stiffness matrix
template <class T>
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
            SparseMatrix_ptr<T> merged(mergeSystemMatrix());
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
    SparseMatrix_ptr<T> mergeSystemMatrix() const;

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

    void MatrixVector(double alpha, const T* in, double beta,
                      T* out) const;

    void MatrixVector_CSR_OFFSET0(double alpha, const double* in, double beta,
                                  double* out) const;

    static SystemMatrix_ptr<double> loadMM_toCSR(const char* filename);

    static SystemMatrix_ptr<double> loadMM_toCSC(const char* filename);

    static int getSystemMatrixTypeId(int solver, int preconditioner,
                                     int package, bool is_complex, bool symmetry,
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
    SparseMatrix_ptr<T> mainBlock;
    /// coupling to neighbouring processors (row - col)
    SparseMatrix_ptr<T> col_coupleBlock;
    /// coupling to neighbouring processors (col - row)
    SparseMatrix_ptr<T> row_coupleBlock;
    /// coupling of rows-cols on neighbouring processors (may not be valid)
    SparseMatrix_ptr<T> remote_coupleBlock;

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

    void solve(T* out, T* in, Options* options) const;
};


void RHS_loadMM_toCSR(const char* filename, double* b, dim_t size);


} // namespace paso

#include "Options.h"
#include "Solver.h"

#include <escript/Data.h>

namespace paso {

template <>
SparseMatrix_ptr<double> PASO_DLL_API SystemMatrix<double>::mergeSystemMatrix() const;
template <>
SparseMatrix_ptr<cplx_t> PASO_DLL_API SystemMatrix<cplx_t>::mergeSystemMatrix() const;
template <>
void PASO_DLL_API SystemMatrix<double>::MatrixVector(double alpha, const double* in, double beta,
                                double* out) const;
template <>
void PASO_DLL_API SystemMatrix<cplx_t>::MatrixVector(double alpha, const cplx_t* in, double beta,
                                cplx_t* out) const;
template <>
void PASO_DLL_API SystemMatrix<double>::solve(double* out, double* in, Options* options) const;
template <>
void PASO_DLL_API SystemMatrix<cplx_t>::solve(cplx_t* out, cplx_t* in, Options* options) const;

template <class T>
SystemMatrix<T>::SystemMatrix()
{
    throw PasoException("SystemMatrix: Illegal to generate default SystemMatrix.");
}

/// Allocates a SystemMatrix of given type using the given matrix pattern.
/// Values are initialized with zero.
/// If patternIsUnrolled and type & MATRIX_FORMAT_BLK1, it is assumed
/// that the pattern is already unrolled to match the requested block size
/// and offsets. Otherwise unrolling and offset adjustment will be performed.
template <class T>
SystemMatrix<T>::SystemMatrix(SystemMatrixType ntype,
                           SystemMatrixPattern_ptr npattern, dim_t rowBlockSize,
                           dim_t colBlockSize, bool patternIsUnrolled,
                           const escript::FunctionSpace& rowFS,
                           const escript::FunctionSpace& colFS) :
    escript::AbstractSystemMatrix(rowBlockSize, rowFS, colBlockSize, colFS),
    type(ntype),
    logical_row_block_size(rowBlockSize),
    logical_col_block_size(colBlockSize),
    is_balanced(false),
    balance_vector(NULL),
    global_id(NULL),
    solver_package(PASO_PASO),
    solver_p(NULL)
{
    if (patternIsUnrolled) {
        if ((ntype & MATRIX_FORMAT_OFFSET1) != (npattern->type & MATRIX_FORMAT_OFFSET1)) {
            throw PasoException("SystemMatrix: requested offset and pattern offset do not match.");
        }
    }
    // do we need to apply unrolling?
    bool unroll
          // we don't like non-square blocks
        = (rowBlockSize != colBlockSize)
#ifndef ESYS_HAVE_LAPACK
          // or any block size bigger than 3
          || (colBlockSize > 3)
#endif
          // or if block size one requested and the block size is not 1
          || ((ntype & MATRIX_FORMAT_BLK1) && colBlockSize > 1)
          // or the offsets don't match
          || ((ntype & MATRIX_FORMAT_OFFSET1) != (npattern->type & MATRIX_FORMAT_OFFSET1));

    SystemMatrixType pattern_format_out = (ntype & MATRIX_FORMAT_OFFSET1)
                             ? MATRIX_FORMAT_OFFSET1 : MATRIX_FORMAT_DEFAULT;

    mpi_info = npattern->mpi_info;

    if (ntype & MATRIX_FORMAT_CSC) {
        if (unroll) {
            if (patternIsUnrolled) {
                pattern=npattern;
            } else {
                pattern = npattern->unrollBlocks(pattern_format_out,
                                                 colBlockSize, rowBlockSize);
            }
            row_block_size = 1;
            col_block_size = 1;
        } else {
            pattern = npattern->unrollBlocks(pattern_format_out, 1, 1);
            row_block_size = rowBlockSize;
            col_block_size = colBlockSize;
        }
        row_distribution = pattern->input_distribution;
        col_distribution = pattern->output_distribution;
    } else {
        if (unroll) {
            if (patternIsUnrolled) {
                pattern = npattern;
            } else {
                pattern = npattern->unrollBlocks(pattern_format_out,
                                                 rowBlockSize, colBlockSize);
            }
            row_block_size = 1;
            col_block_size = 1;
        } else {
            pattern = npattern->unrollBlocks(pattern_format_out, 1, 1);
            row_block_size = rowBlockSize;
            col_block_size = colBlockSize;
        }
        row_distribution = pattern->output_distribution;
        col_distribution = pattern->input_distribution;
    }
    if (ntype & MATRIX_FORMAT_DIAGONAL_BLOCK) {
        block_size = std::min(row_block_size, col_block_size);
    } else {
        block_size = row_block_size*col_block_size;
    }
    col_coupler.reset(new Coupler<real_t>(pattern->col_connector, col_block_size, mpi_info));
    row_coupler.reset(new Coupler<real_t>(pattern->row_connector, row_block_size, mpi_info));
    mainBlock.reset(new SparseMatrix<T>(type, pattern->mainPattern, row_block_size, col_block_size, true));
    col_coupleBlock.reset(new SparseMatrix<T>(type, pattern->col_couplePattern, row_block_size, col_block_size, true));
    row_coupleBlock.reset(new SparseMatrix<T>(type, pattern->row_couplePattern, row_block_size, col_block_size, true));
    const dim_t n_norm = std::max(mainBlock->numCols*col_block_size, mainBlock->numRows*row_block_size);
    balance_vector = new double[n_norm];
#pragma omp parallel for
    for (dim_t i=0; i<n_norm; ++i)
        balance_vector[i] = 1.;
}

// deallocates a SystemMatrix
template <class T>
SystemMatrix<T>::~SystemMatrix()
{
    solve_free(this);
    delete[] balance_vector;
    delete[] global_id;
}

template <class T>
int SystemMatrix<T>::getSystemMatrixTypeId(int solver, int preconditioner,
                                        int package, bool is_complex, bool symmetry,
                                        const escript::JMPI& mpi_info)
{
    int out = -1;
    int true_package = Options::getPackage(Options::mapEscriptOption(solver),
                                           Options::mapEscriptOption(package),
                                           symmetry, mpi_info);

    switch(true_package) {
        case PASO_PASO:
            out = MATRIX_FORMAT_DEFAULT;
        break;

        case PASO_MKL:
            out = MATRIX_FORMAT_BLK1 | MATRIX_FORMAT_OFFSET1;
        break;

        case PASO_UMFPACK:
            if (mpi_info->size > 1) {
                throw PasoException("The selected solver UMFPACK "
                        "requires CSC format which is not supported with "
                        "more than one rank.");
            } else {
                out = MATRIX_FORMAT_CSC | MATRIX_FORMAT_BLK1;
            }
        break;

        case PASO_MUMPS:
            out = MATRIX_FORMAT_BLK1 | MATRIX_FORMAT_OFFSET1;
        break;

        default:
            throw PasoException("unknown package code");
    }
    if (out > 0 && is_complex)
        out |= MATRIX_FORMAT_COMPLEX;
    return out;
}

template <class T>
void SystemMatrix<T>::nullifyRowsAndCols(escript::Data& row_q,
                                      escript::Data& col_q,
                                      double main_diagonal_value)
{
    if (row_q.isComplex() || col_q.isComplex())
    {
        throw PasoException("SystemMatrix::nullifyRowsAndCols: complex arguments not supported");      
    }
    if (col_q.getDataPointSize() != getColumnBlockSize()) {
        throw PasoException("nullifyRowsAndCols: column block size does not match the number of components of column mask.");
    } else if (row_q.getDataPointSize() != getRowBlockSize()) {
        throw PasoException("nullifyRowsAndCols: row block size does not match the number of components of row mask.");
    } else if (col_q.getFunctionSpace() != getColumnFunctionSpace()) {
        throw PasoException("nullifyRowsAndCols: column function space and function space of column mask don't match.");
    } else if (row_q.getFunctionSpace() != getRowFunctionSpace()) {
        throw PasoException("nullifyRowsAndCols: row function space and function space of row mask don't match.");
    }
    row_q.expand();
    col_q.expand();
    row_q.requireWrite();
    col_q.requireWrite();
    double* mask_row = row_q.getExpandedVectorReference(static_cast<escript::DataTypes::real_t>(0)).data();
    double* mask_col = col_q.getExpandedVectorReference(static_cast<escript::DataTypes::real_t>(0)).data();

    if (mpi_info->size > 1) {
        if (type & MATRIX_FORMAT_CSC) {
            throw PasoException("SystemMatrix::nullifyRowsAndCols: "
                                "CSC is not supported with MPI.");
        }

        startColCollect(mask_col);
        startRowCollect(mask_row);
        if (col_block_size==1 && row_block_size==1) {
            mainBlock->nullifyRowsAndCols_CSR_BLK1(mask_row, mask_col, main_diagonal_value);
            double* remote_values = finishColCollect();
            col_coupleBlock->nullifyRowsAndCols_CSR_BLK1(mask_row, remote_values, 0.);
            remote_values = finishRowCollect();
            row_coupleBlock->nullifyRowsAndCols_CSR_BLK1(remote_values, mask_col, 0.);
        } else {
            mainBlock->nullifyRowsAndCols_CSR(mask_row, mask_col, main_diagonal_value);
            double* remote_values = finishColCollect();
            col_coupleBlock->nullifyRowsAndCols_CSR(mask_row, remote_values, 0.);
            remote_values = finishRowCollect();
            row_coupleBlock->nullifyRowsAndCols_CSR(remote_values, mask_col, 0.);
        }
    } else {
        if (col_block_size==1 && row_block_size==1) {
            if (type & MATRIX_FORMAT_CSC) {
                mainBlock->nullifyRowsAndCols_CSC_BLK1(mask_row, mask_col, main_diagonal_value);
            } else {
                mainBlock->nullifyRowsAndCols_CSR_BLK1(mask_row, mask_col, main_diagonal_value);
            }
        } else {
            if (type & MATRIX_FORMAT_CSC) {
                mainBlock->nullifyRowsAndCols_CSC(mask_row, mask_col, main_diagonal_value);
            } else {
                mainBlock->nullifyRowsAndCols_CSR(mask_row, mask_col, main_diagonal_value);
            }
        }
    }
}

template <class T>
void SystemMatrix<T>::resetValues(bool preserveSolverData)
{
    setValues(0.);
    if (!preserveSolverData)
        solve_free(this);
}

template <class T>
void SystemMatrix<T>::setToSolution(escript::Data& out, escript::Data& in,
                                 boost::python::object& options) const
{
#if !defined(ESYS_HAVE_MUMPS)
    if (in.isComplex() || out.isComplex())
    {
        throw PasoException("SystemMatrix::setToSolution: complex arguments not supported.");
    }
#endif
    options.attr("resetDiagnostics")();
    Options paso_options(options);
    if (out.getDataPointSize() != getColumnBlockSize()) {
        throw PasoException("solve: column block size does not match the number of components of solution.");
    } else if (in.getDataPointSize() != getRowBlockSize()) {
        throw PasoException("solve: row block size does not match the number of components of  right hand side.");
    } else if (out.getFunctionSpace() != getColumnFunctionSpace()) {
        throw PasoException("solve: column function space and function space of solution don't match.");
    } else if (in.getFunctionSpace() != getRowFunctionSpace()) {
        throw PasoException("solve: row function space and function space of right hand side don't match.");
    }
    out.expand();
    in.expand();
    out.requireWrite();
    in.requireWrite();
    T* out_dp = out.getExpandedVectorReference(static_cast<T>(0)).data();
    T* in_dp = in.getExpandedVectorReference(static_cast<T>(0)).data();
    solve(out_dp, in_dp, &paso_options);
    paso_options.updateEscriptDiagnostics(options);
}

template <class T>
void SystemMatrix<T>::ypAx(escript::Data& y, escript::Data& x) const 
{
#if !defined(ESYS_HAVE_MUMPS)
    if (x.isComplex() || y.isComplex())
    {
        throw PasoException("SystemMatrix::ypAx: complex arguments not supported.");
    }  
#endif
    if (x.getDataPointSize() != getColumnBlockSize()) {
        throw PasoException("matrix vector product: column block size does not match the number of components in input.");
    } else if (y.getDataPointSize() != getRowBlockSize()) {
        throw PasoException("matrix vector product: row block size does not match the number of components in output.");
    } else if (x.getFunctionSpace() != getColumnFunctionSpace()) {
        throw PasoException("matrix vector product: column function space and function space of input don't match.");
    } else if (y.getFunctionSpace() != getRowFunctionSpace()) {
        throw PasoException("matrix vector product: row function space and function space of output don't match.");
    }
    x.expand();
    y.expand();
    x.requireWrite();
    y.requireWrite();
    T* x_dp = x.getExpandedVectorReference(static_cast<T>(0)).data();
    T* y_dp = y.getExpandedVectorReference(static_cast<T>(0)).data();
    MatrixVector(1., x_dp, 1., y_dp);
}

} // namespace paso

#endif // __PASO_SYSTEMMATRIX_H__

