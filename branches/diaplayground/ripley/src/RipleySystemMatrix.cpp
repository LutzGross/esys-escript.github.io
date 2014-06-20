
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

#include "RipleySystemMatrix.h" 
#include "RipleyException.h" 

#include <esysUtils/index.h>
#include <escript/Data.h>
#include <escript/SolverOptions.h>

#include <cusp/multiply.h>
#include <cusp/krylov/bicgstab.h>
#include <cusp/krylov/cg.h>
#include <cusp/krylov/gmres.h>

#define BLOCK_SIZE 128

namespace ripley {

double gettime()
{
#ifdef _OPENMP
    return omp_get_wtime();
#else
    struct timeval tv;
    gettimeofday(&tv, NULL);
    suseconds_t ret = tv.tv_usec + tv.tv_sec*1e6;
    return 1e-6*(double)ret;
#endif
}

void list_devices(void)
{
#ifdef USE_CUDA
    int deviceCount;
    CUDA_SAFE_CALL(cudaGetDeviceCount(&deviceCount));
    if (deviceCount == 0)
        std::cout << "There is no device supporting CUDA" << std::endl;

    for (int dev = 0; dev < deviceCount; ++dev) {
        cudaDeviceProp deviceProp;
        CUDA_SAFE_CALL(cudaGetDeviceProperties(&deviceProp, dev));

        if (dev == 0) {
            if (deviceProp.major == 9999 && deviceProp.minor == 9999)
                std::cout << "There is no device supporting CUDA." << std::endl;
            else if (deviceCount == 1)
                std::cout << "There is 1 device supporting CUDA" << std:: endl;
            else
                std::cout << "There are " << deviceCount <<  " devices supporting CUDA" << std:: endl;
        }

        std::cout << "\nDevice " << dev << ": \"" << deviceProp.name << "\"" << std::endl;
        std::cout << "  Major revision number:                         " << deviceProp.major << std::endl;
        std::cout << "  Minor revision number:                         " << deviceProp.minor << std::endl;
        std::cout << "  Total amount of global memory:                 " << deviceProp.totalGlobalMem << " bytes" << std::endl;
    }
    std::cout << std::endl;
#endif
}

SystemMatrix::SystemMatrix()
{
    throw RipleyException("Illegal to instantiate SystemMatrix without arguments.");
}

SystemMatrix::SystemMatrix(int blocksize, const escript::FunctionSpace& fs,
                           int nRows, const IndexVector& diagonalOffsets) :
    AbstractSystemMatrix(blocksize, fs, blocksize, fs)
{
    // count nonzero entries
    int numEntries = 0;
    for (size_t i = 0; i < diagonalOffsets.size(); i++) {
        numEntries += blocksize*blocksize *
            (nRows-std::abs(diagonalOffsets[i])*blocksize) / blocksize;
    }

    list_devices();
#ifdef USE_CUDA
    //TODO: give users options...
    cudaSetDevice(0);
#endif

    mat.resize(nRows*blocksize, nRows*blocksize, numEntries, diagonalOffsets.size()*blocksize);
    mat.diagonal_offsets.assign(diagonalOffsets.begin(), diagonalOffsets.end());
}

void SystemMatrix::add(const IndexVector& rowIdx,
                       const std::vector<double>& array)
{
    const int blockSize = getBlockSize();
    const int emSize = rowIdx.size();
    for (int i=0; i<emSize; i++) {
        for (int j=0; j<emSize; j++) {
            const int revi = emSize-1-i;
            const int diag = j%2 + revi%2 + 3*(j/2+revi/2 + j/4+revi/4);
            for (int k=0; k<blockSize; k++) {
                for (int m=0; m<blockSize; m++) {
                    const int row = rowIdx[i]*blockSize + k;
                    const int d = diag*blockSize + m;
                    const int srcIdx = INDEX4(k, m, i, j, blockSize, blockSize, emSize);
                    mat.values(row, d) += array[srcIdx];
                }
            }
        }
    }
}

void SystemMatrix::ypAx(escript::Data& y, escript::Data& x) const
{
    if (x.getDataPointSize() != getBlockSize()) {
        throw RipleyException("matrix vector product: block size does not match the number of components in input.");
    } else if (y.getDataPointSize() != getBlockSize()) {
        throw RipleyException("matrix vector product: block size does not match the number of components in output.");
    } else if (x.getFunctionSpace() != getColumnFunctionSpace()) {
        throw RipleyException("matrix vector product: matrix function space and function space of input don't match.");
    } else if (y.getFunctionSpace() != getRowFunctionSpace()) {
        throw RipleyException("matrix vector product: matrix function space and function space of output don't match.");
    }

    std::cerr << "SystemMatrix::ypAx" << std::endl;
    // expand data object if necessary to be able to grab the whole data
    const_cast<escript::Data*>(&x)->expand();
    y.expand();
    y.requireWrite();
    const double* x_dp = x.getSampleDataRO(0);
    double* y_dp = y.getSampleDataRW(0);
    double T0 = gettime();
    HostVectorType xx(x_dp, x_dp+mat.num_rows);
    HostVectorType yy(mat.num_rows, 0.);
    cusp::multiply(mat, xx, yy);
    thrust::copy(yy.begin(), yy.end(), y_dp);
    std::cerr << "ypAx: " << gettime()-T0 << " seconds." << std::endl;
}

// y = beta y + Ax
void SystemMatrix::matrixVector(const double* x, double beta, double* y) const
{
    //cusp::multiply(mat, x, y);

    const int blockSize = getBlockSize();
    if (blockSize == 1) {
#pragma omp parallel for
        for (size_t ch=0; ch<mat.num_rows; ch+=BLOCK_SIZE) {
            // initialize chunk
            if (std::abs(beta) > 0.) {
                if (beta != 1.) {
                    for (int row=ch; row<std::min(ch+BLOCK_SIZE, mat.num_rows); row++) {
                        y[row] *= beta;
                    }
                }
            } else {
                for (size_t row=ch; row<std::min(ch+BLOCK_SIZE, mat.num_rows); row++) {
                    y[row] = 0.;
                }
            }

            // for each diagonal
            for (size_t d=0; d<mat.diagonal_offsets.size(); d++) {
                for (int row=ch; row<std::min(ch+BLOCK_SIZE, mat.num_rows); row++) {
                    const int col = row + mat.diagonal_offsets[d];
                    if (col >= 0 && col < mat.num_rows) {
                        y[row] += mat.values(row, d) * x[col];
                    }
                }
            }
        }
    } else {
#pragma omp parallel for
        for (size_t ch=0; ch < mat.num_rows; ch+=BLOCK_SIZE) {
            // initialize chunk
            if (std::abs(beta) > 0.) {
                if (beta != 1.) {
                    for (int row=ch; row<std::min(ch+BLOCK_SIZE, mat.num_rows); row++) {
                        y[row] *= beta;
                    }
                }
            } else {
                for (int row=ch; row<std::min(ch+BLOCK_SIZE, mat.num_rows); row++) {
                    y[row] = 0.;
                }
            }

            // for each diagonal block
            for (size_t d=0; d < mat.diagonal_offsets.size(); d++) {
                const int k = mat.diagonal_offsets[d]*blockSize;
                for (int row=ch; row<std::min(ch+BLOCK_SIZE, mat.num_rows); row++) {
                    const int col = blockSize*(row/blockSize) + k;
                    if (col >= 0 && col <= mat.num_rows-blockSize) {
                        // for each column in block
                        for (int i=0; i<blockSize; i++) {
                            const double Aij = mat.values(row, d*blockSize+i);
                            y[row] += Aij * x[col+i];
                        }
                    }
                }
            }
        }
    }
}

template<class LinearOperator,
         class Vector>
void SystemMatrix::runSolver(LinearOperator& A, Vector& x, Vector& b,
                             escript::SolverBuddy& sb) const
{
    typedef typename LinearOperator::memory_space MemorySpace;
    //sb.isVerbose()
    //sb.getPreconditioner()
    cusp::default_monitor<double> monitor(b, sb.getIterMax(), sb.getTolerance(), sb.getAbsoluteTolerance());
    cusp::identity_operator<double, MemorySpace> M(mat.num_rows, mat.num_rows);

    double T0 = gettime();
    switch (sb.getSolverMethod()) {
        case escript::ESCRIPT_DEFAULT:
        case escript::ESCRIPT_PCG:
            cusp::krylov::cg(A, x, b, monitor, M);
            break;
        case escript::ESCRIPT_BICGSTAB:
            cusp::krylov::bicgstab(A, x, b, monitor, M);
            break;
        case escript::ESCRIPT_GMRES:
            {
                const int restart = sb._getRestartForC();
                if (restart < 1)
                    throw RipleyException("Invalid restart parameter for GMRES");
                cusp::krylov::gmres(A, x, b, restart, monitor, M);
            }
            break;
        case escript::ESCRIPT_PRES20:
            {
                const int restart = 20;
                cusp::krylov::gmres(A, x, b, restart, monitor, M);
            }
            break;
        default:
            throw RipleyException("Unsupported solver.");
    }
    double solvertime = gettime()-T0;

    if (monitor.converged()) {
        std::cout << "Solver converged to " << monitor.relative_tolerance()
            << " relative tolerance after " << monitor.iteration_count()
            << " iterations and " << solvertime << " seconds."<< std::endl;
    } else {
        std::cout << "Solver reached iteration limit "
            << monitor.iteration_limit() << " before converging"
            << " to " << monitor.relative_tolerance() << " rel. tolerance."
            << std::endl;
    }
}

void SystemMatrix::setToSolution(escript::Data& out, escript::Data& in,
                                 boost::python::object& options) const
{
    if (out.getDataPointSize() != getBlockSize()) {
        throw RipleyException("solve: block size does not match the number of components of solution.");
    } else if (in.getDataPointSize() != getBlockSize()) {
        throw RipleyException("solve: block size does not match the number of components of right hand side.");
    } else if (out.getFunctionSpace() != getColumnFunctionSpace()) {
        throw RipleyException("solve: matrix function space and function space of solution don't match.");
    } else if (in.getFunctionSpace() != getRowFunctionSpace()) {
        throw RipleyException("solve: matrix function space and function space of right hand side don't match.");
    }

    std::cerr << "SystemMatrix::setToSolution" << std::endl;
    options.attr("resetDiagnostics")();
    escript::SolverBuddy sb = boost::python::extract<escript::SolverBuddy>(options);
    out.expand();
    in.expand();

    double* out_dp = out.getSampleDataRW(0);
    const double* in_dp = in.getSampleDataRO(0);
    double T0;

    if (sb.getSolverTarget() == escript::ESCRIPT_TARGET_GPU) {
#ifdef USE_CUDA
        T0 = gettime();
        DeviceVectorType b(in_dp, in_dp+mat.num_rows);
        double host2dev = gettime()-T0;
        std::cerr << "Copy of b: " << host2dev << " seconds." << std::endl;
        T0 = gettime();
        DeviceMatrixType dmat = mat;
        host2dev = gettime()-T0;
        std::cerr << "Copy of A: " << host2dev << " seconds." << std::endl;
        DeviceVectorType x(mat.num_rows, 0.);
        std::cerr << "Solving on CUDA device" << std::endl;
        runSolver(dmat, x, b, sb);
        T0 = gettime();
        thrust::copy(x.begin(), x.end(), out_dp);
        const double copyTime = gettime()-T0;
        std::cerr << "Copy of x: " << copyTime << " seconds." << std::endl;
#else
        throw RipleyException("solve: GPU-based solver requested but escript "
                              "not compiled with CUDA.");
#endif
    } else { // CPU
        T0 = gettime();
        HostVectorType b(in_dp, in_dp+mat.num_rows);
        double copytime = gettime()-T0;
        std::cerr << "Copy of b: " << copytime << " seconds." << std::endl;
        HostVectorType x(mat.num_rows, 0.);
        std::cerr << "Solving on host" << std::endl;
        runSolver(mat, x, b, sb);
        T0 = gettime();
        thrust::copy(x.begin(), x.end(), out_dp);
        const double copyTime = gettime()-T0;
        std::cerr << "Copy of x: " << copyTime << " seconds." << std::endl;
    }

}

void SystemMatrix::nullifyRowsAndCols(escript::Data& row_q,
                                      escript::Data& col_q,
                                      double mdv)
{
    double T0 = gettime();
    if (col_q.getDataPointSize() != getColumnBlockSize()) {
        throw RipleyException("nullifyRowsAndCols: column block size does not match the number of components of column mask.");
    } else if (row_q.getDataPointSize() != getRowBlockSize()) {
        throw RipleyException("nullifyRowsAndCols: row block size does not match the number of components of row mask.");
    } else if (col_q.getFunctionSpace() != getColumnFunctionSpace()) {
        throw RipleyException("nullifyRowsAndCols: column function space and function space of column mask don't match.");
    } else if (row_q.getFunctionSpace() != getRowFunctionSpace()) {
        throw RipleyException("nullifyRowsAndCols: row function space and function space of row mask don't match.");
    }

    std::cerr << "SystemMatrix::nullifyRowsAndCols mdv=" << mdv << std::endl;
    row_q.expand();
    col_q.expand();
    const double* rowMask = row_q.getSampleDataRO(0);
    const double* colMask = col_q.getSampleDataRO(0);
    const int blockSize = getBlockSize();
#pragma omp parallel for
    for (int row=0; row < mat.num_rows; row++) {
        for (int diag=0; diag < mat.diagonal_offsets.size(); diag++) {
            const int col = blockSize*(row/blockSize)+mat.diagonal_offsets[diag]*blockSize;
            if (col >= 0 && col <= mat.num_rows-blockSize) {
                for (int i=0; i<blockSize; i++) {
                    if (rowMask[row] > 0. || colMask[col+i] > 0.) {
                        mat.values(row, diag*blockSize+i) =
                                                        (row==col+i ? mdv : 0);
                    }
                }
            }
        }
    }
    std::cerr << "nullifyRowsAndCols: " << gettime()-T0 << " seconds." << std::endl;
}

void SystemMatrix::saveMM(const std::string& filename) const
{
    const int blockSize = getBlockSize();

    std::ofstream f(filename.c_str());
    f << "%%%%MatrixMarket matrix coordinate real general" << std::endl;
    f << mat.num_rows << " " << mat.num_rows << " " << mat.num_entries << std::endl;
    for (int row=0; row < mat.num_rows; row++) {
        for (int diag=0; diag < mat.diagonal_offsets.size(); diag++) {
            const int col = blockSize*(row/blockSize)+mat.diagonal_offsets[diag]*blockSize;
            if (col >= 0 && col <= mat.num_rows-blockSize) {
                for (int i=0; i<blockSize; i++) {
                    f << row+1 << " " << col+i+1 << " "
                          << mat.values(row, diag*blockSize+i) << std::endl;
                }
            }
        }
    }
}

void SystemMatrix::saveHB(const std::string& filename) const
{
    throw RipleyException("Harwell-Boeing interface not available.");
}

void SystemMatrix::resetValues()
{
    mat.values.values.assign(mat.values.values.size(), 0.);
}

}  // end of namespace
