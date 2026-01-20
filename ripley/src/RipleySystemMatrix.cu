
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/

#include "RipleySystemMatrix.h" 
#include "RipleyException.h" 

#include <escript/index.h>
#include <escript/Data.h>
#include <escript/SolverOptions.h>

#include <cusp/multiply.h>
#include <cusp/krylov/bicgstab.h>
#include <cusp/krylov/cg.h>
#include <cusp/krylov/cgls.h>
#include <cusp/krylov/gmres.h>
#include <cusp/krylov/lsqr.h>
#include <cusp/precond/diagonal.h>

namespace bp = boost::python;

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

std::vector<int> SystemMatrix::cudaDevices;

void SystemMatrix::checkCUDA()
{
#ifdef USE_CUDA
    cudaDevices.clear();
    int deviceCount = 0;
    cudaError err = cudaGetDeviceCount(&deviceCount);
    if (err == cudaErrorInsufficientDriver) {
        std::cout << "Note: CUDA reported a driver error enumerating CUDA devices" << std::endl;
        cudaDevices.push_back(-1);
        return;
    } else if (deviceCount == 0 || err == cudaErrorNoDevice) {
        std::cout << "Note: There is no device supporting CUDA" << std::endl;
        cudaDevices.push_back(-1);
        return;
    } else if (deviceCount == 1) {
        std::cout << "There is 1 GPU device" << std::endl;
    } else {
        std::cout << "There are " << deviceCount << " GPU devices" << std::endl;
    }

    for (int dev = 0; dev < deviceCount; ++dev) {
        cudaDeviceProp deviceProp;
        CUDA_SAFE_CALL(cudaGetDeviceProperties(&deviceProp, dev));

        std::cout << "\nDevice " << dev << ": \"" << deviceProp.name << "\" -- ";
        if (deviceProp.major == 9999 && deviceProp.minor == 9999) {
            std::cout << "No CUDA support." << std::endl;
        } else {
            cudaDevices.push_back(dev);
            std::cout << "Rev. " << deviceProp.major << "." << deviceProp.minor
                      << ", Mem: "
                      << deviceProp.totalGlobalMem << " bytes." << std::endl;
        }
    }
    // neither of the found devices supports CUDA:
    if (cudaDevices.empty()) {
        cudaDevices.push_back(-1);
    }
#endif
}

SystemMatrix::SystemMatrix(escript::JMPI mpiInfo, int blocksize,
                           const escript::FunctionSpace& fs, int nRows,
                           const IndexVector& diagonalOffsets, bool symm) :
    AbstractSystemMatrix(blocksize, fs, blocksize, fs),
    m_mpiInfo(mpiInfo),
    symmetric(symm)
{
    // count nonzero entries
    int numEntries = 0;
    for (size_t i = 0; i < diagonalOffsets.size(); i++) {
        numEntries += blocksize*blocksize*(nRows-std::abs(diagonalOffsets[i]));
    }

    mat.resize(nRows*blocksize, numEntries, diagonalOffsets.size(), blocksize);
    mat.diagonal_offsets.assign(diagonalOffsets.begin(), diagonalOffsets.end());
    mat.symmetric = symmetric;
    matrixAltered = true;
}

void SystemMatrix::add(const IndexVector& rowIdx,
                       const std::vector<double>& array)
{
    const int blockSize = getBlockSize();
    const int emSize = rowIdx.size();
    if (symmetric) {
        const int offset = mat.diagonal_offsets.size()-1; // 4 in 2D, 12 in 3D
        for (int i=0; i<emSize; i++) {
            for (int j=i; j<emSize; j++) {
                const int revi = emSize-1-i;
                // This ugliness is to convert matrix diagonal index to an
                // index in the storage scheme (i.e. main diagonal = 0)
                const int diag = j%2 + revi%2 + 3*(j/2+revi/2 + j/4+revi/4)-offset;
                for (int k=0; k<blockSize; k++) {
                    for (int m=0; m<blockSize; m++) {
                        const int row = rowIdx[i]*blockSize + k;
                        const int d = diag*blockSize + m;
                        const int srcIdx =
                            INDEX4(k, m, i, j, blockSize, blockSize, emSize);
                        mat.values(row, d) += array[srcIdx];
                    }
                }
            }
        }
    } else {
        for (int i=0; i<emSize; i++) {
            for (int j=0; j<emSize; j++) {
                const int revi = emSize-1-i;
                // This ugliness is to convert matrix diagonal index to an
                // index in the storage scheme (i.e. lower-most diagonal = 0)
                const int diag = j%2 + revi%2 + 3*(j/2+revi/2 + j/4+revi/4);
                for (int k=0; k<blockSize; k++) {
                    for (int m=0; m<blockSize; m++) {
                        const int row = rowIdx[i]*blockSize + k;
                        const int d = diag*blockSize + m;
                        const int srcIdx =
                            INDEX4(k, m, i, j, blockSize, blockSize, emSize);
                        mat.values(row, d) += array[srcIdx];
                    }
                }
            }
        }
    }
    matrixAltered = true;
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

    // expand data object if necessary to be able to grab the whole data
    const_cast<escript::Data*>(&x)->expand();
    y.expand();
    y.requireWrite();
    const double* x_dp = x.getSampleDataRO(0);
    double* y_dp = y.getSampleDataRW(0);
    //double T0 = gettime();
    HostVectorType xx(x_dp, x_dp+mat.num_rows);
    HostVectorType yy(mat.num_rows, 0.);
    cusp::multiply(mat, xx, yy);
    /*
    DeviceVectorType xx(x_dp, x_dp+mat.num_rows);
    DeviceVectorType yy(mat.num_rows, 0.);
    copyMatrixToDevice();
    cusp::multiply(dmat, xx, yy);
    */
    thrust::copy(yy.begin(), yy.end(), y_dp);
    //std::cout << "ypAx: " << gettime()-T0 << " seconds." << std::endl;
}

void SystemMatrix::copyMatrixToDevice(bool verbose) const
{
    if (matrixAltered) {
        double T0 = gettime();
        dmat = mat;
        if (verbose) {
            double host2dev = gettime()-T0;
            std::cout << "Copy of A: " << host2dev << " seconds." << std::endl;
        }
        matrixAltered = false;
    }
}

template<class LinearOperator,
         class Vector,
         class Preconditioner>
void SystemMatrix::runSolver(LinearOperator& A, Vector& x, Vector& b,
                             Preconditioner& M, escript::SolverBuddy& sb) const
{
    //cusp::verbose_monitor<double> monitor(b, sb.getIterMax(), sb.getTolerance(), sb.getAbsoluteTolerance());
    cusp::default_monitor<double> monitor(b, sb.getIterMax(), sb.getTolerance(), sb.getAbsoluteTolerance());
    int solver = sb.getSolverMethod();
    if (solver == escript::SO_DEFAULT) {
        if (sb.isSymmetric()) {
            solver = escript::SO_METHOD_PCG;
        } else {
            solver = escript::SO_METHOD_BICGSTAB;
        }
    }

    double T0 = gettime();
    switch (solver) {
        case escript::SO_DEFAULT:
        case escript::SO_METHOD_PCG:
            cusp::krylov::cg(A, x, b, monitor, M);
            break;
        case escript::SO_METHOD_CGLS:
            cusp::krylov::cgls(A, x, b, 0., monitor); //shift=0 for now
            break;
        case escript::SO_METHOD_LSQR:
            cusp::krylov::lsqr(A, x, b, cusp::krylov::lsqr_parameters<double>(), monitor);
            break;
        case escript::SO_METHOD_BICGSTAB:
            cusp::krylov::bicgstab(A, x, b, monitor, M);
            break;
        case escript::SO_METHOD_GMRES:
            {
                const int restart = (sb.getRestart()==0 ? 1000 : sb.getRestart());
                if (restart < 1)
                    throw RipleyException("Invalid restart parameter for GMRES");
                cusp::krylov::gmres(A, x, b, restart, monitor, M);
            }
            break;
        case escript::SO_METHOD_PRES20:
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
        sb.updateDiagnostics("converged", true);
        if (sb.isVerbose()) {
            std::cout << "Solver converged to " << monitor.relative_tolerance()
                << " relative tolerance after " << monitor.iteration_count()
                << " iterations and " << solvertime << " seconds."<< std::endl;
        }
    } else {
        std::cout << "Solver reached iteration limit "
            << monitor.iteration_limit() << " before converging"
            << " to " << monitor.relative_tolerance() << " rel. tolerance."
            << std::endl;
    }
    sb.updateDiagnostics("num_iter", monitor.iteration_count());
    sb.updateDiagnostics("net_time", solvertime);
}

void SystemMatrix::setToSolution(escript::Data& out, escript::Data& in,
                                 bp::object& options) const
{
    if (m_mpiInfo->size > 1) {
        throw RipleyException("solve: ripley's block diagonal matrix "
                              "is incompatible with MPI.");
    }
    if (out.getDataPointSize() != getBlockSize()) {
        throw RipleyException("solve: block size does not match the number of components of solution.");
    } else if (in.getDataPointSize() != getBlockSize()) {
        throw RipleyException("solve: block size does not match the number of components of right hand side.");
    } else if (out.getFunctionSpace() != getColumnFunctionSpace()) {
        throw RipleyException("solve: matrix function space and function space of solution don't match.");
    } else if (in.getFunctionSpace() != getRowFunctionSpace()) {
        throw RipleyException("solve: matrix function space and function space of right hand side don't match.");
    }

    options.attr("resetDiagnostics")();
    escript::SolverBuddy& sb = bp::extract<escript::SolverBuddy&>(options);
    out.expand();
    in.expand();

    if (sb.isVerbose()) {
        std::cout << "Matrix has " << mat.num_entries << " entries." << std::endl;
    }
    double* out_dp = out.getSampleDataRW(0);
    const double* in_dp = in.getSampleDataRO(0);
    double T0;

    if (sb.getSolverTarget() == escript::SO_TARGET_GPU) {
#ifdef USE_CUDA
        if (cudaDevices.empty()) {
            checkCUDA();
        }

        if (cudaDevices[0] == -1) {
            throw RipleyException("solve: GPU-based solver requested but no "
                                  "CUDA compatible device available.");
        }

        //TODO: give users options...
        if (sb.isVerbose()) {
            std::cout << "Using CUDA device " << cudaDevices[0] << std::endl;
        }
        cudaSetDevice(cudaDevices[0]);

        T0 = gettime();
        DeviceVectorType b(in_dp, in_dp+mat.num_rows);
        double host2dev = gettime()-T0;
        if (sb.isVerbose())
            std::cout << "Copy of b: " << host2dev << " seconds." << std::endl;
        copyMatrixToDevice(sb.isVerbose());
        DeviceVectorType x(mat.num_rows, 0.);
        if (sb.isVerbose())
            std::cout << "Solving on CUDA device..." << std::endl;

        if (sb.getPreconditioner() == escript::SO_PRECONDITIONER_NONE) {
            cusp::identity_operator<double, cusp::device_memory> M(mat.num_rows, mat.num_rows);
            runSolver(dmat, x, b, M, sb);
        } else if (sb.getPreconditioner() == escript::SO_PRECONDITIONER_JACOBI) {
            if (sb.isVerbose())
                std::cout << "Using Jacobi preconditioner" << std::endl;
            // TODO: This should be cached as well but that's not supported
            // at the moment.
            cusp::precond::diagonal<double, cusp::device_memory> M(dmat);
            runSolver(dmat, x, b, M, sb);
        } else {
            throw RipleyException("Unsupported preconditioner requested.");
        }

        T0 = gettime();
        thrust::copy(x.begin(), x.end(), out_dp);
        const double copyTime = gettime()-T0;
        if (sb.isVerbose())
            std::cout << "Copy of x: " << copyTime << " seconds." << std::endl;
#else
        throw RipleyException("solve: GPU-based solver requested but escript "
                              "not compiled with CUDA.");
#endif
    } else { // CPU
        T0 = gettime();
        HostVectorType b(in_dp, in_dp+mat.num_rows);
        double copytime = gettime()-T0;
        if (sb.isVerbose()) {
            std::cout << "Copy of b: " << copytime << " seconds." << std::endl;
            std::cout << "Solving on the CPU..." << std::endl;
        }
        HostVectorType x(mat.num_rows, 0.);
        if (sb.getPreconditioner() == escript::SO_PRECONDITIONER_NONE) {
            cusp::identity_operator<double, cusp::host_memory> M(mat.num_rows, mat.num_rows);
            runSolver(mat, x, b, M, sb);
        } else if (sb.getPreconditioner() == escript::SO_PRECONDITIONER_JACOBI) {
            if (sb.isVerbose())
                std::cout << "Using Jacobi preconditioner" << std::endl;
            // TODO: This should be cached as well but that's not supported
            // at the moment.
            cusp::precond::diagonal<double, cusp::host_memory> M(mat);
            runSolver(mat, x, b, M, sb);
        } else {
            throw RipleyException("Unsupported preconditioner requested.");
        }

        T0 = gettime();
        thrust::copy(x.begin(), x.end(), out_dp);
        const double copyTime = gettime()-T0;
        if (sb.isVerbose()) {
            std::cout << "Copy of x: " << copyTime << " seconds." << std::endl;
        }
    }
}

void SystemMatrix::nullifyRowsAndCols(escript::Data& row_q,
                                      escript::Data& col_q,
                                      double mdv)
{
    //double T0 = gettime();
    if (col_q.getDataPointSize() != getColumnBlockSize()) {
        throw RipleyException("nullifyRowsAndCols: column block size does not match the number of components of column mask.");
    } else if (row_q.getDataPointSize() != getRowBlockSize()) {
        throw RipleyException("nullifyRowsAndCols: row block size does not match the number of components of row mask.");
    } else if (col_q.getFunctionSpace() != getColumnFunctionSpace()) {
        throw RipleyException("nullifyRowsAndCols: column function space and function space of column mask don't match.");
    } else if (row_q.getFunctionSpace() != getRowFunctionSpace()) {
        throw RipleyException("nullifyRowsAndCols: row function space and function space of row mask don't match.");
    }

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
    //std::cout << "nullifyRowsAndCols: " << gettime()-T0 << " seconds." << std::endl;
    matrixAltered = true;
}

void SystemMatrix::saveMM(const std::string& filename) const
{
    const int blockSize = getBlockSize();

    std::ofstream f(filename.c_str());
    f << "%%MatrixMarket matrix coordinate real general" << std::endl;
    f << mat.num_rows << " " << mat.num_cols << " ";
    // here we assume the matrix has a main diagonal block!
    if (symmetric)
       f << 2*mat.num_entries-blockSize*mat.num_rows << std::endl;
    else
       f << mat.num_entries << std::endl;

    f.setf(std::ios_base::scientific, std::ios_base::floatfield);
    f.precision(15);
    for (int row=0; row < mat.num_rows; row++) {
        if (symmetric) {
            for (int diag=mat.diagonal_offsets.size()-1; diag >= 0; diag--) {
                // skip main diagonal which is handled below
                if (mat.diagonal_offsets[diag] == 0)
                    continue;
                const int col = blockSize*(row/blockSize)-mat.diagonal_offsets[diag]*blockSize;
                if (col >= 0 && col <= mat.num_rows-blockSize) {
                    for (int i=0; i<blockSize; i++) {
                        f << row+1 << " " << col+i+1 << " "
                          << mat.values(col+i, diag*blockSize+row%blockSize)
                          << std::endl;
                    }
                }
            }
        }
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

void SystemMatrix::resetValues(bool preserveSolverData)
{
    mat.values.values.assign(mat.values.values.size(), 0.);
    matrixAltered = true;
}

}  // end of namespace
