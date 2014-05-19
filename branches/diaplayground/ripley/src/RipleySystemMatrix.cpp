
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

#define BLOCK_SIZE 128

namespace ripley {

SystemMatrix::SystemMatrix()
{
    throw RipleyException("Illegal to instantiate SystemMatrix without arguments.");
}

SystemMatrix::SystemMatrix(int blocksize, const escript::FunctionSpace& fs,
                           int nRows, const IndexVector& diagonalOffsets) :
    AbstractSystemMatrix(blocksize, fs, blocksize, fs),
    numRows(nRows*blocksize),
    offsets(diagonalOffsets)
{
    values.resize(numRows*offsets.size()*blocksize, 0.);
}

void SystemMatrix::add(const IndexVector& rowIdx,
                       const std::vector<double>& array)
{
    const int blockSize = getBlockSize();
    const int emSize = rowIdx.size();
#if 0
    static bool here=false;
    if (here) return;
    here=true;
    double* arr = const_cast<double*>(&array[0]);

    for (int i=0; i<4*blockSize; i++) {
        for (int j=0; j<4*blockSize; j++) {
            arr[i*4*blockSize+j] = i*4*blockSize+j;
        }
    }
#endif

    //for (k in numEq) for (m in numComp)
    //array[k + numEq*m + numEq*numComp*i + numEq*numComp*emSize*j]
    //std::cerr << "SystemMatrix::add" << std::endl;
    for (int i=0; i<emSize; i++) {
        for (int j=0; j<emSize; j++) {
            const int revi = emSize-1-i;
            const int diag = j%2 + revi%2 + 3*(j/2+revi/2 + j/4+revi/4);
            for (int k=0; k<blockSize; k++) {
                for (int m=0; m<blockSize; m++) {
                    const int destIdx = rowIdx[i]*blockSize + k + numRows*(diag*blockSize+m);
                    const int srcIdx = INDEX4(k, m, i, j, blockSize, blockSize, emSize);
                    values[destIdx] += array[srcIdx];
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
    matrixVector(x_dp, 1., y_dp);
}

// y = beta y + Ax
void SystemMatrix::matrixVector(const double* x, double beta, double* y) const
{
    const int blockSize = getBlockSize();
    if (blockSize == 1) {
#pragma omp parallel for
        for (int ch=0; ch<numRows; ch+=BLOCK_SIZE) {
            // initialize chunk
            if (std::abs(beta) > 0.) {
                if (beta != 1.) {
                    for (int row=ch; row<std::min(ch+BLOCK_SIZE,numRows); row++) {
                        y[row] *= beta;
                    }
                }
            } else {
                for (int row=ch; row<std::min(ch+BLOCK_SIZE,numRows); row++) {
                    y[row] = 0.;
                }
            }

            // for each diagonal
            for (size_t d=0; d<offsets.size(); d++) {
                for (int row=ch; row<std::min(ch+BLOCK_SIZE,numRows); row++) {
                    const int col = row + offsets[d];
                    if (col >= 0 && col < numRows) {
                        y[row] += values[row+d*numRows] * x[col];
                    }
                }
            }
        }
    } else {
#pragma omp parallel for
        for (int ch=0; ch<numRows; ch+=BLOCK_SIZE) {
            // initialize chunk
            if (std::abs(beta) > 0.) {
                if (beta != 1.) {
                    for (int row=ch; row<std::min(ch+BLOCK_SIZE,numRows); row++) {
                        y[row] *= beta;
                    }
                }
            } else {
                for (int row=ch; row<std::min(ch+BLOCK_SIZE,numRows); row++) {
                    y[row] = 0.;
                }
            }

            // for each diagonal block
            for (size_t d=0; d<offsets.size(); d++) {
                const int k = offsets[d]*blockSize;
                for (int row=ch; row<std::min(ch+BLOCK_SIZE,numRows); row++) {
                    const int col = blockSize*(row/blockSize) + k;
                    if (col >= 0 && col <= numRows-blockSize) {
                        // for each column in block
                        for (int i=0; i<blockSize; i++) {
                            const double Aij = values[row+(d*blockSize+i)*numRows];
                            y[row] += Aij * x[col+i];
                        }
                    }
                }
            }
        }
    }
}

inline void axpby(size_t N, const double* x, const double* y, double* z,
                  double alpha, double beta)
{
#pragma omp parallel for
    for (size_t i=0; i<N; i++) {
        z[i] = alpha*x[i] + beta*y[i];
    }
}

inline void axpy(size_t N, const double* y, double* x, double alpha)
{
#pragma omp parallel for
    for (size_t i=0; i<N; i++) {
        x[i] += alpha*y[i];
    }
}

inline double dotc(const std::vector<double>& x, const std::vector<double>& y)
{
    double ret = 0.;
    for (size_t i=0; i<x.size(); i++)
        ret += x[i]*y[i];
    return ret;
}

inline double nrm2(size_t N, const double* x)
{
    double ret = 0.;
    for (size_t i=0; i<N; i++)
        ret += x[i]*x[i];
    return std::sqrt(ret);
}

void SystemMatrix::cg(double* x, const double* b) const
{
    const size_t N = numRows;

    // allocate workspace
    std::vector<double> y(N);
    std::vector<double> z(N);
    std::vector<double> r(N);
    std::vector<double> p(N);

#pragma omp parallel for
    for (size_t i=0; i<N; i++)
        x[i] = b[i];

    // y <- Ax
    matrixVector(x, 0., &y[0]);

    // r <- b - A*x
    axpby(N, b, &y[0], &r[0], 1., -1.);
   
    // z <- M*r
    //cusp::multiply(M, r, z);
    z=r;

    // p <- z
    //blas::copy(z, p);
#pragma omp parallel for
    for (size_t i=0; i<N; i++)
        p[i] = z[i];

    // rz = <r^H, z>
    double rz = dotc(r, z);

    const double b_norm = nrm2(N, b);

    int iteration=0;

    while (nrm2(N, &r[0]) > 1e-9*b_norm && iteration<10000)
    {
        //std::cerr << nrm2(N, &r[0])<< " " << nrm2(N, x)<<std::endl;
        iteration++;
        // y <- Ap
        matrixVector(&p[0], 0., &y[0]);

        // alpha <- <r,z>/<y,p>
        const double alpha =  rz / dotc(y, p);

        // x <- x + alpha * p
        axpy(N, &p[0], x, alpha);

        // r <- r - alpha * y		
        axpy(N, &y[0], &r[0], -alpha);

        // z <- M*r
        //cusp::multiply(M, r, z);
#pragma omp parallel for
        for (size_t i=0; i<N; i++)
            z[i] = r[i];
		
        double rz_old = rz;

        // rz = <r^H, z>
        rz = dotc(r, z);

        // beta <- <r_{i+1},r_{i+1}>/<r,r> 
        double beta = rz / rz_old;
		
        // p <- r + beta*p
        axpby(N, &z[0], &p[0], &p[0], 1., beta);
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
    out.expand();
    in.expand();
    double* out_dp = out.getSampleDataRW(0);
    const double* in_dp = in.getSampleDataRO(0);
    //solve(out_dp, in_dp, &paso_options);
    cg(out_dp, in_dp);
    //std::cerr << out.toString() << std::endl;
}

void SystemMatrix::nullifyRowsAndCols(escript::Data& row_q,
                                      escript::Data& col_q,
                                      double mdv)
{
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
    for (int row=0; row < numRows; row++) {
        for (int diag=0; diag < offsets.size(); diag++) {
            const int col = blockSize*(row/blockSize)+offsets[diag]*blockSize;
            if (col >= 0 && col <= numRows-blockSize) {
                for (int i=0; i<blockSize; i++) {
                    if (rowMask[row] > 0. || colMask[col+i] > 0.) {
                        values[row+(diag*blockSize+i)*numRows] =
                                                        (row==col+i ? mdv : 0);
                    }
                }
            }
        }
    }
}

void SystemMatrix::saveMM(const std::string& filename) const
{
    const int blockSize = getBlockSize();

    // count nonzero entries
    int numNZ = 0;
    for (size_t i = 0; i < offsets.size(); i++) {
        numNZ += blockSize*blockSize*(numRows-std::abs(offsets[i])*blockSize) / blockSize;
    }
    
    std::ofstream f(filename.c_str());
    f << "%%%%MatrixMarket matrix coordinate real general" << std::endl;
    f << numRows << " " << numRows << " " << numNZ << std::endl;
    for (int row=0; row<numRows; row++) {
        for (int diag=0; diag<offsets.size(); diag++) {
            const int col = blockSize*(row/blockSize)+offsets[diag]*blockSize;
            if (col >= 0 && col <= numRows-blockSize) {
                for (int i=0; i<blockSize; i++) {
                    f << row+1 << " " << col+i+1 << " "
                          << values[(diag*blockSize+i)*numRows+row] << std::endl;
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
    values.assign(values.size(), 0.);
}

}  // end of namespace
