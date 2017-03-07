
/*****************************************************************************
*
* Copyright (c) 2014-2017 by The University of Queensland
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

#ifndef __RIPLEY_SYSTEMMATRIX_H__
#define __RIPLEY_SYSTEMMATRIX_H__

#include <escript/AbstractSystemMatrix.h>
#include <escript/FunctionSpace.h>

#include <ripley/Ripley.h>

namespace escript {
class SolverBuddy;
}

#include <cusp/cds_matrix.h>

namespace ripley {

typedef cusp::cds_matrix<int, double, cusp::host_memory> HostMatrixType;
typedef cusp::cds_matrix<int, double, cusp::device_memory> DeviceMatrixType;
typedef cusp::array1d<double, cusp::host_memory> HostVectorType;
typedef cusp::array1d<double, cusp::device_memory> DeviceVectorType;

class SystemMatrix : public escript::AbstractSystemMatrix
{
public:
    SystemMatrix(escript::JMPI mpiInfo, int blocksize,
                 const escript::FunctionSpace& fs, int nRows,
                 const IndexVector& diagonalOffsets, bool symmetric);

    virtual ~SystemMatrix() {}

    virtual void nullifyRowsAndCols(escript::Data& row_q,
                                    escript::Data& col_q,
                                    double mdv);  
  
    virtual void saveMM(const std::string& filename) const;

    virtual void saveHB(const std::string& filename) const;

    virtual void resetValues(bool preserveSolverData = false);

    void add(const IndexVector& rowIndex, const std::vector<double>& array);

    inline int getBlockSize() const { return getRowBlockSize(); }

private:
    template<class LinearOperator, class Vector, class Preconditioner>
    void runSolver(LinearOperator& A, Vector& x, Vector& b, Preconditioner& M,
                   escript::SolverBuddy& sb) const;

    /// copies the current matrix stored on host to device *if required*
    void copyMatrixToDevice(bool verbose=false) const;

    virtual void setToSolution(escript::Data& out, escript::Data& in,
                               boost::python::object& options) const;

    virtual void ypAx(escript::Data& y, escript::Data& x) const;

    static void checkCUDA();

    /// GPU device IDs supporting CUDA
    static std::vector<int> cudaDevices;

    escript::JMPI m_mpiInfo;
    HostMatrixType mat;
    mutable DeviceMatrixType dmat;
    mutable bool matrixAltered;
    bool symmetric;
};

} // namespace ripley

#endif // __RIPLEY_SYSTEMMATRIX_H__

