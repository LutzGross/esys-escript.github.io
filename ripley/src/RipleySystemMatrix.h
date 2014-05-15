
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

#ifndef __RIPLEY_SYSTEMMATRIX_H__
#define __RIPLEY_SYSTEMMATRIX_H__

#include <escript/AbstractSystemMatrix.h>
#include <escript/FunctionSpace.h>

#include <ripley/Ripley.h>

namespace ripley {

class SystemMatrix : public escript::AbstractSystemMatrix
{
public:

    SystemMatrix();

    SystemMatrix(int row_blocksize, const escript::FunctionSpace& row_fs,
                 int col_blocksize, const escript::FunctionSpace& col_fs,
                 int nRows, const IndexVector& diagonalOffsets);

    virtual ~SystemMatrix() {}

    virtual void nullifyRowsAndCols(escript::Data& row_q,
                                    escript::Data& col_q,
                                    double mdv);  
  
    virtual void saveMM(const std::string& filename) const;

    virtual void saveHB(const std::string& filename) const;

    virtual void resetValues();

    void add(const IndexVector& rowIndex, const std::vector<double>& array);

    void matrixVector(const double* x, double beta, double* y) const;

private:
    void cg(double* x, const double* b) const;

    virtual void setToSolution(escript::Data& out, escript::Data& in,
                               boost::python::object& options) const;

    virtual void ypAx(escript::Data& y, escript::Data& x) const;

    int numRows;
    std::vector<int> offsets;
    std::vector<double> values;

};

} // namespace ripley

#endif // __RIPLEY_SYSTEMMATRIX_H__

