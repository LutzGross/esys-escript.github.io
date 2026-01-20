
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

#include "AbstractSystemMatrix.h" 
#include "Data.h"
#include "DataException.h"
#include "DataTypes.h"

namespace escript {

AbstractSystemMatrix::AbstractSystemMatrix(int row_blocksize,
                                           const FunctionSpace& row_fs,
                                           int column_blocksize,
                                           const FunctionSpace& column_fs) :
    m_empty(false),
    m_column_blocksize(column_blocksize),
    m_row_blocksize(row_blocksize),
    m_row_functionspace(row_fs),
    m_column_functionspace(column_fs)
{
    if (row_blocksize <= 0) 
        throw DataException("Negative row block size of system matrix.");
    if (column_blocksize <= 0) 
        throw DataException("Negative column block size of system matrix.");

}

ASM_ptr AbstractSystemMatrix::getPtr()
{
    try {
        return shared_from_this();
    } catch (boost::bad_weak_ptr* p) {
        return ASM_ptr(this);
    }
}

const_ASM_ptr AbstractSystemMatrix::getPtr() const 
{  
    try {
        return shared_from_this();
    } catch (boost::bad_weak_ptr* p) {
        return const_ASM_ptr(this);
    }
}

Data operator*(const AbstractSystemMatrix& left, const Data& right)
{
    return left.vectorMultiply(right);
}

Data AbstractSystemMatrix::vectorMultiply(const Data& right) const
{
    if (isEmpty())
        throw SystemMatrixException("Error - Matrix is empty.");
    if (right.getDataPointSize()!=getColumnBlockSize())
        throw SystemMatrixException("Error - column block size and input data size do not match.");
    DataTypes::ShapeType shape;
    if (getRowBlockSize() > 1)
        shape.push_back(getRowBlockSize());

    Data out = right.isComplex() ?
        Data(DataTypes::cplx_t(0), shape, getRowFunctionSpace(), true) :
        Data(0., shape, getRowFunctionSpace(), true);

    Data in(right, getColumnFunctionSpace());
    ypAx(out, in);
    return out;
}

void AbstractSystemMatrix::ypAx(Data& y, Data& x) const
{
    throw SystemMatrixException("ypAx() is not implemented.");
}

Data AbstractSystemMatrix::solve(const Data& in,
                                 boost::python::object& options) const
{
    if (isEmpty())
        throw SystemMatrixException("Matrix is empty.");
    if (in.getFunctionSpace() != getRowFunctionSpace())
        throw SystemMatrixException("row function space and function space of right hand side do not match.");
    if (in.getDataPointSize() != getRowBlockSize())
        throw SystemMatrixException("row block size and right hand side size do not match.");
    DataTypes::ShapeType shape;
    if (getRowBlockSize() > 1)
        shape.push_back(getColumnBlockSize());
    Data out = in.isComplex() ?
        Data(DataTypes::cplx_t(0), shape, getColumnFunctionSpace(), true) :
        Data(0., shape, getColumnFunctionSpace(), true);
    setToSolution(out, *const_cast<Data*>(&in), options);
    return out;
}
void AbstractSystemMatrix::setToSolution(Data& out, Data& in,
                                         boost::python::object& options) const
{
    throw SystemMatrixException("setToSolution() is not implemented");
}

void AbstractSystemMatrix::nullifyRowsAndCols(Data& row_q,
                                              Data& col_q,
                                              double mdv)
{
    throw SystemMatrixException("nullifyRowsAndCols() is not implemented.");
}

void AbstractSystemMatrix::saveMM(const std::string& filename) const
{
    throw SystemMatrixException("Matrix Market interface not available.");
}

void AbstractSystemMatrix::saveHB(const std::string& fileName) const
{
    throw SystemMatrixException("Harwell-Boeing interface not available.");
}

void AbstractSystemMatrix::resetValues(bool preserveSolverData)
{
    throw SystemMatrixException("resetValues() is not implemented.");
}

}  // end of namespace

