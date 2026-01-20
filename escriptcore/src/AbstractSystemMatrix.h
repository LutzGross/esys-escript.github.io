
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

#ifndef __ESCRIPT_ABSTRACTSYSTEMMATRIX_H__
#define __ESCRIPT_ABSTRACTSYSTEMMATRIX_H__

#include "system_dep.h"
#include "FunctionSpace.h"
#include "Pointers.h"
#include "SystemMatrixException.h"

#include <boost/python/object.hpp>

namespace escript {

//
// Forward declaration
class AbstractSystemMatrix;
class Data;

typedef POINTER_WRAPPER_CLASS(AbstractSystemMatrix) ASM_ptr;
typedef POINTER_WRAPPER_CLASS(const AbstractSystemMatrix) const_ASM_ptr;


/**
   \brief
   Base class for escript system matrices.
*/
class ESCRIPT_DLL_API AbstractSystemMatrix: public REFCOUNT_BASE_CLASS(AbstractSystemMatrix)
{
public:

    /**
        \brief
        Default constructor for AbstractSystemMatrix
    */
    AbstractSystemMatrix() : m_empty(true) {}

    AbstractSystemMatrix(int row_blocksize,
                         const FunctionSpace& row_functionspace,
                         int column_blocksize,
                         const FunctionSpace& column_functionspace);

    /**
        \brief
        Destructor.
    */
    virtual ~AbstractSystemMatrix() {}

    /**
        \brief Returns smart pointer which is managing this object.
        If one does not exist yet it creates one.
    */
    ASM_ptr getPtr();

    /**
        \brief Returns smart pointer which is managing this object.
        If one does not exist yet it creates one.
    */
    const_ASM_ptr getPtr() const; 

    /**
        \brief
        returns the matrix-vector product this*right
    */
    Data vectorMultiply(const Data& right) const;

    /**
        \brief
        returns true if the matrix is empty
    */
    bool isEmpty() const { return m_empty; }

    /**
        \brief
        returns the column function space
    */
    inline FunctionSpace getColumnFunctionSpace() const
    {
        if (isEmpty())
            throw SystemMatrixException("Error - Matrix is empty.");
        return m_column_functionspace;
    }

    /**
        \brief
        returns the row function space
    */
    inline FunctionSpace getRowFunctionSpace() const
    {
        if (isEmpty())
            throw SystemMatrixException("Error - Matrix is empty.");
        return m_row_functionspace;
    }

    /**
        \brief
        returns the row block size
    */
    inline int getRowBlockSize() const
    {
        if (isEmpty())
            throw SystemMatrixException("Error - Matrix is empty.");
        return m_row_blocksize;
    }

    /**
        \brief
        returns the column block size
    */
    inline int getColumnBlockSize() const
    {
        if (isEmpty())
            throw SystemMatrixException("Error - Matrix is empty.");
        return m_column_blocksize;
    }

    /**
        \brief
        returns the solution u of the linear system this*u=in
    */
    Data solve(const Data& in, boost::python::object& options) const;
  
    /**
        \brief
        sets matrix entries to zero in specified rows and columns.
        The rows and columns are marked by positive values in row_q and col_q.
        Values on the main diagonal which are marked to set to zero by both
        row_q and col_q are set to mdv (main diagonal value).
    */
    virtual void nullifyRowsAndCols(Data& row_q, Data& col_q, double mdv);  
  

    /**
        \brief writes the matrix to a file using the Matrix Market file format
    */
    virtual void saveMM(const std::string& filename) const;

    /**
        \brief writes the matrix to a file using the Harwell-Boeing file format
    */
    virtual void saveHB(const std::string& filename) const;

    /**
        \brief resets the matrix entries
    */
    virtual void resetValues(bool preserveSolverData = false);

private:

    /**
        \brief
        solves the linear system this*out=in
    */
    virtual void setToSolution(Data& out, Data& in,
                               boost::python::object& options) const;

    /**
        \brief
        performs y+=this*x
    */
    virtual void ypAx(Data& y, Data& x) const;

    bool m_empty;
    int m_column_blocksize;
    int m_row_blocksize;
    FunctionSpace m_row_functionspace;
    FunctionSpace m_column_functionspace;
};

ESCRIPT_DLL_API
Data operator*(const AbstractSystemMatrix& left, const Data& right);

} // end of namespace

#endif // __ESCRIPT_ABSTRACTSYSTEMMATRIX_H__

