
/*******************************************************
*
* Copyright (c) 2003-2008 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


#if !defined  escript_AbstractSystemMatrix_20040628_H
#define escript_AbstractSystemMatrix_20040628_H
#include "system_dep.h"

#include "FunctionSpace.h"
#include "SystemMatrixException.h"

#include <boost/python/dict.hpp>

//
// Forward declaration
class Data;

namespace escript {

/**
   \brief
   Give a short description of what AbstractSystemMatrix does.

   Description:
   Give a detailed description of AbstractSystemMatrix

   Template Parameters:
   For templates describe any conditions that the parameters used in the
   template must satisfy
*/
class AbstractSystemMatrix {

 public:

  /**
     \brief
     Default constructor for AbstractSystemMatrix

     Description:
     Default constructor for AbstractSystemMatrix

     Preconditions:
     Describe any preconditions

     Throws:
     Describe any exceptions thrown
  */
  ESCRIPT_DLL_API
  AbstractSystemMatrix();

  ESCRIPT_DLL_API
  AbstractSystemMatrix(const int row_blocksize,
                       const FunctionSpace& row_functionspace,
                       const int column_blocksize,
                       const FunctionSpace& column_functionspace);
  /**
    \brief
    Destructor.
  */
  ESCRIPT_DLL_API
  virtual ~AbstractSystemMatrix();

  /**
    \brief
    matrix*vector multiplication
  */
  ESCRIPT_DLL_API
  Data vectorMultiply(Data& right) const;

  /**
    \brief
    returns true if the matrix is empty
  */
  ESCRIPT_DLL_API
  int isEmpty() const;

  /**
    \brief
    returns the column function space
  */
  ESCRIPT_DLL_API
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
  ESCRIPT_DLL_API
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
  ESCRIPT_DLL_API
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
  ESCRIPT_DLL_API
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
  ESCRIPT_DLL_API
  Data solve(Data& in,const boost::python::dict& options) const;

  /**
     \brief writes the matrix to a file using the Matrix Market file format
  */
  ESCRIPT_DLL_API
  virtual void saveMM(const std::string& fileName) const;

  /**
     \brief writes the matrix to a file using the Harwell-Boeing file format
  */
  ESCRIPT_DLL_API
  virtual void saveHB(const std::string& fileName) const;

  /**
     \brief resets the matrix entries
  */
  ESCRIPT_DLL_API
  virtual void resetValues() const;

 protected:

 private:

  /**
     \brief
     solves the linear system this*out=in
  */
  ESCRIPT_DLL_API
  virtual void setToSolution(Data& out,Data& in,const boost::python::dict& options) const;

  /**
     \brief
     performs y+=this*x
  */
  ESCRIPT_DLL_API
  virtual void ypAx(Data& y,Data& x) const;

  int m_empty;
  int m_column_blocksize;
  int m_row_blocksize;
  FunctionSpace m_column_functionspace;
  FunctionSpace m_row_functionspace;

};

ESCRIPT_DLL_API Data operator*(const AbstractSystemMatrix& left,const Data& right) ;



} // end of namespace
#endif
