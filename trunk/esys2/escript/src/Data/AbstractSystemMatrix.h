// $Id$
/* 
 ******************************************************************************
 *                                                                            *
 *       COPYRIGHT  ACcESS 2004 -  All Rights Reserved                        *
 *                                                                            *
 * This software is the property of ACcESS. No part of this code              *
 * may be copied in any form or by any means without the expressed written    *
 * consent of ACcESS.  Copying, use or modification of this software          *
 * by any unauthorised person is illegal unless that person has a software    *
 * license agreement with ACcESS.                                             *
 *                                                                            *
 ******************************************************************************
*/
                                                                           
#if !defined  escript_AbstractSystemMatrix_20040628_H
#define escript_AbstractSystemMatrix_20040628_H

#include "escript/Data/FunctionSpace.h"
#include "escript/Data/SystemMatrixException.h"
#include "escript/Data/Data.h"
#include <boost/python/dict.hpp>

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
  AbstractSystemMatrix();

  AbstractSystemMatrix(const int row_blocksize,
                       const FunctionSpace& row_functionspace,
                       const int column_blocksize,
                       const FunctionSpace& column_functionspace);
  /**
    \brief
    Destructor.
  */
  virtual ~AbstractSystemMatrix();

  /**
    \brief
    matrix*vector multiplication
  */
  Data vectorMultiply(const Data& right) const;

  /**
    \brief
    returns true if the matrix is empty
  */
  int isEmpty() const;

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
  Data solve(const Data& in,const boost::python::dict& options) const;

  /**
     \brief writes the matrix to a file using the Matrix Market file format
  */
  virtual void saveMM(const std::string& fileName) const;
                                                                                                                                                     
  /**
     \brief sets the matrix entries to value
  */
  virtual void setValue(const double value) const;

  /**
     \brief cleans any setting, allocations by the solver. 
  */
  virtual void resetSolver() const;


 protected:

 private:

  /**
     \brief
     solves the linear system this*out=in
  */
  virtual void setToSolution(Data& out,const Data& in,const boost::python::dict& options) const;

  /**
     \brief
     performs y+=this*x
  */
  virtual void ypAx(Data& y,const Data& x) const;

  int m_empty;
  int m_column_blocksize;
  int m_row_blocksize;
  FunctionSpace m_column_functionspace;
  FunctionSpace m_row_functionspace;

};

Data operator*(const AbstractSystemMatrix& left, const Data& right) ;



} // end of namespace
#endif
