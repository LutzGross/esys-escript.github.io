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

#include "escript/Data/AbstractSystemMatrix.h" 
#include "escript/Data/FunctionSpace.h" 
#include "escript/Data/DataException.h"
#include "escript/Data/DataArrayView.h"

namespace escript {

AbstractSystemMatrix::AbstractSystemMatrix() {
    //std::cout << "Called default AbstractSystemMatrix constructor" << std::endl;
    m_empty=1;
}

AbstractSystemMatrix::AbstractSystemMatrix(const int row_blocksize,
                                           const FunctionSpace& row_functionspace,
                                           const int column_blocksize,
                                           const FunctionSpace& column_functionspace)
{
  if (row_blocksize<=0) 
     throw DataException("Error - negative row block size of system matrix.");
  if (column_blocksize<=0) 
     throw DataException("Error - negative column block size of system matrix.");

   m_empty=0;
   m_row_blocksize=row_blocksize;
   m_column_blocksize=column_blocksize;
   m_row_functionspace=row_functionspace;
   m_column_functionspace=column_functionspace;
}

AbstractSystemMatrix::~AbstractSystemMatrix() {
}

int AbstractSystemMatrix::isEmpty() const {
   return m_empty;
}

Data operator*(const AbstractSystemMatrix& left, const Data& right)
{
      return left.vectorMultiply(right);
}

Data AbstractSystemMatrix::vectorMultiply(const Data& right) const
{
     if (isEmpty())
          throw SystemMatrixException("Error - Matrix is empty.");
     if (right.getFunctionSpace()!=getColumnFunctionSpace())
          throw SystemMatrixException("Error - column function space and function space of input data do not match.");
     if (right.getDataPointSize()!=getColumnBlockSize())
          throw SystemMatrixException("Error - column block size and input data size do not match.");
     DataArrayView::ShapeType shape;
     if (getRowBlockSize()>1) shape.push_back(getRowBlockSize());

     Data out=Data(0.,shape,getRowFunctionSpace(),true);
     ypAx(out,right);
     return out;
}

void AbstractSystemMatrix::ypAx(Data& y,const Data& x) const
{
    throw SystemMatrixException("Error - ypAx not available");
}

Data AbstractSystemMatrix::solve(const Data& in,const boost::python::dict& options) const
{
     if (isEmpty())
          throw SystemMatrixException("Error - Matrix is empty.");
     if (in.getFunctionSpace()!=getRowFunctionSpace())
          throw SystemMatrixException("Error - row function space and function space of right hand side do not match.");
     if (in.getDataPointSize()!=getRowBlockSize())
          throw SystemMatrixException("Error - row block size and right hand side size do not match.");
     DataArrayView::ShapeType shape;
     if (getRowBlockSize()>1) shape.push_back(getColumnBlockSize());
     Data out=Data(0.,shape,getColumnFunctionSpace(),true);
     setToSolution(out,in,options);
     return out;
}
void AbstractSystemMatrix::setToSolution(Data& out,const Data& in,const boost::python::dict& options) const
{
    throw SystemMatrixException("Error - setToSolution not available");
}
void AbstractSystemMatrix::saveMM(const std::string& fileName) const
{
    throw SystemMatrixException("Error - Matrix Market interface not available.");
}
void AbstractSystemMatrix:: setValue(const double value) const
{
}

}  // end of namespace
