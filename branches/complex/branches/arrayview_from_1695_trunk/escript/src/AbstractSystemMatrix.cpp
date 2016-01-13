
/* $Id$ */

/*******************************************************
 *
 *           Copyright 2003-2007 by ACceSS MNRF
 *       Copyright 2007 by University of Queensland
 *
 *                http://esscc.uq.edu.au
 *        Primary Business: Queensland, Australia
 *  Licensed under the Open Software License version 3.0
 *     http://www.opensource.org/licenses/osl-3.0.php
 *
 *******************************************************/

#include "AbstractSystemMatrix.h" 
#include "DataException.h"
#include "Data.h"
#include "DataTypes.h"

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

Data operator*(const AbstractSystemMatrix& left,const Data& right)
{
      Data tmp=(Data) right;
      return left.vectorMultiply(tmp);
}

Data AbstractSystemMatrix::vectorMultiply(Data& right) const
{
     if (isEmpty())
          throw SystemMatrixException("Error - Matrix is empty.");
     if (right.getDataPointSize()!=getColumnBlockSize())
          throw SystemMatrixException("Error - column block size and input data size do not match.");
     DataTypes::ShapeType shape;
     if (getRowBlockSize()>1) shape.push_back(getRowBlockSize());

     Data out=Data(0.,shape,getRowFunctionSpace(),true);
     Data in=Data(right,getColumnFunctionSpace());
     ypAx(out,in);
     return out;
}

void AbstractSystemMatrix::ypAx(Data& y,Data& x) const
{
    throw SystemMatrixException("Error - ypAx not available");
}

Data AbstractSystemMatrix::solve(Data& in,const boost::python::dict& options) const
{
     if (isEmpty())
          throw SystemMatrixException("Error - Matrix is empty.");
     if (in.getFunctionSpace()!=getRowFunctionSpace())
          throw SystemMatrixException("Error - row function space and function space of right hand side do not match.");
     if (in.getDataPointSize()!=getRowBlockSize())
          throw SystemMatrixException("Error - row block size and right hand side size do not match.");
     DataTypes::ShapeType shape;
     if (getRowBlockSize()>1) shape.push_back(getColumnBlockSize());
     Data out=Data(0.,shape,getColumnFunctionSpace(),true);
     setToSolution(out,in,options);
     return out;
}
void AbstractSystemMatrix::setToSolution(Data& out,Data& in,const boost::python::dict& options) const
{
    throw SystemMatrixException("Error - setToSolution not available");
}
void AbstractSystemMatrix::saveMM(const std::string& fileName) const
{
    throw SystemMatrixException("Error - Matrix Market interface not available.");
}
void AbstractSystemMatrix::saveHB(const std::string& fileName) const
{
    throw SystemMatrixException("Error - Harwell-Boeing interface not available.");
}
void AbstractSystemMatrix::resetValues() const
{
    throw SystemMatrixException("Error - setValue is not implemented.");
}

}  // end of namespace
