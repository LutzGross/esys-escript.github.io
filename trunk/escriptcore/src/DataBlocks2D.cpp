
/*****************************************************************************
*
* Copyright (c) 2003-2015 by The University of Queensland
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


#include "DataBlocks2D.h"

#include "DataException.h"
#include "esysUtils/EsysAssert.h"

using namespace std;

namespace escript {

DataBlocks2D::DataBlocks2D():
  m_numRows(0),
  m_numCols(0),
  m_blockSize(0)
{
}

DataBlocks2D::DataBlocks2D(const DataBlocks2D& other):
  m_numRows(other.m_numRows),
  m_numCols(other.m_numCols),
  m_blockSize(other.m_blockSize)
{
    m_data=other.m_data;
}

DataBlocks2D::DataBlocks2D(int numRows, int numCols, int blockSize):
  m_numRows(numRows),
  m_numCols(numCols),
  m_blockSize(blockSize)
{
    resize(m_numRows,numCols,blockSize);
}

DataBlocks2D::~DataBlocks2D()
{
    m_numRows=-1;
    m_numCols=-1;
    m_blockSize=-1;
}

void
DataBlocks2D::resize(int numRows, int numCols, int blockSize)
{
    if (numRows < 1 || numCols < 1 || blockSize < 1) {
      stringstream temp;
      temp << "DataBlocks2D: Error - Invalid resize parameter. numRows: " << numRows
	   << " numCols: " << numCols << " blockSize: " << blockSize;
      throw DataException(temp.str());
    }
    ValueType::size_type size=numRows*numCols*blockSize;
    m_data.resize(size, 0.0, numCols*blockSize);
    m_numRows=numRows;
    m_numCols=numCols;
    m_blockSize=blockSize;
}

void
DataBlocks2D::Swap(DataBlocks2D& other)
{
    EsysAssert(((m_numRows >= 0) && (m_numCols >= 0) && (m_blockSize >= 0)), "(DataBlocks2D) Invalid object.");
    swap(m_data,other.m_data);
    swap(m_blockSize,other.m_blockSize);
    swap(m_numRows,other.m_numRows);
    swap(m_numCols,other.m_numCols);
}

DataBlocks2D&
DataBlocks2D::operator=(const DataBlocks2D& other)
{
    EsysAssert(((m_numRows >= 0) && (m_numCols >= 0) && (m_blockSize >= 0)), "(DataBlocks2D) Invalid object.");
    DataBlocks2D temp(other);
    Swap(temp);
    return *this;
}
 
}  // end of namespace
