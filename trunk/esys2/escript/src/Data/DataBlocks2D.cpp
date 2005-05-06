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

#include "escript/Data/DataException.h"
#include "escript/Data/DataBlocks2D.h" 
#include "esysUtils/EsysAssert.h"

#include <sstream>
#include <iostream>

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
}

void
DataBlocks2D::resize(int numRows, int numCols, int blockSize)
{
    if (numRows < 0 || numCols < 0 || blockSize < 0) {
      stringstream mess;
      mess << "Error - Invalid resize parameter. numRows: " << numRows
	   << " numCols: " << numCols << " blockSize: " << blockSize;
      throw DataException(mess.str());
    }
    ValueType::size_type size=numRows*numCols*blockSize;
    m_data.resize(size);
    m_numRows=numRows;
    m_numCols=numCols;
    m_blockSize=blockSize;
}

void
DataBlocks2D::Swap(DataBlocks2D& other)
{
    swap(m_data,other.m_data);
    swap(m_blockSize,other.m_blockSize);
    swap(m_numRows,other.m_numRows);
    swap(m_numCols,other.m_numCols);
}

DataBlocks2D&
DataBlocks2D::operator=(const DataBlocks2D& other)
{
    DataBlocks2D temp(other);
    Swap(temp);
    return *this;
}
 
}  // end of namespace
