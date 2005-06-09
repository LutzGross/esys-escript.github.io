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

#include <iostream>
#include <cassert>

#include "escript/Data/DataVector.h"
#include "escript/Data/Taipan.h"
#include "escript/Data/DataException.h"

using namespace std;
using namespace escript;

namespace escript {

Taipan arrayManager;

DataVector::DataVector() :
  m_array_data(0),
  m_size(0),
  m_dim(0),
  m_N(0)
{
}

DataVector::DataVector(const DataVector& other) :
  m_array_data(0),
  m_size(other.m_size),
  m_dim(other.m_dim),
  m_N(other.m_N)
{
  m_array_data = arrayManager.new_array(m_dim,m_N);
  int i;
  #pragma omp parallel for private(i) schedule(static)
  for (i=0; i<m_size; i++) {
    m_array_data[i] = other.m_array_data[i];
  }
}

DataVector::DataVector(const DataVector::size_type size,
                       const DataVector::value_type val,
                       const DataVector::size_type blockSize) :
  m_array_data(0),
  m_size(size),
  m_dim(blockSize)
{
  resize(size, val, blockSize);
}

DataVector::~DataVector()
{
  // dispose of data array
  arrayManager.delete_array(m_array_data);

  // clear data members
  m_size = -1;
  m_dim = -1;
  m_N = -1;
  m_array_data = 0;
}

void
DataVector::resize(const DataVector::size_type newSize,
                   const DataVector::value_type newValue,
                   const DataVector::size_type newBlockSize)
{
  assert(m_size >= 0);

  if ( (newSize % newBlockSize) != 0) {
    throw DataException("DataVector: invalid blockSize specified");
  }

  arrayManager.delete_array(m_array_data);

  m_size = newSize;
  m_dim = newBlockSize;
  m_N = newSize / newBlockSize;
  m_array_data = arrayManager.new_array(m_dim,m_N);

  int i;
  #pragma omp parallel for private(i) schedule(static)
  for (i=0; i<m_size; i++) {
    m_array_data[i] = newValue;
  }
}

DataVector&
DataVector::operator=(const DataVector& other)
{
  assert(m_size >= 0);

  arrayManager.delete_array(m_array_data);

  m_size = other.m_size;
  m_dim = other.m_dim;
  m_N = other.m_N;

  m_array_data = arrayManager.new_array(m_dim,m_N);
  int i;
  #pragma omp parallel for private(i) schedule(static)
  for (i=0; i<m_size; i++) {
    m_array_data[i] = other.m_array_data[i];
  }

  return *this;
}

bool
DataVector::operator==(const DataVector& other) const
{
  assert(m_size >= 0);

  if (m_size!=other.m_size) {
    return false;
  }
  if (m_dim!=other.m_dim) {
    return false;
  }
  if (m_N!=other.m_N) {
    return false;
  }
  for (int i=0; i<m_size; i++) {
    if (m_array_data[i] != other.m_array_data[i]) {
      return false;
    }
  }
  return true;
}

bool
DataVector::operator!=(const DataVector& other) const
{
  return !(*this==other);
}

} // end of namespace
