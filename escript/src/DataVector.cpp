
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

#include "DataVector.h"

#include "Taipan.h"
#include "DataException.h"
#include <boost/python/extract.hpp>
#include "DataTypes.h"

#include <cassert>

using namespace std;
using namespace escript;
using namespace boost::python;

namespace escript {

Taipan arrayManager;

void releaseUnusedMemory()
{
   arrayManager.release_unused_arrays();
}


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
  if (m_array_data!=0) {
    arrayManager.delete_array(m_array_data);
  }

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

  if ( newBlockSize == 0) {
    throw DataException("DataVector: invalid blockSize specified (newBlockSize)");
  }

  if ( (newSize % newBlockSize) != 0) {
    throw DataException("DataVector: invalid blockSize specified");
  }

  if (m_array_data!=0) {
    arrayManager.delete_array(m_array_data);
  }

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

  if (m_array_data!=0) {
    arrayManager.delete_array(m_array_data);
  }

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

int
DataVector::archiveData(ofstream& archiveFile,
                        const size_type noValues) const
{
  //
  // Check number of values expected to be written matches number in this object
  if (noValues != size()) {
    return 2;
  }

  //
  // Write all values in this object out to archiveFile
  for (int i=0; i<size(); i++) {
    archiveFile.write(reinterpret_cast<char *>(&m_array_data[i]),sizeof(double));
  }

  //
  // Check no errors were encountered before returning
  if (!archiveFile.good()) {
    return 1;
  }

  return 0;
}

int
DataVector::extractData(ifstream& archiveFile,
                        const size_type noValues)
{
  //
  // Check number of values expected to be read matches number in this object
  if (noValues != size()) {
    return 2;
  }

  //
  // Read all values in archiveFile back to this object
  for (int i=0; i<size(); i++) {
    archiveFile.read(reinterpret_cast<char *>(&m_array_data[i]),sizeof(double));
  }

  //
  // Check no errors were encountered before returning
  if (!archiveFile.good()) {
    return 1;
  }

  return 0;
}


void
DataVector::copyFromNumArray(const boost::python::numeric::array& value)
{
  using DataTypes::ValueType;
  if (m_array_data!=0) {
    arrayManager.delete_array(m_array_data);
  }


  m_array_data = arrayManager.new_array(1,value.nelements());

      int si=0,sj=0,sk=0,sl=0;		// bounds for each dimension of the shape
      DataTypes::ShapeType tempShape;    
      for (int i=0; i<value.getrank(); i++) {
         tempShape.push_back(extract<int>(value.getshape()[i]));
      }

      if (value.getrank()==0) {
	m_array_data[0]=extract<double>(value[value.getshape()]);
      } else if (value.getrank()==1) {
	 si=tempShape[0];
         for (ValueType::size_type i=0;i<si;i++) {
            m_array_data[i]=extract<double>(value[i]);
         }
      } else if (value.getrank()==2) {
	 si=tempShape[0];
	 sj=tempShape[1];
         for (ValueType::size_type i=0;i<si;i++) {
            for (ValueType::size_type j=0;j<sj;j++) {
               m_array_data[DataTypes::getRelIndex(tempShape,i,j)]=extract<double>(value[i][j]);
            }
         }
      } else if (value.getrank()==3) {
	 si=tempShape[0];
	 sj=tempShape[1];
	 sk=tempShape[2];
         for (ValueType::size_type i=0;i<si;i++) {
            for (ValueType::size_type j=0;j<sj;j++) {
               for (ValueType::size_type k=0;k<sk;k++) {
                  m_array_data[DataTypes::getRelIndex(tempShape,i,j,k)]=extract<double>(value[i][j][k]);
               }
            }
         }
      } else if (value.getrank()==4) {
	 si=tempShape[0];
	 sj=tempShape[1];
	 sk=tempShape[2];
	 sl=tempShape[3];
         for (ValueType::size_type i=0;i<si;i++) {
            for (ValueType::size_type j=0;j<sj;j++) {
               for (ValueType::size_type k=0;k<sk;k++) {
                  for (ValueType::size_type l=0;l<sl;l++) {
                     m_array_data[DataTypes::getRelIndex(tempShape,i,j,k,l)]=extract<double>(value[i][j][k][l]);
                  }
               }
            }
         }
      }
   }
 


} // end of namespace
