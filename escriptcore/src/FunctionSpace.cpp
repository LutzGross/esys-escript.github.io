
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/

#include "FunctionSpace.h" 

#include "Data.h" 
#include "DataFactory.h" 
#include "FunctionSpaceException.h"
#include "NullDomain.h"

#include <iostream>
#include <sstream>
#include <boost/smart_ptr.hpp>

using namespace std;
using namespace boost;

namespace escript {

  bool canInterpolate(FunctionSpace src, FunctionSpace dest)
  {
      return src.getDomain()->probeInterpolationOnDomain(src.getTypeCode(), dest.getTypeCode());
  }


namespace {
//
// Create a null domain for use with any default-constructed function space
// NullDomain const FunctionSpace::nullDomainValue;
const_Domain_ptr nullDomainValue(new NullDomain());
}

FunctionSpace::FunctionSpace():
//   m_domain(static_cast<const AbstractDomain*>(&nullDomainValue)),
  m_domain(nullDomainValue),
  m_functionSpaceType(dynamic_cast<const NullDomain*>(nullDomainValue.get())->getFunctionCode())
{
}

// FunctionSpace::FunctionSpace(const AbstractDomain& domain,
//                              int functionSpaceType):
// /*  m_domain(dynamic_cast<const AbstractDomain*>(&domain)),*/
//   m_domain(domain.getPtr()),
//   m_functionSpaceType(functionSpaceType)
// {
//   if (!m_domain->isValidFunctionSpaceType(functionSpaceType)) {
//     std::stringstream temp;
//     temp << "Invalid function space type: " << functionSpaceType 
// 	 << " for domain: " << m_domain->getDescription();
//     throw FunctionSpaceException(temp.str());
//   }
// }

FunctionSpace::FunctionSpace(const_Domain_ptr domain,
                             int functionSpaceType):
/*  m_domain(dynamic_cast<const AbstractDomain*>(&domain)),*/
  m_domain(domain),
  m_functionSpaceType(functionSpaceType)
{
  if (!m_domain->isValidFunctionSpaceType(functionSpaceType)) {
    std::stringstream temp;
    temp << "Invalid function space type: " << functionSpaceType 
	 << " for domain: " << m_domain->getDescription();
    throw FunctionSpaceException(temp.str());
  }
}

FunctionSpace::FunctionSpace(const FunctionSpace& other)
:m_domain(other.m_domain),
m_functionSpaceType(other.m_functionSpaceType)
{
}

std::pair<int,DataTypes::dim_t>
FunctionSpace::getDataShape() const
{
  return m_domain->getDataShape(m_functionSpaceType);
}

int
FunctionSpace::getTypeCode() const 
{
  return  m_functionSpaceType;
}

// const
// AbstractDomain&
const_Domain_ptr
FunctionSpace::getDomain() const
{
  return m_domain;
}

Domain_ptr
FunctionSpace::getDomainPython() const
{
  // cast away the const-ness because python ignores it anyway
  return const_pointer_cast<AbstractDomain>(m_domain);
}



std::string
FunctionSpace::toString() const
{
  std::stringstream temp;
  temp << m_domain->functionSpaceTypeAsString(m_functionSpaceType)
       << " on " << m_domain->getDescription();

  return temp.str();
}


#ifdef DEBUG_PY_STRINGS
PyObject *
FunctionSpace::toPyString() const
{
  boost::python::to_python_value<const std::string &> cvtr;
  std::stringstream temp;

  temp << m_domain->functionSpaceTypeAsString(m_functionSpaceType)
       << " on " << m_domain->getDescription();

  return cvtr(temp.str());
}
#endif


int
FunctionSpace::getTagFromSampleNo(DataTypes::dim_t sampleNo) const
{
  return m_domain->getTagFromSampleNo(m_functionSpaceType,sampleNo);
}

int
FunctionSpace::getTagFromDataPointNo(DataTypes::dim_t dataPointNo) const
{
  //
  // Get the number of samples and data-points per sample
  int numSamples = getNumSamples();
  int numDataPointsPerSample = getNumDPPSample();
  int numDataPoints = numSamples * numDataPointsPerSample;

  if (numDataPointsPerSample==0) {
    throw DataException("FunctionSpace::getTagFromDataPointNo error: no data-points associated with this object.");
  }

  if (dataPointNo<0 || dataPointNo>=numDataPoints) {
    throw DataException("FunctionSpace::getTagFromDataPointNo error: invalid data-point number supplied.");
  }

  //
  // Determine the sample number which corresponds to this data-point number
  int sampleNo = dataPointNo / numDataPointsPerSample;

  //
  // Determine the tag number which corresponds to this sample number
  int tagNo = getTagFromSampleNo(sampleNo);

  //
  // return the tag number
  return(tagNo);
}

DataTypes::dim_t FunctionSpace::getReferenceIDFromDataPointNo(DataTypes::dim_t dataPointNo) const
{
     //
     // Get the number of samples and data-points per sample
    DataTypes::dim_t numSamples = getNumSamples();
    int numDataPointsPerSample = getNumDPPSample();
    const DataTypes::dim_t* referenceIDs= borrowSampleReferenceIDs();
    DataTypes::dim_t numDataPoints = numSamples * numDataPointsPerSample;

    if (numDataPointsPerSample==0) {
        throw DataException("FunctionSpace::getReferenceIDFromDataPointNo error: no data-points associated with this object.");
    }
    if (dataPointNo<0 || dataPointNo>numDataPoints) {
        throw DataException("FunctionSpace::getReferenceIDFromDataPointNo error: invalid data-point number supplied.");
    }
    DataTypes::dim_t sampleNo = dataPointNo / numDataPointsPerSample;
    return referenceIDs[sampleNo];
}

const DataTypes::dim_t*
FunctionSpace::borrowSampleReferenceIDs() const
{
  return m_domain->borrowSampleReferenceIDs(m_functionSpaceType);
}

// FunctionSpace instances should not be overwritten to point to different domains/types
// The only time this was actually used was in constructors and the copy constructor can deal with that
FunctionSpace&
FunctionSpace::operator=(const FunctionSpace& other)
{
  throw DataException("FunctionSpace::= should not be called. Programming Error.");
  // explicitly defined assignment operator to emphasise pointer copy
/*  m_functionSpaceType=other.m_functionSpaceType;
  m_domain=other.m_domain;
  return *this;*/
}

bool
FunctionSpace::operator==(const FunctionSpace& other) const
{
  return ((*(other.m_domain)==*(m_domain)) && (other.m_functionSpaceType==m_functionSpaceType));
}

bool
FunctionSpace::operator!=(const FunctionSpace& other) const
{
  return !(operator==(other));
}

escript::Data
FunctionSpace::getX() const 
{
  Data out=escript::Vector(0,*this,true);
  getDomain()->setToX(out);
  out.setProtection();
  return out;
}

#ifdef ESYS_HAVE_BOOST_NUMPY
boost::python::numpy::ndarray 
FunctionSpace::getNumpyX() const
{
    // Initialise boost numpy
    boost::python::numpy::initialize();

    // Get the data
    Data data=escript::Vector(0,*this,true);
    getDomain()->setToX(data);
    data.setProtection();

    // This is needed below in getSampleDataRW
    const DataTypes::real_t onlyreal = 0;

    // Work out how many data points there are
    int numDataPoints = data.getNumSamples();
    int dpps = data.getNumDataPointsPerSample();

    // Work out the data point shape
    std::vector<int> shape = data.getDataPointShape();
    
    // Work out how many spatial dimensions there are
    int dimensions = data.getShapeProduct();

    // Initialise the ndarray
    boost::python::tuple arrayshape = boost::python::make_tuple(dimensions, dpps * numDataPoints);
    boost::python::numpy::dtype datatype = boost::python::numpy::dtype::get_builtin<double>();
    boost::python::numpy::ndarray dataArray = boost::python::numpy::zeros(arrayshape, datatype);

    // Initialise variables
    std::string localmsg;
    std::vector<const DataTypes::real_t*> samplesR(1);
    
// #pragma omp parallel for 
    for (int i = 0; i < numDataPoints; ++i) {
        for (int j = 0; j < shape[0]; j++) {
            dataArray[j][i] = *(data.getSampleDataRW(i, onlyreal)+j);
        }
    }

    // Print out the ndarray to the console - used during debugging 
    // std::cout << "Finished array:\n" << bp::extract<char const *>(bp::str(dataArray)) << std::endl;

    return dataArray;
}
#endif

escript::Data
FunctionSpace::getNormal() const
{
  Data out=escript::Vector(0,*this,true);
  getDomain()->setToNormal(out);
  out.setProtection();
  return out;
}

escript::Data
FunctionSpace::getSize() const
{
  Data out=escript::Scalar(0,*this,true);
  getDomain()->setToSize(out);
  out.setProtection();
  return out;
}

void
FunctionSpace::setTags(const int newTag, const escript::Data& mask) const
{
   if (mask.getFunctionSpace()== *this) {
          m_domain->setTags(m_functionSpaceType,newTag,mask);
   } else {
          throw FunctionSpaceException("illegal function space of mask.");
   }
}

void
FunctionSpace::setTagsByString(const std::string& name, const escript::Data& mask) const
{
   int newTag=m_domain->getTag(name);
   if (mask.getFunctionSpace()== *this) {
          m_domain->setTags(m_functionSpaceType,newTag,mask);
   } else {
          throw FunctionSpaceException("illegal function space of mask.");
   }
}

int 
FunctionSpace::getNumberOfTagsInUse() const
{
   return  m_domain->getNumberOfTagsInUse(m_functionSpaceType);
}

const int* FunctionSpace::borrowListOfTagsInUse() const
{
   return m_domain->borrowListOfTagsInUse(m_functionSpaceType);
}

std::list<int> FunctionSpace::getListOfTagsSTL() const 
{
  const int* tags=borrowListOfTagsInUse();
  std::list<int> taglist(tags, tags+getNumberOfTagsInUse());
  return taglist;
}


boost::python::list FunctionSpace::getListOfTags() const 
{
  const int* tags=borrowListOfTagsInUse();
  boost::python::list taglist;
  for (int i=0; i<getNumberOfTagsInUse(); ++i)
      taglist.append(tags[i]);
  return taglist;
}

bool
FunctionSpace::canTag() const
{
  return m_domain->canTag(m_functionSpaceType);
}

int 
FunctionSpace::getApproximationOrder() const
{
   return m_domain->getApproximationOrder(m_functionSpaceType);
}

}  // end of namespace

