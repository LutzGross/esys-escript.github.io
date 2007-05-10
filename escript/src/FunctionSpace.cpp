/* $Id$ */
/*
 ************************************************************
 *          Copyright 2006 by ACcESS MNRF                   *
 *                                                          *
 *              http://www.access.edu.au                    *
 *       Primary Business: Queensland, Australia            *
 *  Licensed under the Open Software License version 3.0    *
 *     http://www.opensource.org/licenses/osl-3.0.php       *
 *                                                          *
 ************************************************************
*/

#include "FunctionSpace.h" 
#include "FunctionSpaceException.h"
#include "Data.h" 
#include "DataFactory.h" 

#include <iostream>
#include <sstream>

using namespace std;

namespace escript {

//
// Create a null domain for use with any default-constructed function space
ESCRIPT_DLL_API NullDomain FunctionSpace::m_nullDomainValue;

FunctionSpace::FunctionSpace():
  m_domain(static_cast<AbstractDomain*>(&m_nullDomainValue)),
  m_functionSpaceType(m_nullDomainValue.getFunctionCode())
{
}

FunctionSpace::FunctionSpace(const AbstractDomain& domain,
                             int functionSpaceType):
  m_domain(dynamic_cast<const AbstractDomain*>(&domain)),
  m_functionSpaceType(functionSpaceType)
{
  if (!m_domain->isValidFunctionSpaceType(functionSpaceType)) {
    std::stringstream temp;
    temp << "Invalid function space type: " << functionSpaceType 
	 << " for domain: " << m_domain->getDescription();
    throw FunctionSpaceException(temp.str());
  }
}

std::pair<int,int>
FunctionSpace::getDataShape() const
{
  return m_domain->getDataShape(m_functionSpaceType);
}

int
FunctionSpace::getTypeCode() const 
{
  return  m_functionSpaceType;
}

const
AbstractDomain&
FunctionSpace::getDomain() const
{
  return *m_domain;
}

std::string
FunctionSpace::toString() const
{
  std::stringstream temp;
  temp << m_domain->functionSpaceTypeAsString(m_functionSpaceType) << " on " << m_domain->getDescription();
  return temp.str();
}

const
boost::python::str
FunctionSpace::str() const
{
  return boost::python::str(toString().c_str());
}


int
FunctionSpace::getTagFromSampleNo(int sampleNo) const
{
  return m_domain->getTagFromSampleNo(m_functionSpaceType,sampleNo);
}

int
FunctionSpace::getTagFromDataPointNo(int dataPointNo) const
{
  //
  // Get the number of samples and data-points per sample
  int numSamples = getNumSamples();
  int numDataPointsPerSample = getNumDPPSample();
  int numDataPoints = numSamples * numDataPointsPerSample;

  if (numDataPointsPerSample==0) {
    throw DataException("FunctionSpace::getTagFromDataPointNo error: no data-points associated with this object.");
  }

  if (dataPointNo<0 || dataPointNo>numDataPoints) {
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

int*
FunctionSpace::borrowSampleReferenceIDs() const
{
  return m_domain->borrowSampleReferenceIDs(m_functionSpaceType);
}

FunctionSpace&
FunctionSpace::operator=(const FunctionSpace& other)
{
  // explicitly defined assignment operator to emphasise pointer copy
  m_nullDomainValue=other.m_nullDomainValue;
  m_functionSpaceType=other.m_functionSpaceType;
  m_domain=other.m_domain;
  return *this;
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
  getDomain().setToX(out);
  out.setProtection();
  return out;
}

escript::Data
FunctionSpace::getNormal() const
{
  Data out=escript::Vector(0,*this,true);
  getDomain().setToNormal(out);
  out.setProtection();
  return out;
}

escript::Data
FunctionSpace::getSize() const
{
  Data out=escript::Scalar(0,*this,true);
  getDomain().setToSize(out);
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

}  // end of namespace
