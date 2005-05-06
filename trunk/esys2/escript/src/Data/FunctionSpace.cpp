/* $Id$ */
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

#include "escript/Data/AbstractContinuousDomain.h"
#include "escript/Data/FunctionSpaceException.h"
#include "escript/Data/Data.h" 
#include "escript/Data/DataFactory.h" 

#include "escript/Data/FunctionSpace.h" 

#include <iostream>
#include <sstream>

using namespace std;
namespace escript {

//
// Create a null domain for use with a default constructed function space
NullDomain FunctionSpace::m_nullDomainValue;

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
	 <<" for domain: " << m_domain->getDescription();
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
cout << "toString" << m_domain->getDescription() << "\n";
  std::stringstream temp;
  temp << "Function space type: " 
       << m_domain->functionSpaceTypeAsString(m_functionSpaceType)
       << " on " << m_domain->getDescription();
  return temp.str();
}

int
FunctionSpace::getTagFromSampleNo(int sampleNo) const
{
  return m_domain->getTagFromSampleNo(m_functionSpaceType,sampleNo);
}

int
FunctionSpace::getReferenceNoFromSampleNo(int sampleNo) const
{
  return m_domain->getReferenceNoFromSampleNo(m_functionSpaceType,sampleNo);
}

FunctionSpace&
FunctionSpace::operator=(const FunctionSpace& other)
{
  //
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
  return out;
}

escript::Data
FunctionSpace::getNormal() const
{
  Data out=escript::Vector(0,*this,true);
  getDomain().setToNormal(out);
  return out;
}

escript::Data
FunctionSpace::getSize() const
{
  Data out=escript::Scalar(0,*this,true);
  getDomain().setToSize(out);
  return out;
}

}  // end of namespace
