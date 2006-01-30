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

#include "AbstractDomain.h" 
#include "DomainException.h"
#include "AbstractSystemMatrix.h"
#include "FunctionSpace.h"

using namespace std;

namespace escript {

AbstractDomain::AbstractDomain() {
}

AbstractDomain::~AbstractDomain() {
}

void AbstractDomain::throwStandardException(const std::string& functionName) const
{
  throw DomainException("Error - Base class function: " + functionName + " should not be called. Programming error.");
}

bool AbstractDomain::isValidFunctionSpaceType(int functionSpaceType) const
{
  throwStandardException("AbstractDomain::isValidFunctionSpaceType");
  return false;
}

std::string AbstractDomain::getDescription() const
{
  throwStandardException("AbstractDomain::getDescription");
  return "";
}

std::string AbstractDomain::functionSpaceTypeAsString(int functionSpaceType) const
{
  throwStandardException("AbstractDomain::functionSpaceTypeAsString");
  return "";
}

int AbstractDomain::getDim() const
{
  throwStandardException("AbstractDomain::getDim");
  return 0;
}

void AbstractDomain::write(const std::string& filename) const
{
  throwStandardException("AbstractDomain::write");
  return;
}

void AbstractDomain::getTagList(int functionSpaceType, int** tagList, int* numTags) const
{
  throwStandardException("AbstractDomain::getTagList");
  return;
}

void AbstractDomain::getReferenceNoList(int functionSpaceType, int** referenceNoList, int* numReferenceNo) const
{
  throwStandardException("AbstractDomain::getTagList");
  return;
}

std::pair<int,int> AbstractDomain::getDataShape(int functionSpaceCode) const
{
  throwStandardException("AbstractDomain::getDataShape");
  return std::pair<int,int>(0,0);
}

int AbstractDomain::getTagFromSampleNo(int functionSpaceType, int sampleNo) const
{
  throwStandardException("AbstractDomain::getTagFromSampleNo");
  return 0;
}

int AbstractDomain::getReferenceNoFromSampleNo(int functionSpaceType, int sampleNo) const
{
  throwStandardException("AbstractDomain::getReferenceNoFromSampleNo");
  return 0;
}

void AbstractDomain::setNewX(const escript::Data& arg)
{
  throwStandardException("AbstractDomain::setNewX");
  return;
}

void AbstractDomain::interpolateOnDomain(escript::Data& target,const escript::Data& source) const
{
  throwStandardException("AbstractDomain::interpolateOnDomain");
  return;
}
void AbstractDomain::interpolateACross(escript::Data& target, const escript::Data& source) const
{
  throwStandardException("AbstractDomain::interpolateACross");
  return;
}

escript::Data AbstractDomain::getX() const
{
  throwStandardException("AbstractDomain::getX");
  return Data();
}

escript::Data AbstractDomain::getNormal() const
{
  throwStandardException("AbstractDomain::getNormal");
  return Data();
}

escript::Data AbstractDomain::getSize() const
{
  throwStandardException("AbstractDomain::getSize");
  return Data();
}

void AbstractDomain::setToX(escript::Data& out) const
{
  throwStandardException("AbstractDomain::setToX");
  return;
}

void AbstractDomain::setToNormal(escript::Data& out) const
{
  throwStandardException("AbstractDomain::setToNormal");
  return;
}

void AbstractDomain::setToSize(escript::Data& out) const
{
  throwStandardException("AbstractDomain::setToSize");
  return;
}

void AbstractDomain::setToGradient(escript::Data& grad, const escript::Data& arg) const
{
  throwStandardException("AbstractDomain::setToGradient");
  return;
}

void AbstractDomain::saveDX(const std::string& filename,const boost::python::dict& arg) const 
{
  throwStandardException("AbstractDomain::saveDX");
  return;
}

void AbstractDomain::saveVTK(const std::string& filename,const boost::python::dict& arg) const 
{
  throwStandardException("AbstractDomain::saveVTK");
  return;
}

bool AbstractDomain::probeInterpolationOnDomain(int functionSpaceType_source,int functionSpaceType_target) const
{
  throwStandardException("AbstractDomain::probeInterpolationOnDomain");
  return false;
}

bool AbstractDomain::probeInterpolationACross(int functionSpaceType_source,const AbstractDomain& targetDomain, int functionSpaceType_target) const 
{
  throwStandardException("AbstractDomain::probeInterpolationACross");
  return false;
}

bool AbstractDomain::isCellOriented(int functionSpaceCode) const
{
  throwStandardException("AbstractDomain::isCellOriented");
  return false;
}

bool AbstractDomain::operator==(const AbstractDomain& other) const
{
  throwStandardException("AbstractDomain::operator==");
  return false;
}
bool AbstractDomain::operator!=(const AbstractDomain& other) const
{
  throwStandardException("AbstractDomain::operator!=");
  return false;
}

} // end of namespace
