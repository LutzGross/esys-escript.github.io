
/*******************************************************
*
* Copyright (c) 2003-2009 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


#include "AbstractDomain.h" 
#include "DomainException.h"
#include "Data.h"
#include "paso/Paso_MPI.h"

using namespace std;

namespace escript {

// Please see the discusion in DataAbstract's version of this method
// and why I squash the exception
Domain_ptr AbstractDomain::getPtr()
{
  try
  {
      return shared_from_this();
  }
  catch (boost::bad_weak_ptr p)		
  {
      return Domain_ptr(this);
  }
}

const_Domain_ptr AbstractDomain::getPtr() const 
{  
  try
  {
      return shared_from_this();
  }
  catch (boost::bad_weak_ptr p)		
  {
      return const_Domain_ptr(this);
  }
}

AbstractDomain::AbstractDomain() {
}

AbstractDomain::~AbstractDomain() {
}

int AbstractDomain::getMPISize() const
{
   return 1;
}
int AbstractDomain::getMPIRank() const
{
   return 0;
}
void AbstractDomain::MPIBarrier() const
{
   throwStandardException("MPIBarrier");
   return;
}
bool AbstractDomain::onMasterProcessor() const
{
   return true;
}

ESCRIPT_DLL_API
#ifdef PASO_MPI
  MPI_Comm
#else
  unsigned int
#endif
AbstractDomain::getMPIComm() const
{
  throwStandardException("AbstractDomain::getMPIComm");
  return false;
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
void AbstractDomain::dump(const std::string& filename) const
{
  throwStandardException("AbstractDomain::dump");
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

const int* AbstractDomain::borrowSampleReferenceIDs(int functionSpaceType) const
{
  throwStandardException("AbstractDomain::borrowSampleReferenceIDs");
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

void AbstractDomain::setTags(const int functionSpaceType, const int newTag, const escript::Data& mask) const
{
  throwStandardException("AbstractDomain::setTags");
  return;
}

void AbstractDomain::saveDX(const std::string& filename,const boost::python::dict& arg) const 
{
  throwStandardException("AbstractDomain::saveDX");
  return;
}

void AbstractDomain::saveVTK(const std::string& filename,const boost::python::dict& arg, const std::string& metadata, const std::string& metadata_schema) const 
{
  throwStandardException("AbstractDomain::saveVTK");
  return;
}

bool AbstractDomain::probeInterpolationOnDomain(int functionSpaceType_source,int functionSpaceType_target) const
{
  throwStandardException("AbstractDomain::probeInterpolationOnDomain");
  return false;
}

bool AbstractDomain::commonFunctionSpace(const std::vector<int>& fs, int& resultcode) const
{
  throwStandardException("AbstractDomain::commonFunctionSpace");
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

AbstractDomain::StatusType AbstractDomain::getStatus() const 
{
  return 0;
}
void AbstractDomain::setTagMap(const std::string& name,  int tag)
{
  throwStandardException("AbstractDomain::set TagMap is not implemented.");
}
int AbstractDomain::getTag(const std::string& name) const
{
  throwStandardException("AbstractDomain::getTag is not implemented.");
  return 0;
}

bool AbstractDomain::isValidTagName(const std::string& name) const
{
  return false;
}

std::string AbstractDomain::showTagNames() const
{
  throwStandardException("AbstractDomain::showTagNames is not implemented.");
  return string();
}

int AbstractDomain::getNumberOfTagsInUse(int functionSpaceCode) const
{
  throwStandardException("AbstractDomain::getNumberOfTagsInUse is not implemented.");
  return 0;
}

const int* AbstractDomain::borrowListOfTagsInUse(int functionSpaceCode) const
{
  throwStandardException("AbstractDomain::borrowListOfTagsInUse is not implemented.");
  return NULL;
}


bool AbstractDomain::canTag(int functionspacecode) const
{
  throwStandardException("AbstractDomain::canTag is not implemented.");
  return false;
}


} // end of namespace
