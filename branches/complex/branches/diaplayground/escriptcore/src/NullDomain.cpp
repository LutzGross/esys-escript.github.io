
/*****************************************************************************
*
* Copyright (c) 2003-2014 by University of Queensland
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


#include "DomainException.h"
#include "NullDomain.h" 
#include "Data.h"

namespace escript {

namespace {
int defaultList[2]={0,1}; // an array to return in borrowListOfTagsInUse();
}

// Null domains only support 1 functionspace type.
// The choice of =1 as the value is arbitrary
int NullDomain::NullDomainFS = 1;
int NullDomain::referenceID = 10; // arbitrary

std::string NullDomain::getDescription() const 
{
  return "NullDomain";
}

std::string NullDomain::functionSpaceTypeAsString(int functionSpaceType) const
{
    return "Default_FunctionSpace";
}

void NullDomain::interpolateOnDomain(Data& target,const Data& source) const
{
   if (source.getFunctionSpace().getDomain().get()!=this)  
      throw DomainException("Error - Illegal domain of interpolant.");
   if (target.getFunctionSpace().getDomain().get()!=this) 
      throw DomainException("Error - Illegal domain of interpolation target.");
   target=source;
}

bool NullDomain::probeInterpolationOnDomain(int functionSpaceType_source,int functionSpaceType_target) const
{
   if ((functionSpaceType_source!=functionSpaceType_target) || (functionSpaceType_target!=NullDomainFS))
   {
    throw DomainException("Error - Illegal function type for NullDomain.");
   }
   return true;
}

void NullDomain::interpolateACross(Data& target, const Data& source) const
{
   throw DomainException("Error - interpolation to the NullDomain not supported.");
}

std::pair<int,dim_t> NullDomain::getDataShape(int functionSpaceCode) const
{
  //
  // return an arbitrary value
  // - I know it says arbitrary but its not a good idea to change it now.
  // - some tests assume that the null domain holds a single value
  return std::pair<int,dim_t>(1,1);
}

bool NullDomain::operator==(const AbstractDomain& other) const
{
  const NullDomain* temp=dynamic_cast<const NullDomain*>(&other);
  if (temp!=0) {
    return true;
  } else {
    return false;
  }
}

const int* NullDomain::borrowListOfTagsInUse(int functionSpaceCode) const
{
  return defaultList;
}

escript::Data NullDomain::randomFill(const DataTypes::ShapeType& shape,
       const FunctionSpace& what, long seed, const boost::python::tuple& filter) const
{
    throw DataException("Attempted randomFill on NullDomain. NullDomains do not store values.");
}
void NullDomain::dump(std::string const&) const
{
    throwStandardException("NullDomain::dump");
}
void NullDomain::write(std::string const&) const
{
    throwStandardException("NullDomain::write");
}
bool NullDomain::commonFunctionSpace(std::vector<int> const&, int&) const
{
    throwStandardException("NullDomain::commonFunctionSpace");
    return false;
}
bool NullDomain::isCellOriented(int) const
{
    throwStandardException("NullDomain::isCellOriented");
    return false;
}
bool NullDomain::ownSample(int, int) const
{
    throwStandardException("NullDomain::ownSample");
    return false;
}
int NullDomain::getApproximationOrder(int) const
{
    throwStandardException("NullDomain::getApproximationOrder");
    return 0;
}
signed char NullDomain::preferredInterpolationOnDomain(int, int) const
{
    throwStandardException("NullDomain::preferredInterpolationOnDomain");
    return 0;
}
std::string NullDomain::showTagNames() const
{
    throwStandardException("NullDomain::showTagNames");
    return std::string();
}
int NullDomain::getTag(std::string const&) const
{
    throwStandardException("NullDomain::getTag");
    return 0;
}
void NullDomain::setTagMap(std::string const&, int)
{
    throwStandardException("NullDomain::setTagMap");
}
void NullDomain::setTags(int, int, escript::Data const&) const
{
    throwStandardException("NullDomain::setTags");
}
void NullDomain::setNewX(escript::Data const&)
{
    throwStandardException("NullDomain::setNewX");
}
escript::Data NullDomain::getNormal() const
{
    throwStandardException("NullDomain::getNormal");
    return escript::Data();
}
void NullDomain::setToNormal(escript::Data&) const
{
    throwStandardException("NullDomain::setToNormal");
}
void NullDomain::setToGradient(escript::Data&, escript::Data const&) const
{
    throwStandardException("NullDomain::setToGradient");
}
escript::Data NullDomain::getSize() const
{
    throwStandardException("NullDomain::getSize");
    return escript::Data();
}
void NullDomain::setToSize(escript::Data&) const
{
    throwStandardException("NullDomain::setToSize");
}
escript::Data NullDomain::getX() const
{
    throwStandardException("NullDomain::getX");
    return escript::Data();
}
void NullDomain::setToX(escript::Data&) const
{
    throwStandardException("NullDomain::setToX");
}


}  // end of namespace
