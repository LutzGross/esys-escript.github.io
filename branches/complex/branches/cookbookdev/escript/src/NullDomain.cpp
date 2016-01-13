
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


#include "DomainException.h"
#include "NullDomain.h" 
#include "Data.h"

namespace escript {

namespace {
int defaultList[2]={0,1};		// an array to return in borrowListOfTagsInUse();
int NullDomainFS=1;		// Null domains only support 1 functionspace type.
			// The choice of =1 as the value is arbitrary

int referenceID=10;	// arbitrary
}


NullDomain::NullDomain() {
}
 
bool NullDomain::isValidFunctionSpaceType(int functionSpaceType) const 
{
   return (functionSpaceType==NullDomainFS);
}

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

bool NullDomain::probeInterpolationACross(int functionSpaceType_source,const AbstractDomain& targetDomain, int functionSpaceType_target) const
{
   return false;
}

int NullDomain::getContinuousFunctionCode() const 
{
  return NullDomainFS;
}
 
int NullDomain::getFunctionCode() const 
{
  return NullDomainFS;
}

int NullDomain::getFunctionOnBoundaryCode() const 
{
  return NullDomainFS;
}
 
int NullDomain::getFunctionOnContactZeroCode() const
{
  return NullDomainFS;
}

int NullDomain::getFunctionOnContactOneCode() const 
{
  return NullDomainFS;
}
 
int NullDomain::getSolutionCode() const 
{
  return NullDomainFS;
}
 
int NullDomain::getReducedSolutionCode() const
{
  return NullDomainFS;
}

int NullDomain::getDiracDeltaFunctionCode() const
{
  return NullDomainFS;
}

std::pair<int,int> NullDomain::getDataShape(int functionSpaceCode) const
{
  //
  // return an arbitary value
  // - I know it says arbitrary but its not a good idea to change it now.
  // - some tests assume that the null domain holds a single value
  return std::pair<int,int>(1,1);
}

int NullDomain::getTagFromSampleNo(int functionSpaceType, int sampleNo) const
{
  //
  // return an arbitary value
  return 1; 
}



const int* NullDomain::borrowSampleReferenceIDs(int functionSpaceType) const
{
  //
  // return an arbitary value
  return &(referenceID);
}

int NullDomain::getDim() const
{
  //
  // return an arbitary value
  return 1; 
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

bool NullDomain::operator!=(const AbstractDomain& other) const
{
  return(!(*this==other));
}



bool NullDomain::canTag(int functionSpaceCode) const
{
  return true;
}

int NullDomain::getNumberOfTagsInUse(int functionSpaceCode) const
{
  return 1;	// this is not arbitrary. It allows us to report that the default tag is in use
}

const int* NullDomain::borrowListOfTagsInUse(int functionSpaceCode) const
{
  return defaultList;
}



}  // end of namespace
