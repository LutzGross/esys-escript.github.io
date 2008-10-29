
/*******************************************************
*
* Copyright (c) 2003-2008 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


#include "NullDomain.h" 

namespace escript {

namespace {
int defaultList[1]={0};		// an array to return in borrowListOfTagsInUse();
int NullDomainFS=1;		// Null domains only support 1 functionspace type.
			// The choice of =1 as the value is arbitrary
}


NullDomain::NullDomain() {
}
 
bool NullDomain::isValidFunctionSpaceType(int functionSpaceType) const 
{
  //
  // allow anything
  return true;
}

std::string NullDomain::getDescription() const 
{
  return "NullDomain";
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
  return std::pair<int,int>(1,1);
}

int NullDomain::getTagFromSampleNo(int functionSpaceType, int sampleNo) const
{
  //
  // return an arbitary value
  return 1; 
}

int referenceID=10;

int* NullDomain::borrowSampleReferenceIDs(int functionSpaceType) const
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

int* NullDomain::borrowListOfTagsInUse(int functionSpaceCode) const
{
  return defaultList;
}



}  // end of namespace
