
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

#include "NullDomain.h" 

namespace escript {

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
  //
  // return an arbitary value
  return 1;
}
 
int NullDomain::getFunctionCode() const 
{
  //
  // return an arbitary value
  return 1;
}

int NullDomain::getFunctionOnBoundaryCode() const 
{
  //
  // return an arbitary value
  return 1;
}
 
int NullDomain::getFunctionOnContactZeroCode() const
{
  //
  // return an arbitary value
  return 1;
}

int NullDomain::getFunctionOnContactOneCode() const 
{
  //
  // return an arbitary value
  return 1;
}
 
int NullDomain::getSolutionCode() const 
{
  //
  // return an arbitary value
  return 1;
}
 
int NullDomain::getReducedSolutionCode() const
{
  //
  // return an arbitary value
  return 1;
}

int NullDomain::getDiracDeltaFunctionCode() const
{
  //
  // return an arbitary value
  return 1;
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

}  // end of namespace
