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

#include "escript/Data/NullDomain.h" 

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

int NullDomain::getReferenceNoFromSampleNo(int functionSpaceType, int sampleNo) const
{
  //
  // return an arbitary value
  return 1; 
}

int NullDomain::getDim() const
{
  //
  // return an arbitary value
  return 1; 
}

}  // end of namespace
