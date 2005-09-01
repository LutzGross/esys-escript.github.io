// $Id$
/*
 ******************************************************************************
 *                                                                            *
 *       COPYRIGHT  ACcESS 2005 -  All Rights Reserved                        *
 *                                                                            *
 * This software is the property of ACcESS. No part of this code              *
 * may be copied in any form or by any means without the expressed written    *
 * consent of ACcESS.  Copying, use or modification of this software          *
 * by any unauthorised person is illegal unless that person has a software    *
 * license agreement with ACcESS.                                             *
 *                                                                            *
 ******************************************************************************
*/

#include "bruce/Bruce/Bruce.h"

using namespace std;
using namespace escript;

namespace bruce {

const int Bruce::Nodes=0;
const int Bruce::Elements=1;

Bruce::Bruce()
{
}

Bruce::Bruce(const Bruce& other)
{
}

Bruce::~Bruce()
{
}

bool
Bruce::isValidFunctionSpaceType(int functionSpaceType) const
{
  return (true);
}

int
Bruce::getDim() const
{
  return 0;
}

pair<int,int>
Bruce::getDataShape(int functionSpaceCode) const
{
  int numDataPointsPerSample=0;
  int numSamples=0;
  return pair<int,int>(numDataPointsPerSample,numSamples);
}

bool
Bruce::operator==(const AbstractDomain& other) const
{
  const Bruce* temp=dynamic_cast<const Bruce*>(&other);
  if (temp!=0) {
    return (true);
  } else {
    return false;
  }
}

bool
Bruce::operator!=(const AbstractDomain& other) const
{
  return !(operator==(other));
}

Data
Bruce::getX() const
{
  return continuousFunction(asAbstractContinuousDomain()).getX();
}

Data
Bruce::getSize() const
{
  return function(asAbstractContinuousDomain()).getSize();
}

}  // end of namespace
