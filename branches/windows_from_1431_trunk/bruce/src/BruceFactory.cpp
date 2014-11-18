
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

#include "Bruce.h"

#include "BruceException.h"
#include "BruceFactory.h"

#include <boost/python/extract.hpp>

#include <iostream>
#include <sstream>

using namespace std;
using namespace escript;

namespace bruce {

AbstractContinuousDomain* brick(int n0, int n1, int n2,
                                double l0, double l1, double l2)
{
  double v0_len, v1_len, v2_len;

  Bruce::DimVec v0;
  Bruce::DimVec v1;
  Bruce::DimVec v2;
  Bruce::DimVec origin;

  origin.push_back(0);
  origin.push_back(0);
  origin.push_back(0);

  if (n0<2) {
    v0_len = 0;
  } else {
    v0_len = l0/(n0-1);
  }

  if (n1<2) {
    v1_len = 0;
  } else {
    v1_len = l1/(n1-1);
  }

  if (n2<2) {
    v2_len = 0;
  } else {
    v2_len = l2/(n2-1);
  }

  v0.push_back(v0_len);
  v0.push_back(0);
  v0.push_back(0);

  v1.push_back(0);
  v1.push_back(v1_len);
  v1.push_back(0);

  v2.push_back(0);
  v2.push_back(0);
  v2.push_back(v2_len);

  AbstractContinuousDomain* temp=new Bruce(v0, v1, v2, n0, n1, n2, origin);
  return temp;
}

AbstractContinuousDomain* rectangle(int n0, int n1,
                                    double l0, double l1)
{
  double v0_len, v1_len;

  Bruce::DimVec v0;
  Bruce::DimVec v1;
  Bruce::DimVec v2;
  Bruce::DimVec origin;

  origin.push_back(0);
  origin.push_back(0);

  if (n0<2) {
    v0_len = 0;
  } else {
    v0_len = l0/(n0-1);
  }

  if (n1<2) {
    v1_len = 0;
  } else {
    v1_len = l1/(n1-1);
  }

  v0.push_back(v0_len);
  v0.push_back(0);

  v1.push_back(0);
  v1.push_back(v1_len);

  AbstractContinuousDomain* temp=new Bruce(v0, v1, v2, n0, n1, 0, origin);
  return temp;
}

} // end of namespace