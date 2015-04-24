
/*****************************************************************************
*
* Copyright (c) 2003-2015 by The University of Queensland
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


#if !defined DataAlgorithmAdapterTestCase_20040715_H
#define DataAlgorithmAdapterTestCase_20040715_H

#include <cppunit/TestFixture.h>
#include <cppunit/TestSuite.h>

#define REL_TOL ((double)1.e-10)

class DataAlgorithmAdapterTestCase : public CppUnit::TestFixture
{
public:
  void testAll();
  void testAlgorithm();
  void testDpAlgorithm();

  static CppUnit::TestSuite* suite();
};

#endif

