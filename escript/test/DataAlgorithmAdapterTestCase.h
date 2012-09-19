
/*******************************************************
*
* Copyright (c) 2003-2012 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


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

