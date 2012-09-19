
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


#if !defined  DataExpandedTestCase_20040413_H
#define  DataExpandedTestCase_20040413_H

#include <cppunit/TestFixture.h>
#include <cppunit/TestSuite.h>

class DataExpandedTestCase : public CppUnit::TestFixture
{
public:

  // General test case
  void testAll();

  // Test cases to test slicing of DataExpanded objects
  void testSlicing();
  void testSlicing2();
  void testSlicing3();
  void testSliceSetting();
  void testSliceSetting2();

  static CppUnit::TestSuite* suite();
};

#endif

