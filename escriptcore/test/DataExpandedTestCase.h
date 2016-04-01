
/*****************************************************************************
*
* Copyright (c) 2003-2016 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/


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

