
/*****************************************************************************
*
* Copyright (c) 2003-2017 by The University of Queensland
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


#if !defined DataLazyTestCase_H
#define DataLazyTestCase_H

#include <cppunit/TestFixture.h>
#include <cppunit/TestSuite.h>

class DataLazyTestCase : public CppUnit::TestFixture
{
public:
  void testLazy1();
  void testLazy2();
  void testLazy2p();
  void testLazy3();
  void testLazy4();

  static CppUnit::TestSuite* suite();
};

#endif

