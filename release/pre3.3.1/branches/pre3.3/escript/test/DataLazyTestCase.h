
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

