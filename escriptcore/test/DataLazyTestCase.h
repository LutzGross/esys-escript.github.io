
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* https://github.com/LutzGross/esys-escript.github.io
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
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

