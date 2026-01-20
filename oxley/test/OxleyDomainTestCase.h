
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
*
*****************************************************************************/


#ifndef __OXLEY_DOMAIN_TESTCASE_H__
#define __OXLEY_DOMAIN_TESTCASE_H__

#include <cppunit/TestFixture.h>
#include <cppunit/TestSuite.h>

class OxleyDomainTestCase : public CppUnit::TestFixture
{
public:
  void testAll();

  static CppUnit::TestSuite* suite();
};

#endif // __OXLEY_DOMAIN_TESTCASE_H__

