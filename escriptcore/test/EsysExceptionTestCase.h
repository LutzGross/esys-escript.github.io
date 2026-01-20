
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


#if !defined ESYSEXCEPTIONTESTCASE_H
#define ESYSEXCEPTIONTESTCASE_H

#include <cppunit/TestFixture.h>
#include <cppunit/TestSuite.h>

class EsysExceptionTestCase : public CppUnit::TestFixture
{
public:
   void testCase1();
   void testCase2();

   static CppUnit::TestSuite* suite();
};

#endif

