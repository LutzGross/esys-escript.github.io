
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/


#if !defined  TaipanTestCase_20050427_H
#define  TaipanTestCase_20050427_H

#include <cppunit/TestFixture.h>
#include <cppunit/TestSuite.h>

class TaipanTestCase : public CppUnit::TestFixture
{
public:
  void testAll();
  void testN1();
  void testN0();

  static CppUnit::TestSuite* suite();
};

#endif

