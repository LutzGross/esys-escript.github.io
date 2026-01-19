
/*****************************************************************************
*
* Copyright (c) 2003-2020 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/


#if !defined SHAREDDATATESTCASE_H
#define SHAREDDATATESTCASE_H

#include <cppunit/TestFixture.h>
#include <cppunit/TestSuite.h>

class SharedDataTestCase : public CppUnit::TestFixture
{
public:

  void testEQ();
  void testCC();
  void testAssign();
  void testSetToZero();
  void setTaggedValueFromCPP();
  void testSetTaggedValueFromCPP();
  void testGetDataAtOffset();
  void testGetDataPoint();
  void testGetSampleRW();

  static CppUnit::TestSuite* suite();
};

#endif

