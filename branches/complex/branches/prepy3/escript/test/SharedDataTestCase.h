
/*******************************************************
*
* Copyright (c) 2003-2010 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


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

