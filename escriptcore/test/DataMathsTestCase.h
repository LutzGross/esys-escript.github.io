
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


#if !defined  DataMathsTestCase_20080828_H
#define  DataMathsTestCase_20080828_H

#include <cppunit/TestFixture.h>
#include <cppunit/TestSuite.h>

#define REL_TOL ((double)1.e-10)

class DataMathsTestCase : public CppUnit::TestFixture
{
public:

  void testAll();
  void testShapeToString();
  void testScalarView();
  void testResultSliceShape();
  void testSlicing();
  void testMatMult();
  void testUnaryOp();
  void testBinaryOp();
  void testReductionOp();

  static CppUnit::TestSuite* suite();
};

#endif

