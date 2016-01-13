
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

