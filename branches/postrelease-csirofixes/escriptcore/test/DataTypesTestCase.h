
/*****************************************************************************
*
* Copyright (c) 2003-2014 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/


#if !defined DataTypesTestCase_20080827_H
#define DataTypesTestCase_20080827_H

#include <cppunit/TestFixture.h>
#include <cppunit/TestSuite.h>


#define REL_TOL ((double)1.e-10)

class DataTypesTestCase : public CppUnit::TestFixture
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
  void testShapeFns();

  static CppUnit::TestSuite* suite();
};

#endif

