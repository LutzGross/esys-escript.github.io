
/* $Id: DataTypesTestCase.h 1312 2008-08-27 06:18:44Z jfenwick $ */

/*******************************************************
 *
 *           Copyright 2003-2007 by ACceSS MNRF
 *       Copyright 2007 by University of Queensland
 *
 *                http://esscc.uq.edu.au
 *        Primary Business: Queensland, Australia
 *  Licensed under the Open Software License version 3.0
 *     http://www.opensource.org/licenses/osl-3.0.php
 *
 *******************************************************/

#if !defined  DataTypesTestCase_20080827_H
#define  DataTypesTestCase_20080827_H

#include "tools/CppUnitTest/TestCase.h"
#include "tools/CppUnitTest/TestSuite.h"
#include "tools/CppUnitTest/TestCaller.h"


#define REL_TOL ((double)1.e-10)

class DataTypesTestCase : public CppUnitTest::TestCase
{
 public:

  //
  // setUp is called before each test method to set up test state
  void setUp();
  //
  // tearDown is called after each test method is called.
  void tearDown(); 

  //
  // A test method must return void and have no arguments
  // DataArrayView class
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

  DataTypesTestCase (std::string name) : TestCase (name) {}
  ~DataTypesTestCase() {}
  //
  //
  // return the suite of tests to perform
  //
  static CppUnitTest::TestSuite* suite ();
};

#endif
