
/* $Id$ */

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

#if !defined  DataArrayTestCase_20041006_H
#define  DataArrayTestCase_20041006_H

#include "tools/CppUnitTest/TestCase.h"
#include "tools/CppUnitTest/TestSuite.h"
#include "tools/CppUnitTest/TestCaller.h"

class DataArrayTestCase : public CppUnitTest::TestCase
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
  // DataArray class
  void testAll();

  DataArrayTestCase (std::string name) : TestCase (name) {}
  ~DataArrayTestCase() {}
  //
  //
  // return the suite of tests to perform
  //
  static CppUnitTest::TestSuite* suite ();
};

#endif
