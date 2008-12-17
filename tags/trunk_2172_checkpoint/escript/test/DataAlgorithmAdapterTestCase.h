
/*******************************************************
*
* Copyright (c) 2003-2008 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


#if !defined DataAlgorithmAdapterTestCase_20040715_H
#define DataAlgorithmAdapterTestCase_20040715_H

#include "tools/CppUnitTest/TestCase.h"
#include "tools/CppUnitTest/TestSuite.h"
#include "tools/CppUnitTest/TestCaller.h"

#define REL_TOL ((double)1.e-10)

class DataAlgorithmAdapterTestCase : public CppUnitTest::TestCase
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
  // DataAlgorithmAdapter class
  void testAll();
  void testAlgorithm();
  void testDpAlgorithm();

  DataAlgorithmAdapterTestCase (std::string name) : TestCase (name) {}
  ~DataAlgorithmAdapterTestCase() {}
  //
  //
  // return the suite of tests to perform
  //
  static CppUnitTest::TestSuite* suite ();
};

#endif
