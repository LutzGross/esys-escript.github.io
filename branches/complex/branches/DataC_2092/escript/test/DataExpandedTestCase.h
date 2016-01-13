
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


#if !defined  DataExpandedTestCase_20040413_H
#define  DataExpandedTestCase_20040413_H

#include "tools/CppUnitTest/TestCase.h"
#include "tools/CppUnitTest/TestSuite.h"
#include "tools/CppUnitTest/TestCaller.h"

class DataExpandedTestCase : public CppUnitTest::TestCase
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
  // DataExpanded class

  //
  // General test case
  void testAll();

  // Test cases to test slicing of DataExpanded objects
  void testSlicing();
  void testSlicing2();
  void testSlicing3();
  void testSliceSetting();
  void testSliceSetting2();

  DataExpandedTestCase (std::string name) : TestCase (name) {}
  ~DataExpandedTestCase() {}

  //
  // return the suite of tests to perform
  static CppUnitTest::TestSuite* suite ();
};

#endif
