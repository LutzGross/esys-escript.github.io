
/*******************************************************
*
* Copyright (c) 2003-2009 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


#if !defined  DataTaggedTestCase_20040616_H
#define  DataTaggedTestCase_20040616_H

#include "tools/CppUnitTest/TestCase.h"
#include "tools/CppUnitTest/TestSuite.h"
#include "tools/CppUnitTest/TestCaller.h"


class DataTaggedTestCase : public CppUnitTest::TestCase
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
  // DataTagged class
  void testAll();
  void testAddTaggedValues();
  void testSetTaggedValue();
  void testCopyConstructors();
  void testOperations();
  void testGetSlice();
  void testSetSlice();
//   void testFunctionSpaces();


  DataTaggedTestCase (std::string name) : TestCase (name) {}
  ~DataTaggedTestCase() {}
  //
  //
  // return the suite of tests to perform
  //
  static CppUnitTest::TestSuite* suite ();
};

#endif
