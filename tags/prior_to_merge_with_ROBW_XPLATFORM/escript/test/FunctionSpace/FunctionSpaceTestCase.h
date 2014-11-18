/* 
 ************************************************************
 *          Copyright 2006 by ACcESS MNRF                   *
 *                                                          *
 *              http://www.access.edu.au                    *
 *       Primary Business: Queensland, Australia            *
 *  Licensed under the Open Software License version 3.0    *
 *     http://www.opensource.org/licenses/osl-3.0.php       *
 *                                                          *
 ************************************************************
*/
#if !defined  FunctionSpaceTestCase_20040903_H
#define  FunctionSpaceTestCase_20040903_H

#include "CppUnitTest/TestCase.h"
#include "CppUnitTest/TestSuite.h"
#include "CppUnitTest/TestCaller.h"

class FunctionSpaceTestCase : public CppUnitTest::TestCase
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
  // FunctionSpace class
  void testAll();

  FunctionSpaceTestCase (std::string name) : TestCase (name) {}
  ~FunctionSpaceTestCase() {}
  //
  //
  // return the suite of tests to perform
  //
  static CppUnitTest::TestSuite* suite ();
};

#endif