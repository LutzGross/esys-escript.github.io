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

#if !defined BruceTestCase_20050829_H
#define BruceTestCase_20050829_H

#include "CppUnitTest/TestCase.h"
#include "CppUnitTest/TestSuite.h"
#include "CppUnitTest/TestCaller.h"

class BruceTestCase : public CppUnitTest::TestCase
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
  void testConstructorException();
  void testNull();
  void testZero();
  void test1d();
  void test2d();
  void test3d();
  void testBig();
  void testSetToXcon();
  void testSetToXfun();
  void testSetToSizecon();
  void testSetToSizefun();
  void testSetToXException();
  void testSetToSizeException();
  void testGetReferenceNoFromSampleNo();

  BruceTestCase (std::string name) : TestCase (name) {}
  ~BruceTestCase() {}
  //
  //
  // return the suite of tests to perform
  //
  static CppUnitTest::TestSuite* suite ();
};

#endif
