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
#if !defined  TaipanTestCase_20050427_H
#define  TaipanTestCase_20050427_H

#include "CppUnitTest/TestCase.h"
#include "CppUnitTest/TestSuite.h"
#include "CppUnitTest/TestCaller.h"

class TaipanTestCase : public CppUnitTest::TestCase
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
  // Taipan class
  void testAll();
  void testN1();
  void testN0();

  TaipanTestCase (std::string name) : TestCase (name) {}
  ~TaipanTestCase() {}
  //
  //
  // return the suite of tests to perform
  //
  static CppUnitTest::TestSuite* suite ();
};

#endif
