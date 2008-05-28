
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

#if !defined  HANDLETESTCASE_H
#define HANDLETESTCASE_H

#include <string>

#include "tools/CppUnitTest/TestCase.h"
#include "tools/CppUnitTest/TestSuite.h"
#include "tools/CppUnitTest/TestCaller.h"

class EsysExceptionTestCase : public CppUnitTest::TestCase
{
 private:
  //
 protected:
   //
   // The following methods are for tests
   //
   // test methods
   //
   void testCase0();
   void testCase1();
   void testCase2();

 public:
  EsysExceptionTestCase (std::string name) : CppUnitTest::TestCase (name) {}
  ~EsysExceptionTestCase() {}
  //
  //
  // return the suite of tests to perform
  //
  static CppUnitTest::TestSuite* suite ();
};

#endif
