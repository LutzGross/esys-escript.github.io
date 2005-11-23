#if !defined  HANDLETESTCASE_H
#define HANDLETESTCASE_H

#include <string>

#include "CppUnitTest/TestCase.h"
#include "CppUnitTest/TestSuite.h"
#include "CppUnitTest/TestCaller.h"

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
