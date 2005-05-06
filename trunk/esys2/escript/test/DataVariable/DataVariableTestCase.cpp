/* 
 *****************************************************************************
 *                                                                           *
 *       COPYRIGHT  ACcESS  -  All Rights Reserved                           *
 *                                                                           *
 * This software is the property of ACcESS. No part of this code             *
 * may be copied in any form or by any means without the expressed written   *
 * consent of ACcESS.  Copying, use or modification of this software         *
 * by any unauthorised person is illegal unless that person has a software   *
 * license agreement with ACcESS.                                            *
 *                                                                           *
 *****************************************************************************
*/

#include "escript/Data/DataVariable.h"
#include "esysUtils/EsysException.h"

#include "DataVariableTestCase.h"

using namespace std;
using namespace CppUnitTest;
using namespace escript;
using namespace esysUtils;

void DataVariableTestCase::setUp() {
  //
  // This is called before each test is run
}

void DataVariableTestCase::tearDown() {
  //
  // This is called after each test has been run
}

void DataVariableTestCase::testAll() {
  //
  // The test code may be entered here
  // There is nothing special about the function name, it may be renamed to
  // something more suitable. 
  // As many test methods as desired may be added.

  DataVariable d;

}

TestSuite* DataVariableTestCase::suite ()
{
  //
  // create the suite of tests to perform.
  TestSuite *testSuite = new TestSuite ("DataVariableTestCase");

  testSuite->addTest (new TestCaller< DataVariableTestCase>("testAll",&DataVariableTestCase::testAll));
  return testSuite;
}
