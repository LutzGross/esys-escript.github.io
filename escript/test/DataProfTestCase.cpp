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

#include <iostream>

#include "escript/DataProf.h"
#include "esysUtils/EsysException.h"

#include "DataProfTestCase.h"

using namespace std;
using namespace CppUnitTest;
using namespace escript;
using namespace esysUtils;

void DataProfTestCase::setUp() {
  //
  // This is called before each test is run
 
}

void DataProfTestCase::tearDown() {
  //
  // This is called after each test has been run
 
}

void DataProfTestCase::testAll() {
  //
  // The test code may be entered here
  // There is nothing special about the function name, it may be renamed to
  // something more suitable. 
  // As many test methods as desired may be added.

  DataProf dataProf;

  profDataEntry* newEntry1 = dataProf.newData();
  profDataEntry* newEntry2 = dataProf.newData();

  assert(newEntry1->interpolate == 0);
  assert(newEntry1->grad == 0);
  assert(newEntry1->integrate == 0);
  assert(newEntry1->where == 0);
  assert(newEntry1->unary == 0);
  assert(newEntry1->binary == 0);
  assert(newEntry1->reduction1 == 0);
  assert(newEntry1->reduction2 == 0);
  assert(newEntry1->slicing == 0);

  assert(newEntry2->interpolate == 0);
  assert(newEntry2->grad == 0);
  assert(newEntry2->integrate == 0);
  assert(newEntry2->where == 0);
  assert(newEntry2->unary == 0);
  assert(newEntry2->binary == 0);
  assert(newEntry2->reduction1 == 0);
  assert(newEntry2->reduction2 == 0);
  assert(newEntry2->slicing == 0);

  newEntry1->interpolate = 1;
  newEntry1->grad = 1;
  newEntry1->integrate = 1;
  newEntry1->where = 1;
  newEntry1->unary = 1;
  newEntry1->binary = 1;
  newEntry1->reduction1 = 1;
  newEntry1->reduction2 = 1;
  newEntry1->slicing = 1;

  newEntry2->interpolate = 2;
  newEntry2->grad = 2;
  newEntry2->integrate = 2;
  newEntry2->where = 2;
  newEntry2->unary = 2;
  newEntry2->binary = 2;
  newEntry2->reduction1 = 2;
  newEntry2->reduction2 = 2;
  newEntry2->slicing = 2;

  assert(newEntry1->interpolate == 1);
  assert(newEntry1->grad == 1);
  assert(newEntry1->integrate == 1);
  assert(newEntry1->where == 1);
  assert(newEntry1->unary == 1);
  assert(newEntry1->binary == 1);
  assert(newEntry1->reduction1 == 1);
  assert(newEntry1->reduction2 == 1);
  assert(newEntry1->slicing == 1);

  assert(newEntry2->interpolate == 2);
  assert(newEntry2->grad == 2);
  assert(newEntry2->integrate == 2);
  assert(newEntry2->where == 2);
  assert(newEntry2->unary == 2);
  assert(newEntry2->binary == 2);
  assert(newEntry2->reduction1 == 2);
  assert(newEntry2->reduction2 == 2);
  assert(newEntry2->slicing == 2);

  profDataEntry* newEntry3 = dataProf.newData();
  profDataEntry* newEntry4 = dataProf.newData();

  cout << endl;

}

TestSuite* DataProfTestCase::suite ()
{
  //
  // create the suite of tests to perform.
  TestSuite *testSuite = new TestSuite ("DataProfTestCase");

  testSuite->addTest (new TestCaller< DataProfTestCase>("testAll",&DataProfTestCase::testAll));
  return testSuite;
}

