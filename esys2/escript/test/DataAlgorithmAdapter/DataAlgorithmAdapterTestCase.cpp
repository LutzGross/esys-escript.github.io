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
#include "escript/Data/DataAlgorithm.h"
#include "DataAlgorithmAdapterTestCase.h"

#include <iostream>
#include <algorithm>
#include <math.h>
#include <limits>

using namespace CppUnitTest;
using namespace std;
using namespace escript;

void DataAlgorithmAdapterTestCase::setUp() {
  //
  // This is called before each test is run
 
}

void DataAlgorithmAdapterTestCase::tearDown() {
  //
  // This is called after each test has been run
 
}

void DataAlgorithmAdapterTestCase::testAll() {
  cout << endl;
  DataAlgorithmAdapter<FMax> myMax(numeric_limits<double>::min());
  myMax(1);
  myMax(2);
  myMax(3);
  assert(myMax.getResult()==3);
  DataAlgorithmAdapter<AbsMax> Lsup(numeric_limits<double>::min());
  Lsup(-2);
  Lsup(2);
  Lsup(5);
  Lsup(-1000);
  assert(Lsup.getResult()==1000);
  DataAlgorithmAdapter<FMin> inf(numeric_limits<double>::max());
  inf(1);
  inf(2);
  assert(inf.getResult() == 1);
  cout << "Size: " << sizeof(DataAlgorithmAdapter<FMin>) << endl;
}

TestSuite* DataAlgorithmAdapterTestCase::suite ()
{
  //
  // create the suite of tests to perform.
  TestSuite *testSuite = new TestSuite ("DataAlgorithmAdapterTestCase");

  testSuite->addTest (new TestCaller< DataAlgorithmAdapterTestCase>("testAll",&DataAlgorithmAdapterTestCase::testAll));
  return testSuite;
}

