// $Id$
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

  cout << "\tTesting FMax." << endl;

  FMax fmax;
  assert(fmax(5,6)==6);
  assert(fmax(5,-6)==5);
  assert(fmax(0,0)==0);
  assert(fmax(15,-96)==15);

  DataAlgorithmAdapter<FMax> sup(numeric_limits<double>::max()*-1);
  sup.resetResult();
  sup(-1);
  sup(-2);
  sup(-14);
  sup(3);
  assert(sup.getResult()==3);

  cout << "\tTesting AbsMax." << endl;

  AbsMax absmax;
  assert(absmax(5,6)==6);
  assert(absmax(5,-6)==6);
  assert(absmax(0,0)==0);
  assert(absmax(15,-96)==96);

  DataAlgorithmAdapter<AbsMax> Lsup(0);
  Lsup.resetResult();
  Lsup(-2);
  Lsup(2);
  Lsup(5);
  Lsup(-10);
  assert(Lsup.getResult()==10);

  cout << "\tTesting FMin." << endl;

  FMin fmin;
  assert(fmin(5,6)==5);
  assert(fmin(5,-6)==-6);
  assert(fmin(0,0)==0);
  assert(fmin(15,-96)==-96);

  DataAlgorithmAdapter<FMin> inf(numeric_limits<double>::max());
  inf.resetResult();
  inf(1);
  inf(12);
  inf(2);
  inf(99);
  assert(inf.getResult()==1);

  cout << "\tTesting Length." << endl;

  Length lngth;
  assert(lngth(5,6)==std::sqrt(61.0));
  assert(lngth(5,-6)==std::sqrt(61.0));
  assert(lngth(0,0)==std::sqrt(0.0));
  assert(lngth(15,-96)==std::sqrt(9441.0));

  DataAlgorithmAdapter<Length> length(0);
  length.resetResult();
  length(2);
  length(4);
  length(6);
  length(8);
  assert(length.getResult()==std::sqrt(120.0));
  length.resetResult();
  length(1.5);
  length(2.5);
  length(3.5);
  length(4.5);
  assert(length.getResult()==std::sqrt(41.0));

  cout << "\tTesting Trace." << endl;

  Trace trce;
  assert(trce(5,6)==11);
  assert(trce(5,-6)==-1);
  assert(trce(0,0)==0);
  assert(trce(15,-96)==-81);

  DataAlgorithmAdapter<Trace> trace(0);
  trace.resetResult();
  trace(1);
  trace(2);
  trace(3);
  trace(4);
  trace(5);
  assert(trace.getResult()==15);
  trace.resetResult();
  trace(1.5);
  trace(2.5);
  trace(3.5);
  trace(4.5);
  trace(5.5);
  assert(trace.getResult()==17.5);

}

TestSuite* DataAlgorithmAdapterTestCase::suite ()
{
  //
  // create the suite of tests to perform.
  TestSuite *testSuite = new TestSuite ("DataAlgorithmAdapterTestCase");

  testSuite->addTest (new TestCaller< DataAlgorithmAdapterTestCase>("testAll",&DataAlgorithmAdapterTestCase::testAll));
  return testSuite;
}

