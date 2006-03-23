//$Id$
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
#include "DataCached.h"
#include "FunctionSpace.h"
#include "EsysException.h"

#include "DataCachedTestCase.h"

using namespace CppUnitTest;
using namespace escript;
using namespace std;
using namespace esysUtils;

void DataCachedTestCase::setUp() {
  //
  // This is called before each test is run
 
}

void DataCachedTestCase::tearDown() {
  //
  // This is called after each test has been run
 
}

void DataCachedTestCase::testAll() {
  //
  // The test code may be entered here
  // There is nothing special about the function name, it may be renamed to
  // something more suitable. 
  // As many test methods as desired may be added.

  cout << endl;

  cout << "\tTest default constructor." << endl;
  DataCached testData;

  cout << "\tTest toString method." << endl;
  assert(testData.toString() == "(Cached Data)");

  try {
    cout << "\tTest getPointOffset." << endl;
    assert(testData.getPointOffset(0,0) == 0);
    assert(false);
  }
  catch (EsysException& e) {
    //cout << e.toString() << endl;
    assert(true);
  }
  
  try {
    cout << "\tTest getDataPoint." << endl;
    // this function also returns a DataArrayView object - should check that
    testData.getDataPoint(0,0);
    assert(false);
  }
  catch (EsysException& e) {
    //cout << e.toString() << endl;
    assert(true);
  }

  cout << "\tTest getLength." << endl;
  assert(testData.getLength() == 0);

  DataArrayView::RegionType region;

  try {
    cout << "\tTest getSlice." << endl;
    assert(testData.getSlice(region) == 0);
    assert(false);
  }
  catch (EsysException& e) {
    //cout << e.toString() << endl;
    assert(true);
  }

  try {
    cout << "\tTest setSlice." << endl;
    testData.setSlice(0,region);
    assert(false);
  }
  catch (EsysException& e) {
    //cout << e.toString() << endl;
    assert(true);
  }

  DataArrayView::ShapeType shape;
  
  try {
    cout << "\tTest reshapeDataPoint." << endl;
    testData.reshapeDataPoint(shape);
    assert(false);
  }
  catch (EsysException& e) {
    //cout << e.toString() << endl;
    assert(true);
  }

}

TestSuite* DataCachedTestCase::suite ()
{
  //
  // create the suite of tests to perform.
  TestSuite *testSuite = new TestSuite ("DataCachedTestCase");

  testSuite->addTest (new TestCaller< DataCachedTestCase>("testAll",&DataCachedTestCase::testAll));
  return testSuite;
}
