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
#include "escript/Data/FunctionSpace.h"
#include "escript/Data/DataExpanded.h"
#include "esysUtils/EsysException.h"
#include "DataExpandedTestCase.h"

#include <iostream>

using namespace CppUnitTest;
using namespace escript;
using namespace std;
using namespace esysUtils;

void DataExpandedTestCase::setUp() {
  //
  // This is called before each test is run
 
}

void DataExpandedTestCase::tearDown() {
  //
  // This is called after each test has been run
 
}

void DataExpandedTestCase::testReshape() {

  cout << endl;

  //
  // Create a scalar pointData
  DataArrayView::ShapeType shape;
  DataArrayView::ValueType data(DataArrayView::noValues(shape),0);
  DataArrayView pointData(data,shape);

  //
  // assign an arbitrary value
  pointData()=1.0;

  int noDataPointsPerSample=10;
  int noSamples=1000;

  //
  // Test construction
  cout << "Test DataExpanded constructor." << endl;
  DataExpanded testData(pointData,FunctionSpace());

  cout << "Test reshapeDataPoint." << endl;
  shape.push_back(2);
  testData.reshapeDataPoint(shape);
  assert(testData.getPointDataView().getRank() == 1);

  //for (int i=0;i<noSamples;++i) {
   // for (int j=0;j<noDataPointsPerSample;++j) {
    //  assert(testData.getDataPoint(i,j)(0) == pointData());
     // assert(testData.getDataPoint(i,j)(1) == pointData());
    //}
  //}

  try {
    cout << "Test illegal reshape." << endl;
    testData.reshapeDataPoint(shape);
    assert(false);
  }
  catch (EsysException& e) {
    cout << e.toString() << endl;
    assert(true);
  }

  cout << "Test toString." << endl;
  cout << testData.toString() << endl;

  cout << "Test DataExpanded destructor." << endl;

}

void DataExpandedTestCase::testAll() {

  cout << endl;

  //
  // Create a scalar pointData
  DataArrayView::ShapeType shape;
  shape.push_back(3);
  DataArrayView::ValueType data(DataArrayView::noValues(shape),0);
  DataArrayView pointData(data,shape);

  //
  // assign an arbitrary value
  pointData(0)=0.0;
  pointData(1)=1.0;
  pointData(2)=2.0;

  int noDataPointsPerSample=10;
  int noSamples=1000;

  //
  // Test construction
  cout << "Test DataExpanded constructor." << endl;
  DataExpanded testData(pointData,FunctionSpace());

  cout << "Test toString." << endl;
  cout << testData.toString() << endl;

  //for (int i=0;i<noSamples;++i) {
    //for (int j=0;j<noDataPointsPerSample;++j) {
      //assert(testData.getDataPoint(i,j) == pointData);
    //}
  //}

  cout << "Test toString." << endl;
  cout << testData.toString() << endl;

  cout << "Test DataExpanded destructor." << endl;

}

TestSuite* DataExpandedTestCase::suite ()
{
  //
  // create the suite of tests to perform.
  TestSuite *testSuite = new TestSuite ("DataExpandedTestCase");

  testSuite->addTest (new TestCaller< DataExpandedTestCase>("testAll",&DataExpandedTestCase::testAll));
  testSuite->addTest (new TestCaller< DataExpandedTestCase>("testReshape",&DataExpandedTestCase::testReshape));
  return testSuite;
}


