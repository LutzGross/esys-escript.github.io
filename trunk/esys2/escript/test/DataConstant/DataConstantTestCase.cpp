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
#include "escript/Data/DataConstant.h"
#include "escript/Data/FunctionSpace.h"
#include "esysUtils/EsysException.h"

#include "DataConstantTestCase.h"

#include <iostream>

using namespace CppUnitTest;
using namespace escript;
using namespace std;
using namespace esysUtils;

void DataConstantTestCase::setUp() {
  //
  // This is called before each test is run
 
}

void DataConstantTestCase::tearDown() {
  //
  // This is called after each test has been run
 
}

void DataConstantTestCase::testAll() {

  cout << endl;

  //
  // Create a scalar pointData
  DataArrayView::ShapeType shape;
  DataArrayView::ValueType data(DataArrayView::noValues(shape),0);
  DataArrayView pointData(data,shape);
  DataArrayView::RegionType region;

  //
  // assign an arbitrary value
  pointData()=1.0;

  //
  // Test construction
  cout << "\tTesting default constructor." << endl;
  DataConstant testData(pointData, FunctionSpace());

  cout << "\tTest reshape." << endl;
  shape.push_back(2);
  shape.push_back(3);
  shape.push_back(21);
  testData.reshapeDataPoint(shape);
  assert(testData.getPointDataView().getRank()==shape.size());

  cout << "\tTest getPointDataView." << endl;
  for (int k=0;k<shape[1];++k) {
    for (int j=0;j<shape[1];++j) {
      for (int i=0;i<shape[0];++i) {
	assert(testData.getPointDataView()(i,j,k)==pointData());
      }
    }
  }

  try {
    cout << "\tTest illegal reshape." << endl;
    testData.reshapeDataPoint(shape);
    assert(false);
  }
  catch (EsysException& e) {
    //cout << e.toString() << endl;
    assert(true);
  }

  cout << "\tTesting copy constructor." << endl;
  DataConstant testData2(testData);

  cout << "\tTest getLength." << endl;
  assert(testData2.getLength() == 126);

}

TestSuite* DataConstantTestCase::suite ()
{
  //
  // create the suite of tests to perform.
  TestSuite *testSuite = new TestSuite ("DataConstantTestCase");

  testSuite->addTest (new TestCaller< DataConstantTestCase>("testAll",&DataConstantTestCase::testAll));
  return testSuite;
}
