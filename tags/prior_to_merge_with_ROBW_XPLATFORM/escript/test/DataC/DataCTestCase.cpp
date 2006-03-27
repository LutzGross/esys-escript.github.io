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
#include "Data.h"
extern "C" {
#include "DataC.h"
#include "CompareFuncs.h"
}
#include "DataCTestCase.h"

#include <iostream>

using namespace std;
using namespace CppUnitTest;
using namespace escript;

void DataCTestCase::setUp() {
  //
  // This is called before each test is run
 
}

void DataCTestCase::tearDown() {
  //
  // This is called after each test has been run
 
}

void DataCTestCase::testAll() {

  cout << endl;

  cout << "\tTest C interface to escript::Data object." << endl;

  Data myData;
  escriptDataC myDataC=myData.getDataC();
  int typeCode=myData.getFunctionSpace().getTypeCode();

  cout << "\tData typeCode: " << typeCode << endl;
  assert(compareTypeCode(&myDataC,typeCode));

  cout << "\tData isEmpty: " << myData.isEmpty() << endl;
  assert(compareIsEmpty(&myDataC,myData.isEmpty()));

  cout << "\tData isExpanded: " << myData.isExpanded() << endl;
  assert(compareIsExpanded(&myDataC,myData.isExpanded()));

  //cout << "Num DataPoints per sample: " << myData.getNumDPPSample() << " num samples: " <<  myData.getNumSamples() << endl;
  //assert(compareNumSamples(&myDataC,myData.getNumDPPSample(),myData.getNumSamples()));

  //DataArrayView::ShapeType tempShape=myData.getPointDataView().getShape();
  //cout << "Data rank: " << tempShape.size() << endl;
  //assert(comparePointShape(&myDataC,tempShape.size(),&tempShape[0]));

  //cout << "Data value: " << myData.getSampleData(0)[0] << endl;
  //cout << "Data value: " << myData.getTaggedSampleData(0)[0] << endl;
  //assert(compareSampleDataWrite(&myDataC,0,myData.getSampleData(0)));
}

TestSuite* DataCTestCase::suite ()
{
  //
  // create the suite of tests to perform.
  TestSuite *testSuite = new TestSuite ("DataCTestCase");

  testSuite->addTest (new TestCaller< DataCTestCase>("testAll",&DataCTestCase::testAll));
  return testSuite;
}
