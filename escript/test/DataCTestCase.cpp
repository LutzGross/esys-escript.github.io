
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

#include "escript/Data.h"
extern "C" {
#include "escript/DataC.h"
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
