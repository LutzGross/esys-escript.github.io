
/*******************************************************
*
* Copyright (c) 2003-2012 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


#include "DataCTestCase.h"
#include "escript/Data.h"
extern "C" {
#include "escript/DataC.h"
#include "CompareFuncs.h"
}

#include <cppunit/TestCaller.h>
#include <iostream>

using namespace std;
using namespace CppUnit;
using namespace escript;

void DataCTestCase::testAll()
{
  cout << endl;
  cout << "\tTest C interface to escript::Data object." << endl;

  Data myData;
  escriptDataC myDataC=myData.getDataC();
  int typeCode=myData.getFunctionSpace().getTypeCode();

  cout << "\tData typeCode: " << typeCode << endl;
  CPPUNIT_ASSERT(compareTypeCode(&myDataC,typeCode));

  cout << "\tData isEmpty: " << myData.isEmpty() << endl;
  CPPUNIT_ASSERT(compareIsEmpty(&myDataC,myData.isEmpty()));

  cout << "\tData isExpanded: " << myData.isExpanded() << endl;
  CPPUNIT_ASSERT(compareIsExpanded(&myDataC,myData.isExpanded()));

  //cout << "Num DataPoints per sample: " << myData.getNumDPPSample() << " num samples: " <<  myData.getNumSamples() << endl;
  //CPPUNIT_ASSERT(compareNumSamples(&myDataC,myData.getNumDPPSample(),myData.getNumSamples()));

  //DataArrayView::ShapeType tempShape=myData.getPointDataView().getShape();
  //cout << "Data rank: " << tempShape.size() << endl;
  //CPPUNIT_ASSERT(comparePointShape(&myDataC,tempShape.size(),&tempShape[0]));

  //cout << "Data value: " << myData.getSampleData(0)[0] << endl;
  //cout << "Data value: " << myData.getTaggedSampleData(0)[0] << endl;
  //CPPUNIT_ASSERT(compareSampleDataWrite(&myDataC,0,myData.getSampleData(0)));
}

TestSuite* DataCTestCase::suite()
{
  TestSuite *testSuite = new TestSuite("DataCTestCase");

  testSuite->addTest(new TestCaller< DataCTestCase>(
              "testAll",&DataCTestCase::testAll));
  return testSuite;
}

