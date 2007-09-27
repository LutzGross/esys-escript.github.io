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
#include "escript/DataConstant.h"
#include "escript/FunctionSpace.h"
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

  //
  // assign an arbitrary value
  pointData()=1.0;

  //
  // Test construction
  cout << "\tTesting default constructor." << endl;
  DataConstant testData(pointData, FunctionSpace());

  cout << "\tTest getLength." << endl;
  assert(testData.getLength()==1);

  shape.push_back(2);
  shape.push_back(3);
  shape.push_back(21);

/*
  cout << "\tTest reshape." << endl;
  testData.reshapeDataPoint(shape);
  assert((unsigned int)testData.getPointDataView().getRank()==shape.size());

  cout << "\tTest getPointDataView." << endl;
  for (int k=0;k<shape[2];k++) {
    for (int j=0;j<shape[1];j++) {
      for (int i=0;i<shape[0];i++) {
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

  cout << "\tVerify data point attributes." << endl;
  DataArrayView dataView=testData.getPointDataView();
  assert(dataView.getRank()==3);
  assert(dataView.noValues()==126);
  assert(dataView.getShape()[0]==2);
  assert(dataView.getShape()[1]==3);
  assert(dataView.getShape()[2]==21);
*/

  cout << "\tTesting alternative constructor." << endl;
  DataArrayView::ValueType data1(DataArrayView::noValues(shape),1.0);
  // do not call the FunctionSpace constructor directly
  // in the argument of DataConstant
  // GCC chokes on it.
  FunctionSpace tmp_fns;
  DataConstant testData1(tmp_fns, shape, data1);

  for (int k=0;k<shape[2];k++) {
    for (int j=0;j<shape[1];j++) {
      for (int i=0;i<shape[0];i++) {
	assert(testData1.getPointDataView()(i,j,k)==1.0);
      }
    }
  }

  cout << "\tTest getLength." << endl;
  assert(testData1.getLength()==126);

  cout << "\tVerify data point attributes." << endl;
  DataArrayView dataView=testData1.getPointDataView();
  assert(dataView.getRank()==3);
  assert(dataView.noValues()==126);
  assert(dataView.getShape()[0]==2);
  assert(dataView.getShape()[1]==3);
  assert(dataView.getShape()[2]==21);

  cout << "\tTesting copy constructor." << endl;
  DataConstant testData2(testData1);

  for (int k=0;k<shape[2];k++) {
    for (int j=0;j<shape[1];j++) {
      for (int i=0;i<shape[0];i++) {
	assert(testData2.getPointDataView()(i,j,k)==pointData());
      }
    }
  }

  cout << "\tTest getLength." << endl;
  assert(testData2.getLength()==126);

  cout << "\tVerify data point attributes." << endl;
  dataView=testData2.getPointDataView();
  assert(dataView.getRank()==3);
  assert(dataView.noValues()==126);
  assert(dataView.getShape()[0]==2);
  assert(dataView.getShape()[1]==3);
  assert(dataView.getShape()[2]==21);

  cout << "\tTest slicing (whole object)." << endl;

  DataArrayView::RegionType region;
  region.push_back(DataArrayView::RegionType::value_type(0,shape[0]));
  region.push_back(DataArrayView::RegionType::value_type(0,shape[1]));
  region.push_back(DataArrayView::RegionType::value_type(0,shape[2]));

  DataAbstract* testData3=testData2.getSlice(region);

  for (int k=0;k<shape[2];k++) {
    for (int j=0;j<shape[1];j++) {
      for (int i=0;i<shape[0];i++) {
	assert(testData3->getPointDataView()(i,j,k)==pointData());
      }
    }
  }

  assert(testData3->getLength()==126);

  cout << "\tVerify data point attributes." << endl;
  dataView=testData3->getPointDataView();
  assert(dataView.getRank()==3);
  assert(dataView.noValues()==126);
  assert(dataView.getShape()[0]==2);
  assert(dataView.getShape()[1]==3);
  assert(dataView.getShape()[2]==21);

  cout << "\tTest slicing (part object)." << endl;

  DataArrayView::RegionType region2;
  region2.push_back(DataArrayView::RegionType::value_type(0,2));
  region2.push_back(DataArrayView::RegionType::value_type(0,2));
  region2.push_back(DataArrayView::RegionType::value_type(0,2));

  DataAbstract* testData4=testData3->getSlice(region2);

  for (int k=0;k<2;k++) {
    for (int j=0;j<2;j++) {
      for (int i=0;i<2;i++) {
	assert(testData4->getPointDataView()(i,j,k)==pointData());
      }
    }
  }

  assert(testData4->getLength()==8);

  cout << "\tVerify data point attributes." << endl;
  dataView=testData4->getPointDataView();
  assert(dataView.getRank()==3);
  assert(dataView.noValues()==8);
  assert(dataView.getShape()[0]==2);
  assert(dataView.getShape()[1]==2);
  assert(dataView.getShape()[2]==2);

  cout << "\tTest slicing (part object)." << endl;

  DataArrayView::RegionType region3;
  region3.push_back(DataArrayView::RegionType::value_type(1,2));
  region3.push_back(DataArrayView::RegionType::value_type(1,3));
  region3.push_back(DataArrayView::RegionType::value_type(5,9));

  DataAbstract* testData5=testData3->getSlice(region3);

  for (int k=0;k<4;k++) {
    for (int j=0;j<2;j++) {
      for (int i=0;i<1;i++) {
	assert(testData5->getPointDataView()(i,j,k)==pointData());
      }
    }
  }

  assert(testData5->getLength()==8);

  cout << "\tVerify data point attributes." << endl;
  dataView=testData5->getPointDataView();
  assert(dataView.getRank()==3);
  assert(dataView.noValues()==8);
  assert(dataView.getShape()[0]==1);
  assert(dataView.getShape()[1]==2);
  assert(dataView.getShape()[2]==4);

  cout << "\tTest slice setting (1)." << endl;

  DataArrayView::RegionType region4;
  region4.push_back(DataArrayView::RegionType::value_type(0,1));
  region4.push_back(DataArrayView::RegionType::value_type(0,2));
  region4.push_back(DataArrayView::RegionType::value_type(0,4));

  DataAbstract* testData6=testData3->getSlice(region3);

  testData5->setSlice(testData6,region4);

  for (int k=0;k<4;k++) {
    for (int j=0;j<2;j++) {
      for (int i=0;i<1;i++) {
	assert(testData5->getPointDataView()(i,j,k)==pointData());
      }
    }
  }

  assert(testData5->getLength()==8);

  cout << "\tVerify data point attributes." << endl;
  dataView=testData5->getPointDataView();
  assert(dataView.getRank()==3);
  assert(dataView.noValues()==8);
  assert(dataView.getShape()[0]==1);
  assert(dataView.getShape()[1]==2);
  assert(dataView.getShape()[2]==4);

  delete testData3;
  delete testData4;
  delete testData5;
  delete testData6;
}

TestSuite* DataConstantTestCase::suite ()
{
  //
  // create the suite of tests to perform.
  TestSuite *testSuite = new TestSuite ("DataConstantTestCase");

  testSuite->addTest (new TestCaller< DataConstantTestCase>("testAll",&DataConstantTestCase::testAll));
  return testSuite;
}
