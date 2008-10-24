
/*******************************************************
*
* Copyright (c) 2003-2008 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


#include "escript/DataConstant.h"
#include "escript/FunctionSpace.h"
#include "esysUtils/EsysException.h"

#include "DataConstantTestCase.h"

#include <iostream>

using namespace CppUnitTest;
using namespace escript;
using namespace std;
using namespace esysUtils;
using namespace escript::DataTypes;

void DataConstantTestCase::setUp() {
  //
  // This is called before each test is run
}

void DataConstantTestCase::tearDown() {
  //
  // This is called after each test has been run
}


namespace
{

ValueType::reference
getRef(DataReady& data,int i, int j, int k)
{
   return data.getVector()[getRelIndex(data.getShape(),i,j,k)];
}


DataReady_ptr
resolveAndDelete(DataAbstract* p)
{
   DataReady_ptr p2=p->resolve();
   if (p!=p2.get())
   {
	delete p;
   }
   return p2;
}

}


void DataConstantTestCase::testAll() {

  cout << endl;

  //
  // Create a scalar pointData
  DataTypes::ShapeType shape;
  DataTypes::ValueType data(DataTypes::noValues(shape),0);
//  DataArrayView pointData(data,shape);

  //
  // assign an arbitrary value
  data[0]=1.0;

  //
  // Test construction
  cout << "\tTesting default constructor." << endl;
  DataConstant testData(FunctionSpace(),shape,data);

  cout << "\tTest getLength." << endl;
  assert(testData.getLength()==1);

  shape.push_back(2);
  shape.push_back(3);
  shape.push_back(21);

  cout << "\tTesting alternative constructor." << endl;
  DataTypes::ValueType data1(DataTypes::noValues(shape),1.0);
  // do not call the FunctionSpace constructor directly
  // in the argument of DataConstant
  // GCC chokes on it.
  FunctionSpace tmp_fns;
  DataConstant testData1(tmp_fns, shape, data1);

  for (int k=0;k<shape[2];k++) {
    for (int j=0;j<shape[1];j++) {
      for (int i=0;i<shape[0];i++) {
	assert(getRef(testData1,i,j,k)==1.0);
      }
    }
  }

  cout << "\tTest getLength." << endl;
  assert(testData1.getLength()==126);

  cout << "\tTesting copy constructor." << endl;
  DataConstant testData2(testData1);

  for (int k=0;k<shape[2];k++) {
    for (int j=0;j<shape[1];j++) {
      for (int i=0;i<shape[0];i++) {
	assert(getRef(testData2,i,j,k)==data[0]);
      }
    }
  }

  cout << "\tTest getLength." << endl;
  assert(testData2.getLength()==126);

  cout << "\tVerify data point attributes." << endl;
  assert(testData2.getRank()==3);
  assert(testData2.getNoValues()==126);
  assert(testData2.getShape()[0]==2);
  assert(testData2.getShape()[1]==3);
  assert(testData2.getShape()[2]==21);

  cout << "\tTest slicing (whole object)." << endl;

  DataTypes::RegionType region;
  region.push_back(DataTypes::RegionType::value_type(0,shape[0]));
  region.push_back(DataTypes::RegionType::value_type(0,shape[1]));
  region.push_back(DataTypes::RegionType::value_type(0,shape[2]));

  DataReady_ptr testData3=resolveAndDelete(testData2.getSlice(region));

  for (int k=0;k<shape[2];k++) {
    for (int j=0;j<shape[1];j++) {
      for (int i=0;i<shape[0];i++) {
	assert(getRef(*testData3,i,j,k)==data[0]);
      }
    }
  }

  assert(testData3->getLength()==126);

  cout << "\tVerify data point attributes." << endl;
  assert(testData3->getRank()==3);
  assert(testData3->getNoValues()==126);
  assert(testData3->getShape()[0]==2);
  assert(testData3->getShape()[1]==3);
  assert(testData3->getShape()[2]==21);

  cout << "\tTest slicing (part object)." << endl;

  DataTypes::RegionType region2;
  region2.push_back(DataTypes::RegionType::value_type(0,2));
  region2.push_back(DataTypes::RegionType::value_type(0,2));
  region2.push_back(DataTypes::RegionType::value_type(0,2));

  DataReady_ptr testData4=resolveAndDelete(testData3->getSlice(region2));

  for (int k=0;k<2;k++) {
    for (int j=0;j<2;j++) {
      for (int i=0;i<2;i++) {
	assert(getRef(*testData4,i,j,k)==data[0]);
      }
    }
  }

  assert(testData4->getLength()==8);

  cout << "\tVerify data point attributes." << endl;
  assert(testData4->getRank()==3);
  assert(testData4->getNoValues()==8);
  assert(testData4->getShape()[0]==2);
  assert(testData4->getShape()[1]==2);
  assert(testData4->getShape()[2]==2);

  cout << "\tTest slicing (part object)." << endl;

  DataTypes::RegionType region3;
  region3.push_back(DataTypes::RegionType::value_type(1,2));
  region3.push_back(DataTypes::RegionType::value_type(1,3));
  region3.push_back(DataTypes::RegionType::value_type(5,9));

  DataReady_ptr testData5=resolveAndDelete(testData3->getSlice(region3));

  for (int k=0;k<4;k++) {
    for (int j=0;j<2;j++) {
      for (int i=0;i<1;i++) {
	assert(getRef(*testData5,i,j,k)==data[0]);
      }
    }
  }

  assert(testData5->getLength()==8);

  cout << "\tVerify data point attributes." << endl;
//   dataView=testData5->getPointDataView();
  assert(testData5->getRank()==3);
  assert(testData5->getNoValues()==8);
  assert(testData5->getShape()[0]==1);
  assert(testData5->getShape()[1]==2);
  assert(testData5->getShape()[2]==4);

  cout << "\tTest slice setting (1)." << endl;

  DataTypes::RegionType region4;
  region4.push_back(DataTypes::RegionType::value_type(0,1));
  region4.push_back(DataTypes::RegionType::value_type(0,2));
  region4.push_back(DataTypes::RegionType::value_type(0,4));

  DataAbstract* testData6=testData3->getSlice(region3);

  testData5->setSlice(testData6,region4);

  for (int k=0;k<4;k++) {
    for (int j=0;j<2;j++) {
      for (int i=0;i<1;i++) {
	assert(getRef(*testData5,i,j,k)==data[0]);
      }
    }
  }

  assert(testData5->getLength()==8);

  cout << "\tVerify data point attributes." << endl;
//   dataView=testData5->getPointDataView();
  assert(testData5->getRank()==3);
  assert(testData5->getNoValues()==8);
  assert(testData5->getShape()[0]==1);
  assert(testData5->getShape()[1]==2);
  assert(testData5->getShape()[2]==4);

//   delete testData3;
//   delete testData4;
//   delete testData5;
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
