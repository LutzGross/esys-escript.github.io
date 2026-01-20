
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* https://github.com/LutzGross/esys-escript.github.io
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/

#include <escript/DataConstant.h>

#include "DataConstantTestCase.h"

#include <escript/FunctionSpace.h>

#include <cppunit/TestCaller.h>
#include <iostream>

using namespace CppUnit;
using namespace escript;
using namespace std;
using namespace escript::DataTypes;

namespace
{

RealVectorType::const_reference
getRefRO(DataReady& data,int i, int j, int k)
{
   return data.getVectorRO()[getRelIndex(data.getShape(),i,j,k)];
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


void DataConstantTestCase::testAll()
{
  cout << endl;

  //
  // Create a scalar pointData
  DataTypes::ShapeType shape;
  DataTypes::RealVectorType data(DataTypes::noValues(shape),0);
//  DataArrayView pointData(data,shape);

  //
  // assign an arbitrary value
  data[0]=1.0;

  //
  // Test construction
  cout << "\tTesting default constructor." << endl;
  DataConstant testData(FunctionSpace(),shape,data);

  cout << "\tTest getLength." << endl;
  CPPUNIT_ASSERT(testData.getLength()==1);

  shape.push_back(2);
  shape.push_back(3);
  shape.push_back(21);

  cout << "\tTesting alternative constructor." << endl;
  DataTypes::RealVectorType data1(DataTypes::noValues(shape),1.0);
  // do not call the FunctionSpace constructor directly
  // in the argument of DataConstant
  // GCC chokes on it.
  FunctionSpace tmp_fns;
  DataConstant testData1(tmp_fns, shape, data1);

  for (int k=0;k<shape[2];k++) {
    for (int j=0;j<shape[1];j++) {
      for (int i=0;i<shape[0];i++) {
        CPPUNIT_ASSERT(getRefRO(testData1,i,j,k)==1.0);
      }
    }
  }

  cout << "\tTest getLength." << endl;
  CPPUNIT_ASSERT(testData1.getLength()==126);

  cout << "\tTesting copy constructor." << endl;
  DataConstant testData2(testData1);

  for (int k=0;k<shape[2];k++) {
    for (int j=0;j<shape[1];j++) {
      for (int i=0;i<shape[0];i++) {
        CPPUNIT_ASSERT(getRefRO(testData2,i,j,k)==data[0]);
      }
    }
  }

  cout << "\tTest getLength." << endl;
  CPPUNIT_ASSERT(testData2.getLength()==126);

  cout << "\tVerify data point attributes." << endl;
  CPPUNIT_ASSERT(testData2.getRank()==3);
  CPPUNIT_ASSERT(testData2.getNoValues()==126);
  CPPUNIT_ASSERT(testData2.getShape()[0]==2);
  CPPUNIT_ASSERT(testData2.getShape()[1]==3);
  CPPUNIT_ASSERT(testData2.getShape()[2]==21);

  cout << "\tTest slicing (whole object)." << endl;

  DataTypes::RegionType region;
  region.push_back(DataTypes::RegionType::value_type(0,shape[0]));
  region.push_back(DataTypes::RegionType::value_type(0,shape[1]));
  region.push_back(DataTypes::RegionType::value_type(0,shape[2]));

  DataReady_ptr testData3=resolveAndDelete(testData2.getSlice(region));

  for (int k=0;k<shape[2];k++) {
    for (int j=0;j<shape[1];j++) {
      for (int i=0;i<shape[0];i++) {
        CPPUNIT_ASSERT(getRefRO(*testData3,i,j,k)==data[0]);
      }
    }
  }

  CPPUNIT_ASSERT(testData3->getLength()==126);

  cout << "\tVerify data point attributes." << endl;
  CPPUNIT_ASSERT(testData3->getRank()==3);
  CPPUNIT_ASSERT(testData3->getNoValues()==126);
  CPPUNIT_ASSERT(testData3->getShape()[0]==2);
  CPPUNIT_ASSERT(testData3->getShape()[1]==3);
  CPPUNIT_ASSERT(testData3->getShape()[2]==21);

  cout << "\tTest slicing (part object)." << endl;

  DataTypes::RegionType region2;
  region2.push_back(DataTypes::RegionType::value_type(0,2));
  region2.push_back(DataTypes::RegionType::value_type(0,2));
  region2.push_back(DataTypes::RegionType::value_type(0,2));

  DataReady_ptr testData4=resolveAndDelete(testData3->getSlice(region2));

  for (int k=0;k<2;k++) {
    for (int j=0;j<2;j++) {
      for (int i=0;i<2;i++) {
        CPPUNIT_ASSERT(getRefRO(*testData4,i,j,k)==data[0]);
      }
    }
  }

  CPPUNIT_ASSERT(testData4->getLength()==8);

  cout << "\tVerify data point attributes." << endl;
  CPPUNIT_ASSERT(testData4->getRank()==3);
  CPPUNIT_ASSERT(testData4->getNoValues()==8);
  CPPUNIT_ASSERT(testData4->getShape()[0]==2);
  CPPUNIT_ASSERT(testData4->getShape()[1]==2);
  CPPUNIT_ASSERT(testData4->getShape()[2]==2);

  cout << "\tTest slicing (part object)." << endl;

  DataTypes::RegionType region3;
  region3.push_back(DataTypes::RegionType::value_type(1,2));
  region3.push_back(DataTypes::RegionType::value_type(1,3));
  region3.push_back(DataTypes::RegionType::value_type(5,9));

  DataReady_ptr testData5=resolveAndDelete(testData3->getSlice(region3));

  for (int k=0;k<4;k++) {
    for (int j=0;j<2;j++) {
      for (int i=0;i<1;i++) {
        CPPUNIT_ASSERT(getRefRO(*testData5,i,j,k)==data[0]);
      }
    }
  }

  CPPUNIT_ASSERT(testData5->getLength()==8);

  cout << "\tVerify data point attributes." << endl;
//   dataView=testData5->getPointDataView();
  CPPUNIT_ASSERT(testData5->getRank()==3);
  CPPUNIT_ASSERT(testData5->getNoValues()==8);
  CPPUNIT_ASSERT(testData5->getShape()[0]==1);
  CPPUNIT_ASSERT(testData5->getShape()[1]==2);
  CPPUNIT_ASSERT(testData5->getShape()[2]==4);

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
        CPPUNIT_ASSERT(getRefRO(*testData5,i,j,k)==data[0]);
      }
    }
  }

  CPPUNIT_ASSERT(testData5->getLength()==8);

  cout << "\tVerify data point attributes." << endl;
//   dataView=testData5->getPointDataView();
  CPPUNIT_ASSERT(testData5->getRank()==3);
  CPPUNIT_ASSERT(testData5->getNoValues()==8);
  CPPUNIT_ASSERT(testData5->getShape()[0]==1);
  CPPUNIT_ASSERT(testData5->getShape()[1]==2);
  CPPUNIT_ASSERT(testData5->getShape()[2]==4);

//   delete testData3;
//   delete testData4;
//   delete testData5;
  delete testData6;
}

TestSuite* DataConstantTestCase::suite()
{
  TestSuite *testSuite = new TestSuite("DataConstantTestCase");

  testSuite->addTest(new TestCaller<DataConstantTestCase>(
              "testAll",&DataConstantTestCase::testAll));
  return testSuite;
}

