
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/

#include <escript/DataExpanded.h>
#include "DataExpandedTestCase.h"
#include <escript/DataReady.h>
#include <escript/EsysException.h>
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
getRefRO(DataReady& data,int i, int j)
{
   return data.getVectorRO()[getRelIndex(data.getShape(),i,j)];
}

RealVectorType::const_reference
getRefRO(DataReady& data,int i, int j,int k)
{
   return data.getVectorRO()[getRelIndex(data.getShape(),i,j,k)];
}

RealVectorType::reference
getDRef(RealVectorType& data,const ShapeType& shape,int i, int j)
{
   return data[getRelIndex(shape,i,j)];
}

RealVectorType::reference
getDRef(RealVectorType& data,const ShapeType& shape,int i, int j, int k)
{
   return data[getRelIndex(shape,i,j,k)];
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


void DataExpandedTestCase::testAll()
{
  cout << endl;

  //
  // Create a rank 1 pointData
  DataTypes::ShapeType shape;
  shape.push_back(3);
  DataTypes::RealVectorType data(DataTypes::noValues(shape),0);
//  DataArrayView pointData(data,shape);

  //
  // Assign an arbitrary value
  data[0]=0.0;
  data[1]=1.0;
  data[2]=2.0;

  //
  // Test constructor
  cout << "\tTest DataExpanded constructor." << endl;
  DataExpanded testData(FunctionSpace(), shape, data);

//   cout << "\tTest getLength." << endl;
//   CPPUNIT_ASSERT(testData.getLength()==pointData.getNoValues());

  cout << "\tTest getDataAtOffset." << endl;
  for (int i=0;i<testData.getShape()[0];i++) {
      CPPUNIT_ASSERT(testData.getDataAtOffsetRO(i)==data[i]);
  }

  cout << "\tVerify data point attributes." << endl;
//  DataArrayView dataView=testData.getPointDataView();
  CPPUNIT_ASSERT(testData.getRank()==shape.size());
  CPPUNIT_ASSERT(testData.getNoValues()==shape[0]*1);
  CPPUNIT_ASSERT(testData.getShape()[0]==shape[0]);
  CPPUNIT_ASSERT(testData.getNumDPPSample()==1);
  CPPUNIT_ASSERT(testData.getNumSamples()==1);
  CPPUNIT_ASSERT(testData.validSamplePointNo(testData.getNumDPPSample()-1));
  CPPUNIT_ASSERT(testData.validSampleNo(testData.getNumSamples()-1));

  //
  // Test alternative constructor
  cout << "\tTest DataExpanded alternative constructor." << endl;
  data[0]=0.0;
  data[1]=1.0;
  data[2]=2.0;
  FunctionSpace tmp_fns;
  DataExpanded testData1(tmp_fns,shape,data);

//   cout << "\tTest getLength." << endl;
//   CPPUNIT_ASSERT(testData1.getLength()==pointData.noValues());

//   cout << "\tTest getPointDataView." << endl;
//   for (int i=0;i<testData1.getPointDataView().getShape()[0];i++) {
//       CPPUNIT_ASSERT(testData1.getPointDataView()(i)==pointData(i));
//   }

  cout << "\tVerify data point attributes." << endl;
//  dataView=testData1.getPointDataView();
  CPPUNIT_ASSERT(testData1.getRank()==shape.size());
  CPPUNIT_ASSERT(testData1.getNoValues()==shape[0]*1);
  CPPUNIT_ASSERT(testData1.getShape()[0]==shape[0]);
  CPPUNIT_ASSERT(testData1.getNumDPPSample()==1);
  CPPUNIT_ASSERT(testData1.getNumSamples()==1);
  CPPUNIT_ASSERT(testData1.validSamplePointNo(testData1.getNumDPPSample()-1));
  CPPUNIT_ASSERT(testData1.validSampleNo(testData1.getNumSamples()-1));

  //
  // Test copy constructor
  cout << "\tTest DataExpanded copy constructor." << endl;
  DataExpanded testData2(testData);

  cout << "\tTest getLength." << endl;
  CPPUNIT_ASSERT(testData2.getLength()==data.size());

  cout << "\tTest getPointDataView." << endl;
  for (int i=0;i<testData2.getShape()[0];i++) {
    CPPUNIT_ASSERT(testData2.getDataAtOffsetRO(i)==data[i]);
  }

  cout << "\tVerify data point attributes." << endl;
//   dataView=testData2.getPointDataView();
  CPPUNIT_ASSERT(testData2.getRank()==shape.size());
  CPPUNIT_ASSERT(testData2.getNoValues()==shape[0]*1);
  CPPUNIT_ASSERT(testData2.getShape()[0]==shape[0]);
  CPPUNIT_ASSERT(testData2.getNumDPPSample()==1);
  CPPUNIT_ASSERT(testData2.getNumSamples()==1);
  CPPUNIT_ASSERT(testData2.validSamplePointNo(testData2.getNumDPPSample()-1));
  CPPUNIT_ASSERT(testData2.validSampleNo(testData2.getNumSamples()-1));

}


void DataExpandedTestCase::testSlicing() {

  cout << endl;

  //
  // Create a rank 1 pointData
  DataTypes::ShapeType shape;
  shape.push_back(3);
  DataTypes::RealVectorType data(DataTypes::noValues(shape),0);
//   DataArrayView pointData(data,shape);

  //
  // Assign an arbitrary value
  data[0]=0.0;
  data[1]=1.0;
  data[2]=2.0;

  //
  // Create object to test
  cout << "\tCreate rank 1 DataExpanded object." << endl;
  DataExpanded testData(FunctionSpace(),shape,data);

  cout << "\tTest slicing (whole object)." << endl;
  DataTypes::RegionType region;
  region.push_back(DataTypes::RegionType::value_type(0,shape[0]));

  DataReady_ptr testData2=resolveAndDelete(testData.getSlice(region));

  //
  // Verify data values
  cout << "\tVerify data point values." << endl;
  for (int i=0;i<testData2->getShape()[0];i++) {
    CPPUNIT_ASSERT(testData2->getDataAtOffsetRO(i)==data[region[0].first+i]);
  }
}

void DataExpandedTestCase::testSlicing2() {

  cout << endl;

  //
  // Create a rank 2 pointData
  DataTypes::ShapeType shape;
  shape.push_back(3);
  shape.push_back(3);
  DataTypes::RealVectorType data(DataTypes::noValues(shape),0);
//  DataArrayView pointData(data,shape);

  //
  // Assign an arbitrary value
  getDRef(data,shape,0,0)=0.0;
  getDRef(data,shape,1,0)=1.0;
  getDRef(data,shape,2,0)=2.0;
  getDRef(data,shape,0,1)=3.0;
  getDRef(data,shape,1,1)=4.0;
  getDRef(data,shape,2,1)=5.0;
  getDRef(data,shape,0,2)=6.0;
  getDRef(data,shape,1,2)=7.0;
  getDRef(data,shape,2,2)=8.0;

  //
  // Create object to test
  cout << "\tCreate rank 2 DataExpanded object." << endl;
  DataExpanded testData(FunctionSpace(),shape,data);

  cout << "\tTest slicing (part object)." << endl;
  DataTypes::RegionType region;
  region.push_back(DataTypes::RegionType::value_type(0,2));
  region.push_back(DataTypes::RegionType::value_type(0,2));
  
  DataReady_ptr testData2=resolveAndDelete(testData.getSlice(region));

  //
  // Verify data values
  cout << "\tVerify data point values." << endl;
  for (int j=0;j<testData2->getShape()[1];j++) {
    for (int i=0;i<testData2->getShape()[0];i++) {
      CPPUNIT_ASSERT(getRefRO(*testData2,i,j)==data[getRelIndex(shape,region[0].first+i,region[1].first+j)]);
    }
  }

  cout << "\tVerify data point attributes." << endl;
//   DataArrayView dataView=testData2->getPointDataView();
  CPPUNIT_ASSERT(testData2->getRank()==region.size());
  CPPUNIT_ASSERT(testData2->getNoValues()==(region[0].second-region[0].first)*(region[1].second-region[1].first));
  CPPUNIT_ASSERT(testData2->getShape()[0]==(region[0].second-region[0].first));
  CPPUNIT_ASSERT(testData2->getShape()[1]==(region[1].second-region[1].first));

  cout << "\tTest slicing (part object)." << endl;
  DataTypes::RegionType region2;
  region2.push_back(DataTypes::RegionType::value_type(1,3));
  region2.push_back(DataTypes::RegionType::value_type(1,3));
  DataReady_ptr testData3=resolveAndDelete(testData.getSlice(region2));
  //
  // Verify data values
  cout << "\tVerify data point values." << endl;
  for (int j=0;j<testData3->getShape()[1];j++) {
    for (int i=0;i<testData3->getShape()[0];i++) {
      CPPUNIT_ASSERT(getRefRO(*testData3,i,j)==getDRef(data,shape,region2[0].first+i,region2[1].first+j));
    }
  }

  cout << "\tVerify data point attributes." << endl;
//   dataView=testData3->getPointDataView();
  CPPUNIT_ASSERT(testData3->getRank()==region2.size());
  CPPUNIT_ASSERT(testData3->getNoValues()==(region2[0].second-region2[0].first)*(region2[1].second-region2[1].first));
  CPPUNIT_ASSERT(testData3->getShape()[0]==(region2[0].second-region2[0].first));
  CPPUNIT_ASSERT(testData3->getShape()[1]==(region2[1].second-region2[1].first));

//   delete testData2;
//   delete testData3;

}

void DataExpandedTestCase::testSlicing3() {

  cout << endl;

  //
  // Create a rank 3 pointData
  DataTypes::ShapeType shape;
  shape.push_back(3);
  shape.push_back(3);
  shape.push_back(3);
  DataTypes::RealVectorType data(DataTypes::noValues(shape),0);
//   DataArrayView pointData(data,shape);

  //
  // Assign an arbitrary value
  getDRef(data,shape,0,0,0)=0.0;
  getDRef(data,shape,1,0,0)=1.0;
  getDRef(data,shape,2,0,0)=2.0;
  getDRef(data,shape,0,1,0)=3.0;
  getDRef(data,shape,1,1,0)=4.0;
  getDRef(data,shape,2,1,0)=5.0;
  getDRef(data,shape,0,2,0)=6.0;
  getDRef(data,shape,1,2,0)=7.0;
  getDRef(data,shape,2,2,0)=8.0;

  getDRef(data,shape,0,0,1)=9.0;
  getDRef(data,shape,1,0,1)=10.0;
  getDRef(data,shape,2,0,1)=11.0;
  getDRef(data,shape,0,1,1)=12.0;
  getDRef(data,shape,1,1,1)=13.0;
  getDRef(data,shape,2,1,1)=14.0;
  getDRef(data,shape,0,2,1)=15.0;
  getDRef(data,shape,1,2,1)=16.0;
  getDRef(data,shape,2,2,1)=17.0;

  getDRef(data,shape,0,0,2)=18.0;
  getDRef(data,shape,1,0,2)=19.0;
  getDRef(data,shape,2,0,2)=20.0;
  getDRef(data,shape,0,1,2)=21.0;
  getDRef(data,shape,1,1,2)=22.0;
  getDRef(data,shape,2,1,2)=23.0;
  getDRef(data,shape,0,2,2)=24.0;
  getDRef(data,shape,1,2,2)=25.0;
  getDRef(data,shape,2,2,2)=26.0;

  //
  // Create object to test
  cout << "\tCreate rank 3 DataExpanded object." << endl;
  DataExpanded testData(FunctionSpace(),shape,data);

  cout << "\tTest slicing (part object)." << endl;
  DataTypes::RegionType region;
  region.push_back(DataTypes::RegionType::value_type(0,2));
  region.push_back(DataTypes::RegionType::value_type(0,2));
  region.push_back(DataTypes::RegionType::value_type(0,2));
  DataReady_ptr testData2=resolveAndDelete(testData.getSlice(region));

  //
  // Verify data values
  cout << "\tVerify data point values." << endl;
  for (int k=0;k<testData2->getShape()[2];k++) {
    for (int j=0;j<testData2->getShape()[1];j++) {
      for (int i=0;i<testData2->getShape()[0];i++) {
        CPPUNIT_ASSERT(getRefRO(*testData2,i,j,k)==getDRef(data,shape,region[0].first+i,
                                                               region[1].first+j,
                                                               region[2].first+k));
      }
    }
  }

  cout << "\tVerify data point attributes." << endl;
//   DataArrayView dataView=testData2->getPointDataView();
  CPPUNIT_ASSERT(testData2->getRank()==region.size());
  CPPUNIT_ASSERT(testData2->getNoValues()==(region[0].second-region[0].first)
                               *(region[1].second-region[1].first)
                               *(region[2].second-region[2].first));
  CPPUNIT_ASSERT(testData2->getShape()[0]==(region[0].second-region[0].first));
  CPPUNIT_ASSERT(testData2->getShape()[1]==(region[1].second-region[1].first));
  CPPUNIT_ASSERT(testData2->getShape()[2]==(region[2].second-region[2].first));

  cout << "\tTest slicing (part object)." << endl;
  DataTypes::RegionType region2;
  region2.push_back(DataTypes::RegionType::value_type(1,3));
  region2.push_back(DataTypes::RegionType::value_type(1,3));
  region2.push_back(DataTypes::RegionType::value_type(1,3));
  DataReady_ptr testData3=resolveAndDelete(testData.getSlice(region2));

  //
  // Verify data values
  cout << "\tVerify data point values." << endl;
  for (int k=0;k<testData3->getShape()[2];k++) {
    for (int j=0;j<testData3->getShape()[1];j++) {
      for (int i=0;i<testData3->getShape()[0];i++) {
        CPPUNIT_ASSERT(getRefRO(*testData3,i,j,k)==getDRef(data,shape,region2[0].first+i,
                                                               region2[1].first+j,
                                                               region2[2].first+k));
      }
    }
  }

  cout << "\tVerify data point attributes." << endl;
//   dataView=testData2->getPointDataView();
  CPPUNIT_ASSERT(testData2->getRank()==region.size());
  CPPUNIT_ASSERT(testData2->getNoValues()==(region[0].second-region[0].first)
                               *(region[1].second-region[1].first)
                               *(region[2].second-region[2].first));
  CPPUNIT_ASSERT(testData2->getShape()[0]==(region[0].second-region[0].first));
  CPPUNIT_ASSERT(testData2->getShape()[1]==(region[1].second-region[1].first));
  CPPUNIT_ASSERT(testData2->getShape()[2]==(region[2].second-region[2].first));

//   delete testData2;
//  delete testData3;

}


void DataExpandedTestCase::testSliceSetting() {

  cout << endl;

  //
  // Create a rank 2 pointData
  DataTypes::ShapeType shape;
  shape.push_back(2);
  shape.push_back(2);
  DataTypes::RealVectorType data(DataTypes::noValues(shape),0);
//   DataArrayView pointData(data,shape);

  //
  // Assign an arbitrary value
  data[getRelIndex(shape,0,0)]=0.0;
  data[getRelIndex(shape,0,1)]=1.0;
  data[getRelIndex(shape,1,0)]=2.0;
  data[getRelIndex(shape,1,1)]=3.0;

  //
  // Create object to test
  cout << "\tCreate rank 2 DataExpanded object." << endl;
  DataExpanded testData(FunctionSpace(),shape,data);

  //
  // Create another rank 2 pointData
  DataTypes::ShapeType shape2;
  shape2.push_back(3);
  shape2.push_back(3);
  DataTypes::RealVectorType data2(DataTypes::noValues(shape2),0);
//   DataArrayView pointData2(data2,shape2);

  //
  // Assign an arbitrary value
  data2[getRelIndex(shape2,0,0)]=0.1;
  data2[getRelIndex(shape2,0,1)]=1.1;
  data2[getRelIndex(shape2,0,2)]=2.1;
  data2[getRelIndex(shape2,1,0)]=3.1;
  data2[getRelIndex(shape2,1,1)]=4.1;
  data2[getRelIndex(shape2,1,2)]=5.1;
  data2[getRelIndex(shape2,2,0)]=6.1;
  data2[getRelIndex(shape2,2,1)]=7.1;
  data2[getRelIndex(shape2,2,2)]=8.1;

  //
  // Create object to test
  cout << "\tCreate second rank 2 DataExpanded object." << endl;
  DataExpanded testData2(FunctionSpace(),shape2,data2);

  cout << "\tTest slice setting (1)." << endl;

  DataTypes::RegionType region;
  region.push_back(DataTypes::RegionType::value_type(0,2));
  region.push_back(DataTypes::RegionType::value_type(0,2));

  DataTypes::RegionType region2;
  region2.push_back(DataTypes::RegionType::value_type(1,3));
  region2.push_back(DataTypes::RegionType::value_type(1,3));

  DataAbstract* testData3=testData.getSlice(region);

  testData2.setSlice(testData3,region2);
  //
  // Verify data values
  cout << "\tVerify data point values." << endl;
  for (int j=region2[1].first;j<region2[1].second;j++) {
    for (int i=region2[0].first;i<region2[0].second;i++) {
      CPPUNIT_ASSERT(getRefRO(testData2,i,j)==data[getRelIndex(shape,i-(region[0].second-1),j-(region[1].second-1))]);
    }
   }

  delete testData3;

}

void DataExpandedTestCase::testSliceSetting2() {

  cout << endl;

  //
  // Create a rank 0 pointData
  DataTypes::ShapeType shape;
  DataTypes::RealVectorType data(DataTypes::noValues(shape),0);
//   DataArrayView pointData(data,shape);

  //
  // Assign an arbitrary value
  data[0]=0.0;

  //
  // Create object to test
  cout << "\tCreate rank 0 DataExpanded object." << endl;
  DataExpanded testData(FunctionSpace(),shape,data);

  //
  // Create a rank 2 pointData
  DataTypes::ShapeType shape2;
  shape2.push_back(3);
  shape2.push_back(3);
  DataTypes::RealVectorType data2(DataTypes::noValues(shape2),0);
//   DataArrayView pointData2(data2,shape2);

  //
  // Assign an arbitrary value
  data2[getRelIndex(shape2,0,0)]=0.1;
  data2[getRelIndex(shape2,0,1)]=1.1;
  data2[getRelIndex(shape2,0,2)]=2.1;
  data2[getRelIndex(shape2,1,0)]=3.1;
  data2[getRelIndex(shape2,1,1)]=4.1;
  data2[getRelIndex(shape2,1,2)]=5.1;
  data2[getRelIndex(shape2,2,0)]=6.1;
  data2[getRelIndex(shape2,2,1)]=7.1;
  data2[getRelIndex(shape2,2,2)]=8.1;

  //
  // Create object to test
  cout << "\tCreate rank 2 DataExpanded object." << endl;
  DataExpanded testData2(FunctionSpace(),shape2, data2);

  cout << "\tTest slice setting (1)." << endl;

  DataTypes::RegionType region;

  DataTypes::RegionType region2;
  region2.push_back(DataTypes::RegionType::value_type(1,1));
  region2.push_back(DataTypes::RegionType::value_type(1,1));

  DataAbstract* testData3=testData.getSlice(region);

  testData2.setSlice(testData3,region2);

  //
  // Verify data values
  cout << "\tVerify data point values." << endl;
  for (int j=region2[1].first;j<region2[1].second;j++) {
    for (int i=region2[0].first;i<region2[0].second;i++) {
      CPPUNIT_ASSERT(getRefRO(testData2,i,j)==data[0]);
    }
   }

  delete testData3;

}


TestSuite* DataExpandedTestCase::suite()
{
  // Create the suite of tests to perform.
  TestSuite *testSuite = new TestSuite("DataExpandedTestCase");
  testSuite->addTest(new TestCaller<DataExpandedTestCase>(
              "testAll",&DataExpandedTestCase::testAll));
  testSuite->addTest(new TestCaller<DataExpandedTestCase>(
              "testSlicing",&DataExpandedTestCase::testSlicing));
  testSuite->addTest(new TestCaller<DataExpandedTestCase>(
              "testSlicing2",&DataExpandedTestCase::testSlicing2));
  testSuite->addTest(new TestCaller<DataExpandedTestCase>(
              "testSlicing3",&DataExpandedTestCase::testSlicing3));
  testSuite->addTest(new TestCaller<DataExpandedTestCase>(
              "testSliceSetting",&DataExpandedTestCase::testSliceSetting));
  testSuite->addTest(new TestCaller<DataExpandedTestCase>(
              "testSliceSetting2",&DataExpandedTestCase::testSliceSetting2));
  return testSuite;
}

