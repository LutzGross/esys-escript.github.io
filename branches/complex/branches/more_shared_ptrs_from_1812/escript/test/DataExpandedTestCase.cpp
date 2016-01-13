
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


#include "escript/FunctionSpace.h"
#include "escript/DataExpanded.h"
#include "esysUtils/EsysException.h"
#include "DataExpandedTestCase.h"

#include <iostream>

using namespace CppUnitTest;
using namespace escript;
using namespace std;
using namespace esysUtils;
using namespace escript::DataTypes;

void DataExpandedTestCase::setUp() {
  //
  // This is called before each test is run
 
}

void DataExpandedTestCase::tearDown() {
  //
  // This is called after each test has been run
 
}

namespace
{

ValueType::reference
getRef(DataAbstract& data,int i, int j)
{
   return data.getVector()[getRelIndex(data.getShape(),i,j)];
}

ValueType::reference
getRef(DataAbstract& data,int i, int j,int k)
{
   return data.getVector()[getRelIndex(data.getShape(),i,j,k)];
}

ValueType::reference
getDRef(ValueType& data,const ShapeType& shape,int i, int j)
{
   return data[getRelIndex(shape,i,j)];
}

ValueType::reference
getDRef(ValueType& data,const ShapeType& shape,int i, int j, int k)
{
   return data[getRelIndex(shape,i,j,k)];
}

}


void DataExpandedTestCase::testAll() {

  cout << endl;

  //
  // Create a rank 1 pointData
  DataTypes::ShapeType shape;
  shape.push_back(3);
  DataTypes::ValueType data(DataTypes::noValues(shape),0);
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
//   assert(testData.getLength()==pointData.getNoValues());

  cout << "\tTest getDataAtOffset." << endl;
  for (int i=0;i<testData.getShape()[0];i++) {
      assert(testData.getDataAtOffset(i)==data[i]);
  }

  cout << "\tVerify data point attributes." << endl;
//  DataArrayView dataView=testData.getPointDataView();
  assert(testData.getRank()==shape.size());
  assert(testData.getNoValues()==shape[0]*1);
  assert(testData.getShape()[0]==shape[0]);
  assert(testData.getNumDPPSample()==1);
  assert(testData.getNumSamples()==1);
  assert(testData.validSamplePointNo(testData.getNumDPPSample()-1));
  assert(testData.validSampleNo(testData.getNumSamples()-1));

  //
  // Test alternative constructor
  cout << "\tTest DataExpanded alternative constructor." << endl;
  data[0]=0.0;
  data[1]=1.0;
  data[2]=2.0;
  FunctionSpace tmp_fns;
  DataExpanded testData1(tmp_fns,shape,data);

//   cout << "\tTest getLength." << endl;
//   assert(testData1.getLength()==pointData.noValues());

//   cout << "\tTest getPointDataView." << endl;
//   for (int i=0;i<testData1.getPointDataView().getShape()[0];i++) {
//       assert(testData1.getPointDataView()(i)==pointData(i));
//   }

  cout << "\tVerify data point attributes." << endl;
//  dataView=testData1.getPointDataView();
  assert(testData1.getRank()==shape.size());
  assert(testData1.getNoValues()==shape[0]*1);
  assert(testData1.getShape()[0]==shape[0]);
  assert(testData1.getNumDPPSample()==1);
  assert(testData1.getNumSamples()==1);
  assert(testData1.validSamplePointNo(testData1.getNumDPPSample()-1));
  assert(testData1.validSampleNo(testData1.getNumSamples()-1));

  //
  // Test copy constructor
  cout << "\tTest DataExpanded copy constructor." << endl;
  DataExpanded testData2(testData);

  cout << "\tTest getLength." << endl;
  assert(testData2.getLength()==data.size());

  cout << "\tTest getPointDataView." << endl;
  for (int i=0;i<testData2.getShape()[0];i++) {
    assert(testData2.getDataAtOffset(i)==data[i]);
  }

  cout << "\tVerify data point attributes." << endl;
//   dataView=testData2.getPointDataView();
  assert(testData2.getRank()==shape.size());
  assert(testData2.getNoValues()==shape[0]*1);
  assert(testData2.getShape()[0]==shape[0]);
  assert(testData2.getNumDPPSample()==1);
  assert(testData2.getNumSamples()==1);
  assert(testData2.validSamplePointNo(testData2.getNumDPPSample()-1));
  assert(testData2.validSampleNo(testData2.getNumSamples()-1));

}


void DataExpandedTestCase::testSlicing() {

  cout << endl;

  //
  // Create a rank 1 pointData
  DataTypes::ShapeType shape;
  shape.push_back(3);
  DataTypes::ValueType data(DataTypes::noValues(shape),0);
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

  DataAbstract* testData2=testData.getSlice(region);

  //
  // Verify data values
  cout << "\tVerify data point values." << endl;
  for (int i=0;i<testData2->getShape()[0];i++) {
    assert(testData2->getDataAtOffset(i)==data[region[0].first+i]);
  }

//   cout << "\tVerify data point attributes." << endl;
//   DataArrayView dataView=testData2->getPointDataView();
//   assert(data.getRank()==shape.size());
//   assert(data.getNoValues()==shape[0]*1);
//   assert(data.getShape()[0]==shape[0]);

  delete testData2;
}

void DataExpandedTestCase::testSlicing2() {

  cout << endl;

  //
  // Create a rank 2 pointData
  DataTypes::ShapeType shape;
  shape.push_back(3);
  shape.push_back(3);
  DataTypes::ValueType data(DataTypes::noValues(shape),0);
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
  DataAbstract* testData2=testData.getSlice(region);

  //
  // Verify data values
  cout << "\tVerify data point values." << endl;
  for (int j=0;j<testData2->getShape()[1];j++) {
    for (int i=0;i<testData2->getShape()[0];i++) {
      assert(getRef(*testData2,i,j)==data[getRelIndex(shape,region[0].first+i,region[1].first+j)]);
    }
  }

  cout << "\tVerify data point attributes." << endl;
//   DataArrayView dataView=testData2->getPointDataView();
  assert(testData2->getRank()==region.size());
  assert(testData2->getNoValues()==(region[0].second-region[0].first)*(region[1].second-region[1].first));
  assert(testData2->getShape()[0]==(region[0].second-region[0].first));
  assert(testData2->getShape()[1]==(region[1].second-region[1].first));

  cout << "\tTest slicing (part object)." << endl;
  DataTypes::RegionType region2;
  region2.push_back(DataTypes::RegionType::value_type(1,3));
  region2.push_back(DataTypes::RegionType::value_type(1,3));
  DataAbstract* testData3=testData.getSlice(region2);

  //
  // Verify data values
  cout << "\tVerify data point values." << endl;
  for (int j=0;j<testData3->getShape()[1];j++) {
    for (int i=0;i<testData3->getShape()[0];i++) {
      assert(getRef(*testData3,i,j)==getDRef(data,shape,region2[0].first+i,region2[1].first+j));
    }
  }

  cout << "\tVerify data point attributes." << endl;
//   dataView=testData3->getPointDataView();
  assert(testData3->getRank()==region2.size());
  assert(testData3->getNoValues()==(region2[0].second-region2[0].first)*(region2[1].second-region2[1].first));
  assert(testData3->getShape()[0]==(region2[0].second-region2[0].first));
  assert(testData3->getShape()[1]==(region2[1].second-region2[1].first));

  delete testData2;
  delete testData3;

}

void DataExpandedTestCase::testSlicing3() {

  cout << endl;

  //
  // Create a rank 3 pointData
  DataTypes::ShapeType shape;
  shape.push_back(3);
  shape.push_back(3);
  shape.push_back(3);
  DataTypes::ValueType data(DataTypes::noValues(shape),0);
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
  DataAbstract* testData2=testData.getSlice(region);

  //
  // Verify data values
  cout << "\tVerify data point values." << endl;
  for (int k=0;k<testData2->getShape()[2];k++) {
    for (int j=0;j<testData2->getShape()[1];j++) {
      for (int i=0;i<testData2->getShape()[0];i++) {
        assert(getRef(*testData2,i,j,k)==getDRef(data,shape,region[0].first+i,
                                                               region[1].first+j,
                                                               region[2].first+k));
      }
    }
  }

  cout << "\tVerify data point attributes." << endl;
//   DataArrayView dataView=testData2->getPointDataView();
  assert(testData2->getRank()==region.size());
  assert(testData2->getNoValues()==(region[0].second-region[0].first)
                               *(region[1].second-region[1].first)
                               *(region[2].second-region[2].first));
  assert(testData2->getShape()[0]==(region[0].second-region[0].first));
  assert(testData2->getShape()[1]==(region[1].second-region[1].first));
  assert(testData2->getShape()[2]==(region[2].second-region[2].first));

  cout << "\tTest slicing (part object)." << endl;
  DataTypes::RegionType region2;
  region2.push_back(DataTypes::RegionType::value_type(1,3));
  region2.push_back(DataTypes::RegionType::value_type(1,3));
  region2.push_back(DataTypes::RegionType::value_type(1,3));
  DataAbstract* testData3=testData.getSlice(region2);

  //
  // Verify data values
  cout << "\tVerify data point values." << endl;
  for (int k=0;k<testData3->getShape()[2];k++) {
    for (int j=0;j<testData3->getShape()[1];j++) {
      for (int i=0;i<testData3->getShape()[0];i++) {
        assert(getRef(*testData3,i,j,k)==getDRef(data,shape,region2[0].first+i,
                                                               region2[1].first+j,
                                                               region2[2].first+k));
      }
    }
  }

  cout << "\tVerify data point attributes." << endl;
//   dataView=testData2->getPointDataView();
  assert(testData2->getRank()==region.size());
  assert(testData2->getNoValues()==(region[0].second-region[0].first)
                               *(region[1].second-region[1].first)
                               *(region[2].second-region[2].first));
  assert(testData2->getShape()[0]==(region[0].second-region[0].first));
  assert(testData2->getShape()[1]==(region[1].second-region[1].first));
  assert(testData2->getShape()[2]==(region[2].second-region[2].first));

  delete testData2;
  delete testData3;

}


void DataExpandedTestCase::testSliceSetting() {

  cout << endl;

  //
  // Create a rank 2 pointData
  DataTypes::ShapeType shape;
  shape.push_back(2);
  shape.push_back(2);
  DataTypes::ValueType data(DataTypes::noValues(shape),0);
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
  DataTypes::ValueType data2(DataTypes::noValues(shape2),0);
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
      assert(getRef(testData2,i,j)==data[getRelIndex(shape,i-(region[0].second-1),j-(region[1].second-1))]);
    }
   }

  delete testData3;

}

void DataExpandedTestCase::testSliceSetting2() {

  cout << endl;

  //
  // Create a rank 0 pointData
  DataTypes::ShapeType shape;
  DataTypes::ValueType data(DataTypes::noValues(shape),0);
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
  DataTypes::ValueType data2(DataTypes::noValues(shape2),0);
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
      assert(getRef(testData2,i,j)==data[0]);
    }
   }

  delete testData3;

}





TestSuite* DataExpandedTestCase::suite ()
{
  //
  // Create the suite of tests to perform.
  TestSuite *testSuite = new TestSuite ("DataExpandedTestCase");
  testSuite->addTest (new TestCaller< DataExpandedTestCase>("testAll",&DataExpandedTestCase::testAll));
  testSuite->addTest (new TestCaller< DataExpandedTestCase>("testSlicing",&DataExpandedTestCase::testSlicing));
  testSuite->addTest (new TestCaller< DataExpandedTestCase>("testSlicing2",&DataExpandedTestCase::testSlicing2));
  testSuite->addTest (new TestCaller< DataExpandedTestCase>("testSlicing3",&DataExpandedTestCase::testSlicing3));
  testSuite->addTest (new TestCaller< DataExpandedTestCase>("testSliceSetting",&DataExpandedTestCase::testSliceSetting));
  testSuite->addTest (new TestCaller< DataExpandedTestCase>("testSliceSetting2",&DataExpandedTestCase::testSliceSetting2));
  return testSuite;
}
