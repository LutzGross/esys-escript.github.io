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

#include "brucecpp/Bruce.h"

#include "escriptcpp/FunctionSpaceFactory.h"

#include "esysUtils/EsysException.h"

#include "BruceTestCase.h"

using namespace CppUnitTest;

using namespace escript;
using namespace esysUtils;
using namespace bruce;

using namespace std;

void BruceTestCase::setUp() {
  //
  // This is called before each test is run
 
}

void BruceTestCase::tearDown() {
  //
  // This is called after each test has been run
 
}

void BruceTestCase::testConstructorException() {

  cout << endl;

  //
  // test constructor throws an exception for invalid arguments

  Bruce::DimVec v0;
  Bruce::DimVec v1;
  Bruce::DimVec v2;
  Bruce::DimVec origin;

  // a 4d origin is illegal
  origin.push_back(0);
  origin.push_back(0);
  origin.push_back(0);
  origin.push_back(0);

  try {
    Bruce testbruce(v0, v1, v2, 0, 0, 0, origin);
    assert(false);
  }
  catch (EsysException& e) {
    assert(true);
  }

}

void BruceTestCase::testNull() {

  cout << endl;

  //
  // test null case

  cout << "\tTest null Bruce" << endl;

  std::pair<int,int> dataShape;

  Bruce testbruce;

  assert(testbruce.getDim()==0);

  AbstractContinuousDomain testbruce_asAbstractContinuousDomain = testbruce.asAbstractContinuousDomain();
  AbstractDomain testbruce_asAbstractDomain = testbruce.asAbstractDomain();

  assert(testbruce.getDescription() == "Bruce");

  assert(testbruce.isValidFunctionSpaceType(0));
  assert(testbruce.isValidFunctionSpaceType(1));
  assert(!testbruce.isValidFunctionSpaceType(3));

  assert(testbruce.functionSpaceTypeAsString(0)=="Bruce_ContinuousFunction");
  assert(testbruce.functionSpaceTypeAsString(1)=="Bruce_Function");
  try {
    testbruce.functionSpaceTypeAsString(3);
    assert(false);
  }
  catch (EsysException& e) {
    assert(true);
  }

  assert(testbruce.getContinuousFunctionCode()==0);
  assert(testbruce.getFunctionCode()==1);

  dataShape = testbruce.getDataShape(0);
  cout << dataShape.first << " " << dataShape.second << "\n";
  assert(dataShape.first==1);
  assert(dataShape.second==0);
  dataShape = testbruce.getDataShape(1);
  assert(dataShape.first==1);
  assert(dataShape.second==0);
  try {
    dataShape = testbruce.getDataShape(3);
    assert(false);
  }
  catch (EsysException& e) {
    assert(true);
  }

}

void BruceTestCase::testZero() {

  cout << endl;

  //
  // test zero case

  std::pair<int,int> dataShape;

  cout << "\tTest zero Bruce" << endl;

  Bruce::DimVec v0;
  Bruce::DimVec v1;
  Bruce::DimVec v2;
  Bruce::DimVec origin;

  Bruce testbruce0(v0, v1, v2, 0, 0, 0, origin);

  assert(testbruce0.getDim()==0);

  AbstractContinuousDomain testbruce0_asAbstractContinuousDomain = testbruce0.asAbstractContinuousDomain();
  AbstractDomain testbruce0_asAbstractDomain = testbruce0.asAbstractDomain();

  Bruce testbruce1(testbruce0);

  assert(testbruce0==testbruce1);

  assert(testbruce1.getDim()==0);

  AbstractContinuousDomain testbruce1_asAbstractContinuousDomain = testbruce1.asAbstractContinuousDomain();
  AbstractDomain testbruce1_asAbstractDomain = testbruce1.asAbstractDomain();

  assert(testbruce1==testbruce0);
  assert(!(testbruce1!=testbruce0));

  dataShape = testbruce1.getDataShape(0);
  assert(dataShape.first==1);
  assert(dataShape.second==1);
  dataShape = testbruce1.getDataShape(1);
  assert(dataShape.first==1);
  assert(dataShape.second==0);

  Bruce testbruce2(v0, v1, v2, 0, 0, 0, origin);

  assert(testbruce0==testbruce2);
  assert(!(testbruce0!=testbruce2));

}

void BruceTestCase::test1d() {

  cout << endl;

  //
  // test 1d case

  cout << "\tTest 1d Bruces" << endl;

  std::pair<int,int> dataShape;

  Bruce::DimVec v0;
  Bruce::DimVec v1;
  Bruce::DimVec v2;
  Bruce::DimVec origin;

  origin.push_back(0);

  Bruce testbruce0(v0, v1, v2, 1, 0, 0, origin);

  assert(testbruce0.getDim()==1);

  AbstractContinuousDomain testbruce0_asAbstractContinuousDomain = testbruce0.asAbstractContinuousDomain();
  AbstractDomain testbruce0_asAbstractDomain = testbruce0.asAbstractDomain();

  dataShape = testbruce0.getDataShape(0);
  assert(dataShape.first==1);
  assert(dataShape.second==1);
  dataShape = testbruce0.getDataShape(1);
  assert(dataShape.first==1);
  assert(dataShape.second==0);

  v0.push_back(1);

  Bruce testbruce1(v0, v1, v2, 2, 0, 0, origin);

  assert(testbruce1.getDim()==1);

  AbstractContinuousDomain testbruce1_asAbstractContinuousDomain = testbruce1.asAbstractContinuousDomain();
  AbstractDomain testbruce1_asAbstractDomain = testbruce1.asAbstractDomain();

  Bruce testbruce2(testbruce1);

  assert(testbruce2==testbruce1);

  assert(testbruce2.getDim()==1);

  AbstractContinuousDomain testbruce2_asAbstractContinuousDomain = testbruce2.asAbstractContinuousDomain();
  AbstractDomain testbruce2_asAbstractDomain = testbruce2.asAbstractDomain();

  assert(testbruce2==testbruce1);
  assert(!(testbruce2!=testbruce1));

  dataShape = testbruce2.getDataShape(0);
  assert(dataShape.first==1);
  assert(dataShape.second==2);
  dataShape = testbruce2.getDataShape(1);
  assert(dataShape.first==1);
  assert(dataShape.second==1);

}

void BruceTestCase::test2d() {

  cout << endl;

  //
  // test 2d case

  std::pair<int,int> dataShape;

  cout << "\tTest 2d Bruces" << endl;

  Bruce::DimVec v0;
  Bruce::DimVec v1;
  Bruce::DimVec v2;
  Bruce::DimVec origin;

  origin.push_back(0);
  origin.push_back(0);

  Bruce testbruce0(v0, v1, v2, 1, 0, 0, origin);

  assert(testbruce0.getDim()==2);

  AbstractContinuousDomain testbruce0_asAbstractContinuousDomain = testbruce0.asAbstractContinuousDomain();
  AbstractDomain testbruce0_asAbstractDomain = testbruce0.asAbstractDomain();

  dataShape = testbruce0.getDataShape(0);
  assert(dataShape.first==1);
  assert(dataShape.second==1);
  dataShape = testbruce0.getDataShape(1);
  assert(dataShape.first==1);
  assert(dataShape.second==0);

  v0.push_back(1);
  v0.push_back(0);

  Bruce testbruce1(v0, v1, v2, 2, 0, 0, origin);

  assert(testbruce1.getDim()==2);

  AbstractContinuousDomain testbruce1_asAbstractContinuousDomain = testbruce1.asAbstractContinuousDomain();
  AbstractDomain testbruce1_asAbstractDomain = testbruce1.asAbstractDomain();

  dataShape = testbruce1.getDataShape(0);
  assert(dataShape.first==1);
  assert(dataShape.second==2);
  dataShape = testbruce1.getDataShape(1);
  assert(dataShape.first==1);
  assert(dataShape.second==1);

  v1.push_back(0);
  v1.push_back(1);

  Bruce testbruce2(v0, v1, v2, 2, 2, 0, origin);

  assert(testbruce2.getDim()==2);

  AbstractContinuousDomain testbruce2_asAbstractContinuousDomain = testbruce2.asAbstractContinuousDomain();
  AbstractDomain testbruce2_asAbstractDomain = testbruce2.asAbstractDomain();

  Bruce testbruce3(testbruce2);

  assert(testbruce2==testbruce3);

  assert(testbruce3.getDim()==2);

  AbstractContinuousDomain testbruce3_asAbstractContinuousDomain = testbruce3.asAbstractContinuousDomain();
  AbstractDomain testbruce3_asAbstractDomain = testbruce3.asAbstractDomain();

  assert(testbruce3==testbruce2);
  assert(!(testbruce3!=testbruce2));

  dataShape = testbruce3.getDataShape(0);
  assert(dataShape.first==1);
  assert(dataShape.second==4);
  dataShape = testbruce3.getDataShape(1);
  assert(dataShape.first==1);
  assert(dataShape.second==1);

}

void BruceTestCase::test3d() {

  cout << endl;

  //
  // test 3d case

  std::pair<int,int> dataShape;

  cout << "\tTest 3d Bruces" << endl;

  Bruce::DimVec v0;
  Bruce::DimVec v1;
  Bruce::DimVec v2;
  Bruce::DimVec origin;

  origin.push_back(0);
  origin.push_back(0);
  origin.push_back(0);

  Bruce testbruce0(v0, v1, v2, 1, 0, 0, origin);

  assert(testbruce0.getDim()==3);

  AbstractContinuousDomain testbruce0_asAbstractContinuousDomain = testbruce0.asAbstractContinuousDomain();
  AbstractDomain testbruce0_asAbstractDomain = testbruce0.asAbstractDomain();

  dataShape = testbruce0.getDataShape(0);
  assert(dataShape.first==1);
  assert(dataShape.second==1);
  dataShape = testbruce0.getDataShape(1);
  assert(dataShape.first==1);
  assert(dataShape.second==0);

  v0.push_back(1);
  v0.push_back(0);
  v0.push_back(0);

  Bruce testbruce1(v0, v1, v2, 2, 0, 0, origin);

  assert(testbruce1.getDim()==3);

  AbstractContinuousDomain testbruce1_asAbstractContinuousDomain = testbruce1.asAbstractContinuousDomain();
  AbstractDomain testbruce1_asAbstractDomain = testbruce1.asAbstractDomain();

  dataShape = testbruce1.getDataShape(0);
  assert(dataShape.first==1);
  assert(dataShape.second==2);
  dataShape = testbruce1.getDataShape(1);
  assert(dataShape.first==1);
  assert(dataShape.second==1);

  v1.push_back(0);
  v1.push_back(1);
  v1.push_back(0);

  Bruce testbruce2(v0, v1, v2, 2, 2, 0, origin);

  assert(testbruce2.getDim()==3);

  AbstractContinuousDomain testbruce2_asAbstractContinuousDomain = testbruce2.asAbstractContinuousDomain();
  AbstractDomain testbruce2_asAbstractDomain = testbruce2.asAbstractDomain();

  dataShape = testbruce2.getDataShape(0);
  assert(dataShape.first==1);
  assert(dataShape.second==4);
  dataShape = testbruce2.getDataShape(1);
  assert(dataShape.first==1);
  assert(dataShape.second==1);

  v2.push_back(0);
  v2.push_back(0);
  v2.push_back(1);

  Bruce testbruce3(v0, v1, v2, 2, 2, 2, origin);

  assert(testbruce3.getDim()==3);

  AbstractContinuousDomain testbruce3_asAbstractContinuousDomain = testbruce3.asAbstractContinuousDomain();
  AbstractDomain testbruce3_asAbstractDomain = testbruce3.asAbstractDomain();

  Bruce testbruce4(testbruce3);

  assert(testbruce4==testbruce3);

  assert(testbruce4.getDim()==3);

  AbstractContinuousDomain testbruce4_asAbstractContinuousDomain = testbruce4.asAbstractContinuousDomain();
  AbstractDomain testbruce4_asAbstractDomain = testbruce4.asAbstractDomain();

  assert(testbruce3==testbruce4);
  assert(!(testbruce3!=testbruce4));

  dataShape = testbruce4.getDataShape(0);
  assert(dataShape.first==1);
  assert(dataShape.second==8);
  dataShape = testbruce4.getDataShape(1);
  assert(dataShape.first==1);
  assert(dataShape.second==1);

}

void BruceTestCase::testBig() {

  cout << endl;

  //
  // test a few larger Bruce objects

  std::pair<int,int> dataShape;

  cout << "\tTest large Bruces" << endl;

  // large 1d Bruce

  Bruce::DimVec v0;
  Bruce::DimVec v1;
  Bruce::DimVec v2;
  Bruce::DimVec origin;

  origin.push_back(0);
  v0.push_back(1);

  Bruce testbruce0(v0, v1, v2, 100, 0, 0, origin);

  assert(testbruce0.getDim()==1);

  dataShape = testbruce0.getDataShape(0);
  assert(dataShape.first==1);
  assert(dataShape.second==100);
  dataShape = testbruce0.getDataShape(1);
  assert(dataShape.first==1);
  assert(dataShape.second==99);

  // large 2d Bruce

  origin.push_back(0);

  v0.push_back(0);
  v1.push_back(0);
  v1.push_back(1);

  Bruce testbruce1(v0, v1, v2, 100, 100, 0, origin);

  assert(testbruce1.getDim()==2);

  dataShape = testbruce1.getDataShape(0);
  assert(dataShape.first==1);
  assert(dataShape.second==10000);
  dataShape = testbruce1.getDataShape(1);
  assert(dataShape.first==1);
  assert(dataShape.second==9801);

  // large 3d Bruce

  origin.push_back(0);

  v0.push_back(0);

  v1.push_back(0);

  v2.push_back(0);
  v2.push_back(0);
  v2.push_back(1);

  Bruce testbruce2(v0, v1, v2, 100, 100, 100, origin);

  assert(testbruce2.getDim()==3);

  dataShape = testbruce2.getDataShape(0);
  assert(dataShape.first==1);
  assert(dataShape.second==1000000);
  dataShape = testbruce2.getDataShape(1);
  assert(dataShape.first==1);
  assert(dataShape.second==970299);

}

void BruceTestCase::testSetToXcon() {

  cout << endl;

  //
  // test setToX method on continuousFunction domains

  cout << "\tTest Bruce::setToX for continuousFunctions" << endl;

  Bruce::DimVec v0;
  Bruce::DimVec v1;
  Bruce::DimVec v2;
  Bruce::DimVec origin;

  DataArrayView::ShapeType shape;

  int sampleNo;
  DataAbstract::ValueType::value_type* sampleData;

  origin.clear();

  cout << "\t0d Bruce" << endl;

  v0.clear();
  v1.clear();
  v2.clear();

  cout << "\t\t0d domain" << endl;

  Bruce testbruce1(v0, v1, v2, 1, 0, 0, origin);
  shape.clear();
  Data data1(0.0,shape,continuousFunction(testbruce1.asAbstractContinuousDomain()));
  testbruce1.setToX(data1);

  sampleNo=0;
  sampleData = data1.getSampleData(sampleNo);
  assert(sampleData[0]==0);
  sampleNo++;

  cout << "\t1d Bruce" << endl;

  v0.clear();
  v1.clear();
  v2.clear();

  origin.push_back(0);

  cout << "\t\t0d domain" << endl;

  Bruce testbruce2(v0, v1, v2, 1, 0, 0, origin);
  shape.clear();
  shape.push_back(1);
  Data data2(0.0,shape,continuousFunction(testbruce2.asAbstractContinuousDomain()),true);
  testbruce2.setToX(data2);

  sampleNo=0;
  sampleData = data2.getSampleData(sampleNo);
  assert(sampleData[0]==0);
  sampleNo++;

  cout << "\t\t1d domain" << endl;

  v0.push_back(1);

  Bruce testbruce3(v0, v1, v2, 10, 0, 0, origin);
  shape.clear();
  shape.push_back(1);
  Data data3(0.0,shape,continuousFunction(testbruce3.asAbstractContinuousDomain()),true);
  testbruce3.setToX(data3);

  sampleNo=0;
  for (int i=0; i<10; i++) {
    DataAbstract::ValueType::value_type* sampleData = data3.getSampleData(sampleNo);
    assert(sampleData[0]==i);
    sampleNo++;
  }

  cout << "\t2d Bruce" << endl;

  v0.clear();
  v1.clear();
  v2.clear();

  origin.push_back(0);

  cout << "\t\t0d domain" << endl;

  Bruce testbruce4(v0, v1, v2, 1, 0, 0, origin);
  shape.clear();
  shape.push_back(2);
  Data data4(0.0,shape,continuousFunction(testbruce4.asAbstractContinuousDomain()),true);
  testbruce4.setToX(data4);

  sampleNo=0;
  sampleData = data4.getSampleData(sampleNo);
  assert(sampleData[0]==0);
  assert(sampleData[1]==0);
  sampleNo++;

  cout << "\t\t1d domain" << endl;

  v0.push_back(1);
  v0.push_back(0);

  Bruce testbruce5(v0, v1, v2, 10, 0, 0, origin);
  shape.clear();
  shape.push_back(2);
  Data data5(0.0,shape,continuousFunction(testbruce5.asAbstractContinuousDomain()),true);
  testbruce5.setToX(data5);

  sampleNo=0;
  for (int i=0; i<10; i++) {
    DataAbstract::ValueType::value_type* sampleData = data5.getSampleData(sampleNo);
    assert(sampleData[0]==i);
    assert(sampleData[1]==0);
    sampleNo++;
  }

  cout << "\t\t2d domain" << endl;

  v1.push_back(0);
  v1.push_back(1);

  Bruce testbruce6(v0, v1, v2, 10, 10, 0, origin);
  shape.clear();
  shape.push_back(2);
  Data data6(0.0,shape,continuousFunction(testbruce6.asAbstractContinuousDomain()),true);
  testbruce6.setToX(data6);

  sampleNo=0;
  for (int i=0; i<10; i++) {
    for (int j=0; j<10; j++) {
      DataAbstract::ValueType::value_type* sampleData = data6.getSampleData(sampleNo);
      assert(sampleData[0]==i);
      assert(sampleData[1]==j);
      sampleNo++;
    }
  }

  cout << "\t3d Bruce" << endl;

  v0.clear();
  v1.clear();
  v2.clear();

  origin.push_back(0);

  cout << "\t\t0d domain" << endl;

  Bruce testbruce7(v0, v1, v2, 1, 0, 0, origin);
  shape.clear();
  shape.push_back(3);
  Data data7(0.0,shape,continuousFunction(testbruce7.asAbstractContinuousDomain()),true);
  testbruce7.setToX(data7);

  sampleNo=0;
  sampleData = data7.getSampleData(sampleNo);
  assert(sampleData[0]==0);
  assert(sampleData[1]==0);
  assert(sampleData[2]==0);
  sampleNo++;

  cout << "\t\t1d domain" << endl;

  v0.push_back(1);
  v0.push_back(0);
  v0.push_back(0);

  Bruce testbruce8(v0, v1, v2, 10, 0, 0, origin);
  shape.clear();
  shape.push_back(3);
  Data data8(0.0,shape,continuousFunction(testbruce8.asAbstractContinuousDomain()),true);
  testbruce8.setToX(data8);

  sampleNo=0;
  for (int i=0; i<10; i++) {
    DataAbstract::ValueType::value_type* sampleData = data8.getSampleData(sampleNo);
    assert(sampleData[0]==i);
    assert(sampleData[1]==0);
    assert(sampleData[2]==0);
    sampleNo++;
  }

  cout << "\t\t2d domain" << endl;

  v1.push_back(0);
  v1.push_back(1);
  v1.push_back(0);

  Bruce testbruce9(v0, v1, v2, 10, 10, 0, origin);
  shape.clear();
  shape.push_back(3);
  Data data9(0.0,shape,continuousFunction(testbruce9.asAbstractContinuousDomain()),true);
  testbruce9.setToX(data9);

  sampleNo=0;
  for (int i=0; i<10; i++) {
    for (int j=0; j<10; j++) {
      DataAbstract::ValueType::value_type* sampleData = data9.getSampleData(sampleNo);
      assert(sampleData[0]==i);
      assert(sampleData[1]==j);
      assert(sampleData[2]==0);
      sampleNo++;
    }
  }

  cout << "\t\t3d domain" << endl;

  v2.push_back(0);
  v2.push_back(0);
  v2.push_back(1);

  Bruce testbruce10(v0, v1, v2, 10, 10, 10, origin);
  shape.clear();
  shape.push_back(3);
  Data data10(0.0,shape,continuousFunction(testbruce10.asAbstractContinuousDomain()),true);
  testbruce10.setToX(data10);

  sampleNo=0;
  for (int i=0; i<10; i++) {
    for (int j=0; j<10; j++) {
      for (int k=0; k<10; k++) {
        DataAbstract::ValueType::value_type* sampleData = data10.getSampleData(sampleNo);
        assert(sampleData[0]==i);
        assert(sampleData[1]==j);
        assert(sampleData[2]==k);
        sampleNo++;
      }
    }
  }

}

void BruceTestCase::testSetToXfun() {

  cout << endl;

  //
  // test setToX method on Function domains

  cout << "\tTest Bruce::setToX for functions" << endl;

  Bruce::DimVec v0;
  Bruce::DimVec v1;
  Bruce::DimVec v2;
  Bruce::DimVec origin;

  DataArrayView::ShapeType shape;

  int sampleNo;
  DataAbstract::ValueType::value_type* sampleData;

  origin.clear();

  cout << "\t1d Bruce" << endl;

  v0.clear();
  v1.clear();
  v2.clear();

  origin.push_back(0);

  cout << "\t\t1d domain" << endl;

  v0.push_back(1);

  Bruce testbruce3(v0, v1, v2, 10, 0, 0, origin);
  shape.clear();
  shape.push_back(1);
  Data data3(0.0,shape,function(testbruce3.asAbstractContinuousDomain()),true);
  testbruce3.setToX(data3);

  sampleNo=0;
  for (int i=0; i<9; i++) {
    DataAbstract::ValueType::value_type* sampleData = data3.getSampleData(sampleNo);
    assert(sampleData[0]==i+0.5);
    sampleNo++;
  }

  cout << "\t2d Bruce" << endl;

  v0.clear();
  v1.clear();
  v2.clear();

  origin.push_back(0);

  cout << "\t\t1d domain" << endl;

  v0.push_back(1);
  v0.push_back(0);

  Bruce testbruce5(v0, v1, v2, 10, 0, 0, origin);
  shape.clear();
  shape.push_back(2);
  Data data5(0.0,shape,function(testbruce5.asAbstractContinuousDomain()),true);
  testbruce5.setToX(data5);

  sampleNo=0;
  for (int i=0; i<9; i++) {
    DataAbstract::ValueType::value_type* sampleData = data5.getSampleData(sampleNo);
    assert(sampleData[0]==i+0.5);
    assert(sampleData[1]==0);
    sampleNo++;
  }

  cout << "\t\t2d domain" << endl;

  v1.push_back(0);
  v1.push_back(1);

  Bruce testbruce6(v0, v1, v2, 10, 10, 0, origin);
  shape.clear();
  shape.push_back(2);
  Data data6(0.0,shape,function(testbruce6.asAbstractContinuousDomain()),true);
  testbruce6.setToX(data6);

  sampleNo=0;
  for (int i=0; i<9; i++) {
    for (int j=0; j<9; j++) {
      DataAbstract::ValueType::value_type* sampleData = data6.getSampleData(sampleNo);
      assert(sampleData[0]==i+0.5);
      assert(sampleData[1]==j+0.5);
      sampleNo++;
    }
  }

  cout << "\t3d Bruce" << endl;

  v0.clear();
  v1.clear();
  v2.clear();

  origin.push_back(0);

  cout << "\t\t1d domain" << endl;

  v0.push_back(1);
  v0.push_back(0);
  v0.push_back(0);

  Bruce testbruce8(v0, v1, v2, 10, 0, 0, origin);
  shape.clear();
  shape.push_back(3);
  Data data8(0.0,shape,function(testbruce8.asAbstractContinuousDomain()),true);
  testbruce8.setToX(data8);

  sampleNo=0;
  for (int i=0; i<9; i++) {
    DataAbstract::ValueType::value_type* sampleData = data8.getSampleData(sampleNo);
    assert(sampleData[0]==i+0.5);
    assert(sampleData[1]==0);
    assert(sampleData[2]==0);
    sampleNo++;
  }

  cout << "\t\t2d domain" << endl;

  v1.push_back(0);
  v1.push_back(1);
  v1.push_back(0);

  Bruce testbruce9(v0, v1, v2, 10, 10, 0, origin);
  shape.clear();
  shape.push_back(3);
  Data data9(0.0,shape,function(testbruce9.asAbstractContinuousDomain()),true);
  testbruce9.setToX(data9);

  sampleNo=0;
  for (int i=0; i<9; i++) {
    for (int j=0; j<9; j++) {
      DataAbstract::ValueType::value_type* sampleData = data9.getSampleData(sampleNo);
      assert(sampleData[0]==i+0.5);
      assert(sampleData[1]==j+0.5);
      assert(sampleData[2]==0);
      sampleNo++;
    }
  }

  cout << "\t\t3d domain" << endl;

  v2.push_back(0);
  v2.push_back(0);
  v2.push_back(1);

  Bruce testbruce10(v0, v1, v2, 10, 10, 10, origin);
  shape.clear();
  shape.push_back(3);
  Data data10(0.0,shape,function(testbruce10.asAbstractContinuousDomain()),true);
  testbruce10.setToX(data10);

  sampleNo=0;
  for (int i=0; i<9; i++) {
    for (int j=0; j<9; j++) {
      for (int k=0; k<9; k++) {
        DataAbstract::ValueType::value_type* sampleData = data10.getSampleData(sampleNo);
        assert(sampleData[0]==i+0.5);
        assert(sampleData[1]==j+0.5);
        assert(sampleData[2]==k+0.5);
        sampleNo++;
      }
    }
  }

}

void BruceTestCase::testSetToSizecon() {

  cout << endl;

  //
  // test setToSize method on continuousFunction domains

  cout << "\tTest Bruce::setToSize for continuousFunctions" << endl;

  Bruce::DimVec v0;
  Bruce::DimVec v1;
  Bruce::DimVec v2;
  Bruce::DimVec origin;

  DataArrayView::ShapeType shape;
  shape.clear();
  shape.push_back(1);

  int sampleNo;
  DataAbstract::ValueType::value_type* sampleData;

  origin.clear();

  cout << "\t0d Bruce" << endl;

  v0.clear();
  v1.clear();
  v2.clear();

  cout << "\t\t0d domain" << endl;

  Bruce testbruce1(v0, v1, v2, 1, 0, 0, origin);
  Data data1(0.0,shape,continuousFunction(testbruce1.asAbstractContinuousDomain()));
  testbruce1.setToSize(data1);

  sampleNo=0;
  sampleData = data1.getSampleData(sampleNo);
  assert(sampleData[0]==0);
  sampleNo++;

  cout << "\t1d Bruce" << endl;

  v0.clear();
  v1.clear();
  v2.clear();

  origin.push_back(0);

  cout << "\t\t0d domain" << endl;

  Bruce testbruce2(v0, v1, v2, 1, 0, 0, origin);
  Data data2(0.0,shape,continuousFunction(testbruce2.asAbstractContinuousDomain()),true);
  testbruce2.setToSize(data2);

  sampleNo=0;
  sampleData = data2.getSampleData(sampleNo);
  assert(sampleData[0]==0);
  sampleNo++;

  cout << "\t\t1d domain" << endl;

  v0.push_back(1);

  Bruce testbruce3(v0, v1, v2, 10, 0, 0, origin);
  Data data3(0.0,shape,continuousFunction(testbruce3.asAbstractContinuousDomain()),true);
  testbruce3.setToSize(data3);

  sampleNo=0;
  for (int i=0; i<10; i++) {
    DataAbstract::ValueType::value_type* sampleData = data3.getSampleData(sampleNo);
    assert(sampleData[0]==1);
    sampleNo++;
  }

  cout << "\t2d Bruce" << endl;

  v0.clear();
  v1.clear();
  v2.clear();

  origin.push_back(0);

  cout << "\t\t0d domain" << endl;

  Bruce testbruce4(v0, v1, v2, 1, 0, 0, origin);
  Data data4(0.0,shape,continuousFunction(testbruce4.asAbstractContinuousDomain()),true);
  testbruce4.setToSize(data4);

  sampleNo=0;
  sampleData = data4.getSampleData(sampleNo);
  assert(sampleData[0]==0);
  sampleNo++;

  cout << "\t\t1d domain" << endl;

  v0.push_back(1);
  v0.push_back(0);

  Bruce testbruce5(v0, v1, v2, 10, 0, 0, origin);
  Data data5(0.0,shape,continuousFunction(testbruce5.asAbstractContinuousDomain()),true);
  testbruce5.setToSize(data5);

  sampleNo=0;
  for (int i=0; i<10; i++) {
    DataAbstract::ValueType::value_type* sampleData = data5.getSampleData(sampleNo);
    assert(sampleData[0]==1);
    sampleNo++;
  }

  cout << "\t\t2d domain" << endl;

  v1.push_back(0);
  v1.push_back(1);

  Bruce testbruce6(v0, v1, v2, 10, 10, 0, origin);
  Data data6(0.0,shape,continuousFunction(testbruce6.asAbstractContinuousDomain()),true);
  testbruce6.setToSize(data6);

  sampleNo=0;
  for (int i=0; i<10; i++) {
    for (int j=0; j<10; j++) {
      DataAbstract::ValueType::value_type* sampleData = data6.getSampleData(sampleNo);
      assert(sampleData[0]==std::sqrt(2.0));
      sampleNo++;
    }
  }

  cout << "\t3d Bruce" << endl;

  v0.clear();
  v1.clear();
  v2.clear();

  origin.push_back(0);

  cout << "\t\t0d domain" << endl;

  Bruce testbruce7(v0, v1, v2, 1, 0, 0, origin);
  Data data7(0.0,shape,continuousFunction(testbruce7.asAbstractContinuousDomain()),true);
  testbruce7.setToSize(data7);

  sampleNo=0;
  sampleData = data7.getSampleData(sampleNo);
  assert(sampleData[0]==0);
  sampleNo++;

  cout << "\t\t1d domain" << endl;

  v0.push_back(1);
  v0.push_back(0);
  v0.push_back(0);

  Bruce testbruce8(v0, v1, v2, 10, 0, 0, origin);
  Data data8(0.0,shape,continuousFunction(testbruce8.asAbstractContinuousDomain()),true);
  testbruce8.setToSize(data8);

  sampleNo=0;
  for (int i=0; i<10; i++) {
    DataAbstract::ValueType::value_type* sampleData = data8.getSampleData(sampleNo);
    assert(sampleData[0]==1);
    sampleNo++;
  }

  cout << "\t\t2d domain" << endl;

  v1.push_back(0);
  v1.push_back(1);
  v1.push_back(0);

  Bruce testbruce9(v0, v1, v2, 10, 10, 0, origin);
  Data data9(0.0,shape,continuousFunction(testbruce9.asAbstractContinuousDomain()),true);
  testbruce9.setToSize(data9);

  sampleNo=0;
  for (int i=0; i<10; i++) {
    for (int j=0; j<10; j++) {
      DataAbstract::ValueType::value_type* sampleData = data9.getSampleData(sampleNo);
      assert(sampleData[0]==std::sqrt(2.0));
      sampleNo++;
    }
  }

  cout << "\t\t3d domain" << endl;

  v2.push_back(0);
  v2.push_back(0);
  v2.push_back(1);

  Bruce testbruce10(v0, v1, v2, 10, 10, 10, origin);
  Data data10(0.0,shape,continuousFunction(testbruce10.asAbstractContinuousDomain()),true);
  testbruce10.setToSize(data10);

  sampleNo=0;
  for (int i=0; i<10; i++) {
    for (int j=0; j<10; j++) {
      for (int k=0; k<10; k++) {
        DataAbstract::ValueType::value_type* sampleData = data10.getSampleData(sampleNo);
        assert(sampleData[0]==std::sqrt(3.0));
        sampleNo++;
      }
    }
  }

}

void BruceTestCase::testSetToSizefun() {

  cout << endl;

  //
  // test setToSize method on Function domains

  cout << "\tTest Bruce::setToSize for functions" << endl;

  Bruce::DimVec v0;
  Bruce::DimVec v1;
  Bruce::DimVec v2;
  Bruce::DimVec origin;

  DataArrayView::ShapeType shape;
  shape.clear();
  shape.push_back(1);

  int sampleNo;
  DataAbstract::ValueType::value_type* sampleData;

  origin.clear();

  cout << "\t1d Bruce" << endl;

  v0.clear();
  v1.clear();
  v2.clear();

  origin.push_back(0);

  cout << "\t\t1d domain" << endl;

  v0.push_back(1);

  Bruce testbruce3(v0, v1, v2, 10, 0, 0, origin);
  Data data3(0.0,shape,function(testbruce3.asAbstractContinuousDomain()),true);
  testbruce3.setToSize(data3);

  sampleNo=0;
  for (int i=0; i<9; i++) {
    DataAbstract::ValueType::value_type* sampleData = data3.getSampleData(sampleNo);
    assert(sampleData[0]==1);
    sampleNo++;
  }

  cout << "\t2d Bruce" << endl;

  v0.clear();
  v1.clear();
  v2.clear();

  origin.push_back(0);

  cout << "\t\t1d domain" << endl;

  v0.push_back(1);
  v0.push_back(0);

  Bruce testbruce5(v0, v1, v2, 10, 0, 0, origin);
  Data data5(0.0,shape,function(testbruce5.asAbstractContinuousDomain()),true);
  testbruce5.setToSize(data5);

  sampleNo=0;
  for (int i=0; i<9; i++) {
    DataAbstract::ValueType::value_type* sampleData = data5.getSampleData(sampleNo);
    assert(sampleData[0]==1);
    sampleNo++;
  }

  cout << "\t\t2d domain" << endl;

  v1.push_back(0);
  v1.push_back(1);

  Bruce testbruce6(v0, v1, v2, 10, 10, 0, origin);
  Data data6(0.0,shape,function(testbruce6.asAbstractContinuousDomain()),true);
  testbruce6.setToSize(data6);

  sampleNo=0;
  for (int i=0; i<9; i++) {
    for (int j=0; j<9; j++) {
      DataAbstract::ValueType::value_type* sampleData = data6.getSampleData(sampleNo);
      assert(sampleData[0]==std::sqrt(2.0));
      sampleNo++;
    }
  }

  cout << "\t3d Bruce" << endl;

  v0.clear();
  v1.clear();
  v2.clear();

  origin.push_back(0);

  cout << "\t\t1d domain" << endl;

  v0.push_back(1);
  v0.push_back(0);
  v0.push_back(0);

  Bruce testbruce8(v0, v1, v2, 10, 0, 0, origin);
  Data data8(0.0,shape,function(testbruce8.asAbstractContinuousDomain()),true);
  testbruce8.setToSize(data8);

  sampleNo=0;
  for (int i=0; i<9; i++) {
    DataAbstract::ValueType::value_type* sampleData = data8.getSampleData(sampleNo);
    assert(sampleData[0]==1);
    sampleNo++;
  }

  cout << "\t\t2d domain" << endl;

  v1.push_back(0);
  v1.push_back(1);
  v1.push_back(0);

  Bruce testbruce9(v0, v1, v2, 10, 10, 0, origin);
  Data data9(0.0,shape,function(testbruce9.asAbstractContinuousDomain()),true);
  testbruce9.setToSize(data9);

  sampleNo=0;
  for (int i=0; i<9; i++) {
    for (int j=0; j<9; j++) {
      DataAbstract::ValueType::value_type* sampleData = data9.getSampleData(sampleNo);
      assert(sampleData[0]==std::sqrt(2.0));
      sampleNo++;
    }
  }

  cout << "\t\t3d domain" << endl;

  v2.push_back(0);
  v2.push_back(0);
  v2.push_back(1);

  Bruce testbruce10(v0, v1, v2, 10, 10, 10, origin);
  Data data10(0.0,shape,function(testbruce10.asAbstractContinuousDomain()),true);
  testbruce10.setToSize(data10);

  sampleNo=0;
  for (int i=0; i<9; i++) {
    for (int j=0; j<9; j++) {
      for (int k=0; k<9; k++) {
        DataAbstract::ValueType::value_type* sampleData = data10.getSampleData(sampleNo);
        assert(sampleData[0]==std::sqrt(3.0));
        sampleNo++;
      }
    }
  }

}

void BruceTestCase::testSetToXException() {

  cout << endl;

  //
  // test setToX method exception on function domains

  Bruce::DimVec v0;
  Bruce::DimVec v1;
  Bruce::DimVec v2;
  Bruce::DimVec origin;

  DataArrayView::ShapeType shape;
  shape.clear();
  shape.push_back(1);

  origin.clear();

  v0.clear();
  v1.clear();
  v2.clear();

  Bruce testbruce1(v0, v1, v2, 1, 0, 0, origin);
  Data data1(0.0,shape,function(testbruce1.asAbstractContinuousDomain()));

  try {
  testbruce1.setToX(data1);
    assert(false);
  }
  catch (EsysException& e) {
    assert(true);
  }

}

void BruceTestCase::testSetToSizeException() {

  cout << endl;

  //
  // test setToSize method exception on function domains

  Bruce::DimVec v0;
  Bruce::DimVec v1;
  Bruce::DimVec v2;
  Bruce::DimVec origin;

  DataArrayView::ShapeType shape;
  shape.clear();
  shape.push_back(1);

  origin.clear();

  v0.clear();
  v1.clear();
  v2.clear();

  Bruce testbruce1(v0, v1, v2, 1, 0, 0, origin);
  Data data1(0.0,shape,function(testbruce1.asAbstractContinuousDomain()));

  try {
  testbruce1.setToSize(data1);
    assert(false);
  }
  catch (EsysException& e) {
    assert(true);
  }

}

void BruceTestCase::testGetReferenceNoFromSampleNo() {

  cout << endl;

  //
  // test getReferenceNoFromSampleNo method

  Bruce::DimVec v0;
  Bruce::DimVec v1;
  Bruce::DimVec v2;
  Bruce::DimVec origin;

  std::pair<int,int> dataShape;
  int numSamples;

  v0.push_back(1);
  v0.push_back(0);
  v0.push_back(0);

  v1.push_back(0);
  v1.push_back(1);
  v1.push_back(0);

  v2.push_back(0);
  v2.push_back(0);
  v2.push_back(1);

  origin.push_back(0);
  origin.push_back(0);
  origin.push_back(0);

  Bruce testbruce(v0, v1, v2, 10, 10, 10, origin);

  cout << "\tTest Bruce::getReferenceNoFromSampleNo for functions" << endl;

  dataShape = testbruce.getDataShape(0);
  numSamples=dataShape.second;

  for (int i=0; i<numSamples; i++) {
    assert(testbruce.getReferenceNoFromSampleNo(0,i)==i);
  }

  cout << "\tTest Bruce::getReferenceNoFromSampleNo for continuousFunctions" << endl;

  dataShape = testbruce.getDataShape(1);
  numSamples=dataShape.second;

  for (int i=0; i<numSamples; i++) {
    assert(testbruce.getReferenceNoFromSampleNo(1,i)==i);
  }

}

TestSuite* BruceTestCase::suite ()
{
  //
  // create the suite of tests to perform.
  TestSuite *testSuite = new TestSuite ("BruceTestCase");

  testSuite->addTest (new TestCaller< BruceTestCase>("testConstructorException",&BruceTestCase::testConstructorException));
  testSuite->addTest (new TestCaller< BruceTestCase>("testNull",&BruceTestCase::testNull));
  testSuite->addTest (new TestCaller< BruceTestCase>("testZero",&BruceTestCase::testZero));
  testSuite->addTest (new TestCaller< BruceTestCase>("test1d",&BruceTestCase::test1d));
  testSuite->addTest (new TestCaller< BruceTestCase>("test2d",&BruceTestCase::test2d));
  testSuite->addTest (new TestCaller< BruceTestCase>("test3d",&BruceTestCase::test3d));
  testSuite->addTest (new TestCaller< BruceTestCase>("testBig",&BruceTestCase::testBig));
  testSuite->addTest (new TestCaller< BruceTestCase>("testSetToXcon",&BruceTestCase::testSetToXcon));
  testSuite->addTest (new TestCaller< BruceTestCase>("testSetToXfun",&BruceTestCase::testSetToXfun));
  testSuite->addTest (new TestCaller< BruceTestCase>("testSetToSizecon",&BruceTestCase::testSetToSizecon));
  testSuite->addTest (new TestCaller< BruceTestCase>("testSetToSizefun",&BruceTestCase::testSetToSizefun));
  testSuite->addTest (new TestCaller< BruceTestCase>("testSetToXException",&BruceTestCase::testSetToXException));
  testSuite->addTest (new TestCaller< BruceTestCase>("testSetToSizeException",&BruceTestCase::testSetToSizeException));
  testSuite->addTest (new TestCaller< BruceTestCase>("testGetReferenceNoFromSampleNo();",&BruceTestCase::testGetReferenceNoFromSampleNo));
  return testSuite;
}
