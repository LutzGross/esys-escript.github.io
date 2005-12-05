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
#include "esysUtils/EsysException.h"
#include "escript/Data/Data.h"
#include "escript/Data/FunctionSpace.h"

#include "DataTestCase.h"

#include <iostream>
#include <math.h>

using namespace std;
using namespace CppUnitTest;
using namespace escript;
using namespace esysUtils;

void DataTestCase::setUp() {
  //
  // This is called before each test is run
}

void DataTestCase::tearDown() {
  //
  // This is called after each test has been run
}

void DataTestCase::testSlicing() {

  cout << endl;

  {
    DataArrayView::ShapeType viewShape;
    //
    // weak tests for slicing DataConstant
    cout << "\tTest slicing DataConstant" << endl;
    viewShape.push_back(2);
    viewShape.push_back(3);
    Data temp(1.3,viewShape,FunctionSpace(),false);
    DataArrayView::RegionType region;
    region.push_back(DataArrayView::RegionType::value_type(0,0));
    region.push_back(DataArrayView::RegionType::value_type(0,0));
    Data slice(temp.getSlice(region));
    assert(slice.getDataPointRank()==0);
    assert(slice.getDataPoint(0,0)()==1.3);
    //
    // try the same but this time to produce a matrix containing one value
    region.clear();
    region.push_back(DataArrayView::RegionType::value_type(0,1));
    region.push_back(DataArrayView::RegionType::value_type(0,1));
    slice=temp.getSlice(region);
    assert(slice.getDataPointRank()==2);
    assert(slice.getDataPoint(0,0)(0,0)==1.3);
  }

  {
    DataArrayView::ShapeType viewShape;
    //
    // weak tests for slicing DataExpanded
    cout << "\tTest slicing DataExpanded" << endl;
    viewShape.push_back(2);
    viewShape.push_back(3);
    Data temp(1.3,viewShape,FunctionSpace(),true);
    temp.getDataPoint(0,0)(0,0)=0.0;
    temp.getDataPoint(0,0)(1,1)=1.0;
    DataArrayView::RegionType region;
    region.push_back(DataArrayView::RegionType::value_type(0,0));
    region.push_back(DataArrayView::RegionType::value_type(0,0));
    Data slice(temp.getSlice(region));
    assert(slice.getDataPointRank()==0);
    assert(slice.getDataPoint(0,0)()==0.0);
    //
    // try the same but this time to produce a matrix containing one value
    region.clear();
    region.push_back(DataArrayView::RegionType::value_type(0,1));
    region.push_back(DataArrayView::RegionType::value_type(0,1));
    slice=temp.getSlice(region);
    assert(slice.getDataPointRank()==2);
    assert(slice.getDataPoint(0,0)(0,0)==0.0);
    region.clear();
    region.push_back(DataArrayView::RegionType::value_type(0,2));
    region.push_back(DataArrayView::RegionType::value_type(0,2));
    slice=temp.getSlice(region);
    assert(slice.getDataPoint(0,0)(0,0)==0.0);
    assert(slice.getDataPoint(0,0)(1,1)==1.0);
  }

  {
    DataArrayView::ShapeType viewShape;
    //
    // weak tests for slicing DataTagged
    cout << "\tTest slicing DataTagged" << endl;
    viewShape.push_back(2);
    viewShape.push_back(3);
    Data temp(1.3,viewShape,FunctionSpace(),false);
    //
    // convert the data to tagged
    temp.tag();
    temp.getDataPoint(0,0)(0,0)=0.0;
    temp.getDataPoint(0,0)(1,1)=1.0;
    DataArrayView::RegionType region;
    region.push_back(DataArrayView::RegionType::value_type(0,0));
    region.push_back(DataArrayView::RegionType::value_type(0,0));
    Data slice(temp.getSlice(region));
    assert(slice.getDataPointRank()==0);
    assert(slice.getDataPoint(0,0)()==0.0);
    //
    // try the same but this time to produce a matrix containing one value
    region.clear();
    region.push_back(DataArrayView::RegionType::value_type(0,1));
    region.push_back(DataArrayView::RegionType::value_type(0,1));
    slice=temp.getSlice(region);
    assert(slice.getDataPointRank()==2);
    assert(slice.getDataPoint(0,0)(0,0)==0.0);
    region.clear();
    region.push_back(DataArrayView::RegionType::value_type(0,2));
    region.push_back(DataArrayView::RegionType::value_type(0,2));
    slice=temp.getSlice(region);
    assert(slice.getDataPoint(0,0)(0,0)==0.0);
    assert(slice.getDataPoint(0,0)(1,1)==1.0);
  }

  {
    DataArrayView::ShapeType viewShape;
    Data source(10.0,viewShape,FunctionSpace(),false);
    //
    // weak tests for setting a slice of DataConstant
    cout << "\tTest slicing DataConstant" << endl;
    viewShape.push_back(2);
    viewShape.push_back(3);
    Data target(1.3,viewShape,FunctionSpace(),false);
    DataArrayView::RegionType region;
    region.push_back(DataArrayView::RegionType::value_type(0,0));
    region.push_back(DataArrayView::RegionType::value_type(0,0));
    target.setSlice(source,region);
    assert(target.getDataPoint(0,0)(0,0)==source.getDataPoint(0,0)());
  }

  {
    DataArrayView::ShapeType viewShape;
    Data source(10.0,viewShape,FunctionSpace(),true);
    //
    // weak tests for setting a slice of DataExpanded
    viewShape.push_back(2);
    viewShape.push_back(3);
    Data target(1.3,viewShape,FunctionSpace(),true);
    DataArrayView::RegionType region;
    region.push_back(DataArrayView::RegionType::value_type(0,0));
    region.push_back(DataArrayView::RegionType::value_type(0,0));
    target.setSlice(source,region);
    assert(target.getDataPoint(0,0)(0,0)==source.getDataPoint(0,0)());
  }

  {
    DataArrayView::ShapeType viewShape;
    Data source(10.0,viewShape,FunctionSpace(),false);
    source.tag();
    //
    // weak tests for slicing DataTagged
    cout << "\tTest slicing DataTagged" << endl;
    viewShape.push_back(2);
    viewShape.push_back(3);
    Data target(1.3,viewShape,FunctionSpace(),false);
    //
    // convert the data to tagged
    target.tag();
    DataArrayView::RegionType region;
    region.push_back(DataArrayView::RegionType::value_type(0,0));
    region.push_back(DataArrayView::RegionType::value_type(0,0));
    target.setSlice(source,region);
    assert(target.getDataPoint(0,0)(0,0)==source.getDataPoint(0,0)());
  }

}

void DataTestCase::testMore() {

  cout << endl;

  cout << "\tCreate a Data object from a DataArrayView" << endl;

  DataArrayView::ShapeType viewShape;
  viewShape.push_back(3);
  DataArrayView::ValueType viewData(3);
  for (int i=0;i<viewShape[0];++i) {
    viewData[i]=i;
  }
  DataArrayView myView(viewData,viewShape);

  bool expanded=true;
  Data exData(myView,FunctionSpace(),expanded);
  Data cData(myView);
  Data result;

  assert(exData.isExpanded());
  assert(cData.isConstant());
  assert(result.isEmpty());

  cout << "\tTest some basic operations" << endl;
  result=exData*cData;
  assert(result.isExpanded());

  assert(result.Lsup()==4);
  assert(result.sup()==4);
  assert(result.inf()==0);

  result=exData+cData;
  result=exData-cData;
  result=exData/cData;

  cout << "\tExercise wherePositive method" << endl;
  assert(!exData.wherePositive().isEmpty());

  cout << "\tExercise copyWithMask method" << endl;
  exData.copyWithMask(result, exData.wherePositive());
  assert(!exData.wherePositive().isEmpty());

}

void DataTestCase::testAll() {

  cout << endl;

  cout << "\tCreate a Data object from a DataArrayView" << endl;

  DataArrayView::ShapeType viewShape;
  viewShape.push_back(3);
  DataArrayView::ValueType viewData(3);
  for (int i=0;i<viewShape[0];++i) {
    viewData[i]=i;
  }
  DataArrayView myView(viewData,viewShape);

  bool expanded=true;

  Data exData(myView,FunctionSpace(),expanded);
  Data cData(myView);
  Data result;

  assert(exData.isExpanded());
  assert(cData.isConstant());
  assert(result.isEmpty());

  cout << "\tTest some basic operations" << endl;
  result=exData*cData;
  assert(result.isExpanded());

}

void DataTestCase::testDataConstant() {

  cout << endl;

  cout << "\tCreate a DataConstant object from a DataArrayView" << endl;

  DataArrayView::ShapeType viewShape;
  viewShape.push_back(2);
  viewShape.push_back(3);
  viewShape.push_back(4);
  DataArrayView::ValueType viewData(2*3*4);
  for (int i=0;i<DataArrayView::noValues(viewShape);++i) {
    viewData[i]=i;
  }
  DataArrayView myView(viewData,viewShape);

  Data left(myView);
  Data right(myView);
  Data result;

  cout << "\tTest some basic operations" << endl;

  result=left-right;

  assert(left.isConstant());
  assert(right.isConstant());
  assert(result.isConstant());

  result=left+right;

  assert(left.isConstant());
  assert(right.isConstant());
  assert(result.isConstant());

  assert(!result.isExpanded());
  assert(!result.isTagged());

}

void DataTestCase::testDataTaggedExceptions() {

  cout << endl;

  cout << "\tTest DataTagged operations exceptions." << endl;

  Data myData;
  DataArrayView myView;

  try {
      myData.getSampleDataByTag(0);;
      assert(false);
  }
  catch (EsysException& e) {
      //cout << e.what() << endl;
      assert(true);
  }

  /*
  try {
      myData.setTaggedValue(0,myView);;
      assert(false);
  }
  catch (EsysException& e) {
      //cout << e.what() << endl;
      assert(true);
  }
  */

}

void DataTestCase::testDataTagged() {

  cout << endl;

  cout << "\tCreate a DataTagged object from a DataArrayView" << endl;

  DataTagged::TagListType keys;
  DataTagged::ValueListType values;
  DataArrayView::ShapeType viewShape;
  viewShape.push_back(3);
  DataArrayView::ValueType viewData(3);
  for (int i=0;i<viewShape[0];++i) {
    viewData[i]=i;
  }
  DataArrayView myView(viewData,viewShape);

  // create tagged data with no tag values just a default
  bool expanded=false;

  Data myData(keys,values,myView,FunctionSpace(),expanded);
  assert(myData.isTagged());

  cout << "\tTest some basic operations" << endl;

  Data myDataCopy(myData);
  myDataCopy.expand();
  assert(myDataCopy.isExpanded());

}

void DataTestCase::testConstructors() {

  cout << endl;

  DataArrayView::ShapeType viewShape;
  {
    cout << "\tCreate an Empty Data object" << endl;
    Data temp(1.3,viewShape,FunctionSpace(),false);
  }
  {
    cout << "\tCreate a rank 2 Data object" << endl;
    viewShape.push_back(2);
    viewShape.push_back(3);
    Data temp(1.3,viewShape,FunctionSpace(),false);
  }
}

void DataTestCase::testOperations() {

  cout << endl;

  // define the shape for the DataArrayView test data
  DataArrayView::ShapeType shape;
  shape.push_back(2);
  shape.push_back(3);

  // allocate the data for the DataArrayView
  DataArrayView::ValueType data(DataArrayView::noValues(shape),0);

  // construct DataArrayView
  DataArrayView dataView(data,shape);

  // assign values to the data
  for (int i=0;i<shape[0];i++) {
    for (int j=0;j<shape[1];j++) {
      dataView(i,j)=dataView.index(i,j);
    }
  }

  Data base(dataView);

  // test unary operations

  cout << "\tTest Data::pow." << endl;
  Data power(3.0,shape,FunctionSpace(),true);
  Data result(base.powD(power));
  for (int i=0;i<shape[0];i++) {
    for (int j=0;j<shape[1];j++) {
      assert(result.getPointDataView()(i,j) == pow(dataView.index(i,j),3.0));
    }
  }

  cout << "\tTest Data::sin." << endl;
  result.copy(base.sin());
  assert(true);

  cout << "\tTest Data::cos." << endl;
  result.copy(base.cos());
  assert(true);

  cout << "\tTest Data::tan." << endl;
  result.copy(base.tan());
  assert(true);

  cout << "\tTest Data::asin." << endl;
  result.copy(base.asin());
  assert(true);

  cout << "\tTest Data::acos." << endl;
  result.copy(base.acos());
  assert(true);

  cout << "\tTest Data::atan." << endl;
  result.copy(base.atan());
  assert(true);

  cout << "\tTest Data::sinh." << endl;
  result.copy(base.sinh());
  assert(true);

  cout << "\tTest Data::cosh." << endl;
  result.copy(base.cosh());
  assert(true);

  cout << "\tTest Data::tanh." << endl;
  result.copy(base.tanh());
  assert(true);

  cout << "\tTest Data::asinh." << endl;
  result.copy(base.asinh());
  assert(true);

  cout << "\tTest Data::acosh." << endl;
  result.copy(base.acosh());
  assert(true);

  cout << "\tTest Data::atanh." << endl;
  result.copy(base.atanh());
  assert(true);

  cout << "\tTest Data::log." << endl;
  result.copy(base.log());
  assert(true);

  //cout << "\tTest Data::ln." << endl;
  //result.copy(base.ln());
  //assert(true);

  cout << "\tTest Data::abs." << endl;
  result.copy(base.abs());
  assert(true);

  cout << "\tTest Data::sign." << endl;
  result.copy(base.sign());
  assert(true);

  cout << "\tTest Data::exp." << endl;
  result.copy(base.exp());
  assert(true);

  cout << "\tTest Data::sqrt." << endl;
  result.copy(base.sqrt());
  assert(true);

  cout << "\tTest Data::neg." << endl;
  result.copy(base.neg());
  assert(true);

  cout << "\tTest Data::pos." << endl;
  result.copy(base.pos());
  for (int i=0;i<shape[0];i++) {
    for (int j=0;j<shape[1];j++) {
      assert(result.getPointDataView()(i,j) == dataView.index(i,j));
    }
  }

  // test reduction operations

  cout << "\tTest Data::Lsup." << endl;
  assert(base.Lsup() == 5);

  cout << "\tTest Data::sup." << endl;
  assert(base.sup() == 5);

  cout << "\tTest Data::inf." << endl;
  assert(base.inf() == 0);

  // test data-point reduction operations

  cout << "\tTest Data::minval." << endl;
  result.copy(base.minval());
  assert(result.getPointDataView()() == 0);

  cout << "\tTest Data::maxval." << endl;
  result.copy(base.maxval());
  assert(result.getPointDataView()() == 5);

  //cout << "\tTest Data::length." << endl;
  //result.copy(base.length());
  //assert(pow(result.getPointDataView()(),2.0) == 55);

  cout << "\tTest Data::trace." << endl;
  result.copy(base.trace());
  assert(result.getPointDataView()() == 15);

  //result.copy(base.transpose(0));
  //assert(true);

}

void DataTestCase::testRefValue() {

  //
  // Note - this test can't be run as boost::python::numeric::array
  // objects can only be created and used from within a python thread!
  //

  cout << endl;

  cout << "\tTest Data object RefValue methods." << endl;

  // Create three Data object - DataExpanded, DataConstant and DataEmpty
  DataArrayView::ShapeType viewShape;
  viewShape.push_back(3);
  DataArrayView::ValueType viewData(3);
  for (int i=0;i<viewShape[0];++i) {
    viewData[i]=i;
  }
  DataArrayView myView(viewData,viewShape);

  bool expanded=true;

  Data expandedData(myView,FunctionSpace(),expanded);
  Data constantData(myView);
  Data emptyData;

  assert(expandedData.isExpanded());
  assert(constantData.isConstant());
  assert(emptyData.isEmpty());

  // Check assertions are thrown for RefValue methods on DataEmpty

  int ref = 0;
  boost::python::numeric::array num_array(1.0);

  try {
      emptyData.getRefValue(ref,num_array);
      assert(false);
  }
  catch (EsysException& e) {
      assert(true);
  }
  try {
      emptyData.setRefValue(ref,num_array);
      assert(false);
  }
  catch (EsysException& e) {
      assert(true);
  }

  // Check assertions are thrown for RefValue methods on DataConstant
  try {
      constantData.getRefValue(ref,num_array);
      assert(false);
  }
  catch (EsysException& e) {
      assert(true);
  }
  try {
      constantData.setRefValue(ref,num_array);
      assert(false);
  }
  catch (EsysException& e) {
      assert(true);
  }

  // Check calls to RefValue methods on DataExpanded
  expandedData.getRefValue(ref,num_array);
  expandedData.setRefValue(ref,num_array);

}

void DataTestCase::testMemAlloc() {

  //
  // Simple little sanity check for the memory allocator

  cout << endl;

  Data *testData;
  for (int i=0; i<1000; i++) {
    testData = new Data(0.0, DataArrayView::ShapeType(), FunctionSpace(), true);
    delete testData;
  }

  DataArrayView::ShapeType viewShape;
  viewShape.push_back(10);
  viewShape.push_back(10);
  viewShape.push_back(10);

  Data *testData2;
  Data *testData3 = new Data(0.0, viewShape, FunctionSpace(), true);
  for (int i=0; i<1000; i++) {
    testData2 = new Data(0.0, viewShape, FunctionSpace(), true);
    delete testData2;
  }
  delete testData3;

}

TestSuite* DataTestCase::suite ()
{
  //
  // create the suite of tests to perform.
  TestSuite *testSuite = new TestSuite ("DataTestCase");

  testSuite->addTest (new TestCaller< DataTestCase>("testAll",&DataTestCase::testAll));
  testSuite->addTest (new TestCaller< DataTestCase>("testMore",&DataTestCase::testMore));
  testSuite->addTest (new TestCaller< DataTestCase>("testDataConstant",&DataTestCase::testDataConstant));
  testSuite->addTest (new TestCaller< DataTestCase>("testDataTagged",&DataTestCase::testDataTagged));
  testSuite->addTest (new TestCaller< DataTestCase>("testDataTaggedExceptions",&DataTestCase::testDataTaggedExceptions));
  testSuite->addTest (new TestCaller< DataTestCase>("testConstructors",&DataTestCase::testConstructors));
  testSuite->addTest (new TestCaller< DataTestCase>("testSlicing",&DataTestCase::testSlicing));
  testSuite->addTest (new TestCaller< DataTestCase>("testOperations",&DataTestCase::testOperations));
  //testSuite->addTest (new TestCaller< DataTestCase>("testRefValue",&DataTestCase::testRefValue));
  testSuite->addTest (new TestCaller< DataTestCase>("testMemAlloc",&DataTestCase::testMemAlloc));

  return testSuite;
}
