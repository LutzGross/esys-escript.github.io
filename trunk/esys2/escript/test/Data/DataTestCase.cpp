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
    // weak tests for setting a slice of DataConstant
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

  DataArrayView::ValueType viewData;
  DataArrayView::ShapeType viewShape;
  viewShape.push_back(3);
  for (int i=0;i<viewShape[0];++i) {
    viewData.push_back(i);
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
  //assert(exData.wherePositive()==exData.wherePositive());

  cout << "\tExercise copyWithMask method" << endl;
  exData.copyWithMask(result, exData.wherePositive());
  assert(!exData.wherePositive().isEmpty());

}

void DataTestCase::testAll() {

  cout << endl;

  cout << "\tCreate a Data object from a DataArrayView" << endl;

  DataArrayView::ValueType viewData;
  DataArrayView::ShapeType viewShape;
  viewShape.push_back(3);
  for (int i=0;i<viewShape[0];++i) {
    viewData.push_back(i);
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

  DataArrayView::ValueType viewData;
  DataArrayView::ShapeType viewShape;
  viewShape.push_back(2);
  viewShape.push_back(3);
  viewShape.push_back(4);
  for (int i=0;i<DataArrayView::noValues(viewShape);++i) {
    viewData.push_back(i);
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
  DataArrayView::ValueType viewData;
  DataArrayView::ShapeType viewShape;
  viewShape.push_back(3);
  for (int i=0;i<viewShape[0];++i) {
    viewData.push_back(i);
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
    cout << "\tDump it toString:" << endl;
    cout << temp.toString() << endl;
  }
}

void DataTestCase::testOperations() {
  cout << endl;

  cout << "\tCreate a rank 2 Data object" << endl;
  DataArrayView::ShapeType viewShape;
  viewShape.push_back(2);
  viewShape.push_back(3);

  Data base(2.0,viewShape,FunctionSpace(),false);
  Data power(3.0,viewShape,FunctionSpace(),false);

  cout << "\tPerform basic exercises of unary operations" << endl;

  Data result(base.powD(power));
  assert(result.getDataPoint(0,0)(0,0) == 8);

  result.copy(base.sin());
  assert(true);

  result.copy(base.cos());
  assert(true);

  result.copy(base.tan());
  assert(true);

  result.copy(base.log());
  assert(true);

  result.copy(base.ln());
  assert(true);

  result.copy(base.abs());
  assert(true);

  result.copy(base.maxval());
  assert(true);

  result.copy(base.minval());
  assert(true);

  result.copy(base.length());
  assert(true);

  result.copy(base.sign());
  assert(true);

  result.copy(base.transpose(0));
  assert(true);

  result.copy(base.trace());
  assert(true);

  result.copy(base.exp());
  assert(true);

  result.copy(base.sqrt());
  assert(true);

  result.copy(base.neg());
  assert(true);

  result.copy(base.pos());
  assert(true);
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

  return testSuite;
}
