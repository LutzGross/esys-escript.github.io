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
#include "escript/Data/DataTagged.h"
#include "escript/Data/BinaryOp.h"
#include "escript/Data/UnaryOp.h"
#include "esysUtils/EsysException.h"

/*
#include "finley/CPPAdapter/MeshAdapter.h"
#include "finley/CPPAdapter/MeshAdapterFactory.h"
#include "escript/Data/AbstractContinuousDomain.h"
*/

#include "escript/Data/FunctionSpaceFactory.h"
#include "escript/Data/DataFactory.h"

#include "DataTaggedTestCase.h"

#include <iostream>
#include <functional>
#include <algorithm>

using namespace CppUnitTest;
using namespace escript;
using namespace esysUtils;
using namespace std;

//using namespace finley;

void DataTaggedTestCase::setUp() {
  //
  // This is called before each test is run
 
}

void DataTaggedTestCase::tearDown() {
  //
  // This is called after each test has been run
 
}

void DataTaggedTestCase::testReshape() {

  cout << endl;

  {
    cout << "\tTest reshape of default constructed DataTagged to rank 1." << endl;
    DataTagged value;
    value.getPointDataView()()=1.0;
    DataArrayView::ShapeType shape;
    shape.push_back(2);
    value.reshapeDataPoint(shape);
    for (int i=0;i<shape[0];++i) {
      assert(value.getDefaultValue()(i)==1);
    }
  }
  {
    cout << "\tTest reshape of default constructed DataTagged to rank 2." << endl;
    DataTagged value;
    value.getPointDataView()()=0.0;
    DataArray vOne(1.0);
    DataArray vTwo(2.0);
    value.addTaggedValue(1,vOne.getView());
    value.addTaggedValue(2,vTwo.getView());
    DataArrayView::ShapeType shape;
    shape.push_back(2);
    shape.push_back(5);
    value.reshapeDataPoint(shape);
    for (int j=0;j<shape[1];++j) {
      for (int i=0;i<shape[0];++i) {
	assert(value.getDefaultValue()(i,j)==0.0);
	assert(value.getDataPointByTag(1)(i,j)==vOne.getView()());
	assert(value.getDataPointByTag(2)(i,j)==vTwo.getView()());
      }
    }
  }
}

void DataTaggedTestCase::testOperations() {

  cout << endl;

  {
    cout << "\tTest default DataTagged contains only a default value." << endl;
    DataTagged left;
    DataTagged right;
    binaryOp(left,right,plus<double>());
    assert(left.getPointDataView()()==0);
    assert(right.getPointDataView()()==0);
    cout << "\tTest addTaggedValue and binaryOp(plus)." << endl;
    DataArray vOne(1.0);
    DataArray vTwo(2.0);
    right.addTaggedValue(1,vOne.getView());
    right.addTaggedValue(2,vTwo.getView());
    binaryOp(left,right,plus<double>());
    assert(left.getPointDataView()()==0);
    assert(left.getDataPointByTag(1)==vOne.getView());
    assert(left.getDataPointByTag(2)==vTwo.getView());
    cout << "\tTest setTaggedValue and binaryOp(multiplies)." << endl;
    DataArray vZero(0.0);
    right.setTaggedValue(1,vZero.getView());
    right.setTaggedValue(2,vZero.getView());
    binaryOp(left,right,multiplies<double>());
    assert(left.getPointDataView()()==0);
    assert(left.getDataPointByTag(1)==vZero.getView());
    assert(left.getDataPointByTag(2)==vZero.getView());
  }
  {
    DataArrayView::ValueType viewData;
    DataArrayView::ShapeType viewShape;
    viewShape.push_back(3);
    DataTagged::TagListType keys;
    DataTagged::ValueListType values;
    for (int i=0;i<viewShape[0];++i) {
      viewData.push_back(i);
    }
    DataArrayView myView(viewData,viewShape);
    cout << "\tCreate tagged data with no tag values just a default." << endl;
    DataTagged left(keys,values,myView,FunctionSpace());
    DataTagged right(keys,values,myView,FunctionSpace());
    binaryOp(left,right,minus<double>());
    for (int i=0;i<viewShape[0];++i) {
      assert(left.getDefaultValue()(i)==0);
    }
    double mVal=10.0;
    for (int i=0;i<viewShape[0];++i) {
      viewData[i]=i*mVal;
    }
    cout << "\tTest binaryOp(minus) with the right hand side as a view." << endl;
    binaryOp(left,myView,minus<double>());
    for (int i=0;i<viewShape[0];++i) {
      assert(left.getDefaultValue()(i)==-(i*mVal));
    }
  }
  {
    cout << "\tTest unaryOp(negate) on default DataTagged." << endl;
    DataTagged data;
    unaryOp(data,negate<double>());
    assert(data.getDefaultValue()()==0);
    DataArray vOne(1);
    binaryOp(data,vOne.getView(),plus<double>());
    assert(data.getDefaultValue()()==1);
    unaryOp(data,negate<double>());
    assert(data.getDefaultValue()()==-1);
  }
  {
    cout << "\tTest unaryOp(negate) on DataTagged with 3 tags." << endl;
    DataArrayView::ShapeType vShape;
    vShape.push_back(3);
    vShape.push_back(2);
    vShape.push_back(1);
    DataArray defData(vShape,0.0);
    DataArrayView& defView=defData.getView();
    DataArray tOneData(vShape,1.0);
    DataArrayView& tOneView=tOneData.getView();
    DataArray tTwoData(vShape,2.0);
    DataArrayView& tTwoView=tTwoData.getView();
    DataArray tThreeData(vShape,3.0);
    DataArrayView& tThreeView=tThreeData.getView();
    DataTagged::TagListType keys;
    DataTagged::ValueListType values;
    keys.push_back(1);
    keys.push_back(2);
    keys.push_back(3);
    values.push_back(tOneView);
    values.push_back(tTwoView);
    values.push_back(tThreeView);
    DataTagged tData(keys,values,defView,FunctionSpace());
    unaryOp(tData,negate<double>());
    unaryOp(tData,negate<double>());
    assert(tData.getDataPointByTag(1)==tOneView);
    assert(tData.getDataPointByTag(2)==tTwoView);
    assert(tData.getDataPointByTag(3)==tThreeView);
  }

}

void DataTaggedTestCase::testAll() {
  cout << endl;
  {
    cout << "\tTest default construction." << endl;
    DataTagged myData;
    assert(myData.getPointDataView()()==0);
    assert(myData.getNumDPPSample()==1);
    assert(myData.getNumSamples()==1);
    cout << "\tTest adding two keys with empty value list." << endl;
    DataTagged::TagListType keys;
    DataTagged::ValueListType values;
    keys.push_back(1);
    keys.push_back(2);
    myData.addTaggedValues(keys,values);
    for (int i=0;i<keys.size();++i) {
      assert(myData.getPointDataView()()==0);
    }
  }
  {
    DataArrayView::ValueType viewData;
    DataArrayView::ShapeType viewShape;
    viewShape.push_back(3);
    DataTagged::TagListType keys;
    DataTagged::ValueListType values;
    for (int i=0;i<viewShape[0];++i) {
      viewData.push_back(0.0);
    }
    DataArrayView myView(viewData,viewShape);
    cout << "\tCreate tagged data with no tag values just a default." << endl;
    DataTagged myData(keys,values,myView,FunctionSpace());
    assert(myData.getNumDPPSample()==1);
    assert(myData.getNumSamples()==1);
    cout << "\tTest non existent tag returns the default value." << endl;
    assert(myData.getDataPointByTag(1)==myView);
    cout << "\tTest adding a single tag value." << endl;
    for (int i=0;i<myView.getShape()[0];++i) {
      myView(i)=i;
    }
    values.push_back(myView);
    keys.push_back(1);
    myData.addTaggedValues(keys,values);
    assert(myData.getDataPointByTag(1)==myView);
    cout << "\tTest addition of further tags." << endl;
    keys.clear();
    keys.push_back(3);
    for (int i=0;i<myView.getShape()[0];++i) {
      myView(i)=i+1.5;
    }
    myData.addTaggedValues(keys,values);
    assert(myData.getDataPointByTag(3)==myView);
    assert(myData.getDataPointByTag(1)!=myView);
    cout << "\tTrigger the size mismatch exception." << endl;
    try {
      values.push_back(myView);
      myData.addTaggedValues(keys,values);
      assert(false);
    }
    catch (EsysException& e) {
      //cout << e.what() << endl;
      assert(true);
    }
  }
  {
    cout << "\tTest creation of tagged data with multiple tags." << endl;
    DataArrayView::ValueType viewData;
    DataArrayView::ShapeType viewShape;
    viewShape.push_back(3);
    DataTagged::TagListType keys;
    DataTagged::ValueListType values;
    for (int i=0;i<viewShape[0];++i) {
      viewData.push_back(0.0);
    }
    DataArrayView myView(viewData,viewShape);
    DataArray eOne(myView);
    DataArray eTwo(myView);
    DataArray eThree(myView);
    for (int i=0;i<eOne.getView().getShape()[0];++i) {
      eOne.getView()(i)=i+1.0;
    }
    for (int i=0;i<eTwo.getView().getShape()[0];++i) {
      eTwo.getView()(i)=i+2.0;
    }
    for (int i=0;i<eThree.getView().getShape()[0];++i) {
      eThree.getView()(i)=i+3.0;
    }
    values.push_back(eOne.getView());
    values.push_back(eTwo.getView());
    values.push_back(eThree.getView());
    keys.push_back(1);
    keys.push_back(2);
    keys.push_back(3);
    DataTagged myData(keys,values,myView,FunctionSpace());
    assert(myData.getDataPointByTag(1)==eOne.getView());
    assert(myData.getDataPointByTag(2)==eTwo.getView());
    assert(myData.getDataPointByTag(3)==eThree.getView());
    cout << "\tTest isCurrentTag function." << endl;
    for (int i=0;i<keys.size();++i) {
      assert(myData.isCurrentTag(keys[i]));
    }
    cout << "\tCheck correct operation for key that doesn't exist." << endl;
    assert(!myData.isCurrentTag(123));
    cout << "\tTrigger bad shape in input values exception." << endl;
    viewShape.clear();
    viewShape.push_back(1);
    keys.clear();
    values.clear();
    viewData.clear();
    for (int i=0;i<viewShape[0];++i) {
      viewData.push_back(0.0);
    }
    DataArrayView myView2(viewData,viewShape);
    try {
      myData.addTaggedValue(5,myView2);
      assert(false);
    }
    catch (EsysException& e) {
      //cout << e.what() << endl;
      assert(true);
    }
    cout << "\tTest setTaggedValues." << endl;
    DataTagged myData2;
    myData2.reshapeDataPoint(myView.getShape());
    keys.clear();
    values.clear();
    keys.push_back(1);
    keys.push_back(2);
    keys.push_back(3);
    values.push_back(eOne.getView());
    values.push_back(eTwo.getView());
    values.push_back(eThree.getView());
    myData2.setTaggedValues(keys,values);
    assert(myData2.getDataPointByTag(1)==eOne.getView());
    assert(myData2.getDataPointByTag(2)==eTwo.getView());
    assert(myData2.getDataPointByTag(3)==eThree.getView());
    cout << "\tTest setTaggedValue." << endl;
    DataTagged myData3;
    myData3.reshapeDataPoint(myView.getShape());
    myData3.setTaggedValue(1,eOne.getView());
    myData3.setTaggedValue(2,eTwo.getView());
    myData3.setTaggedValue(3,eThree.getView());
    assert(myData3.getDataPointByTag(1)==eOne.getView());
    assert(myData3.getDataPointByTag(2)==eTwo.getView());
    assert(myData3.getDataPointByTag(3)==eThree.getView());
  }

}

void DataTaggedTestCase::testSubtraction() {

  // An error with FinleyMesh::getTagList used to cause binary operations
  // between DataExpanded and DataTagged objects to seg-fault. This test
  // case will provoke this error if it arises again.

  // This test requires the version of setTaggedData which takes a Data
  // object as an argument. This version is not currently available, so
  // this test is disabled for now

  /*

  cout << endl;

  cout << "\tCreate domain and function-space." << endl;
  AbstractContinuousDomain* myDomain = rectangle(10,10);
  FunctionSpace f = functionOnBoundary(*myDomain);

  cout << "\tCreate two vectors, one being DataExpanded." << endl;
  Data A = Vector(0,f);
  Data B = Vector(0,f,true);

  cout << "\tCreate some tags and values to add to the other." << endl;
  DataArrayView::ValueType viewData;
  DataArrayView::ShapeType viewShape;
  viewShape.push_back(2);
  for (int i=0;i<viewShape[0];++i) {
    viewData.push_back(0.0);
  }
  DataArrayView myView(viewData,viewShape);
  DataArray eOne(myView);
  for (int i=0;i<eOne.getView().getShape()[0];++i) {
    eOne.getView()(i)=i+1.0;
  }

  A.setTaggedValue(2,eOne.getView());

  //cout << A.toString() << endl;

  cout << "\tCalculate difference." << endl;
  Data difference = B - A;

  cout << "\tCalculate other binaryOps just to be sure." << endl;
  Data sum = B + A;
  Data product = B * A;
  Data dividend = B / A;

  // If we get here, subtraction operation did not seg-fault.
  assert(true);

  */

}

TestSuite* DataTaggedTestCase::suite ()
{
  //
  // create the suite of tests to perform.
  TestSuite *testSuite = new TestSuite ("DataTaggedTestCase");
  testSuite->addTest (new TestCaller< DataTaggedTestCase>("testAll",&DataTaggedTestCase::testAll));
  testSuite->addTest (new TestCaller< DataTaggedTestCase>("testOperations",&DataTaggedTestCase::testOperations));
  testSuite->addTest (new TestCaller< DataTaggedTestCase>("testReshape",&DataTaggedTestCase::testReshape));
  //testSuite->addTest (new TestCaller< DataTaggedTestCase>("testSubtraction",&DataTaggedTestCase::testSubtraction));
  return testSuite;
}
