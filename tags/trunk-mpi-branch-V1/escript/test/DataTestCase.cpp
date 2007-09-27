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
#include <iostream>
#if (defined _WIN32) && (defined __INTEL_COMPILER)
#include <mathimf.h>
#else
#include <math.h>
#endif

#include "DataTestCase.h"

#include "escript/FunctionSpace.h"
#include "esysUtils/EsysException.h"

#include "escript/Data.h"


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

    cout << "\tTest get-slicing DataConstant" << endl;

    DataArrayView::ShapeType viewShape;
    viewShape.push_back(2);
    viewShape.push_back(3);
    Data data(1.3,viewShape,FunctionSpace(),false);

    //cout << data.toString() << endl;

    DataArrayView::RegionType region;
    region.push_back(DataArrayView::RegionType::value_type(0,0));
    region.push_back(DataArrayView::RegionType::value_type(0,0));

    Data slice1(data.getSlice(region));

    //cout << slice1.toString() << endl;

    assert(slice1.getDataPointRank()==0);
    assert(slice1.getDataPoint(0,0)()==1.3);

    region.clear();
    region.push_back(DataArrayView::RegionType::value_type(0,1));
    region.push_back(DataArrayView::RegionType::value_type(0,1));

    Data slice2(data.getSlice(region));

    //cout << slice2.toString() << endl;

    assert(slice2.getDataPointRank()==2);
    assert(slice2.getDataPoint(0,0)(0,0)==1.3);

    region.clear();
    region.push_back(DataArrayView::RegionType::value_type(0,1));
    region.push_back(DataArrayView::RegionType::value_type(0,2));

    Data slice3(data.getSlice(region));

    //cout << slice3.toString() << endl;

    assert(slice3.getDataPointRank()==2);
    assert(slice3.getDataPoint(0,0)(0,0)==1.3);
    assert(slice3.getDataPoint(0,0)(0,1)==1.3);

  }

  {

    cout << "\tTest set-slicing DataConstant" << endl;

    DataArrayView::ShapeType viewShape;
    Data source(10.0,viewShape,FunctionSpace(),false);

    //cout << source.toString() << endl;

    viewShape.push_back(2);
    viewShape.push_back(3);
    Data target(1.3,viewShape,FunctionSpace(),false);

    //cout << target.toString() << endl;

    DataArrayView::RegionType region;
    region.push_back(DataArrayView::RegionType::value_type(0,0));
    region.push_back(DataArrayView::RegionType::value_type(0,0));

    target.setSlice(source,region);

    //cout << target.toString() << endl;

    assert(target.getDataPoint(0,0)(0,0)==source.getDataPoint(0,0)());

  }

  {

    cout << "\tTest get-slicing DataTagged" << endl;

    //
    // create a DataTagged with a default value only

    DataArrayView::ShapeType viewShape;
    viewShape.push_back(2);
    viewShape.push_back(3);
    Data data(1.3,viewShape,FunctionSpace(),false);
    data.tag();
    data.getDataPoint(0,0)(0,0)=1.0;
    data.getDataPoint(0,0)(1,1)=2.0;

    //cout << data.toString() << endl;

    //
    // create a scalar slice

    DataArrayView::RegionType region;
    region.push_back(DataArrayView::RegionType::value_type(0,0));
    region.push_back(DataArrayView::RegionType::value_type(0,0));

    Data slice1(data.getSlice(region));

    //cout << slice1.toString() << endl;

    assert(slice1.isTagged());
    assert(slice1.getDataPointRank()==0);
    assert(slice1.getDataPoint(0,0)()==1.0);

    //
    // create a rank 2 slice with one value

    region.clear();
    region.push_back(DataArrayView::RegionType::value_type(0,1));
    region.push_back(DataArrayView::RegionType::value_type(0,1));

    Data slice2(data.getSlice(region));

    //cout << slice2.toString() << endl;

    assert(slice2.isTagged());
    assert(slice2.getDataPointRank()==2);
    assert(slice2.getDataPoint(0,0)(0,0)==1.0);

    //
    // create a rank 2 slice with four values

    region.clear();
    region.push_back(DataArrayView::RegionType::value_type(0,2));
    region.push_back(DataArrayView::RegionType::value_type(0,2));

    Data slice3(data.getSlice(region));

    //cout << slice3.toString() << endl;

    assert(slice3.isTagged());
    assert(slice3.getDataPointRank()==2);
    assert(slice3.getDataPoint(0,0)(0,0)==1.0);
    assert(slice3.getDataPoint(0,0)(0,1)==1.3);
    assert(slice3.getDataPoint(0,0)(1,0)==1.3);
    assert(slice3.getDataPoint(0,0)(1,1)==2.0);

    //
    // add a value for tag "1"

    DataArrayView::ValueType viewData(6);
    for (int i=0;i<viewData.size();i++) {
      viewData[i]=i;
    }
    DataArrayView dataView(viewData,viewShape);

    data.setTaggedValueFromCPP(1, dataView);

    //
    // create a full slice

    region.clear();
    region.push_back(DataArrayView::RegionType::value_type(0,2));
    region.push_back(DataArrayView::RegionType::value_type(0,3));

    Data slice4(data.getSlice(region));

    //cout << slice4.toString() << endl;

    assert(slice4.isTagged());
    assert(slice4.getDataPointRank()==2);
    assert(slice4.getDataPoint(0,0)(0,0)==0);
    assert(slice4.getDataPoint(0,0)(0,1)==2);
    assert(slice4.getDataPoint(0,0)(0,2)==4);
    assert(slice4.getDataPoint(0,0)(1,0)==1);
    assert(slice4.getDataPoint(0,0)(1,1)==3);
    assert(slice4.getDataPoint(0,0)(1,2)==5);

  }

  {

    cout << "\tTest set-slicing DataTagged" << endl;

    //
    // create a source DataTagged with a scalar default value only

    DataArrayView::ShapeType viewShape;
    Data source(10.0,viewShape,FunctionSpace(),false);
    source.tag();

    //cout << "source:\n" << source.toString() << endl;

    //
    // create a target DataTagged with a rank 2 default value only

    viewShape.push_back(2);
    viewShape.push_back(3);
    Data target(1.3,viewShape,FunctionSpace(),false);
    target.tag();

    //cout << "target:\n" << target.toString() << endl;

    //
    // set a slice in target from source

    DataArrayView::RegionType region;
    region.push_back(DataArrayView::RegionType::value_type(1,1));
    region.push_back(DataArrayView::RegionType::value_type(1,1));

    target.setSlice(source,region);

    //cout << "target:\n" << target.toString() << endl;

    assert(target.isTagged());
    assert(target.getDataPointRank()==2);
    assert(target.getDataPoint(0,0)(0,0)==1.3);
    assert(target.getDataPoint(0,0)(0,1)==1.3);
    assert(target.getDataPoint(0,0)(0,2)==1.3);
    assert(target.getDataPoint(0,0)(1,0)==1.3);
    assert(target.getDataPoint(0,0)(1,1)==source.getDataPoint(0,0)());
    assert(target.getDataPoint(0,0)(1,2)==1.3);

    //
    // add a value for tag "1" to target

    DataArrayView::ValueType viewData(6);
    for (int i=0;i<viewData.size();i++) {
      viewData[i]=i;
    }
    DataArrayView dataView(viewData,viewShape);

    target.setTaggedValueFromCPP(1, dataView);

    //cout << "target:\n" << target.toString() << endl;

    //
    // set a slice in target from source

    region.clear();
    region.push_back(DataArrayView::RegionType::value_type(0,0));
    region.push_back(DataArrayView::RegionType::value_type(1,1));

    target.setSlice(source,region);

    //cout << "target:\n" << target.toString() << endl;

    assert(target.isTagged());
    assert(target.getDataPointRank()==2);
    assert(target.getDataPoint(0,0)(0,0)==0);
    assert(target.getDataPoint(0,0)(0,1)==source.getDataPoint(0,0)());
    assert(target.getDataPoint(0,0)(0,2)==4);
    assert(target.getDataPoint(0,0)(1,0)==1);
    assert(target.getDataPoint(0,0)(1,1)==3);
    assert(target.getDataPoint(0,0)(1,2)==5);

    //
    // add a value for tag "2" to source

    DataArrayView::ShapeType viewShape2;
    DataArrayView::ValueType viewData2(1);
    viewData2[0]=6;
    DataArrayView dataView2(viewData2,viewShape2);

    source.setTaggedValueFromCPP(2, dataView2);

    //cout << "source:\n" << source.toString() << endl;

    //
    // set a slice in target from source

    region.clear();
    region.push_back(DataArrayView::RegionType::value_type(0,0));
    region.push_back(DataArrayView::RegionType::value_type(1,1));

    target.setSlice(source,region);

    //cout << "target:\n" << target.toString() << endl;

    assert(target.isTagged());
    assert(target.getDataPointRank()==2);

    // use a non-existant tag so we get a pointer to the default value
    // ie: the first element in the data array
    DataAbstract::ValueType::value_type* targetData=target.getSampleDataByTag(9);
    for (int i=0; i<target.getLength(); i++) {
      assert(targetData[i]>=0);
    }
    assert(targetData[0]==1.3);
    assert(targetData[1]==1.3);
    assert(targetData[2]==10);
    assert(targetData[3]==10);
    assert(targetData[4]==1.3);
    assert(targetData[5]==1.3);
    assert(targetData[6]==0);
    assert(targetData[7]==1);
    assert(targetData[8]==10);
    assert(targetData[9]==3);
    assert(targetData[10]==4);
    assert(targetData[11]==5);
    assert(targetData[12]==1.3);
    assert(targetData[13]==1.3);
    assert(targetData[14]==6);
    assert(targetData[15]==10);
    assert(targetData[16]==1.3);
    assert(targetData[17]==1.3);

  }

  {

    cout << "\tTest get-slicing DataExpanded" << endl;

    DataArrayView::ShapeType viewShape;
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

    cout << "\tTest set-slicing DataExpanded" << endl;

    DataArrayView::ShapeType viewShape;
    Data source(10.0,viewShape,FunctionSpace(),true);

    viewShape.push_back(2);
    viewShape.push_back(3);
    Data target(1.3,viewShape,FunctionSpace(),true);

    DataArrayView::RegionType region;
    region.push_back(DataArrayView::RegionType::value_type(0,0));
    region.push_back(DataArrayView::RegionType::value_type(0,0));

    target.setSlice(source,region);

    assert(target.getDataPoint(0,0)(0,0)==source.getDataPoint(0,0)());

  }

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

void DataTestCase::testDataTagged() {

  cout << endl;

  {

    cout << "\tCreate a DataTagged object with a default value only." << endl;

    DataTagged::TagListType keys;

    DataTagged::ValueListType values;

    DataArrayView::ShapeType viewShape;
    viewShape.push_back(3);

    DataArrayView::ValueType viewData(3);
    for (int i=0;i<viewShape[0];i++) {
      viewData[i]=i;
    }
    DataArrayView defaultValue(viewData,viewShape);

    bool expanded=false;
  
    Data myData(keys,values,defaultValue,FunctionSpace(),expanded);

    // cout << myData.toString() << endl;

    assert(!myData.isEmpty());
    assert(myData.isTagged());
    assert(myData.getTagNumber(0)==1);
    assert(myData.getDataPointRank()==1);
    assert(myData.getLength()==3);

    DataArrayView myDataView = myData.getPointDataView();
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==0);
    assert(myDataView.getRank()==1);
    assert(myDataView.noValues()==3);
    assert(myDataView.getShape().size()==1);
    assert(myDataView(0)==0.0);
    assert(myDataView(1)==1.0);
    assert(myDataView(2)==2.0);

    myDataView = myData.getDataPoint(0,0);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==0);
    assert(myDataView.getRank()==1);
    assert(myDataView.noValues()==3);
    assert(myDataView.getShape().size()==1);
    assert(myDataView(0)==0.0);
    assert(myDataView(1)==1.0);
    assert(myDataView(2)==2.0);

    double* sampleData=myData.getSampleData(0);
    for (int i=0; i<myDataView.noValues(); i++) {
      assert(sampleData[i]==i);
    }
    // use a non-existent tag so we get a pointer to
    // the first element of the data array
    sampleData=myData.getSampleDataByTag(9);
    for (int i=0; i<myData.getLength(); i++) {
      assert(sampleData[i]==i);
    }

    cout << "\tTest setting of a tag and associated value." << endl;

    // value for tag "1"
    DataArray eTwo(defaultValue);
    for (int i=0;i<eTwo.getView().getShape()[0];i++) {
      eTwo.getView()(i)=i+2.0;
    }

    myData.setTaggedValueFromCPP(1,eTwo.getView());

    assert(myData.getLength()==6);

    myDataView = myData.getDataPoint(0,0);
    assert(myDataView==eTwo.getView());
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==3);
    assert(myDataView.getRank()==1);
    assert(myDataView.noValues()==3);
    assert(myDataView.getShape().size()==1);
    assert(myDataView(0)==2);
    assert(myDataView(1)==3);
    assert(myDataView(2)==4);

    sampleData=myData.getSampleDataByTag(1);
    for (int i=0; i<myDataView.noValues(); i++) {
      assert(sampleData[i]==i+2);
    }

  }

  {

    cout << "\tCreate a DataTagged object via tag() method." << endl;

    DataArrayView::ShapeType viewShape;
    viewShape.push_back(2);
    viewShape.push_back(3);
    Data myData(1.3,viewShape,FunctionSpace(),false);
    myData.tag();

    //cout << myData.toString() << endl;

    assert(!myData.isEmpty());
    assert(myData.isTagged());
    assert(myData.getTagNumber(0)==1);
    assert(myData.getDataPointRank()==2);
    assert(myData.getLength()==6);

    // check default value
    DataArrayView myDataView = myData.getPointDataView();
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==0);
    assert(myDataView.getRank()==2);
    assert(myDataView.noValues()==6);
    assert(myDataView.getShape().size()==2);
    assert(myDataView(0,0)==1.3);
    assert(myDataView(0,1)==1.3);
    assert(myDataView(0,2)==1.3);
    assert(myDataView(1,0)==1.3);
    assert(myDataView(1,1)==1.3);
    assert(myDataView(1,2)==1.3);

    // check value for data-point (0,0).
    myDataView = myData.getDataPoint(0,0);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==0);
    assert(myDataView.getRank()==2);
    assert(myDataView.noValues()==6);
    assert(myDataView.getShape().size()==2);
    assert(myDataView(0,0)==1.3);
    assert(myDataView(0,1)==1.3);
    assert(myDataView(0,2)==1.3);
    assert(myDataView(1,0)==1.3);
    assert(myDataView(1,1)==1.3);
    assert(myDataView(1,2)==1.3);

  }

}

void DataTestCase::testDataTaggedExceptions() {

  cout << endl;

  cout << "\tTest DataTagged exceptions." << endl;

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

  try {
      myData.setTaggedValueFromCPP(0,myView);;
      assert(false);
  }
  catch (EsysException& e) {
      //cout << e.what() << endl;
      assert(true);
  }

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

  Data baseEx(dataView,FunctionSpace(),true);
  Data baseCon(dataView,FunctionSpace(),false);
  Data baseTag(dataView,FunctionSpace(),false);
  baseTag.tag();

  assert(baseEx.isExpanded());
  assert(baseCon.isConstant());
  assert(baseTag.isTagged());

  Data resultEx;
  Data resultCon;
  Data resultTag;

  // test unary operations

  double tmp;
  cout << "\tTest Data::pow." << endl;
  Data power(3.0,shape,FunctionSpace(),true);
  resultEx.copy(baseEx.powD(power));
  resultCon.copy(baseCon.powD(power));
  resultTag.copy(baseTag.powD(power));
  for (int i=0;i<shape[0];i++) {
    for (int j=0;j<shape[1];j++) {
      tmp=pow(dataView.index(i,j),3.0);
      assert(std::abs(resultEx.getPointDataView()(i,j) - tmp) <= REL_TOL*std::abs(tmp));
      assert(std::abs(resultCon.getPointDataView()(i,j)- tmp) <= REL_TOL*std::abs(tmp));
      assert(std::abs(resultTag.getPointDataView()(i,j)- tmp) <= REL_TOL*std::abs(tmp));
    }
  }

  cout << "\tTest Data::sin." << endl;
  resultEx.copy(baseEx.sin());
  resultCon.copy(baseCon.sin());
  resultTag.copy(baseTag.sin());
  for (int i=0;i<shape[0];i++) {
    for (int j=0;j<shape[1];j++) {
      tmp=sin((double)dataView.index(i,j));
      assert(std::abs(resultEx.getPointDataView()(i,j) - tmp) <= REL_TOL*std::abs(tmp));
      assert(std::abs(resultCon.getPointDataView()(i,j)- tmp) <= REL_TOL*std::abs(tmp));
      assert(std::abs(resultTag.getPointDataView()(i,j)- tmp) <= REL_TOL*std::abs(tmp));
    }
  }

  cout << "\tTest Data::cos." << endl;
  resultEx.copy(baseEx.cos());
  resultCon.copy(baseCon.cos());
  resultTag.copy(baseTag.cos());
  for (int i=0;i<shape[0];i++) {
    for (int j=0;j<shape[1];j++) {
      tmp=cos((double)dataView.index(i,j));
      assert(std::abs(resultEx.getPointDataView()(i,j) - tmp) <= REL_TOL*std::abs(tmp));
      assert(std::abs(resultCon.getPointDataView()(i,j)- tmp) <= REL_TOL*std::abs(tmp));
      assert(std::abs(resultTag.getPointDataView()(i,j)- tmp) <= REL_TOL*std::abs(tmp));
    }
  }

  cout << "\tTest Data::tan." << endl;
  resultEx.copy(baseEx.tan());
  resultCon.copy(baseCon.tan());
  resultTag.copy(baseTag.tan());
  for (int i=0;i<shape[0];i++) {
    for (int j=0;j<shape[1];j++) {
      tmp=tan((double)dataView.index(i,j));
      assert(std::abs(resultEx.getPointDataView()(i,j) - tmp) <= REL_TOL*std::abs(tmp));
      assert(std::abs(resultCon.getPointDataView()(i,j)- tmp) <= REL_TOL*std::abs(tmp));
      assert(std::abs(resultTag.getPointDataView()(i,j)- tmp) <= REL_TOL*std::abs(tmp));
    }
  }

  cout << "\tTest Data::asin." << endl;
  resultEx.copy(baseEx.asin());
  resultCon.copy(baseCon.asin());
  resultTag.copy(baseTag.asin());
  assert(true);

  cout << "\tTest Data::acos." << endl;
  resultEx.copy(baseEx.acos());
  resultCon.copy(baseCon.acos());
  resultTag.copy(baseTag.acos());
  assert(true);

  cout << "\tTest Data::atan." << endl;
  resultEx.copy(baseEx.atan());
  resultCon.copy(baseCon.atan());
  resultTag.copy(baseTag.atan());
  for (int i=0;i<shape[0];i++) {
    for (int j=0;j<shape[1];j++) {
      tmp=atan((double)dataView.index(i,j));
      assert(std::abs(resultEx.getPointDataView()(i,j) - tmp) <= REL_TOL*std::abs(tmp));
      assert(std::abs(resultCon.getPointDataView()(i,j)- tmp) <= REL_TOL*std::abs(tmp));
      assert(std::abs(resultTag.getPointDataView()(i,j)- tmp) <= REL_TOL*std::abs(tmp));
    }
  }

  cout << "\tTest Data::sinh." << endl;
  resultEx.copy(baseEx.sinh());
  resultCon.copy(baseCon.sinh());
  resultTag.copy(baseTag.sinh());
  for (int i=0;i<shape[0];i++) {
    for (int j=0;j<shape[1];j++) {
      tmp=sinh((double)dataView.index(i,j));
      assert(std::abs(resultEx.getPointDataView()(i,j) - tmp) <= REL_TOL*std::abs(tmp));
      assert(std::abs(resultCon.getPointDataView()(i,j)- tmp) <= REL_TOL*std::abs(tmp));
      assert(std::abs(resultTag.getPointDataView()(i,j)- tmp) <= REL_TOL*std::abs(tmp));
    }
  }

  cout << "\tTest Data::cosh." << endl;
  resultEx.copy(baseEx.cosh());
  resultCon.copy(baseCon.cosh());
  resultTag.copy(baseTag.cosh());
  for (int i=0;i<shape[0];i++) {
    for (int j=0;j<shape[1];j++) {
      tmp=cosh((double)dataView.index(i,j));
      assert(std::abs(resultEx.getPointDataView()(i,j) - tmp) <= REL_TOL*std::abs(tmp));
      assert(std::abs(resultCon.getPointDataView()(i,j)- tmp) <= REL_TOL*std::abs(tmp));
      assert(std::abs(resultTag.getPointDataView()(i,j)- tmp) <= REL_TOL*std::abs(tmp));
    }
  }

  cout << "\tTest Data::tanh." << endl;
  resultEx.copy(baseEx.tanh());
  resultCon.copy(baseCon.tanh());
  resultTag.copy(baseTag.tanh());
  for (int i=0;i<shape[0];i++) {
    for (int j=0;j<shape[1];j++) {
      tmp=tanh((double)dataView.index(i,j));
      assert(std::abs(resultEx.getPointDataView()(i,j) - tmp) <= REL_TOL*std::abs(tmp));
      assert(std::abs(resultCon.getPointDataView()(i,j)- tmp) <= REL_TOL*std::abs(tmp));
      assert(std::abs(resultTag.getPointDataView()(i,j)- tmp) <= REL_TOL*std::abs(tmp));
    }
  }

  cout << "\tTest Data::asinh." << endl;
  resultEx.copy(baseEx.asinh());
  resultCon.copy(baseCon.asinh());
  resultTag.copy(baseTag.asinh());
  assert(true);

  cout << "\tTest Data::acosh." << endl;
  resultEx.copy(baseEx.acosh());
  resultCon.copy(baseCon.acosh());
  resultTag.copy(baseTag.acosh());
  assert(true);

  cout << "\tTest Data::atanh." << endl;
  resultEx.copy(baseEx.atanh());
  resultCon.copy(baseCon.atanh());
  resultTag.copy(baseTag.atanh());
  assert(true);

  cout << "\tTest Data::log." << endl;
  resultEx.copy(baseEx.log());
  resultCon.copy(baseCon.log());
  resultTag.copy(baseTag.log());
  assert(true);

  cout << "\tTest Data::abs." << endl;
  resultEx.copy(baseEx.abs());
  resultCon.copy(baseCon.abs());
  resultTag.copy(baseTag.abs());
  for (int i=0;i<shape[0];i++) {
    for (int j=0;j<shape[1];j++) {
      tmp=abs((double)dataView.index(i,j));
      assert(std::abs(resultEx.getPointDataView()(i,j) - tmp) <= REL_TOL*std::abs(tmp));
      assert(std::abs(resultCon.getPointDataView()(i,j)- tmp) <= REL_TOL*std::abs(tmp));
      assert(std::abs(resultTag.getPointDataView()(i,j)- tmp) <= REL_TOL*std::abs(tmp));
    }
  }

  cout << "\tTest Data::sign." << endl;
  resultEx.copy(baseEx.sign());
  resultCon.copy(baseCon.sign());
  resultTag.copy(baseTag.sign());
  assert(true);

  cout << "\tTest Data::exp." << endl;
  resultEx.copy(baseEx.exp());
  resultCon.copy(baseCon.exp());
  resultTag.copy(baseTag.exp());
  for (int i=0;i<shape[0];i++) {
    for (int j=0;j<shape[1];j++) {
      tmp=exp((double)dataView.index(i,j));
      assert(std::abs(resultEx.getPointDataView()(i,j) - tmp) <= REL_TOL*std::abs(tmp));
      assert(std::abs(resultCon.getPointDataView()(i,j)- tmp) <= REL_TOL*std::abs(tmp));
      assert(std::abs(resultTag.getPointDataView()(i,j)- tmp) <= REL_TOL*std::abs(tmp));
    }
  }

  cout << "\tTest Data::sqrt." << endl;
  resultEx.copy(baseEx.sqrt());
  resultCon.copy(baseCon.sqrt());
  resultTag.copy(baseTag.sqrt());
  for (int i=0;i<shape[0];i++) {
    for (int j=0;j<shape[1];j++) {
      tmp=sqrt((double)dataView.index(i,j));
      assert(std::abs(resultEx.getPointDataView()(i,j) - tmp) <= REL_TOL*std::abs(tmp));
      assert(std::abs(resultCon.getPointDataView()(i,j)- tmp) <= REL_TOL*std::abs(tmp));
      assert(std::abs(resultTag.getPointDataView()(i,j)- tmp) <= REL_TOL*std::abs(tmp));
    }
  }

  cout << "\tTest Data::neg." << endl;
  resultEx.copy(baseEx.neg());
  resultCon.copy(baseCon.neg());
  resultTag.copy(baseTag.neg());
  assert(true);

  cout << "\tTest Data::pos." << endl;
  resultEx.copy(baseEx.pos());
  resultCon.copy(baseCon.pos());
  resultTag.copy(baseTag.pos());
  for (int i=0;i<shape[0];i++) {
    for (int j=0;j<shape[1];j++) {
      assert(std::abs(resultEx.getPointDataView()(i,j) - dataView.index(i,j)) <= REL_TOL*std::abs(dataView.index(i,j)));
      assert(std::abs(resultCon.getPointDataView()(i,j) - dataView.index(i,j)) <= REL_TOL*std::abs(dataView.index(i,j)));
      assert(std::abs(resultTag.getPointDataView()(i,j) - dataView.index(i,j)) <= REL_TOL*std::abs(dataView.index(i,j)));
    }
  }

  // test reduction operations

  cout << "\tTest Data::Lsup." << endl;
  assert(std::abs(baseEx.Lsup() - 5) <= REL_TOL*5);
  assert(std::abs(baseCon.Lsup() - 5) <= REL_TOL*5);
  assert(std::abs(baseTag.Lsup() - 5) <= REL_TOL*5);

  cout << "\tTest Data::sup." << endl;
  assert(std::abs(baseEx.sup() - 5) <= REL_TOL*5);
  assert(std::abs(baseCon.sup() - 5) <= REL_TOL*5);
  assert(std::abs(baseTag.sup() - 5) <= REL_TOL*5);

  cout << "\tTest Data::inf." << endl;
  assert(std::abs(baseEx.inf() - 0) <= REL_TOL*0);
  assert(std::abs(baseCon.inf() - 0) <= REL_TOL*0);
  assert(std::abs(baseTag.inf() - 0) <= REL_TOL*0);

  // test data-point reduction operations

  cout << "\tTest Data::minval." << endl;
  resultEx.copy(baseEx.minval());
  resultCon.copy(baseCon.minval());
  resultTag.copy(baseTag.minval());
  assert(std::abs(resultEx.getPointDataView()() - 0) <= REL_TOL*0);
  assert(std::abs(resultCon.getPointDataView()() - 0) <= REL_TOL*0);
  assert(std::abs(resultTag.getPointDataView()() - 0) <= REL_TOL*0);

  cout << "\tTest Data::maxval." << endl;
  resultEx.copy(baseEx.maxval());
  resultCon.copy(baseCon.maxval());
  resultTag.copy(baseTag.maxval());
  assert(std::abs(resultEx.getPointDataView()() - 5) <= REL_TOL*5);
  assert(std::abs(resultCon.getPointDataView()() - 5) <= REL_TOL*5);
  assert(std::abs(resultTag.getPointDataView()() - 5) <= REL_TOL*5);

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
