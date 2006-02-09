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

#include "EsysException.h"

#include "DataTagged.h"
#include "DataConstant.h"

#include "DataTaggedTestCase.h"

#include "BinaryOp.h"
#include "UnaryOp.h"
#include "FunctionSpaceFactory.h"
#include "DataFactory.h"

#include <iostream>
#include <functional>
#include <algorithm>

using namespace CppUnitTest;
using namespace escript;
using namespace esysUtils;
using namespace std;

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
    cout << "\tTest binaryOp addition of two default DataTagged objects." << endl;

    DataTagged myData;
    DataTagged right;

    binaryOp(myData,right,plus<double>());

    //cout << myData.toString() << endl;

    assert(myData.getNumSamples()==1);
    assert(myData.getNumDPPSample()==1);

    assert(myData.validSamplePointNo(0));
    assert(myData.validSampleNo(0));
    assert(!myData.validSamplePointNo(1));
    assert(!myData.validSampleNo(1));

    // data-point 0 has tag number 1 by default
    assert(myData.getTagNumber(0)==1);

    assert(!myData.isCurrentTag(1));

    assert(myData.getTagLookup().size()==0);

    assert(myData.getLength()==1);

    assert(myData.getPointOffset(0,0)==0);

    DataArrayView myDataView = myData.getDataPoint(0,0);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==0);
    assert(myDataView.getRank()==0);
    assert(myDataView.noValues()==1);
    assert(myDataView.getShape().size()==0);
    assert(myDataView()==0.0);

    // Test non-existent tag returns the default value.
    myDataView = myData.getDataPointByTag(1);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==0);
    assert(myDataView.getRank()==0);
    assert(myDataView.noValues()==1);
    assert(myDataView.getShape().size()==0);
    assert(myDataView()==0.0);

    myDataView = myData.getDefaultValue();
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==0);
    assert(myDataView.getRank()==0);
    assert(myDataView.noValues()==1);
    assert(myDataView.getShape().size()==0);
    assert(myDataView()==0.0);

    // use a non-existent tag so we get a pointer to
    // the first element of the data array
    double* sampleData=myData.getSampleDataByTag(9);
    for (int i=0; i<myData.getLength(); i++) {
      assert(sampleData[i]==i);
    }

  }

  {
    cout << "\tTest binaryOp addition of two DataTagged objects with default values only." << endl;

    DataArrayView::ShapeType viewShape;
    viewShape.push_back(3);

    DataTagged::TagListType keys;

    DataTagged::ValueListType values;

    DataArrayView::ValueType viewData(3);
    for (int i=0;i<viewShape[0];i++) {
      viewData[i]=i;
    }
    DataArrayView myView(viewData,viewShape);

    DataTagged myData(keys,values,myView,FunctionSpace());
    DataTagged right(keys,values,myView,FunctionSpace());

    binaryOp(myData,right,plus<double>());

    //cout << myData.toString() << endl;

    assert(myData.getNumSamples()==1);
    assert(myData.getNumDPPSample()==1);

    assert(myData.validSamplePointNo(0));
    assert(myData.validSampleNo(0));
    assert(!myData.validSamplePointNo(1));
    assert(!myData.validSampleNo(1));

    // data-point 0 has tag number 1 by default
    assert(myData.getTagNumber(0)==1);

    assert(!myData.isCurrentTag(1));

    assert(myData.getTagLookup().size()==0);

    assert(myData.getLength()==3);

    assert(myData.getPointOffset(0,0)==0);

    DataArrayView myDataView = myData.getDefaultValue();
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==0);
    assert(myDataView.getRank()==1);
    assert(myDataView.noValues()==3);
    assert(myDataView.getShape().size()==1);
    assert(myDataView(0)==0);
    assert(myDataView(1)==2);
    assert(myDataView(2)==4);

    // use a non-existent tag so we get a pointer to
    // the first element of the data array
    double* sampleData=myData.getSampleDataByTag(9);
    for (int i=0; i<myData.getLength(); i++) {
      assert(sampleData[i]==i*2);
    }

  }

  {
    cout << "\tTest binaryOp addition of two DataTagged objects with one identical tag each." << endl;

    DataTagged myData;
    DataTagged right;

    DataArray vOne(1.0);
    myData.addTaggedValue(1,vOne.getView());
    right.addTaggedValue(1,vOne.getView());

    binaryOp(myData,right,plus<double>());

    assert(myData.getNumSamples()==1);
    assert(myData.getNumDPPSample()==1);

    assert(myData.validSamplePointNo(0));
    assert(myData.validSampleNo(0));
    assert(!myData.validSamplePointNo(1));
    assert(!myData.validSampleNo(1));

    // data-point 0 has tag number 1 by default
    assert(myData.getTagNumber(0)==1);

    assert(myData.isCurrentTag(1));

    assert(myData.getTagLookup().size()==1);

    assert(myData.getLength()==2);

    assert(myData.getPointOffset(0,0)==1);

    // check result value for tag "1"
    DataArrayView myDataView = myData.getDataPointByTag(1);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==1);
    assert(myDataView.getRank()==0);
    assert(myDataView.noValues()==1);
    assert(myDataView.getShape().size()==0);
    assert(myDataView()==2.0);

    // check result for default value
    myDataView = myData.getDefaultValue();
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==0);
    assert(myDataView.getRank()==0);
    assert(myDataView.noValues()==1);
    assert(myDataView.getShape().size()==0);
    assert(myDataView()==0.0);

    // use a non-existent tag so we get a pointer to
    // the first element of the data array
    double* sampleData=myData.getSampleDataByTag(9);
    for (int i=0; i<myData.getLength(); i++) {
      assert(sampleData[i]==i*2);
    }

  }

  {
    cout << "\tTest binaryOp addition of two DataTagged objects with one different tag each." << endl;

    DataTagged myData;
    DataTagged right;

    DataArray vOne(1.0);
    DataArray vTwo(2.0);
    myData.addTaggedValue(1,vOne.getView());
    right.addTaggedValue(2,vTwo.getView());

    binaryOp(myData,right,plus<double>());

    assert(myData.getNumSamples()==1);
    assert(myData.getNumDPPSample()==1);

    assert(myData.validSamplePointNo(0));
    assert(myData.validSampleNo(0));
    assert(!myData.validSamplePointNo(1));
    assert(!myData.validSampleNo(1));

    // data-point 0 has tag number 1 by default
    assert(myData.getTagNumber(0)==1);

    assert(myData.isCurrentTag(1));
    assert(myData.isCurrentTag(2));

    assert(myData.getTagLookup().size()==2);

    assert(myData.getLength()==3);

    assert(myData.getPointOffset(0,0)==1);

    // check result value for tag "1"
    DataArrayView myDataView = myData.getDataPointByTag(1);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==1);
    assert(myDataView.getRank()==0);
    assert(myDataView.noValues()==1);
    assert(myDataView.getShape().size()==0);
    assert(myDataView()==1.0);

    // check result value for tag "2"
    myDataView = myData.getDataPointByTag(2);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==2);
    assert(myDataView.getRank()==0);
    assert(myDataView.noValues()==1);
    assert(myDataView.getShape().size()==0);
    assert(myDataView()==2.0);

    // check result for default value
    myDataView = myData.getDefaultValue();
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==0);
    assert(myDataView.getRank()==0);
    assert(myDataView.noValues()==1);
    assert(myDataView.getShape().size()==0);
    assert(myDataView()==0.0);

    // use a non-existent tag so we get a pointer to
    // the first element of the data array
    double* sampleData=myData.getSampleDataByTag(9);
    for (int i=0; i<myData.getLength(); i++) {
      assert(sampleData[i]==i);
    }

  }

  {
    cout << "\tTest binaryOp addition of two DataTagged objects with overlapping tag sets." << endl;

    DataTagged myData;
    DataTagged right;

    DataArray vOne(1.0);
    myData.addTaggedValue(1,vOne.getView());
    myData.addTaggedValue(2,vOne.getView());
    right.addTaggedValue(2,vOne.getView());
    right.addTaggedValue(3,vOne.getView());

    binaryOp(myData,right,plus<double>());

    assert(myData.getNumSamples()==1);
    assert(myData.getNumDPPSample()==1);

    assert(myData.validSamplePointNo(0));
    assert(myData.validSampleNo(0));
    assert(!myData.validSamplePointNo(1));
    assert(!myData.validSampleNo(1));

    // data-point 0 has tag number 1 by default
    assert(myData.getTagNumber(0)==1);

    assert(myData.isCurrentTag(1));
    assert(myData.isCurrentTag(2));
    assert(myData.isCurrentTag(3));

    assert(myData.getTagLookup().size()==3);

    assert(myData.getLength()==4);

    assert(myData.getPointOffset(0,0)==1);

    // check result value for tag "1"
    DataArrayView myDataView = myData.getDataPointByTag(1);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==1);
    assert(myDataView.getRank()==0);
    assert(myDataView.noValues()==1);
    assert(myDataView.getShape().size()==0);
    assert(myDataView()==1.0);

    // check result value for tag "2"
    myDataView = myData.getDataPointByTag(2);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==2);
    assert(myDataView.getRank()==0);
    assert(myDataView.noValues()==1);
    assert(myDataView.getShape().size()==0);
    assert(myDataView()==2.0);

    // check result value for tag "3"
    myDataView = myData.getDataPointByTag(3);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==3);
    assert(myDataView.getRank()==0);
    assert(myDataView.noValues()==1);
    assert(myDataView.getShape().size()==0);
    assert(myDataView()==1.0);

    // check result for default value
    myDataView = myData.getDefaultValue();
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==0);
    assert(myDataView.getRank()==0);
    assert(myDataView.noValues()==1);
    assert(myDataView.getShape().size()==0);
    assert(myDataView()==0.0);

    // use a non-existent tag so we get a pointer to
    // the first element of the data array
    double* sampleData=myData.getSampleDataByTag(9);
    for (int i=0; i<myData.getLength(); i++) {
      if (i<3) {
        assert(sampleData[i]==i);
      }
      if (i>=3) {
        assert(sampleData[i]==i-2);
      }
    }

  }

  {
    cout << "\tTest binaryOp multiplication of two DataTagged objects with default values only." << endl;

    DataArrayView::ShapeType viewShape;
    viewShape.push_back(3);

    DataTagged::TagListType keys;

    DataTagged::ValueListType values;

    DataArrayView::ValueType viewData(3);
    for (int i=0;i<viewShape[0];i++) {
      viewData[i]=i;
    }
    DataArrayView myView(viewData,viewShape);

    DataTagged myData(keys,values,myView,FunctionSpace());
    DataTagged right(keys,values,myView,FunctionSpace());

    binaryOp(myData,right,multiplies<double>());

    //cout << myData.toString() << endl;

    assert(myData.getNumSamples()==1);
    assert(myData.getNumDPPSample()==1);

    assert(myData.validSamplePointNo(0));
    assert(myData.validSampleNo(0));
    assert(!myData.validSamplePointNo(1));
    assert(!myData.validSampleNo(1));

    // data-point 0 has tag number 1 by default
    assert(myData.getTagNumber(0)==1);

    assert(!myData.isCurrentTag(1));

    assert(myData.getTagLookup().size()==0);

    assert(myData.getLength()==3);

    assert(myData.getPointOffset(0,0)==0);

    DataArrayView myDataView = myData.getDefaultValue();
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==0);
    assert(myDataView.getRank()==1);
    assert(myDataView.noValues()==3);
    assert(myDataView.getShape().size()==1);
    assert(myDataView(0)==0);
    assert(myDataView(1)==1);
    assert(myDataView(2)==4);

    // use a non-existent tag so we get a pointer to
    // the first element of the data array
    double* sampleData=myData.getSampleDataByTag(9);
    for (int i=0; i<myData.getLength(); i++) {
      assert(sampleData[i]==i*i);
    }

  }

  {
    cout << "\tTest binaryOp multiplication of two DataTagged objects with overlapping tag sets." << endl;

    DataTagged myData;
    DataTagged right;

    DataArray vOne(1.0);
    DataArray vTwo(2.0);
    myData.addTaggedValue(1,vOne.getView());
    myData.addTaggedValue(2,vOne.getView());
    right.addTaggedValue(2,vTwo.getView());
    right.addTaggedValue(3,vTwo.getView());

    binaryOp(myData,right,multiplies<double>());

    assert(myData.getNumSamples()==1);
    assert(myData.getNumDPPSample()==1);

    assert(myData.validSamplePointNo(0));
    assert(myData.validSampleNo(0));
    assert(!myData.validSamplePointNo(1));
    assert(!myData.validSampleNo(1));

    // data-point 0 has tag number 1 by default
    assert(myData.getTagNumber(0)==1);

    assert(myData.isCurrentTag(1));
    assert(myData.isCurrentTag(2));
    assert(myData.isCurrentTag(3));

    assert(myData.getTagLookup().size()==3);

    assert(myData.getLength()==4);

    assert(myData.getPointOffset(0,0)==1);

    // check result value for tag "1"
    DataArrayView myDataView = myData.getDataPointByTag(1);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==1);
    assert(myDataView.getRank()==0);
    assert(myDataView.noValues()==1);
    assert(myDataView.getShape().size()==0);
    assert(myDataView()==0.0);

    // check result value for tag "2"
    myDataView = myData.getDataPointByTag(2);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==2);
    assert(myDataView.getRank()==0);
    assert(myDataView.noValues()==1);
    assert(myDataView.getShape().size()==0);
    assert(myDataView()==2.0);

    // check result value for tag "3"
    myDataView = myData.getDataPointByTag(3);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==3);
    assert(myDataView.getRank()==0);
    assert(myDataView.noValues()==1);
    assert(myDataView.getShape().size()==0);
    assert(myDataView()==0.0);

    // check result for default value
    myDataView = myData.getDefaultValue();
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==0);
    assert(myDataView.getRank()==0);
    assert(myDataView.noValues()==1);
    assert(myDataView.getShape().size()==0);
    assert(myDataView()==0.0);

    // use a non-existent tag so we get a pointer to
    // the first element of the data array
    double* sampleData=myData.getSampleDataByTag(9);
    for (int i=0; i<myData.getLength(); i++) {
      if (i==2) {
        assert(sampleData[i]==2);
      } else {
        assert(sampleData[i]==0);
      }
    }

  }

  {
    cout << "\tTest unaryOp negate on default DataTagged object." << endl;

    DataTagged myData;

    unaryOp(myData,negate<double>());

    //cout << myData.toString() << endl;

    assert(myData.getNumSamples()==1);
    assert(myData.getNumDPPSample()==1);

    assert(myData.validSamplePointNo(0));
    assert(myData.validSampleNo(0));
    assert(!myData.validSamplePointNo(1));
    assert(!myData.validSampleNo(1));

    // data-point 0 has tag number 1 by default
    assert(myData.getTagNumber(0)==1);

    assert(!myData.isCurrentTag(1));

    assert(myData.getTagLookup().size()==0);

    assert(myData.getLength()==1);

    assert(myData.getPointOffset(0,0)==0);

    DataArrayView myDataView = myData.getDataPoint(0,0);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==0);
    assert(myDataView.getRank()==0);
    assert(myDataView.noValues()==1);
    assert(myDataView.getShape().size()==0);
    assert(myDataView()==0.0);

    // Test non-existent tag returns the default value.
    myDataView = myData.getDataPointByTag(1);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==0);
    assert(myDataView.getRank()==0);
    assert(myDataView.noValues()==1);
    assert(myDataView.getShape().size()==0);
    assert(myDataView()==0.0);

    myDataView = myData.getDefaultValue();
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==0);
    assert(myDataView.getRank()==0);
    assert(myDataView.noValues()==1);
    assert(myDataView.getShape().size()==0);
    assert(myDataView()==0.0);

    // use a non-existent tag so we get a pointer to
    // the first element of the data array
    double* sampleData=myData.getSampleDataByTag(9);
    for (int i=0; i<myData.getLength(); i++) {
      assert(sampleData[i]==i);
    }

  }

  {
    cout << "\tTest unaryOp negate on DataTagged object with default value only." << endl;

    DataArrayView::ShapeType viewShape;
    viewShape.push_back(3);

    DataTagged::TagListType keys;

    DataTagged::ValueListType values;

    DataArrayView::ValueType viewData(3);
    for (int i=0;i<viewShape[0];i++) {
      viewData[i]=i;
    }
    DataArrayView myView(viewData,viewShape);

    DataTagged myData(keys,values,myView,FunctionSpace());

    unaryOp(myData,negate<double>());

    //cout << myData.toString() << endl;

    assert(myData.getNumSamples()==1);
    assert(myData.getNumDPPSample()==1);

    assert(myData.validSamplePointNo(0));
    assert(myData.validSampleNo(0));
    assert(!myData.validSamplePointNo(1));
    assert(!myData.validSampleNo(1));

    // data-point 0 has tag number 1 by default
    assert(myData.getTagNumber(0)==1);

    assert(!myData.isCurrentTag(1));

    assert(myData.getTagLookup().size()==0);

    assert(myData.getLength()==3);

    assert(myData.getPointOffset(0,0)==0);

    DataArrayView myDataView = myData.getDefaultValue();
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==0);
    assert(myDataView.getRank()==1);
    assert(myDataView.noValues()==3);
    assert(myDataView.getShape().size()==1);
    assert(myDataView(0)==0);
    assert(myDataView(1)==-1);
    assert(myDataView(2)==-2);

    // use a non-existent tag so we get a pointer to
    // the first element of the data array
    double* sampleData=myData.getSampleDataByTag(9);
    for (int i=0; i<myData.getLength(); i++) {
      assert(sampleData[i]==0-i);
    }

  }

  {
    cout << "\tTest unaryOp negate on DataTagged object with two tags." << endl;

    DataTagged myData;

    DataArray vOne(1.0);
    DataArray vTwo(2.0);
    myData.addTaggedValue(1,vOne.getView());
    myData.addTaggedValue(2,vTwo.getView());

    unaryOp(myData,negate<double>());

    assert(myData.getNumSamples()==1);
    assert(myData.getNumDPPSample()==1);

    assert(myData.validSamplePointNo(0));
    assert(myData.validSampleNo(0));
    assert(!myData.validSamplePointNo(1));
    assert(!myData.validSampleNo(1));

    // data-point 0 has tag number 1 by default
    assert(myData.getTagNumber(0)==1);

    assert(myData.isCurrentTag(1));
    assert(myData.isCurrentTag(2));

    assert(myData.getTagLookup().size()==2);

    assert(myData.getLength()==3);

    assert(myData.getPointOffset(0,0)==1);

    // check result value for tag "1"
    DataArrayView myDataView = myData.getDataPointByTag(1);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==1);
    assert(myDataView.getRank()==0);
    assert(myDataView.noValues()==1);
    assert(myDataView.getShape().size()==0);
    assert(myDataView()==-1.0);

    // check result value for tag "2"
    myDataView = myData.getDataPointByTag(2);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==2);
    assert(myDataView.getRank()==0);
    assert(myDataView.noValues()==1);
    assert(myDataView.getShape().size()==0);
    assert(myDataView()==-2.0);

    // check result for default value
    myDataView = myData.getDefaultValue();
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==0);
    assert(myDataView.getRank()==0);
    assert(myDataView.noValues()==1);
    assert(myDataView.getShape().size()==0);
    assert(myDataView()==0.0);

    // use a non-existent tag so we get a pointer to
    // the first element of the data array
    double* sampleData=myData.getSampleDataByTag(9);
    for (int i=0; i<myData.getLength(); i++) {
      assert(sampleData[i]==0-i);
    }

  }

}

void DataTaggedTestCase::testAddTaggedValues() {

  cout << endl;

  {

    cout << "\tTest adding one key with empty value list to default DataTagged." << endl;
    DataTagged myData;

    DataTagged::TagListType keys;
    keys.push_back(1);

    DataTagged::ValueListType values;

    myData.addTaggedValues(keys,values);

    assert(myData.isCurrentTag(1));

    assert(myData.getTagLookup().size()==1);

    assert(myData.getLength()==2);

    // data-point 0 has tag number 1 by default
    assert(myData.getTagNumber(0)==1);

    assert(myData.getPointOffset(0,0)==1);

    DataArrayView myDataView = myData.getDataPoint(0,0);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==1);
    assert(myDataView.getRank()==0);
    assert(myDataView.noValues()==1);
    assert(myDataView.getShape().size()==0);
    assert(myDataView()==0.0);

    myDataView = myData.getDataPointByTag(1);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==1);
    assert(myDataView.getRank()==0);
    assert(myDataView.noValues()==1);
    assert(myDataView.getShape().size()==0);
    assert(myDataView()==0.0);

    myDataView = myData.getDefaultValue();
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==0);
    assert(myDataView.getRank()==0);
    assert(myDataView.noValues()==1);
    assert(myDataView.getShape().size()==0);
    assert(myDataView()==0.0);

    // use a non-existent tag so we get a pointer to
    // the first element of the data array
    double* sampleData=myData.getSampleDataByTag(9);
    for (int i=0; i<myData.getLength(); i++) {
      assert(sampleData[i]==0);
    }

  }

  {

    cout << "\tTest adding one key with one value to default DataTagged." << endl;
    DataTagged myData;

    DataTagged::TagListType keys;
    keys.push_back(1);

    DataTagged::ValueListType values;

    DataArrayView::ShapeType viewShape;
    DataArrayView::ValueType viewData(1);
    viewData[0]=1.0;
    DataArrayView myView(viewData,viewShape);
    values.push_back(myView);

    myData.addTaggedValues(keys,values);

    assert(myData.isCurrentTag(1));

    assert(myData.getTagLookup().size()==1);

    assert(myData.getLength()==2);

    // data-point 0 has tag number 1 by default
    assert(myData.getTagNumber(0)==1);

    assert(myData.getPointOffset(0,0)==1);

    DataArrayView myDataView = myData.getDataPoint(0,0);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==1);
    assert(myDataView.getRank()==0);
    assert(myDataView.noValues()==1);
    assert(myDataView.getShape().size()==0);
    assert(myDataView()==1.0);

    myDataView = myData.getDataPointByTag(1);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==1);
    assert(myDataView.getRank()==0);
    assert(myDataView.noValues()==1);
    assert(myDataView.getShape().size()==0);
    assert(myDataView()==1.0);

    myDataView = myData.getDefaultValue();
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==0);
    assert(myDataView.getRank()==0);
    assert(myDataView.noValues()==1);
    assert(myDataView.getShape().size()==0);
    assert(myDataView()==0.0);

    // use a non-existent tag so we get a pointer to
    // the first element of the data array
    double* sampleData=myData.getSampleDataByTag(9);
    for (int i=0; i<myData.getLength(); i++) {
      assert(sampleData[i]==i);
    }

  }

  {

    cout << "\tTest adding three keys with one value to default DataTagged." << endl;
    DataTagged myData;

    DataTagged::TagListType keys;
    keys.push_back(1);
    keys.push_back(2);
    keys.push_back(3);

    DataTagged::ValueListType values;

    DataArrayView::ShapeType viewShape;
    DataArrayView::ValueType viewData(1);
    viewData[0]=1.0;
    DataArrayView myView(viewData,viewShape);
    values.push_back(myView);

    myData.addTaggedValues(keys,values);

    assert(myData.isCurrentTag(1));
    assert(myData.isCurrentTag(2));
    assert(myData.isCurrentTag(3));

    assert(myData.getTagLookup().size()==3);

    assert(myData.getLength()==4);

    // data-point 0 has tag number 1 by default
    assert(myData.getTagNumber(0)==1);

    assert(myData.getPointOffset(0,0)==1);

    DataArrayView myDataView = myData.getDataPoint(0,0);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==1);
    assert(myDataView.getRank()==0);
    assert(myDataView.noValues()==1);
    assert(myDataView.getShape().size()==0);
    assert(myDataView()==1.0);

    myDataView = myData.getDataPointByTag(1);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==1);
    assert(myDataView.getRank()==0);
    assert(myDataView.noValues()==1);
    assert(myDataView.getShape().size()==0);
    assert(myDataView()==1.0);

    myDataView = myData.getDataPointByTag(2);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==2);
    assert(myDataView.getRank()==0);
    assert(myDataView.noValues()==1);
    assert(myDataView.getShape().size()==0);
    assert(myDataView()==1.0);

    myDataView = myData.getDataPointByTag(3);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==3);
    assert(myDataView.getRank()==0);
    assert(myDataView.noValues()==1);
    assert(myDataView.getShape().size()==0);
    assert(myDataView()==1.0);

    myDataView = myData.getDefaultValue();
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==0);
    assert(myDataView.getRank()==0);
    assert(myDataView.noValues()==1);
    assert(myDataView.getShape().size()==0);
    assert(myDataView()==0.0);

    // use a non-existent tag so we get a pointer to
    // the first element of the data array
    double* sampleData=myData.getSampleDataByTag(9);
    for (int i=0; i<myData.getLength(); i++) {
      if (i==0) {
        assert(sampleData[i]==0);
      } else {
        assert(sampleData[i]==1);
      }
    }

  }

  {

    cout << "\tTest adding three keys with three values to default DataTagged." << endl;
    DataTagged myData;

    DataTagged::TagListType keys;
    keys.push_back(1);
    keys.push_back(2);
    keys.push_back(3);

    DataTagged::ValueListType values;

    DataArrayView::ShapeType viewShape;
    DataArrayView::ValueType viewData1(1);
    viewData1[0]=1.0;
    DataArrayView::ValueType viewData2(1);
    viewData2[0]=2.0;
    DataArrayView::ValueType viewData3(1);
    viewData3[0]=3.0;
    DataArrayView myView1(viewData1,viewShape);
    DataArrayView myView2(viewData2,viewShape);
    DataArrayView myView3(viewData3,viewShape);
    values.push_back(myView1);
    values.push_back(myView2);
    values.push_back(myView3);

    myData.addTaggedValues(keys,values);

    assert(myData.isCurrentTag(1));
    assert(myData.isCurrentTag(2));
    assert(myData.isCurrentTag(3));

    assert(myData.getTagLookup().size()==3);

    assert(myData.getLength()==4);

    // data-point 0 has tag number 1 by default
    assert(myData.getTagNumber(0)==1);

    assert(myData.getPointOffset(0,0)==1);

    DataArrayView myDataView = myData.getDataPoint(0,0);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==1);
    assert(myDataView.getRank()==0);
    assert(myDataView.noValues()==1);
    assert(myDataView.getShape().size()==0);
    assert(myDataView()==1.0);

    myDataView = myData.getDataPointByTag(1);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==1);
    assert(myDataView.getRank()==0);
    assert(myDataView.noValues()==1);
    assert(myDataView.getShape().size()==0);
    assert(myDataView()==1.0);

    myDataView = myData.getDataPointByTag(2);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==2);
    assert(myDataView.getRank()==0);
    assert(myDataView.noValues()==1);
    assert(myDataView.getShape().size()==0);
    assert(myDataView()==2.0);

    myDataView = myData.getDataPointByTag(3);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==3);
    assert(myDataView.getRank()==0);
    assert(myDataView.noValues()==1);
    assert(myDataView.getShape().size()==0);
    assert(myDataView()==3.0);

    myDataView = myData.getDefaultValue();
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==0);
    assert(myDataView.getRank()==0);
    assert(myDataView.noValues()==1);
    assert(myDataView.getShape().size()==0);
    assert(myDataView()==0.0);

    // use a non-existent tag so we get a pointer to
    // the first element of the data array
    double* sampleData=myData.getSampleDataByTag(9);
    for (int i=0; i<myData.getLength(); i++) {
      assert(sampleData[i]==i);
    }

  }

  {

    cout << "\tTest adding one key with empty value list to DataTagged with default value only." << endl;

    DataArrayView::ShapeType viewShape;
    viewShape.push_back(3);

    DataTagged::TagListType keys;

    DataTagged::ValueListType values;

    DataArrayView::ValueType viewData(3);
    for (int i=0;i<viewShape[0];i++) {
      viewData[i]=i;
    }
    DataArrayView myView(viewData,viewShape);

    DataTagged myData(keys,values,myView,FunctionSpace());

    keys.push_back(1);
    values.clear();

    myData.addTaggedValues(keys,values);

    assert(myData.isCurrentTag(1));

    assert(myData.getTagLookup().size()==1);

    assert(myData.getLength()==6);

    // data-point 0 has tag number 1 by default
    assert(myData.getTagNumber(0)==1);

    assert(myData.getPointOffset(0,0)==3);

    DataArrayView myDataView = myData.getDataPoint(0,0);
    assert(myDataView==myView);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==3);
    assert(myDataView.getRank()==1);
    assert(myDataView.noValues()==3);
    assert(myDataView.getShape().size()==1);
    assert(myDataView(0)==0);
    assert(myDataView(1)==1);
    assert(myDataView(2)==2);

    myDataView = myData.getDataPointByTag(1);
    assert(myDataView==myView);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==3);
    assert(myDataView.getRank()==1);
    assert(myDataView.noValues()==3);
    assert(myDataView.getShape().size()==1);
    assert(myDataView(0)==0);
    assert(myDataView(1)==1);
    assert(myDataView(2)==2);

    myDataView = myData.getDefaultValue();
    assert(myDataView==myView);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==0);
    assert(myDataView.getRank()==1);
    assert(myDataView.noValues()==3);
    assert(myDataView.getShape().size()==1);
    assert(myDataView(0)==0);
    assert(myDataView(1)==1);
    assert(myDataView(2)==2);

    // use a non-existent tag so we get a pointer to
    // the first element of the data array
    double* sampleData=myData.getSampleDataByTag(9);
    for (int i=0; i<myData.getLength(); i++) {
      assert(sampleData[i]==i%3);
    }

  }

  {

    cout << "\tTest adding one key with one value to DataTagged with default value only." << endl;

    DataArrayView::ShapeType viewShape;
    viewShape.push_back(3);

    DataTagged::TagListType keys;

    DataTagged::ValueListType values;

    DataArrayView::ValueType viewData(3);
    for (int i=0;i<viewShape[0];i++) {
      viewData[i]=i;
    }
    DataArrayView myView(viewData,viewShape);

    DataTagged myData(keys,values,myView,FunctionSpace());

    keys.push_back(1);

    DataArrayView::ValueType viewData1(3);
    for (int i=0;i<viewShape[0];i++) {
      viewData1[i]=i+3;
    }
    DataArrayView myView1(viewData1,viewShape);
    values.push_back(myView1);

    myData.addTaggedValues(keys,values);

    assert(myData.isCurrentTag(1));

    assert(myData.getTagLookup().size()==1);

    assert(myData.getLength()==6);

    // data-point 0 has tag number 1 by default
    assert(myData.getTagNumber(0)==1);

    assert(myData.getPointOffset(0,0)==3);

    DataArrayView myDataView = myData.getDataPoint(0,0);
    assert(myDataView==myView1);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==3);
    assert(myDataView.getRank()==1);
    assert(myDataView.noValues()==3);
    assert(myDataView.getShape().size()==1);
    assert(myDataView(0)==3);
    assert(myDataView(1)==4);
    assert(myDataView(2)==5);

    myDataView = myData.getDataPointByTag(1);
    assert(myDataView==myView1);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==3);
    assert(myDataView.getRank()==1);
    assert(myDataView.noValues()==3);
    assert(myDataView.getShape().size()==1);
    assert(myDataView(0)==3);
    assert(myDataView(1)==4);
    assert(myDataView(2)==5);

    myDataView = myData.getDefaultValue();
    assert(myDataView==myView);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==0);
    assert(myDataView.getRank()==1);
    assert(myDataView.noValues()==3);
    assert(myDataView.getShape().size()==1);
    assert(myDataView(0)==0);
    assert(myDataView(1)==1);
    assert(myDataView(2)==2);

    // use a non-existent tag so we get a pointer to
    // the first element of the data array
    double* sampleData=myData.getSampleDataByTag(9);
    for (int i=0; i<myData.getLength(); i++) {
      assert(sampleData[i]==i);
    }

  }

  {

    cout << "\tTest adding three keys with one value to DataTagged with default value only." << endl;

    DataArrayView::ShapeType viewShape;
    viewShape.push_back(3);

    DataTagged::TagListType keys;

    DataTagged::ValueListType values;

    DataArrayView::ValueType viewData(3);
    for (int i=0;i<viewShape[0];i++) {
      viewData[i]=i;
    }
    DataArrayView myView(viewData,viewShape);

    DataTagged myData(keys,values,myView,FunctionSpace());

    keys.push_back(1);
    keys.push_back(2);
    keys.push_back(3);

    DataArrayView::ValueType viewData1(3);
    for (int i=0;i<viewShape[0];i++) {
      viewData1[i]=3;
    }
    DataArrayView myView1(viewData1,viewShape);
    values.push_back(myView1);

    myData.addTaggedValues(keys,values);

    assert(myData.isCurrentTag(1));
    assert(myData.isCurrentTag(2));
    assert(myData.isCurrentTag(3));

    assert(myData.getTagLookup().size()==3);

    assert(myData.getLength()==12);

    // data-point 0 has tag number 1 by default
    assert(myData.getTagNumber(0)==1);

    assert(myData.getPointOffset(0,0)==3);

    DataArrayView myDataView = myData.getDataPoint(0,0);
    assert(myDataView==myView1);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==3);
    assert(myDataView.getRank()==1);
    assert(myDataView.noValues()==3);
    assert(myDataView.getShape().size()==1);
    assert(myDataView(0)==3);
    assert(myDataView(1)==3);
    assert(myDataView(2)==3);

    myDataView = myData.getDataPointByTag(1);
    assert(myDataView==myView1);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==3);
    assert(myDataView.getRank()==1);
    assert(myDataView.noValues()==3);
    assert(myDataView.getShape().size()==1);
    assert(myDataView(0)==3);
    assert(myDataView(1)==3);
    assert(myDataView(2)==3);

    myDataView = myData.getDataPointByTag(2);
    assert(myDataView==myView1);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==6);
    assert(myDataView.getRank()==1);
    assert(myDataView.noValues()==3);
    assert(myDataView.getShape().size()==1);
    assert(myDataView(0)==3);
    assert(myDataView(1)==3);
    assert(myDataView(2)==3);

    myDataView = myData.getDataPointByTag(3);
    assert(myDataView==myView1);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==9);
    assert(myDataView.getRank()==1);
    assert(myDataView.noValues()==3);
    assert(myDataView.getShape().size()==1);
    assert(myDataView(0)==3);
    assert(myDataView(1)==3);
    assert(myDataView(2)==3);

    myDataView = myData.getDefaultValue();
    assert(myDataView==myView);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==0);
    assert(myDataView.getRank()==1);
    assert(myDataView.noValues()==3);
    assert(myDataView.getShape().size()==1);
    assert(myDataView(0)==0);
    assert(myDataView(1)==1);
    assert(myDataView(2)==2);

    // use a non-existent tag so we get a pointer to
    // the first element of the data array
    double* sampleData=myData.getSampleDataByTag(9);
    for (int i=0; i<myData.getLength(); i++) {
      if (i<3) {
        assert(sampleData[i]==i);
      } else {
        assert(sampleData[i]==3);
      }
    }

  }

  {

    cout << "\tTest adding three keys with three values to DataTagged with default value only." << endl;

    DataArrayView::ShapeType viewShape;
    viewShape.push_back(3);

    DataTagged::TagListType keys;

    DataTagged::ValueListType values;

    DataArrayView::ValueType viewData(3);
    for (int i=0;i<viewShape[0];i++) {
      viewData[i]=i;
    }
    DataArrayView myView(viewData,viewShape);

    DataTagged myData(keys,values,myView,FunctionSpace());

    keys.push_back(1);
    keys.push_back(2);
    keys.push_back(3);

    DataArrayView::ValueType viewData1(3);
    for (int i=0;i<viewShape[0];i++) {
      viewData1[i]=i+1;
    }
    DataArrayView myView1(viewData1,viewShape);
    values.push_back(myView1);

    DataArrayView::ValueType viewData2(3);
    for (int i=0;i<viewShape[0];i++) {
      viewData2[i]=i+2;
    }
    DataArrayView myView2(viewData2,viewShape);
    values.push_back(myView2);

    DataArrayView::ValueType viewData3(3);
    for (int i=0;i<viewShape[0];i++) {
      viewData3[i]=i+3;
    }
    DataArrayView myView3(viewData3,viewShape);
    values.push_back(myView3);

    myData.addTaggedValues(keys,values);

    assert(myData.isCurrentTag(1));
    assert(myData.isCurrentTag(2));
    assert(myData.isCurrentTag(3));

    assert(myData.getTagLookup().size()==3);

    assert(myData.getLength()==12);

    // data-point 0 has tag number 1 by default
    assert(myData.getTagNumber(0)==1);

    assert(myData.getPointOffset(0,0)==3);

    DataArrayView myDataView = myData.getDataPoint(0,0);
    assert(myDataView==myView1);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==3);
    assert(myDataView.getRank()==1);
    assert(myDataView.noValues()==3);
    assert(myDataView.getShape().size()==1);
    assert(myDataView(0)==1);
    assert(myDataView(1)==2);
    assert(myDataView(2)==3);

    myDataView = myData.getDataPointByTag(1);
    assert(myDataView==myView1);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==3);
    assert(myDataView.getRank()==1);
    assert(myDataView.noValues()==3);
    assert(myDataView.getShape().size()==1);
    assert(myDataView(0)==1);
    assert(myDataView(1)==2);
    assert(myDataView(2)==3);

    myDataView = myData.getDataPointByTag(2);
    assert(myDataView==myView2);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==6);
    assert(myDataView.getRank()==1);
    assert(myDataView.noValues()==3);
    assert(myDataView.getShape().size()==1);
    assert(myDataView(0)==2);
    assert(myDataView(1)==3);
    assert(myDataView(2)==4);

    myDataView = myData.getDataPointByTag(3);
    assert(myDataView==myView3);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==9);
    assert(myDataView.getRank()==1);
    assert(myDataView.noValues()==3);
    assert(myDataView.getShape().size()==1);
    assert(myDataView(0)==3);
    assert(myDataView(1)==4);
    assert(myDataView(2)==5);

    myDataView = myData.getDefaultValue();
    assert(myDataView==myView);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==0);
    assert(myDataView.getRank()==1);
    assert(myDataView.noValues()==3);
    assert(myDataView.getShape().size()==1);
    assert(myDataView(0)==0);
    assert(myDataView(1)==1);
    assert(myDataView(2)==2);

    // use a non-existent tag so we get a pointer to
    // the first element of the data array
    double* sampleData=myData.getSampleDataByTag(9);
    for (int i=0; i<myData.getLength(); i++) {
      if (i<3) {
        assert(sampleData[i]==i);
      }
      if ((i>=3) && (i<6)) {
        assert(sampleData[i]==i-2);
      }
      if ((i>=6) && (i<9)) {
        assert(sampleData[i]==i-4);
      }
      if (i>=9) {
        assert(sampleData[i]==i-6);
      }
    }

  }

  {

    cout << "\tTest adding one key with empty value list to DataTagged with three tags." << endl;

    DataTagged::TagListType keys;
    keys.push_back(1);
    keys.push_back(2);
    keys.push_back(3);

    DataTagged::ValueListType values;

    DataArrayView::ShapeType viewShape;
    viewShape.push_back(3);

    // default value
    DataArrayView::ValueType viewData(3);
    for (int i=0;i<viewShape[0];i++) {
      viewData[i]=i;
    }
    DataArrayView myView(viewData,viewShape);

    // value for tag "1"
    DataArray eOne(myView);
    for (int i=0;i<eOne.getView().getShape()[0];i++) {
      eOne.getView()(i)=i+1.0;
    }
    values.push_back(eOne.getView());

    // value for tag "2"
    DataArray eTwo(myView);
    for (int i=0;i<eTwo.getView().getShape()[0];i++) {
      eTwo.getView()(i)=i+2.0;
    }
    values.push_back(eTwo.getView());

    // value for tag "3"
    DataArray eThree(myView);
    for (int i=0;i<eThree.getView().getShape()[0];i++) {
      eThree.getView()(i)=i+3.0;
    }
    values.push_back(eThree.getView());

    DataTagged myData(keys,values,myView,FunctionSpace());

    keys.clear();
    keys.push_back(4);
    values.clear();

    myData.addTaggedValues(keys,values);

    assert(myData.isCurrentTag(4));

    assert(myData.getTagLookup().size()==4);

    assert(myData.getLength()==15);

    DataArrayView myDataView = myData.getDataPointByTag(4);
    assert(myDataView==myView);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==12);
    assert(myDataView.getRank()==1);
    assert(myDataView.noValues()==3);
    assert(myDataView.getShape().size()==1);
    assert(myDataView(0)==0);
    assert(myDataView(1)==1);
    assert(myDataView(2)==2);

    // use a non-existent tag so we get a pointer to
    // the first element of the data array
    double* sampleData=myData.getSampleDataByTag(9);
    for (int i=0; i<myData.getLength(); i++) {
      if (i<3) {
        assert(sampleData[i]==i);
      }
      if ((i>=3) && (i<6)) {
        assert(sampleData[i]==i-2);
      }
      if ((i>=6) && (i<9)) {
        assert(sampleData[i]==i-4);
      }
      if ((i>=9) && (i<12)) {
        assert(sampleData[i]==i-6);
      }
      if (i>=12) {
        assert(sampleData[i]==i-12);
      }
    }

  }

  {

    cout << "\tTest adding one key with one value to DataTagged with three tags." << endl;

    DataTagged::TagListType keys;
    keys.push_back(1);
    keys.push_back(2);
    keys.push_back(3);

    DataTagged::ValueListType values;

    DataArrayView::ShapeType viewShape;
    viewShape.push_back(3);

    // default value
    DataArrayView::ValueType viewData(3);
    for (int i=0;i<viewShape[0];i++) {
      viewData[i]=i;
    }
    DataArrayView myView(viewData,viewShape);

    // value for tag "1"
    DataArray eOne(myView);
    for (int i=0;i<eOne.getView().getShape()[0];i++) {
      eOne.getView()(i)=i+1.0;
    }
    values.push_back(eOne.getView());

    // value for tag "2"
    DataArray eTwo(myView);
    for (int i=0;i<eTwo.getView().getShape()[0];i++) {
      eTwo.getView()(i)=i+2.0;
    }
    values.push_back(eTwo.getView());

    // value for tag "3"
    DataArray eThree(myView);
    for (int i=0;i<eThree.getView().getShape()[0];i++) {
      eThree.getView()(i)=i+3.0;
    }
    values.push_back(eThree.getView());

    DataTagged myData(keys,values,myView,FunctionSpace());

    keys.clear();
    keys.push_back(4);

    values.clear();
    // value for tag "4"
    DataArray eFour(myView);
    for (int i=0;i<eFour.getView().getShape()[0];i++) {
      eFour.getView()(i)=i+4.0;
    }
    values.push_back(eFour.getView());

    myData.addTaggedValues(keys,values);

    assert(myData.isCurrentTag(4));

    assert(myData.getTagLookup().size()==4);

    assert(myData.getLength()==15);

    DataArrayView myDataView = myData.getDataPointByTag(4);
    assert(myDataView==eFour.getView());
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==12);
    assert(myDataView.getRank()==1);
    assert(myDataView.noValues()==3);
    assert(myDataView.getShape().size()==1);
    assert(myDataView(0)==4);
    assert(myDataView(1)==5);
    assert(myDataView(2)==6);

    // use a non-existent tag so we get a pointer to
    // the first element of the data array
    double* sampleData=myData.getSampleDataByTag(9);
    for (int i=0; i<myData.getLength(); i++) {
      if (i<3) {
        assert(sampleData[i]==i);
      }
      if ((i>=3) && (i<6)) {
        assert(sampleData[i]==i-2);
      }
      if ((i>=6) && (i<9)) {
        assert(sampleData[i]==i-4);
      }
      if ((i>=9) && (i<12)) {
        assert(sampleData[i]==i-6);
      }
      if (i>=12) {
        assert(sampleData[i]==i-8);
      }
    }

  }

  {

    cout << "\tTest adding three keys with one value to DataTagged with three tags." << endl;

    DataTagged::TagListType keys;
    keys.push_back(1);
    keys.push_back(2);
    keys.push_back(3);

    DataTagged::ValueListType values;

    DataArrayView::ShapeType viewShape;
    viewShape.push_back(3);

    // default value
    DataArrayView::ValueType viewData(3);
    for (int i=0;i<viewShape[0];i++) {
      viewData[i]=i;
    }
    DataArrayView myView(viewData,viewShape);

    // value for tag "1"
    DataArray eOne(myView);
    for (int i=0;i<eOne.getView().getShape()[0];i++) {
      eOne.getView()(i)=i+1.0;
    }
    values.push_back(eOne.getView());

    // value for tag "2"
    DataArray eTwo(myView);
    for (int i=0;i<eTwo.getView().getShape()[0];i++) {
      eTwo.getView()(i)=i+2.0;
    }
    values.push_back(eTwo.getView());

    // value for tag "3"
    DataArray eThree(myView);
    for (int i=0;i<eThree.getView().getShape()[0];i++) {
      eThree.getView()(i)=i+3.0;
    }
    values.push_back(eThree.getView());

    DataTagged myData(keys,values,myView,FunctionSpace());

    keys.clear();
    keys.push_back(4);
    keys.push_back(5);
    keys.push_back(6);

    values.clear();
    // value for tags "4", "5" and "6"
    DataArray eFour(myView);
    for (int i=0;i<eFour.getView().getShape()[0];i++) {
      eFour.getView()(i)=i+4.0;
    }
    values.push_back(eFour.getView());

    myData.addTaggedValues(keys,values);

    assert(myData.isCurrentTag(4));
    assert(myData.isCurrentTag(5));
    assert(myData.isCurrentTag(6));

    assert(myData.getTagLookup().size()==6);

    assert(myData.getLength()==21);

    DataArrayView myDataView = myData.getDataPointByTag(4);
    assert(myDataView==eFour.getView());
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==12);
    assert(myDataView.getRank()==1);
    assert(myDataView.noValues()==3);
    assert(myDataView.getShape().size()==1);
    assert(myDataView(0)==4);
    assert(myDataView(1)==5);
    assert(myDataView(2)==6);

    myDataView = myData.getDataPointByTag(5);
    assert(myDataView==eFour.getView());
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==15);
    assert(myDataView.getRank()==1);
    assert(myDataView.noValues()==3);
    assert(myDataView.getShape().size()==1);
    assert(myDataView(0)==4);
    assert(myDataView(1)==5);
    assert(myDataView(2)==6);

    myDataView = myData.getDataPointByTag(6);
    assert(myDataView==eFour.getView());
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==18);
    assert(myDataView.getRank()==1);
    assert(myDataView.noValues()==3);
    assert(myDataView.getShape().size()==1);
    assert(myDataView(0)==4);
    assert(myDataView(1)==5);
    assert(myDataView(2)==6);

    // use a non-existent tag so we get a pointer to
    // the first element of the data array
    double* sampleData=myData.getSampleDataByTag(9);
    for (int i=0; i<myData.getLength(); i++) {
      if (i<3) {
        assert(sampleData[i]==i);
      }
      if ((i>=3) && (i<6)) {
        assert(sampleData[i]==i-2);
      }
      if ((i>=6) && (i<9)) {
        assert(sampleData[i]==i-4);
      }
      if ((i>=9) && (i<12)) {
        assert(sampleData[i]==i-6);
      }
      if ((i>=12) && (i<15)) {
        assert(sampleData[i]==i-8);
      }
      if ((i>=15) && (i<18)) {
        assert(sampleData[i]==i-11);
      }
      if (i>=18) {
        assert(sampleData[i]==i-14);
      }
    }

  }

  {

    cout << "\tTest adding three keys with three values to DataTagged with three tags." << endl;

    DataTagged::TagListType keys;
    keys.push_back(1);
    keys.push_back(2);
    keys.push_back(3);

    DataTagged::ValueListType values;

    DataArrayView::ShapeType viewShape;
    viewShape.push_back(3);

    // default value
    DataArrayView::ValueType viewData(3);
    for (int i=0;i<viewShape[0];i++) {
      viewData[i]=i;
    }
    DataArrayView myView(viewData,viewShape);

    // value for tag "1"
    DataArray eOne(myView);
    for (int i=0;i<eOne.getView().getShape()[0];i++) {
      eOne.getView()(i)=i+1.0;
    }
    values.push_back(eOne.getView());

    // value for tag "2"
    DataArray eTwo(myView);
    for (int i=0;i<eTwo.getView().getShape()[0];i++) {
      eTwo.getView()(i)=i+2.0;
    }
    values.push_back(eTwo.getView());

    // value for tag "3"
    DataArray eThree(myView);
    for (int i=0;i<eThree.getView().getShape()[0];i++) {
      eThree.getView()(i)=i+3.0;
    }
    values.push_back(eThree.getView());

    DataTagged myData(keys,values,myView,FunctionSpace());

    keys.clear();
    keys.push_back(4);
    keys.push_back(5);
    keys.push_back(6);

    values.clear();

    // value for tag "4"
    DataArray eFour(myView);
    for (int i=0;i<eFour.getView().getShape()[0];i++) {
      eFour.getView()(i)=i+4.0;
    }
    values.push_back(eFour.getView());

    // value for tag "5"
    DataArray eFive(myView);
    for (int i=0;i<eFive.getView().getShape()[0];i++) {
      eFive.getView()(i)=i+5.0;
    }
    values.push_back(eFive.getView());

    // value for tag "6"
    DataArray eSix(myView);
    for (int i=0;i<eSix.getView().getShape()[0];i++) {
      eSix.getView()(i)=i+6.0;
    }
    values.push_back(eSix.getView());

    myData.addTaggedValues(keys,values);

    assert(myData.isCurrentTag(4));
    assert(myData.isCurrentTag(5));
    assert(myData.isCurrentTag(6));

    assert(myData.getTagLookup().size()==6);

    assert(myData.getLength()==21);

    DataArrayView myDataView = myData.getDataPointByTag(4);
    assert(myDataView==eFour.getView());
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==12);
    assert(myDataView.getRank()==1);
    assert(myDataView.noValues()==3);
    assert(myDataView.getShape().size()==1);
    assert(myDataView(0)==4);
    assert(myDataView(1)==5);
    assert(myDataView(2)==6);

    myDataView = myData.getDataPointByTag(5);
    assert(myDataView==eFive.getView());
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==15);
    assert(myDataView.getRank()==1);
    assert(myDataView.noValues()==3);
    assert(myDataView.getShape().size()==1);
    assert(myDataView(0)==5);
    assert(myDataView(1)==6);
    assert(myDataView(2)==7);

    myDataView = myData.getDataPointByTag(6);
    assert(myDataView==eSix.getView());
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==18);
    assert(myDataView.getRank()==1);
    assert(myDataView.noValues()==3);
    assert(myDataView.getShape().size()==1);
    assert(myDataView(0)==6);
    assert(myDataView(1)==7);
    assert(myDataView(2)==8);

    // use a non-existent tag so we get a pointer to
    // the first element of the data array
    double* sampleData=myData.getSampleDataByTag(9);
    for (int i=0; i<myData.getLength(); i++) {
      if (i<3) {
        assert(sampleData[i]==i);
      }
      if ((i>=3) && (i<6)) {
        assert(sampleData[i]==i-2);
      }
      if ((i>=6) && (i<9)) {
        assert(sampleData[i]==i-4);
      }
      if ((i>=9) && (i<12)) {
        assert(sampleData[i]==i-6);
      }
      if ((i>=12) && (i<15)) {
        assert(sampleData[i]==i-8);
      }
      if ((i>=15) && (i<18)) {
        assert(sampleData[i]==i-10);
      }
      if (i>=18) {
        assert(sampleData[i]==i-12);
      }
    }

  }

}

void DataTaggedTestCase::testSetTaggedValue() {

  cout << endl;

  {

    cout << "\tTest setting key in DataTagged with three tags." << endl;

    DataTagged::TagListType keys;
    keys.push_back(1);
    keys.push_back(2);
    keys.push_back(3);

    DataTagged::ValueListType values;

    DataArrayView::ShapeType viewShape;
    viewShape.push_back(3);

    // default value
    DataArrayView::ValueType viewData(3);
    for (int i=0;i<viewShape[0];i++) {
      viewData[i]=i;
    }
    DataArrayView myView(viewData,viewShape);

    // value for tag "1"
    DataArray eOne(myView);
    for (int i=0;i<eOne.getView().getShape()[0];i++) {
      eOne.getView()(i)=i+1.0;
    }
    values.push_back(eOne.getView());

    // value for tag "2"
    DataArray eTwo(myView);
    for (int i=0;i<eTwo.getView().getShape()[0];i++) {
      eTwo.getView()(i)=i+2.0;
    }
    values.push_back(eTwo.getView());

    // value for tag "3"
    DataArray eThree(myView);
    for (int i=0;i<eThree.getView().getShape()[0];i++) {
      eThree.getView()(i)=i+3.0;
    }
    values.push_back(eThree.getView());

    DataTagged myData(keys,values,myView,FunctionSpace());

    // new value for tag "2"
    for (int i=0;i<eTwo.getView().getShape()[0];i++) {
      eTwo.getView()(i)=i+5.0;
    }

    myData.setTaggedValue(2,eTwo.getView());

    assert(myData.isCurrentTag(2));

    assert(myData.getTagLookup().size()==3);

    assert(myData.getLength()==12);

    DataArrayView myDataView = myData.getDataPointByTag(2);
    assert(myDataView==eTwo.getView());
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==6);
    assert(myDataView.getRank()==1);
    assert(myDataView.noValues()==3);
    assert(myDataView.getShape().size()==1);
    assert(myDataView(0)==5);
    assert(myDataView(1)==6);
    assert(myDataView(2)==7);

    // use a non-existent tag so we get a pointer to
    // the first element of the data array
    double* sampleData=myData.getSampleDataByTag(9);
    for (int i=0; i<myData.getLength(); i++) {
      if (i<3) {
        assert(sampleData[i]==i);
      }
      if ((i>=3) && (i<6)) {
        assert(sampleData[i]==i-2);
      }
      if ((i>=6) && (i<9)) {
        assert(sampleData[i]==i-1);
      }
      if ((i>=9) && (i<12)) {
        assert(sampleData[i]==i-6);
      }
    }

  }

}

void DataTaggedTestCase::testAll() {

  cout << endl;

  {

    cout << "\tTest default DataTagged." << endl;
    DataTagged myData;

    //cout << myData.toString() << endl;

    assert(myData.getNumSamples()==1);
    assert(myData.getNumDPPSample()==1);

    assert(myData.validSamplePointNo(0));
    assert(myData.validSampleNo(0));
    assert(!myData.validSamplePointNo(1));
    assert(!myData.validSampleNo(1));

    // data-point 0 has tag number 1 by default
    assert(myData.getTagNumber(0)==1);

    assert(!myData.isCurrentTag(1));

    assert(myData.getTagLookup().size()==0);

    assert(myData.getLength()==1);

    assert(myData.getPointOffset(0,0)==0);

    DataArrayView myDataView = myData.getDataPoint(0,0);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==0);
    assert(myDataView.getRank()==0);
    assert(myDataView.noValues()==1);
    assert(myDataView.getShape().size()==0);
    assert(myDataView()==0.0);

    // Test non-existent tag returns the default value.
    myDataView = myData.getDataPointByTag(1);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==0);
    assert(myDataView.getRank()==0);
    assert(myDataView.noValues()==1);
    assert(myDataView.getShape().size()==0);
    assert(myDataView()==0.0);

    myDataView = myData.getDefaultValue();
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==0);
    assert(myDataView.getRank()==0);
    assert(myDataView.noValues()==1);
    assert(myDataView.getShape().size()==0);
    assert(myDataView()==0.0);

    // use a non-existent tag so we get a pointer to
    // the first element of the data array
    double* sampleData=myData.getSampleDataByTag(9);
    for (int i=0; i<myData.getLength(); i++) {
      assert(sampleData[i]==i);
    }

  }

  {

    cout << "\tTest DataTagged with default value only." << endl;

    DataArrayView::ShapeType viewShape;
    viewShape.push_back(3);

    DataTagged::TagListType keys;

    DataTagged::ValueListType values;

    DataArrayView::ValueType viewData(3);
    for (int i=0;i<viewShape[0];i++) {
      viewData[i]=i;
    }
    DataArrayView myView(viewData,viewShape);

    DataTagged myData(keys,values,myView,FunctionSpace());

    //cout << myData.toString() << endl;

    assert(myData.getNumSamples()==1);
    assert(myData.getNumDPPSample()==1);

    assert(myData.validSamplePointNo(0));
    assert(myData.validSampleNo(0));
    assert(!myData.validSamplePointNo(1));
    assert(!myData.validSampleNo(1));

    // data-point 0 has tag number 1 by default
    assert(myData.getTagNumber(0)==1);

    assert(!myData.isCurrentTag(1));

    assert(myData.getTagLookup().size()==0);

    assert(myData.getLength()==3);

    assert(myData.getPointOffset(0,0)==0);

    DataArrayView myDataView = myData.getDataPoint(0,0);
    assert(myDataView==myView);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==0);
    assert(myDataView.getRank()==1);
    assert(myDataView.noValues()==3);
    assert(myDataView.getShape().size()==1);
    assert(myDataView(0)==0);
    assert(myDataView(1)==1);
    assert(myDataView(2)==2);

    // Test non-existent tag returns the default value.
    myDataView = myData.getDataPointByTag(1);
    assert(myDataView==myView);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==0);
    assert(myDataView.getRank()==1);
    assert(myDataView.noValues()==3);
    assert(myDataView.getShape().size()==1);
    assert(myDataView(0)==0);
    assert(myDataView(1)==1);
    assert(myDataView(2)==2);

    myDataView = myData.getDefaultValue();
    assert(myDataView==myView);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==0);
    assert(myDataView.getRank()==1);
    assert(myDataView.noValues()==3);
    assert(myDataView.getShape().size()==1);
    assert(myDataView(0)==0);
    assert(myDataView(1)==1);
    assert(myDataView(2)==2);

    // use a non-existent tag so we get a pointer to
    // the first element of the data array
    double* sampleData=myData.getSampleDataByTag(9);
    for (int i=0; i<myData.getLength(); i++) {
      assert(sampleData[i]==i);
    }

  }

  {

    cout << "\tTest DataTagged with one tag." << endl;

    // the one data-point has tag value "1"

    DataTagged::TagListType keys;
    keys.push_back(1);

    DataTagged::ValueListType values;

    DataArrayView::ShapeType viewShape;
    viewShape.push_back(3);

    // default value
    DataArrayView::ValueType viewData(3);
    for (int i=0;i<viewShape[0];i++) {
      viewData[i]=i;
    }
    DataArrayView myView(viewData,viewShape);

    // value for tag "1"
    DataArray eOne(myView);
    for (int i=0;i<eOne.getView().getShape()[0];i++) {
      eOne.getView()(i)=i+1.0;
    }
    values.push_back(eOne.getView());

    DataTagged myData(keys,values,myView,FunctionSpace());

    //cout << myData.toString() << endl;

    assert(myData.getNumSamples()==1);
    assert(myData.getNumDPPSample()==1);

    assert(myData.validSamplePointNo(0));
    assert(myData.validSampleNo(0));
    assert(!myData.validSamplePointNo(1));
    assert(!myData.validSampleNo(1));

    // data-point 0 has tag number 1 by default
    assert(myData.getTagNumber(0)==1);

    assert(!myData.isCurrentTag(0));
    assert(myData.isCurrentTag(1));

    assert(myData.getTagLookup().size()==1);

    assert(myData.getLength()==6);

    assert(myData.getPointOffset(0,0)==3);

    DataArrayView myDataView = myData.getDataPoint(0,0);
    assert(myDataView==eOne.getView());
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==3);
    assert(myDataView.getRank()==1);
    assert(myDataView.noValues()==3);
    assert(myDataView.getShape().size()==1);
    assert(myDataView(0)==1);
    assert(myDataView(1)==2);
    assert(myDataView(2)==3);

    myDataView = myData.getDataPointByTag(1);
    assert(myDataView==eOne.getView());
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==3);
    assert(myDataView.getRank()==1);
    assert(myDataView.noValues()==3);
    assert(myDataView.getShape().size()==1);
    assert(myDataView(0)==1);
    assert(myDataView(1)==2);
    assert(myDataView(2)==3);

    // Test non-existent tag returns the default value.
    myDataView = myData.getDataPointByTag(0);
    assert(myDataView==myView);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==0);
    assert(myDataView.getRank()==1);
    assert(myDataView.noValues()==3);
    assert(myDataView.getShape().size()==1);
    assert(myDataView(0)==0);
    assert(myDataView(1)==1);
    assert(myDataView(2)==2);

    myDataView = myData.getDefaultValue();
    assert(myDataView==myView);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==0);
    assert(myDataView.getRank()==1);
    assert(myDataView.noValues()==3);
    assert(myDataView.getShape().size()==1);
    assert(myDataView(0)==0);
    assert(myDataView(1)==1);
    assert(myDataView(2)==2);

    // use a non-existent tag so we get a pointer to
    // the first element of the data array
    double* sampleData=myData.getSampleDataByTag(9);
    for (int i=0; i<myData.getLength(); i++) {
      if (i<3) {
        assert(sampleData[i]==i);
      }
      if (i>=3) {
        assert(sampleData[i]==i-2);
      }
    }

  }

  {

    cout << "\tTest DataTagged with multiple tags." << endl;

    // the one data-point has tag value "1"

    DataTagged::TagListType keys;
    keys.push_back(1);
    keys.push_back(2);
    keys.push_back(3);

    DataTagged::ValueListType values;

    DataArrayView::ShapeType viewShape;
    viewShape.push_back(3);

    // default value
    DataArrayView::ValueType viewData(3);
    for (int i=0;i<viewShape[0];i++) {
      viewData[i]=i;
    }
    DataArrayView myView(viewData,viewShape);

    // value for tag "1"
    DataArray eOne(myView);
    for (int i=0;i<eOne.getView().getShape()[0];i++) {
      eOne.getView()(i)=i+1.0;
    }
    values.push_back(eOne.getView());

    // value for tag "2"
    DataArray eTwo(myView);
    for (int i=0;i<eTwo.getView().getShape()[0];i++) {
      eTwo.getView()(i)=i+2.0;
    }
    values.push_back(eTwo.getView());

    // value for tag "3"
    DataArray eThree(myView);
    for (int i=0;i<eThree.getView().getShape()[0];i++) {
      eThree.getView()(i)=i+3.0;
    }
    values.push_back(eThree.getView());

    DataTagged myData(keys,values,myView,FunctionSpace());

    //cout << myData.toString() << endl;

    assert(myData.getNumSamples()==1);
    assert(myData.getNumDPPSample()==1);

    assert(myData.validSamplePointNo(0));
    assert(myData.validSampleNo(0));
    assert(!myData.validSamplePointNo(1));
    assert(!myData.validSampleNo(1));

    // data-point 0 has tag number 1 by default
    assert(myData.getTagNumber(0)==1);

    assert(!myData.isCurrentTag(0));
    assert(myData.isCurrentTag(1));
    assert(myData.isCurrentTag(2));
    assert(myData.isCurrentTag(3));

    assert(myData.getTagLookup().size()==3);

    assert(myData.getLength()==12);

    assert(myData.getPointOffset(0,0)==3);

    DataArrayView myDataView = myData.getDataPoint(0,0);
    assert(myDataView==eOne.getView());
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==3);
    assert(myDataView.getRank()==1);
    assert(myDataView.noValues()==3);
    assert(myDataView.getShape().size()==1);
    assert(myDataView(0)==1);
    assert(myDataView(1)==2);
    assert(myDataView(2)==3);

    myDataView = myData.getDataPointByTag(1);
    assert(myDataView==eOne.getView());
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==3);
    assert(myDataView.getRank()==1);
    assert(myDataView.noValues()==3);
    assert(myDataView.getShape().size()==1);
    assert(myDataView(0)==1);
    assert(myDataView(1)==2);
    assert(myDataView(2)==3);

    // Test non-existent tag returns the default value.
    myDataView = myData.getDataPointByTag(0);
    assert(myDataView==myView);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==0);
    assert(myDataView.getRank()==1);
    assert(myDataView.noValues()==3);
    assert(myDataView.getShape().size()==1);
    assert(myDataView(0)==0);
    assert(myDataView(1)==1);
    assert(myDataView(2)==2);

    myDataView = myData.getDefaultValue();
    assert(myDataView==myView);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==0);
    assert(myDataView.getRank()==1);
    assert(myDataView.noValues()==3);
    assert(myDataView.getShape().size()==1);
    assert(myDataView(0)==0);
    assert(myDataView(1)==1);
    assert(myDataView(2)==2);

    // Test data-points held for remaining tags
    myDataView = myData.getDataPointByTag(2);
    assert(myDataView==eTwo.getView());
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==6);
    assert(myDataView.getRank()==1);
    assert(myDataView.noValues()==3);
    assert(myDataView.getShape().size()==1);
    assert(myDataView(0)==2);
    assert(myDataView(1)==3);
    assert(myDataView(2)==4);

    myDataView = myData.getDataPointByTag(3);
    assert(myDataView==eThree.getView());
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==9);
    assert(myDataView.getRank()==1);
    assert(myDataView.noValues()==3);
    assert(myDataView.getShape().size()==1);
    assert(myDataView(0)==3);
    assert(myDataView(1)==4);
    assert(myDataView(2)==5);

    // use a non-existent tag so we get a pointer to
    // the first element of the data array
    double* sampleData=myData.getSampleDataByTag(9);
    for (int i=0; i<myData.getLength(); i++) {
      if (i<3) {
        assert(sampleData[i]==i);
      }
      if ((i>=3) && (i<6)) {
        assert(sampleData[i]==i-2);
      }
      if ((i>=6) && (i<9)) {
        assert(sampleData[i]==i-4);
      }
      if (i>=9) {
        assert(sampleData[i]==i-6);
      }
    }

  }

}

void DataTaggedTestCase::testCopyConstructors() {

  cout << endl;

  {

    cout << "\tTest DataTagged copy constructor for DataTagged with multiple tags." << endl;

    // the one data-point has tag value "1"

    DataTagged::TagListType keys;
    keys.push_back(1);
    keys.push_back(2);
    keys.push_back(3);

    DataTagged::ValueListType values;

    DataArrayView::ShapeType viewShape;
    viewShape.push_back(3);

    // default value
    DataArrayView::ValueType viewData(3);
    for (int i=0;i<viewShape[0];i++) {
      viewData[i]=i;
    }
    DataArrayView myView(viewData,viewShape);

    // value for tag "1"
    DataArray eOne(myView);
    for (int i=0;i<eOne.getView().getShape()[0];i++) {
      eOne.getView()(i)=i+1.0;
    }
    values.push_back(eOne.getView());

    // value for tag "2"
    DataArray eTwo(myView);
    for (int i=0;i<eTwo.getView().getShape()[0];i++) {
      eTwo.getView()(i)=i+2.0;
    }
    values.push_back(eTwo.getView());

    // value for tag "3"
    DataArray eThree(myView);
    for (int i=0;i<eThree.getView().getShape()[0];i++) {
      eThree.getView()(i)=i+3.0;
    }
    values.push_back(eThree.getView());

    DataTagged myData(keys,values,myView,FunctionSpace());

    DataTagged myDataCopy(myData);

    //cout << myDataCopy.toString() << endl;

    assert(myDataCopy.getNumSamples()==1);
    assert(myDataCopy.getNumDPPSample()==1);

    assert(myDataCopy.validSamplePointNo(0));
    assert(myDataCopy.validSampleNo(0));
    assert(!myDataCopy.validSamplePointNo(1));
    assert(!myDataCopy.validSampleNo(1));

    // data-point 0 has tag number 1 by default
    assert(myDataCopy.getTagNumber(0)==1);

    assert(!myDataCopy.isCurrentTag(0));
    assert(myDataCopy.isCurrentTag(1));
    assert(myDataCopy.isCurrentTag(2));
    assert(myDataCopy.isCurrentTag(3));

    assert(myDataCopy.getTagLookup().size()==3);

    assert(myDataCopy.getLength()==12);

    assert(myDataCopy.getPointOffset(0,0)==3);

    DataArrayView myDataView = myDataCopy.getDataPoint(0,0);
    assert(myDataView==eOne.getView());
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==3);
    assert(myDataView.getRank()==1);
    assert(myDataView.noValues()==3);
    assert(myDataView.getShape().size()==1);
    assert(myDataView(0)==1);
    assert(myDataView(1)==2);
    assert(myDataView(2)==3);

    myDataView = myDataCopy.getDataPointByTag(1);
    assert(myDataView==eOne.getView());
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==3);
    assert(myDataView.getRank()==1);
    assert(myDataView.noValues()==3);
    assert(myDataView.getShape().size()==1);
    assert(myDataView(0)==1);
    assert(myDataView(1)==2);
    assert(myDataView(2)==3);

    // Test non-existent tag returns the default value.
    myDataView = myDataCopy.getDataPointByTag(0);
    assert(myDataView==myView);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==0);
    assert(myDataView.getRank()==1);
    assert(myDataView.noValues()==3);
    assert(myDataView.getShape().size()==1);
    assert(myDataView(0)==0);
    assert(myDataView(1)==1);
    assert(myDataView(2)==2);

    myDataView = myDataCopy.getDefaultValue();
    assert(myDataView==myView);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==0);
    assert(myDataView.getRank()==1);
    assert(myDataView.noValues()==3);
    assert(myDataView.getShape().size()==1);
    assert(myDataView(0)==0);
    assert(myDataView(1)==1);
    assert(myDataView(2)==2);

    // Test data-points held for remaining tags
    myDataView = myDataCopy.getDataPointByTag(2);
    assert(myDataView==eTwo.getView());
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==6);
    assert(myDataView.getRank()==1);
    assert(myDataView.noValues()==3);
    assert(myDataView.getShape().size()==1);
    assert(myDataView(0)==2);
    assert(myDataView(1)==3);
    assert(myDataView(2)==4);

    myDataView = myDataCopy.getDataPointByTag(3);
    assert(myDataView==eThree.getView());
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==9);
    assert(myDataView.getRank()==1);
    assert(myDataView.noValues()==3);
    assert(myDataView.getShape().size()==1);
    assert(myDataView(0)==3);
    assert(myDataView(1)==4);
    assert(myDataView(2)==5);

    // use a non-existent tag so we get a pointer to
    // the first element of the data array
    double* sampleData=myDataCopy.getSampleDataByTag(9);
    for (int i=0; i<myData.getLength(); i++) {
      if (i<3) {
        assert(sampleData[i]==i);
      }
      if ((i>=3) && (i<6)) {
        assert(sampleData[i]==i-2);
      }
      if ((i>=6) && (i<9)) {
        assert(sampleData[i]==i-4);
      }
      if (i>=9) {
        assert(sampleData[i]==i-6);
      }
    }

  }

  {

    cout << "\tTest DataTagged copy constructor for DataConstant." << endl;

    // Create a DataConstant
    DataArrayView::ShapeType shape;
    DataArrayView::ValueType data(DataArrayView::noValues(shape),0);
    DataArrayView pointData(data,shape);
    pointData()=1.0;
    DataConstant myConstantData(pointData, FunctionSpace());

    // use this DataConstant to initialise a DataTagged
    DataTagged myData(myConstantData);

    //cout << myData.toString() << endl;

    assert(myData.getNumSamples()==1);
    assert(myData.getNumDPPSample()==1);

    assert(myData.validSamplePointNo(0));
    assert(myData.validSampleNo(0));
    assert(!myData.validSamplePointNo(1));
    assert(!myData.validSampleNo(1));

    // data-point 0 has tag number 1 by default
    assert(myData.getTagNumber(0)==1);

    assert(!myData.isCurrentTag(1));

    assert(myData.getTagLookup().size()==0);

    assert(myData.getLength()==1);

    assert(myData.getPointOffset(0,0)==0);

    DataArrayView myDataView = myData.getDataPoint(0,0);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==0);
    assert(myDataView.getRank()==0);
    assert(myDataView.noValues()==1);
    assert(myDataView.getShape().size()==0);
    assert(myDataView()==1.0);

    // Test non-existent tag returns the default value.
    myDataView = myData.getDataPointByTag(1);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==0);
    assert(myDataView.getRank()==0);
    assert(myDataView.noValues()==1);
    assert(myDataView.getShape().size()==0);
    assert(myDataView()==1.0);

    myDataView = myData.getDefaultValue();
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==0);
    assert(myDataView.getRank()==0);
    assert(myDataView.noValues()==1);
    assert(myDataView.getShape().size()==0);
    assert(myDataView()==1.0);

    // use a non-existent tag so we get a pointer to
    // the first element of the data array
    double* sampleData=myData.getSampleDataByTag(9);
    for (int i=0; i<myData.getLength(); i++) {
      assert(sampleData[i]==i+1);
    }

  }

}

TestSuite* DataTaggedTestCase::suite ()
{
  //
  // create the suite of tests to perform.
  TestSuite *testSuite = new TestSuite ("DataTaggedTestCase");
  testSuite->addTest (new TestCaller< DataTaggedTestCase>("testAll",&DataTaggedTestCase::testAll));
  testSuite->addTest (new TestCaller< DataTaggedTestCase>("testAddTaggedValues",&DataTaggedTestCase::testAddTaggedValues));
  testSuite->addTest (new TestCaller< DataTaggedTestCase>("testSetTaggedValue",&DataTaggedTestCase::testSetTaggedValue));
  testSuite->addTest (new TestCaller< DataTaggedTestCase>("testCopyConstructors",&DataTaggedTestCase::testCopyConstructors));
  testSuite->addTest (new TestCaller< DataTaggedTestCase>("testOperations",&DataTaggedTestCase::testOperations));
//  testSuite->addTest (new TestCaller< DataTaggedTestCase>("testReshape",&DataTaggedTestCase::testReshape));
  return testSuite;
}
