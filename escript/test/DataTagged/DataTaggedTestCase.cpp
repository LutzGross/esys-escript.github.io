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

    cout << "\tTest rank 1 reshape of default DataTagged." << endl;

    DataTagged myData;
    myData.getPointDataView()()=1.0;

    DataArrayView::ShapeType shape;
    shape.push_back(2);

    myData.reshapeDataPoint(shape);

    for (int i=0;i<shape[0];i++) {
      assert(myData.getDefaultValue()(i)==1);
    }

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

    assert(myData.getLength()==2);

    assert(myData.getPointOffset(0,0)==0);

    DataArrayView myDataView = myData.getDataPoint(0,0);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==0);
    assert(myDataView.getRank()==1);
    assert(myDataView.noValues()==2);
    assert(myDataView.getShape().size()==1);
    assert(myDataView(0)==1.0);
    assert(myDataView(1)==1.0);

    // Test non-existent tag returns the default value.
    myDataView = myData.getDataPointByTag(1);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==0);
    assert(myDataView.getRank()==1);
    assert(myDataView.noValues()==2);
    assert(myDataView.getShape().size()==1);
    assert(myDataView(0)==1.0);
    assert(myDataView(1)==1.0);

    myDataView = myData.getDefaultValue();
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==0);
    assert(myDataView.getRank()==1);
    assert(myDataView.noValues()==2);
    assert(myDataView.getShape().size()==1);
    assert(myDataView(0)==1.0);
    assert(myDataView(1)==1.0);

    // use a non-existent tag so we get a pointer to
    // the first element of the data array
    double* sampleData=myData.getSampleDataByTag(9);
    for (int i=0; i<myData.getLength(); i++) {
      assert(sampleData[i]==1.0);
    }

  }

  {

    cout << "\tTest rank 2 reshape of DataTagged with one tag." << endl;

    DataTagged::TagListType keys;
    keys.push_back(1);

    DataTagged::ValueListType values;

    DataArrayView::ShapeType viewShape;

    // default value
    DataArrayView::ValueType viewData(1);
    viewData[0]=1.0;
    DataArrayView myView(viewData,viewShape);

    // value for tag "1"
    DataArray eOne(myView);
    eOne.getView()()=2.0;
    values.push_back(eOne.getView());

    DataTagged myData(keys,values,myView,FunctionSpace());

    DataArrayView::ShapeType shape;
    shape.push_back(2);
    shape.push_back(5);

    myData.reshapeDataPoint(shape);

    //cout << myData.toString() << endl;

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

    assert(myData.getLength()==20);

    assert(myData.getPointOffset(0,0)==10);

    DataArrayView myDataView = myData.getDataPoint(0,0);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==10);
    assert(myDataView.getRank()==2);
    assert(myDataView.noValues()==10);
    assert(myDataView.getShape().size()==2);
    for (int j=0;j<shape[1];j++) {
      for (int i=0;i<shape[0];i++) {
        assert(myDataView(i,j)==2.0);
      }
    }

    myDataView = myData.getDataPointByTag(1);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==10);
    assert(myDataView.getRank()==2);
    assert(myDataView.noValues()==10);
    assert(myDataView.getShape().size()==2);
    for (int j=0;j<shape[1];j++) {
      for (int i=0;i<shape[0];i++) {
        assert(myDataView(i,j)==2.0);
      }
    }

    myDataView = myData.getDefaultValue();
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==0);
    assert(myDataView.getRank()==2);
    assert(myDataView.noValues()==10);
    assert(myDataView.getShape().size()==2);
    for (int j=0;j<shape[1];j++) {
      for (int i=0;i<shape[0];i++) {
        assert(myDataView(i,j)==1.0);
      }
    }

    // use a non-existent tag so we get a pointer to
    // the first element of the data array
    double* sampleData=myData.getSampleDataByTag(9);
    for (int i=0; i<myData.getLength(); i++) {
      if (i<10) {
        assert(sampleData[i]==1.0);
      } else {
        assert(sampleData[i]==2.0);
      }
    }

  }

  {

    cout << "\tTest rank 3 reshape of DataTagged with three tags." << endl;

    DataTagged::TagListType keys;
    keys.push_back(1);
    keys.push_back(2);
    keys.push_back(3);

    DataTagged::ValueListType values;

    DataArrayView::ShapeType viewShape;

    // default value
    DataArrayView::ValueType viewData(1);
    viewData[0]=0.0;
    DataArrayView myView(viewData,viewShape);

    // value for tag "1"
    DataArray eOne(myView);
    eOne.getView()()=1.0;
    values.push_back(eOne.getView());

    // value for tag "2"
    DataArray eTwo(myView);
    eTwo.getView()()=2.0;
    values.push_back(eTwo.getView());

    // value for tag "3"
    DataArray eThree(myView);
    eThree.getView()()=3.0;
    values.push_back(eThree.getView());

    DataTagged myData(keys,values,myView,FunctionSpace());

    DataArrayView::ShapeType shape;
    shape.push_back(2);
    shape.push_back(2);
    shape.push_back(2);

    myData.reshapeDataPoint(shape);

    //cout << myData.toString() << endl;

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

    assert(myData.getLength()==32);

    assert(myData.getPointOffset(0,0)==8);

    DataArrayView myDataView = myData.getDataPoint(0,0);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==8);
    assert(myDataView.getRank()==3);
    assert(myDataView.noValues()==8);
    assert(myDataView.getShape().size()==3);
    for (int k=0;k<shape[2];k++) {
      for (int j=0;j<shape[1];j++) {
        for (int i=0;i<shape[0];i++) {
          assert(myDataView(i,j,k)==1.0);
        }
      }
    }

    myDataView = myData.getDataPointByTag(1);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==8);
    assert(myDataView.getRank()==3);
    assert(myDataView.noValues()==8);
    assert(myDataView.getShape().size()==3);
    for (int k=0;k<shape[2];k++) {
      for (int j=0;j<shape[1];j++) {
        for (int i=0;i<shape[0];i++) {
          assert(myDataView(i,j,k)==1.0);
        }
      }
    }

    myDataView = myData.getDataPointByTag(2);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==16);
    assert(myDataView.getRank()==3);
    assert(myDataView.noValues()==8);
    assert(myDataView.getShape().size()==3);
    for (int k=0;k<shape[2];k++) {
      for (int j=0;j<shape[1];j++) {
        for (int i=0;i<shape[0];i++) {
          assert(myDataView(i,j,k)==2.0);
        }
      }
    }

    myDataView = myData.getDataPointByTag(3);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==24);
    assert(myDataView.getRank()==3);
    assert(myDataView.noValues()==8);
    assert(myDataView.getShape().size()==3);
    for (int k=0;k<shape[2];k++) {
      for (int j=0;j<shape[1];j++) {
        for (int i=0;i<shape[0];i++) {
          assert(myDataView(i,j,k)==3.0);
        }
      }
    }

    myDataView = myData.getDefaultValue();
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==0);
    assert(myDataView.getRank()==3);
    assert(myDataView.noValues()==8);
    assert(myDataView.getShape().size()==3);
    for (int k=0;k<shape[2];k++) {
      for (int j=0;j<shape[1];j++) {
        for (int i=0;i<shape[0];i++) {
          assert(myDataView(i,j,k)==0.0);
        }
      }
    }

    // use a non-existent tag so we get a pointer to
    // the first element of the data array
    double* sampleData=myData.getSampleDataByTag(9);
    for (int i=0; i<myData.getLength(); i++) {
      if (i<8) {
        assert(sampleData[i]==0.0);
      } else if ((i>=8) && (i<16)) {
        assert(sampleData[i]==1.0);
      } else if ((i>=16) && (i<24)) {
        assert(sampleData[i]==2.0);
      } else {
        assert(sampleData[i]==3.0);
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

    // it's important that default values are different, as we need to be able to
    // verify that the tag values in each object are being added to the correct
    // default values - since the tag lists don't match, the default values will
    // be used for missing tags in each object
    myData.getDefaultValue()()=1.0;
    right.getDefaultValue()()=2.0;

    DataArray vOne(3.0);
    DataArray vTwo(4.0);
    myData.addTaggedValue(1,vOne.getView());
    right.addTaggedValue(2,vTwo.getView());

    //cout << myData.toString() << endl;
    //cout << right.toString() << endl;

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
    assert(myDataView()==5.0);

    // check result value for tag "2"
    myDataView = myData.getDataPointByTag(2);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==2);
    assert(myDataView.getRank()==0);
    assert(myDataView.noValues()==1);
    assert(myDataView.getShape().size()==0);
    assert(myDataView()==5.0);

    // check result for default value
    myDataView = myData.getDefaultValue();
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==0);
    assert(myDataView.getRank()==0);
    assert(myDataView.noValues()==1);
    assert(myDataView.getShape().size()==0);
    assert(myDataView()==3.0);

    // use a non-existent tag so we get a pointer to
    // the first element of the data array
    double* sampleData=myData.getSampleDataByTag(9);
    assert(sampleData[0]==3);
    assert(sampleData[1]==5);
    assert(sampleData[2]==5);

  }

  {
    cout << "\tTest binaryOp addition of two DataTagged objects with overlapping tag sets." << endl;

    DataTagged myData;
    DataTagged right;

    // it's important that default values are different, as we need to be able to
    // verify that the tag values in each object are being added to the correct
    // default values - since the tag lists don't match, the default values will
    // be used for missing tags in each object
    myData.getDefaultValue()()=2.0;
    right.getDefaultValue()()=3.0;

    DataArray vOne(1.0);
    myData.addTaggedValue(1,vOne.getView());
    myData.addTaggedValue(2,vOne.getView());
    right.addTaggedValue(2,vOne.getView());
    right.addTaggedValue(3,vOne.getView());

    //cout << myData.toString() << endl;
    //cout << right.toString() << endl;

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
    assert(myDataView()==4.0);

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
    assert(myDataView()==3.0);

    // check result for default value
    myDataView = myData.getDefaultValue();
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==0);
    assert(myDataView.getRank()==0);
    assert(myDataView.noValues()==1);
    assert(myDataView.getShape().size()==0);
    assert(myDataView()==5.0);

    // use a non-existent tag so we get a pointer to
    // the first element of the data array
    double* sampleData=myData.getSampleDataByTag(9);
    assert(sampleData[0]==5);
    assert(sampleData[1]==4);
    assert(sampleData[2]==2);
    assert(sampleData[3]==3);

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

    // it's important that default values are different, as we need to be able to
    // verify that the tag values in each object are being added to the correct
    // default values - since the tag lists don't match, the default values will
    // be used for missing tags in each object
    myData.getDefaultValue()()=2.0;
    right.getDefaultValue()()=3.0;

    DataArray vOne(1.0);
    DataArray vTwo(2.0);
    myData.addTaggedValue(1,vOne.getView());
    myData.addTaggedValue(2,vOne.getView());
    right.addTaggedValue(2,vTwo.getView());
    right.addTaggedValue(3,vTwo.getView());

    //cout << myData.toString() << endl;
    //cout << right.toString() << endl;

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
    assert(myDataView()==3.0);

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
    assert(myDataView()==4.0);

    // check result for default value
    myDataView = myData.getDefaultValue();
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==0);
    assert(myDataView.getRank()==0);
    assert(myDataView.noValues()==1);
    assert(myDataView.getShape().size()==0);
    assert(myDataView()==6.0);

    // use a non-existent tag so we get a pointer to
    // the first element of the data array
    double* sampleData=myData.getSampleDataByTag(9);
    assert(sampleData[0]==6);
    assert(sampleData[1]==3);
    assert(sampleData[2]==2);
    assert(sampleData[3]==4);

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
      } else if ((i>=3) && (i<6)) {
        assert(sampleData[i]==i-2);
      } else if ((i>=6) && (i<9)) {
        assert(sampleData[i]==i-4);
      } else  {
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
      } else if ((i>=3) && (i<6)) {
        assert(sampleData[i]==i-2);
      } else if ((i>=6) && (i<9)) {
        assert(sampleData[i]==i-4);
      } else if ((i>=9) && (i<12)) {
        assert(sampleData[i]==i-6);
      } else {
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
      } else if ((i>=3) && (i<6)) {
        assert(sampleData[i]==i-2);
      } else if ((i>=6) && (i<9)) {
        assert(sampleData[i]==i-4);
      } else if ((i>=9) && (i<12)) {
        assert(sampleData[i]==i-6);
      } else {
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
      } else if ((i>=3) && (i<6)) {
        assert(sampleData[i]==i-2);
      } else if ((i>=6) && (i<9)) {
        assert(sampleData[i]==i-4);
      } else if ((i>=9) && (i<12)) {
        assert(sampleData[i]==i-6);
      } else if ((i>=12) && (i<15)) {
        assert(sampleData[i]==i-8);
      } else if ((i>=15) && (i<18)) {
        assert(sampleData[i]==i-11);
      } else {
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
      } else if ((i>=3) && (i<6)) {
        assert(sampleData[i]==i-2);
      } else if ((i>=6) && (i<9)) {
        assert(sampleData[i]==i-4);
      } else if ((i>=9) && (i<12)) {
        assert(sampleData[i]==i-6);
      } else if ((i>=12) && (i<15)) {
        assert(sampleData[i]==i-8);
      } else if ((i>=15) && (i<18)) {
        assert(sampleData[i]==i-10);
      } else {
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
      } else if ((i>=3) && (i<6)) {
        assert(sampleData[i]==i-2);
      } else if ((i>=6) && (i<9)) {
        assert(sampleData[i]==i-1);
      } else {
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
    sampleData=myData.getSampleData(0);
    for (int i=0; i<myDataView.noValues(); i++) {
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
    sampleData=myData.getSampleDataByTag(0);
    for (int i=0; i<myDataView.noValues(); i++) {
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
    myDataView = myData.getDataPointByTag(9);
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
      } else {
        assert(sampleData[i]==i-2);
      }
    }
    sampleData=myData.getSampleData(0);
    for (int i=0; i<myDataView.noValues(); i++) {
      assert(sampleData[i]==i+1);
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
      } else if ((i>=3) && (i<6)) {
        assert(sampleData[i]==i-2);
      } else if ((i>=6) && (i<9)) {
        assert(sampleData[i]==i-4);
      } else {
        assert(sampleData[i]==i-6);
      }
    }
    sampleData=myData.getSampleData(0);
    for (int i=0; i<myDataView.noValues(); i++) {
      assert(sampleData[i]==i+1);
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
      } else if ((i>=3) && (i<6)) {
        assert(sampleData[i]==i-2);
      } else if ((i>=6) && (i<9)) {
        assert(sampleData[i]==i-4);
      } else {
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

void DataTaggedTestCase::testGetSlice() {

  cout << endl;

  {

    cout << "\tTest slicing default DataTagged." << endl;

    DataTagged myData;

    DataArrayView::RegionType region;

    DataAbstract* slicedDefault = myData.getSlice(region);

    // cout << slicedDefault->toString() << endl;

    const DataTagged* myDataSliced=dynamic_cast<const DataTagged*>(slicedDefault);

    assert(myDataSliced->getTagLookup().size()==0);

    assert(myDataSliced->getLength()==1);

    DataArrayView myDataView = myDataSliced->getDefaultValue();
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==0);
    assert(myDataView.getRank()==0);
    assert(myDataView.noValues()==1);
    assert(myDataView.getShape().size()==0);
    assert(myDataView()==0.0);

  }

  {

    cout << "\tTest slicing DataTagged with rank 1 default value only." << endl;

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

    // full slice

    std::pair<int, int> region_element;
    region_element.first=0;
    region_element.second=3;
    DataArrayView::RegionType region;
    region.push_back(region_element);

    DataAbstract* slicedDefault = myData.getSlice(region);

    //cout << slicedDefault->toString() << endl;

    const DataTagged* myDataSliced=dynamic_cast<const DataTagged*>(slicedDefault);

    assert(myDataSliced->getTagLookup().size()==0);

    assert(myDataSliced->getLength()==3);

    DataArrayView myDataView = myDataSliced->getDefaultValue();
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==0);
    assert(myDataView.getRank()==1);
    assert(myDataView.noValues()==3);
    assert(myDataView.getShape().size()==1);
    assert(myDataView(0)==0.0);
    assert(myDataView(1)==1.0);
    assert(myDataView(2)==2.0);

    // scalar slice

    region.clear();
    region_element.first=0;
    region_element.second=0;
    region.push_back(region_element);

    slicedDefault = myData.getSlice(region);

    //cout << slicedDefault->toString() << endl;

    myDataSliced=dynamic_cast<const DataTagged*>(slicedDefault);

    assert(myDataSliced->getTagLookup().size()==0);

    assert(myDataSliced->getLength()==1);

    myDataView = myDataSliced->getDefaultValue();
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==0);
    assert(myDataView.getRank()==0);
    assert(myDataView.noValues()==1);
    assert(myDataView.getShape().size()==0);
    assert(myDataView()==0.0);

  }

  {

    cout << "\tTest slicing DataTagged with rank 3 default value only." << endl;

    DataArrayView::ShapeType viewShape;
    viewShape.push_back(3);
    viewShape.push_back(3);
    viewShape.push_back(3);

    DataTagged::TagListType keys;

    DataTagged::ValueListType values;

    DataArrayView::ValueType viewData(27);
    for (int i=0;i<viewData.size();i++) {
      viewData[i]=i;
    }
    DataArrayView myView(viewData,viewShape);

    DataTagged myData(keys,values,myView,FunctionSpace());

    //cout << myData.toString() << endl;

    // full slice

    std::pair<int, int> region_element;
    region_element.first=0;
    region_element.second=3;
    DataArrayView::RegionType region;
    region.push_back(region_element);
    region.push_back(region_element);
    region.push_back(region_element);

    DataAbstract* slicedDefault = myData.getSlice(region);

    //cout << slicedDefault->toString() << endl;

    const DataTagged* myDataSliced=dynamic_cast<const DataTagged*>(slicedDefault);

    assert(myDataSliced->getTagLookup().size()==0);

    assert(myDataSliced->getLength()==27);

    DataArrayView myDataView = myDataSliced->getDefaultValue();
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==0);
    assert(myDataView.getRank()==3);
    assert(myDataView.noValues()==27);
    assert(myDataView.getShape().size()==3);

    // rank 1 slice

    region.clear();
    region.push_back(region_element);
    region_element.second=0;
    region.push_back(region_element);
    region.push_back(region_element);

    slicedDefault = myData.getSlice(region);

    //cout << slicedDefault->toString() << endl;

    myDataSliced=dynamic_cast<const DataTagged*>(slicedDefault);

    assert(myDataSliced->getTagLookup().size()==0);

    assert(myDataSliced->getLength()==3);

    myDataView = myDataSliced->getDefaultValue();
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==0);
    assert(myDataView.getRank()==1);
    assert(myDataView.noValues()==3);
    assert(myDataView.getShape().size()==1);
    assert(myDataView(0)==0.0);
    assert(myDataView(1)==1.0);
    assert(myDataView(2)==2.0);

    // scalar slice

    region.clear();
    region_element.first=2;
    region_element.second=2;
    region.push_back(region_element);
    region.push_back(region_element);
    region.push_back(region_element);

    slicedDefault = myData.getSlice(region);

    //cout << slicedDefault->toString() << endl;

    myDataSliced=dynamic_cast<const DataTagged*>(slicedDefault);

    assert(myDataSliced->getTagLookup().size()==0);

    assert(myDataSliced->getLength()==1);

    myDataView = myDataSliced->getDefaultValue();
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==0);
    assert(myDataView.getRank()==0);
    assert(myDataView.noValues()==1);
    assert(myDataView.getShape().size()==0);
    assert(myDataView()==26);

  }

  {

    cout << "\tTest slicing DataTagged with scalar values and one tag." << endl;

    DataTagged::TagListType keys;
    keys.push_back(1);

    DataTagged::ValueListType values;

    DataArrayView::ShapeType viewShape;

    // default value
    DataArrayView::ValueType viewData(1);
    viewData[0]=0.0;
    DataArrayView myView(viewData,viewShape);

    // value for tag "1"
    DataArray eOne(myView);
    eOne.getView()()=1.0;
    values.push_back(eOne.getView());

    DataTagged myData(keys,values,myView,FunctionSpace());

    //cout << myData.toString() << endl;

    // full slice

    DataArrayView::RegionType region;

    DataAbstract* slicedDefault = myData.getSlice(region);

    //cout << slicedDefault->toString() << endl;

    const DataTagged* myDataSliced=dynamic_cast<const DataTagged*>(slicedDefault);

    assert(myDataSliced->getTagLookup().size()==1);

    assert(myDataSliced->getLength()==2);

    DataArrayView myDataView = myDataSliced->getDefaultValue();
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==0);
    assert(myDataView.getRank()==0);
    assert(myDataView.noValues()==1);
    assert(myDataView.getShape().size()==0);
    assert(myDataView()==0);

    myDataView = myDataSliced->getDataPointByTag(1);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==1);
    assert(myDataView.getRank()==0);
    assert(myDataView.noValues()==1);
    assert(myDataView.getShape().size()==0);
    assert(myDataView()==1);

  }

  {

    cout << "\tTest slicing DataTagged with rank 1 values and one tag." << endl;

    DataArrayView::ShapeType viewShape;
    viewShape.push_back(3);

    DataTagged::TagListType keys;
    keys.push_back(1);

    DataTagged::ValueListType values;

    // default value
    DataArrayView::ValueType viewData(3);
    for (int i=0;i<viewShape[0];i++) {
      viewData[i]=i;
    }
    DataArrayView myView(viewData,viewShape);

    // value for tag "1"
    DataArray eOne(myView);
    for (int i=0;i<eOne.getView().getShape()[0];i++) {
      eOne.getView()(i)=i+3.0;
    }
    values.push_back(eOne.getView());

    DataTagged myData(keys,values,myView,FunctionSpace());

    //cout << myData.toString() << endl;

    // full slice

    std::pair<int, int> region_element;
    region_element.first=0;
    region_element.second=3;
    DataArrayView::RegionType region;
    region.push_back(region_element);

    DataAbstract* slicedDefault = myData.getSlice(region);

    //cout << slicedDefault->toString() << endl;

    const DataTagged* myDataSliced=dynamic_cast<const DataTagged*>(slicedDefault);

    assert(myDataSliced->getTagLookup().size()==1);

    assert(myDataSliced->getLength()==6);

    DataArrayView myDataView = myDataSliced->getDefaultValue();
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==0);
    assert(myDataView.getRank()==1);
    assert(myDataView.noValues()==3);
    assert(myDataView.getShape().size()==1);
    assert(myDataView(0)==0);
    assert(myDataView(1)==1);
    assert(myDataView(2)==2);

    myDataView = myDataSliced->getDataPointByTag(1);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==3);
    assert(myDataView.getRank()==1);
    assert(myDataView.noValues()==3);
    assert(myDataView.getShape().size()==1);
    assert(myDataView(0)==3);
    assert(myDataView(1)==4);
    assert(myDataView(2)==5);

    // scalar slice

    region_element.first=1;
    region_element.second=1;
    region.clear();
    region.push_back(region_element);

    slicedDefault = myData.getSlice(region);

    //cout << slicedDefault->toString() << endl;

    myDataSliced=dynamic_cast<const DataTagged*>(slicedDefault);

    assert(myDataSliced->getTagLookup().size()==1);

    assert(myDataSliced->getLength()==2);

    myDataView = myDataSliced->getDefaultValue();
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==0);
    assert(myDataView.getRank()==0);
    assert(myDataView.noValues()==1);
    assert(myDataView.getShape().size()==0);
    assert(myDataView()==1);

    myDataView = myDataSliced->getDataPointByTag(1);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==1);
    assert(myDataView.getRank()==0);
    assert(myDataView.noValues()==1);
    assert(myDataView.getShape().size()==0);
    assert(myDataView()==4);

  }

  {

    cout << "\tTest slicing DataTagged with rank 3 values and one tag." << endl;

    DataArrayView::ShapeType viewShape;
    viewShape.push_back(3);
    viewShape.push_back(3);
    viewShape.push_back(3);

    DataTagged::TagListType keys;
    keys.push_back(1);

    DataTagged::ValueListType values;

    // default value
    DataArrayView::ValueType viewData(27);
    for (int i=0;i<viewData.size();i++) {
      viewData[i]=i;
    }
    DataArrayView myView(viewData,viewShape);

    // value for tag "1"
    DataArrayView::ValueType viewData1(27);
    for (int i=0;i<viewData1.size();i++) {
      viewData1[i]=i+27.0;
    }
    DataArrayView myView1(viewData1,viewShape);
    values.push_back(myView1);

    DataTagged myData(keys,values,myView,FunctionSpace());

    //cout << myData.toString() << endl;

    // full slice

    std::pair<int, int> region_element;
    region_element.first=0;
    region_element.second=3;
    DataArrayView::RegionType region;
    region.push_back(region_element);
    region.push_back(region_element);
    region.push_back(region_element);

    DataAbstract* slicedDefault = myData.getSlice(region);

    //cout << slicedDefault->toString() << endl;

    const DataTagged* myDataSliced=dynamic_cast<const DataTagged*>(slicedDefault);

    assert(myDataSliced->getTagLookup().size()==1);

    assert(myDataSliced->getLength()==54);

    DataArrayView myDataView = myDataSliced->getDefaultValue();
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==0);
    assert(myDataView.getRank()==3);
    assert(myDataView.noValues()==27);
    assert(myDataView.getShape().size()==3);

    myDataView = myDataSliced->getDataPointByTag(1);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==27);
    assert(myDataView.getRank()==3);
    assert(myDataView.noValues()==27);
    assert(myDataView.getShape().size()==3);

    // rank 1 slice

    region.clear();
    region.push_back(region_element);
    region_element.second=0;
    region.push_back(region_element);
    region.push_back(region_element);

    slicedDefault = myData.getSlice(region);

    //cout << slicedDefault->toString() << endl;

    myDataSliced=dynamic_cast<const DataTagged*>(slicedDefault);

    assert(myDataSliced->getTagLookup().size()==1);

    assert(myDataSliced->getLength()==6);

    myDataView = myDataSliced->getDefaultValue();
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==0);
    assert(myDataView.getRank()==1);
    assert(myDataView.noValues()==3);
    assert(myDataView.getShape().size()==1);
    assert(myDataView(0)==0);
    assert(myDataView(1)==1);
    assert(myDataView(2)==2);

    myDataView = myDataSliced->getDataPointByTag(1);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==3);
    assert(myDataView.getRank()==1);
    assert(myDataView.noValues()==3);
    assert(myDataView.getShape().size()==1);
    assert(myDataView(0)==27);
    assert(myDataView(1)==28);
    assert(myDataView(2)==29);

    // scalar slice

    region_element.first=1;
    region_element.second=1;
    region.clear();
    region.push_back(region_element);
    region.push_back(region_element);
    region.push_back(region_element);

    slicedDefault = myData.getSlice(region);

    //cout << slicedDefault->toString() << endl;

    myDataSliced=dynamic_cast<const DataTagged*>(slicedDefault);

    assert(myDataSliced->getTagLookup().size()==1);

    assert(myDataSliced->getLength()==2);

    myDataView = myDataSliced->getDefaultValue();
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==0);
    assert(myDataView.getRank()==0);
    assert(myDataView.noValues()==1);
    assert(myDataView.getShape().size()==0);
    assert(myDataView()==13);

    myDataView = myDataSliced->getDataPointByTag(1);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==1);
    assert(myDataView.getRank()==0);
    assert(myDataView.noValues()==1);
    assert(myDataView.getShape().size()==0);
    assert(myDataView()==40);

  }

  {

    cout << "\tTest slicing DataTagged with scalar values and three tags." << endl;

    DataTagged::TagListType keys;
    keys.push_back(1);
    keys.push_back(2);
    keys.push_back(3);

    DataTagged::ValueListType values;

    DataArrayView::ShapeType viewShape;

    // default value
    DataArrayView::ValueType viewData(1);
    viewData[0]=0.0;
    DataArrayView myView(viewData,viewShape);

    // value for tag "1"
    DataArray eOne(myView);
    eOne.getView()()=1.0;
    values.push_back(eOne.getView());

    // value for tag "2"
    DataArray eTwo(myView);
    eTwo.getView()()=2.0;
    values.push_back(eTwo.getView());

    // value for tag "3"
    DataArray eThree(myView);
    eThree.getView()()=3.0;
    values.push_back(eThree.getView());

    DataTagged myData(keys,values,myView,FunctionSpace());

    // cout << myData.toString() << endl;

    // full slice

    DataArrayView::RegionType region;

    DataAbstract* slicedDefault = myData.getSlice(region);

    //cout << slicedDefault->toString() << endl;

    const DataTagged* myDataSliced=dynamic_cast<const DataTagged*>(slicedDefault);

    assert(myDataSliced->getTagLookup().size()==3);

    assert(myDataSliced->getLength()==4);

    DataArrayView myDataView = myDataSliced->getDefaultValue();
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==0);
    assert(myDataView.getRank()==0);
    assert(myDataView.noValues()==1);
    assert(myDataView.getShape().size()==0);
    assert(myDataView()==0);

    myDataView = myDataSliced->getDataPointByTag(1);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==1);
    assert(myDataView.getRank()==0);
    assert(myDataView.noValues()==1);
    assert(myDataView.getShape().size()==0);
    assert(myDataView()==1);

    myDataView = myDataSliced->getDataPointByTag(2);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==2);
    assert(myDataView.getRank()==0);
    assert(myDataView.noValues()==1);
    assert(myDataView.getShape().size()==0);
    assert(myDataView()==2);

    myDataView = myDataSliced->getDataPointByTag(3);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==3);
    assert(myDataView.getRank()==0);
    assert(myDataView.noValues()==1);
    assert(myDataView.getShape().size()==0);
    assert(myDataView()==3);

  }

  {

    cout << "\tTest slicing DataTagged with rank 1 values and three tags." << endl;

    DataArrayView::ShapeType viewShape;
    viewShape.push_back(3);

    DataTagged::TagListType keys;
    keys.push_back(1);
    keys.push_back(2);
    keys.push_back(3);

    DataTagged::ValueListType values;

    // default value
    DataArrayView::ValueType viewData(3);
    for (int i=0;i<viewShape[0];i++) {
      viewData[i]=i;
    }
    DataArrayView myView(viewData,viewShape);

    // value for tag "1"
    DataArray eOne(myView);
    for (int i=0;i<eOne.getView().getShape()[0];i++) {
      eOne.getView()(i)=i+3.0;
    }
    values.push_back(eOne.getView());

    // value for tag "2"
    DataArray eTwo(myView);
    for (int i=0;i<eTwo.getView().getShape()[0];i++) {
      eTwo.getView()(i)=i+6.0;
    }
    values.push_back(eTwo.getView());

    // value for tag "3"
    DataArray eThree(myView);
    for (int i=0;i<eThree.getView().getShape()[0];i++) {
      eThree.getView()(i)=i+9.0;
    }
    values.push_back(eThree.getView());

    DataTagged myData(keys,values,myView,FunctionSpace());

    //cout << myData.toString() << endl;

    // full slice

    std::pair<int, int> region_element;
    region_element.first=0;
    region_element.second=3;
    DataArrayView::RegionType region;
    region.push_back(region_element);

    DataAbstract* slicedDefault = myData.getSlice(region);

    //cout << slicedDefault->toString() << endl;

    const DataTagged* myDataSliced=dynamic_cast<const DataTagged*>(slicedDefault);

    assert(myDataSliced->getTagLookup().size()==3);

    assert(myDataSliced->getLength()==12);

    DataArrayView myDataView = myDataSliced->getDefaultValue();
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==0);
    assert(myDataView.getRank()==1);
    assert(myDataView.noValues()==3);
    assert(myDataView.getShape().size()==1);
    assert(myDataView(0)==0);
    assert(myDataView(1)==1);
    assert(myDataView(2)==2);

    myDataView = myDataSliced->getDataPointByTag(1);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==3);
    assert(myDataView.getRank()==1);
    assert(myDataView.noValues()==3);
    assert(myDataView.getShape().size()==1);
    assert(myDataView(0)==3);
    assert(myDataView(1)==4);
    assert(myDataView(2)==5);

    myDataView = myDataSliced->getDataPointByTag(2);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==6);
    assert(myDataView.getRank()==1);
    assert(myDataView.noValues()==3);
    assert(myDataView.getShape().size()==1);
    assert(myDataView(0)==6);
    assert(myDataView(1)==7);
    assert(myDataView(2)==8);

    myDataView = myDataSliced->getDataPointByTag(3);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==9);
    assert(myDataView.getRank()==1);
    assert(myDataView.noValues()==3);
    assert(myDataView.getShape().size()==1);
    assert(myDataView(0)==9);
    assert(myDataView(1)==10);
    assert(myDataView(2)==11);

    // scalar slice

    region.clear();
    region_element.first=1;
    region_element.second=1;
    region.push_back(region_element);

    slicedDefault = myData.getSlice(region);

    //cout << slicedDefault->toString() << endl;

    myDataSliced=dynamic_cast<const DataTagged*>(slicedDefault);

    assert(myDataSliced->getTagLookup().size()==3);

    assert(myDataSliced->getLength()==4);

    myDataView = myDataSliced->getDefaultValue();
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==0);
    assert(myDataView.getRank()==0);
    assert(myDataView.noValues()==1);
    assert(myDataView.getShape().size()==0);
    assert(myDataView()==1);

    myDataView = myDataSliced->getDataPointByTag(1);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==1);
    assert(myDataView.getRank()==0);
    assert(myDataView.noValues()==1);
    assert(myDataView.getShape().size()==0);
    assert(myDataView()==4);

    myDataView = myDataSliced->getDataPointByTag(2);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==2);
    assert(myDataView.getRank()==0);
    assert(myDataView.noValues()==1);
    assert(myDataView.getShape().size()==0);
    assert(myDataView()==7);

    myDataView = myDataSliced->getDataPointByTag(3);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==3);
    assert(myDataView.getRank()==0);
    assert(myDataView.noValues()==1);
    assert(myDataView.getShape().size()==0);
    assert(myDataView()==10);

  }

  {

    cout << "\tTest slicing DataTagged with rank 3 values and three tags." << endl;

    DataArrayView::ShapeType viewShape;
    viewShape.push_back(3);
    viewShape.push_back(3);
    viewShape.push_back(3);

    DataTagged::TagListType keys;
    keys.push_back(1);
    keys.push_back(2);
    keys.push_back(3);

    DataTagged::ValueListType values;

    // default value
    DataArrayView::ValueType viewData(27);
    for (int i=0;i<viewData.size();i++) {
      viewData[i]=i;
    }
    DataArrayView myView(viewData,viewShape);

    // value for tag "1"
    DataArrayView::ValueType viewData1(27);
    for (int i=0;i<viewData1.size();i++) {
      viewData1[i]=i+27.0;
    }
    DataArrayView myView1(viewData1,viewShape);
    values.push_back(myView1);

    // value for tag "2"
    DataArrayView::ValueType viewData2(27);
    for (int i=0;i<viewData2.size();i++) {
      viewData2[i]=i+54.0;
    }
    DataArrayView myView2(viewData2,viewShape);
    values.push_back(myView2);

    // value for tag "3"
    DataArrayView::ValueType viewData3(27);
    for (int i=0;i<viewData3.size();i++) {
      viewData3[i]=i+81.0;
    }
    DataArrayView myView3(viewData3,viewShape);
    values.push_back(myView3);

    DataTagged myData(keys,values,myView,FunctionSpace());

    //cout << myData.toString() << endl;

    // full slice

    std::pair<int, int> region_element;
    region_element.first=0;
    region_element.second=3;
    DataArrayView::RegionType region;
    region.push_back(region_element);
    region.push_back(region_element);
    region.push_back(region_element);

    DataAbstract* slicedDefault = myData.getSlice(region);

    //cout << slicedDefault->toString() << endl;

    const DataTagged* myDataSliced=dynamic_cast<const DataTagged*>(slicedDefault);

    assert(myDataSliced->getTagLookup().size()==3);

    assert(myDataSliced->getLength()==108);

    DataArrayView myDataView = myDataSliced->getDefaultValue();
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==0);
    assert(myDataView.getRank()==3);
    assert(myDataView.noValues()==27);
    assert(myDataView.getShape().size()==3);

    myDataView = myDataSliced->getDataPointByTag(1);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==27);
    assert(myDataView.getRank()==3);
    assert(myDataView.noValues()==27);
    assert(myDataView.getShape().size()==3);

    myDataView = myDataSliced->getDataPointByTag(2);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==54);
    assert(myDataView.getRank()==3);
    assert(myDataView.noValues()==27);
    assert(myDataView.getShape().size()==3);

    myDataView = myDataSliced->getDataPointByTag(3);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==81);
    assert(myDataView.getRank()==3);
    assert(myDataView.noValues()==27);
    assert(myDataView.getShape().size()==3);

    // rank 1 slice

    region.clear();
    region.push_back(region_element);
    region_element.second=0;
    region.push_back(region_element);
    region.push_back(region_element);

    slicedDefault = myData.getSlice(region);

    // cout << slicedDefault->toString() << endl;

    myDataSliced=dynamic_cast<const DataTagged*>(slicedDefault);

    assert(myDataSliced->getTagLookup().size()==3);

    assert(myDataSliced->getLength()==12);

    myDataView = myDataSliced->getDefaultValue();
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==0);
    assert(myDataView.getRank()==1);
    assert(myDataView.noValues()==3);
    assert(myDataView.getShape().size()==1);
    assert(myDataView(0)==0);
    assert(myDataView(1)==1);
    assert(myDataView(2)==2);

    myDataView = myDataSliced->getDataPointByTag(1);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==3);
    assert(myDataView.getRank()==1);
    assert(myDataView.noValues()==3);
    assert(myDataView.getShape().size()==1);
    assert(myDataView(0)==27);
    assert(myDataView(1)==28);
    assert(myDataView(2)==29);

    myDataView = myDataSliced->getDataPointByTag(2);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==6);
    assert(myDataView.getRank()==1);
    assert(myDataView.noValues()==3);
    assert(myDataView.getShape().size()==1);
    assert(myDataView(0)==54);
    assert(myDataView(1)==55);
    assert(myDataView(2)==56);

    myDataView = myDataSliced->getDataPointByTag(3);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==9);
    assert(myDataView.getRank()==1);
    assert(myDataView.noValues()==3);
    assert(myDataView.getShape().size()==1);
    assert(myDataView(0)==81);
    assert(myDataView(1)==82);
    assert(myDataView(2)==83);

    // scalar slice

    region_element.first=1;
    region_element.second=1;
    region.clear();
    region.push_back(region_element);
    region.push_back(region_element);
    region.push_back(region_element);

    slicedDefault = myData.getSlice(region);

    //cout << slicedDefault->toString() << endl;

    myDataSliced=dynamic_cast<const DataTagged*>(slicedDefault);

    assert(myDataSliced->getTagLookup().size()==3);

    assert(myDataSliced->getLength()==4);

    myDataView = myDataSliced->getDefaultValue();
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==0);
    assert(myDataView.getRank()==0);
    assert(myDataView.noValues()==1);
    assert(myDataView.getShape().size()==0);
    assert(myDataView()==13);

    myDataView = myDataSliced->getDataPointByTag(1);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==1);
    assert(myDataView.getRank()==0);
    assert(myDataView.noValues()==1);
    assert(myDataView.getShape().size()==0);
    assert(myDataView()==40);

    myDataView = myDataSliced->getDataPointByTag(2);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==2);
    assert(myDataView.getRank()==0);
    assert(myDataView.noValues()==1);
    assert(myDataView.getShape().size()==0);
    assert(myDataView()==67);

    myDataView = myDataSliced->getDataPointByTag(3);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==3);
    assert(myDataView.getRank()==0);
    assert(myDataView.noValues()==1);
    assert(myDataView.getShape().size()==0);
    assert(myDataView()==94);

  }

}

void DataTaggedTestCase::testSetSlice() {

  cout << endl;

  {

    cout << "\tTest slicing default DataTagged." << endl;

    DataTagged myData1;
    DataTagged myData2;

    DataArrayView::RegionType region;

    myData2.getDefaultValue()()=1.0;

    myData1.setSlice(&myData2, region);

    //cout << myData1.toString() << endl;

    assert(myData1.getTagLookup().size()==0);

    assert(myData1.getLength()==1);

    DataArrayView myDataView = myData1.getDefaultValue();
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==0);
    assert(myDataView.getRank()==0);
    assert(myDataView.noValues()==1);
    assert(myDataView.getShape().size()==0);
    assert(myDataView()==1.0);

  }

  {

    cout << "\tTest slicing DataTagged with rank 1 default value only." << endl;

    DataTagged::TagListType keys;

    DataTagged::ValueListType values;

    DataArrayView::ShapeType viewShape;
    viewShape.push_back(3);

    DataArrayView::ValueType viewData1(3);
    for (int i=0;i<viewShape[0];i++) {
      viewData1[i]=i;
    }
    DataArrayView myView1(viewData1,viewShape);
    DataTagged myData1(keys,values,myView1,FunctionSpace());

    DataArrayView::ValueType viewData2(3);
    for (int i=0;i<viewShape[0];i++) {
      viewData2[i]=i+3;
    }
    DataArrayView myView2(viewData2,viewShape);
    DataTagged myData2(keys,values,myView2,FunctionSpace());

    // full slice

    std::pair<int, int> region_element;
    region_element.first=0;
    region_element.second=3;
    DataArrayView::RegionType region;
    region.push_back(region_element);

    myData1.setSlice(&myData2, region);

    //cout << myData1.toString() << endl;

    assert(myData1.getTagLookup().size()==0);

    assert(myData1.getLength()==3);

    DataArrayView myDataView = myData1.getDefaultValue();
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==0);
    assert(myDataView.getRank()==1);
    assert(myDataView.noValues()==3);
    assert(myDataView.getShape().size()==1);
    assert(myDataView(0)==3.0);
    assert(myDataView(1)==4.0);
    assert(myDataView(2)==5.0);

    // rank 1 slice

    viewShape.clear();
    viewShape.push_back(1);

    DataArrayView::ValueType viewData3(1);
    viewData3[0]=6.0;
    DataArrayView myView3(viewData3,viewShape);
    DataTagged myData3(keys,values,myView3,FunctionSpace());

    region.clear();
    region_element.first=1;
    region_element.second=2;
    region.push_back(region_element);

    myData1.setSlice(&myData3, region);

    //cout << myData1.toString() << endl;

    assert(myData1.getTagLookup().size()==0);

    assert(myData1.getLength()==3);

    myDataView = myData1.getDefaultValue();
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==0);
    assert(myDataView.getRank()==1);
    assert(myDataView.noValues()==3);
    assert(myDataView.getShape().size()==1);
    assert(myDataView(0)==3.0);
    assert(myDataView(1)==6.0);
    assert(myDataView(2)==5.0);

    // scalar slice

    region.clear();
    region_element.first=0;
    region_element.second=0;
    region.push_back(region_element);

    DataTagged myData4;
    myData4.getDefaultValue()()=7.0;

    myData1.setSlice(&myData4, region);

    //cout << myData3.toString() << endl;

    assert(myData1.getTagLookup().size()==0);

    assert(myData1.getLength()==3);

    myDataView = myData1.getDefaultValue();
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==0);
    assert(myDataView.getRank()==1);
    assert(myDataView.noValues()==3);
    assert(myDataView.getShape().size()==1);
    assert(myDataView(0)==7.0);
    assert(myDataView(1)==6.0);
    assert(myDataView(2)==5.0);

  }

  {

    cout << "\tTest slicing DataTagged with rank 3 default value only." << endl;

    DataTagged::TagListType keys;

    DataTagged::ValueListType values;

    DataArrayView::ShapeType viewShape;
    viewShape.push_back(3);
    viewShape.push_back(3);
    viewShape.push_back(3);

    DataArrayView::ValueType viewData1(27);
    for (int i=0;i<viewData1.size();i++) {
      viewData1[i]=i;
    }
    DataArrayView myView1(viewData1,viewShape);
    DataTagged myData1(keys,values,myView1,FunctionSpace());

    DataArrayView::ValueType viewData2(27);
    for (int i=0;i<viewData2.size();i++) {
      viewData2[i]=i+27;
    }
    DataArrayView myView2(viewData2,viewShape);
    DataTagged myData2(keys,values,myView2,FunctionSpace());

    // full slice

    std::pair<int, int> region_element;
    region_element.first=0;
    region_element.second=3;
    DataArrayView::RegionType region;
    region.push_back(region_element);
    region.push_back(region_element);
    region.push_back(region_element);

    myData1.setSlice(&myData2, region);

    //cout << myData1.toString() << endl;

    assert(myData1.getTagLookup().size()==0);

    assert(myData1.getLength()==27);

    DataArrayView myDataView = myData1.getDefaultValue();
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==0);
    assert(myDataView.getRank()==3);
    assert(myDataView.noValues()==27);
    assert(myDataView.getShape().size()==3);

    // rank 1 slice

    viewShape.clear();
    viewShape.push_back(3);

    DataArrayView::ValueType viewData3(3);
    for (int i=0;i<viewData3.size();i++) {
      viewData3[i]=i+60;
    }
    DataArrayView myView3(viewData3,viewShape);
    DataTagged myData3(keys,values,myView3,FunctionSpace());

    region.clear();
    region.push_back(region_element);
    region_element.first=0;
    region_element.second=0;
    region.push_back(region_element);
    region.push_back(region_element);

    myData1.setSlice(&myData3, region);

    //cout << myData1.toString() << endl;

    assert(myData1.getTagLookup().size()==0);

    assert(myData1.getLength()==27);

    myDataView = myData1.getDefaultValue();
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==0);
    assert(myDataView.getRank()==3);
    assert(myDataView.noValues()==27);
    assert(myDataView.getShape().size()==3);
    assert(myDataView(0,0,0)==60.0);
    assert(myDataView(1,0,0)==61.0);
    assert(myDataView(2,0,0)==62.0);

    // scalar slice

    region.clear();
    region_element.first=0;
    region_element.second=0;
    region.push_back(region_element);
    region.push_back(region_element);
    region.push_back(region_element);

    DataTagged myData4;
    myData4.getDefaultValue()()=70.0;

    myData1.setSlice(&myData4, region);

    //cout << myData1.toString() << endl;

    assert(myData1.getTagLookup().size()==0);

    assert(myData1.getLength()==27);

    myDataView = myData1.getDefaultValue();
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==0);
    assert(myDataView.getRank()==3);
    assert(myDataView.noValues()==27);
    assert(myDataView.getShape().size()==3);
    assert(myDataView(0,0,0)==70.0);

  }

  {

    cout << "\tTest slicing DataTagged with scalar values and one tag." << endl;

    DataTagged::TagListType keys;
    keys.push_back(1);

    DataTagged::ValueListType values;

    DataArrayView::ShapeType viewShape;

    // default value for Data1
    DataArrayView::ValueType viewData1(1);
    viewData1[0]=0.0;
    DataArrayView myView1(viewData1,viewShape);

    // value for tag "1" for Data1
    DataArrayView::ValueType viewData2(1);
    viewData2[0]=0.0;
    DataArrayView myView2(viewData2,viewShape);
    values.push_back(myView2);

    DataTagged myData1(keys,values,myView1,FunctionSpace());

    values.clear();

    // default value for Data2
    DataArrayView::ValueType viewData3(1);
    viewData3[0]=1.0;
    DataArrayView myView3(viewData3,viewShape);

    // value for tag "1" for Data2
    DataArrayView::ValueType viewData4(1);
    viewData4[0]=2.0;
    DataArrayView myView4(viewData4,viewShape);
    values.push_back(myView4);

    DataTagged myData2(keys,values,myView3,FunctionSpace());

    // full slice

    DataArrayView::RegionType region;

    myData1.setSlice(&myData2, region);

    //cout << myData1.toString() << endl;

    assert(myData1.getTagLookup().size()==1);

    assert(myData1.getLength()==2);

    DataArrayView myDataView = myData1.getDefaultValue();
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==0);
    assert(myDataView.getRank()==0);
    assert(myDataView.noValues()==1);
    assert(myDataView.getShape().size()==0);
    assert(myDataView()==1.0);

    myDataView = myData1.getDataPointByTag(1);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==1);
    assert(myDataView.getRank()==0);
    assert(myDataView.noValues()==1);
    assert(myDataView.getShape().size()==0);
    assert(myDataView()==2.0);

  }

  {

    cout << "\tTest slicing DataTagged with rank 1 values and one tag." << endl;

    DataTagged::TagListType keys;
    keys.push_back(1);

    DataTagged::ValueListType values;

    DataArrayView::ShapeType viewShape;
    viewShape.push_back(3);

    // default value for Data1
    DataArrayView::ValueType viewData1(3);
    for (int i=0;i<viewData1.size();i++) {
      viewData1[i]=0.0;
    }
    DataArrayView myView1(viewData1,viewShape);

    // value for tag "1" for Data1
    DataArrayView::ValueType viewData2(3);
    for (int i=0;i<viewData2.size();i++) {
      viewData2[i]=0.0;
    }
    DataArrayView myView2(viewData2,viewShape);
    values.push_back(myView2);

    DataTagged myData1(keys,values,myView1,FunctionSpace());

    values.clear();

    // default value for Data2
    DataArrayView::ValueType viewData3(3);
    for (int i=0;i<viewData3.size();i++) {
      viewData3[i]=1.0;
    }
    DataArrayView myView3(viewData3,viewShape);

    // value for tag "1" for Data2
    DataArrayView::ValueType viewData4(3);
    for (int i=0;i<viewData4.size();i++) {
      viewData4[i]=2.0;
    }
    DataArrayView myView4(viewData4,viewShape);
    values.push_back(myView4);

    DataTagged myData2(keys,values,myView3,FunctionSpace());

    // full slice

    std::pair<int, int> region_element;
    region_element.first=0;
    region_element.second=3;
    DataArrayView::RegionType region;
    region.push_back(region_element);

    myData1.setSlice(&myData2, region);

    //cout << myData1.toString() << endl;

    assert(myData1.getTagLookup().size()==1);

    assert(myData1.getLength()==6);

    DataArrayView myDataView = myData1.getDefaultValue();
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==0);
    assert(myDataView.getRank()==1);
    assert(myDataView.noValues()==3);
    assert(myDataView.getShape().size()==1);
    assert(myDataView(0)==1.0);
    assert(myDataView(1)==1.0);
    assert(myDataView(2)==1.0);

    myDataView = myData1.getDataPointByTag(1);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==3);
    assert(myDataView.getRank()==1);
    assert(myDataView.noValues()==3);
    assert(myDataView.getShape().size()==1);
    assert(myDataView(0)==2.0);
    assert(myDataView(1)==2.0);
    assert(myDataView(2)==2.0);

    // rank 1 slice

    viewShape.clear();
    viewShape.push_back(1);

    DataArrayView::ValueType viewData5(1);
    viewData5[0]=3.0;
    DataArrayView myView5(viewData5,viewShape);

    values.clear();

    DataArrayView::ValueType viewData6(1);
    viewData6[0]=4.0;
    DataArrayView myView6(viewData6,viewShape);
    values.push_back(myView6);

    DataTagged myData3(keys,values,myView5,FunctionSpace());

    region.clear();
    region_element.first=1;
    region_element.second=2;
    region.push_back(region_element);

    myData1.setSlice(&myData3, region);

    //cout << myData1.toString() << endl;

    assert(myData1.getTagLookup().size()==1);

    assert(myData1.getLength()==6);

    myDataView = myData1.getDefaultValue();
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==0);
    assert(myDataView.getRank()==1);
    assert(myDataView.noValues()==3);
    assert(myDataView.getShape().size()==1);
    assert(myDataView(0)==1.0);
    assert(myDataView(1)==3.0);
    assert(myDataView(2)==1.0);

    myDataView = myData1.getDataPointByTag(1);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==3);
    assert(myDataView.getRank()==1);
    assert(myDataView.noValues()==3);
    assert(myDataView.getShape().size()==1);
    assert(myDataView(0)==2.0);
    assert(myDataView(1)==4.0);
    assert(myDataView(2)==2.0);

    // scalar slice

    viewShape.clear();

    DataArrayView::ValueType viewData7(1);
    viewData7[0]=5.0;
    DataArrayView myView7(viewData7,viewShape);

    values.clear();

    DataArrayView::ValueType viewData8(1);
    viewData8[0]=6.0;
    DataArrayView myView8(viewData8,viewShape);
    values.push_back(myView8);

    DataTagged myData4(keys,values,myView7,FunctionSpace());

    region.clear();
    region_element.first=0;
    region_element.second=0;
    region.push_back(region_element);

    myData1.setSlice(&myData4, region);

    //cout << myData1.toString() << endl;

    myDataView = myData1.getDefaultValue();
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==0);
    assert(myDataView.getRank()==1);
    assert(myDataView.noValues()==3);
    assert(myDataView.getShape().size()==1);
    assert(myDataView(0)==5.0);
    assert(myDataView(1)==3.0);
    assert(myDataView(2)==1.0);

    myDataView = myData1.getDataPointByTag(1);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==3);
    assert(myDataView.getRank()==1);
    assert(myDataView.noValues()==3);
    assert(myDataView.getShape().size()==1);
    assert(myDataView(0)==6.0);
    assert(myDataView(1)==4.0);
    assert(myDataView(2)==2.0);

  }

  {

    cout << "\tTest slicing DataTagged with rank 3 values and one tag." << endl;

    DataTagged::TagListType keys;
    keys.push_back(1);

    DataTagged::ValueListType values;

    DataArrayView::ShapeType viewShape;
    viewShape.push_back(3);
    viewShape.push_back(3);
    viewShape.push_back(3);

    // default value for Data1
    DataArrayView::ValueType viewData1(27);
    for (int i=0;i<viewData1.size();i++) {
      viewData1[i]=0.0;
    }
    DataArrayView myView1(viewData1,viewShape);

    // value for tag "1" for Data1
    DataArrayView::ValueType viewData2(27);
    for (int i=0;i<viewData2.size();i++) {
      viewData2[i]=0.0;
    }
    DataArrayView myView2(viewData2,viewShape);
    values.push_back(myView2);

    DataTagged myData1(keys,values,myView1,FunctionSpace());

    values.clear();

    // default value for Data2
    DataArrayView::ValueType viewData3(27);
    for (int i=0;i<viewData3.size();i++) {
      viewData3[i]=1.0;
    }
    DataArrayView myView3(viewData3,viewShape);

    // value for tag "1" for Data2
    DataArrayView::ValueType viewData4(27);
    for (int i=0;i<viewData4.size();i++) {
      viewData4[i]=2.0;
    }
    DataArrayView myView4(viewData4,viewShape);
    values.push_back(myView4);

    DataTagged myData2(keys,values,myView3,FunctionSpace());

    // full slice

    std::pair<int, int> region_element;
    region_element.first=0;
    region_element.second=3;
    DataArrayView::RegionType region;
    region.push_back(region_element);
    region.push_back(region_element);
    region.push_back(region_element);

    myData1.setSlice(&myData2, region);

    //cout << myData1.toString() << endl;

    assert(myData1.getTagLookup().size()==1);

    assert(myData1.getLength()==54);

    DataArrayView myDataView = myData1.getDefaultValue();
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==0);
    assert(myDataView.getRank()==3);
    assert(myDataView.noValues()==27);
    assert(myDataView.getShape().size()==3);
    assert(myDataView(0,0,0)==1.0);
    assert(myDataView(1,1,1)==1.0);
    assert(myDataView(2,2,2)==1.0);

    myDataView = myData1.getDataPointByTag(1);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==27);
    assert(myDataView.getRank()==3);
    assert(myDataView.noValues()==27);
    assert(myDataView.getShape().size()==3);
    assert(myDataView(0,0,0)==2.0);
    assert(myDataView(1,1,1)==2.0);
    assert(myDataView(2,2,2)==2.0);

    // rank 1 slice

    viewShape.clear();
    viewShape.push_back(3);

    DataArrayView::ValueType viewData5(3);
    for (int i=0;i<viewData5.size();i++) {
      viewData5[i]=3.0;
    }
    DataArrayView myView5(viewData5,viewShape);

    values.clear();

    DataArrayView::ValueType viewData6(3);
    for (int i=0;i<viewData6.size();i++) {
      viewData6[i]=4.0;
    }
    DataArrayView myView6(viewData6,viewShape);
    values.push_back(myView6);

    DataTagged myData3(keys,values,myView5,FunctionSpace());

    region.clear();
    region.push_back(region_element);
    region_element.first=0;
    region_element.second=0;
    region.push_back(region_element);
    region.push_back(region_element);

    myData1.setSlice(&myData3, region);

    //cout << myData1.toString() << endl;

    assert(myData1.getTagLookup().size()==1);

    assert(myData1.getLength()==54);

    myDataView = myData1.getDefaultValue();
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==0);
    assert(myDataView.getRank()==3);
    assert(myDataView.noValues()==27);
    assert(myDataView.getShape().size()==3);
    assert(myDataView(0,0,0)==3.0);
    assert(myDataView(1,0,0)==3.0);
    assert(myDataView(2,0,0)==3.0);

    myDataView = myData1.getDataPointByTag(1);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==27);
    assert(myDataView.getRank()==3);
    assert(myDataView.noValues()==27);
    assert(myDataView.getShape().size()==3);
    assert(myDataView(0,0,0)==4.0);
    assert(myDataView(1,0,0)==4.0);
    assert(myDataView(2,0,0)==4.0);

    // scalar slice

    viewShape.clear();

    DataArrayView::ValueType viewData7(1);
    viewData7[0]=5.0;
    DataArrayView myView7(viewData7,viewShape);

    values.clear();

    DataArrayView::ValueType viewData8(1);
    viewData8[0]=6.0;
    DataArrayView myView8(viewData8,viewShape);
    values.push_back(myView8);

    DataTagged myData4(keys,values,myView7,FunctionSpace());

    region.clear();
    region_element.first=0;
    region_element.second=0;
    region.push_back(region_element);
    region.push_back(region_element);
    region.push_back(region_element);

    myData1.setSlice(&myData4, region);

    //cout << myData1.toString() << endl;

    myDataView = myData1.getDefaultValue();
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==0);
    assert(myDataView.getRank()==3);
    assert(myDataView.noValues()==27);
    assert(myDataView.getShape().size()==3);
    assert(myDataView(0,0,0)==5.0);

    myDataView = myData1.getDataPointByTag(1);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==27);
    assert(myDataView.getRank()==3);
    assert(myDataView.noValues()==27);
    assert(myDataView.getShape().size()==3);
    assert(myDataView(0,0,0)==6.0);

  }

  {

    cout << "\tTest slicing DataTagged with scalar values and three tags." << endl;

    DataTagged::TagListType keys;
    keys.push_back(1);
    keys.push_back(2);
    keys.push_back(3);

    DataTagged::ValueListType values;

    DataArrayView::ShapeType viewShape;

    // default value for Data1
    DataArrayView::ValueType viewData1(1);
    viewData1[0]=0.0;
    DataArrayView myView1(viewData1,viewShape);

    // value for tag "1" for Data1
    DataArrayView::ValueType viewData2(1);
    viewData2[0]=0.0;
    DataArrayView myView2(viewData2,viewShape);
    values.push_back(myView2);

    // value for tag "2" for Data1
    DataArrayView::ValueType viewData5(1);
    viewData5[0]=0.0;
    DataArrayView myView5(viewData5,viewShape);
    values.push_back(myView5);

    // value for tag "3" for Data1
    DataArrayView::ValueType viewData6(1);
    viewData6[0]=0.0;
    DataArrayView myView6(viewData6,viewShape);
    values.push_back(myView6);

    DataTagged myData1(keys,values,myView1,FunctionSpace());

    values.clear();

    // default value for Data2
    DataArrayView::ValueType viewData3(1);
    viewData3[0]=1.0;
    DataArrayView myView3(viewData3,viewShape);

    // value for tag "1" for Data2
    DataArrayView::ValueType viewData4(1);
    viewData4[0]=2.0;
    DataArrayView myView4(viewData4,viewShape);
    values.push_back(myView4);

    // value for tag "2" for Data2
    DataArrayView::ValueType viewData7(1);
    viewData7[0]=3.0;
    DataArrayView myView7(viewData7,viewShape);
    values.push_back(myView7);

    // value for tag "3" for Data2
    DataArrayView::ValueType viewData8(1);
    viewData8[0]=4.0;
    DataArrayView myView8(viewData8,viewShape);
    values.push_back(myView8);

    DataTagged myData2(keys,values,myView3,FunctionSpace());

    // full slice

    DataArrayView::RegionType region;

    myData1.setSlice(&myData2, region);

    //cout << myData1.toString() << endl;

    assert(myData1.getTagLookup().size()==3);

    assert(myData1.getLength()==4);

    DataArrayView myDataView = myData1.getDefaultValue();
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==0);
    assert(myDataView.getRank()==0);
    assert(myDataView.noValues()==1);
    assert(myDataView.getShape().size()==0);
    assert(myDataView()==1.0);

    myDataView = myData1.getDataPointByTag(1);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==1);
    assert(myDataView.getRank()==0);
    assert(myDataView.noValues()==1);
    assert(myDataView.getShape().size()==0);
    assert(myDataView()==2.0);

    myDataView = myData1.getDataPointByTag(2);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==2);
    assert(myDataView.getRank()==0);
    assert(myDataView.noValues()==1);
    assert(myDataView.getShape().size()==0);
    assert(myDataView()==3.0);

    myDataView = myData1.getDataPointByTag(3);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==3);
    assert(myDataView.getRank()==0);
    assert(myDataView.noValues()==1);
    assert(myDataView.getShape().size()==0);
    assert(myDataView()==4.0);

  }

  {

    cout << "\tTest slicing DataTagged with rank 1 values and three tags." << endl;

    DataTagged::TagListType keys;
    keys.push_back(1);
    keys.push_back(2);
    keys.push_back(3);

    DataTagged::ValueListType values;

    DataArrayView::ShapeType viewShape;
    viewShape.push_back(3);

    // default value for Data1
    DataArrayView::ValueType viewData1(3);
    for (int i=0;i<viewData1.size();i++) {
      viewData1[i]=0.0;
    }
    DataArrayView myView1(viewData1,viewShape);

    // value for tag "1" for Data1
    DataArrayView::ValueType viewData2(3);
    for (int i=0;i<viewData2.size();i++) {
      viewData2[i]=0.0;
    }
    DataArrayView myView2(viewData2,viewShape);
    values.push_back(myView2);

    // value for tag "2" for Data1
    DataArrayView::ValueType viewData3(3);
    for (int i=0;i<viewData3.size();i++) {
      viewData3[i]=0.0;
    }
    DataArrayView myView3(viewData3,viewShape);
    values.push_back(myView3);

    // value for tag "3" for Data1
    DataArrayView::ValueType viewData4(3);
    for (int i=0;i<viewData4.size();i++) {
      viewData4[i]=0.0;
    }
    DataArrayView myView4(viewData4,viewShape);
    values.push_back(myView4);

    DataTagged myData1(keys,values,myView1,FunctionSpace());

    values.clear();

    // default value for Data2
    DataArrayView::ValueType viewData5(3);
    for (int i=0;i<viewData5.size();i++) {
      viewData5[i]=1.0;
    }
    DataArrayView myView5(viewData5,viewShape);

    // value for tag "1" for Data2
    DataArrayView::ValueType viewData6(3);
    for (int i=0;i<viewData6.size();i++) {
      viewData6[i]=2.0;
    }
    DataArrayView myView6(viewData6,viewShape);
    values.push_back(myView6);

    // value for tag "2" for Data2
    DataArrayView::ValueType viewData7(3);
    for (int i=0;i<viewData7.size();i++) {
      viewData7[i]=3.0;
    }
    DataArrayView myView7(viewData7,viewShape);
    values.push_back(myView7);

    // value for tag "3" for Data2
    DataArrayView::ValueType viewData8(3);
    for (int i=0;i<viewData8.size();i++) {
      viewData8[i]=4.0;
    }
    DataArrayView myView8(viewData8,viewShape);
    values.push_back(myView8);

    DataTagged myData2(keys,values,myView5,FunctionSpace());

    // full slice

    std::pair<int, int> region_element;
    region_element.first=0;
    region_element.second=3;
    DataArrayView::RegionType region;
    region.push_back(region_element);

    myData1.setSlice(&myData2, region);

    //cout << myData1.toString() << endl;

    assert(myData1.getTagLookup().size()==3);

    assert(myData1.getLength()==12);

    DataArrayView myDataView = myData1.getDefaultValue();
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==0);
    assert(myDataView.getRank()==1);
    assert(myDataView.noValues()==3);
    assert(myDataView.getShape().size()==1);
    assert(myDataView(0)==1.0);
    assert(myDataView(1)==1.0);
    assert(myDataView(2)==1.0);

    myDataView = myData1.getDataPointByTag(1);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==3);
    assert(myDataView.getRank()==1);
    assert(myDataView.noValues()==3);
    assert(myDataView.getShape().size()==1);
    assert(myDataView(0)==2.0);
    assert(myDataView(1)==2.0);
    assert(myDataView(2)==2.0);

    myDataView = myData1.getDataPointByTag(2);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==6);
    assert(myDataView.getRank()==1);
    assert(myDataView.noValues()==3);
    assert(myDataView.getShape().size()==1);
    assert(myDataView(0)==3.0);
    assert(myDataView(1)==3.0);
    assert(myDataView(2)==3.0);

    myDataView = myData1.getDataPointByTag(3);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==9);
    assert(myDataView.getRank()==1);
    assert(myDataView.noValues()==3);
    assert(myDataView.getShape().size()==1);
    assert(myDataView(0)==4.0);
    assert(myDataView(1)==4.0);
    assert(myDataView(2)==4.0);

    // rank 1 slice

    viewShape.clear();
    viewShape.push_back(1);

    DataArrayView::ValueType viewData9(1);
    viewData9[0]=6.0;
    DataArrayView myView9(viewData9,viewShape);

    values.clear();

    DataArrayView::ValueType viewData10(1);
    viewData10[0]=7.0;
    DataArrayView myView10(viewData10,viewShape);
    values.push_back(myView10);

    DataArrayView::ValueType viewData11(1);
    viewData11[0]=8.0;
    DataArrayView myView11(viewData11,viewShape);
    values.push_back(myView11);

    DataArrayView::ValueType viewData12(1);
    viewData12[0]=9.0;
    DataArrayView myView12(viewData12,viewShape);
    values.push_back(myView12);

    DataTagged myData3(keys,values,myView9,FunctionSpace());

    region.clear();
    region_element.first=1;
    region_element.second=2;
    region.push_back(region_element);

    myData1.setSlice(&myData3, region);

    //cout << myData1.toString() << endl;

    assert(myData1.getTagLookup().size()==3);

    assert(myData1.getLength()==12);

    myDataView = myData1.getDefaultValue();
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==0);
    assert(myDataView.getRank()==1);
    assert(myDataView.noValues()==3);
    assert(myDataView.getShape().size()==1);
    assert(myDataView(0)==1.0);
    assert(myDataView(1)==6.0);
    assert(myDataView(2)==1.0);

    myDataView = myData1.getDataPointByTag(1);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==3);
    assert(myDataView.getRank()==1);
    assert(myDataView.noValues()==3);
    assert(myDataView.getShape().size()==1);
    assert(myDataView(0)==2.0);
    assert(myDataView(1)==7.0);
    assert(myDataView(2)==2.0);

    myDataView = myData1.getDataPointByTag(2);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==6);
    assert(myDataView.getRank()==1);
    assert(myDataView.noValues()==3);
    assert(myDataView.getShape().size()==1);
    assert(myDataView(0)==3.0);
    assert(myDataView(1)==8.0);
    assert(myDataView(2)==3.0);

    myDataView = myData1.getDataPointByTag(3);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==9);
    assert(myDataView.getRank()==1);
    assert(myDataView.noValues()==3);
    assert(myDataView.getShape().size()==1);
    assert(myDataView(0)==4.0);
    assert(myDataView(1)==9.0);
    assert(myDataView(2)==4.0);

    // scalar slice

    viewShape.clear();

    DataArrayView::ValueType viewData13(1);
    viewData13[0]=10.0;
    DataArrayView myView13(viewData13,viewShape);

    values.clear();

    DataArrayView::ValueType viewData14(1);
    viewData14[0]=11.0;
    DataArrayView myView14(viewData14,viewShape);
    values.push_back(myView14);

    DataArrayView::ValueType viewData15(2);
    viewData15[0]=12.0;
    DataArrayView myView15(viewData15,viewShape);
    values.push_back(myView15);

    DataArrayView::ValueType viewData16(3);
    viewData16[0]=13.0;
    DataArrayView myView16(viewData16,viewShape);
    values.push_back(myView16);

    DataTagged myData4(keys,values,myView13,FunctionSpace());

    region.clear();
    region_element.first=0;
    region_element.second=0;
    region.push_back(region_element);

    myData1.setSlice(&myData4, region);

    //cout << myData1.toString() << endl;

    myDataView = myData1.getDefaultValue();
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==0);
    assert(myDataView.getRank()==1);
    assert(myDataView.noValues()==3);
    assert(myDataView.getShape().size()==1);
    assert(myDataView(0)==10.0);
    assert(myDataView(1)==6.0);
    assert(myDataView(2)==1.0);

    myDataView = myData1.getDataPointByTag(1);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==3);
    assert(myDataView.getRank()==1);
    assert(myDataView.noValues()==3);
    assert(myDataView.getShape().size()==1);
    assert(myDataView(0)==11.0);
    assert(myDataView(1)==7.0);
    assert(myDataView(2)==2.0);

    myDataView = myData1.getDataPointByTag(2);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==6);
    assert(myDataView.getRank()==1);
    assert(myDataView.noValues()==3);
    assert(myDataView.getShape().size()==1);
    assert(myDataView(0)==12.0);
    assert(myDataView(1)==8.0);
    assert(myDataView(2)==3.0);

    myDataView = myData1.getDataPointByTag(3);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==9);
    assert(myDataView.getRank()==1);
    assert(myDataView.noValues()==3);
    assert(myDataView.getShape().size()==1);
    assert(myDataView(0)==13.0);
    assert(myDataView(1)==9.0);
    assert(myDataView(2)==4.0);

  }

  {

    cout << "\tTest slicing DataTagged with rank 3 values and three tags." << endl;

    DataTagged::TagListType keys;
    keys.push_back(1);
    keys.push_back(2);
    keys.push_back(3);

    DataTagged::ValueListType values;

    DataArrayView::ShapeType viewShape;
    viewShape.push_back(3);
    viewShape.push_back(3);
    viewShape.push_back(3);

    // default value for Data1
    DataArrayView::ValueType viewData1(27);
    for (int i=0;i<viewData1.size();i++) {
      viewData1[i]=0.0;
    }
    DataArrayView myView1(viewData1,viewShape);

    // value for tag "1" for Data1
    DataArrayView::ValueType viewData2(27);
    for (int i=0;i<viewData2.size();i++) {
      viewData2[i]=0.0;
    }
    DataArrayView myView2(viewData2,viewShape);
    values.push_back(myView2);

    // value for tag "2" for Data1
    DataArrayView::ValueType viewData3(27);
    for (int i=0;i<viewData3.size();i++) {
      viewData3[i]=0.0;
    }
    DataArrayView myView3(viewData3,viewShape);
    values.push_back(myView3);

    // value for tag "3" for Data1
    DataArrayView::ValueType viewData4(27);
    for (int i=0;i<viewData4.size();i++) {
      viewData4[i]=0.0;
    }
    DataArrayView myView4(viewData4,viewShape);
    values.push_back(myView4);

    DataTagged myData1(keys,values,myView1,FunctionSpace());

    values.clear();

    // default value for Data2
    DataArrayView::ValueType viewData5(27);
    for (int i=0;i<viewData5.size();i++) {
      viewData5[i]=1.0;
    }
    DataArrayView myView5(viewData5,viewShape);

    // value for tag "1" for Data2
    DataArrayView::ValueType viewData6(27);
    for (int i=0;i<viewData6.size();i++) {
      viewData6[i]=2.0;
    }
    DataArrayView myView6(viewData6,viewShape);
    values.push_back(myView6);

    // value for tag "2" for Data2
    DataArrayView::ValueType viewData7(27);
    for (int i=0;i<viewData7.size();i++) {
      viewData7[i]=3.0;
    }
    DataArrayView myView7(viewData7,viewShape);
    values.push_back(myView7);

    // value for tag "3" for Data2
    DataArrayView::ValueType viewData8(27);
    for (int i=0;i<viewData8.size();i++) {
      viewData8[i]=4.0;
    }
    DataArrayView myView8(viewData8,viewShape);
    values.push_back(myView8);

    DataTagged myData2(keys,values,myView5,FunctionSpace());

    // full slice

    std::pair<int, int> region_element;
    region_element.first=0;
    region_element.second=3;
    DataArrayView::RegionType region;
    region.push_back(region_element);
    region.push_back(region_element);
    region.push_back(region_element);

    myData1.setSlice(&myData2, region);

    //cout << myData1.toString() << endl;

    assert(myData1.getTagLookup().size()==3);

    assert(myData1.getLength()==108);

    DataArrayView myDataView = myData1.getDefaultValue();
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==0);
    assert(myDataView.getRank()==3);
    assert(myDataView.noValues()==27);
    assert(myDataView.getShape().size()==3);
    assert(myDataView(0,0,0)==1.0);
    assert(myDataView(1,0,0)==1.0);
    assert(myDataView(2,0,0)==1.0);

    myDataView = myData1.getDataPointByTag(1);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==27);
    assert(myDataView.getRank()==3);
    assert(myDataView.noValues()==27);
    assert(myDataView.getShape().size()==3);
    assert(myDataView(0,0,0)==2.0);
    assert(myDataView(1,0,0)==2.0);
    assert(myDataView(2,0,0)==2.0);

    myDataView = myData1.getDataPointByTag(2);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==54);
    assert(myDataView.getRank()==3);
    assert(myDataView.noValues()==27);
    assert(myDataView.getShape().size()==3);
    assert(myDataView(0,0,0)==3.0);
    assert(myDataView(1,0,0)==3.0);
    assert(myDataView(2,0,0)==3.0);

    myDataView = myData1.getDataPointByTag(3);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==81);
    assert(myDataView.getRank()==3);
    assert(myDataView.noValues()==27);
    assert(myDataView.getShape().size()==3);
    assert(myDataView(0,0,0)==4.0);
    assert(myDataView(1,0,0)==4.0);
    assert(myDataView(2,0,0)==4.0);

    // rank 1 slice

    viewShape.clear();
    viewShape.push_back(3);

    DataArrayView::ValueType viewData9(3);
    for (int i=0;i<viewData9.size();i++) {
      viewData9[i]=6.0;
    }
    DataArrayView myView9(viewData9,viewShape);

    values.clear();

    DataArrayView::ValueType viewData10(3);
    for (int i=0;i<viewData10.size();i++) {
      viewData10[i]=7.0;
    }
    DataArrayView myView10(viewData10,viewShape);
    values.push_back(myView10);

    DataArrayView::ValueType viewData11(3);
    for (int i=0;i<viewData11.size();i++) {
      viewData11[i]=8.0;
    }
    DataArrayView myView11(viewData11,viewShape);
    values.push_back(myView11);

    DataArrayView::ValueType viewData12(3);
    for (int i=0;i<viewData12.size();i++) {
      viewData12[i]=9.0;
    }
    DataArrayView myView12(viewData12,viewShape);
    values.push_back(myView12);

    DataTagged myData3(keys,values,myView9,FunctionSpace());

    region.clear();
    region_element.first=0;
    region_element.second=3;
    region.push_back(region_element);
    region_element.second=0;
    region.push_back(region_element);
    region.push_back(region_element);

    myData1.setSlice(&myData3, region);

    //cout << myData1.toString() << endl;

    assert(myData1.getTagLookup().size()==3);

    assert(myData1.getLength()==108);

    myDataView = myData1.getDefaultValue();
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==0);
    assert(myDataView.getRank()==3);
    assert(myDataView.noValues()==27);
    assert(myDataView.getShape().size()==3);
    assert(myDataView(0,0,0)==6.0);
    assert(myDataView(1,0,0)==6.0);
    assert(myDataView(2,0,0)==6.0);

    myDataView = myData1.getDataPointByTag(1);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==27);
    assert(myDataView.getRank()==3);
    assert(myDataView.noValues()==27);
    assert(myDataView.getShape().size()==3);
    assert(myDataView(0,0,0)==7.0);
    assert(myDataView(1,0,0)==7.0);
    assert(myDataView(2,0,0)==7.0);

    myDataView = myData1.getDataPointByTag(2);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==54);
    assert(myDataView.getRank()==3);
    assert(myDataView.noValues()==27);
    assert(myDataView.getShape().size()==3);
    assert(myDataView(0,0,0)==8.0);
    assert(myDataView(1,0,0)==8.0);
    assert(myDataView(2,0,0)==8.0);

    myDataView = myData1.getDataPointByTag(3);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==81);
    assert(myDataView.getRank()==3);
    assert(myDataView.noValues()==27);
    assert(myDataView.getShape().size()==3);
    assert(myDataView(0,0,0)==9.0);
    assert(myDataView(1,0,0)==9.0);
    assert(myDataView(2,0,0)==9.0);

    // scalar slice

    viewShape.clear();

    DataArrayView::ValueType viewData13(1);
    viewData13[0]=10.0;
    DataArrayView myView13(viewData13,viewShape);

    values.clear();

    DataArrayView::ValueType viewData14(1);
    viewData14[0]=11.0;
    DataArrayView myView14(viewData14,viewShape);
    values.push_back(myView14);

    DataArrayView::ValueType viewData15(2);
    viewData15[0]=12.0;
    DataArrayView myView15(viewData15,viewShape);
    values.push_back(myView15);

    DataArrayView::ValueType viewData16(3);
    viewData16[0]=13.0;
    DataArrayView myView16(viewData16,viewShape);
    values.push_back(myView16);

    DataTagged myData4(keys,values,myView13,FunctionSpace());

    region.clear();
    region_element.first=0;
    region_element.second=0;
    region.push_back(region_element);
    region.push_back(region_element);
    region.push_back(region_element);

    myData1.setSlice(&myData4, region);

    //cout << myData1.toString() << endl;

    myDataView = myData1.getDefaultValue();
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==0);
    assert(myDataView.getRank()==3);
    assert(myDataView.noValues()==27);
    assert(myDataView.getShape().size()==3);
    assert(myDataView(0,0,0)==10.0);

    myDataView = myData1.getDataPointByTag(1);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==27);
    assert(myDataView.getRank()==3);
    assert(myDataView.noValues()==27);
    assert(myDataView.getShape().size()==3);
    assert(myDataView(0,0,0)==11.0);

    myDataView = myData1.getDataPointByTag(2);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==54);
    assert(myDataView.getRank()==3);
    assert(myDataView.noValues()==27);
    assert(myDataView.getShape().size()==3);
    assert(myDataView(0,0,0)==12.0);

    myDataView = myData1.getDataPointByTag(3);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==81);
    assert(myDataView.getRank()==3);
    assert(myDataView.noValues()==27);
    assert(myDataView.getShape().size()==3);
    assert(myDataView(0,0,0)==13.0);

  }

  {

    cout << "\tTest slicing DataTagged with scalar values and three *mismatched* tags." << endl;

    DataTagged::TagListType keys1;
    keys1.push_back(1);
    keys1.push_back(2);
    keys1.push_back(3);

    DataTagged::TagListType keys2;
    keys2.push_back(3);
    keys2.push_back(4);
    keys2.push_back(5);

    DataTagged::ValueListType values;

    DataArrayView::ShapeType viewShape;

    // default value for Data1
    DataArrayView::ValueType viewData1(1);
    viewData1[0]=0.0;
    DataArrayView myView1(viewData1,viewShape);

    // value for tag "1" for Data1
    DataArrayView::ValueType viewData2(1);
    viewData2[0]=0.0;
    DataArrayView myView2(viewData2,viewShape);
    values.push_back(myView2);

    // value for tag "2" for Data1
    DataArrayView::ValueType viewData5(1);
    viewData5[0]=0.0;
    DataArrayView myView5(viewData5,viewShape);
    values.push_back(myView5);

    // value for tag "3" for Data1
    DataArrayView::ValueType viewData6(1);
    viewData6[0]=0.0;
    DataArrayView myView6(viewData6,viewShape);
    values.push_back(myView6);

    DataTagged myData1(keys1,values,myView1,FunctionSpace());

    values.clear();

    // default value for Data2
    DataArrayView::ValueType viewData3(1);
    viewData3[0]=1.0;
    DataArrayView myView3(viewData3,viewShape);

    // value for tag "3" for Data2
    DataArrayView::ValueType viewData4(1);
    viewData4[0]=2.0;
    DataArrayView myView4(viewData4,viewShape);
    values.push_back(myView4);

    // value for tag "4" for Data2
    DataArrayView::ValueType viewData7(1);
    viewData7[0]=3.0;
    DataArrayView myView7(viewData7,viewShape);
    values.push_back(myView7);

    // value for tag "5" for Data2
    DataArrayView::ValueType viewData8(1);
    viewData8[0]=4.0;
    DataArrayView myView8(viewData8,viewShape);
    values.push_back(myView8);

    DataTagged myData2(keys2,values,myView3,FunctionSpace());

    //cout << myData1.toString() << endl;
    //cout << myData2.toString() << endl;

    // full slice

    DataArrayView::RegionType region;

    myData1.setSlice(&myData2, region);

    //cout << myData1.toString() << endl;

    assert(myData1.getTagLookup().size()==5);

    assert(myData1.getLength()==6);

    DataArrayView myDataView = myData1.getDefaultValue();
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==0);
    assert(myDataView.getRank()==0);
    assert(myDataView.noValues()==1);
    assert(myDataView.getShape().size()==0);
    assert(myDataView()==1.0);

    myDataView = myData1.getDataPointByTag(1);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==1);
    assert(myDataView.getRank()==0);
    assert(myDataView.noValues()==1);
    assert(myDataView.getShape().size()==0);
    assert(myDataView()==1.0);

    myDataView = myData1.getDataPointByTag(2);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==2);
    assert(myDataView.getRank()==0);
    assert(myDataView.noValues()==1);
    assert(myDataView.getShape().size()==0);
    assert(myDataView()==1.0);

    myDataView = myData1.getDataPointByTag(3);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==3);
    assert(myDataView.getRank()==0);
    assert(myDataView.noValues()==1);
    assert(myDataView.getShape().size()==0);
    assert(myDataView()==2.0);

    myDataView = myData1.getDataPointByTag(4);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==4);
    assert(myDataView.getRank()==0);
    assert(myDataView.noValues()==1);
    assert(myDataView.getShape().size()==0);
    assert(myDataView()==3.0);

    myDataView = myData1.getDataPointByTag(5);
    assert(!myDataView.isEmpty());
    assert(myDataView.getOffset()==5);
    assert(myDataView.getRank()==0);
    assert(myDataView.noValues()==1);
    assert(myDataView.getShape().size()==0);
    assert(myDataView()==4.0);

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
  testSuite->addTest (new TestCaller< DataTaggedTestCase>("testReshape",&DataTaggedTestCase::testReshape));
  testSuite->addTest (new TestCaller< DataTaggedTestCase>("testGetSlice",&DataTaggedTestCase::testGetSlice));
  testSuite->addTest (new TestCaller< DataTaggedTestCase>("testSetSlice",&DataTaggedTestCase::testSetSlice));
  return testSuite;
}
