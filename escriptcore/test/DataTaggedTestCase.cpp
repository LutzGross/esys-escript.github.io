
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

#include <escript/DataTypes.h>

#include "DataTaggedTestCase.h"

#include <escript/BinaryDataReadyOps.h>
#include <escript/DataConstant.h>
#include <escript/DataFactory.h>
#include <escript/DataTagged.h>
#include <escript/DataVector.h>
#include <escript/EsysException.h>
#include <escript/FunctionSpace.h>
#include <escript/FunctionSpaceFactory.h>

#include <cppunit/TestCaller.h>

#include <algorithm>
#include <functional>
#include <iostream>


using namespace CppUnit;
using namespace escript;
using namespace std;
using namespace escript::DataTypes;

namespace {

RealVectorType::const_reference
getRefRO(DataTagged& data,int offset, int i, int j, int k)
{
   return data.getVectorRO()[offset+getRelIndex(data.getShape(),i,j,k)];
}

RealVectorType::const_reference
getRefRO(const DataTagged& data,int offset, int i)
{
   return data.getVectorRO()[offset+getRelIndex(data.getShape(),i)];
}

}

namespace
{
    DataTagged makeTagged()
    {
        int a[1]={0};
	DataTypes::RealVectorType v;
	v.resize(1, 0.0,1);
        return DataTagged(FunctionSpace(), DataTypes::scalarShape, a, v);
    }
}


void DataTaggedTestCase::testAddTaggedValues() {

  cout << endl;

  {

    cout << "\tTest adding one key with empty value list to default DataTagged." << endl;
    DataTagged myData=makeTagged();

    DataTagged::TagListType keys;
    keys.push_back(1);

    DataTagged::FloatBatchType values;

    myData.addTaggedValues(keys,values,DataTypes::scalarShape);

    CPPUNIT_ASSERT(myData.isCurrentTag(1));

    CPPUNIT_ASSERT(myData.getTagLookup().size()==1);

    CPPUNIT_ASSERT(myData.getLength()==2);

    // data-point 0 has tag number 1 by default
    CPPUNIT_ASSERT(myData.getTagNumber(0)==1);

    CPPUNIT_ASSERT(myData.getPointOffset(0,0)==1);

    CPPUNIT_ASSERT(myData.getRank()==0);
    CPPUNIT_ASSERT(myData.getNoValues()==1);
    CPPUNIT_ASSERT(myData.getShape().size()==0);


    int offset=myData.getPointOffset(0,0);
    CPPUNIT_ASSERT(offset==1);
    CPPUNIT_ASSERT(myData.getDataAtOffsetRO(offset)==0.0);

    offset=myData.getOffsetForTag(1);
    CPPUNIT_ASSERT(offset==1);
    CPPUNIT_ASSERT(myData.getDataAtOffsetRO(offset)==0.0);

    offset=myData.getDefaultOffset();
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(myData.getDataAtOffsetRO(offset)==0.0);

    // use a non-existent tag so we get a pointer to
    // the first element of the data array
    double* sampleData=myData.getSampleDataByTag(9);
    for (int i=0; i<myData.getLength(); i++) {
      CPPUNIT_ASSERT(sampleData[i]==0);
    }

  }

  {

    cout << "\tTest adding one key with one value to default DataTagged." << endl;
    DataTagged myData=makeTagged();

    DataTagged::TagListType keys;
    keys.push_back(1);

    DataTagged::FloatBatchType values;

    DataTypes::ShapeType viewShape;
    values.push_back(1.0);

    myData.addTaggedValues(keys,values,viewShape);

    CPPUNIT_ASSERT(myData.isCurrentTag(1));

    CPPUNIT_ASSERT(myData.getTagLookup().size()==1);

    CPPUNIT_ASSERT(myData.getLength()==2);

    // data-point 0 has tag number 1 by default
    CPPUNIT_ASSERT(myData.getTagNumber(0)==1);

    CPPUNIT_ASSERT(myData.getPointOffset(0,0)==1);

    CPPUNIT_ASSERT(myData.getRank()==0);
    CPPUNIT_ASSERT(myData.getNoValues()==1);
    CPPUNIT_ASSERT(myData.getShape().size()==0);


    int offset=myData.getPointOffset(0,0);
    CPPUNIT_ASSERT(offset==1);
    CPPUNIT_ASSERT(myData.getDataAtOffsetRO(offset)==1.0);

    offset=myData.getOffsetForTag(1);
    CPPUNIT_ASSERT(offset==1);
    CPPUNIT_ASSERT(myData.getDataAtOffsetRO(offset)==1.0);

    offset=myData.getDefaultOffset();
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(myData.getDataAtOffsetRO(offset)==0.0);

    // use a non-existent tag so we get a pointer to
    // the first element of the data array
    double* sampleData=myData.getSampleDataByTag(9);
    for (int i=0; i<myData.getLength(); i++) {
      CPPUNIT_ASSERT(sampleData[i]==i);
    }

  }

  {

    cout << "\tTest adding three keys with one value to default DataTagged." << endl;
    DataTagged myData=makeTagged();

    DataTagged::TagListType keys;
    keys.push_back(1);
    keys.push_back(2);
    keys.push_back(3);

    DataTagged::FloatBatchType values;

    DataTypes::ShapeType viewShape;
/*    DataTypes::RealVectorType viewData(1);
    viewData[0]=1.0;
    DataArrayView myView(viewData,viewShape);*/
    values.push_back(1.0);

    myData.addTaggedValues(keys,values,viewShape);

    CPPUNIT_ASSERT(myData.isCurrentTag(1));
    CPPUNIT_ASSERT(myData.isCurrentTag(2));
    CPPUNIT_ASSERT(myData.isCurrentTag(3));

    CPPUNIT_ASSERT(myData.getTagLookup().size()==3);

    CPPUNIT_ASSERT(myData.getLength()==4);

    // data-point 0 has tag number 1 by default
    CPPUNIT_ASSERT(myData.getTagNumber(0)==1);

    CPPUNIT_ASSERT(myData.getPointOffset(0,0)==1);

    CPPUNIT_ASSERT(myData.getRank()==0);
    CPPUNIT_ASSERT(myData.getNoValues()==1);
    CPPUNIT_ASSERT(myData.getShape().size()==0);


    int offset=myData.getPointOffset(0,0);
    CPPUNIT_ASSERT(offset==1);
    CPPUNIT_ASSERT(myData.getDataAtOffsetRO(offset)==1.0);

    offset=myData.getOffsetForTag(1);
    CPPUNIT_ASSERT(offset==1);
    CPPUNIT_ASSERT(myData.getDataAtOffsetRO(offset)==1.0);

    offset=myData.getOffsetForTag(2);
    CPPUNIT_ASSERT(offset==2);
    CPPUNIT_ASSERT(myData.getDataAtOffsetRO(offset)==1.0);

    offset=myData.getOffsetForTag(3);
    CPPUNIT_ASSERT(offset==3);
    CPPUNIT_ASSERT(myData.getDataAtOffsetRO(offset)==1.0);

    offset=myData.getDefaultOffset();
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(myData.getDataAtOffsetRO(offset)==0.0);

    // use a non-existent tag so we get a pointer to
    // the first element of the data array
    double* sampleData=myData.getSampleDataByTag(9);
    for (int i=0; i<myData.getLength(); i++) {
      if (i==0) {
        CPPUNIT_ASSERT(sampleData[i]==0);
      } else {
        CPPUNIT_ASSERT(sampleData[i]==1);
      }
    }

  }

  {

    cout << "\tTest adding three keys with three values to default DataTagged." << endl;
    DataTagged myData=makeTagged();

    DataTagged::TagListType keys;
    keys.push_back(1);
    keys.push_back(2);
    keys.push_back(3);

    DataTagged::FloatBatchType values;

    DataTypes::ShapeType viewShape;
/*    DataTypes::RealVectorType viewData1(1);
    viewData1[0]=1.0;
    DataTypes::RealVectorType viewData2(1);
    viewData2[0]=2.0;
    DataTypes::RealVectorType viewData3(1);
    viewData3[0]=3.0;
    DataArrayView myView1(viewData1,viewShape);
    DataArrayView myView2(viewData2,viewShape);
    DataArrayView myView3(viewData3,viewShape);*/
    values.push_back(1.0);
    values.push_back(2.0);
    values.push_back(3.0);

    myData.addTaggedValues(keys,values,viewShape);

    CPPUNIT_ASSERT(myData.isCurrentTag(1));
    CPPUNIT_ASSERT(myData.isCurrentTag(2));
    CPPUNIT_ASSERT(myData.isCurrentTag(3));

    CPPUNIT_ASSERT(myData.getTagLookup().size()==3);

    CPPUNIT_ASSERT(myData.getLength()==4);

    // data-point 0 has tag number 1 by default
    CPPUNIT_ASSERT(myData.getTagNumber(0)==1);

    CPPUNIT_ASSERT(myData.getPointOffset(0,0)==1);

    CPPUNIT_ASSERT(myData.getRank()==0);
    CPPUNIT_ASSERT(myData.getNoValues()==1);
    CPPUNIT_ASSERT(myData.getShape().size()==0);


    int offset=myData.getPointOffset(0,0);
    CPPUNIT_ASSERT(offset==1);
    CPPUNIT_ASSERT(myData.getDataAtOffsetRO(offset)==1.0);

    offset=myData.getOffsetForTag(1);
    CPPUNIT_ASSERT(offset==1);
    CPPUNIT_ASSERT(myData.getDataAtOffsetRO(offset)==1.0);

    offset=myData.getOffsetForTag(2);
    CPPUNIT_ASSERT(offset==2);
    CPPUNIT_ASSERT(myData.getDataAtOffsetRO(offset)==2.0);

    offset=myData.getOffsetForTag(3);
    CPPUNIT_ASSERT(offset==3);
    CPPUNIT_ASSERT(myData.getDataAtOffsetRO(offset)==3.0);

    offset=myData.getDefaultOffset();
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(myData.getDataAtOffsetRO(offset)==0.0);

    // use a non-existent tag so we get a pointer to
    // the first element of the data array
    double* sampleData=myData.getSampleDataByTag(9);
    for (int i=0; i<myData.getLength(); i++) {
      CPPUNIT_ASSERT(sampleData[i]==i);
    }

  }

  {

    cout << "\tTest adding one key with empty value list to DataTagged with default value only." << endl;

    DataTypes::ShapeType viewShape;
    viewShape.push_back(3);

    DataTagged::TagListType keys;

    DataTagged::FloatBatchType values;

    DataTypes::RealVectorType viewData(3);
    for (int i=0;i<viewShape[0];i++) {
      viewData[i]=i;
    }

    DataTagged myData(FunctionSpace(),viewShape,viewData);

    keys.push_back(1);
    values.clear();

    myData.addTaggedValues(keys,values,viewShape);

    CPPUNIT_ASSERT(myData.isCurrentTag(1));

    CPPUNIT_ASSERT(myData.getTagLookup().size()==1);

    CPPUNIT_ASSERT(myData.getLength()==6);

    // data-point 0 has tag number 1 by default
    CPPUNIT_ASSERT(myData.getTagNumber(0)==1);

    CPPUNIT_ASSERT(myData.getPointOffset(0,0)==3);

    CPPUNIT_ASSERT(myData.getRank()==1);
    CPPUNIT_ASSERT(myData.getNoValues()==3);
    CPPUNIT_ASSERT(myData.getShape().size()==1);


    int offset=myData.getPointOffset(0,0);
    CPPUNIT_ASSERT(offset==3);
    CPPUNIT_ASSERT(getRefRO(myData,offset,0)==0);
    CPPUNIT_ASSERT(getRefRO(myData,offset,1)==1);
    CPPUNIT_ASSERT(getRefRO(myData,offset,2)==2);

    offset=myData.getOffsetForTag(1);
    CPPUNIT_ASSERT(offset==3);
    CPPUNIT_ASSERT(getRefRO(myData,offset,0)==0);
    CPPUNIT_ASSERT(getRefRO(myData,offset,1)==1);
    CPPUNIT_ASSERT(getRefRO(myData,offset,2)==2);

    offset=myData.getDefaultOffset();
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(getRefRO(myData,offset,0)==0);
    CPPUNIT_ASSERT(getRefRO(myData,offset,1)==1);
    CPPUNIT_ASSERT(getRefRO(myData,offset,2)==2);

    // use a non-existent tag so we get a pointer to
    // the first element of the data array
    double* sampleData=myData.getSampleDataByTag(9);
    for (int i=0; i<myData.getLength(); i++) {
      CPPUNIT_ASSERT(sampleData[i]==i%3);
    }

  }

  {

    cout << "\tTest adding one key with one value to DataTagged with default value only." << endl;

    DataTypes::ShapeType viewShape;
    viewShape.push_back(3);

    DataTagged::TagListType keys;

    DataTagged::FloatBatchType values;

    DataTypes::RealVectorType viewData(3);
    for (int i=0;i<viewShape[0];i++) {
      viewData[i]=i;
    }

    DataTagged myData(FunctionSpace(),viewShape,viewData);

    keys.push_back(1);

    for (int i=0;i<viewShape[0];i++) {
	values.push_back(i+3);
    }

    myData.addTaggedValues(keys,values,viewShape);

    CPPUNIT_ASSERT(myData.isCurrentTag(1));

    CPPUNIT_ASSERT(myData.getTagLookup().size()==1);

    CPPUNIT_ASSERT(myData.getLength()==6);

    // data-point 0 has tag number 1 by default
    CPPUNIT_ASSERT(myData.getTagNumber(0)==1);

    CPPUNIT_ASSERT(myData.getPointOffset(0,0)==3);

    CPPUNIT_ASSERT(myData.getRank()==1);
    CPPUNIT_ASSERT(myData.getNoValues()==3);
    CPPUNIT_ASSERT(myData.getShape().size()==1);


    int offset=myData.getPointOffset(0,0);
    CPPUNIT_ASSERT(offset==3);
    CPPUNIT_ASSERT(getRefRO(myData,offset,0)==3);
    CPPUNIT_ASSERT(getRefRO(myData,offset,1)==4);
    CPPUNIT_ASSERT(getRefRO(myData,offset,2)==5);

    offset=myData.getOffsetForTag(1);
    CPPUNIT_ASSERT(offset==3);
    CPPUNIT_ASSERT(getRefRO(myData,offset,0)==3);
    CPPUNIT_ASSERT(getRefRO(myData,offset,1)==4);
    CPPUNIT_ASSERT(getRefRO(myData,offset,2)==5);

    offset=myData.getDefaultOffset();
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(getRefRO(myData,offset,0)==0);
    CPPUNIT_ASSERT(getRefRO(myData,offset,1)==1);
    CPPUNIT_ASSERT(getRefRO(myData,offset,2)==2);

    // use a non-existent tag so we get a pointer to
    // the first element of the data array
    double* sampleData=myData.getSampleDataByTag(9);
    for (int i=0; i<myData.getLength(); i++) {
      CPPUNIT_ASSERT(sampleData[i]==i);
    }

  }

  {

    cout << "\tTest adding three keys with one value to DataTagged with default value only." << endl;

    DataTypes::ShapeType viewShape;
    viewShape.push_back(3);

    DataTagged::TagListType keys;

    DataTagged::FloatBatchType values;

    DataTypes::RealVectorType viewData(3);
    for (int i=0;i<viewShape[0];i++) {
      viewData[i]=i;
    }

    DataTagged myData(FunctionSpace(),viewShape,viewData);

    keys.push_back(1);
    keys.push_back(2);
    keys.push_back(3);

    for (int i=0;i<viewShape[0];i++) {
	values.push_back(3);
    }

    myData.addTaggedValues(keys,values,viewShape);

    CPPUNIT_ASSERT(myData.isCurrentTag(1));
    CPPUNIT_ASSERT(myData.isCurrentTag(2));
    CPPUNIT_ASSERT(myData.isCurrentTag(3));

    CPPUNIT_ASSERT(myData.getTagLookup().size()==3);

    CPPUNIT_ASSERT(myData.getLength()==12);

    // data-point 0 has tag number 1 by default
    CPPUNIT_ASSERT(myData.getTagNumber(0)==1);

    CPPUNIT_ASSERT(myData.getPointOffset(0,0)==3);

    CPPUNIT_ASSERT(myData.getRank()==1);
    CPPUNIT_ASSERT(myData.getNoValues()==3);
    CPPUNIT_ASSERT(myData.getShape().size()==1);

    int offset=myData.getPointOffset(0,0);
    CPPUNIT_ASSERT(offset==3);
    CPPUNIT_ASSERT(getRefRO(myData,offset,0)==3);
    CPPUNIT_ASSERT(getRefRO(myData,offset,1)==3);
    CPPUNIT_ASSERT(getRefRO(myData,offset,2)==3);

    offset=myData.getOffsetForTag(1);
    CPPUNIT_ASSERT(offset==3);
    CPPUNIT_ASSERT(getRefRO(myData,offset,0)==3);
    CPPUNIT_ASSERT(getRefRO(myData,offset,1)==3);
    CPPUNIT_ASSERT(getRefRO(myData,offset,2)==3);

    offset=myData.getOffsetForTag(2);
    CPPUNIT_ASSERT(offset==6);
    CPPUNIT_ASSERT(getRefRO(myData,offset,0)==3);
    CPPUNIT_ASSERT(getRefRO(myData,offset,1)==3);
    CPPUNIT_ASSERT(getRefRO(myData,offset,2)==3);

    offset=myData.getOffsetForTag(3);
    CPPUNIT_ASSERT(offset==9);
    CPPUNIT_ASSERT(getRefRO(myData,offset,0)==3);
    CPPUNIT_ASSERT(getRefRO(myData,offset,1)==3);
    CPPUNIT_ASSERT(getRefRO(myData,offset,2)==3);

    offset=myData.getDefaultOffset();
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(getRefRO(myData,offset,0)==0);
    CPPUNIT_ASSERT(getRefRO(myData,offset,1)==1);
    CPPUNIT_ASSERT(getRefRO(myData,offset,2)==2);

    // use a non-existent tag so we get a pointer to
    // the first element of the data array
    double* sampleData=myData.getSampleDataByTag(9);
    for (int i=0; i<myData.getLength(); i++) {
      if (i<3) {
        CPPUNIT_ASSERT(sampleData[i]==i);
      } else {
        CPPUNIT_ASSERT(sampleData[i]==3);
      }
    }

  }

  {

    cout << "\tTest adding three keys with three values to DataTagged with default value only." << endl;

    DataTypes::ShapeType viewShape;
    viewShape.push_back(3);

    DataTagged::TagListType keys;

    DataTagged::FloatBatchType values;

    DataTypes::RealVectorType viewData(3);
    for (int i=0;i<viewShape[0];i++) {
      viewData[i]=i;
    }

    DataTagged myData(FunctionSpace(),viewShape,viewData);

    keys.push_back(1);
    keys.push_back(2);
    keys.push_back(3);

    for (int i=0;i<viewShape[0];i++) {
	values.push_back(i+1);
    }

    for (int i=0;i<viewShape[0];i++) {
	values.push_back(i+2);
    }

    for (int i=0;i<viewShape[0];i++) {
	values.push_back(i+3);
    }

    myData.addTaggedValues(keys,values,viewShape);

    CPPUNIT_ASSERT(myData.isCurrentTag(1));
    CPPUNIT_ASSERT(myData.isCurrentTag(2));
    CPPUNIT_ASSERT(myData.isCurrentTag(3));

    CPPUNIT_ASSERT(myData.getTagLookup().size()==3);

    CPPUNIT_ASSERT(myData.getLength()==12);

    // data-point 0 has tag number 1 by default
    CPPUNIT_ASSERT(myData.getTagNumber(0)==1);

    CPPUNIT_ASSERT(myData.getPointOffset(0,0)==3);

    CPPUNIT_ASSERT(myData.getRank()==1);
    CPPUNIT_ASSERT(myData.getNoValues()==3);
    CPPUNIT_ASSERT(myData.getShape().size()==1);


    int offset=myData.getPointOffset(0,0);
    CPPUNIT_ASSERT(offset==3);
    CPPUNIT_ASSERT(getRefRO(myData,offset,0)==1);
    CPPUNIT_ASSERT(getRefRO(myData,offset,1)==2);
    CPPUNIT_ASSERT(getRefRO(myData,offset,2)==3);

    offset=myData.getOffsetForTag(1);
    CPPUNIT_ASSERT(offset==3);
    CPPUNIT_ASSERT(getRefRO(myData,offset,0)==1);
    CPPUNIT_ASSERT(getRefRO(myData,offset,1)==2);
    CPPUNIT_ASSERT(getRefRO(myData,offset,2)==3);

    offset=myData.getOffsetForTag(2);
    CPPUNIT_ASSERT(offset==6);
    CPPUNIT_ASSERT(getRefRO(myData,offset,0)==2);
    CPPUNIT_ASSERT(getRefRO(myData,offset,1)==3);
    CPPUNIT_ASSERT(getRefRO(myData,offset,2)==4);

    offset=myData.getOffsetForTag(3);
    CPPUNIT_ASSERT(offset==9);
    CPPUNIT_ASSERT(getRefRO(myData,offset,0)==3);
    CPPUNIT_ASSERT(getRefRO(myData,offset,1)==4);
    CPPUNIT_ASSERT(getRefRO(myData,offset,2)==5);

    offset=myData.getDefaultOffset();
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(getRefRO(myData,offset,0)==0);
    CPPUNIT_ASSERT(getRefRO(myData,offset,1)==1);
    CPPUNIT_ASSERT(getRefRO(myData,offset,2)==2);

    // use a non-existent tag so we get a pointer to
    // the first element of the data array
    double* sampleData=myData.getSampleDataByTag(9);
    for (int i=0; i<myData.getLength(); i++) {
      if (i<3) {
        CPPUNIT_ASSERT(sampleData[i]==i);
      } else if ((i>=3) && (i<6)) {
        CPPUNIT_ASSERT(sampleData[i]==i-2);
      } else if ((i>=6) && (i<9)) {
        CPPUNIT_ASSERT(sampleData[i]==i-4);
      } else  {
        CPPUNIT_ASSERT(sampleData[i]==i-6);
      }
    }

  }

  {

    cout << "\tTest adding one key with empty value list to DataTagged with three tags." << endl;

    DataTagged::TagListType keys;
    keys.push_back(1);
    keys.push_back(2);
    keys.push_back(3);

    DataTagged::FloatBatchType values;

    DataTypes::ShapeType viewShape;
    viewShape.push_back(3);

    // default value
    DataTypes::RealVectorType viewData(3*4);
    for (int i=0;i<viewShape[0];i++) {
      viewData[i]=i;
    }

    // value for tag "1"
    for (int i=0;i<viewShape[0];i++) {
      viewData[viewShape[0]+i]=i+1.0;
    }

    // value for tag "2"
    for (int i=0;i<viewShape[0];i++) {
	viewData[2*viewShape[0]+i]=i+2.0;
    }

    // value for tag "3"
    for (int i=0;i<viewShape[0];i++) {
	viewData[3*viewShape[0]+i]=i+3.0;
    }

    DataTagged myData(FunctionSpace(),viewShape,keys,viewData);

    keys.clear();
    keys.push_back(4);
    values.clear();

    myData.addTaggedValues(keys,values,viewShape);

    CPPUNIT_ASSERT(myData.isCurrentTag(4));

    CPPUNIT_ASSERT(myData.getTagLookup().size()==4);

    CPPUNIT_ASSERT(myData.getLength()==15);

    CPPUNIT_ASSERT(myData.getRank()==1);
    CPPUNIT_ASSERT(myData.getNoValues()==3);
    CPPUNIT_ASSERT(myData.getShape().size()==1);

    int offset=myData.getOffsetForTag(4);
    CPPUNIT_ASSERT(offset==12);
    CPPUNIT_ASSERT(getRefRO(myData,offset,0)==0);
    CPPUNIT_ASSERT(getRefRO(myData,offset,1)==1);
    CPPUNIT_ASSERT(getRefRO(myData,offset,2)==2);

    // use a non-existent tag so we get a pointer to
    // the first element of the data array
    double* sampleData=myData.getSampleDataByTag(9);
    for (int i=0; i<myData.getLength(); i++) {
      if (i<3) {
        CPPUNIT_ASSERT(sampleData[i]==i);
      } else if ((i>=3) && (i<6)) {
        CPPUNIT_ASSERT(sampleData[i]==i-2);
      } else if ((i>=6) && (i<9)) {
        CPPUNIT_ASSERT(sampleData[i]==i-4);
      } else if ((i>=9) && (i<12)) {
        CPPUNIT_ASSERT(sampleData[i]==i-6);
      } else {
        CPPUNIT_ASSERT(sampleData[i]==i-12);
      }
    }

  }

  {

    cout << "\tTest adding one key with one value to DataTagged with three tags." << endl;

    DataTagged::TagListType keys;
    keys.push_back(1);
    keys.push_back(2);
    keys.push_back(3);

    DataTagged::FloatBatchType values;

    DataTypes::ShapeType viewShape;
    viewShape.push_back(3);

    // default value
    DataTypes::RealVectorType viewData(3*4);
    for (int i=0;i<viewShape[0];i++) {
      viewData[i]=i;
    }

    // value for tag "1"
    for (int i=0;i<viewShape[0];i++) {
	viewData[viewShape[0]+i]=i+1.0;
    }

    // value for tag "2"
    for (int i=0;i<viewShape[0];i++) {
	viewData[2*viewShape[0]+i]=i+2.0;
    }

    // value for tag "3"
    for (int i=0;i<viewShape[0];i++) {
	viewData[3*viewShape[0]+i]=i+3.0;
    }

    DataTagged myData(FunctionSpace(),viewShape,keys,viewData);

    keys.clear();
    keys.push_back(4);

    values.clear();
    // value for tag "4"
    for (int i=0;i<viewShape[0];i++) {
      values.push_back(i+4.0);
    }

    myData.addTaggedValues(keys,values,viewShape);

    CPPUNIT_ASSERT(myData.isCurrentTag(4));

    CPPUNIT_ASSERT(myData.getTagLookup().size()==4);

    CPPUNIT_ASSERT(myData.getLength()==15);

    CPPUNIT_ASSERT(myData.getRank()==1);
    CPPUNIT_ASSERT(myData.getNoValues()==3);
    CPPUNIT_ASSERT(myData.getShape().size()==1);


    int offset=myData.getOffsetForTag(4);
    CPPUNIT_ASSERT(offset==12);
    CPPUNIT_ASSERT(getRefRO(myData,offset,0)==4);
    CPPUNIT_ASSERT(getRefRO(myData,offset,1)==5);
    CPPUNIT_ASSERT(getRefRO(myData,offset,2)==6);

    // use a non-existent tag so we get a pointer to
    // the first element of the data array
    double* sampleData=myData.getSampleDataByTag(9);
    for (int i=0; i<myData.getLength(); i++) {
      if (i<3) {
        CPPUNIT_ASSERT(sampleData[i]==i);
      } else if ((i>=3) && (i<6)) {
        CPPUNIT_ASSERT(sampleData[i]==i-2);
      } else if ((i>=6) && (i<9)) {
        CPPUNIT_ASSERT(sampleData[i]==i-4);
      } else if ((i>=9) && (i<12)) {
        CPPUNIT_ASSERT(sampleData[i]==i-6);
      } else {
        CPPUNIT_ASSERT(sampleData[i]==i-8);
      }
    }

  }

  {

    cout << "\tTest adding three keys with one value to DataTagged with three tags." << endl;

    DataTagged::TagListType keys;
    keys.push_back(1);
    keys.push_back(2);
    keys.push_back(3);

    DataTagged::FloatBatchType values;

    DataTypes::ShapeType viewShape;
    viewShape.push_back(3);

    // default value
    DataTypes::RealVectorType viewData(3*4);
    for (int i=0;i<viewShape[0];i++) {
      viewData[i]=i;
    }

    // value for tag "1"
    for (int i=0;i<viewShape[0];i++) {
      viewData[viewShape[0]+i]=i+1.0;
    }

    // value for tag "2"
    for (int i=0;i<viewShape[0];i++) {
	viewData[2*viewShape[0]+i]=i+2.0;
    }

    // value for tag "3"
    for (int i=0;i<viewShape[0];i++) {
	viewData[3*viewShape[0]+i]=i+3.0;
    }

    DataTagged myData(FunctionSpace(),viewShape,keys,viewData);

    keys.clear();
    keys.push_back(4);
    keys.push_back(5);
    keys.push_back(6);

    values.clear();
    // value for tags "4", "5" and "6"
    for (int i=0;i<viewShape[0];i++) {
	values.push_back(i+4.0);
    }

    myData.addTaggedValues(keys,values,viewShape);

    CPPUNIT_ASSERT(myData.isCurrentTag(4));
    CPPUNIT_ASSERT(myData.isCurrentTag(5));
    CPPUNIT_ASSERT(myData.isCurrentTag(6));

    CPPUNIT_ASSERT(myData.getTagLookup().size()==6);

    CPPUNIT_ASSERT(myData.getLength()==21);

    CPPUNIT_ASSERT(myData.getRank()==1);
    CPPUNIT_ASSERT(myData.getNoValues()==3);
    CPPUNIT_ASSERT(myData.getShape().size()==1);

    int offset=myData.getOffsetForTag(4);
    CPPUNIT_ASSERT(offset==12);
    CPPUNIT_ASSERT(getRefRO(myData,offset,0)==4);
    CPPUNIT_ASSERT(getRefRO(myData,offset,1)==5);
    CPPUNIT_ASSERT(getRefRO(myData,offset,2)==6);

    offset=myData.getOffsetForTag(5);
    CPPUNIT_ASSERT(offset==15);
    CPPUNIT_ASSERT(getRefRO(myData,offset,0)==4);
    CPPUNIT_ASSERT(getRefRO(myData,offset,1)==5);
    CPPUNIT_ASSERT(getRefRO(myData,offset,2)==6);

    offset=myData.getOffsetForTag(6);
    CPPUNIT_ASSERT(offset==18);
    CPPUNIT_ASSERT(getRefRO(myData,offset,0)==4);
    CPPUNIT_ASSERT(getRefRO(myData,offset,1)==5);
    CPPUNIT_ASSERT(getRefRO(myData,offset,2)==6);

    // use a non-existent tag so we get a pointer to
    // the first element of the data array
    double* sampleData=myData.getSampleDataByTag(9);
    for (int i=0; i<myData.getLength(); i++) {
      if (i<3) {
        CPPUNIT_ASSERT(sampleData[i]==i);
      } else if ((i>=3) && (i<6)) {
        CPPUNIT_ASSERT(sampleData[i]==i-2);
      } else if ((i>=6) && (i<9)) {
        CPPUNIT_ASSERT(sampleData[i]==i-4);
      } else if ((i>=9) && (i<12)) {
        CPPUNIT_ASSERT(sampleData[i]==i-6);
      } else if ((i>=12) && (i<15)) {
        CPPUNIT_ASSERT(sampleData[i]==i-8);
      } else if ((i>=15) && (i<18)) {
        CPPUNIT_ASSERT(sampleData[i]==i-11);
      } else {
        CPPUNIT_ASSERT(sampleData[i]==i-14);
      }
    }

  }

  {

    cout << "\tTest adding three keys with three values to DataTagged with three tags." << endl;

    DataTagged::TagListType keys;
    keys.push_back(1);
    keys.push_back(2);
    keys.push_back(3);

    DataTagged::FloatBatchType values;

    DataTypes::ShapeType viewShape;
    viewShape.push_back(3);

    // default value
    DataTypes::RealVectorType viewData(3*4);
    for (int i=0;i<viewShape[0];i++) {
      viewData[i]=i;
    }

    // value for tag "1"
    for (int i=0;i<viewShape[0];i++) {
      viewData[viewShape[0]+i]=i+1.0;
    }

    // value for tag "2"
    for (int i=0;i<viewShape[0];i++) {
	viewData[2*viewShape[0]+i]=i+2.0;
    }

    // value for tag "3"
    for (int i=0;i<viewShape[0];i++) {
      viewData[3*viewShape[0]+i]=i+3.0;
    }

    DataTagged myData(FunctionSpace(),viewShape,keys,viewData);

    keys.clear();
    keys.push_back(4);
    keys.push_back(5);
    keys.push_back(6);

    values.clear();

    // value for tag "4"
    for (int i=0;i<viewShape[0];i++) {
      values.push_back(i+4.0);
    }

    // value for tag "5"
    for (int i=0;i<viewShape[0];i++) {
      values.push_back(i+5.0);
    }

    // value for tag "6"
    for (int i=0;i<viewShape[0];i++) {
	values.push_back(i+6.0);
    }

    myData.addTaggedValues(keys,values,viewShape);

    CPPUNIT_ASSERT(myData.isCurrentTag(4));
    CPPUNIT_ASSERT(myData.isCurrentTag(5));
    CPPUNIT_ASSERT(myData.isCurrentTag(6));

    CPPUNIT_ASSERT(myData.getTagLookup().size()==6);

    CPPUNIT_ASSERT(myData.getLength()==21);

    CPPUNIT_ASSERT(myData.getRank()==1);
    CPPUNIT_ASSERT(myData.getNoValues()==3);
    CPPUNIT_ASSERT(myData.getShape().size()==1);

    int offset=myData.getOffsetForTag(4);
    CPPUNIT_ASSERT(offset==12);
    CPPUNIT_ASSERT(getRefRO(myData,offset,0)==4);
    CPPUNIT_ASSERT(getRefRO(myData,offset,1)==5);
    CPPUNIT_ASSERT(getRefRO(myData,offset,2)==6);

    offset=myData.getOffsetForTag(5);
    CPPUNIT_ASSERT(offset==15);
    CPPUNIT_ASSERT(getRefRO(myData,offset,0)==5);
    CPPUNIT_ASSERT(getRefRO(myData,offset,1)==6);
    CPPUNIT_ASSERT(getRefRO(myData,offset,2)==7);

    offset=myData.getOffsetForTag(6);
    CPPUNIT_ASSERT(offset==18);
    CPPUNIT_ASSERT(getRefRO(myData,offset,0)==6);
    CPPUNIT_ASSERT(getRefRO(myData,offset,1)==7);
    CPPUNIT_ASSERT(getRefRO(myData,offset,2)==8);

    // use a non-existent tag so we get a pointer to
    // the first element of the data array
    double* sampleData=myData.getSampleDataByTag(9);
    for (int i=0; i<myData.getLength(); i++) {
      if (i<3) {
        CPPUNIT_ASSERT(sampleData[i]==i);
      } else if ((i>=3) && (i<6)) {
        CPPUNIT_ASSERT(sampleData[i]==i-2);
      } else if ((i>=6) && (i<9)) {
        CPPUNIT_ASSERT(sampleData[i]==i-4);
      } else if ((i>=9) && (i<12)) {
        CPPUNIT_ASSERT(sampleData[i]==i-6);
      } else if ((i>=12) && (i<15)) {
        CPPUNIT_ASSERT(sampleData[i]==i-8);
      } else if ((i>=15) && (i<18)) {
        CPPUNIT_ASSERT(sampleData[i]==i-10);
      } else {
        CPPUNIT_ASSERT(sampleData[i]==i-12);
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

    DataTagged::FloatBatchType values;

    DataTypes::ShapeType viewShape;
    viewShape.push_back(3);

    // default value
    DataTypes::RealVectorType viewData(3*4);
    for (int i=0;i<viewShape[0];i++) {
      viewData[i]=i;
    }

    // value for tag "1"
    for (int i=0;i<viewShape[0];i++) {
      viewData[viewShape[0]+i]=i+1.0;
    }

    // value for tag "2"
    for (int i=0;i<viewShape[0];i++) {
      viewData[2*viewShape[0]+i]=i+2.0;
    }
    // value for tag "3"
    for (int i=0;i<viewShape[0];i++) {
	viewData[3*viewShape[0]+i]=i+3.0;
    }

    DataTagged myData(FunctionSpace(),viewShape,keys,viewData);

    // new value for tag "2"
    RealVectorType tmp(viewShape[0]);
    for (int i=0;i<viewShape[0];i++) {
      tmp[i]=i+5.0;
    }

    myData.setTaggedValue(2,viewShape,tmp);

    CPPUNIT_ASSERT(myData.isCurrentTag(2));

    CPPUNIT_ASSERT(myData.getTagLookup().size()==3);

    CPPUNIT_ASSERT(myData.getLength()==12);

    CPPUNIT_ASSERT(myData.getRank()==1);
    CPPUNIT_ASSERT(myData.getNoValues()==3);
    CPPUNIT_ASSERT(myData.getShape().size()==1);
    int offset=myData.getOffsetForTag(2);
    CPPUNIT_ASSERT(offset==6);
    CPPUNIT_ASSERT(getRefRO(myData,offset,0)==5);
    CPPUNIT_ASSERT(getRefRO(myData,offset,1)==6);
    CPPUNIT_ASSERT(getRefRO(myData,offset,2)==7);

    // use a non-existent tag so we get a pointer to
    // the first element of the data array
    double* sampleData=myData.getSampleDataByTag(9);
    for (int i=0; i<myData.getLength(); i++) {
      if (i<3) {
        CPPUNIT_ASSERT(sampleData[i]==i);
      } else if ((i>=3) && (i<6)) {
        CPPUNIT_ASSERT(sampleData[i]==i-2);
      } else if ((i>=6) && (i<9)) {
        CPPUNIT_ASSERT(sampleData[i]==i-1);
      } else {
        CPPUNIT_ASSERT(sampleData[i]==i-6);
      }
    }

  }

}

void DataTaggedTestCase::testAll() {

  cout << endl;

  {

    cout << "\tTest default DataTagged." << endl;
    DataTagged myData=makeTagged();


    CPPUNIT_ASSERT(myData.getNumSamples()==1);
    CPPUNIT_ASSERT(myData.getNumDPPSample()==1);

    CPPUNIT_ASSERT(myData.validSamplePointNo(0));
    CPPUNIT_ASSERT(myData.validSampleNo(0));
    CPPUNIT_ASSERT(!myData.validSamplePointNo(1));
    CPPUNIT_ASSERT(!myData.validSampleNo(1));

    // data-point 0 has tag number 1 by default
    CPPUNIT_ASSERT(myData.getTagNumber(0)==1);

    CPPUNIT_ASSERT(!myData.isCurrentTag(1));

    CPPUNIT_ASSERT(myData.getTagLookup().size()==0);

    int ll=myData.getLength();
    cout << "\t" << ll << endl;
    
    
    CPPUNIT_ASSERT(myData.getLength()==1);

    CPPUNIT_ASSERT(myData.getPointOffset(0,0)==0);

    CPPUNIT_ASSERT(myData.getRank()==0);
    CPPUNIT_ASSERT(myData.getNoValues()==1);
    CPPUNIT_ASSERT(myData.getShape().size()==0);


    int offset=myData.getPointOffset(0,0);
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(myData.getDataAtOffsetRO(offset)==0.0);

    // Test non-existent tag returns the default value.
    offset=myData.getOffsetForTag(1);
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(myData.getDataAtOffsetRO(offset)==0.0);

    offset=myData.getDefaultOffset();
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(myData.getDataAtOffsetRO(offset)==0.0);

    // use a non-existent tag so we get a pointer to
    // the first element of the data array
    const double* sampleData=myData.getSampleDataByTag(9);
    for (int i=0; i<myData.getLength(); i++) {
      CPPUNIT_ASSERT(sampleData[i]==i);
    }
    sampleData=myData.getSampleDataRO(0);
    for (unsigned int i=0; i<myData.getNoValues(); i++) {
      CPPUNIT_ASSERT(sampleData[i]==i);
    }

  }

  {

    cout << "\tTest DataTagged with default value only." << endl;

    DataTypes::ShapeType viewShape;
    viewShape.push_back(3);

    DataTagged::TagListType keys;

    DataTagged::FloatBatchType values;

    DataTypes::RealVectorType viewData(3);
    for (int i=0;i<viewShape[0];i++) {
      viewData[i]=i;
    }
    DataTagged myData(FunctionSpace(),viewShape, viewData);

    CPPUNIT_ASSERT(myData.getNumSamples()==1);
    CPPUNIT_ASSERT(myData.getNumDPPSample()==1);

    CPPUNIT_ASSERT(myData.validSamplePointNo(0));
    CPPUNIT_ASSERT(myData.validSampleNo(0));
    CPPUNIT_ASSERT(!myData.validSamplePointNo(1));
    CPPUNIT_ASSERT(!myData.validSampleNo(1));

    // data-point 0 has tag number 1 by default
    CPPUNIT_ASSERT(myData.getTagNumber(0)==1);

    CPPUNIT_ASSERT(!myData.isCurrentTag(1));

    CPPUNIT_ASSERT(myData.getTagLookup().size()==0);

    CPPUNIT_ASSERT(myData.getLength()==3);

    CPPUNIT_ASSERT(myData.getPointOffset(0,0)==0);

    CPPUNIT_ASSERT(myData.getRank()==1);
    CPPUNIT_ASSERT(myData.getNoValues()==3);
    CPPUNIT_ASSERT(myData.getShape().size()==1);

    int offset=myData.getPointOffset(0,0);
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(getRefRO(myData,offset,0)==0);
    CPPUNIT_ASSERT(getRefRO(myData,offset,1)==1);
    CPPUNIT_ASSERT(getRefRO(myData,offset,2)==2);

    // Test non-existent tag returns the default value.
    offset=myData.getOffsetForTag(1);
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(getRefRO(myData,offset,0)==0);
    CPPUNIT_ASSERT(getRefRO(myData,offset,1)==1);
    CPPUNIT_ASSERT(getRefRO(myData,offset,2)==2);

    offset=myData.getDefaultOffset();
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(getRefRO(myData,offset,0)==0);
    CPPUNIT_ASSERT(getRefRO(myData,offset,1)==1);
    CPPUNIT_ASSERT(getRefRO(myData,offset,2)==2);

    // use a non-existent tag so we get a pointer to
    // the first element of the data array
    double* sampleData=myData.getSampleDataByTag(9);
    for (int i=0; i<myData.getLength(); i++) {
      CPPUNIT_ASSERT(sampleData[i]==i);
    }
    sampleData=myData.getSampleDataByTag(0);
    for (unsigned int i=0; i<myData.getNoValues(); i++) {
      CPPUNIT_ASSERT(sampleData[i]==i);
    }

  }

  {

    cout << "\tTest DataTagged with one tag." << endl;

    // the one data-point has tag value "1"

    DataTagged::TagListType keys;
    keys.push_back(1);

    DataTagged::FloatBatchType values;

    DataTypes::ShapeType viewShape;
    viewShape.push_back(3);

    // default value
    DataTypes::RealVectorType viewData(3*2);
    for (int i=0;i<viewShape[0];i++) {
      viewData[i]=i;
    }

    // value for tag "1"
    for (int i=0;i<viewShape[0];i++) {
	viewData[viewShape[0]+i]=i+1.0;
    }
    DataTagged myData(FunctionSpace(),viewShape,keys,viewData);

    CPPUNIT_ASSERT(myData.getNumSamples()==1);
    CPPUNIT_ASSERT(myData.getNumDPPSample()==1);

    CPPUNIT_ASSERT(myData.validSamplePointNo(0));
    CPPUNIT_ASSERT(myData.validSampleNo(0));
    CPPUNIT_ASSERT(!myData.validSamplePointNo(1));
    CPPUNIT_ASSERT(!myData.validSampleNo(1));

    // data-point 0 has tag number 1 by default
    CPPUNIT_ASSERT(myData.getTagNumber(0)==1);

    CPPUNIT_ASSERT(!myData.isCurrentTag(0));
    CPPUNIT_ASSERT(myData.isCurrentTag(1));

    CPPUNIT_ASSERT(myData.getTagLookup().size()==1);

    CPPUNIT_ASSERT(myData.getLength()==6);

    CPPUNIT_ASSERT(myData.getPointOffset(0,0)==3);

    CPPUNIT_ASSERT(myData.getRank()==1);
    CPPUNIT_ASSERT(myData.getNoValues()==3);
    CPPUNIT_ASSERT(myData.getShape().size()==1);


    int offset=myData.getPointOffset(0,0);
    CPPUNIT_ASSERT(offset==3);
    CPPUNIT_ASSERT(getRefRO(myData,offset,0)==1);
    CPPUNIT_ASSERT(getRefRO(myData,offset,1)==2);
    CPPUNIT_ASSERT(getRefRO(myData,offset,2)==3);

    offset=myData.getOffsetForTag(1);
    CPPUNIT_ASSERT(offset==3);
    CPPUNIT_ASSERT(getRefRO(myData,offset,0)==1);
    CPPUNIT_ASSERT(getRefRO(myData,offset,1)==2);
    CPPUNIT_ASSERT(getRefRO(myData,offset,2)==3);

    // Test non-existent tag returns the default value.
    offset=myData.getOffsetForTag(9);
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(getRefRO(myData,offset,0)==0);
    CPPUNIT_ASSERT(getRefRO(myData,offset,1)==1);
    CPPUNIT_ASSERT(getRefRO(myData,offset,2)==2);

    offset=myData.getDefaultOffset();
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(getRefRO(myData,offset,0)==0);
    CPPUNIT_ASSERT(getRefRO(myData,offset,1)==1);
    CPPUNIT_ASSERT(getRefRO(myData,offset,2)==2);

    // use a non-existent tag so we get a pointer to
    // the first element of the data array
    const double* sampleData=myData.getSampleDataByTag(9);
    for (int i=0; i<myData.getLength(); i++) {
      if (i<3) {
        CPPUNIT_ASSERT(sampleData[i]==i);
      } else {
        CPPUNIT_ASSERT(sampleData[i]==i-2);
      }
    }
    sampleData=myData.getSampleDataRO(0);
    for (unsigned int i=0; i<myData.getNoValues(); i++) {
      CPPUNIT_ASSERT(sampleData[i]==i+1);
    }

  }

  {

    cout << "\tTest DataTagged with multiple tags." << endl;

    // the one data-point has tag value "1"

    DataTagged::TagListType keys;
    keys.push_back(1);
    keys.push_back(2);
    keys.push_back(3);

    DataTagged::FloatBatchType values;

    DataTypes::ShapeType viewShape;
    viewShape.push_back(3);

    // default value
    DataTypes::RealVectorType viewData(3*4);
    for (int i=0;i<viewShape[0];i++) {
      viewData[i]=i;
    }

    // value for tag "1"
    for (int i=0;i<viewShape[0];i++) {
      viewData[viewShape[0]+i]=i+1.0;
    }

    // value for tag "2"
    for (int i=0;i<viewShape[0];i++) {
      viewData[2*viewShape[0]+i]=i+2.0;
    }

    // value for tag "3"
    for (int i=0;i<viewShape[0];i++) {
      viewData[3*viewShape[0]+i]=i+3.0;
    }

    DataTagged myData(FunctionSpace(),viewShape,keys,viewData);

    CPPUNIT_ASSERT(myData.getNumSamples()==1);
    CPPUNIT_ASSERT(myData.getNumDPPSample()==1);

    CPPUNIT_ASSERT(myData.validSamplePointNo(0));
    CPPUNIT_ASSERT(myData.validSampleNo(0));
    CPPUNIT_ASSERT(!myData.validSamplePointNo(1));
    CPPUNIT_ASSERT(!myData.validSampleNo(1));

    // data-point 0 has tag number 1 by default
    CPPUNIT_ASSERT(myData.getTagNumber(0)==1);

    CPPUNIT_ASSERT(!myData.isCurrentTag(0));
    CPPUNIT_ASSERT(myData.isCurrentTag(1));
    CPPUNIT_ASSERT(myData.isCurrentTag(2));
    CPPUNIT_ASSERT(myData.isCurrentTag(3));

    CPPUNIT_ASSERT(myData.getTagLookup().size()==3);

    CPPUNIT_ASSERT(myData.getLength()==12);

    CPPUNIT_ASSERT(myData.getPointOffset(0,0)==3);

    CPPUNIT_ASSERT(myData.getRank()==1);
    CPPUNIT_ASSERT(myData.getNoValues()==3);
    CPPUNIT_ASSERT(myData.getShape().size()==1);


    int offset=myData.getPointOffset(0,0);
    CPPUNIT_ASSERT(offset==3);
    CPPUNIT_ASSERT(getRefRO(myData,offset,0)==1);
    CPPUNIT_ASSERT(getRefRO(myData,offset,1)==2);
    CPPUNIT_ASSERT(getRefRO(myData,offset,2)==3);

    offset=myData.getOffsetForTag(1);
    CPPUNIT_ASSERT(offset==3);
    CPPUNIT_ASSERT(getRefRO(myData,offset,0)==1);
    CPPUNIT_ASSERT(getRefRO(myData,offset,1)==2);
    CPPUNIT_ASSERT(getRefRO(myData,offset,2)==3);

    // Test non-existent tag returns the default value.
    offset=myData.getOffsetForTag(0);
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(getRefRO(myData,offset,0)==0);
    CPPUNIT_ASSERT(getRefRO(myData,offset,1)==1);
    CPPUNIT_ASSERT(getRefRO(myData,offset,2)==2);

    offset=myData.getDefaultOffset();
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(getRefRO(myData,offset,0)==0);
    CPPUNIT_ASSERT(getRefRO(myData,offset,1)==1);
    CPPUNIT_ASSERT(getRefRO(myData,offset,2)==2);

    // Test data-points held for remaining tags
    offset=myData.getOffsetForTag(2);
    CPPUNIT_ASSERT(offset==6);
    CPPUNIT_ASSERT(getRefRO(myData,offset,0)==2);
    CPPUNIT_ASSERT(getRefRO(myData,offset,1)==3);
    CPPUNIT_ASSERT(getRefRO(myData,offset,2)==4);

    offset=myData.getOffsetForTag(3);
    CPPUNIT_ASSERT(offset==9);
    CPPUNIT_ASSERT(getRefRO(myData,offset,0)==3);
    CPPUNIT_ASSERT(getRefRO(myData,offset,1)==4);
    CPPUNIT_ASSERT(getRefRO(myData,offset,2)==5);

    // use a non-existent tag so we get a pointer to
    // the first element of the data array
    const double* sampleData=myData.getSampleDataByTag(9);
    for (int i=0; i<myData.getLength(); i++) {
      if (i<3) {
        CPPUNIT_ASSERT(sampleData[i]==i);
      } else if ((i>=3) && (i<6)) {
        CPPUNIT_ASSERT(sampleData[i]==i-2);
      } else if ((i>=6) && (i<9)) {
        CPPUNIT_ASSERT(sampleData[i]==i-4);
      } else {
        CPPUNIT_ASSERT(sampleData[i]==i-6);
      }
    }
    sampleData=myData.getSampleDataRO(0);
    for (unsigned int i=0; i<myData.getNoValues(); i++) {
      CPPUNIT_ASSERT(sampleData[i]==i+1);
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

    DataTagged::FloatBatchType values;

    DataTypes::ShapeType viewShape;
    viewShape.push_back(3);

    // default value
    DataTypes::RealVectorType viewData(3*4);
    for (int i=0;i<viewShape[0];i++) {
      viewData[i]=i;
    }

    // value for tag "1"
    for (int i=0;i<viewShape[0];i++) {
      viewData[viewShape[0]+i]=i+1.0;
    }

    // value for tag "2"
    for (int i=0;i<viewShape[0];i++) {
	viewData[2*viewShape[0]+i]=i+2.0;
    }

    // value for tag "3"
    for (int i=0;i<viewShape[0];i++) {
	viewData[3*viewShape[0]+i]=i+3.0;
    }

    DataTagged myData(FunctionSpace(),viewShape,keys,viewData);

    DataTagged myDataCopy(myData);

    CPPUNIT_ASSERT(myDataCopy.getNumSamples()==1);
    CPPUNIT_ASSERT(myDataCopy.getNumDPPSample()==1);

    CPPUNIT_ASSERT(myDataCopy.validSamplePointNo(0));
    CPPUNIT_ASSERT(myDataCopy.validSampleNo(0));
    CPPUNIT_ASSERT(!myDataCopy.validSamplePointNo(1));
    CPPUNIT_ASSERT(!myDataCopy.validSampleNo(1));

    // data-point 0 has tag number 1 by default
    CPPUNIT_ASSERT(myDataCopy.getTagNumber(0)==1);

    CPPUNIT_ASSERT(!myDataCopy.isCurrentTag(0));
    CPPUNIT_ASSERT(myDataCopy.isCurrentTag(1));
    CPPUNIT_ASSERT(myDataCopy.isCurrentTag(2));
    CPPUNIT_ASSERT(myDataCopy.isCurrentTag(3));

    CPPUNIT_ASSERT(myDataCopy.getTagLookup().size()==3);

    CPPUNIT_ASSERT(myDataCopy.getLength()==12);

    CPPUNIT_ASSERT(myDataCopy.getPointOffset(0,0)==3);

    CPPUNIT_ASSERT(myDataCopy.getRank()==1);
    CPPUNIT_ASSERT(myDataCopy.getNoValues()==3);
    CPPUNIT_ASSERT(myDataCopy.getShape().size()==1);

    int offset=myDataCopy.getPointOffset(0,0);
    CPPUNIT_ASSERT(offset==3);
    CPPUNIT_ASSERT(getRefRO(myDataCopy,offset,0)==1);
    CPPUNIT_ASSERT(getRefRO(myDataCopy,offset,1)==2);
    CPPUNIT_ASSERT(getRefRO(myDataCopy,offset,2)==3);

    offset=myDataCopy.getOffsetForTag(1);
    CPPUNIT_ASSERT(offset==3);
    CPPUNIT_ASSERT(getRefRO(myDataCopy,offset,0)==1);
    CPPUNIT_ASSERT(getRefRO(myDataCopy,offset,1)==2);
    CPPUNIT_ASSERT(getRefRO(myDataCopy,offset,2)==3);

    // Test non-existent tag returns the default value.
    offset=myDataCopy.getOffsetForTag(0);
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(getRefRO(myDataCopy,offset,0)==0);
    CPPUNIT_ASSERT(getRefRO(myDataCopy,offset,1)==1);
    CPPUNIT_ASSERT(getRefRO(myDataCopy,offset,2)==2);

    offset=myDataCopy.getDefaultOffset();
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(getRefRO(myDataCopy,offset,0)==0);
    CPPUNIT_ASSERT(getRefRO(myDataCopy,offset,1)==1);
    CPPUNIT_ASSERT(getRefRO(myDataCopy,offset,2)==2);

    // Test data-points held for remaining tags
    offset=myDataCopy.getOffsetForTag(2);
    CPPUNIT_ASSERT(offset==6);
    CPPUNIT_ASSERT(getRefRO(myDataCopy,offset,0)==2);
    CPPUNIT_ASSERT(getRefRO(myDataCopy,offset,1)==3);
    CPPUNIT_ASSERT(getRefRO(myDataCopy,offset,2)==4);

    offset=myDataCopy.getOffsetForTag(3);
    CPPUNIT_ASSERT(offset==9);
    CPPUNIT_ASSERT(getRefRO(myDataCopy,offset,0)==3);
    CPPUNIT_ASSERT(getRefRO(myDataCopy,offset,1)==4);
    CPPUNIT_ASSERT(getRefRO(myDataCopy,offset,2)==5);

    // use a non-existent tag so we get a pointer to
    // the first element of the data array
    double* sampleData=myDataCopy.getSampleDataByTag(9);
    for (int i=0; i<myData.getLength(); i++) {
      if (i<3) {
        CPPUNIT_ASSERT(sampleData[i]==i);
      } else if ((i>=3) && (i<6)) {
        CPPUNIT_ASSERT(sampleData[i]==i-2);
      } else if ((i>=6) && (i<9)) {
        CPPUNIT_ASSERT(sampleData[i]==i-4);
      } else {
        CPPUNIT_ASSERT(sampleData[i]==i-6);
      }
    }

  }

  {

    cout << "\tTest DataTagged copy constructor for DataConstant." << endl;

    // Create a DataConstant
    DataTypes::ShapeType shape;
    DataTypes::RealVectorType data(DataTypes::noValues(shape),0);
    data[0]=1.0;
    DataConstant myConstantData(FunctionSpace(),shape,data);

    // use this DataConstant to initialise a DataTagged
    DataTagged myData(myConstantData);

    //cout << myData.toString() << endl;

    CPPUNIT_ASSERT(myData.getNumSamples()==1);
    CPPUNIT_ASSERT(myData.getNumDPPSample()==1);

    CPPUNIT_ASSERT(myData.validSamplePointNo(0));
    CPPUNIT_ASSERT(myData.validSampleNo(0));
    CPPUNIT_ASSERT(!myData.validSamplePointNo(1));
    CPPUNIT_ASSERT(!myData.validSampleNo(1));

    // data-point 0 has tag number 1 by default
    CPPUNIT_ASSERT(myData.getTagNumber(0)==1);

    CPPUNIT_ASSERT(!myData.isCurrentTag(1));

    CPPUNIT_ASSERT(myData.getTagLookup().size()==0);

    CPPUNIT_ASSERT(myData.getLength()==1);

    CPPUNIT_ASSERT(myData.getPointOffset(0,0)==0);

    CPPUNIT_ASSERT(myData.getRank()==0);
    CPPUNIT_ASSERT(myData.getNoValues()==1);
    CPPUNIT_ASSERT(myData.getShape().size()==0);


    int offset=myData.getPointOffset(0,0);
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(myData.getDataAtOffsetRO(offset)==1.0);

    // Test non-existent tag returns the default value.
    offset=myData.getOffsetForTag(1);
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(myData.getDataAtOffsetRO(offset)==1.0);

    offset=myData.getDefaultOffset();
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(myData.getDataAtOffsetRO(offset)==1.0);

    // use a non-existent tag so we get a pointer to
    // the first element of the data array
    double* sampleData=myData.getSampleDataByTag(9);
    for (int i=0; i<myData.getLength(); i++) {
      CPPUNIT_ASSERT(sampleData[i]==i+1);
    }

  }

}

void DataTaggedTestCase::testGetSlice() {

  cout << endl;

  {

    cout << "\tTest slicing default DataTagged." << endl;

    DataTagged myData=makeTagged();

    DataTypes::RegionType region;

    DataAbstract* slicedDefault = myData.getSlice(region);

    const DataTagged* myDataSliced=dynamic_cast<const DataTagged*>(slicedDefault);

    CPPUNIT_ASSERT(myDataSliced->getTagLookup().size()==0);

    CPPUNIT_ASSERT(myDataSliced->getLength()==1);

    int offset=myDataSliced->getDefaultOffset();
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(myDataSliced->getRank()==0);
    CPPUNIT_ASSERT(myDataSliced->getNoValues()==1);
    CPPUNIT_ASSERT(myDataSliced->getShape().size()==0);
    CPPUNIT_ASSERT(myDataSliced->getVectorRO()[0]==0.0);

    delete slicedDefault;
  }

  {

    cout << "\tTest slicing DataTagged with rank 1 default value only." << endl;

    DataTypes::ShapeType viewShape;
    viewShape.push_back(3);

    DataTagged::TagListType keys;

    DataTagged::FloatBatchType values;

    DataTypes::RealVectorType viewData(3);
    for (int i=0;i<viewShape[0];i++) {
      viewData[i]=i;
    }

    DataTagged myData(FunctionSpace(),viewShape,viewData);

    // full slice

    std::pair<int, int> region_element;
    region_element.first=0;
    region_element.second=3;
    DataTypes::RegionType region;
    region.push_back(region_element);

    DataAbstract* slicedDefault = myData.getSlice(region);

    const DataTagged* myDataSliced=dynamic_cast<const DataTagged*>(slicedDefault);

    CPPUNIT_ASSERT(myDataSliced->getTagLookup().size()==0);

    CPPUNIT_ASSERT(myDataSliced->getLength()==3);

    CPPUNIT_ASSERT(myDataSliced->getRank()==1);
    CPPUNIT_ASSERT(myDataSliced->getNoValues()==3);
    CPPUNIT_ASSERT(myDataSliced->getShape().size()==1);

    int offset=myDataSliced->getDefaultOffset();
    CPPUNIT_ASSERT(offset==0);

    CPPUNIT_ASSERT(getRefRO(*myDataSliced,offset,0)==0.0);
    CPPUNIT_ASSERT(getRefRO(*myDataSliced,offset,1)==1.0);
    CPPUNIT_ASSERT(getRefRO(*myDataSliced,offset,2)==2.0);

    // scalar slice

    region.clear();
    region_element.first=0;
    region_element.second=0;
    region.push_back(region_element);

    delete slicedDefault;

    slicedDefault = myData.getSlice(region);


    myDataSliced=dynamic_cast<const DataTagged*>(slicedDefault);

    CPPUNIT_ASSERT(myDataSliced->getTagLookup().size()==0);

    CPPUNIT_ASSERT(myDataSliced->getLength()==1);

    offset=myDataSliced->getDefaultOffset();
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(myDataSliced->getRank()==0);
    CPPUNIT_ASSERT(myDataSliced->getNoValues()==1);
    CPPUNIT_ASSERT(myDataSliced->getShape().size()==0);
    CPPUNIT_ASSERT(myDataSliced->getVectorRO()[0]==0.0);

    delete slicedDefault;
  }

  {

    cout << "\tTest slicing DataTagged with rank 3 default value only." << endl;

    DataTypes::ShapeType viewShape;
    viewShape.push_back(3);
    viewShape.push_back(3);
    viewShape.push_back(3);

    DataTagged::TagListType keys;

    DataTagged::FloatBatchType values;

    DataTypes::RealVectorType viewData(27);
    for (int i=0;i<viewData.size();i++) {
      viewData[i]=i;
    }

    DataTagged myData(FunctionSpace(),viewShape,viewData);


    // full slice

    std::pair<int, int> region_element;
    region_element.first=0;
    region_element.second=3;
    DataTypes::RegionType region;
    region.push_back(region_element);
    region.push_back(region_element);
    region.push_back(region_element);

    DataAbstract* slicedDefault = myData.getSlice(region);

    const DataTagged* myDataSliced=dynamic_cast<const DataTagged*>(slicedDefault);

    CPPUNIT_ASSERT(myDataSliced->getTagLookup().size()==0);

    CPPUNIT_ASSERT(myDataSliced->getLength()==27);

    int offset=myDataSliced->getDefaultOffset();
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(myDataSliced->getRank()==3);
    CPPUNIT_ASSERT(myDataSliced->getNoValues()==27);
    CPPUNIT_ASSERT(myDataSliced->getShape().size()==3);

    // rank 1 slice

    region.clear();
    region.push_back(region_element);
    region_element.second=0;
    region.push_back(region_element);
    region.push_back(region_element);

    delete slicedDefault;

    slicedDefault = myData.getSlice(region);


    myDataSliced=dynamic_cast<const DataTagged*>(slicedDefault);

    CPPUNIT_ASSERT(myDataSliced->getTagLookup().size()==0);

    CPPUNIT_ASSERT(myDataSliced->getLength()==3);

    CPPUNIT_ASSERT(myDataSliced->getRank()==1);
    CPPUNIT_ASSERT(myDataSliced->getNoValues()==3);
    CPPUNIT_ASSERT(myDataSliced->getShape().size()==1);

    offset=myDataSliced->getDefaultOffset();
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(getRefRO(*myDataSliced,offset,0)==0.0);
    CPPUNIT_ASSERT(getRefRO(*myDataSliced,offset,1)==1.0);
    CPPUNIT_ASSERT(getRefRO(*myDataSliced,offset,2)==2.0);

    // scalar slice

    region.clear();
    region_element.first=2;
    region_element.second=2;
    region.push_back(region_element);
    region.push_back(region_element);
    region.push_back(region_element);

    delete slicedDefault;

    slicedDefault = myData.getSlice(region);


    myDataSliced=dynamic_cast<const DataTagged*>(slicedDefault);

    CPPUNIT_ASSERT(myDataSliced->getTagLookup().size()==0);

    CPPUNIT_ASSERT(myDataSliced->getLength()==1);

    CPPUNIT_ASSERT(myDataSliced->getRank()==0);
    CPPUNIT_ASSERT(myDataSliced->getNoValues()==1);
    CPPUNIT_ASSERT(myDataSliced->getShape().size()==0);


    offset=myDataSliced->getDefaultOffset();
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(myDataSliced->getVectorRO()[0]==26);
    delete slicedDefault;
  }

  {

    cout << "\tTest slicing DataTagged with scalar values and one tag." << endl;

    DataTagged::TagListType keys;
    keys.push_back(1);

    DataTagged::FloatBatchType values;

    DataTypes::ShapeType viewShape;

    // default value
    DataTypes::RealVectorType viewData(1*2);
    viewData[0]=0.0;

    // value for tag "1"
    viewData[1]=1.0;

    DataTagged myData(FunctionSpace(),viewShape,keys, viewData);


    // full slice

    DataTypes::RegionType region;

    DataAbstract* slicedDefault = myData.getSlice(region);


    const DataTagged* myDataSliced=dynamic_cast<const DataTagged*>(slicedDefault);

    CPPUNIT_ASSERT(myDataSliced->getTagLookup().size()==1);

    CPPUNIT_ASSERT(myDataSliced->getLength()==2);

    int offset=myDataSliced->getDefaultOffset();
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(myDataSliced->getRank()==0);
    CPPUNIT_ASSERT(myDataSliced->getNoValues()==1);
    CPPUNIT_ASSERT(myDataSliced->getShape().size()==0);
    CPPUNIT_ASSERT(myDataSliced->getVectorRO()[offset]==0);

    offset=myDataSliced->getOffsetForTag(1);
    CPPUNIT_ASSERT(offset==1);
    CPPUNIT_ASSERT(myDataSliced->getVectorRO()[offset]==1);

    delete slicedDefault;

  }

  {

    cout << "\tTest slicing DataTagged with rank 1 values and one tag." << endl;

    DataTypes::ShapeType viewShape;
    viewShape.push_back(3);

    DataTagged::TagListType keys;
    keys.push_back(1);

    DataTagged::FloatBatchType values;

    // default value
    DataTypes::RealVectorType viewData(3*2);
    for (int i=0;i<viewShape[0];i++) {
      viewData[i]=i;
    }

    // value for tag "1"
    for (int i=0;i<viewShape[0];i++) {
       viewData[viewShape[0]+i]=i+3.0;
    }

    DataTagged myData(FunctionSpace(),viewShape,keys,viewData);

    // full slice

    std::pair<int, int> region_element;
    region_element.first=0;
    region_element.second=3;
    DataTypes::RegionType region;
    region.push_back(region_element);

    DataAbstract* slicedDefault = myData.getSlice(region);

    const DataTagged* myDataSliced=dynamic_cast<const DataTagged*>(slicedDefault);

    CPPUNIT_ASSERT(myDataSliced->getTagLookup().size()==1);

    CPPUNIT_ASSERT(myDataSliced->getLength()==6);
    CPPUNIT_ASSERT(myDataSliced->getRank()==1);
    CPPUNIT_ASSERT(myDataSliced->getNoValues()==3);
    CPPUNIT_ASSERT(myDataSliced->getShape().size()==1);
    int offset=myDataSliced->getDefaultOffset();
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(getRefRO(*myDataSliced,offset,0)==0);
    CPPUNIT_ASSERT(getRefRO(*myDataSliced,offset,1)==1);
    CPPUNIT_ASSERT(getRefRO(*myDataSliced,offset,2)==2);

    offset=myDataSliced->getOffsetForTag(1);
    CPPUNIT_ASSERT(offset==3);
    CPPUNIT_ASSERT(getRefRO(*myDataSliced,offset,0)==3);
    CPPUNIT_ASSERT(getRefRO(*myDataSliced,offset,1)==4);
    CPPUNIT_ASSERT(getRefRO(*myDataSliced,offset,2)==5);

    // scalar slice

    region_element.first=1;
    region_element.second=1;
    region.clear();
    region.push_back(region_element);

    delete slicedDefault;

    slicedDefault = myData.getSlice(region);

    myDataSliced=dynamic_cast<const DataTagged*>(slicedDefault);

    CPPUNIT_ASSERT(myDataSliced->getTagLookup().size()==1);

    CPPUNIT_ASSERT(myDataSliced->getLength()==2);

    CPPUNIT_ASSERT(myDataSliced->getRank()==0);
    CPPUNIT_ASSERT(myDataSliced->getNoValues()==1);
    CPPUNIT_ASSERT(myDataSliced->getShape().size()==0);


    offset=myDataSliced->getDefaultOffset();
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(myDataSliced->getVectorRO()[offset]==1);

    offset=myDataSliced->getOffsetForTag(1);
    CPPUNIT_ASSERT(offset==1);
    CPPUNIT_ASSERT(myDataSliced->getVectorRO()[offset]==4);
    delete slicedDefault;

  }

  {

    cout << "\tTest slicing DataTagged with rank 3 values and one tag." << endl;
    DataTypes::ShapeType viewShape;
    viewShape.push_back(3);
    viewShape.push_back(3);
    viewShape.push_back(3);

    DataTagged::TagListType keys;
    keys.push_back(1);

    DataTagged::FloatBatchType values;

    // default value
    DataTypes::RealVectorType viewData(27*2);
    for (int i=0;i<noValues(viewShape);i++) {
      viewData[i]=i;
    }

    // value for tag "1"
    for (int i=0;i<noValues(viewShape);i++) {
      viewData[noValues(viewShape)+i]=i+27.0;
    }

    DataTagged myData(FunctionSpace(),viewShape,keys,viewData);

    // full slice

    std::pair<int, int> region_element;
    region_element.first=0;
    region_element.second=3;
    DataTypes::RegionType region;
    region.push_back(region_element);
    region.push_back(region_element);
    region.push_back(region_element);

    DataAbstract* slicedDefault = myData.getSlice(region);


    const DataTagged* myDataSliced=dynamic_cast<const DataTagged*>(slicedDefault);

    CPPUNIT_ASSERT(myDataSliced->getTagLookup().size()==1);

    CPPUNIT_ASSERT(myDataSliced->getLength()==54);

    CPPUNIT_ASSERT(myDataSliced->getRank()==3);
    CPPUNIT_ASSERT(myDataSliced->getNoValues()==27);
    CPPUNIT_ASSERT(myDataSliced->getShape().size()==3);

    int offset=myDataSliced->getDefaultOffset();
    CPPUNIT_ASSERT(offset==0);

    offset=myDataSliced->getOffsetForTag(1);
    CPPUNIT_ASSERT(offset==27);

    // rank 1 slice

    region.clear();
    region.push_back(region_element);
    region_element.second=0;
    region.push_back(region_element);
    region.push_back(region_element);

    delete slicedDefault;

    slicedDefault = myData.getSlice(region);

    myDataSliced=dynamic_cast<const DataTagged*>(slicedDefault);

    CPPUNIT_ASSERT(myDataSliced->getTagLookup().size()==1);

    CPPUNIT_ASSERT(myDataSliced->getLength()==6);

    CPPUNIT_ASSERT(myDataSliced->getRank()==1);
    CPPUNIT_ASSERT(myDataSliced->getNoValues()==3);
    CPPUNIT_ASSERT(myDataSliced->getShape().size()==1);

    offset=myDataSliced->getDefaultOffset();
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(getRefRO(*myDataSliced,offset,0)==0);
    CPPUNIT_ASSERT(getRefRO(*myDataSliced,offset,1)==1);
    CPPUNIT_ASSERT(getRefRO(*myDataSliced,offset,2)==2);

    offset=myDataSliced->getOffsetForTag(1);
    CPPUNIT_ASSERT(offset==3);
    CPPUNIT_ASSERT(getRefRO(*myDataSliced,offset,0)==27);
    CPPUNIT_ASSERT(getRefRO(*myDataSliced,offset,1)==28);
    CPPUNIT_ASSERT(getRefRO(*myDataSliced,offset,2)==29);
    // scalar slice

    region_element.first=1;
    region_element.second=1;
    region.clear();
    region.push_back(region_element);
    region.push_back(region_element);
    region.push_back(region_element);

    delete slicedDefault;

    slicedDefault = myData.getSlice(region);


    myDataSliced=dynamic_cast<const DataTagged*>(slicedDefault);

    CPPUNIT_ASSERT(myDataSliced->getTagLookup().size()==1);

    CPPUNIT_ASSERT(myDataSliced->getLength()==2);

    CPPUNIT_ASSERT(myDataSliced->getRank()==0);
    CPPUNIT_ASSERT(myDataSliced->getNoValues()==1);
    CPPUNIT_ASSERT(myDataSliced->getShape().size()==0);

    offset=myDataSliced->getDefaultOffset();
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(myDataSliced->getVectorRO()[offset]==13);

    offset=myDataSliced->getOffsetForTag(1);
    CPPUNIT_ASSERT(offset==1);
    CPPUNIT_ASSERT(myDataSliced->getVectorRO()[offset]==40);

    delete slicedDefault;
  }

  {

    cout << "\tTest slicing DataTagged with scalar values and three tags." << endl;

    DataTagged::TagListType keys;
    keys.push_back(1);
    keys.push_back(2);
    keys.push_back(3);

    DataTagged::FloatBatchType values;

    DataTypes::ShapeType viewShape;

    // default value
    DataTypes::RealVectorType viewData(1*4);
    viewData[0]=0.0;

    // value for tag "1"
    viewData[1]=1.0;

    // value for tag "2"
    viewData[2]=2.0;

    // value for tag "3"
    viewData[3]=3.0;

    DataTagged myData(FunctionSpace(),viewShape,keys,viewData);

    // full slice

    DataTypes::RegionType region;

    DataAbstract* slicedDefault = myData.getSlice(region);

    const DataTagged* myDataSliced=dynamic_cast<const DataTagged*>(slicedDefault);

    CPPUNIT_ASSERT(myDataSliced->getTagLookup().size()==3);

    CPPUNIT_ASSERT(myDataSliced->getLength()==4);

    CPPUNIT_ASSERT(myDataSliced->getRank()==0);
    CPPUNIT_ASSERT(myDataSliced->getNoValues()==1);
    CPPUNIT_ASSERT(myDataSliced->getShape().size()==0);

    int offset=myDataSliced->getDefaultOffset();
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(myDataSliced->getVectorRO()[offset]==0);

    offset=myDataSliced->getOffsetForTag(1);
    CPPUNIT_ASSERT(offset==1);
    CPPUNIT_ASSERT(myDataSliced->getVectorRO()[offset]==1);

    offset=myDataSliced->getOffsetForTag(2);
    CPPUNIT_ASSERT(offset==2);
    CPPUNIT_ASSERT(myDataSliced->getVectorRO()[offset]==2);

    offset=myDataSliced->getOffsetForTag(3);
    CPPUNIT_ASSERT(offset==3);
    CPPUNIT_ASSERT(myDataSliced->getVectorRO()[offset]==3);

    delete slicedDefault;
  }

  {

    cout << "\tTest slicing DataTagged with rank 1 values and three tags." << endl;

    DataTypes::ShapeType viewShape;
    viewShape.push_back(3);

    DataTagged::TagListType keys;
    keys.push_back(1);
    keys.push_back(2);
    keys.push_back(3);

    DataTagged::FloatBatchType values;

    // default value
    DataTypes::RealVectorType viewData(3*4);
    for (int i=0;i<viewShape[0];i++) {
      viewData[i]=i;
    }

    // value for tag "1"
    for (int i=0;i<viewShape[0];i++) {
	viewData[viewShape[0]+i]=i+3.0;
    }

    // value for tag "2"
    for (int i=0;i<viewShape[0];i++) {
	viewData[2*viewShape[0]+i]=i+6.0;
    }

    // value for tag "3"
    for (int i=0;i<viewShape[0];i++) {
	viewData[3*viewShape[0]+i]=i+9.0;
    }

    DataTagged myData(FunctionSpace(),viewShape, keys, viewData);


    // full slice

    std::pair<int, int> region_element;
    region_element.first=0;
    region_element.second=3;
    DataTypes::RegionType region;
    region.push_back(region_element);

    DataAbstract* slicedDefault = myData.getSlice(region);

    const DataTagged* myDataSliced=dynamic_cast<const DataTagged*>(slicedDefault);

    CPPUNIT_ASSERT(myDataSliced->getTagLookup().size()==3);

    CPPUNIT_ASSERT(myDataSliced->getLength()==12);

    CPPUNIT_ASSERT(myDataSliced->getRank()==1);
    CPPUNIT_ASSERT(myDataSliced->getNoValues()==3);
    CPPUNIT_ASSERT(myDataSliced->getShape().size()==1);


    int offset=myDataSliced->getDefaultOffset();
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(getRefRO(*myDataSliced,offset,0)==0);
    CPPUNIT_ASSERT(getRefRO(*myDataSliced,offset,1)==1);
    CPPUNIT_ASSERT(getRefRO(*myDataSliced,offset,2)==2);

    offset=myDataSliced->getOffsetForTag(1);
    CPPUNIT_ASSERT(offset==3);
    CPPUNIT_ASSERT(getRefRO(*myDataSliced,offset,0)==3);
    CPPUNIT_ASSERT(getRefRO(*myDataSliced,offset,1)==4);
    CPPUNIT_ASSERT(getRefRO(*myDataSliced,offset,2)==5);

    offset=myDataSliced->getOffsetForTag(2);
    CPPUNIT_ASSERT(offset==6);
    CPPUNIT_ASSERT(getRefRO(*myDataSliced,offset,0)==6);
    CPPUNIT_ASSERT(getRefRO(*myDataSliced,offset,1)==7);
    CPPUNIT_ASSERT(getRefRO(*myDataSliced,offset,2)==8);

    offset=myDataSliced->getOffsetForTag(3);
    CPPUNIT_ASSERT(offset==9);
    CPPUNIT_ASSERT(getRefRO(*myDataSliced,offset,0)==9);
    CPPUNIT_ASSERT(getRefRO(*myDataSliced,offset,1)==10);
    CPPUNIT_ASSERT(getRefRO(*myDataSliced,offset,2)==11);

    // scalar slice

    region.clear();
    region_element.first=1;
    region_element.second=1;
    region.push_back(region_element);

    delete slicedDefault;

    slicedDefault = myData.getSlice(region);


    myDataSliced=dynamic_cast<const DataTagged*>(slicedDefault);

    CPPUNIT_ASSERT(myDataSliced->getTagLookup().size()==3);

    CPPUNIT_ASSERT(myDataSliced->getLength()==4);

    CPPUNIT_ASSERT(myDataSliced->getRank()==0);
    CPPUNIT_ASSERT(myDataSliced->getNoValues()==1);
    CPPUNIT_ASSERT(myDataSliced->getShape().size()==0);

    offset=myDataSliced->getDefaultOffset();
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(myDataSliced->getVectorRO()[offset]==1);

    offset=myDataSliced->getOffsetForTag(1);
    CPPUNIT_ASSERT(offset==1);
    CPPUNIT_ASSERT(myDataSliced->getVectorRO()[offset]==4);

    offset=myDataSliced->getOffsetForTag(2);
    CPPUNIT_ASSERT(offset==2);
    CPPUNIT_ASSERT(myDataSliced->getVectorRO()[offset]==7);

    offset=myDataSliced->getOffsetForTag(3);
    CPPUNIT_ASSERT(offset==3);
    CPPUNIT_ASSERT(myDataSliced->getVectorRO()[offset]==10);
    
    delete slicedDefault;
  }

  {

    cout << "\tTest slicing DataTagged with rank 3 values and three tags." << endl;

    DataTypes::ShapeType viewShape;
    viewShape.push_back(3);
    viewShape.push_back(3);
    viewShape.push_back(3);

    DataTagged::TagListType keys;
    keys.push_back(1);
    keys.push_back(2);
    keys.push_back(3);

    DataTagged::FloatBatchType values;

    int nvals=27;
    // default value
    DataTypes::RealVectorType viewData(27*4);
    for (int i=0;i<nvals;i++) {
      viewData[i]=i;
    }

    // value for tag "1"
    for (int i=0;i<nvals;i++) {
      viewData[nvals+i]=i+27.0;
    }

    // value for tag "2"
    for (int i=0;i<nvals;i++) {
      viewData[2*nvals+i]=i+54.0;
    }

    // value for tag "3"
    for (int i=0;i<nvals;i++) {
      viewData[3*nvals+i]=i+81.0;
    }

    DataTagged myData(FunctionSpace(),viewShape,keys,viewData);


    // full slice

    std::pair<int, int> region_element;
    region_element.first=0;
    region_element.second=3;
    DataTypes::RegionType region;
    region.push_back(region_element);
    region.push_back(region_element);
    region.push_back(region_element);

    DataAbstract* slicedDefault = myData.getSlice(region);


    const DataTagged* myDataSliced=dynamic_cast<const DataTagged*>(slicedDefault);

    CPPUNIT_ASSERT(myDataSliced->getTagLookup().size()==3);

    CPPUNIT_ASSERT(myDataSliced->getLength()==108);

    CPPUNIT_ASSERT(myDataSliced->getRank()==3);
    CPPUNIT_ASSERT(myDataSliced->getNoValues()==27);
    CPPUNIT_ASSERT(myDataSliced->getShape().size()==3);

    int offset=myDataSliced->getDefaultOffset();
    CPPUNIT_ASSERT(offset==0);

    offset=myDataSliced->getOffsetForTag(1);
    CPPUNIT_ASSERT(offset==27);

    offset=myDataSliced->getOffsetForTag(2);
    CPPUNIT_ASSERT(offset==54);

    offset=myDataSliced->getOffsetForTag(3);
    CPPUNIT_ASSERT(offset==81);

    // rank 1 slice

    region.clear();
    region.push_back(region_element);
    region_element.second=0;
    region.push_back(region_element);
    region.push_back(region_element);

    delete slicedDefault;

    slicedDefault = myData.getSlice(region);


    myDataSliced=dynamic_cast<const DataTagged*>(slicedDefault);

    CPPUNIT_ASSERT(myDataSliced->getTagLookup().size()==3);

    CPPUNIT_ASSERT(myDataSliced->getLength()==12);

    CPPUNIT_ASSERT(myDataSliced->getRank()==1);
    CPPUNIT_ASSERT(myDataSliced->getNoValues()==3);
    CPPUNIT_ASSERT(myDataSliced->getShape().size()==1);

    offset=myDataSliced->getDefaultOffset();
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(getRefRO(*myDataSliced,offset,0)==0);
    CPPUNIT_ASSERT(getRefRO(*myDataSliced,offset,1)==1);
    CPPUNIT_ASSERT(getRefRO(*myDataSliced,offset,2)==2);

    offset=myDataSliced->getOffsetForTag(1);
    CPPUNIT_ASSERT(offset==3);
    CPPUNIT_ASSERT(getRefRO(*myDataSliced,offset,0)==27);
    CPPUNIT_ASSERT(getRefRO(*myDataSliced,offset,1)==28);
    CPPUNIT_ASSERT(getRefRO(*myDataSliced,offset,2)==29);

    offset=myDataSliced->getOffsetForTag(2);
    CPPUNIT_ASSERT(offset==6);
    CPPUNIT_ASSERT(getRefRO(*myDataSliced,offset,0)==54);
    CPPUNIT_ASSERT(getRefRO(*myDataSliced,offset,1)==55);
    CPPUNIT_ASSERT(getRefRO(*myDataSliced,offset,2)==56);

    offset=myDataSliced->getOffsetForTag(3);
    CPPUNIT_ASSERT(offset==9);
    CPPUNIT_ASSERT(getRefRO(*myDataSliced,offset,0)==81);
    CPPUNIT_ASSERT(getRefRO(*myDataSliced,offset,1)==82);
    CPPUNIT_ASSERT(getRefRO(*myDataSliced,offset,2)==83);

    // scalar slice

    region_element.first=1;
    region_element.second=1;
    region.clear();
    region.push_back(region_element);
    region.push_back(region_element);
    region.push_back(region_element);

    delete myDataSliced;

    slicedDefault = myData.getSlice(region);


    myDataSliced=dynamic_cast<const DataTagged*>(slicedDefault);

    CPPUNIT_ASSERT(myDataSliced->getTagLookup().size()==3);

    CPPUNIT_ASSERT(myDataSliced->getLength()==4);
    CPPUNIT_ASSERT(myDataSliced->getRank()==0);
    CPPUNIT_ASSERT(myDataSliced->getNoValues()==1);
    CPPUNIT_ASSERT(myDataSliced->getShape().size()==0);

    offset=myDataSliced->getDefaultOffset();
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(myDataSliced->getVectorRO()[offset]==13);

    offset=myDataSliced->getOffsetForTag(1);
    CPPUNIT_ASSERT(offset==1);
    CPPUNIT_ASSERT(myDataSliced->getVectorRO()[offset]==40);

    offset=myDataSliced->getOffsetForTag(2);
    CPPUNIT_ASSERT(offset==2);
    CPPUNIT_ASSERT(myDataSliced->getVectorRO()[offset]==67);

    offset=myDataSliced->getOffsetForTag(3);
    CPPUNIT_ASSERT(offset==3);
    CPPUNIT_ASSERT(myDataSliced->getVectorRO()[offset]==94);

    delete slicedDefault;
  }

}

void DataTaggedTestCase::testSetSlice() {

  cout << endl;

  {

    cout << "\tTest slicing default DataTagged." << endl;

    DataTagged myData1=makeTagged();
    DataTagged myData2=makeTagged();

    DataTypes::RegionType region;

    myData2.getDataAtOffsetRW(myData2.getDefaultOffset())=1.0;
    myData1.setSlice(&myData2, region);

    CPPUNIT_ASSERT(myData1.getTagLookup().size()==0);

    CPPUNIT_ASSERT(myData1.getLength()==1);
    CPPUNIT_ASSERT(myData1.getRank()==0);
    CPPUNIT_ASSERT(myData1.getNoValues()==1);
    CPPUNIT_ASSERT(myData1.getShape().size()==0);

    int offset=myData1.getDefaultOffset();
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(myData1.getDataAtOffsetRW(offset)==1.0);

  }

  {

    cout << "\tTest slicing DataTagged with rank 1 default value only." << endl;

    DataTagged::TagListType keys;

    DataTagged::FloatBatchType values;

    DataTypes::ShapeType viewShape;
    viewShape.push_back(3);

    DataTypes::RealVectorType viewData1(3);
    for (int i=0;i<viewShape[0];i++) {
      viewData1[i]=i;
    }
    DataTagged myData1(FunctionSpace(),viewShape,viewData1);

    DataTypes::RealVectorType viewData2(3);
    for (int i=0;i<viewShape[0];i++) {
      viewData2[i]=i+3;
    }
    DataTagged myData2(FunctionSpace(),viewShape,viewData2);

    // full slice

    std::pair<int, int> region_element;
    region_element.first=0;
    region_element.second=3;
    DataTypes::RegionType region;
    region.push_back(region_element);

    myData1.setSlice(&myData2, region);

    CPPUNIT_ASSERT(myData1.getTagLookup().size()==0);

    CPPUNIT_ASSERT(myData1.getLength()==3);

    CPPUNIT_ASSERT(myData1.getRank()==1);
    CPPUNIT_ASSERT(myData1.getNoValues()==3);
    CPPUNIT_ASSERT(myData1.getShape().size()==1);

    int offset=myData1.getDefaultOffset();
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,0)==3.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,1)==4.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,2)==5.0);

    // rank 1 slice

    viewShape.clear();
    viewShape.push_back(1);

    DataTypes::RealVectorType viewData3(1);
    viewData3[0]=6.0;
    DataTagged myData3(FunctionSpace(),viewShape,viewData3);

    region.clear();
    region_element.first=1;
    region_element.second=2;
    region.push_back(region_element);

    myData1.setSlice(&myData3, region);

    CPPUNIT_ASSERT(myData1.getTagLookup().size()==0);

    CPPUNIT_ASSERT(myData1.getLength()==3);
    CPPUNIT_ASSERT(myData1.getRank()==1);
    CPPUNIT_ASSERT(myData1.getNoValues()==3);
    CPPUNIT_ASSERT(myData1.getShape().size()==1);
    offset=myData1.getDefaultOffset();
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,0)==3.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,1)==6.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,2)==5.0);

    // scalar slice

    region.clear();
    region_element.first=0;
    region_element.second=0;
    region.push_back(region_element);

    DataTagged myData4=makeTagged();
    myData4.getDataAtOffsetRW(myData4.getDefaultOffset())=7.0;

    myData1.setSlice(&myData4, region);

    CPPUNIT_ASSERT(myData1.getTagLookup().size()==0);

    CPPUNIT_ASSERT(myData1.getLength()==3);
    CPPUNIT_ASSERT(myData1.getRank()==1);
    CPPUNIT_ASSERT(myData1.getNoValues()==3);
    CPPUNIT_ASSERT(myData1.getShape().size()==1);

    offset=myData1.getDefaultOffset();
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,0)==7.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,1)==6.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,2)==5.0);

  }

  {

    cout << "\tTest slicing DataTagged with rank 3 default value only." << endl;

    DataTagged::TagListType keys;

    DataTagged::FloatBatchType values;

    DataTypes::ShapeType viewShape;
    viewShape.push_back(3);
    viewShape.push_back(3);
    viewShape.push_back(3);

    DataTypes::RealVectorType viewData1(27);
    for (int i=0;i<viewData1.size();i++) {
      viewData1[i]=i;
    }
    DataTagged myData1(FunctionSpace(),viewShape,viewData1);

    DataTypes::RealVectorType viewData2(27);
    for (int i=0;i<viewData2.size();i++) {
      viewData2[i]=i+27;
    }
    DataTagged myData2(FunctionSpace(),viewShape,viewData2);

    // full slice

    std::pair<int, int> region_element;
    region_element.first=0;
    region_element.second=3;
    DataTypes::RegionType region;
    region.push_back(region_element);
    region.push_back(region_element);
    region.push_back(region_element);

    myData1.setSlice(&myData2, region);


    CPPUNIT_ASSERT(myData1.getTagLookup().size()==0);

    CPPUNIT_ASSERT(myData1.getLength()==27);
    CPPUNIT_ASSERT(myData1.getRank()==3);
    CPPUNIT_ASSERT(myData1.getNoValues()==27);
    CPPUNIT_ASSERT(myData1.getShape().size()==3);

    int offset=myData1.getDefaultOffset();
    CPPUNIT_ASSERT(offset==0);


    // rank 1 slice

    viewShape.clear();
    viewShape.push_back(3);

    DataTypes::RealVectorType viewData3(3);
    for (int i=0;i<viewData3.size();i++) {
      viewData3[i]=i+60;
    }
    DataTagged myData3(FunctionSpace(),viewShape,viewData3);

    region.clear();
    region.push_back(region_element);
    region_element.first=0;
    region_element.second=0;
    region.push_back(region_element);
    region.push_back(region_element);

    myData1.setSlice(&myData3, region);


    CPPUNIT_ASSERT(myData1.getTagLookup().size()==0);

    CPPUNIT_ASSERT(myData1.getLength()==27);
    CPPUNIT_ASSERT(myData1.getRank()==3);
    CPPUNIT_ASSERT(myData1.getNoValues()==27);
    CPPUNIT_ASSERT(myData1.getShape().size()==3);
    offset=myData1.getDefaultOffset();
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,0,0,0)==60.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,1,0,0)==61.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,2,0,0)==62.0);

    // scalar slice

    region.clear();
    region_element.first=0;
    region_element.second=0;
    region.push_back(region_element);
    region.push_back(region_element);
    region.push_back(region_element);

    DataTagged myData4=makeTagged();
    myData4.getDataAtOffsetRW(myData4.getDefaultOffset())=70.0;

    myData1.setSlice(&myData4, region);

    CPPUNIT_ASSERT(myData1.getTagLookup().size()==0);

    CPPUNIT_ASSERT(myData1.getLength()==27);

    CPPUNIT_ASSERT(myData1.getRank()==3);
    CPPUNIT_ASSERT(myData1.getNoValues()==27);
    CPPUNIT_ASSERT(myData1.getShape().size()==3);

    offset=myData1.getDefaultOffset();
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,0,0,0)==70.0);

  }

  {

    cout << "\tTest slicing DataTagged with scalar values and one tag." << endl;

    DataTagged::TagListType keys;
    keys.push_back(1);

    DataTagged::FloatBatchType values;

    DataTypes::ShapeType viewShape;

    // default value for Data1
    DataTypes::RealVectorType viewData1(1*2);
    viewData1[0]=0.0;

    // value for tag "1" for Data1
    viewData1[1]=0.0;

    DataTagged myData1(FunctionSpace(),viewShape,keys,viewData1);

    values.clear();

    // default value for Data2
    DataTypes::RealVectorType viewData3(1*2);
    viewData3[0]=1.0;

    // value for tag "1" for Data2
    viewData3[1]=2.0;

    DataTagged myData2(FunctionSpace(),viewShape,keys,viewData3);

    // full slice

    DataTypes::RegionType region;

    myData1.setSlice(&myData2, region);

    CPPUNIT_ASSERT(myData1.getTagLookup().size()==1);

    CPPUNIT_ASSERT(myData1.getLength()==2);
    CPPUNIT_ASSERT(myData1.getRank()==0);
    CPPUNIT_ASSERT(myData1.getNoValues()==1);
    CPPUNIT_ASSERT(myData1.getShape().size()==0);
    int offset=myData1.getDefaultOffset();
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(myData1.getVectorRO()[offset]==1.0);

    offset=myData1.getOffsetForTag(1);
    CPPUNIT_ASSERT(offset==1);
    CPPUNIT_ASSERT(myData1.getVectorRO()[offset]==2.0);
  }

  {

    cout << "\tTest slicing DataTagged with rank 1 values and one tag." << endl;
    DataTagged::TagListType keys;
    keys.push_back(1);

    DataTagged::FloatBatchType values;

    DataTypes::ShapeType viewShape;
    viewShape.push_back(3);

    int nvals=3;
    // default value for Data1
    DataTypes::RealVectorType viewData1(3*2);
    for (int i=0;i<nvals;i++) {
      viewData1[i]=0.0;
    }

    // value for tag "1" for Data1
    for (int i=0;i<nvals;i++) {
      viewData1[nvals+i]=0.0;
    }

    DataTagged myData1(FunctionSpace(),viewShape,keys,viewData1);
    values.clear();

    // default value for Data2
    DataTypes::RealVectorType viewData3(3*2);
    for (int i=0;i<nvals;i++) {
      viewData3[i]=1.0;
    }

    // value for tag "1" for Data2
    for (int i=0;i<nvals;i++) {
      viewData3[nvals+i]=2.0;
    }

    DataTagged myData2(FunctionSpace(),viewShape,keys,viewData3);

    // full slice

    std::pair<int, int> region_element;
    region_element.first=0;
    region_element.second=3;
    DataTypes::RegionType region;
    region.push_back(region_element);

    myData1.setSlice(&myData2, region);

    CPPUNIT_ASSERT(myData1.getTagLookup().size()==1);

    CPPUNIT_ASSERT(myData1.getLength()==6);
    CPPUNIT_ASSERT(myData1.getRank()==1);
    CPPUNIT_ASSERT(myData1.getNoValues()==3);
    CPPUNIT_ASSERT(myData1.getShape().size()==1);

    int offset=myData1.getDefaultOffset();
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,0)==1.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,1)==1.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,2)==1.0);

    offset=myData1.getOffsetForTag(1);
    CPPUNIT_ASSERT(offset==3);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,0)==2.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,1)==2.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,2)==2.0);

    // rank 1 slice

    viewShape.clear();
    viewShape.push_back(1);

    DataTypes::RealVectorType viewData5(1*2);
    viewData5[0]=3.0;

    values.clear();

    viewData5[1]=4.0;

    DataTagged myData3(FunctionSpace(),viewShape,keys,viewData5);

    region.clear();
    region_element.first=1;
    region_element.second=2;
    region.push_back(region_element);

    myData1.setSlice(&myData3, region);


    CPPUNIT_ASSERT(myData1.getTagLookup().size()==1);

    CPPUNIT_ASSERT(myData1.getLength()==6);
    CPPUNIT_ASSERT(myData1.getRank()==1);
    CPPUNIT_ASSERT(myData1.getNoValues()==3);
    CPPUNIT_ASSERT(myData1.getShape().size()==1);
    offset=myData1.getDefaultOffset();
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,0)==1.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,1)==3.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,2)==1.0);

    offset=myData1.getOffsetForTag(1);
    CPPUNIT_ASSERT(offset==3);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,0)==2.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,1)==4.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,2)==2.0);

    // scalar slice

    viewShape.clear();

    DataTypes::RealVectorType viewData7(1*2);
    viewData7[0]=5.0;

    values.clear();

    viewData7[1]=6.0;

    DataTagged myData4(FunctionSpace(),viewShape,keys,viewData7);

    region.clear();
    region_element.first=0;
    region_element.second=0;
    region.push_back(region_element);

    myData1.setSlice(&myData4, region);

    CPPUNIT_ASSERT(myData1.getRank()==1);
    CPPUNIT_ASSERT(myData1.getNoValues()==3);
    CPPUNIT_ASSERT(myData1.getShape().size()==1);

    offset=myData1.getDefaultOffset();
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,0)==5.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,1)==3.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,2)==1.0);

    offset=myData1.getOffsetForTag(1);
    CPPUNIT_ASSERT(offset==3);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,0)==6.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,1)==4.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,2)==2.0);

  }

  {

    cout << "\tTest slicing DataTagged with rank 3 values and one tag." << endl;
    DataTagged::TagListType keys;
    keys.push_back(1);

    DataTagged::FloatBatchType values;

    DataTypes::ShapeType viewShape;
    viewShape.push_back(3);
    viewShape.push_back(3);
    viewShape.push_back(3);

    int nvals=27;
    // default value for Data1
    DataTypes::RealVectorType viewData1(27*2);
    for (int i=0;i<nvals;i++) {
      viewData1[i]=0.0;
    }

    // value for tag "1" for Data1
    for (int i=0;i<nvals;i++) {
      viewData1[nvals+i]=0.0;
    }

    DataTagged myData1(FunctionSpace(),viewShape,keys,viewData1);

    values.clear();

    // default value for Data2
    DataTypes::RealVectorType viewData3(27*2);
    for (int i=0;i<nvals;i++) {
      viewData3[i]=1.0;
    }

    // value for tag "1" for Data2
    for (int i=0;i<nvals;i++) {
      viewData3[nvals+i]=2.0;
    }

    DataTagged myData2(FunctionSpace(),viewShape,keys,viewData3);

    // full slice

    std::pair<int, int> region_element;
    region_element.first=0;
    region_element.second=3;
    DataTypes::RegionType region;
    region.push_back(region_element);
    region.push_back(region_element);
    region.push_back(region_element);

    myData1.setSlice(&myData2, region);

    CPPUNIT_ASSERT(myData1.getTagLookup().size()==1);

    CPPUNIT_ASSERT(myData1.getLength()==54);

    CPPUNIT_ASSERT(myData1.getRank()==3);
    CPPUNIT_ASSERT(myData1.getNoValues()==27);
    CPPUNIT_ASSERT(myData1.getShape().size()==3);


    int offset=myData1.getDefaultOffset();
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,0,0,0)==1.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,1,1,1)==1.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,2,2,2)==1.0);

    offset=myData1.getOffsetForTag(1);
    CPPUNIT_ASSERT(offset==27);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,0,0,0)==2.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,1,1,1)==2.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,2,2,2)==2.0);

    // rank 1 slice

    viewShape.clear();
    viewShape.push_back(3);
  
    nvals=3;

    DataTypes::RealVectorType viewData5(3*2);
    for (int i=0;i<nvals;i++) {
      viewData5[i]=3.0;
    }

    values.clear();

    for (int i=0;i<nvals;i++) {
      viewData5[nvals+i]=4.0;
    }

    DataTagged myData3(FunctionSpace(),viewShape,keys,viewData5);

    region.clear();
    region.push_back(region_element);
    region_element.first=0;
    region_element.second=0;
    region.push_back(region_element);
    region.push_back(region_element);

    myData1.setSlice(&myData3, region);

    CPPUNIT_ASSERT(myData1.getTagLookup().size()==1);

    CPPUNIT_ASSERT(myData1.getLength()==54);

    CPPUNIT_ASSERT(myData1.getRank()==3);
    CPPUNIT_ASSERT(myData1.getNoValues()==27);
    CPPUNIT_ASSERT(myData1.getShape().size()==3);


    offset=myData1.getDefaultOffset();
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,0,0,0)==3.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,1,0,0)==3.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,2,0,0)==3.0);

    offset=myData1.getOffsetForTag(1);
    CPPUNIT_ASSERT(offset==27);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,0,0,0)==4.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,1,0,0)==4.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,2,0,0)==4.0);

    // scalar slice

    viewShape.clear();

    DataTypes::RealVectorType viewData7(1*2);
    viewData7[0]=5.0;

    values.clear();

    viewData7[1]=6.0;

    DataTagged myData4(FunctionSpace(),viewShape,keys,viewData7);

    region.clear();
    region_element.first=0;
    region_element.second=0;
    region.push_back(region_element);
    region.push_back(region_element);
    region.push_back(region_element);

    myData1.setSlice(&myData4, region);

    CPPUNIT_ASSERT(myData1.getRank()==3);
    CPPUNIT_ASSERT(myData1.getNoValues()==27);
    CPPUNIT_ASSERT(myData1.getShape().size()==3);


    offset=myData1.getDefaultOffset();
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,0,0,0)==5.0);

    offset=myData1.getOffsetForTag(1);
    CPPUNIT_ASSERT(offset==27);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,0,0,0)==6.0);

  }

  {

    cout << "\tTest slicing DataTagged with scalar values and three tags." << endl;
    DataTagged::TagListType keys;
    keys.push_back(1);
    keys.push_back(2);
    keys.push_back(3);

    DataTagged::FloatBatchType values;

    DataTypes::ShapeType viewShape;

    // default value for Data1
    DataTypes::RealVectorType viewData1(1*4);
    viewData1[0]=0.0;

    // value for tag "1" for Data1
    viewData1[1]=0.0;

    // value for tag "2" for Data1
    viewData1[2]=0.0;

    // value for tag "3" for Data1
    viewData1[3]=0.0;

    DataTagged myData1(FunctionSpace(),viewShape,keys,viewData1);

    values.clear();

    // default value for Data2
    DataTypes::RealVectorType viewData3(1*4);
    viewData3[0]=1.0;

    // value for tag "1" for Data2
    viewData3[1]=2.0;

    // value for tag "2" for Data2
    viewData3[2]=3.0;

    // value for tag "3" for Data2
    viewData3[3]=4.0;

    DataTagged myData2(FunctionSpace(),viewShape,keys,viewData3);

    // full slice

    DataTypes::RegionType region;

    myData1.setSlice(&myData2, region);

    CPPUNIT_ASSERT(myData1.getTagLookup().size()==3);

    CPPUNIT_ASSERT(myData1.getLength()==4);

    CPPUNIT_ASSERT(myData1.getRank()==0);
    CPPUNIT_ASSERT(myData1.getNoValues()==1);
    CPPUNIT_ASSERT(myData1.getShape().size()==0);

    int offset=myData1.getDefaultOffset();
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(myData1.getDataAtOffsetRO(offset)==1.0);

    offset=myData1.getOffsetForTag(1);
    CPPUNIT_ASSERT(offset==1);
    CPPUNIT_ASSERT(myData1.getDataAtOffsetRO(offset)==2.0);

    offset=myData1.getOffsetForTag(2);
    CPPUNIT_ASSERT(offset==2);

    CPPUNIT_ASSERT(myData1.getDataAtOffsetRO(offset)==3.0);

    offset=myData1.getOffsetForTag(3);
    CPPUNIT_ASSERT(offset==3);
    CPPUNIT_ASSERT(myData1.getDataAtOffsetRO(offset)==4.0);

  }

  {

    cout << "\tTest slicing DataTagged with rank 1 values and three tags." << endl;
    DataTagged::TagListType keys;
    keys.push_back(1);
    keys.push_back(2);
    keys.push_back(3);

    DataTagged::FloatBatchType values;

    DataTypes::ShapeType viewShape;
    viewShape.push_back(3);

    int nvals=3;

    // default value for Data1
    DataTypes::RealVectorType viewData1(3*4);
    for (int i=0;i<viewData1.size();i++) {
      viewData1[i]=0.0;
    }

    // value for tag "1" for Data1
    for (int i=0;i<nvals;i++) {
      viewData1[nvals+i]=0.0;
    }

    // value for tag "2" for Data1
    for (int i=0;i<nvals;i++) {
      viewData1[2*nvals+i]=0.0;
    }

    // value for tag "3" for Data1
    for (int i=0;i<nvals;i++) {
      viewData1[3*nvals+i]=0.0;
    }

    DataTagged myData1(FunctionSpace(),viewShape,keys,viewData1);

    values.clear();

    nvals=3;

    // default value for Data2
    DataTypes::RealVectorType viewData5(3*4);
    for (int i=0;i<nvals;i++) {
      viewData5[i]=1.0;
    }

    // value for tag "1" for Data2
    for (int i=0;i<nvals;i++) {
      viewData5[nvals+i]=2.0;
    }

    // value for tag "2" for Data2
    for (int i=0;i<nvals;i++) {
      viewData5[2*nvals+i]=3.0;
    }

    // value for tag "3" for Data2
    for (int i=0;i<nvals;i++) {
      viewData5[3*nvals+i]=4.0;
    }
    DataTagged myData2(FunctionSpace(),viewShape,keys,viewData5);

    // full slice

    std::pair<int, int> region_element;
    region_element.first=0;
    region_element.second=3;
    DataTypes::RegionType region;
    region.push_back(region_element);

    myData1.setSlice(&myData2, region);

    CPPUNIT_ASSERT(myData1.getTagLookup().size()==3);

    CPPUNIT_ASSERT(myData1.getLength()==12);

    CPPUNIT_ASSERT(myData1.getRank()==1);
    CPPUNIT_ASSERT(myData1.getNoValues()==3);
    CPPUNIT_ASSERT(myData1.getShape().size()==1);

    int offset=myData1.getDefaultOffset();
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,0)==1.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,1)==1.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,2)==1.0);

    offset=myData1.getOffsetForTag(1);
    CPPUNIT_ASSERT(offset==3);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,0)==2.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,1)==2.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,2)==2.0);

    offset=myData1.getOffsetForTag(2);
    CPPUNIT_ASSERT(offset==6);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,0)==3.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,1)==3.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,2)==3.0);

    offset=myData1.getOffsetForTag(3);
    CPPUNIT_ASSERT(offset==9);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,0)==4.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,1)==4.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,2)==4.0);

    // rank 1 slice

    viewShape.clear();
    viewShape.push_back(1);

    DataTypes::RealVectorType viewData9(1*4);
    viewData9[0]=6.0;

    values.clear();

    viewData9[1]=7.0;
    viewData9[2]=8.0;
    viewData9[3]=9.0;

    DataTagged myData3(FunctionSpace(),viewShape, keys, viewData9);

    region.clear();
    region_element.first=1;
    region_element.second=2;
    region.push_back(region_element);

    myData1.setSlice(&myData3, region);

    CPPUNIT_ASSERT(myData1.getTagLookup().size()==3);

    CPPUNIT_ASSERT(myData1.getLength()==12);

    CPPUNIT_ASSERT(myData1.getRank()==1);
    CPPUNIT_ASSERT(myData1.getNoValues()==3);
    CPPUNIT_ASSERT(myData1.getShape().size()==1);


    offset=myData1.getDefaultOffset();
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,0)==1.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,1)==6.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,2)==1.0);

    offset=myData1.getOffsetForTag(1);
    CPPUNIT_ASSERT(offset==3);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,0)==2.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,1)==7.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,2)==2.0);

    offset=myData1.getOffsetForTag(2);
    CPPUNIT_ASSERT(offset==6);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,0)==3.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,1)==8.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,2)==3.0);

    offset=myData1.getOffsetForTag(3);
    CPPUNIT_ASSERT(offset==9);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,0)==4.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,1)==9.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,2)==4.0);

    // scalar slice

    viewShape.clear();

    DataTypes::RealVectorType viewData13(1*4);
    viewData13[0]=10.0;

    values.clear();

    viewData13[1]=11.0;
    viewData13[2]=12.0;
    viewData13[3]=13.0;

    DataTagged myData4(FunctionSpace(),viewShape,keys,viewData13);

    region.clear();
    region_element.first=0;
    region_element.second=0;
    region.push_back(region_element);

    myData1.setSlice(&myData4, region);

    CPPUNIT_ASSERT(myData1.getRank()==1);
    CPPUNIT_ASSERT(myData1.getNoValues()==3);
    CPPUNIT_ASSERT(myData1.getShape().size()==1);


    offset=myData1.getDefaultOffset();
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,0)==10.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,1)==6.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,2)==1.0);

    offset=myData1.getOffsetForTag(1);
    CPPUNIT_ASSERT(offset==3);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,0)==11.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,1)==7.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,2)==2.0);

    offset=myData1.getOffsetForTag(2);
    CPPUNIT_ASSERT(offset==6);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,0)==12.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,1)==8.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,2)==3.0);

    offset=myData1.getOffsetForTag(3);
    CPPUNIT_ASSERT(offset==9);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,0)==13.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,1)==9.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,2)==4.0);

  }

  {

    cout << "\tTest slicing DataTagged with rank 3 values and three tags." << endl;

    DataTagged::TagListType keys;
    keys.push_back(1);
    keys.push_back(2);
    keys.push_back(3);

    DataTagged::FloatBatchType values;

    DataTypes::ShapeType viewShape;
    viewShape.push_back(3);
    viewShape.push_back(3);
    viewShape.push_back(3);
    int nvals=noValues(viewShape);

    // default value for Data1
    DataTypes::RealVectorType viewData1(27*4);
    for (int i=0;i<nvals;i++) {
      viewData1[i]=0.0;
    }

    // value for tag "1" for Data1
    for (int i=0;i<nvals;i++) {
      viewData1[nvals+i]=0.0;
    }

    // value for tag "2" for Data1
    for (int i=0;i<nvals;i++) {
      viewData1[2*nvals+i]=0.0;
    }

    // value for tag "3" for Data1
    for (int i=0;i<nvals;i++) {
      viewData1[3*nvals+i]=0.0;
    }

    DataTagged myData1(FunctionSpace(),viewShape,keys,viewData1);

    values.clear();

    // default value for Data2
    DataTypes::RealVectorType viewData5(27*4);
    for (int i=0;i<nvals;i++) {
      viewData5[i]=1.0;
    }

    // value for tag "1" for Data2
    for (int i=0;i<nvals;i++) {
      viewData5[nvals+i]=2.0;
    }

    // value for tag "2" for Data2
    for (int i=0;i<nvals;i++) {
      viewData5[2*nvals+i]=3.0;
    }
    // value for tag "3" for Data2
    DataTypes::RealVectorType viewData8(27);
    for (int i=0;i<nvals;i++) {
      viewData5[3*nvals+i]=4.0;
    }
    DataTagged myData2(FunctionSpace(),viewShape,keys,viewData5);

    // full slice

    std::pair<int, int> region_element;
    region_element.first=0;
    region_element.second=3;
    DataTypes::RegionType region;
    region.push_back(region_element);
    region.push_back(region_element);
    region.push_back(region_element);

    myData1.setSlice(&myData2, region);

    CPPUNIT_ASSERT(myData1.getTagLookup().size()==3);

    CPPUNIT_ASSERT(myData1.getLength()==108);

    CPPUNIT_ASSERT(myData1.getRank()==3);
    CPPUNIT_ASSERT(myData1.getNoValues()==27);
    CPPUNIT_ASSERT(myData1.getShape().size()==3);

    int offset=myData1.getDefaultOffset();
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,0,0,0)==1.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,1,0,0)==1.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,2,0,0)==1.0);

    offset=myData1.getOffsetForTag(1);
    CPPUNIT_ASSERT(offset==27);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,0,0,0)==2.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,1,0,0)==2.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,2,0,0)==2.0);

    offset=myData1.getOffsetForTag(2);
    CPPUNIT_ASSERT(offset==54);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,0,0,0)==3.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,1,0,0)==3.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,2,0,0)==3.0);

    offset=myData1.getOffsetForTag(3);
    CPPUNIT_ASSERT(offset==81);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,0,0,0)==4.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,1,0,0)==4.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,2,0,0)==4.0);

    // rank 1 slice

    viewShape.clear();
    viewShape.push_back(3);

    nvals=3;

    DataTypes::RealVectorType viewData9(3*4);
    for (int i=0;i<nvals;i++) {
      viewData9[i]=6.0;
    }
    values.clear();

    for (int i=0;i<nvals;i++) {
      viewData9[nvals+i]=7.0;
    }
    for (int i=0;i<nvals;i++) {
      viewData9[2*nvals+i]=8.0;
    }
    for (int i=0;i<nvals;i++) {
      viewData9[3*nvals+i]=9.0;
    }

    DataTagged myData3(FunctionSpace(),viewShape,keys,viewData9);

    region.clear();
    region_element.first=0;
    region_element.second=3;
    region.push_back(region_element);
    region_element.second=0;
    region.push_back(region_element);
    region.push_back(region_element);

    myData1.setSlice(&myData3, region);

    CPPUNIT_ASSERT(myData1.getTagLookup().size()==3);

    CPPUNIT_ASSERT(myData1.getLength()==108);

    CPPUNIT_ASSERT(myData1.getRank()==3);
    CPPUNIT_ASSERT(myData1.getNoValues()==27);
    CPPUNIT_ASSERT(myData1.getShape().size()==3);

    offset=myData1.getDefaultOffset();
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,0,0,0)==6.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,1,0,0)==6.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,2,0,0)==6.0);

    offset=myData1.getOffsetForTag(1);
    CPPUNIT_ASSERT(offset==27);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,0,0,0)==7.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,1,0,0)==7.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,2,0,0)==7.0);

    offset=myData1.getOffsetForTag(2);
    CPPUNIT_ASSERT(offset==54);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,0,0,0)==8.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,1,0,0)==8.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,2,0,0)==8.0);

    offset=myData1.getOffsetForTag(3);
    CPPUNIT_ASSERT(offset==81);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,0,0,0)==9.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,1,0,0)==9.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,2,0,0)==9.0);

    // scalar slice

    viewShape.clear();

    DataTypes::RealVectorType viewData13(1*4);
    viewData13[0]=10.0;

    values.clear();

    viewData13[1]=11.0;
    viewData13[2]=12.0;
    viewData13[3]=13.0;

    DataTagged myData4(FunctionSpace(),viewShape,keys,viewData13);

    region.clear();
    region_element.first=0;
    region_element.second=0;
    region.push_back(region_element);
    region.push_back(region_element);
    region.push_back(region_element);

    myData1.setSlice(&myData4, region);
    offset=myData1.getDefaultOffset();
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,0,0,0)==10.0);
    offset=myData1.getOffsetForTag(1);
    CPPUNIT_ASSERT(offset==27);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,0,0,0)==11.0);

    offset=myData1.getOffsetForTag(2);
    CPPUNIT_ASSERT(offset==54);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,0,0,0)==12.0);

    offset=myData1.getOffsetForTag(3);
    CPPUNIT_ASSERT(offset==81);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,0,0,0)==13.0);

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

    DataTagged::FloatBatchType values;

    DataTypes::ShapeType viewShape;

    // default value for Data1
    DataTypes::RealVectorType viewData1(1*4);
    viewData1[0]=0.0;
    viewData1[1]=0.0;
    viewData1[2]=0.0;
    viewData1[3]=0.0;

    DataTagged myData1(FunctionSpace(),viewShape,keys1,viewData1);

    values.clear();

    // default value for Data2
    DataTypes::RealVectorType viewData3(1*4);
    viewData3[0]=1.0;
    viewData3[1]=2.0;
    viewData3[2]=3.0;
    viewData3[3]=4.0;

    DataTagged myData2(FunctionSpace(),viewShape,keys2,viewData3);

    // full slice

    DataTypes::RegionType region;

    myData1.setSlice(&myData2, region);

    CPPUNIT_ASSERT(myData1.getTagLookup().size()==5);

    CPPUNIT_ASSERT(myData1.getLength()==6);

    CPPUNIT_ASSERT(myData1.getRank()==0);
    CPPUNIT_ASSERT(myData1.getNoValues()==1);
    CPPUNIT_ASSERT(myData1.getShape().size()==0);
    CPPUNIT_ASSERT(myData1.getDataAtOffsetRO(0)==1.0);


    CPPUNIT_ASSERT(myData1.getDefaultOffset()==0);


    CPPUNIT_ASSERT(myData1.getOffsetForTag(1)==1);
    CPPUNIT_ASSERT(myData1.getDataAtOffsetRO(1)==1.0);

    CPPUNIT_ASSERT(myData1.getOffsetForTag(2)==2);
    CPPUNIT_ASSERT(myData1.getDataAtOffsetRO(2)==1.0);

    CPPUNIT_ASSERT(myData1.getOffsetForTag(3)==3);
    CPPUNIT_ASSERT(myData1.getDataAtOffsetRO(3)==2.0);

    CPPUNIT_ASSERT(myData1.getOffsetForTag(4)==4);
    CPPUNIT_ASSERT(myData1.getDataAtOffsetRO(4)==3.0);

    CPPUNIT_ASSERT(myData1.getOffsetForTag(5)==5);
    CPPUNIT_ASSERT(myData1.getDataAtOffsetRO(5)==4.0);

  }

}


TestSuite* DataTaggedTestCase::suite()
{
  // create the suite of tests to perform.
  TestSuite *testSuite = new TestSuite("DataTaggedTestCase");
  testSuite->addTest(new TestCaller<DataTaggedTestCase>(
              "testAll",&DataTaggedTestCase::testAll));
  testSuite->addTest(new TestCaller<DataTaggedTestCase>(
              "testAddTaggedValues",&DataTaggedTestCase::testAddTaggedValues));
  testSuite->addTest(new TestCaller<DataTaggedTestCase>(
              "testSetTaggedValue",&DataTaggedTestCase::testSetTaggedValue));
  testSuite->addTest(new TestCaller<DataTaggedTestCase>(
              "testCopyConstructors",&DataTaggedTestCase::testCopyConstructors));
  testSuite->addTest(new TestCaller<DataTaggedTestCase>(
              "testGetSlice",&DataTaggedTestCase::testGetSlice));
  testSuite->addTest(new TestCaller<DataTaggedTestCase>(
              "testSetSlice",&DataTaggedTestCase::testSetSlice));
  return testSuite;
}

