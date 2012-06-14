
/*******************************************************
*
* Copyright (c) 2003-2012 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


#include "DataTaggedTestCase.h"

#include "esysUtils/EsysException.h"

#include "escript/BinaryOp.h"
#include "escript/DataConstant.h"
#include "escript/DataFactory.h"
#include "escript/DataTagged.h"
#include "escript/DataTypes.h"
#include "escript/DataVector.h"
#include "escript/FunctionSpace.h"
#include "escript/FunctionSpaceFactory.h"
#include "escript/UnaryOp.h"

#include <cppunit/TestCaller.h>
#include <iostream>
#include <functional>
#include <algorithm>


using namespace CppUnit;
using namespace escript;
using namespace esysUtils;
using namespace std;
using namespace escript::DataTypes;

// namespace {
// std::string constr(FunctionSpace& fs)
// {
//    
//    try
//    {
// 	int t[1];
// 	DataTagged dt(fs,DataTypes::scalarShape,t,DataTypes::ValueType());
// 	
// 	return "DataTagged(const FunctionSpace& what, const DataTypes::ShapeType &shape, const int tags[], const ValueType& data) was supposed to throw.";
//    } catch (DataException d){}
//    try
//    {
// 	DataTagged t(fs,DataTypes::scalarShape,DataTagged::TagListType(),DataTypes::ValueType());
// 	return "DataTagged(const FunctionSpace& what, const DataTypes::ShapeType &shape, const TagListType& tags, const ValueType& data) was supposed to throw.";
//    } catch (DataException d){}
//    try
//    {
// 	DataTagged t(fs,DataTypes::scalarShape,DataTypes::ValueType());
// 	return "  DataTagged(const FunctionSpace& what, const DataTypes::ShapeType& shape, const DataTypes::ValueType& defaultvalue, const DataTagged* tagsource=0) was supposed to throw.";
//    } catch (DataException d){}
//    try
//    {
//     	DataTypes::ValueType viewData1(1);
//     	viewData1[0]=0.0;
// 	DataConstant c(fs,DataTypes::scalarShape, viewData1);
// 	DataTagged t(c);
// 	return "DataTagged(const DataConstant& other) was supposed to throw.";
//    } catch (DataException d){}
// 
// }
// 
// }

namespace {

ValueType::const_reference
getRefRO(DataTagged& data,int offset, int i, int j, int k)
{
   return data.getVectorRO()[offset+getRelIndex(data.getShape(),i,j,k)];
}

ValueType::const_reference
getRefRO(DataTagged& data,int offset, int i, int j, int k, int l)
{
   return data.getVectorRO()[offset+getRelIndex(data.getShape(),i,j,k,l)];
}

ValueType::const_reference
getRefRO(DataTagged& data,int offset, int i, int j)
{
   return data.getVectorRO()[offset+getRelIndex(data.getShape(),i,j)];
}

ValueType::const_reference
getRefRO(const DataTagged& data,int offset, int i)
{
   return data.getVectorRO()[offset+getRelIndex(data.getShape(),i)];
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

//     DataArrayView myDataView = myData.getDataPoint(0,0);
//     CPPUNIT_ASSERT(!myDataView.isEmpty());
    CPPUNIT_ASSERT(myData.getPointOffset(0,0)==0);
    CPPUNIT_ASSERT(myData.getRank()==0);
    CPPUNIT_ASSERT(myData.getNoValues()==1);
    CPPUNIT_ASSERT(myData.getShape().size()==0);
    CPPUNIT_ASSERT(myData.getDataAtOffsetRW(0)==0.0);

    // Test non-existent tag returns the default value.
//     myDataView = myData.getDataPointByTag(1);
//     CPPUNIT_ASSERT(!myDataView.isEmpty());
    CPPUNIT_ASSERT(myData.getOffsetForTag(1)==0);
    CPPUNIT_ASSERT(myData.getRank()==0);
    CPPUNIT_ASSERT(myData.getNoValues()==1);
    CPPUNIT_ASSERT(myData.getShape().size()==0);
    CPPUNIT_ASSERT(myData.getDataAtOffsetRW(0)==0.0);

//     myDataView = myData.getDefaultValue();
//     CPPUNIT_ASSERT(!myDataView.isEmpty());
    CPPUNIT_ASSERT(myData.getDefaultOffset()==0);
//     CPPUNIT_ASSERT(myDataView.getRank()==0);		// there is no point in testing this again
//     CPPUNIT_ASSERT(myDataView.noValues()==1);	// since we are not building DataArrayViews
//     CPPUNIT_ASSERT(myDataView.getShape().size()==0);
//     CPPUNIT_ASSERT(myDataView()==0.0);

    // use a non-existent tag so we get a pointer to
    // the first element of the data array
    double* sampleData=myData.getSampleDataByTag(9);
    for (int i=0; i<myData.getLength(); i++) {
      CPPUNIT_ASSERT(sampleData[i]==i);
    }

  }

  {
    cout << "\tTest binaryOp addition of two DataTagged objects with default values only." << endl;

    DataTypes::ShapeType viewShape;
    viewShape.push_back(3);

//     DataTagged::TagListType keys;
// 
//     DataTagged::ValueListType values;

    DataTypes::ValueType viewData(3);
    for (int i=0;i<viewShape[0];i++) {
      viewData[i]=i;
    }

//     DataTagged myData(keys,values,myView,FunctionSpace());
//     DataTagged right(keys,values,myView,FunctionSpace());
    DataTagged myData(FunctionSpace(),viewShape,viewData);
    DataTagged right(FunctionSpace(),viewShape,viewData);


    binaryOp(myData,right,plus<double>());

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

    CPPUNIT_ASSERT(myData.getLength()==3);

    CPPUNIT_ASSERT(myData.getPointOffset(0,0)==0);

    CPPUNIT_ASSERT(myData.getRank()==1);
    CPPUNIT_ASSERT(myData.getNoValues()==3);
    CPPUNIT_ASSERT(myData.getShape().size()==1);


    int offset=myData.getDefaultOffset();
//     DataArrayView myDataView = myData.getDefaultValue();
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(getRefRO(myData,offset,0)==0);
    CPPUNIT_ASSERT(getRefRO(myData,offset,1)==2);
    CPPUNIT_ASSERT(getRefRO(myData,offset,2)==4);

    // use a non-existent tag so we get a pointer to
    // the first element of the data array
    double* sampleData=myData.getSampleDataByTag(9);
    for (int i=0; i<myData.getLength(); i++) {
      CPPUNIT_ASSERT(sampleData[i]==i*2);
    }

  }

  {
    cout << "\tTest binaryOp addition of two DataTagged objects with one identical tag each." << endl;

    DataTagged myData;
    DataTagged right;

    DataVector vOneData(1, 1.0 ,1);
    // create a view with an empty shape, a scalar.
//     DataArrayView vOneView(vOneData,DataTypes::ShapeType());

    myData.addTaggedValue(1,DataTypes::scalarShape,vOneData);
    right.addTaggedValue(1,DataTypes::scalarShape,vOneData);

    binaryOp(myData,right,plus<double>());

    CPPUNIT_ASSERT(myData.getNumSamples()==1);
    CPPUNIT_ASSERT(myData.getNumDPPSample()==1);

    CPPUNIT_ASSERT(myData.validSamplePointNo(0));
    CPPUNIT_ASSERT(myData.validSampleNo(0));
    CPPUNIT_ASSERT(!myData.validSamplePointNo(1));
    CPPUNIT_ASSERT(!myData.validSampleNo(1));

    // data-point 0 has tag number 1 by default
    CPPUNIT_ASSERT(myData.getTagNumber(0)==1);

    CPPUNIT_ASSERT(myData.isCurrentTag(1));

    CPPUNIT_ASSERT(myData.getTagLookup().size()==1);

    CPPUNIT_ASSERT(myData.getLength()==2);

    CPPUNIT_ASSERT(myData.getPointOffset(0,0)==1);


    CPPUNIT_ASSERT(myData.getRank()==0);
    CPPUNIT_ASSERT(myData.getNoValues()==1);
    CPPUNIT_ASSERT(myData.getShape().size()==0);



    // check result value for tag "1"
//     DataArrayView myDataView = myData.getDataPointByTag(1);
    int offset=myData.getOffsetForTag(1);
    CPPUNIT_ASSERT(offset==1);
    CPPUNIT_ASSERT(myData.getVectorRO()[offset]==2.0);

    // check result for default value
//     myDataView = myData.getDefaultValue();
    offset=myData.getDefaultOffset();
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(myData.getVectorRO()[offset]==0.0);

    // use a non-existent tag so we get a pointer to
    // the first element of the data array
    double* sampleData=myData.getSampleDataByTag(9);
    for (int i=0; i<myData.getLength(); i++) {
      CPPUNIT_ASSERT(sampleData[i]==i*2);
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
//     myData.getDefaultValue()()=1.0;
//     right.getDefaultValue()()=2.0;
    myData.getVectorRW()[myData.getDefaultOffset()]=1.0;
    right.getVectorRW()[right.getDefaultOffset()]=2.0;

    DataVector vOneData(1, 3.0 ,1);
    // create a view with an empty shape, a scalar.
//     DataArrayView vOneView(vOneData,DataTypes::ShapeType());

    DataVector vTwoData(1, 4.0 ,1);
    // create a view with an empty shape, a scalar.
//     DataArrayView vTwoView(vTwoData,DataTypes::ShapeType());

    myData.addTaggedValue(1,DataTypes::scalarShape,vOneData);
    right.addTaggedValue(2,DataTypes::scalarShape,vTwoData);

    //cout << myData.toString() << endl;
    //cout << right.toString() << endl;

    binaryOp(myData,right,plus<double>());

    //cout << myData.toString() << endl;

    CPPUNIT_ASSERT(myData.getNumSamples()==1);
    CPPUNIT_ASSERT(myData.getNumDPPSample()==1);

    CPPUNIT_ASSERT(myData.validSamplePointNo(0));
    CPPUNIT_ASSERT(myData.validSampleNo(0));
    CPPUNIT_ASSERT(!myData.validSamplePointNo(1));
    CPPUNIT_ASSERT(!myData.validSampleNo(1));

    // data-point 0 has tag number 1 by default
    CPPUNIT_ASSERT(myData.getTagNumber(0)==1);

    CPPUNIT_ASSERT(myData.isCurrentTag(1));
    CPPUNIT_ASSERT(myData.isCurrentTag(2));

    CPPUNIT_ASSERT(myData.getTagLookup().size()==2);

    CPPUNIT_ASSERT(myData.getLength()==3);

    CPPUNIT_ASSERT(myData.getPointOffset(0,0)==1);

    CPPUNIT_ASSERT(myData.getRank()==0);
    CPPUNIT_ASSERT(myData.getNoValues()==1);
    CPPUNIT_ASSERT(myData.getShape().size()==0);


    // check result value for tag "1"
//     DataArrayView myDataView = myData.getDataPointByTag(1);
    int offset=myData.getOffsetForTag(1);
    CPPUNIT_ASSERT(offset==1);
    CPPUNIT_ASSERT(myData.getDataAtOffsetRO(offset)==5.0);

    // check result value for tag "2"
//     myDataView = myData.getDataPointByTag(2);
    offset=myData.getOffsetForTag(2);
    CPPUNIT_ASSERT(offset==2);
    CPPUNIT_ASSERT(myData.getDataAtOffsetRO(offset)==5.0);

    // check result for default value
//     myDataView = myData.getDefaultValue();
    offset=myData.getDefaultOffset();
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(myData.getDataAtOffsetRO(offset)==3.0);

    // use a non-existent tag so we get a pointer to
    // the first element of the data array
    double* sampleData=myData.getSampleDataByTag(9);
    CPPUNIT_ASSERT(sampleData[0]==3);
    CPPUNIT_ASSERT(sampleData[1]==5);
    CPPUNIT_ASSERT(sampleData[2]==5);

  }

  {
    cout << "\tTest binaryOp addition of two DataTagged objects with overlapping tag sets." << endl;

    DataTagged myData;
    DataTagged right;

    // it's important that default values are different, as we need to be able to
    // verify that the tag values in each object are being added to the correct
    // default values - since the tag lists don't match, the default values will
    // be used for missing tags in each object
/*    myData.getDefaultValue()()=2.0;
    right.getDefaultValue()()=3.0;*/
    myData.getVectorRW()[myData.getDefaultOffset()]=2.0;
    right.getVectorRW()[right.getDefaultOffset()]=3.0;


    DataVector vOneData(1, 1.0 ,1);
    // create a view with an empty shape, a scalar.
//     DataArrayView vOneView(vOneData,DataTypes::ShapeType());

    myData.addTaggedValue(1,DataTypes::scalarShape,vOneData);
    myData.addTaggedValue(2,DataTypes::scalarShape,vOneData);
    right.addTaggedValue(2,DataTypes::scalarShape,vOneData);
    right.addTaggedValue(3,DataTypes::scalarShape,vOneData);

    //cout << myData.toString() << endl;
    //cout << right.toString() << endl;

    binaryOp(myData,right,plus<double>());

    //cout << myData.toString() << endl;

    CPPUNIT_ASSERT(myData.getNumSamples()==1);
    CPPUNIT_ASSERT(myData.getNumDPPSample()==1);

    CPPUNIT_ASSERT(myData.validSamplePointNo(0));
    CPPUNIT_ASSERT(myData.validSampleNo(0));
    CPPUNIT_ASSERT(!myData.validSamplePointNo(1));
    CPPUNIT_ASSERT(!myData.validSampleNo(1));

    // data-point 0 has tag number 1 by default
    CPPUNIT_ASSERT(myData.getTagNumber(0)==1);

    CPPUNIT_ASSERT(myData.isCurrentTag(1));
    CPPUNIT_ASSERT(myData.isCurrentTag(2));
    CPPUNIT_ASSERT(myData.isCurrentTag(3));

    CPPUNIT_ASSERT(myData.getTagLookup().size()==3);

    CPPUNIT_ASSERT(myData.getLength()==4);

    CPPUNIT_ASSERT(myData.getPointOffset(0,0)==1);

    CPPUNIT_ASSERT(myData.getRank()==0);
    CPPUNIT_ASSERT(myData.getNoValues()==1);
    CPPUNIT_ASSERT(myData.getShape().size()==0);


    // check result value for tag "1"
//     DataArrayView myDataView = myData.getDataPointByTag(1);
    int offset=myData.getOffsetForTag(1);
    CPPUNIT_ASSERT(offset==1);
    CPPUNIT_ASSERT(myData.getDataAtOffsetRO(offset)==4.0);

    // check result value for tag "2"
//     myDataView = myData.getDataPointByTag(2);
    offset=myData.getOffsetForTag(2);
    CPPUNIT_ASSERT(offset==2);
    CPPUNIT_ASSERT(myData.getDataAtOffsetRO(offset)==2.0);

    // check result value for tag "3"
//     myDataView = myData.getDataPointByTag(3);
    offset=myData.getOffsetForTag(3);
    CPPUNIT_ASSERT(offset==3);
    CPPUNIT_ASSERT(myData.getDataAtOffsetRO(offset)==3.0);

    // check result for default value
//     myDataView = myData.getDefaultValue();
    offset=myData.getDefaultOffset();
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(myData.getDataAtOffsetRO(offset)==5.0);

    // use a non-existent tag so we get a pointer to
    // the first element of the data array
    double* sampleData=myData.getSampleDataByTag(9);
    CPPUNIT_ASSERT(sampleData[0]==5);
    CPPUNIT_ASSERT(sampleData[1]==4);
    CPPUNIT_ASSERT(sampleData[2]==2);
    CPPUNIT_ASSERT(sampleData[3]==3);

  }

  {
    cout << "\tTest binaryOp multiplication of two DataTagged objects with default values only." << endl;

    DataTypes::ShapeType viewShape;
    viewShape.push_back(3);

//     DataTagged::TagListType keys;

//     DataTagged::ValueListType values;

    DataTypes::ValueType viewData(3);
    for (int i=0;i<viewShape[0];i++) {
      viewData[i]=i;
    }
//     DataArrayView myView(viewData,viewShape);

    DataTagged myData(FunctionSpace(),viewShape,viewData);
    DataTagged right(FunctionSpace(),viewShape,viewData);

    binaryOp(myData,right,multiplies<double>());

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

    CPPUNIT_ASSERT(myData.getLength()==3);

    CPPUNIT_ASSERT(myData.getPointOffset(0,0)==0);

    CPPUNIT_ASSERT(myData.getRank()==1);
    CPPUNIT_ASSERT(myData.getNoValues()==3);
    CPPUNIT_ASSERT(myData.getShape().size()==1);


//     DataArrayView myDataView = myData.getDefaultValue();
    int offset=myData.getDefaultOffset();
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(getRefRO(myData,offset,0)==0);
    CPPUNIT_ASSERT(getRefRO(myData,offset,1)==1);
    CPPUNIT_ASSERT(getRefRO(myData,offset,2)==4);

    // use a non-existent tag so we get a pointer to
    // the first element of the data array
    double* sampleData=myData.getSampleDataByTag(9);
    for (int i=0; i<myData.getLength(); i++) {
      CPPUNIT_ASSERT(sampleData[i]==i*i);
    }

  }

  {

    cout << "\tTest binaryOp multiplication of DataTagged object with a scalar." << endl;

    DataTagged myData;

    DataVector vOneData(1, 1.0 ,1);
    // create a view with an empty shape, a scalar.
//     DataArrayView vOneView(vOneData,DataTypes::ShapeType());

    DataVector vTwoData(1, 2.0 ,1);
    // create a view with an empty shape, a scalar.
//     DataArrayView vTwoView(vTwoData,DataTypes::ShapeType());

    myData.addTaggedValue(1,DataTypes::scalarShape,vOneData);
    myData.addTaggedValue(2,DataTypes::scalarShape,vTwoData);

    DataVector vThreeData(1, 3.0 ,1);
    // create a view with an empty shape, a scalar.
//     DataArrayView vThreeView(vThreeData,DataTypes::ShapeType());

//     DataArrayView right=vThreeView;

    //cout << myData.toString() << endl;
    //cout << right.toString() << endl;

    binaryOp(myData,vThreeData, DataTypes::scalarShape,multiplies<double>());

    //cout << myData.toString() << endl;

    CPPUNIT_ASSERT(myData.getNumSamples()==1);
    CPPUNIT_ASSERT(myData.getNumDPPSample()==1);

    CPPUNIT_ASSERT(myData.validSamplePointNo(0));
    CPPUNIT_ASSERT(myData.validSampleNo(0));
    CPPUNIT_ASSERT(!myData.validSamplePointNo(1));
    CPPUNIT_ASSERT(!myData.validSampleNo(1));

    // data-point 0 has tag number 1 by default
    CPPUNIT_ASSERT(myData.getTagNumber(0)==1);

    CPPUNIT_ASSERT(myData.isCurrentTag(1));
    CPPUNIT_ASSERT(myData.isCurrentTag(2));

    CPPUNIT_ASSERT(myData.getTagLookup().size()==2);

    CPPUNIT_ASSERT(myData.getLength()==3);

    CPPUNIT_ASSERT(myData.getPointOffset(0,0)==1);

    CPPUNIT_ASSERT(myData.getRank()==0);
    CPPUNIT_ASSERT(myData.getNoValues()==1);
    CPPUNIT_ASSERT(myData.getShape().size()==0);

    // check result value for tag "1"
//     DataArrayView myDataView = myData.getDataPointByTag(1);
    int offset=myData.getOffsetForTag(1);
    CPPUNIT_ASSERT(offset==1);
    CPPUNIT_ASSERT(myData.getDataAtOffsetRO(offset)==3.0);

    // check result value for tag "2"
//     myDataView = myData.getDataPointByTag(2);
    offset=myData.getOffsetForTag(2);
    CPPUNIT_ASSERT(offset==2);
    CPPUNIT_ASSERT(myData.getDataAtOffsetRO(offset)==6.0);

    // check result for default value
//     myDataView = myData.getDefaultValue();
    offset=myData.getDefaultOffset();
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(myData.getDataAtOffsetRO(offset)==0.0);

    // use a non-existent tag so we get a pointer to
    // the first element of the data array
    double* sampleData=myData.getSampleDataByTag(9);
    CPPUNIT_ASSERT(sampleData[0]==0);
    CPPUNIT_ASSERT(sampleData[1]==3);
    CPPUNIT_ASSERT(sampleData[2]==6);

  }

  {
    cout << "\tTest binaryOp multiplication of two DataTagged objects with overlapping tag sets." << endl;

    DataTagged myData;
    DataTagged right;

    // it's important that default values are different, as we need to be able to
    // verify that the tag values in each object are being added to the correct
    // default values - since the tag lists don't match, the default values will
    // be used for missing tags in each object
//     myData.getDefaultValue()()=2.0;
//     right.getDefaultValue()()=3.0;
    myData.getVectorRW()[myData.getDefaultOffset()]=2.0;
    right.getVectorRW()[right.getDefaultOffset()]=3.0;

    DataVector vOneData(1, 1.0 ,1);
    // create a view with an empty shape, a scalar.
//     DataArrayView vOneView(vOneData,DataTypes::ShapeType());

    DataVector vTwoData(1, 2.0 ,1);
    // create a view with an empty shape, a scalar.
//     DataArrayView vTwoView(vTwoData,DataTypes::ShapeType());

    myData.addTaggedValue(1,DataTypes::scalarShape,vOneData);
    myData.addTaggedValue(2,DataTypes::scalarShape,vOneData);
    right.addTaggedValue(2,DataTypes::scalarShape,vTwoData);
    right.addTaggedValue(3,DataTypes::scalarShape,vTwoData);

    //cout << myData.toString() << endl;
    //cout << right.toString() << endl;

    binaryOp(myData,right,multiplies<double>());

    //cout << myData.toString() << endl;

    CPPUNIT_ASSERT(myData.getNumSamples()==1);
    CPPUNIT_ASSERT(myData.getNumDPPSample()==1);

    CPPUNIT_ASSERT(myData.validSamplePointNo(0));
    CPPUNIT_ASSERT(myData.validSampleNo(0));
    CPPUNIT_ASSERT(!myData.validSamplePointNo(1));
    CPPUNIT_ASSERT(!myData.validSampleNo(1));

    // data-point 0 has tag number 1 by default
    CPPUNIT_ASSERT(myData.getTagNumber(0)==1);

    CPPUNIT_ASSERT(myData.isCurrentTag(1));
    CPPUNIT_ASSERT(myData.isCurrentTag(2));
    CPPUNIT_ASSERT(myData.isCurrentTag(3));

    CPPUNIT_ASSERT(myData.getTagLookup().size()==3);

    CPPUNIT_ASSERT(myData.getLength()==4);

    CPPUNIT_ASSERT(myData.getPointOffset(0,0)==1);

    CPPUNIT_ASSERT(myData.getRank()==0);
    CPPUNIT_ASSERT(myData.getNoValues()==1);
    CPPUNIT_ASSERT(myData.getShape().size()==0);


    // check result value for tag "1"
//     DataArrayView myDataView = myData.getDataPointByTag(1);
    int offset=myData.getOffsetForTag(1);
    CPPUNIT_ASSERT(offset==1);
    CPPUNIT_ASSERT(myData.getDataAtOffsetRO(offset)==3.0);

    // check result value for tag "2"
//     myDataView = myData.getDataPointByTag(2);
    offset=myData.getOffsetForTag(2);
    CPPUNIT_ASSERT(offset==2);
    CPPUNIT_ASSERT(myData.getDataAtOffsetRO(offset)==2.0);

    // check result value for tag "3"
//     myDataView = myData.getDataPointByTag(3);
    offset=myData.getOffsetForTag(3);
    CPPUNIT_ASSERT(offset==3);
    CPPUNIT_ASSERT(myData.getDataAtOffsetRO(offset)==4.0);

    // check result for default value
//     myDataView = myData.getDefaultValue();
    offset=myData.getDefaultOffset();
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(myData.getDataAtOffsetRO(offset)==6.0);

    // use a non-existent tag so we get a pointer to
    // the first element of the data array
    double* sampleData=myData.getSampleDataByTag(9);
    CPPUNIT_ASSERT(sampleData[0]==6);
    CPPUNIT_ASSERT(sampleData[1]==3);
    CPPUNIT_ASSERT(sampleData[2]==2);
    CPPUNIT_ASSERT(sampleData[3]==4);

  }

  {
    cout << "\tTest unaryOp negate on default DataTagged object." << endl;

    DataTagged myData;

    unaryOp(myData,negate<double>());

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


//     DataArrayView myDataView = myData.getDataPoint(0,0);
    int offset=myData.getPointOffset(0,0);
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(myData.getDataAtOffsetRO(offset)==0.0);

    // Test non-existent tag returns the default value.
//     myDataView = myData.getDataPointByTag(1);
    offset=myData.getOffsetForTag(1);
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(myData.getDataAtOffsetRO(offset)==0.0);

//     myDataView = myData.getDefaultValue();
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
    cout << "\tTest unaryOp negate on DataTagged object with default value only." << endl;

    DataTypes::ShapeType viewShape;
    viewShape.push_back(3);

//     DataTagged::TagListType keys;

//     DataTagged::ValueListType values;

    DataTypes::ValueType viewData(3);
    for (int i=0;i<viewShape[0];i++) {
      viewData[i]=i;
    }
//     DataArrayView myView(viewData,viewShape);

    DataTagged myData(FunctionSpace(),viewShape,viewData);

    unaryOp(myData,negate<double>());

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

    CPPUNIT_ASSERT(myData.getLength()==3);

    CPPUNIT_ASSERT(myData.getPointOffset(0,0)==0);

    CPPUNIT_ASSERT(myData.getRank()==1);
    CPPUNIT_ASSERT(myData.getNoValues()==3);
    CPPUNIT_ASSERT(myData.getShape().size()==1);


    int offset=myData.getDefaultOffset();
//     DataArrayView myDataView = myData.getDefaultValue();
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(getRefRO(myData,offset,0)==0);
    CPPUNIT_ASSERT(getRefRO(myData,offset,1)==-1);
    CPPUNIT_ASSERT(getRefRO(myData,offset,2)==-2);

    // use a non-existent tag so we get a pointer to
    // the first element of the data array
    double* sampleData=myData.getSampleDataByTag(9);
    for (int i=0; i<myData.getLength(); i++) {
      CPPUNIT_ASSERT(sampleData[i]==0-i);
    }

  }

  {
    cout << "\tTest unnaryOp negate on DataTagged object with two tags." << endl;

    DataTagged myData;

    DataVector vOneData(1, 1.0 ,1);
    // create a view with an empty shape, a scalar.
//     DataArrayView vOneView(vOneData,DataTypes::ShapeType());

    DataVector vTwoData(1, 2.0 ,1);
    // create a view with an empty shape, a scalar.
//     DataArrayView vTwoView(vTwoData,DataTypes::ShapeType());

    myData.addTaggedValue(1,DataTypes::scalarShape,vOneData);
    myData.addTaggedValue(2,DataTypes::scalarShape,vTwoData);

    unaryOp(myData,negate<double>());

    CPPUNIT_ASSERT(myData.getNumSamples()==1);
    CPPUNIT_ASSERT(myData.getNumDPPSample()==1);

    CPPUNIT_ASSERT(myData.validSamplePointNo(0));
    CPPUNIT_ASSERT(myData.validSampleNo(0));
    CPPUNIT_ASSERT(!myData.validSamplePointNo(1));
    CPPUNIT_ASSERT(!myData.validSampleNo(1));

    // data-point 0 has tag number 1 by default
    CPPUNIT_ASSERT(myData.getTagNumber(0)==1);

    CPPUNIT_ASSERT(myData.isCurrentTag(1));
    CPPUNIT_ASSERT(myData.isCurrentTag(2));

    CPPUNIT_ASSERT(myData.getTagLookup().size()==2);

    CPPUNIT_ASSERT(myData.getLength()==3);

    CPPUNIT_ASSERT(myData.getPointOffset(0,0)==1);

    CPPUNIT_ASSERT(myData.getRank()==0);
    CPPUNIT_ASSERT(myData.getNoValues()==1);
    CPPUNIT_ASSERT(myData.getShape().size()==0);


    // check result value for tag "1"
//     DataArrayView myDataView = myData.getDataPointByTag(1);
    int offset=myData.getOffsetForTag(1);
    CPPUNIT_ASSERT(offset==1);
    CPPUNIT_ASSERT(myData.getDataAtOffsetRO(offset)==-1.0);

    // check result value for tag "2"
//     myDataView = myData.getDataPointByTag(2);
    offset=myData.getOffsetForTag(2);
    CPPUNIT_ASSERT(offset==2);
    CPPUNIT_ASSERT(myData.getDataAtOffsetRO(offset)==-2.0);

    // check result for default value
//     myDataView = myData.getDefaultValue();
    offset=myData.getDefaultOffset();
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(myData.getDataAtOffsetRO(offset)==0.0);

    // use a non-existent tag so we get a pointer to
    // the first element of the data array
    double* sampleData=myData.getSampleDataByTag(9);
    for (int i=0; i<myData.getLength(); i++) {
      CPPUNIT_ASSERT(sampleData[i]==0-i);
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

    DataTagged::ValueBatchType values;

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


//     DataArrayView myDataView = myData.getDataPoint(0,0);
    int offset=myData.getPointOffset(0,0);
    CPPUNIT_ASSERT(offset==1);
    CPPUNIT_ASSERT(myData.getDataAtOffsetRO(offset)==0.0);

//     myDataView = myData.getDataPointByTag(1);
    offset=myData.getOffsetForTag(1);
    CPPUNIT_ASSERT(offset==1);
    CPPUNIT_ASSERT(myData.getDataAtOffsetRO(offset)==0.0);

//     myDataView = myData.getDefaultValue();
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
    DataTagged myData;

    DataTagged::TagListType keys;
    keys.push_back(1);

    DataTagged::ValueBatchType values;

    DataTypes::ShapeType viewShape;
/*    DataTypes::ValueType viewData(1);
    viewData[0]=1.0;*/
//     DataArrayView myView(viewData,viewShape);
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


//     DataArrayView myDataView = myData.getDataPoint(0,0);
    int offset=myData.getPointOffset(0,0);
    CPPUNIT_ASSERT(offset==1);
    CPPUNIT_ASSERT(myData.getDataAtOffsetRO(offset)==1.0);

//     myDataView = myData.getDataPointByTag(1);
    offset=myData.getOffsetForTag(1);
    CPPUNIT_ASSERT(offset==1);
    CPPUNIT_ASSERT(myData.getDataAtOffsetRO(offset)==1.0);

//     myDataView = myData.getDefaultValue();
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
    DataTagged myData;

    DataTagged::TagListType keys;
    keys.push_back(1);
    keys.push_back(2);
    keys.push_back(3);

    DataTagged::ValueBatchType values;

    DataTypes::ShapeType viewShape;
/*    DataTypes::ValueType viewData(1);
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


//     DataArrayView myDataView = myData.getDataPoint(0,0);
    int offset=myData.getPointOffset(0,0);
    CPPUNIT_ASSERT(offset==1);
    CPPUNIT_ASSERT(myData.getDataAtOffsetRO(offset)==1.0);

//     myDataView = myData.getDataPointByTag(1);
    offset=myData.getOffsetForTag(1);
    CPPUNIT_ASSERT(offset==1);
    CPPUNIT_ASSERT(myData.getDataAtOffsetRO(offset)==1.0);

//     myDataView = myData.getDataPointByTag(2);
    offset=myData.getOffsetForTag(2);
    CPPUNIT_ASSERT(offset==2);
    CPPUNIT_ASSERT(myData.getDataAtOffsetRO(offset)==1.0);

//     myDataView = myData.getDataPointByTag(3);
    offset=myData.getOffsetForTag(3);
    CPPUNIT_ASSERT(offset==3);
    CPPUNIT_ASSERT(myData.getDataAtOffsetRO(offset)==1.0);

//     myDataView = myData.getDefaultValue();
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
    DataTagged myData;

    DataTagged::TagListType keys;
    keys.push_back(1);
    keys.push_back(2);
    keys.push_back(3);

    DataTagged::ValueBatchType values;

    DataTypes::ShapeType viewShape;
/*    DataTypes::ValueType viewData1(1);
    viewData1[0]=1.0;
    DataTypes::ValueType viewData2(1);
    viewData2[0]=2.0;
    DataTypes::ValueType viewData3(1);
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


//     DataArrayView myDataView = myData.getDataPoint(0,0);
    int offset=myData.getPointOffset(0,0);
    CPPUNIT_ASSERT(offset==1);
    CPPUNIT_ASSERT(myData.getDataAtOffsetRO(offset)==1.0);

//     myDataView = myData.getDataPointByTag(1);
    offset=myData.getOffsetForTag(1);
    CPPUNIT_ASSERT(offset==1);
    CPPUNIT_ASSERT(myData.getDataAtOffsetRO(offset)==1.0);

//     myDataView = myData.getDataPointByTag(2);
    offset=myData.getOffsetForTag(2);
    CPPUNIT_ASSERT(offset==2);
    CPPUNIT_ASSERT(myData.getDataAtOffsetRO(offset)==2.0);

//     myDataView = myData.getDataPointByTag(3);
    offset=myData.getOffsetForTag(3);
    CPPUNIT_ASSERT(offset==3);
    CPPUNIT_ASSERT(myData.getDataAtOffsetRO(offset)==3.0);

//     myDataView = myData.getDefaultValue();
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

    DataTagged::ValueBatchType values;

    DataTypes::ValueType viewData(3);
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


//     DataArrayView myDataView = myData.getDataPoint(0,0);
    int offset=myData.getPointOffset(0,0);
    CPPUNIT_ASSERT(offset==3);
    CPPUNIT_ASSERT(getRefRO(myData,offset,0)==0);
    CPPUNIT_ASSERT(getRefRO(myData,offset,1)==1);
    CPPUNIT_ASSERT(getRefRO(myData,offset,2)==2);

//     myDataView = myData.getDataPointByTag(1);
    offset=myData.getOffsetForTag(1);
    CPPUNIT_ASSERT(offset==3);
    CPPUNIT_ASSERT(getRefRO(myData,offset,0)==0);
    CPPUNIT_ASSERT(getRefRO(myData,offset,1)==1);
    CPPUNIT_ASSERT(getRefRO(myData,offset,2)==2);

//     myDataView = myData.getDefaultValue();
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

    DataTagged::ValueBatchType values;

    DataTypes::ValueType viewData(3);
    for (int i=0;i<viewShape[0];i++) {
      viewData[i]=i;
    }
//     DataArrayView myView(viewData,viewShape);

    DataTagged myData(FunctionSpace(),viewShape,viewData);

    keys.push_back(1);

//     DataTypes::ValueType viewData1(3);
    for (int i=0;i<viewShape[0];i++) {
//       viewData1[i]=i+3;
	values.push_back(i+3);
    }
//     DataArrayView myView1(viewData1,viewShape);
//     values.push_back(myView1);

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
//     DataArrayView myDataView = myData.getDataPoint(0,0);
    CPPUNIT_ASSERT(offset==3);
    CPPUNIT_ASSERT(getRefRO(myData,offset,0)==3);
    CPPUNIT_ASSERT(getRefRO(myData,offset,1)==4);
    CPPUNIT_ASSERT(getRefRO(myData,offset,2)==5);

//     myDataView = myData.getDataPointByTag(1);
    offset=myData.getOffsetForTag(1);
    CPPUNIT_ASSERT(offset==3);
    CPPUNIT_ASSERT(getRefRO(myData,offset,0)==3);
    CPPUNIT_ASSERT(getRefRO(myData,offset,1)==4);
    CPPUNIT_ASSERT(getRefRO(myData,offset,2)==5);

//     myDataView = myData.getDefaultValue();
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

    DataTagged::ValueBatchType values;

    DataTypes::ValueType viewData(3);
    for (int i=0;i<viewShape[0];i++) {
      viewData[i]=i;
    }
//     DataArrayView myView(viewData,viewShape);

    DataTagged myData(FunctionSpace(),viewShape,viewData);

    keys.push_back(1);
    keys.push_back(2);
    keys.push_back(3);

//     DataTypes::ValueType viewData1(3);
    for (int i=0;i<viewShape[0];i++) {
//       viewData1[i]=3;
	values.push_back(3);
    }
//     DataArrayView myView1(viewData1,viewShape);
//     values.push_back(myView1);

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

//     DataArrayView myDataView = myData.getDataPoint(0,0);
    int offset=myData.getPointOffset(0,0);
    CPPUNIT_ASSERT(offset==3);
    CPPUNIT_ASSERT(getRefRO(myData,offset,0)==3);
    CPPUNIT_ASSERT(getRefRO(myData,offset,1)==3);
    CPPUNIT_ASSERT(getRefRO(myData,offset,2)==3);

//     myDataView = myData.getDataPointByTag(1);
    offset=myData.getOffsetForTag(1);
    CPPUNIT_ASSERT(offset==3);
    CPPUNIT_ASSERT(getRefRO(myData,offset,0)==3);
    CPPUNIT_ASSERT(getRefRO(myData,offset,1)==3);
    CPPUNIT_ASSERT(getRefRO(myData,offset,2)==3);

//     myDataView = myData.getDataPointByTag(2);
    offset=myData.getOffsetForTag(2);
    CPPUNIT_ASSERT(offset==6);
    CPPUNIT_ASSERT(getRefRO(myData,offset,0)==3);
    CPPUNIT_ASSERT(getRefRO(myData,offset,1)==3);
    CPPUNIT_ASSERT(getRefRO(myData,offset,2)==3);

//     myDataView = myData.getDataPointByTag(3);
    offset=myData.getOffsetForTag(3);
    CPPUNIT_ASSERT(offset==9);
    CPPUNIT_ASSERT(getRefRO(myData,offset,0)==3);
    CPPUNIT_ASSERT(getRefRO(myData,offset,1)==3);
    CPPUNIT_ASSERT(getRefRO(myData,offset,2)==3);

//     myDataView = myData.getDefaultValue();
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

    DataTagged::ValueBatchType values;

    DataTypes::ValueType viewData(3);
    for (int i=0;i<viewShape[0];i++) {
      viewData[i]=i;
    }
//     DataArrayView myView(viewData,viewShape);

    DataTagged myData(FunctionSpace(),viewShape,viewData);

    keys.push_back(1);
    keys.push_back(2);
    keys.push_back(3);

//     DataTypes::ValueType viewData1(3);
    for (int i=0;i<viewShape[0];i++) {
//       viewData1[i]=i+1;
	values.push_back(i+1);
    }
//     DataArrayView myView1(viewData1,viewShape);
//     values.push_back(myView1);

//     DataTypes::ValueType viewData2(3);
    for (int i=0;i<viewShape[0];i++) {
//       viewData2[i]=i+2;
	values.push_back(i+2);
    }
//     DataArrayView myView2(viewData2,viewShape);
//     values.push_back(myView2);

//     DataTypes::ValueType viewData3(3);
    for (int i=0;i<viewShape[0];i++) {
//       viewData3[i]=i+3;
	values.push_back(i+3);
    }
//     DataArrayView myView3(viewData3,viewShape);
//     values.push_back(myView3);

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


//     DataArrayView myDataView = myData.getDataPoint(0,0);
    int offset=myData.getPointOffset(0,0);
    CPPUNIT_ASSERT(offset==3);
    CPPUNIT_ASSERT(getRefRO(myData,offset,0)==1);
    CPPUNIT_ASSERT(getRefRO(myData,offset,1)==2);
    CPPUNIT_ASSERT(getRefRO(myData,offset,2)==3);

//     myDataView = myData.getDataPointByTag(1);
    offset=myData.getOffsetForTag(1);
    CPPUNIT_ASSERT(offset==3);
    CPPUNIT_ASSERT(getRefRO(myData,offset,0)==1);
    CPPUNIT_ASSERT(getRefRO(myData,offset,1)==2);
    CPPUNIT_ASSERT(getRefRO(myData,offset,2)==3);

//     myDataView = myData.getDataPointByTag(2);
    offset=myData.getOffsetForTag(2);
    CPPUNIT_ASSERT(offset==6);
    CPPUNIT_ASSERT(getRefRO(myData,offset,0)==2);
    CPPUNIT_ASSERT(getRefRO(myData,offset,1)==3);
    CPPUNIT_ASSERT(getRefRO(myData,offset,2)==4);

//     myDataView = myData.getDataPointByTag(3);
    offset=myData.getOffsetForTag(3);
    CPPUNIT_ASSERT(offset==9);
    CPPUNIT_ASSERT(getRefRO(myData,offset,0)==3);
    CPPUNIT_ASSERT(getRefRO(myData,offset,1)==4);
    CPPUNIT_ASSERT(getRefRO(myData,offset,2)==5);

//     myDataView = myData.getDefaultValue();
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

    DataTagged::ValueBatchType values;

    DataTypes::ShapeType viewShape;
    viewShape.push_back(3);

    // default value
    DataTypes::ValueType viewData(3*4);
    for (int i=0;i<viewShape[0];i++) {
      viewData[i]=i;
    }
//     DataArrayView myView(viewData,viewShape);

    // value for tag "1"
//     DataTypes::ValueType eOneData(viewData);
//     DataArrayView eOneView(eOneData, viewShape);
    for (int i=0;i<viewShape[0];i++) {
//       eOneView(i)=i+1.0;
      viewData[viewShape[0]+i]=i+1.0;
    }
//     values.push_back(eOneView);

    // value for tag "2"
//     DataTypes::ValueType eTwoData(viewData);
//     DataArrayView eTwoView(eTwoData, viewShape);
    for (int i=0;i<viewShape[0];i++) {
//       eTwoView(i)=i+2.0;
	viewData[2*viewShape[0]+i]=i+2.0;
    }
//     values.push_back(eTwoView);

    // value for tag "3"
//     DataTypes::ValueType eThreeData(viewData);
//     DataArrayView eThreeView(eThreeData, viewShape);
    for (int i=0;i<viewShape[0];i++) {
//       eThreeView(i)=i+3.0;
	viewData[3*viewShape[0]+i]=i+3.0;
    }
//     values.push_back(eThreeView);

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

//     DataArrayView myDataView = myData.getDataPointByTag(4);
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

    DataTagged::ValueBatchType values;

    DataTypes::ShapeType viewShape;
    viewShape.push_back(3);

    // default value
    DataTypes::ValueType viewData(3*4);
    for (int i=0;i<viewShape[0];i++) {
      viewData[i]=i;
    }
//     DataArrayView myView(viewData,viewShape);

    // value for tag "1"
//     DataTypes::ValueType eOneData(viewData);
//     DataArrayView eOneView(eOneData, viewShape);
    for (int i=0;i<viewShape[0];i++) {
//       eOneView(i)=i+1.0;
	viewData[viewShape[0]+i]=i+1.0;
    }
//     values.push_back(eOneView);

    // value for tag "2"
//     DataTypes::ValueType eTwoData(viewData);
//     DataArrayView eTwoView(eTwoData, viewShape);
    for (int i=0;i<viewShape[0];i++) {
//       eTwoView(i)=i+2.0;
	viewData[2*viewShape[0]+i]=i+2.0;
    }
//     values.push_back(eTwoView);

    // value for tag "3"
//     DataTypes::ValueType eThreeData(viewData);
//     DataArrayView eThreeView(eThreeData, viewShape);
    for (int i=0;i<viewShape[0];i++) {
//       eThreeView(i)=i+3.0;
	viewData[3*viewShape[0]+i]=i+3.0;
    }
//     values.push_back(eThreeView);

    DataTagged myData(FunctionSpace(),viewShape,keys,viewData);

    keys.clear();
    keys.push_back(4);

    values.clear();
    // value for tag "4"
//     DataTypes::ValueType eFourData(viewData);
//     DataArrayView eFourView(eFourData, viewShape);
    for (int i=0;i<viewShape[0];i++) {
      values.push_back(i+4.0);
    }
//     values.push_back(eFourView);

    myData.addTaggedValues(keys,values,viewShape);

    CPPUNIT_ASSERT(myData.isCurrentTag(4));

    CPPUNIT_ASSERT(myData.getTagLookup().size()==4);

    CPPUNIT_ASSERT(myData.getLength()==15);

    CPPUNIT_ASSERT(myData.getRank()==1);
    CPPUNIT_ASSERT(myData.getNoValues()==3);
    CPPUNIT_ASSERT(myData.getShape().size()==1);


//     DataArrayView myDataView = myData.getDataPointByTag(4);
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

    DataTagged::ValueBatchType values;

    DataTypes::ShapeType viewShape;
    viewShape.push_back(3);

    // default value
    DataTypes::ValueType viewData(3*4);
    for (int i=0;i<viewShape[0];i++) {
      viewData[i]=i;
    }
//     DataArrayView myView(viewData,viewShape);

    // value for tag "1"
//     DataTypes::ValueType eOneData(viewData);
//     DataArrayView eOneView(eOneData, viewShape);
    for (int i=0;i<viewShape[0];i++) {
      viewData[viewShape[0]+i]=i+1.0;
    }
//     values.push_back(eOneView);

    // value for tag "2"
//     DataTypes::ValueType eTwoData(viewData);
//     DataArrayView eTwoView(eTwoData, viewShape);
    for (int i=0;i<viewShape[0];i++) {
//       eTwoView(i)=i+2.0;
	viewData[2*viewShape[0]+i]=i+2.0;
    }
//     values.push_back(eTwoView);

    // value for tag "3"
//     DataTypes::ValueType eThreeData(viewData);
//     DataArrayView eThreeView(eThreeData, viewShape);
    for (int i=0;i<viewShape[0];i++) {
//       eThreeView(i)=i+3.0;
	viewData[3*viewShape[0]+i]=i+3.0;
    }
//     values.push_back(eThreeView);

    DataTagged myData(FunctionSpace(),viewShape,keys,viewData);

    keys.clear();
    keys.push_back(4);
    keys.push_back(5);
    keys.push_back(6);

    values.clear();
    // value for tags "4", "5" and "6"
//     DataTypes::ValueType eFourData(viewData);
//     DataArrayView eFourView(eFourData, viewShape);
    for (int i=0;i<viewShape[0];i++) {
//       eFourView(i)=i+4.0;
	values.push_back(i+4.0);
    }
//     values.push_back(eFourView);

    myData.addTaggedValues(keys,values,viewShape);

    CPPUNIT_ASSERT(myData.isCurrentTag(4));
    CPPUNIT_ASSERT(myData.isCurrentTag(5));
    CPPUNIT_ASSERT(myData.isCurrentTag(6));

    CPPUNIT_ASSERT(myData.getTagLookup().size()==6);

    CPPUNIT_ASSERT(myData.getLength()==21);

    CPPUNIT_ASSERT(myData.getRank()==1);
    CPPUNIT_ASSERT(myData.getNoValues()==3);
    CPPUNIT_ASSERT(myData.getShape().size()==1);

//     DataArrayView myDataView = myData.getDataPointByTag(4);
    int offset=myData.getOffsetForTag(4);
    CPPUNIT_ASSERT(offset==12);
    CPPUNIT_ASSERT(getRefRO(myData,offset,0)==4);
    CPPUNIT_ASSERT(getRefRO(myData,offset,1)==5);
    CPPUNIT_ASSERT(getRefRO(myData,offset,2)==6);

//     myDataView = myData.getDataPointByTag(5);
    offset=myData.getOffsetForTag(5);
    CPPUNIT_ASSERT(offset==15);
    CPPUNIT_ASSERT(getRefRO(myData,offset,0)==4);
    CPPUNIT_ASSERT(getRefRO(myData,offset,1)==5);
    CPPUNIT_ASSERT(getRefRO(myData,offset,2)==6);

//     myDataView = myData.getDataPointByTag(6);
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

    DataTagged::ValueBatchType values;

    DataTypes::ShapeType viewShape;
    viewShape.push_back(3);

    // default value
    DataTypes::ValueType viewData(3*4);
    for (int i=0;i<viewShape[0];i++) {
      viewData[i]=i;
    }
//     DataArrayView myView(viewData,viewShape);

    // value for tag "1"
//     DataTypes::ValueType eOneData(viewData);
//     DataArrayView eOneView(eOneData, viewShape);
    for (int i=0;i<viewShape[0];i++) {
      viewData[viewShape[0]+i]=i+1.0;
    }
//     values.push_back(eOneView);

    // value for tag "2"
//     DataTypes::ValueType eTwoData(viewData);
//     DataArrayView eTwoView(eTwoData, viewShape);
    for (int i=0;i<viewShape[0];i++) {
//       eTwoView(i)=i+2.0;
	viewData[2*viewShape[0]+i]=i+2.0;
    }
//     values.push_back(eTwoView);

    // value for tag "3"
//     DataTypes::ValueType eThreeData(viewData);
//     DataArrayView eThreeView(eThreeData, viewShape);
    for (int i=0;i<viewShape[0];i++) {
      viewData[3*viewShape[0]+i]=i+3.0;
    }
//     values.push_back(eThreeView);

    DataTagged myData(FunctionSpace(),viewShape,keys,viewData);

    keys.clear();
    keys.push_back(4);
    keys.push_back(5);
    keys.push_back(6);

    values.clear();

    // value for tag "4"
//     DataTypes::ValueType eFourData(viewData);
//     DataArrayView eFourView(eFourData, viewShape);
    for (int i=0;i<viewShape[0];i++) {
      values.push_back(i+4.0);
    }
//     values.push_back(eFourView);

    // value for tag "5"
//     DataTypes::ValueType eFiveData(viewData);
//     DataArrayView eFiveView(eFiveData, viewShape);
    for (int i=0;i<viewShape[0];i++) {
      values.push_back(i+5.0);
    }
//     values.push_back(eFiveView);

    // value for tag "6"
//     DataTypes::ValueType eSixData(viewData);
//     DataArrayView eSixView(eSixData, viewShape);
    for (int i=0;i<viewShape[0];i++) {
//       eSixView(i)=i+6.0;
	values.push_back(i+6.0);
    }
//     values.push_back(eSixView);

    myData.addTaggedValues(keys,values,viewShape);

    CPPUNIT_ASSERT(myData.isCurrentTag(4));
    CPPUNIT_ASSERT(myData.isCurrentTag(5));
    CPPUNIT_ASSERT(myData.isCurrentTag(6));

    CPPUNIT_ASSERT(myData.getTagLookup().size()==6);

    CPPUNIT_ASSERT(myData.getLength()==21);

    CPPUNIT_ASSERT(myData.getRank()==1);
    CPPUNIT_ASSERT(myData.getNoValues()==3);
    CPPUNIT_ASSERT(myData.getShape().size()==1);

//     DataArrayView myDataView = myData.getDataPointByTag(4);
    int offset=myData.getOffsetForTag(4);
    CPPUNIT_ASSERT(offset==12);
    CPPUNIT_ASSERT(getRefRO(myData,offset,0)==4);
    CPPUNIT_ASSERT(getRefRO(myData,offset,1)==5);
    CPPUNIT_ASSERT(getRefRO(myData,offset,2)==6);

//     myDataView = myData.getDataPointByTag(5);
    offset=myData.getOffsetForTag(5);
    CPPUNIT_ASSERT(offset==15);
    CPPUNIT_ASSERT(getRefRO(myData,offset,0)==5);
    CPPUNIT_ASSERT(getRefRO(myData,offset,1)==6);
    CPPUNIT_ASSERT(getRefRO(myData,offset,2)==7);

//     myDataView = myData.getDataPointByTag(6);
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

    DataTagged::ValueBatchType values;

    DataTypes::ShapeType viewShape;
    viewShape.push_back(3);

    // default value
    DataTypes::ValueType viewData(3*4);
    for (int i=0;i<viewShape[0];i++) {
      viewData[i]=i;
    }
//     DataArrayView myView(viewData,viewShape);

    // value for tag "1"
//     DataTypes::ValueType eOneData(viewData);
//     DataArrayView eOneView(eOneData, viewShape);
    for (int i=0;i<viewShape[0];i++) {
      viewData[viewShape[0]+i]=i+1.0;
    }
//     values.push_back(eOneView);

    // value for tag "2"
//     DataTypes::ValueType eTwoData(viewData);
//     DataArrayView eTwoView(eTwoData, viewShape);
    for (int i=0;i<viewShape[0];i++) {
//       eTwoView(i)=i+2.0;
      viewData[2*viewShape[0]+i]=i+2.0;
    }
//     values.push_back(eTwoView);

    // value for tag "3"
//     DataTypes::ValueType eThreeData(viewData);
//     DataArrayView eThreeView(eThreeData, viewShape);
    for (int i=0;i<viewShape[0];i++) {
//       eThreeView(i)=i+3.0;
	viewData[3*viewShape[0]+i]=i+3.0;
    }
//     values.push_back(eThreeView);

    DataTagged myData(FunctionSpace(),viewShape,keys,viewData);

    // new value for tag "2"
    ValueType tmp(viewShape[0]);
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
//     DataArrayView myDataView = myData.getDataPointByTag(2);
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
    DataTagged myData;

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


//     DataArrayView myDataView = myData.getDataPoint(0,0);
    int offset=myData.getPointOffset(0,0);
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(myData.getDataAtOffsetRO(offset)==0.0);

    // Test non-existent tag returns the default value.
//     myDataView = myData.getDataPointByTag(1);
    offset=myData.getOffsetForTag(1);
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(myData.getDataAtOffsetRO(offset)==0.0);

//     myDataView = myData.getDefaultValue();
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

    DataTagged::ValueBatchType values;

    DataTypes::ValueType viewData(3);
    for (int i=0;i<viewShape[0];i++) {
      viewData[i]=i;
    }
//     DataArrayView myView(viewData,viewShape);
    DataTagged myData(FunctionSpace(),viewShape, viewData);
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

    CPPUNIT_ASSERT(myData.getLength()==3);

    CPPUNIT_ASSERT(myData.getPointOffset(0,0)==0);

    CPPUNIT_ASSERT(myData.getRank()==1);
    CPPUNIT_ASSERT(myData.getNoValues()==3);
    CPPUNIT_ASSERT(myData.getShape().size()==1);

//    DataArrayView myDataView = myData.getDataPoint(0,0);
    int offset=myData.getPointOffset(0,0);
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(getRefRO(myData,offset,0)==0);
    CPPUNIT_ASSERT(getRefRO(myData,offset,1)==1);
    CPPUNIT_ASSERT(getRefRO(myData,offset,2)==2);

    // Test non-existent tag returns the default value.
//     myDataView = myData.getDataPointByTag(1);
    offset=myData.getOffsetForTag(1);
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(getRefRO(myData,offset,0)==0);
    CPPUNIT_ASSERT(getRefRO(myData,offset,1)==1);
    CPPUNIT_ASSERT(getRefRO(myData,offset,2)==2);

//     myDataView = myData.getDefaultValue();
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

    DataTagged::ValueBatchType values;

    DataTypes::ShapeType viewShape;
    viewShape.push_back(3);

    // default value
    DataTypes::ValueType viewData(3*2);
    for (int i=0;i<viewShape[0];i++) {
      viewData[i]=i;
    }
//     DataArrayView myView(viewData,viewShape);

    // value for tag "1"
//     DataTypes::ValueType eOneData(viewData);
//     DataArrayView eOneView(eOneData, viewShape);
    for (int i=0;i<viewShape[0];i++) {
//       eOneView(i)=i+1.0;
	viewData[viewShape[0]+i]=i+1.0;
    }
//     values.push_back(eOneView);
    DataTagged myData(FunctionSpace(),viewShape,keys,viewData);

    //cout << myData.toString() << endl;

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
//     DataArrayView myDataView = myData.getDataPoint(0,0);
    CPPUNIT_ASSERT(offset==3);
    CPPUNIT_ASSERT(getRefRO(myData,offset,0)==1);
    CPPUNIT_ASSERT(getRefRO(myData,offset,1)==2);
    CPPUNIT_ASSERT(getRefRO(myData,offset,2)==3);

//     myDataView = myData.getDataPointByTag(1);
    offset=myData.getOffsetForTag(1);
    CPPUNIT_ASSERT(offset==3);
    CPPUNIT_ASSERT(getRefRO(myData,offset,0)==1);
    CPPUNIT_ASSERT(getRefRO(myData,offset,1)==2);
    CPPUNIT_ASSERT(getRefRO(myData,offset,2)==3);

    // Test non-existent tag returns the default value.
//     myDataView = myData.getDataPointByTag(9);
    offset=myData.getOffsetForTag(9);
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(getRefRO(myData,offset,0)==0);
    CPPUNIT_ASSERT(getRefRO(myData,offset,1)==1);
    CPPUNIT_ASSERT(getRefRO(myData,offset,2)==2);

//     myDataView = myData.getDefaultValue();
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

    DataTagged::ValueBatchType values;

    DataTypes::ShapeType viewShape;
    viewShape.push_back(3);

    // default value
    DataTypes::ValueType viewData(3*4);
    for (int i=0;i<viewShape[0];i++) {
      viewData[i]=i;
    }
//     DataArrayView myView(viewData,viewShape);

    // value for tag "1"
//     DataTypes::ValueType eOneData(viewData);
//     DataArrayView eOneView(eOneData, viewShape);
    for (int i=0;i<viewShape[0];i++) {
/*      eOneView(i)=i+1.0;*/
      viewData[viewShape[0]+i]=i+1.0;
    }
//     values.push_back(eOneView);

    // value for tag "2"
//     DataTypes::ValueType eTwoData(viewData);
//     DataArrayView eTwoView(eTwoData, viewShape);
    for (int i=0;i<viewShape[0];i++) {
      viewData[2*viewShape[0]+i]=i+2.0;
    }
//     values.push_back(eTwoView);

    // value for tag "3"
//     DataTypes::ValueType eThreeData(viewData);
//     DataArrayView eThreeView(eThreeData, viewShape);
    for (int i=0;i<viewShape[0];i++) {
/*      eThreeView(i)=i+3.0;*/
      viewData[3*viewShape[0]+i]=i+3.0;
    }
//     values.push_back(eThreeView);

    DataTagged myData(FunctionSpace(),viewShape,keys,viewData);

    //cout << myData.toString() << endl;

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
//     DataArrayView myDataView = myData.getDataPoint(0,0);
    CPPUNIT_ASSERT(offset==3);
    CPPUNIT_ASSERT(getRefRO(myData,offset,0)==1);
    CPPUNIT_ASSERT(getRefRO(myData,offset,1)==2);
    CPPUNIT_ASSERT(getRefRO(myData,offset,2)==3);

//     myDataView = myData.getDataPointByTag(1);
    offset=myData.getOffsetForTag(1);
    CPPUNIT_ASSERT(offset==3);
    CPPUNIT_ASSERT(getRefRO(myData,offset,0)==1);
    CPPUNIT_ASSERT(getRefRO(myData,offset,1)==2);
    CPPUNIT_ASSERT(getRefRO(myData,offset,2)==3);

    // Test non-existent tag returns the default value.
//     myDataView = myData.getDataPointByTag(0);
    offset=myData.getOffsetForTag(0);
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(getRefRO(myData,offset,0)==0);
    CPPUNIT_ASSERT(getRefRO(myData,offset,1)==1);
    CPPUNIT_ASSERT(getRefRO(myData,offset,2)==2);

//     myDataView = myData.getDefaultValue();
    offset=myData.getDefaultOffset();
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(getRefRO(myData,offset,0)==0);
    CPPUNIT_ASSERT(getRefRO(myData,offset,1)==1);
    CPPUNIT_ASSERT(getRefRO(myData,offset,2)==2);

    // Test data-points held for remaining tags
//     myDataView = myData.getDataPointByTag(2);
    offset=myData.getOffsetForTag(2);
    CPPUNIT_ASSERT(offset==6);
    CPPUNIT_ASSERT(getRefRO(myData,offset,0)==2);
    CPPUNIT_ASSERT(getRefRO(myData,offset,1)==3);
    CPPUNIT_ASSERT(getRefRO(myData,offset,2)==4);

//     myDataView = myData.getDataPointByTag(3);
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

    DataTagged::ValueBatchType values;

    DataTypes::ShapeType viewShape;
    viewShape.push_back(3);

    // default value
    DataTypes::ValueType viewData(3*4);
    for (int i=0;i<viewShape[0];i++) {
      viewData[i]=i;
    }
//     DataArrayView myView(viewData,viewShape);

    // value for tag "1"
//     DataTypes::ValueType eOneData(viewData);
//     DataArrayView eOneView(eOneData, viewShape);
    for (int i=0;i<viewShape[0];i++) {
      viewData[viewShape[0]+i]=i+1.0;
    }
//     values.push_back(eOneView);

    // value for tag "2"
//     DataTypes::ValueType eTwoData(viewData);
//     DataArrayView eTwoView(eTwoData, viewShape);
    for (int i=0;i<viewShape[0];i++) {
//       eTwoView(i)=i+2.0;
	viewData[2*viewShape[0]+i]=i+2.0;
    }
//     values.push_back(eTwoView);

    // value for tag "3"
//     DataTypes::ValueType eThreeData(viewData);
//     DataArrayView eThreeView(eThreeData, viewShape);
    for (int i=0;i<viewShape[0];i++) {
//       eThreeView(i)=i+3.0;
	viewData[3*viewShape[0]+i]=i+3.0;
    }
//     values.push_back(eThreeView);

    DataTagged myData(FunctionSpace(),viewShape,keys,viewData);

    DataTagged myDataCopy(myData);

    //cout << myDataCopy.toString() << endl;

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
//     DataArrayView myDataView = myDataCopy.getDataPoint(0,0);
    CPPUNIT_ASSERT(offset==3);
    CPPUNIT_ASSERT(getRefRO(myDataCopy,offset,0)==1);
    CPPUNIT_ASSERT(getRefRO(myDataCopy,offset,1)==2);
    CPPUNIT_ASSERT(getRefRO(myDataCopy,offset,2)==3);

//     myDataView = myDataCopy.getDataPointByTag(1);
    offset=myDataCopy.getOffsetForTag(1);
    CPPUNIT_ASSERT(offset==3);
    CPPUNIT_ASSERT(getRefRO(myDataCopy,offset,0)==1);
    CPPUNIT_ASSERT(getRefRO(myDataCopy,offset,1)==2);
    CPPUNIT_ASSERT(getRefRO(myDataCopy,offset,2)==3);

    // Test non-existent tag returns the default value.
//     myDataView = myDataCopy.getDataPointByTag(0);
    offset=myDataCopy.getOffsetForTag(0);
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(getRefRO(myDataCopy,offset,0)==0);
    CPPUNIT_ASSERT(getRefRO(myDataCopy,offset,1)==1);
    CPPUNIT_ASSERT(getRefRO(myDataCopy,offset,2)==2);

    //myDataView = myDataCopy.getDefaultValue();
    offset=myDataCopy.getDefaultOffset();
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(getRefRO(myDataCopy,offset,0)==0);
    CPPUNIT_ASSERT(getRefRO(myDataCopy,offset,1)==1);
    CPPUNIT_ASSERT(getRefRO(myDataCopy,offset,2)==2);

    // Test data-points held for remaining tags
//     myDataView = myDataCopy.getDataPointByTag(2);
    offset=myDataCopy.getOffsetForTag(2);
    CPPUNIT_ASSERT(offset==6);
    CPPUNIT_ASSERT(getRefRO(myDataCopy,offset,0)==2);
    CPPUNIT_ASSERT(getRefRO(myDataCopy,offset,1)==3);
    CPPUNIT_ASSERT(getRefRO(myDataCopy,offset,2)==4);

//     myDataView = myDataCopy.getDataPointByTag(3);
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
    DataTypes::ValueType data(DataTypes::noValues(shape),0);
//     DataArrayView pointData(data,shape);
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


//     DataArrayView myDataView = myData.getDataPoint(0,0);
    int offset=myData.getPointOffset(0,0);
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(myData.getDataAtOffsetRO(offset)==1.0);

    // Test non-existent tag returns the default value.
//     myDataView = myData.getDataPointByTag(1);
    offset=myData.getOffsetForTag(1);
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(myData.getDataAtOffsetRO(offset)==1.0);

//     myDataView = myData.getDefaultValue();
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

    DataTagged myData;

    DataTypes::RegionType region;

    DataAbstract* slicedDefault = myData.getSlice(region);

    // cout << slicedDefault->toString() << endl;

    const DataTagged* myDataSliced=dynamic_cast<const DataTagged*>(slicedDefault);

    CPPUNIT_ASSERT(myDataSliced->getTagLookup().size()==0);

    CPPUNIT_ASSERT(myDataSliced->getLength()==1);

//     DataArrayView myDataView = myDataSliced->getDefaultValue();
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

    DataTagged::ValueBatchType values;

    DataTypes::ValueType viewData(3);
    for (int i=0;i<viewShape[0];i++) {
      viewData[i]=i;
    }
//     DataArrayView myView(viewData,viewShape);

    DataTagged myData(FunctionSpace(),viewShape,viewData);

    // full slice

    std::pair<int, int> region_element;
    region_element.first=0;
    region_element.second=3;
    DataTypes::RegionType region;
    region.push_back(region_element);

    DataAbstract* slicedDefault = myData.getSlice(region);

    //cout << slicedDefault->toString() << endl;

    const DataTagged* myDataSliced=dynamic_cast<const DataTagged*>(slicedDefault);

    CPPUNIT_ASSERT(myDataSliced->getTagLookup().size()==0);

    CPPUNIT_ASSERT(myDataSliced->getLength()==3);

    CPPUNIT_ASSERT(myDataSliced->getRank()==1);
    CPPUNIT_ASSERT(myDataSliced->getNoValues()==3);
    CPPUNIT_ASSERT(myDataSliced->getShape().size()==1);

//     DataArrayView myDataView = myDataSliced->getDefaultValue();
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

    //cout << slicedDefault->toString() << endl;

    myDataSliced=dynamic_cast<const DataTagged*>(slicedDefault);

    CPPUNIT_ASSERT(myDataSliced->getTagLookup().size()==0);

    CPPUNIT_ASSERT(myDataSliced->getLength()==1);

//     myDataView = myDataSliced->getDefaultValue();
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

    DataTagged::ValueBatchType values;

    DataTypes::ValueType viewData(27);
    for (int i=0;i<viewData.size();i++) {
      viewData[i]=i;
    }
//     DataArrayView myView(viewData,viewShape);

    DataTagged myData(FunctionSpace(),viewShape,viewData);

    //cout << myData.toString() << endl;

    // full slice

    std::pair<int, int> region_element;
    region_element.first=0;
    region_element.second=3;
    DataTypes::RegionType region;
    region.push_back(region_element);
    region.push_back(region_element);
    region.push_back(region_element);

    DataAbstract* slicedDefault = myData.getSlice(region);

    //cout << slicedDefault->toString() << endl;

    const DataTagged* myDataSliced=dynamic_cast<const DataTagged*>(slicedDefault);

    CPPUNIT_ASSERT(myDataSliced->getTagLookup().size()==0);

    CPPUNIT_ASSERT(myDataSliced->getLength()==27);

//     DataArrayView myDataView = myDataSliced->getDefaultValue();
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

    //cout << slicedDefault->toString() << endl;

    myDataSliced=dynamic_cast<const DataTagged*>(slicedDefault);

    CPPUNIT_ASSERT(myDataSliced->getTagLookup().size()==0);

    CPPUNIT_ASSERT(myDataSliced->getLength()==3);

    CPPUNIT_ASSERT(myDataSliced->getRank()==1);
    CPPUNIT_ASSERT(myDataSliced->getNoValues()==3);
    CPPUNIT_ASSERT(myDataSliced->getShape().size()==1);

//     myDataView = myDataSliced->getDefaultValue();
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

    //cout << slicedDefault->toString() << endl;

    myDataSliced=dynamic_cast<const DataTagged*>(slicedDefault);

    CPPUNIT_ASSERT(myDataSliced->getTagLookup().size()==0);

    CPPUNIT_ASSERT(myDataSliced->getLength()==1);

    CPPUNIT_ASSERT(myDataSliced->getRank()==0);
    CPPUNIT_ASSERT(myDataSliced->getNoValues()==1);
    CPPUNIT_ASSERT(myDataSliced->getShape().size()==0);


//     myDataView = myDataSliced->getDefaultValue();
    offset=myDataSliced->getDefaultOffset();
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(myDataSliced->getVectorRO()[0]==26);
    delete slicedDefault;
  }

  {

    cout << "\tTest slicing DataTagged with scalar values and one tag." << endl;

    DataTagged::TagListType keys;
    keys.push_back(1);

    DataTagged::ValueBatchType values;

    DataTypes::ShapeType viewShape;

    // default value
    DataTypes::ValueType viewData(1*2);
    viewData[0]=0.0;
//     DataArrayView myView(viewData,viewShape);

    // value for tag "1"
//     DataTypes::ValueType eOneData(viewData);
//     DataArrayView eOneView(eOneData, viewShape);
    viewData[1]=1.0;
//     values.push_back(eOneView);

    DataTagged myData(FunctionSpace(),viewShape,keys, viewData);

    //cout << myData.toString() << endl;

    // full slice

    DataTypes::RegionType region;

    DataAbstract* slicedDefault = myData.getSlice(region);

    //cout << slicedDefault->toString() << endl;

    const DataTagged* myDataSliced=dynamic_cast<const DataTagged*>(slicedDefault);

    CPPUNIT_ASSERT(myDataSliced->getTagLookup().size()==1);

    CPPUNIT_ASSERT(myDataSliced->getLength()==2);

//     DataArrayView myDataView = myDataSliced->getDefaultValue();
    int offset=myDataSliced->getDefaultOffset();
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(myDataSliced->getRank()==0);
    CPPUNIT_ASSERT(myDataSliced->getNoValues()==1);
    CPPUNIT_ASSERT(myDataSliced->getShape().size()==0);
    CPPUNIT_ASSERT(myDataSliced->getVectorRO()[offset]==0);

//     myDataView = myDataSliced->getDataPointByTag(1);
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

    DataTagged::ValueBatchType values;

    // default value
    DataTypes::ValueType viewData(3*2);
    for (int i=0;i<viewShape[0];i++) {
      viewData[i]=i;
    }
//     DataArrayView myView(viewData,viewShape);

    // value for tag "1"
//     DataTypes::ValueType eOneData(viewData);
//     DataArrayView eOneView(eOneData, viewShape);
    for (int i=0;i<viewShape[0];i++) {
//       eOneView(i)=i+3.0;
       viewData[viewShape[0]+i]=i+3.0;
    }
//     values.push_back(eOneView);

    DataTagged myData(FunctionSpace(),viewShape,keys,viewData);

    //cout << myData.toString() << endl;

    // full slice

    std::pair<int, int> region_element;
    region_element.first=0;
    region_element.second=3;
    DataTypes::RegionType region;
    region.push_back(region_element);

    DataAbstract* slicedDefault = myData.getSlice(region);

    //cout << slicedDefault->toString() << endl;

    const DataTagged* myDataSliced=dynamic_cast<const DataTagged*>(slicedDefault);

    CPPUNIT_ASSERT(myDataSliced->getTagLookup().size()==1);

    CPPUNIT_ASSERT(myDataSliced->getLength()==6);
    CPPUNIT_ASSERT(myDataSliced->getRank()==1);
    CPPUNIT_ASSERT(myDataSliced->getNoValues()==3);
    CPPUNIT_ASSERT(myDataSliced->getShape().size()==1);
//     DataArrayView myDataView = myDataSliced->getDefaultValue();
    int offset=myDataSliced->getDefaultOffset();
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(getRefRO(*myDataSliced,offset,0)==0);
    CPPUNIT_ASSERT(getRefRO(*myDataSliced,offset,1)==1);
    CPPUNIT_ASSERT(getRefRO(*myDataSliced,offset,2)==2);

//     myDataView = myDataSliced->getDataPointByTag(1);
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

    //cout << slicedDefault->toString() << endl;

    myDataSliced=dynamic_cast<const DataTagged*>(slicedDefault);

    CPPUNIT_ASSERT(myDataSliced->getTagLookup().size()==1);

    CPPUNIT_ASSERT(myDataSliced->getLength()==2);

    CPPUNIT_ASSERT(myDataSliced->getRank()==0);
    CPPUNIT_ASSERT(myDataSliced->getNoValues()==1);
    CPPUNIT_ASSERT(myDataSliced->getShape().size()==0);


//     myDataView = myDataSliced->getDefaultValue();
    offset=myDataSliced->getDefaultOffset();
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(myDataSliced->getVectorRO()[offset]==1);

//     myDataView = myDataSliced->getDataPointByTag(1);
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

    DataTagged::ValueBatchType values;

    // default value
    DataTypes::ValueType viewData(27*2);
    for (int i=0;i<noValues(viewShape);i++) {
      viewData[i]=i;
    }
//     DataArrayView myView(viewData,viewShape);

    // value for tag "1"
//     DataTypes::ValueType viewData1(27);
    for (int i=0;i<noValues(viewShape);i++) {
      viewData[noValues(viewShape)+i]=i+27.0;
    }
//     DataArrayView myView1(viewData1,viewShape);
//     values.push_back(myView1);

    DataTagged myData(FunctionSpace(),viewShape,keys,viewData);

    //cout << myData.toString() << endl;

    // full slice

    std::pair<int, int> region_element;
    region_element.first=0;
    region_element.second=3;
    DataTypes::RegionType region;
    region.push_back(region_element);
    region.push_back(region_element);
    region.push_back(region_element);

    DataAbstract* slicedDefault = myData.getSlice(region);

    //cout << slicedDefault->toString() << endl;

    const DataTagged* myDataSliced=dynamic_cast<const DataTagged*>(slicedDefault);

    CPPUNIT_ASSERT(myDataSliced->getTagLookup().size()==1);

    CPPUNIT_ASSERT(myDataSliced->getLength()==54);

    CPPUNIT_ASSERT(myDataSliced->getRank()==3);
    CPPUNIT_ASSERT(myDataSliced->getNoValues()==27);
    CPPUNIT_ASSERT(myDataSliced->getShape().size()==3);

//     DataArrayView myDataView = myDataSliced->getDefaultValue();
    int offset=myDataSliced->getDefaultOffset();
    CPPUNIT_ASSERT(offset==0);

//     myDataView = myDataSliced->getDataPointByTag(1);
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
    //cout << slicedDefault->toString() << endl;

    myDataSliced=dynamic_cast<const DataTagged*>(slicedDefault);

    CPPUNIT_ASSERT(myDataSliced->getTagLookup().size()==1);

    CPPUNIT_ASSERT(myDataSliced->getLength()==6);

    CPPUNIT_ASSERT(myDataSliced->getRank()==1);
    CPPUNIT_ASSERT(myDataSliced->getNoValues()==3);
    CPPUNIT_ASSERT(myDataSliced->getShape().size()==1);

//     myDataView = myDataSliced->getDefaultValue();
    offset=myDataSliced->getDefaultOffset();
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(getRefRO(*myDataSliced,offset,0)==0);
    CPPUNIT_ASSERT(getRefRO(*myDataSliced,offset,1)==1);
    CPPUNIT_ASSERT(getRefRO(*myDataSliced,offset,2)==2);

//     myDataView = myDataSliced->getDataPointByTag(1);
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

    //cout << slicedDefault->toString() << endl;

    myDataSliced=dynamic_cast<const DataTagged*>(slicedDefault);

    CPPUNIT_ASSERT(myDataSliced->getTagLookup().size()==1);

    CPPUNIT_ASSERT(myDataSliced->getLength()==2);

    CPPUNIT_ASSERT(myDataSliced->getRank()==0);
    CPPUNIT_ASSERT(myDataSliced->getNoValues()==1);
    CPPUNIT_ASSERT(myDataSliced->getShape().size()==0);

//     myDataView = myDataSliced->getDefaultValue();
    offset=myDataSliced->getDefaultOffset();
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(myDataSliced->getVectorRO()[offset]==13);

//     myDataView = myDataSliced->getDataPointByTag(1);
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

    DataTagged::ValueBatchType values;

    DataTypes::ShapeType viewShape;

    // default value
    DataTypes::ValueType viewData(1*4);
    viewData[0]=0.0;
//     DataArrayView myView(viewData,viewShape);

    // value for tag "1"
//     DataTypes::ValueType eOneData(viewData);
//     DataArrayView eOneView(eOneData, viewShape);
//     eOneView()=1.0;
    viewData[1]=1.0;
//     values.push_back(eOneView);

    // value for tag "2"
//     DataTypes::ValueType eTwoData(viewData);
//     DataArrayView eTwoView(eTwoData, viewShape);
//     eTwoView()=2.0;
    viewData[2]=2.0;
//     values.push_back(eTwoView);

    // value for tag "3"
//     DataTypes::ValueType eThreeData(viewData);
//     DataArrayView eThreeView(eThreeData, viewShape);
//     eThreeView()=3.0;
    viewData[3]=3.0;
//     values.push_back(eThreeView);

    DataTagged myData(FunctionSpace(),viewShape,keys,viewData);

    // cout << myData.toString() << endl;

    // full slice

    DataTypes::RegionType region;

    DataAbstract* slicedDefault = myData.getSlice(region);

    //cout << slicedDefault->toString() << endl;

    const DataTagged* myDataSliced=dynamic_cast<const DataTagged*>(slicedDefault);

    CPPUNIT_ASSERT(myDataSliced->getTagLookup().size()==3);

    CPPUNIT_ASSERT(myDataSliced->getLength()==4);

    CPPUNIT_ASSERT(myDataSliced->getRank()==0);
    CPPUNIT_ASSERT(myDataSliced->getNoValues()==1);
    CPPUNIT_ASSERT(myDataSliced->getShape().size()==0);

//     DataArrayView myDataView = myDataSliced->getDefaultValue();
    int offset=myDataSliced->getDefaultOffset();
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(myDataSliced->getVectorRO()[offset]==0);

//     myDataView = myDataSliced->getDataPointByTag(1);
    offset=myDataSliced->getOffsetForTag(1);
    CPPUNIT_ASSERT(offset==1);
    CPPUNIT_ASSERT(myDataSliced->getVectorRO()[offset]==1);

//     myDataView = myDataSliced->getDataPointByTag(2);
    offset=myDataSliced->getOffsetForTag(2);
    CPPUNIT_ASSERT(offset==2);
    CPPUNIT_ASSERT(myDataSliced->getVectorRO()[offset]==2);

//     myDataView = myDataSliced->getDataPointByTag(3);
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

    DataTagged::ValueBatchType values;

    // default value
    DataTypes::ValueType viewData(3*4);
    for (int i=0;i<viewShape[0];i++) {
      viewData[i]=i;
    }
//     DataArrayView myView(viewData,viewShape);

    // value for tag "1"
//     DataTypes::ValueType eOneData(viewData);
//     DataArrayView eOneView(eOneData, viewShape);
    for (int i=0;i<viewShape[0];i++) {
//       eOneView(i)=i+3.0;
	viewData[viewShape[0]+i]=i+3.0;
    }
//     values.push_back(eOneView);

    // value for tag "2"
//     DataTypes::ValueType eTwoData(viewData);
//     DataArrayView eTwoView(eTwoData, viewShape);
    for (int i=0;i<viewShape[0];i++) {
//       eTwoView(i)=i+6.0;
	viewData[2*viewShape[0]+i]=i+6.0;
    }
//     values.push_back(eTwoView);

    // value for tag "3"
//     DataTypes::ValueType eThreeData(viewData);
//     DataArrayView eThreeView(eThreeData, viewShape);
    for (int i=0;i<viewShape[0];i++) {
//       eThreeView(i)=i+9.0;
	viewData[3*viewShape[0]+i]=i+9.0;
    }
//     values.push_back(eThreeView);

    DataTagged myData(FunctionSpace(),viewShape, keys, viewData);

    //cout << myData.toString() << endl;

    // full slice

    std::pair<int, int> region_element;
    region_element.first=0;
    region_element.second=3;
    DataTypes::RegionType region;
    region.push_back(region_element);

    DataAbstract* slicedDefault = myData.getSlice(region);

    //cout << slicedDefault->toString() << endl;

    const DataTagged* myDataSliced=dynamic_cast<const DataTagged*>(slicedDefault);

    CPPUNIT_ASSERT(myDataSliced->getTagLookup().size()==3);

    CPPUNIT_ASSERT(myDataSliced->getLength()==12);

    CPPUNIT_ASSERT(myDataSliced->getRank()==1);
    CPPUNIT_ASSERT(myDataSliced->getNoValues()==3);
    CPPUNIT_ASSERT(myDataSliced->getShape().size()==1);


//     DataArrayView myDataView = myDataSliced->getDefaultValue();
    int offset=myDataSliced->getDefaultOffset();
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(getRefRO(*myDataSliced,offset,0)==0);
    CPPUNIT_ASSERT(getRefRO(*myDataSliced,offset,1)==1);
    CPPUNIT_ASSERT(getRefRO(*myDataSliced,offset,2)==2);

//     myDataView = myDataSliced->getDataPointByTag(1);
    offset=myDataSliced->getOffsetForTag(1);
    CPPUNIT_ASSERT(offset==3);
    CPPUNIT_ASSERT(getRefRO(*myDataSliced,offset,0)==3);
    CPPUNIT_ASSERT(getRefRO(*myDataSliced,offset,1)==4);
    CPPUNIT_ASSERT(getRefRO(*myDataSliced,offset,2)==5);

//     myDataView = myDataSliced->getDataPointByTag(2);
    offset=myDataSliced->getOffsetForTag(2);
    CPPUNIT_ASSERT(offset==6);
    CPPUNIT_ASSERT(getRefRO(*myDataSliced,offset,0)==6);
    CPPUNIT_ASSERT(getRefRO(*myDataSliced,offset,1)==7);
    CPPUNIT_ASSERT(getRefRO(*myDataSliced,offset,2)==8);

//     myDataView = myDataSliced->getDataPointByTag(3);
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

    //cout << slicedDefault->toString() << endl;

    myDataSliced=dynamic_cast<const DataTagged*>(slicedDefault);

    CPPUNIT_ASSERT(myDataSliced->getTagLookup().size()==3);

    CPPUNIT_ASSERT(myDataSliced->getLength()==4);

    CPPUNIT_ASSERT(myDataSliced->getRank()==0);
    CPPUNIT_ASSERT(myDataSliced->getNoValues()==1);
    CPPUNIT_ASSERT(myDataSliced->getShape().size()==0);

//     myDataView = myDataSliced->getDefaultValue();
    offset=myDataSliced->getDefaultOffset();
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(myDataSliced->getVectorRO()[offset]==1);

//     myDataView = myDataSliced->getDataPointByTag(1);
    offset=myDataSliced->getOffsetForTag(1);
    CPPUNIT_ASSERT(offset==1);
    CPPUNIT_ASSERT(myDataSliced->getVectorRO()[offset]==4);

//     myDataView = myDataSliced->getDataPointByTag(2);
    offset=myDataSliced->getOffsetForTag(2);
    CPPUNIT_ASSERT(offset==2);
    CPPUNIT_ASSERT(myDataSliced->getVectorRO()[offset]==7);

//     myDataView = myDataSliced->getDataPointByTag(3);
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

    DataTagged::ValueBatchType values;

    int nvals=27;
    // default value
    DataTypes::ValueType viewData(27*4);
    for (int i=0;i<nvals;i++) {
      viewData[i]=i;
    }
//     DataArrayView myView(viewData,viewShape);

    // value for tag "1"
//     DataTypes::ValueType viewData1(27);
    for (int i=0;i<nvals;i++) {
      viewData[nvals+i]=i+27.0;
    }
//     DataArrayView myView1(viewData1,viewShape);
//     values.push_back(myView1);

    // value for tag "2"
//     DataTypes::ValueType viewData2(27);
    for (int i=0;i<nvals;i++) {
      viewData[2*nvals+i]=i+54.0;
    }
//     DataArrayView myView2(viewData2,viewShape);
//     values.push_back(myView2);

    // value for tag "3"
//     DataTypes::ValueType viewData3(27);
    for (int i=0;i<nvals;i++) {
      viewData[3*nvals+i]=i+81.0;
    }
//     DataArrayView myView3(viewData3,viewShape);
//     values.push_back(myView3);

    DataTagged myData(FunctionSpace(),viewShape,keys,viewData);

    //cout << myData.toString() << endl;

    // full slice

    std::pair<int, int> region_element;
    region_element.first=0;
    region_element.second=3;
    DataTypes::RegionType region;
    region.push_back(region_element);
    region.push_back(region_element);
    region.push_back(region_element);

    DataAbstract* slicedDefault = myData.getSlice(region);

    //cout << slicedDefault->toString() << endl;

    const DataTagged* myDataSliced=dynamic_cast<const DataTagged*>(slicedDefault);

    CPPUNIT_ASSERT(myDataSliced->getTagLookup().size()==3);

    CPPUNIT_ASSERT(myDataSliced->getLength()==108);

    CPPUNIT_ASSERT(myDataSliced->getRank()==3);
    CPPUNIT_ASSERT(myDataSliced->getNoValues()==27);
    CPPUNIT_ASSERT(myDataSliced->getShape().size()==3);

//     DataArrayView myDataView = myDataSliced->getDefaultValue();
    int offset=myDataSliced->getDefaultOffset();
    CPPUNIT_ASSERT(offset==0);

//     myDataView = myDataSliced->getDataPointByTag(1);
    offset=myDataSliced->getOffsetForTag(1);
    CPPUNIT_ASSERT(offset==27);

//     myDataView = myDataSliced->getDataPointByTag(2);
    offset=myDataSliced->getOffsetForTag(2);
    CPPUNIT_ASSERT(offset==54);

//     myDataView = myDataSliced->getDataPointByTag(3);
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

    // cout << slicedDefault->toString() << endl;

    myDataSliced=dynamic_cast<const DataTagged*>(slicedDefault);

    CPPUNIT_ASSERT(myDataSliced->getTagLookup().size()==3);

    CPPUNIT_ASSERT(myDataSliced->getLength()==12);

    CPPUNIT_ASSERT(myDataSliced->getRank()==1);
    CPPUNIT_ASSERT(myDataSliced->getNoValues()==3);
    CPPUNIT_ASSERT(myDataSliced->getShape().size()==1);

//     myDataView = myDataSliced->getDefaultValue();
    offset=myDataSliced->getDefaultOffset();
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(getRefRO(*myDataSliced,offset,0)==0);
    CPPUNIT_ASSERT(getRefRO(*myDataSliced,offset,1)==1);
    CPPUNIT_ASSERT(getRefRO(*myDataSliced,offset,2)==2);

//     myDataView = myDataSliced->getDataPointByTag(1);
    offset=myDataSliced->getOffsetForTag(1);
    CPPUNIT_ASSERT(offset==3);
    CPPUNIT_ASSERT(getRefRO(*myDataSliced,offset,0)==27);
    CPPUNIT_ASSERT(getRefRO(*myDataSliced,offset,1)==28);
    CPPUNIT_ASSERT(getRefRO(*myDataSliced,offset,2)==29);

//     myDataView = myDataSliced->getDataPointByTag(2);
    offset=myDataSliced->getOffsetForTag(2);
    CPPUNIT_ASSERT(offset==6);
    CPPUNIT_ASSERT(getRefRO(*myDataSliced,offset,0)==54);
    CPPUNIT_ASSERT(getRefRO(*myDataSliced,offset,1)==55);
    CPPUNIT_ASSERT(getRefRO(*myDataSliced,offset,2)==56);

//     myDataView = myDataSliced->getDataPointByTag(3);
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

    //cout << slicedDefault->toString() << endl;

    myDataSliced=dynamic_cast<const DataTagged*>(slicedDefault);

    CPPUNIT_ASSERT(myDataSliced->getTagLookup().size()==3);

    CPPUNIT_ASSERT(myDataSliced->getLength()==4);
    CPPUNIT_ASSERT(myDataSliced->getRank()==0);
    CPPUNIT_ASSERT(myDataSliced->getNoValues()==1);
    CPPUNIT_ASSERT(myDataSliced->getShape().size()==0);

//     myDataView = myDataSliced->getDefaultValue();
    offset=myDataSliced->getDefaultOffset();
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(myDataSliced->getVectorRO()[offset]==13);

//     myDataView = myDataSliced->getDataPointByTag(1);
    offset=myDataSliced->getOffsetForTag(1);
    CPPUNIT_ASSERT(offset==1);
    CPPUNIT_ASSERT(myDataSliced->getVectorRO()[offset]==40);

//     myDataView = myDataSliced->getDataPointByTag(2);
    offset=myDataSliced->getOffsetForTag(2);
    CPPUNIT_ASSERT(offset==2);
    CPPUNIT_ASSERT(myDataSliced->getVectorRO()[offset]==67);

//     myDataView = myDataSliced->getDataPointByTag(3);
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

    DataTagged myData1;
    DataTagged myData2;

    DataTypes::RegionType region;

    myData2.getDataAtOffsetRW(myData2.getDefaultOffset())=1.0;
    myData1.setSlice(&myData2, region);
    //cout << myData1.toString() << endl;

    CPPUNIT_ASSERT(myData1.getTagLookup().size()==0);

    CPPUNIT_ASSERT(myData1.getLength()==1);
    CPPUNIT_ASSERT(myData1.getRank()==0);
    CPPUNIT_ASSERT(myData1.getNoValues()==1);
    CPPUNIT_ASSERT(myData1.getShape().size()==0);

//     DataArrayView myDataView = myData1.getDefaultValue();
    int offset=myData1.getDefaultOffset();
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(myData1.getDataAtOffsetRW(offset)==1.0);

  }

  {

    cout << "\tTest slicing DataTagged with rank 1 default value only." << endl;

    DataTagged::TagListType keys;

    DataTagged::ValueBatchType values;

    DataTypes::ShapeType viewShape;
    viewShape.push_back(3);

    DataTypes::ValueType viewData1(3);
    for (int i=0;i<viewShape[0];i++) {
      viewData1[i]=i;
    }
//     DataArrayView myView1(viewData1,viewShape);
    DataTagged myData1(FunctionSpace(),viewShape,viewData1);

    DataTypes::ValueType viewData2(3);
    for (int i=0;i<viewShape[0];i++) {
      viewData2[i]=i+3;
    }
//     DataArrayView myView2(viewData2,viewShape);
    DataTagged myData2(FunctionSpace(),viewShape,viewData2);

    // full slice

    std::pair<int, int> region_element;
    region_element.first=0;
    region_element.second=3;
    DataTypes::RegionType region;
    region.push_back(region_element);

    myData1.setSlice(&myData2, region);
    //cout << myData1.toString() << endl;

    CPPUNIT_ASSERT(myData1.getTagLookup().size()==0);

    CPPUNIT_ASSERT(myData1.getLength()==3);

    CPPUNIT_ASSERT(myData1.getRank()==1);
    CPPUNIT_ASSERT(myData1.getNoValues()==3);
    CPPUNIT_ASSERT(myData1.getShape().size()==1);

//     DataArrayView myDataView = myData1.getDefaultValue();
    int offset=myData1.getDefaultOffset();
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,0)==3.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,1)==4.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,2)==5.0);

    // rank 1 slice

    viewShape.clear();
    viewShape.push_back(1);

    DataTypes::ValueType viewData3(1);
    viewData3[0]=6.0;
//     DataArrayView myView3(viewData3,viewShape);
    DataTagged myData3(FunctionSpace(),viewShape,viewData3);

    region.clear();
    region_element.first=1;
    region_element.second=2;
    region.push_back(region_element);

    myData1.setSlice(&myData3, region);
    //cout << myData1.toString() << endl;

    CPPUNIT_ASSERT(myData1.getTagLookup().size()==0);

    CPPUNIT_ASSERT(myData1.getLength()==3);
    CPPUNIT_ASSERT(myData1.getRank()==1);
    CPPUNIT_ASSERT(myData1.getNoValues()==3);
    CPPUNIT_ASSERT(myData1.getShape().size()==1);
//     myDataView = myData1.getDefaultValue();
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

    DataTagged myData4;
    myData4.getDataAtOffsetRW(myData4.getDefaultOffset())=7.0;
//     myData4.getDefaultValue()()=7.0;

    myData1.setSlice(&myData4, region);

    //cout << myData3.toString() << endl;

    CPPUNIT_ASSERT(myData1.getTagLookup().size()==0);

    CPPUNIT_ASSERT(myData1.getLength()==3);
    CPPUNIT_ASSERT(myData1.getRank()==1);
    CPPUNIT_ASSERT(myData1.getNoValues()==3);
    CPPUNIT_ASSERT(myData1.getShape().size()==1);

//     myDataView = myData1.getDefaultValue();
    offset=myData1.getDefaultOffset();
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,0)==7.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,1)==6.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,2)==5.0);

  }

  {

    cout << "\tTest slicing DataTagged with rank 3 default value only." << endl;

    DataTagged::TagListType keys;

    DataTagged::ValueBatchType values;

    DataTypes::ShapeType viewShape;
    viewShape.push_back(3);
    viewShape.push_back(3);
    viewShape.push_back(3);

    DataTypes::ValueType viewData1(27);
    for (int i=0;i<viewData1.size();i++) {
      viewData1[i]=i;
    }
//     DataArrayView myView1(viewData1,viewShape);
    DataTagged myData1(FunctionSpace(),viewShape,viewData1);

    DataTypes::ValueType viewData2(27);
    for (int i=0;i<viewData2.size();i++) {
      viewData2[i]=i+27;
    }
//     DataArrayView myView2(viewData2,viewShape);
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

    //cout << myData1.toString() << endl;

    CPPUNIT_ASSERT(myData1.getTagLookup().size()==0);

    CPPUNIT_ASSERT(myData1.getLength()==27);
    CPPUNIT_ASSERT(myData1.getRank()==3);
    CPPUNIT_ASSERT(myData1.getNoValues()==27);
    CPPUNIT_ASSERT(myData1.getShape().size()==3);

//     DataArrayView myDataView = myData1.getDefaultValue();
    int offset=myData1.getDefaultOffset();
    CPPUNIT_ASSERT(offset==0);


    // rank 1 slice

    viewShape.clear();
    viewShape.push_back(3);

    DataTypes::ValueType viewData3(3);
    for (int i=0;i<viewData3.size();i++) {
      viewData3[i]=i+60;
    }
//     DataArrayView myView3(viewData3,viewShape);
    DataTagged myData3(FunctionSpace(),viewShape,viewData3);

    region.clear();
    region.push_back(region_element);
    region_element.first=0;
    region_element.second=0;
    region.push_back(region_element);
    region.push_back(region_element);

    myData1.setSlice(&myData3, region);

    //cout << myData1.toString() << endl;

    CPPUNIT_ASSERT(myData1.getTagLookup().size()==0);

    CPPUNIT_ASSERT(myData1.getLength()==27);
    CPPUNIT_ASSERT(myData1.getRank()==3);
    CPPUNIT_ASSERT(myData1.getNoValues()==27);
    CPPUNIT_ASSERT(myData1.getShape().size()==3);
//     myDataView = myData1.getDefaultValue();
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

    DataTagged myData4;
    myData4.getDataAtOffsetRW(myData4.getDefaultOffset())=70.0;

    myData1.setSlice(&myData4, region);

    //cout << myData1.toString() << endl;

    CPPUNIT_ASSERT(myData1.getTagLookup().size()==0);

    CPPUNIT_ASSERT(myData1.getLength()==27);

    CPPUNIT_ASSERT(myData1.getRank()==3);
    CPPUNIT_ASSERT(myData1.getNoValues()==27);
    CPPUNIT_ASSERT(myData1.getShape().size()==3);

//     myDataView = myData1.getDefaultValue();
    offset=myData1.getDefaultOffset();
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,0,0,0)==70.0);

  }

  {

    cout << "\tTest slicing DataTagged with scalar values and one tag." << endl;

    DataTagged::TagListType keys;
    keys.push_back(1);

    DataTagged::ValueBatchType values;

    DataTypes::ShapeType viewShape;

    // default value for Data1
    DataTypes::ValueType viewData1(1*2);
    viewData1[0]=0.0;
//     DataArrayView myView1(viewData1,viewShape);

    // value for tag "1" for Data1
//     DataTypes::ValueType viewData2(1);
    viewData1[1]=0.0;
//     DataArrayView myView2(viewData2,viewShape);
//     values.push_back(myView2);

    DataTagged myData1(FunctionSpace(),viewShape,keys,viewData1);

    values.clear();

    // default value for Data2
    DataTypes::ValueType viewData3(1*2);
    viewData3[0]=1.0;
//     DataArrayView myView3(viewData3,viewShape);

    // value for tag "1" for Data2
//     DataTypes::ValueType viewData4(1);
    viewData3[1]=2.0;
//     DataArrayView myView4(viewData4,viewShape);
//     values.push_back(myView4);

    DataTagged myData2(FunctionSpace(),viewShape,keys,viewData3);

    // full slice

    DataTypes::RegionType region;

    myData1.setSlice(&myData2, region);

    //cout << myData1.toString() << endl;

    CPPUNIT_ASSERT(myData1.getTagLookup().size()==1);

    CPPUNIT_ASSERT(myData1.getLength()==2);
    CPPUNIT_ASSERT(myData1.getRank()==0);
    CPPUNIT_ASSERT(myData1.getNoValues()==1);
    CPPUNIT_ASSERT(myData1.getShape().size()==0);
//     DataArrayView myDataView = myData1.getDefaultValue();
    int offset=myData1.getDefaultOffset();
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(myData1.getVectorRO()[offset]==1.0);

//     myDataView = myData1.getDataPointByTag(1);
    offset=myData1.getOffsetForTag(1);
    CPPUNIT_ASSERT(offset==1);
    CPPUNIT_ASSERT(myData1.getVectorRO()[offset]==2.0);
  }

  {

    cout << "\tTest slicing DataTagged with rank 1 values and one tag." << endl;
    DataTagged::TagListType keys;
    keys.push_back(1);

    DataTagged::ValueBatchType values;

    DataTypes::ShapeType viewShape;
    viewShape.push_back(3);

    int nvals=3;
    // default value for Data1
    DataTypes::ValueType viewData1(3*2);
    for (int i=0;i<nvals;i++) {
      viewData1[i]=0.0;
    }
//     DataArrayView myView1(viewData1,viewShape);

    // value for tag "1" for Data1
//     DataTypes::ValueType viewData2(3);
    for (int i=0;i<nvals;i++) {
      viewData1[nvals+i]=0.0;
    }
//     DataArrayView myView2(viewData2,viewShape);
//     values.push_back(myView2);

    DataTagged myData1(FunctionSpace(),viewShape,keys,viewData1);
    values.clear();

    // default value for Data2
    DataTypes::ValueType viewData3(3*2);
    for (int i=0;i<nvals;i++) {
      viewData3[i]=1.0;
    }
//     DataArrayView myView3(viewData3,viewShape);

    // value for tag "1" for Data2
//     DataTypes::ValueType viewData4(3);
    for (int i=0;i<nvals;i++) {
      viewData3[nvals+i]=2.0;
    }
//     DataArrayView myView4(viewData4,viewShape);
//     values.push_back(myView4);

    DataTagged myData2(FunctionSpace(),viewShape,keys,viewData3);

    // full slice

    std::pair<int, int> region_element;
    region_element.first=0;
    region_element.second=3;
    DataTypes::RegionType region;
    region.push_back(region_element);

    myData1.setSlice(&myData2, region);

    //cout << myData1.toString() << endl;

    CPPUNIT_ASSERT(myData1.getTagLookup().size()==1);

    CPPUNIT_ASSERT(myData1.getLength()==6);
    CPPUNIT_ASSERT(myData1.getRank()==1);
    CPPUNIT_ASSERT(myData1.getNoValues()==3);
    CPPUNIT_ASSERT(myData1.getShape().size()==1);

//     DataArrayView myDataView = myData1.getDefaultValue();
    int offset=myData1.getDefaultOffset();
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,0)==1.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,1)==1.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,2)==1.0);

//     myDataView = myData1.getDataPointByTag(1);
    offset=myData1.getOffsetForTag(1);
    CPPUNIT_ASSERT(offset==3);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,0)==2.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,1)==2.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,2)==2.0);

    // rank 1 slice

    viewShape.clear();
    viewShape.push_back(1);

    DataTypes::ValueType viewData5(1*2);
    viewData5[0]=3.0;
//     DataArrayView myView5(viewData5,viewShape);

    values.clear();

//     DataTypes::ValueType viewData6(1);
    viewData5[1]=4.0;
//     DataArrayView myView6(viewData6,viewShape);
//     values.push_back(myView6);

    DataTagged myData3(FunctionSpace(),viewShape,keys,viewData5);

    region.clear();
    region_element.first=1;
    region_element.second=2;
    region.push_back(region_element);

    myData1.setSlice(&myData3, region);

    //cout << myData1.toString() << endl;

    CPPUNIT_ASSERT(myData1.getTagLookup().size()==1);

    CPPUNIT_ASSERT(myData1.getLength()==6);
    CPPUNIT_ASSERT(myData1.getRank()==1);
    CPPUNIT_ASSERT(myData1.getNoValues()==3);
    CPPUNIT_ASSERT(myData1.getShape().size()==1);
//     myDataView = myData1.getDefaultValue();
    offset=myData1.getDefaultOffset();
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,0)==1.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,1)==3.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,2)==1.0);

//     myDataView = myData1.getDataPointByTag(1);
    offset=myData1.getOffsetForTag(1);
    CPPUNIT_ASSERT(offset==3);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,0)==2.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,1)==4.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,2)==2.0);

    // scalar slice

    viewShape.clear();

    DataTypes::ValueType viewData7(1*2);
    viewData7[0]=5.0;
//     DataArrayView myView7(viewData7,viewShape);

    values.clear();

//     DataTypes::ValueType viewData8(1);
    viewData7[1]=6.0;
//     DataArrayView myView8(viewData8,viewShape);
//     values.push_back(myView8);

    DataTagged myData4(FunctionSpace(),viewShape,keys,viewData7);

    region.clear();
    region_element.first=0;
    region_element.second=0;
    region.push_back(region_element);

    myData1.setSlice(&myData4, region);

    //cout << myData1.toString() << endl;
    CPPUNIT_ASSERT(myData1.getRank()==1);
    CPPUNIT_ASSERT(myData1.getNoValues()==3);
    CPPUNIT_ASSERT(myData1.getShape().size()==1);

//     myDataView = myData1.getDefaultValue();
    offset=myData1.getDefaultOffset();
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,0)==5.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,1)==3.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,2)==1.0);

//     myDataView = myData1.getDataPointByTag(1);
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

    DataTagged::ValueBatchType values;

    DataTypes::ShapeType viewShape;
    viewShape.push_back(3);
    viewShape.push_back(3);
    viewShape.push_back(3);

    int nvals=27;
    // default value for Data1
    DataTypes::ValueType viewData1(27*2);
    for (int i=0;i<nvals;i++) {
      viewData1[i]=0.0;
    }
//     DataArrayView myView1(viewData1,viewShape);

    // value for tag "1" for Data1
//     DataTypes::ValueType viewData2(27);
    for (int i=0;i<nvals;i++) {
      viewData1[nvals+i]=0.0;
    }
//     DataArrayView myView2(viewData2,viewShape);
//     values.push_back(myView2);

    DataTagged myData1(FunctionSpace(),viewShape,keys,viewData1);

    values.clear();

    // default value for Data2
    DataTypes::ValueType viewData3(27*2);
    for (int i=0;i<nvals;i++) {
      viewData3[i]=1.0;
    }
//     DataArrayView myView3(viewData3,viewShape);

    // value for tag "1" for Data2
//     DataTypes::ValueType viewData4(27);
    for (int i=0;i<nvals;i++) {
      viewData3[nvals+i]=2.0;
    }
//     DataArrayView myView4(viewData4,viewShape);
//     values.push_back(myView4);

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

    //cout << myData1.toString() << endl;

    CPPUNIT_ASSERT(myData1.getTagLookup().size()==1);

    CPPUNIT_ASSERT(myData1.getLength()==54);

    CPPUNIT_ASSERT(myData1.getRank()==3);
    CPPUNIT_ASSERT(myData1.getNoValues()==27);
    CPPUNIT_ASSERT(myData1.getShape().size()==3);


    int offset=myData1.getDefaultOffset();
//     DataArrayView myDataView = myData1.getDefaultValue();
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,0,0,0)==1.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,1,1,1)==1.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,2,2,2)==1.0);

//     myDataView = myData1.getDataPointByTag(1);
    offset=myData1.getOffsetForTag(1);
    CPPUNIT_ASSERT(offset==27);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,0,0,0)==2.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,1,1,1)==2.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,2,2,2)==2.0);

    // rank 1 slice

    viewShape.clear();
    viewShape.push_back(3);
  
    nvals=3;

    DataTypes::ValueType viewData5(3*2);
    for (int i=0;i<nvals;i++) {
      viewData5[i]=3.0;
    }
//     DataArrayView myView5(viewData5,viewShape);

    values.clear();

//     DataTypes::ValueType viewData6(3);
    for (int i=0;i<nvals;i++) {
      viewData5[nvals+i]=4.0;
    }
//     DataArrayView myView6(viewData6,viewShape);
//     values.push_back(myView6);

    DataTagged myData3(FunctionSpace(),viewShape,keys,viewData5);

    region.clear();
    region.push_back(region_element);
    region_element.first=0;
    region_element.second=0;
    region.push_back(region_element);
    region.push_back(region_element);

    myData1.setSlice(&myData3, region);

    //cout << myData1.toString() << endl;

    CPPUNIT_ASSERT(myData1.getTagLookup().size()==1);

    CPPUNIT_ASSERT(myData1.getLength()==54);

    CPPUNIT_ASSERT(myData1.getRank()==3);
    CPPUNIT_ASSERT(myData1.getNoValues()==27);
    CPPUNIT_ASSERT(myData1.getShape().size()==3);


    offset=myData1.getDefaultOffset();
//     myDataView = myData1.getDefslicing DataTagged with rank 3 values and one tagaultValue();
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,0,0,0)==3.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,1,0,0)==3.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,2,0,0)==3.0);

//     myDataView = myData1.getDataPointByTag(1);
    offset=myData1.getOffsetForTag(1);
    CPPUNIT_ASSERT(offset==27);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,0,0,0)==4.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,1,0,0)==4.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,2,0,0)==4.0);

    // scalar slice

    viewShape.clear();

    DataTypes::ValueType viewData7(1*2);
    viewData7[0]=5.0;
//     DataArrayView myView7(viewData7,viewShape);

    values.clear();

//     DataTypes::ValueType viewData8(1);
    viewData7[1]=6.0;
//     DataArrayView myView8(viewData8,viewShape);
//     values.push_back(myView8);

    DataTagged myData4(FunctionSpace(),viewShape,keys,viewData7);

    region.clear();
    region_element.first=0;
    region_element.second=0;
    region.push_back(region_element);
    region.push_back(region_element);
    region.push_back(region_element);

    myData1.setSlice(&myData4, region);

    //cout << myData1.toString() << endl;

    CPPUNIT_ASSERT(myData1.getRank()==3);
    CPPUNIT_ASSERT(myData1.getNoValues()==27);
    CPPUNIT_ASSERT(myData1.getShape().size()==3);


    offset=myData1.getDefaultOffset();
//     myDataView = myData1.getDefaultValue();
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,0,0,0)==5.0);

//     myDataView = myData1.getDataPointByTag(1);
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

    DataTagged::ValueBatchType values;

    DataTypes::ShapeType viewShape;

    // default value for Data1
    DataTypes::ValueType viewData1(1*4);
    viewData1[0]=0.0;
//     DataArrayView myView1(viewData1,viewShape);

    // value for tag "1" for Data1
//     DataTypes::ValueType viewData2(1);
    viewData1[1]=0.0;
//     DataArrayView myView2(viewData2,viewShape);
//     values.push_back(myView2);

    // value for tag "2" for Data1
//     DataTypes::ValueType viewData5(1);
    viewData1[2]=0.0;
//     DataArrayView myView5(viewData5,viewShape);
//     values.push_back(myView5);

    // value for tag "3" for Data1
//     DataTypes::ValueType viewData6(1);
    viewData1[3]=0.0;
//     DataArrayView myView6(viewData6,viewShape);
//     values.push_back(myView6);

    DataTagged myData1(FunctionSpace(),viewShape,keys,viewData1);

    values.clear();

    // default value for Data2
    DataTypes::ValueType viewData3(1*4);
    viewData3[0]=1.0;
//     DataArrayView myView3(viewData3,viewShape);

    // value for tag "1" for Data2
//     DataTypes::ValueType viewData4(1);
    viewData3[1]=2.0;
//     DataArrayView myView4(viewData4,viewShape);
//     values.push_back(myView4);

    // value for tag "2" for Data2
//     DataTypes::ValueType viewData7(1);
    viewData3[2]=3.0;
//     DataArrayView myView7(viewData7,viewShape);
//     values.push_back(myView7);

    // value for tag "3" for Data2
//     DataTypes::ValueType viewData8(1);
    viewData3[3]=4.0;
//     DataArrayView myView8(viewData8,viewShape);
//     values.push_back(myView8);

    DataTagged myData2(FunctionSpace(),viewShape,keys,viewData3);

    // full slice

    DataTypes::RegionType region;

    myData1.setSlice(&myData2, region);

    //cout << myData1.toString() << endl;

    CPPUNIT_ASSERT(myData1.getTagLookup().size()==3);

    CPPUNIT_ASSERT(myData1.getLength()==4);

    CPPUNIT_ASSERT(myData1.getRank()==0);
    CPPUNIT_ASSERT(myData1.getNoValues()==1);
    CPPUNIT_ASSERT(myData1.getShape().size()==0);

    int offset=myData1.getDefaultOffset();
//     DataArrayView myDataView = myData1.getDefaultValue();
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(myData1.getDataAtOffsetRO(offset)==1.0);

    offset=myData1.getOffsetForTag(1);
//     myDataView = myData1.getDataPointByTag(1);
    CPPUNIT_ASSERT(offset==1);
    CPPUNIT_ASSERT(myData1.getDataAtOffsetRO(offset)==2.0);

    offset=myData1.getOffsetForTag(2);
//     myDataView = myData1.getDataPointByTag(2);
    CPPUNIT_ASSERT(offset==2);

    CPPUNIT_ASSERT(myData1.getDataAtOffsetRO(offset)==3.0);

//     myDataView = myData1.getDataPointByTag(3);
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

    DataTagged::ValueBatchType values;

    DataTypes::ShapeType viewShape;
    viewShape.push_back(3);

    int nvals=3;

    // default value for Data1
    DataTypes::ValueType viewData1(3*4);
    for (int i=0;i<viewData1.size();i++) {
      viewData1[i]=0.0;
    }
//     DataArrayView myView1(viewData1,viewShape);

    // value for tag "1" for Data1
//     DataTypes::ValueType viewData2(3);
    for (int i=0;i<nvals;i++) {
      viewData1[nvals+i]=0.0;
    }
//     DataArrayView myView2(viewData2,viewShape);
//     values.push_back(myView2);

    // value for tag "2" for Data1
//     DataTypes::ValueType viewData3(3);
    for (int i=0;i<nvals;i++) {
      viewData1[2*nvals+i]=0.0;
    }
//     DataArrayView myView3(viewData3,viewShape);
//     values.push_back(myView3);

    // value for tag "3" for Data1
//     DataTypes::ValueType viewData4(3);
    for (int i=0;i<nvals;i++) {
      viewData1[3*nvals+i]=0.0;
    }
//     DataArrayView myView4(viewData4,viewShape);
//     values.push_back(myView4);

    DataTagged myData1(FunctionSpace(),viewShape,keys,viewData1);

    values.clear();

    nvals=3;

    // default value for Data2
    DataTypes::ValueType viewData5(3*4);
    for (int i=0;i<nvals;i++) {
      viewData5[i]=1.0;
    }
//     DataArrayView myView5(viewData5,viewShape);

    // value for tag "1" for Data2
//     DataTypes::ValueType viewData6(3);
    for (int i=0;i<nvals;i++) {
      viewData5[nvals+i]=2.0;
    }
//     DataArrayView myView6(viewData6,viewShape);
//     values.push_back(myView6);

    // value for tag "2" for Data2
//     DataTypes::ValueType viewData7(3);
    for (int i=0;i<nvals;i++) {
      viewData5[2*nvals+i]=3.0;
    }
//     DataArrayView myView7(viewData7,viewShape);
//     values.push_back(myView7);

    // value for tag "3" for Data2
//     DataTypes::ValueType viewData8(3);
    for (int i=0;i<nvals;i++) {
      viewData5[3*nvals+i]=4.0;
    }
//     DataArrayView myView8(viewData8,viewShape);
//     values.push_back(myView8);

    DataTagged myData2(FunctionSpace(),viewShape,keys,viewData5);

    // full slice

    std::pair<int, int> region_element;
    region_element.first=0;
    region_element.second=3;
    DataTypes::RegionType region;
    region.push_back(region_element);

    myData1.setSlice(&myData2, region);

    //cout << myData1.toString() << endl;

    CPPUNIT_ASSERT(myData1.getTagLookup().size()==3);

    CPPUNIT_ASSERT(myData1.getLength()==12);

    CPPUNIT_ASSERT(myData1.getRank()==1);
    CPPUNIT_ASSERT(myData1.getNoValues()==3);
    CPPUNIT_ASSERT(myData1.getShape().size()==1);


//     DataArrayView myDataView = myData1.getDefaultValue();
    int offset=myData1.getDefaultOffset();
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,0)==1.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,1)==1.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,2)==1.0);

//     myDataView = myData1.getDataPointByTag(1);
    offset=myData1.getOffsetForTag(1);
    CPPUNIT_ASSERT(offset==3);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,0)==2.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,1)==2.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,2)==2.0);

//     myDataView = myData1.getDataPointByTag(2);
    offset=myData1.getOffsetForTag(2);
    CPPUNIT_ASSERT(offset==6);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,0)==3.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,1)==3.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,2)==3.0);

//     myDataView = myData1.getDataPointByTag(3);
    offset=myData1.getOffsetForTag(3);
    CPPUNIT_ASSERT(offset==9);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,0)==4.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,1)==4.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,2)==4.0);

    // rank 1 slice

    viewShape.clear();
    viewShape.push_back(1);

    DataTypes::ValueType viewData9(1*4);
    viewData9[0]=6.0;
//     DataArrayView myView9(viewData9,viewShape);

    values.clear();

//     DataTypes::ValueType viewData10(1);
    viewData9[1]=7.0;
//     DataArrayView myView10(viewData10,viewShape);
//     values.push_back(myView10);

//     DataTypes::ValueType viewData11(1);
    viewData9[2]=8.0;
//     DataArrayView myView11(viewData11,viewShape);
//     values.push_back(myView11);

//     DataTypes::ValueType viewData12(1);
    viewData9[3]=9.0;
//     DataArrayView myView12(viewData12,viewShape);
//     values.push_back(myView12);

    DataTagged myData3(FunctionSpace(),viewShape, keys, viewData9);

    region.clear();
    region_element.first=1;
    region_element.second=2;
    region.push_back(region_element);

    myData1.setSlice(&myData3, region);

    //cout << myData1.toString() << endl;

    CPPUNIT_ASSERT(myData1.getTagLookup().size()==3);

    CPPUNIT_ASSERT(myData1.getLength()==12);

    CPPUNIT_ASSERT(myData1.getRank()==1);
    CPPUNIT_ASSERT(myData1.getNoValues()==3);
    CPPUNIT_ASSERT(myData1.getShape().size()==1);


    offset=myData1.getDefaultOffset();
//     myDataView = myData1.getDefaultValue();
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,0)==1.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,1)==6.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,2)==1.0);

//     myDataView = myData1.getDataPointByTag(1);
    offset=myData1.getOffsetForTag(1);
    CPPUNIT_ASSERT(offset==3);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,0)==2.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,1)==7.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,2)==2.0);

//     myDataView = myData1.getDataPointByTag(2);
    offset=myData1.getOffsetForTag(2);
    CPPUNIT_ASSERT(offset==6);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,0)==3.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,1)==8.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,2)==3.0);

//     myDataView = myData1.getDataPointByTag(3);
    offset=myData1.getOffsetForTag(3);
    CPPUNIT_ASSERT(offset==9);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,0)==4.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,1)==9.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,2)==4.0);

    // scalar slice

    viewShape.clear();

    DataTypes::ValueType viewData13(1*4);
    viewData13[0]=10.0;
//     DataArrayView myView13(viewData13,viewShape);

    values.clear();

//     DataTypes::ValueType viewData14(1);
    viewData13[1]=11.0;
//     DataArrayView myView14(viewData14,viewShape);
//     values.push_back(myView14);

//     DataTypes::ValueType viewData15(2);
    viewData13[2]=12.0;
//     DataArrayView myView15(viewData15,viewShape);
//     values.push_back(myView15);

//     DataTypes::ValueType viewData16(3);
    viewData13[3]=13.0;
//     DataArrayView myView16(viewData16,viewShape);
//     values.push_back(myView16);

    DataTagged myData4(FunctionSpace(),viewShape,keys,viewData13);

    region.clear();
    region_element.first=0;
    region_element.second=0;
    region.push_back(region_element);

    myData1.setSlice(&myData4, region);

    //cout << myData1.toString() << endl;

    CPPUNIT_ASSERT(myData1.getRank()==1);
    CPPUNIT_ASSERT(myData1.getNoValues()==3);
    CPPUNIT_ASSERT(myData1.getShape().size()==1);


    offset=myData1.getDefaultOffset();
/*    myDataView = myData1.getDefaultValue();*/
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,0)==10.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,1)==6.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,2)==1.0);

//     myDataView = myData1.getDataPointByTag(1);
    offset=myData1.getOffsetForTag(1);
    CPPUNIT_ASSERT(offset==3);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,0)==11.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,1)==7.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,2)==2.0);

//     myDataView = myData1.getDataPointByTag(2);
    offset=myData1.getOffsetForTag(2);
    CPPUNIT_ASSERT(offset==6);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,0)==12.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,1)==8.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,2)==3.0);

//     myDataView = myData1.getDataPointByTag(3);
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

    DataTagged::ValueBatchType values;

    DataTypes::ShapeType viewShape;
    viewShape.push_back(3);
    viewShape.push_back(3);
    viewShape.push_back(3);
    int nvals=noValues(viewShape);

    // default value for Data1
    DataTypes::ValueType viewData1(27*4);
    for (int i=0;i<nvals;i++) {
      viewData1[i]=0.0;
    }
//     DataArrayView myView1(viewData1,viewShape);

    // value for tag "1" for Data1
//     DataTypes::ValueType viewData2(27);
    for (int i=0;i<nvals;i++) {
      viewData1[nvals+i]=0.0;
    }
//     DataArrayView myView2(viewData2,viewShape);
//     values.push_back(myView2);

    // value for tag "2" for Data1
//     DataTypes::ValueType viewData3(27);
    for (int i=0;i<nvals;i++) {
      viewData1[2*nvals+i]=0.0;
    }
//     DataArrayView myView3(viewData3,viewShape);
//     values.push_back(myView3);

    // value for tag "3" for Data1
//     DataTypes::ValueType viewData4(27);
    for (int i=0;i<nvals;i++) {
      viewData1[3*nvals+i]=0.0;
    }
//     DataArrayView myView4(viewData4,viewShape);
//     values.push_back(myView4);

    DataTagged myData1(FunctionSpace(),viewShape,keys,viewData1);

    values.clear();

    // default value for Data2
    DataTypes::ValueType viewData5(27*4);
    for (int i=0;i<nvals;i++) {
      viewData5[i]=1.0;
    }
//     DataArrayView myView5(viewData5,viewShape);

    // value for tag "1" for Data2
//     DataTypes::ValueType viewData6(27);
    for (int i=0;i<nvals;i++) {
      viewData5[nvals+i]=2.0;
    }
//     DataArrayView myView6(viewData6,viewShape);
//     values.push_back(myView6);

    // value for tag "2" for Data2
//     DataTypes::ValueType viewData7(27);
    for (int i=0;i<nvals;i++) {
      viewData5[2*nvals+i]=3.0;
    }
//     DataArrayView myView7(viewData7,viewShape);
//     values.push_back(myView7);

    // value for tag "3" for Data2
    DataTypes::ValueType viewData8(27);
    for (int i=0;i<nvals;i++) {
      viewData5[3*nvals+i]=4.0;
    }
//     DataArrayView myView8(viewData8,viewShape);
//     values.push_back(myView8);

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

    //cout << myData1.toString() << endl;

    CPPUNIT_ASSERT(myData1.getTagLookup().size()==3);

    CPPUNIT_ASSERT(myData1.getLength()==108);

    CPPUNIT_ASSERT(myData1.getRank()==3);
    CPPUNIT_ASSERT(myData1.getNoValues()==27);
    CPPUNIT_ASSERT(myData1.getShape().size()==3);

    int offset=myData1.getDefaultOffset();
//     DataArrayView myDataView = myData1.getDefaultValue();
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,0,0,0)==1.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,1,0,0)==1.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,2,0,0)==1.0);

//     myDataView = myData1.getDataPointByTag(1);
    offset=myData1.getOffsetForTag(1);
    CPPUNIT_ASSERT(offset==27);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,0,0,0)==2.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,1,0,0)==2.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,2,0,0)==2.0);

//     myDataView = myData1.getDataPointByTag(2);
    offset=myData1.getOffsetForTag(2);
    CPPUNIT_ASSERT(offset==54);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,0,0,0)==3.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,1,0,0)==3.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,2,0,0)==3.0);

//     myDataView = myData1.getDataPointByTag(3);
    offset=myData1.getOffsetForTag(3);
    CPPUNIT_ASSERT(offset==81);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,0,0,0)==4.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,1,0,0)==4.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,2,0,0)==4.0);

    // rank 1 slice

    viewShape.clear();
    viewShape.push_back(3);

    nvals=3;

    DataTypes::ValueType viewData9(3*4);
    for (int i=0;i<nvals;i++) {
      viewData9[i]=6.0;
    }
//     DataArrayView myView9(viewData9,viewShape);

    values.clear();

//     DataTypes::ValueType viewData10(3);
    for (int i=0;i<nvals;i++) {
      viewData9[nvals+i]=7.0;
    }
//     DataArrayView myView10(viewData10,viewShape);
//     values.push_back(myView10);

//     DataTypes::ValueType viewData11(3);
    for (int i=0;i<nvals;i++) {
      viewData9[2*nvals+i]=8.0;
    }
//     DataArrayView myView11(viewData11,viewShape);
//     values.push_back(myView11);

//     DataTypes::ValueType viewData12(3);
    for (int i=0;i<nvals;i++) {
      viewData9[3*nvals+i]=9.0;
    }
//     DataArrayView myView12(viewData12,viewShape);
//     values.push_back(myView12);

    DataTagged myData3(FunctionSpace(),viewShape,keys,viewData9);

    region.clear();
    region_element.first=0;
    region_element.second=3;
    region.push_back(region_element);
    region_element.second=0;
    region.push_back(region_element);
    region.push_back(region_element);

    myData1.setSlice(&myData3, region);

    //cout << myData1.toString() << endl;

    CPPUNIT_ASSERT(myData1.getTagLookup().size()==3);

    CPPUNIT_ASSERT(myData1.getLength()==108);

    CPPUNIT_ASSERT(myData1.getRank()==3);
    CPPUNIT_ASSERT(myData1.getNoValues()==27);
    CPPUNIT_ASSERT(myData1.getShape().size()==3);

//    myDataView = myData1.getDefaultValue();
//     CPPUNIT_ASSERT(!myDataView.isEmpty());
    offset=myData1.getDefaultOffset();
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,0,0,0)==6.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,1,0,0)==6.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,2,0,0)==6.0);

//     myDataView = myData1.getDataPointByTag(1);
    offset=myData1.getOffsetForTag(1);
    CPPUNIT_ASSERT(offset==27);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,0,0,0)==7.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,1,0,0)==7.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,2,0,0)==7.0);

    offset=myData1.getOffsetForTag(2);
//     myDataView = myData1.getDataPointByTag(2);
//     CPPUNIT_ASSERT(!myDataView.isEmpty());
    CPPUNIT_ASSERT(offset==54);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,0,0,0)==8.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,1,0,0)==8.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,2,0,0)==8.0);

//     myDataView = myData1.getDataPointByTag(3);
    offset=myData1.getOffsetForTag(3);
    CPPUNIT_ASSERT(offset==81);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,0,0,0)==9.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,1,0,0)==9.0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,2,0,0)==9.0);

    // scalar slice

    viewShape.clear();

    DataTypes::ValueType viewData13(1*4);
    viewData13[0]=10.0;
//     DataArrayView myView13(viewData13,viewShape);

    values.clear();

//     DataTypes::ValueType viewData14(1);
    viewData13[1]=11.0;
//     DataArrayView myView14(viewData14,viewShape);
//     values.push_back(myView14);

//     DataTypes::ValueType viewData15(2);
    viewData13[2]=12.0;
//     DataArrayView myView15(viewData15,viewShape);
//     values.push_back(myView15);

//     DataTypes::ValueType viewData16(3);
    viewData13[3]=13.0;
//     DataArrayView myView16(viewData16,viewShape);
//     values.push_back(myView16);

    DataTagged myData4(FunctionSpace(),viewShape,keys,viewData13);

    region.clear();
    region_element.first=0;
    region_element.second=0;
    region.push_back(region_element);
    region.push_back(region_element);
    region.push_back(region_element);

    myData1.setSlice(&myData4, region);

    //cout << myData1.toString() << endl;

//     myDataView = myData1.getDefaultValue();
    offset=myData1.getDefaultOffset();
    CPPUNIT_ASSERT(offset==0);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,0,0,0)==10.0);
//     myDataView = myData1.getDataPointByTag(1);
    offset=myData1.getOffsetForTag(1);
    CPPUNIT_ASSERT(offset==27);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,0,0,0)==11.0);

//     myDataView = myData1.getDataPointByTag(2);
    offset=myData1.getOffsetForTag(2);
    CPPUNIT_ASSERT(offset==54);
    CPPUNIT_ASSERT(getRefRO(myData1,offset,0,0,0)==12.0);

//     myDataView = myData1.getDataPointByTag(3);
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

    DataTagged::ValueBatchType values;

    DataTypes::ShapeType viewShape;

    // default value for Data1
    DataTypes::ValueType viewData1(1*4);
    viewData1[0]=0.0;
//     DataArrayView myView1(viewData1,viewShape);

    // value for tag "1" for Data1
//     DataTypes::ValueType viewData2(1);
    viewData1[1]=0.0;
//     DataArrayView myView2(viewData2,viewShape);
//     values.push_back(myView2);

    // value for tag "2" for Data1
//     DataTypes::ValueType viewData5(1);
    viewData1[2]=0.0;
//     DataArrayView myView5(viewData5,viewShape);
//     values.push_back(myView5);

    // value for tag "3" for Data1
//     DataTypes::ValueType viewData6(1);
    viewData1[3]=0.0;
//     DataArrayView myView6(viewData6,viewShape);
//     values.push_back(myView6);

    DataTagged myData1(FunctionSpace(),viewShape,keys1,viewData1);

    values.clear();

    // default value for Data2
    DataTypes::ValueType viewData3(1*4);
    viewData3[0]=1.0;
//     DataArrayView myView3(viewData3,viewShape);

    // value for tag "3" for Data2
//     DataTypes::ValueType viewData4(1);
    viewData3[1]=2.0;
//     DataArrayView myView4(viewData4,viewShape);
//     values.push_back(myView4);

    // value for tag "4" for Data2
//     DataTypes::ValueType viewData7(1);
    viewData3[2]=3.0;
//     DataArrayView myView7(viewData7,viewShape);
//     values.push_back(myView7);

    // value for tag "5" for Data2
//     DataTypes::ValueType viewData8(1);
    viewData3[3]=4.0;
//     DataArrayView myView8(viewData8,viewShape);
//     values.push_back(myView8);

    DataTagged myData2(FunctionSpace(),viewShape,keys2,viewData3);

    //cout << myData1.toString() << endl;
    //cout << myData2.toString() << endl;

    // full slice

    DataTypes::RegionType region;

    myData1.setSlice(&myData2, region);

    //cout << myData1.toString() << endl;

    CPPUNIT_ASSERT(myData1.getTagLookup().size()==5);

    CPPUNIT_ASSERT(myData1.getLength()==6);

    CPPUNIT_ASSERT(myData1.getRank()==0);
    CPPUNIT_ASSERT(myData1.getNoValues()==1);
    CPPUNIT_ASSERT(myData1.getShape().size()==0);
    CPPUNIT_ASSERT(myData1.getDataAtOffsetRO(0)==1.0);


//     DataArrayView myDataView = myData1.getDefaultValue();
//     CPPUNIT_ASSERT(!myDataView.isEmpty());
    CPPUNIT_ASSERT(myData1.getDefaultOffset()==0);
//     CPPUNIT_ASSERT(myDataView.getOffset()==0);


//     myDataView = myData1.getDataPointByTag(1);
//     CPPUNIT_ASSERT(!myDataView.isEmpty());
    CPPUNIT_ASSERT(myData1.getOffsetForTag(1)==1);
//     CPPUNIT_ASSERT(myDataView.getRank()==0);
//     CPPUNIT_ASSERT(myDataView.noValues()==1);
//     CPPUNIT_ASSERT(myDataView.getShape().size()==0);
    CPPUNIT_ASSERT(myData1.getDataAtOffsetRO(1)==1.0);

//     myDataView = myData1.getDataPointByTag(2);
//     CPPUNIT_ASSERT(!myDataView.isEmpty());
    CPPUNIT_ASSERT(myData1.getOffsetForTag(2)==2);
//     CPPUNIT_ASSERT(myDataView.getRank()==0);
//     CPPUNIT_ASSERT(myDataView.noValues()==1);
//     CPPUNIT_ASSERT(myDataView.getShape().size()==0);
    CPPUNIT_ASSERT(myData1.getDataAtOffsetRO(2)==1.0);

//     myDataView = myData1.getDataPointByTag(3);
//     CPPUNIT_ASSERT(!myDataView.isEmpty());
    CPPUNIT_ASSERT(myData1.getOffsetForTag(3)==3);
/*    CPPUNIT_ASSERT(myDataView.getRank()==0);
    CPPUNIT_ASSERT(myDataView.noValues()==1);
    CPPUNIT_ASSERT(myDataView.getShape().size()==0);*/
    CPPUNIT_ASSERT(myData1.getDataAtOffsetRO(3)==2.0);

//     myDataView = myData1.getDataPointByTag(4);
//     CPPUNIT_ASSERT(!myDataView.isEmpty());
    CPPUNIT_ASSERT(myData1.getOffsetForTag(4)==4);
//     CPPUNIT_ASSERT(myDataView.getRank()==0);
//     CPPUNIT_ASSERT(myDataView.noValues()==1);
//     CPPUNIT_ASSERT(myDataView.getShape().size()==0);
    CPPUNIT_ASSERT(myData1.getDataAtOffsetRO(4)==3.0);

//     myDataView = myData1.getDataPointByTag(5);
//     CPPUNIT_ASSERT(!myDataView.isEmpty());
    CPPUNIT_ASSERT(myData1.getOffsetForTag(5)==5);
/*    CPPUNIT_ASSERT(myData1.getRank()==0);
    CPPUNIT_ASSERT(myData1.noValues()==1);
    CPPUNIT_ASSERT(myData1.getShape().size()==0);*/
    CPPUNIT_ASSERT(myData1.getDataAtOffsetRO(5)==4.0);

  }

}

/*
// Testing to see if FunctionSpaces are checked for taggability before use
void DataTaggedTestCase::testFunctionSpaces()
{
   {
	cout << "\tTest Non-Taggable Degrees Of Freedom." << endl;
	MeshAdapter d;
	FunctionSpace fs=solution(d);
	std::string res=constr(fs);
	if (!res.empty())
	{
		cout << "\t\t" << res << endl;
		CPPUNIT_ASSERT(false);
	}
   }
   {
	cout << "\tTest Non-Taggable Reduced Degrees Of Freedom." << endl;
	MeshAdapter d;
	FunctionSpace fs=reducedSolution(d);
	std::string res=constr(fs);
	if (!res.empty())
	{
		cout << "\t\t" << res << endl;
		CPPUNIT_ASSERT(false);
	}
   }
   {
	cout << "\tTest Non-Taggable Reduced Degrees Of Freedom." << endl;
	MeshAdapter d;
	FunctionSpace fs(d,MeshAdapter::ReducedNodes);
	std::string res=constr(fs);
	if (!res.empty())
	{
		cout << "\t\t" << res << endl;
		CPPUNIT_ASSERT(false);
	}
   }
}*/

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
              "testOperations",&DataTaggedTestCase::testOperations));
//   testSuite->addTest(new TestCaller<DataTaggedTestCase>(
//              "testFunctionSpaces",&DataTaggedTestCase::testFunctionSpaces));
  testSuite->addTest(new TestCaller<DataTaggedTestCase>(
              "testGetSlice",&DataTaggedTestCase::testGetSlice));
  testSuite->addTest(new TestCaller<DataTaggedTestCase>(
              "testSetSlice",&DataTaggedTestCase::testSetSlice));
  return testSuite;
}

