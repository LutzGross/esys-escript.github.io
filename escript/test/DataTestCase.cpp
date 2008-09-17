
/* $Id$ */

/*******************************************************
 *
 *           Copyright 2003-2007 by ACceSS MNRF
 *       Copyright 2007 by University of Queensland
 *
 *                http://esscc.uq.edu.au
 *        Primary Business: Queensland, Australia
 *  Licensed under the Open Software License version 3.0
 *     http://www.opensource.org/licenses/osl-3.0.php
 *
 *******************************************************/

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
using namespace escript::DataTypes;

void DataTestCase::setUp() {
  //
  // This is called before each test is run
}

void DataTestCase::tearDown() {
  //
  // This is called after each test has been run
}



namespace
{

inline
DataTypes::ValueType::reference
getRef(Data& d,int s1, int p1, int x, int y)
{
	return d.getDataAtOffset(d.getDataOffset(s1,p1)+getRelIndex(d.getDataPointShape(),x,y));
}

inline
DataTypes::ValueType::reference
getRef(Data& d, int x, int y)
{
	return d.getDataAtOffset(getRelIndex(d.getDataPointShape(),x,y));
}

}

// This is to test new copy routines, existing tests should remain where they are
void DataTestCase::testCopying()
{
  using namespace escript::DataTypes;
  cout << endl;
  {
	// first we test the deep copy
	cout << "\tTest deep copy DataConstant" << endl;
	DataTypes::ShapeType shape;
	shape.push_back(2);
	shape.push_back(3);
	DataTypes::ValueType data(DataTypes::noValues(shape),1);
	DataConstant* dc=new DataConstant(FunctionSpace(),shape,data);
	Data d(dc);
	
	Data* deep=d.copySelf();	// test self copy
	for (int i=0;i<DataTypes::noValues(shape);++i)
	{
	   if (d.getDataAtOffset(i)!=deep->getDataAtOffset(i))
		assert(false);
	}
	d.setToZero();
	for (int i=0;i<DataTypes::noValues(shape);++i)
	{
	   if (d.getDataAtOffset(i)==deep->getDataAtOffset(i))
		assert(false);
	}
	d.copy(*deep);			// test copy from object
	for (int i=0;i<DataTypes::noValues(shape);++i)
	{
	   if (d.getDataAtOffset(i)!=deep->getDataAtOffset(i))
		assert(false);
	}
	d.setToZero();
	for (int i=0;i<DataTypes::noValues(shape);++i)
	{
	   if (d.getDataAtOffset(i)==deep->getDataAtOffset(i))
		assert(false);
	}
	delete deep;
  }

  {
	// first we test the deep copy
	cout << "\tTest deep copy DataExpanded" << endl;
	DataTypes::ShapeType shape;
	shape.push_back(2);
	shape.push_back(3);
	DataTypes::ValueType data(DataTypes::noValues(shape),1);
	DataExpanded* dc=new DataExpanded(FunctionSpace(),shape,data);
	Data d(dc);
	
	Data* deep=d.copySelf();	// test self copy
	for (int i=0;i<DataTypes::noValues(shape);++i)
	{
	   if (d.getDataAtOffset(i)!=deep->getDataAtOffset(i))
		assert(false);
	}
	d.setToZero();
	for (int i=0;i<DataTypes::noValues(shape);++i)
	{
	   if (d.getDataAtOffset(i)==deep->getDataAtOffset(i))
		assert(false);
	}
	d.copy(*deep);			// test copy from object
	for (int i=0;i<DataTypes::noValues(shape);++i)
	{
	   if (d.getDataAtOffset(i)!=deep->getDataAtOffset(i))
		assert(false);
	}
	d.setToZero();
	for (int i=0;i<DataTypes::noValues(shape);++i)
	{
	   if (d.getDataAtOffset(i)==deep->getDataAtOffset(i))
		assert(false);
	}
	delete deep;
  }
  {
	// first we test the deep copy
	cout << "\tTest deep copy DataTagged" << endl;
	DataTypes::ShapeType shape;
	shape.push_back(2);
	shape.push_back(3);
	DataTypes::ValueType data(DataTypes::noValues(shape),1);
	DataTagged* dc=new DataTagged(FunctionSpace(),shape,data);
	Data d(dc);
	
	Data* deep=d.copySelf();	// test self copy
	for (int i=0;i<DataTypes::noValues(shape);++i)
	{
	   if (d.getDataAtOffset(i)!=deep->getDataAtOffset(i))
		assert(false);
	}
	d.setToZero();
	for (int i=0;i<DataTypes::noValues(shape);++i)
	{
	   if (d.getDataAtOffset(i)==deep->getDataAtOffset(i))
		assert(false);
	}
	d.copy(*deep);			// test copy from object
	for (int i=0;i<DataTypes::noValues(shape);++i)
	{
	   if (d.getDataAtOffset(i)!=deep->getDataAtOffset(i))
		assert(false);
	}
	d.setToZero();
	for (int i=0;i<DataTypes::noValues(shape);++i)
	{
	   if (d.getDataAtOffset(i)==deep->getDataAtOffset(i))
		assert(false);
	}
	delete deep;
  }

}

void DataTestCase::testSlicing() {

  using namespace escript::DataTypes;
  cout << endl;

  {

    cout << "\tTest get-slicing DataConstant" << endl;

    DataTypes::ShapeType viewShape;
    viewShape.push_back(2);
    viewShape.push_back(3);
    Data data(1.3,viewShape,FunctionSpace(),false);

    //cout << data.toString() << endl;

    DataTypes::RegionType region;
    region.push_back(DataTypes::RegionType::value_type(0,0));
    region.push_back(DataTypes::RegionType::value_type(0,0));

    Data slice1(data.getSlice(region));

    //cout << slice1.toString() << endl;

    assert(slice1.getDataPointRank()==0);
    assert(slice1.getDataPoint(0,0)==1.3);

    region.clear();
    region.push_back(DataTypes::RegionType::value_type(0,1));
    region.push_back(DataTypes::RegionType::value_type(0,1));

    Data slice2(data.getSlice(region));

    //cout << slice2.toString() << endl;

    assert(slice2.getDataPointRank()==2);
    int off1=slice2.getDataOffset(0,0);
//     assert(slice2.getDataPoint(0,0)(0,0)==1.3);
    assert(slice2.getDataAtOffset(off1+getRelIndex(slice2.getDataPointShape(),0,0))==1.3);

    region.clear();
    region.push_back(DataTypes::RegionType::value_type(0,1));
    region.push_back(DataTypes::RegionType::value_type(0,2));

    Data slice3(data.getSlice(region));

    //cout << slice3.toString() << endl;

    assert(slice3.getDataPointRank()==2);
    off1=slice3.getDataOffset(0,0);
//     assert(slice3.getDataPoint(0,0)(0,0)==1.3);
//     assert(slice3.getDataPoint(0,0)(0,1)==1.3);
    assert(slice3.getDataAtOffset(off1+getRelIndex(slice3.getDataPointShape(),0,0))==1.3);
    assert(slice3.getDataAtOffset(off1+getRelIndex(slice3.getDataPointShape(),0,1))==1.3);

  }

  {

    cout << "\tTest set-slicing DataConstant" << endl;

    DataTypes::ShapeType viewShape;
    Data source(10.0,viewShape,FunctionSpace(),false);

    //cout << source.toString() << endl;

    viewShape.push_back(2);
    viewShape.push_back(3);
    Data target(1.3,viewShape,FunctionSpace(),false);

    //cout << target.toString() << endl;

    DataTypes::RegionType region;
    region.push_back(DataTypes::RegionType::value_type(0,0));
    region.push_back(DataTypes::RegionType::value_type(0,0));

    target.setSlice(source,region);

    //cout << target.toString() << endl;

    int off1=target.getDataOffset(0,0);
    assert(target.getDataAtOffset(off1+getRelIndex(target.getDataPointShape(),0,0)==source.getDataPoint(0,0)));

  }

  {

    cout << "\tTest get-slicing DataTagged" << endl;
    //
    // create a DataTagged with a default value only

    DataTypes::ShapeType viewShape;
    viewShape.push_back(2);
    viewShape.push_back(3);
    Data data(1.3,viewShape,FunctionSpace(),false);
    data.tag();
    data.getDataAtOffset(data.getDataOffset(0,0)+getRelIndex(viewShape,0,0))=1.0;
    data.getDataAtOffset(data.getDataOffset(0,0)+getRelIndex(viewShape,1,1))=2.0;   
//     data.getDataPoint(0,0)(0,0)=1.0;
//     data.getDataPoint(0,0)(1,1)=2.0;

    //cout << data.toString() << endl;
    //
    // create a scalar slice

    DataTypes::RegionType region;
    region.push_back(DataTypes::RegionType::value_type(0,0));
    region.push_back(DataTypes::RegionType::value_type(0,0));

    Data slice1(data.getSlice(region));
    //cout << slice1.toString() << endl;

    assert(slice1.isTagged());
    assert(slice1.getDataPointRank()==0);
    assert(slice1.getDataPoint(0,0)==1.0);
    //
    // create a rank 2 slice with one value

    region.clear();
    region.push_back(DataTypes::RegionType::value_type(0,1));
    region.push_back(DataTypes::RegionType::value_type(0,1));

    Data slice2(data.getSlice(region));

    //cout << slice2.toString() << endl;

    assert(slice2.isTagged());
    assert(slice2.getDataPointRank()==2);

    assert(slice2.getDataAtOffset(slice2.getDataOffset(0,0)+getRelIndex(slice2.getDataPointShape(),0,0))==1.0);

    //
    // create a rank 2 slice with four values

    region.clear();
    region.push_back(DataTypes::RegionType::value_type(0,2));
    region.push_back(DataTypes::RegionType::value_type(0,2));

    Data slice3(data.getSlice(region));

    //cout << slice3.toString() << endl;

    assert(slice3.isTagged());
    assert(slice3.getDataPointRank()==2);
    assert(getRef(slice3,0,0,0,0)==1.0);
    assert(getRef(slice3,0,0,0,1)==1.3);
    assert(getRef(slice3,0,0,1,0)==1.3);
    assert(getRef(slice3,0,0,1,1)==2.0);

    //
    // add a value for tag "1"

    DataTypes::ValueType viewData(6);
    for (int i=0;i<viewData.size();i++) {
      viewData[i]=i;
    }
//    DataArrayView dataView(viewData,viewShape);

//     data.setTaggedValueFromCPP(1, dataView);
    data.setTaggedValueFromCPP(1, viewShape, viewData);


    //
    // create a full slice

    region.clear();
    region.push_back(DataTypes::RegionType::value_type(0,2));
    region.push_back(DataTypes::RegionType::value_type(0,3));

    Data slice4(data.getSlice(region));

    //cout << slice4.toString() << endl;

    assert(slice4.isTagged());
    assert(slice4.getDataPointRank()==2);
    assert(getRef(slice4,0,0,0,0)==0);
    assert(getRef(slice4,0,0,0,1)==2);
    assert(getRef(slice4,0,0,0,2)==4);
    assert(getRef(slice4,0,0,1,0)==1);
    assert(getRef(slice4,0,0,1,1)==3);
    assert(getRef(slice4,0,0,1,2)==5);

  }

  {

    cout << "\tTest set-slicing DataTagged" << endl;

    //
    // create a source DataTagged with a scalar default value only

    DataTypes::ShapeType viewShape;
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

    DataTypes::RegionType region;
    region.push_back(DataTypes::RegionType::value_type(1,1));
    region.push_back(DataTypes::RegionType::value_type(1,1));

    target.setSlice(source,region);

    //cout << "target:\n" << target.toString() << endl;

    assert(target.isTagged());
    assert(target.getDataPointRank()==2);
    assert(getRef(target,0,0,0,0)==1.3);
    assert(getRef(target,0,0,0,1)==1.3);
    assert(getRef(target,0,0,0,2)==1.3);
    assert(getRef(target,0,0,1,0)==1.3);
    assert(getRef(target,0,0,1,1)==source.getDataPoint(0,0));
    assert(getRef(target,0,0,1,2)==1.3);

    //
    // add a value for tag "1" to target

    DataTypes::ValueType viewData(6);
    for (int i=0;i<viewData.size();i++) {
      viewData[i]=i;
    }
//     DataArrayView dataView(viewData,viewShape);
// 
//     target.setTaggedValueFromCPP(1, dataView);
    target.setTaggedValueFromCPP(1, viewShape, viewData);

    //cout << "target:\n" << target.toString() << endl;

    //
    // set a slice in target from source

    region.clear();
    region.push_back(DataTypes::RegionType::value_type(0,0));
    region.push_back(DataTypes::RegionType::value_type(1,1));

    target.setSlice(source,region);

    //cout << "target:\n" << target.toString() << endl;

    assert(target.isTagged());
    assert(target.getDataPointRank()==2);
    assert(getRef(target,0,0,0,0)==0);
    assert(getRef(target,0,0,0,1)==source.getDataPoint(0,0));
    assert(getRef(target,0,0,0,2)==4);
    assert(getRef(target,0,0,1,0)==1);
    assert(getRef(target,0,0,1,1)==3);
    assert(getRef(target,0,0,1,2)==5);

    //
    // add a value for tag "2" to source

    DataTypes::ShapeType viewShape2;
    DataTypes::ValueType viewData2(1);
    viewData2[0]=6;
//     DataArrayView dataView2(viewData2,viewShape2);
// 
//     source.setTaggedValueFromCPP(2, dataView2);
    source.setTaggedValueFromCPP(2, viewShape2, viewData2);

    //cout << "source:\n" << source.toString() << endl;

    //
    // set a slice in target from source

    region.clear();
    region.push_back(DataTypes::RegionType::value_type(0,0));
    region.push_back(DataTypes::RegionType::value_type(1,1));

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

    DataTypes::ShapeType viewShape;
    viewShape.push_back(2);
    viewShape.push_back(3);
    Data temp(1.3,viewShape,FunctionSpace(),true);

    getRef(temp,0,0,0,0)=0.0;
    getRef(temp,0,0,1,1)=1.0;

    DataTypes::RegionType region;
    region.push_back(DataTypes::RegionType::value_type(0,0));
    region.push_back(DataTypes::RegionType::value_type(0,0));

    Data slice(temp.getSlice(region));

    assert(slice.getDataPointRank()==0);
    assert(slice.getDataPoint(0,0)==0.0);

    region.clear();
    region.push_back(DataTypes::RegionType::value_type(0,1));
    region.push_back(DataTypes::RegionType::value_type(0,1));

    slice=temp.getSlice(region);

    assert(slice.getDataPointRank()==2);
    assert(getRef(slice,0,0,0,0)==0.0);

    region.clear();
    region.push_back(DataTypes::RegionType::value_type(0,2));
    region.push_back(DataTypes::RegionType::value_type(0,2));

    slice=temp.getSlice(region);

    assert(getRef(slice,0,0,0,0)==0.0);
    assert(getRef(slice,0,0,1,1)==1.0);

  }

  {

    cout << "\tTest set-slicing DataExpanded" << endl;

    DataTypes::ShapeType viewShape;
    Data source(10.0,viewShape,FunctionSpace(),true);

    viewShape.push_back(2);
    viewShape.push_back(3);
    Data target(1.3,viewShape,FunctionSpace(),true);

    DataTypes::RegionType region;
    region.push_back(DataTypes::RegionType::value_type(0,0));
    region.push_back(DataTypes::RegionType::value_type(0,0));

    target.setSlice(source,region);

    assert(getRef(target,0,0,0,0)==source.getDataPoint(0,0));

  }

}

void DataTestCase::testAll() {

  cout << endl;

  cout << "\tCreate a Data object from a DataArrayView" << endl;

  DataTypes::ShapeType viewShape;
  viewShape.push_back(3);
  DataTypes::ValueType viewData(3);
  for (int i=0;i<viewShape[0];++i) {
    viewData[i]=i;
  }
//   DataArrayView myView(viewData,viewShape);

  bool expanded=true;
  Data exData(viewData,viewShape,FunctionSpace(),expanded);
//   Data cData(myView);
  Data cData(viewData,viewShape,FunctionSpace());
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

  DataTypes::ShapeType viewShape;
  viewShape.push_back(3);
  DataTypes::ValueType viewData(3);
  for (int i=0;i<viewShape[0];++i) {
    viewData[i]=i;
  }
//   DataArrayView myView(viewData,viewShape);

  bool expanded=true;
  Data exData(viewData,viewShape,FunctionSpace(),expanded);
  Data cData(viewData,viewShape);
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

  DataTypes::ShapeType viewShape;
  viewShape.push_back(2);
  viewShape.push_back(3);
  viewShape.push_back(4);
  DataTypes::ValueType viewData(2*3*4);
  for (int i=0;i<DataTypes::noValues(viewShape);++i) {
    viewData[i]=i;
  }
//   DataArrayView myView(viewData,viewShape);

  Data left(viewData,viewShape);
  Data right(viewData,viewShape);
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

//    DataTagged::TagListType keys;

//    DataTagged::ValueListType values;

    DataTypes::ShapeType viewShape;
    viewShape.push_back(3);

    DataTypes::ValueType viewData(3);
    for (int i=0;i<viewShape[0];i++) {
      viewData[i]=i;
    }
/*    DataArrayView defaultValue(viewData,viewShape);

    bool expanded=false;

    Data myData(keys,values,defaultValue,FunctionSpace(),expanded);*/
    int arr[1]={1};		// iso c++ does not like empty arrays
    DataTagged* dt=new DataTagged(FunctionSpace(),viewShape,arr,viewData); 
    Data myData(dt);

    // cout << myData.toString() << endl;

    assert(!myData.isEmpty());
    assert(myData.isTagged());
    assert(myData.getTagNumber(0)==1);
    assert(myData.getDataPointRank()==1);
    assert(myData.getLength()==3);
    
    assert(myData.getNoValues()==3);
    assert(myData.getDataAtOffset(0)==0.0);
    assert(myData.getDataAtOffset(1)==1.0);
    assert(myData.getDataAtOffset(2)==2.0);

//     DataArrayView myDataView = myData.getPointDataView();
//     assert(!myDataView.isEmpty());
//     assert(myDataView.getOffset()==0);
//     assert(myDataView.getRank()==1);
//     assert(myDataView.noValues()==3);
//     assert(myDataView.getShape().size()==1);
//     assert(myDataView(0)==0.0);
//     assert(myDataView(1)==1.0);
//     assert(myDataView(2)==2.0);

//     myDataView = myData.getDataPoint(0,0);
//     assert(!myDataView.isEmpty());
//     assert(myDataView.getOffset()==0);
//     assert(myDataView.getRank()==1);
//     assert(myDataView.noValues()==3);
//     assert(myDataView.getShape().size()==1);
//     assert(myDataView(0)==0.0);
//     assert(myDataView(1)==1.0);
//     assert(myDataView(2)==2.0);

    double* sampleData=myData.getSampleData(0);
    for (int i=0; i<myData.getNoValues(); i++) {
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
    DataTypes::ValueType eTwoData(viewData);
 //   DataArrayView eTwoView(eTwoData, viewShape);
    for (int i=0;i<viewShape[0];i++) {
      eTwoData[i]=i+2.0;
    }

    myData.setTaggedValueFromCPP(1,viewShape, eTwoData);

    assert(myData.getLength()==6);

    int offset=myData.getDataOffset(0,0);
//    myDataView = myData.getDataPoint(0,0);
//     assert(myDataView==eTwoView);
//     assert(!myDataView.isEmpty());
    assert(offset==3);
    assert(myData.getDataPointRank()==1);
    assert(myData.getNoValues()==3);
//    assert(myDataView.getShape().size()==1);

    assert(myData.getDataAtOffset(offset+0)==2);
    assert(myData.getDataAtOffset(offset+1)==3);
    assert(myData.getDataAtOffset(offset+2)==4);

    sampleData=myData.getSampleDataByTag(1);
    for (int i=0; i<myData.getNoValues(); i++) {
      assert(sampleData[i]==i+2);
    }

  }

  {

    cout << "\tCreate a DataTagged object via tag() method." << endl;

    DataTypes::ShapeType viewShape;
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
//     DataArrayView myDataView = myData.getPointDataView();
    assert(!myData.isEmpty());
//     assert(myDataView.getOffset()==0);
    assert(myData.getDataPointRank()==2);
    assert(myData.getNoValues()==6);
    assert(myData.getDataPointShape().size()==2);
    assert(getRef(myData,0,0)==1.3);
    assert(getRef(myData,0,1)==1.3);
    assert(getRef(myData,0,2)==1.3);
    assert(getRef(myData,1,0)==1.3);
    assert(getRef(myData,1,1)==1.3);
    assert(getRef(myData,1,2)==1.3);

    // check value for data-point (0,0).
//     myDataView = myData.getDataPoint(0,0);
    assert(!myData.isEmpty());
//     assert(myDataView.getOffset()==0);
    assert(myData.getDataPointRank()==2);
    assert(myData.getNoValues()==6);
    assert(myData.getDataPointShape().size()==2);
    assert(getRef(myData,0,0)==1.3);
    assert(getRef(myData,0,1)==1.3);
    assert(getRef(myData,0,2)==1.3);
    assert(getRef(myData,1,0)==1.3);
    assert(getRef(myData,1,1)==1.3);
    assert(getRef(myData,1,2)==1.3);

  }

}

void DataTestCase::testDataTaggedExceptions() {

  cout << endl;

  cout << "\tTest DataTagged exceptions." << endl;

  Data myData;

  try {
      myData.getSampleDataByTag(0);;
      assert(false);
  }
  catch (EsysException&) {
      //cout << e.what() << endl;
      assert(true);
  }

  try {
      myData.setTaggedValueFromCPP(0,DataTypes::ShapeType(), DataTypes::ValueType());;
      assert(false);
  }
  catch (EsysException&) {
      //cout << e.what() << endl;
      assert(true);
  }

}

void DataTestCase::testConstructors() {

  cout << endl;

  DataTypes::ShapeType viewShape;
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
  DataTypes::ShapeType shape;
  shape.push_back(2);
  shape.push_back(3);

  // allocate the data for the DataArrayView
  DataTypes::ValueType data(DataTypes::noValues(shape),0);

  // construct DataArrayView
//   DataArrayView dataView(data,shape);

  // assign values to the data
  for (int i=0;i<shape[0];i++) {
    for (int j=0;j<shape[1];j++) {
      data[getRelIndex(shape,i,j)]=getRelIndex(shape,i,j);
    }
  }

  Data baseEx(data,shape,FunctionSpace(),true);
  Data baseCon(data,shape,FunctionSpace(),false);
  Data baseTag(data,shape,FunctionSpace(),false);
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
      tmp=pow((double)data[getRelIndex(shape,i,j)],(double)3.0);
      assert(std::abs(getRef(resultEx,i,j) - tmp) <= REL_TOL*std::abs(tmp));
      assert(std::abs(getRef(resultCon,i,j) - tmp) <= REL_TOL*std::abs(tmp));
      assert(std::abs(getRef(resultTag,i,j) - tmp) <= REL_TOL*std::abs(tmp));
    }
  }

  cout << "\tTest Data::sin." << endl;
  resultEx.copy(baseEx.sin());
  resultCon.copy(baseCon.sin());
  resultTag.copy(baseTag.sin());
  for (int i=0;i<shape[0];i++) {
    for (int j=0;j<shape[1];j++) {
      tmp=sin((double)data[getRelIndex(shape,i,j)]);
      assert(std::abs(getRef(resultEx,i,j) - tmp) <= REL_TOL*std::abs(tmp));
      assert(std::abs(getRef(resultCon,i,j) - tmp) <= REL_TOL*std::abs(tmp));
      assert(std::abs(getRef(resultTag,i,j) - tmp) <= REL_TOL*std::abs(tmp));
    }
  }

  cout << "\tTest Data::cos." << endl;
  resultEx.copy(baseEx.cos());
  resultCon.copy(baseCon.cos());
  resultTag.copy(baseTag.cos());
  for (int i=0;i<shape[0];i++) {
    for (int j=0;j<shape[1];j++) {
      tmp=cos((double)data[getRelIndex(shape,i,j)]);
      assert(std::abs(getRef(resultEx,i,j) - tmp) <= REL_TOL*std::abs(tmp));
      assert(std::abs(getRef(resultCon,i,j) - tmp) <= REL_TOL*std::abs(tmp));
      assert(std::abs(getRef(resultTag,i,j) - tmp) <= REL_TOL*std::abs(tmp));
    }
  }

  cout << "\tTest Data::tan." << endl;
  resultEx.copy(baseEx.tan());
  resultCon.copy(baseCon.tan());
  resultTag.copy(baseTag.tan());
  for (int i=0;i<shape[0];i++) {
    for (int j=0;j<shape[1];j++) {
      tmp=tan((double)data[getRelIndex(shape,i,j)]);
      assert(std::abs(getRef(resultEx,i,j) - tmp) <= REL_TOL*std::abs(tmp));
      assert(std::abs(getRef(resultCon,i,j)- tmp) <= REL_TOL*std::abs(tmp));
      assert(std::abs(getRef(resultTag,i,j)- tmp) <= REL_TOL*std::abs(tmp));
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
      tmp=atan((double)data[getRelIndex(shape,i,j)]);
      assert(std::abs(getRef(resultEx,i,j) - tmp) <= REL_TOL*std::abs(tmp));
      assert(std::abs(getRef(resultCon,i,j)- tmp) <= REL_TOL*std::abs(tmp));
      assert(std::abs(getRef(resultTag,i,j)- tmp) <= REL_TOL*std::abs(tmp));
    }
  }

  cout << "\tTest Data::sinh." << endl;
  resultEx.copy(baseEx.sinh());
  resultCon.copy(baseCon.sinh());
  resultTag.copy(baseTag.sinh());
  for (int i=0;i<shape[0];i++) {
    for (int j=0;j<shape[1];j++) {
      tmp=sinh((double)data[getRelIndex(shape,i,j)]);
      assert(std::abs(getRef(resultEx,i,j) - tmp) <= REL_TOL*std::abs(tmp));
      assert(std::abs(getRef(resultCon,i,j)- tmp) <= REL_TOL*std::abs(tmp));
      assert(std::abs(getRef(resultTag,i,j)- tmp) <= REL_TOL*std::abs(tmp));
    }
  }

  cout << "\tTest Data::cosh." << endl;
  resultEx.copy(baseEx.cosh());
  resultCon.copy(baseCon.cosh());
  resultTag.copy(baseTag.cosh());
  for (int i=0;i<shape[0];i++) {
    for (int j=0;j<shape[1];j++) {
      tmp=cosh((double)data[getRelIndex(shape,i,j)]);
      assert(std::abs(getRef(resultEx,i,j) - tmp) <= REL_TOL*std::abs(tmp));
      assert(std::abs(getRef(resultCon,i,j)- tmp) <= REL_TOL*std::abs(tmp));
      assert(std::abs(getRef(resultTag,i,j)- tmp) <= REL_TOL*std::abs(tmp));
    }
  }

  cout << "\tTest Data::tanh." << endl;
  resultEx.copy(baseEx.tanh());
  resultCon.copy(baseCon.tanh());
  resultTag.copy(baseTag.tanh());
  for (int i=0;i<shape[0];i++) {
    for (int j=0;j<shape[1];j++) {
      tmp=tanh((double)data[getRelIndex(shape,i,j)]);
      assert(std::abs(getRef(resultEx,i,j) - tmp) <= REL_TOL*std::abs(tmp));
      assert(std::abs(getRef(resultCon,i,j)- tmp) <= REL_TOL*std::abs(tmp));
      assert(std::abs(getRef(resultTag,i,j)- tmp) <= REL_TOL*std::abs(tmp));
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
      tmp=abs((double)data[getRelIndex(shape,i,j)]);
      assert(std::abs(getRef(resultEx,i,j) - tmp) <= REL_TOL*std::abs(tmp));
      assert(std::abs(getRef(resultCon,i,j)- tmp) <= REL_TOL*std::abs(tmp));
      assert(std::abs(getRef(resultTag,i,j)- tmp) <= REL_TOL*std::abs(tmp));
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
      tmp=exp((double)data[getRelIndex(shape,i,j)]);
      assert(std::abs(getRef(resultEx,i,j) - tmp) <= REL_TOL*std::abs(tmp));
      assert(std::abs(getRef(resultCon,i,j)- tmp) <= REL_TOL*std::abs(tmp));
      assert(std::abs(getRef(resultTag,i,j)- tmp) <= REL_TOL*std::abs(tmp));
    }
  }

  cout << "\tTest Data::sqrt." << endl;
  resultEx.copy(baseEx.sqrt());
  resultCon.copy(baseCon.sqrt());
  resultTag.copy(baseTag.sqrt());
  for (int i=0;i<shape[0];i++) {
    for (int j=0;j<shape[1];j++) {
      tmp=sqrt((double)data[getRelIndex(shape,i,j)]);
      assert(std::abs(getRef(resultEx,i,j) - tmp) <= REL_TOL*std::abs(tmp));
      assert(std::abs(getRef(resultCon,i,j)- tmp) <= REL_TOL*std::abs(tmp));
      assert(std::abs(getRef(resultTag,i,j)- tmp) <= REL_TOL*std::abs(tmp));
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
      assert(std::abs(getRef(resultEx,i,j) - getRelIndex(shape,i,j)) <= REL_TOL*std::abs(data[getRelIndex(shape,i,j)]));
      assert(std::abs(getRef(resultCon,i,j) - getRelIndex(shape,i,j)) <= REL_TOL*std::abs(data[getRelIndex(shape,i,j)]));
      assert(std::abs(getRef(resultTag,i,j) - getRelIndex(shape,i,j)) <= REL_TOL*std::abs(data[getRelIndex(shape,i,j)]));
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
//   assert(std::abs(resultEx.getPointDataView()() - 0) <= REL_TOL*0);
//   assert(std::abs(resultCon.getPointDataView()() - 0) <= REL_TOL*0);
//   assert(std::abs(resultTag.getPointDataView()() - 0) <= REL_TOL*0);
  assert(std::abs(resultEx.getDataAtOffset(0) - 0) <= REL_TOL*0);
  assert(std::abs(resultCon.getDataAtOffset(0) - 0) <= REL_TOL*0);
  assert(std::abs(resultTag.getDataAtOffset(0) - 0) <= REL_TOL*0);

  cout << "\tTest Data::maxval." << endl;
  resultEx.copy(baseEx.maxval());
  resultCon.copy(baseCon.maxval());
  resultTag.copy(baseTag.maxval());
  assert(std::abs(resultEx.getDataAtOffset(0) - 5) <= REL_TOL*5);
  assert(std::abs(resultCon.getDataAtOffset(0) - 5) <= REL_TOL*5);
  assert(std::abs(resultTag.getDataAtOffset(0) - 5) <= REL_TOL*5);

}


void DataTestCase::testMemAlloc() {

  //
  // Simple little sanity check for the memory allocator

  cout << endl;

  Data *testData;
  for (int i=0; i<1000; i++) {
    testData = new Data(0.0, DataTypes::ShapeType(), FunctionSpace(), true);
    delete testData;
  }

  DataTypes::ShapeType viewShape;
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
  testSuite->addTest (new TestCaller< DataTestCase>("testCopying",&DataTestCase::testCopying));
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
