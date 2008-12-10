
/*******************************************************
*
* Copyright (c) 2003-2008 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


#include <iostream>
#if defined (_WIN32) && defined (__INTEL_COMPILER)
#include <mathimf.h>
#else
#include <math.h>
#endif

#include "DataTestCase.h"

#include "escript/FunctionSpace.h"
#include "esysUtils/EsysException.h"

#include "escript/Data.h"
#include "escript/DataLazy.h"
#include "escript/EscriptParams.h"

#define AUTOLAZYON setEscriptParamInt("AUTOLAZY",1);
#define AUTOLAZYOFF setEscriptParamInt("AUTOLAZY",0);
#define CHECKAUTOLAZY (getEscriptParamInt("AUTOLAZY")==1)
#define SAVELAZYSTATE int LAZYSTATE=getEscriptParamInt("AUTOLAZY");
#define RESTORELAZYSTATE setEscriptParamInt("AUTOLAZY",LAZYSTATE);



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


void DataTestCase::testCopyingWorker(bool delayed)
{

  using namespace escript::DataTypes;
  cout << endl;

  DataTypes::ShapeType shape;
  shape.push_back(2);
  shape.push_back(3);
  DataTypes::ValueType data(DataTypes::noValues(shape),1);
  const int NUMDATS=3;
  Data* dats[NUMDATS];
  const char* strs[]={"DataConstant", "DataTagged", "DataExpanded"};
  dats[0]=new Data(new DataConstant(FunctionSpace(),shape,data));
  dats[1]=new Data(new DataTagged(FunctionSpace(),shape,data));
  dats[2]=new Data(new DataExpanded(FunctionSpace(),shape,data));
  if (delayed)
  {
    dats[0]->delaySelf();
    dats[1]->delaySelf();
    dats[2]->delaySelf();
  }

  for (int k=0;k<NUMDATS;++k)
  {
	cout << "\tTest deep copy " << strs[k] << endl;
	Data* d=dats[k];
	Data deep=d->copySelf();	// test self copy
	if (delayed)
	{
	  assert(deep.isLazy());
	}
	for (int i=0;i<DataTypes::noValues(shape);++i)
	{
	if (d->getDataAtOffset(i)!=deep.getDataAtOffset(i))
		assert(false);
	}
	if (delayed)
	{
	   d->delaySelf();
	}
	d->setToZero();
	if (delayed)
	{
	  assert(d->isLazy());
	}
	for (int i=0;i<DataTypes::noValues(shape);++i)
	{
	if (d->getDataAtOffset(i)==deep.getDataAtOffset(i))
		assert(false);
	}
        if (delayed)
	{
	   d->delaySelf();
	   deep.delaySelf();
	}
	d->copy(deep);			// test copy from object
	if (delayed)
	{
	  assert(d->isLazy());
	}
	for (int i=0;i<DataTypes::noValues(shape);++i)
	{
	if (d->getDataAtOffset(i)!=deep.getDataAtOffset(i))
		assert(false);
	}
	d->setToZero();
	for (int i=0;i<DataTypes::noValues(shape);++i)
	{
	if (d->getDataAtOffset(i)==deep.getDataAtOffset(i))
		assert(false);
	}
	delete dats[k];
  }






}


void DataTestCase::testSlicingWorker(bool delayed)
{

  using namespace escript::DataTypes;
  cout << endl;
  {
   DataTypes::ShapeType viewShape;
   viewShape.push_back(2);
   viewShape.push_back(3);

   const int NUMDATS=3;
   const char* strs[]={"DataConstant", "DataTagged","DataExpanded"};
   bool tags[]={false,true,false};	// is the slice of this data supposed to be tagged
   Data* dats[NUMDATS];
   for (int k=0;k<NUMDATS;++k)
   {
    	dats[k]=new Data(1.3, viewShape);
   }
   dats[1]->tag();
   dats[2]->expand();
   for (int k=0;k<NUMDATS;++k)
   {
	Data* temp=dats[k];
	dats[k]=new Data(dats[k]->delay());
	delete temp;
   }
   for (int k=0;k<NUMDATS;++k)
   {
	cout << "\t\tTest get-slicing " << strs[k] << endl;
    	dats[k]->getDataAtOffset(dats[k]->getDataOffset(0,0)+getRelIndex(viewShape,0,0))=1.0;
    	dats[k]->getDataAtOffset(dats[k]->getDataOffset(0,0)+getRelIndex(viewShape,1,1))=2.0;

    	DataTypes::RegionType region;
    	region.push_back(DataTypes::RegionType::value_type(0,0));
    	region.push_back(DataTypes::RegionType::value_type(0,0));

    	Data slice1(dats[k]->getSlice(region));

    	if (tags[k]) {assert(slice1.isTagged());}
    	assert(slice1.getDataPointRank()==0);
    	assert(slice1.getDataPoint(0,0)==1.0);

	//
	// create a rank 2 slice with one value
	
	region.clear();
	region.push_back(DataTypes::RegionType::value_type(0,1));
	region.push_back(DataTypes::RegionType::value_type(0,1));
	
	Data slice2(dats[k]->getSlice(region));
	
	//cout << slice2.toString() << endl;
	
	if (tags[k]) {assert(slice2.isTagged());}
	assert(slice2.getDataPointRank()==2);
	
	assert(slice2.getDataAtOffset(slice2.getDataOffset(0,0)+getRelIndex(slice2.getDataPointShape(),0,0))==1.0);

	//
	// create a rank 2 slice with four values
	
	region.clear();
	region.push_back(DataTypes::RegionType::value_type(0,2));
	region.push_back(DataTypes::RegionType::value_type(0,2));
	
	Data slice3(dats[k]->getSlice(region));
	
	//cout << slice3.toString() << endl;
	
	if (tags[k]) {assert(slice3.isTagged());}
	assert(slice3.getDataPointRank()==2);
	assert(getRef(slice3,0,0,0,0)==1.0);
	assert(getRef(slice3,0,0,0,1)==1.3);
	assert(getRef(slice3,0,0,1,0)==1.3);
	assert(getRef(slice3,0,0,1,1)==2.0);
   }

   // now some extra tests for tagged data (dats[1])

   //
   // add a value for tag "1"

   DataTypes::ValueType viewData(6);
   for (int i=0;i<viewData.size();i++) {
    viewData[i]=i;
   }
   dats[1]->setTaggedValueFromCPP(1, viewShape, viewData);

    //
    // create a full slice

   DataTypes::RegionType region;
   region.push_back(DataTypes::RegionType::value_type(0,2));
   region.push_back(DataTypes::RegionType::value_type(0,3));

   Data slice4(dats[1]->getSlice(region));

   assert(slice4.isTagged());
   assert(slice4.getDataPointRank()==2);
   assert(getRef(slice4,0,0,0,0)==0);
   assert(getRef(slice4,0,0,0,1)==2);
   assert(getRef(slice4,0,0,0,2)==4);
   assert(getRef(slice4,0,0,1,0)==1);
   assert(getRef(slice4,0,0,1,1)==3);
   assert(getRef(slice4,0,0,1,2)==5);

   for (int k=0;k<NUMDATS;++k)
   {
	delete dats[k];
   }
 }

 {
  DataTypes::ShapeType viewShape;
  viewShape.push_back(2);
  viewShape.push_back(3);

  const int NUMDATS=3;
  const char* strs[]={"DataConstant", "DataTagged","DataExpanded"};
  Data* dats[NUMDATS];
  Data* src[NUMDATS];
  for (int k=0;k<NUMDATS;++k)
  {
 	dats[k]=new Data(1.3, viewShape);
    	src[k]=new Data(10,DataTypes::scalarShape);
  }
  dats[1]->tag();
  src[1]->tag();
  dats[2]->expand();
  src[2]->expand();
  if (delayed)
  {
    for(int k=0;k<NUMDATS;++k)
    {
	if (delayed)
	{
	  Data* temp=dats[k];
	  dats[k]=new Data(dats[k]->delay());	// coz delay returns an object not a pointer
	  delete temp;
	  temp=src[k];
	  src[k]=new Data(src[k]->delay());
	  delete temp;
	}
    }
  }
  for (int k=0;k<NUMDATS;++k)
  {
	cout << "\t\tTest set-slicing " << strs[k] << endl;
	Data target(1.3,viewShape);
	if (k==2) {target.expand();}
	DataTypes::RegionType region;
	region.push_back(DataTypes::RegionType::value_type(1,1));
	region.push_back(DataTypes::RegionType::value_type(1,1));
	target.setSlice(*(src[k]),region);
	assert(getRef(target,0,0,1,1)==src[k]->getDataPoint(0,0));
  }
  
  // some extra tests on tagged data

  //
  // add a value for tag "1" to target

  DataTypes::ValueType viewData(6);
  for (int i=0;i<viewData.size();i++) {
	viewData[i]=i;
  }

  Data target(1.3,viewShape,FunctionSpace(),false);
  target.tag();
  target.setTaggedValueFromCPP(1, viewShape, viewData);

    //cout << "target:\n" << target.toString() << endl;

    //
    // set a slice in target from source

  DataTypes::RegionType region;
  region.push_back(DataTypes::RegionType::value_type(0,0));
  region.push_back(DataTypes::RegionType::value_type(1,1));

  target.setSlice(*src[1],region);

  assert(target.isTagged());
  assert(target.getDataPointRank()==2);
  assert(getRef(target,0,0,0,0)==0);
  assert(getRef(target,0,0,0,1)==src[1]->getDataPoint(0,0));
  assert(getRef(target,0,0,0,2)==4);
  assert(getRef(target,0,0,1,0)==1);
  assert(getRef(target,0,0,1,1)==3);
  assert(getRef(target,0,0,1,2)==5);

  //
  // add a value for tag "2" to source

  DataTypes::ShapeType viewShape2;
  DataTypes::ValueType viewData2(1);
  viewData2[0]=6;
  src[1]->setTaggedValueFromCPP(2, viewShape2, viewData2);

  region.clear();
  region.push_back(DataTypes::RegionType::value_type(0,0));
  region.push_back(DataTypes::RegionType::value_type(1,1));

  target.setSlice(*src[1],region);

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
  assert(targetData[3]==1.3);
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
  assert(targetData[15]==1.3);
  assert(targetData[16]==1.3);
  assert(targetData[17]==1.3);


  for (int k=0;k<NUMDATS;++k)
  {
	delete dats[k];
	delete src[k];
  }

 }

}

// This is to test new copy routines, existing tests should remain where they are
void DataTestCase::testCopying()
{
  cout << "\n\tReadyData." << endl;
  testCopyingWorker(false);
  cout << "\n\tLazyData." << endl;
  testCopyingWorker(true);
}

void DataTestCase::testSlicing() {
  cout << "\n\tReadyData." << endl;
  testSlicingWorker(false);
  cout << "\n\tLazyData." << endl;
  testSlicingWorker(true);
}


void DataTestCase::testSomeDriver(bool autolazy)
{
  cout << endl;
  SAVELAZYSTATE
  if (autolazy)
  {
	AUTOLAZYON
	cout << "\tNow testing using autolazy." << endl;
  }
  cout << "\tCreate a Data object." << endl;

  DataTypes::ShapeType viewShape;
  viewShape.push_back(3);
  DataTypes::ValueType viewData(3);
  for (int i=0;i<viewShape[0];++i) {
    viewData[i]=i;
  }

  bool expanded=true;
  Data exData(viewData,viewShape,FunctionSpace(),expanded);
  Data cData(viewData,viewShape);
  Data result;

  assert(exData.isExpanded());
  assert(cData.isConstant());
  assert(result.isEmpty());

  cout << "\tTest some basic operations" << endl;
  result=exData*cData;
  cout << CHECKAUTOLAZY << " " << result.isLazy() << " " << result.isExpanded()<< endl;
  assert(CHECKAUTOLAZY?result.isLazy():result.isExpanded());

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
  RESTORELAZYSTATE

}

void DataTestCase::testSome() {
  testSomeDriver(false);
  testSomeDriver(true);
}



// This method tests to see if resolve() produces results of the correct type
void DataTestCase::testResolveType()
{
  cout << endl;
  cout << "\tTesting resolve()\n";
  DataTypes::ShapeType viewShape;
  viewShape.push_back(2);
  viewShape.push_back(3);
  viewShape.push_back(4);
  DataTypes::ValueType viewData(2*3*4);
  for (int i=0;i<DataTypes::noValues(viewShape);++i) {
    viewData[i]=i;
  }
  Data c1(viewData,viewShape);
  Data t1(viewData,viewShape);
  Data e1(viewData,viewShape);
  t1.tag();
  e1.expand();
  c1.delaySelf();
  t1.delaySelf();
  e1.delaySelf();
  Data d1=c1+c1;
  assert(d1.isLazy());
  assert((d1.resolve(),d1.isConstant()));
  d1=c1+t1;
  assert(d1.isLazy());
  assert((d1.resolve(),d1.isTagged()));
  d1=t1+c1;
  assert(d1.isLazy());
  assert((d1.resolve(),d1.isTagged()));
  d1=t1+t1;
  assert(d1.isLazy());
  assert((d1.resolve(),d1.isTagged()));
  d1=c1+e1;
  assert(d1.isLazy());
  assert((d1.resolve(),d1.isExpanded()));
  d1=e1+c1;
  assert(d1.isLazy());
  assert((d1.resolve(),d1.isExpanded()));
  d1=e1+t1;
  assert(d1.isLazy());
  assert((d1.resolve(),d1.isExpanded()));
  d1=t1+e1;
  assert(d1.isLazy());
  assert((d1.resolve(),d1.isExpanded()));
  d1=e1+e1;
  assert(d1.isLazy());
  assert((d1.resolve(),d1.isExpanded()));
  cout << "\tTesting tag()\n";
  c1.tag();
  assert(c1.isTagged());
  t1.tag();
  assert(t1.isTagged());
  try
  {
	e1.tag();
	assert(false);		// this should have thrown
  } catch(...) {}
  cout << "\tTesting expand()\n";
  Data c2(viewData,viewShape);
  Data t2(viewData,viewShape);
  Data e2(viewData,viewShape);
  t2.tag();
  e2.expand();
  c2.delaySelf();
  t2.delaySelf();
  e2.delaySelf();
  c2.expand();
  assert(c2.isExpanded());
  t2.expand();
  assert(t2.isExpanded());
  e2.expand();
  assert(e2.isExpanded());
}

void DataTestCase::testDataConstant() {

  cout << endl;

  cout << "\tCreate a DataConstant object." << endl;

  DataTypes::ShapeType viewShape;
  viewShape.push_back(2);
  viewShape.push_back(3);
  viewShape.push_back(4);
  DataTypes::ValueType viewData(2*3*4);
  for (int i=0;i<DataTypes::noValues(viewShape);++i) {
    viewData[i]=i;
  }

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

    DataTypes::ShapeType viewShape;
    viewShape.push_back(3);

    DataTypes::ValueType viewData(3);
    for (int i=0;i<viewShape[0];i++) {
      viewData[i]=i;
    }
    int arr[1]={1};		// iso c++ does not like empty arrays
    DataTagged* dt=new DataTagged(FunctionSpace(),viewShape,arr,viewData); 
    Data myData(dt);

    assert(!myData.isEmpty());
    assert(myData.isTagged());
    assert(myData.getTagNumber(0)==1);
    assert(myData.getDataPointRank()==1);
    assert(myData.getLength()==3);
    
    assert(myData.getNoValues()==3);
    assert(myData.getDataAtOffset(0)==0.0);
    assert(myData.getDataAtOffset(1)==1.0);
    assert(myData.getDataAtOffset(2)==2.0);

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
    assert(offset==3);
    assert(myData.getDataPointRank()==1);
    assert(myData.getNoValues()==3);

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

    assert(!myData.isEmpty());
    assert(myData.isTagged());
    assert(myData.getTagNumber(0)==1);
    assert(myData.getDataPointRank()==2);
    assert(myData.getLength()==6);

    // check default value
    assert(!myData.isEmpty());
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


void DataTestCase::testOperations()
{

  cout << endl;

  // define the shape for the test data
  DataTypes::ShapeType shape;
  shape.push_back(2);
  shape.push_back(3);

  // allocate the data 
  DataTypes::ValueType data(DataTypes::noValues(shape),0);

  // assign values to the data
  for (int i=0;i<shape[0];i++) {
    for (int j=0;j<shape[1];j++) {
      data[getRelIndex(shape,i,j)]=getRelIndex(shape,i,j);
    }
  }



  Data dats[]={Data(data,shape,FunctionSpace(),false),
		Data(data,shape,FunctionSpace(),false),
		Data(data,shape,FunctionSpace(),true),
		Data(data,shape,FunctionSpace(),false),
		Data(data,shape,FunctionSpace(),false),
		Data(data,shape,FunctionSpace(),true)};
  const int NUMDATS=6;
  const int LAZY=3;		// where do the lazy objects start?

//   Data baseEx(data,shape,FunctionSpace(),true);
//   Data baseCon(data,shape,FunctionSpace(),false);
//   Data baseTag(data,shape,FunctionSpace(),false);
  Data& baseCon=dats[0];
  Data& baseTag=dats[1];
  Data& baseEx=dats[2];
  baseTag.tag();
  dats[4].tag();
  dats[3].delaySelf();
  dats[4].delaySelf();
  dats[5].delaySelf();

  assert(baseEx.isExpanded());
  assert(baseCon.isConstant());
  assert(baseTag.isTagged());

  Data results[NUMDATS];
//   Data& resultEx=results[0];
//   Data& resultCon=results[1];
//   Data& resultTag=results[2];

  // create 0 <= smalldata <= 1 for testing trig functions

  DataTypes::ValueType smalldata(DataTypes::noValues(shape),0);

  // assign values to the data
  for (int i=0;i<shape[0];i++) {
    for (int j=0;j<shape[1];j++) {
      smalldata[getRelIndex(shape,i,j)]=(i==0 && j==0)?0:1.0/(getRelIndex(shape,i,j)+1);
    }
  }
  Data sdats[]={Data(smalldata,shape,FunctionSpace(),false),
		Data(smalldata,shape,FunctionSpace(),false),
		Data(smalldata,shape,FunctionSpace(),true),
		Data(smalldata,shape,FunctionSpace(),false),
		Data(smalldata,shape,FunctionSpace(),false),
		Data(smalldata,shape,FunctionSpace(),true)};
  sdats[1].tag();
  sdats[4].tag();
  sdats[3].delaySelf();
  sdats[4].delaySelf();
  sdats[5].delaySelf();



  // test unary operations

  double tmp;
  cout << "\tTest Data::pow." << endl;
  Data power(3.0,shape,FunctionSpace(),true);
  for (int z=0;z<NUMDATS;++z)
  {
    results[z].copy(dats[z].powD(power));
    if (z>=LAZY)
    {
	assert(results[z].isLazy());
    }
  }
  for (int i=0;i<shape[0];i++) {
    for (int j=0;j<shape[1];j++) {
      tmp=pow((double)data[getRelIndex(shape,i,j)],(double)3.0);
      for (int z=0;z<NUMDATS;++z)
      {
	assert(std::abs(getRef(results[z],i,j) - tmp) <= REL_TOL*std::abs(tmp));
      }
    }
  }

  cout << "\tTest Data::sin." << endl;
  for (int z=0;z<NUMDATS;++z)
  {
    results[z].copy(dats[z].sin());
    if (z>=LAZY)
    {
	assert(results[z].isLazy());
    }
  }
  for (int i=0;i<shape[0];i++) {
    for (int j=0;j<shape[1];j++) {
      tmp=sin((double)data[getRelIndex(shape,i,j)]);
      for (int z=0;z<NUMDATS;++z)
      {
	assert(std::abs(getRef(results[z],i,j) - tmp) <= REL_TOL*std::abs(tmp));
      }
    }
  }

  cout << "\tTest Data::cos." << endl;
  for (int z=0;z<NUMDATS;++z)
  {
    results[z].copy(dats[z].cos());
    if (z>=LAZY)
    {
	assert(results[z].isLazy());
    }
  }
  for (int i=0;i<shape[0];i++) {
    for (int j=0;j<shape[1];j++) {
      tmp=cos((double)data[getRelIndex(shape,i,j)]);
      for (int z=0;z<NUMDATS;++z)
      {
	assert(std::abs(getRef(results[z],i,j) - tmp) <= REL_TOL*std::abs(tmp));
      }
    }
  }

  cout << "\tTest Data::tan." << endl;
  for (int z=0;z<NUMDATS;++z)
  {
    results[z].copy(dats[z].tan());
    if (z>=LAZY)
    {
	assert(results[z].isLazy());
    }
  }
  for (int i=0;i<shape[0];i++) {
    for (int j=0;j<shape[1];j++) {
      tmp=tan((double)data[getRelIndex(shape,i,j)]);
      for (int z=0;z<NUMDATS;++z)
      {
	assert(std::abs(getRef(results[z],i,j) - tmp) <= REL_TOL*std::abs(tmp));
      }
    }
  }

  cout << "\tTest Data::asin." << endl;
  for (int z=0;z<NUMDATS;++z)
  {
    results[z].copy(sdats[z].asin());
    if (z>=LAZY)
    {
	assert(results[z].isLazy());
    }
  }
  for (int i=0;i<shape[0];i++) {
    for (int j=0;j<shape[1];j++) {
      tmp=asin((double)smalldata[getRelIndex(shape,i,j)]);
      for (int z=0;z<NUMDATS;++z)
      {
	assert(std::abs(getRef(results[z],i,j) - tmp) <= REL_TOL*std::abs(tmp));
      }
    }
  }

  cout << "\tTest Data::acos." << endl;
  for (int z=0;z<NUMDATS;++z)
  {
    results[z].copy(sdats[z].acos());
    if (z>=LAZY)
    {
	assert(results[z].isLazy());
    }
  }
  for (int i=0;i<shape[0];i++) {
    for (int j=0;j<shape[1];j++) {
      tmp=acos((double)smalldata[getRelIndex(shape,i,j)]);
      for (int z=0;z<NUMDATS;++z)
      {
	assert(std::abs(getRef(results[z],i,j) - tmp) <= REL_TOL*std::abs(tmp));
      }
    }
  }

  cout << "\tTest Data::atan." << endl;
  for (int z=0;z<NUMDATS;++z)
  {
    results[z].copy(sdats[z].atan());
    if (z>=LAZY)
    {
	assert(results[z].isLazy());
    }
  }
  for (int i=0;i<shape[0];i++) {
    for (int j=0;j<shape[1];j++) {
      tmp=atan((double)smalldata[getRelIndex(shape,i,j)]);
      for (int z=0;z<NUMDATS;++z)
      {
	assert(std::abs(getRef(results[z],i,j) - tmp) <= REL_TOL*std::abs(tmp));
      }
    }
  }

  cout << "\tTest Data::sinh." << endl;
  for (int z=0;z<NUMDATS;++z)
  {
    results[z].copy(dats[z].sinh());
    if (z>=LAZY)
    {
	assert(results[z].isLazy());
    }
  }
  for (int i=0;i<shape[0];i++) {
    for (int j=0;j<shape[1];j++) {
      tmp=sinh((double)data[getRelIndex(shape,i,j)]);
      for (int z=0;z<NUMDATS;++z)
      {
	assert(std::abs(getRef(results[z],i,j) - tmp) <= REL_TOL*std::abs(tmp));
      }
    }
  }

  cout << "\tTest Data::cosh." << endl;
  for (int z=0;z<NUMDATS;++z)
  {
    results[z].copy(dats[z].cosh());
    if (z>=LAZY)
    {
	assert(results[z].isLazy());
    }
  }
  for (int i=0;i<shape[0];i++) {
    for (int j=0;j<shape[1];j++) {
      tmp=cosh((double)data[getRelIndex(shape,i,j)]);
      for (int z=0;z<NUMDATS;++z)
      {
	assert(std::abs(getRef(results[z],i,j) - tmp) <= REL_TOL*std::abs(tmp));
      }
    }
  }

  cout << "\tTest Data::tanh." << endl;
  for (int z=0;z<NUMDATS;++z)
  {
    results[z].copy(dats[z].tanh());
    if (z>=LAZY)
    {
	assert(results[z].isLazy());
    }
  }
  for (int i=0;i<shape[0];i++) {
    for (int j=0;j<shape[1];j++) {
      tmp=tanh((double)data[getRelIndex(shape,i,j)]);
      for (int z=0;z<NUMDATS;++z)
      {
	assert(std::abs(getRef(results[z],i,j) - tmp) <= REL_TOL*std::abs(tmp));
      }
    }
  }

  // rather than accomodate the different windows operations directly I'll just use inverse functions
  cout << "\tTest Data::asinh." << endl;
  for (int z=0;z<NUMDATS;++z)
  {
    results[z].copy(dats[z].asinh().sinh());
    if (z>=LAZY)
    {
	assert(results[z].isLazy());
    }
  }
  for (int i=0;i<shape[0];i++) {
    for (int j=0;j<shape[1];j++) {
      tmp=data[getRelIndex(shape,i,j)];
      for (int z=0;z<NUMDATS;++z)
      {
	assert(std::abs(getRef(results[z],i,j) - tmp) <= REL_TOL*std::abs(tmp));
      }
    }
  }

  cout << "\tTest Data::acosh." << endl;
  for (int z=0;z<NUMDATS;++z)
  {
    results[z].copy(dats[z].acosh().cosh());
    if (z>=LAZY)
    {
	assert(results[z].isLazy());
    }
  }
  for (int i=0;i<shape[0];i++) {
    for (int j=0;j<shape[1];j++) {
      if (i==0 && j==0) break;
      tmp=data[getRelIndex(shape,i,j)];
      for (int z=0;z<NUMDATS;++z)
      {
	assert(std::abs(getRef(results[z],i,j) - tmp) <= REL_TOL*std::abs(tmp));
      }
    }
  }

  cout << "\tTest Data::atanh." << endl;
  for (int z=0;z<NUMDATS;++z)
  {
    results[z].copy(dats[z].tanh().atanh());		// if these are the other way around the results are
    if (z>=LAZY)					// undefined
    {
	assert(results[z].isLazy());
    }
  }
  for (int i=0;i<shape[0];i++) {
    for (int j=0;j<shape[1];j++) {
      tmp=data[getRelIndex(shape,i,j)];
      for (int z=0;z<NUMDATS;++z)
      {
	assert(std::abs(getRef(results[z],i,j) - tmp) <= REL_TOL*std::abs(tmp));
      }
    }
  }

  cout << "\tTest Data::log." << endl;
  for (int z=0;z<NUMDATS;++z)
  {
    results[z].copy(dats[z].log());
    if (z>=LAZY)
    {
	assert(results[z].isLazy());
    }
  }
  for (int i=0;i<shape[0];i++) {
    for (int j=0;j<shape[1];j++) {
      if (i==0 && j==0) break; 
      tmp=log((double)data[getRelIndex(shape,i,j)]);
      for (int z=0;z<NUMDATS;++z)
      {
	assert(std::abs(getRef(results[z],i,j) - tmp) <= REL_TOL*std::abs(tmp));
      }
    }
  }

  cout << "\tTest Data::log10." << endl;
  for (int z=0;z<NUMDATS;++z)
  {
    results[z].copy(dats[z].log10());
    if (z>=LAZY)
    {
	assert(results[z].isLazy());
    }
  }
  for (int i=0;i<shape[0];i++) {
    for (int j=0;j<shape[1];j++) {
      if (i==0 && j==0) break; 
      tmp=log10((double)data[getRelIndex(shape,i,j)]);
      for (int z=0;z<NUMDATS;++z)
      {
	assert(std::abs(getRef(results[z],i,j) - tmp) <= REL_TOL*std::abs(tmp));
      }
    }
  }
#if defined (_WIN32) && !defined (__INTEL_COMPILER)
  cout << "\tSkip test Data::erf on windows with MSVC compiler." << endl;
#else
  cout << "\tTest Data::erf." << endl;
  for (int z=0;z<NUMDATS;++z)
  {
    results[z].copy(dats[z].erf());
    if (z>=LAZY)
    {
	assert(results[z].isLazy());
    }
  }
  for (int i=0;i<shape[0];i++) {
    for (int j=0;j<shape[1];j++) {
      if (i==0 && j==0) break; 
      tmp=erf((double)data[getRelIndex(shape,i,j)]);
      for (int z=0;z<NUMDATS;++z)
      {
	assert(std::abs(getRef(results[z],i,j) - tmp) <= REL_TOL*std::abs(tmp));
      }
    }
  }
#endif


  cout << "\tTest Data::abs." << endl;
  for (int z=0;z<NUMDATS;++z)
  {
    results[z].copy(dats[z].abs());
    if (z>=LAZY)
    {
	assert(results[z].isLazy());
    }
  }
  for (int i=0;i<shape[0];i++) {
    for (int j=0;j<shape[1];j++) {
      tmp=abs((double)data[getRelIndex(shape,i,j)]);
      for (int z=0;z<NUMDATS;++z)
      {
	assert(std::abs(getRef(results[z],i,j) - tmp) <= REL_TOL*std::abs(tmp));
      }
    }
  }

  cout << "\tTest Data::sign (positive)." << endl;
  for (int z=0;z<NUMDATS;++z)
  {
    results[z].copy(dats[z].sign());
    if (z>=LAZY)
    {
	assert(results[z].isLazy());
    }
  }
  for (int i=0;i<shape[0];i++) {
    for (int j=0;j<shape[1];j++) {
      tmp=(i==0 && j==0)?0:1;
      for (int z=0;z<NUMDATS;++z)
      {
	assert(std::abs(getRef(results[z],i,j) - tmp) <= REL_TOL*std::abs(tmp));
      }
    }
  } 

  cout << "\tTest Data::sign (negative)." << endl;
  for (int z=0;z<NUMDATS;++z)
  {
    results[z].copy(dats[z].neg().sign());
    if (z>=LAZY)
    {
	assert(results[z].isLazy());
    }
  }
  for (int i=0;i<shape[0];i++) {
    for (int j=0;j<shape[1];j++) {
      tmp=(i==0 && j==0)?0:-1;
      for (int z=0;z<NUMDATS;++z)
      {
	assert(std::abs(getRef(results[z],i,j) - tmp) <= REL_TOL*std::abs(tmp));
      }
    }
  } 


  cout << "\tTest Data::exp." << endl;
  for (int z=0;z<NUMDATS;++z)
  {
    results[z].copy(dats[z].exp());
    if (z>=LAZY)
    {
	assert(results[z].isLazy());
    }
  }
  for (int i=0;i<shape[0];i++) {
    for (int j=0;j<shape[1];j++) {
      tmp=exp((double)data[getRelIndex(shape,i,j)]);
      for (int z=0;z<NUMDATS;++z)
      {
	assert(std::abs(getRef(results[z],i,j) - tmp) <= REL_TOL*std::abs(tmp));
      }
    }
  }

  cout << "\tTest Data::sqrt." << endl;
  for (int z=0;z<NUMDATS;++z)
  {
    results[z].copy(dats[z].sqrt());
    if (z>=LAZY)
    {
	assert(results[z].isLazy());
    }
  }
  for (int i=0;i<shape[0];i++) {
    for (int j=0;j<shape[1];j++) {
      tmp=sqrt((double)data[getRelIndex(shape,i,j)]);
      for (int z=0;z<NUMDATS;++z)
      {
	assert(std::abs(getRef(results[z],i,j) - tmp) <= REL_TOL*std::abs(tmp));
      }
    }
  }

  cout << "\tTest Data::neg." << endl;
  for (int z=0;z<NUMDATS;++z)
  {
    results[z].copy(dats[z].neg());
    if (z>=LAZY)
    {
	assert(results[z].isLazy());
    }
  }
  for (int i=0;i<shape[0];i++) {
    for (int j=0;j<shape[1];j++) {
      tmp=-data[getRelIndex(shape,i,j)];
      for (int z=0;z<NUMDATS;++z)
      {
	assert(std::abs(getRef(results[z],i,j) - tmp) <= REL_TOL*std::abs(tmp));
      }
    }
  }

  cout << "\tTest Data::pos." << endl;
  for (int z=0;z<NUMDATS;++z)
  {
    results[z].copy(dats[z].pos());
    if (z>=LAZY)
    {
	assert(results[z].isLazy());
    }
  }
  for (int i=0;i<shape[0];i++) {
    for (int j=0;j<shape[1];j++) {
      for (int z=0;z<NUMDATS;++z)
      {
	assert(std::abs(getRef(results[z],i,j) - getRelIndex(shape,i,j)) <= REL_TOL*std::abs(data[getRelIndex(shape,i,j)]));
      }
    }
  }

  // test reduction operations

  cout << "\tTest Data::Lsup." << endl;
  for (int z=0;z<NUMDATS;++z)
  {
    assert(std::abs(dats[z].Lsup() - 5) <= REL_TOL*5);
  }

  cout << "\tTest Data::sup." << endl;
  for (int z=0;z<NUMDATS;++z)
  {
    assert(std::abs(dats[z].sup() - 5) <= REL_TOL*5);
  }

  cout << "\tTest Data::inf." << endl;
  for (int z=0;z<NUMDATS;++z)
  {
    assert(std::abs(dats[z].inf() - 0) <= REL_TOL*0);
  }

  // test data-point reduction operations

  cout << "\tTest Data::minval." << endl;
  for (int z=0;z<NUMDATS;++z)
  {
    results[z].copy(dats[z].minval());
  }
  for (int z=0;z<NUMDATS;++z)
  {
    assert(std::abs(results[z].getDataAtOffset(0) - 0) <= REL_TOL*0);
  }
  

  cout << "\tTest Data::maxval." << endl;
  for (int z=0;z<NUMDATS;++z)
  {
    results[z].copy(dats[z].maxval());
  }
  for (int z=0;z<NUMDATS;++z)
  {
    assert(std::abs(results[z].getDataAtOffset(0) - 5) <= REL_TOL*5);
  }

}


// Here we test the binary operators in complex expressions
void DataTestCase::testBinary()
{

  cout << endl;

  // define the shape for the test data
  DataTypes::ShapeType shape;
  shape.push_back(2);
  shape.push_back(3);

  // allocate the data 
  DataTypes::ValueType data(DataTypes::noValues(shape),0);

  // assign values to the data
  for (int i=0;i<shape[0];i++) {
    for (int j=0;j<shape[1];j++) {
      data[getRelIndex(shape,i,j)]=getRelIndex(shape,i,j)+2;	// so we get no zeros
    }
  }


  Data one(1.0,DataTypes::scalarShape,FunctionSpace());
  Data two(2.0,DataTypes::scalarShape,FunctionSpace());
  Data dats[]={Data(data,shape,FunctionSpace(),false),
		Data(data,shape,FunctionSpace(),false),
		Data(data,shape,FunctionSpace(),true),
		Data(data,shape,FunctionSpace(),false),
		Data(data,shape,FunctionSpace(),false),
		Data(data,shape,FunctionSpace(),true)};
  dats[1].tag();
  dats[4].tag();
  const int NUMDATS=6;
  const int LAZY=3;
  dats[3].delaySelf();
  dats[4].delaySelf();
  dats[5].delaySelf();
  for (int z=0;z<NUMDATS;++z)
  {
	Data& a=dats[z];
	Data r1=(((a+a)/two)+a-a)*one;	// scalar/*, matrix+-
	Data r2=(((a*a)/a+one)-one);	// scalar+-, matrix*/
	Data r3=(a.powD(two)/a.powD(one)); // scalar power
	Data r4=a.powD(a);		// matrix power
	if (z>LAZY)
	{
	  assert(r1.isLazy() && r2.isLazy() && r3.isLazy() && r4.isLazy());
	}
	for (int i=0;i<DataTypes::noValues(shape);++i)
	{
	  assert(std::abs(r1.getDataAtOffset(i)-data[i]) <= REL_TOL*data[i]);
	  assert(std::abs(r2.getDataAtOffset(i)-data[i]) <= REL_TOL*data[i]);
	  assert(std::abs(r3.getDataAtOffset(i)-data[i]) <= REL_TOL*data[i]);
	  assert(std::abs(r4.getDataAtOffset(i)-pow(data[i],i)) <=REL_TOL*pow(data[i],i));
	}
  }
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
  testSuite->addTest (new TestCaller< DataTestCase>("testSome",&DataTestCase::testSome));
  testSuite->addTest (new TestCaller< DataTestCase>("testDataConstant",&DataTestCase::testDataConstant));
  testSuite->addTest (new TestCaller< DataTestCase>("testDataTagged",&DataTestCase::testDataTagged));
  testSuite->addTest (new TestCaller< DataTestCase>("testDataTaggedExceptions",&DataTestCase::testDataTaggedExceptions));
  testSuite->addTest (new TestCaller< DataTestCase>("testConstructors",&DataTestCase::testConstructors));
  testSuite->addTest (new TestCaller< DataTestCase>("testSlicing",&DataTestCase::testSlicing));
  testSuite->addTest (new TestCaller< DataTestCase>("testOperations",&DataTestCase::testOperations));
  //testSuite->addTest (new TestCaller< DataTestCase>("testRefValue",&DataTestCase::testRefValue));
  testSuite->addTest (new TestCaller< DataTestCase>("testMemAlloc",&DataTestCase::testMemAlloc));
  testSuite->addTest (new TestCaller< DataTestCase>("Resolving",&DataTestCase::testResolveType));

  return testSuite;
}
