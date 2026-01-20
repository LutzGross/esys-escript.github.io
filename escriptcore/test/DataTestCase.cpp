
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* https://github.com/LutzGross/esys-escript.github.io
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/

#include <escript/Data.h>

#include "DataTestCase.h"

#include <escript/DataLazy.h>
#include <escript/EscriptParams.h>
#include <escript/EsysException.h>
#include <escript/FunctionSpace.h>
#include <escript/TestDomain.h>

#include <cmath>
#include <cppunit/TestCaller.h>

#define AUTOLAZYON setEscriptParamInt("AUTOLAZY",1);
#define AUTOLAZYOFF setEscriptParamInt("AUTOLAZY",0);
#define CHECKAUTOLAZY (getEscriptParamInt("AUTOLAZY")==1)
#define SAVELAZYSTATE int LAZYSTATE=getEscriptParamInt("AUTOLAZY");
#define RESTORELAZYSTATE setEscriptParamInt("AUTOLAZY",LAZYSTATE);


using namespace std;
using namespace CppUnit;
using namespace escript;
using namespace escript::DataTypes;

namespace
{

inline
DataTypes::RealVectorType::const_reference
getRef(Data& d,int s1, int p1, int x, int y)
{
	return d.getDataAtOffsetRO(d.getDataOffset(s1,p1)+getRelIndex(d.getDataPointShape(),x,y), static_cast<DataTypes::real_t>(0));
}

inline
DataTypes::RealVectorType::const_reference
getRef(Data& d, int x, int y)
{
	return d.getDataAtOffsetRO(getRelIndex(d.getDataPointShape(),x,y), static_cast<DataTypes::real_t>(0));
}

}


void DataTestCase::testCopyingWorker(bool delayed)
{

  using namespace escript::DataTypes;
  cout << endl;
  TestDomain* tdp=new TestDomain(2,3,2);	// 2 points per sample, 3 samples, 2D coords
  Domain_ptr p(tdp);
  FunctionSpace fs=FunctionSpace(p, tdp->getContinuousFunctionCode());  
  DataTypes::ShapeType shape;
  shape.push_back(2);
  shape.push_back(3);
  DataTypes::RealVectorType data(DataTypes::noValues(shape),1);
  const int NUMDATS=3;
  Data* dats[NUMDATS];
  const char* strs[]={"DataConstant", "DataTagged", "DataExpanded"};
  dats[0]=new Data(new DataConstant(fs,shape,data));
  dats[1]=new Data(new DataTagged(fs,shape,data));
  dats[2]=new Data(new DataExpanded(fs,shape,data));
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
	  CPPUNIT_ASSERT(deep.isLazy());
	}
	if (!d->hasNoSamples())
	{
	    for (int i=0;i<DataTypes::noValues(shape);++i)
	    {
	      CPPUNIT_ASSERT(d->getDataAtOffsetRO(i, static_cast<DataTypes::real_t>(0))==deep.getDataAtOffsetRO(i, static_cast<DataTypes::real_t>(0)));
	    }
	}
	if (delayed)
	{
	   d->delaySelf();
	}
	d->setToZero();
	if (delayed)
	{
	  CPPUNIT_ASSERT(d->isLazy());
	}
	if (!d->hasNoSamples())
	{	
	    for (int i=0;i<DataTypes::noValues(shape);++i)
	    {
	      CPPUNIT_ASSERT(d->getDataAtOffsetRO(i, static_cast<DataTypes::real_t>(0))!=deep.getDataAtOffsetRO(i, static_cast<DataTypes::real_t>(0)));
	    }
	}
	if (delayed)
	{
	   d->delaySelf();
	   deep.delaySelf();
	}
	d->copy(deep);			// test copy from object
	if (delayed)
	{
	  CPPUNIT_ASSERT(d->isLazy());
	}
	if (!d->hasNoSamples())
	{	
	    for (int i=0;i<DataTypes::noValues(shape);++i)
	    {
	      CPPUNIT_ASSERT(d->getDataAtOffsetRO(i, static_cast<DataTypes::real_t>(0))==deep.getDataAtOffsetRO(i, static_cast<DataTypes::real_t>(0)));
	    }
	}
	d->setToZero();
	if (!d->hasNoSamples())
	{	
	    for (int i=0;i<DataTypes::noValues(shape);++i)
	    {
	      CPPUNIT_ASSERT(d->getDataAtOffsetRO(i, static_cast<DataTypes::real_t>(0))!=deep.getDataAtOffsetRO(i, static_cast<DataTypes::real_t>(0)));
	    }
	}
	delete dats[k];
  }
}


void DataTestCase::testSlicingWorker(bool delayed)
{

  using namespace escript::DataTypes;
  TestDomain* tdp=new TestDomain(2,3,2);	// 2 points per sample, 3 samples, 2D coords
  Domain_ptr p(tdp);
  FunctionSpace fs=FunctionSpace(p, tdp->getContinuousFunctionCode());    
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
    	dats[k]=new Data(1.3, viewShape, fs, false);
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
	dats[k]->requireWrite();
	if (!dats[k]->hasNoSamples())
	{
	    dats[k]->getDataAtOffsetRW(dats[k]->getDataOffset(0,0)+getRelIndex(viewShape,0,0), static_cast<DataTypes::real_t>(0))=1.0;
	    dats[k]->getDataAtOffsetRW(dats[k]->getDataOffset(0,0)+getRelIndex(viewShape,1,1), static_cast<DataTypes::real_t>(0))=2.0;
	}
    	DataTypes::RegionType region;
    	region.push_back(DataTypes::RegionType::value_type(0,0));
    	region.push_back(DataTypes::RegionType::value_type(0,0));

    	Data slice1(dats[k]->getSlice(region));

    	if (tags[k]) { CPPUNIT_ASSERT(slice1.isTagged()); }
	if (!slice1.hasNoSamples())
	{    	
	    CPPUNIT_ASSERT(slice1.getDataPointRank()==0);
	    CPPUNIT_ASSERT(slice1.getDataPointRO(0,0)==1.0);
	}    

	//
	// create a rank 2 slice with one value
	
	region.clear();
	region.push_back(DataTypes::RegionType::value_type(0,1));
	region.push_back(DataTypes::RegionType::value_type(0,1));
	
	Data slice2(dats[k]->getSlice(region));
	
	//cout << slice2.toString() << endl;
	
	if (tags[k]) {CPPUNIT_ASSERT(slice2.isTagged());}
	CPPUNIT_ASSERT(slice2.getDataPointRank()==2);
	
	if (!dats[k]->hasNoSamples())
	{
	    CPPUNIT_ASSERT(slice2.getDataAtOffsetRO(slice2.getDataOffset(0,0)+getRelIndex(slice2.getDataPointShape(),0,0), static_cast<DataTypes::real_t>(0))==1.0);
	}
	//
	// create a rank 2 slice with four values
	
	region.clear();
	region.push_back(DataTypes::RegionType::value_type(0,2));
	region.push_back(DataTypes::RegionType::value_type(0,2));
	
	Data slice3(dats[k]->getSlice(region));
	
	//cout << slice3.toString() << endl;
	
	if (tags[k]) {CPPUNIT_ASSERT(slice3.isTagged());}
	CPPUNIT_ASSERT(slice3.getDataPointRank()==2);
	if (!slice3.hasNoSamples())
	{
	    CPPUNIT_ASSERT(getRef(slice3,0,0,0,0)==1.0);
	    CPPUNIT_ASSERT(getRef(slice3,0,0,0,1)==1.3);
	    CPPUNIT_ASSERT(getRef(slice3,0,0,1,0)==1.3);
	    CPPUNIT_ASSERT(getRef(slice3,0,0,1,1)==2.0);
	}
   }

   // now some extra tests for tagged data (dats[1])

   //
   // add a value for tag "1"

   DataTypes::RealVectorType viewData(6);
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
   tdp->addUsedTag(1);
   vector<int> ts(3,1);
   tdp->assignTags(ts);
   
   CPPUNIT_ASSERT(slice4.isTagged());
   CPPUNIT_ASSERT(slice4.getDataPointRank()==2);
    if (!slice4.hasNoSamples())
    {
	CPPUNIT_ASSERT(getRef(slice4,0,0,0,0)==0);
	CPPUNIT_ASSERT(getRef(slice4,0,0,0,1)==2);
	CPPUNIT_ASSERT(getRef(slice4,0,0,0,2)==4);
	CPPUNIT_ASSERT(getRef(slice4,0,0,1,0)==1);
	CPPUNIT_ASSERT(getRef(slice4,0,0,1,1)==3);
	CPPUNIT_ASSERT(getRef(slice4,0,0,1,2)==5);
    }
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
 	dats[k]=new Data(1.3, viewShape,fs,false);
    	src[k]=new Data(10,DataTypes::scalarShape,fs,false);
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
	Data target(1.3,viewShape,fs,false);
	if (k==2) {target.expand();}
	DataTypes::RegionType region;
	region.push_back(DataTypes::RegionType::value_type(1,1));
	region.push_back(DataTypes::RegionType::value_type(1,1));
	target.setSlice(*(src[k]),region);
	if (!target.hasNoSamples())
	{
	    CPPUNIT_ASSERT(getRef(target,0,0,1,1)==src[k]->getDataPointRO(0,0));
	}
  }
  
  // some extra tests on tagged data

  //
  // add a value for tag "1" to target

  DataTypes::RealVectorType viewData(6);
  for (int i=0;i<viewData.size();i++) {
	viewData[i]=i;
  }

  Data target(1.3,viewShape,fs,false);
  target.tag();
  target.setTaggedValueFromCPP(1, viewShape, viewData);

    //cout << "target:\n" << target.toString() << endl;

    //
    // set a slice in target from source

  DataTypes::RegionType region;
  region.push_back(DataTypes::RegionType::value_type(0,0));
  region.push_back(DataTypes::RegionType::value_type(1,1));

  target.setSlice(*src[1],region);

  CPPUNIT_ASSERT(target.isTagged());
  CPPUNIT_ASSERT(target.getDataPointRank()==2);
  if (!target.hasNoSamples())
  {
      CPPUNIT_ASSERT(getRef(target,0,0,0,0)==0);
      CPPUNIT_ASSERT(getRef(target,0,0,0,1)==src[1]->getDataPointRO(0,0));
      CPPUNIT_ASSERT(getRef(target,0,0,0,2)==4);
      CPPUNIT_ASSERT(getRef(target,0,0,1,0)==1);
      CPPUNIT_ASSERT(getRef(target,0,0,1,1)==3);
      CPPUNIT_ASSERT(getRef(target,0,0,1,2)==5);
  }
  //
  // add a value for tag "2" to source

  DataTypes::ShapeType viewShape2;
  DataTypes::RealVectorType viewData2(1);
  viewData2[0]=6;
  src[1]->setTaggedValueFromCPP(2, viewShape2, viewData2);

  region.clear();
  region.push_back(DataTypes::RegionType::value_type(0,0));
  region.push_back(DataTypes::RegionType::value_type(1,1));

  target.setSlice(*src[1],region);

  CPPUNIT_ASSERT(target.isTagged());
  CPPUNIT_ASSERT(target.getDataPointRank()==2);

    // use a non-existent tag so we get a pointer to the default value
    // i.e.: the first element in the data array
  DataTypes::real_t* targetData=target.getSampleDataByTag(9);
  for (int i=0; i<target.getLength(); i++) {
      CPPUNIT_ASSERT(targetData[i]>=0);
  }
  CPPUNIT_ASSERT(targetData[0]==1.3);
  CPPUNIT_ASSERT(targetData[1]==1.3);
  CPPUNIT_ASSERT(targetData[2]==10);
  CPPUNIT_ASSERT(targetData[3]==1.3);
  CPPUNIT_ASSERT(targetData[4]==1.3);
  CPPUNIT_ASSERT(targetData[5]==1.3);
  CPPUNIT_ASSERT(targetData[6]==0);
  CPPUNIT_ASSERT(targetData[7]==1);
  CPPUNIT_ASSERT(targetData[8]==10);
  CPPUNIT_ASSERT(targetData[9]==3);
  CPPUNIT_ASSERT(targetData[10]==4);
  CPPUNIT_ASSERT(targetData[11]==5);
  CPPUNIT_ASSERT(targetData[12]==1.3);
  CPPUNIT_ASSERT(targetData[13]==1.3);
  CPPUNIT_ASSERT(targetData[14]==6);
  CPPUNIT_ASSERT(targetData[15]==1.3);
  CPPUNIT_ASSERT(targetData[16]==1.3);
  CPPUNIT_ASSERT(targetData[17]==1.3);


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

void DataTestCase::testSlicing()
{
  cout << "\n\tReadyData." << endl;
  testSlicingWorker(false);
  cout << "\n\tLazyData." << endl;
  testSlicingWorker(true);
}


void DataTestCase::testSomeDriver(bool autolazy)
{
  cout << endl;
  TestDomain* tdp=new TestDomain(2,3,2);	// 2 points per sample, 3 samples, 2D coords
  Domain_ptr p(tdp);
  FunctionSpace fs=FunctionSpace(p, tdp->getContinuousFunctionCode());    
  SAVELAZYSTATE
  if (autolazy)
  {
	AUTOLAZYON
	cout << "\tNow testing using autolazy." << endl;
  }
  cout << "\tCreate a Data object." << endl;

  DataTypes::ShapeType viewShape;
  viewShape.push_back(3);
  DataTypes::RealVectorType viewData(3);
  for (int i=0;i<viewShape[0];++i) {
    viewData[i]=i;
  }

  bool expanded=true;
  Data exData(viewData,viewShape, fs,expanded);
  Data cData(viewData,viewShape, fs,false);
  Data result;

  CPPUNIT_ASSERT(exData.isExpanded());
  CPPUNIT_ASSERT(cData.isConstant());
  CPPUNIT_ASSERT(result.isEmpty());

  cout << "\tTest some basic operations" << endl;
  result=exData*cData;
  cout << CHECKAUTOLAZY << " " << result.isLazy() << " " << result.isExpanded()<< endl;
  CPPUNIT_ASSERT(CHECKAUTOLAZY?result.isLazy():result.isExpanded());

  CPPUNIT_ASSERT(result.Lsup()==4);
  CPPUNIT_ASSERT(result.sup()==4);
  CPPUNIT_ASSERT(result.inf()==0);

  result=exData+cData;
  result=exData-cData;
  result=exData/cData;

  cout << "\tExercise wherePositive method" << endl;
  CPPUNIT_ASSERT(!exData.wherePositive().isEmpty());

  cout << "\tExercise copyWithMask method" << endl;
  exData.copyWithMask(result, exData.wherePositive());
  CPPUNIT_ASSERT(!exData.wherePositive().isEmpty());
  RESTORELAZYSTATE

}

void DataTestCase::testSome()
{
  testSomeDriver(false);
  testSomeDriver(true);
}



// This method tests to see if resolve() produces results of the correct type
void DataTestCase::testResolveType()
{
  cout << endl;
  cout << "\tTesting resolve()\n";
  TestDomain* tdp=new TestDomain(2,3,2);	// 2 points per sample, 3 samples, 2D coords
  Domain_ptr p(tdp);
  FunctionSpace fs=FunctionSpace(p, tdp->getContinuousFunctionCode());    
  DataTypes::ShapeType viewShape;
  viewShape.push_back(2);
  viewShape.push_back(3);
  viewShape.push_back(4);
  DataTypes::RealVectorType viewData(2*3*4);
  for (int i=0;i<DataTypes::noValues(viewShape);++i) {
    viewData[i]=i;
  }
  Data c1(viewData,viewShape,fs,false);
  Data t1(viewData,viewShape,fs,false);
  Data e1(viewData,viewShape,fs,false);
  t1.tag();
  e1.expand();
  c1.delaySelf();
  t1.delaySelf();
  e1.delaySelf();
  Data d1=c1+c1;
  CPPUNIT_ASSERT(d1.isLazy());
  CPPUNIT_ASSERT((d1.resolve(),d1.isConstant()));
  d1=c1+t1;
  CPPUNIT_ASSERT(d1.isLazy());
  CPPUNIT_ASSERT((d1.resolve(),d1.isTagged()));
  d1=t1+c1;
  CPPUNIT_ASSERT(d1.isLazy());
  CPPUNIT_ASSERT((d1.resolve(),d1.isTagged()));
  d1=t1+t1;
  CPPUNIT_ASSERT(d1.isLazy());
  CPPUNIT_ASSERT((d1.resolve(),d1.isTagged()));
  d1=c1+e1;
  CPPUNIT_ASSERT(d1.isLazy());
  CPPUNIT_ASSERT((d1.resolve(),d1.isExpanded()));
  d1=e1+c1;
  CPPUNIT_ASSERT(d1.isLazy());
  CPPUNIT_ASSERT((d1.resolve(),d1.isExpanded()));
  d1=e1+t1;
  CPPUNIT_ASSERT(d1.isLazy());
  CPPUNIT_ASSERT((d1.resolve(),d1.isExpanded()));
  d1=t1+e1;
  CPPUNIT_ASSERT(d1.isLazy());
  CPPUNIT_ASSERT((d1.resolve(),d1.isExpanded()));
  d1=e1+e1;
  CPPUNIT_ASSERT(d1.isLazy());
  CPPUNIT_ASSERT((d1.resolve(),d1.isExpanded()));
  cout << "\tTesting tag()\n";
  c1.tag();
  CPPUNIT_ASSERT(c1.isTagged());
  t1.tag();
  CPPUNIT_ASSERT(t1.isTagged());
  CPPUNIT_ASSERT_THROW(e1.tag(), DataException);
  cout << "\tTesting expand()\n";
  Data c2(viewData,viewShape,fs,false);
  Data t2(viewData,viewShape,fs,false);
  Data e2(viewData,viewShape,fs,false);
  t2.tag();
  e2.expand();
  c2.delaySelf();
  t2.delaySelf();
  e2.delaySelf();
  c2.expand();
  CPPUNIT_ASSERT(c2.isExpanded());
  t2.expand();
  CPPUNIT_ASSERT(t2.isExpanded());
  e2.expand();
  CPPUNIT_ASSERT(e2.isExpanded());
}

void DataTestCase::testDataConstant()
{
  cout << endl;
  cout << "\tCreate a DataConstant object." << endl;

  DataTypes::ShapeType viewShape;
  viewShape.push_back(2);
  viewShape.push_back(3);
  viewShape.push_back(4);
  DataTypes::RealVectorType viewData(2*3*4);
  for (int i=0;i<DataTypes::noValues(viewShape);++i) {
    viewData[i]=i;
  }

  Data left(viewData,viewShape,FunctionSpace(),false);
  Data right(viewData,viewShape,FunctionSpace(),false);
  Data result;

  cout << "\tTest some basic operations" << endl;

  result=left-right;

  CPPUNIT_ASSERT(left.isConstant());
  CPPUNIT_ASSERT(right.isConstant());
  CPPUNIT_ASSERT(result.isConstant());

  result=left+right;

  CPPUNIT_ASSERT(left.isConstant());
  CPPUNIT_ASSERT(right.isConstant());
  CPPUNIT_ASSERT(result.isConstant());

  CPPUNIT_ASSERT(!result.isExpanded());
  CPPUNIT_ASSERT(!result.isTagged());

}

void DataTestCase::testDataTagged()
{
  cout << endl;
  DataTypes::real_t dummyr=0;

  {

    cout << "\tCreate a DataTagged object with a default value only." << endl;

    DataTypes::ShapeType viewShape;
    viewShape.push_back(3);

    DataTypes::RealVectorType viewData(3);
    for (int i=0;i<viewShape[0];i++) {
      viewData[i]=i;
    }
    int arr[1]={1};		// iso c++ does not like empty arrays
    DataTagged* dt=new DataTagged(FunctionSpace(),viewShape,arr,viewData); 
    Data myData(dt);

    CPPUNIT_ASSERT(!myData.isEmpty());
    CPPUNIT_ASSERT(myData.isTagged());
    CPPUNIT_ASSERT(myData.getTagNumber(0)==1);
    CPPUNIT_ASSERT(myData.getDataPointRank()==1);
    

    cerr << "\n\n\n\n" << myData.getLength() << endl;
    cout << "\n\n\n\n" << myData.getLength() << endl;
    cout.flush();
    
    CPPUNIT_ASSERT(myData.getLength()==3);
    
    CPPUNIT_ASSERT(myData.getNoValues()==3);
    CPPUNIT_ASSERT(myData.getDataAtOffsetRO(0, dummyr)==0.0);
    CPPUNIT_ASSERT(myData.getDataAtOffsetRO(1, dummyr)==1.0);
    CPPUNIT_ASSERT(myData.getDataAtOffsetRO(2, dummyr)==2.0);

#ifdef EXWRITECHK		
    myData.requireWrite();
#endif	    
    double* sampleData=myData.getSampleDataRW(0, dummyr);
    for (int i=0; i<myData.getNoValues(); i++) {
      CPPUNIT_ASSERT(sampleData[i]==i);
    }
    // use a non-existent tag so we get a pointer to
    // the first element of the data array
    sampleData=myData.getSampleDataByTag(9);
    for (int i=0; i<myData.getLength(); i++) {
      CPPUNIT_ASSERT(sampleData[i]==i);
    }

    cout << "\tTest setting of a tag and associated value." << endl;

    // value for tag "1"
    DataTypes::RealVectorType eTwoData(viewData);
 //   DataArrayView eTwoView(eTwoData, viewShape);
    for (int i=0;i<viewShape[0];i++) {
      eTwoData[i]=i+2.0;
    }

    myData.setTaggedValueFromCPP(1,viewShape, eTwoData);

    CPPUNIT_ASSERT(myData.getLength()==6);

    int offset=myData.getDataOffset(0,0);
    CPPUNIT_ASSERT(offset==3);
    CPPUNIT_ASSERT(myData.getDataPointRank()==1);
    CPPUNIT_ASSERT(myData.getNoValues()==3);

    CPPUNIT_ASSERT(myData.getDataAtOffsetRO(offset+0, dummyr)==2);
    CPPUNIT_ASSERT(myData.getDataAtOffsetRO(offset+1, dummyr)==3);
    CPPUNIT_ASSERT(myData.getDataAtOffsetRO(offset+2, dummyr)==4);

    sampleData=myData.getSampleDataByTag(1);
    for (int i=0; i<myData.getNoValues(); i++) {
      CPPUNIT_ASSERT(sampleData[i]==i+2);
    }

  }

  {

    cout << "\tCreate a DataTagged object via tag() method." << endl;

    DataTypes::ShapeType viewShape;
    viewShape.push_back(2);
    viewShape.push_back(3);
    Data myData(1.3,viewShape,FunctionSpace(),false);
    myData.tag();

    CPPUNIT_ASSERT(!myData.isEmpty());
    CPPUNIT_ASSERT(myData.isTagged());
    CPPUNIT_ASSERT(myData.getTagNumber(0)==1);
    CPPUNIT_ASSERT(myData.getDataPointRank()==2);
    CPPUNIT_ASSERT(myData.getLength()==6);

    // check default value
    CPPUNIT_ASSERT(!myData.isEmpty());
    CPPUNIT_ASSERT(myData.getDataPointRank()==2);
    CPPUNIT_ASSERT(myData.getNoValues()==6);
    CPPUNIT_ASSERT(myData.getDataPointShape().size()==2);
    CPPUNIT_ASSERT(getRef(myData,0,0)==1.3);
    CPPUNIT_ASSERT(getRef(myData,0,1)==1.3);
    CPPUNIT_ASSERT(getRef(myData,0,2)==1.3);
    CPPUNIT_ASSERT(getRef(myData,1,0)==1.3);
    CPPUNIT_ASSERT(getRef(myData,1,1)==1.3);
    CPPUNIT_ASSERT(getRef(myData,1,2)==1.3);

    // check value for data-point (0,0).
//     myDataView = myData.getDataPoint(0,0);
    CPPUNIT_ASSERT(!myData.isEmpty());
//     CPPUNIT_ASSERT(myDataView.getOffset()==0);
    CPPUNIT_ASSERT(myData.getDataPointRank()==2);
    CPPUNIT_ASSERT(myData.getNoValues()==6);
    CPPUNIT_ASSERT(myData.getDataPointShape().size()==2);
    CPPUNIT_ASSERT(getRef(myData,0,0)==1.3);
    CPPUNIT_ASSERT(getRef(myData,0,1)==1.3);
    CPPUNIT_ASSERT(getRef(myData,0,2)==1.3);
    CPPUNIT_ASSERT(getRef(myData,1,0)==1.3);
    CPPUNIT_ASSERT(getRef(myData,1,1)==1.3);
    CPPUNIT_ASSERT(getRef(myData,1,2)==1.3);

  }

}

void DataTestCase::testDataTaggedExceptions()
{
  cout << endl;
  cout << "\tTest DataTagged exceptions." << endl;

  Data myData;

  CPPUNIT_ASSERT_THROW(myData.getSampleDataByTag(0), EsysException);
  CPPUNIT_ASSERT_THROW(myData.setTaggedValueFromCPP(0,DataTypes::ShapeType(), DataTypes::RealVectorType()), EsysException);
}

void DataTestCase::testConstructors()
{
  cout << endl;
  TestDomain* tdp=new TestDomain(2,3,2);	// 2 points per sample, 3 samples, 2D coords
  Domain_ptr p(tdp);
  FunctionSpace fs=FunctionSpace(p, tdp->getContinuousFunctionCode());    
  DataTypes::ShapeType viewShape;
  {
    cout << "\tCreate an Empty Data object" << endl;
    Data temp(1.3,viewShape,fs,false);
  }
  {
    cout << "\tCreate a rank 2 Data object" << endl;
    viewShape.push_back(2);
    viewShape.push_back(3);
    Data temp(1.3,viewShape,fs,false);
  }
}

void DataTestCase::testMoreOperations()
{
   cout << endl;
  TestDomain* tdp=new TestDomain(2,3,2);	// 2 points per sample, 3 samples, 2D coords
  Domain_ptr p(tdp);
  FunctionSpace fs=FunctionSpace(p, tdp->getContinuousFunctionCode());    
   
   DataTypes::ShapeType shape;
   shape.push_back(3);
   shape.push_back(3);

  // allocate the data 
  DataTypes::RealVectorType data(DataTypes::noValues(shape),0);

  // assign values to the data
  for (int i=0;i<shape[0];i++) {
    for (int j=0;j<shape[1];j++) {
      data[getRelIndex(shape,i,j)]=getRelIndex(shape,i,j);
    }
  }



  Data dats[]={Data(data,shape,fs,false),
		Data(data,shape,fs,false),
		Data(data,shape,fs,true),
		Data(data,shape,fs,false),
		Data(data,shape,fs,false),
		Data(data,shape,fs,true)};
  const int NUMDATS=6;
//  const int LAZY=3;		// where do the lazy objects start?

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

  CPPUNIT_ASSERT(baseEx.isExpanded());
  CPPUNIT_ASSERT(baseCon.isConstant());
  CPPUNIT_ASSERT(baseTag.isTagged());

  Data results[NUMDATS];
  double tmp;
  cout << "\tTest Data::trace(0)." << endl;
  for (int z=0;z<NUMDATS;++z)
  {
    results[z].copy(dats[z].trace(0));
  }
  for (int z=0;z<NUMDATS;++z)
  {
	tmp=0;
	if (!dats[z].hasNoSamples())
	{	
	    for (int i=0;i<shape[0];++i)
	    {
	      tmp+=getRef(dats[z],i,i);
	    }
	    CPPUNIT_ASSERT(std::abs(results[z].getDataAtOffsetRO(0, static_cast<DataTypes::real_t>(0)) - tmp) <= REL_TOL*std::abs(tmp));
	}
  }


}

void DataTestCase::testOperations()
{
  DataTypes::real_t dummyr=0;
  cout << endl;
  TestDomain* tdp=new TestDomain(1,1,1);	// 1 points per sample, 1 samples, 1D coords
  Domain_ptr p(tdp);
  FunctionSpace fs=FunctionSpace(p, tdp->getContinuousFunctionCode()); 
  
  // define the shape for the test data
  DataTypes::ShapeType shape;
  shape.push_back(2);
  shape.push_back(3);

  // allocate the data 
  DataTypes::RealVectorType data(DataTypes::noValues(shape),0);

  // assign values to the data
  for (int i=0;i<shape[0];i++) {
    for (int j=0;j<shape[1];j++) {
      data[getRelIndex(shape,i,j)]=getRelIndex(shape,i,j);
    }
  }

  DataTypes::RealVectorType data2(DataTypes::noValues(shape),0);
  // assign values to the data
  for (int i=0;i<shape[0];i++) {
    for (int j=0;j<shape[1];j++) {
      data2[getRelIndex(shape,i,j)]=12345678.8976+getRelIndex(shape,i,j);
    }
  }


  Data dats[]={Data(data,shape,fs,false),
		Data(data,shape,fs,false),
		Data(data,shape,fs,true),
		Data(data,shape,fs,false),
		Data(data,shape,fs,false),
		Data(data,shape,fs,true)
  };
  const int NUMDATS=6;
  const int LAZY=3;		// where do the lazy objects start?

  Data& baseCon=dats[0];
  Data& baseTag=dats[1];
  Data& baseEx=dats[2];
  baseTag.tag();
  dats[4].tag();
  dats[3].delaySelf();
  dats[4].delaySelf();
  dats[5].delaySelf();

  CPPUNIT_ASSERT(baseEx.isExpanded());
  CPPUNIT_ASSERT(baseCon.isConstant());
  CPPUNIT_ASSERT(baseTag.isTagged());

  Data results[NUMDATS];

  // create 0 <= smalldata <= 1 for testing trig functions

  DataTypes::RealVectorType smalldata(DataTypes::noValues(shape),0);

  // assign values to the data
  for (int i=0;i<shape[0];i++) {
    for (int j=0;j<shape[1];j++) {
      smalldata[getRelIndex(shape,i,j)]=(i==0 && j==0)?0:1.0/(getRelIndex(shape,i,j)+1);
    }
  }
  Data sdats[]={Data(smalldata,shape,fs,false),
		Data(smalldata,shape,fs,false),
		Data(smalldata,shape,fs,true),
		Data(smalldata,shape,fs,false),
		Data(smalldata,shape,fs,false),
		Data(smalldata,shape,fs,true)};
  sdats[1].tag();
  sdats[4].tag();
  sdats[3].delaySelf();		// 3 is a lazy constant
  sdats[4].delaySelf();		// 4 is a lazy tagged
  sdats[5].delaySelf();		// 5 is a lazy expanded



  // test unary operations

  double tmp;
  cout << "\tTest Data::pow." << endl;
  Data power(3.0,shape,fs,true);
  for (int z=0;z<NUMDATS;++z)
  {
    results[z].copy(dats[z].powD(power));
    if (z>=LAZY)
    {
	CPPUNIT_ASSERT(results[z].isLazy());
    }
  }
  for (int i=0;i<shape[0];i++) {
    for (int j=0;j<shape[1];j++) {
      tmp=pow((double)data[getRelIndex(shape,i,j)],(double)3.0);
cerr << tmp << endl;      
      for (int z=0;z<NUMDATS;++z)
      {
	results[z].resolve();
	if (!results[z].hasNoSamples())
	{
	    CPPUNIT_ASSERT(std::abs(getRef(results[z],i,j) - tmp) <= (REL_TOL*std::abs(tmp)+0.00000001));
	}
      }
    }
  }
cerr << "Ending pow" << endl;
  cout << "\tTest Data::sin." << endl;
  for (int z=0;z<NUMDATS;++z)
  {
    results[z].copy(dats[z].sin());
    if (z>=LAZY)
    {
	CPPUNIT_ASSERT(results[z].isLazy());
    }
  }
cerr << "Ending sin" << endl;  
  for (int i=0;i<shape[0];i++) {
    for (int j=0;j<shape[1];j++) {
      tmp=sin((double)data[getRelIndex(shape,i,j)]);
      for (int z=0;z<NUMDATS;++z)
      {
	results[z].resolve();	
	if (!results[z].hasNoSamples())
	{
	    CPPUNIT_ASSERT(std::abs(getRef(results[z],i,j) - tmp) <= REL_TOL*std::abs(tmp));
	}
      }
    }
  }

  cout << "\tTest Data::cos." << endl;
  for (int z=0;z<NUMDATS;++z)
  {
    results[z].copy(dats[z].cos());
    if (z>=LAZY)
    {
	CPPUNIT_ASSERT(results[z].isLazy());
    }
  }
  for (int i=0;i<shape[0];i++) {
    for (int j=0;j<shape[1];j++) {
      tmp=cos((double)data[getRelIndex(shape,i,j)]);
      for (int z=0;z<NUMDATS;++z)
      {
	results[z].resolve();	
	if (!results[z].hasNoSamples())
	{
	    CPPUNIT_ASSERT(std::abs(getRef(results[z],i,j) - tmp) <= REL_TOL*std::abs(tmp));
	}
      }
    }
  }

  cout << "\tTest Data::tan." << endl;
  for (int z=0;z<NUMDATS;++z)
  {
    results[z].copy(dats[z].tan());
    if (z>=LAZY)
    {
	CPPUNIT_ASSERT(results[z].isLazy());
    }
  }
  for (int i=0;i<shape[0];i++) {
    for (int j=0;j<shape[1];j++) {
      tmp=tan((double)data[getRelIndex(shape,i,j)]);
      for (int z=0;z<NUMDATS;++z)
      {
	results[z].resolve();
	if (!results[z].hasNoSamples())
	{
	    CPPUNIT_ASSERT(std::abs(getRef(results[z],i,j) - tmp) <= REL_TOL*std::abs(tmp));
	}
      }
    }
  }

  cout << "\tTest Data::asin." << endl;
  for (int z=0;z<NUMDATS;++z)
  {
    results[z].copy(sdats[z].asin());
    if (z>=LAZY)
    {
	CPPUNIT_ASSERT(results[z].isLazy());
    }
  }
  for (int i=0;i<shape[0];i++) {
    for (int j=0;j<shape[1];j++) {
      tmp=asin((double)smalldata[getRelIndex(shape,i,j)]);
      for (int z=0;z<NUMDATS;++z)
      {
	results[z].resolve();
	if (!results[z].hasNoSamples())
	{
	    CPPUNIT_ASSERT(std::abs(getRef(results[z],i,j) - tmp) <= REL_TOL*std::abs(tmp));
	}
      }
    }
  }

  cout << "\tTest Data::acos." << endl;
  for (int z=0;z<NUMDATS;++z)
  {
    results[z].copy(sdats[z].acos());
    if (z>=LAZY)
    {
	CPPUNIT_ASSERT(results[z].isLazy());
    }
  }
  for (int i=0;i<shape[0];i++) {
    for (int j=0;j<shape[1];j++) {
      tmp=acos((double)smalldata[getRelIndex(shape,i,j)]);
      for (int z=0;z<NUMDATS;++z)
      {
	results[z].resolve();	
	if (!results[z].hasNoSamples())
	{	
	    CPPUNIT_ASSERT(std::abs(getRef(results[z],i,j) - tmp) <= REL_TOL*std::abs(tmp));
	}
      }
    }
  }

  cout << "\tTest Data::atan." << endl;
  for (int z=0;z<NUMDATS;++z)
  {
    results[z].copy(sdats[z].atan());
    if (z>=LAZY)
    {
	CPPUNIT_ASSERT(results[z].isLazy());
    }
  }
  for (int i=0;i<shape[0];i++) {
    for (int j=0;j<shape[1];j++) {
      tmp=atan((double)smalldata[getRelIndex(shape,i,j)]);
      for (int z=0;z<NUMDATS;++z)
      {
	results[z].resolve();	
	if (!results[z].hasNoSamples())
	{	
	    CPPUNIT_ASSERT(std::abs(getRef(results[z],i,j) - tmp) <= REL_TOL*std::abs(tmp));
	}
      }
    }
  }

  cout << "\tTest Data::sinh." << endl;
  for (int z=0;z<NUMDATS;++z)
  {
    results[z].copy(dats[z].sinh());
    if (z>=LAZY)
    {
	CPPUNIT_ASSERT(results[z].isLazy());
    }
  }
  for (int i=0;i<shape[0];i++) {
    for (int j=0;j<shape[1];j++) {
      tmp=sinh((double)data[getRelIndex(shape,i,j)]);
      for (int z=0;z<NUMDATS;++z)
      {
	results[z].resolve();	
	if (!results[z].hasNoSamples())
	{
	    CPPUNIT_ASSERT(std::abs(getRef(results[z],i,j) - tmp) <= REL_TOL*std::abs(tmp));
	}
      }
    }
  }

  cout << "\tTest Data::cosh." << endl;
  for (int z=0;z<NUMDATS;++z)
  {
    results[z].copy(dats[z].cosh());
    if (z>=LAZY)
    {
	CPPUNIT_ASSERT(results[z].isLazy());
    }
  }
  for (int i=0;i<shape[0];i++) {
    for (int j=0;j<shape[1];j++) {
      tmp=cosh((double)data[getRelIndex(shape,i,j)]);
      for (int z=0;z<NUMDATS;++z)
      {
	results[z].resolve();	
	if (!results[z].hasNoSamples())
	{	
	    CPPUNIT_ASSERT(std::abs(getRef(results[z],i,j) - tmp) <= REL_TOL*std::abs(tmp));
	}
      }
    }
  }

  cout << "\tTest Data::tanh." << endl;
  for (int z=0;z<NUMDATS;++z)
  {
    results[z].copy(dats[z].tanh());
    if (z>=LAZY)
    {
	CPPUNIT_ASSERT(results[z].isLazy());
    }
  }
  for (int i=0;i<shape[0];i++) {
    for (int j=0;j<shape[1];j++) {
      tmp=tanh((double)data[getRelIndex(shape,i,j)]);
      for (int z=0;z<NUMDATS;++z)
      {
	results[z].resolve();	
	if (!results[z].hasNoSamples())
	{	
	    CPPUNIT_ASSERT(std::abs(getRef(results[z],i,j) - tmp) <= REL_TOL*std::abs(tmp));
	}    
      }
    }
  }

  // rather than accommodate the different windows operations directly I'll just use inverse functions
  cout << "\tTest Data::asinh." << endl;
  for (int z=0;z<NUMDATS;++z)
  {
    results[z].copy(dats[z].asinh().sinh());
    if (z>=LAZY)
    {
	CPPUNIT_ASSERT(results[z].isLazy());
    }
  }
  for (int i=0;i<shape[0];i++) {
    for (int j=0;j<shape[1];j++) {
      tmp=data[getRelIndex(shape,i,j)];
      for (int z=0;z<NUMDATS;++z)
      {
	results[z].resolve();	
	if (!results[z].hasNoSamples())
	{
	    CPPUNIT_ASSERT(std::abs(getRef(results[z],i,j) - tmp) <= REL_TOL*std::abs(tmp));
	}
      }
    }
  }

  cout << "\tTest Data::acosh." << endl;
  for (int z=0;z<NUMDATS;++z)
  {
    results[z].copy(dats[z].acosh().cosh());
    if (z>=LAZY)
    {
	CPPUNIT_ASSERT(results[z].isLazy());
    }
  }
  for (int i=0;i<shape[0];i++) {
    for (int j=0;j<shape[1];j++) {
      if (i==0 && j==0) break;
      tmp=data[getRelIndex(shape,i,j)];
      for (int z=0;z<NUMDATS;++z)
      {
	results[z].resolve();	
	if (!results[z].hasNoSamples())
	{	
	    CPPUNIT_ASSERT(std::abs(getRef(results[z],i,j) - tmp) <= REL_TOL*std::abs(tmp));
	}
      }
    }
  }

  cout << "\tTest Data::atanh." << endl;
  for (int z=0;z<NUMDATS;++z)
  {
    results[z].copy(dats[z].tanh().atanh());		// if these are the other way around the results are
    if (z>=LAZY)					// undefined
    {
	CPPUNIT_ASSERT(results[z].isLazy());
    }
  }
  for (int i=0;i<shape[0];i++) {
    for (int j=0;j<shape[1];j++) {
      tmp=data[getRelIndex(shape,i,j)];
      for (int z=0;z<NUMDATS;++z)
      {
	results[z].resolve();	
	if (!results[z].hasNoSamples())
	{	
	    CPPUNIT_ASSERT(std::abs(getRef(results[z],i,j) - tmp) <= REL_TOL*std::abs(tmp));
	}
      }
    }
  }

  cout << "\tTest Data::log." << endl;
  for (int z=0;z<NUMDATS;++z)
  {
    results[z].copy(dats[z].log());
    if (z>=LAZY)
    {
	CPPUNIT_ASSERT(results[z].isLazy());
    }
  }
  for (int i=0;i<shape[0];i++) {
    for (int j=0;j<shape[1];j++) {
      if (i==0 && j==0) break; 
      tmp=log((double)data[getRelIndex(shape,i,j)]);
      for (int z=0;z<NUMDATS;++z)
      {
	results[z].resolve();	
	if (!results[z].hasNoSamples())
	{	
	    CPPUNIT_ASSERT(std::abs(getRef(results[z],i,j) - tmp) <= REL_TOL*std::abs(tmp));
	}
      }
    }
  }

  cout << "\tTest Data::log10." << endl;
  for (int z=0;z<NUMDATS;++z)
  {
    results[z].copy(dats[z].log10());
    if (z>=LAZY)
    {
	CPPUNIT_ASSERT(results[z].isLazy());
    }
  }
  for (int i=0;i<shape[0];i++) {
    for (int j=0;j<shape[1];j++) {
      if (i==0 && j==0) break; 
      tmp=log10((double)data[getRelIndex(shape,i,j)]);
      for (int z=0;z<NUMDATS;++z)
      {
	results[z].resolve();	
	if (!results[z].hasNoSamples())
	{	
	    CPPUNIT_ASSERT(std::abs(getRef(results[z],i,j) - tmp) <= REL_TOL*std::abs(tmp));
	}
      }
    }
  }
  cout << "\tTest Data::erf." << endl;
  for (int z=0;z<NUMDATS;++z)
  {
    results[z].copy(dats[z].erf());
    if (z>=LAZY)
    {
	CPPUNIT_ASSERT(results[z].isLazy());
    }
  }
  for (int i=0;i<shape[0];i++) {
    for (int j=0;j<shape[1];j++) {
      if (i==0 && j==0) break; 
      tmp=erf((double)data[getRelIndex(shape,i,j)]);
      for (int z=0;z<NUMDATS;++z)
      {
	results[z].resolve();	
	if (!results[z].hasNoSamples())
	{	
	    CPPUNIT_ASSERT(std::abs(getRef(results[z],i,j) - tmp) <= REL_TOL*std::abs(tmp));
	}
      }
    }
  }


  cout << "\tTest Data::abs." << endl;
  for (int z=0;z<NUMDATS;++z)
  {
    results[z].copy(dats[z].abs());
    if (z>=LAZY)
    {
	CPPUNIT_ASSERT(results[z].isLazy());
    }
  }
  for (int i=0;i<shape[0];i++) {
    for (int j=0;j<shape[1];j++) {
      tmp=abs((double)data[getRelIndex(shape,i,j)]);
      for (int z=0;z<NUMDATS;++z)
      {
	results[z].resolve();	
	if (!results[z].hasNoSamples())
	{	
	    CPPUNIT_ASSERT(std::abs(getRef(results[z],i,j) - tmp) <= REL_TOL*std::abs(tmp));
	}
      }
    }
  }
  {
      Data inp(data2, shape, fs, true);
      Data res=inp.abs();
      CPPUNIT_ASSERT(res.inf()>12345678);
  }
  cout << "\tTest Data::sign (positive)." << endl;
  for (int z=0;z<NUMDATS;++z)
  {
    results[z].copy(dats[z].sign());
    if (z>=LAZY)
    {
	CPPUNIT_ASSERT(results[z].isLazy());
    }
  }
  for (int i=0;i<shape[0];i++) {
    for (int j=0;j<shape[1];j++) {
      tmp=(i==0 && j==0)?0:1;
      for (int z=0;z<NUMDATS;++z)
      {
	results[z].resolve();	
	if (!results[z].hasNoSamples())
	{	
	    CPPUNIT_ASSERT(std::abs(getRef(results[z],i,j) - tmp) <= REL_TOL*std::abs(tmp));
	}
      }
    }
  } 

  cout << "\tTest Data::sign (negative)." << endl;
  for (int z=0;z<NUMDATS;++z)
  {
    results[z].copy(dats[z].neg().sign());
    if (z>=LAZY)
    {
	CPPUNIT_ASSERT(results[z].isLazy());
    }
  }
  for (int i=0;i<shape[0];i++) {
    for (int j=0;j<shape[1];j++) {
      tmp=(i==0 && j==0)?0:-1;
      for (int z=0;z<NUMDATS;++z)
      {
	results[z].resolve();		
	if (!results[z].hasNoSamples())
	{	
	    CPPUNIT_ASSERT(std::abs(getRef(results[z],i,j) - tmp) <= REL_TOL*std::abs(tmp));
	}
      }
    }
  } 


  cout << "\tTest Data::exp." << endl;
  for (int z=0;z<NUMDATS;++z)
  {
    results[z].copy(dats[z].exp());
    if (z>=LAZY)
    {
	CPPUNIT_ASSERT(results[z].isLazy());
    }
  }
  for (int i=0;i<shape[0];i++) {
    for (int j=0;j<shape[1];j++) {
      tmp=exp((double)data[getRelIndex(shape,i,j)]);
      for (int z=0;z<NUMDATS;++z)
      {
	results[z].resolve();		
	if (!results[z].hasNoSamples())
	{	
	    CPPUNIT_ASSERT(std::abs(getRef(results[z],i,j) - tmp) <= REL_TOL*std::abs(tmp));
	}
      }
    }
  }

  cout << "\tTest Data::sqrt." << endl;
  for (int z=0;z<NUMDATS;++z)
  {
    results[z].copy(dats[z].sqrt());
    if (z>=LAZY)
    {
	CPPUNIT_ASSERT(results[z].isLazy());
    }
  }
  for (int i=0;i<shape[0];i++) {
    for (int j=0;j<shape[1];j++) {
      tmp=sqrt((double)data[getRelIndex(shape,i,j)]);
      for (int z=0;z<NUMDATS;++z)
      {
	results[z].resolve();		
	if (!results[z].hasNoSamples())
	{	
	    CPPUNIT_ASSERT(std::abs(getRef(results[z],i,j) - tmp) <= REL_TOL*std::abs(tmp));
	}
      }
    }
  }

  cout << "\tTest Data::neg." << endl;
  for (int z=0;z<NUMDATS;++z)
  {
    results[z].copy(dats[z].neg());
    if (z>=LAZY)
    {
	CPPUNIT_ASSERT(results[z].isLazy());
    }
  }
  for (int i=0;i<shape[0];i++) {
    for (int j=0;j<shape[1];j++) {
      tmp=-data[getRelIndex(shape,i,j)];
      for (int z=0;z<NUMDATS;++z)
      {
	results[z].resolve();		
	if (!results[z].hasNoSamples())
	{	
	    CPPUNIT_ASSERT(std::abs(getRef(results[z],i,j) - tmp) <= REL_TOL*std::abs(tmp));
	}
      }
    }
  }

  cout << "\tTest Data::pos." << endl;
  for (int z=0;z<NUMDATS;++z)
  {
      if (!dats[z].isLazy())     // python handles this differently
      {
            results[z].copy(dats[z].pos());
            if (z>=LAZY)
            {
            CPPUNIT_ASSERT(results[z].isLazy());
            }
      }
  }
  for (int i=0;i<shape[0];i++) {
    for (int j=0;j<shape[1];j++) {
      for (int z=0;z<NUMDATS;++z)
      {
          if (!dats[z].isLazy())
          {
                results[z].resolve();		
                if (!results[z].hasNoSamples())
                {	
                    CPPUNIT_ASSERT(std::abs(getRef(results[z],i,j) - getRelIndex(shape,i,j)) <= REL_TOL*std::abs(data[getRelIndex(shape,i,j)]));
                }
          }
      }
    }
  }

  // test reduction operations

  cout << "\tTest Data::Lsup." << endl;
  for (int z=0;z<NUMDATS;++z)
  {
    CPPUNIT_ASSERT(std::abs(dats[z].Lsup() - 5) <= REL_TOL*5);
  }

  cout << "\tTest Data::sup." << endl;
  for (int z=0;z<NUMDATS;++z)
  {
    CPPUNIT_ASSERT(std::abs(dats[z].sup() - 5) <= REL_TOL*5);
  }

  cout << "\tTest Data::inf." << endl;
  for (int z=0;z<NUMDATS;++z)
  {
    CPPUNIT_ASSERT(std::abs(dats[z].inf() - 0) <= REL_TOL*0);
  }

  // test data-point reduction operations

  cout << "\tTest Data::minval." << endl;
  for (int z=0;z<NUMDATS;++z)
  {
    results[z].copy(dats[z].minval());
  }
  for (int z=0;z<NUMDATS;++z)
  {
      results[z].resolve();		
      if (!results[z].hasNoSamples())
      {    
	  CPPUNIT_ASSERT(std::abs(results[z].getDataAtOffsetRO(0, dummyr) - 0) <= REL_TOL*0); 
      }
  }
  

  cout << "\tTest Data::maxval." << endl;
  for (int z=0;z<NUMDATS;++z)
  {
    results[z].copy(dats[z].maxval());
  }
  for (int z=0;z<NUMDATS;++z)
  {
    results[z].resolve();		
    if (!results[z].hasNoSamples())
    {    
	CPPUNIT_ASSERT(std::abs(results[z].getDataAtOffsetRO(0, dummyr) - 5) <= REL_TOL*5);
    }
  }

  cout << "\tTest Data::whereZero." << endl;
  for (int z=0;z<NUMDATS;++z)
  {
    results[z].copy(dats[z].whereZero(2));
  }
  for (int i=0;i<shape[0];i++) {
    for (int j=0;j<shape[1];j++) {
      tmp=(getRelIndex(shape,i,j)<=2);
      for (int z=0;z<NUMDATS;++z)
      {
	results[z].resolve();		
	if (!results[z].hasNoSamples())
	{	
	    CPPUNIT_ASSERT(std::abs(getRef(results[z],i,j) - tmp) <= REL_TOL*std::abs(tmp));
	}
      }
    }
  } 

  cout << "\tTest Data::whereNonZero." << endl;
  for (int z=0;z<NUMDATS;++z)
  {
    results[z].copy(dats[z].whereNonZero(2));
  }
  for (int i=0;i<shape[0];i++) {
    for (int j=0;j<shape[1];j++) {
      tmp=!(getRelIndex(shape,i,j)<=2);
      for (int z=0;z<NUMDATS;++z)
      {
	results[z].resolve();		
	if (!results[z].hasNoSamples())
	{	
	    CPPUNIT_ASSERT(std::abs(getRef(results[z],i,j) - tmp) <= REL_TOL*std::abs(tmp));
	}
      }
    }
  } 

  cout << "\tTest Data::transpose(1)." << endl;
  for (int z=0;z<NUMDATS;++z)
  {
    results[z].copy(dats[z].transpose(1));
  }
  for (int i=0;i<shape[0];i++) {
    for (int j=0;j<shape[1];j++) {
     for (int z=0;z<NUMDATS;++z)
     {
	results[z].resolve();		
	if (!results[z].hasNoSamples())
	{       
	    tmp=getRef(dats[z],i,j);
	    CPPUNIT_ASSERT(std::abs(getRef(results[z],j,i) - tmp) <= REL_TOL*std::abs(tmp));
	}
     }
    }
  } 

  cout << "\tTest Data::swapaxes(0,1)." << endl;
  for (int z=0;z<NUMDATS;++z)
  {
    results[z].copy(dats[z].swapaxes(0,1));
  }
  for (int i=0;i<shape[0];i++) {
    for (int j=0;j<shape[1];j++) {
     for (int z=0;z<NUMDATS;++z)
     {
	results[z].resolve();		
	if (!results[z].hasNoSamples())
	{       
	    tmp=getRef(dats[z],i,j);
	    CPPUNIT_ASSERT(std::abs(getRef(results[z],j,i) - tmp) <= REL_TOL*std::abs(tmp));
	}
     }
    }
  } 
}


// Here we test the binary operators in complex expressions
void DataTestCase::testBinary()
{
  DataTypes::real_t dummyr=0;
  cout << endl;

  // define the shape for the test data
  DataTypes::ShapeType shape;
  shape.push_back(2);
  shape.push_back(3);

  // allocate the data 
  DataTypes::RealVectorType data(DataTypes::noValues(shape),0);

  // assign values to the data
  for (int i=0;i<shape[0];i++) {
    for (int j=0;j<shape[1];j++) {
      data[getRelIndex(shape,i,j)]=getRelIndex(shape,i,j)+2;	// so we get no zeros
    }
  }


  Data one(1.0,DataTypes::scalarShape,FunctionSpace(),false);
  Data two(2.0,DataTypes::scalarShape,FunctionSpace(),false);
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
	  CPPUNIT_ASSERT(r1.isLazy() && r2.isLazy() && r3.isLazy() && r4.isLazy());
	}
	for (int i=0;i<DataTypes::noValues(shape);++i)
	{
	  CPPUNIT_ASSERT(std::abs(r1.getDataAtOffsetRO(i, dummyr)-data[i]) <= REL_TOL*data[i]);
	  CPPUNIT_ASSERT(std::abs(r2.getDataAtOffsetRO(i, dummyr)-data[i]) <= REL_TOL*data[i]);
	  CPPUNIT_ASSERT(std::abs(r3.getDataAtOffsetRO(i, dummyr)-data[i]) <= REL_TOL*data[i]);
	  CPPUNIT_ASSERT(std::abs(r4.getDataAtOffsetRO(i, dummyr)-pow(data[i],i)) <=REL_TOL*pow(data[i],i));
	}
  }
}


void DataTestCase::testComplexSamples()
{
    FunctionSpace fs=getTestDomainFunctionSpace(4,1,1);	// 4 points per sample, there is one sample and each point has one value in it
    Data x(5, DataTypes::scalarShape, fs, false);
    x.complicate();
    
    const DataTypes::cplx_t* r=x.getSampleDataRO(0, DataTypes::cplx_t(0));
    CPPUNIT_ASSERT(r[0]==DataTypes::cplx_t(5,0));
    
    RealVectorType v(1);
    Data t(0,DataTypes::scalarShape, fs, false);
    t.tag();
    for (int i=1;i<5;++i)
    {
        v[0]=i;
        t.setTaggedValueFromCPP(i,DataTypes::scalarShape, v);
    }	
    t.complicate();
    for (int i=1;i<5;++i)
    {
	CPPUNIT_ASSERT(t.getSampleDataByTag(i,DataTypes::cplx_t(0))[0]==DataTypes::cplx_t(i,0));
    }
}

void DataTestCase::testMemAlloc()
{
  //
  // Simple little sanity check for the memory allocator

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

TestSuite* DataTestCase::suite()
{
  // create the suite of tests to perform.
  TestSuite *testSuite = new TestSuite("DataTestCase");
  testSuite->addTest(new TestCaller<DataTestCase>(
              "testComplexSamples",&DataTestCase::testComplexSamples));
  testSuite->addTest(new TestCaller<DataTestCase>(
              "testCopying",&DataTestCase::testCopying));
  testSuite->addTest(new TestCaller<DataTestCase>(
              "testSome",&DataTestCase::testSome));
  testSuite->addTest(new TestCaller<DataTestCase>(
              "testDataConstant",&DataTestCase::testDataConstant));
  testSuite->addTest(new TestCaller<DataTestCase>(
              "testDataTagged",&DataTestCase::testDataTagged));
  testSuite->addTest(new TestCaller<DataTestCase>(
              "testDataTaggedExceptions",&DataTestCase::testDataTaggedExceptions));
  testSuite->addTest(new TestCaller<DataTestCase>(
              "testConstructors",&DataTestCase::testConstructors));
  testSuite->addTest(new TestCaller<DataTestCase>(
              "testSlicing",&DataTestCase::testSlicing));
  testSuite->addTest(new TestCaller<DataTestCase>(
              "testOperations",&DataTestCase::testOperations));
  testSuite->addTest(new TestCaller<DataTestCase>(
              "testMoreOperations",&DataTestCase::testMoreOperations));
  testSuite->addTest(new TestCaller<DataTestCase>(
              "testMemAlloc",&DataTestCase::testMemAlloc));
  testSuite->addTest(new TestCaller<DataTestCase>(
              "Resolving",&DataTestCase::testResolveType));
  
  return testSuite;
}

