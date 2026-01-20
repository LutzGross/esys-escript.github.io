
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

#include <escript/DataTypes.h>
#include "DataCombinationsTestCase.h"

#include <escript/Data.h>
#include <escript/TestDomain.h>
#include <cppunit/TestCaller.h>

#include <escript/Utils.h>

using namespace escript;
using namespace DataTypes;
using namespace std;

namespace
{
  
inline
DataTypes::RealVectorType::const_reference
getRef(Data& d, int x, int y)
{
	return d.getDataAtOffsetRO(getRelIndex(d.getDataPointShape(),x,y), static_cast<DataTypes::real_t>(0));
}  


DataTypes::RealVectorType getVector(const ShapeType& shape, double seed)
{
    DataTypes::RealVectorType data(DataTypes::noValues(shape),0);

    // assign values to the data
    for (int i=0;i<shape[0];i++) {
      for (int j=0;j<shape[1];j++) {
	data[getRelIndex(shape,i,j)]=getRelIndex(shape,i,j)+seed;	// so we get no zeros
      }
    }
    return data;
}  
  
  
Data getConstant(FunctionSpace fs, bool rank0, double seed)
{
    if (rank0)
    {
        return Data(seed,  DataTypes::ShapeType(),  fs,false);  
    }
    DataTypes::ShapeType shape;
    shape.push_back(2);
    shape.push_back(3);    

    DataTypes::RealVectorType data(DataTypes::noValues(shape),0);

    // assign values to the data
    for (int i=0;i<shape[0];i++) {
      for (int j=0;j<shape[1];j++) {
	data[getRelIndex(shape,i,j)]=getRelIndex(shape,i,j)+seed;	// so we get no zeros
      }
    }
    return Data(data, shape, fs, false);
}


Data getTagged(FunctionSpace fs, bool rank0, double seed, int tag1, int tag2, int tag3)
{
    if (rank0)
    {
	Data d(seed, DataTypes::ShapeType(), fs, false);
	d.tag();
	DataTypes::RealVectorType data(1,0);
	data[0]=seed*2;
	d.setTaggedValueFromCPP(tag1, DataTypes::ShapeType(), data);
	data[0]=seed*4;	
	d.setTaggedValueFromCPP(tag2, DataTypes::ShapeType(), data);
	data[0]=seed*8;		
	d.setTaggedValueFromCPP(tag3, DataTypes::ShapeType(), data);
	return d;
    }
    DataTypes::ShapeType shape;
    shape.push_back(2);
    shape.push_back(3);    

    DataTypes::RealVectorType data(DataTypes::noValues(shape),0);
    // assign values to the data
    for (int i=0;i<shape[0];i++) {
      for (int j=0;j<shape[1];j++) {
	data[getRelIndex(shape,i,j)]=getRelIndex(shape,i,j)+seed;	// so we get no zeros
      }
    }
    Data d(data, shape, fs, false);
    d.tag();
    for (int i=0;i<data.size();++i)
    {
	data[i]=data[i]*2;
    }
    d.setTaggedValueFromCPP(tag1,
		      shape,
		      data);
    for (int i=0;i<data.size();++i)
    {
	data[i]=data[i]*2;
    }
    d.setTaggedValueFromCPP(tag2,
		      shape,
		      data);
    for (int i=0;i<data.size();++i)
    {
	data[i]=data[i]*2;
    }
    d.setTaggedValueFromCPP(tag3,
		      shape,
		      data);
    return d;
}

Data getExpanded(FunctionSpace fs, bool rank0, double seed)
{
    if (rank0)
    {
	Data z(seed, DataTypes::ShapeType(), fs,false);
	Data d=fs.getDomain()->getX()*z;
	return d;     
    }
    DataTypes::ShapeType shape;
    shape.push_back(2);
    shape.push_back(3);    

    DataTypes::RealVectorType data(DataTypes::noValues(shape),0);
    // assign values to the data
    for (int i=0;i<shape[0];i++) {
      for (int j=0;j<shape[1];j++) {
	data[getRelIndex(shape,i,j)]=getRelIndex(shape,i,j)+seed;	// so we get no zeros
      }
    }
    Data d=fs.getDomain()->getX()*Data(data, shape, fs, true);
    return d;
}

void createConsts(FunctionSpace& fs, Data& c0s, Data& c4s, Data& c1, Data& c5)
{
    c0s=getConstant(fs, true, 0);
    c4s=getConstant(fs, true, 4);
  
    c1=getConstant(fs, false, 1);
    c5=getConstant(fs, false, 5); 
}

void createTagged(FunctionSpace& fs, Data& t1s, Data& t2s, Data& t3s, Data& t4s, Data& t5s, Data& t1, Data& t2, Data& t3, Data& t4, Data& t5)
{
  t1s=getTagged(fs, true, 1, 1, 2, 3);	// to check for problems with tag ordering
  t2s=getTagged(fs, true, 5, 1, 3, 2);
  t3s=getTagged(fs, true, 7, 3, 2, 1);
  t4s=getTagged(fs, true, 3, 2, 1, 3); 
  t5s=getTagged(fs, true, 8, 2, 4, 5);     // ensure we test mismatching tags
  
  t1=getTagged(fs, false, 1, 1, 2, 3);	// to check for problems with tag ordering
  t2=getTagged(fs, false, 5, 1, 2, 3);
  t3=getTagged(fs, false, 7, 3, 2, 1);
  t4=getTagged(fs, false, 3, 2, 1, 3);  
  t5=getTagged(fs, false, 9, 5, 2, 4);	// ensure we test mismatching tags  
}

void createExpand(FunctionSpace& fs, Data& es1, Data& es2, Data& e1, Data& e2)
{
    es1=getExpanded(fs, true, 1);
    es2=getExpanded(fs, true, 2);
  
    e1=getExpanded(fs, false, 1);
    e2=getExpanded(fs, false, 2);
}


}


void DataCombinationsTestCase::testUpdate()
{
  cout << endl;
  TestDomain* tdp=new TestDomain(2,3,1);	// 2 points per sample, 3 samples, 1D coords
  Domain_ptr p(tdp);
  FunctionSpace fs=FunctionSpace(p, tdp->getContinuousFunctionCode());  
  
  Data c0s;
  Data c4s;
  
  Data c1;
  Data c5;
  
  
  Data t1s;
  Data t2s;
  Data t3s;
  Data t4s; 
  Data t5s;
  
  Data t1;
  Data t2;
  Data t3;
  Data t4;  
  Data t5;
  
  Data es1;
  Data es2;
  
  Data e1;
  Data e2;  
  
  
  createConsts(fs, c0s, c4s, c1, c5);
  
  cout << "Constants (self):\n";
  c0s+=c4s;
  CPPUNIT_ASSERT(fabs(c0s.Lsup()-4)<0.01);
  
  c1+=c5;
  CPPUNIT_ASSERT(fabs(c1.Lsup()-16)<0.01);  
  
  createConsts(fs, c0s, c4s, c1, c5);

  c1+=c4s;
  CPPUNIT_ASSERT(fabs(c1.Lsup()-10)<0.01);  

  cout << "Tagged (self):\n";  
  createTagged(fs, t1s, t2s, t3s, t4s, t5s, t1, t2, t3, t4, t5); 
  
  t1+=t2;
  CPPUNIT_ASSERT(fabs(t1.inf()-6)<0.01);
  CPPUNIT_ASSERT(fabs(t1.Lsup()-16)<0.01);
  
  tdp->addUsedTag(1);
  CPPUNIT_ASSERT(fabs(t1.inf()-6)<0.01);
  CPPUNIT_ASSERT(fabs(t1.Lsup()-32)<0.01);
  
  tdp->addUsedTag(2);
  CPPUNIT_ASSERT(fabs(t1.inf()-6)<0.01);
  CPPUNIT_ASSERT(fabs(t1.Lsup()-64)<0.01);

  tdp->addUsedTag(3);
  CPPUNIT_ASSERT(fabs(t1.inf()-6)<0.01);
  CPPUNIT_ASSERT(fabs(t1.Lsup()-128)<0.01);
 
  tdp->clearUsedTags();
  
  t2s+=t1s;
  CPPUNIT_ASSERT(fabs(t2s.Lsup()-6)<0.01);
  
  tdp->addUsedTag(1);
  CPPUNIT_ASSERT(fabs(t2s.inf()-6)<0.01);
  CPPUNIT_ASSERT(fabs(t2s.Lsup()-12)<0.01);
  
  tdp->addUsedTag(2);
  CPPUNIT_ASSERT(fabs(t2s.Lsup()-44)<0.01);

  tdp->clearUsedTags();
  tdp->addUsedTag(3);
  CPPUNIT_ASSERT(fabs(t2s.Lsup()-28)<0.01);  
  
  createTagged(fs, t1s, t2s, t3s, t4s, t5s, t1, t2, t3, t4, t5); 

  tdp->clearUsedTags();
  t5+=t2s;

  CPPUNIT_ASSERT(fabs(t5.Lsup()-19)<0.01);
  CPPUNIT_ASSERT(fabs(t5.inf()-14)<0.01);
  tdp->addUsedTag(5);
  CPPUNIT_ASSERT(fabs(t5.Lsup()-33)<0.01);
  CPPUNIT_ASSERT(fabs(t5.inf()-14)<0.01);
  tdp->addUsedTag(3);
  
  CPPUNIT_ASSERT(fabs(t5.Lsup()-34)<0.01);
  CPPUNIT_ASSERT(fabs(t5.inf()-14)<0.01);  
  tdp->addUsedTag(4);  
  CPPUNIT_ASSERT(fabs(t5.Lsup()-117)<0.01);
  CPPUNIT_ASSERT(fabs(t5.inf()-14)<0.01);    
  tdp->clearUsedTags();  
  
  cout << "Expanded (self):\n";
  createExpand(fs, es1, es2, e1, e2);
  
  es2+=es1;
  CPPUNIT_ASSERT(fabs(es2.inf())<0.01);    
  CPPUNIT_ASSERT(fabs(es2.Lsup()-7.5)<0.01);
  
  e2+=e1;
  
  CPPUNIT_ASSERT(fabs(e2.inf())<0.01);    
  CPPUNIT_ASSERT(fabs(e2.Lsup()-32.5)<0.01);
  
  createExpand(fs, es1, es2, e1, e2);

  e1+=es2;
  CPPUNIT_ASSERT(fabs(e1.inf())<0.01);    
  CPPUNIT_ASSERT(fabs(e1.Lsup()-20)<0.01);
  
  cout << "Constant (self update by others):\n";
  createConsts(fs, c0s, c4s, c1, c5);
  createTagged(fs, t1s, t2s, t3s, t4s, t5s, t1, t2, t3, t4, t5);
  createExpand(fs, es1, es2, e1, e2);
  
  c4s+=t1s;
  CPPUNIT_ASSERT(fabs(c4s.Lsup()-5)<0.01);
  createConsts(fs, c0s, c4s, c1, c5);

  createConsts(fs, c0s, c4s, c1, c5);
  c5+=t1s;
  CPPUNIT_ASSERT(fabs(c5.Lsup()-11)<0.01);  
  
  createConsts(fs, c0s, c4s, c1, c5);
  c5+=t3;
  CPPUNIT_ASSERT(fabs(c5.Lsup()-22)<0.01);  

  createExpand(fs, es1, es2, e1, e2);
  c4s+=es1;
  CPPUNIT_ASSERT(fabs(c4s.Lsup()-6.5)<0.01);
  
  createConsts(fs, c0s, c4s, c1, c5);
  c5+=es1;
  CPPUNIT_ASSERT(fabs(c5.Lsup()-12.5)<0.01);
  
  createConsts(fs, c0s, c4s, c1, c5);
  c5+=e2;
  CPPUNIT_ASSERT(fabs(c5.Lsup()-27.5)<0.01);
  
  cout << "Tagged (self update by others):\n";
  createConsts(fs, c0s, c4s, c1, c5);
  createTagged(fs, t1s, t2s, t3s, t4s, t5s, t1, t2, t3, t4, t5);
  createExpand(fs, es1, es2, e1, e2);
  
  tdp->clearUsedTags();
  t3s+=c4s;
  CPPUNIT_ASSERT(fabs(t3s.Lsup()-11)<0.01); 
  tdp->addUsedTag(1);
  tdp->addUsedTag(2);
  tdp->addUsedTag(3);  
  CPPUNIT_ASSERT(fabs(t3s.Lsup()-60)<0.01);
  tdp->clearUsedTags();
  
  createTagged(fs, t1s, t2s, t3s, t4s, t5s, t1, t2, t3, t4, t5);
  t2+=c4s;
  CPPUNIT_ASSERT(fabs(t2.Lsup()-14)<0.01);  
  tdp->addUsedTag(1);
  tdp->addUsedTag(2);
  tdp->addUsedTag(3);  
  CPPUNIT_ASSERT(fabs(t2.Lsup()-84)<0.01); 
  tdp->clearUsedTags();
  
  
  t4+=c4s;
  CPPUNIT_ASSERT(fabs(t4.Lsup()-12)<0.01);  
  tdp->addUsedTag(1);
  tdp->addUsedTag(2);
  tdp->addUsedTag(3);  
  CPPUNIT_ASSERT(fabs(t4.Lsup()-68)<0.01);  
  tdp->clearUsedTags();
  
  
  createTagged(fs, t1s, t2s, t3s, t4s, t5s, t1, t2, t3, t4, t5);
  t3s+=es1;
  CPPUNIT_ASSERT(fabs(t3s.Lsup()-9.5)<0.01);  
  tdp->addUsedTag(1);
  tdp->addUsedTag(2);
  tdp->addUsedTag(3);  
  CPPUNIT_ASSERT(fabs(t3s.Lsup()-9.5)<0.01);  
  tdp->clearUsedTags();
  

  createTagged(fs, t1s, t2s, t3s, t4s, t5s, t1, t2, t3, t4, t5);
  t3+=es1;
  
  CPPUNIT_ASSERT(fabs(t3.Lsup()-14.5)<0.01);
  tdp->addUsedTag(1);
  tdp->addUsedTag(2);
  tdp->addUsedTag(3);  
  CPPUNIT_ASSERT(fabs(t3.Lsup()-14.5)<0.01);
  tdp->clearUsedTags();
  

  t5+=es2;
  CPPUNIT_ASSERT(fabs(t5.Lsup()-19)<0.01);  
  tdp->addUsedTag(1);
  tdp->addUsedTag(2);
  tdp->addUsedTag(3);  
  CPPUNIT_ASSERT(fabs(t5.Lsup()-19)<0.01);
  tdp->clearUsedTags();
  
  
  cout << "Expanded (self update by others):\n";
  createConsts(fs, c0s, c4s, c1, c5);
  createTagged(fs, t1s, t2s, t3s, t4s, t5s, t1, t2, t3, t4, t5);
  createExpand(fs, es1, es2, e1, e2);


  es1+=c4s;
  
  CPPUNIT_ASSERT(fabs(es1.Lsup()-6.5)<0.01);

  es2+=t4s;
  createExpand(fs, es1, es2, e1, e2);
  
  std::vector<int> t;
  t.push_back(1);
  t.push_back(1);
  t.push_back(0);  
  tdp->assignTags(t);
  es2+=t4s;
  CPPUNIT_ASSERT(fabs(es2.Lsup()-15)<0.01);   

  createExpand(fs, es1, es2, e1, e2);
  t[0]=2;
  t[1]=3;
  t[2]=1;
  tdp->assignTags(t);    
  es2+=t4s;
  CPPUNIT_ASSERT(fabs(es2.Lsup()-27)<0.01);   
  
  createExpand(fs, es1, es2, e1, e2);
  
  e2+=c4s;
  CPPUNIT_ASSERT(fabs(e2.Lsup()-21.5)<0.01);  
  
  e1+=c5;
  
  CPPUNIT_ASSERT(fabs(e1.Lsup()-25)<0.01);

  
  createExpand(fs, es1, es2, e1, e2);

  e2+=t3;

  CPPUNIT_ASSERT(fabs(e2.Lsup()-113.5)<0.01);
  

  t[0]=0;
  t[1]=1;
  t[2]=1;   
  tdp->assignTags(t);  
  e2+=t3;
  CPPUNIT_ASSERT(fabs(e2.Lsup()-209.5)<0.01);
 
}

// The purpose of this test is to check all the various combinations 
// of DataReady interactions
// These will test only the (Data + Data)
// hopefully it includes all relevant combinations
void DataCombinationsTestCase::testNonUpdate()
{

  cout << endl;
  TestDomain* tdp=new TestDomain(2,3,1);	// 2 points per sample, 3 samples, 1D coords
  Domain_ptr p(tdp);
  FunctionSpace fs=FunctionSpace(p, tdp->getContinuousFunctionCode());
  
    
  Data c0s=getConstant(fs, true, 0);
  Data c4s=getConstant(fs, true, 4);
  
  Data c1=getConstant(fs, false, 1);
  Data c5=getConstant(fs, false, 5);
  
  
  Data t1s=getTagged(fs, true, 1, 1, 2, 3);	// to check for problems with tag ordering
  Data t2s=getTagged(fs, true, 5, 1, 3, 2);
  Data t3s=getTagged(fs, true, 7, 3, 2, 1);
  Data t4s=getTagged(fs, true, 3, 2, 1, 3); 
  Data t5s=getTagged(fs, true, 8, 2, 4, 5);     // ensure we test mismatching tags
  
  Data t1=getTagged(fs, false, 1, 1, 2, 3);	// to check for problems with tag ordering
  Data t2=getTagged(fs, false, 5, 1, 2, 3);
  Data t3=getTagged(fs, false, 7, 3, 2, 1);
  Data t4=getTagged(fs, false, 3, 2, 1, 3);  
  Data t5=getTagged(fs, false, 9, 5, 2, 4);	// ensure we test mismatching tags
  
  Data es1=getExpanded(fs, true, 1);
  Data es2=getExpanded(fs, true, 2);
  
  Data e1=getExpanded(fs, false, 1);
  Data e2=getExpanded(fs, false, 2);
  
  real_t rr=0;
  Data res, res1, res2;
  cout << "Identical adds:\n";		// strictly, this only tests for inverse relationship between + and -
 
  rr=(c0s+c0s).Lsup();
  CPPUNIT_ASSERT(rr<0.01);
  
  res=c4s+c4s;
  rr=(res-c4s-c4s).Lsup();
  CPPUNIT_ASSERT(rr<0.01);
  
  res=(c1+c1);
  rr=(res-c1-c1).Lsup();
  CPPUNIT_ASSERT(rr<0.01);
  res=(c5+c5);
  rr=(res-c5-c5).Lsup();;
  CPPUNIT_ASSERT(rr<0.01);
  
  res=(t1s+t1s);
  rr=(res-t1s-t1s).Lsup();;
  CPPUNIT_ASSERT(rr<0.01);

  res=(t2s+t2s);
  rr=(res-t2s-t2s).Lsup();;
  CPPUNIT_ASSERT(rr<0.01);

  res=(t3s+t3s);
  rr=(res-t3s-t3s).Lsup();;
  CPPUNIT_ASSERT(rr<0.01);

  res=(t4s+t4s);
  rr=(res-t4s-t4s).Lsup();;
  CPPUNIT_ASSERT(rr<0.01);
  
  res=(es1+es1);
  rr=(res-es1-es1).Lsup();;
  CPPUNIT_ASSERT(rr<0.01);

  res=(es2+es2);
  rr=(res-es2-es2).Lsup();;
  CPPUNIT_ASSERT(rr<0.01);

  res=(e1+e1);
  rr=(res-e1-e1).Lsup();;
  CPPUNIT_ASSERT(rr<0.01);
  cout << "Rank0 constant adds:\n";
  res1=(c4s+c5);
  res2=(c5+c4s);
  rr=(res1-res2).Lsup();
  
  CPPUNIT_ASSERT(rr<0.01);

    // so the answers are the same, are they correct?
  
  DataTypes::ShapeType shape;
  shape.push_back(2);
  shape.push_back(3);        
  
  DataTypes::RealVectorType dat=getVector(shape, 9);  
  bool mismatch=false;
  for (int i=0;i<2;++i)
  {
      for (int j=0;j<3;++j)
      {
	  if (!res1.hasNoSamples())
	  {	
	      if (getRef(res1,i,j)!=dat[getRelIndex(shape,i,j)]) {
		  cout << "Mismatch at " << i << ',' << j << "::" << getRef(res1,i,j) << dat[getRelIndex(shape,i,j)] << endl;
		  mismatch=true;
	      }
	  }
      }
  }
  // need to do a global check to see if anyone mismatched
  assert(getMPIWorldMax(mismatch)==0);  
  
  res1=(c1+c5);
  res2=(c5+c1);
  rr=(res1-res2).Lsup();
  CPPUNIT_ASSERT(rr<0.01);
  dat=getVector(shape, 1);
  auto dat2=getVector(shape,5);
  mismatch=false;
  for (int i=0;i<2;++i)
  {
      for (int j=0;j<3;++j)
      {
	  if (!res1.hasNoSamples())
	  {
	      if (getRef(res1,i,j)!=(dat[getRelIndex(shape,i,j)]+dat2[getRelIndex(shape,i,j)])) {
		  cout << "Mismatch at " << i << ',' << j << endl;
		  mismatch=true;
	      }
	  }
      }
  }
  // need to do a global check to see if anyone mismatched
  assert(getMPIWorldMax(mismatch)==0);
  
  cout << "CCC looks good\n";

  // will test EEE next
  
  cout << "EEE\n";
  res1=e1+e2;
  res2=e2+e1;
  rr=(res1-res2).Lsup();
  CPPUNIT_ASSERT(rr<0.01);
  CPPUNIT_ASSERT(fabs(res1.Lsup()-32.5)<0.01);

  
  
  res1=es1+es2;
  res2=es2+es1;
  rr=(res1-res2).Lsup();
  CPPUNIT_ASSERT(rr<0.01);
  CPPUNIT_ASSERT(fabs(res1.Lsup()-7.5)<0.01);  

  res1=e1+es2;
  res2=es2+e1;
  rr=(res1-res2).Lsup();
  CPPUNIT_ASSERT(rr<0.01);
  CPPUNIT_ASSERT(fabs(res1.Lsup()-20)<0.01);    

  res1=es1+e2;
  res2=e2+es1;
  rr=(res1-res2).Lsup();
  CPPUNIT_ASSERT(rr<0.01);
  CPPUNIT_ASSERT(fabs(res1.Lsup()-20)<0.01);
  CPPUNIT_ASSERT(fabs(res1.inf())<0.000001);
  cout << "EEE looks good\n";

  cout << "E and C\n";
  res1=c5+es2;
  res2=es2+c5;
  rr=(res1-res2).Lsup();
  CPPUNIT_ASSERT(rr<0.01);
  CPPUNIT_ASSERT(fabs(res1.Lsup()-15)<0.01);
  CPPUNIT_ASSERT(fabs(res1.inf()-5)<0.000001);  

  res1=c4s+e2;
  res2=e2+c4s;
  rr=(res1-res2).Lsup();
  CPPUNIT_ASSERT(rr<0.01);
  CPPUNIT_ASSERT(fabs(res1.Lsup()-21.5)<0.01);
  CPPUNIT_ASSERT(fabs(res1.inf()-4)<0.000001);    
  
  res1=c5+e2;
  res2=e2+c5;
  rr=(res1-res2).Lsup();
  CPPUNIT_ASSERT(rr<0.01);
  CPPUNIT_ASSERT(fabs(res1.Lsup()-27.5)<0.01);
  CPPUNIT_ASSERT(fabs(res1.inf()-5)<0.000001); 
 
  cout << "TTT\n";
  
  res1=t1+t2;
  res2=t2+t1;
  rr=(res1-res2).Lsup();
  CPPUNIT_ASSERT(rr<0.01);
  CPPUNIT_ASSERT(fabs(res1.inf()-6)<0.001);
  CPPUNIT_ASSERT(fabs(res1.Lsup()-16)<0.001);  
  tdp->addUsedTag(2);
  CPPUNIT_ASSERT(fabs(res1.inf()-6)<0.001);
  CPPUNIT_ASSERT(fabs(res1.Lsup()-64)<0.001);
  tdp->addUsedTag(3);
  CPPUNIT_ASSERT(fabs(res1.inf()-6)<0.001);
  CPPUNIT_ASSERT(fabs(res1.Lsup()-128)<0.001);
  tdp->clearUsedTags();
  
  res1=t2+t3;		// tags in different orders
  res2=t3+t2;
  rr=(res1-res2).Lsup();
  CPPUNIT_ASSERT(rr<0.01);
  tdp->addUsedTag(3);
  CPPUNIT_ASSERT(fabs(res1.inf()-12)<0.001);
  CPPUNIT_ASSERT(fabs(res1.Lsup()-104)<0.001);
  tdp->clearUsedTags();

  res1=t5s+t4;
  res2=t4+t5s;
  rr=(res1-res2).Lsup();
  CPPUNIT_ASSERT(rr<0.01);
  tdp->addUsedTag(2);
  CPPUNIT_ASSERT(fabs(res1.inf()-11)<0.001);
  CPPUNIT_ASSERT(fabs(res1.Lsup()-32)<0.001);  
  tdp->clearUsedTags();
  tdp->addUsedTag(1);
  CPPUNIT_ASSERT(fabs(res1.inf()-11)<0.001);
  CPPUNIT_ASSERT(fabs(res1.Lsup()-40)<0.001);  
  tdp->clearUsedTags();
  tdp->addUsedTag(5);
  CPPUNIT_ASSERT(fabs(res1.inf()-11)<0.001);
  CPPUNIT_ASSERT(fabs(res1.Lsup()-72)<0.001);  
  tdp->clearUsedTags();
  
  res1=t2s+t5s;
  res2=t5s+t2s;
  rr=(res1-res2).Lsup();
  tdp->addUsedTag(1);
  CPPUNIT_ASSERT(fabs(res1.inf()-13)<0.001);
  CPPUNIT_ASSERT(fabs(res1.Lsup()-18)<0.001);    
  tdp->addUsedTag(3);
  CPPUNIT_ASSERT(fabs(res1.inf()-13)<0.001);
  CPPUNIT_ASSERT(fabs(res1.Lsup()-28)<0.001);      
  tdp->addUsedTag(5);
  CPPUNIT_ASSERT(fabs(res1.inf()-13)<0.001);
  CPPUNIT_ASSERT(fabs(res1.Lsup()-69)<0.001);    
  tdp->clearUsedTags();  
  
  cout << "TTT looks ok\n";
  cout << "T and C\n";
  
  
  res1=c4s+t4s;
  res2=t4s+c4s;
  CPPUNIT_ASSERT((res1-res2).Lsup()<0.01);
  CPPUNIT_ASSERT(fabs(res1.inf()-7)<0.001);
  CPPUNIT_ASSERT(fabs(res1.Lsup()-7)<0.001);    
  tdp->addUsedTag(2);
  CPPUNIT_ASSERT(fabs(res1.inf()-7)<0.001);
  CPPUNIT_ASSERT(fabs(res1.Lsup()-10)<0.001);      
  tdp->addUsedTag(3);
  CPPUNIT_ASSERT(fabs(res1.inf()-7)<0.001);
  CPPUNIT_ASSERT(fabs(res1.Lsup()-28)<0.001);     
  tdp->clearUsedTags();    
  
  res1=t4+c5;
  res2=c5+t4;
  CPPUNIT_ASSERT((res1-res2).Lsup()<0.01);
  CPPUNIT_ASSERT(fabs(res1.inf()-8)<0.001);
  CPPUNIT_ASSERT(fabs(res1.Lsup()-18)<0.001);    
  tdp->addUsedTag(2);
  CPPUNIT_ASSERT(fabs(res1.inf()-8)<0.001);
  CPPUNIT_ASSERT(fabs(res1.Lsup()-26)<0.001);      
  tdp->addUsedTag(3);
  CPPUNIT_ASSERT(fabs(res1.inf()-8)<0.001);
  CPPUNIT_ASSERT(fabs(res1.Lsup()-74)<0.001);     
  tdp->clearUsedTags();  
  
  res1=c4s+t4;
  res2=t4+c4s;
  CPPUNIT_ASSERT((res1-res2).Lsup()<0.01);
 
  CPPUNIT_ASSERT(fabs(res1.inf()-7)<0.001);
  CPPUNIT_ASSERT(fabs(res1.Lsup()-12)<0.001);    
  tdp->addUsedTag(1);
  
  CPPUNIT_ASSERT(fabs(res1.inf()-7)<0.001);
  CPPUNIT_ASSERT(fabs(res1.Lsup()-36)<0.001);      
  tdp->addUsedTag(3);
  CPPUNIT_ASSERT(fabs(res1.inf()-7)<0.001);
  CPPUNIT_ASSERT(fabs(res1.Lsup()-68)<0.001);     
  tdp->clearUsedTags();    
  
  res1=t4s+c5;
  res2=c5+t4s;
  CPPUNIT_ASSERT((res1-res2).Lsup()<0.01);
  CPPUNIT_ASSERT(fabs(res1.inf()-8)<0.001);
  CPPUNIT_ASSERT(fabs(res1.Lsup()-13)<0.001);       
  tdp->addUsedTag(3);
  CPPUNIT_ASSERT(fabs(res1.inf()-8)<0.001);
  CPPUNIT_ASSERT(fabs(res1.Lsup()-34)<0.001);     
  tdp->clearUsedTags();
  
  cout << "T and C look good\n";
  cout << "T and E\n";
  
  // We need to try various combinations of tags
  // with each combination of input types
  
  std::vector<int> t;
  t.push_back(0);
  t.push_back(1);
  t.push_back(0);
  
  tdp->assignTags(t);

  res1=t1+es2;
  res2=es2+t1;
  CPPUNIT_ASSERT((res1-res2).Lsup()<0.01);
  CPPUNIT_ASSERT(fabs(res1.inf()-1)<0.001);
  CPPUNIT_ASSERT(fabs(res1.Lsup()-15)<0.001);       

  res1=t1+e2;
  res2=e2+t1;  
  
  CPPUNIT_ASSERT((res1-res2).Lsup()<0.01);
  CPPUNIT_ASSERT(fabs(res1.inf()-1)<0.001);
  CPPUNIT_ASSERT(fabs(res1.Lsup()-23.5)<0.001);   
  
  res1=t1s+e2;
  res2=e2+t1s;
  CPPUNIT_ASSERT((res1-res2).Lsup()<0.01);
  CPPUNIT_ASSERT(fabs(res1.inf()-1)<0.001);
  CPPUNIT_ASSERT(fabs(res1.Lsup()-18.5)<0.001);     
  
  res1=t1s+es2;
  res2=es2+t1s;
  CPPUNIT_ASSERT((res1-res2).Lsup()<0.01);
  CPPUNIT_ASSERT(fabs(res1.inf()-1)<0.001);
  CPPUNIT_ASSERT(fabs(res1.Lsup()-6)<0.001);     

  cout << " round 0,1,0 OK\n";
  
  // now repeat with a different tag combination
  t[0]=2;
  t[1]=3;
  t[2]=1;
  tdp->assignTags(t);
  
  res1=t3+es2;
  res2=es2+t3;
  CPPUNIT_ASSERT((res1-res2).Lsup()<0.01);
  CPPUNIT_ASSERT(fabs(res1.inf()-16)<0.001);
  CPPUNIT_ASSERT(fabs(res1.Lsup()-101)<0.001);       
  
  res1=t3+e2;
  res2=e2+t3;  
  
  CPPUNIT_ASSERT((res1-res2).Lsup()<0.01);
  CPPUNIT_ASSERT(fabs(res1.inf()-16)<0.001);
  CPPUNIT_ASSERT(fabs(res1.Lsup()-113.5)<0.001);  
  
  res1=t3s+e2;
  res2=e2+t3s;
  CPPUNIT_ASSERT((res1-res2).Lsup()<0.01);
  CPPUNIT_ASSERT(fabs(res1.inf()-16)<0.001);
  CPPUNIT_ASSERT(fabs(res1.Lsup()-73.5)<0.001); 
  
  res1=t3s+es2;
  res2=es2+t3s;
  CPPUNIT_ASSERT((res1-res2).Lsup()<0.01);
  CPPUNIT_ASSERT(fabs(res1.inf()-16)<0.001);
  CPPUNIT_ASSERT(fabs(res1.Lsup()-61)<0.001);     

  cout << " round 2,3,1 OK\n";  
  
  t[0]=3;
  t[1]=1;
  t[2]=0;
  tdp->assignTags(t);
  
  res1=t3+es2;
  res2=es2+t3;
  CPPUNIT_ASSERT((res1-res2).Lsup()<0.01);
  CPPUNIT_ASSERT(fabs(res1.inf()-11)<0.001);
  CPPUNIT_ASSERT(fabs(res1.Lsup()-99)<0.001);       
  
  res1=t3+e2;
  res2=e2+t3;  
  CPPUNIT_ASSERT((res1-res2).Lsup()<0.01);
  CPPUNIT_ASSERT(fabs(res1.inf()-11)<0.001);
  CPPUNIT_ASSERT(fabs(res1.Lsup()-106.5)<0.001); 
  
  res1=t3s+e2;
  res2=e2+t3s;
  CPPUNIT_ASSERT((res1-res2).Lsup()<0.01);
  CPPUNIT_ASSERT(fabs(res1.inf()-11)<0.001);
  CPPUNIT_ASSERT(fabs(res1.Lsup()-66.5)<0.001); 
  
  res1=t3s+es2;
  res2=es2+t3s;
  CPPUNIT_ASSERT((res1-res2).Lsup()<0.01);
  CPPUNIT_ASSERT(fabs(res1.inf()-11)<0.001);
  CPPUNIT_ASSERT(fabs(res1.Lsup()-59)<0.001);     
  cout << " round 3,1,0 OK\n";   

  
}

using namespace CppUnit;

CppUnit::TestSuite* DataCombinationsTestCase::suite()
{
  // create the suite of tests to perform.
  CppUnit::TestSuite *testSuite = new TestSuite("DataCombinationsTestCase");
  testSuite->addTest(new TestCaller<DataCombinationsTestCase>(
              "testNonUpdate",&DataCombinationsTestCase::testNonUpdate));
  testSuite->addTest(new TestCaller<DataCombinationsTestCase>(
              "testUpdate",&DataCombinationsTestCase::testUpdate));
  
  return testSuite;
}
