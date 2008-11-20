
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


#include "escript/DataConstant.h"
#include "escript/DataLazy.h"
#include "escript/FunctionSpace.h"
#include "esysUtils/EsysException.h"

#include "DataLazyTestCase.h"

#include <iostream>
#include <boost/shared_ptr.hpp>	// for the cast operator

using namespace CppUnitTest;
using namespace escript;
using namespace std;
using namespace esysUtils;
using namespace escript::DataTypes;
using namespace boost;

void DataLazyTestCase::setUp() {
  //
  // This is called before each test is run
}

void DataLazyTestCase::tearDown() {
  //
  // This is called after each test has been run
}


// This test file checks the basic properties of lazy data.
// It does not check the correctness of particular operations 


namespace
{

ValueType::reference
getRef(DataReady& data,int i, int j, int k)
{
   return data.getVector()[getRelIndex(data.getShape(),i,j,k)];
}


DataReady_ptr
resolveAndDelete(DataAbstract* p)
{
   DataReady_ptr p2=p->resolve();
   if (p!=p2.get())
   {
	delete p;
   }
   return p2;
}


DataAbstract_ptr
getLazy(DataTypes::ShapeType& shape,bool minus=false)
{
  int pts=DataTypes::noValues(shape);
  DataTypes::ValueType data(pts,0);
  for (int i=0;i<pts;++i)
  {
	data[i]=minus?-(i+1):i+1;
  }
  DataConstant* p=new DataConstant(FunctionSpace(),shape,data);
  DataAbstract_ptr pp(p);
  DataLazy* l=new DataLazy(pp);
  return DataAbstract_ptr(l);
}

DataAbstract_ptr
getLazyU(DataTypes::ShapeType& shape, ES_optype typ)
{
  int pts=DataTypes::noValues(shape);
  DataTypes::ValueType data(pts,0);
  for (int i=0;i<pts;++i)
  {
	data[i]=(i+1);
  }
  DataConstant* p=new DataConstant(FunctionSpace(),shape,data);
  DataAbstract_ptr pp(p);
  DataLazy* l=new DataLazy(pp,typ);
  return DataAbstract_ptr(l);
}


DataAbstract_ptr
getLazyB(DataTypes::ShapeType& shape, ES_optype typ)
{
  int pts=DataTypes::noValues(shape);
  DataTypes::ValueType data(pts,0);
  DataTypes::ValueType data2(pts,0);
  for (int i=0;i<pts;++i)
  {
	data[i]=(i+1);
	data2[i]=-(i+1);
  }
  DataConstant* p=new DataConstant(FunctionSpace(),shape,data);
  DataConstant* p2=new DataConstant(FunctionSpace(),shape,data2);
  DataAbstract_ptr pp(p);
  DataAbstract_ptr pp2(p2);
  DataLazy* l=new DataLazy(pp,pp2,typ);
  return DataAbstract_ptr(l);
}

#define TESTOP(X,V) { DataAbstract_ptr d1=getLazy(shape); assert(d1->X()==V); assert(d1->isLazy());}

}

// This method tests the identity constructor
void DataLazyTestCase::testLazy1()
{
  cout << endl;
  cout << "\tTesting IDENTITY constructor\n";

  DataTypes::ShapeType shape;
  DataAbstract_ptr d1=getLazy(shape);
  assert(d1->isLazy());
  assert(dynamic_pointer_cast<DataLazy>(d1)->getBuffsRequired()==1);

  for (int i=0;i<5;++i)
  {
    TESTOP(getRank,i);
    TESTOP(getNoValues,DataTypes::noValues(shape));
    TESTOP(getShape,shape);
    TESTOP(getNumDPPSample,1);
    TESTOP(getNumSamples,1);
    shape.push_back(3);
  }
}

#define TESTOPU(X,V,O) { DataAbstract_ptr d1=getLazyU(shape,O); assert(dynamic_pointer_cast<DataLazy>(d1)->X()==V); assert(d1->isLazy());}
// This method tests the unary op  constructor
// We aren't checking the correctness of the results here, just that they have the right properties
void DataLazyTestCase::testLazy2()
{
  cout << endl;
  cout << "\tTesting UNARY constructor (basic checks only)\n";

  DataTypes::ShapeType shape;
  DataAbstract_ptr d1=getLazyU(shape,LOG);
  assert(d1->isLazy());

  for (int j=SIN;j<=LEZ;++j)
  {
    shape=DataTypes::scalarShape;
    ES_optype op=(ES_optype)(j);			// not even reinterpret_cast works here
					// if other compilers object I'll write a switch 
    cout << "\t" << opToString(op) << endl;
    for (int i=0;i<5;++i)
    {
	TESTOPU(getRank,i,op);
    	TESTOPU(getNoValues,DataTypes::noValues(shape),op);
    	TESTOPU(getShape,shape,op);
    	TESTOPU(getNumDPPSample,1,op);
    	TESTOPU(getNumSamples,1,op);
	TESTOPU(getBuffsRequired,1,op);
    	shape.push_back(3);
    }
  }
}

#define TESTOPB(X,V,O) { DataAbstract_ptr d1=getLazyB(shape,O); assert(dynamic_pointer_cast<DataLazy>(d1)->X()==V); assert(dynamic_pointer_cast<DataLazy>(d1)->isLazy());}
// This method tests the unary op  constructor
// We aren't checking the correctness of the results here, just that they have the right properties
void DataLazyTestCase::testLazy3()
{
  cout << endl;
  cout << "\tTesting BINARY constructor (basic checks only)\n";

  DataTypes::ShapeType shape;
  DataAbstract_ptr d1=getLazyB(shape,ADD);
  assert(d1->isLazy());

  for (int j=ADD;j<=POW;++j)
  {
    shape=DataTypes::scalarShape;
    ES_optype op=(ES_optype)(j);			// not even reinterpret_cast works here
					// if other compilers object I'll write a switch 
    cout << "\t" << opToString(op) << endl;
    for (int i=0;i<5;++i)
    {
	TESTOPB(getRank,i,op);
    	TESTOPB(getNoValues,DataTypes::noValues(shape),op);
    	TESTOPB(getShape,shape,op);
    	TESTOPB(getNumDPPSample,1,op);
    	TESTOPB(getNumSamples,1,op);
	TESTOPB(getBuffsRequired,2,op);
    	shape.push_back(3);
    }
  }
}

void DataLazyTestCase::testBuffers()
{
  cout << endl;
  cout << "\tTesting Buffs required\n";
  DataTypes::ShapeType shape;
  DataAbstract_ptr p=(new DataLazy(getLazy(shape),getLazy(shape),ADD))->getPtr();
  DataAbstract_ptr p2=(new DataLazy(p,SIN))->getPtr();
  DataAbstract_ptr p3=(new DataLazy(p2,COS))->getPtr();
  DataAbstract_ptr p4=(new DataLazy(p3,getLazy(shape),ADD))->getPtr();
  assert(dynamic_pointer_cast<DataLazy>(p4)->getBuffsRequired()==2);
  DataAbstract_ptr p5=(new DataLazy(p2,p4,ADD))->getPtr();
  assert(dynamic_pointer_cast<DataLazy>(p5)->getBuffsRequired()==3);
}


TestSuite* DataLazyTestCase::suite ()
{
  //
  // create the suite of tests to perform.
  TestSuite *testSuite = new TestSuite ("DataLazyTestCase");

  testSuite->addTest (new TestCaller< DataLazyTestCase>("Identity",&DataLazyTestCase::testLazy1));
  testSuite->addTest (new TestCaller< DataLazyTestCase>("Unary",&DataLazyTestCase::testLazy2));
  testSuite->addTest (new TestCaller< DataLazyTestCase>("Binary",&DataLazyTestCase::testLazy3));
  testSuite->addTest (new TestCaller< DataLazyTestCase>("Buffers",&DataLazyTestCase::testBuffers));
  return testSuite;
}
