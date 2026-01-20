
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

#include <escript/DataConstant.h>
#include "DataLazyTestCase.h"

#include <escript/DataLazy.h>
#include <escript/FunctionSpace.h>

#include <iostream>
#include <cppunit/TestCaller.h>
#include <boost/shared_ptr.hpp>	// for the cast operator

using namespace CppUnit;
using namespace escript;
using namespace std;
using namespace escript::DataTypes;
using namespace boost;

// This test file checks the basic properties of lazy data.
// It does not check the correctness of particular operations 


namespace
{

//DataReady_ptr
//resolveAndDelete(DataAbstract* p)
//{
//   DataReady_ptr p2=p->resolve();
//   if (p!=p2.get())
//   {
//	delete p;
//   }
//   return p2;
//}


DataAbstract_ptr
getLazy(DataTypes::ShapeType& shape,bool minus=false)
{
  int pts=DataTypes::noValues(shape);
  DataTypes::RealVectorType data(pts,0);
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
  DataTypes::RealVectorType data(pts,0);
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
getLazyUP(DataTypes::ShapeType& shape, ES_optype typ, int par)
{
  int pts=DataTypes::noValues(shape);
  DataTypes::RealVectorType data(pts,0);
  for (int i=0;i<pts;++i)
  {
	data[i]=(i+1);
  }
  DataConstant* p=new DataConstant(FunctionSpace(),shape,data);
  DataAbstract_ptr pp(p);
  DataLazy* l=new DataLazy(pp,typ,par);
  return DataAbstract_ptr(l);
}


DataAbstract_ptr
getLazyB(DataTypes::ShapeType& shape, ES_optype typ)
{
  int pts=DataTypes::noValues(shape);
  DataTypes::RealVectorType data(pts,0);
  DataTypes::RealVectorType data2(pts,0);
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

DataAbstract_ptr
getLazyGTP(DataTypes::ShapeType& shape, ES_optype typ, int ax, int tr)
{
  int pts=DataTypes::noValues(shape);
  DataTypes::RealVectorType data(pts,0);
  DataTypes::RealVectorType data2(pts,0);
  for (int i=0;i<pts;++i)
  {
	data[i]=(i+1);
	data2[i]=-(i+1);
  }
  DataConstant* p=new DataConstant(FunctionSpace(),shape,data);
  DataConstant* p2=new DataConstant(FunctionSpace(),shape,data2);
  DataAbstract_ptr pp(p);
  DataAbstract_ptr pp2(p2);
  DataLazy* l=new DataLazy(pp,pp2,typ,ax,tr);
  return DataAbstract_ptr(l);
}



#define TESTOP(X,V) { DataAbstract_ptr d1=getLazy(shape); CPPUNIT_ASSERT(d1->X()==V); CPPUNIT_ASSERT(d1->isLazy());}

}

// This method tests the identity constructor
void DataLazyTestCase::testLazy1()
{
  cout << endl;
  cout << "\tTesting IDENTITY constructor\n";

  DataTypes::ShapeType shape;
  DataAbstract_ptr d1=getLazy(shape);
  CPPUNIT_ASSERT(d1->isLazy());

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

#define TESTOPU(X,V,O) { DataAbstract_ptr d1=getLazyU(shape,O); CPPUNIT_ASSERT(dynamic_pointer_cast<DataLazy>(d1)->X()==V); CPPUNIT_ASSERT(d1->isLazy());}
// This method tests the unary op  constructor
// We aren't checking the correctness of the results here, just that they have the right properties
void DataLazyTestCase::testLazy2()
{
  cout << endl;
  cout << "\tTesting UNARY constructor (basic checks only)\n";

  DataTypes::ShapeType shape;
  DataAbstract_ptr d1=getLazyU(shape,LOG);
  CPPUNIT_ASSERT(d1->isLazy());

  for (int j=SIN;j<=LEZ;++j)
  {
    shape=DataTypes::scalarShape;
    ES_optype op=(ES_optype)(j);			// not even reinterpret_cast works here
					// if other compilers object I'll write a switch 
    if (op==POS)
    {
        continue;       // not testing this, python handles it differently
    }
    cout << "\t" << opToString(op) << endl;
    for (int i=0;i<5;++i)
    {
	TESTOPU(getRank,i,op);
    	TESTOPU(getNoValues,DataTypes::noValues(shape),op);
    	TESTOPU(getShape,shape,op);
    	TESTOPU(getNumDPPSample,1,op);
    	TESTOPU(getNumSamples,1,op);
    	shape.push_back(3);
    }
  }
}

#define TESTOPUP(X,V,O) { DataAbstract_ptr d1=getLazyUP(shape,O,0); CPPUNIT_ASSERT(dynamic_pointer_cast<DataLazy>(d1)->X()==V); CPPUNIT_ASSERT(d1->isLazy());}
// This method tests the unary op  constructor
// We aren't checking the correctness of the results here, just that they have the right properties
void DataLazyTestCase::testLazy2p()
{
  cout << endl;
  cout << "\tTesting UNARY (with arg) constructor (basic checks only)\n";

  DataTypes::ShapeType shape;
  DataTypes::ShapeType traceshape;
  DataAbstract_ptr d1=getLazyUP(shape,TRANS,0);
  CPPUNIT_ASSERT(d1->isLazy());

  for (int j=TRANS;j<=TRACE;++j)
  {
    shape=DataTypes::scalarShape;	// traceshape is only used once so not initialised here
    
    ES_optype op=(ES_optype)(j);			// not even reinterpret_cast works here
					// if other compilers object I'll write a switch 
    cout << "\t" << opToString(op) << endl;
    for (int i=0;i<5;++i)
    {
	if (op==TRACE)
	{
	   if (i>1)	// trace only works 2 and up
	   {
	      TESTOPUP(getRank,i-2,op);
    	      TESTOPUP(getNoValues, DataTypes::noValues(traceshape),op);
    	      TESTOPUP(getShape,traceshape,op);

    	      TESTOPUP(getNumDPPSample,1,op);
    	      TESTOPUP(getNumSamples,1,op);
	      traceshape.push_back(3);
	   }
	}
	else
	{
	   TESTOPUP(getRank,i,op);
    	   TESTOPUP(getNoValues,DataTypes::noValues(shape),op);
    	   TESTOPUP(getShape,shape,op);
    	   TESTOPUP(getNumDPPSample,1,op);
    	   TESTOPUP(getNumSamples,1,op);
	}
    	shape.push_back(3);
    }
  }
}

#define TESTOPB(X,V,O) { DataAbstract_ptr d1=getLazyB(shape,O); CPPUNIT_ASSERT(dynamic_pointer_cast<DataLazy>(d1)->X()==V); CPPUNIT_ASSERT(dynamic_pointer_cast<DataLazy>(d1)->isLazy());}
// This method tests the binary op  constructor
// We aren't checking the correctness of the results here, just that they have the right properties
void DataLazyTestCase::testLazy3()
{
  cout << endl;
  cout << "\tTesting BINARY constructor (basic checks only)\n";

  DataTypes::ShapeType shape;
  DataAbstract_ptr d1=getLazyB(shape,ADD);
  CPPUNIT_ASSERT(d1->isLazy());

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
    	shape.push_back(3);
    }
  }
}




#define TESTOPGTP(X,V,O) { DataAbstract_ptr d1=getLazyGTP(shape,O,0,0); CPPUNIT_ASSERT(dynamic_pointer_cast<DataLazy>(d1)->X()==V); CPPUNIT_ASSERT(dynamic_pointer_cast<DataLazy>(d1)->isLazy());}

// This method tests the GeneralTensorproduct  constructor
// We aren't checking the correctness of the results here, just that they have the right properties
void DataLazyTestCase::testLazy4()
{
  cout << endl;
  cout << "\tTesting GTP constructor (basic checks only)\n";

  DataTypes::ShapeType shape;
  DataTypes::ShapeType prodshape;

  DataAbstract_ptr d1=getLazyGTP(shape,PROD,0,0);
  CPPUNIT_ASSERT(d1->isLazy());

  for (int j=PROD;j<=PROD;++j)
  {
    shape=DataTypes::scalarShape;
    ES_optype op=(ES_optype)(j);			// not even reinterpret_cast works here
					// if other compilers object I'll write a switch 
    cout << "\t" << opToString(op) << endl;
    for (int i=0;i<3;++i)
    {
	ShapeType ns;
	for (int k=0;k<i;++k)
	{
	  ns.push_back(3);
	  ns.push_back(3);
	}
	TESTOPGTP(getRank,i*2,op);
    	TESTOPGTP(getNoValues,DataTypes::noValues(ns),op);
    	TESTOPGTP(getShape,ns,op);
    	TESTOPGTP(getNumDPPSample,1,op);
    	TESTOPGTP(getNumSamples,1,op);
    	shape.push_back(3);
    }
  }
}



TestSuite* DataLazyTestCase::suite()
{
  // create the suite of tests to perform.
  TestSuite *testSuite = new TestSuite("DataLazyTestCase");

  testSuite->addTest(new TestCaller<DataLazyTestCase>(
              "Identity",&DataLazyTestCase::testLazy1));
  testSuite->addTest(new TestCaller<DataLazyTestCase>(
              "Unary",&DataLazyTestCase::testLazy2));
  testSuite->addTest(new TestCaller<DataLazyTestCase>(
              "Unary (params)",&DataLazyTestCase::testLazy2p));
  testSuite->addTest(new TestCaller<DataLazyTestCase>(
              "Binary",&DataLazyTestCase::testLazy3));
  testSuite->addTest(new TestCaller<DataLazyTestCase>(
              "GTP",&DataLazyTestCase::testLazy4));
  return testSuite;
}

