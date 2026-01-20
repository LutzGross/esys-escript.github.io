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

// The purpose of these tests is to check for unwanted sharing of between Data objects

#include <escript/Data.h>
#include <escript/TestDomain.h>

#include "SharedDataTestCase.h"
#include <escript/EscriptParams.h>

#include <cppunit/TestCaller.h>
#include <iostream>

using namespace escript;
using namespace std;
using namespace CppUnit;
using namespace escript::DataTypes;

FunctionSpace getSharedFs()
{
    static FunctionSpace fs=getTestDomainFunctionSpace(1,50,1);
    return fs;
  
}


// Create a data, involve it in a lazy expression. Then modify the original
// and see if the value of the lazy is affected.
#define TESTEQOP(OP) { Data d((double)42,DataTypes::scalarShape, getSharedFs(),false); Data L=d.delay(); L-=Data((double)42,DataTypes::scalarShape, getSharedFs(),false); d OP Data(2,DataTypes::scalarShape, getSharedFs(),false); CPPUNIT_ASSERT(L.Lsup()<0.001);} 

// Test if the copy constructor shares a DataAbstract with its originator
void SharedDataTestCase::testEQ()
{
  cout << endl << "Testing +=" << flush;
  TESTEQOP(+=)
  cout << "\tOK" << endl << "Testing -=";
  TESTEQOP(-=)
  cout << "\tOK" << endl << "Testing *=";
  TESTEQOP(*=)
  cout << "\tOK" << endl << "Testing /=";
  TESTEQOP(/=)
  cout << "\tOK" << endl;
}

// Test for shared data caused by using a copy constructor
void SharedDataTestCase::testCC()
{
  cout << endl;
  Data d(42, DataTypes::scalarShape, getSharedFs(),false);
  Data shared(d);
  d+=Data(20,DataTypes::scalarShape, getSharedFs(),false);
  shared-=Data(42,DataTypes::scalarShape, getSharedFs(),false);
  CPPUNIT_ASSERT(shared.Lsup()<0.001);
}

// Test for shared data caused by using = operator
void SharedDataTestCase::testAssign()
{
  cout << endl;
  Data d(42, DataTypes::scalarShape, getSharedFs(),false);
  Data shared=d;
  d+=Data(20,DataTypes::scalarShape, getSharedFs(),false);
  shared-=Data(42,DataTypes::scalarShape, getSharedFs(),false);
  CPPUNIT_ASSERT(shared.Lsup()<0.001);
}

void SharedDataTestCase::testSetToZero()
{
  Data d((double)42,DataTypes::scalarShape, getSharedFs(),false); 
  Data L=d.delay(); 
  L-=Data((double)42,DataTypes::scalarShape, getSharedFs(),false);
  d.setToZero();
  CPPUNIT_ASSERT(L.Lsup()<0.001);
}

void SharedDataTestCase::testSetTaggedValueFromCPP()
{
  Data d((double)42,DataTypes::scalarShape, getSharedFs(),false);
  d.tag(); 
  Data L=d.delay();
  RealVectorType v(1,17);
  d.setTaggedValueFromCPP(1,DataTypes::scalarShape,v);
  L.resolve();
  // at this point, d should have a tag and L should not
  // unfortunately its a little tricky to find out what tags a Data object has so I'll use strings
  string s=L.toString();
  CPPUNIT_ASSERT(s.find("Tag(1)")==string::npos);		// if the tag shows up we have shared data
}

void SharedDataTestCase::testGetDataAtOffset()
{
  Data d((double)42,DataTypes::scalarShape, getSharedFs(),false);
  Data L=d.delay();
  // now change the data directly
  d.requireWrite();
  d.getDataAtOffsetRW(0, static_cast<DataTypes::real_t>(0))=17;
  CPPUNIT_ASSERT(L.getDataAtOffsetRO(0, static_cast<DataTypes::real_t>(0))==42);
}

void SharedDataTestCase::testGetDataPoint()
{
  Data d((double)42,DataTypes::scalarShape, getSharedFs(),false);
  Data L=d.delay();
  // now change the data directly
  d.requireWrite();
  d.getDataPointRW(0,0)=17;
  CPPUNIT_ASSERT(L.getDataPointRO(0,0)==42);
}

void SharedDataTestCase::testGetSampleRW()
{
  Data d((double)42,DataTypes::scalarShape, getSharedFs(),false);
  Data L=d.delay();
  
  std::cerr << "Please ignore the shared object message. We are testing the error checking.\n";
#ifdef SLOWSHARECHECK
  // now change the data directly
  CPPUNIT_ASSERT_THROW(*d.getSampleDataRW(0)=17, DataException);
#endif  
  std::cerr << "End ignore message\n";
  // Now try again properly 
  d.requireWrite();
  *d.getSampleDataRW(0)=17;
  L.resolve();
  CPPUNIT_ASSERT(*L.getSampleDataRO(0)==42);
}

TestSuite* SharedDataTestCase::suite()
{
  // create the suite of tests to perform.
  TestSuite *testSuite = new TestSuite("SharedDataTestCase");

  testSuite->addTest(new TestCaller<SharedDataTestCase>(
              "Arithmetic Assignment operators",&SharedDataTestCase::testEQ));
  testSuite->addTest(new TestCaller<SharedDataTestCase>(
              "Copy Constructor",&SharedDataTestCase::testCC));
  testSuite->addTest(new TestCaller<SharedDataTestCase>(
              "Assignment operator",&SharedDataTestCase::testAssign));
  testSuite->addTest(new TestCaller<SharedDataTestCase>(
              "setToZero",&SharedDataTestCase::testSetToZero));
  testSuite->addTest(new TestCaller<SharedDataTestCase>(
              "setTaggedValueFromCPP",&SharedDataTestCase::testSetTaggedValueFromCPP));
  testSuite->addTest(new TestCaller<SharedDataTestCase>(
              "getDataAtOffset",&SharedDataTestCase::testGetDataAtOffset));
  testSuite->addTest(new TestCaller<SharedDataTestCase>(
              "getDataPoint",&SharedDataTestCase::testGetDataPoint));
  testSuite->addTest(new TestCaller<SharedDataTestCase>(
              "getSampleRW",&SharedDataTestCase::testGetSampleRW));
  return testSuite;
}

