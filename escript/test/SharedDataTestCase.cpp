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

// The purpose of these tests is to check for unwanted sharing of between Data objects


#include "SharedDataTestCase.h"
#include "escript/Data.h"
#include "escript/EscriptParams.h"

#include <iostream>

using namespace escript;
using namespace std;
using namespace CppUnitTest;
using namespace escript::DataTypes;

void SharedDataTestCase::setUp() 
{
  //
  // This is called before each test is run
}

void SharedDataTestCase::tearDown() 
{
  //
  // This is called after each test has been run
}

// Create a data, involve it in a lazy expression. Then modify the original and see if the value of the lazy is affected.
#define TESTEQOP(OP) { Data d((double)42,DataTypes::scalarShape); Data L=d.delay(); L-=Data((double)42,DataTypes::scalarShape); d OP Data(2,DataTypes::scalarShape); assert(L.Lsup()<0.001);} 

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
}

// Test for shared data caused by using a copy constructor
void SharedDataTestCase::testCC()
{
  cout << endl;
  Data d(42, DataTypes::scalarShape);
  Data shared(d);
  d+=Data(20,DataTypes::scalarShape);
  shared-=Data(42,DataTypes::scalarShape);
  assert(shared.Lsup()<0.001);
}

// Test for shared data caused by using = operator
void SharedDataTestCase::testAssign()
{
  cout << endl;
  Data d(42, DataTypes::scalarShape);
  Data shared=d;
  d+=Data(20,DataTypes::scalarShape);
  shared-=Data(42,DataTypes::scalarShape);
  assert(shared.Lsup()<0.001);
}

void SharedDataTestCase::testSetToZero()
{
  Data d((double)42,DataTypes::scalarShape); 
  Data L=d.delay(); 
  L-=Data((double)42,DataTypes::scalarShape);
  d.setToZero();
  assert(L.Lsup()<0.001);
}

void SharedDataTestCase::testSetTaggedValueFromCPP()
{
  Data d((double)42,DataTypes::scalarShape);
  d.tag(); 
  Data L=d.delay();
  ValueType v(1,17);
  d.setTaggedValueFromCPP(1,DataTypes::scalarShape,v);
  L.resolve();
  // at this point, d should have a tag and L should not
  // unfortunately its a little tricky to find out what tags a Data object has so I'll use strings
  string s=L.toString();
  assert(s.find("Tag(1)")==string::npos);		// if the tag shows up we have shared data
}

void SharedDataTestCase::testGetDataAtOffset()
{
  Data d((double)42,DataTypes::scalarShape);
  Data L=d.delay();
  // now change the data directly
  d.getDataAtOffset(0)=17;
  assert(L.getDataAtOffset(0)==42);
}

void SharedDataTestCase::testGetDataPoint()
{
  Data d((double)42,DataTypes::scalarShape);
  Data L=d.delay();
  // now change the data directly
  d.getDataPoint(0,0)=17;
  assert(L.getDataPoint(0,0)==42);
}

void SharedDataTestCase::testGetSampleRW()
{
  Data d((double)42,DataTypes::scalarShape);
  Data L=d.delay();
  // now change the data directly
  try
  {
  	*d.getSampleDataRW(0)=17;
	assert(false);			// should have thrown 
  } catch (DataException e)
  {
  }
  // Now try again properly 
  d.requireWrite();
  *d.getSampleDataRW(0)=17;
  L.resolve();
  assert(*L.getSampleDataRO(0)==42);
}

TestSuite* SharedDataTestCase::suite ()
{
  //
  // create the suite of tests to perform.
  TestSuite *testSuite = new TestSuite ("SharedDataTestCase");

  testSuite->addTest (new TestCaller< SharedDataTestCase>("Arithmetic Assignment operators",&SharedDataTestCase::testEQ));
  testSuite->addTest (new TestCaller< SharedDataTestCase>("Copy Constructor",&SharedDataTestCase::testCC));
  testSuite->addTest (new TestCaller< SharedDataTestCase>("Assignment operator",&SharedDataTestCase::testAssign));
  testSuite->addTest (new TestCaller< SharedDataTestCase>("setToZero",&SharedDataTestCase::testSetToZero));
  testSuite->addTest (new TestCaller< SharedDataTestCase>("setTaggedValueFromCPP",&SharedDataTestCase::testSetTaggedValueFromCPP));
  testSuite->addTest (new TestCaller< SharedDataTestCase>("getDataAtOffset",&SharedDataTestCase::testGetDataAtOffset));
  testSuite->addTest (new TestCaller< SharedDataTestCase>("getDataPoint",&SharedDataTestCase::testGetDataPoint));
  testSuite->addTest (new TestCaller< SharedDataTestCase>("getSampleRW",&SharedDataTestCase::testGetSampleRW));
  return testSuite;
}