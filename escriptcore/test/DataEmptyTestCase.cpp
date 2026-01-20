
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

#include <escript/DataEmpty.h>

#include "DataEmptyTestCase.h"

#include <escript/FunctionSpace.h>
#include <escript/EsysException.h>

#include <cppunit/TestCaller.h>

using namespace CppUnit;
using namespace escript;
using namespace std;

void DataEmptyTestCase::testAll()
{
  cout << endl;
  cout << "\tTest default constructor." << endl;
  DataEmpty testData;

  cout << "\tTest toString method." << endl;
  CPPUNIT_ASSERT(testData.toString() == "(Empty Data)");

  cout << "\tTest getPointOffset." << endl;
  CPPUNIT_ASSERT_THROW(testData.getPointOffset(0,0), EsysException);
  
  cout << "\tTest getDataPoint." << endl;
  CPPUNIT_ASSERT_THROW(testData.getPointOffset(0,0), EsysException);

  cout << "\tTest getLength." << endl;
  CPPUNIT_ASSERT(testData.getLength() == 0);

  DataTypes::RegionType region;

  cout << "\tTest getSlice." << endl;
  CPPUNIT_ASSERT_THROW(testData.getSlice(region), EsysException);

  cout << "\tTest setSlice." << endl;
  CPPUNIT_ASSERT_THROW(testData.setSlice(0,region), EsysException);
}

TestSuite* DataEmptyTestCase::suite()
{
  TestSuite *testSuite = new TestSuite("DataEmptyTestCase");

  testSuite->addTest(new TestCaller<DataEmptyTestCase>(
              "testAll",&DataEmptyTestCase::testAll));
  return testSuite;
}

