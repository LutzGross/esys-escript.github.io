
/*****************************************************************************
*
* Copyright (c) 2003-2016 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/

#include <escript/FunctionSpace.h>

#include "FunctionSpaceTestCase.h"

#include <escript/NullDomain.h>

#include <cppunit/TestCaller.h>
#include <iostream>
#include <sstream>

using namespace CppUnit;
using namespace escript;
using namespace std;


void FunctionSpaceTestCase::testAll()
{
  cout << endl;
  cout << "\tTest default FunctionSpace constructor." << endl;

  // Test Default constructor
  FunctionSpace testFunctionSpace0;
  FunctionSpace testFunctionSpace1;

  CPPUNIT_ASSERT(testFunctionSpace0==testFunctionSpace1);

  cout << "\tTest FunctionSpace constructor." << endl;

  // Test constructor
  NullDomain* nullDomain=new NullDomain();	// the shared ptr will deal with it
  int testfunctionSpaceType = nullDomain->getFunctionCode();
  Domain_ptr nulldom(nullDomain);

  FunctionSpace testFunctionSpace2(nulldom, testfunctionSpaceType);
  
  CPPUNIT_ASSERT(testFunctionSpace1.getTypeCode()==testfunctionSpaceType);
  CPPUNIT_ASSERT(*(testFunctionSpace2.getDomain())==*nullDomain);
  CPPUNIT_ASSERT(testFunctionSpace1.getDim()==1);
  CPPUNIT_ASSERT(testFunctionSpace1==testFunctionSpace1);
  CPPUNIT_ASSERT(!(testFunctionSpace1!=testFunctionSpace1));

  FunctionSpace testFunctionSpace3=testFunctionSpace2;	// test copy constructor
cout << "Testing equality\n";
  CPPUNIT_ASSERT(testFunctionSpace3==testFunctionSpace2);
}

TestSuite* FunctionSpaceTestCase::suite()
{
  TestSuite *testSuite = new TestSuite("FunctionSpaceTestCase");

  testSuite->addTest(new TestCaller<FunctionSpaceTestCase>(
              "testAll",&FunctionSpaceTestCase::testAll));
  return testSuite;
}

