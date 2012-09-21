
/*****************************************************************************
*
* Copyright (c) 2003-2012 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development since 2012 by School of Earth Sciences
*
*****************************************************************************/

#ifdef BADPYTHONMACROS
// This hack is required for BSD/OSX builds with python 2.7
// (and possibly others).  It must be the first include.
// From bug reports online it seems that python redefines
// some c macros that are functions in c++.
// c++ doesn't like that!
#include <Python.h>
#undef BADPYTHONMACROS
#endif

#include "FunctionSpaceTestCase.h"

#include "escript/FunctionSpace.h"
#include "escript/NullDomain.h"

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

