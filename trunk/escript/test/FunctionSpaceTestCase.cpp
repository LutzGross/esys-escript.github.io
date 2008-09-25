
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


#include "escript/FunctionSpace.h"
#include "escript/NullDomain.h"

#include "FunctionSpaceTestCase.h"

#include <iostream>
#include <sstream>

using namespace CppUnitTest;
using namespace escript;
using namespace std;

void FunctionSpaceTestCase::setUp() {
  //
  // This is called before each test is run
 
}

void FunctionSpaceTestCase::tearDown() {
  //
  // This is called after each test has been run
 
}

void FunctionSpaceTestCase::testAll() {
  //
  // The test code may be entered here
  // There is nothing special about the function name, it may be renamed to
  // something more suitable. 
  // As many test methods as desired may be added.

  cout << endl;

  cout << "\tTest default FunctionSpace constructor." << endl;

  // Test Default constructor
  FunctionSpace testFunctionSpace0;
  FunctionSpace testFunctionSpace1;

  assert(testFunctionSpace0==testFunctionSpace1);

  cout << "\tTest FunctionSpace constructor." << endl;

  // Test constructor
  NullDomain nullDomain;
  int testfunctionSpaceType = nullDomain.getFunctionCode();

  FunctionSpace testFunctionSpace2(nullDomain, testfunctionSpaceType);
  
  assert(testFunctionSpace1.getTypeCode()==testfunctionSpaceType);
  assert(testFunctionSpace2.getDomain()==nullDomain);
  assert(testFunctionSpace1.getDim()==1);
  assert(testFunctionSpace1==testFunctionSpace1);
  assert(!(testFunctionSpace1!=testFunctionSpace1));

  testFunctionSpace1=testFunctionSpace2;
  assert(testFunctionSpace1==testFunctionSpace2);

}

TestSuite* FunctionSpaceTestCase::suite ()
{
  //
  // create the suite of tests to perform.
  TestSuite *testSuite = new TestSuite ("FunctionSpaceTestCase");

  testSuite->addTest (new TestCaller< FunctionSpaceTestCase>("testAll",&FunctionSpaceTestCase::testAll));
  return testSuite;
}
