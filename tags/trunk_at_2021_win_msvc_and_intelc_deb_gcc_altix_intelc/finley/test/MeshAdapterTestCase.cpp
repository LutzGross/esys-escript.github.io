
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


#include "finley/CppAdapter/MeshAdapter.h"
#include "finley/CppAdapter/MeshAdapterFactory.h"

#include "escript/AbstractContinuousDomain.h"

#include "MeshAdapterTestCase.h"

#include <boost/scoped_ptr.hpp>

using namespace escript;
using namespace finley;
using namespace CppUnitTest;

void MeshAdapterTestCase::setUp() {
  //
  // This is called before each test is run
 
}

void MeshAdapterTestCase::tearDown() {
  //
  // This is called after each test has been run
 
}

void MeshAdapterTestCase::testAll() {
  //
  // test construction of a mesh using the brick factory method
//   boost::scoped_ptr<AbstractContinuousDomain> myMesh(brick());
	brick();		// brick now returns a Domain_ptr which will auto delete
}

TestSuite* MeshAdapterTestCase::suite ()
{
  //
  // create the suite of tests to perform.
  TestSuite *testSuite = new TestSuite ("MeshAdapterTestCase");

  testSuite->addTest (new TestCaller< MeshAdapterTestCase>("testAll",&MeshAdapterTestCase::testAll));
  return testSuite;
}

