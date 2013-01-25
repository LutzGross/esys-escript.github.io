
/* $Id$ */

/*******************************************************
 *
 *           Copyright 2003-2007 by ACceSS MNRF
 *       Copyright 2007 by University of Queensland
 *
 *                http://esscc.uq.edu.au
 *        Primary Business: Queensland, Australia
 *  Licensed under the Open Software License version 3.0
 *     http://www.opensource.org/licenses/osl-3.0.php
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
  boost::scoped_ptr<AbstractContinuousDomain> myMesh(brick());
}

TestSuite* MeshAdapterTestCase::suite ()
{
  //
  // create the suite of tests to perform.
  TestSuite *testSuite = new TestSuite ("MeshAdapterTestCase");

  testSuite->addTest (new TestCaller< MeshAdapterTestCase>("testAll",&MeshAdapterTestCase::testAll));
  return testSuite;
}

