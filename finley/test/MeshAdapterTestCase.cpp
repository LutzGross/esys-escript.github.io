/* 
 *****************************************************************************
 *                                                                           *
 *       COPYRIGHT  ACcESS  -  All Rights Reserved                           *
 *                                                                           *
 * This software is the property of ACcESS. No part of this code             *
 * may be copied in any form or by any means without the expressed written   *
 * consent of ACcESS.  Copying, use or modification of this software         *
 * by any unauthorised person is illegal unless that person has a software   *
 * license agreement with ACcESS.                                            *
 *                                                                           *
 *****************************************************************************
*/

#include "finleycpp/CppAdapter/MeshAdapter.h"
#include "finleycpp/CppAdapter/MeshAdapterFactory.h"

#include "escriptcpp/AbstractContinuousDomain.h"

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

