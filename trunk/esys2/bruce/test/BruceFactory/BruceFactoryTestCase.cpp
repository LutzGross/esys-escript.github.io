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

#include "escript/Data/AbstractContinuousDomain.h"
#include "bruce/Bruce/BruceFactory.h"
#include "BruceFactoryTestCase.h"

using namespace CppUnitTest;

using namespace escript;
using namespace bruce;

void BruceFactoryTestCase::setUp() {
  //
  // This is called before each test is run
 
}

void BruceFactoryTestCase::tearDown() {
  //
  // This is called after each test has been run
 
}

void BruceFactoryTestCase::testAll() {
  //
  // The test code may be entered here
  // There is nothing special about the function name, it may be renamed to
  // something more suitable. 
  // As many test methods as desired may be added.

  AbstractContinuousDomain* brick;
  AbstractContinuousDomain* rectangle;
  
}

TestSuite* BruceFactoryTestCase::suite ()
{
  //
  // create the suite of tests to perform.
  TestSuite *testSuite = new TestSuite ("BruceFactoryTestCase");

  testSuite->addTest (new TestCaller< BruceFactoryTestCase>("testAll",&BruceFactoryTestCase::testAll));
  return testSuite;
}
