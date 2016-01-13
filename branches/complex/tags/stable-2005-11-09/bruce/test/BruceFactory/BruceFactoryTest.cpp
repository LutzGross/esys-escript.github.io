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

#include "BruceFactoryTestCase.h"
#include "CppUnitTest/TestRunner.h"

using namespace CppUnitTest;

int main(int argc, char **argv) {
  //
  // object which runs all of the tests
  TestRunner runner;
  //
  // add the RangeTestCase suite of tests to the runner
  runner.addTest ("BruceFactory", BruceFactoryTestCase::suite());

  // actually run the unit tests.
  runner.run (argc, argv);
  return 0;
}
