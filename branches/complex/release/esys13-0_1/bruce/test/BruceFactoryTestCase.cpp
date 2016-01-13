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


#include "escript/AbstractContinuousDomain.h"

#include "bruce/BruceFactory.h"

#include "BruceFactoryTestCase.h"

#include <iostream>

using namespace CppUnitTest;

using namespace escript;
using namespace bruce;

using namespace std;

void BruceFactoryTestCase::setUp() {
  //
  // This is called before each test is run
 
}

void BruceFactoryTestCase::tearDown() {
  //
  // This is called after each test has been run
 
}

void BruceFactoryTestCase::testAll() {

  cout << endl;

  cout << "\tTest brick factory." << endl;

  AbstractContinuousDomain* b = brick();

  assert(b->getDim()==3);

  std::pair<int,int> b_shape = b->getDataShape(0);

  assert(b_shape.first==1);
  assert(b_shape.second==8);

  delete b;

  b = brick(11,11,11,10,10,10);

  assert(b->getDim()==3);

  b_shape = b->getDataShape(0);

  assert(b_shape.first==1);
  assert(b_shape.second==1331);

  delete b;


  cout << "\tTest rectangle factory." << endl;

  AbstractContinuousDomain* r = rectangle();

  assert(r->getDim()==2);

  std::pair<int,int> r_shape = r->getDataShape(0);
  
  assert(r_shape.first==1);
  assert(r_shape.second==4);

  delete r;

  r = rectangle(11,11,10,10);

  assert(r->getDim()==2);

  r_shape = r->getDataShape(0);
  
  assert(r_shape.first==1);
  assert(r_shape.second==121);

  delete r;

}

TestSuite* BruceFactoryTestCase::suite ()
{
  //
  // create the suite of tests to perform.
  TestSuite *testSuite = new TestSuite ("BruceFactoryTestCase");

  testSuite->addTest (new TestCaller< BruceFactoryTestCase>("testAll",&BruceFactoryTestCase::testAll));
  return testSuite;
}
