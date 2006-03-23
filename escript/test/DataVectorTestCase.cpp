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
#include "escriptcpp/DataVector.h"
#include "esysUtils/EsysException.h"

#include "DataVectorTestCase.h"

#include <iostream>

using namespace std;
using namespace CppUnitTest;
using namespace escript;
using namespace esysUtils;

void DataVectorTestCase::setUp() {
  //
  // This is called before each test is run
 
}

void DataVectorTestCase::tearDown() {
  //
  // This is called after each test has been run
 
}

void DataVectorTestCase::testAll() {

  cout << endl;

  {
    cout << "\tCreate and check an empty DataVector object." << endl;

    DataVector vec;
    assert(vec.size() == 0);
  }
 
  {
    cout << "\tCheck DataVector resize operation." << endl;

    DataVector vec;
    assert(vec.size() == 0);

    vec.resize(1,0,1);
    assert(vec.size() == 1);

    vec.resize(1000,0,1);
    assert(vec.size() == 1000);

    vec.resize(0,0,1);
    assert(vec.size() == 0);
  }

  {
    cout << "\tCreate and check DataVector objects of various sizes." << endl;

    DataVector vec1(0,0,1);
    assert(vec1.size() == 0);

    DataVector vec2(1,0,1);
    assert(vec2.size() == 1);

    DataVector vec3(1000,0,1);
    assert(vec3.size() == 1000);
  }

  {
    cout << "\tAssign and check various elements to a DataVector." << endl;

    DataVector vec(1000,0,1);

    for (int i=0; i < 1000; i++) {
      vec[i] = i;
    }

    for (int i=0; i < 1000; i++) {
      assert(vec[i] == i);
    }

    for (int i=0; i < 1000; i++) {
      vec[i] = i/1000;
    }

    for (int i=0; i < 1000; i++) {
      assert(vec[i] == i/1000);
    }
  }

  {
    cout << "\tCheck DataVector copy constructor." << endl;

    DataVector vec1(1000,0,1);

    for (int i=0; i < 1000; i++) {
      vec1[i] = i;
    }

    DataVector vec2(vec1);

    assert(vec1.size() == vec2.size());

    for (int i=0; i < 1000; i++) {
      assert(vec2[i] == i);
    }
  }
 
  {
    cout << "\tCheck DataVector = operator." << endl;

    DataVector vec1(1000,0,1);

    for (int i=0; i < 1000; i++) {
      vec1[i] = i;
    }

    DataVector vec2;

    vec2 = vec1;

    assert(vec1.size() == vec2.size());

    for (int i=0; i < 1000; i++) {
      assert(vec2[i] == i);
    }
  }
 
  {
    cout << "\tCheck DataVector == operator." << endl;

    DataVector vec1(1000,0,1);

    for (int i=0; i < 1000; i++) {
      vec1[i] = i;
    }

    DataVector vec2;

    vec2 = vec1;

    assert(vec1 == vec2);
  }
 
  {
    cout << "\tCheck DataVector != operator." << endl;

    DataVector vec1(1000,0,1);

    for (int i=0; i < 1000; i++) {
      vec1[i] = i;
    }

    DataVector vec2;

    assert(vec1 != vec2);
  }
 
  {
    cout << "\tCheck DataVector index exception." << endl;

    DataVector vec(1000,0,1);

    try {
      double x = vec[1001];
      assert(false);
    }

    catch (EsysException& e) {
      //cout << e.toString() << endl;
      assert(true);
    }

  }

}

TestSuite* DataVectorTestCase::suite ()
{
  //
  // create the suite of tests to perform.
  TestSuite *testSuite = new TestSuite ("DataVectorTestCase");

  testSuite->addTest (new TestCaller< DataVectorTestCase>("testAll",&DataVectorTestCase::testAll));
  return testSuite;
}

