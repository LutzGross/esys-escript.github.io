
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

#define ESNEEDPYTHON
#include "esysUtils/first.h"


#include "DataVectorTestCase.h"

#include "escript/DataVector.h"
#include "escript/EsysException.h"

#include <cppunit/TestCaller.h>
#include <iostream>

using namespace std;
using namespace CppUnit;
using namespace escript;
using namespace escript::DataTypes;


void DataVectorTestCase::testAll()
{
  cout << endl;

  {
    cout << "\tCreate and check an empty DataVector object." << endl;

    RealVectorType vec;
    CPPUNIT_ASSERT(vec.size() == 0);
  }
 
  {
    cout << "\tCheck DataVector resize operation." << endl;

    RealVectorType vec;
    CPPUNIT_ASSERT(vec.size() == 0);

    vec.resize(1,0,1);
    CPPUNIT_ASSERT(vec.size() == 1);

    vec.resize(1000,0,1);
    CPPUNIT_ASSERT(vec.size() == 1000);

    vec.resize(0,0,1);
    CPPUNIT_ASSERT(vec.size() == 0);
  }

  {
    cout << "\tCreate and check DataVector objects of various sizes." << endl;

    RealVectorType vec1(0,0,1);
    CPPUNIT_ASSERT(vec1.size() == 0);

    RealVectorType vec2(1,0,1);
    CPPUNIT_ASSERT(vec2.size() == 1);

    RealVectorType vec3(1000,0,1);
    CPPUNIT_ASSERT(vec3.size() == 1000);
  }

  {
    cout << "\tAssign and check various elements to a DataVector." << endl;

    RealVectorType vec(1000,0,1);

    for (int i=0; i < 1000; i++) {
      vec[i] = i;
    }

    for (int i=0; i < 1000; i++) {
      CPPUNIT_ASSERT(vec[i] == i);
    }

    for (int i=0; i < 1000; i++) {
      vec[i] = i/1000;
    }

    for (int i=0; i < 1000; i++) {
      CPPUNIT_ASSERT(vec[i] == i/1000);
    }
  }

  {
    cout << "\tCheck DataVector copy constructor." << endl;

    RealVectorType vec1(1000,0,1);

    for (int i=0; i < 1000; i++) {
      vec1[i] = i;
    }

    RealVectorType vec2(vec1);

    CPPUNIT_ASSERT(vec1.size() == vec2.size());

    for (int i=0; i < 1000; i++) {
      CPPUNIT_ASSERT(vec2[i] == i);
    }
  }
 
  {
    cout << "\tCheck DataVector = operator." << endl;

    RealVectorType vec1(1000,0,1);

    for (int i=0; i < 1000; i++) {
      vec1[i] = i;
    }

    RealVectorType vec2;

    vec2 = vec1;

    CPPUNIT_ASSERT(vec1.size() == vec2.size());

    for (int i=0; i < 1000; i++) {
      CPPUNIT_ASSERT(vec2[i] == i);
    }
  }
 
  {
    cout << "\tCheck DataVector == operator." << endl;

    RealVectorType vec1(1000,0,1);

    for (int i=0; i < 1000; i++) {
      vec1[i] = i;
    }

    RealVectorType vec2;

    vec2 = vec1;

    CPPUNIT_ASSERT(vec1 == vec2);
  }
 
  {
    cout << "\tCheck DataVector != operator." << endl;

    RealVectorType vec1(1000,0,1);

    for (int i=0; i < 1000; i++) {
      vec1[i] = i;
    }

    RealVectorType vec2;

    CPPUNIT_ASSERT(vec1 != vec2);
  }
  #if defined DOASSERT
  {
    cout << "\tCheck DataVector index exception." << endl;

    RealVectorType vec(1000,0,1);

    CPPUNIT_ASSERT_THROW( (void) vec[1001],  EsysException);
  }
  #endif

}

TestSuite* DataVectorTestCase::suite()
{
  // create the suite of tests to perform.
  TestSuite *testSuite = new TestSuite("DataVectorTestCase");

  testSuite->addTest(new TestCaller<DataVectorTestCase>(
              "testAll",&DataVectorTestCase::testAll));
  return testSuite;
}

