
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/

#include <escript/DataTypes.h>

#include "TaipanTestCase.h"

#include <escript/Taipan.h>

#include <cppunit/TestCaller.h>
#include <iostream>

using namespace std;
using namespace escript;
using namespace CppUnit;

void TaipanTestCase::testN1() {

  // Arrays with N=1 are handled differently by the Taipan memory manager, these
  // are never disposed of to maximise reusability, so they are tested separately.

  cout << endl;

  Taipan t;

  const int dim = 1;

  double* array;
  double* arraY[10];

  cout << "\tTest Taipan memTable with one array of size 1 under management." << endl;

  array = t.new_array(dim,1);

  CPPUNIT_ASSERT(array[0] == 0.0);

  CPPUNIT_ASSERT(t.num_arrays() == 1);
  CPPUNIT_ASSERT(t.num_arrays(1) == 1);
  CPPUNIT_ASSERT(t.num_free(1) == 0);
  CPPUNIT_ASSERT(t.num_elements() == 1);

  t.delete_array(array);

  CPPUNIT_ASSERT(t.num_arrays() == 1);
  CPPUNIT_ASSERT(t.num_elements() == 1);

  cout << "\tTest Taipan memTable with ten arrays of size 1 under management." << endl;

  // allocate all ten arrays
  // the number of arrays under management will increase with each allocation
//  arraY[10];
  for (int i=0; i<10; i++) {
    arraY[i] = t.new_array(dim,1);
    CPPUNIT_ASSERT(t.num_arrays() == i+1);
    CPPUNIT_ASSERT(t.num_arrays(1) == i+1);
    CPPUNIT_ASSERT(t.num_free(1) == 0);
    CPPUNIT_ASSERT(t.num_elements() == i+1);
  }

  for (int i=0; i<10; i++) {
    arraY[i][0] = i;
  }
  for (int i=0; i<10; i++) {
    CPPUNIT_ASSERT(arraY[i][0] == i);
  }

  // delete all but the last array under management
  // the number of arrays under management does not change
  // but the number of free arrays will increase with each deallocation
  for (int i=0; i<9; i++) {
    t.delete_array(arraY[i]);
    CPPUNIT_ASSERT(t.num_arrays() == 10);
    CPPUNIT_ASSERT(t.num_arrays(1) == 10);
    CPPUNIT_ASSERT(t.num_free(1) == i+1);
    CPPUNIT_ASSERT(t.num_elements() == 10);
  }
  CPPUNIT_ASSERT(arraY[9][0] == 9);

  // now reallocate the 9 free arrays
  // the preserved arrays will be reused
  // and the number of free arrays will decrease with each allocation
  for (int i=0; i<9; i++) {
    arraY[i] = t.new_array(dim,1);
    CPPUNIT_ASSERT(t.num_arrays() == 10);
    CPPUNIT_ASSERT(t.num_arrays(1) == 10);
    CPPUNIT_ASSERT(t.num_free(1) == 8-i);
    CPPUNIT_ASSERT(t.num_elements() == 10);
  }
  for (int i=0; i<10; i++) {
    CPPUNIT_ASSERT(arraY[i][0] == i);
  }

  // delete all but the last array under management
  // the number of arrays under management does not change
  // but the number of free arrays will increase with each deallocation
  for (int i=0; i<9; i++) {
    t.delete_array(arraY[i]);
    CPPUNIT_ASSERT(t.num_arrays() == 10);
    CPPUNIT_ASSERT(t.num_arrays(1) == 10);
    CPPUNIT_ASSERT(t.num_free(1) == i+1);
    CPPUNIT_ASSERT(t.num_elements() == 10);
  }
  CPPUNIT_ASSERT(arraY[9][0] == 9);

  // finally delete the last array
  t.delete_array(arraY[9]);
  CPPUNIT_ASSERT(t.num_arrays() == 10);
  CPPUNIT_ASSERT(t.num_elements() == 10);
}

void TaipanTestCase::testN0() {

  // Arrays with N=0 are handled differently by the Taipan memory manager, these
  // are never disposed of to maximise reusability, so they are tested separately.

  cout << endl;

  Taipan t;

  const int dim = 1;

  double* array;
  double* arraY[10];

  cout << "\tTest Taipan memTable with one array of size 0 under management." << endl;

  array = t.new_array(dim,0);

  CPPUNIT_ASSERT(t.num_arrays() == 1);
  CPPUNIT_ASSERT(t.num_arrays(0) == 1);
  CPPUNIT_ASSERT(t.num_free(0) == 0);
  CPPUNIT_ASSERT(t.num_elements() == 0);

  t.delete_array(array);

  CPPUNIT_ASSERT(t.num_arrays() == 1);
  CPPUNIT_ASSERT(t.num_elements() == 0);

  cout << "\tTest Taipan memTable with ten arrays of size 0 under management." << endl;

  // allocate all ten arrays
  // the number of arrays under management will increase with each allocation
  for (int i=0; i<10; i++) {
    arraY[i] = t.new_array(dim,0);
    CPPUNIT_ASSERT(t.num_arrays() == i+1);
    CPPUNIT_ASSERT(t.num_arrays(0) == i+1);
    CPPUNIT_ASSERT(t.num_free(0) == 0);
    CPPUNIT_ASSERT(t.num_elements() == 0);
  }

  // delete all but the last array under management
  // the number of arrays under management does not change
  // but the number of free arrays will increase with each deallocation
  for (int i=0; i<9; i++) {
    t.delete_array(arraY[i]);
    CPPUNIT_ASSERT(t.num_arrays() == 10);
    CPPUNIT_ASSERT(t.num_arrays(0) == 10);
    CPPUNIT_ASSERT(t.num_free(0) == i+1);
    CPPUNIT_ASSERT(t.num_elements() == 0);
  }

  // now reallocate the 9 free arrays
  // the preserved arrays will be reused
  // and the number of free arrays will decrease with each allocation
  for (int i=0; i<9; i++) {
    arraY[i] = t.new_array(dim,0);
    CPPUNIT_ASSERT(t.num_arrays() == 10);
    CPPUNIT_ASSERT(t.num_arrays(0) == 10);
    CPPUNIT_ASSERT(t.num_free(0) == 8-i);
    CPPUNIT_ASSERT(t.num_elements() == 0);
  }

  // delete all but the last array under management
  // the number of arrays under management does not change
  // but the number of free arrays will increase with each deallocation
  for (int i=0; i<9; i++) {
    t.delete_array(arraY[i]);
    CPPUNIT_ASSERT(t.num_arrays() == 10);
    CPPUNIT_ASSERT(t.num_arrays(0) == 10);
    CPPUNIT_ASSERT(t.num_free(0) == i+1);
    CPPUNIT_ASSERT(t.num_elements() == 0);
  }

  // finally delete the last array
  t.delete_array(arraY[9]);
  CPPUNIT_ASSERT(t.num_arrays() == 10);
  CPPUNIT_ASSERT(t.num_elements() == 0);
}

void TaipanTestCase::testAll() {

  cout << endl;

  Taipan t;

  double* array;
  double* arraY[10];

  const int dim = 1;

  cout << "\tTest empty Taipan memTable." << endl;

  CPPUNIT_ASSERT(t.num_arrays() == 0);
  CPPUNIT_ASSERT(t.num_elements() == 0);

  cout << "\tTest Taipan memTable with one array of size 10 under management." << endl;

  array = t.new_array(dim,10);

  for (int i=0; i<10; i++) {
    CPPUNIT_ASSERT(array[i] == 0.0);
  }

  for (int i=0; i<10; i++) {
    array[i] = i;
  }
  for (int i=0; i<10; i++) {
    CPPUNIT_ASSERT(array[i] == i);
  }

  CPPUNIT_ASSERT(t.num_arrays() == 1);
  CPPUNIT_ASSERT(t.num_arrays(10) == 1);
  CPPUNIT_ASSERT(t.num_free(10) == 0);
  CPPUNIT_ASSERT(t.num_elements() == 10);

  t.delete_array(array);

  CPPUNIT_ASSERT(t.num_arrays() == 0);
  CPPUNIT_ASSERT(t.num_elements() == 0);
 
  cout << "\tTest Taipan memTable with ten arrays of size 10 under management." << endl;

  // allocate all ten arrays
  // the number of arrays under management will increase with each allocation
//  arraY[10];
  for (int i=0; i<10; i++) {
    arraY[i] = t.new_array(dim,10);
    CPPUNIT_ASSERT(t.num_arrays() == i+1);
    CPPUNIT_ASSERT(t.num_arrays(10) == i+1);
    CPPUNIT_ASSERT(t.num_free(10) == 0);
    CPPUNIT_ASSERT(t.num_elements() == (i+1)*10);
  }

  int val = 0;
  for (int i=0; i<10; i++) {
    for (int j=0; j<10; j++) {
      arraY[i][j] = val;
      val++;
    }
  }
  val = 0;
  for (int i=0; i<10; i++) {
    for (int j=0; j<10; j++) {
      CPPUNIT_ASSERT(arraY[i][j] == val);
      val++;
    }
  }

  // delete all but the last array under management
  // the number of arrays under management does not change
  // but the number of free arrays will increase with each deallocation
  for (int i=0; i<9; i++) {
    t.delete_array(arraY[i]);
    CPPUNIT_ASSERT(t.num_arrays() == 10);
    CPPUNIT_ASSERT(t.num_arrays(10) == 10);
    CPPUNIT_ASSERT(t.num_free(10) == i+1);
    CPPUNIT_ASSERT(t.num_elements() == 100);
  }
  val = 90;
  for (int j=0; j<10; j++) {
    CPPUNIT_ASSERT(arraY[9][j] == val);
    val++;
  }

  // now reallocate the 9 free arrays
  // the preserved arrays will be reused
  // and the number of free arrays will decrease with each allocation
  for (int i=0; i<9; i++) {
    arraY[i] = t.new_array(dim,10);
    CPPUNIT_ASSERT(t.num_arrays() == 10);
    CPPUNIT_ASSERT(t.num_arrays(10) == 10);
    CPPUNIT_ASSERT(t.num_free(10) == 8-i);
    CPPUNIT_ASSERT(t.num_elements() == 100);
  }
  val = 0;
  for (int i=0; i<10; i++) {
    for (int j=0; j<10; j++) {
      CPPUNIT_ASSERT(arraY[i][j] == val);
      val++;
    }
  }

  // delete all but the last array under management
  // the number of arrays under management does not change
  // but the number of free arrays will increase with each deallocation
  for (int i=0; i<9; i++) {
    t.delete_array(arraY[i]);
    CPPUNIT_ASSERT(t.num_arrays() == 10);
    CPPUNIT_ASSERT(t.num_arrays(10) == 10);
    CPPUNIT_ASSERT(t.num_free(10) == i+1);
    CPPUNIT_ASSERT(t.num_elements() == 100);
  }
  val = 90;
  for (int j=0; j<10; j++) {
    CPPUNIT_ASSERT(arraY[9][j] == val);
    val++;
  }

  // finally delete the last array
  // all arrays under management are deleted
  t.delete_array(arraY[9]);
  CPPUNIT_ASSERT(t.num_arrays() == 0);
  CPPUNIT_ASSERT(t.num_elements() == 0);

  cout << "\tTest Taipan memTable with ten arrays of various sizes under management." << endl;

  // allocate all ten arrays
  // the number of arrays under management will increase with each allocation
  arraY[0] = t.new_array(dim,2);
  arraY[1] = t.new_array(dim,2);
  arraY[2] = t.new_array(dim,2);
  arraY[3] = t.new_array(dim,10);
  arraY[4] = t.new_array(dim,10);
  arraY[5] = t.new_array(dim,10);
  arraY[6] = t.new_array(dim,100);
  arraY[7] = t.new_array(dim,100);
  arraY[8] = t.new_array(dim,100);
  arraY[9] = t.new_array(dim,100);

  CPPUNIT_ASSERT(t.num_arrays() == 10);
  CPPUNIT_ASSERT(t.num_arrays(2) == 3);
  CPPUNIT_ASSERT(t.num_arrays(10) == 3);
  CPPUNIT_ASSERT(t.num_arrays(100) == 4);
  CPPUNIT_ASSERT(t.num_free(2) == 0);
  CPPUNIT_ASSERT(t.num_free(10) == 0);
  CPPUNIT_ASSERT(t.num_free(100) == 0);
  CPPUNIT_ASSERT(t.num_elements() == 436);

  // delete all but the first array of each size under management
  // the number of arrays under management does not change
  // but the number of free arrays will increase with each deallocation
  t.delete_array(arraY[1]);
  t.delete_array(arraY[2]);
  t.delete_array(arraY[4]);
  t.delete_array(arraY[5]);
  t.delete_array(arraY[7]);
  t.delete_array(arraY[8]);
  t.delete_array(arraY[9]);

  CPPUNIT_ASSERT(t.num_arrays() == 10);
  CPPUNIT_ASSERT(t.num_arrays(2) == 3);
  CPPUNIT_ASSERT(t.num_arrays(10) == 3);
  CPPUNIT_ASSERT(t.num_arrays(100) == 4);
  CPPUNIT_ASSERT(t.num_free(2) == 2);
  CPPUNIT_ASSERT(t.num_free(10) == 2);
  CPPUNIT_ASSERT(t.num_free(100) == 3);
  CPPUNIT_ASSERT(t.num_elements() == 436);

  // now reallocate the free arrays
  // the preserved arrays will be reused
  // and the number of free arrays will decrease with each allocation
  arraY[1] = t.new_array(dim,2);
  arraY[2] = t.new_array(dim,2);
  arraY[4] = t.new_array(dim,10);
  arraY[5] = t.new_array(dim,10);
  arraY[7] = t.new_array(dim,100);
  arraY[8] = t.new_array(dim,100);
  arraY[9] = t.new_array(dim,100);

  CPPUNIT_ASSERT(t.num_arrays() == 10);
  CPPUNIT_ASSERT(t.num_arrays(2) == 3);
  CPPUNIT_ASSERT(t.num_arrays(10) == 3);
  CPPUNIT_ASSERT(t.num_arrays(100) == 4);
  CPPUNIT_ASSERT(t.num_free(2) == 0);
  CPPUNIT_ASSERT(t.num_free(10) == 0);
  CPPUNIT_ASSERT(t.num_free(100) == 0);
  CPPUNIT_ASSERT(t.num_elements() == 436);

  // delete all but the last array of each size under management
  // the number of arrays under management does not change
  // but the number of free arrays will increase with each deallocation
  t.delete_array(arraY[1]);
  t.delete_array(arraY[2]);
  t.delete_array(arraY[4]);
  t.delete_array(arraY[5]);
  t.delete_array(arraY[7]);
  t.delete_array(arraY[8]);
  t.delete_array(arraY[9]);

  CPPUNIT_ASSERT(t.num_arrays() == 10);
  CPPUNIT_ASSERT(t.num_arrays(2) == 3);
  CPPUNIT_ASSERT(t.num_arrays(10) == 3);
  CPPUNIT_ASSERT(t.num_arrays(100) == 4);
  CPPUNIT_ASSERT(t.num_free(2) == 2);
  CPPUNIT_ASSERT(t.num_free(10) == 2);
  CPPUNIT_ASSERT(t.num_free(100) == 3);
  CPPUNIT_ASSERT(t.num_elements() == 436);

  // deallocate the last of the arrays
  t.delete_array(arraY[0]);
  t.delete_array(arraY[3]);
  t.delete_array(arraY[6]);

  CPPUNIT_ASSERT(t.num_arrays() == 0);
  CPPUNIT_ASSERT(t.num_elements() == 0);
}

TestSuite* TaipanTestCase::suite()
{
  //
  // create the suite of tests to perform.
  TestSuite *testSuite = new TestSuite("TaipanTestCase");

  testSuite->addTest(new TestCaller<TaipanTestCase>(
              "testAll",&TaipanTestCase::testAll));
  testSuite->addTest(new TestCaller<TaipanTestCase>(
              "testN1",&TaipanTestCase::testN1));
  testSuite->addTest(new TestCaller<TaipanTestCase>(
              "testN0",&TaipanTestCase::testN0));
  return testSuite;
}

