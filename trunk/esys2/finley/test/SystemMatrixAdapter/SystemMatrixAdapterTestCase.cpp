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
#include "finley/CPPAdapter/SystemMatrixAdapter.h"
#include "finley/CPPAdapter/FinleyAdapterException.h"
#include "finley/CPPAdapter/FinleyError.h"

#include "SystemMatrixAdapterTestCase.h"

using namespace std;

using namespace CppUnitTest;

using namespace escript;
using namespace finley;

static Finley_Mesh *mesh;
static Finley_SystemMatrix *system_matrix;

static Finley_SystemMatrixType type;

static int symmetric;

static int row_blocksize;
static int column_blocksize;

static int reduce_row_order;
static int reduce_col_order;

static FunctionSpace row_functionspace;
static FunctionSpace colum_functionspace;

void SystemMatrixAdapterTestCase::setUp() {
  //
  // This is called before each test is run

  mesh = Finley_Mesh_alloc("foo", 2, 1);

  type = UNKNOWN;

  symmetric = 0;

  row_blocksize = 10;
  column_blocksize = 10;

  reduce_row_order = 0;
  reduce_col_order = 0;

  system_matrix = Finley_SystemMatrix_alloc(mesh, type, symmetric, row_blocksize, reduce_row_order, column_blocksize, reduce_col_order);

}

void SystemMatrixAdapterTestCase::tearDown() {
  //
  // This is called after each test has been run

  Finley_SystemMatrix_dealloc(system_matrix);

}

void SystemMatrixAdapterTestCase::testAll() {
  //
  // The test code may be entered here
  // There is nothing special about the function name, it may be renamed to
  // something more suitable. 
  // As many test methods as desired may be added.

  cout << endl;

  try {
    cout << "\tTest illegal default construction." << endl;
    SystemMatrixAdapter system_matrix_adapter;
    assert(false);
  }
  catch (FinleyAdapterException& e) {
    cout << "\t" << e.toString() << endl;
    assert(true);
  }

  cout << "\tTest constructor." << endl;

  SystemMatrixAdapter system_matrix_adapter(system_matrix, row_blocksize, row_functionspace, column_blocksize, colum_functionspace);
  Finley_SystemMatrix* adapter_system_matrix = system_matrix_adapter.getFinley_SystemMatrix();
  assert(adapter_system_matrix == system_matrix);

  cout << "\tExercise nullifyRowsAndCols." << endl;

  Data row_q;
  Data col_q;
  double mdv = 1.0;

  system_matrix_adapter.nullifyRowsAndCols(row_q, col_q, mdv);

}

TestSuite* SystemMatrixAdapterTestCase::suite ()
{
  //
  // create the suite of tests to perform.
  TestSuite *testSuite = new TestSuite ("SystemMatrixAdapterTestCase");

  testSuite->addTest (new TestCaller< SystemMatrixAdapterTestCase>("testAll",&SystemMatrixAdapterTestCase::testAll));
  return testSuite;
}
