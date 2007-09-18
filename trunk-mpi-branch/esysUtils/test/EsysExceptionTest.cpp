
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

#include "EsysExceptionTestCase.h"
#include "tools/CppUnitTest/TestRunner.h"

using namespace CppUnitTest;

extern "C"{
#include "paso/Paso_MPI.h"
}

int main(int argc, char **argv) {

#ifdef PASO_MPI
  int status = MPI_Init(&argc, &argv);
  if (status != MPI_SUCCESS) {
    std::cerr << argv[0] << ": MPI_Init failed, exiting." << std::endl;
    return status;
  }
#endif

  //
  // object which runs all of the tests
  TestRunner runner;
  //
  // add the RangeTestCase suite of tests to the runner
  runner.addTest ("EsysException", EsysExceptionTestCase::suite());

  // actually run the unit tests.
  runner.run (argc, argv);

#ifdef PASO_MPI
  MPI_Finalize();
#endif

  return 0;
}
