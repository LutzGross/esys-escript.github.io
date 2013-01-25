
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

#include <iostream>

#include "DataEmptyTestCase.h"
#include "DataConstantTestCase.h"
#include "DataTaggedTestCase.h"
#include "DataExpandedTestCase.h"
#include "DataFactoryTestCase.h"
#include "DataArrayTestCase.h"
#include "DataArrayViewTestCase.h"
#include "DataBlocks2DTestCase.h"
#include "DataVectorTestCase.h"
#include "TaipanTestCase.h"
#include "DataCTestCase.h"
#include "DataAlgorithmAdapterTestCase.h"
#include "FunctionSpaceTestCase.h"
#include "DataTestCase.h"

#include "tools/CppUnitTest/TestRunner.h"

using namespace CppUnitTest;

extern "C"{
#include "paso/Paso_MPI.h"
}

int main(int argc, char* argv[])
{ 
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
	runner.addTest ("DataEmpty", DataEmptyTestCase::suite());
	runner.addTest ("DataConstant", DataConstantTestCase::suite());
	runner.addTest ("DataTagged", DataTaggedTestCase::suite());
	runner.addTest ("DataExpanded", DataExpandedTestCase::suite());
	runner.addTest ("DataFactory", DataFactoryTestCase::suite());
	runner.addTest ("DataArray", DataArrayTestCase::suite());
	runner.addTest ("DataArrayView", DataArrayViewTestCase::suite());
	runner.addTest ("DataBlocks2D", DataBlocks2DTestCase::suite());
	runner.addTest ("DataVector", DataVectorTestCase::suite());
	runner.addTest ("Taipan", TaipanTestCase::suite());
	runner.addTest ("DataC", DataCTestCase::suite());
	runner.addTest ("DataAlgorithmAdapter", DataAlgorithmAdapterTestCase::suite());
	runner.addTest ("FunctionSpace", FunctionSpaceTestCase::suite());
	runner.addTest ("Data", DataTestCase::suite());

	// actually run the unit tests.
	runner.run (argc, argv);

#ifdef PASO_MPI
        MPI_Finalize();
#endif

	return 0;
}


