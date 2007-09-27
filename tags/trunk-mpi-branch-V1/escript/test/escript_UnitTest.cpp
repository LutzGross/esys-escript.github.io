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
#include "DataProfTestCase.h"
#include "DataTestCase.h"

#include "tools/CppUnitTest/TestRunner.h"

using namespace CppUnitTest;

int main(int argc, char* argv[])
{ 
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
	runner.addTest ("DataProf", DataProfTestCase::suite());
	runner.addTest ("Data", DataTestCase::suite());

	// actually run the unit tests.
	runner.run (argc, argv);
	return 0;
}


