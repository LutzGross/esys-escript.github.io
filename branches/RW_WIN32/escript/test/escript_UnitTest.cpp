// escript_UnitTest.cpp : Defines the entry point for the console application.
//


#include <iostream>

// TODO : Can't we just use int main(int, char* []) on all platforms?
#ifdef MSVC
#include <tchar.h>
#else
#define _TCHAR char
#define _tmain main
#endif

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
#include "DataVariableTestCase.h"
#include "DataCTestCase.h"
#include "DataAlgorithmAdapterTestCase.h"
#include "FunctionSpaceTestCase.h"
#include "DataProfTestCase.h"
#include "DataTestCase.h"

#include "CppUnitTest/TestRunner.h"

using namespace CppUnitTest;

int _tmain(int argc, _TCHAR* argv[])
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
	// TODO: put this back RW runner.addTest ("DataArrayView", DataArrayViewTestCase::suite());
	runner.addTest ("DataBlocks2D", DataBlocks2DTestCase::suite());
	runner.addTest ("DataVector", DataVectorTestCase::suite());
	runner.addTest ("Taipan", TaipanTestCase::suite());
	runner.addTest ("DataVariable", DataVariableTestCase::suite());
	runner.addTest ("DataC", DataCTestCase::suite());
	runner.addTest ("DataAlgorithmAdapter", DataAlgorithmAdapterTestCase::suite());
	runner.addTest ("FunctionSpace", FunctionSpaceTestCase::suite());
	runner.addTest ("DataProf", DataProfTestCase::suite());
	runner.addTest ("Data", DataTestCase::suite());

	// actually run the unit tests.
	runner.run (argc, argv);
	return 0;
}

