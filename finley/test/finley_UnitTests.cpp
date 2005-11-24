// finley_UnitTests.cpp : Defines the entry point for the console application.
//

#include <iostream>
#include <tchar.h>

#include "MeshAdapterTestCase.h"

#include "CppUnitTest/TestRunner.h"

using namespace CppUnitTest;

int _tmain(int argc, _TCHAR* argv[])
{
	//
	// object which runs all of the tests
	TestRunner runner;
	//
	// add the RangeTestCase suite of tests to the runner
	runner.addTest ("MeshAdapter", MeshAdapterTestCase::suite());
   
	// actually run the unit tests.
	runner.run (argc, argv);

	return 0;
}

