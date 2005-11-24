// bruce_UnitTests.cpp : Defines the entry point for the console application.
//

#include <iostream>
#include <tchar.h>

#include "BruceTestCase.h"
#include "BruceFactoryTestCase.h"

#include "CppUnitTest/TestRunner.h"

using namespace CppUnitTest;

int _tmain(int argc, _TCHAR* argv[])
{
	//
	// object which runs all of the tests
	TestRunner runner;
	//
	// add the RangeTestCase suite of tests to the runner
	runner.addTest ("Bruce", BruceTestCase::suite());
	runner.addTest ("BruceFactory", BruceFactoryTestCase::suite());

	// actually run the unit tests.
	runner.run (argc, argv);
	return 0;
}

