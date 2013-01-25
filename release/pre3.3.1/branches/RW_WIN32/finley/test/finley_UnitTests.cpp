// finley_UnitTests.cpp : Defines the entry point for the console application.
//

#include <iostream>

// TODO : Can't we just use int main(int, char* []) on all platforms?
#ifdef MSVC
#include <tchar.h>
#else
#define _TCHAR char
#define _tmain main
#endif

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

