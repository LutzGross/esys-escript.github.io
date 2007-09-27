#include "BruceTestCase.h"
#include "BruceFactoryTestCase.h"

#include "tools/CppUnitTest/TestRunner.h"

using namespace CppUnitTest;

int main(int argc, char* argv[])
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

