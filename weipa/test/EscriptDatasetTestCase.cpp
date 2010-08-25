
/*******************************************************
*
* Copyright (c) 2003-2010 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


#include "weipa/EscriptDataset.h"
#include "EscriptDatasetTestCase.h"

using namespace CppUnitTest;
using namespace weipa;
using namespace std;

void EscriptDatasetTestCase::setUp()
{
    // This is called before each test is run
}

void EscriptDatasetTestCase::tearDown()
{
    // This is called after each test has been run
}

void EscriptDatasetTestCase::testAll()
{
    cout << endl;
    cout << "\tTest default constructor." << endl;
    EscriptDataset dataset;

    cout << "\tTest saveSilo without data." << endl;
    assert(dataset.saveSilo("dummy") == false);

    cout << "\tTest saveVTK without data." << endl;
    assert(dataset.saveVTK("dummy") == false);

    cout << "\tTest getConvertedDomain without data." << endl;
    assert(dataset.getConvertedDomain().size() == 0);

    cout << "\tTest getVariables without data." << endl;
    assert(dataset.getVariables().size() == 0);

    cout << "\tTest getMeshVariables without data." << endl;
    assert(dataset.getMeshVariables().size() == 0);
}

TestSuite* EscriptDatasetTestCase::suite()
{
    //
    // create the suite of tests to perform.
    TestSuite *testSuite = new TestSuite("EscriptDatasetTestCase");

    testSuite->addTest(new TestCaller<EscriptDatasetTestCase>(
                "testAll",&EscriptDatasetTestCase::testAll));
    return testSuite;
}

