
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


#include "escript/DataFactory.h"
#include "finley/CppAdapter/MeshAdapterFactory.h"
#include "weipa/EscriptDataset.h"
#include "EscriptDatasetTestCase.h"

using namespace CppUnitTest;
using namespace escript;
using namespace finley;
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
    EscriptDataset_ptr dataset(new EscriptDataset());

    cout << "\tTest saveSilo without data." << endl;
    assert(dataset->saveSilo("dummy") == false);

    cout << "\tTest saveVTK without data." << endl;
    assert(dataset->saveVTK("dummy") == false);

    cout << "\tTest getConvertedDomain without data." << endl;
    assert(dataset->getConvertedDomain().size() == 0);

    cout << "\tTest getVariables without data." << endl;
    assert(dataset->getVariables().size() == 0);

    cout << "\tTest getMeshVariables without data." << endl;
    assert(dataset->getMeshVariables().size() == 0);

    // instantiate a domain and data
    Domain_ptr dom2d(rectangle());
    //Domain_ptr dom3d(brick());
    escript::Data dataOn2d = Scalar(0.0, continuousFunction(*dom2d), true);
    //escript::Data dataOn3d = Scalar(0.0, continuousFunction(*dom3d), true);

    cout << "\tTest addData with NULL domain." << endl;
    assert(dataset->addData(dataOn2d, "foo", "bar") == false);

    cout << "\tTest setDomain." << endl;
    assert(dataset->setDomain(dom2d.get()) == true);

    cout << "\tTest bogus setDomain call." << endl;
    assert(dataset->setDomain(dom2d.get()) == false);

    cout << "\tTest getConvertedDomain." << endl;
    assert(dataset->getConvertedDomain().size() > 0);

    cout << "\tTest addData with valid data." << endl;
    assert(dataset->addData(dataOn2d, "test2d", "") == true);
    assert(dataset->getVariables().size() == 1);
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

