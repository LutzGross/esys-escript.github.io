/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* https://github.com/LutzGross/esys-escript.github.io
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
*
*****************************************************************************/

#include <escript/AbstractContinuousDomain.h>

#include "OxleyDomainTestCase.h"

#include <cppunit/TestCaller.h>

#include <boost/scoped_ptr.hpp>

using namespace escript;
// using namespace oxley;
using namespace CppUnit;

void OxleyDomainTestCase::testAll()
{
 //    JMPI info = makeInfo(MPI_COMM_WORLD);
	// Domain_ptr dom(brick(info));
 //    CPPUNIT_ASSERT(dom->getDim() == 3);
}

TestSuite* OxleyDomainTestCase::suite()
{
    TestSuite *testSuite = new TestSuite("OxleyDomainTestCase");

    testSuite->addTest(new TestCaller<OxleyDomainTestCase>(
                "testAll",&OxleyDomainTestCase::testAll));
    return testSuite;
}

