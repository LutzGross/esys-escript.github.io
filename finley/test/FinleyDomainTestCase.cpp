
/*****************************************************************************
*
* Copyright (c) 2003-2020 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/

#include <escript/AbstractContinuousDomain.h>

#include "FinleyDomainTestCase.h"

#include <finley/DomainFactory.h>

#include <cppunit/TestCaller.h>
#include <boost/scoped_ptr.hpp>

using namespace escript;
using namespace finley;
using namespace CppUnit;

void FinleyDomainTestCase::testAll()
{
    JMPI info = makeInfo(MPI_COMM_WORLD);
	Domain_ptr dom(brick(info));
    CPPUNIT_ASSERT(dom->getDim() == 3);
}

TestSuite* FinleyDomainTestCase::suite()
{
    TestSuite *testSuite = new TestSuite("FinleyDomainTestCase");

    testSuite->addTest(new TestCaller<FinleyDomainTestCase>(
                "testAll",&FinleyDomainTestCase::testAll));
    return testSuite;
}

