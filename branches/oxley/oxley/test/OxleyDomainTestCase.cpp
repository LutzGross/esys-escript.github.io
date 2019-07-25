/*****************************************************************************
*
* Copyright (c) 2003-2019 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/

#include <escript/AbstractContinuousDomain.h>

#include "OxleyDomainTestCase.h"

#include <cppunit/TestCaller.h>

#include <boost/scoped_ptr.hpp>

using namespace escript;
using namespace oxley;
using namespace CppUnit;

void OxleyDomainTestCase::testAll()
{
    JMPI info = makeInfo(MPI_COMM_WORLD);
	Domain_ptr dom(brick(info));
    CPPUNIT_ASSERT(dom->getDim() == 3);
}

TestSuite* OxleyDomainTestCase::suite()
{
    TestSuite *testSuite = new TestSuite("OxleyDomainTestCase");

    testSuite->addTest(new TestCaller<OxleyDomainTestCase>(
                "testAll",&OxleyDomainTestCase::testAll));
    return testSuite;
}

