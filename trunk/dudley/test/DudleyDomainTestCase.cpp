
/*****************************************************************************
*
* Copyright (c) 2003-2017 by The University of Queensland
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

#include "DudleyDomainTestCase.h"

#include <dudley/DomainFactory.h>

#include <cppunit/TestCaller.h>
#include <boost/scoped_ptr.hpp>

using namespace escript;
using namespace dudley;
using namespace CppUnit;

void DudleyDomainTestCase::testAll()
{
    JMPI info = makeInfo(MPI_COMM_WORLD);
	Domain_ptr dom(brick(info));
    CPPUNIT_ASSERT(dom->getDim() == 3);
}

TestSuite* DudleyDomainTestCase::suite()
{
    TestSuite *testSuite = new TestSuite("DudleyDomainTestCase");

    testSuite->addTest(new TestCaller<DudleyDomainTestCase>(
                "testAll", &DudleyDomainTestCase::testAll));
    return testSuite;
}

