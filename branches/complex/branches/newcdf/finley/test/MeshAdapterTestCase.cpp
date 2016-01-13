
/*****************************************************************************
*
* Copyright (c) 2003-2014 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/


#include "MeshAdapterTestCase.h"

#include "finley/CppAdapter/MeshAdapter.h"
#include "finley/CppAdapter/MeshAdapterFactory.h"

#include "escript/AbstractContinuousDomain.h"

#include <cppunit/TestCaller.h>
#include <boost/scoped_ptr.hpp>

using namespace escript;
using namespace finley;
using namespace CppUnit;

void MeshAdapterTestCase::testAll()
{
    // test construction of a mesh using the brick factory method
    // boost::scoped_ptr<AbstractContinuousDomain> myMesh(brick());
	esysUtils::JMPI info=esysUtils::makeInfo(MPI_COMM_WORLD);
	brick(info);	// brick now returns a Domain_ptr which will auto delete
}

TestSuite* MeshAdapterTestCase::suite()
{
    TestSuite *testSuite = new TestSuite("MeshAdapterTestCase");

    testSuite->addTest(new TestCaller<MeshAdapterTestCase>(
                "testAll",&MeshAdapterTestCase::testAll));
    return testSuite;
}

