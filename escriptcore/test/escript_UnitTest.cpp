
/*****************************************************************************
*
* Copyright (c) 2003-2016 by The University of Queensland
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

#include <escript/EsysMPI.h>

#include "DataAlgorithmAdapterTestCase.h"
#include "DataConstantTestCase.h"
#include "DataEmptyTestCase.h"
#include "DataExpandedTestCase.h"
#include "DataFactoryTestCase.h"
#include "DataLazyTestCase.h"
#include "DataMathsTestCase.h"
#include "DataTaggedTestCase.h"
#include "DataTestCase.h"
#include "DataTypesTestCase.h"
#include "DataVectorTestCase.h"
#include "EsysExceptionTestCase.h"
#include "FileWriterTestCase.h"
#include "FunctionSpaceTestCase.h"
#include "SharedDataTestCase.h"
#include "TaipanTestCase.h"
#include "DataCombinationsTestCase.h"

#include <iostream>

#include <cppunit/CompilerOutputter.h>
#include <cppunit/TestResult.h>
#include <cppunit/TestResultCollector.h>
#include <cppunit/TestRunner.h>

using namespace CppUnit;

int main(int argc, char* argv[])
{
#ifdef ESYS_MPI
    int status = MPI_Init(&argc, &argv);
    if (status != MPI_SUCCESS) {
        std::cerr << argv[0] << ": MPI_Init failed, exiting." << std::endl;
        return status;
    }
#endif
    TestResult controller;
    TestResultCollector result;
    controller.addListener(&result);
	TestRunner runner;
    runner.addTest(EsysExceptionTestCase::suite());
	runner.addTest(SharedDataTestCase::suite());
	runner.addTest(DataTypesTestCase::suite());
	runner.addTest(DataFactoryTestCase::suite());
	runner.addTest(DataEmptyTestCase::suite());
	runner.addTest(DataConstantTestCase::suite());
 	runner.addTest(DataTaggedTestCase::suite());
	runner.addTest(DataExpandedTestCase::suite());
	runner.addTest(DataVectorTestCase::suite());
	runner.addTest(DataMathsTestCase::suite());
 	runner.addTest(DataAlgorithmAdapterTestCase::suite());
	runner.addTest(FileWriterTestCase::suite());
	runner.addTest(TaipanTestCase::suite());
	runner.addTest(FunctionSpaceTestCase::suite());
	runner.addTest(DataTestCase::suite());
	runner.addTest(DataLazyTestCase::suite());
	runner.addTest(DataCombinationsTestCase::suite());

	runner.run(controller);
    CompilerOutputter outputter( &result, std::cerr );
    outputter.write();
#ifdef ESYS_MPI
        MPI_Finalize();
#endif
    return result.wasSuccessful() ? 0 : 1;
}

