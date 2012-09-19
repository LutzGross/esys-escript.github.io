
/*******************************************************
*
* Copyright (c) 2003-2012 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


#include <iostream>

#include "DataEmptyTestCase.h"
#include "DataConstantTestCase.h"
#include "DataTaggedTestCase.h"
#include "DataExpandedTestCase.h"
#include "DataFactoryTestCase.h"
#include "DataBlocks2DTestCase.h"
#include "DataVectorTestCase.h"
#include "TaipanTestCase.h"
#include "DataCTestCase.h"
#include "DataAlgorithmAdapterTestCase.h"
#include "FunctionSpaceTestCase.h"
#include "DataTestCase.h"
#include "DataMathsTestCase.h"
#include "DataTypesTestCase.h"
#include "DataLazyTestCase.h"
#include "SharedDataTestCase.h"

#include <cppunit/CompilerOutputter.h>
#include <cppunit/TestResult.h>
#include <cppunit/TestResultCollector.h>
#include <cppunit/TestRunner.h>

using namespace CppUnit;

extern "C"{
#include "esysUtils/Esys_MPI.h"
}

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
	runner.addTest(SharedDataTestCase::suite());
	runner.addTest(DataTypesTestCase::suite());
	runner.addTest(DataMathsTestCase::suite());
	runner.addTest(DataEmptyTestCase::suite());
	runner.addTest(DataConstantTestCase::suite());
 	runner.addTest(DataTaggedTestCase::suite());
	runner.addTest(DataExpandedTestCase::suite());
	runner.addTest(DataFactoryTestCase::suite());
	runner.addTest(DataBlocks2DTestCase::suite());
	runner.addTest(DataVectorTestCase::suite());
	runner.addTest(TaipanTestCase::suite());
	runner.addTest(DataCTestCase::suite());
 	runner.addTest(DataAlgorithmAdapterTestCase::suite());
	runner.addTest(FunctionSpaceTestCase::suite());
	runner.addTest(DataTestCase::suite());
	runner.addTest(DataLazyTestCase::suite());

	runner.run(controller);
    CompilerOutputter outputter( &result, std::cerr );
    outputter.write();
#ifdef ESYS_MPI
        MPI_Finalize();
#endif
    return result.wasSuccessful() ? 0 : 1;
}

