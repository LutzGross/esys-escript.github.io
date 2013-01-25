
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

#include <iostream>

#include "EscriptDatasetTestCase.h"
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
    runner.addTest(EscriptDatasetTestCase::suite());
    runner.run(controller);
    CompilerOutputter outputter( &result, std::cerr );
    outputter.write();
#ifdef ESYS_MPI
    MPI_Finalize();
#endif
    return result.wasSuccessful() ? 0 : 1;
}

