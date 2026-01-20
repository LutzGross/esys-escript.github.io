
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/

#include <escript/EsysMPI.h>

#include "SystemMatrixTestCase.h"

#include <cppunit/CompilerOutputter.h>
#include <cppunit/TestResult.h>
#include <cppunit/TestResultCollector.h>
#include <cppunit/TestRunner.h>

#include <iostream>

using namespace CppUnit;


int main(int argc, char* argv[])
{
    int mpiRank = 0;
    int mpiSize = 1;
#ifdef ESYS_MPI
    int status = MPI_Init(&argc, &argv);
    if (status != MPI_SUCCESS) {
        std::cerr << argv[0] << ": MPI_Init failed, exiting." << std::endl;
        return status;
    }
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
#endif
    TestResult controller;
    TestResultCollector result;
    controller.addListener(&result);
    TestRunner runner;
    if (mpiSize == 1) {
        runner.addTest(SystemMatrixTestCase::suite());
    } else {
        if (mpiRank == 0)
            std::cout << "Skipping SystemMatrixTestCase with more than one rank."
                      << std::endl;
    }
    runner.run(controller);
    CompilerOutputter outputter( &result, std::cerr );
    if (mpiRank == 0)
        outputter.write();
#ifdef ESYS_MPI
    MPI_Finalize();
#endif
    return result.wasSuccessful() ? 0 : 1;
}

