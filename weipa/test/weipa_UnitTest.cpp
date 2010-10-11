
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
#include "tools/CppUnitTest/TestRunner.h"

using namespace CppUnitTest;

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
    // object which runs all of the tests
    TestRunner runner;
    runner.addTest("EscriptDataset", EscriptDatasetTestCase::suite());

    // actually run the unit tests.
    runner.run(argc, argv);

#ifdef ESYS_MPI
    MPI_Finalize();
#endif

    return 0;
}

