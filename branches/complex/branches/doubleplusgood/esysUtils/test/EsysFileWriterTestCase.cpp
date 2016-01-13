
/*****************************************************************************
*
* Copyright (c) 2003-2013 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development since 2012 by School of Earth Sciences
*
*****************************************************************************/


#include "EsysFileWriterTestCase.h"
#include "esysUtils/esysFileWriter.h"
#include <cppunit/TestCaller.h>
#include <fstream>
#include <sstream>

extern "C"{
#include "esysUtils/Esys_MPI.h"
}

using namespace CppUnit;
using namespace std;

using esysUtils::FileWriter;

void EsysFileWriterTestCase::testAll()
{
    const string filename("fwtest_file");
    int mpisize=1, mpirank=0;
    FileWriter* fw;
#ifdef ESYS_MPI
    MPI_Comm_size(MPI_COMM_WORLD, &mpisize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
    fw = new FileWriter(MPI_COMM_WORLD);
#else
    fw = new FileWriter();
#endif

    const char crank = (const char)mpirank;
    const char data[] = {crank,crank,crank,crank};
    ostringstream oss;
    oss.write(data, 4);

    cout << endl;
    cout << "\tTest open file." << endl;
    CPPUNIT_ASSERT(fw->openFile(filename) == true);
    cout << "\tTest writeOrdered." << endl;
    CPPUNIT_ASSERT(fw->writeOrdered(oss) == true);
    CPPUNIT_ASSERT(oss.str().length() == 0);
    fw->close();
    CPPUNIT_ASSERT(fileSize(filename) == 4*mpisize);

    cout << "\tTest open file with initial size." << endl;
    CPPUNIT_ASSERT(fw->openFile(filename, 100) == true);
    fw->close();
    CPPUNIT_ASSERT(fileSize(filename) == 100);

    CPPUNIT_ASSERT(fw->openFile(filename) == true);
    oss.write(data, 4);
    cout << "\tTest writeShared." << endl;
    CPPUNIT_ASSERT(fw->writeShared(oss) == true);
    CPPUNIT_ASSERT(oss.str().length() == 0);
    fw->close();
    CPPUNIT_ASSERT(fileSize(filename) == 4*mpisize);

    CPPUNIT_ASSERT(fw->openFile(filename) == true);
    oss.write(data, 4);
    cout << "\tTest writeAt." << endl;
    CPPUNIT_ASSERT(fw->writeAt(oss, 16*(mpirank+1)) == true);
    CPPUNIT_ASSERT(oss.str().length() == 0);
    fw->close();
    CPPUNIT_ASSERT(fileSize(filename) == 16*mpisize+4);
}

long EsysFileWriterTestCase::fileSize(string filename)
{
    ifstream f(filename.c_str());
    f.seekg(0, f.end);
    long pos = f.tellg();
    f.close();
#ifdef ESYS_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    return pos;
}

TestSuite* EsysFileWriterTestCase::suite()
{
    TestSuite *testSuite = new TestSuite("EsysFileWriterTestCase");
    testSuite->addTest(new TestCaller<EsysFileWriterTestCase>(
                "testAll",&EsysFileWriterTestCase::testAll));
    return testSuite;
}

