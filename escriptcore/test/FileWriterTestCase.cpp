
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

#include <escript/FileWriter.h>

#include "FileWriterTestCase.h"

#include <cppunit/TestCaller.h>
#include <fstream>
#include <sstream>


using namespace CppUnit;
using namespace std;

using escript::FileWriter;

void FileWriterTestCase::testAll()
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

    const char crank = static_cast<char>(mpirank<128?mpirank:128);
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
    cout << "\tTest writeShared." << endl;
    if (mpirank == mpisize-1) {
        oss.write(data, 4);
        CPPUNIT_ASSERT(fw->writeShared(oss) == true);
        CPPUNIT_ASSERT(oss.str().length() == 0);
    }
    fw->close();
#ifdef ESYS_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    CPPUNIT_ASSERT_EQUAL(fileSize(filename), 4L);

    CPPUNIT_ASSERT(fw->openFile(filename) == true);
    oss.write(data, 4);
    cout << "\tTest writeAt." << endl;
    CPPUNIT_ASSERT(fw->writeAt(oss, 16*(mpirank+1)) == true);
    CPPUNIT_ASSERT(oss.str().length() == 0);
    fw->close();
#ifdef ESYS_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    CPPUNIT_ASSERT(fileSize(filename) == 16*mpisize+4);
    delete fw;
}

long FileWriterTestCase::fileSize(string filename)
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

TestSuite* FileWriterTestCase::suite()
{
    TestSuite *testSuite = new TestSuite("FileWriterTestCase");
    testSuite->addTest(new TestCaller<FileWriterTestCase>(
                "testAll",&FileWriterTestCase::testAll));
    return testSuite;
}

