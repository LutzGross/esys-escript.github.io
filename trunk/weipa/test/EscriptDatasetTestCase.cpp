
/*****************************************************************************
*
* Copyright (c) 2003-2016 by The University of Queensland
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

#include "EscriptDatasetTestCase.h"

#include <weipa/EscriptDataset.h>

#include <escript/DataFactory.h>
#include <escript/FunctionSpaceFactory.h>

#if USE_DUDLEY
#include <dudley/DomainFactory.h>
#endif
#if USE_FINLEY
#include <finley/CppAdapter/MeshAdapterFactory.h>
#endif
#if USE_RIPLEY
#include <ripley/Brick.h>
#endif
#if USE_SPECKLEY
#include <speckley/Brick.h>
#endif

#include <cppunit/TestCaller.h>

using namespace CppUnit;
using namespace escript;
using namespace weipa;
using namespace std;

TestSuite* EscriptDatasetTestCase::suite()
{
    TestSuite *testSuite = new TestSuite("EscriptDatasetTestCase");
    testSuite->addTest(new TestCaller<EscriptDatasetTestCase>(
                "testBase",&EscriptDatasetTestCase::testBase));
#if USE_DUDLEY
    testSuite->addTest(new TestCaller<EscriptDatasetTestCase>(
                "testDudley",&EscriptDatasetTestCase::testDudley));
#endif
#if USE_FINLEY
    testSuite->addTest(new TestCaller<EscriptDatasetTestCase>(
                "testFinley",&EscriptDatasetTestCase::testFinley));
#endif
#if USE_RIPLEY
    testSuite->addTest(new TestCaller<EscriptDatasetTestCase>(
                "testRipley",&EscriptDatasetTestCase::testRipley));
#endif
#if USE_SPECKLEY
    testSuite->addTest(new TestCaller<EscriptDatasetTestCase>(
                "testSpeckley",&EscriptDatasetTestCase::testSpeckley));
#endif
    return testSuite;
}

void EscriptDatasetTestCase::testBase()
{
    cout << endl;
    cout << "\tTest default constructor." << endl;
    EscriptDataset_ptr dataset(new EscriptDataset());

    cout << "\tTest saveSilo without data." << endl;
    CPPUNIT_ASSERT(dataset->saveSilo("dummy") == false);

    cout << "\tTest saveVTK without data." << endl;
    CPPUNIT_ASSERT(dataset->saveVTK("dummy") == false);

    cout << "\tTest getConvertedDomain without data." << endl;
    CPPUNIT_ASSERT(dataset->getConvertedDomain().size() == 0);

    cout << "\tTest getVariables without data." << endl;
    CPPUNIT_ASSERT(dataset->getVariables().size() == 0);

    cout << "\tTest getMeshVariables without data." << endl;
    CPPUNIT_ASSERT(dataset->getMeshVariables().size() == 0);
}

#if USE_DUDLEY
void EscriptDatasetTestCase::testDudley()
{
    JMPI info=makeInfo(MPI_COMM_WORLD);
    Domain_ptr dom(dudley::brick(info));
    cout << "Running Dudley tests..." << endl;
    runDomainTests(dom);
}
#endif

#if USE_FINLEY
void EscriptDatasetTestCase::testFinley()
{
    JMPI info=makeInfo(MPI_COMM_WORLD);
    Domain_ptr dom(finley::brick(info));
    cout << "Running Finley tests..." << endl;
    runDomainTests(dom);
}
#endif

#if USE_RIPLEY
void EscriptDatasetTestCase::testRipley()
{
    Domain_ptr dom(new ripley::Brick(5,4,3, 0,0,0, 1,1,1));
    cout << "Running Ripley tests..." << endl;
    runDomainTests(dom);
}
#endif

#if USE_SPECKLEY
void EscriptDatasetTestCase::testSpeckley()
{
    for (int i = 2; i < 11; i++) {
        Domain_ptr dom(new speckley::Brick(i, 5,4,3, 0,0,0, 1,1,1));
        cout << "Running Speckley order " << i << " tests..." << endl;
        runDomainTests(dom);
    }
}
#endif

void EscriptDatasetTestCase::runDomainTests(Domain_ptr dom)
{
    EscriptDataset_ptr dataset(new EscriptDataset());
    Data data = Scalar(0.0, continuousFunction(*dom), true);

    cout << "\tTest addData with NULL domain." << endl;
    CPPUNIT_ASSERT(dataset->addData(data, "foo", "bar") == false);

    cout << "\tTest setDomain." << endl;
    CPPUNIT_ASSERT(dataset->setDomain(dom.get()) == true);
    CPPUNIT_ASSERT(dataset->getMeshVariables().size() > 0);

    cout << "\tTest bogus setDomain call." << endl;
    CPPUNIT_ASSERT(dataset->setDomain(dom.get()) == false);

    cout << "\tTest getConvertedDomain." << endl;
    DomainChunks chunks = dataset->getConvertedDomain();
    CPPUNIT_ASSERT(chunks.size() > 0);

    StringVec varfiles, varnames;
    varfiles.push_back("testvar%04d.nc");
    varnames.push_back("testvar");
    cout << "\tTest bogus loadNetCDF call 1." << endl;
    CPPUNIT_ASSERT(dataset->loadNetCDF("mesh%04d.nc", varfiles, varnames, 1) == false);

    cout << "\tTest bogus loadNetCDF call 2." << endl;
    CPPUNIT_ASSERT(dataset->loadNetCDF(chunks, varfiles, varnames) == false);

    cout << "\tTest addData with valid data." << endl;
    CPPUNIT_ASSERT(dataset->addData(data, "testvar", "cm") == true);
    CPPUNIT_ASSERT(dataset->getVariables().size() == 1);

    cout << "\tTest set/getCycleAndTime." << endl;
    dataset->setCycleAndTime(42, 3.1415);
    CPPUNIT_ASSERT(dataset->getCycle() == 42);
    CPPUNIT_ASSERT(dataset->getTime()-3.1415 < 0.001);

    dataset->setMetadataSchemaString("xmlns:test=\"http://myschema.com/test\"",
            "<MyValue>4711</MyValue>");
    dataset->setMeshLabels("x-axis", "y-axis", "z-axis");
    dataset->setMeshUnits("km", "cm", "mm");
    dataset->setSaveMeshData(true);

#ifdef ESYS_HAVE_SILO
    cout << "\tTest saveSilo." << endl;
    CPPUNIT_ASSERT(dataset->saveSilo("domaintest.silo") == true);
    ifstream f("domaintest.silo");
    CPPUNIT_ASSERT(f.is_open());
    f.close();
#endif

    cout << "\tTest saveVTK." << endl;
    CPPUNIT_ASSERT(dataset->saveVTK("domaintest.vtu") == true);
    checkVTKfile("domaintest.vtu");
}

int EscriptDatasetTestCase::getDataArrayLength(std::istream& is)
{
    int length=0;
    char line[256];
    while (is.good()) {
        is.getline(line, 256);
        string s(line);
        if (s.find("</DataArray") != 0)
            length++;
        else
            break;
    }
    return length;
}

void EscriptDatasetTestCase::checkVTKfile(std::string filename)
{
    ifstream f(filename.c_str());
    CPPUNIT_ASSERT(f.is_open());

    char line[256];
    int numPoints=0, numCells=0;
    while (f.good()) {
        f.getline(line, 256);
        string s(line);
        size_t pp = s.find("NumberOfPoints=");
        size_t cp = s.find("NumberOfCells=");
        if (pp!=s.npos && cp!=s.npos) {
            stringstream ss;
            string tmp(s.substr(pp+16));
            ss.str(tmp);
            ss >> numPoints;
            tmp = s.substr(cp+15);
            ss.str(tmp);
            ss >> numCells;
            break;
        }
    }
    CPPUNIT_ASSERT(numPoints>0);
    CPPUNIT_ASSERT(numCells>0);

    bool pointsFound=false, cellsFound=false;
    int numPointData=0, numCellData=0;

    while (f.good()) {
        f.getline(line, 256);
        string s(line);
        if (s.compare("<Points>") == 0) {
            pointsFound=true;
            // check node coordinates
            while (f.good() && s.find("</Points>") != 0) {
                f.getline(line, 256);
                s = line;
                if (s.find("<DataArray") == 0) {
                    CPPUNIT_ASSERT_EQUAL(numPoints, getDataArrayLength(f));
                }
            }
        } else if (s.find("<Cells>") == 0) {
            cellsFound=true;
            // check cell info (connectivity, offsets, types)
            while (f.good() && s.find("</Cells>") != 0) {
                f.getline(line, 256);
                s = line;
                if (s.find("<DataArray") == 0) {
                    CPPUNIT_ASSERT_EQUAL(numCells, getDataArrayLength(f));
                }
            }
        } else if (s.compare("<PointData>") == 0) {
            // check nodal data
            while (f.good() && s.find("</PointData>") != 0) {
                f.getline(line, 256);
                s = line;
                if (s.find("<DataArray") == 0) {
                    numPointData++;
                    CPPUNIT_ASSERT_EQUAL(numPoints, getDataArrayLength(f));
                }
            }
        } else if (s.find("<CellData>") == 0) {
            // check cell data
            while (f.good() && s.find("</CellData>") != 0) {
                f.getline(line, 256);
                s = line;
                if (s.find("<DataArray") == 0) {
                    numCellData++;
                    CPPUNIT_ASSERT_EQUAL(numCells, getDataArrayLength(f));
                }
            }
        }
    }

    CPPUNIT_ASSERT(pointsFound);
    CPPUNIT_ASSERT(cellsFound);
    CPPUNIT_ASSERT(numPointData>0);
    CPPUNIT_ASSERT(numCellData>0);
}

