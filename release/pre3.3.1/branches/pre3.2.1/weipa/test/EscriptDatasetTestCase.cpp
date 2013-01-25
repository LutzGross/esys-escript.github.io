
/*******************************************************
*
* Copyright (c) 2003-2011 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


#include "escript/DataFactory.h"
#include "finley/CppAdapter/MeshAdapterFactory.h"
#include "weipa/EscriptDataset.h"
#include "EscriptDatasetTestCase.h"

using namespace CppUnitTest;
using namespace escript;
using namespace weipa;
using namespace std;

void EscriptDatasetTestCase::setUp()
{
    // This is called before each test is run
}

void EscriptDatasetTestCase::tearDown()
{
    // This is called after each test has been run
}

void EscriptDatasetTestCase::testAll()
{
    cout << endl;
    cout << "\tTest default constructor." << endl;
    EscriptDataset_ptr dataset(new EscriptDataset());

    cout << "\tTest saveSilo without data." << endl;
    assert(dataset->saveSilo("dummy") == false);

    cout << "\tTest saveVTK without data." << endl;
    assert(dataset->saveVTK("dummy") == false);

    cout << "\tTest getConvertedDomain without data." << endl;
    assert(dataset->getConvertedDomain().size() == 0);

    cout << "\tTest getVariables without data." << endl;
    assert(dataset->getVariables().size() == 0);

    cout << "\tTest getMeshVariables without data." << endl;
    assert(dataset->getMeshVariables().size() == 0);

    // instantiate a domain and data
    Domain_ptr dom(finley::brick());
    escript::Data data = Scalar(0.0, continuousFunction(*dom), true);

    cout << "\tTest addData with NULL domain." << endl;
    assert(dataset->addData(data, "foo", "bar") == false);

    cout << "\tTest setDomain." << endl;
    assert(dataset->setDomain(dom.get()) == true);
    assert(dataset->getMeshVariables().size() > 0);

    cout << "\tTest bogus setDomain call." << endl;
    assert(dataset->setDomain(dom.get()) == false);

    cout << "\tTest getConvertedDomain." << endl;
    DomainChunks chunks = dataset->getConvertedDomain();
    assert(chunks.size() > 0);

    StringVec varfiles, varnames;
    varfiles.push_back("testvar%04d.nc");
    varnames.push_back("testvar");
    cout << "\tTest bogus loadNetCDF call 1." << endl;
    assert(dataset->loadNetCDF("mesh%04d.nc", varfiles, varnames, 1) == false);

    cout << "\tTest bogus loadNetCDF call 2." << endl;
    assert(dataset->loadNetCDF(chunks, varfiles, varnames) == false);

    cout << "\tTest addData with valid data." << endl;
    assert(dataset->addData(data, "testvar", "cm") == true);
    assert(dataset->getVariables().size() == 1);

    cout << "\tTest set/getCycleAndTime." << endl;
    dataset->setCycleAndTime(42, 3.1415);
    assert(dataset->getCycle() == 42);
    assert(dataset->getTime()-3.1415 < 0.001);

    dataset->setMetadataSchemaString("xmlns:test=\"http://myschema.com/test\"",
            "<MyValue>4711</MyValue>");
    dataset->setMeshLabels("x-axis", "y-axis", "z-axis");
    dataset->setMeshUnits("km", "cm", "mm");

#if USE_SILO
    cout << "\tTest saveSilo." << endl;
    assert(dataset->saveSilo("weipatest.silo") == true);
    ifstream f("weipatest.silo");
    assert(f.is_open());
    f.close();
#endif

    cout << "\tTest saveVTK." << endl;
    assert(dataset->saveVTK("weipatest.vtu") == true);
    checkVTKfile("weipatest.vtu");

    //varnames.push_back("dummy");
    //cout << "\tTest loadNetCDF with invalid params." << endl;
    //assert(dataset->loadNetCDF(blocks, varfiles, varnames) == false);
}

TestSuite* EscriptDatasetTestCase::suite()
{
    //
    // create the suite of tests to perform.
    TestSuite *testSuite = new TestSuite("EscriptDatasetTestCase");

    testSuite->addTest(new TestCaller<EscriptDatasetTestCase>(
                "testAll",&EscriptDatasetTestCase::testAll));
    return testSuite;
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
    assert(f.is_open());

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
    assert(numPoints>0);
    assert(numCells>0);

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
                    assertLongsEqual(numPoints, getDataArrayLength(f));
                }
            }
        } else if (s.find("<Cells>") == 0) {
            cellsFound=true;
            // check cell info (connectivity, offsets, types)
            while (f.good() && s.find("</Cells>") != 0) {
                f.getline(line, 256);
                s = line;
                if (s.find("<DataArray") == 0) {
                    assertLongsEqual(numCells, getDataArrayLength(f));
                }
            }
        } else if (s.compare("<PointData>") == 0) {
            // check nodal data
            while (f.good() && s.find("</PointData>") != 0) {
                f.getline(line, 256);
                s = line;
                if (s.find("<DataArray") == 0) {
                    numPointData++;
                    assertLongsEqual(numPoints, getDataArrayLength(f));
                }
            }
        } else if (s.find("<CellData>") == 0) {
            // check cell data
            while (f.good() && s.find("</CellData>") != 0) {
                f.getline(line, 256);
                s = line;
                if (s.find("<DataArray") == 0) {
                    numCellData++;
                    assertLongsEqual(numCells, getDataArrayLength(f));
                }
            }
        }
    }

    assert(pointsFound);
    assert(cellsFound);
    assert(numPointData>0);
    assert(numCellData>0);
}

