
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

#include <escriptexport/EscriptDataset.h>
#include <escriptexport/FinleyMesh.h>
#include <escriptexport/DataVar.h>

#if USE_SILO
#include <silo.h>
#endif

#include <cstring>
#include <fstream>
#include <iostream>

using namespace std;
using namespace escriptexport;

string insertTimestep(const string& fString, int timeStep, int tsMultiplier)
{
    string s(fString);
    size_t pos;
    if ((pos = s.find("%")) != s.npos) {
        size_t endpos = pos+1;
        while (endpos<s.length() && s[endpos] != 'd')
            endpos++;
        string fmtStr = s.substr(pos, endpos-pos+1);
        char ts[255];
        snprintf(ts, 255, fmtStr.c_str(), timeStep*tsMultiplier);
        s.replace(pos, endpos-pos+1, ts);
    }
    return s;
}

int main(int argc, char** argv)
{
#if HAVE_MPI
    MPI_Init(&argc, &argv);
#endif

#if USE_SILO
    // turn off for debugging purposes
    bool writeMultiMesh = true;

    // whether time-varying datasets should use the same mesh (from T=0)
    // TODO: Add a command line option for this
    bool writeMeshOnce = true;

    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <file.esd>" << endl; fflush(stderr);
        return -1;
    }

    string esdFile(argv[1]);

    ifstream in(esdFile.c_str());
    if (!in.is_open()) {
        cerr << "Could not open " << esdFile << "." << endl;
        return -1;
    }

    // first line (header) should be "#escript datafile Vx.y"
    char line[256];
    in.getline(line, 256);
    int major, minor;
    if (sscanf(line, "#escript datafile V%d.%d", &major, &minor) != 2) {
        cerr << esdFile << " is not a valid escript datafile." << endl;
        in.close();
        return -1;
    }

    int nParts=0, nTimesteps=1, tsMultiplier=1;
    string meshFile;
    StringVec varFiles;
    StringVec varNames;

    while (in.good()) {
        in.getline(line, 256);
        if (line[0] == '#' || strlen(line) == 0)
            continue;
        int iVal;
        char sVal[256];
        if (sscanf(line, "N=%d", &iVal) == 1)
            nParts = iVal;
        else if (sscanf(line, "T=%d", &iVal) == 1)
            nTimesteps = iVal;
        else if (sscanf(line, "DT=%d", &iVal) == 1)
            tsMultiplier = iVal;
        else if (sscanf(line, "M=%s", sVal) == 1)
            meshFile = sVal;
        else if (sscanf(line, "V=%s", sVal) == 1 && strchr(sVal, ':')) {
            // split into filename and variable name
            char* colon = strchr(sVal, ':');
            *colon = 0;
            varFiles.push_back(sVal);
            varNames.push_back(colon+1);
        } else {
            cerr << esdFile << " is not a valid escript datafile." << endl;
            in.close();
            return -1;
        }
    }

    in.close();
    
    if (nParts < 1 || meshFile == "" || nTimesteps < 1 || tsMultiplier < 1) {
        cerr << esdFile << " is not a valid escript datafile." << endl;
        return -1;
    }

    cout << "Converting " << esdFile << "..." << endl;

    MeshBlocks meshFromTzero;

    // load and save all timesteps
    for (int timeStep = 0; timeStep < nTimesteps; timeStep++) {
        StringVec varFilesTS;
        StringVec::const_iterator it;

        // look for "%d" in filename and replace by timestep*multiplier if found
        for (it=varFiles.begin(); it!=varFiles.end(); it++) {
            string v = insertTimestep(*it, timeStep, tsMultiplier);
            if (nParts > 1)
                v.append(".nc.%04d");
            else
                v.append(".nc");
            varFilesTS.push_back(v);
        }

        if (nTimesteps > 1)
            cout << "T = " << timeStep << endl;

        EscriptDataset* ds;
#if HAVE_MPI
        ds = new EscriptDataset(MPI_COMM_WORLD);
#else
        ds = new EscriptDataset();
#endif

        if (writeMeshOnce && timeStep > 0) {
            if (!ds->loadNetCDF(meshFromTzero, varFilesTS, varNames)) {
                delete ds;
                break;
            }
        } else {
            string meshTS = insertTimestep(meshFile, timeStep, tsMultiplier);
            if (nParts > 1)
                meshTS.append(".nc.%04d");
            else
                meshTS.append(".nc");

            if (!ds->loadNetCDF(meshTS, varFilesTS, varNames, nParts)) {
                delete ds;
                break;
            }
        }

        string baseName(esdFile);
        size_t dot = baseName.rfind('.');
        if (dot != baseName.npos)
            baseName.erase(dot, baseName.length()-dot);

        ostringstream siloFile;
        siloFile << baseName;
        if (nTimesteps > 1)
            siloFile << "." << timeStep;
        siloFile << ".silo";

        ds->setCycleAndTime(timeStep, (double)timeStep);
        ds->saveSilo(siloFile.str(), writeMultiMesh);

        // keep mesh from first timestep if it should be reused
        if (writeMeshOnce && nTimesteps > 1 && timeStep == 0) {
            meshFromTzero = ds->extractMesh();
            meshFile = siloFile.str();
            MeshBlocks::iterator meshIt;
            for (meshIt = meshFromTzero.begin(); meshIt != meshFromTzero.end();
                    meshIt++)
            {
                // Prepend Silo mesh paths with the filename of the mesh to be
                // used
                string fullSiloPath = meshFile + string(":");
                fullSiloPath += (*meshIt)->getSiloPath();
                (*meshIt)->setSiloPath(fullSiloPath);
            }
        }

        delete ds;
    }
    
    // clean up
    MeshBlocks::iterator meshIt;

    cout << "All done." << endl;
#else // !USE_SILO
    cerr << "Error: escriptconvert was compiled without Silo support! Exiting."
        << endl;
#endif

#if HAVE_MPI
    MPI_Finalize();
#endif

    return 0;
}

