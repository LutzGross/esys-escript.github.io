
/*******************************************************
*
* Copyright (c) 2003-2008 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/

//
// MPDataSet.cpp
//
#include <escriptreader/MPDataSet.h>
#include <escriptreader/MeshWithElements.h>
#include <escriptreader/DataVar.h>
#include <cstring>
#include <netcdf.hh>
#if HAVE_SILO
#include <silo.h>
#endif

using namespace std;

//
// Constructor
//
MPDataSet::MPDataSet() : numParts(0), keepMesh(false)
{
}

//
// Destructor
//
MPDataSet::~MPDataSet()
{
    VarVector::iterator viIt;
    for (viIt = variables.begin(); viIt != variables.end(); viIt++)
        delete (*viIt).dataVar;

    for (viIt = meshVariables.begin(); viIt != meshVariables.end(); viIt++)
        delete (*viIt).dataVar;

    if (!keepMesh) {
        MeshBlocks::iterator meshIt;
        for (meshIt = meshBlocks.begin(); meshIt != meshBlocks.end(); meshIt++)
            delete *meshIt;
    }
}

//
//
//
bool MPDataSet::load(const string meshFile, const StringVec& varFiles,
                     const StringVec& varNames, int nBlocks)
{
    numParts = nBlocks;
    meshFmt = meshFile;
    
    // initialize variables
    StringVec::const_iterator fileIt = varFiles.begin();
    StringVec::const_iterator nameIt = varNames.begin();
    for (; fileIt != varFiles.end(); fileIt++, nameIt++) {
        VarInfo vi;
        vi.fileName = *fileIt;
        vi.varName = *nameIt;
        vi.valid = true;
        variables.push_back(vi);
    }

    // Load the mesh files
    if (!readMeshes()) {
        cerr << "Reading the domain failed." << endl;
        return false;
    }

    // Retrieve the mesh variables
    convertMeshVariables();

    // Load the variables
    if (!readVariables()) {
        cerr << "Reading the variables failed." << endl;
        return false;
    }

    return true;
}

//
// Load only variables using provided mesh
//
bool MPDataSet::load(const MeshBlocks& m, const string siloFile,
                     const StringVec& varFiles, const StringVec& varNames)
{
    siloMeshFile = siloFile;
    meshBlocks = m;
    numParts = meshBlocks.size();
    keepMesh = true;

    // initialize variables
    StringVec::const_iterator fileIt = varFiles.begin();
    StringVec::const_iterator nameIt = varNames.begin();
    for (; fileIt != varFiles.end(); fileIt++, nameIt++) {
        VarInfo vi;
        vi.fileName = *fileIt;
        vi.varName = *nameIt;
        vi.valid = true;
        variables.push_back(vi);
    }

    // Load the variables
    if (!readVariables()) {
        cerr << "Reading the variables failed." << endl;
        return false;
    }

    return true;
}


//
//
//
bool MPDataSet::readMeshes()
{
    bool ok = true;
    char* str = new char[meshFmt.length()+10];
    for (int idx=0; idx < numParts; idx++) {
        MeshWithElements* meshPart = new MeshWithElements();
        sprintf(str, meshFmt.c_str(), idx);
        string meshfile = str;
        if (meshPart->readFromNc(meshfile)) {
            if (numParts > 1)
                meshPart->handleGhostZones(idx);
            meshBlocks.push_back(meshPart);
        } else {
            delete meshPart;
            ok = false;
            break;
        }
    }
    delete[] str;
    return ok;
}

void MPDataSet::convertMeshVariables()
{
    const StringVec& varNames = meshBlocks[0]->getVarNames();
    StringVec::const_iterator it;
    for (it = varNames.begin(); it != varNames.end(); it++) {
        VarInfo vi;
        vi.varName = *it;
        vi.valid = true;
        DataParts parts;
        // get all parts of current variable
        for (int idx=0; idx < numParts; idx++) {
            const IntVec& data = meshBlocks[idx]->getVarDataByName(*it);
            DataVar* var = new DataVar(*it, data, meshBlocks[idx]);
            parts.push_back(var);
        }
        // if we have more than one part, assemble variable
        if (numParts > 1) {
            assembleVariable(vi, parts);
            for (size_t i=0; i<parts.size(); i++)
                delete parts[i];
        } else
            vi.dataVar = parts.back();

        meshVariables.push_back(vi);
    }
}

//
//
//
bool MPDataSet::readVariables()
{
    if (variables.size() == 0)
        return true;

    VarVector::iterator it;
    bool ok = false;
    for (it = variables.begin(); it != variables.end(); it++) {
        bool curVarOk = true;
        VarInfo& vi = (*it);
        DataParts parts;
        char* str = new char[vi.fileName.length()+10];
        // read all parts of current variable
        for (int idx=0; idx < numParts; idx++) {
            sprintf(str, vi.fileName.c_str(), idx);
            string dfile = str;
            DataVar* var = new DataVar(vi.varName);
            if (var->readFromNc(dfile))
                parts.push_back(var);
            else {
                delete var;
                curVarOk = false;
                break;
            }
        }
        delete[] str;
        if (curVarOk) {
            // at least one variable was read without problems
            ok = true;

            // if we have more than one part, assemble variable
            if (numParts > 1) {
                assembleVariable(vi, parts);
                for (size_t i=0; i<parts.size(); i++)
                    delete parts[i];
            } else
                vi.dataVar = parts.back();
        } else {
            vi.dataVar = NULL;
            vi.valid = false;
        }
    }
    return ok;
}

//
//
//
void MPDataSet::assembleVariable(VarInfo& vi, const DataParts& parts)
{
    DataVar* var = new DataVar(*parts[0]);
    vi.dataVar = var;

    DataParts::const_iterator it;
    for (it = parts.begin()+1; it != parts.end(); it++) {
        var->append(*(*it));
    }
}

bool MPDataSet::saveAsSilo(string siloFile, bool useMultiMesh)
{
#if HAVE_SILO
    if (numParts == 0)
        return false;

    DBfile* dbfile;
    dbfile = DBCreate(siloFile.c_str(), DB_CLOBBER, DB_LOCAL, "Data", DB_PDB);
    if (!dbfile) {
        cerr << "Could not create Silo file." << endl;
        return false;
    }

    MeshBlocks::iterator meshIt;
    VarVector::iterator viIt;
    int idx = 0;
    for (meshIt = meshBlocks.begin(); meshIt != meshBlocks.end(); meshIt++, idx++) {
        string siloPath("");
        //if (numParts > 1) {
            char str[64];
            snprintf(str, 64, "/block%04d", idx);
            siloPath = str;
            DBMkdir(dbfile, siloPath.c_str());
        //}
        // write block of the mesh if we don't use an external mesh
        if (!siloMeshFile.length()) {
            if (! (*meshIt)->writeToSilo(dbfile, siloPath)) {
                cerr << "Error writing mesh of block " << idx << " to Silo!\n";
                break;
            }
        }

        // write variables for current mesh block
        for (viIt = variables.begin(); viIt != variables.end(); viIt++) {
            // do not attempt to write this variable if previous steps failed
            if (!(*viIt).valid) continue;
            DataVar* var = (*viIt).dataVar;
            if (!var->setMesh(*meshIt) || !var->writeToSilo(dbfile, siloPath)) {
                cerr << "Error writing block " << idx << " of "
                    << var->getName() << " to Silo!" << endl;
                (*viIt).valid = false;
            }
        }
    }

    if (useMultiMesh) {
        const StringVec& meshNames = meshBlocks[0]->getMeshNames();
        StringVec::const_iterator it;
        for (it = meshNames.begin(); it != meshNames.end(); it++)
            putSiloMultiMesh(dbfile, *it);

        const StringVec& varNames = meshBlocks[0]->getVarNames();
        for (it = varNames.begin(); it != varNames.end(); it++)
            putSiloMultiVar(dbfile, *it, true);

        for (viIt = variables.begin(); viIt != variables.end(); viIt++) {
            if (!(*viIt).valid) continue;
            DataVar* var = (*viIt).dataVar;
            if (var->getRank() < 2)
                putSiloMultiVar(dbfile, var->getName());
            else
                putSiloMultiTensor(dbfile, var);
        }
    }

    vector<char*> tensorNames;
    vector<string> tensorDefStrings;
    vector<char*> tensorDefs;

    // collect tensors for their Silo definitions
    for (viIt = variables.begin(); viIt != variables.end(); viIt++) {
        if (!(*viIt).valid) continue;
        DataVar* var = (*viIt).dataVar;
        if (var->getRank() == 2) {
            tensorDefStrings.push_back(var->getTensorDef());
            tensorDefs.push_back((char*)tensorDefStrings.back().c_str());
            tensorNames.push_back((char*)var->getName().c_str());
        }
    }

    if (tensorDefs.size()) {
        vector<int> defTypes(tensorDefs.size(), DB_VARTYPE_TENSOR);
        DBPutDefvars(dbfile, "tensors", tensorDefs.size(), &tensorNames[0],
                &defTypes[0], &tensorDefs[0], NULL);
    }

    DBClose(dbfile);
    return true;

#else // !HAVE_SILO
    return false;
#endif
}

void MPDataSet::putSiloMultiMesh(DBfile* dbfile, string meshName)
{
#if HAVE_SILO
    vector<int> meshtypes(meshBlocks.size(), DB_UCDMESH);
    vector<string> tempstrings;
    vector<char*> meshnames;
    vector<Mesh*>::iterator it;

    for (size_t idx = 0; idx < meshBlocks.size(); idx++) {
        string siloPath = meshBlocks[idx]->getSiloPath();
        tempstrings.push_back(siloPath + string("/") + meshName);
        meshnames.push_back((char*)tempstrings.back().c_str());
    }
    DBPutMultimesh(dbfile, meshName.c_str(), meshBlocks.size(), &meshnames[0],
            &meshtypes[0], NULL);
#endif
}

void MPDataSet::putSiloMultiVar(DBfile* dbfile, string varName,
                                bool useMeshFile)
{
#if HAVE_SILO
    vector<int> vartypes(meshBlocks.size(), DB_UCDVAR);
    vector<string> tempstrings;
    vector<char*> varnames;
    for (size_t idx = 0; idx < meshBlocks.size(); idx++) {
        string siloPath;
        if (useMeshFile)
            siloPath = meshBlocks[idx]->getSiloPath();
        else {
            char str[64];
            snprintf(str, 64, "/block%04d", static_cast<int>(idx));
            siloPath = str;
        }
        tempstrings.push_back(siloPath + string("/") + varName);
        varnames.push_back((char*)tempstrings.back().c_str());
    }
    DBPutMultivar(dbfile, varName.c_str(), meshBlocks.size(), &varnames[0],
            &vartypes[0], NULL);
#endif
}

void MPDataSet::putSiloMultiTensor(DBfile* dbfile, const DataVar* var)
{
#if HAVE_SILO
    string tensorDir = var->getName()+string("_comps/");
    DBSetDir(dbfile, "/");
    DBMkdir(dbfile, tensorDir.c_str());
    int one = 1;
    DBoptlist* optList = DBMakeOptlist(1);
    DBAddOption(optList, DBOPT_HIDE_FROM_GUI, &one);
    vector<int> vartypes(meshBlocks.size(), DB_UCDVAR);
    const IntVec& shape = var->getShape();
    for (int i=0; i<shape[1]; i++) {
        for (int j=0; j<shape[0]; j++) {
            vector<string> tempstrings;
            vector<char*> varnames;
            char comp[255];
            snprintf(comp, 255, "%s_comps/a_%d%d", var->getName().c_str(), i,j);
            for (size_t idx = 0; idx < meshBlocks.size(); idx++) {
                string siloPath;
                char str[64];
                snprintf(str, 64, "/block%04d", static_cast<int>(idx));
                siloPath = str;
                tempstrings.push_back(siloPath + string("/") + string(comp));
                varnames.push_back((char*)tempstrings.back().c_str());
            }
            DBPutMultivar(dbfile, comp, meshBlocks.size(), &varnames[0],
                    &vartypes[0], optList);
        }
    }
    DBFreeOptlist(optList);
#endif
}

