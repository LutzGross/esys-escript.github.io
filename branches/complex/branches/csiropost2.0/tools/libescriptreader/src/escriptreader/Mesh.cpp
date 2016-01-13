
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
// Mesh.cpp
//
#include <escriptreader/Mesh.h>
#include <escriptreader/ElementData.h>
#include <netcdf.hh>
#if HAVE_SILO
#include <silo.h>
#endif

using namespace std;

namespace EscriptReader {

//
//
//
Mesh::Mesh(CoordArray c, int nDims, int nNodes) :
    coords(c), numDims(nDims), numNodes(nNodes), name("Mesh")
{
}

//
//
//
Mesh::Mesh(const Mesh& m)
{
    numDims = m.numDims;
    numNodes = m.numNodes;
    nodeID = m.nodeID;
    nodeID2idx = m.nodeID2idx;
    name = m.name;
    siloPath = m.siloPath;
    for (int i=0; i<numDims; i++) {
        float* c = new float[numNodes];
        copy(m.coords[i], m.coords[i]+numNodes, c);
        coords.push_back(c);
    }
}

//
//
//
Mesh::~Mesh()
{
    CoordArray::iterator it;
    for (it = coords.begin(); it != coords.end(); it++)
        delete[] *it;
}

//
//
//
bool Mesh::readFromNc(const string& ncFile)
{
    NcError ncerr(NcError::silent_nonfatal);
    NcFile* input;
    NcAtt* att;
    NcVar* var;
 
    input = new NcFile(ncFile.c_str());
    if (!input->is_valid()) {
        cerr << "Could not open input file " << ncFile.c_str() << ".\n";
        delete input;
        return false;
    }

    att = input->get_att("numDim");
    numDims = att->as_int(0);

    att = input->get_att("numNodes");
    numNodes = att->as_int(0);

    coords.clear();
    var = input->get_var("Nodes_Coordinates");
    for (int i=0; i<numDims; i++) {
        float* c = new float[numNodes];
        var->set_cur(0, i);
        var->get(c, numNodes, 1);
        coords.push_back(c);
    }

    nodeID.clear();
    nodeID.insert(nodeID.end(), numNodes, 0);
    var = input->get_var("Nodes_Id");
    var->get(&nodeID[0], numNodes);

    buildIndexMap();

    delete input;
    return true;
}

//
//
//
bool Mesh::writeToSilo(DBfile* dbfile, const string& pathInSilo)
{
#if HAVE_SILO
    int ret;

    if (pathInSilo != "") {
        siloPath = pathInSilo;
        ret = DBSetDir(dbfile, siloPath.c_str());
        if (ret != 0)
            return false;
    } else {
        siloPath = "/";
    }
   
    // Write node-centered variable
    ret = DBPutUcdvar1(dbfile, "Nodes_Id", getFullSiloName().c_str(),
            (float*)&nodeID[0], numNodes, NULL, 0, DB_INT, DB_NODECENT, NULL);

    DBSetDir(dbfile, "/");
    return (ret == 0);

#else // !HAVE_SILO
    return false;
#endif
}

} // namespace EscriptReader

