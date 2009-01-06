
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
// MeshWithElements.cpp
//
#include <escriptreader/MeshWithElements.h>
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
MeshWithElements::MeshWithElements() : Mesh(),
    cells(NULL),
    faces(NULL),
    contacts(NULL),
    points(NULL)
{
    name = "Elements";
}

//
//
//
MeshWithElements::MeshWithElements(const MeshWithElements& m) :
    Mesh(m)
{
    nodeTag = m.nodeTag;
    nodeGDOF = m.nodeGDOF;
    nodeGNI = m.nodeGNI;
    nodeGRDFI = m.nodeGRDFI;
    nodeGRNI = m.nodeGRNI;
    cells = new ElementData(*m.cells);
    faces = new ElementData(*m.faces);
    contacts = new ElementData(*m.contacts);
    points = new ElementData(*m.points);
}

//
//
//
MeshWithElements::~MeshWithElements()
{
    delete cells;
    delete faces;
    delete contacts;
    delete points;
}

//
// Returns a vector of strings containing Silo mesh names that have been
// written
//
StringVec MeshWithElements::getMeshNames() const
{
    StringVec res;
    res.push_back(name);
    StringVec tmpVec;
    tmpVec = cells->getMeshNames();
    res.insert(res.end(), tmpVec.begin(), tmpVec.end());
    tmpVec = faces->getMeshNames();
    res.insert(res.end(), tmpVec.begin(), tmpVec.end());
    tmpVec = contacts->getMeshNames();
    res.insert(res.end(), tmpVec.begin(), tmpVec.end());
    tmpVec = points->getMeshNames();
    res.insert(res.end(), tmpVec.begin(), tmpVec.end());
    return res;
}

Mesh* MeshWithElements::getMeshByName(const string name) const
{
    Mesh* ret = NULL;
    if (name == "Elements")
        ret = cells->fullMesh;
    else if (name == "ReducedElements")
        ret = cells->reducedMesh;
    else if (name == "FaceElements")
        ret = faces->fullMesh;
    else if (name == "ReducedFaceElements")
        ret = faces->reducedMesh;
    else if (name == "ContactElements")
        ret = contacts->fullMesh;
    else if (name == "ReducedContactElements")
        ret = contacts->reducedMesh;
    else if (name == "Points")
        ret = points->fullMesh;

    return ret;
}

//
// Returns a vector of strings containing Silo variable names that have
// been written
//
StringVec MeshWithElements::getVarNames() const
{
    StringVec res;
    res.push_back("Nodes_Id");
    res.push_back("Nodes_Tag");
    res.push_back("Nodes_gDOF");
    res.push_back("Nodes_gNI");
    res.push_back("Nodes_grDfI");
    res.push_back("Nodes_grNI");

    StringVec tmpVec;
    tmpVec = cells->getVarNames();
    res.insert(res.end(), tmpVec.begin(), tmpVec.end());
    tmpVec = faces->getVarNames();
    res.insert(res.end(), tmpVec.begin(), tmpVec.end());
    tmpVec = contacts->getVarNames();
    res.insert(res.end(), tmpVec.begin(), tmpVec.end());
    tmpVec = points->getVarNames();
    res.insert(res.end(), tmpVec.begin(), tmpVec.end());

    return res;
}

const IntVec& MeshWithElements::getVarDataByName(const string name) const
{
    if (name == "Nodes_Id")
        return nodeID;
    else if (name == "Nodes_Tag")
        return nodeTag;
    else if (name == "Nodes_gDOF")
        return nodeGDOF;
    else if (name == "Nodes_gNI")
        return nodeGNI;
    else if (name == "Nodes_grDfI")
        return nodeGRDFI;
    else if (name == "Nodes_grNI")
        return nodeGRNI;
    else if (name.compare(0, 9, "Elements_") == 0)
        return cells->getVarDataByName(name);
    else if (name.compare(0, 13, "FaceElements_") == 0)
        return faces->getVarDataByName(name);
    else if (name.compare(0, 16, "ContactElements_") == 0)
        return contacts->getVarDataByName(name);
    else if (name.compare(0, 7, "Points_") == 0)
        return points->getVarDataByName(name);
    else
        return *(IntVec*)(NULL);
}

ElementData* MeshWithElements::getElementsByName(const string name) const
{
    ElementData* ret = NULL;
    if (name == "Elements" || name == "ReducedElements")
        ret = cells;
    else if (name == "FaceElements" || name == "ReducedFaceElements")
        ret = faces;
    else if (name == "ContactElements" || name == "ReducedContactElements")
        ret = contacts;
    else if (name == "Points")
        ret = points;

    return ret;
}

//
// Reads mesh and element data from NetCDF file with given name
//
bool MeshWithElements::readFromNc(const string& filename)
{
    if (!Mesh::readFromNc(filename))
        return false;

    NcError ncerr(NcError::silent_nonfatal);
    NcFile* input;
    NcVar* var;
 
    input = new NcFile(filename.c_str());
    if (!input->is_valid()) {
        cerr << "Could not open input file " << filename << "." << endl;
        delete input;
        return false;
    }

    nodeTag.clear();
    nodeTag.insert(nodeTag.end(), numNodes, 0);
    var = input->get_var("Nodes_Tag");
    var->get(&nodeTag[0], numNodes);

    nodeGDOF.clear();
    nodeGDOF.insert(nodeGDOF.end(), numNodes, 0);
    var = input->get_var("Nodes_gDOF");
    var->get(&nodeGDOF[0], numNodes);

    nodeGNI.clear();
    nodeGNI.insert(nodeGNI.end(), numNodes, 0);
    var = input->get_var("Nodes_gNI");
    var->get(&nodeGNI[0], numNodes);

    nodeGRDFI.clear();
    nodeGRDFI.insert(nodeGRDFI.end(), numNodes, 0);
    var = input->get_var("Nodes_grDfI");
    var->get(&nodeGRDFI[0], numNodes);

    nodeGRNI.clear();
    nodeGRNI.insert(nodeGRNI.end(), numNodes, 0);
    var = input->get_var("Nodes_grNI");
    var->get(&nodeGRNI[0], numNodes);

    // Read all element types
    cells = new ElementData("Elements", this);
    cells->readFromNc(input);
    faces = new ElementData("FaceElements", this);
    faces->readFromNc(input);
    contacts = new ElementData("ContactElements", this);
    contacts->readFromNc(input);
    points = new ElementData("Points", this);
    points->readFromNc(input);

    delete input;
    return true;
}

//
//
//
void MeshWithElements::handleGhostZones(int ownIndex)
{
    cells->handleGhostZones(ownIndex);
    faces->handleGhostZones(ownIndex);
    contacts->handleGhostZones(ownIndex);
    points->handleGhostZones(ownIndex);
#ifdef _DEBUG
    cout << "block " << ownIndex << " has " << cells->getGhostCount()
         << " (" << cells->getReducedGhostCount() << ") ghost zones," << endl;
    cout << faces->getGhostCount() << " (" << faces->getReducedGhostCount()
         << ") ghost faces," << endl;
    cout << contacts->getGhostCount() << " ("
         << contacts->getReducedGhostCount() << ") ghost contacts," << endl;
    cout << points->getGhostCount() << " (" << points->getReducedGhostCount()
         << ") ghost points." << endl;
#endif
}

//
//
//
void MeshWithElements::removeGhostZones()
{
    cells->removeGhostZones();
    faces->removeGhostZones();
    contacts->removeGhostZones();
    points->removeGhostZones();
#ifdef _DEBUG
    cout << "After removing ghost zones there are" << endl;
    cout << "    " << numNodes << " Nodes, ";
    cout << cells->getCount() << " Elements, ";
    cout << faces->getCount() << " Face elements, ";
    cout << contacts->getCount() << " Contact elements, ";
    cout << points->getCount() << " Points left." << endl;
#endif
}

//
//
//
bool MeshWithElements::writeToSilo(DBfile* dbfile, const string& pathInSilo)
{
#if HAVE_SILO
    int ret;

    // Write meshes and zone-centered variables
    if (!cells->writeToSilo(dbfile, pathInSilo) ||
            !faces->writeToSilo(dbfile, pathInSilo) ||
            !contacts->writeToSilo(dbfile, pathInSilo) ||
            !points->writeToSilo(dbfile, pathInSilo))
        return false;

    if (!Mesh::writeToSilo(dbfile, pathInSilo))
        return false;
    
    if (pathInSilo != "") {
        ret = DBSetDir(dbfile, pathInSilo.c_str());
        if (ret != 0)
            return false;
    }

    string siloMeshName = getFullSiloName();

    // Write node-centered variables
    ret = DBPutUcdvar1(dbfile, "Nodes_Tag", siloMeshName.c_str(),
                (float*)&nodeTag[0], numNodes, NULL, 0, DB_INT,
                DB_NODECENT, NULL);
    if (ret == 0)
        ret = DBPutUcdvar1(dbfile, "Nodes_gDOF", siloMeshName.c_str(),
                (float*)&nodeGDOF[0], numNodes, NULL, 0, DB_INT,
                DB_NODECENT, NULL);
    if (ret == 0)
        ret = DBPutUcdvar1(dbfile, "Nodes_gNI", siloMeshName.c_str(),
                (float*)&nodeGNI[0], numNodes, NULL, 0, DB_INT,
                DB_NODECENT, NULL);
    if (ret == 0)
        ret = DBPutUcdvar1(dbfile, "Nodes_grDfI", siloMeshName.c_str(),
                (float*)&nodeGRDFI[0], numNodes, NULL, 0, DB_INT,
                DB_NODECENT, NULL);
    if (ret == 0)
        ret = DBPutUcdvar1(dbfile, "Nodes_grNI", siloMeshName.c_str(),
                (float*)&nodeGRNI[0], numNodes, NULL, 0, DB_INT,
                DB_NODECENT, NULL);

    DBSetDir(dbfile, "/");
    return (ret == 0);

#else // !HAVE_SILO
    return false;
#endif
}

} // namespace EscriptReader

