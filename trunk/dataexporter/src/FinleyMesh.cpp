
/*******************************************************
*
* Copyright (c) 2003-2009 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/

#include <escriptexport/FinleyMesh.h>
#include <escriptexport/ElementData.h>
#include <escriptexport/NodeData.h>

#include <finley/CppAdapter/MeshAdapter.h>
extern "C" {
#include <finley/Mesh.h>
}

#include <iostream>

#if USE_NETCDF
#include <netcdf.hh>
#endif

#if USE_SILO
#include <silo.h>
#endif

using namespace std;

namespace escriptexport {

//
//
//
FinleyMesh::FinleyMesh() :
    initialized(false)
{
}

//
//
//
FinleyMesh::FinleyMesh(const FinleyMesh& m)
{
    nodes = NodeData_ptr(new NodeData(*m.nodes));
    cells = ElementData_ptr(new ElementData(*m.cells));
    faces = ElementData_ptr(new ElementData(*m.faces));
    contacts = ElementData_ptr(new ElementData(*m.contacts));
    initialized = m.initialized;
}

//
//
//
FinleyMesh::~FinleyMesh()
{
    cleanup();
}

//
//
//
void FinleyMesh::cleanup()
{
    nodes.reset();
    cells.reset();
    faces.reset();
    contacts.reset();
    initialized = false;
}

//
//
//
bool FinleyMesh::initFromEscript(escript::const_Domain_ptr escriptDomain)
{
    cleanup();

    finleyMesh = dynamic_cast<const finley::MeshAdapter*>(escriptDomain.get())->getFinley_Mesh();
    if (!finleyMesh) {
        return false;
    }

    nodes = NodeData_ptr(new NodeData("Elements"));
    cells = ElementData_ptr(new ElementData("Elements", nodes));
    faces = ElementData_ptr(new ElementData("FaceElements", nodes));
    contacts = ElementData_ptr(new ElementData("ContactElements", nodes));

    if (nodes->initFromFinley(finleyMesh->Nodes) &&
            cells->initFromFinley(finleyMesh->Elements) &&
            faces->initFromFinley(finleyMesh->FaceElements) &&
            contacts->initFromFinley(finleyMesh->ContactElements)) {
        initialized = true;
    }

    return initialized;
}

//
// Reads mesh and element data from NetCDF file with given name
//
bool FinleyMesh::initFromNetCDF(const string& filename)
{
    cleanup();
    
#if USE_NETCDF
    NcError ncerr(NcError::silent_nonfatal);
    NcFile* input;
 
    input = new NcFile(filename.c_str());
    if (!input->is_valid()) {
        cerr << "Could not open input file " << filename << "." << endl;
        delete input;
        return false;
    }

    nodes = NodeData_ptr(new NodeData("Elements"));
    if (!nodes->readFromNc(input))
        return false;

    // Read all element types
    cells = ElementData_ptr(new ElementData("Elements", nodes));
    cells->readFromNc(input);
    faces = ElementData_ptr(new ElementData("FaceElements", nodes));
    faces->readFromNc(input);
    contacts = ElementData_ptr(new ElementData("ContactElements", nodes));
    contacts->readFromNc(input);

    delete input;
    initialized = true;
#endif

    return initialized;
}

//
//
//
NodeData_ptr FinleyMesh::getMeshForFinleyFS(int functionSpace) const
{
    NodeData_ptr result;

    if (!initialized)
        return result;

    ElementData_ptr elements = getElementsForFinleyFS(functionSpace);
    if (elements != NULL) {
        if (functionSpace == FINLEY_NODES)
            result = elements->getNodeMesh();
        else {
            if (elements->getReducedNumElements() > 0)
                result = elements->getReducedNodeMesh();
            else
                result = elements->getNodeMesh();
        }
    }
 
    return result;
}

//
//
//
ElementData_ptr FinleyMesh::getElementsForFinleyFS(int functionSpace) const
{
    ElementData_ptr result;

    if (!initialized) {
        return result;
    }

    switch (functionSpace) {
        case FINLEY_NODES:
        case FINLEY_REDUCED_NODES:
        case FINLEY_REDUCED_ELEMENTS:
        case FINLEY_ELEMENTS:
            result = cells;
            break;

        case FINLEY_REDUCED_FACE_ELEMENTS:
        case FINLEY_FACE_ELEMENTS:
            result = faces;
            break;

        case FINLEY_REDUCED_CONTACT_ELEMENTS_1:
        case FINLEY_REDUCED_CONTACT_ELEMENTS_2:
        case FINLEY_CONTACT_ELEMENTS_1:
        case FINLEY_CONTACT_ELEMENTS_2:
            result = contacts;
            break;

        default:
            cerr << "Unsupported function space type " << functionSpace
                << "!" << endl;
    }
    return result;
}

//
// Returns a vector of strings containing Silo mesh names that have been
// written
//
StringVec FinleyMesh::getMeshNames() const
{
    StringVec res;
    if (initialized) {
        StringVec tmpVec;
        tmpVec = cells->getMeshNames();
        res.insert(res.end(), tmpVec.begin(), tmpVec.end());
        tmpVec = faces->getMeshNames();
        res.insert(res.end(), tmpVec.begin(), tmpVec.end());
        tmpVec = contacts->getMeshNames();
        res.insert(res.end(), tmpVec.begin(), tmpVec.end());
    }
    return res;
}

//
//
//
NodeData_ptr FinleyMesh::getMeshByName(const string& name) const
{
    NodeData_ptr ret;
    if (!initialized) {
        return ret;
    }
    if (name == "Elements")
        ret = cells->getNodeMesh();
    else if (name == "ReducedElements")
        ret = cells->getReducedNodeMesh();
    else if (name == "FaceElements")
        ret = faces->getNodeMesh();
    else if (name == "ReducedFaceElements")
        ret = faces->getReducedNodeMesh();
    else if (name == "ContactElements")
        ret = contacts->getNodeMesh();
    else if (name == "ReducedContactElements")
        ret = contacts->getReducedNodeMesh();

    return ret;
}

//
// Returns a vector of strings containing Silo variable names that have
// been written
//
StringVec FinleyMesh::getVarNames() const
{
    StringVec res;
 
    if (initialized) {
        res = nodes->getVarNames();
        StringVec tmpVec = cells->getVarNames();
        res.insert(res.end(), tmpVec.begin(), tmpVec.end());
        tmpVec = faces->getVarNames();
        res.insert(res.end(), tmpVec.begin(), tmpVec.end());
        tmpVec = contacts->getVarNames();
        res.insert(res.end(), tmpVec.begin(), tmpVec.end());
    }

    return res;
}

//
//
//
const IntVec& FinleyMesh::getVarDataByName(const string& name) const
{
    if (!initialized) {
        throw "Mesh not initialized";
    }
    
    if (name.compare(0, 6, "Nodes_") == 0)
        return nodes->getVarDataByName(name);
    else if (name.compare(0, 9, "Elements_") == 0)
        return cells->getVarDataByName(name);
    else if (name.compare(0, 13, "FaceElements_") == 0)
        return faces->getVarDataByName(name);
    else if (name.compare(0, 16, "ContactElements_") == 0)
        return contacts->getVarDataByName(name);
    else
        throw "Invalid variable name";
}

//
//
//
ElementData_ptr FinleyMesh::getElementsByName(const string& name) const
{
    ElementData_ptr ret;
    if (name == "Elements" || name == "ReducedElements")
        ret = cells;
    else if (name == "FaceElements" || name == "ReducedFaceElements")
        ret = faces;
    else if (name == "ContactElements" || name == "ReducedContactElements")
        ret = contacts;

    return ret;
}

//
//
//
void FinleyMesh::reorderGhostZones(int ownIndex)
{
    if (initialized) {
        cells->reorderGhostZones(ownIndex);
        faces->reorderGhostZones(ownIndex);
        contacts->reorderGhostZones(ownIndex);
#ifdef _DEBUG
        cout << "block " << ownIndex << " has " << cells->getGhostCount()
             << " (reduced=" << cells->getReducedGhostCount()
             << ") ghost zones," << endl;
        cout << "\t" << faces->getGhostCount() << " ("
             << faces->getReducedGhostCount() << ") ghost faces," << endl;
        cout << "\t" << contacts->getGhostCount() << " ("
             << contacts->getReducedGhostCount() << ") ghost contacts." << endl;
#endif
    }
}

//
//
//
void FinleyMesh::removeGhostZones(int ownIndex)
{
    if (initialized) {
        cells->removeGhostZones(ownIndex);
        faces->removeGhostZones(ownIndex);
        contacts->removeGhostZones(ownIndex);
#ifdef _DEBUG
        cout << "After removing ghost zones there are" << endl;
        cout << "    " << nodes->getNumNodes() << " Nodes, ";
        cout << cells->getCount() << " Elements, ";
        cout << faces->getCount() << " Face elements, ";
        cout << contacts->getCount() << " Contact elements left." << endl;
#endif
    }
}

//
//
//
bool FinleyMesh::writeToSilo(DBfile* dbfile, const string& pathInSilo)
{
#if USE_SILO
    // Write nodes, elements and mesh variables
    if (!initialized ||
            !cells->writeToSilo(dbfile, pathInSilo) ||
            !faces->writeToSilo(dbfile, pathInSilo) ||
            !contacts->writeToSilo(dbfile, pathInSilo))
        return false;

    siloPath = pathInSilo;
    return true;

#else // !USE_SILO
    return false;
#endif
}

} // namespace escriptexport

