
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

#include <weipa/FinleyDomain.h>
#include <weipa/FinleyNodes.h>
#include <weipa/DataVar.h>

#ifndef VISIT_PLUGIN
#include <dudley/CppAdapter/MeshAdapter.h>
#include <finley/CppAdapter/MeshAdapter.h>
#include <ripley/RipleyDomain.h>
extern "C" {
#include <dudley/Mesh.h>
#include <finley/Mesh.h>
}
#endif

#include <iostream>

#if USE_NETCDF
#include <netcdfcpp.h>
#endif

#if USE_SILO
#include <silo.h>
#endif

using namespace std;

namespace weipa {

//
//
//
FinleyDomain::FinleyDomain() :
    initialized(false)
{
}

//
//
//
FinleyDomain::FinleyDomain(const FinleyDomain& m) :
    boost::enable_shared_from_this<FinleyDomain>()
{
    nodes = FinleyNodes_ptr(new FinleyNodes(*m.nodes));
    cells = FinleyElements_ptr(new FinleyElements(*m.cells));
    faces = FinleyElements_ptr(new FinleyElements(*m.faces));
    contacts = FinleyElements_ptr(new FinleyElements(*m.contacts));
    initialized = m.initialized;
}

//
//
//
FinleyDomain::~FinleyDomain()
{
    cleanup();
}

//
//
//
void FinleyDomain::cleanup()
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
bool FinleyDomain::initFromEscript(const escript::AbstractDomain* escriptDomain)
{
#ifndef VISIT_PLUGIN
    cleanup();

    if (dynamic_cast<const finley::MeshAdapter*>(escriptDomain)) {
        const Finley_Mesh* finleyMesh =
            dynamic_cast<const finley::MeshAdapter*>(escriptDomain)
                ->getFinley_Mesh();

        nodes = FinleyNodes_ptr(new FinleyNodes("Elements"));
        cells = FinleyElements_ptr(new FinleyElements("Elements", nodes));
        faces = FinleyElements_ptr(new FinleyElements("FaceElements", nodes));
        contacts = FinleyElements_ptr(new FinleyElements("ContactElements", nodes));

        if (nodes->initFromFinley(finleyMesh->Nodes) &&
                cells->initFromFinley(finleyMesh->Elements) &&
                faces->initFromFinley(finleyMesh->FaceElements) &&
                contacts->initFromFinley(finleyMesh->ContactElements)) {
            initialized = true;
        }
    } else if (dynamic_cast<const dudley::MeshAdapter*>(escriptDomain)) {
        const Dudley_Mesh* dudleyMesh =
            dynamic_cast<const dudley::MeshAdapter*>(escriptDomain)
                ->getDudley_Mesh();

        nodes = FinleyNodes_ptr(new FinleyNodes("Elements"));
        cells = FinleyElements_ptr(new FinleyElements("Elements", nodes));
        faces = FinleyElements_ptr(new FinleyElements("FaceElements", nodes));
        contacts = FinleyElements_ptr(new FinleyElements("ContactElements", nodes));

        if (nodes->initFromDudley(dudleyMesh->Nodes) &&
                cells->initFromDudley(dudleyMesh->Elements) &&
                faces->initFromDudley(dudleyMesh->FaceElements)) {
            initialized = true;
        }
    } else if (dynamic_cast<const ripley::RipleyDomain*>(escriptDomain)) {
        const ripley::RipleyDomain* ripleyMesh =
            dynamic_cast<const ripley::RipleyDomain*>(escriptDomain);

        nodes = FinleyNodes_ptr(new FinleyNodes("Elements"));
        cells = FinleyElements_ptr(new FinleyElements("Elements", nodes));
        faces = FinleyElements_ptr(new FinleyElements("FaceElements", nodes));
        contacts = FinleyElements_ptr(new FinleyElements("ContactElements", nodes));

        if (nodes->initFromRipley(ripleyMesh->getNodes().get()) &&
                cells->initFromRipley(ripleyMesh->getElements()) &&
                faces->initFromRipley(ripleyMesh->getFaceElements())) {
            initialized = true;
        }
    }

    return initialized;
#else // VISIT_PLUGIN
    return false;
#endif
}

//
// Reads mesh and element data from NetCDF file with given name
//
bool FinleyDomain::initFromFile(const string& filename)
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

    nodes = FinleyNodes_ptr(new FinleyNodes("Elements"));
    if (!nodes->readFromNc(input))
        return false;

    // Read all element types
    cells = FinleyElements_ptr(new FinleyElements("Elements", nodes));
    cells->readFromNc(input);
    faces = FinleyElements_ptr(new FinleyElements("FaceElements", nodes));
    faces->readFromNc(input);
    contacts = FinleyElements_ptr(new FinleyElements("ContactElements", nodes));
    contacts->readFromNc(input);

    delete input;
    initialized = true;
#endif

    return initialized;
}

Centering FinleyDomain::getCenteringForFunctionSpace(int fsCode) const
{
    return (fsCode==FINLEY_REDUCED_NODES || fsCode==FINLEY_NODES ?
        NODE_CENTERED : ZONE_CENTERED);
}

//
//
//
NodeData_ptr FinleyDomain::getMeshForFunctionSpace(int fsCode) const
{
    NodeData_ptr result;

    if (!initialized)
        return result;

    ElementData_ptr elements = getElementsForFunctionSpace(fsCode);
    if (elements != NULL)
        result = elements->getNodes();
 
    return result;
}

//
//
//
ElementData_ptr FinleyDomain::getElementsForFunctionSpace(int fsCode) const
{
    ElementData_ptr result;

    if (!initialized) {
        return result;
    }

    if (fsCode == FINLEY_NODES) {
        result = cells;
    } else if (fsCode == FINLEY_REDUCED_NODES) {
        result = cells->getReducedElements();
        if (!result)
            result = cells;
    } else {
        switch (fsCode) {
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

            default: {
                cerr << "Unsupported function space type " << fsCode
                    << "!" << endl;
                return result;
            }
        }
        int typeId = static_cast<FinleyElements*>(result.get())
            ->getFinleyTypeId();
        if (typeId != Finley_Line3Macro &&
                typeId != Finley_Rec9 && typeId != Finley_Rec9Macro &&
                typeId != Finley_Hex27 &&
                typeId != Finley_Hex27Macro && typeId != Finley_Tri6 &&
                typeId != Finley_Tri6Macro && typeId != Finley_Tet10 &&
                typeId != Finley_Tet10Macro) {
            if (result->getReducedElements())
                result = result->getReducedElements();
        }
    }

    return result;
}

//
// Returns a vector of strings containing mesh names for this domain
//
StringVec FinleyDomain::getMeshNames() const
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
// Returns a vector of strings containing mesh variable names for this domain
//
StringVec FinleyDomain::getVarNames() const
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
DataVar_ptr FinleyDomain::getDataVarByName(const string& name) const
{
    if (!initialized) {
        throw "Domain not initialized";
    }

    DataVar_ptr var(new DataVar(name));
    if (name.find("ContactElements_") != name.npos) {
        const IntVec& data = contacts->getVarDataByName(name);
        string elementName = name.substr(0, name.find('_'));
        ElementData_ptr elements = getElementsByName(elementName);
        var->initFromMeshData(shared_from_this(), data,
                FINLEY_CONTACT_ELEMENTS_1, ZONE_CENTERED, elements->getNodes(),
                elements->getIDs());
    } else if (name.find("FaceElements_") != name.npos) {
        const IntVec& data =  faces->getVarDataByName(name);
        string elementName = name.substr(0, name.find('_'));
        ElementData_ptr elements = getElementsByName(elementName);
        var->initFromMeshData(shared_from_this(), data,
                FINLEY_FACE_ELEMENTS, ZONE_CENTERED, elements->getNodes(),
                elements->getIDs());
    } else if (name.find("Elements_") != name.npos) {
        const IntVec& data =  cells->getVarDataByName(name);
        string elementName = name.substr(0, name.find('_'));
        ElementData_ptr elements = getElementsByName(elementName);
        var->initFromMeshData(shared_from_this(), data, FINLEY_ELEMENTS,
                ZONE_CENTERED, elements->getNodes(), elements->getIDs());
    } else if (name.find("Nodes_") != name.npos) {
        const IntVec& data =  nodes->getVarDataByName(name);
        var->initFromMeshData(shared_from_this(), data, FINLEY_NODES,
                NODE_CENTERED, getNodes(), getNodes()->getNodeIDs());
    } else {
        cerr << "WARNING: Unrecognized domain variable '" << name << "'\n";
        return DataVar_ptr();
    }

    return var;
}

//
//
//
ElementData_ptr FinleyDomain::getElementsByName(const string& name) const
{
    ElementData_ptr ret;
    if (name == "Elements")
        ret = cells;
    else if (name == "ReducedElements")
        ret = cells->getReducedElements();
    else if (name == "FaceElements")
        ret = faces;
    else if (name == "ReducedFaceElements")
        ret = faces->getReducedElements();
    else if (name == "ContactElements")
        ret = contacts;
    else if (name == "ReducedContactElements")
        ret = contacts->getReducedElements();

    return ret;
}

//
//
//
NodeData_ptr FinleyDomain::getMeshByName(const string& name) const
{
    NodeData_ptr ret;
    if (initialized) {
        ElementData_ptr els = getElementsByName(name);
        if (els)
            ret = els->getNodes();
    }

    return ret;
}

//
//
//
void FinleyDomain::reorderGhostZones(int ownIndex)
{
    if (initialized) {
        cells->reorderGhostZones(ownIndex);
        faces->reorderGhostZones(ownIndex);
        contacts->reorderGhostZones(ownIndex);
#ifdef _DEBUG
        cout << "block " << ownIndex << " has " << cells->getGhostCount()
             << " ghost zones," << endl;
        cout << "\t" << faces->getGhostCount() << " ghost faces," << endl;
        cout << "\t" << contacts->getGhostCount() << " ghost contacts." << endl;
#endif
    }
}

//
//
//
void FinleyDomain::removeGhostZones(int ownIndex)
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
bool FinleyDomain::writeToSilo(DBfile* dbfile, const string& pathInSilo,
                               const StringVec& labels, const StringVec& units,
                               bool writeMeshData)
{
#if USE_SILO
    // Write nodes, elements and mesh variables
    if (!initialized ||
            !cells->writeToSilo(dbfile, pathInSilo, labels, units, writeMeshData) ||
            !faces->writeToSilo(dbfile, pathInSilo, labels, units, writeMeshData) ||
            !contacts->writeToSilo(dbfile, pathInSilo, labels, units, writeMeshData))
        return false;

    siloPath = pathInSilo;
    return true;

#else // !USE_SILO
    return false;
#endif
}

} // namespace weipa

