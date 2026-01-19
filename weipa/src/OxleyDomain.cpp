
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

#include <oxley/OxleyDomain.h>

#include <weipa/OxleyDomain.h>
#include <weipa/OxleyNodes.h>
#include <weipa/DataVar.h>

#include <iostream>

#ifdef ESYS_HAVE_SILO
#include <silo.h>
#endif

using namespace std;

namespace weipa {

OxleyDomain::OxleyDomain() :
    initialized(false)
{
}

OxleyDomain::OxleyDomain(const OxleyDomain& m) :
    boost::enable_shared_from_this<OxleyDomain>()
{
    nodes = OxleyNodes_ptr(new OxleyNodes(*m.nodes));
    cells = OxleyElements_ptr(new OxleyElements(*m.cells));
    faces = OxleyElements_ptr(new OxleyElements(*m.faces));
    initialized = m.initialized;
}

bool OxleyDomain::initFromEscript(const escript::AbstractDomain* escriptDomain)
{
    initialized = false;
#ifndef VISIT_PLUGIN
     const oxley::OxleyDomain* dom = dynamic_cast<const oxley::OxleyDomain*>(escriptDomain);
     if (dom) {
         nodes = OxleyNodes_ptr(new OxleyNodes("Elements"));
         cells = OxleyElements_ptr(new OxleyElements("Elements", nodes));
         faces = OxleyElements_ptr(new OxleyElements("FaceElements", nodes));

         if (nodes->initFromOxley(dom) &&
                 cells->initFromOxley(dom, oxley::Elements) &&
                 faces->initFromOxley(dom, oxley::FaceElements)) {
             initialized = true;
         }
     }
#endif // VISIT_PLUGIN

    return initialized;
}

bool OxleyDomain::initFromFile(const string& filename)
{
    // not supported yet
    return false;
}

Centering OxleyDomain::getCenteringForFunctionSpace(int fsCode) const
{
    return (fsCode==oxley::ReducedNodes || fsCode==oxley::Nodes ?
        NODE_CENTERED : ZONE_CENTERED);
    // return NODE_CENTERED;
}

NodeData_ptr OxleyDomain::getMeshForFunctionSpace(int fsCode) const
{
    NodeData_ptr result;

    if (!initialized)
        return result;

    ElementData_ptr elements = getElementsForFunctionSpace(fsCode);
    if (elements)
        result = elements->getNodes();
 
    return result;
}

ElementData_ptr OxleyDomain::getElementsForFunctionSpace(int fsCode) const
{
    ElementData_ptr result;

    if (!initialized)
        return result;

    switch (fsCode) {
        case oxley::Nodes:
        case oxley::ReducedNodes: // FIXME: reduced
        case oxley::ReducedElements:
        case oxley::Elements:
            result = cells;
            break;

        case oxley::ReducedFaceElements:
        case oxley::FaceElements:
            result = faces;
            break;

        default: {
            cerr << "Unsupported function space type " << fsCode
                << "!" << endl;
            return result;
        }
    }

    return result;
}

//
// Returns a vector of strings containing mesh names for this domain
//
StringVec OxleyDomain::getMeshNames() const
{
    StringVec res;
    if (initialized) {
        StringVec tmpVec;
        tmpVec = cells->getMeshNames();
        res.insert(res.end(), tmpVec.begin(), tmpVec.end());
        tmpVec = faces->getMeshNames();
        res.insert(res.end(), tmpVec.begin(), tmpVec.end());
    }
    return res;
}

//
// Returns a vector of strings containing mesh variable names for this domain
//
StringVec OxleyDomain::getVarNames() const
{
    StringVec res;
 
    if (initialized) {
        res = nodes->getVarNames();
        StringVec tmpVec = cells->getVarNames();
        res.insert(res.end(), tmpVec.begin(), tmpVec.end());
        tmpVec = faces->getVarNames();
        res.insert(res.end(), tmpVec.begin(), tmpVec.end());
    }

    return res;
}

DataVar_ptr OxleyDomain::getDataVarByName(const string& name) const
{
    if (!initialized) {
        throw "Domain not initialized";
    }

    DataVar_ptr var(new DataVar(name));
    if (name.find("FaceElements_") != name.npos) {
        const IntVec& data =  faces->getVarDataByName(name);
        string elementName = name.substr(0, name.find('_'));
        ElementData_ptr elements = getElementsByName(elementName);
        var->initFromMeshData(shared_from_this(), data,
                oxley::FaceElements, ZONE_CENTERED, elements->getNodes(),
                elements->getIDs());
    } else if (name.find("Elements_") != name.npos) {
        const IntVec& data =  cells->getVarDataByName(name);
        string elementName = name.substr(0, name.find('_'));
        ElementData_ptr elements = getElementsByName(elementName);
        var->initFromMeshData(shared_from_this(), data, oxley::Elements,
                ZONE_CENTERED, elements->getNodes(), elements->getIDs());
    } else if (name.find("Nodes_") != name.npos) {
        const IntVec& data =  nodes->getVarDataByName(name);
        var->initFromMeshData(shared_from_this(), data, oxley::Nodes,
                NODE_CENTERED, getNodes(), getNodes()->getNodeIDs());
    } else {
        cerr << "WARNING: Unrecognized domain variable '" << name << "'\n";
        return DataVar_ptr();
    }

    return var;
}

ElementData_ptr OxleyDomain::getElementsByName(const string& name) const
{
    ElementData_ptr ret;
    if (name == "Elements")
        ret = cells;
    else if (name == "FaceElements")
        ret = faces;

    return ret;
}

NodeData_ptr OxleyDomain::getMeshByName(const string& name) const
{
    NodeData_ptr ret;
    if (initialized) {
        ElementData_ptr els = getElementsByName(name);
        if (els)
            ret = els->getNodes();
    }

    return ret;
}

void OxleyDomain::reorderGhostZones(int ownIndex)
{
    if (initialized) {
        cells->reorderGhostZones(ownIndex);
        faces->reorderGhostZones(ownIndex);
#ifdef _DEBUG
        cout << "block " << ownIndex << " has " << cells->getGhostCount()
             << " ghost zones," << endl;
        cout << "\t" << faces->getGhostCount() << " ghost faces." << endl;
#endif
    }
}

void OxleyDomain::removeGhostZones(int ownIndex)
{
    if (initialized) {
        cells->removeGhostZones(ownIndex);
        faces->removeGhostZones(ownIndex);
#ifdef _DEBUG
        cout << "After removing ghost zones there are" << endl;
        cout << "    " << nodes->getNumNodes() << " Nodes, ";
        cout << cells->getCount() << " Elements, ";
        cout << faces->getCount() << " Face elements left." << endl;
#endif
    }
}

bool OxleyDomain::writeToSilo(DBfile* dbfile, const string& pathInSilo,
                               const StringVec& labels, const StringVec& units,
                               bool writeMeshData)
{
#ifdef ESYS_HAVE_SILO
    // Write nodes, elements and mesh variables
    if (!initialized
            || !cells->writeToSilo(dbfile, pathInSilo, labels, units, writeMeshData)
            || !faces->writeToSilo(dbfile, pathInSilo, labels, units, writeMeshData))
        return false;

    siloPath = pathInSilo;
    return true;

#else // !ESYS_HAVE_SILO
    return false;
#endif
}

} // namespace weipa

