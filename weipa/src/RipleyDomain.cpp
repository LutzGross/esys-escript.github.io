
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/

#include <weipa/RipleyDomain.h>
#include <weipa/RipleyNodes.h>
#include <weipa/DataVar.h>

#ifndef VISIT_PLUGIN
#include <ripley/RipleyDomain.h>
#endif

#include <iostream>

#ifdef ESYS_HAVE_SILO
#include <silo.h>
#endif

using namespace std;

namespace weipa {

RipleyDomain::RipleyDomain() :
    initialized(false)
{
}

RipleyDomain::RipleyDomain(const RipleyDomain& m) :
    boost::enable_shared_from_this<RipleyDomain>()
{
    nodes = RipleyNodes_ptr(new RipleyNodes(*m.nodes));
    cells = RipleyElements_ptr(new RipleyElements(*m.cells));
    faces = RipleyElements_ptr(new RipleyElements(*m.faces));
    initialized = m.initialized;
}

bool RipleyDomain::initFromEscript(const escript::AbstractDomain* escriptDomain)
{
    initialized = false;
#ifndef VISIT_PLUGIN
    const ripley::RipleyDomain* dom = dynamic_cast<const ripley::RipleyDomain*>(escriptDomain);
    if (dom) {
        nodes = RipleyNodes_ptr(new RipleyNodes("Elements"));
        cells = RipleyElements_ptr(new RipleyElements("Elements", nodes));
        faces = RipleyElements_ptr(new RipleyElements("FaceElements", nodes));

        if (nodes->initFromRipley(dom) &&
                cells->initFromRipley(dom, ripley::Elements) &&
                faces->initFromRipley(dom, ripley::FaceElements)) {
            initialized = true;
        }
    }
#endif // VISIT_PLUGIN

    return initialized;
}

bool RipleyDomain::initFromFile(const string& filename)
{
    // not supported yet
    return false;
}

Centering RipleyDomain::getCenteringForFunctionSpace(int fsCode) const
{
    return (fsCode==ripley::ReducedNodes || fsCode==ripley::Nodes ?
        NODE_CENTERED : ZONE_CENTERED);
}

NodeData_ptr RipleyDomain::getMeshForFunctionSpace(int fsCode) const
{
    NodeData_ptr result;

    if (!initialized)
        return result;

    ElementData_ptr elements = getElementsForFunctionSpace(fsCode);
    if (elements)
        result = elements->getNodes();
 
    return result;
}

ElementData_ptr RipleyDomain::getElementsForFunctionSpace(int fsCode) const
{
    ElementData_ptr result;

    if (!initialized)
        return result;

    switch (fsCode) {
        case ripley::Nodes:
        case ripley::ReducedNodes: // FIXME: reduced
        case ripley::ReducedElements:
        case ripley::Elements:
            result = cells;
            break;

        case ripley::ReducedFaceElements:
        case ripley::FaceElements:
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
StringVec RipleyDomain::getMeshNames() const
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
StringVec RipleyDomain::getVarNames() const
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

DataVar_ptr RipleyDomain::getDataVarByName(const string& name) const
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
                ripley::FaceElements, ZONE_CENTERED, elements->getNodes(),
                elements->getIDs());
    } else if (name.find("Elements_") != name.npos) {
        const IntVec& data =  cells->getVarDataByName(name);
        string elementName = name.substr(0, name.find('_'));
        ElementData_ptr elements = getElementsByName(elementName);
        var->initFromMeshData(shared_from_this(), data, ripley::Elements,
                ZONE_CENTERED, elements->getNodes(), elements->getIDs());
    } else if (name.find("Nodes_") != name.npos) {
        const IntVec& data =  nodes->getVarDataByName(name);
        var->initFromMeshData(shared_from_this(), data, ripley::Nodes,
                NODE_CENTERED, getNodes(), getNodes()->getNodeIDs());
    } else {
        cerr << "WARNING: Unrecognized domain variable '" << name << "'\n";
        return DataVar_ptr();
    }

    return var;
}

ElementData_ptr RipleyDomain::getElementsByName(const string& name) const
{
    ElementData_ptr ret;
    if (name == "Elements")
        ret = cells;
    else if (name == "FaceElements")
        ret = faces;

    return ret;
}

NodeData_ptr RipleyDomain::getMeshByName(const string& name) const
{
    NodeData_ptr ret;
    if (initialized) {
        ElementData_ptr els = getElementsByName(name);
        if (els)
            ret = els->getNodes();
    }

    return ret;
}

void RipleyDomain::reorderGhostZones(int ownIndex)
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

void RipleyDomain::removeGhostZones(int ownIndex)
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

bool RipleyDomain::writeToSilo(DBfile* dbfile, const string& pathInSilo,
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

