
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* https://github.com/LutzGross/esys-escript.github.io
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/

#include <weipa/SpeckleyDomain.h>
#include <weipa/SpeckleyNodes.h>
#include <weipa/DataVar.h>

#ifndef VISIT_PLUGIN
#include <speckley/SpeckleyDomain.h>
#endif

#include <iostream>

#ifdef ESYS_HAVE_SILO
#include <silo.h>
#endif

using namespace std;

namespace weipa {

SpeckleyDomain::SpeckleyDomain() :
    initialized(false)
{
}

SpeckleyDomain::SpeckleyDomain(const SpeckleyDomain& m) :
    boost::enable_shared_from_this<SpeckleyDomain>()
{
    nodes = SpeckleyNodes_ptr(new SpeckleyNodes(*m.nodes));
    cells = SpeckleyElements_ptr(new SpeckleyElements(*m.cells));
    faces = SpeckleyElements_ptr(new SpeckleyElements(*m.faces));
    initialized = m.initialized;
}

bool SpeckleyDomain::initFromEscript(const escript::AbstractDomain* escriptDomain)
{
    initialized = false;
#ifndef VISIT_PLUGIN
    const speckley::SpeckleyDomain* dom = dynamic_cast<const speckley::SpeckleyDomain*>(escriptDomain);
    if (dom) {
        nodes = SpeckleyNodes_ptr(new SpeckleyNodes("Elements"));
        cells = SpeckleyElements_ptr(new SpeckleyElements("Elements", nodes));
        faces = SpeckleyElements_ptr(new SpeckleyElements("FaceElements", nodes));

        if (nodes->initFromSpeckley(dom) &&
                cells->initFromSpeckley(dom, speckley::Elements)) {
            initialized = true;
        }
    }
#endif // VISIT_PLUGIN
    return initialized;
}

bool SpeckleyDomain::initFromFile(const string& filename)
{
    // not supported yet
    return false;
}

Centering SpeckleyDomain::getCenteringForFunctionSpace(int fsCode) const
{
    return (fsCode==speckley::ReducedNodes || fsCode==speckley::Nodes ?
        NODE_CENTERED : ZONE_CENTERED);
}

NodeData_ptr SpeckleyDomain::getMeshForFunctionSpace(int fsCode) const
{
    NodeData_ptr result;

    if (!initialized) {
        std::cerr << "uninitialised skipping getElementsForFunctionSpace\n";
        return result;
    }

    ElementData_ptr elements = getElementsForFunctionSpace(fsCode);
    if (elements)
        result = elements->getNodes();
 
    return result;
}

ElementData_ptr SpeckleyDomain::getElementsForFunctionSpace(int fsCode) const
{
    ElementData_ptr result;

    if (!initialized) {
        std::cerr << "uninitialised skipping getElementsForFunctionSpace\n";
        return result;
    }
    switch (fsCode) {
        case speckley::Nodes:
//        case speckley::ReducedNodes: // FIXME: reduced
//        case speckley::ReducedElements:
//        case speckley::Elements:
            result = cells;
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
StringVec SpeckleyDomain::getMeshNames() const
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
StringVec SpeckleyDomain::getVarNames() const
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

DataVar_ptr SpeckleyDomain::getDataVarByName(const string& name) const
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
                speckley::FaceElements, ZONE_CENTERED, elements->getNodes(),
                elements->getIDs());
    } else if (name.find("Elements_") != name.npos) {
        const IntVec& data =  cells->getVarDataByName(name);
        string elementName = name.substr(0, name.find('_'));
        ElementData_ptr elements = getElementsByName(elementName);
        var->initFromMeshData(shared_from_this(), data, speckley::Elements,
                ZONE_CENTERED, elements->getNodes(), elements->getIDs());
    } else if (name.find("Nodes_") != name.npos) {
        const IntVec& data =  nodes->getVarDataByName(name);
        var->initFromMeshData(shared_from_this(), data, speckley::Nodes,
                NODE_CENTERED, getNodes(), getNodes()->getNodeIDs());
    } else {
        cerr << "WARNING: Unrecognized domain variable '" << name << "'\n";
        return DataVar_ptr();
    }

    return var;
}

ElementData_ptr SpeckleyDomain::getElementsByName(const string& name) const
{
    ElementData_ptr ret;
    if (name == "Elements")
        ret = cells;
    else if (name == "FaceElements")
        ret = faces;

    return ret;
}

NodeData_ptr SpeckleyDomain::getMeshByName(const string& name) const
{
    NodeData_ptr ret;
    if (initialized) {
        ElementData_ptr els = getElementsByName(name);
        if (els)
            ret = els->getNodes();
    }

    return ret;
}

void SpeckleyDomain::reorderGhostZones(int ownIndex)
{
    return;
}

void SpeckleyDomain::removeGhostZones(int ownIndex)
{
    return;
}

bool SpeckleyDomain::writeToSilo(DBfile* dbfile, const string& pathInSilo,
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
    std::cerr << "skipping writeToSilo, not built with Silo support\n";
    return false;
#endif
}

} // namespace weipa

