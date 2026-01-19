
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

#include <weipa/VisItData.h>
#include <weipa/DataVar.h>
#include <weipa/ElementData.h>
#include <weipa/NodeData.h>

#include <VisItDataInterface_V2.h>

#include <boost/python/dict.hpp>
#include <boost/python/extract.hpp>

#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <sstream>

using std::map;
using std::set;
using std::string;
using std::vector;

namespace weipa {

//
// Returns simulation metadata from this dataset.
//
visit_handle VisItData::getSimMetaData()
{
    visit_handle mdh = VISIT_INVALID_HANDLE;
    VisIt_SimulationMetaData_alloc(&mdh);
    vector<string>::const_iterator it;
    for (it=cmdNames.begin(); it!=cmdNames.end(); it++) {
        visit_handle cmd = VISIT_INVALID_HANDLE;
        if (VisIt_CommandMetaData_alloc(&cmd) == VISIT_OKAY) {
            VisIt_CommandMetaData_setName(cmd, it->c_str());
            VisIt_SimulationMetaData_addGenericCommand(mdh, cmd);
        }
    }
    VisIt_SimulationMetaData_setMode(mdh, runFlag ?
        VISIT_SIMMODE_RUNNING : VISIT_SIMMODE_STOPPED);

    if (dataset.get() == NULL) {
        return mdh;
    }

    VisIt_SimulationMetaData_setCycleTime(
            mdh, dataset->getCycle(), dataset->getTime());

    set<string> usedMeshes;

    // add "special" mesh variable metadata
    const VarVector& meshVars = dataset->getMeshVariables();
    VarVector::const_iterator vIt;
    for (vIt = meshVars.begin(); vIt != meshVars.end(); vIt++) {
        // skip empty variable
        int numSamples = accumulate(vIt->sampleDistribution.begin(),
                vIt->sampleDistribution.end(), 0);
        if (numSamples > 0 && vIt->valid) {
            DataVar_ptr readerVar = vIt->dataChunks[0];
            string meshName = readerVar->getMeshName();
            usedMeshes.insert(meshName);
            int centering = readerVar->isNodeCentered() ?
                VISIT_VARCENTERING_NODE : VISIT_VARCENTERING_ZONE;
            string varPath = "mesh_vars/"+vIt->varName;
            variables[varPath] = readerVar;
            addVariableMetadata(mdh, varPath, meshName, centering, 0);
        }
    }

    // add other variables
    const VarVector& vars = dataset->getVariables();
    for (vIt = vars.begin(); vIt != vars.end(); vIt++) {
        // skip empty variable
        int numSamples = accumulate(vIt->sampleDistribution.begin(),
                vIt->sampleDistribution.end(), 0);
        if (numSamples == 0 || !vIt->valid) {
            continue;
        }

        const string& varName = vIt->varName;
        DataVar_ptr readerVar = vIt->dataChunks[0];
        string meshName = readerVar->getMeshName();
        usedMeshes.insert(meshName);
        int centering = readerVar->isNodeCentered() ?
            VISIT_VARCENTERING_NODE : VISIT_VARCENTERING_ZONE;
        int rank = readerVar->getRank();
        addVariableMetadata(mdh, varName, meshName, centering, rank);
        variables[varName] = readerVar;
    }

    // add all meshes
    int mpiSize=1;
#ifdef WEIPA_HAVE_MPI
    MPI_Comm comm = dataset->getMPIComm();
    MPI_Comm_size(comm, &mpiSize);
#endif
    set<string>::const_iterator sIt;
    int dim = dataset->getConvertedDomain()[0]->getNodes()->getNumDims();
    for (sIt = usedMeshes.begin(); sIt != usedMeshes.end(); sIt++) {
        addMeshMetadata(mdh, *sIt, dim, mpiSize);
    }

    return mdh;
}

///
/// Returns the domain list
///
visit_handle VisItData::getDomainList()
{
    // each MPI rank serves exactly one domain (chunk) of the data
    visit_handle domainList = VISIT_INVALID_HANDLE;
    if (VisIt_DomainList_alloc(&domainList) == VISIT_OKAY) {
        int mpiRank=0, mpiSize=1;
#ifdef WEIPA_HAVE_MPI
        MPI_Comm comm = dataset->getMPIComm();
        MPI_Comm_rank(comm, &mpiRank);
        MPI_Comm_size(comm, &mpiSize);
#endif
        visit_handle hdl;
        int* rank = (int*)malloc(sizeof(int));
        *rank = mpiRank;
        VisIt_VariableData_alloc(&hdl);
        VisIt_VariableData_setDataI(hdl, VISIT_OWNER_VISIT, 1, 1, rank);
        VisIt_DomainList_setDomains(domainList, mpiSize, hdl);
    }
    return domainList;
}

///
/// Returns mesh data for the specified mesh
///
visit_handle VisItData::getMesh(const char* name)
{
    visit_handle hmesh = VISIT_INVALID_HANDLE;
    DomainChunk_ptr dom = dataset->getConvertedDomain()[0];
    NodeData_ptr nodes = dom->getMeshByName(name);
    ElementData_ptr elements = dom->getElementsByName(name);

    if (nodes && elements
            && VisIt_UnstructuredMesh_alloc(&hmesh)) {

        // set node coordinates
        visit_handle hx, hy;
        VisIt_VariableData_alloc(&hx);
        VisIt_VariableData_alloc(&hy);
        VisIt_VariableData_setDataF(hx, VISIT_OWNER_SIM, 1,
                nodes->getNumNodes(), nodes->getCoords()[0]);
        VisIt_VariableData_setDataF(hy, VISIT_OWNER_SIM, 1,
                nodes->getNumNodes(), nodes->getCoords()[1]);

        if (nodes->getNumDims() > 2) {
            visit_handle hz;
            VisIt_VariableData_alloc(&hz);
            VisIt_VariableData_setDataF(hz, VISIT_OWNER_SIM, 1,
                    nodes->getNumNodes(), nodes->getCoords()[2]);
            VisIt_UnstructuredMesh_setCoordsXYZ(hmesh, hx, hy, hz);
        } else {
            VisIt_UnstructuredMesh_setCoordsXY(hmesh, hx, hy);
        }

        // set connectivity
        size_t sz = (elements->getNodesPerElement()+1) *
                     elements->getNumElements() * sizeof(int);
        int* conn = (int*)malloc(sz);
        int* cptr = conn;
        const int* nodeList = &elements->getNodeList()[0];
        int type = 0;
        memset(conn, 0, sz);
        switch (elements->getType()) {
            case weipa::ZONETYPE_BEAM: type = VISIT_CELL_BEAM; break;
            case weipa::ZONETYPE_HEX: type = VISIT_CELL_HEX; break;
            case weipa::ZONETYPE_POLYGON: type = VISIT_CELL_POINT; break;
            case weipa::ZONETYPE_QUAD: type = VISIT_CELL_QUAD; break;
            case weipa::ZONETYPE_TET: type = VISIT_CELL_TET; break;
            case weipa::ZONETYPE_TRIANGLE: type = VISIT_CELL_TRI; break;
            default:
                std::cerr << "Unhandled zone type " << elements->getType()
                    << std::endl;
                break;
        }
        for (size_t i=0; i<elements->getNumElements(); i++) {
            *cptr++ = type;
            for (size_t j=0; j<elements->getNodesPerElement(); j++) {
                *cptr++ = *nodeList++;
            }
        }

        visit_handle hc;
        VisIt_VariableData_alloc(&hc);
        VisIt_VariableData_setDataI(hc, VISIT_OWNER_VISIT, 1,
                sz / sizeof(int), conn);
        VisIt_UnstructuredMesh_setConnectivity(hmesh,
                elements->getNumElements(), hc);

        // set ghost information
        VisIt_UnstructuredMesh_setRealIndices(hmesh, 0,
            elements->getNumElements()-elements->getGhostCount()-1);
    }

    return hmesh;
}

///
/// Returns variable data for the specified variable
///
visit_handle VisItData::getVariable(const char* name)
{
    visit_handle hvar = VISIT_INVALID_HANDLE;
    map<string, DataVar_ptr>::const_iterator it;
    it = variables.find(name);
    if (it != variables.end()) {
        DataVar_ptr var = it->second;
        VisIt_VariableData_alloc(&hvar);
        VisIt_VariableData_setDataF(hvar, VISIT_OWNER_VISIT,
            var->getNumberOfComponents(),
            var->getNumberOfSamples(),
            var->getDataFlat());
    }

    return hvar;
}

///
/// Adds an expression to the simulation metadata smd
///
void VisItData::addExpressionMetadata(visit_handle smd, const string& name,
                                      const string& def, int type)
{
    visit_handle emd;
    if (VisIt_ExpressionMetaData_alloc(&emd) == VISIT_OKAY) {
        VisIt_ExpressionMetaData_setName(emd, name.c_str());
        VisIt_ExpressionMetaData_setDefinition(emd, def.c_str());
        VisIt_ExpressionMetaData_setType(emd, type);
        VisIt_SimulationMetaData_addExpression(smd, emd);
    }
}

///
/// Adds a mesh to the simulation metadata smd
///
void VisItData::addMeshMetadata(visit_handle smd, const string& name, int dim,
                                int numDoms)
{
    visit_handle mmd;
    if (VisIt_MeshMetaData_alloc(&mmd) == VISIT_OKAY) {
        VisIt_MeshMetaData_setName(mmd, name.c_str());
        VisIt_MeshMetaData_setMeshType(mmd, VISIT_MESHTYPE_UNSTRUCTURED);
        VisIt_MeshMetaData_setTopologicalDimension(mmd, dim);
        VisIt_MeshMetaData_setSpatialDimension(mmd, dim);
        VisIt_MeshMetaData_setNumDomains(mmd, numDoms);
        VisIt_MeshMetaData_setDomainTitle(mmd, "domains");
        VisIt_SimulationMetaData_addMesh(smd, mmd);
    }
}

///
/// Adds a variable to the simulation metadata smd
///
void VisItData::addVariableMetadata(visit_handle smd, const string& name,
                                    const string& meshName, int centering,
                                    int rank)
{
    visit_handle vmd;
    if (VisIt_VariableMetaData_alloc(&vmd) == VISIT_OKAY) {
        VisIt_VariableMetaData_setName(vmd, name.c_str());
        VisIt_VariableMetaData_setMeshName(vmd, meshName.c_str());
        if (rank == 0) {
            VisIt_VariableMetaData_setType(vmd, VISIT_VARTYPE_SCALAR);
        } else if (rank == 1) {
            VisIt_VariableMetaData_setType(vmd, VISIT_VARTYPE_VECTOR);
        } else {
            VisIt_VariableMetaData_setType(vmd, VISIT_VARTYPE_TENSOR);
        }
        VisIt_VariableMetaData_setCentering(vmd, centering);
        VisIt_SimulationMetaData_addVariable(smd, vmd);
    }
}

} // namespace weipa

