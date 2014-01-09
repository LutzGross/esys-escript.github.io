
/*****************************************************************************
*
* Copyright (c) 2003-2013 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development since 2012 by School of Earth Sciences
*
*****************************************************************************/

#include <pasowrap/PasoException.h>
#include <pasowrap/TransportProblemAdapter.h>
#include "MeshAdapter.h"
#include "escript/Data.h"
#include "escript/DataFactory.h"
#include "esysUtils/blocktimer.h"

#include <boost/python/import.hpp>
#ifdef USE_NETCDF
#include <netcdfcpp.h>
#endif

using namespace std;
using namespace paso;

namespace finley {

// define the static constants
MeshAdapter::FunctionSpaceNamesMapType MeshAdapter::m_functionSpaceTypeNames;
const int MeshAdapter::DegreesOfFreedom=FINLEY_DEGREES_OF_FREEDOM;
const int MeshAdapter::ReducedDegreesOfFreedom=FINLEY_REDUCED_DEGREES_OF_FREEDOM;
const int MeshAdapter::Nodes=FINLEY_NODES;
const int MeshAdapter::ReducedNodes=FINLEY_REDUCED_NODES;
const int MeshAdapter::Elements=FINLEY_ELEMENTS;
const int MeshAdapter::ReducedElements=FINLEY_REDUCED_ELEMENTS;
const int MeshAdapter::FaceElements=FINLEY_FACE_ELEMENTS;
const int MeshAdapter::ReducedFaceElements=FINLEY_REDUCED_FACE_ELEMENTS;
const int MeshAdapter::Points=FINLEY_POINTS;
const int MeshAdapter::ContactElementsZero=FINLEY_CONTACT_ELEMENTS_1;
const int MeshAdapter::ReducedContactElementsZero=FINLEY_REDUCED_CONTACT_ELEMENTS_1;
const int MeshAdapter::ContactElementsOne=FINLEY_CONTACT_ELEMENTS_2;
const int MeshAdapter::ReducedContactElementsOne=FINLEY_REDUCED_CONTACT_ELEMENTS_2;

MeshAdapter::MeshAdapter(Mesh* finleyMesh)
{
    setFunctionSpaceTypeNames();
    // need to use a null_deleter as Finley_Mesh_free deletes the pointer
    // for us.
    m_finleyMesh.reset(finleyMesh, null_deleter());
}

//
// The copy constructor should just increment the use count
MeshAdapter::MeshAdapter(const MeshAdapter& in) :
    m_finleyMesh(in.m_finleyMesh)
{
    setFunctionSpaceTypeNames();
}

MeshAdapter::~MeshAdapter()
{
    // I hope the case for the pointer being zero has been taken care of
    if (m_finleyMesh.unique()) {
        delete m_finleyMesh.get();
    }
}

int MeshAdapter::getMPISize() const
{
    return m_finleyMesh.get()->MPIInfo->size;
}

int MeshAdapter::getMPIRank() const
{
    return m_finleyMesh.get()->MPIInfo->rank;
}

void MeshAdapter::MPIBarrier() const
{
#ifdef ESYS_MPI
    MPI_Barrier(m_finleyMesh.get()->MPIInfo->comm);
#endif
}

bool MeshAdapter::onMasterProcessor() const
{
    return m_finleyMesh.get()->MPIInfo->rank == 0;
}

#ifdef ESYS_MPI
MPI_Comm
#else
unsigned int
#endif
MeshAdapter::getMPIComm() const
{
#ifdef ESYS_MPI
    return m_finleyMesh->MPIInfo->comm;
#else
    return 0;
#endif
}

Mesh* MeshAdapter::getFinley_Mesh() const
{
    return m_finleyMesh.get();
}

void MeshAdapter::write(const string& fileName) const
{
    m_finleyMesh.get()->write(fileName);
    checkFinleyError();
}

void MeshAdapter::Print_Mesh_Info(const bool full) const
{
    m_finleyMesh.get()->printInfo(full);
}

void MeshAdapter::dump(const string& fileName) const
{
#ifdef USE_NETCDF
    const NcDim* ncdims[12] = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
    NcVar *ids;
    int *int_ptr;
    Mesh *mesh = m_finleyMesh.get();
    int num_Tags = 0;
    int mpi_size                         = mesh->MPIInfo->size;
    int mpi_rank                         = mesh->MPIInfo->rank;
    int numDim                           = mesh->Nodes->numDim;
    int numNodes                         = mesh->Nodes->numNodes;
    int num_Elements                     = mesh->Elements->numElements;
    int num_FaceElements                 = mesh->FaceElements->numElements;
    int num_ContactElements              = mesh->ContactElements->numElements;
    int num_Points                       = mesh->Points->numElements;
    int num_Elements_numNodes            = mesh->Elements->numNodes;
    int num_FaceElements_numNodes        = mesh->FaceElements->numNodes;
    int num_ContactElements_numNodes     = mesh->ContactElements->numNodes;
#ifdef ESYS_MPI
    MPI_Status status;
#endif

/* Incoming token indicates it's my turn to write */
#ifdef ESYS_MPI
    if (mpi_rank>0) {
        MPI_Recv(&num_Tags, 0, MPI_INT, mpi_rank-1, 81800, mesh->MPIInfo->comm, &status);
    }
#endif

    const std::string newFileName(esysUtils::appendRankToFileName(fileName,
                                                     mpi_size, mpi_rank));

    // Figure out how much storage is required for tags
    num_Tags = mesh->tagMap.size();

    // NetCDF error handler
    NcError err(NcError::verbose_nonfatal);
    // Create the file.
    NcFile dataFile(newFileName.c_str(), NcFile::Replace);
    string msgPrefix("Error in MeshAdapter::dump: NetCDF operation failed - ");
    // check if writing was successful
    if (!dataFile.is_valid())
        throw FinleyAdapterException(msgPrefix+"Open file for output");

    // Define dimensions (num_Elements and dim_Elements are identical,
    // dim_Elements only appears if > 0)
    if (! (ncdims[0] = dataFile.add_dim("numNodes", numNodes)) )
        throw FinleyAdapterException(msgPrefix+"add_dim(numNodes)");
    if (! (ncdims[1] = dataFile.add_dim("numDim", numDim)) )
        throw FinleyAdapterException(msgPrefix+"add_dim(numDim)");
    if (! (ncdims[2] = dataFile.add_dim("mpi_size_plus_1", mpi_size+1)) )
        throw FinleyAdapterException(msgPrefix+"add_dim(mpi_size)");
    if (num_Elements>0)
        if (! (ncdims[3] = dataFile.add_dim("dim_Elements", num_Elements)) )
            throw FinleyAdapterException(msgPrefix+"add_dim(dim_Elements)");
    if (num_FaceElements>0)
        if (! (ncdims[4] = dataFile.add_dim("dim_FaceElements", num_FaceElements)) )
            throw FinleyAdapterException(msgPrefix+"add_dim(dim_FaceElements)");
    if (num_ContactElements>0)
        if (! (ncdims[5] = dataFile.add_dim("dim_ContactElements", num_ContactElements)) )
            throw FinleyAdapterException(msgPrefix+"add_dim(dim_ContactElements)");
    if (num_Points>0)
        if (! (ncdims[6] = dataFile.add_dim("dim_Points", num_Points)) )
            throw FinleyAdapterException(msgPrefix+"add_dim(dim_Points)");
    if (num_Elements>0)
        if (! (ncdims[7] = dataFile.add_dim("dim_Elements_Nodes", num_Elements_numNodes)) )
            throw FinleyAdapterException(msgPrefix+"add_dim(dim_Elements_Nodes)");
    if (num_FaceElements>0)
        if (! (ncdims[8] = dataFile.add_dim("dim_FaceElements_numNodes", num_FaceElements_numNodes)) )
            throw FinleyAdapterException(msgPrefix+"add_dim(dim_FaceElements_numNodes)");
    if (num_ContactElements>0)
        if (! (ncdims[9] = dataFile.add_dim("dim_ContactElements_numNodes", num_ContactElements_numNodes)) )
            throw FinleyAdapterException(msgPrefix+"add_dim(dim_ContactElements_numNodes)");
    if (num_Tags>0)
        if (! (ncdims[10] = dataFile.add_dim("dim_Tags", num_Tags)) )
            throw FinleyAdapterException(msgPrefix+"add_dim(dim_Tags)");

    // Attributes: MPI size, MPI rank, Name, order, reduced_order
    if (!dataFile.add_att("mpi_size", mpi_size))
        throw FinleyAdapterException(msgPrefix+"add_att(mpi_size)");
    if (!dataFile.add_att("mpi_rank", mpi_rank))
        throw FinleyAdapterException(msgPrefix+"add_att(mpi_rank)");
    if (!dataFile.add_att("Name",mesh->m_name.c_str()))
        throw FinleyAdapterException(msgPrefix+"add_att(Name)");
    if (!dataFile.add_att("numDim",numDim))
        throw FinleyAdapterException(msgPrefix+"add_att(order)");
    if (!dataFile.add_att("order",mesh->integrationOrder))
        throw FinleyAdapterException(msgPrefix+"add_att(order)");
    if (!dataFile.add_att("reduced_order",mesh->reducedIntegrationOrder))
        throw FinleyAdapterException(msgPrefix+"add_att(reduced_order)");
    if (!dataFile.add_att("numNodes",numNodes))
        throw FinleyAdapterException(msgPrefix+"add_att(numNodes)");
    if (!dataFile.add_att("num_Elements",num_Elements))
        throw FinleyAdapterException(msgPrefix+"add_att(num_Elements)");
    if (!dataFile.add_att("num_FaceElements",num_FaceElements))
        throw FinleyAdapterException(msgPrefix+"add_att(num_FaceElements)");
    if (!dataFile.add_att("num_ContactElements",num_ContactElements))
        throw FinleyAdapterException(msgPrefix+"add_att(num_ContactElements)");
    if (!dataFile.add_att("num_Points",num_Points))
        throw FinleyAdapterException(msgPrefix+"add_att(num_Points)");
    if (!dataFile.add_att("num_Elements_numNodes",num_Elements_numNodes))
        throw FinleyAdapterException(msgPrefix+"add_att(num_Elements_numNodes)");
    if (!dataFile.add_att("num_FaceElements_numNodes",num_FaceElements_numNodes) )
        throw FinleyAdapterException(msgPrefix+"add_att(num_FaceElements_numNodes)");
    if (!dataFile.add_att("num_ContactElements_numNodes",num_ContactElements_numNodes) )
        throw FinleyAdapterException(msgPrefix+"add_att(num_ContactElements_numNodes)");
    if (!dataFile.add_att("Elements_TypeId", mesh->Elements->referenceElementSet->referenceElement->Type->TypeId) )
        throw FinleyAdapterException(msgPrefix+"add_att(Elements_TypeId)");
    if (!dataFile.add_att("FaceElements_TypeId", mesh->FaceElements->referenceElementSet->referenceElement->Type->TypeId) )
        throw FinleyAdapterException(msgPrefix+"add_att(FaceElements_TypeId)");
    if (!dataFile.add_att("ContactElements_TypeId", mesh->ContactElements->referenceElementSet->referenceElement->Type->TypeId) )
        throw FinleyAdapterException(msgPrefix+"add_att(ContactElements_TypeId)");
    if (!dataFile.add_att("Points_TypeId", mesh->Points->referenceElementSet->referenceElement->Type->TypeId) )
        throw FinleyAdapterException(msgPrefix+"add_att(Points_TypeId)");
    if (!dataFile.add_att("num_Tags", num_Tags) )
        throw FinleyAdapterException(msgPrefix+"add_att(num_Tags)");

    // // // // // Nodes // // // // //

    // Nodes nodeDistribution
    if (! (ids = dataFile.add_var("Nodes_NodeDistribution", ncInt, ncdims[2])) )
        throw FinleyAdapterException(msgPrefix+"add_var(Nodes_NodeDistribution)");
    int_ptr = &mesh->Nodes->nodesDistribution->first_component[0];
    if (! (ids->put(int_ptr, mpi_size+1)) )
        throw FinleyAdapterException(msgPrefix+"put(Nodes_NodeDistribution)");

    // Nodes degreesOfFreedomDistribution
    if (! ( ids = dataFile.add_var("Nodes_DofDistribution", ncInt, ncdims[2])) )
        throw FinleyAdapterException(msgPrefix+"add_var(Nodes_DofDistribution)");
    int_ptr = &mesh->Nodes->degreesOfFreedomDistribution->first_component[0];
    if (! (ids->put(int_ptr, mpi_size+1)) )
        throw FinleyAdapterException(msgPrefix+"put(Nodes_DofDistribution)");

    // Only write nodes if non-empty because NetCDF doesn't like empty arrays
    // (it treats them as NC_UNLIMITED)
    if (numNodes>0) {
        // Nodes Id
        if (! ( ids = dataFile.add_var("Nodes_Id", ncInt, ncdims[0])) )
            throw FinleyAdapterException(msgPrefix+"add_var(Nodes_Id)");
        int_ptr = &mesh->Nodes->Id[0];
        if (! (ids->put(int_ptr, numNodes)) )
            throw FinleyAdapterException(msgPrefix+"put(Nodes_Id)");

        // Nodes Tag
        if (! ( ids = dataFile.add_var("Nodes_Tag", ncInt, ncdims[0])) )
            throw FinleyAdapterException(msgPrefix+"add_var(Nodes_Tag)");
        int_ptr = &mesh->Nodes->Tag[0];
        if (! (ids->put(int_ptr, numNodes)) )
            throw FinleyAdapterException(msgPrefix+"put(Nodes_Tag)");

        // Nodes gDOF
        if (! ( ids = dataFile.add_var("Nodes_gDOF", ncInt, ncdims[0])) )
            throw FinleyAdapterException(msgPrefix+"add_var(Nodes_gDOF)");
        int_ptr = &mesh->Nodes->globalDegreesOfFreedom[0];
        if (! (ids->put(int_ptr, numNodes)) )
            throw FinleyAdapterException(msgPrefix+"put(Nodes_gDOF)");

        // Nodes global node index
        if (! ( ids = dataFile.add_var("Nodes_gNI", ncInt, ncdims[0])) )
            throw FinleyAdapterException(msgPrefix+"add_var(Nodes_gNI)");
        int_ptr = &mesh->Nodes->globalNodesIndex[0];
        if (! (ids->put(int_ptr, numNodes)) )
            throw FinleyAdapterException(msgPrefix+"put(Nodes_gNI)");

        // Nodes grDof
        if (! ( ids = dataFile.add_var("Nodes_grDfI", ncInt, ncdims[0])) )
            throw FinleyAdapterException(msgPrefix+"add_var(Nodes_grDfI)");
        int_ptr = &mesh->Nodes->globalReducedDOFIndex[0];
        if (! (ids->put(int_ptr, numNodes)) )
            throw FinleyAdapterException(msgPrefix+"put(Nodes_grDfI)");

        // Nodes grNI
        if (! ( ids = dataFile.add_var("Nodes_grNI", ncInt, ncdims[0])) )
            throw FinleyAdapterException(msgPrefix+"add_var(Nodes_grNI)");
        int_ptr = &mesh->Nodes->globalReducedNodesIndex[0];
        if (! (ids->put(int_ptr, numNodes)) )
            throw FinleyAdapterException(msgPrefix+"put(Nodes_grNI)");

        // Nodes Coordinates
        if (! ( ids = dataFile.add_var("Nodes_Coordinates", ncDouble, ncdims[0], ncdims[1]) ) )
            throw FinleyAdapterException(msgPrefix+"add_var(Nodes_Coordinates)");
        if (! (ids->put(&(mesh->Nodes->Coordinates[INDEX2(0,0,numDim)]), numNodes, numDim)) )
            throw FinleyAdapterException(msgPrefix+"put(Nodes_Coordinates)");
    }

    // // // // // Elements // // // // //
    if (num_Elements>0) {
        // Elements_Id
        if (! ( ids = dataFile.add_var("Elements_Id", ncInt, ncdims[3])) )
            throw FinleyAdapterException(msgPrefix+"add_var(Elements_Id)");
        int_ptr = &mesh->Elements->Id[0];
        if (! (ids->put(int_ptr, num_Elements)) )
            throw FinleyAdapterException(msgPrefix+"put(Elements_Id)");

        // Elements_Tag
        if (! ( ids = dataFile.add_var("Elements_Tag", ncInt, ncdims[3])) )
            throw FinleyAdapterException(msgPrefix+"add_var(Elements_Tag)");
        int_ptr = &mesh->Elements->Tag[0];
        if (! (ids->put(int_ptr, num_Elements)) )
            throw FinleyAdapterException(msgPrefix+"put(Elements_Tag)");

        // Elements_Owner
        if (! ( ids = dataFile.add_var("Elements_Owner", ncInt, ncdims[3])) )
            throw FinleyAdapterException(msgPrefix+"add_var(Elements_Owner)");
        int_ptr = &mesh->Elements->Owner[0];
        if (! (ids->put(int_ptr, num_Elements)) )
            throw FinleyAdapterException(msgPrefix+"put(Elements_Owner)");

        // Elements_Color
        if (! ( ids = dataFile.add_var("Elements_Color", ncInt, ncdims[3])) )
            throw FinleyAdapterException(msgPrefix+"add_var(Elements_Color)");
        int_ptr = &mesh->Elements->Color[0];
        if (! (ids->put(int_ptr, num_Elements)) )
            throw FinleyAdapterException(msgPrefix+"put(Elements_Color)");

        // Elements_Nodes
        if (! ( ids = dataFile.add_var("Elements_Nodes", ncInt, ncdims[3], ncdims[7]) ) )
            throw FinleyAdapterException(msgPrefix+"add_var(Elements_Nodes)");
        if (! (ids->put(&(mesh->Elements->Nodes[0]), num_Elements, num_Elements_numNodes)) )
            throw FinleyAdapterException(msgPrefix+"put(Elements_Nodes)");
    }

    // // // // // Face_Elements // // // // //
    if (num_FaceElements>0) {
        // FaceElements_Id
        if (! ( ids = dataFile.add_var("FaceElements_Id", ncInt, ncdims[4])) )
            throw FinleyAdapterException(msgPrefix+"add_var(FaceElements_Id)");
        int_ptr = &mesh->FaceElements->Id[0];
        if (! (ids->put(int_ptr, num_FaceElements)) )
            throw FinleyAdapterException(msgPrefix+"put(FaceElements_Id)");

        // FaceElements_Tag
        if (! ( ids = dataFile.add_var("FaceElements_Tag", ncInt, ncdims[4])) )
            throw FinleyAdapterException(msgPrefix+"add_var(FaceElements_Tag)");
        int_ptr = &mesh->FaceElements->Tag[0];
        if (! (ids->put(int_ptr, num_FaceElements)) )
            throw FinleyAdapterException(msgPrefix+"put(FaceElements_Tag)");

        // FaceElements_Owner
        if (! ( ids = dataFile.add_var("FaceElements_Owner", ncInt, ncdims[4])) )
            throw FinleyAdapterException(msgPrefix+"add_var(FaceElements_Owner)");
        int_ptr = &mesh->FaceElements->Owner[0];
        if (! (ids->put(int_ptr, num_FaceElements)) )
            throw FinleyAdapterException(msgPrefix+"put(FaceElements_Owner)");

        // FaceElements_Color
        if (! ( ids = dataFile.add_var("FaceElements_Color", ncInt, ncdims[4])) )
            throw FinleyAdapterException(msgPrefix+"add_var(FaceElements_Color)");
        int_ptr = &mesh->FaceElements->Color[0];
        if (! (ids->put(int_ptr, num_FaceElements)) )
            throw FinleyAdapterException(msgPrefix+"put(FaceElements_Color)");

        // FaceElements_Nodes
        if (! ( ids = dataFile.add_var("FaceElements_Nodes", ncInt, ncdims[4], ncdims[8]) ) )
            throw FinleyAdapterException(msgPrefix+"add_var(FaceElements_Nodes)");
        if (! (ids->put(&(mesh->FaceElements->Nodes[0]), num_FaceElements, num_FaceElements_numNodes)) )
            throw FinleyAdapterException(msgPrefix+"put(FaceElements_Nodes)");
    }

    // // // // // Contact_Elements // // // // //
    if (num_ContactElements>0) {

        // ContactElements_Id
        if (! ( ids = dataFile.add_var("ContactElements_Id", ncInt, ncdims[5])) )
            throw FinleyAdapterException(msgPrefix+"add_var(ContactElements_Id)");
        int_ptr = &mesh->ContactElements->Id[0];
        if (! (ids->put(int_ptr, num_ContactElements)) )
            throw FinleyAdapterException(msgPrefix+"put(ContactElements_Id)");

        // ContactElements_Tag
        if (! ( ids = dataFile.add_var("ContactElements_Tag", ncInt, ncdims[5])) )
            throw FinleyAdapterException(msgPrefix+"add_var(ContactElements_Tag)");
        int_ptr = &mesh->ContactElements->Tag[0];
        if (! (ids->put(int_ptr, num_ContactElements)) )
            throw FinleyAdapterException(msgPrefix+"put(ContactElements_Tag)");

        // ContactElements_Owner
        if (! ( ids = dataFile.add_var("ContactElements_Owner", ncInt, ncdims[5])) )
            throw FinleyAdapterException(msgPrefix+"add_var(ContactElements_Owner)");
        int_ptr = &mesh->ContactElements->Owner[0];
        if (! (ids->put(int_ptr, num_ContactElements)) )
            throw FinleyAdapterException(msgPrefix+"put(ContactElements_Owner)");

        // ContactElements_Color
        if (! ( ids = dataFile.add_var("ContactElements_Color", ncInt, ncdims[5])) )
            throw FinleyAdapterException(msgPrefix+"add_var(ContactElements_Color)");
        int_ptr = &mesh->ContactElements->Color[0];
        if (! (ids->put(int_ptr, num_ContactElements)) )
            throw FinleyAdapterException(msgPrefix+"put(ContactElements_Color)");

        // ContactElements_Nodes
        if (! ( ids = dataFile.add_var("ContactElements_Nodes", ncInt, ncdims[5], ncdims[9]) ) )
            throw FinleyAdapterException(msgPrefix+"add_var(ContactElements_Nodes)");
        if (! (ids->put(&(mesh->ContactElements->Nodes[0]), num_ContactElements, num_ContactElements_numNodes)) )
            throw FinleyAdapterException(msgPrefix+"put(ContactElements_Nodes)");
    }

    // // // // // Points // // // // //
    if (num_Points>0) {
        // Points_Id
        if (! ( ids = dataFile.add_var("Points_Id", ncInt, ncdims[6])) )
            throw FinleyAdapterException(msgPrefix+"add_var(Points_Id)");
        int_ptr = &mesh->Points->Id[0];
        if (! (ids->put(int_ptr, num_Points)) )
            throw FinleyAdapterException(msgPrefix+"put(Points_Id)");

        // Points_Tag
        if (! ( ids = dataFile.add_var("Points_Tag", ncInt, ncdims[6])) )
            throw FinleyAdapterException(msgPrefix+"add_var(Points_Tag)");
        int_ptr = &mesh->Points->Tag[0];
        if (! (ids->put(int_ptr, num_Points)) )
            throw FinleyAdapterException(msgPrefix+"put(Points_Tag)");

        // Points_Owner
        if (! ( ids = dataFile.add_var("Points_Owner", ncInt, ncdims[6])) )
            throw FinleyAdapterException(msgPrefix+"add_var(Points_Owner)");
        int_ptr = &mesh->Points->Owner[0];
        if (! (ids->put(int_ptr, num_Points)) )
            throw FinleyAdapterException(msgPrefix+"put(Points_Owner)");

        // Points_Color
        if (! ( ids = dataFile.add_var("Points_Color", ncInt, ncdims[6])) )
            throw FinleyAdapterException(msgPrefix+"add_var(Points_Color)");
        int_ptr = &mesh->Points->Color[0];
        if (! (ids->put(int_ptr, num_Points)) )
            throw FinleyAdapterException(msgPrefix+"put(Points_Color)");

        // Points_Nodes
        // mesh->Nodes->Id[mesh->Points->Nodes[INDEX2(0,i,1)]]
        if (! ( ids = dataFile.add_var("Points_Nodes", ncInt, ncdims[6]) ) )
            throw FinleyAdapterException(msgPrefix+"add_var(Points_Nodes)");
        if (! (ids->put(&(mesh->Points->Nodes[0]), num_Points)) )
            throw FinleyAdapterException(msgPrefix+"put(Points_Nodes)");
    }

    // // // // // TagMap // // // // //
    if (num_Tags>0) {
        // Temp storage to gather node IDs
        vector<int> Tags_keys;

        // Copy tag data into temp arrays
        TagMap::const_iterator it;
        for (it=mesh->tagMap.begin(); it!=mesh->tagMap.end(); it++) {
            Tags_keys.push_back(it->second);
        }

        // Tags_keys
        if (! (ids = dataFile.add_var("Tags_keys", ncInt, ncdims[10])) )
            throw FinleyAdapterException(msgPrefix+"add_var(Tags_keys)");
        int_ptr = &Tags_keys[0];
        if (! (ids->put(int_ptr, num_Tags)) )
            throw FinleyAdapterException(msgPrefix+"put(Tags_keys)");

        // Tags_names_*
        // This is an array of strings, it should be stored as an array but
        // instead I have hacked in one attribute per string because the NetCDF
        // manual doesn't tell how to do an array of strings
        int i = 0;
        for (it=mesh->tagMap.begin(); it!=mesh->tagMap.end(); it++, i++) {
            stringstream tagnamestream;
            tagnamestream << "Tags_name_" << i;
            const string tagname = tagnamestream.str();
            if (!dataFile.add_att(tagname.c_str(), it->first.c_str()))
                throw FinleyAdapterException(msgPrefix+"add_att(Tags_names_X)");
        }
    }

    // Send token to next MPI process so he can take his turn
#ifdef ESYS_MPI
    if (mpi_rank<mpi_size-1) {
        MPI_Send(&num_Tags, 0, MPI_INT, mpi_rank+1, 81800, mesh->MPIInfo->comm);
    }
#endif

   // NetCDF file is closed by destructor of NcFile object

#else // USE_NETCDF
    setError(IO_ERROR, "MeshAdapter::dump: not configured with netCDF. Please contact your installation manager.");
#endif // USE_NETCDF
    checkFinleyError();
}

string MeshAdapter::getDescription() const
{
    return "FinleyMesh";
}

string MeshAdapter::functionSpaceTypeAsString(int functionSpaceType) const
{
    FunctionSpaceNamesMapType::iterator loc;
    loc=m_functionSpaceTypeNames.find(functionSpaceType);
    if (loc==m_functionSpaceTypeNames.end()) {
        return "Invalid function space type code.";
    } else {
        return loc->second;
    }
}

bool MeshAdapter::isValidFunctionSpaceType(int functionSpaceType) const
{
    FunctionSpaceNamesMapType::iterator loc;
    loc=m_functionSpaceTypeNames.find(functionSpaceType);
    return (loc!=m_functionSpaceTypeNames.end());
}

void MeshAdapter::setFunctionSpaceTypeNames()
{
    m_functionSpaceTypeNames.insert(FunctionSpaceNamesMapType::value_type(
                DegreesOfFreedom,"Finley_DegreesOfFreedom"));
    m_functionSpaceTypeNames.insert(FunctionSpaceNamesMapType::value_type(
                ReducedDegreesOfFreedom,"Finley_ReducedDegreesOfFreedom"));
    m_functionSpaceTypeNames.insert(FunctionSpaceNamesMapType::value_type(
                Nodes,"Finley_Nodes"));
    m_functionSpaceTypeNames.insert(FunctionSpaceNamesMapType::value_type(
                ReducedNodes,"Finley_Reduced_Nodes"));
    m_functionSpaceTypeNames.insert(FunctionSpaceNamesMapType::value_type(
                Elements,"Finley_Elements"));
    m_functionSpaceTypeNames.insert(FunctionSpaceNamesMapType::value_type(
                ReducedElements,"Finley_Reduced_Elements"));
    m_functionSpaceTypeNames.insert(FunctionSpaceNamesMapType::value_type(
                FaceElements,"Finley_Face_Elements"));
    m_functionSpaceTypeNames.insert(FunctionSpaceNamesMapType::value_type(
                ReducedFaceElements,"Finley_Reduced_Face_Elements"));
    m_functionSpaceTypeNames.insert(FunctionSpaceNamesMapType::value_type(
                Points,"Finley_Points"));
    m_functionSpaceTypeNames.insert(FunctionSpaceNamesMapType::value_type(
                ContactElementsZero,"Finley_Contact_Elements_0"));
    m_functionSpaceTypeNames.insert(FunctionSpaceNamesMapType::value_type(
                ReducedContactElementsZero,"Finley_Reduced_Contact_Elements_0"));
    m_functionSpaceTypeNames.insert(FunctionSpaceNamesMapType::value_type(
                ContactElementsOne,"Finley_Contact_Elements_1"));
    m_functionSpaceTypeNames.insert(FunctionSpaceNamesMapType::value_type(
                ReducedContactElementsOne,"Finley_Reduced_Contact_Elements_1"));
}

int MeshAdapter::getContinuousFunctionCode() const
{
    return Nodes;
}

int MeshAdapter::getReducedContinuousFunctionCode() const
{
    return ReducedNodes;
}

int MeshAdapter::getFunctionCode() const
{
    return Elements;
}

int MeshAdapter::getReducedFunctionCode() const
{
    return ReducedElements;
}

int MeshAdapter::getFunctionOnBoundaryCode() const
{
    return FaceElements;
}

int MeshAdapter::getReducedFunctionOnBoundaryCode() const
{
    return ReducedFaceElements;
}

int MeshAdapter::getFunctionOnContactZeroCode() const
{
    return ContactElementsZero;
}

int MeshAdapter::getReducedFunctionOnContactZeroCode() const
{
    return ReducedContactElementsZero;
}

int MeshAdapter::getFunctionOnContactOneCode() const
{
    return ContactElementsOne;
}

int MeshAdapter::getReducedFunctionOnContactOneCode() const
{
    return ReducedContactElementsOne;
}

int MeshAdapter::getSolutionCode() const
{
    return DegreesOfFreedom;
}

int MeshAdapter::getReducedSolutionCode() const
{
    return ReducedDegreesOfFreedom;
}

int MeshAdapter::getDiracDeltaFunctionsCode() const
{
    return Points;
}

//
// return the spatial dimension of the Mesh:
//
int MeshAdapter::getDim() const
{
    Mesh* mesh=m_finleyMesh.get();
    return mesh->getDim();
}

//
// Return the number of data points summed across all MPI processes
//
int MeshAdapter::getNumDataPointsGlobal() const
{
    return m_finleyMesh.get()->Nodes->getGlobalNumNodes();
}

//
// return the number of data points per sample and the number of samples
// needed to represent data on a parts of the mesh.
//
pair<int,int> MeshAdapter::getDataShape(int functionSpaceCode) const
{
    int numDataPointsPerSample=0;
    int numSamples=0;
    Mesh* mesh=m_finleyMesh.get();
    switch (functionSpaceCode) {
        case Nodes:
            numDataPointsPerSample=1;
            numSamples=mesh->Nodes->getNumNodes();
            break;
        case ReducedNodes:
            numDataPointsPerSample=1;
            numSamples=mesh->Nodes->getNumReducedNodes();
            break;
        case Elements:
            if (mesh->Elements!=NULL) {
                numSamples=mesh->Elements->numElements;
                numDataPointsPerSample=mesh->Elements->referenceElementSet->referenceElement->Parametrization->numQuadNodes;
            }
            break;
        case ReducedElements:
            if (mesh->Elements!=NULL) {
                numSamples=mesh->Elements->numElements;
                numDataPointsPerSample=mesh->Elements->referenceElementSet->referenceElementReducedQuadrature->Parametrization->numQuadNodes;
            }
            break;
        case FaceElements:
            if (mesh->FaceElements!=NULL) {
                numDataPointsPerSample=mesh->FaceElements->referenceElementSet->referenceElement->Parametrization->numQuadNodes;
                numSamples=mesh->FaceElements->numElements;
            }
            break;
        case ReducedFaceElements:
            if (mesh->FaceElements!=NULL) {
                numDataPointsPerSample=mesh->FaceElements->referenceElementSet->referenceElementReducedQuadrature->Parametrization->numQuadNodes;
                numSamples=mesh->FaceElements->numElements;
            }
            break;
        case Points:
            if (mesh->Points!=NULL) {
                numDataPointsPerSample=1;
                numSamples=mesh->Points->numElements;
            }
            break;
        case ContactElementsZero:
            if (mesh->ContactElements!=NULL) {
                numDataPointsPerSample=mesh->ContactElements->referenceElementSet->referenceElement->Parametrization->numQuadNodes;
                numSamples=mesh->ContactElements->numElements;
            }
            break;
        case ReducedContactElementsZero:
            if (mesh->ContactElements!=NULL) {
                numDataPointsPerSample=mesh->ContactElements->referenceElementSet->referenceElementReducedQuadrature->Parametrization->numQuadNodes;
                numSamples=mesh->ContactElements->numElements;
            }
            break;
        case ContactElementsOne:
            if (mesh->ContactElements!=NULL) {
                numDataPointsPerSample=mesh->ContactElements->referenceElementSet->referenceElement->Parametrization->numQuadNodes;
                numSamples=mesh->ContactElements->numElements;
            }
            break;
        case ReducedContactElementsOne:
            if (mesh->ContactElements!=NULL) {
                numDataPointsPerSample=mesh->ContactElements->referenceElementSet->referenceElementReducedQuadrature->Parametrization->numQuadNodes;
                numSamples=mesh->ContactElements->numElements;
            }
            break;
        case DegreesOfFreedom:
            if (mesh->Nodes!=NULL) {
                numDataPointsPerSample=1;
                numSamples=mesh->Nodes->getNumDegreesOfFreedom();
            }
            break;
        case ReducedDegreesOfFreedom:
            if (mesh->Nodes!=NULL) {
                numDataPointsPerSample=1;
                numSamples=mesh->Nodes->getNumReducedDegreesOfFreedom();
            }
            break;
        default:
            stringstream temp;
            temp << "Error - Invalid function space type: " << functionSpaceCode << " for domain: " << getDescription();
            throw FinleyAdapterException(temp.str());
            break;
    }
    return pair<int,int>(numDataPointsPerSample,numSamples);
}

//
// adds linear PDE of second order into a given stiffness matrix and right
// hand side:
//
void MeshAdapter::addPDEToSystem(
        escript::AbstractSystemMatrix& mat, escript::Data& rhs,
        const escript::Data& A, const escript::Data& B, const escript::Data& C,
        const escript::Data& D, const escript::Data& X, const escript::Data& Y,
        const escript::Data& d, const escript::Data& y,
        const escript::Data& d_contact, const escript::Data& y_contact,
        const escript::Data& d_dirac, const escript::Data& y_dirac) const
{
    SystemMatrixAdapter* smat=dynamic_cast<SystemMatrixAdapter*>(&mat);
    if (!smat)
        throw FinleyAdapterException("finley only supports Paso system matrices.");

    Mesh* mesh=m_finleyMesh.get();
    Paso_SystemMatrix* S = smat->getPaso_SystemMatrix();
    Assemble_PDE(mesh->Nodes, mesh->Elements, S, rhs, A, B, C, D, X, Y);
    checkFinleyError();


    Assemble_PDE(mesh->Nodes, mesh->FaceElements, S, rhs,
            escript::Data(), escript::Data(), escript::Data(), d,
            escript::Data(), y);
    checkFinleyError();

    Assemble_PDE(mesh->Nodes, mesh->ContactElements, S, rhs,
            escript::Data(), escript::Data(), escript::Data(), d_contact,
            escript::Data(), y_contact);
    checkFinleyError();

    Assemble_PDE(mesh->Nodes, mesh->Points, S, rhs, escript::Data(),
            escript::Data(), escript::Data(), d_dirac, escript::Data(), y_dirac);
    checkFinleyError();
}

void MeshAdapter::addPDEToLumpedSystem(escript::Data& mat,
        const escript::Data& D, const escript::Data& d,
        const escript::Data& d_dirac, const bool useHRZ) const
{
    Mesh* mesh=m_finleyMesh.get();
    Assemble_LumpedSystem(mesh->Nodes, mesh->Elements, mat, D, useHRZ);
    checkFinleyError();

    Assemble_LumpedSystem(mesh->Nodes, mesh->FaceElements, mat, d, useHRZ);
    checkFinleyError();

    Assemble_LumpedSystem(mesh->Nodes, mesh->Points, mat, d_dirac, useHRZ);
    checkFinleyError();
}

//
// adds linear PDE of second order into the right hand side only
//
void MeshAdapter::addPDEToRHS(escript::Data& rhs, const escript::Data& X,
        const escript::Data& Y, const escript::Data& y,
        const escript::Data& y_contact, const escript::Data& y_dirac) const
{
    Mesh* mesh=m_finleyMesh.get();
    Assemble_PDE(mesh->Nodes, mesh->Elements, 0, rhs,
            escript::Data(), escript::Data(), escript::Data(), escript::Data(),
            X, Y);
    checkFinleyError();

    Assemble_PDE(mesh->Nodes, mesh->FaceElements, 0, rhs,
            escript::Data(), escript::Data(), escript::Data(), escript::Data(),
            escript::Data(), y);
    checkFinleyError();

    Assemble_PDE(mesh->Nodes, mesh->ContactElements, 0, rhs,
            escript::Data(), escript::Data(), escript::Data(),
            escript::Data(), escript::Data(), y_contact);
    checkFinleyError();

    Assemble_PDE(mesh->Nodes, mesh->Points, 0, rhs,
            escript::Data(), escript::Data(), escript::Data(), escript::Data(),
            escript::Data(), y_dirac);
    checkFinleyError();
}

//
// adds PDE of second order into a transport problem
//
void MeshAdapter::addPDEToTransportProblem(
        escript::AbstractTransportProblem& tp, escript::Data& source,
        const escript::Data& M, const escript::Data& A, const escript::Data& B,
        const escript::Data& C, const escript::Data& D, const escript::Data& X,
        const escript::Data& Y, const escript::Data& d, const escript::Data& y,
        const escript::Data& d_contact, const escript::Data& y_contact,
        const escript::Data& d_dirac, const escript::Data& y_dirac) const
{
    TransportProblemAdapter* tpa=dynamic_cast<TransportProblemAdapter*>(&tp);
    if (!tpa)
        throw FinleyAdapterException("finley only supports Paso transport problems.");

    source.expand();

    Mesh* mesh=m_finleyMesh.get();
    Paso_TransportProblem* _tp = tpa->getPaso_TransportProblem();

    Assemble_PDE(mesh->Nodes, mesh->Elements, _tp->mass_matrix, source,
                        escript::Data(), escript::Data(), escript::Data(),
                        M, escript::Data(), escript::Data());
    checkFinleyError();

    Assemble_PDE(mesh->Nodes, mesh->Elements, _tp->transport_matrix,
                        source, A, B, C, D, X, Y);
    checkFinleyError();

    Assemble_PDE(mesh->Nodes, mesh->FaceElements, _tp->transport_matrix,
                        source, escript::Data(), escript::Data(),
                        escript::Data(), d, escript::Data(), y);
    checkFinleyError();

    Assemble_PDE(mesh->Nodes, mesh->ContactElements,
                        _tp->transport_matrix, source, escript::Data(),
                        escript::Data(), escript::Data(), d_contact,
                        escript::Data(), y_contact);
    checkFinleyError();

    Assemble_PDE(mesh->Nodes, mesh->Points, _tp->transport_matrix,
                        source, escript::Data(), escript::Data(),
                        escript::Data(), d_dirac, escript::Data(), y_dirac);
    checkFinleyError();
}

//
// interpolates data between different function spaces
//
void MeshAdapter::interpolateOnDomain(escript::Data& target, const escript::Data& in) const
{
    const MeshAdapter& inDomain=dynamic_cast<const MeshAdapter&>(*(in.getFunctionSpace().getDomain()));
    const MeshAdapter& targetDomain=dynamic_cast<const MeshAdapter&>(*(target.getFunctionSpace().getDomain()));
    if (inDomain!=*this)
        throw FinleyAdapterException("Error - Illegal domain of interpolant.");
    if (targetDomain!=*this)
        throw FinleyAdapterException("Error - Illegal domain of interpolation target.");

    Mesh* mesh=m_finleyMesh.get();
    switch(in.getFunctionSpace().getTypeCode()) {
        case Nodes:
            switch(target.getFunctionSpace().getTypeCode()) {
                case Nodes:
                case ReducedNodes:
                case DegreesOfFreedom:
                case ReducedDegreesOfFreedom:
                    Assemble_CopyNodalData(mesh->Nodes, target, in);
                    break;
                case Elements:
                case ReducedElements:
                    Assemble_interpolate(mesh->Nodes, mesh->Elements, in,target);
                    break;
                case FaceElements:
                case ReducedFaceElements:
                    Assemble_interpolate(mesh->Nodes, mesh->FaceElements, in, target);
                    break;
                case Points:
                    Assemble_interpolate(mesh->Nodes, mesh->Points, in, target);
                    break;
                case ContactElementsZero:
                case ReducedContactElementsZero:
                case ContactElementsOne:
                case ReducedContactElementsOne:
                    Assemble_interpolate(mesh->Nodes, mesh->ContactElements, in, target);
                    break;
                default:
                    stringstream temp;
                    temp << "Error - Interpolation on Domain: Finley does not know anything about function space type " << target.getFunctionSpace().getTypeCode();
                    throw FinleyAdapterException(temp.str());
                    break;
            }
            break;
        case ReducedNodes:
            switch(target.getFunctionSpace().getTypeCode()) {
                case Nodes:
                case ReducedNodes:
                case DegreesOfFreedom:
                case ReducedDegreesOfFreedom:
                    Assemble_CopyNodalData(mesh->Nodes, target, in);
                    break;
                case Elements:
                case ReducedElements:
                    Assemble_interpolate(mesh->Nodes, mesh->Elements, in, target);
                    break;
                case FaceElements:
                case ReducedFaceElements:
                    Assemble_interpolate(mesh->Nodes, mesh->FaceElements, in, target);
                    break;
                case Points:
                    Assemble_interpolate(mesh->Nodes, mesh->Points, in, target);
                    break;
                case ContactElementsZero:
                case ReducedContactElementsZero:
                case ContactElementsOne:
                case ReducedContactElementsOne:
                    Assemble_interpolate(mesh->Nodes, mesh->ContactElements, in, target);
                    break;
                default:
                    stringstream temp;
                    temp << "Error - Interpolation on Domain: Finley does not know anything about function space type " << target.getFunctionSpace().getTypeCode();
                    throw FinleyAdapterException(temp.str());
                    break;
            }
            break;
        case Elements:
            if (target.getFunctionSpace().getTypeCode()==Elements) {
                Assemble_CopyElementData(mesh->Elements, target, in);
            } else if (target.getFunctionSpace().getTypeCode()==ReducedElements) {
                Assemble_AverageElementData(mesh->Elements, target, in);
            } else {
                throw FinleyAdapterException("Error - No interpolation with data on elements possible.");
            }
            break;
        case ReducedElements:
            if (target.getFunctionSpace().getTypeCode()==ReducedElements) {
                Assemble_CopyElementData(mesh->Elements, target, in);
            } else {
                throw FinleyAdapterException("Error - No interpolation with data on elements with reduced integration order possible.");
            }
            break;
        case FaceElements:
            if (target.getFunctionSpace().getTypeCode()==FaceElements) {
                Assemble_CopyElementData(mesh->FaceElements, target, in);
            } else if (target.getFunctionSpace().getTypeCode()==ReducedFaceElements) {
                Assemble_AverageElementData(mesh->FaceElements, target, in);
            } else {
                throw FinleyAdapterException("Error - No interpolation with data on face elements possible.");
            }
            break;
        case ReducedFaceElements:
            if (target.getFunctionSpace().getTypeCode()==ReducedFaceElements) {
                Assemble_CopyElementData(mesh->FaceElements, target, in);
            } else {
                throw FinleyAdapterException("Error - No interpolation with data on face elements with reduced integration order possible.");
            }
            break;
        case Points:
            if (target.getFunctionSpace().getTypeCode()==Points) {
                Assemble_CopyElementData(mesh->Points, target, in);
            } else {
                throw FinleyAdapterException("Error - No interpolation with data on points possible.");
            }
            break;
        case ContactElementsZero:
        case ContactElementsOne:
            if (target.getFunctionSpace().getTypeCode()==ContactElementsZero || target.getFunctionSpace().getTypeCode()==ContactElementsOne) {
                Assemble_CopyElementData(mesh->ContactElements, target, in);
            } else if (target.getFunctionSpace().getTypeCode()==ReducedContactElementsZero || target.getFunctionSpace().getTypeCode()==ReducedContactElementsOne) {
                Assemble_AverageElementData(mesh->ContactElements, target, in);
            } else {
                throw FinleyAdapterException("Error - No interpolation with data on contact elements possible.");
            }
            break;
        case ReducedContactElementsZero:
        case ReducedContactElementsOne:
            if (target.getFunctionSpace().getTypeCode()==ReducedContactElementsZero || target.getFunctionSpace().getTypeCode()==ReducedContactElementsOne) {
                Assemble_CopyElementData(mesh->ContactElements, target, in);
            } else {
                throw FinleyAdapterException("Error - No interpolation with data on contact elements with reduced integration order possible.");
            }
            break;
        case DegreesOfFreedom:
            switch(target.getFunctionSpace().getTypeCode()) {
                case ReducedDegreesOfFreedom:
                case DegreesOfFreedom:
                    Assemble_CopyNodalData(mesh->Nodes, target, in);
                    break;

                case Nodes:
                case ReducedNodes:
                    if (getMPISize()>1) {
                        escript::Data in2=escript::Data(in);
                        in2.expand();
                        Assemble_CopyNodalData(mesh->Nodes, target, in2);
                    } else {
                        Assemble_CopyNodalData(mesh->Nodes, target, in);
                    }
                    break;
                case Elements:
                case ReducedElements:
                    if (getMPISize()>1) {
                        escript::Data in2=escript::Data(in, continuousFunction(*this));
                        Assemble_interpolate(mesh->Nodes, mesh->Elements, in2, target);
                    } else {
                        Assemble_interpolate(mesh->Nodes, mesh->Elements, in, target);
                    }
                    break;
                case FaceElements:
                case ReducedFaceElements:
                    if (getMPISize()>1) {
                        escript::Data in2=escript::Data(in, continuousFunction(*this));
                        Assemble_interpolate(mesh->Nodes, mesh->FaceElements, in2, target);
                    } else {
                        Assemble_interpolate(mesh->Nodes, mesh->FaceElements, in, target);
                    }
                    break;
                case Points:
                    if (getMPISize()>1) {
                        //escript::Data in2=escript::Data(in, continuousFunction(*this) );
                    } else {
                        Assemble_interpolate(mesh->Nodes, mesh->Points, in, target);
                    }
                    break;
                case ContactElementsZero:
                case ContactElementsOne:
                case ReducedContactElementsZero:
                case ReducedContactElementsOne:
                    if (getMPISize()>1) {
                        escript::Data in2=escript::Data(in, continuousFunction(*this));
                        Assemble_interpolate(mesh->Nodes, mesh->ContactElements, in2, target);
                    } else {
                        Assemble_interpolate(mesh->Nodes, mesh->ContactElements, in, target);
                    }
                    break;
                default:
                    stringstream temp;
                    temp << "Error - Interpolation On Domain: Finley does not know anything about function space type " << target.getFunctionSpace().getTypeCode();
                    throw FinleyAdapterException(temp.str());
            }
            break;
        case ReducedDegreesOfFreedom:
            switch(target.getFunctionSpace().getTypeCode()) {
                case Nodes:
                    throw FinleyAdapterException("Error - Finley does not support interpolation from reduced degrees of freedom to mesh nodes.");
                case ReducedNodes:
                    if (getMPISize()>1) {
                        escript::Data in2=escript::Data(in);
                        in2.expand();
                        Assemble_CopyNodalData(mesh->Nodes, target, in2);
                    } else {
                        Assemble_CopyNodalData(mesh->Nodes, target, in);
                    }
                    break;
                case DegreesOfFreedom:
                    throw FinleyAdapterException("Error - Finley does not support interpolation from reduced degrees of freedom to degrees of freedom");
                    break;
                case ReducedDegreesOfFreedom:
                    Assemble_CopyNodalData(mesh->Nodes, target, in);
                    break;
                case Elements:
                case ReducedElements:
                    if (getMPISize()>1) {
                        escript::Data in2=escript::Data(in, reducedContinuousFunction(*this) );
                        Assemble_interpolate(mesh->Nodes, mesh->Elements, in2, target);
                    } else {
                        Assemble_interpolate(mesh->Nodes, mesh->Elements, in, target);
                    }
                    break;
                case FaceElements:
                case ReducedFaceElements:
                    if (getMPISize()>1) {
                        escript::Data in2=escript::Data(in, reducedContinuousFunction(*this) );
                        Assemble_interpolate(mesh->Nodes, mesh->FaceElements, in2, target);
                    } else {
                        Assemble_interpolate(mesh->Nodes, mesh->FaceElements, in, target);
                    }
                    break;
                case Points:
                    if (getMPISize()>1) {
                        escript::Data in2=escript::Data(in, reducedContinuousFunction(*this));
                        Assemble_interpolate(mesh->Nodes, mesh->Points, in2, target);
                    } else {
                        Assemble_interpolate(mesh->Nodes, mesh->Points, in, target);
                    }
                    break;
                case ContactElementsZero:
                case ContactElementsOne:
                case ReducedContactElementsZero:
                case ReducedContactElementsOne:
                    if (getMPISize()>1) {
                        escript::Data in2=escript::Data(in, reducedContinuousFunction(*this));
                        Assemble_interpolate(mesh->Nodes, mesh->ContactElements, in2, target);
                    } else {
                        Assemble_interpolate(mesh->Nodes, mesh->ContactElements, in, target);
                    }
                    break;
                default:
                    stringstream temp;
                    temp << "Error - Interpolation On Domain: Finley does not know anything about function space type " << target.getFunctionSpace().getTypeCode();
                    throw FinleyAdapterException(temp.str());
                    break;
            }
            break;
        default:
            stringstream temp;
            temp << "Error - Interpolation On Domain: Finley does not know anything about function space type %d" << in.getFunctionSpace().getTypeCode();
            throw FinleyAdapterException(temp.str());
            break;
    }
    checkFinleyError();
}

//
// copies the locations of sample points into x
//
void MeshAdapter::setToX(escript::Data& arg) const
{
    const MeshAdapter& argDomain=dynamic_cast<const MeshAdapter&>(*(arg.getFunctionSpace().getDomain()));
    if (argDomain!=*this)
        throw FinleyAdapterException("Error - Illegal domain of data point locations");
    Mesh* mesh=m_finleyMesh.get();
    // in case of appropriate function space we can do the job directly:
    if (arg.getFunctionSpace().getTypeCode()==Nodes) {
        Assemble_NodeCoordinates(mesh->Nodes, arg);
    } else {
        escript::Data tmp_data=Vector(0., continuousFunction(*this), true);
        Assemble_NodeCoordinates(mesh->Nodes, tmp_data);
        // this is then interpolated onto arg:
        interpolateOnDomain(arg, tmp_data);
    }
    checkFinleyError();
}

//
// return the normal vectors at the location of data points as a Data object
//
void MeshAdapter::setToNormal(escript::Data& normal) const
{
    const MeshAdapter& normalDomain=dynamic_cast<const MeshAdapter&>(*(normal.getFunctionSpace().getDomain()));
    if (normalDomain!=*this)
        throw FinleyAdapterException("Error - Illegal domain of normal locations");
    Mesh* mesh=m_finleyMesh.get();
    switch(normal.getFunctionSpace().getTypeCode()) {
        case Nodes:
            throw FinleyAdapterException("Error - Finley does not support surface normal vectors for nodes");
            break;
        case ReducedNodes:
            throw FinleyAdapterException("Error - Finley does not support surface normal vectors for reduced nodes");
            break;
        case Elements:
            throw FinleyAdapterException("Error - Finley does not support surface normal vectors for elements");
            break;
        case ReducedElements:
            throw FinleyAdapterException("Error - Finley does not support surface normal vectors for elements with reduced integration order");
            break;
        case FaceElements:
        case ReducedFaceElements:
            Assemble_getNormal(mesh->Nodes, mesh->FaceElements, normal);
            break;
        case Points:
            throw FinleyAdapterException("Error - Finley does not support surface normal vectors for point elements");
            break;
        case ContactElementsOne:
        case ContactElementsZero:
        case ReducedContactElementsOne:
        case ReducedContactElementsZero:
            Assemble_getNormal(mesh->Nodes, mesh->ContactElements, normal);
            break;
        case DegreesOfFreedom:
            throw FinleyAdapterException("Error - Finley does not support surface normal vectors for degrees of freedom.");
            break;
        case ReducedDegreesOfFreedom:
            throw FinleyAdapterException("Error - Finley does not support surface normal vectors for reduced degrees of freedom.");
            break;
        default:
            stringstream temp;
            temp << "Error - Normal Vectors: Finley does not know anything about function space type " << normal.getFunctionSpace().getTypeCode();
            throw FinleyAdapterException(temp.str());
            break;
    }
    checkFinleyError();
}

//
// interpolates data to other domain
//
void MeshAdapter::interpolateACross(escript::Data& target, const escript::Data& source) const
{
    throw FinleyAdapterException("Error - Finley does not allow interpolation across domains.");
}

//
// calculates the integral of a function defined on arg
//
void MeshAdapter::setToIntegrals(vector<double>& integrals, const escript::Data& arg) const
{
    const MeshAdapter& argDomain=dynamic_cast<const MeshAdapter&>(*(arg.getFunctionSpace().getDomain()));
    if (argDomain!=*this)
        throw FinleyAdapterException("Error - Illegal domain of integration kernel");

    double blocktimer_start = blocktimer_time();
    Mesh* mesh=m_finleyMesh.get();
    switch(arg.getFunctionSpace().getTypeCode()) {
        case Nodes:
            {
                escript::Data temp(arg, escript::function(*this));
                Assemble_integrate(mesh->Nodes, mesh->Elements, temp, &integrals[0]);
            }
            break;
        case ReducedNodes:
            {
                escript::Data temp(arg, escript::function(*this));
                Assemble_integrate(mesh->Nodes, mesh->Elements, temp, &integrals[0]);
            }
            break;
        case Elements:
        case ReducedElements:
            Assemble_integrate(mesh->Nodes, mesh->Elements, arg, &integrals[0]);
            break;
        case FaceElements:
        case ReducedFaceElements:
            Assemble_integrate(mesh->Nodes, mesh->FaceElements, arg, &integrals[0]);
            break;
        case Points:
            throw FinleyAdapterException("Error - Integral of data on points is not supported.");
            break;
        case ContactElementsZero:
        case ReducedContactElementsZero:
        case ContactElementsOne:
        case ReducedContactElementsOne:
            Assemble_integrate(mesh->Nodes, mesh->ContactElements, arg, &integrals[0]);
            break;
        case DegreesOfFreedom:
        case ReducedDegreesOfFreedom:
            {
                escript::Data temp(arg, escript::function(*this));
                Assemble_integrate(mesh->Nodes, mesh->Elements, temp, &integrals[0]);
            }
            break;
        default:
            stringstream temp;
            temp << "Error - Integrals: Finley does not know anything about function space type " << arg.getFunctionSpace().getTypeCode();
            throw FinleyAdapterException(temp.str());
            break;
    }
    checkFinleyError();
    blocktimer_increment("integrate()", blocktimer_start);
}

//
// calculates the gradient of arg
//
void MeshAdapter::setToGradient(escript::Data& grad, const escript::Data& arg) const
{
    const MeshAdapter& argDomain=dynamic_cast<const MeshAdapter&>(*(arg.getFunctionSpace().getDomain()));
    if (argDomain!=*this)
        throw FinleyAdapterException("Error - Illegal domain of gradient argument");
    const MeshAdapter& gradDomain=dynamic_cast<const MeshAdapter&>(*(grad.getFunctionSpace().getDomain()));
    if (gradDomain!=*this)
        throw FinleyAdapterException("Error - Illegal domain of gradient");

    Mesh* mesh=m_finleyMesh.get();
    escript::Data nodeData;
    if (getMPISize()>1) {
        if (arg.getFunctionSpace().getTypeCode() == DegreesOfFreedom) {
            nodeData=escript::Data(arg, continuousFunction(*this));
        } else if(arg.getFunctionSpace().getTypeCode() == ReducedDegreesOfFreedom) {
            nodeData=escript::Data(arg, reducedContinuousFunction(*this));
        } else {
            nodeData = arg;
        }
    } else {
        nodeData = arg;
    }
    switch(grad.getFunctionSpace().getTypeCode()) {
        case Nodes:
            throw FinleyAdapterException("Error - Gradient at nodes is not supported.");
            break;
        case ReducedNodes:
            throw FinleyAdapterException("Error - Gradient at reduced nodes is not supported.");
            break;
        case Elements:
        case ReducedElements:
            Assemble_gradient(mesh->Nodes, mesh->Elements, grad, nodeData);
            break;
        case FaceElements:
        case ReducedFaceElements:
            Assemble_gradient(mesh->Nodes, mesh->FaceElements, grad, nodeData);
            break;
        case Points:
            throw FinleyAdapterException("Error - Gradient at points is not supported.");
            break;
        case ContactElementsZero:
        case ReducedContactElementsZero:
        case ContactElementsOne:
        case ReducedContactElementsOne:
            Assemble_gradient(mesh->Nodes, mesh->ContactElements, grad, nodeData);
            break;
        case DegreesOfFreedom:
            throw FinleyAdapterException("Error - Gradient at degrees of freedom is not supported.");
            break;
        case ReducedDegreesOfFreedom:
            throw FinleyAdapterException("Error - Gradient at reduced degrees of freedom is not supported.");
            break;
        default:
            stringstream temp;
            temp << "Error - Gradient: Finley does not know anything about function space type " << arg.getFunctionSpace().getTypeCode();
            throw FinleyAdapterException(temp.str());
            break;
    }
    checkFinleyError();
}

//
// returns the size of elements
//
void MeshAdapter::setToSize(escript::Data& size) const
{
    Mesh* mesh=m_finleyMesh.get();
    switch(size.getFunctionSpace().getTypeCode()) {
        case Nodes:
            throw FinleyAdapterException("Error - Size of nodes is not supported.");
            break;
        case ReducedNodes:
            throw FinleyAdapterException("Error - Size of reduced nodes is not supported.");
            break;
        case Elements:
        case ReducedElements:
            Assemble_getSize(mesh->Nodes, mesh->Elements, size);
            break;
        case FaceElements:
        case ReducedFaceElements:
            Assemble_getSize(mesh->Nodes, mesh->FaceElements, size);
            break;
        case Points:
            throw FinleyAdapterException("Error - Size of point elements is not supported.");
            break;
        case ContactElementsZero:
        case ContactElementsOne:
        case ReducedContactElementsZero:
        case ReducedContactElementsOne:
            Assemble_getSize(mesh->Nodes,mesh->ContactElements,size);
            break;
        case DegreesOfFreedom:
            throw FinleyAdapterException("Error - Size of degrees of freedom is not supported.");
            break;
        case ReducedDegreesOfFreedom:
            throw FinleyAdapterException("Error - Size of reduced degrees of freedom is not supported.");
            break;
        default:
            stringstream temp;
            temp << "Error - Element size: Finley does not know anything about function space type " << size.getFunctionSpace().getTypeCode();
            throw FinleyAdapterException(temp.str());
            break;
    }
    checkFinleyError();
}

//
// sets the location of nodes
//
void MeshAdapter::setNewX(const escript::Data& new_x)
{
    Mesh* mesh=m_finleyMesh.get();
    const MeshAdapter& newDomain=dynamic_cast<const MeshAdapter&>(*(new_x.getFunctionSpace().getDomain()));
    if (newDomain!=*this)
        throw FinleyAdapterException("Error - Illegal domain of new point locations");
    if (new_x.getFunctionSpace() == continuousFunction(*this)) {
        mesh->setCoordinates(new_x);
    } else {
        throw FinleyAdapterException("As of escript version 3.3 SetX() only accepts ContinuousFunction arguments. Please interpolate.");
    }
    checkFinleyError();
}

bool MeshAdapter::ownSample(int fs_code, index_t id) const
{
    if (getMPISize() > 1) {
#ifdef ESYS_MPI
        int myFirstNode=0, myLastNode=0, k=0;
        int* globalNodeIndex=0;
        Mesh* mesh_p=m_finleyMesh.get();
        /*
         * this method is only used by saveDataCSV which would use the returned
         * values for reduced nodes wrongly so this case is disabled for now
        if (fs_code == FINLEY_REDUCED_NODES) {
            myFirstNode = NodeFile_getFirstReducedNode(mesh_p->Nodes);
            myLastNode = NodeFile_getLastReducedNode(mesh_p->Nodes);
            globalNodeIndex = NodeFile_borrowGlobalReducedNodesIndex(mesh_p->Nodes);
        } else
        */
        if (fs_code == FINLEY_NODES) {
            myFirstNode = mesh_p->Nodes->getFirstNode();
            myLastNode = mesh_p->Nodes->getLastNode();
            globalNodeIndex = mesh_p->Nodes->borrowGlobalNodesIndex();
        } else {
            throw FinleyAdapterException("Unsupported function space type for ownSample()");
        }

        k=globalNodeIndex[id];
        return ((myFirstNode <= k) && (k < myLastNode));
#endif
    }
    return true;
}


//
// creates a SystemMatrixAdapter stiffness matrix an initializes it with zeros
//
escript::ASM_ptr MeshAdapter::newSystemMatrix(const int row_blocksize,
                            const escript::FunctionSpace& row_functionspace,
                            const int column_blocksize,
                            const escript::FunctionSpace& column_functionspace,
                            const int type) const
{
    // is the domain right?
    const MeshAdapter& row_domain=dynamic_cast<const MeshAdapter&>(*(row_functionspace.getDomain()));
    if (row_domain!=*this)
        throw FinleyAdapterException("Error - domain of row function space does not match the domain of matrix generator.");
    const MeshAdapter& col_domain=dynamic_cast<const MeshAdapter&>(*(column_functionspace.getDomain()));
    if (col_domain!=*this)
        throw FinleyAdapterException("Error - domain of column function space does not match the domain of matrix generator.");

    int reduceRowOrder=0;
    int reduceColOrder=0;
    // is the function space type right?
    if (row_functionspace.getTypeCode() == ReducedDegreesOfFreedom) {
        reduceRowOrder=1;
    } else if (row_functionspace.getTypeCode() != DegreesOfFreedom) {
        throw FinleyAdapterException("Error - illegal function space type for system matrix rows.");
    }
    if (column_functionspace.getTypeCode() == ReducedDegreesOfFreedom) {
        reduceColOrder=1;
    } else if (column_functionspace.getTypeCode() != DegreesOfFreedom) {
        throw FinleyAdapterException("Error - illegal function space type for system matrix columns.");
    }

    // generate matrix:
    Paso_SystemMatrixPattern* fsystemMatrixPattern=
        getFinley_Mesh()->getPattern(reduceRowOrder, reduceColOrder);
    checkFinleyError();
    Paso_SystemMatrix* fsystemMatrix;
    const int trilinos = 0;
    if (trilinos) {
#ifdef TRILINOS
        // FIXME: Allocation Epetra_VrbMatrix here...
#endif
    } else {
        fsystemMatrix=Paso_SystemMatrix_alloc(type, fsystemMatrixPattern,
                row_blocksize, column_blocksize, FALSE);
    }
    checkPasoError();
    Paso_SystemMatrixPattern_free(fsystemMatrixPattern);
    SystemMatrixAdapter* sma=new SystemMatrixAdapter(fsystemMatrix, row_blocksize, row_functionspace, column_blocksize, column_functionspace);
    return escript::ASM_ptr(sma);
}

//
// creates a TransportProblemAdapter
//
escript::ATP_ptr MeshAdapter::newTransportProblem(const int blocksize,
        const escript::FunctionSpace& functionspace, const int type) const
{
    // is the domain right?
    const MeshAdapter& domain=dynamic_cast<const MeshAdapter&>(*(functionspace.getDomain()));
    if (domain!=*this)
        throw FinleyAdapterException("Error - domain of function space does not match the domain of transport problem generator.");

    // is the function space type right?
    int reduceOrder=0;
    if (functionspace.getTypeCode() == ReducedDegreesOfFreedom) {
        reduceOrder=1;
    } else if (functionspace.getTypeCode() != DegreesOfFreedom) {
        throw FinleyAdapterException("Error - illegal function space type for transport problem.");
    }

    // generate transport problem:
    Paso_SystemMatrixPattern* fsystemMatrixPattern=
        getFinley_Mesh()->getPattern(reduceOrder, reduceOrder);
    checkFinleyError();
    Paso_TransportProblem* transportProblem;
    transportProblem=Paso_TransportProblem_alloc(fsystemMatrixPattern, blocksize);
    checkPasoError();
    Paso_SystemMatrixPattern_free(fsystemMatrixPattern);
    TransportProblemAdapter* tpa=new TransportProblemAdapter(
            transportProblem, blocksize, functionspace);
    return escript::ATP_ptr(tpa);
}

//
// returns true if data on functionSpaceCode is considered as being cell centered
//
bool MeshAdapter::isCellOriented(int functionSpaceCode) const
{
    switch(functionSpaceCode) {
        case Nodes:
        case DegreesOfFreedom:
        case ReducedDegreesOfFreedom:
            return false;
        case Elements:
        case FaceElements:
        case Points:
        case ContactElementsZero:
        case ContactElementsOne:
        case ReducedElements:
        case ReducedFaceElements:
        case ReducedContactElementsZero:
        case ReducedContactElementsOne:
            return true;
        default:
            stringstream temp;
            temp << "Error - Cell: Finley does not know anything about function space type " << functionSpaceCode;
            throw FinleyAdapterException(temp.str());
            break;
    }
    return false;
}

bool
MeshAdapter::commonFunctionSpace(const vector<int>& fs, int& resultcode) const
{
    /* The idea is to use equivalence classes. [Types which can be interpolated
       back and forth]:
        class 1: DOF <-> Nodes
        class 2: ReducedDOF <-> ReducedNodes
        class 3: Points
        class 4: Elements
        class 5: ReducedElements
        class 6: FaceElements
        class 7: ReducedFaceElements
        class 8: ContactElementZero <-> ContactElementOne
        class 9: ReducedContactElementZero <-> ReducedContactElementOne

    There is also a set of lines. Interpolation is possible down a line but
    not between lines.
    class 1 and 2 belong to all lines so aren't considered.
        line 0: class 3
        line 1: class 4,5
        line 2: class 6,7
        line 3: class 8,9

    For classes with multiple members (eg class 2) we have vars to record if
    there is at least one instance.
    e.g. hasnodes is true if we have at least one instance of Nodes.
    */
    if (fs.empty())
        return false;

    vector<int> hasclass(10);
    vector<int> hasline(4);     
    bool hasnodes=false;
    bool hasrednodes=false;
    bool hascez=false;
    bool hasrcez=false;
    for (int i=0;i<fs.size();++i) {
        switch(fs[i]) {
            case Nodes:
                hasnodes=true; // fall through
            case DegreesOfFreedom:
                hasclass[1]=1;
                break;
            case ReducedNodes:
                hasrednodes=true; // fall through
            case ReducedDegreesOfFreedom:
                hasclass[2]=1;
                break;
            case Points:
                hasline[0]=1;
                hasclass[3]=1;
                break;
            case Elements:
                hasclass[4]=1;
                hasline[1]=1;
                break;
            case ReducedElements:
                hasclass[5]=1;
                hasline[1]=1;
                break;
            case FaceElements:
                hasclass[6]=1;
                hasline[2]=1;
                break;
            case ReducedFaceElements:
                hasclass[7]=1;
                hasline[2]=1;
                break;
            case ContactElementsZero:
                hascez=true; // fall through
            case ContactElementsOne:
                hasclass[8]=1;
                hasline[3]=1;
                break;
            case ReducedContactElementsZero:
                hasrcez=true; // fall through
            case ReducedContactElementsOne:
                hasclass[9]=1;
                hasline[3]=1;
                break;
            default:
                return false;
        }
    }
    int totlines=hasline[0]+hasline[1]+hasline[2]+hasline[3];

    // fail if we have more than one leaf group
    if (totlines>1)
        return false; // there are at least two branches we can't interpolate between
    else if (totlines==1) {
        if (hasline[0]==1)              // we have points
            resultcode=Points;
        else if (hasline[1]==1) {
            if (hasclass[5]==1)
                resultcode=ReducedElements;
            else
                resultcode=Elements;
        } else if (hasline[2]==1) {
            if (hasclass[7]==1)
                resultcode=ReducedFaceElements;
            else
                resultcode=FaceElements;
        } else {   // so we must be in line3
            if (hasclass[9]==1) {
                // need something from class 9
                resultcode=(hasrcez ? ReducedContactElementsZero : ReducedContactElementsOne);
            } else {
                // something from class 8
                resultcode=(hascez?ContactElementsZero:ContactElementsOne);
            }
        }
    } else { // totlines==0
        if (hasclass[2]==1) {
            // something from class 2
            resultcode=(hasrednodes ? ReducedNodes : ReducedDegreesOfFreedom);
        } else { 
            // something from class 1
            resultcode=(hasnodes ? Nodes : DegreesOfFreedom);
        }
    }
    return true;
}

bool MeshAdapter::probeInterpolationOnDomain(int functionSpaceType_source,
                                             int functionSpaceType_target) const
{
    switch(functionSpaceType_source) {
        case Nodes:
            switch(functionSpaceType_target) {
                case Nodes:
                case ReducedNodes:
                case ReducedDegreesOfFreedom:
                case DegreesOfFreedom:
                case Elements:
                case ReducedElements:
                case FaceElements:
                case ReducedFaceElements:
                case Points:
                case ContactElementsZero:
                case ReducedContactElementsZero:
                case ContactElementsOne:
                case ReducedContactElementsOne:
                    return true;
                default:
                    stringstream temp;
                    temp << "Error - Interpolation On Domain: Finley does not know anything about function space type " << functionSpaceType_target;
                    throw FinleyAdapterException(temp.str());
            }
            break;
        case ReducedNodes:
            switch(functionSpaceType_target) {
                case ReducedNodes:
                case ReducedDegreesOfFreedom:
                case Elements:
                case ReducedElements:
                case FaceElements:
                case ReducedFaceElements:
                case Points:
                case ContactElementsZero:
                case ReducedContactElementsZero:
                case ContactElementsOne:
                case ReducedContactElementsOne:
                    return true;
                case Nodes:
                case DegreesOfFreedom:
                    return false;
                default:
                    stringstream temp;
                    temp << "Error - Interpolation On Domain: Finley does not know anything about function space type " << functionSpaceType_target;
                    throw FinleyAdapterException(temp.str());
            }
            break;
        case Elements:
            if (functionSpaceType_target==Elements) {
                return true;
            } else if (functionSpaceType_target==ReducedElements) {
                return true;
            } else {
                return false;
            }
        case ReducedElements:
            if (functionSpaceType_target==ReducedElements) {
                return true;
            } else {
                return false;
            }
        case FaceElements:
            if (functionSpaceType_target==FaceElements) {
                return true;
            } else if (functionSpaceType_target==ReducedFaceElements) {
                return true;
            } else {
                return false;
            }
        case ReducedFaceElements:
            if (functionSpaceType_target==ReducedFaceElements) {
                return true;
            } else {
                return false;
            }
        case Points:
            if (functionSpaceType_target==Points) {
                return true;
            } else {
                return false;
            }
        case ContactElementsZero:
        case ContactElementsOne:
            if (functionSpaceType_target==ContactElementsZero || functionSpaceType_target==ContactElementsOne) {
                return true;
            } else if (functionSpaceType_target==ReducedContactElementsZero || functionSpaceType_target==ReducedContactElementsOne) {
                return true;
            } else {
                return false;
            }
        case ReducedContactElementsZero:
        case ReducedContactElementsOne:
            if (functionSpaceType_target==ReducedContactElementsZero || functionSpaceType_target==ReducedContactElementsOne) {
                return true;
            } else {
                return false;
            }
        case DegreesOfFreedom:
            switch(functionSpaceType_target) {
                case ReducedDegreesOfFreedom:
                case DegreesOfFreedom:
                case Nodes:
                case ReducedNodes:
                case Elements:
                case ReducedElements:
                case Points:
                case FaceElements:
                case ReducedFaceElements:
                case ContactElementsZero:
                case ReducedContactElementsZero:
                case ContactElementsOne:
                case ReducedContactElementsOne:
                    return true;
                default:
                    stringstream temp;
                    temp << "Error - Interpolation On Domain: Finley does not know anything about function space type " << functionSpaceType_target;
                    throw FinleyAdapterException(temp.str());
            }
            break;
        case ReducedDegreesOfFreedom:
            switch(functionSpaceType_target) {
                case ReducedDegreesOfFreedom:
                case ReducedNodes:
                case Elements:
                case ReducedElements:
                case FaceElements:
                case ReducedFaceElements:
                case Points:
                case ContactElementsZero:
                case ReducedContactElementsZero:
                case ContactElementsOne:
                case ReducedContactElementsOne:
                    return true;
                case Nodes:
                case DegreesOfFreedom:
                    return false;
                default:
                    stringstream temp;
                    temp << "Error - Interpolation On Domain: Finley does not know anything about function space type " << functionSpaceType_target;
                    throw FinleyAdapterException(temp.str());
            }
            break;
        default:
            stringstream temp;
            temp << "Error - Interpolation On Domain: Finley does not know anything about function space type " << functionSpaceType_source;
            throw FinleyAdapterException(temp.str());
            break;
    }
    return false;
}

signed char MeshAdapter::preferredInterpolationOnDomain(int functionSpaceType_source, int functionSpaceType_target) const
{
    if (probeInterpolationOnDomain(functionSpaceType_source, functionSpaceType_target))
        return 1;

    if (probeInterpolationOnDomain(functionSpaceType_target, functionSpaceType_source))
        return -1;

    return 0;
}

bool MeshAdapter::probeInterpolationACross(int functionSpaceType_source,
        const escript::AbstractDomain& targetDomain,
        int functionSpaceType_target) const
{
    return false;
}

bool MeshAdapter::operator==(const escript::AbstractDomain& other) const
{
    const MeshAdapter* temp=dynamic_cast<const MeshAdapter*>(&other);
    if (temp) {
        return (m_finleyMesh==temp->m_finleyMesh);
    } else {
        return false;
    }
}

bool MeshAdapter::operator!=(const escript::AbstractDomain& other) const
{
    return !(operator==(other));
}

int MeshAdapter::getSystemMatrixTypeId(const int solver, const int preconditioner, const int package, const bool symmetry) const
{
    Mesh* mesh=m_finleyMesh.get();
    return SystemMatrixAdapter::getSystemMatrixTypeId(solver, preconditioner,
               package, symmetry, mesh->MPIInfo);
}

int MeshAdapter::getTransportTypeId(const int solver, const int preconditioner, const int package, const bool symmetry) const
{
    Mesh* mesh=m_finleyMesh.get();
    return TransportProblemAdapter::getTransportTypeId(solver, preconditioner,
                package, symmetry, mesh->MPIInfo);
}

escript::Data MeshAdapter::getX() const
{
    return continuousFunction(*this).getX();
}

escript::Data MeshAdapter::getNormal() const
{
    return functionOnBoundary(*this).getNormal();
}

escript::Data MeshAdapter::getSize() const
{
    return escript::function(*this).getSize();
}

const int* MeshAdapter::borrowSampleReferenceIDs(int functionSpaceType) const
{
    int *out = NULL;
    Mesh* mesh=m_finleyMesh.get();
    switch (functionSpaceType) {
        case Nodes:
            out=mesh->Nodes->Id;
            break;
        case ReducedNodes:
            out=mesh->Nodes->reducedNodesId;
            break;
        case Elements:
        case ReducedElements:
            out=mesh->Elements->Id;
            break;
        case FaceElements:
        case ReducedFaceElements:
            out=mesh->FaceElements->Id;
            break;
        case Points:
            out=mesh->Points->Id;
            break;
        case ContactElementsZero:
        case ReducedContactElementsZero:
        case ContactElementsOne:
        case ReducedContactElementsOne:
            out=mesh->ContactElements->Id;
            break;
        case DegreesOfFreedom:
            out=mesh->Nodes->degreesOfFreedomId;
            break;
        case ReducedDegreesOfFreedom:
            out=mesh->Nodes->reducedDegreesOfFreedomId;
            break;
        default:
            stringstream temp;
            temp << "Error - Invalid function space type: " << functionSpaceType << " for domain: " << getDescription();
            throw FinleyAdapterException(temp.str());
            break;
    }
    return out;
}
int MeshAdapter::getTagFromSampleNo(int functionSpaceType, int sampleNo) const
{
    int out=0;
    Mesh* mesh=m_finleyMesh.get();
    switch (functionSpaceType) {
        case Nodes:
            out=mesh->Nodes->Tag[sampleNo];
            break;
        case ReducedNodes:
            throw FinleyAdapterException(" Error - ReducedNodes does not support tags.");
            break;
        case Elements:
        case ReducedElements:
            out=mesh->Elements->Tag[sampleNo];
            break;
        case FaceElements:
        case ReducedFaceElements:
            out=mesh->FaceElements->Tag[sampleNo];
            break;
        case Points:
            out=mesh->Points->Tag[sampleNo];
            break;
        case ContactElementsZero:
        case ReducedContactElementsZero:
        case ContactElementsOne:
        case ReducedContactElementsOne:
            out=mesh->ContactElements->Tag[sampleNo];
            break;
        case DegreesOfFreedom:
            throw FinleyAdapterException(" Error - DegreesOfFreedom does not support tags.");
            break;
        case ReducedDegreesOfFreedom:
            throw FinleyAdapterException(" Error - ReducedDegreesOfFreedom does not support tags.");
            break;
        default:
            stringstream temp;
            temp << "Error - Invalid function space type: " << functionSpaceType << " for domain: " << getDescription();
            throw FinleyAdapterException(temp.str());
            break;
    }
    return out;
}


void MeshAdapter::setTags(const int functionSpaceType, const int newTag, const escript::Data& mask) const
{
    Mesh* mesh=m_finleyMesh.get();
    switch(functionSpaceType) {
        case Nodes:
            mesh->Nodes->setTags(newTag, mask);
            break;
        case ReducedNodes:
            throw FinleyAdapterException("Error - ReducedNodes does not support tags");
        case DegreesOfFreedom:
            throw FinleyAdapterException("Error - DegreesOfFreedom does not support tags");
        case ReducedDegreesOfFreedom:
            throw FinleyAdapterException("Error - ReducedDegreesOfFreedom does not support tags");
        case Elements:
        case ReducedElements:
            mesh->Elements->setTags(newTag, mask);
            break;
        case FaceElements:
        case ReducedFaceElements:
            mesh->FaceElements->setTags(newTag, mask);
            break;
        case Points:
            mesh->Points->setTags(newTag, mask);
            break;
        case ContactElementsZero:
        case ReducedContactElementsZero:
        case ContactElementsOne:
        case ReducedContactElementsOne:
            mesh->ContactElements->setTags(newTag, mask);
            break;
        default:
            stringstream temp;
            temp << "Error - Finley does not know anything about function space type " << functionSpaceType;
            throw FinleyAdapterException(temp.str());
    }
    checkFinleyError();
}

void MeshAdapter::setTagMap(const string& name, int tag)
{
    Mesh* mesh=m_finleyMesh.get();
    mesh->addTagMap(name.c_str(), tag);
    checkFinleyError();
}

int MeshAdapter::getTag(const string& name) const
{
    Mesh* mesh=m_finleyMesh.get();
    int tag = mesh->getTag(name.c_str());
    checkFinleyError();
    return tag;
}

bool MeshAdapter::isValidTagName(const string& name) const
{
    Mesh* mesh=m_finleyMesh.get();
    return mesh->isValidTagName(name.c_str());
}

string MeshAdapter::showTagNames() const
{
    stringstream temp;
    Mesh* mesh=m_finleyMesh.get();
    TagMap::const_iterator it = mesh->tagMap.begin();
    while (it != mesh->tagMap.end()) {
        temp << it->first;
        ++it;
        if (it != mesh->tagMap.end())
            temp << ", ";
    }
    return temp.str();
}

int MeshAdapter::getNumberOfTagsInUse(int functionSpaceCode) const
{
    Mesh* mesh=m_finleyMesh.get();
    switch(functionSpaceCode) {
        case Nodes:
            return mesh->Nodes->tagsInUse.size();
        case ReducedNodes:
            throw FinleyAdapterException("Error - ReducedNodes does not support tags");
        case DegreesOfFreedom:
            throw FinleyAdapterException("Error - DegreesOfFreedom does not support tags");
        case ReducedDegreesOfFreedom:
            throw FinleyAdapterException("Error - ReducedDegreesOfFreedom does not support tags");
        case Elements:
        case ReducedElements:
            return mesh->Elements->tagsInUse.size();
        case FaceElements:
        case ReducedFaceElements:
            return mesh->FaceElements->tagsInUse.size();
        case Points:
            return mesh->Points->tagsInUse.size();
        case ContactElementsZero:
        case ReducedContactElementsZero:
        case ContactElementsOne:
        case ReducedContactElementsOne:
            return mesh->ContactElements->tagsInUse.size();
        default:
            stringstream ss;
            ss << "Finley does not know anything about function space type "
                 << functionSpaceCode;
            throw FinleyAdapterException(ss.str());
    }
    return 0;
}

const int* MeshAdapter::borrowListOfTagsInUse(int functionSpaceCode) const
{
    Mesh* mesh=m_finleyMesh.get();
    switch(functionSpaceCode) {
        case Nodes:
            return &mesh->Nodes->tagsInUse[0];
        case ReducedNodes:
            throw FinleyAdapterException("Error - ReducedNodes does not support tags");
        case DegreesOfFreedom:
            throw FinleyAdapterException("Error - DegreesOfFreedom does not support tags");
        case ReducedDegreesOfFreedom:
            throw FinleyAdapterException("Error - ReducedDegreesOfFreedom does not support tags");
        case Elements:
        case ReducedElements:
            return &mesh->Elements->tagsInUse[0];
        case FaceElements:
        case ReducedFaceElements:
            return &mesh->FaceElements->tagsInUse[0];
        case Points:
            return &mesh->Points->tagsInUse[0];
        case ContactElementsZero:
        case ReducedContactElementsZero:
        case ContactElementsOne:
        case ReducedContactElementsOne:
            return &mesh->ContactElements->tagsInUse[0];
        default:
            stringstream temp;
            temp << "Error - Finley does not know anything about function space type " << functionSpaceCode;
            throw FinleyAdapterException(temp.str());
    }
    return NULL;
}


bool MeshAdapter::canTag(int functionSpaceCode) const
{
    switch(functionSpaceCode) {
        case Nodes:
        case Elements:
        case ReducedElements:
        case FaceElements:
        case ReducedFaceElements:
        case Points:
        case ContactElementsZero:
        case ReducedContactElementsZero:
        case ContactElementsOne:
        case ReducedContactElementsOne:
            return true;
        default:
            return false;
    }
}

escript::AbstractDomain::StatusType MeshAdapter::getStatus() const
{
    Mesh* mesh=m_finleyMesh.get();
    return mesh->getStatus();
}

int MeshAdapter::getApproximationOrder(const int functionSpaceCode) const
{
    Mesh* mesh=m_finleyMesh.get();
    int order =-1;
    switch(functionSpaceCode) {
        case Nodes:
        case DegreesOfFreedom:
            order=mesh->approximationOrder;
            break;
        case ReducedNodes:
        case ReducedDegreesOfFreedom:
            order=mesh->reducedApproximationOrder;
            break;
        case Elements:
        case FaceElements:
        case Points:
        case ContactElementsZero:
        case ContactElementsOne:
            order=mesh->integrationOrder;
            break;
        case ReducedElements:
        case ReducedFaceElements:
        case ReducedContactElementsZero:
        case ReducedContactElementsOne:
            order=mesh->reducedIntegrationOrder;
            break;
        default:
            stringstream temp;
            temp << "Error - Finley does not know anything about function space type " << functionSpaceCode;
            throw FinleyAdapterException(temp.str());
    }
    return order;
}

bool MeshAdapter::supportsContactElements() const
{
    return true;
}

void MeshAdapter::addDiracPoints(const vector<double>& points,
                                 const vector<int>& tags) const
{
    // points will be flattened
    const int dim = getDim();
    int numPoints=points.size()/dim;
    int numTags=tags.size();
    Mesh* mesh=m_finleyMesh.get();

    if ( points.size() % dim != 0 ) {
        throw FinleyAdapterException("Error - number of coords does not appear to be a multiple of dimension.");
    }

    if ((numTags > 0) && (numPoints != numTags))
        throw FinleyAdapterException("Error - if tags are given number of tags and points must match.");

    mesh->addPoints(numPoints, &points[0], &tags[0]);
    checkFinleyError();
}

// void MeshAdapter::addDiracPoints(const boost::python::list& points, const boost::python::list& tags) const
// {
//       const int dim = getDim();
//       int numPoints=boost::python::extract<int>(points.attr("__len__")());
//       int numTags=boost::python::extract<int>(tags.attr("__len__")());
//       Mesh* mesh=m_finleyMesh.get();
//
//       if  ( (numTags > 0) && ( numPoints !=  numTags ) )
//       throw FinleyAdapterException("Error - if tags are given number of tags and points must match.");
//
//       double* points_ptr=TMPMEMALLOC(numPoints * dim, double);
//       int*    tags_ptr= TMPMEMALLOC(numPoints, int);
//
//       for (int i=0;i<numPoints;++i) {
//         int tag_id=-1;
//         int numComps=boost::python::extract<int>(points[i].attr("__len__")());
//         if  ( numComps !=   dim ) {
//                stringstream temp;
//                temp << "Error - illegal number of components " << numComps << " for point " << i;
//                throw FinleyAdapterException(temp.str());
//         }
//         points_ptr[ i * dim     ] = boost::python::extract<double>(points[i][0]);
//         if ( dim > 1 ) points_ptr[ i * dim + 1 ] = boost::python::extract<double>(points[i][1]);
//         if ( dim > 2 ) points_ptr[ i * dim + 2 ] = boost::python::extract<double>(points[i][2]);
//
//         if ( numTags > 0) {
//                boost::python::extract<string> ex_str(tags[i]);
//                if  ( ex_str.check() ) {
//                    tag_id=getTag( ex_str());
//                } else {
//                     boost::python::extract<int> ex_int(tags[i]);
//                     if ( ex_int.check() ) {
//                         tag_id=ex_int();
//                     } else {
//                          stringstream temp;
//                          temp << "Error - unable to extract tag for point " << i;
//                          throw FinleyAdapterException(temp.str());
//                    }
//                }
//         }
//            tags_ptr[i]=tag_id;
//       }
//
//       Finley_Mesh_addPoints(mesh, numPoints, points_ptr, tags_ptr);
//       checkPasoError();
//
//       TMPMEMFREE(points_ptr);
//       TMPMEMFREE(tags_ptr);
// }

/*
void MeshAdapter:: addDiracPoint( const boost::python::list& point, const int tag) const
{
    boost::python::list points =  boost::python::list();
    boost::python::list tags =  boost::python::list();
    points.append(point);
    tags.append(tag);
    addDiracPoints(points, tags);
}
*/

/*
void MeshAdapter:: addDiracPointWithTagName( const boost::python::list& point, const string& tag) const
{
    boost::python::list points =   boost::python::list();
    boost::python::list tags =   boost::python::list();
    points.append(point);
    tags.append(tag);
    addDiracPoints(points, tags);
}
*/
}  // end of namespace

