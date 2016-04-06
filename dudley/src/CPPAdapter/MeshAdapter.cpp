
/*****************************************************************************
*
* Copyright (c) 2003-2016 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/

#include "MeshAdapter.h"
#include <dudley/Assemble.h>
#include <dudley/DudleyException.h>

#include <escript/Data.h>
#include <escript/DataFactory.h>
#include <escript/Random.h>
#include <escript/SolverOptions.h>

#include <paso/SystemMatrix.h>
#include <paso/Transport.h>

#ifdef ESYS_HAVE_NETCDF
#include <netcdfcpp.h>
#endif

#ifdef ESYS_HAVE_TRILINOS
#include <trilinoswrap/TrilinosMatrixAdapter.h>

using esys_trilinos::TrilinosMatrixAdapter;
using esys_trilinos::const_TrilinosGraph_ptr;
#endif

using namespace std;
namespace bp = boost::python;
using escript::ValueError;

namespace dudley {

// define the static constants
MeshAdapter::FunctionSpaceNamesMapType MeshAdapter::m_functionSpaceTypeNames;
// dudley only supports single approximation order 
const int MeshAdapter::DegreesOfFreedom=DUDLEY_DEGREES_OF_FREEDOM;
const int MeshAdapter::Nodes=DUDLEY_NODES;
const int MeshAdapter::Elements=DUDLEY_ELEMENTS;
const int MeshAdapter::ReducedElements=DUDLEY_REDUCED_ELEMENTS;
const int MeshAdapter::FaceElements=DUDLEY_FACE_ELEMENTS;
const int MeshAdapter::ReducedFaceElements=DUDLEY_REDUCED_FACE_ELEMENTS;
const int MeshAdapter::Points=DUDLEY_POINTS;

MeshAdapter::MeshAdapter(Mesh* dudleyMesh) :
    m_dudleyMesh(dudleyMesh)
{
    setFunctionSpaceTypeNames();
}

// The copy constructor should just increment the use count
MeshAdapter::MeshAdapter(const MeshAdapter& in) :
    m_dudleyMesh(in.m_dudleyMesh)
{
    setFunctionSpaceTypeNames();
}

MeshAdapter::~MeshAdapter()
{
}

escript::JMPI MeshAdapter::getMPI() const
{
    return m_dudleyMesh->MPIInfo;
}

int MeshAdapter::getMPISize() const
{
    return getMPI()->size;
}

int MeshAdapter::getMPIRank() const
{
    return getMPI()->rank;
}

void MeshAdapter::MPIBarrier() const
{
#ifdef ESYS_MPI
    MPI_Barrier(getMPIComm());
#endif
}

bool MeshAdapter::onMasterProcessor() const
{
    return getMPIRank() == 0;
}

MPI_Comm MeshAdapter::getMPIComm() const
{
    return getMPI()->comm;
}

Mesh* MeshAdapter::getMesh() const
{
    return m_dudleyMesh.get();
}

void MeshAdapter::write(const string& fileName) const
{
    m_dudleyMesh->write(fileName);
}

void MeshAdapter::Print_Mesh_Info(bool full) const
{
    m_dudleyMesh->printInfo(full);
}

void MeshAdapter::dump(const string& fileName) const
{
#ifdef ESYS_HAVE_NETCDF
   const NcDim* ncdims[12] = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
   NcVar* ids;
   index_t* index_ptr;
#ifdef ESYS_INDEXTYPE_LONG
    NcType ncIdxType = ncLong;
#else
    NcType ncIdxType = ncInt;
#endif
   Mesh* mesh = m_dudleyMesh.get();
   int num_Tags = 0;
   int mpi_size                  = getMPISize();
   int mpi_rank                  = getMPIRank();
   int numDim                    = mesh->Nodes->numDim;
   dim_t numNodes                = mesh->Nodes->getNumNodes();
   dim_t num_Elements            = mesh->Elements->numElements;
   dim_t num_FaceElements        = mesh->FaceElements->numElements;
   dim_t num_Points              = mesh->Points->numElements;
   int num_Elements_numNodes     = mesh->Elements->numNodes;
   int num_FaceElements_numNodes = mesh->FaceElements->numNodes;
#ifdef ESYS_MPI
   MPI_Status status;
#endif

    // Incoming token indicates it's my turn to write
#ifdef ESYS_MPI
    if (mpi_rank > 0)
        MPI_Recv(&num_Tags, 0, MPI_INT, mpi_rank-1, 81800, getMPIComm(), &status);
#endif

    const string newFileName(mesh->MPIInfo->appendRankToFileName(fileName));

    // Figure out how much storage is required for tags
    num_Tags = mesh->tagMap.size();

    // NetCDF error handler
    NcError err(NcError::verbose_nonfatal);
    // Create the file
    NcFile dataFile(newFileName.c_str(), NcFile::Replace);
    string msgPrefix("Error in MeshAdapter::dump: NetCDF operation failed - ");
    // check if writing was successful
    if (!dataFile.is_valid())
        throw DudleyException(msgPrefix + "Open file for output");

    // Define dimensions (num_Elements and dim_Elements are identical,
    // dim_Elements only appears if > 0)
    if (! (ncdims[0] = dataFile.add_dim("numNodes", numNodes)) )
        throw DudleyException(msgPrefix+"add_dim(numNodes)");
    if (! (ncdims[1] = dataFile.add_dim("numDim", numDim)) )
        throw DudleyException(msgPrefix+"add_dim(numDim)");
    if (! (ncdims[2] = dataFile.add_dim("mpi_size_plus_1", mpi_size+1)) )
        throw DudleyException(msgPrefix+"add_dim(mpi_size)");
    if (num_Elements > 0)
        if (! (ncdims[3] = dataFile.add_dim("dim_Elements", num_Elements)) )
            throw DudleyException(msgPrefix+"add_dim(dim_Elements)");
    if (num_FaceElements > 0)
        if (! (ncdims[4] = dataFile.add_dim("dim_FaceElements", num_FaceElements)) )
         throw DudleyException(msgPrefix+"add_dim(dim_FaceElements)");
    if (num_Points > 0)
        if (! (ncdims[6] = dataFile.add_dim("dim_Points", num_Points)) )
            throw DudleyException(msgPrefix+"add_dim(dim_Points)");
    if (num_Elements > 0)
        if (! (ncdims[7] = dataFile.add_dim("dim_Elements_Nodes", num_Elements_numNodes)) )
            throw DudleyException(msgPrefix+"add_dim(dim_Elements_Nodes)");
    if (num_FaceElements > 0)
        if (! (ncdims[8] = dataFile.add_dim("dim_FaceElements_numNodes", num_FaceElements_numNodes)) )
            throw DudleyException(msgPrefix+"add_dim(dim_FaceElements_numNodes)");
    if (num_Tags > 0)
        if (! (ncdims[10] = dataFile.add_dim("dim_Tags", num_Tags)) )
            throw DudleyException(msgPrefix+"add_dim(dim_Tags)");

    // Attributes: MPI size, MPI rank, Name, order, reduced_order
    if (!dataFile.add_att("index_size", (int)sizeof(index_t)))
        throw DudleyException(msgPrefix+"add_att(index_size)");
    if (!dataFile.add_att("mpi_size", mpi_size) )
        throw DudleyException(msgPrefix+"add_att(mpi_size)");
    if (!dataFile.add_att("mpi_rank", mpi_rank) )
        throw DudleyException(msgPrefix+"add_att(mpi_rank)");
    if (!dataFile.add_att("Name",mesh->m_name.c_str()) )
        throw DudleyException(msgPrefix+"add_att(Name)");
    if (!dataFile.add_att("numDim",numDim) )
        throw DudleyException(msgPrefix+"add_att(order)");
    if (!dataFile.add_att("order",mesh->integrationOrder) )
        throw DudleyException(msgPrefix+"add_att(order)");
    if (!dataFile.add_att("reduced_order",mesh->reducedIntegrationOrder) )
        throw DudleyException(msgPrefix+"add_att(reduced_order)");
    if (!dataFile.add_att("numNodes",numNodes) )
        throw DudleyException(msgPrefix+"add_att(numNodes)");
    if (!dataFile.add_att("num_Elements",num_Elements) )
        throw DudleyException(msgPrefix+"add_att(num_Elements)");
    if (!dataFile.add_att("num_FaceElements",num_FaceElements) )
        throw DudleyException(msgPrefix+"add_att(num_FaceElements)");
    if (!dataFile.add_att("num_Points",num_Points) )
        throw DudleyException(msgPrefix+"add_att(num_Points)");
    if (!dataFile.add_att("num_Elements_numNodes",num_Elements_numNodes) )
        throw DudleyException(msgPrefix+"add_att(num_Elements_numNodes)");
    if (!dataFile.add_att("num_FaceElements_numNodes",num_FaceElements_numNodes) )
        throw DudleyException(msgPrefix+"add_att(num_FaceElements_numNodes)");
    if (!dataFile.add_att("Elements_TypeId", mesh->Elements->etype) )
        throw DudleyException(msgPrefix+"add_att(Elements_TypeId)");
    if (!dataFile.add_att("FaceElements_TypeId", mesh->FaceElements->etype) )
        throw DudleyException(msgPrefix+"add_att(FaceElements_TypeId)");
    if (!dataFile.add_att("Points_TypeId", mesh->Points->etype) )
        throw DudleyException(msgPrefix+"add_att(Points_TypeId)");
    if (!dataFile.add_att("num_Tags", num_Tags) )
        throw DudleyException(msgPrefix+"add_att(num_Tags)");

    // // // // // Nodes // // // // //

    // Nodes nodeDistribution
    if (! ( ids = dataFile.add_var("Nodes_NodeDistribution", ncIdxType, ncdims[2])) )
        throw DudleyException(msgPrefix+"add_var(Nodes_NodeDistribution)");
    index_ptr = &mesh->Nodes->nodesDistribution->first_component[0];
    if (! (ids->put(index_ptr, mpi_size+1)) )
        throw DudleyException(msgPrefix+"put(Nodes_NodeDistribution)");

    // Nodes degreesOfFreedomDistribution
    if (! ( ids = dataFile.add_var("Nodes_DofDistribution", ncIdxType, ncdims[2])) )
        throw DudleyException(msgPrefix+"add_var(Nodes_DofDistribution)");
    index_ptr = &mesh->Nodes->dofDistribution->first_component[0];
    if (! (ids->put(index_ptr, mpi_size+1)) )
        throw DudleyException(msgPrefix+"put(Nodes_DofDistribution)");

    // Only write nodes if non-empty because NetCDF doesn't like empty arrays
    // (it treats them as NC_UNLIMITED)
    if (numNodes > 0) {
        // Nodes Id
        if (! ( ids = dataFile.add_var("Nodes_Id", ncIdxType, ncdims[0])) )
            throw DudleyException(msgPrefix+"add_var(Nodes_Id)");
        if (! (ids->put(mesh->Nodes->Id, numNodes)) )
            throw DudleyException(msgPrefix+"put(Nodes_Id)");

        // Nodes Tag
        if (! ( ids = dataFile.add_var("Nodes_Tag", ncInt, ncdims[0])) )
            throw DudleyException(msgPrefix+"add_var(Nodes_Tag)");
        if (! (ids->put(mesh->Nodes->Tag, numNodes)) )
            throw DudleyException(msgPrefix+"put(Nodes_Tag)");

        // Nodes gDOF
        if (! ( ids = dataFile.add_var("Nodes_gDOF", ncIdxType, ncdims[0])) )
            throw DudleyException(msgPrefix+"add_var(Nodes_gDOF)");
        if (! (ids->put(mesh->Nodes->globalDegreesOfFreedom, numNodes)) )
            throw DudleyException(msgPrefix+"put(Nodes_gDOF)");

        // Nodes global node index
        if (! ( ids = dataFile.add_var("Nodes_gNI", ncIdxType, ncdims[0])) )
            throw DudleyException(msgPrefix+"add_var(Nodes_gNI)");
        if (! (ids->put(mesh->Nodes->globalNodesIndex, numNodes)) )
            throw DudleyException(msgPrefix+"put(Nodes_gNI)");

        // Nodes Coordinates
        if (! ( ids = dataFile.add_var("Nodes_Coordinates", ncDouble, ncdims[0], ncdims[1]) ) )
            throw DudleyException(msgPrefix+"add_var(Nodes_Coordinates)");
        if (! (ids->put(mesh->Nodes->Coordinates, numNodes, numDim)) )
            throw DudleyException(msgPrefix+"put(Nodes_Coordinates)");
    }

    // // // // // Elements // // // // //
    if (num_Elements > 0) {
        // Elements_Id
        if (! ( ids = dataFile.add_var("Elements_Id", ncIdxType, ncdims[3])) )
            throw DudleyException(msgPrefix+"add_var(Elements_Id)");
        if (! (ids->put(mesh->Elements->Id, num_Elements)) )
            throw DudleyException(msgPrefix+"put(Elements_Id)");

        // Elements_Tag
        if (! ( ids = dataFile.add_var("Elements_Tag", ncInt, ncdims[3])) )
            throw DudleyException(msgPrefix+"add_var(Elements_Tag)");
        if (! (ids->put(mesh->Elements->Tag, num_Elements)) )
            throw DudleyException(msgPrefix+"put(Elements_Tag)");

        // Elements_Owner
        if (! ( ids = dataFile.add_var("Elements_Owner", ncInt, ncdims[3])) )
            throw DudleyException(msgPrefix+"add_var(Elements_Owner)");
        if (! (ids->put(mesh->Elements->Owner, num_Elements)) )
            throw DudleyException(msgPrefix+"put(Elements_Owner)");

        // Elements_Color
        if (! ( ids = dataFile.add_var("Elements_Color", ncIdxType, ncdims[3])) )
            throw DudleyException(msgPrefix+"add_var(Elements_Color)");
        if (! (ids->put(mesh->Elements->Color, num_Elements)) )
            throw DudleyException(msgPrefix+"put(Elements_Color)");

        // Elements_Nodes
        if (! ( ids = dataFile.add_var("Elements_Nodes", ncIdxType, ncdims[3], ncdims[7]) ) )
            throw DudleyException(msgPrefix+"add_var(Elements_Nodes)");
        if (! (ids->put(mesh->Elements->Nodes, num_Elements, num_Elements_numNodes)) )
            throw DudleyException(msgPrefix+"put(Elements_Nodes)");
    }

    // // // // // Face_Elements // // // // //
    if (num_FaceElements > 0) {
        // FaceElements_Id
        if (!(ids = dataFile.add_var("FaceElements_Id", ncIdxType, ncdims[4])))
            throw DudleyException(msgPrefix+"add_var(FaceElements_Id)");
        if (!(ids->put(mesh->FaceElements->Id, num_FaceElements)))
            throw DudleyException(msgPrefix+"put(FaceElements_Id)");

        // FaceElements_Tag
        if (!(ids = dataFile.add_var("FaceElements_Tag", ncInt, ncdims[4])))
            throw DudleyException(msgPrefix+"add_var(FaceElements_Tag)");
        if (!(ids->put(mesh->FaceElements->Tag, num_FaceElements)))
            throw DudleyException(msgPrefix+"put(FaceElements_Tag)");

        // FaceElements_Owner
        if (!(ids = dataFile.add_var("FaceElements_Owner", ncInt, ncdims[4])))
            throw DudleyException(msgPrefix+"add_var(FaceElements_Owner)");
        if (!(ids->put(mesh->FaceElements->Owner, num_FaceElements)))
            throw DudleyException(msgPrefix+"put(FaceElements_Owner)");

        // FaceElements_Color
        if (!(ids = dataFile.add_var("FaceElements_Color", ncIdxType, ncdims[4])))
            throw DudleyException(msgPrefix+"add_var(FaceElements_Color)");
        if (!(ids->put(mesh->FaceElements->Color, num_FaceElements)))
            throw DudleyException(msgPrefix+"put(FaceElements_Color)");

        // FaceElements_Nodes
        if (!(ids = dataFile.add_var("FaceElements_Nodes", ncIdxType, ncdims[4], ncdims[8])))
            throw DudleyException(msgPrefix+"add_var(FaceElements_Nodes)");
        if (!(ids->put(mesh->FaceElements->Nodes, num_FaceElements, num_FaceElements_numNodes)))
            throw DudleyException(msgPrefix+"put(FaceElements_Nodes)");
    }

    // // // // // Points // // // // //
    if (num_Points > 0) {
        // Points_Id
        if (!(ids = dataFile.add_var("Points_Id", ncIdxType, ncdims[6])))
            throw DudleyException(msgPrefix+"add_var(Points_Id)");
        if (!(ids->put(mesh->Points->Id, num_Points)))
            throw DudleyException(msgPrefix+"put(Points_Id)");

        // Points_Tag
        if (!(ids = dataFile.add_var("Points_Tag", ncInt, ncdims[6])))
            throw DudleyException(msgPrefix+"add_var(Points_Tag)");
        if (!(ids->put(mesh->Points->Tag, num_Points)))
            throw DudleyException(msgPrefix+"put(Points_Tag)");

        // Points_Owner
        if (!(ids = dataFile.add_var("Points_Owner", ncInt, ncdims[6])))
            throw DudleyException(msgPrefix+"add_var(Points_Owner)");
        if (!(ids->put(mesh->Points->Owner, num_Points)))
            throw DudleyException(msgPrefix+"put(Points_Owner)");

        // Points_Color
        if (!(ids = dataFile.add_var("Points_Color", ncIdxType, ncdims[6])))
            throw DudleyException(msgPrefix+"add_var(Points_Color)");
        if (!(ids->put(mesh->Points->Color, num_Points)))
            throw DudleyException(msgPrefix+"put(Points_Color)");

        // Points_Nodes
        if (!(ids = dataFile.add_var("Points_Nodes", ncIdxType, ncdims[6])))
            throw DudleyException(msgPrefix+"add_var(Points_Nodes)");
        if (!(ids->put(mesh->Points->Nodes, num_Points)))
            throw DudleyException(msgPrefix+"put(Points_Nodes)");
    }

    // // // // // TagMap // // // // //
    if (num_Tags > 0) {
        // Temp storage to gather node IDs
        vector<int> Tags_keys;

        // Copy tag data into temp arrays
        TagMap::const_iterator it;
        for (it = mesh->tagMap.begin(); it != mesh->tagMap.end(); it++) {
            Tags_keys.push_back(it->second);
        }

        // Tags_keys
        if (!(ids = dataFile.add_var("Tags_keys", ncInt, ncdims[10])))
            throw DudleyException(msgPrefix+"add_var(Tags_keys)");
        if (!(ids->put(&Tags_keys[0], num_Tags)))
            throw DudleyException(msgPrefix+"put(Tags_keys)");

        // Tags_names_*
        // This is an array of strings, it should be stored as an array but
        // instead I have hacked in one attribute per string because the NetCDF
        // manual doesn't tell how to do an array of strings
        int i = 0;
        for (it = mesh->tagMap.begin(); it != mesh->tagMap.end(); it++, i++) {
            stringstream ss;
            ss << "Tags_name_" << i;
            const string name(ss.str());
            if (!dataFile.add_att(name.c_str(), it->first.c_str()))
                throw DudleyException(msgPrefix+"add_att(Tags_names_XX)");
        }
    }

    // Send token to next MPI process so he can take his turn
#ifdef ESYS_MPI
    if (mpi_rank < mpi_size-1)
        MPI_Send(&num_Tags, 0, MPI_INT, mpi_rank+1, 81800, getMPIComm());
#endif

    // NetCDF file is closed by destructor of NcFile object

#else
    throw DudleyException("MeshAdapter::dump: not configured with netCDF. "
                          "Please contact your installation manager.");
#endif // ESYS_HAVE_NETCDF
}

string MeshAdapter::getDescription() const
{
    return "DudleyMesh";
}

string MeshAdapter::functionSpaceTypeAsString(int functionSpaceType) const
{
    FunctionSpaceNamesMapType::iterator loc;
    loc = m_functionSpaceTypeNames.find(functionSpaceType);
    if (loc == m_functionSpaceTypeNames.end()) {
        return "Invalid function space type code.";
    } else {
        return loc->second;
    }
}

bool MeshAdapter::isValidFunctionSpaceType(int functionSpaceType) const
{
    FunctionSpaceNamesMapType::iterator loc;
    loc = m_functionSpaceTypeNames.find(functionSpaceType);
    return (loc != m_functionSpaceTypeNames.end());
}

void MeshAdapter::setFunctionSpaceTypeNames()
{
    m_functionSpaceTypeNames.insert(FunctionSpaceNamesMapType::value_type(
                DegreesOfFreedom,"Dudley_DegreesOfFreedom [Solution(domain)]"));
    m_functionSpaceTypeNames.insert(FunctionSpaceNamesMapType::value_type(
                Nodes,"Dudley_Nodes [ContinuousFunction(domain)]"));
    m_functionSpaceTypeNames.insert(FunctionSpaceNamesMapType::value_type(
                Elements,"Dudley_Elements [Function(domain)]"));
    m_functionSpaceTypeNames.insert(FunctionSpaceNamesMapType::value_type(
                ReducedElements,"Dudley_Reduced_Elements [ReducedFunction(domain)]"));
    m_functionSpaceTypeNames.insert(FunctionSpaceNamesMapType::value_type(
                FaceElements,"Dudley_Face_Elements [FunctionOnBoundary(domain)]"));
    m_functionSpaceTypeNames.insert(FunctionSpaceNamesMapType::value_type(
                ReducedFaceElements,"Dudley_Reduced_Face_Elements [ReducedFunctionOnBoundary(domain)]"));
    m_functionSpaceTypeNames.insert(FunctionSpaceNamesMapType::value_type(
                Points,"Dudley_Points [DiracDeltaFunctions(domain)]"));
}

int MeshAdapter::getContinuousFunctionCode() const
{
    return Nodes;
}

int MeshAdapter::getReducedContinuousFunctionCode() const
{
    return Nodes;
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
    throw DudleyException("Dudley does not support contact elements.");
}

int MeshAdapter::getReducedFunctionOnContactZeroCode() const
{
    throw DudleyException("Dudley does not support contact elements.");
}

int MeshAdapter::getFunctionOnContactOneCode() const
{
    throw DudleyException("Dudley does not support contact elements.");
}

int MeshAdapter::getReducedFunctionOnContactOneCode() const
{
    throw DudleyException("Dudley does not support contact elements.");
}

int MeshAdapter::getSolutionCode() const
{
    return DegreesOfFreedom;
}

int MeshAdapter::getReducedSolutionCode() const
{
    return DegreesOfFreedom;
}

int MeshAdapter::getDiracDeltaFunctionsCode() const
{
    return Points;
}

int MeshAdapter::getDim() const
{
    return m_dudleyMesh->getDim();
}

//
// Return the number of data points summed across all MPI processes
//
dim_t MeshAdapter::getNumDataPointsGlobal() const
{
    return m_dudleyMesh->Nodes->getGlobalNumNodes();
}

//
// return the number of data points per sample and the number of samples
// needed to represent data on a parts of the mesh.
//
pair<int,dim_t> MeshAdapter::getDataShape(int functionSpaceCode) const
{
    int numDataPointsPerSample = 0;
    dim_t numSamples = 0;
    Mesh* mesh = getMesh();
    switch (functionSpaceCode) {
        case Nodes:
            numDataPointsPerSample = 1;
            numSamples = mesh->Nodes->getNumNodes();
        break;
        case Elements:
            if (mesh->Elements) {
                numSamples = mesh->Elements->numElements;
                numDataPointsPerSample = mesh->Elements->numLocalDim+1;
            }
        break;
        case ReducedElements:
            if (mesh->Elements) {
                numSamples = mesh->Elements->numElements;
                numDataPointsPerSample =(mesh->Elements->numLocalDim==0)?0:1;
            }
        break;
        case FaceElements:
            if (mesh->FaceElements) {
                numDataPointsPerSample = mesh->FaceElements->numLocalDim+1;
                numSamples = mesh->FaceElements->numElements;
            }
        break;
        case ReducedFaceElements:
            if (mesh->FaceElements) {
                numDataPointsPerSample = (mesh->FaceElements->numLocalDim==0)?0:1;
                numSamples = mesh->FaceElements->numElements;
            }
        break;
        case Points:
            if (mesh->Points) {
                numDataPointsPerSample = 1;
                numSamples = mesh->Points->numElements;
            }
        break;
        case DegreesOfFreedom:
            if (mesh->Nodes) {
                numDataPointsPerSample = 1;
                numSamples = mesh->Nodes->getNumDegreesOfFreedom();
            }
        break;
        default:
            stringstream ss;
            ss << "Invalid function space type: " << functionSpaceCode
                << " for domain " << getDescription();
            throw DudleyException(ss.str());
            break;
    }
    return pair<int,dim_t>(numDataPointsPerSample,numSamples);
}

//
// adds linear PDE of second order into a given stiffness matrix and right hand side:
//
void MeshAdapter::addPDEToSystem(
        escript::AbstractSystemMatrix& mat, escript::Data& rhs,
        const escript::Data& A, const escript::Data& B, const escript::Data& C,
        const escript::Data& D, const escript::Data& X, const escript::Data& Y,
        const escript::Data& d, const escript::Data& y,
        const escript::Data& d_contact, const escript::Data& y_contact,
        const escript::Data& d_dirac, const escript::Data& y_dirac) const
{
    if (!d_contact.isEmpty() || !y_contact.isEmpty())
        throw DudleyException("Dudley does not support contact elements");

#ifdef ESYS_HAVE_TRILINOS
    TrilinosMatrixAdapter* tm = dynamic_cast<TrilinosMatrixAdapter*>(&mat);
    if (tm) {
        tm->resumeFill();
    }
#endif

    Mesh* mesh = m_dudleyMesh.get();
    Assemble_PDE(mesh->Nodes, mesh->Elements, mat.getPtr(), rhs,
                 A, B, C, D, X, Y);
    Assemble_PDE(mesh->Nodes, mesh->FaceElements, mat.getPtr(), rhs,
                 escript::Data(), escript::Data(), escript::Data(), d,
                 escript::Data(), y);
    Assemble_PDE(mesh->Nodes, mesh->Points, mat.getPtr(), rhs, escript::Data(),
                 escript::Data(), escript::Data(), d_dirac,
                 escript::Data(), y_dirac);

#ifdef ESYS_HAVE_TRILINOS
    if (tm) {
        tm->fillComplete(true);
    }
#endif
}

void MeshAdapter::addPDEToLumpedSystem(escript::Data& mat,
                                       const escript::Data& D,
                                       const escript::Data& d,
                                       const escript::Data& d_dirac,
                                       bool useHRZ) const
{
    Mesh* mesh = m_dudleyMesh.get();
    Assemble_LumpedSystem(mesh->Nodes, mesh->Elements, mat, D, useHRZ);
    Assemble_LumpedSystem(mesh->Nodes, mesh->FaceElements, mat, d, useHRZ);
    Assemble_LumpedSystem(mesh->Nodes, mesh->Points, mat, d_dirac, useHRZ);
}

//
// adds linear PDE of second order into the right hand side only
//
void MeshAdapter::addPDEToRHS(escript::Data& rhs, const escript::Data& X,
          const escript::Data& Y, const escript::Data& y,
          const escript::Data& y_contact, const escript::Data& y_dirac) const
{
    if (!y_contact.isEmpty())
        throw DudleyException("Dudley does not support y_contact");

    Mesh* mesh=m_dudleyMesh.get();

    Assemble_PDE(mesh->Nodes, mesh->Elements, escript::ASM_ptr(), rhs,
                 escript::Data(), escript::Data(), escript::Data(),
                 escript::Data(), X, Y);

    Assemble_PDE(mesh->Nodes, mesh->FaceElements, escript::ASM_ptr(), rhs,
                 escript::Data(), escript::Data(), escript::Data(),
                 escript::Data(), escript::Data(), y);

    Assemble_PDE(mesh->Nodes, mesh->Points, escript::ASM_ptr(), rhs,
                 escript::Data(), escript::Data(), escript::Data(),
                 escript::Data(), escript::Data(), y_dirac);
}

//
// adds PDE of second order into a transport problem
//
void MeshAdapter::addPDEToTransportProblem(
        escript::AbstractTransportProblem& tp, escript::Data& source,
        const escript::Data& M, const escript::Data& A, const escript::Data& B,
        const escript::Data& C, const  escript::Data& D,const escript::Data& X,
        const escript::Data& Y, const escript::Data& d, const escript::Data& y,
        const escript::Data& d_contact,const escript::Data& y_contact,
        const escript::Data& d_dirac, const escript::Data& y_dirac) const
{
    if (!d_contact.isEmpty())
        throw DudleyException("Dudley does not support d_contact");
    if (!y_contact.isEmpty())
        throw DudleyException("Dudley does not support y_contact");

    paso::TransportProblem* ptp = dynamic_cast<paso::TransportProblem*>(&tp);
    if (!ptp)
        throw DudleyException("Dudley only accepts Paso transport problems");

    source.expand();

    Mesh* mesh = m_dudleyMesh.get();

    Assemble_PDE(mesh->Nodes, mesh->Elements, ptp->borrowMassMatrix(), source,
                 escript::Data(), escript::Data(), escript::Data(), M,
                 escript::Data(), escript::Data());

    Assemble_PDE(mesh->Nodes, mesh->Elements, ptp->borrowTransportMatrix(),
                 source, A, B, C, D, X, Y);

    Assemble_PDE(mesh->Nodes, mesh->FaceElements, ptp->borrowTransportMatrix(),
                 source, escript::Data(), escript::Data(), escript::Data(), d,
                 escript::Data(), y);

    Assemble_PDE(mesh->Nodes, mesh->Points, ptp->borrowTransportMatrix(),
                 source, escript::Data(), escript::Data(), escript::Data(),
                 d_dirac, escript::Data(), y_dirac);
}

//
// interpolates data between different function spaces:
//
void MeshAdapter::interpolateOnDomain(escript::Data& target,
                                      const escript::Data& in) const
{
    if (*in.getFunctionSpace().getDomain() != *this)  
        throw DudleyException("Illegal domain of interpolant.");
    if (*target.getFunctionSpace().getDomain() != *this) 
        throw DudleyException("Illegal domain of interpolation target.");

    Mesh* mesh = m_dudleyMesh.get();
    switch (in.getFunctionSpace().getTypeCode()) {
        case Nodes:
            switch (target.getFunctionSpace().getTypeCode()) {
                case Nodes:
                case DegreesOfFreedom:
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
                default:
                    stringstream ss;
                    ss << "interpolateOnDomain: Dudley does not know anything "
                          "about function space type "
                          << target.getFunctionSpace().getTypeCode();
                    throw DudleyException(ss.str());
                    break;
            }
        break;
        case Elements:
            if (target.getFunctionSpace().getTypeCode() == Elements) {
                Assemble_CopyElementData(mesh->Elements, target, in);
            } else if (target.getFunctionSpace().getTypeCode()==ReducedElements) {
                Assemble_AverageElementData(mesh->Elements, target, in);
            } else {
                throw DudleyException("No interpolation with data on elements possible.");
            }
            break;
        case ReducedElements:
            if (target.getFunctionSpace().getTypeCode() == ReducedElements) {
                Assemble_CopyElementData(mesh->Elements, target, in);
            } else {
                throw DudleyException("No interpolation with data on elements "
                                   "with reduced integration order possible.");
            }
            break;
        case FaceElements:
            if (target.getFunctionSpace().getTypeCode() == FaceElements) {
                Assemble_CopyElementData(mesh->FaceElements, target, in);
            } else if (target.getFunctionSpace().getTypeCode() == ReducedFaceElements) {
                Assemble_AverageElementData(mesh->FaceElements, target, in);
            } else {
                throw DudleyException("No interpolation with data on face elements possible.");
            }
            break;
        case ReducedFaceElements:
            if (target.getFunctionSpace().getTypeCode() == ReducedFaceElements) {
                Assemble_CopyElementData(mesh->FaceElements, target, in);
            } else {
                throw DudleyException("No interpolation with data on face "
                          "elements with reduced integration order possible.");
            }
            break;
        case Points:
            if (target.getFunctionSpace().getTypeCode() == Points) {
                Assemble_CopyElementData(mesh->Points, target, in);
            } else {
                throw DudleyException("No interpolation with data on points possible.");
            }
            break;
        case DegreesOfFreedom:
            switch (target.getFunctionSpace().getTypeCode()) {
                case DegreesOfFreedom:
                Assemble_CopyNodalData(mesh->Nodes, target, in);
                break;
            
                case Nodes:
                if (getMPISize() > 1) {
                    escript::Data temp = escript::Data(in);
                    temp.expand();
                    Assemble_CopyNodalData(mesh->Nodes, target, temp);
                } else {
                    Assemble_CopyNodalData(mesh->Nodes, target, in);
                }
                break;
                case Elements:
                case ReducedElements:
                if (getMPISize() > 1) {
                    escript::Data temp = escript::Data(in, continuousFunction(*this) );
                    Assemble_interpolate(mesh->Nodes, mesh->Elements, temp, target);
                } else {
                    Assemble_interpolate(mesh->Nodes, mesh->Elements, in, target);
                }
                break;
                case FaceElements:
                case ReducedFaceElements:
                if (getMPISize() > 1) {
                    escript::Data temp = escript::Data(in, continuousFunction(*this) );
                    Assemble_interpolate(mesh->Nodes, mesh->FaceElements, temp, target);
                } else {
                    Assemble_interpolate(mesh->Nodes, mesh->FaceElements, in, target);
                }
                break;
                case Points:
                if (getMPISize() > 1) {
                    //escript::Data temp=escript::Data(in, continuousFunction(*this) );
                } else {
                    Assemble_interpolate(mesh->Nodes, mesh->Points, in, target);
                }
                break;
                default:
                    stringstream ss;
                    ss << "interpolateOnDomain: Dudley does not know anything "
                          "about function space type "
                       << target.getFunctionSpace().getTypeCode();
                    throw DudleyException(ss.str());
                    break;
            }
            break;
       default:
          stringstream ss;
          ss << "interpolateOnDomain: Dudley does not know anything about "
                "function space type " << in.getFunctionSpace().getTypeCode();
          throw DudleyException(ss.str());
          break;
    }
}

//
// copies the locations of sample points into x:
//
void MeshAdapter::setToX(escript::Data& arg) const
{
    if (*arg.getFunctionSpace().getDomain() != *this) 
        throw DudleyException("setToX: Illegal domain of data point locations");

    Mesh* mesh = m_dudleyMesh.get();
    // in case of appropriate function space we can do the job directly:
    if (arg.getFunctionSpace().getTypeCode() == Nodes) {
        Assemble_NodeCoordinates(mesh->Nodes, arg);
    } else {
        escript::Data tmp_data = Vector(0., continuousFunction(*this), true);
        Assemble_NodeCoordinates(mesh->Nodes, tmp_data);
        // this is then interpolated onto arg:
        interpolateOnDomain(arg, tmp_data);
    }
}

//
// return the normal vectors at the location of data points as a Data object:
//
void MeshAdapter::setToNormal(escript::Data& normal) const
{
    if (*normal.getFunctionSpace().getDomain() != *this) 
        throw ValueError("setToNormal: Illegal domain of normal locations");

    Mesh* mesh=m_dudleyMesh.get();
    if (normal.getFunctionSpace().getTypeCode() == FaceElements ||
            normal.getFunctionSpace().getTypeCode() == ReducedFaceElements) {
        Assemble_getNormal(mesh->Nodes, mesh->FaceElements, normal);
    } else {
        stringstream ss;
        ss << "setToNormal: Illegal function space type "
           << normal.getFunctionSpace().getTypeCode();
        throw ValueError(ss.str());
    }
}

//
// interpolates data to other domain
//
void MeshAdapter::interpolateAcross(escript::Data& target,
                                    const escript::Data& source) const
{
    throw escript::NotImplementedError("Dudley does not allow interpolation "
                                       "across domains.");
}

//
// calculates the integral of a function defined of arg:
//
void MeshAdapter::setToIntegrals(vector<double>& integrals,
                                 const escript::Data& arg) const
{
    if (*arg.getFunctionSpace().getDomain() != *this) 
        throw ValueError("setToIntegrals: Illegal domain of integration kernel");

    Mesh* mesh = m_dudleyMesh.get();
    switch (arg.getFunctionSpace().getTypeCode()) {
        case Nodes: // fall through
        case DegreesOfFreedom:
        {
            escript::Data temp(arg, escript::function(*this));
            Assemble_integrate(mesh->Nodes, mesh->Elements, temp, integrals);
        }
        break;
        case Elements: // fall through
        case ReducedElements:
            Assemble_integrate(mesh->Nodes,mesh->Elements, arg, integrals);
        break;
        case FaceElements: // fall through
        case ReducedFaceElements:
            Assemble_integrate(mesh->Nodes,mesh->FaceElements, arg, integrals);
        break;
        case Points:
            throw ValueError("Integral of data on points is not supported.");
        break;
        default:
            stringstream ss;
            ss << "setToIntegrals: Dudley does not know anything about "
                "function space type " << arg.getFunctionSpace().getTypeCode();
            throw DudleyException(ss.str());
    }
}

//
// calculates the gradient of arg:
//
void MeshAdapter::setToGradient(escript::Data& grad, const escript::Data& arg) const
{
    if (*arg.getFunctionSpace().getDomain() != *this)
        throw ValueError("setToGradient: Illegal domain of gradient argument");
    if (*grad.getFunctionSpace().getDomain() != *this)
        throw ValueError("setToGradient: Illegal domain of gradient");

    Mesh* mesh = m_dudleyMesh.get();
    escript::Data nodeData;
    if (getMPISize() > 1) {
        if (arg.getFunctionSpace().getTypeCode() == DegreesOfFreedom) {
            nodeData = escript::Data(arg, continuousFunction(*this));
        } else {
            nodeData = arg;
        }
    } else {
        nodeData = arg;
    }
    switch (grad.getFunctionSpace().getTypeCode()) {
        case Nodes:
            throw DudleyException("Gradient at nodes is not supported.");
        break;
        case Elements:
            Assemble_gradient(mesh->Nodes, mesh->Elements, grad, nodeData);
        break;
        case ReducedElements:
            Assemble_gradient(mesh->Nodes, mesh->Elements, grad, nodeData);
        break;
        case FaceElements:
            Assemble_gradient(mesh->Nodes,mesh->FaceElements, grad, nodeData);
        break;
        case ReducedFaceElements:
            Assemble_gradient(mesh->Nodes, mesh->FaceElements, grad, nodeData);
        break;
        case Points:
            throw DudleyException("Gradient at points is not supported.");
        break;
        case DegreesOfFreedom:
            throw DudleyException("Gradient at degrees of freedom is not supported.");
        break;
        default:
            stringstream ss;
            ss << "Gradient: Dudley does not know anything about function space type " << arg.getFunctionSpace().getTypeCode();
            throw DudleyException(ss.str());
    }
}

//
// returns the size of elements
//
void MeshAdapter::setToSize(escript::Data& size) const
{
    Mesh* mesh=m_dudleyMesh.get();
    switch (size.getFunctionSpace().getTypeCode()) {
        case Nodes:
            throw DudleyException("Size of nodes is not supported.");
        break;
        case Elements:
            Assemble_getSize(mesh->Nodes, mesh->Elements, size);
        break;
        case ReducedElements:
            Assemble_getSize(mesh->Nodes, mesh->Elements, size);
        break;
        case FaceElements:
            Assemble_getSize(mesh->Nodes, mesh->FaceElements, size);
        break;
        case ReducedFaceElements:
            Assemble_getSize(mesh->Nodes, mesh->FaceElements, size);
        break;
        case Points:
            throw DudleyException("Size of point elements is not supported.");
        break;
        case DegreesOfFreedom:
            throw DudleyException("Size of degrees of freedom is not supported.");
        break;
        default:
            stringstream ss;
            ss << "setToSize: Dudley does not know anything about function "
                  "space type " << size.getFunctionSpace().getTypeCode();
            throw ValueError(ss.str());
    }
}

//
// sets the location of nodes
//
void MeshAdapter::setNewX(const escript::Data& newX)
{
    if (*newX.getFunctionSpace().getDomain() != *this) 
        throw DudleyException("Illegal domain of new point locations");

    if (newX.getFunctionSpace() == continuousFunction(*this)) {
        m_dudleyMesh->setCoordinates(newX);
    } else {
        throw DudleyException("As of escript version 3.3 - setNewX only "
                "accepts ContinuousFunction arguments. Please interpolate.");
    }
}

bool MeshAdapter::ownSample(int fs_code, index_t id) const
{
#ifdef ESYS_MPI
    if (getMPISize() > 1) {
        Mesh* mesh = m_dudleyMesh.get();
        if (fs_code == DUDLEY_NODES) {
            const index_t myFirstNode = mesh->Nodes->getFirstNode();
            const index_t myLastNode = mesh->Nodes->getLastNode();
            const index_t k = mesh->Nodes->borrowGlobalNodesIndex()[id];
            return (myFirstNode <= k && k < myLastNode);
        } else {
            throw ValueError("ownSample: unsupported function space type");
        }
    }
#endif
    return true;
}

#ifdef ESYS_HAVE_TRILINOS
const_TrilinosGraph_ptr MeshAdapter::getTrilinosGraph() const
{
    if (m_graph.is_null()) {
        m_graph = m_dudleyMesh->createTrilinosGraph();
    }
    return m_graph;
}
#endif

//
// creates a stiffness matrix and initializes it with zeros
//
escript::ASM_ptr MeshAdapter::newSystemMatrix(int row_blocksize,
                            const escript::FunctionSpace& row_functionspace,
                            int column_blocksize,
                            const escript::FunctionSpace& column_functionspace,
                            int type) const
{
    // is the domain right?
    if (*row_functionspace.getDomain() != *this) 
        throw ValueError("domain of row function space does not match the domain of matrix generator.");
    if (*column_functionspace.getDomain() != *this) 
        throw DudleyException("domain of column function space does not match the domain of matrix generator.");

    // is the function space type right?
    if (row_functionspace.getTypeCode() != DegreesOfFreedom) {
        throw DudleyException("illegal function space type for system matrix rows.");
    }
    if (column_functionspace.getTypeCode() != DegreesOfFreedom) {
        throw DudleyException("illegal function space type for system matrix columns.");
    }

    // generate matrix
    if (type & (int)SMT_TRILINOS) {
#ifdef ESYS_HAVE_TRILINOS
        const_TrilinosGraph_ptr graph(getTrilinosGraph());
        escript::ASM_ptr sm(new TrilinosMatrixAdapter(m_dudleyMesh->MPIInfo,
                    row_blocksize, row_functionspace, graph));
        return sm;
#else
        throw DudleyException("newSystemMatrix: dudley was not compiled "
                "with Trilinos support so the Trilinos solver stack cannot be "
                "used.");
#endif
    } else if (type & (int)SMT_PASO) {
        paso::SystemMatrixPattern_ptr pattern(getMesh()->getPasoPattern());
        paso::SystemMatrix_ptr sm(new paso::SystemMatrix(type, pattern,
                  row_blocksize, column_blocksize, false, row_functionspace,
                  column_functionspace));
        return sm;
    } else {
        throw DudleyException("newSystemMatrix: unknown matrix type ID");
    }
}

//
// creates a TransportProblem
//
escript::ATP_ptr MeshAdapter::newTransportProblem(int blocksize,
                                             const escript::FunctionSpace& fs,
                                             int type) const
{
    // is the domain right?
    if (*fs.getDomain() != *this) 
        throw DudleyException("domain of function space does not match the domain of transport problem generator.");
    // is the function space type right 
    if (fs.getTypeCode() != DegreesOfFreedom) {
        throw DudleyException("illegal function space type for system matrix rows.");
    }

    // generate matrix
    paso::SystemMatrixPattern_ptr pattern(getMesh()->getPasoPattern());
    paso::TransportProblem_ptr transportProblem(new paso::TransportProblem(
                                              pattern, blocksize, fs));
    return transportProblem;
}

//
// returns true if data at the atom_type is considered as being cell centered:
bool MeshAdapter::isCellOriented(int functionSpaceCode) const
{
    switch (functionSpaceCode) {
        case Nodes:
        case DegreesOfFreedom:
            return false;
            break;
        case Elements:
        case FaceElements:
        case Points:
        case ReducedElements:
        case ReducedFaceElements:
            return true;
            break;
        default:
            stringstream ss;
            ss << "isCellOriented: Dudley does not know anything about "
                  "function space type " << functionSpaceCode;
            throw ValueError(ss.str());
    }
    return false;
}

bool
MeshAdapter::commonFunctionSpace(const vector<int>& fs, int& resultcode) const
{
    if (fs.empty())
        return false;
    // The idea is to use equivalence classes, i.e. types which can be
    // interpolated back and forth
    //    class 1: DOF <-> Nodes
    //    class 3: Points
    //    class 4: Elements
    //    class 5: ReducedElements
    //    class 6: FaceElements
    //    class 7: ReducedFaceElements

    // There is also a set of lines. Interpolation is possible down a line but
    // not between lines.
    // class 1 and 2 belong to all lines so aren't considered.
    //    line 0: class 3
    //    line 1: class 4,5
    //    line 2: class 6,7

    // For classes with multiple members (class 1) we have vars to record
    // if there is at least one instance -> hasnodes is true if we have at
    // least one instance of Nodes.
    vector<int> hasclass(8);
    vector<int> hasline(3);
    bool hasnodes = false;
    for (int i = 0; i < fs.size(); ++i) {
        switch (fs[i]) {
            case Nodes:
                hasnodes = true; // fall through
            case DegreesOfFreedom:
                hasclass[1] = 1;
                break;
            case Points:
                hasline[0] = 1;
                hasclass[3] = 1;
                break;
            case Elements:
                hasclass[4] = 1;
                hasline[1] = 1;
                break;
            case ReducedElements:
                hasclass[5] = 1;
                hasline[1] = 1;
                break;
            case FaceElements:
                hasclass[6] = 1;
                hasline[2] = 1;
                break;
            case ReducedFaceElements:
                hasclass[7] = 1;
                hasline[2] = 1;
                break;
            default:
                return false;
        }
    }
    int totlines = hasline[0]+hasline[1]+hasline[2];
    // fail if we have more than one leaf group
    if (totlines > 1)
        // there are at least two branches we can't interpolate between
        return false;

    if (totlines == 1) {
        if (hasline[0] == 1) // we have points
            resultcode = Points;
        else if (hasline[1] == 1) {
            if (hasclass[5]==1)
                resultcode=ReducedElements;
            else
                resultcode=Elements;
        } else if (hasline[2]==1) {
            if (hasclass[7]==1)
                resultcode=ReducedFaceElements;
            else
                resultcode=FaceElements;
        }
    } else { // totlines==0
        // something from class 1
        resultcode = (hasnodes ? Nodes : DegreesOfFreedom);
    }
    return true;
}

bool MeshAdapter::probeInterpolationOnDomain(int functionSpaceType_source,
                                             int functionSpaceType_target) const
{
    switch(functionSpaceType_source) {
        case Nodes:
            switch (functionSpaceType_target) {
                case Nodes:
                case DegreesOfFreedom:
                case Elements:
                case ReducedElements:
                case FaceElements:
                case ReducedFaceElements:
                case Points:
                    return true;
                default:
                    stringstream ss;
                    ss << "Interpolation On Domain: Dudley does not know "
                        "anything about function space type "
                        << functionSpaceType_target;
                    throw ValueError(ss.str());
            }
        break;
        case Elements:
            return (functionSpaceType_target == Elements ||
                    functionSpaceType_target == ReducedElements);
        case ReducedElements:
            return (functionSpaceType_target == ReducedElements);
        case FaceElements:
            return (functionSpaceType_target == FaceElements ||
                    functionSpaceType_target == ReducedFaceElements);
        case ReducedFaceElements:
            return (functionSpaceType_target == ReducedFaceElements);
        case Points:
            return (functionSpaceType_target == Points);
        case DegreesOfFreedom:
            switch (functionSpaceType_target) {
                case DegreesOfFreedom:
                case Nodes:
                case Elements:
                case ReducedElements:
                case Points:
                case FaceElements:
                case ReducedFaceElements:
                    return true;
                default:
                    stringstream ss;
                    ss << "Interpolation On Domain: Dudley does not know "
                          "anything about function space type "
                       << functionSpaceType_target;
                    throw DudleyException(ss.str());
            }
            break;
        default:
            stringstream ss;
            ss << "Interpolation On Domain: Dudley does not know anything "
                  "about function space type " << functionSpaceType_source;
            throw DudleyException(ss.str());
    }
    return false;
}

signed char MeshAdapter::preferredInterpolationOnDomain(
        int functionSpaceType_source,int functionSpaceType_target) const
{
    if (probeInterpolationOnDomain(functionSpaceType_source, functionSpaceType_target))
        return 1;
    else if (probeInterpolationOnDomain(functionSpaceType_target, functionSpaceType_source))
        return -1;
    return 0;
}

bool MeshAdapter::probeInterpolationAcross(int functionSpaceType_source,
        const AbstractDomain& targetDomain, int functionSpaceType_target) const
{
    return false;
}

bool MeshAdapter::operator==(const AbstractDomain& other) const
{
    const MeshAdapter* temp = dynamic_cast<const MeshAdapter*>(&other);
    if (temp) {
        return (m_dudleyMesh == temp->m_dudleyMesh);
    }
    return false;
}

bool MeshAdapter::operator!=(const AbstractDomain& other) const
{
    return !(operator==(other));
}

int MeshAdapter::getSystemMatrixTypeId(const bp::object& options) const
{
    const escript::SolverBuddy& sb = bp::extract<escript::SolverBuddy>(options);

    int package = sb.getPackage();
    if (package == escript::SO_PACKAGE_TRILINOS) {
#ifdef ESYS_HAVE_TRILINOS
        return (int)SMT_TRILINOS;
#else
        throw DudleyException("Trilinos requested but not built with Trilinos.");       
#endif
    }
    return (int)SMT_PASO | paso::SystemMatrix::getSystemMatrixTypeId(
                sb.getSolverMethod(), sb.getPreconditioner(), sb.getPackage(),
                sb.isSymmetric(), m_dudleyMesh->MPIInfo);
}

int MeshAdapter::getTransportTypeId(int solver, int preconditioner,
                                    int package, bool symmetry) const
{
    return paso::TransportProblem::getTypeId(solver, preconditioner, package,
                                             symmetry, getMPI());
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

const index_t* MeshAdapter::borrowSampleReferenceIDs(int functionSpaceType) const
{
    index_t* out = NULL;
    switch (functionSpaceType) {
        case Nodes:
            out = getMesh()->Nodes->Id;
        break;
        case Elements:
            out = getMesh()->Elements->Id;
        break;
        case ReducedElements:
            out = getMesh()->Elements->Id;
        break;
        case FaceElements:
            out = getMesh()->FaceElements->Id;
        break;
        case ReducedFaceElements:
            out = getMesh()->FaceElements->Id;
        break;
        case Points:
            out = getMesh()->Points->Id;
        break;
        case DegreesOfFreedom:
            out = getMesh()->Nodes->degreesOfFreedomId;
        break;
        default:
            stringstream ss;
            ss << "Invalid function space type: " << functionSpaceType
               << " for domain: " << getDescription();
            throw ValueError(ss.str());
    }
    return out;
}

int MeshAdapter::getTagFromSampleNo(int functionSpaceType, index_t sampleNo) const
{
    int out = 0;
    switch (functionSpaceType) {
        case Nodes:
            out = getMesh()->Nodes->Tag[sampleNo];
        break;
        case Elements:
            out = getMesh()->Elements->Tag[sampleNo];
        break;
        case ReducedElements:
            out = getMesh()->Elements->Tag[sampleNo];
        break;
        case FaceElements:
            out = getMesh()->FaceElements->Tag[sampleNo];
        break;
        case ReducedFaceElements:
            out = getMesh()->FaceElements->Tag[sampleNo];
        break;
        case Points:
            out = getMesh()->Points->Tag[sampleNo];
        break;
        case DegreesOfFreedom:
            throw DudleyException("DegreesOfFreedom does not support tags.");
        break;
        default:
            stringstream ss;
            ss << "Invalid function space type: " << functionSpaceType
               << " for domain: " << getDescription();
            throw DudleyException(ss.str());
    }
    return out;
}


void MeshAdapter::setTags(int functionSpaceType, int newTag, const escript::Data& mask) const
{
    switch (functionSpaceType) {
        case Nodes:
            getMesh()->Nodes->setTags(newTag, mask);
            break;
        case DegreesOfFreedom:
            throw DudleyException("DegreesOfFreedom does not support tags");
            break;
        case Elements: // fall through
        case ReducedElements:
            getMesh()->Elements->setTags(newTag, mask);
            break;
        case FaceElements:
        case ReducedFaceElements:
            getMesh()->FaceElements->setTags(newTag, mask);
            break;
        case Points:
            getMesh()->Points->setTags(newTag, mask);
            break;
        default:
            stringstream ss;
            ss << "Dudley does not know anything about function space type "
               << functionSpaceType;
            throw ValueError(ss.str());
    }
}

void MeshAdapter::setTagMap(const string& name,  int tag)
{
    getMesh()->addTagMap(name, tag);
}

int MeshAdapter::getTag(const string& name) const
{
    return getMesh()->getTag(name);
}

bool MeshAdapter::isValidTagName(const string& name) const
{
    return getMesh()->isValidTagName(name);
}

string MeshAdapter::showTagNames() const
{
    stringstream ss;
    TagMap::const_iterator it = getMesh()->tagMap.begin();
    while (it != getMesh()->tagMap.end()) {
        ss << it->first;
        ++it;
        if (it != getMesh()->tagMap.end())
            ss << ", ";
    }
    return ss.str();
}

int MeshAdapter::getNumberOfTagsInUse(int functionSpaceCode) const
{
    switch (functionSpaceCode) {
        case Nodes:
            return getMesh()->Nodes->tagsInUse.size();
        case DegreesOfFreedom:
            throw ValueError("DegreesOfFreedom does not support tags");
        case Elements: // fall through
        case ReducedElements:
            return getMesh()->Elements->tagsInUse.size();
        case FaceElements: // fall through
        case ReducedFaceElements:
            return getMesh()->FaceElements->tagsInUse.size();
        case Points:
            return getMesh()->Points->tagsInUse.size();
        default:
            stringstream ss;
            ss << "Dudley does not know anything about function space type "
               << functionSpaceCode;
            throw ValueError(ss.str());
    }
    return 0;
}

const int* MeshAdapter::borrowListOfTagsInUse(int functionSpaceCode) const
{
    switch (functionSpaceCode) {
        case Nodes:
            if (getMesh()->Nodes->tagsInUse.empty())
                return NULL;
            else
                return &getMesh()->Nodes->tagsInUse[0];
        case DegreesOfFreedom:
            throw DudleyException("DegreesOfFreedom does not support tags");
        case Elements: // fall through
        case ReducedElements:
            if (getMesh()->Elements->tagsInUse.empty())
                return NULL;
            else
                return &getMesh()->Elements->tagsInUse[0];
        case FaceElements: // fall through
        case ReducedFaceElements:
            if (getMesh()->FaceElements->tagsInUse.empty())
                return NULL;
            else
                return &getMesh()->FaceElements->tagsInUse[0];
        case Points:
            if (getMesh()->Points->tagsInUse.empty())
                return NULL;
            else
                return &getMesh()->Points->tagsInUse[0];
        default:
            stringstream ss;
            ss << "Dudley does not know anything about function space type "
               << functionSpaceCode;
            throw DudleyException(ss.str());
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
            return true;
        default:
            return false;
    }
}

MeshAdapter::StatusType MeshAdapter::getStatus() const
{
    return getMesh()->getStatus();
}

int MeshAdapter::getApproximationOrder(int functionSpaceCode) const
{
    switch (functionSpaceCode) {
        case Nodes:
        case DegreesOfFreedom:
            return getMesh()->approximationOrder;
        case Elements:
        case FaceElements:
        case Points:
            return getMesh()->integrationOrder;
        case ReducedElements:
        case ReducedFaceElements:
            return getMesh()->reducedIntegrationOrder;
        default:
            stringstream ss;
            ss << "Dudley does not know anything about function space type "
               << functionSpaceCode;
            throw ValueError(ss.str());
    }
    return 0;
}

bool MeshAdapter::supportsContactElements() const
{
    return false;
}

escript::Data MeshAdapter::randomFill(
        const escript::DataTypes::ShapeType& shape,
        const escript::FunctionSpace& what, long seed,
        const bp::tuple& filter) const
{
    escript::Data towipe(0, shape, what, true);
    // since we just made this object, no sharing is possible and we don't
    // need to check for exclusive write
    escript::DataTypes::RealVectorType& dv(towipe.getExpandedVectorReference());
    escript::randomFillArray(seed, &dv[0], dv.size());
    return towipe;       
}

} // end of namespace

