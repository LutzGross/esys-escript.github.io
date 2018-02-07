
/*****************************************************************************
*
* Copyright (c) 2003-2018 by The University of Queensland
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

#include <dudley/DudleyDomain.h>
#include <dudley/Assemble.h>
#include <dudley/DudleyException.h>
#include <dudley/IndexList.h>

#include <escript/Data.h>
#include <escript/DataFactory.h>
#include <escript/Random.h>
#include <escript/SolverOptions.h>

#ifdef ESYS_HAVE_PASO
#include <paso/SystemMatrix.h>
#include <paso/Transport.h>
#endif

#ifdef ESYS_HAVE_TRILINOS
#include <trilinoswrap/TrilinosMatrixAdapter.h>

using esys_trilinos::TrilinosMatrixAdapter;
using esys_trilinos::const_TrilinosGraph_ptr;
#endif

#include <boost/scoped_array.hpp>

#ifdef ESYS_HAVE_NETCDF
 #ifdef NETCDF4
  #include <ncVar.h>
  #include <ncDim.h>
  #include <escript/NCHelper.h>

 #else
   #include <netcdfcpp.h>
 #endif
#endif

using namespace std;
namespace bp = boost::python;
using escript::NotImplementedError;
using escript::ValueError;

#ifdef NETCDF4
using namespace netCDF;
#endif

namespace dudley {

using escript::DataTypes::real_t;
using escript::DataTypes::cplx_t;

DudleyDomain::FunctionSpaceNamesMapType DudleyDomain::m_functionSpaceTypeNames;

DudleyDomain::DudleyDomain(const string& name, int numDim, escript::JMPI jmpi) :
    m_mpiInfo(jmpi),
    m_name(name),
    m_elements(NULL),
    m_faceElements(NULL),
    m_points(NULL)
{
    // allocate node table
    m_nodes = new NodeFile(numDim, m_mpiInfo);
    setFunctionSpaceTypeNames();
}

DudleyDomain::DudleyDomain(const DudleyDomain& in) :
    m_mpiInfo(in.m_mpiInfo),
    m_name(in.m_name),
    m_nodes(in.m_nodes),
    m_elements(in.m_elements),
    m_faceElements(in.m_faceElements),
    m_points(in.m_points)
{
    setFunctionSpaceTypeNames();
}

DudleyDomain::~DudleyDomain()
{
    delete m_nodes;
    delete m_elements;
    delete m_faceElements;
    delete m_points;
}

void DudleyDomain::MPIBarrier() const
{
#ifdef ESYS_MPI
    MPI_Barrier(getMPIComm());
#endif
}

void DudleyDomain::setElements(ElementFile* elements)
{
    delete m_elements;
    m_elements = elements;
}

void DudleyDomain::setFaceElements(ElementFile* elements)
{
    delete m_faceElements;
    m_faceElements = elements;
}

void DudleyDomain::setPoints(ElementFile* elements)
{
    delete m_points;
    m_points = elements;
}

void DudleyDomain::createMappings(const IndexVector& dofDist,
                                  const IndexVector& nodeDist)
{
    m_nodes->createNodeMappings(dofDist, nodeDist);

#ifdef ESYS_HAVE_TRILINOS
    // TODO?: the following block should probably go into prepare() but
    // Domain::load() only calls createMappings which is why it's here...
    // make sure trilinos distribution graph is available for matrix building
    // and interpolation
    const index_t numTargets = m_nodes->getNumDegreesOfFreedomTargets();
    const index_t* target = m_nodes->borrowTargetDegreesOfFreedom();
    boost::scoped_array<IndexList> indexList(new IndexList[numTargets]);

#pragma omp parallel
    {
        // insert contributions from element matrices into columns in
        // index list
        IndexList_insertElements(indexList.get(), m_elements, target);
        IndexList_insertElements(indexList.get(), m_faceElements, target);
        IndexList_insertElements(indexList.get(), m_points, target);
    }
    m_nodes->createTrilinosGraph(indexList.get());
#endif
}

void DudleyDomain::markNodes(vector<short>& mask, index_t offset) const
{
    m_elements->markNodes(mask, offset);
    m_faceElements->markNodes(mask, offset);
    m_points->markNodes(mask, offset);
}

void DudleyDomain::relabelElementNodes(const index_t* newNode, index_t offset)
{
    m_elements->relabelNodes(newNode, offset);
    m_faceElements->relabelNodes(newNode, offset);
    m_points->relabelNodes(newNode, offset);
}

#ifdef NETCDF4

void DudleyDomain::dump(const string& fileName) const
{
#ifdef ESYS_HAVE_NETCDF
    NcDim ncdims[12];
    NcVar ids;
    index_t* index_ptr;
#ifdef ESYS_INDEXTYPE_LONG
    NcType ncIdxType = ncLong;
#else
    NcType ncIdxType = ncInt;
#endif
    int num_Tags = 0;
    int mpi_size                  = getMPISize();
    int mpi_rank                  = getMPIRank();
    int numDim                    = m_nodes->numDim;
    dim_t numNodes                = m_nodes->getNumNodes();
    dim_t num_Elements            = m_elements->numElements;
    dim_t num_FaceElements        = m_faceElements->numElements;
    dim_t num_Points              = m_points->numElements;
    int num_Elements_numNodes     = m_elements->numNodes;
    int num_FaceElements_numNodes = m_faceElements->numNodes;
#ifdef ESYS_MPI
    MPI_Status status;
#endif

    // Incoming token indicates it's my turn to write
#ifdef ESYS_MPI
    if (mpi_rank > 0)
        MPI_Recv(&num_Tags, 0, MPI_INT, mpi_rank-1, 81800, getMPIComm(), &status);
#endif

    const string newFileName(m_mpiInfo->appendRankToFileName(fileName));

    // Figure out how much storage is required for tags
    num_Tags = m_tagMap.size();

    // Create the file
    NcFile dataFile;
    try
    {
        dataFile.open(newFileName.c_str(), NcFile::FileMode::replace,   NcFile::FileFormat::classic64);
    }
    catch (exceptions::NcException e)
    {
        throw DudleyException("Error - DudleyDomain:: opening of netCDF file for output failed.");
    }        

    string msgPrefix("Error in DudleyDomain::dump: NetCDF operation failed - ");
    // Define dimensions (num_Elements and dim_Elements are identical,
    // dim_Elements only appears if > 0)
    if ((ncdims[0] = dataFile.addDim("numNodes", numNodes)).isNull() )
        throw DudleyException(msgPrefix+"add_dim(numNodes)");
    if ((ncdims[1] = dataFile.addDim("numDim", numDim)).isNull() )
        throw DudleyException(msgPrefix+"add_dim(numDim)");
    if ((ncdims[2] = dataFile.addDim("mpi_size_plus_1", mpi_size+1)).isNull() )
        throw DudleyException(msgPrefix+"add_dim(mpi_size)");
    if (num_Elements > 0)
        if ((ncdims[3] = dataFile.addDim("dim_Elements", num_Elements)).isNull() )
            throw DudleyException(msgPrefix+"add_dim(dim_Elements)");
    if (num_FaceElements > 0)
        if ((ncdims[4] = dataFile.addDim("dim_FaceElements", num_FaceElements)).isNull() )
         throw DudleyException(msgPrefix+"add_dim(dim_FaceElements)");
    if (num_Points > 0)
        if ((ncdims[6] = dataFile.addDim("dim_Points", num_Points)).isNull() )
            throw DudleyException(msgPrefix+"add_dim(dim_Points)");
    if (num_Elements > 0)
        if ((ncdims[7] = dataFile.addDim("dim_Elements_Nodes", num_Elements_numNodes)).isNull() )
            throw DudleyException(msgPrefix+"add_dim(dim_Elements_Nodes)");
    if (num_FaceElements > 0)
        if ((ncdims[8] = dataFile.addDim("dim_FaceElements_numNodes", num_FaceElements_numNodes)).isNull() )
            throw DudleyException(msgPrefix+"add_dim(dim_FaceElements_numNodes)");
    if (num_Tags > 0)
        if ((ncdims[10] = dataFile.addDim("dim_Tags", num_Tags)).isNull() )
            throw DudleyException(msgPrefix+"add_dim(dim_Tags)");

    // Attributes: MPI size, MPI rank, Name, order, reduced_order
    NcInt ni;        
    if (dataFile.putAtt("index_size", ni, (int)sizeof(index_t)).isNull())
        throw DudleyException(msgPrefix+"add_att(index_size)");
    if (dataFile.putAtt("mpi_size", ni, mpi_size).isNull())
        throw DudleyException(msgPrefix+"add_att(mpi_size)");
    if (dataFile.putAtt("mpi_rank", ni, mpi_rank).isNull())
        throw DudleyException(msgPrefix+"add_att(mpi_rank)");
    if (dataFile.putAtt("Name", m_name.c_str()).isNull())
        throw DudleyException(msgPrefix+"add_att(Name)");
    if (dataFile.putAtt("numDim", ni, numDim).isNull())
        throw DudleyException(msgPrefix+"add_att(numDim)");
    if (dataFile.putAtt("order", ni, 2).isNull())
        throw DudleyException(msgPrefix+"add_att(order)");
    if (dataFile.putAtt("reduced_order", ni, 0).isNull())
        throw DudleyException(msgPrefix+"add_att(reduced_order)");
    if (dataFile.putAtt("numNodes", ni, numNodes).isNull())
        throw DudleyException(msgPrefix+"add_att(numNodes)");
    if (dataFile.putAtt("num_Elements", ni, num_Elements).isNull())
        throw DudleyException(msgPrefix+"add_att(num_Elements)");
    if (dataFile.putAtt("num_FaceElements", ni, num_FaceElements).isNull())
        throw DudleyException(msgPrefix+"add_att(num_FaceElements)");
    if (dataFile.putAtt("num_Points", ni, num_Points).isNull())
        throw DudleyException(msgPrefix+"add_att(num_Points)");
    if (dataFile.putAtt("num_Elements_numNodes", ni, num_Elements_numNodes).isNull())
        throw DudleyException(msgPrefix+"add_att(num_Elements_numNodes)");
    if (dataFile.putAtt("num_FaceElements_numNodes", ni, num_FaceElements_numNodes).isNull())
        throw DudleyException(msgPrefix+"add_att(num_FaceElements_numNodes)");
    if (dataFile.putAtt("Elements_TypeId", ni, m_elements->etype).isNull())
        throw DudleyException(msgPrefix+"add_att(Elements_TypeId)");
    if (dataFile.putAtt("FaceElements_TypeId", ni, m_faceElements->etype).isNull())
        throw DudleyException(msgPrefix+"add_att(FaceElements_TypeId)");
    if (dataFile.putAtt("Points_TypeId", ni, m_points->etype).isNull())
        throw DudleyException(msgPrefix+"add_att(Points_TypeId)");
    if (dataFile.putAtt("num_Tags", ni, num_Tags).isNull())
        throw DudleyException(msgPrefix+"add_att(num_Tags)");

    // // // // // Nodes // // // // //

    // Nodes nodeDistribution
    if ((ids = dataFile.addVar("Nodes_NodeDistribution", ncIdxType, ncdims[2])).isNull() )
        throw DudleyException(msgPrefix+"add_var(Nodes_NodeDistribution)");
    index_ptr = &m_nodes->nodesDistribution->first_component[0];
    ids.putVar(index_ptr);  // mpi_size+1

    // Nodes degreesOfFreedomDistribution
    if (( ids = dataFile.addVar("Nodes_DofDistribution", ncIdxType, ncdims[2])).isNull() )
        throw DudleyException(msgPrefix+"add_var(Nodes_DofDistribution)");
    index_ptr = &m_nodes->dofDistribution->first_component[0];
    ids.putVar(index_ptr);  // mpi_size+1

    // Only write nodes if non-empty because NetCDF doesn't like empty arrays
    // (it treats them as NC_UNLIMITED)
    if (numNodes > 0) {
        // Nodes Id
        if (( ids = dataFile.addVar("Nodes_Id", ncIdxType, ncdims[0])).isNull() )
            throw DudleyException(msgPrefix+"add_var(Nodes_Id)");
        ids.putVar(m_nodes->Id);    // numNodes

        // Nodes Tag
        if (( ids = dataFile.addVar("Nodes_Tag", ncInt, ncdims[0])).isNull() )
            throw DudleyException(msgPrefix+"add_var(Nodes_Tag)");
        ids.putVar(m_nodes->Tag);  // numNodes

        // Nodes gDOF
        if (( ids = dataFile.addVar("Nodes_gDOF", ncIdxType, ncdims[0])).isNull() )
            throw DudleyException(msgPrefix+"add_var(Nodes_gDOF)");
        ids.putVar(m_nodes->globalDegreesOfFreedom);  // numNodes

        // Nodes global node index
        if (( ids = dataFile.addVar("Nodes_gNI", ncIdxType, ncdims[0])).isNull() )
            throw DudleyException(msgPrefix+"add_var(Nodes_gNI)");
        ids.putVar(m_nodes->globalNodesIndex); // numNodes

        // Nodes Coordinates
        vector<NcDim> ncds;
        ncds.push_back(ncdims[0]);
        ncds.push_back(ncdims[1]);
        if (( ids = dataFile.addVar("Nodes_Coordinates", ncDouble, ncds)).isNull() )
            throw DudleyException(msgPrefix+"add_var(Nodes_Coordinates)");
        ids.putVar(m_nodes->Coordinates);  // should be (numNodes, numDim) values written        
        
        
    }

    // // // // // Elements // // // // //
    if (num_Elements > 0) {
        // Elements_Id
        if (( ids = dataFile.addVar("Elements_Id", ncIdxType, ncdims[3])).isNull() )
            throw DudleyException(msgPrefix+"add_var(Elements_Id)");
        ids.putVar(m_elements->Id);    // num_Elements

        // Elements_Tag
        if (( ids = dataFile.addVar("Elements_Tag", ncInt, ncdims[3])).isNull() )
            throw DudleyException(msgPrefix+"add_var(Elements_Tag)");
        ids.putVar(m_elements->Tag);   // num_Elements

        // Elements_Owner
        if (( ids = dataFile.addVar("Elements_Owner", ncInt, ncdims[3])).isNull() )
            throw DudleyException(msgPrefix+"add_var(Elements_Owner)");
        ids.putVar(m_elements->Owner); // num_Elements

        // Elements_Color
        if (( ids = dataFile.addVar("Elements_Color", ncIdxType, ncdims[3])).isNull() )
            throw DudleyException(msgPrefix+"add_var(Elements_Color)");
        ids.putVar(m_elements->Color); // num_Elements

        // Elements_Nodes
        vector<NcDim> dv;
        dv.push_back(ncdims[3]);
        dv.push_back(ncdims[7]);
        if (( ids = dataFile.addVar("Elements_Nodes", ncIdxType, dv) ).isNull() )
            throw DudleyException(msgPrefix+"add_var(Elements_Nodes)");
        ids.putVar(m_elements->Nodes); //(, num_Elements, num_Elements_numNodes) values written        
        
        
    }

    // // // // // Face_Elements // // // // //
    if (num_FaceElements > 0) {
        // FaceElements_Id
        if ((ids = dataFile.addVar("FaceElements_Id", ncIdxType, ncdims[4])).isNull())
            throw DudleyException(msgPrefix+"add_var(FaceElements_Id)");
        ids.putVar(m_faceElements->Id); // num_FaceElements

        // FaceElements_Tag
        if ((ids = dataFile.addVar("FaceElements_Tag", ncInt, ncdims[4])).isNull())
            throw DudleyException(msgPrefix+"add_var(FaceElements_Tag)");
        ids.putVar(m_faceElements->Tag);    // num_FaceElements

        // FaceElements_Owner
        if ((ids = dataFile.addVar("FaceElements_Owner", ncInt, ncdims[4])).isNull())
            throw DudleyException(msgPrefix+"add_var(FaceElements_Owner)");
        ids.putVar(m_faceElements->Owner);  // num_FaceElements

        // FaceElements_Color
        if ((ids = dataFile.addVar("FaceElements_Color", ncIdxType, ncdims[4])).isNull())
            throw DudleyException(msgPrefix+"add_var(FaceElements_Color)");
        ids.putVar(m_faceElements->Color);  // num_FaceElements

        // FaceElements_Nodes
        vector<NcDim> dv;
        dv.push_back(ncdims[4]);
        dv.push_back(ncdims[8]);
        if ((ids = dataFile.addVar("FaceElements_Nodes", ncIdxType, dv)).isNull())
            throw DudleyException(msgPrefix+"add_var(FaceElements_Nodes)");
        ids.putVar(m_faceElements->Nodes);  // num_FaceElements, num_FaceElements_numNodes
    }

    // // // // // Points // // // // //
    if (num_Points > 0) {
        // Points_Id
        if ((ids = dataFile.addVar("Points_Id", ncIdxType, ncdims[6])).isNull())
            throw DudleyException(msgPrefix+"add_var(Points_Id)");
        ids.putVar(m_points->Id);   // num_Points

        // Points_Tag
        if ((ids = dataFile.addVar("Points_Tag", ncInt, ncdims[6])).isNull())
            throw DudleyException(msgPrefix+"add_var(Points_Tag)");
        ids.putVar(m_points->Tag);  // num_Points

        // Points_Owner
        if ((ids = dataFile.addVar("Points_Owner", ncInt, ncdims[6])).isNull())
            throw DudleyException(msgPrefix+"add_var(Points_Owner)");
        ids.putVar(m_points->Owner);    // num_Points

        // Points_Color
        if ((ids = dataFile.addVar("Points_Color", ncIdxType, ncdims[6])).isNull())
            throw DudleyException(msgPrefix+"add_var(Points_Color)");
        ids.putVar(m_points->Color);    // num_Points

        // Points_Nodes
        if ((ids = dataFile.addVar("Points_Nodes", ncIdxType, ncdims[6])).isNull())
            throw DudleyException(msgPrefix+"add_var(Points_Nodes)");
        ids.putVar(m_points->Nodes);    // num_Points
    }

    // // // // // TagMap // // // // //
    if (num_Tags > 0) {
        // Temp storage to gather node IDs
        vector<int> Tags_keys;

        // Copy tag data into temp arrays
        TagMap::const_iterator it;
        for (it = m_tagMap.begin(); it != m_tagMap.end(); it++) {
            Tags_keys.push_back(it->second);
        }

        // Tags_keys
        if ((ids = dataFile.addVar("Tags_keys", ncInt, ncdims[10])).isNull())
            throw DudleyException(msgPrefix+"add_var(Tags_keys)");
        ids.putVar(&Tags_keys[0]);  // num_Tags

        // Tags_names_*
        // This is an array of strings, it should be stored as an array but
        // instead I have hacked in one attribute per string because the NetCDF
        // manual doesn't tell how to do an array of strings
        int i = 0;
        for (it = m_tagMap.begin(); it != m_tagMap.end(); it++, i++) {
            stringstream ss;
            ss << "Tags_name_" << i;
            const string name(ss.str());
            if (dataFile.putAtt(name.c_str(), it->first.c_str()).isNull())
                throw DudleyException(msgPrefix+"add_att(Tags_names_X)");
        }
    }

    // Send token to next MPI process so he can take his turn
#ifdef ESYS_MPI
    if (mpi_rank < mpi_size-1)
        MPI_Send(&num_Tags, 0, MPI_INT, mpi_rank+1, 81800, getMPIComm());
#endif

    // NetCDF file is closed by destructor of NcFile object

#else
    throw DudleyException("DudleyDomain::dump: not configured with netCDF. "
                          "Please contact your installation manager.");
#endif // ESYS_HAVE_NETCDF
}

#else

void DudleyDomain::dump(const string& fileName) const
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
    int num_Tags = 0;
    int mpi_size                  = getMPISize();
    int mpi_rank                  = getMPIRank();
    int numDim                    = m_nodes->numDim;
    dim_t numNodes                = m_nodes->getNumNodes();
    dim_t num_Elements            = m_elements->numElements;
    dim_t num_FaceElements        = m_faceElements->numElements;
    dim_t num_Points              = m_points->numElements;
    int num_Elements_numNodes     = m_elements->numNodes;
    int num_FaceElements_numNodes = m_faceElements->numNodes;
#ifdef ESYS_MPI
    MPI_Status status;
#endif

    // Incoming token indicates it's my turn to write
#ifdef ESYS_MPI
    if (mpi_rank > 0)
        MPI_Recv(&num_Tags, 0, MPI_INT, mpi_rank-1, 81800, getMPIComm(), &status);
#endif

    const string newFileName(m_mpiInfo->appendRankToFileName(fileName));

    // Figure out how much storage is required for tags
    num_Tags = m_tagMap.size();

    // NetCDF error handler
    NcError err(NcError::verbose_nonfatal);
    // Create the file
    NcFile dataFile(newFileName.c_str(), NcFile::Replace);
    string msgPrefix("Error in DudleyDomain::dump: NetCDF operation failed - ");
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
    if (!dataFile.add_att("mpi_size", mpi_size))
        throw DudleyException(msgPrefix+"add_att(mpi_size)");
    if (!dataFile.add_att("mpi_rank", mpi_rank))
        throw DudleyException(msgPrefix+"add_att(mpi_rank)");
    if (!dataFile.add_att("Name", m_name.c_str()))
        throw DudleyException(msgPrefix+"add_att(Name)");
    if (!dataFile.add_att("numDim", numDim))
        throw DudleyException(msgPrefix+"add_att(numDim)");
    if (!dataFile.add_att("order", 2))
        throw DudleyException(msgPrefix+"add_att(order)");
    if (!dataFile.add_att("reduced_order", 0))
        throw DudleyException(msgPrefix+"add_att(reduced_order)");
    if (!dataFile.add_att("numNodes", numNodes))
        throw DudleyException(msgPrefix+"add_att(numNodes)");
    if (!dataFile.add_att("num_Elements", num_Elements))
        throw DudleyException(msgPrefix+"add_att(num_Elements)");
    if (!dataFile.add_att("num_FaceElements", num_FaceElements))
        throw DudleyException(msgPrefix+"add_att(num_FaceElements)");
    if (!dataFile.add_att("num_Points", num_Points))
        throw DudleyException(msgPrefix+"add_att(num_Points)");
    if (!dataFile.add_att("num_Elements_numNodes", num_Elements_numNodes))
        throw DudleyException(msgPrefix+"add_att(num_Elements_numNodes)");
    if (!dataFile.add_att("num_FaceElements_numNodes", num_FaceElements_numNodes))
        throw DudleyException(msgPrefix+"add_att(num_FaceElements_numNodes)");
    if (!dataFile.add_att("Elements_TypeId", m_elements->etype))
        throw DudleyException(msgPrefix+"add_att(Elements_TypeId)");
    if (!dataFile.add_att("FaceElements_TypeId", m_faceElements->etype))
        throw DudleyException(msgPrefix+"add_att(FaceElements_TypeId)");
    if (!dataFile.add_att("Points_TypeId", m_points->etype))
        throw DudleyException(msgPrefix+"add_att(Points_TypeId)");
    if (!dataFile.add_att("num_Tags", num_Tags))
        throw DudleyException(msgPrefix+"add_att(num_Tags)");

    // // // // // Nodes // // // // //

    // Nodes nodeDistribution
    if (! (ids = dataFile.add_var("Nodes_NodeDistribution", ncIdxType, ncdims[2])) )
        throw DudleyException(msgPrefix+"add_var(Nodes_NodeDistribution)");
    index_ptr = &m_nodes->nodesDistribution->first_component[0];
    if (! (ids->put(index_ptr, mpi_size+1)) )
        throw DudleyException(msgPrefix+"put(Nodes_NodeDistribution)");

    // Nodes degreesOfFreedomDistribution
    if (! ( ids = dataFile.add_var("Nodes_DofDistribution", ncIdxType, ncdims[2])) )
        throw DudleyException(msgPrefix+"add_var(Nodes_DofDistribution)");
    index_ptr = &m_nodes->dofDistribution->first_component[0];
    if (! (ids->put(index_ptr, mpi_size+1)) )
        throw DudleyException(msgPrefix+"put(Nodes_DofDistribution)");

    // Only write nodes if non-empty because NetCDF doesn't like empty arrays
    // (it treats them as NC_UNLIMITED)
    if (numNodes > 0) {
        // Nodes Id
        if (! ( ids = dataFile.add_var("Nodes_Id", ncIdxType, ncdims[0])) )
            throw DudleyException(msgPrefix+"add_var(Nodes_Id)");
        if (! (ids->put(m_nodes->Id, numNodes)) )
            throw DudleyException(msgPrefix+"put(Nodes_Id)");

        // Nodes Tag
        if (! ( ids = dataFile.add_var("Nodes_Tag", ncInt, ncdims[0])) )
            throw DudleyException(msgPrefix+"add_var(Nodes_Tag)");
        if (! (ids->put(m_nodes->Tag, numNodes)) )
            throw DudleyException(msgPrefix+"put(Nodes_Tag)");

        // Nodes gDOF
        if (! ( ids = dataFile.add_var("Nodes_gDOF", ncIdxType, ncdims[0])) )
            throw DudleyException(msgPrefix+"add_var(Nodes_gDOF)");
        if (! (ids->put(m_nodes->globalDegreesOfFreedom, numNodes)) )
            throw DudleyException(msgPrefix+"put(Nodes_gDOF)");

        // Nodes global node index
        if (! ( ids = dataFile.add_var("Nodes_gNI", ncIdxType, ncdims[0])) )
            throw DudleyException(msgPrefix+"add_var(Nodes_gNI)");
        if (! (ids->put(m_nodes->globalNodesIndex, numNodes)) )
            throw DudleyException(msgPrefix+"put(Nodes_gNI)");

        // Nodes Coordinates
        if (! ( ids = dataFile.add_var("Nodes_Coordinates", ncDouble, ncdims[0], ncdims[1]) ) )
            throw DudleyException(msgPrefix+"add_var(Nodes_Coordinates)");
        if (! (ids->put(m_nodes->Coordinates, numNodes, numDim)) )
            throw DudleyException(msgPrefix+"put(Nodes_Coordinates)");
    }

    // // // // // Elements // // // // //
    if (num_Elements > 0) {
        // Elements_Id
        if (! ( ids = dataFile.add_var("Elements_Id", ncIdxType, ncdims[3])) )
            throw DudleyException(msgPrefix+"add_var(Elements_Id)");
        if (! (ids->put(m_elements->Id, num_Elements)) )
            throw DudleyException(msgPrefix+"put(Elements_Id)");

        // Elements_Tag
        if (! ( ids = dataFile.add_var("Elements_Tag", ncInt, ncdims[3])) )
            throw DudleyException(msgPrefix+"add_var(Elements_Tag)");
        if (! (ids->put(m_elements->Tag, num_Elements)) )
            throw DudleyException(msgPrefix+"put(Elements_Tag)");

        // Elements_Owner
        if (! ( ids = dataFile.add_var("Elements_Owner", ncInt, ncdims[3])) )
            throw DudleyException(msgPrefix+"add_var(Elements_Owner)");
        if (! (ids->put(m_elements->Owner, num_Elements)) )
            throw DudleyException(msgPrefix+"put(Elements_Owner)");

        // Elements_Color
        if (! ( ids = dataFile.add_var("Elements_Color", ncIdxType, ncdims[3])) )
            throw DudleyException(msgPrefix+"add_var(Elements_Color)");
        if (! (ids->put(m_elements->Color, num_Elements)) )
            throw DudleyException(msgPrefix+"put(Elements_Color)");

        // Elements_Nodes
        if (! ( ids = dataFile.add_var("Elements_Nodes", ncIdxType, ncdims[3], ncdims[7]) ) )
            throw DudleyException(msgPrefix+"add_var(Elements_Nodes)");
        if (! (ids->put(m_elements->Nodes, num_Elements, num_Elements_numNodes)) )
            throw DudleyException(msgPrefix+"put(Elements_Nodes)");
    }

    // // // // // Face_Elements // // // // //
    if (num_FaceElements > 0) {
        // FaceElements_Id
        if (!(ids = dataFile.add_var("FaceElements_Id", ncIdxType, ncdims[4])))
            throw DudleyException(msgPrefix+"add_var(FaceElements_Id)");
        if (!(ids->put(m_faceElements->Id, num_FaceElements)))
            throw DudleyException(msgPrefix+"put(FaceElements_Id)");

        // FaceElements_Tag
        if (!(ids = dataFile.add_var("FaceElements_Tag", ncInt, ncdims[4])))
            throw DudleyException(msgPrefix+"add_var(FaceElements_Tag)");
        if (!(ids->put(m_faceElements->Tag, num_FaceElements)))
            throw DudleyException(msgPrefix+"put(FaceElements_Tag)");

        // FaceElements_Owner
        if (!(ids = dataFile.add_var("FaceElements_Owner", ncInt, ncdims[4])))
            throw DudleyException(msgPrefix+"add_var(FaceElements_Owner)");
        if (!(ids->put(m_faceElements->Owner, num_FaceElements)))
            throw DudleyException(msgPrefix+"put(FaceElements_Owner)");

        // FaceElements_Color
        if (!(ids = dataFile.add_var("FaceElements_Color", ncIdxType, ncdims[4])))
            throw DudleyException(msgPrefix+"add_var(FaceElements_Color)");
        if (!(ids->put(m_faceElements->Color, num_FaceElements)))
            throw DudleyException(msgPrefix+"put(FaceElements_Color)");

        // FaceElements_Nodes
        if (!(ids = dataFile.add_var("FaceElements_Nodes", ncIdxType, ncdims[4], ncdims[8])))
            throw DudleyException(msgPrefix+"add_var(FaceElements_Nodes)");
        if (!(ids->put(m_faceElements->Nodes, num_FaceElements, num_FaceElements_numNodes)))
            throw DudleyException(msgPrefix+"put(FaceElements_Nodes)");
    }

    // // // // // Points // // // // //
    if (num_Points > 0) {
        // Points_Id
        if (!(ids = dataFile.add_var("Points_Id", ncIdxType, ncdims[6])))
            throw DudleyException(msgPrefix+"add_var(Points_Id)");
        if (!(ids->put(m_points->Id, num_Points)))
            throw DudleyException(msgPrefix+"put(Points_Id)");

        // Points_Tag
        if (!(ids = dataFile.add_var("Points_Tag", ncInt, ncdims[6])))
            throw DudleyException(msgPrefix+"add_var(Points_Tag)");
        if (!(ids->put(m_points->Tag, num_Points)))
            throw DudleyException(msgPrefix+"put(Points_Tag)");

        // Points_Owner
        if (!(ids = dataFile.add_var("Points_Owner", ncInt, ncdims[6])))
            throw DudleyException(msgPrefix+"add_var(Points_Owner)");
        if (!(ids->put(m_points->Owner, num_Points)))
            throw DudleyException(msgPrefix+"put(Points_Owner)");

        // Points_Color
        if (!(ids = dataFile.add_var("Points_Color", ncIdxType, ncdims[6])))
            throw DudleyException(msgPrefix+"add_var(Points_Color)");
        if (!(ids->put(m_points->Color, num_Points)))
            throw DudleyException(msgPrefix+"put(Points_Color)");

        // Points_Nodes
        if (!(ids = dataFile.add_var("Points_Nodes", ncIdxType, ncdims[6])))
            throw DudleyException(msgPrefix+"add_var(Points_Nodes)");
        if (!(ids->put(m_points->Nodes, num_Points)))
            throw DudleyException(msgPrefix+"put(Points_Nodes)");
    }

    // // // // // TagMap // // // // //
    if (num_Tags > 0) {
        // Temp storage to gather node IDs
        vector<int> Tags_keys;

        // Copy tag data into temp arrays
        TagMap::const_iterator it;
        for (it = m_tagMap.begin(); it != m_tagMap.end(); it++) {
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
        for (it = m_tagMap.begin(); it != m_tagMap.end(); it++, i++) {
            stringstream ss;
            ss << "Tags_name_" << i;
            const string name(ss.str());
            if (!dataFile.add_att(name.c_str(), it->first.c_str()))
                throw DudleyException(msgPrefix+"add_att(Tags_names_X)");
        }
    }

    // Send token to next MPI process so he can take his turn
#ifdef ESYS_MPI
    if (mpi_rank < mpi_size-1)
        MPI_Send(&num_Tags, 0, MPI_INT, mpi_rank+1, 81800, getMPIComm());
#endif

    // NetCDF file is closed by destructor of NcFile object

#else
    throw DudleyException("DudleyDomain::dump: not configured with netCDF. "
                          "Please contact your installation manager.");
#endif // ESYS_HAVE_NETCDF
}

#endif

string DudleyDomain::getDescription() const
{
    return "DudleyMesh";
}

string DudleyDomain::functionSpaceTypeAsString(int functionSpaceType) const
{
    FunctionSpaceNamesMapType::iterator loc;
    loc = m_functionSpaceTypeNames.find(functionSpaceType);
    if (loc == m_functionSpaceTypeNames.end()) {
        return "Invalid function space type code.";
    } else {
        return loc->second;
    }
}

bool DudleyDomain::isValidFunctionSpaceType(int functionSpaceType) const
{
    FunctionSpaceNamesMapType::iterator loc;
    loc = m_functionSpaceTypeNames.find(functionSpaceType);
    return (loc != m_functionSpaceTypeNames.end());
}

void DudleyDomain::setFunctionSpaceTypeNames()
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

int DudleyDomain::getContinuousFunctionCode() const
{
    return Nodes;
}

int DudleyDomain::getReducedContinuousFunctionCode() const
{
    return Nodes;
}

int DudleyDomain::getFunctionCode() const
{
    return Elements;
}

int DudleyDomain::getReducedFunctionCode() const
{
    return ReducedElements;
}

int DudleyDomain::getFunctionOnBoundaryCode() const
{
    return FaceElements;
}

int DudleyDomain::getReducedFunctionOnBoundaryCode() const
{
    return ReducedFaceElements;
}

int DudleyDomain::getFunctionOnContactZeroCode() const
{
    throw DudleyException("Dudley does not support contact elements.");
}

int DudleyDomain::getReducedFunctionOnContactZeroCode() const
{
    throw DudleyException("Dudley does not support contact elements.");
}

int DudleyDomain::getFunctionOnContactOneCode() const
{
    throw DudleyException("Dudley does not support contact elements.");
}

int DudleyDomain::getReducedFunctionOnContactOneCode() const
{
    throw DudleyException("Dudley does not support contact elements.");
}

int DudleyDomain::getSolutionCode() const
{
    return DegreesOfFreedom;
}

int DudleyDomain::getReducedSolutionCode() const
{
    return DegreesOfFreedom;
}

int DudleyDomain::getDiracDeltaFunctionsCode() const
{
    return Points;
}

//
// Return the number of data points summed across all MPI processes
//
dim_t DudleyDomain::getNumDataPointsGlobal() const
{
    return m_nodes->getGlobalNumNodes();
}

//
// return the number of data points per sample and the number of samples
// needed to represent data on a parts of the mesh.
//
pair<int,dim_t> DudleyDomain::getDataShape(int functionSpaceCode) const
{
    int numDataPointsPerSample = 0;
    dim_t numSamples = 0;
    switch (functionSpaceCode) {
        case Nodes:
            numDataPointsPerSample = 1;
            numSamples = m_nodes->getNumNodes();
        break;
        case Elements:
            if (m_elements) {
                numSamples = m_elements->numElements;
                numDataPointsPerSample = m_elements->numLocalDim + 1;
            }
        break;
        case ReducedElements:
            if (m_elements) {
                numSamples = m_elements->numElements;
                numDataPointsPerSample = (m_elements->numLocalDim==0) ? 0 : 1;
            }
        break;
        case FaceElements:
            if (m_faceElements) {
                numSamples = m_faceElements->numElements;
                numDataPointsPerSample = m_faceElements->numLocalDim+1;
            }
        break;
        case ReducedFaceElements:
            if (m_faceElements) {
                numSamples = m_faceElements->numElements;
                numDataPointsPerSample = (m_faceElements->numLocalDim==0)? 0:1;
            }
        break;
        case Points:
            if (m_points) {
                numSamples = m_points->numElements;
                numDataPointsPerSample = 1;
            }
        break;
        case DegreesOfFreedom:
            if (m_nodes) {
                numSamples = m_nodes->getNumDegreesOfFreedom();
                numDataPointsPerSample = 1;
            }
        break;
        default:
            stringstream ss;
            ss << "Invalid function space type: " << functionSpaceCode
                << " for domain " << getDescription();
            throw ValueError(ss.str());
    }
    return pair<int,dim_t>(numDataPointsPerSample, numSamples);
}

//
// adds linear PDE of second order into a given stiffness matrix and
// right hand side
//
void DudleyDomain::addPDEToSystem(
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

    Assemble_PDE(m_nodes, m_elements, mat.getPtr(), rhs, A, B, C, D, X, Y);
    Assemble_PDE(m_nodes, m_faceElements, mat.getPtr(), rhs,
                 escript::Data(), escript::Data(), escript::Data(), d,
                 escript::Data(), y);
    Assemble_PDE(m_nodes, m_points, mat.getPtr(), rhs, escript::Data(),
                 escript::Data(), escript::Data(), d_dirac,
                 escript::Data(), y_dirac);

#ifdef ESYS_HAVE_TRILINOS
    if (tm) {
        tm->fillComplete(true);
    }
#endif
}

void DudleyDomain::addPDEToLumpedSystem(escript::Data& mat,
                                        const escript::Data& D,
                                        const escript::Data& d,
                                        const escript::Data& d_dirac,
                                        bool useHRZ) const
{
    Assemble_LumpedSystem(m_nodes, m_elements, mat, D, useHRZ);
    Assemble_LumpedSystem(m_nodes, m_faceElements, mat, d, useHRZ);
    Assemble_LumpedSystem(m_nodes, m_points, mat, d_dirac, useHRZ);
}

//
// adds linear PDE of second order into the right hand side only
//
void DudleyDomain::addPDEToRHS(escript::Data& rhs, const escript::Data& X,
          const escript::Data& Y, const escript::Data& y,
          const escript::Data& y_contact, const escript::Data& y_dirac) const
{
    if (!y_contact.isEmpty())
        throw DudleyException("Dudley does not support y_contact");

    Assemble_PDE(m_nodes, m_elements, escript::ASM_ptr(), rhs,
                 escript::Data(), escript::Data(), escript::Data(),
                 escript::Data(), X, Y);

    Assemble_PDE(m_nodes, m_faceElements, escript::ASM_ptr(), rhs,
                 escript::Data(), escript::Data(), escript::Data(),
                 escript::Data(), escript::Data(), y);

    Assemble_PDE(m_nodes, m_points, escript::ASM_ptr(), rhs,
                 escript::Data(), escript::Data(), escript::Data(),
                 escript::Data(), escript::Data(), y_dirac);
}

//
// adds PDE of second order into a transport problem
//
void DudleyDomain::addPDEToTransportProblem(
        escript::AbstractTransportProblem& tp, escript::Data& source,
        const escript::Data& M, const escript::Data& A, const escript::Data& B,
        const escript::Data& C, const escript::Data& D, const escript::Data& X,
        const escript::Data& Y, const escript::Data& d, const escript::Data& y,
        const escript::Data& d_contact, const escript::Data& y_contact,
        const escript::Data& d_dirac, const escript::Data& y_dirac) const
{
    if (!d_contact.isEmpty())
        throw DudleyException("Dudley does not support d_contact");
    if (!y_contact.isEmpty())
        throw DudleyException("Dudley does not support y_contact");

#ifdef ESYS_HAVE_PASO
    paso::TransportProblem* ptp = dynamic_cast<paso::TransportProblem*>(&tp);
    if (!ptp)
        throw ValueError("Dudley only supports Paso transport problems.");

    source.expand();

    escript::ASM_ptr mm(boost::static_pointer_cast<escript::AbstractSystemMatrix>(
                ptp->borrowMassMatrix()));
    escript::ASM_ptr tm(boost::static_pointer_cast<escript::AbstractSystemMatrix>(
                ptp->borrowTransportMatrix()));

    Assemble_PDE(m_nodes, m_elements, mm, source, escript::Data(),
                 escript::Data(), escript::Data(), M, escript::Data(),
                 escript::Data());
    Assemble_PDE(m_nodes, m_elements, tm, source, A, B, C, D, X, Y);
    Assemble_PDE(m_nodes, m_faceElements, tm, source, escript::Data(),
                 escript::Data(), escript::Data(), d, escript::Data(), y);
    Assemble_PDE(m_nodes, m_points, tm, source, escript::Data(),
                 escript::Data(), escript::Data(), d_dirac, escript::Data(),
                 y_dirac);
#else
    throw DudleyException("Transport problems require the Paso library which "
                          "is not available.");
#endif
}

//
// interpolates data between different function spaces
//
void DudleyDomain::interpolateOnDomain(escript::Data& target,
                                      const escript::Data& in) const
{
    if (*in.getFunctionSpace().getDomain() != *this)
        throw ValueError("Illegal domain of interpolant.");
    if (*target.getFunctionSpace().getDomain() != *this)
        throw ValueError("Illegal domain of interpolation target.");

    switch (in.getFunctionSpace().getTypeCode()) {
        case Nodes:
            switch (target.getFunctionSpace().getTypeCode()) {
                case Nodes:
                case DegreesOfFreedom:
                    if (in.isComplex())
                        Assemble_CopyNodalData<cplx_t>(m_nodes, target, in);
                    else
                        Assemble_CopyNodalData<real_t>(m_nodes, target, in);
                break;
                case Elements:
                case ReducedElements:
                    if (in.isComplex())
                        Assemble_interpolate<cplx_t>(m_nodes, m_elements, in, target);
                    else
                        Assemble_interpolate<real_t>(m_nodes, m_elements, in, target);
                break;
                case FaceElements:
                case ReducedFaceElements:
                    if (in.isComplex())
                        Assemble_interpolate<cplx_t>(m_nodes, m_faceElements, in, target);
                    else
                        Assemble_interpolate<real_t>(m_nodes, m_faceElements, in, target);
                break;
                case Points:
                    if (in.isComplex())
                        Assemble_interpolate<cplx_t>(m_nodes, m_points, in, target);
                    else
                        Assemble_interpolate<real_t>(m_nodes, m_points, in, target);
                break;
                default:
                    stringstream ss;
                    ss << "interpolateOnDomain: Dudley does not know anything "
                          "about function space type "
                          << target.getFunctionSpace().getTypeCode();
                    throw ValueError(ss.str());
            }
        break;
        case Elements:
            if (target.getFunctionSpace().getTypeCode() == Elements) {
                if (in.isComplex())
                    Assemble_CopyElementData<cplx_t>(m_elements, target, in);
                else
                    Assemble_CopyElementData<real_t>(m_elements, target, in);
            } else if (target.getFunctionSpace().getTypeCode()==ReducedElements) {
                if (in.isComplex())
                    Assemble_AverageElementData<cplx_t>(m_elements, target, in);
                else
                    Assemble_AverageElementData<real_t>(m_elements, target, in);
            } else {
                throw ValueError("No interpolation with data on elements possible.");
            }
            break;
        case ReducedElements:
            if (target.getFunctionSpace().getTypeCode() == ReducedElements) {
                if (in.isComplex())
                    Assemble_CopyElementData<cplx_t>(m_elements, target, in);
                else
                    Assemble_CopyElementData<real_t>(m_elements, target, in);
            } else {
                throw ValueError("No interpolation with data on elements "
                                 "with reduced integration order possible.");
            }
            break;
        case FaceElements:
            if (target.getFunctionSpace().getTypeCode() == FaceElements) {
                if (in.isComplex())
                    Assemble_CopyElementData<cplx_t>(m_faceElements, target, in);
                else
                    Assemble_CopyElementData<real_t>(m_faceElements, target, in);
            } else if (target.getFunctionSpace().getTypeCode() == ReducedFaceElements) {
                if (in.isComplex())
                    Assemble_AverageElementData<cplx_t>(m_faceElements, target, in);
                else
                    Assemble_AverageElementData<real_t>(m_faceElements, target, in);
            } else {
                throw ValueError("No interpolation with data on face elements possible.");
            }
            break;
        case ReducedFaceElements:
            if (target.getFunctionSpace().getTypeCode() == ReducedFaceElements) {
                if (in.isComplex())
                    Assemble_CopyElementData<cplx_t>(m_faceElements, target, in);
                else
                    Assemble_CopyElementData<real_t>(m_faceElements, target, in);
            } else {
                throw ValueError("No interpolation with data on face "
                         "elements with reduced integration order possible.");
            }
            break;
        case Points:
            if (target.getFunctionSpace().getTypeCode() == Points) {
                if (in.isComplex())
                    Assemble_CopyElementData<cplx_t>(m_points, target, in);
                else
                    Assemble_CopyElementData<real_t>(m_points, target, in);
            } else {
                throw ValueError("No interpolation with data on points possible.");
            }
            break;
        case DegreesOfFreedom:
            switch (target.getFunctionSpace().getTypeCode()) {
                case DegreesOfFreedom:
                    if (in.isComplex())
                        Assemble_CopyNodalData<cplx_t>(m_nodes, target, in);
                    else
                        Assemble_CopyNodalData<real_t>(m_nodes, target, in);
                break;

                case Nodes:
                    if (getMPISize() > 1) {
                        escript::Data temp(in);
                        temp.expand();
                        if (in.isComplex())
                            Assemble_CopyNodalData<cplx_t>(m_nodes, target, temp);
                        else
                            Assemble_CopyNodalData<real_t>(m_nodes, target, temp);
                    } else {
                        if (in.isComplex())
                            Assemble_CopyNodalData<cplx_t>(m_nodes, target, in);
                        else
                            Assemble_CopyNodalData<real_t>(m_nodes, target, in);
                    }
                break;
                case Elements:
                case ReducedElements:
                    if (getMPISize() > 1) {
                        escript::Data temp(in, continuousFunction(*this));
                        if (in.isComplex())
                            Assemble_interpolate<cplx_t>(m_nodes, m_elements, temp, target);
                        else
                            Assemble_interpolate<real_t>(m_nodes, m_elements, temp, target);
                    } else {
                        if (in.isComplex())
                            Assemble_interpolate<cplx_t>(m_nodes, m_elements, in, target);
                        else
                            Assemble_interpolate<real_t>(m_nodes, m_elements, in, target);
                    }
                break;
                case FaceElements:
                case ReducedFaceElements:
                    if (getMPISize() > 1) {
                        escript::Data temp(in, continuousFunction(*this));
                        if (in.isComplex())
                            Assemble_interpolate<cplx_t>(m_nodes, m_faceElements, temp, target);
                        else
                            Assemble_interpolate<real_t>(m_nodes, m_faceElements, temp, target);
                    } else {
                        if (in.isComplex())
                            Assemble_interpolate<cplx_t>(m_nodes, m_faceElements, in, target);
                        else
                            Assemble_interpolate<real_t>(m_nodes, m_faceElements, in, target);
                    }
                break;
                case Points:
                    if (getMPISize() > 1) {
                        //escript::Data temp(in, continuousFunction(*this));
                    } else {
                        if (in.isComplex())
                            Assemble_interpolate<cplx_t>(m_nodes, m_points, in, target);
                        else
                            Assemble_interpolate<real_t>(m_nodes, m_points, in, target);
                    }
                break;
                default:
                    stringstream ss;
                    ss << "interpolateOnDomain: Dudley does not know anything "
                          "about function space type "
                       << target.getFunctionSpace().getTypeCode();
                    throw ValueError(ss.str());
            }
            break;
        default:
            stringstream ss;
            ss << "interpolateOnDomain: Dudley does not know anything about "
                "function space type " << in.getFunctionSpace().getTypeCode();
            throw ValueError(ss.str());
    }
}

//
// copies the locations of sample points into x
//
void DudleyDomain::setToX(escript::Data& arg) const
{
    if (*arg.getFunctionSpace().getDomain() != *this)
        throw ValueError("setToX: Illegal domain of data point locations");

    // in case of appropriate function space we can do the job directly:
    if (arg.getFunctionSpace().getTypeCode() == Nodes) {
        Assemble_NodeCoordinates(m_nodes, arg);
    } else {
        escript::Data tmp_data = Vector(0., continuousFunction(*this), true);
        Assemble_NodeCoordinates(m_nodes, tmp_data);
        // this is then interpolated onto arg:
        interpolateOnDomain(arg, tmp_data);
    }
}

//
// return the normal vectors at the location of data points as a Data object
//
void DudleyDomain::setToNormal(escript::Data& normal) const
{
    if (*normal.getFunctionSpace().getDomain() != *this)
        throw ValueError("setToNormal: Illegal domain of normal locations");

    if (normal.getFunctionSpace().getTypeCode() == FaceElements ||
            normal.getFunctionSpace().getTypeCode() == ReducedFaceElements) {
        Assemble_getNormal(m_nodes, m_faceElements, normal);
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
void DudleyDomain::interpolateAcross(escript::Data& /*target*/,
                                    const escript::Data& /*source*/) const
{
    throw escript::NotImplementedError("Dudley does not allow interpolation "
                                       "across domains.");
}

//
// calculates the integral of a function defined on arg
//
void DudleyDomain::setToIntegrals(vector<real_t>& integrals,
                                  const escript::Data& arg) const
{
    setToIntegralsWorker<real_t>(integrals, arg);
}

void DudleyDomain::setToIntegrals(vector<cplx_t>& integrals,
                                  const escript::Data& arg) const
{
    setToIntegralsWorker<cplx_t>(integrals, arg);
}

template<typename Scalar>
void DudleyDomain::setToIntegralsWorker(vector<Scalar>& integrals,
                                        const escript::Data& arg) const
{
    if (*arg.getFunctionSpace().getDomain() != *this)
        throw ValueError("setToIntegrals: Illegal domain of integration kernel");

    switch (arg.getFunctionSpace().getTypeCode()) {
        case Nodes: // fall through
        case DegreesOfFreedom:
        {
            escript::Data temp(arg, escript::function(*this));
            Assemble_integrate(m_nodes, m_elements, temp, integrals);
        }
        break;
        case Elements: // fall through
        case ReducedElements:
            Assemble_integrate(m_nodes, m_elements, arg, integrals);
        break;
        case FaceElements: // fall through
        case ReducedFaceElements:
            Assemble_integrate(m_nodes, m_faceElements, arg, integrals);
        break;
        case Points:
            throw ValueError("Integral of data on points is not supported.");
        default:
            stringstream ss;
            ss << "setToIntegrals: Dudley does not know anything about "
                "function space type " << arg.getFunctionSpace().getTypeCode();
            throw ValueError(ss.str());
    }
}

//
// calculates the gradient of arg
//
void DudleyDomain::setToGradient(escript::Data& grad, const escript::Data& arg) const
{
    if (*arg.getFunctionSpace().getDomain() != *this)
        throw ValueError("setToGradient: Illegal domain of gradient argument");
    if (*grad.getFunctionSpace().getDomain() != *this)
        throw ValueError("setToGradient: Illegal domain of gradient");
    if (grad.isComplex() != arg.isComplex())
        throw ValueError("setToGradient: Complexity of input and output must match");

    escript::Data nodeData;
    if (getMPISize() > 1 && arg.getFunctionSpace().getTypeCode() == DegreesOfFreedom) {
        nodeData = escript::Data(arg, continuousFunction(*this));
    } else {
        nodeData = arg;
    }
    switch (grad.getFunctionSpace().getTypeCode()) {
        case Nodes:
            throw ValueError("Gradient at nodes is not supported.");
        case Elements:
        case ReducedElements:
            if (arg.isComplex())
                Assemble_gradient<cplx_t>(m_nodes, m_elements, grad, nodeData);
            else
                Assemble_gradient<real_t>(m_nodes, m_elements, grad, nodeData);
            break;
        case FaceElements:
        case ReducedFaceElements:
            if (arg.isComplex())
                Assemble_gradient<cplx_t>(m_nodes, m_faceElements, grad, nodeData);
            else
                Assemble_gradient<real_t>(m_nodes, m_faceElements, grad, nodeData);
            break;
        case Points:
            throw ValueError("Gradient at points is not supported.");
        case DegreesOfFreedom:
            throw ValueError("Gradient at degrees of freedom is not supported.");
        default:
            stringstream ss;
            ss << "Gradient: Dudley does not know anything about function "
                  "space type " << arg.getFunctionSpace().getTypeCode();
            throw ValueError(ss.str());
    }
}

//
// returns the size of elements
//
void DudleyDomain::setToSize(escript::Data& size) const
{
    switch (size.getFunctionSpace().getTypeCode()) {
        case Nodes:
            throw ValueError("Size of nodes is not supported.");
        case Elements:
        case ReducedElements:
            Assemble_getSize(m_nodes, m_elements, size);
            break;
        case FaceElements:
        case ReducedFaceElements:
            Assemble_getSize(m_nodes, m_faceElements, size);
            break;
        case Points:
            throw ValueError("Size of point elements is not supported.");
        case DegreesOfFreedom:
            throw ValueError("Size of degrees of freedom is not supported.");
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
void DudleyDomain::setNewX(const escript::Data& newX)
{
    if (*newX.getFunctionSpace().getDomain() != *this)
        throw DudleyException("Illegal domain of new point locations");

    if (newX.getFunctionSpace() == continuousFunction(*this)) {
        m_nodes->setCoordinates(newX);
    } else {
        throw ValueError("As of escript version 3.3 setNewX only accepts "
                         "ContinuousFunction arguments. Please interpolate.");
    }
}

bool DudleyDomain::ownSample(int fs_code, index_t id) const
{
#ifdef ESYS_MPI
    if (getMPISize() > 1) {
        if (fs_code == Nodes) {
            const index_t myFirstNode = m_nodes->getFirstNode();
            const index_t myLastNode = m_nodes->getLastNode();
            const index_t k = m_nodes->borrowGlobalNodesIndex()[id];
            return (myFirstNode <= k && k < myLastNode);
        } else {
            std::ostringstream oss;
            oss << "ownSample: unsupported function space type (" << fs_code <<")";
            throw ValueError(oss.str());            
        }
    }
#endif
    return true;
}

//
// creates a stiffness matrix and initializes it with zeros
//
escript::ASM_ptr DudleyDomain::newSystemMatrix(int row_blocksize,
                            const escript::FunctionSpace& row_functionspace,
                            int column_blocksize,
                            const escript::FunctionSpace& column_functionspace,
                            int type) const
{
    // is the domain right?
    if (*row_functionspace.getDomain() != *this)
        throw ValueError("domain of row function space does not match the domain of matrix generator.");
    if (*column_functionspace.getDomain() != *this)
        throw ValueError("domain of column function space does not match the domain of matrix generator.");

    // is the function space type right?
    if (row_functionspace.getTypeCode() != DegreesOfFreedom) {
        throw ValueError("illegal function space type for system matrix rows.");
    }
    if (column_functionspace.getTypeCode() != DegreesOfFreedom) {
        throw ValueError("illegal function space type for system matrix columns.");
    }

    // generate matrix
    if (type & (int)SMT_TRILINOS) {
#ifdef ESYS_HAVE_TRILINOS
        const_TrilinosGraph_ptr graph(getTrilinosGraph());
        bool isComplex = (type & (int)SMT_COMPLEX);
        bool unroll = (type & (int)SMT_UNROLL);
        escript::ASM_ptr sm(new TrilinosMatrixAdapter(m_mpiInfo, row_blocksize,
                    row_functionspace, graph, isComplex, unroll));
        return sm;
#else
        throw DudleyException("newSystemMatrix: dudley was not compiled "
                "with Trilinos support so the Trilinos solver stack cannot be "
                "used.");
#endif
    } else if (type & (int)SMT_PASO) {
#ifdef ESYS_HAVE_PASO
        paso::SystemMatrixPattern_ptr pattern(getPasoPattern());
        paso::SystemMatrix_ptr sm(new paso::SystemMatrix(type, pattern,
                  row_blocksize, column_blocksize, false, row_functionspace,
                  column_functionspace));
        return sm;
#else
        throw DudleyException("newSystemMatrix: dudley was not compiled "
                "with Paso support so the Paso solver stack cannot be used.");
#endif
    } else {
        throw DudleyException("newSystemMatrix: unknown matrix type ID");
    }
}

//
// creates a TransportProblem
//
escript::ATP_ptr DudleyDomain::newTransportProblem(int blocksize,
                                             const escript::FunctionSpace& fs,
                                             int type) const
{
    // is the domain right?
    if (*fs.getDomain() != *this)
        throw ValueError("domain of function space does not match the domain of transport problem generator.");
    // is the function space type right
    if (fs.getTypeCode() != DegreesOfFreedom) {
        throw ValueError("illegal function space type for transport problem.");
    }

#ifdef ESYS_HAVE_PASO
    // generate transport problem
    paso::SystemMatrixPattern_ptr pattern(getPasoPattern());
    paso::TransportProblem_ptr transportProblem(new paso::TransportProblem(
                                              pattern, blocksize, fs));
    return transportProblem;
#else
    throw DudleyException("Transport problems require the Paso library which "
                          "is not available.");
#endif
}

//
// returns true if data on functionSpaceCode is considered as being cell centered
//
bool DudleyDomain::isCellOriented(int functionSpaceCode) const
{
    switch (functionSpaceCode) {
        case Nodes:
        case DegreesOfFreedom:
            return false;
        case Elements:
        case FaceElements:
        case Points:
        case ReducedElements:
        case ReducedFaceElements:
            return true;
    }
    stringstream ss;
    ss << "isCellOriented: Dudley does not know anything about "
          "function space type " << functionSpaceCode;
    throw ValueError(ss.str());
}

bool
DudleyDomain::commonFunctionSpace(const vector<int>& fs, int& resultcode) const
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
            if (hasclass[5] == 1)
                resultcode=ReducedElements;
            else
                resultcode=Elements;
        } else if (hasline[2] == 1) {
            if (hasclass[7] == 1)
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

bool DudleyDomain::probeInterpolationOnDomain(int functionSpaceType_source,
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
                    throw ValueError(ss.str());
            }
            break;
    }
    stringstream ss;
    ss << "Interpolation On Domain: Dudley does not know anything "
          "about function space type " << functionSpaceType_source;
    throw ValueError(ss.str());
}

signed char DudleyDomain::preferredInterpolationOnDomain(
        int functionSpaceType_source, int functionSpaceType_target) const
{
    if (probeInterpolationOnDomain(functionSpaceType_source, functionSpaceType_target))
        return 1;
    if (probeInterpolationOnDomain(functionSpaceType_target, functionSpaceType_source))
        return -1;

    return 0;
}

bool DudleyDomain::probeInterpolationAcross(int /*source*/,
        const AbstractDomain& /*targetDomain*/, int /*target*/) const
{
    return false;
}

bool DudleyDomain::operator==(const AbstractDomain& other) const
{
    const DudleyDomain* temp = dynamic_cast<const DudleyDomain*>(&other);
    if (temp) {
        return (m_nodes == temp->m_nodes &&
                m_elements == temp->m_elements &&
                m_faceElements == temp->m_faceElements &&
                m_points == temp->m_points);
    }
    return false;
}

bool DudleyDomain::operator!=(const AbstractDomain& other) const
{
    return !(operator==(other));
}

int DudleyDomain::getSystemMatrixTypeId(const bp::object& options) const
{
    const escript::SolverBuddy& sb = bp::extract<escript::SolverBuddy>(options);

    int package = sb.getPackage();
    escript::SolverOptions method = sb.getSolverMethod();
#ifdef ESYS_HAVE_TRILINOS
    bool isDirect = escript::isDirectSolver(method);
#endif

    // the configuration of dudley should have taken care that we have either
    // paso or trilinos so here's how we prioritize
#if defined(ESYS_HAVE_PASO) && defined(ESYS_HAVE_TRILINOS)
    // we have Paso & Trilinos so use Trilinos for parallel direct solvers and
    // for complex problems
    if (package == escript::SO_DEFAULT) {
        if ((method == escript::SO_METHOD_DIRECT && getMPISize() > 1)
                || isDirect
                || sb.isComplex()) {
            package = escript::SO_PACKAGE_TRILINOS;
        }
    }
#endif
#ifdef ESYS_HAVE_PASO
    if (package == escript::SO_DEFAULT)
        package = escript::SO_PACKAGE_PASO;
#endif
#ifdef ESYS_HAVE_TRILINOS
    if (package == escript::SO_DEFAULT)
        package = escript::SO_PACKAGE_TRILINOS;
#endif
    if (package == escript::SO_PACKAGE_TRILINOS) {
#ifdef ESYS_HAVE_TRILINOS
        int type = (int)SMT_TRILINOS;
        if (sb.isComplex())
            type |= (int)SMT_COMPLEX;
        // This is required because MueLu (AMG) and Amesos2 (direct) do not
        // support block matrices at this point. Remove if they ever do...
        if (sb.getPreconditioner() == escript::SO_PRECONDITIONER_AMG ||
                sb.getPreconditioner() == escript::SO_PRECONDITIONER_ILUT ||
                isDirect) {
            type |= (int)SMT_UNROLL;
        }
        return type;
#else
        throw DudleyException("Trilinos requested but not built with Trilinos.");
#endif
    }
#ifdef ESYS_HAVE_PASO
    if (sb.isComplex()) {
        throw NotImplementedError("Paso does not support complex-valued matrices");
    }
    return (int)SMT_PASO | paso::SystemMatrix::getSystemMatrixTypeId(
                method, sb.getPreconditioner(), sb.getPackage(),
                sb.isSymmetric(), m_mpiInfo);
#else
    throw DudleyException("Unable to find a working solver library!");
#endif
}

int DudleyDomain::getTransportTypeId(int solver, int preconditioner,
                                    int package, bool symmetry) const
{
#ifdef ESYS_HAVE_PASO
    return paso::TransportProblem::getTypeId(solver, preconditioner, package,
                                             symmetry, getMPI());
#else
    throw DudleyException("Transport solvers require Paso but dudley was not "
                          "compiled with Paso!");
#endif
}

escript::Data DudleyDomain::getX() const
{
    return continuousFunction(*this).getX();
}

escript::Data DudleyDomain::getNormal() const
{
    return functionOnBoundary(*this).getNormal();
}

escript::Data DudleyDomain::getSize() const
{
    return escript::function(*this).getSize();
}

const index_t* DudleyDomain::borrowSampleReferenceIDs(int functionSpaceType) const
{
    index_t* out = NULL;
    switch (functionSpaceType) {
        case Nodes:
            out = m_nodes->Id;
            break;
        case Elements:
            out = m_elements->Id;
            break;
        case ReducedElements:
            out = m_elements->Id;
            break;
        case FaceElements:
            out = m_faceElements->Id;
            break;
        case ReducedFaceElements:
            out = m_faceElements->Id;
            break;
        case Points:
            out = m_points->Id;
            break;
        case DegreesOfFreedom:
            out = m_nodes->degreesOfFreedomId;
            break;
        default:
            stringstream ss;
            ss << "Invalid function space type: " << functionSpaceType
               << " for domain: " << getDescription();
            throw ValueError(ss.str());
    }
    return out;
}

int DudleyDomain::getTagFromSampleNo(int functionSpaceType, index_t sampleNo) const
{
    int out = 0;
    switch (functionSpaceType) {
        case Nodes:
            out = m_nodes->Tag[sampleNo];
            break;
        case Elements:
        case ReducedElements:
            out = m_elements->Tag[sampleNo];
            break;
        case FaceElements:
        case ReducedFaceElements:
            out = m_faceElements->Tag[sampleNo];
            break;
        case Points:
            out = m_points->Tag[sampleNo];
            break;
        case DegreesOfFreedom:
            throw ValueError("DegreesOfFreedom does not support tags.");
        default:
            stringstream ss;
            ss << "Invalid function space type: " << functionSpaceType
               << " for domain: " << getDescription();
            throw ValueError(ss.str());
    }
    return out;
}


void DudleyDomain::setTags(int functionSpaceType, int newTag,
                           const escript::Data& mask) const
{
    switch (functionSpaceType) {
        case Nodes:
            m_nodes->setTags(newTag, mask);
            break;
        case DegreesOfFreedom:
            throw ValueError("DegreesOfFreedom does not support tags");
        case Elements: // fall through
        case ReducedElements:
            m_elements->setTags(newTag, mask);
            break;
        case FaceElements:
        case ReducedFaceElements:
            m_faceElements->setTags(newTag, mask);
            break;
        case Points:
            m_points->setTags(newTag, mask);
            break;
        default:
            stringstream ss;
            ss << "Dudley does not know anything about function space type "
               << functionSpaceType;
            throw ValueError(ss.str());
    }
}

void DudleyDomain::setTagMap(const string& name, int tag)
{
    m_tagMap[name] = tag;
}

int DudleyDomain::getTag(const string& name) const
{
    TagMap::const_iterator it = m_tagMap.find(name);
    if (it == m_tagMap.end()) {
        stringstream ss;
        ss << "getTag: unknown tag name " << name << ".";
        throw escript::ValueError(ss.str());
    }
    return it->second;
}

bool DudleyDomain::isValidTagName(const string& name) const
{
    return (m_tagMap.count(name) > 0);
}

string DudleyDomain::showTagNames() const
{
    stringstream ss;
    TagMap::const_iterator it = m_tagMap.begin();
    while (it != m_tagMap.end()) {
        ss << it->first;
        ++it;
        if (it != m_tagMap.end())
            ss << ", ";
    }
    return ss.str();
}

int DudleyDomain::getNumberOfTagsInUse(int functionSpaceCode) const
{
    switch (functionSpaceCode) {
        case Nodes:
            return m_nodes->tagsInUse.size();
        case DegreesOfFreedom:
            throw ValueError("DegreesOfFreedom does not support tags");
        case Elements: // fall through
        case ReducedElements:
            return m_elements->tagsInUse.size();
        case FaceElements: // fall through
        case ReducedFaceElements:
            return m_faceElements->tagsInUse.size();
        case Points:
            return m_points->tagsInUse.size();
    }
    stringstream ss;
    ss << "Dudley does not know anything about function space type "
       << functionSpaceCode;
    throw ValueError(ss.str());
}

const int* DudleyDomain::borrowListOfTagsInUse(int functionSpaceCode) const
{
    switch (functionSpaceCode) {
        case Nodes:
            if (m_nodes->tagsInUse.empty())
                return NULL;
            else
                return &m_nodes->tagsInUse[0];
        case DegreesOfFreedom:
            throw DudleyException("DegreesOfFreedom does not support tags");
        case Elements: // fall through
        case ReducedElements:
            if (m_elements->tagsInUse.empty())
                return NULL;
            else
                return &m_elements->tagsInUse[0];
        case FaceElements: // fall through
        case ReducedFaceElements:
            if (m_faceElements->tagsInUse.empty())
                return NULL;
            else
                return &m_faceElements->tagsInUse[0];
        case Points:
            if (m_points->tagsInUse.empty())
                return NULL;
            else
                return &m_points->tagsInUse[0];
    }
    stringstream ss;
    ss << "Dudley does not know anything about function space type "
       << functionSpaceCode;
    throw ValueError(ss.str());
}


bool DudleyDomain::canTag(int functionSpaceCode) const
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

DudleyDomain::StatusType DudleyDomain::getStatus() const
{
    return m_nodes->status;
}

int DudleyDomain::getApproximationOrder(int functionSpaceCode) const
{
    switch (functionSpaceCode) {
        case Nodes:
        case DegreesOfFreedom:
            return 1;
        case Elements:
        case FaceElements:
        case Points:
            return 2;
        case ReducedElements:
        case ReducedFaceElements:
            return 0;
    }
    stringstream ss;
    ss << "Dudley does not know anything about function space type "
       << functionSpaceCode;
    throw ValueError(ss.str());
}

escript::Data DudleyDomain::randomFill(
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

/// prepares the mesh for further use
void DudleyDomain::prepare(bool optimize)
{
    // first step is to distribute the elements according to a global
    // distribution of DOF
    IndexVector distribution(m_mpiInfo->size + 1);

    // first we create dense labeling for the DOFs
    dim_t newGlobalNumDOFs = m_nodes->createDenseDOFLabeling();

    // create a distribution of the global DOFs and determine the MPI rank
    // controlling the DOFs on this processor
    m_mpiInfo->setDistribution(0, newGlobalNumDOFs - 1, &distribution[0]);

    // now the mesh is re-distributed according to the distribution vector
    // this will redistribute the Nodes and Elements including overlap and
    // will create an element colouring but will not create any mappings
    // (see later in this function)
    distributeByRankOfDOF(distribution);

    // at this stage we are able to start an optimization of the DOF
    // distribution using ParaMetis. On return distribution is altered and
    // new DOF IDs have been assigned
    if (optimize && m_mpiInfo->size > 1) {
        optimizeDOFDistribution(distribution);
        distributeByRankOfDOF(distribution);
    }
    // the local labelling of the degrees of freedom is optimized
    if (optimize) {
        optimizeDOFLabeling(distribution);
    }

    // rearrange elements with the aim of bringing elements closer to memory
    // locations of the nodes (distributed shared memory!):
    optimizeElementOrdering();

    // create the global indices
    IndexVector nodeDistribution(m_mpiInfo->size + 1);

    m_nodes->createDenseNodeLabeling(nodeDistribution, distribution);
    // create the missing mappings
    createMappings(distribution, nodeDistribution);

    updateTagList();
}

/// tries to reduce the number of colours for all element files
void DudleyDomain::createColoring(const index_t* node_localDOF_map)
{
    m_elements->createColoring(m_nodes->getNumNodes(), node_localDOF_map);
    m_faceElements->createColoring(m_nodes->getNumNodes(), node_localDOF_map);
    m_points->createColoring(m_nodes->getNumNodes(), node_localDOF_map);
}

/// redistributes elements to minimize communication during assemblage
void DudleyDomain::optimizeElementOrdering()
{
    m_elements->optimizeOrdering();
    m_faceElements->optimizeOrdering();
    m_points->optimizeOrdering();
}

/// regenerates list of tags in use for node file and element files
void DudleyDomain::updateTagList()
{
    m_nodes->updateTagList();
    m_elements->updateTagList();
    m_faceElements->updateTagList();
    m_points->updateTagList();
}

} // end of namespace

