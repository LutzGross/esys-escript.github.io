
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

#include "FinleyDomain.h"
#include "Assemble.h"
#include "FinleyException.h"
#include "IndexList.h"
#include "ReferenceElements.h"

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

#ifdef ESYS_HAVE_HDF5
  #include <H5Cpp.h>
#endif


using namespace std;
namespace bp = boost::python;
using escript::NotImplementedError;
using escript::ValueError;

namespace finley {

using escript::DataTypes::real_t;
using escript::DataTypes::cplx_t;

// define the static constants
FinleyDomain::FunctionSpaceNamesMapType FinleyDomain::m_functionSpaceTypeNames;

FinleyDomain::FinleyDomain(const string& name, int numDim, escript::JMPI jmpi) :
    escript::AbstractContinuousDomain(jmpi),
    m_name(name),
    approximationOrder(-1),
    reducedApproximationOrder(-1),
    integrationOrder(-1),
    reducedIntegrationOrder(-1),
    m_elements(NULL),
    m_faceElements(NULL),
    m_contactElements(NULL),
    m_points(NULL)
{
    // allocate node table
    m_nodes = new NodeFile(numDim, m_mpiInfo);
    setFunctionSpaceTypeNames();
}

FinleyDomain::FinleyDomain(const FinleyDomain& in) :
    escript::AbstractContinuousDomain(in.m_mpiInfo),
    m_name(in.m_name),
    approximationOrder(in.approximationOrder),
    reducedApproximationOrder(in.reducedApproximationOrder),
    integrationOrder(in.integrationOrder),
    reducedIntegrationOrder(in.reducedIntegrationOrder),
    m_nodes(in.m_nodes),
    m_elements(in.m_elements),
    m_faceElements(in.m_faceElements),
    m_contactElements(in.m_contactElements),
    m_points(in.m_points)
{
    setFunctionSpaceTypeNames();
}

FinleyDomain::~FinleyDomain()
{
    delete m_nodes;
    delete m_elements;
    delete m_faceElements;
    delete m_contactElements;
    delete m_points;
}

void FinleyDomain::setElements(ElementFile* elements)
{
    delete m_elements;
    m_elements = elements;
}

void FinleyDomain::setFaceElements(ElementFile* elements)
{
    delete m_faceElements;
    m_faceElements = elements;
}

void FinleyDomain::setContactElements(ElementFile* elements)
{
    delete m_contactElements;
    m_contactElements = elements;
}

void FinleyDomain::setPoints(ElementFile* elements)
{
    delete m_points;
    m_points = elements;
}

void FinleyDomain::setOrders()
{
    const int ORDER_MAX = 9999999;
    int locals[4] = { ORDER_MAX, ORDER_MAX, ORDER_MAX, ORDER_MAX };

    if (m_elements != NULL && m_elements->numElements > 0) {
        locals[0] = std::min(locals[0], m_elements->referenceElementSet->referenceElement->BasisFunctions->Type->numOrder);
        locals[1] = std::min(locals[1], m_elements->referenceElementSet->referenceElement->LinearBasisFunctions->Type->numOrder);
        locals[2] = std::min(locals[2], m_elements->referenceElementSet->referenceElement->integrationOrder);
        locals[3] = std::min(locals[3], m_elements->referenceElementSet->referenceElementReducedQuadrature->integrationOrder);
    }
    if (m_faceElements != NULL && m_faceElements->numElements > 0) {
        locals[0] = std::min(locals[0], m_faceElements->referenceElementSet->referenceElement->BasisFunctions->Type->numOrder);
        locals[1] = std::min(locals[1], m_faceElements->referenceElementSet->referenceElement->LinearBasisFunctions->Type->numOrder);
        locals[2] = std::min(locals[2], m_faceElements->referenceElementSet->referenceElement->integrationOrder);
        locals[3] = std::min(locals[3], m_faceElements->referenceElementSet->referenceElementReducedQuadrature->integrationOrder);
    }
    if (m_contactElements != NULL && m_contactElements->numElements > 0) {
        locals[0] = std::min(locals[0], m_contactElements->referenceElementSet->referenceElement->BasisFunctions->Type->numOrder);
        locals[1] = std::min(locals[1], m_contactElements->referenceElementSet->referenceElement->LinearBasisFunctions->Type->numOrder);
        locals[2] = std::min(locals[2], m_contactElements->referenceElementSet->referenceElement->integrationOrder);
        locals[3] = std::min(locals[3], m_contactElements->referenceElementSet->referenceElementReducedQuadrature->integrationOrder);
    }

#ifdef ESYS_MPI
    int globals[4];
    MPI_Allreduce(locals, globals, 4, MPI_INT, MPI_MIN, m_mpiInfo->comm);
    approximationOrder = (globals[0] < ORDER_MAX ? globals[0] : -1);
    reducedApproximationOrder = (globals[1] < ORDER_MAX ? globals[1] : -1);
    integrationOrder = (globals[2] < ORDER_MAX ? globals[2] : -1);
    reducedIntegrationOrder = (globals[3] < ORDER_MAX ? globals[3] : -1);
#else
    approximationOrder = (locals[0] < ORDER_MAX ? locals[0] : -1);
    reducedApproximationOrder = (locals[1] < ORDER_MAX ? locals[1] : -1);
    integrationOrder = (locals[2] < ORDER_MAX ? locals[2] : -1);
    reducedIntegrationOrder = (locals[3] < ORDER_MAX ? locals[3] : -1);
#endif
}

void FinleyDomain::createMappings(const IndexVector& dofDist,
                                  const IndexVector& nodeDist)
{
    std::vector<short> maskReducedNodes(m_nodes->getNumNodes(), -1);
    markNodes(maskReducedNodes, 0, true);
    IndexVector indexReducedNodes = util::packMask(maskReducedNodes);
    m_nodes->createNodeMappings(indexReducedNodes, dofDist, nodeDist);
}

void FinleyDomain::markNodes(vector<short>& mask, index_t offset,
                             bool useLinear) const
{
    m_elements->markNodes(mask, offset, useLinear);
    m_faceElements->markNodes(mask, offset, useLinear);
    m_contactElements->markNodes(mask, offset, useLinear);
    m_points->markNodes(mask, offset, useLinear);
}

void FinleyDomain::relabelElementNodes(const IndexVector& newNode, index_t offset)
{
    m_elements->relabelNodes(newNode, offset);
    m_faceElements->relabelNodes(newNode, offset);
    m_contactElements->relabelNodes(newNode, offset);
    m_points->relabelNodes(newNode, offset);
}

void FinleyDomain::dump(const string& fileName) const
{
#ifdef ESYS_HAVE_HDF5
    int mpi_size = getMPISize();
    int mpi_rank = getMPIRank();
    #ifdef ESYS_INDEXTYPE_LONG
        H5::DataType h5_type_index = H5::PredType::NATIVE_LONG;
    #else
        H5::DataType h5_type_index = H5::PredType::NATIVE_INT;
    #endif
    H5::DataType h5_type_double = H5::PredType::NATIVE_DOUBLE;
    #ifdef ESYS_MPI
    MPI_Status status;
    int dummy = 0;
    #endif
    #ifdef ESYS_MPI
        if (mpi_rank > 0)
            MPI_Recv(&dummy, 1, MPI_INT, mpi_rank-1, 81800, getMPIComm(), &status);
    #endif
    try
    {
        const string newFileName(m_mpiInfo->appendRankToFileName(fileName));
        // Incoming token indicates it's my turn to write

        // Open HDF5 file :
        H5::H5File h5_file(newFileName, H5F_ACC_TRUNC);
        // Save same meta-data
        H5::Group h5_grp_meta = h5_file.createGroup("Meta");

        hsize_t h5_dims_mpi[1] = { 2 };
        long long h5_values_mpi[2] = { mpi_size,  mpi_rank };
        H5::Attribute h5_attr_mpi = h5_grp_meta.createAttribute("mpi", H5::PredType::NATIVE_LLONG, H5::DataSpace(1, h5_dims_mpi ) );
        h5_attr_mpi.write( H5::PredType::NATIVE_LLONG, h5_values_mpi );

        hsize_t h5_dims_order[1] = { 1 };
        H5::Attribute h5_attr_approximationOrder = h5_grp_meta.createAttribute("approximationOrder", H5::PredType::NATIVE_INT, H5::DataSpace(1, h5_dims_order ) );
        h5_attr_approximationOrder.write( H5::PredType::NATIVE_INT, &approximationOrder);
        H5::Attribute h5_attr_reducedApproximationOrder = h5_grp_meta.createAttribute("reducedApproximationOrder", H5::PredType::NATIVE_INT, H5::DataSpace(1, h5_dims_order ) );
        h5_attr_reducedApproximationOrder.write( H5::PredType::NATIVE_INT, &reducedApproximationOrder );
        H5::Attribute h5_attr_integrationOrder = h5_grp_meta.createAttribute("integrationOrder", H5::PredType::NATIVE_INT, H5::DataSpace(1, h5_dims_order ) );
        h5_attr_integrationOrder.write( H5::PredType::NATIVE_INT, &integrationOrder );
        H5::Attribute h5_attr_reducedIntegrationOrder = h5_grp_meta.createAttribute("reducedIntegrationOrder", H5::PredType::NATIVE_INT, H5::DataSpace(1, h5_dims_order ) );
        h5_attr_reducedIntegrationOrder.write( H5::PredType::NATIVE_INT, &reducedIntegrationOrder );

        hsize_t h5_dims_name[1] = {1};
        H5::StrType h5_type_name(0, H5T_VARIABLE);
        H5::Attribute h5_attr_name = h5_grp_meta.createAttribute("name", h5_type_name, H5::DataSpace(1, h5_dims_name));
        h5_attr_name.write(h5_type_name, &m_name);


        // Node File
        dim_t numNodes = m_nodes->getNumNodes();
        int numDim = m_nodes->numDim;
        H5::Group h5_grp_nodes = h5_file.createGroup("Nodes");

        hsize_t h5_dims_numDim[1] = { 1 };
        uint h5_values_numDim[1] = { static_cast<uint>(numDim) };
        H5::Attribute h5_attr_numDim = h5_grp_nodes.createAttribute("numDim", H5::PredType::NATIVE_UINT, H5::DataSpace(1, h5_dims_numDim ) );
        h5_attr_numDim.write( H5::PredType::NATIVE_UINT, h5_values_numDim );

        hsize_t h5_dims_numNodes[1] = { 1 };
        long h5_values_numNodes[1] = { static_cast<long>(numNodes) };
        H5::Attribute h5_attr_numNodes = h5_grp_nodes.createAttribute("numNodes", H5::PredType::NATIVE_LONG, H5::DataSpace(1, h5_dims_numNodes ) );
        h5_attr_numNodes.write( H5::PredType::NATIVE_LONG, h5_values_numNodes );

        hsize_t h5_dims_nodes[1] = { static_cast<hsize_t>(numNodes * numDim) };
        H5::DataSet h5_ds_nodes = h5_grp_nodes.createDataSet("Coordinates",  H5::DataType(h5_type_double), H5::DataSpace(1, h5_dims_nodes ) );
        h5_ds_nodes.write( &m_nodes->Coordinates[0],  H5::DataType(h5_type_double) );

        hsize_t h5_dims_node_distr[1] = { static_cast<hsize_t>(mpi_size + 1) };
        H5::DataSet h5_ds_node_distr = h5_grp_nodes.createDataSet("NodeDistribution",  H5::DataType(h5_type_index), H5::DataSpace(1, h5_dims_node_distr ) );
        h5_ds_node_distr.write( &m_nodes->nodesDistribution->first_component[0],  H5::DataType(h5_type_index));
        H5::DataSet h5_ds_dof_distr = h5_grp_nodes.createDataSet("DoFDistribution", H5::DataType(h5_type_index), H5::DataSpace(1, h5_dims_node_distr ) );
        h5_ds_dof_distr.write( &m_nodes->degreesOfFreedomDistribution->first_component[0], H5::DataType(h5_type_index));

        hsize_t h5_dims_ids[1] = { static_cast<hsize_t>(numNodes)  };
        H5::DataSet h5_ds_node_ids = h5_grp_nodes.createDataSet("Ids", H5::DataType(h5_type_index), H5::DataSpace(1, h5_dims_ids ) );
        h5_ds_node_ids.write( &m_nodes->Id[0], H5::DataType(h5_type_index));

        H5::DataSet h5_ds_node_tags = h5_grp_nodes.createDataSet("Tags", H5::PredType::NATIVE_INT, H5::DataSpace(1, h5_dims_ids ) );
        h5_ds_node_tags.write(&m_nodes->Tag[0], H5::PredType::NATIVE_INT);

        H5::DataSet h5_ds_node_gdof = h5_grp_nodes.createDataSet("globalDegreesOfFreedom", H5::DataType(h5_type_index), H5::DataSpace(1, h5_dims_ids ) );
        h5_ds_node_gdof.write( &m_nodes->globalDegreesOfFreedom[0], H5::DataType(h5_type_index));

        H5::DataSet h5_ds_node_gnidx = h5_grp_nodes.createDataSet("globalNodesIndex", H5::DataType(h5_type_index), H5::DataSpace(1, h5_dims_ids ) );
        h5_ds_node_gnidx.write( &m_nodes->globalNodesIndex[0], H5::DataType(h5_type_index));

        H5::DataSet h5_ds_node_rgdof = h5_grp_nodes.createDataSet("globalReducedDOFIndex", H5::DataType(h5_type_index), H5::DataSpace(1, h5_dims_ids ) );
        h5_ds_node_rgdof.write( &m_nodes->globalReducedDOFIndex[0], H5::DataType(h5_type_index));

        H5::DataSet h5_ds_node_rgnidx = h5_grp_nodes.createDataSet("globalReducedNodesIndex", H5::DataType(h5_type_index), H5::DataSpace(1, h5_dims_ids ) );
        h5_ds_node_rgnidx.write( &m_nodes->globalReducedNodesIndex[0], H5::DataType(h5_type_index));

        // dump elements
        m_elements->dump(h5_file.createGroup("Elements"));
        m_faceElements->dump(h5_file.createGroup("FaceElements"));
        m_contactElements->dump(h5_file.createGroup("ContactElements"));
        m_points->dump(h5_file.createGroup("Points"));
        // save tags:
        uint numTags = static_cast<uint>(m_tagMap.size());
        H5::Group h5_grp_tags = h5_file.createGroup("Tags");
        hsize_t h5_dims_numTags[1] = { 1 };
        uint h5_values_numTags[1] = { numTags };
        H5::Attribute h5_attr_numTags = h5_grp_tags.createAttribute("numTags", H5::PredType::NATIVE_UINT, H5::DataSpace(1, h5_dims_numTags ) );
        h5_attr_numTags.write( H5::PredType::NATIVE_UINT, h5_values_numTags );

        vector<int> Tag_keys;
        vector<const char*> Tag_names;
        TagMap::const_iterator it;
        for (it = m_tagMap.begin(); it != m_tagMap.end(); it++) {
                Tag_keys.push_back(it->second);
                Tag_names.push_back(it->first.c_str());
        }

        hsize_t h5_dims_tags[1] = { static_cast<hsize_t>(numTags) };
        H5::DataSet h5_ds_tags_ids = h5_grp_tags.createDataSet("TagIds", H5::PredType::NATIVE_INT, H5::DataSpace(1, h5_dims_tags ) );
        h5_ds_tags_ids.write( &Tag_keys[0], H5::PredType::NATIVE_INT);

        H5::StrType h5_type_tagnames(0, H5T_VARIABLE);
        H5::DataSet h5_ds_tagnames = h5_grp_tags.createDataSet("TagNames", h5_type_tagnames, H5::DataSpace(1, h5_dims_tags));
        h5_ds_tagnames.write( &Tag_names[0], h5_type_tagnames);

        h5_file.close();
    }
    // catch failure caused by the H5File operations
    catch (H5::Exception& error)
    {
        #ifdef ESYS_MPI
        if (mpi_rank < mpi_size-1)
            MPI_Send(&dummy, 1, MPI_INT, mpi_rank+1, 81800, getMPIComm());
        #endif
        error.printErrorStack();
        throw FinleyException("Error - DataConstant:: creating HDF5 file failed.");
    }
    // pass the MPI token:
    #ifdef ESYS_MPI
        if (mpi_rank < mpi_size-1)
            MPI_Send(&dummy, 1, MPI_INT, mpi_rank+1, 81800, getMPIComm());
    #endif
#else
    throw FinleyException("FinleyDomain::dump: not configured with NDF5. Please contact your installation manager.");
#endif
}

string FinleyDomain::getDescription() const
{
    return "FinleyMesh";
}

string FinleyDomain::functionSpaceTypeAsString(int functionSpaceType) const
{
    FunctionSpaceNamesMapType::iterator loc;
    loc = m_functionSpaceTypeNames.find(functionSpaceType);
    if (loc == m_functionSpaceTypeNames.end()) {
        return "Invalid function space type code.";
    } else {
        return loc->second;
    }
}

bool FinleyDomain::isValidFunctionSpaceType(int functionSpaceType) const
{
    FunctionSpaceNamesMapType::iterator loc;
    loc = m_functionSpaceTypeNames.find(functionSpaceType);
    return (loc != m_functionSpaceTypeNames.end());
}

void FinleyDomain::setFunctionSpaceTypeNames()
{
    m_functionSpaceTypeNames.insert(FunctionSpaceNamesMapType::value_type(
                DegreesOfFreedom,"Finley_DegreesOfFreedom [Solution(domain)]"));
    m_functionSpaceTypeNames.insert(FunctionSpaceNamesMapType::value_type(
                ReducedDegreesOfFreedom,"Finley_ReducedDegreesOfFreedom [ReducedSolution(domain)]"));
    m_functionSpaceTypeNames.insert(FunctionSpaceNamesMapType::value_type(
                Nodes,"Finley_Nodes [ContinuousFunction(domain)]"));
    m_functionSpaceTypeNames.insert(FunctionSpaceNamesMapType::value_type(
                ReducedNodes,"Finley_Reduced_Nodes [ReducedContinuousFunction(domain)]"));
    m_functionSpaceTypeNames.insert(FunctionSpaceNamesMapType::value_type(
                Elements,"Finley_Elements [Function(domain)]"));
    m_functionSpaceTypeNames.insert(FunctionSpaceNamesMapType::value_type(
                ReducedElements,"Finley_Reduced_Elements [ReducedFunction(domain)]"));
    m_functionSpaceTypeNames.insert(FunctionSpaceNamesMapType::value_type(
                FaceElements,"Finley_Face_Elements [FunctionOnBoundary(domain)]"));
    m_functionSpaceTypeNames.insert(FunctionSpaceNamesMapType::value_type(
                ReducedFaceElements,"Finley_Reduced_Face_Elements [ReducedFunctionOnBoundary(domain)]"));
    m_functionSpaceTypeNames.insert(FunctionSpaceNamesMapType::value_type(
                Points,"Finley_Points [DiracDeltaFunctions(domain)]"));
    m_functionSpaceTypeNames.insert(FunctionSpaceNamesMapType::value_type(
                ContactElementsZero,"Finley_Contact_Elements_0 [FunctionOnContactZero(domain)]"));
    m_functionSpaceTypeNames.insert(FunctionSpaceNamesMapType::value_type(
                ReducedContactElementsZero,"Finley_Reduced_Contact_Elements_0 [ReducedFunctionOnContactZero(domain)]"));
    m_functionSpaceTypeNames.insert(FunctionSpaceNamesMapType::value_type(
                ContactElementsOne,"Finley_Contact_Elements_1 [FunctionOnContactOne(domain)]"));
    m_functionSpaceTypeNames.insert(FunctionSpaceNamesMapType::value_type(
                ReducedContactElementsOne,"Finley_Reduced_Contact_Elements_1 [ReducedFunctionOnContactOne(domain)]"));
}

int FinleyDomain::getContinuousFunctionCode() const
{
    return Nodes;
}

int FinleyDomain::getReducedContinuousFunctionCode() const
{
    return ReducedNodes;
}

int FinleyDomain::getFunctionCode() const
{
    return Elements;
}

int FinleyDomain::getReducedFunctionCode() const
{
    return ReducedElements;
}

int FinleyDomain::getFunctionOnBoundaryCode() const
{
    return FaceElements;
}

int FinleyDomain::getReducedFunctionOnBoundaryCode() const
{
    return ReducedFaceElements;
}

int FinleyDomain::getFunctionOnContactZeroCode() const
{
    return ContactElementsZero;
}

int FinleyDomain::getReducedFunctionOnContactZeroCode() const
{
    return ReducedContactElementsZero;
}

int FinleyDomain::getFunctionOnContactOneCode() const
{
    return ContactElementsOne;
}

int FinleyDomain::getReducedFunctionOnContactOneCode() const
{
    return ReducedContactElementsOne;
}

int FinleyDomain::getSolutionCode() const
{
    return DegreesOfFreedom;
}

int FinleyDomain::getReducedSolutionCode() const
{
    return ReducedDegreesOfFreedom;
}

int FinleyDomain::getDiracDeltaFunctionsCode() const
{
    return Points;
}

//
// Return the number of data points summed across all MPI processes
//
dim_t FinleyDomain::getNumDataPointsGlobal() const
{
    return m_nodes->getGlobalNumNodes();
}

//
// return the number of data points per sample and the number of samples
// needed to represent data on a parts of the mesh.
//
pair<int,dim_t> FinleyDomain::getDataShape(int functionSpaceCode) const
{
    int numDataPointsPerSample = 0;
    dim_t numSamples = 0;
    switch (functionSpaceCode) {
        case Nodes:
            numDataPointsPerSample = 1;
            numSamples = m_nodes->getNumNodes();
        break;
        case ReducedNodes:
            numDataPointsPerSample = 1;
            numSamples = m_nodes->getNumReducedNodes();
        break;
        case Elements:
            if (m_elements) {
                numSamples = m_elements->numElements;
                numDataPointsPerSample = m_elements->referenceElementSet->referenceElement->Parametrization->numQuadNodes;
            }
        break;
        case ReducedElements:
            if (m_elements) {
                numSamples = m_elements->numElements;
                numDataPointsPerSample = m_elements->referenceElementSet->referenceElementReducedQuadrature->Parametrization->numQuadNodes;
            }
        break;
        case FaceElements:
            if (m_faceElements) {
                numSamples = m_faceElements->numElements;
                numDataPointsPerSample = m_faceElements->referenceElementSet->referenceElement->Parametrization->numQuadNodes;
            }
        break;
        case ReducedFaceElements:
            if (m_faceElements) {
                numSamples = m_faceElements->numElements;
                numDataPointsPerSample = m_faceElements->referenceElementSet->referenceElementReducedQuadrature->Parametrization->numQuadNodes;
            }
        break;
        case Points:
            if (m_points) {
                numSamples = m_points->numElements;
                numDataPointsPerSample = 1;
            }
        break;
        case ContactElementsZero:
        case ContactElementsOne:
            if (m_contactElements) {
                numSamples = m_contactElements->numElements;
                numDataPointsPerSample = m_contactElements->referenceElementSet->referenceElement->Parametrization->numQuadNodes;
            }
            break;
        case ReducedContactElementsZero:
        case ReducedContactElementsOne:
            if (m_contactElements) {
                numSamples = m_contactElements->numElements;
                numDataPointsPerSample = m_contactElements->referenceElementSet->referenceElementReducedQuadrature->Parametrization->numQuadNodes;
            }
            break;
        case DegreesOfFreedom:
            if (m_nodes) {
                numSamples = m_nodes->getNumDegreesOfFreedom();
                numDataPointsPerSample = 1;
            }
        break;
        case ReducedDegreesOfFreedom:
            if (m_nodes) {
                numSamples = m_nodes->getNumReducedDegreesOfFreedom();
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
void FinleyDomain::addPDEToSystem(
        escript::AbstractSystemMatrix& mat, escript::Data& rhs,
        const escript::Data& A, const escript::Data& B, const escript::Data& C,
        const escript::Data& D, const escript::Data& X, const escript::Data& Y,
        const escript::Data& d, const escript::Data& y,
        const escript::Data& d_contact, const escript::Data& y_contact,
        const escript::Data& d_dirac, const escript::Data& y_dirac) const
{
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
    Assemble_PDE(m_nodes, m_contactElements, mat.getPtr(), rhs,
                 escript::Data(), escript::Data(), escript::Data(), d_contact,
                 escript::Data(), y_contact);
    Assemble_PDE(m_nodes, m_points, mat.getPtr(), rhs, escript::Data(),
                 escript::Data(), escript::Data(), d_dirac,
                 escript::Data(), y_dirac);

#ifdef ESYS_HAVE_TRILINOS
    if (tm) {
        tm->fillComplete(true);
    }
#endif
}

void FinleyDomain::addPDEToLumpedSystem(escript::Data& mat,
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
void FinleyDomain::addPDEToRHS(escript::Data& rhs, const escript::Data& X,
          const escript::Data& Y, const escript::Data& y,
          const escript::Data& y_contact, const escript::Data& y_dirac) const
{
    Assemble_PDE(m_nodes, m_elements, escript::ASM_ptr(), rhs,
                 escript::Data(), escript::Data(), escript::Data(),
                 escript::Data(), X, Y);

    Assemble_PDE(m_nodes, m_faceElements, escript::ASM_ptr(), rhs,
                 escript::Data(), escript::Data(), escript::Data(),
                 escript::Data(), escript::Data(), y);

    Assemble_PDE(m_nodes, m_contactElements, escript::ASM_ptr(),
                 rhs, escript::Data(), escript::Data(), escript::Data(),
                 escript::Data(), escript::Data(), y_contact);

    Assemble_PDE(m_nodes, m_points, escript::ASM_ptr(), rhs,
                 escript::Data(), escript::Data(), escript::Data(),
                 escript::Data(), escript::Data(), y_dirac);
}

//
// adds PDE of second order into a transport problem
//
void FinleyDomain::addPDEToTransportProblem(
        escript::AbstractTransportProblem& tp, escript::Data& source,
        const escript::Data& M, const escript::Data& A, const escript::Data& B,
        const escript::Data& C, const escript::Data& D, const escript::Data& X,
        const escript::Data& Y, const escript::Data& d, const escript::Data& y,
        const escript::Data& d_contact, const escript::Data& y_contact,
        const escript::Data& d_dirac, const escript::Data& y_dirac) const
{
#ifdef ESYS_HAVE_PASO
    paso::TransportProblem* ptp = dynamic_cast<paso::TransportProblem*>(&tp);
    if (!ptp)
        throw ValueError("Finley only supports Paso transport problems.");

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
    Assemble_PDE(m_nodes, m_contactElements, tm, source,
                 escript::Data(), escript::Data(), escript::Data(), d_contact,
                 escript::Data(), y_contact);
    Assemble_PDE(m_nodes, m_points, tm, source, escript::Data(),
                 escript::Data(), escript::Data(), d_dirac, escript::Data(),
                 y_dirac);
#else
    throw FinleyException("Transport problems require the Paso library which "
                          "is not available.");
#endif
}

//
// interpolates data between different function spaces
//
void FinleyDomain::interpolateOnDomain(escript::Data& target,
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
                case ReducedNodes:
                case DegreesOfFreedom:
                case ReducedDegreesOfFreedom:
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
                case ContactElementsZero:
                case ReducedContactElementsZero:
                case ContactElementsOne:
                case ReducedContactElementsOne:
                    if (in.isComplex())
                        Assemble_interpolate<cplx_t>(m_nodes, m_contactElements, in, target);
                    else
                        Assemble_interpolate<real_t>(m_nodes, m_contactElements, in, target);
                break;
                default:
                    stringstream ss;
                    ss << "interpolateOnDomain: Finley does not know anything "
                          "about function space type "
                          << target.getFunctionSpace().getTypeCode();
                    throw ValueError(ss.str());
            }
        break;
        case ReducedNodes:
            switch(target.getFunctionSpace().getTypeCode()) {
                case Nodes:
                case ReducedNodes:
                case DegreesOfFreedom:
                case ReducedDegreesOfFreedom:
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
                case ContactElementsZero:
                case ReducedContactElementsZero:
                case ContactElementsOne:
                case ReducedContactElementsOne:
                    if (in.isComplex())
                        Assemble_interpolate<cplx_t>(m_nodes, m_contactElements, in, target);
                    else
                        Assemble_interpolate<real_t>(m_nodes, m_contactElements, in, target);
                break;
                default:
                    stringstream ss;
                    ss << "interpolateOnDomain: Finley does not know anything "
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
            } else if (target.getFunctionSpace().getTypeCode() == Elements) {
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
        case ContactElementsZero:
        case ContactElementsOne:
            if (target.getFunctionSpace().getTypeCode()==ContactElementsZero || target.getFunctionSpace().getTypeCode()==ContactElementsOne) {
                if (in.isComplex())
                    Assemble_CopyElementData<cplx_t>(m_contactElements, target, in);
                else
                    Assemble_CopyElementData<real_t>(m_contactElements, target, in);
            } else if (target.getFunctionSpace().getTypeCode()==ReducedContactElementsZero || target.getFunctionSpace().getTypeCode()==ReducedContactElementsOne) {
                if (in.isComplex())
                    Assemble_AverageElementData<cplx_t>(m_contactElements, target, in);
                else
                    Assemble_AverageElementData<real_t>(m_contactElements, target, in);
            } else {
                throw ValueError("No interpolation with data on contact elements possible.");
            }
            break;
        case ReducedContactElementsZero:
        case ReducedContactElementsOne:
            if (target.getFunctionSpace().getTypeCode()==ReducedContactElementsZero || target.getFunctionSpace().getTypeCode()==ReducedContactElementsOne) {
                if (in.isComplex())
                    Assemble_CopyElementData<cplx_t>(m_contactElements, target, in);
                else
                    Assemble_CopyElementData<real_t>(m_contactElements, target, in);
            } else {
                throw ValueError("No interpolation with data on contact elements with reduced integration order possible.");
            }
            break;
        case DegreesOfFreedom:
            switch (target.getFunctionSpace().getTypeCode()) {
                case ReducedDegreesOfFreedom:
                case DegreesOfFreedom:
                    if (in.isComplex())
                        Assemble_CopyNodalData<cplx_t>(m_nodes, target, in);
                    else
                        Assemble_CopyNodalData<real_t>(m_nodes, target, in);
                break;

                case Nodes:
                case ReducedNodes:
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
                        escript::Data temp(in, continuousFunction(*this));
                        if (in.isComplex())
                            Assemble_interpolate<cplx_t>(m_nodes, m_points, temp, target);
                        else
                            Assemble_interpolate<real_t>(m_nodes, m_points, temp, target);
                    } else {
                        if (in.isComplex())
                            Assemble_interpolate<cplx_t>(m_nodes, m_points, in, target);
                        else
                            Assemble_interpolate<real_t>(m_nodes, m_points, in, target);
                    }
                break;
                case ContactElementsZero:
                case ContactElementsOne:
                case ReducedContactElementsZero:
                case ReducedContactElementsOne:
                    if (getMPISize() > 1) {
                        escript::Data temp(in, continuousFunction(*this));
                        if (in.isComplex())
                            Assemble_interpolate<cplx_t>(m_nodes, m_contactElements, temp, target);
                        else
                            Assemble_interpolate<real_t>(m_nodes, m_contactElements, temp, target);
                    } else {
                        if (in.isComplex())
                            Assemble_interpolate<cplx_t>(m_nodes, m_contactElements, in, target);
                        else
                            Assemble_interpolate<real_t>(m_nodes, m_contactElements, in, target);
                    }
                    break;
                default:
                    stringstream ss;
                    ss << "interpolateOnDomain: Finley does not know anything "
                          "about function space type "
                       << target.getFunctionSpace().getTypeCode();
                    throw ValueError(ss.str());
            }
            break;
        case ReducedDegreesOfFreedom:
            switch (target.getFunctionSpace().getTypeCode()) {
                case Nodes:
                    throw ValueError("Finley does not support interpolation from reduced degrees of freedom to mesh nodes.");
                case ReducedNodes:
                    if (getMPISize() > 1) {
                        escript::Data in2(in);
                        in2.expand();
                        if (in.isComplex())
                            Assemble_CopyNodalData<cplx_t>(m_nodes, target, in2);
                        else
                            Assemble_CopyNodalData<real_t>(m_nodes, target, in2);
                    } else {
                        if (in.isComplex())
                            Assemble_CopyNodalData<cplx_t>(m_nodes, target, in);
                        else
                            Assemble_CopyNodalData<real_t>(m_nodes, target, in);
                    }
                    break;
                case DegreesOfFreedom:
                    throw ValueError("Finley does not support interpolation from reduced degrees of freedom to degrees of freedom");
                case ReducedDegreesOfFreedom:
                    if (in.isComplex())
                        Assemble_CopyNodalData<cplx_t>(m_nodes, target, in);
                    else
                        Assemble_CopyNodalData<real_t>(m_nodes, target, in);
                    break;
                case Elements:
                case ReducedElements:
                    if (getMPISize() > 1) {
                        escript::Data in2(in, reducedContinuousFunction(*this));
                        if (in.isComplex())
                            Assemble_interpolate<cplx_t>(m_nodes, m_elements, in2, target);
                        else
                            Assemble_interpolate<real_t>(m_nodes, m_elements, in2, target);
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
                        escript::Data in2(in, reducedContinuousFunction(*this));
                        if (in.isComplex())
                            Assemble_interpolate<cplx_t>(m_nodes, m_faceElements, in2, target);
                        else
                            Assemble_interpolate<real_t>(m_nodes, m_faceElements, in2, target);
                    } else {
                        if (in.isComplex())
                            Assemble_interpolate<cplx_t>(m_nodes, m_faceElements, in, target);
                        else
                            Assemble_interpolate<real_t>(m_nodes, m_faceElements, in, target);
                    }
                    break;
                case Points:
                    if (getMPISize() > 1) {
                        escript::Data in2(in, reducedContinuousFunction(*this));
                        if (in.isComplex())
                            Assemble_interpolate<cplx_t>(m_nodes, m_points, in2, target);
                        else
                            Assemble_interpolate<real_t>(m_nodes, m_points, in2, target);
                    } else {
                        if (in.isComplex())
                            Assemble_interpolate<cplx_t>(m_nodes, m_points, in, target);
                        else
                            Assemble_interpolate<real_t>(m_nodes, m_points, in, target);
                    }
                    break;
                case ContactElementsZero:
                case ContactElementsOne:
                case ReducedContactElementsZero:
                case ReducedContactElementsOne:
                    if (getMPISize()>1) {
                        escript::Data in2(in, reducedContinuousFunction(*this));
                        if (in.isComplex())
                            Assemble_interpolate<cplx_t>(m_nodes, m_contactElements, in2, target);
                        else
                            Assemble_interpolate<real_t>(m_nodes, m_contactElements, in2, target);
                    } else {
                        if (in.isComplex())
                            Assemble_interpolate<cplx_t>(m_nodes, m_contactElements, in, target);
                        else
                            Assemble_interpolate<real_t>(m_nodes, m_contactElements, in, target);
                    }
                    break;
                default:
                    stringstream ss;
                    ss << "interpolateOnDomain: Finley does not know anything about function space type " << target.getFunctionSpace().getTypeCode();
                    throw ValueError(ss.str());
            }
            break;
        default:
            stringstream ss;
            ss << "interpolateOnDomain: Finley does not know anything about "
                "function space type " << in.getFunctionSpace().getTypeCode();
            throw ValueError(ss.str());
    }
}

//
// copies the locations of sample points into x
//
void FinleyDomain::setToX(escript::Data& arg) const
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
void FinleyDomain::setToNormal(escript::Data& normal) const
{
    if (*normal.getFunctionSpace().getDomain() != *this)
        throw ValueError("setToNormal: Illegal domain of normal locations");

    if (normal.getFunctionSpace().getTypeCode() == FaceElements ||
            normal.getFunctionSpace().getTypeCode() == ReducedFaceElements) {
        Assemble_getNormal(m_nodes, m_faceElements, normal);
    } else if (normal.getFunctionSpace().getTypeCode() == ContactElementsOne ||
            normal.getFunctionSpace().getTypeCode() == ContactElementsZero ||
            normal.getFunctionSpace().getTypeCode() == ReducedContactElementsOne ||
            normal.getFunctionSpace().getTypeCode() == ReducedContactElementsZero) {
        Assemble_getNormal(m_nodes, m_contactElements, normal);
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
void FinleyDomain::interpolateAcross(escript::Data& /*target*/,
                                    const escript::Data& /*source*/) const
{
    throw NotImplementedError("Finley does not allow interpolation across "
                              "domains.");
}

//
// calculates the integral of a function defined on arg
//
void FinleyDomain::setToIntegrals(vector<real_t>& integrals,
                                  const escript::Data& arg) const
{
    setToIntegralsWorker<real_t>(integrals, arg);
}

void FinleyDomain::setToIntegrals(vector<cplx_t>& integrals,
                                  const escript::Data& arg) const
{
    setToIntegralsWorker<cplx_t>(integrals, arg);
}

template<typename Scalar>
void FinleyDomain::setToIntegralsWorker(vector<Scalar>& integrals,
                                        const escript::Data& arg) const
{
    if (*arg.getFunctionSpace().getDomain() != *this)
        throw ValueError("setToIntegrals: Illegal domain of integration kernel");

    switch (arg.getFunctionSpace().getTypeCode()) {
        case Nodes:
        case ReducedNodes:
        case DegreesOfFreedom:
        case ReducedDegreesOfFreedom:
        {
            escript::Data temp(arg, escript::function(*this));
            Assemble_integrate(m_nodes, m_elements, temp, &integrals[0]);
        }
        break;
        case Points:
            Assemble_integrate_points(m_points, arg, &integrals[0]);
        break;
        case Elements:
        case ReducedElements:
            Assemble_integrate(m_nodes, m_elements, arg, &integrals[0]);
        break;
        case FaceElements:
        case ReducedFaceElements:
            Assemble_integrate(m_nodes, m_faceElements, arg, &integrals[0]);
        break;
        case ContactElementsZero:
        case ReducedContactElementsZero:
        case ContactElementsOne:
        case ReducedContactElementsOne:
            Assemble_integrate(m_nodes, m_contactElements, arg, &integrals[0]);
        break;
        default:
            stringstream ss;
            ss << "setToIntegrals: Finley does not know anything about "
                "function space type " << arg.getFunctionSpace().getTypeCode();
            throw ValueError(ss.str());
    }
}

//
// calculates the gradient of arg
//
void FinleyDomain::setToGradient(escript::Data& grad, const escript::Data& arg) const
{
    if (*arg.getFunctionSpace().getDomain() != *this)
        throw ValueError("setToGradient: Illegal domain of gradient argument");
    if (*grad.getFunctionSpace().getDomain() != *this)
        throw ValueError("setToGradient: Illegal domain of gradient");
    if (grad.isComplex() != arg.isComplex())
        throw ValueError("setToGradient: Complexity of input and output must match");

    escript::Data nodeData;
    if (getMPISize() > 1) {
        if (arg.getFunctionSpace().getTypeCode() == DegreesOfFreedom) {
            nodeData = escript::Data(arg, continuousFunction(*this));
        } else if (arg.getFunctionSpace().getTypeCode() == ReducedDegreesOfFreedom) {
            nodeData = escript::Data(arg, reducedContinuousFunction(*this));
        } else {
            nodeData = arg;
        }
    } else {
        nodeData = arg;
    }
    switch (grad.getFunctionSpace().getTypeCode()) {
        case Nodes:
            throw ValueError("Gradient at nodes is not supported.");
        case ReducedNodes:
            throw ValueError("Gradient at reduced nodes is not supported.");
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
        case ContactElementsZero:
        case ReducedContactElementsZero:
        case ContactElementsOne:
        case ReducedContactElementsOne:
            if (arg.isComplex())
                Assemble_gradient<cplx_t>(m_nodes, m_contactElements, grad, nodeData);
            else
                Assemble_gradient<real_t>(m_nodes, m_contactElements, grad, nodeData);
        break;
        case Points:
            throw ValueError("Gradient at points is not supported.");
        case DegreesOfFreedom:
            throw ValueError("Gradient at degrees of freedom is not supported.");
        case ReducedDegreesOfFreedom:
            throw ValueError("Gradient at reduced degrees of freedom is not supported.");
        default:
            stringstream ss;
            ss << "Gradient: Finley does not know anything about function "
                  "space type " << arg.getFunctionSpace().getTypeCode();
            throw ValueError(ss.str());
    }
}

//
// returns the size of elements
//
void FinleyDomain::setToSize(escript::Data& size) const
{
    switch (size.getFunctionSpace().getTypeCode()) {
        case Nodes:
            throw ValueError("Size of nodes is not supported.");
        case ReducedNodes:
            throw ValueError("Size of reduced nodes is not supported.");
        case Elements:
        case ReducedElements:
            Assemble_getSize(m_nodes, m_elements, size);
            break;
        case FaceElements:
        case ReducedFaceElements:
            Assemble_getSize(m_nodes, m_faceElements, size);
            break;
        case ContactElementsZero:
        case ContactElementsOne:
        case ReducedContactElementsZero:
        case ReducedContactElementsOne:
            Assemble_getSize(m_nodes, m_contactElements, size);
            break;
        case Points:
            throw ValueError("Size of point elements is not supported.");
        case DegreesOfFreedom:
            throw ValueError("Size of degrees of freedom is not supported.");
        case ReducedDegreesOfFreedom:
            throw ValueError("Size of reduced degrees of freedom is not supported.");
        default:
            stringstream ss;
            ss << "setToSize: Finley does not know anything about function "
                  "space type " << size.getFunctionSpace().getTypeCode();
            throw ValueError(ss.str());
    }
}

//
// sets the location of nodes
//
void FinleyDomain::setNewX(const escript::Data& newX)
{
    if (*newX.getFunctionSpace().getDomain() != *this)
        throw ValueError("Illegal domain of new point locations");

    if (newX.getFunctionSpace() == continuousFunction(*this)) {
        m_nodes->setCoordinates(newX);
    } else {
        throw ValueError("As of escript version 3.3 setNewX only accepts "
                         "ContinuousFunction arguments. Please interpolate.");
    }
}

bool FinleyDomain::ownSample(int fs_code, index_t id) const
{

#ifdef ESYS_MPI
    if (getMPISize() > 1 && fs_code != FINLEY_DEGREES_OF_FREEDOM &&
            fs_code != FINLEY_REDUCED_DEGREES_OF_FREEDOM) {
        if (fs_code == Nodes || fs_code == Elements ||
            fs_code == ReducedNodes || fs_code == ReducedElements ||
            fs_code == FaceElements || fs_code == ReducedFaceElements) {
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
escript::ASM_ptr FinleyDomain::newSystemMatrix(int row_blocksize,
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

    bool reduceRowOrder = false;
    bool reduceColOrder = false;
    // is the function space type right?
    if (row_functionspace.getTypeCode() == ReducedDegreesOfFreedom) {
        reduceRowOrder = true;
    } else if (row_functionspace.getTypeCode() != DegreesOfFreedom) {
        throw ValueError("illegal function space type for system matrix rows.");
    }
    if (column_functionspace.getTypeCode() == ReducedDegreesOfFreedom) {
        reduceColOrder = true;
    } else if (column_functionspace.getTypeCode() != DegreesOfFreedom) {
        throw ValueError("illegal function space type for system matrix columns.");
    }

    // generate matrix
    if (type & (int)SMT_TRILINOS) {
#ifdef ESYS_HAVE_TRILINOS
        if (reduceRowOrder != reduceColOrder)
            throw ValueError("element order of matrix rows and columns must "
                             "match when using Trilinos");
        const_TrilinosGraph_ptr graph(getTrilinosGraph(reduceRowOrder));
        bool isComplex = (type & (int)SMT_COMPLEX);
        bool unroll = (type & (int)SMT_UNROLL);
        escript::ASM_ptr sm(new TrilinosMatrixAdapter(m_mpiInfo, row_blocksize,
                    row_functionspace, graph, isComplex, unroll));
        return sm;
#else
        throw FinleyException("newSystemMatrix: finley was not compiled "
                "with Trilinos support so the Trilinos solver stack cannot be "
                "used.");
#endif
    } else if (type & (int)SMT_PASO) {
#ifdef ESYS_HAVE_PASO
        paso::SystemMatrixPattern_ptr pattern(getPasoPattern(
                                            reduceRowOrder, reduceColOrder));
        if (type & (int)MATRIX_FORMAT_COMPLEX) {
            paso::SystemMatrix_ptr<cplx_t> sm(new paso::SystemMatrix<cplx_t>(type, pattern,
                      row_blocksize, column_blocksize, false, row_functionspace,
                      column_functionspace));
            return sm;
        } else {
            paso::SystemMatrix_ptr<double> sm(new paso::SystemMatrix<double>(type, pattern,
                      row_blocksize, column_blocksize, false, row_functionspace,
                      column_functionspace));
            return sm;
        }
#else
        throw FinleyException("newSystemMatrix: finley was not compiled "
                "with Paso support so the Paso solver stack cannot be used.");
#endif
    } else {
        throw FinleyException("newSystemMatrix: unknown matrix type ID");
    }
}

//
// creates a TransportProblem
//
escript::ATP_ptr FinleyDomain::newTransportProblem(int blocksize,
                                             const escript::FunctionSpace& fs,
                                             int type) const
{
    // is the domain right?
    if (*fs.getDomain() != *this)
        throw ValueError("domain of function space does not match the domain of transport problem generator.");

#ifdef ESYS_HAVE_PASO
    // is the function space type right
    bool reduceOrder = false;
    if (fs.getTypeCode() == ReducedDegreesOfFreedom) {
        reduceOrder = true;
    } else if (fs.getTypeCode() != DegreesOfFreedom) {
        throw ValueError("illegal function space type for transport problem.");
    }

    // generate transport problem
    paso::SystemMatrixPattern_ptr pattern(getPasoPattern(
                                                  reduceOrder, reduceOrder));
    paso::TransportProblem_ptr transportProblem(new paso::TransportProblem(
                                              pattern, blocksize, fs));
    return transportProblem;
#else
    throw FinleyException("Transport problems require the Paso library which "
                          "is not available.");
#endif
}

//
// returns true if data on functionSpaceCode is considered as being cell centered
//
bool FinleyDomain::isCellOriented(int functionSpaceCode) const
{
    switch (functionSpaceCode) {
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
    }
    stringstream ss;
    ss << "isCellOriented: Finley does not know anything about "
          "function space type " << functionSpaceCode;
    throw ValueError(ss.str());
}

bool
FinleyDomain::commonFunctionSpace(const vector<int>& fs, int& resultcode) const
{
    if (fs.empty())
        return false;
    // The idea is to use equivalence classes, i.e. types which can be
    // interpolated back and forth
    //    class 1: DOF <-> Nodes
    //    class 2: ReducedDOF <-> ReducedNodes
    //    class 3: Points
    //    class 4: Elements
    //    class 5: ReducedElements
    //    class 6: FaceElements
    //    class 7: ReducedFaceElements
    //    class 8: ContactElementZero <-> ContactElementOne
    //    class 9: ReducedContactElementZero <-> ReducedContactElementOne

    // There is also a set of lines. Interpolation is possible down a line but
    // not between lines.
    // class 1 and 2 belong to all lines so aren't considered.
    //    line 0: class 3
    //    line 1: class 4,5
    //    line 2: class 6,7
    //    line 3: class 8,9

    // For classes with multiple members (e.g. class 2) we have vars to record
    // if there is at least one instance.
    // e.g. hasnodes is true if we have at least one instance of Nodes.
    vector<int> hasclass(10);
    vector<int> hasline(4);
    bool hasnodes = false;
    bool hasrednodes = false;
    bool hascez = false;
    bool hasrcez = false;
    for (int i = 0; i < fs.size(); ++i) {
        switch (fs[i]) {
            case Nodes:
                hasnodes = true; // fall through
            case DegreesOfFreedom:
                hasclass[1] = 1;
                break;
            case ReducedNodes:
                hasrednodes = true; // fall through
            case ReducedDegreesOfFreedom:
                hasclass[2] = 1;
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
            case ContactElementsZero:
                hascez = true; // fall through
            case ContactElementsOne:
                hasclass[8] = 1;
                hasline[3] = 1;
                break;
            case ReducedContactElementsZero:
                hasrcez = true; // fall through
            case ReducedContactElementsOne:
                hasclass[9] = 1;
                hasline[3] = 1;
                break;
            default:
                return false;
        }
    }
    int totlines = hasline[0]+hasline[1]+hasline[2]+hasline[3];

    // fail if we have more than one leaf group
    if (totlines > 1)
        // there are at least two branches we can't interpolate between
        return false;
    else if (totlines == 1) {
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
        } else { // so we must be in line3
            if (hasclass[9] == 1) {
                // need something from class 9
                resultcode = (hasrcez ? ReducedContactElementsZero : ReducedContactElementsOne);
            } else {
                // something from class 8
                resultcode = (hascez ? ContactElementsZero : ContactElementsOne);
            }
        }
    } else { // totlines==0
        if (hasclass[2] == 1) {
            // something from class 2
            resultcode = (hasrednodes ? ReducedNodes : ReducedDegreesOfFreedom);
        } else {
            // something from class 1
            resultcode = (hasnodes ? Nodes : DegreesOfFreedom);
        }
    }
    return true;
}

bool FinleyDomain::probeInterpolationOnDomain(int functionSpaceType_source,
                                              int functionSpaceType_target) const
{
    switch(functionSpaceType_source) {
        case Nodes:
            switch (functionSpaceType_target) {
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
                    stringstream ss;
                    ss << "Interpolation On Domain: Finley does not know "
                          "anything about function space type "
                       << functionSpaceType_target;
                    throw ValueError(ss.str());
            }
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
                    stringstream ss;
                    ss << "Interpolation On Domain: Finley does not know "
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
        case ContactElementsZero:
        case ContactElementsOne:
            return (functionSpaceType_target == ContactElementsZero ||
                    functionSpaceType_target == ContactElementsOne ||
                    functionSpaceType_target == ReducedContactElementsZero ||
                    functionSpaceType_target == ReducedContactElementsOne);
        case ReducedContactElementsZero:
        case ReducedContactElementsOne:
            return (functionSpaceType_target == ReducedContactElementsZero ||
                    functionSpaceType_target == ReducedContactElementsOne);
        case DegreesOfFreedom:
            switch (functionSpaceType_target) {
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
                    stringstream ss;
                    ss << "Interpolation On Domain: Finley does not know "
                          "anything about function space type "
                       << functionSpaceType_target;
                    throw ValueError(ss.str());
            }
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
                    stringstream ss;
                    ss << "Interpolation On Domain: Finley does not know "
                          "anything about function space type "
                       << functionSpaceType_target;
                    throw ValueError(ss.str());
            }
    }
    stringstream ss;
    ss << "Interpolation On Domain: Finley does not know anything "
          "about function space type " << functionSpaceType_source;
    throw ValueError(ss.str());
}

signed char FinleyDomain::preferredInterpolationOnDomain(
        int functionSpaceType_source, int functionSpaceType_target) const
{
    if (probeInterpolationOnDomain(functionSpaceType_source, functionSpaceType_target))
        return 1;
    if (probeInterpolationOnDomain(functionSpaceType_target, functionSpaceType_source))
        return -1;

    return 0;
}

bool FinleyDomain::probeInterpolationAcross(int /*source*/,
        const AbstractDomain& /*targetDomain*/, int /*target*/) const
{
    return false;
}

bool FinleyDomain::operator==(const AbstractDomain& other) const
{
    const FinleyDomain* temp = dynamic_cast<const FinleyDomain*>(&other);
    if (temp) {
        return (m_nodes == temp->m_nodes &&
                m_elements == temp->m_elements &&
                m_faceElements == temp->m_faceElements &&
                m_contactElements == temp->m_contactElements &&
                m_points == temp->m_points);
    }
    return false;
}

bool FinleyDomain::operator!=(const AbstractDomain& other) const
{
    return !(operator==(other));
}

int FinleyDomain::getSystemMatrixTypeId(const bp::object& options) const
{
    const escript::SolverBuddy& sb = bp::extract<escript::SolverBuddy>(options);

    int package = sb.getPackage();
    escript::SolverOptions method = sb.getSolverMethod();
#ifdef ESYS_HAVE_TRILINOS
    bool isDirect = escript::isDirectSolver(method);
#endif

    // the configuration of finley should have taken care that we have either
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
        throw FinleyException("Trilinos requested but not built with Trilinos.");
#endif
    }
#ifdef ESYS_HAVE_PASO
    if (sb.isComplex()) {
#ifdef ESYS_HAVE_MUMPS
        return (int)SMT_PASO | paso::SystemMatrix<cplx_t>::getSystemMatrixTypeId(
                    method, sb.getPreconditioner(), sb.getPackage(),
                    sb.isComplex(), sb.isSymmetric(), m_mpiInfo);
#else
        throw NotImplementedError("Paso requires MUMPS for complex-valued matrices.");
#endif
    } else {
        return (int)SMT_PASO | paso::SystemMatrix<double>::getSystemMatrixTypeId(
                    method, sb.getPreconditioner(), sb.getPackage(),
                    sb.isComplex(), sb.isSymmetric(), m_mpiInfo);
    }
#else
    throw FinleyException("Unable to find a working solver library!");
#endif
}

int FinleyDomain::getTransportTypeId(int solver, int preconditioner,
                                    int package, bool symmetry) const
{
#ifdef ESYS_HAVE_PASO
    return paso::TransportProblem::getTypeId(solver, preconditioner, package,
                                             symmetry, getMPI());
#else
    throw FinleyException("Transport solvers require Paso but finley was not "
                          "compiled with Paso!");
#endif
}

escript::Data FinleyDomain::getX() const
{
    return continuousFunction(*this).getX();
}

#ifdef ESYS_HAVE_BOOST_NUMPY
boost::python::numpy::ndarray FinleyDomain::getNumpyX() const
{
    return continuousFunction(*this).getNumpyX();
}

boost::python::numpy::ndarray FinleyDomain::getConnectivityInfo() const
{
    // Initialise boost numpy
    boost::python::numpy::initialize();

    // Get the node information
    escript::DataTypes::index_t * nodedata = m_elements->Nodes;

    // Work out how many elements there are
    int numberOfElements = m_elements->numElements;

    // Work out how many data points there are per element
    int numNodesPerElement = m_elements->numNodes;

    // Initialise the ndarray
    boost::python::tuple arrayshape = boost::python::make_tuple(numberOfElements, numNodesPerElement);
    boost::python::numpy::dtype datatype = boost::python::numpy::dtype::get_builtin<double>();
    boost::python::numpy::ndarray dataArray = boost::python::numpy::zeros(arrayshape, datatype);

    // Initialise variables
    std::string localmsg;
    std::vector<const escript::DataTypes::real_t*> samplesR(1);

    // Copy the information over
// #pragma omp parallel for
    for (int k = 0; k < numberOfElements; k++) {
        for (int j = 0; j < numNodesPerElement; j++) {
            dataArray[k][j] = nodedata[j+k*numNodesPerElement];
        }
    }

    // Print out the ndarray to the console - used during debugging
    // std::cout << "Finished array:\n" << bp::extract<char const *>(bp::str(dataArray)) << std::endl;

    return dataArray;
}
#endif

int FinleyDomain::getVTKElementType() const
{
    const_ReferenceElementSet_ptr refElement = m_elements->referenceElementSet;
    const_ReferenceElement_ptr borrowedRefElement = refElement->borrowReferenceElement(false);
    const ReferenceElementInfo* type = borrowedRefElement->Type;

    // From vtkCellType.h
    // #define VTK_EMPTY_CELL 0
    // #define VTK_VERTEX 1
    // #define VTK_POLY_VERTEX 2
    // #define VTK_LINE 3
    // #define VTK_POLY_LINE 4
    // #define VTK_TRIANGLE 5
    // #define VTK_TRIANGLE_STRIP 6
    // #define VTK_POLYGON 7
    // #define VTK_PIXEL 8
    // #define VTK_QUAD 9
    // #define VTK_TETRA 10
    // #define VTK_VOXEL 11
    // #define VTK_HEXAHEDRON 12
    // #define VTK_WEDGE 13
    // #define VTK_PYRAMID 14

    if(std::strcmp(type->Name, "Tri3")
        || std::strcmp(type->Name, "Tri6")
        || std::strcmp(type->Name, "Tri9")
        || std::strcmp(type->Name, "Tri10"))
    {
        return 5; //VTK_TRIANGLE
    }
    else if (std::strcmp(type->Name, "Rec4")
        || std::strcmp(type->Name, "Rec8")
        || std::strcmp(type->Name, "Rec9")
        || std::strcmp(type->Name, "Rec12")
        || std::strcmp(type->Name, "Rec16"))
    {
        return 8; // VTK_PIXEL
    }
    else if (std::strcmp(type->Name, "Tet4")
        || std::strcmp(type->Name, "Tet10")
        || std::strcmp(type->Name, "Tet16"))
    {
        return 10; // VTK_TETRA
    }
    else if (std::strcmp(type->Name, "Hex8")
        || std::strcmp(type->Name, "Hex20")
        || std::strcmp(type->Name, "Hex27")
        || std::strcmp(type->Name, "Hex32"))
    {
        return 11; // VTK_VOXEL
    }
    else
    {
        throw FinleyException("getVTKElementType: unknown element type");
    }
}

escript::Data FinleyDomain::getNormal() const
{
    return functionOnBoundary(*this).getNormal();
}

escript::Data FinleyDomain::getSize() const
{
    return escript::function(*this).getSize();
}

const index_t* FinleyDomain::borrowSampleReferenceIDs(int functionSpaceType) const
{
    index_t* out = NULL;
    switch (functionSpaceType) {
        case Nodes:
            out = m_nodes->Id;
            break;
        case ReducedNodes:
            out = m_nodes->reducedNodesId;
            break;
        case Elements:
        case ReducedElements:
            out = m_elements->Id;
            break;
        case FaceElements:
        case ReducedFaceElements:
            out = m_faceElements->Id;
            break;
        case Points:
            out = m_points->Id;
            break;
        case ContactElementsZero:
        case ReducedContactElementsZero:
        case ContactElementsOne:
        case ReducedContactElementsOne:
            out = m_contactElements->Id;
            break;
        case DegreesOfFreedom:
            out = m_nodes->degreesOfFreedomId;
            break;
        case ReducedDegreesOfFreedom:
            out = m_nodes->reducedDegreesOfFreedomId;
            break;
        default:
            stringstream ss;
            ss << "Invalid function space type: " << functionSpaceType
               << " for domain: " << getDescription();
            throw ValueError(ss.str());
    }
    return out;
}
int FinleyDomain::getTagFromSampleNo(int functionSpaceType, index_t sampleNo) const
{
    int out = 0;
    switch (functionSpaceType) {
        case Nodes:
            out = m_nodes->Tag[sampleNo];
            break;
        case ReducedNodes:
            throw ValueError("ReducedNodes does not support tags.");
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
        case ContactElementsZero:
        case ReducedContactElementsZero:
        case ContactElementsOne:
        case ReducedContactElementsOne:
            out = m_contactElements->Tag[sampleNo];
            break;
        case DegreesOfFreedom:
            throw ValueError("DegreesOfFreedom does not support tags.");
        case ReducedDegreesOfFreedom:
            throw ValueError("ReducedDegreesOfFreedom does not support tags.");
        default:
            stringstream ss;
            ss << "Invalid function space type: " << functionSpaceType
               << " for domain: " << getDescription();
            throw ValueError(ss.str());
    }
    return out;
}


void FinleyDomain::setTags(int functionSpaceType, int newTag, const escript::Data& mask) const
{
    switch (functionSpaceType) {
        case Nodes:
            m_nodes->setTags(newTag, mask);
            break;
        case ReducedNodes:
            throw ValueError("ReducedNodes does not support tags");
        case DegreesOfFreedom:
            throw ValueError("DegreesOfFreedom does not support tags");
        case ReducedDegreesOfFreedom:
            throw ValueError("ReducedDegreesOfFreedom does not support tags");
        case Elements:
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
        case ContactElementsZero:
        case ReducedContactElementsZero:
        case ContactElementsOne:
        case ReducedContactElementsOne:
            m_contactElements->setTags(newTag, mask);
            break;
        default:
            stringstream ss;
            ss << "Finley does not know anything about function space type "
               << functionSpaceType;
            throw ValueError(ss.str());
    }
}

void FinleyDomain::setTagMap(const string& name, int tag)
{
    m_tagMap[name] = tag;
}

int FinleyDomain::getTag(const string& name) const
{
    TagMap::const_iterator it = m_tagMap.find(name);
    if (it == m_tagMap.end()) {
        stringstream ss;
        ss << "getTag: unknown tag name " << name << ".";
        throw escript::ValueError(ss.str());
    }
    return it->second;
}

bool FinleyDomain::isValidTagName(const string& name) const
{
    return (m_tagMap.count(name) > 0);
}

string FinleyDomain::showTagNames() const
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

int FinleyDomain::getNumberOfTagsInUse(int functionSpaceCode) const
{
    switch (functionSpaceCode) {
        case Nodes:
            return m_nodes->tagsInUse.size();
        case ReducedNodes:
            throw ValueError("ReducedNodes does not support tags");
        case DegreesOfFreedom:
            throw ValueError("DegreesOfFreedom does not support tags");
        case ReducedDegreesOfFreedom:
            throw ValueError("ReducedDegreesOfFreedom does not support tags");
        case Elements:
        case ReducedElements:
            return m_elements->tagsInUse.size();
        case FaceElements:
        case ReducedFaceElements:
            return m_faceElements->tagsInUse.size();
        case Points:
            return m_points->tagsInUse.size();
        case ContactElementsZero:
        case ReducedContactElementsZero:
        case ContactElementsOne:
        case ReducedContactElementsOne:
            return m_contactElements->tagsInUse.size();
    }
    stringstream ss;
    ss << "Finley does not know anything about function space type "
       << functionSpaceCode;
    throw ValueError(ss.str());
}

const int* FinleyDomain::borrowListOfTagsInUse(int functionSpaceCode) const
{
    switch (functionSpaceCode) {
        case Nodes:
            if (m_nodes->tagsInUse.empty())
                return NULL;
            else
                return &m_nodes->tagsInUse[0];
        case ReducedNodes:
            throw ValueError("ReducedNodes does not support tags");
        case DegreesOfFreedom:
            throw ValueError("DegreesOfFreedom does not support tags");
        case ReducedDegreesOfFreedom:
            throw ValueError("ReducedDegreesOfFreedom does not support tags");
        case Elements:
        case ReducedElements:
            if (m_elements->tagsInUse.empty())
                return NULL;
            else
                return &m_elements->tagsInUse[0];
        case FaceElements:
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
        case ContactElementsZero:
        case ReducedContactElementsZero:
        case ContactElementsOne:
        case ReducedContactElementsOne:
            if (m_contactElements->tagsInUse.empty())
                return NULL;
            else
                return &m_contactElements->tagsInUse[0];
    }
    stringstream ss;
    ss << "Finley does not know anything about function space type "
       << functionSpaceCode;
    throw ValueError(ss.str());
}

bool FinleyDomain::canTag(int functionSpaceCode) const
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

FinleyDomain::StatusType FinleyDomain::getStatus() const
{
    return m_nodes->status;
}

int FinleyDomain::getApproximationOrder(int functionSpaceCode) const
{
    switch(functionSpaceCode) {
        case Nodes:
        case DegreesOfFreedom:
            return approximationOrder;
        case ReducedNodes:
        case ReducedDegreesOfFreedom:
            return reducedApproximationOrder;
        case Elements:
        case FaceElements:
        case Points:
        case ContactElementsZero:
        case ContactElementsOne:
            return integrationOrder;
        case ReducedElements:
        case ReducedFaceElements:
        case ReducedContactElementsZero:
        case ReducedContactElementsOne:
            return reducedIntegrationOrder;
    }
    stringstream ss;
    ss << "Finley does not know anything about function space type "
       << functionSpaceCode;
    throw ValueError(ss.str());
}

escript::Data FinleyDomain::randomFill(
                                const escript::DataTypes::ShapeType& shape,
                                const escript::FunctionSpace& what, long seed,
                                const bp::tuple& filter) const
{
    escript::Data towipe(0, shape, what, true);
    // since we just made this object, no sharing is possible and we don't
    // need to check for exclusive write
    escript::DataTypes::RealVectorType& dv(towipe.getExpandedVectorReference());
    escript::randomFillArray(seed, &dv[0], dv.size(), m_mpiInfo);
    return towipe;
}

/// prepares the mesh for further use
void FinleyDomain::prepare(bool optimize)
{
    setOrders();

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
    std::vector<short> maskReducedNodes(m_nodes->getNumNodes(), -1);
    IndexVector nodeDistribution(m_mpiInfo->size + 1);
    markNodes(maskReducedNodes, 0, true);
    IndexVector indexReducedNodes = util::packMask(maskReducedNodes);

    m_nodes->createDenseNodeLabeling(nodeDistribution, distribution);
    // created reduced DOF labeling
    m_nodes->createDenseReducedLabeling(maskReducedNodes, false);
    // created reduced node labeling
    m_nodes->createDenseReducedLabeling(maskReducedNodes, true);
    // create the missing mappings
    m_nodes->createNodeMappings(indexReducedNodes, distribution, nodeDistribution);

    updateTagList();
}

/// redistributes the Nodes and Elements including overlap
/// according to the DOF distribution. It will create an element colouring
/// but will not create any mappings.
void FinleyDomain::distributeByRankOfDOF(const std::vector<index_t>& dof_distribution)
{
    std::vector<int> mpiRankOfDOF(m_nodes->getNumNodes());
    m_nodes->assignMPIRankToDOFs(mpiRankOfDOF, dof_distribution);

    // first, the elements are redistributed according to mpiRankOfDOF
    // at the input the Node tables refer to the local labeling of the nodes
    // while at the output they refer to the global labeling which is rectified
    // in the next step
    m_elements->distributeByRankOfDOF(mpiRankOfDOF, m_nodes->Id);
    m_faceElements->distributeByRankOfDOF(mpiRankOfDOF, m_nodes->Id);
    m_contactElements->distributeByRankOfDOF(mpiRankOfDOF, m_nodes->Id);
    m_points->distributeByRankOfDOF(mpiRankOfDOF, m_nodes->Id);

    resolveNodeIds();

    // create a local labeling of the DOFs
    const std::pair<index_t,index_t> dof_range(m_nodes->getDOFRange());
    const index_t len = dof_range.second-dof_range.first+1;
    // local mask for used nodes
    std::vector<index_t> localDOF_mask(len, -1);
    std::vector<index_t> localDOF_map(m_nodes->getNumNodes(), -1);

#pragma omp parallel for
    for (index_t n = 0; n < m_nodes->getNumNodes(); n++) {
#ifdef BOUNDS_CHECK
        ESYS_ASSERT(m_nodes->globalDegreesOfFreedom[n]-dof_range.first < len, "BOUNDS_CHECK");
        ESYS_ASSERT(m_nodes->globalDegreesOfFreedom[n]-dof_range.first >= 0, "BOUNDS_CHECK");
#endif
        localDOF_mask[m_nodes->globalDegreesOfFreedom[n]-dof_range.first] = n;
    }

    index_t numDOFs = 0;
    for (index_t n = 0; n < len; n++) {
        const index_t k = localDOF_mask[n];
        if (k >= 0) {
             localDOF_mask[n] = numDOFs;
             numDOFs++;
          }
    }
#pragma omp parallel for
    for (index_t n = 0; n < m_nodes->getNumNodes(); n++) {
        const index_t k = localDOF_mask[m_nodes->globalDegreesOfFreedom[n]-dof_range.first];
        localDOF_map[n] = k;
    }
    // create element coloring
    createColoring(localDOF_map);
}

/// optimizes the labeling of the DOFs on each processor
void FinleyDomain::optimizeDOFLabeling(const IndexVector& distribution)
{
    // this method relies on Pattern::reduceBandwidth so requires PASO
    // at the moment
#ifdef ESYS_HAVE_PASO
    const int myRank = getMPIRank();
    const int mpiSize = getMPISize();
    const index_t myFirstVertex = distribution[myRank];
    const index_t myLastVertex = distribution[myRank+1];
    const dim_t myNumVertices = myLastVertex-myFirstVertex;
    dim_t len = 0;
    for (int p = 0; p < mpiSize; ++p)
        len=std::max(len, distribution[p+1]-distribution[p]);

    boost::scoped_array<IndexList> index_list(new IndexList[myNumVertices]);
    boost::scoped_array<index_t> newGlobalDOFID(new index_t[len]);

    // create the adjacency structure xadj and adjncy
#pragma omp parallel
    {
        // insert contributions from element matrices into columns index
        IndexList_insertElementsWithRowRangeNoMainDiagonal(index_list.get(),
                myFirstVertex, myLastVertex, m_elements,
                m_nodes->globalDegreesOfFreedom,
                m_nodes->globalDegreesOfFreedom);
        IndexList_insertElementsWithRowRangeNoMainDiagonal(index_list.get(),
                myFirstVertex, myLastVertex, m_faceElements,
                m_nodes->globalDegreesOfFreedom,
                m_nodes->globalDegreesOfFreedom);
        IndexList_insertElementsWithRowRangeNoMainDiagonal(index_list.get(),
                myFirstVertex, myLastVertex, m_contactElements,
                m_nodes->globalDegreesOfFreedom,
                m_nodes->globalDegreesOfFreedom);
        IndexList_insertElementsWithRowRangeNoMainDiagonal(index_list.get(),
                myFirstVertex, myLastVertex, m_points,
                m_nodes->globalDegreesOfFreedom,
                m_nodes->globalDegreesOfFreedom);
    }
    // create the local matrix pattern
    paso::Pattern_ptr pattern = paso::Pattern::fromIndexListArray(0,
            myNumVertices, index_list.get(), myFirstVertex, myLastVertex,
            -myFirstVertex);

    pattern->reduceBandwidth(&newGlobalDOFID[0]);

    // shift new labeling to create a global id
#pragma omp parallel for
    for (index_t i = 0; i < myNumVertices; ++i)
        newGlobalDOFID[i] += myFirstVertex;

    // distribute new labeling to other processors
#ifdef ESYS_MPI
    const int dest = m_mpiInfo->mod_rank(myRank + 1);
    const int source = m_mpiInfo->mod_rank(myRank - 1);
#endif
    int current_rank = myRank;
    for (int p = 0; p < mpiSize; ++p) {
        const index_t firstVertex = distribution[current_rank];
        const index_t lastVertex = distribution[current_rank + 1];
#pragma omp parallel for
        for (index_t i = 0; i < m_nodes->getNumNodes(); ++i) {
            const index_t k = m_nodes->globalDegreesOfFreedom[i];
            if (firstVertex <= k && k < lastVertex) {
                m_nodes->globalDegreesOfFreedom[i] = newGlobalDOFID[k-firstVertex];
            }
        }

        if (p < mpiSize - 1) { // the final send can be skipped
#ifdef ESYS_MPI
            MPI_Status status;
            MPI_Sendrecv_replace(&newGlobalDOFID[0], len, MPI_DIM_T,
                                 dest, m_mpiInfo->counter(), source,
                                 m_mpiInfo->counter(), m_mpiInfo->comm, &status);
            m_mpiInfo->incCounter();
#endif
            current_rank = m_mpiInfo->mod_rank(current_rank - 1);
        }
    }
#endif // ESYS_HAVE_PASO
}

void FinleyDomain::resolveNodeIds()
{
    // find the minimum and maximum id used by elements
    index_t min_id = escript::DataTypes::index_t_max();
    index_t max_id = -escript::DataTypes::index_t_max();
    std::pair<index_t,index_t> range(m_elements->getNodeRange());
    max_id = std::max(max_id, range.second);
    min_id = std::min(min_id, range.first);
    range = m_faceElements->getNodeRange();
    max_id = std::max(max_id, range.second);
    min_id = std::min(min_id, range.first);
    range = m_contactElements->getNodeRange();
    max_id = std::max(max_id, range.second);
    min_id = std::min(min_id, range.first);
    range = m_points->getNodeRange();
    max_id = std::max(max_id, range.second);
    min_id = std::min(min_id, range.first);
#ifdef Finley_TRACE
    index_t global_min_id, global_max_id;
#ifdef ESYS_MPI
    index_t id_range[2], global_id_range[2];
    id_range[0] = -min_id;
    id_range[1] = max_id;
    MPI_Allreduce(id_range, global_id_range, 2, MPI_DIM_T, MPI_MAX, m_mpiInfo->comm);
    global_min_id = -global_id_range[0];
    global_max_id = global_id_range[1];
#else
    global_min_id = min_id;
    global_max_id = max_id;
#endif
    printf("Node id range used by elements is %d:%d\n", global_min_id, global_max_id);
#endif
    if (min_id > max_id) {
        max_id = -1;
        min_id = 0;
    }

    // allocate mappings for new local node labeling to global node labeling
    // (newLocalToGlobalNodeLabels) and global node labeling to the new local
    // node labeling (globalToNewLocalNodeLabels[i-min_id] is the new local id
    // of global node i)
    index_t len = (max_id >= min_id) ? max_id - min_id + 1 : 0;

    // mark the nodes referred by elements in usedMask
    std::vector<short> usedMask(len, -1);
    markNodes(usedMask, min_id, false);

    // create a local labeling newLocalToGlobalNodeLabels of the local nodes
    // by packing the mask usedMask
    std::vector<index_t> newLocalToGlobalNodeLabels =  util::packMask(usedMask);
    const dim_t newNumNodes = newLocalToGlobalNodeLabels.size();

    usedMask.clear();

    // invert the new labeling and shift the index newLocalToGlobalNodeLabels
    // to global node IDs
    std::vector<index_t> globalToNewLocalNodeLabels(len, -1);

#pragma omp parallel for
    for (index_t n = 0; n < newNumNodes; n++) {
#ifdef BOUNDS_CHECK
        ESYS_ASSERT(newLocalToGlobalNodeLabels[n] < len, "BOUNDS_CHECK");
        ESYS_ASSERT(newLocalToGlobalNodeLabels[n] >= 0, "BOUNDS_CHECK");
#endif
        globalToNewLocalNodeLabels[newLocalToGlobalNodeLabels[n]] = n;
        newLocalToGlobalNodeLabels[n] += min_id;
    }
    // create a new node file
    NodeFile* newNodeFile = new NodeFile(getDim(), m_mpiInfo);
    newNodeFile->allocTable(newNumNodes);
    if (len)
        newNodeFile->gather_global(&newLocalToGlobalNodeLabels[0], m_nodes);
    else
        newNodeFile->gather_global(NULL, m_nodes);

    delete m_nodes;
    m_nodes = newNodeFile;
    // relabel nodes of the elements
    relabelElementNodes(globalToNewLocalNodeLabels, min_id);
}

/// tries to reduce the number of colours for all element files
void FinleyDomain::createColoring(const IndexVector& dofMap)
{
    m_elements->createColoring(dofMap);
    m_faceElements->createColoring(dofMap);
    m_points->createColoring(dofMap);
    m_contactElements->createColoring(dofMap);
}

/// redistributes elements to minimize communication during assemblage
void FinleyDomain::optimizeElementOrdering()
{
    m_elements->optimizeOrdering();
    m_faceElements->optimizeOrdering();
    m_points->optimizeOrdering();
    m_contactElements->optimizeOrdering();
}

/// regenerates list of tags in use for node file and element files
void FinleyDomain::updateTagList()
{
    m_nodes->updateTagList();
    m_elements->updateTagList();
    m_faceElements->updateTagList();
    m_points->updateTagList();
    m_contactElements->updateTagList();
}


}  // end of namespace

