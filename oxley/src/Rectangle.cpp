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
*
*****************************************************************************/

#include <algorithm>
#include <cmath>
#include <ctime>
#include <exception>
#include <random>
#include <string>
#include <vector>

#include <escript/Assert.h>
#include <escript/Data.h>
#include <escript/DataFactory.h>
#include <escript/FunctionSpaceFactory.h>
#include <escript/Random.h>
#include <escript/Utils.h>

#include <oxley/AbstractAssembler.h>
#include <oxley/DefaultAssembler2D.h>
#include <oxley/domainhelpers.h>
#include <oxley/InitAlgorithms.h>
#include <oxley/Oxley.h>
#include <oxley/OxleyData.h>
#include <oxley/Rectangle.h>
#include <oxley/RefinementAlgorithms.h>
#include <oxley/RefinementType.h>
#include <oxley/RefinementZone.h>

// p4est headers will include MPI via sc.h when SC_ENABLE_MPI is defined
#include <p4est.h>
#include <p4est_algorithms.h>
#include <p4est_bits.h>
#include <p4est_connectivity.h>
#include <p4est_communication.h>
#include <p4est_extended.h>
#include <p4est_io.h>
#include <p4est_iterate.h>
#include <p4est_lnodes.h>
#include <p4est_ghost.h>
#include <p4est_vtk.h>

#include <unordered_map>
#include <array>

#include <sc_mpi.h>

// Include after p4est to get MPI_COMM_WORLD

#ifdef ESYS_HAVE_SILO
#include <silo.h>
#ifdef ESYS_MPI
#include <pmpio.h>
#endif
#endif

#ifdef ENABLE_OPENMP
#include <omp.h>
#endif

namespace bp = boost::python;

namespace oxley {


/**
   \brief
   Constructor with custom MPI communicator
*/
Rectangle::Rectangle(escript::JMPI jmpi, int order,
    dim_t n0, dim_t n1,
    double x0, double y0,
    double x1, double y1,
    const std::vector<double>& points,
    const std::vector<int>& tags,
    const TagMap& tagnamestonums,
    const std::vector<int>& refine_level):
    OxleyDomain(2, order, jmpi){

    // MPI communicator passed to base class constructor
    // Caller is responsible for ensuring MPI is initialized

    // n0/n1 are the number of BLOCKS (p4est trees) per axis; refine_level is the
    // subdivision applied to each block. It is either a single value (uniform, so
    // the base mesh has n0*2^L by n1*2^L elements) or one value per block (row-
    // major over (n0,n1), index i*n1+j), which lets blocks carry different levels
    // and so introduces hanging nodes at the block seams.
    if(n0 <= 0 || n1 <= 0)
        throw OxleyException("Number of blocks in each spatial dimension must be positive");
    const dim_t num_trees = n0 * n1;
    if(refine_level.empty())
        throw OxleyException("refine_level must not be empty");
    if(refine_level.size() != 1 && (dim_t) refine_level.size() != num_trees)
        throw OxleyException("refine_level must be a single value or one value per block (n0*n1)");
    for(size_t i = 0; i < refine_level.size(); i++)
        if(refine_level[i] < 0)
            throw OxleyException("refine_level must be non-negative");
    const int min_level = *std::min_element(refine_level.begin(), refine_level.end());
    const int max_level = *std::max_element(refine_level.begin(), refine_level.end());

    // Domain decomposition across MPI ranks is handled by p4est (see
    // p4est_partition below), not by a Cartesian d0 x d1 block grid.

    connectivity = new_rectangle_connectivity(n0, n1, false, false, x0, y0, x1, y1);

#ifdef OXLEY_ENABLE_DEBUG_CHECKS //These checks are turned off by default as they can be very timeconsuming
    std::cout << "In Rectangle() constructor..." << std::endl;
    std::cout << "Checking connectivity ... ";
    if(!p4est_connectivity_is_valid(connectivity))
        std::cout << "broken" << std::endl;
    else
        std::cout << "OK" << std::endl;
#endif

    // Create a p4est - use the custom communicator.
    // fill_uniform + min_level builds a uniform base mesh where every block is
    // subdivided min_level times. min_quadrants MUST be 0: it is p4est's
    // PER-PROCESSOR minimum, so a positive value forces extra refinement under
    // MPI and makes the mesh depend on the rank count. (A6.)
    p4est_locidx_t min_quadrants = 0;
    int fill_uniform = 1;

    p4est = p4est_new_ext(m_mpiInfo->comm, connectivity, min_quadrants,
            min_level, fill_uniform, sizeof(quadrantData), init_rectangle_data, (void *) &forestData);

    // If the blocks do not all share the same level, refine each block up to its
    // own target level and 2:1-balance the base mesh. block_levels is indexed by
    // p4est tree id, so map each tree to its block (i,j) via its lower-left
    // connectivity vertex (the trees are Morton-ordered, not row-major).
    if(max_level > min_level) {
        const double dxb = (x1-x0)/n0, dyb = (y1-y0)/n1;
        forestData.block_levels.assign(num_trees, min_level);
        for(p4est_topidx_t t = 0; t < num_trees; t++) {
            const p4est_topidx_t v = connectivity->tree_to_vertex[P4EST_CHILDREN*t + 0];
            const double vx = connectivity->vertices[3*v + 0];
            const double vy = connectivity->vertices[3*v + 1];
            int bi = (int) std::lround((vx - x0)/dxb);
            int bj = (int) std::lround((vy - y0)/dyb);
            if(bi < 0) bi = 0; else if(bi >= n0) bi = n0-1;
            if(bj < 0) bj = 0; else if(bj >= n1) bj = n1-1;
            forestData.block_levels[t] = refine_level[(size_t) bi*n1 + bj];
        }
        const int recursive = 1;
        p4est_refine_ext(p4est, recursive, max_level, refine_to_block_level,
                         init_rectangle_data, NULL);
        p4est_balance_ext(p4est, P4EST_CONNECT_FULL, init_rectangle_data, NULL);
    }

#ifdef OXLEY_ENABLE_DEBUG_CHECKS //These checks are turned off by default as they can be very timeconsuming
    std::cout << "Checking p4est ... ";
    if(!p4est_is_valid(p4est))
        std::cout << "broken" << std::endl;
    else
        std::cout << "OK" << std::endl;
#endif

    // Distribute the p4est across the processors. This MUST happen before the
    // lnodes are built: p4est_partition moves quadrants between ranks, and an
    // lnodes built beforehand still describes the old partition (its
    // num_local_elements and element_nodes refer to quadrants this rank no longer
    // owns), so every later leaf walk runs off the end of element_nodes. It went
    // unnoticed while every mesh was uniform, because then the partition is
    // already balanced and p4est_partition is a no-op. (Milestone B.)
    int allow_coarsening = 0;
    p4est_partition(p4est, allow_coarsening, NULL);

    // Nodes numbering
    p4est_ghost_t * ghost = p4est_ghost_new(p4est, P4EST_CONNECT_FULL);
    nodes = p4est_lnodes_new(p4est, ghost, 1);
    p4est_ghost_destroy(ghost);

    // This information is needed by the assembler
    m_NE[0] = n0;
    m_NE[1] = n1;
    m_NX[0] = (x1-x0)/n0;
    m_NX[1] = (y1-y0)/n1;
    m_NN[0] = n0;
    m_NN[1] = n1;

    // Record the physical dimensions of the domain and the location of the origin
    forestData.m_origin[0] = x0;
    forestData.m_origin[1] = y0;
    forestData.m_lxy[0] = x1;
    forestData.m_lxy[1] = y1;
    forestData.m_length[0] = x1-x0;
    forestData.m_length[1] = y1-y0;
    forestData.m_NX[0] = (x1-x0)/n0;
    forestData.m_NX[1] = (y1-y0)/n1;

    // Periodic boundaries are not implemented: the connectivity above is built
    // non-periodic. forestData.periodic stays at its default (false).

    // Find the grid spacing for each level of refinement in the mesh
#pragma omp parallel for
    for(int i = 0; i<=P4EST_MAXLEVEL; i++){
        double numberOfSubDivisions = (p4est_qcoord_t) (1 << (P4EST_MAXLEVEL - i));
        forestData.m_dx[0][i] = forestData.m_NX[0] / (numberOfSubDivisions);
        forestData.m_dx[1][i] = forestData.m_NX[1] / (numberOfSubDivisions);
    }

    // max levels of refinement
    forestData.max_levels_refinement = MAXREFINEMENTLEVELS;

    // element order
    m_order = order;

    // Number of dimensions
    m_numDim=2;

    // (p4est_partition happens above, before the lnodes are built)

    // Number the nodes
    updateNodeIncrements();
    renumberNodes();
    updateRowsColumns();
    updateNodeDistribution();
    updateElementIds();
    updateFaceOffset();
    updateFaceElementCount();
    updateQuadrantIDinformation();

    // Tags
    populateSampleIds();
    for (TagMap::const_iterator i = tagnamestonums.begin(); i != tagnamestonums.end(); i++) {
        setTagMap(i->first, i->second);
    }

    // Dirac points and tags
    addPoints(points, tags);

    // To prevent segmentation faults when using numpy ndarray
#ifdef ESYS_HAVE_BOOST_NUMPY
    Py_Initialize();
    boost::python::numpy::initialize();
#endif

#ifdef ESYS_HAVE_PASO

    /// local array length shared
    dim_t local_length = 0;

    /// list of the processors sharing values with this processor
    std::vector<int> neighbour = {};

    /// offsetInShared[i] points to the first input value in array shared
    /// for processor i. Has length numNeighbors+1
    std::vector<index_t> offsetInShared = {0};

    /// list of the (local) components which are shared with other processors.
    /// Has length numSharedComponents
    index_t* shared = {};

    /// = offsetInShared[numNeighbours]
    dim_t numSharedComponents = 0;

    IndexVector sendShared, recvShared;

    createPasoConnector(neighbour, offsetInShared, offsetInShared, sendShared, recvShared);

#endif

    oxleytimer.toc("Class initialised");
}

Rectangle::Rectangle(const oxley::Rectangle& R, int order):
    OxleyDomain(2, order){

    m_mpiInfo = R.m_mpiInfo;

    p4est=p4est_copy(R.p4est,1);

#ifdef OXLEY_ENABLE_DEBUG_CHECKS //These checks are turned off by default as they can be very timeconsuming
    std::cout << "In Rectangle() constructor..." << std::endl;
    std::cout << "Checking connectivity ... ";
    if(!p4est_connectivity_is_valid(connectivity))
        std::cout << "broken" << std::endl;
    else
        std::cout << "OK" << std::endl;
    std::cout << "Checking p4est ... ";
    if(!p4est_is_valid(p4est))
        std::cout << "broken" << std::endl;
    else
        std::cout << "OK" << std::endl;
#endif

    nodes=R.nodes;

    // This information is needed by the assembler
    m_NE[0] = R.m_NE[0];
    m_NE[1] = R.m_NE[1];
    m_NX[0] = R.m_NX[0];
    m_NX[1] = R.m_NX[1];
    m_NN[0] = R.m_NN[0];
    m_NN[1] = R.m_NN[1];

    forestData.m_origin[0] = R.forestData.m_origin[0];
    forestData.m_origin[1] = R.forestData.m_origin[1];
    forestData.m_lxy[0] = R.forestData.m_lxy[0];
    forestData.m_lxy[1] = R.forestData.m_lxy[1];
    forestData.m_length[0] = R.forestData.m_length[0];
    forestData.m_length[1] = R.forestData.m_length[1];
    forestData.m_NX[0] = R.forestData.m_NX[0];
    forestData.m_NX[1] = R.forestData.m_NX[1];

    connectivity=new_rectangle_connectivity(m_NE[0], m_NE[1], false, false, 
                                            forestData.m_origin[0], forestData.m_origin[1], 
                                            forestData.m_lxy[0], forestData.m_lxy[1]);

    // Whether or not we have periodic boundaries
    forestData.periodic[0] = R.forestData.periodic[0];
    forestData.periodic[1] = R.forestData.periodic[1];

    // Find the grid spacing for each level of refinement in the mesh
#pragma omp parallel for
    for(int i = 0; i<=P4EST_MAXLEVEL; i++){
        double numberOfSubDivisions = (p4est_qcoord_t) (1 << (P4EST_MAXLEVEL - i));
        forestData.m_dx[0][i] = forestData.m_NX[0] / (numberOfSubDivisions);
        forestData.m_dx[1][i] = forestData.m_NX[1] / (numberOfSubDivisions);
    }

    // max levels of refinement
    forestData.max_levels_refinement = MAXREFINEMENTLEVELS;

    // Update the user_data pointer in p4est 
    p4est->user_pointer=&forestData;

    // element order
    m_order = R.m_order;

    // Number of dimensions
    m_numDim=2;

    // Distribute the p4est across the processors. MUST precede the lnodes: see
    // the note in the main constructor -- an lnodes built before the partition
    // describes the old one.
    int allow_coarsening = 0;
    p4est_partition(p4est, allow_coarsening, NULL);

    // lnodes
    p4est_ghost_t * ghost = p4est_ghost_new(p4est, P4EST_CONNECT_FULL);
    nodes = p4est_lnodes_new(p4est, ghost, 1);
    p4est_ghost_destroy(ghost);

    // Number the nodes
    updateNodeIncrements();
    renumberNodes();
    updateRowsColumns();
    updateNodeDistribution();
    updateElementIds();
    updateFaceOffset();
    updateFaceElementCount();
updateQuadrantIDinformation();

    // Tags
    populateSampleIds();
    m_tagMap=R.m_tagMap;
    
    // Dirac points and tags
    m_diracPoints=R.m_diracPoints;

    // To prevent segmentation faults when using numpy ndarray
#ifdef ESYS_HAVE_BOOST_NUMPY
    Py_Initialize();
    boost::python::numpy::initialize();
#endif

    oxleytimer.toc("Class initialised");
}

/**
   \brief
   Destructor.
*/
Rectangle::~Rectangle(){
#ifdef OXLEY_ENABLE_DEBUG_CHECKS
    std::cout << "In Rectangle() destructor" << std::endl;
    std::cout << "checking p4est ... ";
    if(!p4est_is_valid(p4est))
        std::cout << "broken" << std::endl;
    else
        std::cout << "OK" << std::endl;
    std::cout << "checking connectivity ... ";
    if(!p4est_connectivity_is_valid(connectivity))
        std::cout << "broken" << std::endl;
    else
        std::cout << "OK" << std::endl;
#endif
    if (m_ghost) { p4est_ghost_destroy(m_ghost); m_ghost = nullptr; }
}

/**
   \brief
   returns a description for this domain
*/
std::string Rectangle::getDescription() const{
    return "oxley::rectangle";
}


/**
   \brief
   writes the current mesh to a file with the given name
   \param filename The name of the file to write to
*/
void Rectangle::write(const std::string& filename) const
{
    throw OxleyException("write: not supported");
}

void Rectangle::interpolateAcross(escript::Data& target, const escript::Data& source) const
{
    const Rectangle *other = dynamic_cast<const Rectangle *>(target.getDomain().get());
    if (other == NULL)
        throw OxleyException("Invalid interpolation: Domains must both be instances of Rectangle");
    //shouldn't ever happen, but I want to know if it does
    if (other == this)
        throw OxleyException("interpolateAcross: this domain is the target");
        
    validateInterpolationAcross(source.getFunctionSpace().getTypeCode(),
            *(target.getDomain().get()), target.getFunctionSpace().getTypeCode());
    int fsSource = source.getFunctionSpace().getTypeCode();
    int fsTarget = target.getFunctionSpace().getTypeCode();

    std::stringstream msg;
    msg << "Invalid interpolation: interpolation not implemented for function space "
        << functionSpaceTypeAsString(fsSource)
        << " -> "
        << functionSpaceTypeAsString(fsTarget);

#ifdef OXLEY_ENABLE_DEBUG_INTERPOLATE_ACROSS
    std::cout << "InterpolateAcross" << std::endl;
    std::cout << "Doing " << functionSpaceTypeAsString(fsSource)
                          << " -> "
                          << functionSpaceTypeAsString(fsTarget)
                          << std::endl;
#endif

    switch (fsSource) {
        case Nodes:
            switch (fsTarget) {
                case Nodes:
                case ReducedNodes:
                case DegreesOfFreedom:
                case ReducedDegreesOfFreedom:
                    interpolateNodesToNodesFiner(source, target, *other);
                    return;
                case Elements:
                    interpolateNodesToElementsFiner(source, target, *other);
                    return;
                default:
                    throw OxleyException(msg.str());
            }
            break;
        case Elements:
            switch (fsTarget) {
                case Elements:
                    interpolateElementsToElementsFiner(source, target, *other);
                    return;
                default:
                    throw OxleyException(msg.str());
            }
            break;
        case ReducedElements:
            switch (fsTarget) {
                case Elements:
                    interpolateReducedToElementsFiner(source, target, *other);
                    return;
                default:
                    throw OxleyException(msg.str());
            }
            break;
        case DegreesOfFreedom:
            switch (fsTarget) {
                case Nodes:
                case ReducedNodes:
                case DegreesOfFreedom:
                case ReducedDegreesOfFreedom:
                    interpolateNodesToNodesFiner(source, target, *other);
                    return;
                case Elements:
                    interpolateNodesToElementsFiner(source, target, *other);
                    return;
                default:
                    throw OxleyException(msg.str());
            }
            break;
        default:
            throw OxleyException(msg.str());
    }
}

void Rectangle::validateInterpolationAcross(int fsType_source, const escript::AbstractDomain& domain, int fsType_target) const
{
    const Rectangle *other = dynamic_cast<const Rectangle *>(&domain);
    if (other == NULL)
        throw OxleyException("Invalid interpolation: domains must both be instances of oxley::Rectangle");

    // TODO
    // if(!p4est_is_equal(borrow_p4est, other->borrow_p4est(), 0))
    //     throw OxleyException("Invalid interpolation: domains have different p4ests");

    // if(!p4est_connectivity_is_equivalent(borrow_connectivity(),other->borrow_connectivity()))
    //     throw OxleyException("Invalid interpolation: domains have different connectivities");


}

// downward interpolation
void Rectangle::interpolateNodesToNodesCoarser(const escript::Data& source, escript::Data& target, const Rectangle& other) const
{

}

void Rectangle::interpolateNodesToNodesFiner(const escript::Data& source, escript::Data& target, const Rectangle& other) const
{
    if(source.isComplex())
        interpolateNodesToNodesWorker<cplx_t>(source, target, other);
    else
        interpolateNodesToNodesWorker<real_t>(source, target, other);
}

// upward interpolation
template <typename S>
void Rectangle::interpolateNodesToNodesWorker(const escript::Data& source, escript::Data& target, const Rectangle& other) const
{
    // Cross-domain (refinement) interpolation is not yet reimplemented on the
    // lnodes numbering; probeInterpolationAcross() returns false so this is
    // never reached. (Milestone B.)
    throw OxleyException("interpolateNodesToNodesWorker: not implemented for the lnodes numbering yet");
}

void Rectangle::interpolateNodesToElementsFiner(const escript::Data& source, escript::Data& target, const Rectangle& other)  const
{
    throw OxleyException("Not yet implemented"); //TODO
}

void Rectangle::interpolateElementsToElementsCoarser(const escript::Data& source, escript::Data& target, const Rectangle& other)  const
{
    
}

void Rectangle::interpolateElementsToElementsFiner(const escript::Data& source, escript::Data& target, const Rectangle& other)  const
{
    if(source.isComplex())
        interpolateElementsToElementsWorker<cplx_t>(source, target, other);
    else
        interpolateElementsToElementsWorker<real_t>(source, target, other);
}

template <typename S>
void Rectangle::interpolateElementsToElementsWorker(const escript::Data& source, escript::Data& target, const Rectangle& other)  const
{
    // Cross-domain (refinement) interpolation is not yet reimplemented on the
    // lnodes numbering; probeInterpolationAcross() returns false so this is
    // never reached. (Milestone B.)
    throw OxleyException("interpolateElementsToElementsWorker: not implemented for the lnodes numbering yet");
}

void Rectangle::interpolateReducedToElementsFiner(const escript::Data& source, escript::Data& target, const Rectangle& other)  const
{

}

void Rectangle::interpolateReducedToReducedFiner(const escript::Data& source, escript::Data& target, const Rectangle& other)  const
{

}

void Rectangle::setToNormal(escript::Data& out) const
{
    if (out.getFunctionSpace().getTypeCode() == FaceElements) {
        out.requireWrite();
#pragma omp parallel
        {
            if (m_faceOffset[0] > -1) {
#pragma omp for nowait
                for (index_t k=0; k<NodeIDsLeft.size(); k++) {
                    double* o = out.getSampleDataRW(m_faceOffset[0]+k);
                    // set vector at two quadrature points
                    *o++ = -1.;
                    *o++ = 0.;
                    *o++ = -1.;
                    *o = 0.;
                }
            }

            if (m_faceOffset[1] > -1) {
#pragma omp for nowait
                for (index_t k=0; k<NodeIDsRight.size(); k++) {
                    double* o = out.getSampleDataRW(m_faceOffset[1]+k);
                    // set vector at two quadrature points
                    *o++ = 1.;
                    *o++ = 0.;
                    *o++ = 1.;
                    *o = 0.;
                }
            }

            if (m_faceOffset[2] > -1) {
#pragma omp for nowait
                for (index_t k=0; k<NodeIDsBottom.size(); k++) {
                    double* o = out.getSampleDataRW(m_faceOffset[2]+k);
                    // set vector at two quadrature points
                    *o++ = 0.;
                    *o++ = -1.;
                    *o++ = 0.;
                    *o = -1.;
                }
            }

            if (m_faceOffset[3] > -1) {
#pragma omp for nowait
                for (index_t k=0; k<NodeIDsTop.size(); k++) {
                    double* o = out.getSampleDataRW(m_faceOffset[3]+k);
                    // set vector at two quadrature points
                    *o++ = 0.;
                    *o++ = 1.;
                    *o++ = 0.;
                    *o = 1.;
                }
            }
        } // end of parallel section
    } else if (out.getFunctionSpace().getTypeCode() == ReducedFaceElements) {
        out.requireWrite();
#pragma omp parallel
        {
            if (m_faceOffset[0] > -1) {
#pragma omp for nowait
                for (index_t k=0; k<NodeIDsLeft.size(); k++) {
                    double* o = out.getSampleDataRW(m_faceOffset[0]+k);
                    *o++ = -1.;
                    *o = 0.;
                }
            }

            if (m_faceOffset[1] > -1) {
#pragma omp for nowait
                for (index_t k=0; k<NodeIDsRight.size(); k++) {
                    double* o = out.getSampleDataRW(m_faceOffset[1]+k);
                    *o++ = 1.;
                    *o = 0.;
                }
            }

            if (m_faceOffset[2] > -1) {
#pragma omp for nowait
                for (index_t k=0; k<NodeIDsBottom.size(); k++) {
                    double* o = out.getSampleDataRW(m_faceOffset[2]+k);
                    *o++ = 0.;
                    *o = -1.;
                }
            }

            if (m_faceOffset[3] > -1) {
#pragma omp for nowait
                for (index_t k=0; k<NodeIDsTop.size(); k++) {
                    double* o = out.getSampleDataRW(m_faceOffset[3]+k);
                    *o++ = 0.;
                    *o = 1.;
                }
            }
        } // end of parallel section

    } else {
        std::stringstream msg;
        msg << "setToNormal: invalid function space type "
            << out.getFunctionSpace().getTypeCode();
        throw ValueError(msg.str());
    }

    #ifdef OXLEY_ENABLE_DEBUG_SETTONORMAL
        std::cout << "setToNormal:" << std::endl;
        out.print();
    #endif
}

void Rectangle::setToSize(escript::Data& out) const
{
    if (out.getFunctionSpace().getTypeCode() == Elements
        || out.getFunctionSpace().getTypeCode() == ReducedElements)
    {
        out.requireWrite();

        // Find the maximum level of refinement in the mesh
        int max_level = 0;
        for(p4est_topidx_t tree = p4est->first_local_tree; tree <= p4est->last_local_tree; tree++) {
            p4est_tree_t * tree_t = p4est_tree_array_index(p4est->trees, tree);
            max_level = tree_t->maxlevel > max_level ? tree_t->maxlevel : max_level;
        }
        // Work out the size at each level
        std::vector<double> size_vect(max_level+1, -1.0);
        for(int i = 0 ; i <= max_level ; i++)
        {
            size_vect[i] = sqrt((forestData.m_dx[0][P4EST_MAXLEVEL-i]*forestData.m_dx[0][P4EST_MAXLEVEL-i]
                                                    +forestData.m_dx[1][P4EST_MAXLEVEL-i]*forestData.m_dx[1][P4EST_MAXLEVEL-i]));
        }

        const dim_t numQuad = out.getNumDataPointsPerSample();
        long id = 0;   // running local leaf index (lnodes / element sample order)
        for (p4est_topidx_t t = p4est->first_local_tree; t <= p4est->last_local_tree; t++)
        {
            p4est_tree_t * currenttree = p4est_tree_array_index(p4est->trees, t);
            sc_array_t * tquadrants = &currenttree->quadrants;
            p4est_locidx_t Q = (p4est_locidx_t) tquadrants->elem_count;
            for (int q = 0; q < Q; ++q, ++id)
            {
                p4est_quadrant_t * quad = p4est_quadrant_array_index(tquadrants, q);
                int l = quad->level;
                const double size = size_vect[l];
                double* o = out.getSampleDataRW(id);
                std::fill(o, o+numQuad, size);
            }
        }
    }
    else if (out.getFunctionSpace().getTypeCode() == FaceElements
            || out.getFunctionSpace().getTypeCode() == ReducedFaceElements) 
    {
        out.requireWrite();
        const dim_t numQuad=out.getNumDataPointsPerSample();

        if (m_faceOffset[0] > -1) {
            for (index_t k=0; k<NodeIDsLeft.size(); k++) {
                borderNodeInfo tmp = NodeIDsLeft[k];

                double* o = out.getSampleDataRW(m_faceOffset[0]+k);
                std::fill(o, o+numQuad, forestData.m_dx[1][P4EST_MAXLEVEL-tmp.level]);
            }
        }

        if (m_faceOffset[1] > -1) {
            for (index_t k=0; k<NodeIDsRight.size(); k++) {
                borderNodeInfo tmp = NodeIDsRight[k];
                double* o = out.getSampleDataRW(m_faceOffset[1]+k);
                std::fill(o, o+numQuad, forestData.m_dx[1][P4EST_MAXLEVEL-tmp.level]);
            }
        }

        if (m_faceOffset[2] > -1) {
            for (index_t k=0; k<NodeIDsBottom.size(); k++) {
                borderNodeInfo tmp = NodeIDsBottom[k];
                double* o = out.getSampleDataRW(m_faceOffset[2]+k);
                std::fill(o, o+numQuad, forestData.m_dx[0][P4EST_MAXLEVEL-tmp.level]);
            }
        }

        if (m_faceOffset[3] > -1) {
            for (index_t k=0; k<NodeIDsTop.size(); k++) {
                borderNodeInfo tmp = NodeIDsTop[k];
                double* o = out.getSampleDataRW(m_faceOffset[3]+k);
                std::fill(o, o+numQuad, forestData.m_dx[0][P4EST_MAXLEVEL-tmp.level]);
            }
        }
    } else {
        std::stringstream msg;
        msg << "setToSize: invalid function space type "
            << out.getFunctionSpace().getTypeCode();
        throw ValueError(msg.str());
    }
}

bool Rectangle::ownSample(int fsType, index_t id) const
{
    if (getMPISize()==1)
        return true;

    switch (fsType) {
        case Nodes:
        case ReducedNodes: // FIXME: reduced
            return (m_dofMap[id] < getNumDOF());
        case DegreesOfFreedom:
        case ReducedDegreesOfFreedom:
            return true;
        case Elements:
        case ReducedElements:
        case FaceElements:
        case ReducedFaceElements:
        case Points:
            // p4est partitions leaves (and hence boundary faces) uniquely across
            // ranks, and Dirac points are claimed by a single owner (addPoints),
            // so every local element/face/point sample is owned by this rank.
            return true;
        default:
            break;
    }

    std::stringstream msg;
    msg << "ownSample: invalid function space type " << fsType;
    throw ValueError(msg.str());
}


dim_t Rectangle::getNumDataPointsGlobal() const
{
    // total number of (owned) nodes across all ranks
    if(!nodes) return 0;
    dim_t total = 0;
    for(int r = 0; r < m_mpiInfo->size; ++r)
        total += (dim_t) nodes->global_owned_count[r];
    return total;
}

void Rectangle::dump(const std::string& fileName) const
{
#ifdef ESYS_HAVE_SILO
    
    // Add the suffix to the filename if required
    std::string fn = fileName;
    if (fileName.length() < 6 || fileName.compare(fileName.length()-5, 5, ".silo") != 0) {
        fn+=".silo";
    }

    // int driver=DB_HDF5;

    // Silo file pointer
    DBfile* dbfile = NULL; 

    // The coordinate arrays
    float *pNodex = nullptr;
    float *pNodey = nullptr;
    long int *pNode_ids = nullptr;
    double * pValues = nullptr;

    // node coordinates and element connectivity come from the lnodes-based
    // mesh-access view (no coordinate hashing).
    const MeshAccess m = getMeshAccess();
    const int V = m.nodesPerElement;

    pNodex = new float[m.numNodes];
    pNodey = new float[m.numNodes];
    pNode_ids = new long int [m.numNodes];

    for(long i = 0; i < m.numNodes; ++i)
    {
        pNodex[i]    = (float) m.nodeCoords[(size_t) i*m.numDim + 0];
        pNodey[i]    = (float) m.nodeCoords[(size_t) i*m.numDim + 1];
        pNode_ids[i] = m.nodeLnodesId[i];
    }

    // Array of the coordinate arrays
    float * pCoordinates[2];
    pCoordinates[0]=pNodex;
    pCoordinates[1]=pNodey;

    // Create the file
    dbfile = DBCreate(fn.c_str(), DB_CLOBBER, DB_LOCAL, getDescription().c_str(), DB_HDF5);
    if (!dbfile)
        throw escript::IOError("dump: Could not create Silo file");

    // create the nodelist (Silo quad winding: z-order 0,2,3,1)
    std::vector<int> nodelist;
    for(long e = 0; e < m.numElements; ++e)
    {
        nodelist.push_back((int) m.elementNodes[(size_t) e*V + 0]);
        nodelist.push_back((int) m.elementNodes[(size_t) e*V + 2]);
        nodelist.push_back((int) m.elementNodes[(size_t) e*V + 3]);
        nodelist.push_back((int) m.elementNodes[(size_t) e*V + 1]);
    }

    int* nodelistarray = &nodelist[0];

    // write mesh
    int lnodelist = nodelist.size();
    int shapesize[] = {4};
    int shapecounts[] = {lnodelist/4};
    int nshapetypes = 1;
    int shapetype[1] = {DB_ZONETYPE_QUAD};

    // This is deprecated
    // DBPutZonelist(dbfile, "quads", getNumElements(), 2, nodelistarray, lnodelist, 0, 
    //                 shapesize, shapecounts, nshapetypes);
    DBPutZonelist2(dbfile, "quads", getNumElements(), 2, nodelistarray, lnodelist, 0,
                0, 0, shapetype, shapesize, shapecounts, nshapetypes, NULL);
        

    DBPutUcdmesh(dbfile, "mesh", 2, NULL, pCoordinates, getNumNodes(), getNumElements(), 
                    "quads", NULL, DB_FLOAT, NULL);

    // Coordinates
    DBPutPointmesh(dbfile, "nodes", 2, pCoordinates, getNumNodes(), DB_FLOAT, NULL) ;

    // Node IDs
    DBPutPointvar1(dbfile, "id", "nodes", pNode_ids, getNumNodes(), DB_LONG, NULL);

    DBClose(dbfile);

    delete [] pNodex;
    delete [] pNodey;
    delete [] pNode_ids;

#else // ESYS_HAVE_SILO
    throw OxleyException("dump: escript was not compiled with Silo enabled");
#endif
}

const dim_t* Rectangle::borrowSampleReferenceIDs(int fsType) const
{
    switch (fsType) {
        case Nodes:
        case ReducedNodes:
            return &m_nodeId[0];
        case DegreesOfFreedom:
        case ReducedDegreesOfFreedom:
            return &myRows[0];
        case Elements:
        case ReducedElements:
            return &m_elementId[0];
        case FaceElements:
        case ReducedFaceElements:
            return &m_faceId[0];
        case Points:
            return &m_diracPointNodeIDs[0];
        default:
            std::stringstream msg;
            msg << "borrowSampleReferenceIDs: invalid function space type " << fsType;
            throw ValueError(msg.str());
    }    
}

void Rectangle::writeToVTK(std::string filename, bool writeMesh) const
{
    // Write to file
    const char * name = filename.c_str();
    if(writeMesh)
    {
        p4est_vtk_write_file(p4est, NULL, name);
    }
    else
    {
        // Create the context for the VTK file
        p4est_vtk_context_t * context = p4est_vtk_context_new(p4est, name);

        // Continuous point data
        p4est_vtk_context_set_continuous(context, true);

        // Set the scale
        p4est_vtk_context_set_scale(context, 1.0);

        // Write the header
        context = p4est_vtk_write_header(context);

        // Get the point and cell data together
        p4est_locidx_t numquads = p4est->local_num_quadrants;

        //  Info
        sc_array_t * quadTag = sc_array_new_count(sizeof(double), numquads);
        p4est_iterate(p4est, NULL, (void *) quadTag, getQuadTagVector, NULL, NULL);
        sc_array_t * xcoord = sc_array_new_count(sizeof(double), numquads);
        p4est_iterate(p4est, NULL, (void *) xcoord, getXCoordVector, NULL, NULL);
        sc_array_t * ycoord = sc_array_new_count(sizeof(double), numquads);
        p4est_iterate(p4est, NULL, (void *) ycoord, getYCoordVector, NULL, NULL);
        // sc_array_t * NodeNumber = sc_array_new_count(sizeof(double), numquads);
        // p4est_iterate(p4est, NULL, (void *) NodeNumber, getNodeNumber, NULL, NULL);

        // Write the cell Data
#ifdef OXLEY_ENABLE_DEBUG
        context = p4est_vtk_write_cell_dataf(context,1,1,0,0,3,0,"tag",quadTag,"x",xcoord,"y",ycoord,context);
#else
        context = p4est_vtk_write_cell_dataf(context,0,0,0,0,3,0,"tag",quadTag,"x",xcoord,"y",ycoord,context);
#endif
        if(context == NULL)
            throw OxleyException("Error writing cell data");

        // Write the point Data
        context = p4est_vtk_write_point_dataf(context, 0, 0, context);
        if(context == NULL)
            throw OxleyException("Error writing point data");

        // Write the footer
        if(p4est_vtk_write_footer(context)) // The context is destroyed by this function
                throw OxleyException("Error writing footer.");

        // Cleanup
        sc_array_reset(quadTag);
        sc_array_destroy(quadTag);
        sc_array_reset(xcoord);
        sc_array_destroy(xcoord);
        sc_array_reset(ycoord);
        sc_array_destroy(ycoord);
    }
}

#ifdef ESYS_HAVE_TRILINOS
void Rectangle::saveMesh(std::string filename) 
{
    std::string fnames=filename+".p4est";
    std::string cnames=filename+".conn";

    const char * fname=fnames.c_str();
    const char * cname=cnames.c_str();

    p4est_deflate_quadrants(p4est, NULL);

#ifdef ESYS_MPI
    if(escript::getMPIRankWorld()==0)
    {
#endif
        int retval = p4est_connectivity_save(cname, connectivity)==0;
        ESYS_ASSERT(retval!=0,"Failed to save connectivity");
        int save_partition = 0;
        int save_data = 1;
        p4est_save_ext(fname, p4est, save_data, save_partition); // Should abort on file error
        // p4est_save(fname,p4est,1);
#ifdef ESYS_MPI
    }
#endif
}

void Rectangle::loadMesh(std::string filename) 
{
    std::string fnames=filename+".p4est";
    std::string cnames=filename+".conn";

    const char * fname=fnames.c_str();
    const char * cname=cnames.c_str();

    int load_data=true;
    int autopartition=true;
    int broadcasthead=false;

    // Delete the old structure
    p4est_connectivity_destroy(connectivity);
    p4est_destroy(p4est);

    // Load the new information
    // connectivity=p4est_connectivity_load(cname, NULL);
    // ESYS_ASSERT(p4est_connectivity_is_valid(connectivity), "Invalid connectivity file");
    p4est=p4est_load_ext(fname, m_mpiInfo->comm, sizeof(quadrantData), load_data, 
                    autopartition, broadcasthead, &forestData, &connectivity);
    // p4est = p4est_load(fname, m_mpiInfo->comm, sizeof(quadrantData), 1, NULL, &connectivity);
    ESYS_ASSERT(p4est_is_valid(p4est),"Invalid p4est file");

    // Update the nodes
    p4est_lnodes_destroy(nodes);
    p4est_ghost_t * ghost = p4est_ghost_new(p4est, P4EST_CONNECT_FULL);
    nodes = p4est_lnodes_new(p4est, ghost, 1);
    p4est_ghost_destroy(ghost);

    // Update rectangle
    if(autoMeshUpdates)
        updateMesh();

    // Need to update these now that the mesh has changed
    z_needs_update=true;
    iz_needs_update=true;
}

void Rectangle::refineMesh(std::string algorithmname)
{
    oxleytimer.toc("refineMesh...");

    z_needs_update=true;
    iz_needs_update=true;

    p4estData * pForestData;
    pForestData = &forestData;
    p4est->user_pointer = pForestData;

    if(!algorithmname.compare("uniform"))
    {
        p4est_refine_ext(p4est, true, -1, refine_uniform, init_rectangle_data, NULL);
        p4est_balance_ext(p4est, P4EST_CONNECT_FULL, init_rectangle_data, NULL);
    }
    else if(!algorithmname.compare("MARE2DEM") || !algorithmname.compare("mare2dem"))
    {
        if(adaptive_refinement == true)
        {
            p4est_refine_ext(p4est, true, -1, refine_mare2dem, init_rectangle_data, NULL);
            p4est_balance_ext(p4est, P4EST_CONNECT_FULL, init_rectangle_data, NULL);
        }
        else
        {
#ifdef OXLEY_ENABLE_DEBUG
            std::cout << "Warning: Adaptive mesh refinement is disabled." << std::endl;
#endif
        }
    }
    else {
        throw OxleyException("Unknown refinement algorithm name.");
    }

    // Make sure that nothing went wrong
#ifdef OXLEY_ENABLE_DEBUG
    if(!p4est_is_valid(p4est))
        throw OxleyException("p4est broke during refinement");
    if(!p4est_connectivity_is_valid(connectivity))
        throw OxleyException("connectivity broke during refinement");
#endif

    bool partition_for_coarsening = true;
    p4est_partition_ext(p4est, partition_for_coarsening, NULL);

    // Update the nodes
    p4est_lnodes_destroy(nodes);
    p4est_ghost_t * ghost = p4est_ghost_new(p4est, P4EST_CONNECT_FULL);
    nodes = p4est_lnodes_new(p4est, ghost, 1);
    p4est_ghost_destroy(ghost);

    // Update
    if(autoMeshUpdates)
        updateMesh();

    oxleytimer.toc("refineMesh...done");
}

void Rectangle::refineBoundary(std::string boundaryname, double dx)
{
    oxleytimer.toc("refineBoundary...");

    z_needs_update=true;
    iz_needs_update=true;

    forestData.refinement_depth = dx;

    if(!boundaryname.compare("top") || !boundaryname.compare("Top")
        || !boundaryname.compare("t") || !boundaryname.compare("T")
        || !boundaryname.compare("TOP"))
    {
        p4est_refine_ext(p4est, true, -1, refine_north, init_rectangle_data, NULL);
        p4est_balance_ext(p4est, P4EST_CONNECT_FULL, init_rectangle_data, NULL);
    } 
    else if(!boundaryname.compare("bottom") || !boundaryname.compare("Bottom")
        || !boundaryname.compare("b") || !boundaryname.compare("B")
        || !boundaryname.compare("BOTTOM"))
    {
        p4est_refine_ext(p4est, true, -1, refine_south, init_rectangle_data, NULL);
        p4est_balance_ext(p4est, P4EST_CONNECT_FULL, init_rectangle_data, NULL);
    }
    else if(!boundaryname.compare("left") || !boundaryname.compare("Left")
        || !boundaryname.compare("l") || !boundaryname.compare("L")
        || !boundaryname.compare("LEFT"))
    {
        p4est_refine_ext(p4est, true, -1, refine_west, init_rectangle_data, NULL);
        p4est_balance_ext(p4est, P4EST_CONNECT_FULL, init_rectangle_data, NULL);
    }
    else if(!boundaryname.compare("right") || !boundaryname.compare("Right")
        || !boundaryname.compare("r") || !boundaryname.compare("R")
        || !boundaryname.compare("RIGHT"))
    {
        p4est_refine_ext(p4est, true, -1, refine_east, init_rectangle_data, NULL);
        p4est_balance_ext(p4est, P4EST_CONNECT_FULL, init_rectangle_data, NULL);  
    }
    else {
        throw OxleyException("Unknown boundary name. Please try 'top', 'bottom', 'left' or 'right'.");
    }

    // Make sure that nothing went wrong
#ifdef OXLEY_ENABLE_DEBUG
    if(!p4est_is_valid(p4est))
        throw OxleyException("p4est broke during refinement");
    if(!p4est_connectivity_is_valid(connectivity))
        throw OxleyException("connectivity broke during refinement");
#endif

    bool partition_for_coarsening = true;
    p4est_partition_ext(p4est, partition_for_coarsening, NULL);

    // Update the nodes
    p4est_lnodes_destroy(nodes);
    p4est_ghost_t * ghost = p4est_ghost_new(p4est, P4EST_CONNECT_FULL);
    nodes = p4est_lnodes_new(p4est, ghost, 1);
    p4est_ghost_destroy(ghost);

    // Update
    if(autoMeshUpdates)
        updateMesh();

    oxleytimer.toc("refineBoundary...Done");
}

void Rectangle::refineRegion(double x0, double x1, double y0, double y1)
{
    oxleytimer.toc("refineRegion...");

    z_needs_update=true;
    iz_needs_update=true;

    // If the boundaries were not specified by the user, default to the border of the domain
    forestData.refinement_boundaries[0] = x0 == -1 ? forestData.m_origin[0] : x0; 
    forestData.refinement_boundaries[1] = x1 == -1 ? forestData.m_origin[1] : x1;
    forestData.refinement_boundaries[2] = y0 == -1 ? forestData.m_lxy[0] : y0;
    forestData.refinement_boundaries[3] = y1 == -1 ? forestData.m_lxy[1] : y1;

#ifdef OXLEY_ENABLE_DEBUG_REFINE_REGION
    std::cout << "Rectangle::refineRegion" << std::endl;
    std::cout << "Region boundaries = " << x0 << ", " << y0 << " and " << x1 << ", " << y1 << std::endl;
#endif

    p4est_refine_ext(p4est, true, -1, refine_region, init_rectangle_data, NULL);
    p4est_balance_ext(p4est, P4EST_CONNECT_FULL, init_rectangle_data, NULL);

    // Make sure that nothing went wrong
#ifdef OXLEY_ENABLE_DEBUG
    if(!p4est_is_valid(p4est))
        throw OxleyException("p4est broke during refinement");
    if(!p4est_connectivity_is_valid(connectivity))
        throw OxleyException("connectivity broke during refinement");
#endif

    bool partition_for_coarsening = true;
    p4est_partition_ext(p4est, partition_for_coarsening, NULL);

    // Update the nodes
    p4est_lnodes_destroy(nodes);
    p4est_ghost_t * ghost = p4est_ghost_new(p4est, P4EST_CONNECT_FULL);
    nodes = p4est_lnodes_new(p4est, ghost, 1);
    p4est_ghost_destroy(ghost);

    // Update
    if(autoMeshUpdates)
        updateMesh();

    oxleytimer.toc("refineRegion...Done");
}

void Rectangle::refinePoint(double x0, double y0)
{
    oxleytimer.toc("refinePoint...");

    z_needs_update=true;
    iz_needs_update=true;

    // Check that the point is inside the domain
    if(x0 < forestData.m_origin[0] || x0 > forestData.m_length[0] 
        || y0 < forestData.m_origin[1] || y0 > forestData.m_length[1] )
    {
        throw OxleyException("Coordinates lie outside the domain.");
    }

    // If the boundaries were not specified by the user, default to the border of the domain
    forestData.refinement_boundaries[0] = x0;
    forestData.refinement_boundaries[1] = y0;
    p4est_refine_ext(p4est, true, -1, refine_point, init_rectangle_data, NULL);
    p4est_balance_ext(p4est, P4EST_CONNECT_FULL, init_rectangle_data, NULL);

    // Make sure that nothing went wrong
#ifdef OXLEY_ENABLE_DEBUG
    if(!p4est_is_valid(p4est))
        throw OxleyException("p4est broke during refinement");
    if(!p4est_connectivity_is_valid(connectivity))
        throw OxleyException("connectivity broke during refinement");
#endif

    bool partition_for_coarsening = true;
    p4est_partition_ext(p4est, partition_for_coarsening, NULL);

    // Update
    if(autoMeshUpdates)
        updateMesh();

    oxleytimer.toc("refinePoint...Done");
}

void Rectangle::refineCircle(double x0, double y0, double r)
{
    oxleytimer.toc("refineCircle...");

    z_needs_update=true;
    iz_needs_update=true;
    
    // Check that the point is inside the domain
    if(x0 < forestData.m_origin[0] || x0 > forestData.m_lxy[0] 
        || y0 < forestData.m_origin[1] || y0 > forestData.m_lxy[1] )
    {
        throw OxleyException("Coordinates lie outside the domain.");
    }

    // If the boundaries were not specified by the user, default to the border of the domain
    forestData.refinement_boundaries[0] = x0;
    forestData.refinement_boundaries[1] = y0;
    forestData.refinement_boundaries[2] = r;
    p4est_refine_ext(p4est, true, -1, refine_circle, init_rectangle_data, NULL);
    p4est_balance_ext(p4est, P4EST_CONNECT_FULL, init_rectangle_data, NULL);

    // Make sure that nothing went wrong
#ifdef OXLEY_ENABLE_DEBUG
    if(!p4est_is_valid(p4est))
        throw OxleyException("p4est broke during refinement");
    if(!p4est_connectivity_is_valid(connectivity))
        throw OxleyException("connectivity broke during refinement");
#endif

    bool partition_for_coarsening = true;
    p4est_partition_ext(p4est, partition_for_coarsening, NULL);

    // Update the nodes
    p4est_lnodes_destroy(nodes);
    p4est_ghost_t * ghost = p4est_ghost_new(p4est, P4EST_CONNECT_FULL);
    nodes = p4est_lnodes_new(p4est, ghost, 1);
    p4est_ghost_destroy(ghost);

    // Update
    if(autoMeshUpdates)
        updateMesh();

    oxleytimer.toc("refineCircle...Done");
}
#endif //ESYS_HAVE_TRILINOS

void Rectangle::refineMask(escript::Data mask)
{
    oxleytimer.toc("refineCircle...");

    z_needs_update=true;
    iz_needs_update=true;

    // update the quadrant id information
    updateQuadrantIDinformation();

    // If the boundaries were not specified by the user, default to the border of the domain
    forestData.mask = mask;
    bool refine_recursively = false;
    p4est_refine_ext(p4est, refine_recursively, -1, 
                refine_mask, init_rectangle_data, NULL);
    p4est_balance_ext(p4est, P4EST_CONNECT_FULL, init_rectangle_data, NULL);

    // Make sure that nothing went wrong
#ifdef OXLEY_ENABLE_DEBUG
    if(!p4est_is_valid(p4est))
        throw OxleyException("p4est broke during refinement");
    if(!p4est_connectivity_is_valid(connectivity))
        throw OxleyException("connectivity broke during refinement");
#endif

    bool partition_for_coarsening = true;
    p4est_partition_ext(p4est, partition_for_coarsening, NULL);

    // Update the nodes
    p4est_lnodes_destroy(nodes);
    p4est_ghost_t * ghost = p4est_ghost_new(p4est, P4EST_CONNECT_FULL);
    nodes = p4est_lnodes_new(p4est, ghost, 1);
    p4est_ghost_destroy(ghost);

    // Update
    if(autoMeshUpdates)
        updateMesh();

    oxleytimer.toc("refineCircle...Done");
}

void Rectangle::updateQuadrantIDinformation()
{
    oxleytimer.toc("updateQuadrantIDinformation");
    // quadrant IDs
    // TODO speed up
    // for(p4est_topidx_t treeid = p4est->first_local_tree; treeid <= p4est->last_local_tree; ++treeid) {
    //     p4est_tree_t * tree = p4est_tree_array_index(p4est->trees, treeid);
    //     sc_array_t * tquadrants = &tree->quadrants;
    //     p4est_locidx_t Q = (p4est_locidx_t) tquadrants->elem_count;
    //     for(int q = 0; q < Q; ++q) { 
    //         p4est_quadrant_t * quad = p4est_quadrant_array_index(tquadrants, q);
    //         quadrantData * quadData = (quadrantData *) quad->p.user_data;
    //         double xy[3];
    //         p4est_qcoord_to_vertex(p4est->connectivity, treeid, quad->x, quad->y, xy);
    //         long nodeid = NodeIDs.find(std::make_pair(xy[0],xy[1]))->second;
    //         quadData->nodeid=nodeid;
    //         quadData->treeid=treeid;
    //         quadData->xy[0]=xy[0];
    //         quadData->xy[1]=xy[1];
    //     }
    // }
    oxleytimer.toc("done");
}

escript::Data Rectangle::getX() const
{
    escript::Data out=escript::Vector(0,escript::continuousFunction(*this),true);
    setToX(out);
    out.setProtection();
    return out;
}

void Rectangle::print_debug_report(std::string locat)
{
    std::cout << "report for " <<  locat << std::endl;
    std::cout << "p4est = " << &p4est << std::endl;
    if(!p4est_is_valid(p4est))
        std::cout << "WARNING: p4est is invalid" << std::endl;
    std::cout << "forestData = " << &forestData << std::endl;
    std::cout << "connectivity = " << &connectivity << std::endl;
    if(!p4est_connectivity_is_valid(connectivity))
        std::cout << "WARNING: connectivity is invalid" << std::endl;
    std::cout << "temp_data = " << &temp_data << std::endl;

}

Assembler_ptr Rectangle::createAssembler(std::string type, const DataMap& constants) const
{
    bool isComplex = false;
    DataMap::const_iterator it;
    for(it = constants.begin(); it != constants.end(); it++) {
        if(!it->second.isEmpty() && it->second.isComplex()) {
            isComplex = true;
            break;
        }
    }

    if(type.compare("DefaultAssembler") == 0) {
        if(isComplex) {
            return Assembler_ptr(new DefaultAssembler2D<cplx_t>(shared_from_this()));
        } else {
            return Assembler_ptr(new DefaultAssembler2D<real_t>(shared_from_this()));
        }
    } 
    throw escript::NotImplementedError("oxley::rectangle does not support the requested assembler");
}

// return True for a boundary node and False for an internal node
bool Rectangle::isBoundaryNode(p4est_quadrant_t * quad, int n, p4est_topidx_t treeid, p4est_qcoord_t length) const
{
    double lx = length * ((int) (n % 2) == 1);
    double ly = length * ((int) (n / 2) == 1);
    double xy[3];
    p4est_qcoord_to_vertex(p4est->connectivity, treeid, quad->x+lx, quad->y+ly, xy);
    return (xy[0] == forestData.m_origin[0]) || (xy[0] == forestData.m_lxy[0]) || (xy[1] == forestData.m_origin[1]) || (xy[1] == forestData.m_lxy[1]);
}

// returns True for a boundary node on the north or east of the domain
bool Rectangle::isUpperBoundaryNode(p4est_quadrant_t * quad, int n, p4est_topidx_t treeid, p4est_qcoord_t length) const
{
    double lx = length * ((int) (n % 2) == 1);
    double ly = length * ((int) (n / 2) == 1);
    double xy[3];
    p4est_qcoord_to_vertex(p4est->connectivity, treeid, quad->x+lx, quad->y+ly, xy);
    return (xy[0] == forestData.m_lxy[0]) || (xy[1] == forestData.m_lxy[1]);
}

// returns True for a boundary node on the south or west of the domain
bool Rectangle::isLowerBoundaryNode(p4est_quadrant_t * quad, int n, p4est_topidx_t treeid, p4est_qcoord_t length) const
{
    double lx = length * ((int) (n % 2) == 1);
    double ly = length * ((int) (n / 2) == 1);
    double xy[3];
    p4est_qcoord_to_vertex(p4est->connectivity, treeid, quad->x+lx, quad->y+ly, xy);
    return (xy[0] == forestData.m_origin[0]) || (xy[1] == forestData.m_origin[1]);
}

// returns True for a boundary node on the left boundary
bool Rectangle::isLeftBoundaryNode(p4est_quadrant_t * quad, int n, p4est_topidx_t treeid, p4est_qcoord_t length) const
{
    double lx = length * ((int) (n % 2) == 1);
    double ly = length * ((int) (n / 2) == 1);
    double xy[3];
    p4est_qcoord_to_vertex(p4est->connectivity, treeid, quad->x+lx, quad->y+ly, xy);
    return (xy[0] == forestData.m_origin[0]);
}

// returns True for a boundary node on the right boundary
bool Rectangle::isRightBoundaryNode(p4est_quadrant_t * quad, int n, p4est_topidx_t treeid, p4est_qcoord_t length) const
{
    double lx = length * ((int) (n % 2) == 1);
    double ly = length * ((int) (n / 2) == 1);
    double xy[3];
    p4est_qcoord_to_vertex(p4est->connectivity, treeid, quad->x+lx, quad->y+ly, xy);
    return (xy[0] == forestData.m_lxy[0]);
}

// returns True for a boundary node on the bottom boundary
bool Rectangle::isBottomBoundaryNode(p4est_quadrant_t * quad, int n, p4est_topidx_t treeid, p4est_qcoord_t length) const
{
    double lx = length * ((int) (n % 2) == 1);
    double ly = length * ((int) (n / 2) == 1);
    double xy[3];
    p4est_qcoord_to_vertex(p4est->connectivity, treeid, quad->x+lx, quad->y+ly, xy);
    return (xy[1] == forestData.m_origin[1]);
}

// returns True for a boundary node on the top boundary
bool Rectangle::isTopBoundaryNode(p4est_quadrant_t * quad, int n, p4est_topidx_t treeid, p4est_qcoord_t length) const
{
    double lx = length * ((int) (n % 2) == 1);
    double ly = length * ((int) (n / 2) == 1);
    double xy[3];
    p4est_qcoord_to_vertex(p4est->connectivity, treeid, quad->x+lx, quad->y+ly, xy);
    return (xy[1] == forestData.m_lxy[1]);
}

// return True for a hanging node and False for an non-hanging node
bool Rectangle::isHangingNode(p4est_lnodes_code_t face_code, int n) const
{
    if(face_code == 0)
    {
        return false;
    }
    else
    {
        int8_t c = face_code & 1;
        int8_t ishanging0 = (face_code >> 2) & 1;
        int8_t ishanging1 = (face_code >> 3) & 1;

        int8_t f0 = p4est_corner_faces[c][0];
        int8_t f1 = p4est_corner_faces[c][1];

        // int d0 = c / 2 == 0;
        // int d1 = c % 2 == 0;
        
        return ((ishanging0 == 1) && (p4est_corner_face_corners[c][f0] == 1)) 
            || ((ishanging1 == 1) && (p4est_corner_face_corners[c][f1] == 1));
    }
}

bool Rectangle::getHangingNodes(p4est_lnodes_code_t face_code, int hanging_corner[P4EST_CHILDREN]) const
{
    static const int ones = P4EST_CHILDREN - 1;

#pragma omp parallel for
    for(int i =0;i<P4EST_CHILDREN;i++)
        hanging_corner[i]=-1;

    if (face_code) {
        const int c = (int) (face_code & ones);
        int work = (int) (face_code >> P4EST_DIM);

        /* These two corners are never hanging by construction. */
        hanging_corner[c] = hanging_corner[c ^ ones] = -1;
        for (int i = 0; i < P4EST_DIM; ++i) {
            /* Process face hanging corners. */
           int h = c ^ (1 << i);
            hanging_corner[h ^ ones] = (work & 1) ? c : -1;
            work >>= 1;
        }
        return 1;
    }
    else
    {
        return 0;
    }
}

//protected
void Rectangle::updateNodeIncrements()
{
    nodeIncrements[0] = 1;
    for(p4est_topidx_t treeid = p4est->first_local_tree+1, k=1; treeid <= p4est->last_local_tree; treeid++, k++) 
    {
        p4est_tree_t * tree = p4est_tree_array_index(p4est->trees, treeid);
        sc_array_t * tquadrants = &tree->quadrants;
        p4est_locidx_t Q = (p4est_locidx_t) tquadrants->elem_count;
        nodeIncrements[k] = nodeIncrements[k-1] + Q;
    }
}

// void Rectangle::renumberHangingNodes()
// {
//     hangingNodeIDs.clear();
//     long numNodes = getNumNodes();

// #pragma omp for
//     for(p4est_topidx_t treeid = p4est->first_local_tree; treeid <= p4est->last_local_tree; ++treeid) {
//         p4est_tree_t * tree = p4est_tree_array_index(p4est->trees, treeid);
//         sc_array_t * tquadrants = &tree->quadrants;
//         p4est_locidx_t Q = (p4est_locidx_t) tquadrants->elem_count;
//         // This records which faces have hanging nodes on them
//         bool hanging[4] = {false};

//         for(int q = 0; q < Q; ++q) { 
//             int k = q - Q + nodeIncrements[treeid - p4est->first_local_tree];
   
//             // If there are no hanging nodes here, skip to the next quadrant
//             if(nodes->face_code[k] == 0)
//             {
//                 continue;
//             }

//             // Decode the info
//             int8_t face_code = nodes->face_code[k];
//             int8_t c = face_code & 0x03;
//             int8_t ishanging0 = (face_code >> 2) & 0x01;
//             int8_t ishanging1 = face_code >> 3;

//             int8_t f0 = p4est_corner_faces[c][0];
//             int8_t f1 = p4est_corner_faces[c][1];

//             int d0 = c / 2 == 0;
//             int d1 = c % 2 == 0;

//             int tmp0 = p4est_face_corners[f0][d0];
//             int tmp1 = p4est_face_corners[f1][d1];
            
//             // Record which nodes are hanging
//             if(ishanging0 || ishanging1)
//             {
//                 p4est_quadrant_t * quad = p4est_quadrant_array_index(tquadrants, q);
//                 p4est_qcoord_t length = P4EST_QUADRANT_LEN(quad->level);
//                 double xy[3];

//                 if(ishanging0)
//                 {                   
//                     double lx = length * ((int) (tmp0 % 2) == 1);
//                     double ly = length * ((int) (tmp0 / 2) == 1);
//                     p4est_qcoord_to_vertex(p4est->connectivity, treeid, quad->x+lx, quad->y+ly, xy);
//                     if(!hangingNodeIDs.count(std::make_pair(xy[0],xy[1])))
//                         hangingNodeIDs[std::make_pair(xy[0],xy[1])]=hangingNodeIDs.size()+numNodes;
//                 }
//                 if(ishanging1)
//                 {
//                     double lx = length * ((int) (tmp1 % 2) == 1);
//                     double ly = length * ((int) (tmp1 / 2) == 1);
//                     p4est_qcoord_to_vertex(p4est->connectivity, treeid, quad->x+lx, quad->y+ly, xy);
//                     if(!hangingNodeIDs.count(std::make_pair(xy[0],xy[1])))
//                         hangingNodeIDs[std::make_pair(xy[0],xy[1])]=hangingNodeIDs.size()+numNodes;
//                 }
//             }   
//         }
//     }
// }

void Rectangle::renumberNodes()
{
    oxleytimer.toc("renumberNodes...");

    // The global node numbering now comes directly from p4est_lnodes (built in
    // the constructor / after refinement); this routine only derives m_nodeId,
    // the global id of each local node, from the lnodes owned/ghost partition.
    // The legacy coordinate-hash containers are retired.
    quadrantInfo.clear();
    hanging_face_orientation.clear();

    const long nOwned = (long) nodes->owned_count;
    const long nLocal = (long) nodes->num_local_nodes;
    m_nodeId.clear();
    m_nodeId.resize(nLocal);
    for(long i = 0; i < nOwned; ++i)
        m_nodeId[i] = (long) nodes->global_offset + i;
    for(long i = nOwned; i < nLocal; ++i)
        m_nodeId[i] = (long) nodes->nonlocal_nodes[i - nOwned];
    m_nodeId.shrink_to_fit();

    // Trilinos map inputs: row map = owned global ids; col map = all local
    // global ids (owned first, then ghost -- the lnodes local ordering).
    myColumns.assign(m_nodeId.begin(), m_nodeId.end());
    myRows.assign(m_nodeId.begin(), m_nodeId.begin() + nOwned);

    // MPI: build the ghost element halo and extend myColumns with any 2nd-layer
    // ghost nodes (nodes that appear only on ghost elements). (A6.)
    buildParallelOverlap();

    oxleytimer.toc("renumberNodes...Done");
}

//protected
void Rectangle::buildParallelOverlap()
{
    if (m_ghost) { p4est_ghost_destroy(m_ghost); m_ghost = nullptr; }
    m_ghostElemNodes.clear();

    // Serial: owned == all local nodes, no halo needed.
    if (m_mpiInfo->size <= 1)
        return;

    const int V = nodes->vnodes;                       // 4 corners (degree-1)
    const long nLocal = (long) nodes->num_local_nodes;
    const long nLocalElem = (long) nodes->num_local_elements;

    // global node id -> local (column) index for all lnodes-local nodes
    std::unordered_map<long,long> g2l;
    g2l.reserve((size_t) nLocal * 2);
    for (long i = 0; i < nLocal; ++i)
        g2l[(long) m_nodeId[i]] = i;

    // FULL (face+corner) ghost layer of the current forest -- this is exactly
    // the one-element halo incident to the owned nodes.
    m_ghost = p4est_ghost_new(p4est, P4EST_CONNECT_FULL);
    const long nGhost  = (long) m_ghost->ghosts.elem_count;
    const long nMirror = (long) m_ghost->mirrors.elem_count;

    // For each local element, its V corner GLOBAL node ids (the data we mirror
    // to the ranks that hold that element as a ghost).
    std::vector<p4est_gloidx_t> localElemGN((size_t) nLocalElem * V);
    for (long e = 0; e < nLocalElem; ++e)
        for (int c = 0; c < V; ++c)
            localElemGN[(size_t) e*V + c] =
                (p4est_gloidx_t) m_nodeId[ nodes->element_nodes[(size_t) e*V + c] ];

    // mirror_data[m] -> the V global ids of the local element that is mirror m
    std::vector<void*> mirror_data((size_t) nMirror, nullptr);
    for (long m = 0; m < nMirror; ++m) {
        p4est_quadrant_t* mq = p4est_quadrant_array_index(&m_ghost->mirrors, m);
        const long le = (long) mq->p.piggy3.local_num;   // cumulative local elem id
        mirror_data[m] = (void*) &localElemGN[(size_t) le*V];
    }

    // Receive, per ghost quadrant, its V corner global node ids.
    std::vector<p4est_gloidx_t> ghostElemGN((size_t) nGhost * V);
    p4est_ghost_exchange_custom(p4est, m_ghost,
                                (size_t) V * sizeof(p4est_gloidx_t),
                                mirror_data.data(), ghostElemGN.data());

    // Map ghost element corners to extended local column indices; nodes not
    // already local (2nd layer) get a fresh column index appended to myColumns.
    m_ghostElemNodes.resize((size_t) nGhost * V);
    long nextLocal = nLocal;
    for (long g = 0; g < nGhost; ++g) {
        for (int c = 0; c < V; ++c) {
            const long gid = (long) ghostElemGN[(size_t) g*V + c];
            auto it = g2l.find(gid);
            long lidx;
            if (it != g2l.end()) {
                lidx = it->second;
            } else {
                lidx = nextLocal++;
                g2l[gid] = lidx;
                myColumns.push_back((index_t) gid);
            }
            m_ghostElemNodes[(size_t) g*V + c] = (index_t) lidx;
        }
    }
}

//protected
void Rectangle::assembleCoordinates(escript::Data& arg) const
{
    
    if (!arg.isDataPointShapeEqual(1, &m_numDim))
        throw ValueError("assembleCoordinates: Invalid Data object shape");
    if (!arg.numSamplesEqual(1, getNumNodes()))
        throw ValueError("assembleCoordinates: Illegal number of samples in Data object");
    arg.requireWrite();

#ifdef OXLEY_ENABLE_DEBUG_ASSEMBLE_COORDINATES
    std::cout << "assemble coordinates " << std::endl;
    float bounds[4]={0.0};
#endif

    std::vector<bool> duplicates(getNumNodes(),false);

    const int V = nodes->vnodes;   // 4 corners (degree-1 lnodes)
    long e = 0;                    // running local leaf index (lnodes order)
    for(p4est_topidx_t treeid = p4est->first_local_tree; treeid <= p4est->last_local_tree; ++treeid) {
        p4est_tree_t * tree = p4est_tree_array_index(p4est->trees, treeid);
        sc_array_t * tquadrants = &tree->quadrants;
        p4est_locidx_t Q = (p4est_locidx_t) tquadrants->elem_count;

        for(int q = 0; q < Q; ++q, ++e) { // Loop over the elements attached to the tree
            p4est_quadrant_t * quad = p4est_quadrant_array_index(tquadrants, q);
            p4est_qcoord_t length = P4EST_QUADRANT_LEN(quad->level);

            // A HANGING corner has no node of its own: the slot holds a MASTER,
            // which is a node of the coarse neighbour and lies OUTSIDE this
            // element. Writing this corner's position into it would move the
            // master to the hanging position, so skip those slots. Usually the
            // master gets its coordinate from an element where it is a real
            // corner; under MPI it need not have one here, which the second pass
            // below deals with.
            int hangingCorner[P4EST_CHILDREN];
            const bool anyHanging = getHangingNodes(nodes->face_code[e],
                                                    hangingCorner);

            // Loop over the four corners of the quadrant (z-order matches lnodes)
            for(int n = 0; n < 4; ++n){
                if(anyHanging && hangingCorner[n] >= 0)
                    continue;

                double lx = length * ((int) (n % 2) == 1);
                double ly = length * ((int) (n / 2) == 1);
                double xy[3];
                p4est_qcoord_to_vertex(p4est->connectivity, treeid, quad->x+lx, quad->y+ly, xy);

                // lnodes local node id for this corner (no coordinate hashing)
                long lni = (long) nodes->element_nodes[(size_t) e * V + n];

                if(duplicates[lni] == true)
                    continue;
                else
                    duplicates[lni] = true;

                double * point = arg.getSampleDataRW(lni);
                point[0] = xy[0];
                point[1] = xy[1];
#ifdef OXLEY_ENABLE_DEBUG_ASSEMBLE_COORDINATES
                std::cout<<"lni="<<lni<<"\tCorner=("<<xy[0]<<","<<xy[1]<<"),"<<std::endl;
                if(xy[0] < bounds[0]) bounds[0]=xy[0];
                if(xy[0] > bounds[1]) bounds[1]=xy[0];
                if(xy[1] < bounds[2]) bounds[2]=xy[1];
                if(xy[1] > bounds[3]) bounds[3]=xy[1];
#endif
                // }
            }
        }
    }
    // Fill in any node our elements only ever reference as a hanging slot; it
    // has no coordinate yet and would silently stay at the origin. (MPI only.)
    {
        std::vector<long> fmIds;
        std::vector<double> fmXY;
        farMasterCoords(fmIds, fmXY);
        for (size_t i = 0; i < fmIds.size(); ++i) {
            const long lni = fmIds[i];
            if (lni < 0 || lni >= (long) getNumNodes() || duplicates[lni])
                continue;
            duplicates[lni] = true;
            double * point = arg.getSampleDataRW(lni);
            point[0] = fmXY[2*i];
            point[1] = fmXY[2*i+1];
        }
    }

#ifdef OXLEY_ENABLE_DEBUG_ASSEMBLE_COORDINATES_POINTS
    std::cout << "assembleCoordinates new points are..." << std::endl;
    for(int i = 0; i < getNumNodes() ; i++)
    {
        double * point = arg.getSampleDataRW(i);
        std::cout << i << ": " << point[0] << ", " << point[1] << std::endl;
    }
#endif
#ifdef OXLEY_ENABLE_DEBUG_ASSEMBLE_COORDINATES
    std::cout << "bounds " << bounds[0] << ", " << bounds[1] 
        << " and " << bounds[2] << ", " << bounds[3] << std::endl;
#endif
}

// //private
//void Rectangle::populateDofMap()
//{
//     // const dim_t nDOF0 = getNumDOFInAxis(0);
//     // const dim_t nDOF1 = getNumDOFInAxis(1);
//     // const index_t left = getFirstInDim(0);
//     // const index_t bottom = getFirstInDim(1);

//     // populate node->DOF mapping with own degrees of freedom.
//     // The rest is assigned in the loop further down
//     m_dofMap.assign(getNumNodes(), 0);
// #pragma omp parallel for
//     for(index_t i=bottom; i<bottom+nDOF1; i++) {
//         for(index_t j=left; j<left+nDOF0; j++) {
//             m_dofMap[i*m_NN[0]+j]=(i-bottom)*nDOF0+j-left;
//         }
//     }

//     // build list of shared components and neighbours by looping through
//     // all potential neighbouring ranks and checking if positions are
//     // within bounds
//     const dim_t numDOF=nDOF0*nDOF1;
//     RankVector neighbour;
//     IndexVector offsetInShared(1,0);
//     IndexVector sendShared, recvShared;
//     const int x=m_mpiInfo->rank%m_NX[0];
//     const int y=m_mpiInfo->rank/m_NX[0];
//     // numShared will contain the number of shared DOFs after the following
//     // blocks
//     dim_t numShared=0;
//     // sharing bottom edge
//     if (y > 0) {
//         neighbour.push_back((y-1)*m_NX[0] + x);
//         const dim_t num = nDOF0;
//         offsetInShared.push_back(offsetInShared.back()+num);
//         for(dim_t i=0; i<num; i++, numShared++) {
//             sendShared.push_back(i);
//             recvShared.push_back(numDOF+numShared);
//             m_dofMap[left+i]=numDOF+numShared;
//         }
//     }
//     // sharing top edge
//     if (y < m_NX[1] - 1) {
//         neighbour.push_back((y+1)*m_NX[0] + x);
//         const dim_t num = nDOF0;
//         offsetInShared.push_back(offsetInShared.back()+num);
//         for(dim_t i=0; i<num; i++, numShared++) {
//             sendShared.push_back(numDOF-num+i);
//             recvShared.push_back(numDOF+numShared);
//             m_dofMap[m_NN[0]*(m_NN[1]-1)+left+i]=numDOF+numShared;
//         }
//     }
//     // sharing left edge
//     if (x > 0) {
//         neighbour.push_back(y*m_NX[0] + x-1);
//         const dim_t num = nDOF1;
//         offsetInShared.push_back(offsetInShared.back()+num);
//         for(dim_t i=0; i<num; i++, numShared++) {
//             sendShared.push_back(i*nDOF0);
//             recvShared.push_back(numDOF+numShared);
//             m_dofMap[(bottom+i)*m_NN[0]]=numDOF+numShared;
//         }
//     }
//     // sharing right edge
//     if (x < m_NX[0] - 1) {
//         neighbour.push_back(y*m_NX[0] + x+1);
//         const dim_t num = nDOF1;
//         offsetInShared.push_back(offsetInShared.back()+num);
//         for(dim_t i=0; i<num; i++, numShared++) {
//             sendShared.push_back((i+1)*nDOF0-1);
//             recvShared.push_back(numDOF+numShared);
//             m_dofMap[(bottom+1+i)*m_NN[0]-1]=numDOF+numShared;
//         }
//     }
//     // sharing bottom-left node
//     if (x > 0 && y > 0) {
//         neighbour.push_back((y-1)*m_NX[0] + x-1);
//         // sharing a node
//         offsetInShared.push_back(offsetInShared.back()+1);
//         sendShared.push_back(0);
//         recvShared.push_back(numDOF+numShared);
//         m_dofMap[0]=numDOF+numShared;
//         ++numShared;
//     }
//     // sharing top-left node
//     if (x > 0 && y < m_NX[1]-1) {
//         neighbour.push_back((y+1)*m_NX[0] + x-1);
//         offsetInShared.push_back(offsetInShared.back()+1);
//         sendShared.push_back(numDOF-nDOF0);
//         recvShared.push_back(numDOF+numShared);
//         m_dofMap[m_NN[0]*(m_NN[1]-1)]=numDOF+numShared;
//         ++numShared;
//     }
//     // sharing bottom-right node
//     if (x < m_NX[0]-1 && y > 0) {
//         neighbour.push_back((y-1)*m_NX[0] + x+1);
//         offsetInShared.push_back(offsetInShared.back()+1);
//         sendShared.push_back(nDOF0-1);
//         recvShared.push_back(numDOF+numShared);
//         m_dofMap[m_NN[0]-1]=numDOF+numShared;
//         ++numShared;
//     }
//     // sharing top-right node
//     if (x < m_NX[0]-1 && y < m_NX[1]-1) {
//         neighbour.push_back((y+1)*m_NX[0] + x+1);
//         offsetInShared.push_back(offsetInShared.back()+1);
//         sendShared.push_back(numDOF-1);
//         recvShared.push_back(numDOF+numShared);
//         m_dofMap[m_NN[0]*m_NN[1]-1]=numDOF+numShared;
//         ++numShared;
//     }

// #ifdef ESYS_HAVE_PASO
//     createPasoConnector(neighbour, offsetInShared, offsetInShared, sendShared,
//                         recvShared);
// #endif

    // useful debug output
    /*
    std::cout << "--- rcv_shcomp ---" << std::endl;
    std::cout << "numDOF=" << numDOF << ", numNeighbors=" << neighbour.size() << std::endl;
    for(size_t i=0; i<neighbour.size(); i++) {
        std::cout << "neighbor[" << i << "]=" << neighbour[i]
            << " offsetInShared[" << i+1 << "]=" << offsetInShared[i+1] << std::endl;
    }
    for(size_t i=0; i<recvShared.size(); i++) {
        std::cout << "shared[" << i << "]=" << recvShared[i] << std::endl;
    }
    std::cout << "--- snd_shcomp ---" << std::endl;
    for(size_t i=0; i<sendShared.size(); i++) {
        std::cout << "shared[" << i << "]=" << sendShared[i] << std::endl;
    }
    std::cout << "--- dofMap ---" << std::endl;
    for(size_t i=0; i<m_dofMap.size(); i++) {
        std::cout << "m_dofMap[" << i << "]=" << m_dofMap[i] << std::endl;
    }
    */
//}


//private
template<typename Scalar>
void Rectangle::addToMatrixAndRHS(escript::AbstractSystemMatrix* S, escript::Data& F,
         const std::vector<Scalar>& EM_S, const std::vector<Scalar>& EM_F, 
         bool addS, bool addF, index_t e, index_t t, int nEq, int nComp) const
{    
    IndexVector rowIndex(4);
    p4est_tree_t * currenttree = p4est_tree_array_index(p4est->trees, t);
    // global local-leaf index in lnodes order; quadrants_offset is the cumulative
    // number of local quadrants in the trees before t.
    const long g = (long) currenttree->quadrants_offset + (long) e;
    const int V = nodes->vnodes;   // 4 corners (z-order matches lxy above)
    for(int i = 0; i < 4; i++)
        rowIndex[i] = (index_t) nodes->element_nodes[(size_t) g * V + i];

    if(addF)
    {
        Scalar* F_p = F.getSampleDataRW(0, static_cast<Scalar>(0));
        for(index_t i=0; i<rowIndex.size(); i++) {
            if (rowIndex[i]<getNumDOF()) {
                for(int eq=0; eq<nEq; eq++) {
                    F_p[INDEX2(eq, rowIndex[i], nEq)]+=EM_F[INDEX2(eq,i,nEq)];
                }
            }
        }
    }
    if(addS)
    {
        addToSystemMatrix<Scalar>(S, rowIndex, nEq, EM_S);
    }
}

template<typename Scalar>
void Rectangle::addToMatrixAndRHS(escript::AbstractSystemMatrix* S, escript::Data& F,
         const std::vector<Scalar>& EM_S, const std::vector<Scalar>& EM_F, 
         bool addS, bool addF, borderNodeInfo quad, int nEq, int nComp) const
{
    // the 4 quad corners (lnodes, z-order) were stored when the boundary lists
    // were built (updateFaceOffset), so no coordinate-hash lookup is needed.
    long rowIndex[4] = { quad.neighbours[0], quad.neighbours[1],
                         quad.neighbours[2], quad.neighbours[3] };
    if(addF)
    {
        Scalar* F_p = F.getSampleDataRW(0, static_cast<Scalar>(0));
        for(index_t i=0; i<4; i++) {
            if (rowIndex[i]<getNumDOF()) {
                for(int eq=0; eq<nEq; eq++) {
                    F_p[INDEX2(eq, rowIndex[i], nEq)]+=EM_F[INDEX2(eq,i,nEq)];
                }
            }
        }
    }
    if(addS)
    {
        IndexVector rowInd(4);
        for(int i = 0; i < 4; i++)
            rowInd[i]=rowIndex[i];
        addToSystemMatrix<Scalar>(S, rowInd, nEq, EM_S);
    }
}

template
void Rectangle::addToMatrixAndRHS<real_t>(escript::AbstractSystemMatrix* S, escript::Data& F,
         const std::vector<real_t>& EM_S, const std::vector<real_t>& EM_F, 
         bool addS, bool addF, borderNodeInfo firstNode, int nEq, int nComp) const;

template
void Rectangle::addToMatrixAndRHS<cplx_t>(escript::AbstractSystemMatrix* S, escript::Data& F,
         const std::vector<cplx_t>& EM_S, const std::vector<cplx_t>& EM_F, 
         bool addS, bool addF, borderNodeInfo firstNode, int nEq, int nComp) const;

template
void Rectangle::addToMatrixAndRHS<real_t>(escript::AbstractSystemMatrix* S, escript::Data& F,
         const std::vector<real_t>& EM_S, const std::vector<real_t>& EM_F, 
         bool addS, bool addF, index_t e, index_t t, int nEq, int nComp) const;

template
void Rectangle::addToMatrixAndRHS<cplx_t>(escript::AbstractSystemMatrix* S, escript::Data& F,
         const std::vector<cplx_t>& EM_S, const std::vector<cplx_t>& EM_F,
         bool addS, bool addF, index_t e, index_t t, int nEq, int nComp) const;

//protected
template<typename Scalar>
void Rectangle::addToMatrixAndRHSGhost(escript::AbstractSystemMatrix* S, escript::Data& F,
         const std::vector<Scalar>& EM_S, const std::vector<Scalar>& EM_F,
         bool addS, bool addF, const index_t* rowIndex, int nEq, int nComp) const
{
    // rowIndex are extended-local corner node ids of a ghost (halo) element.
    // Only OWNED rows (< getNumDOF()) are kept; the matrix wrapper likewise
    // drops non-owned rows. Columns may be 2nd-layer ghost nodes (valid colMap).
    if(addF)
    {
        Scalar* F_p = F.getSampleDataRW(0, static_cast<Scalar>(0));
        for(int i=0; i<4; i++) {
            if (rowIndex[i]<getNumDOF()) {
                for(int eq=0; eq<nEq; eq++) {
                    F_p[INDEX2(eq, rowIndex[i], nEq)]+=EM_F[INDEX2(eq,i,nEq)];
                }
            }
        }
    }
    if(addS)
    {
        IndexVector rowInd(rowIndex, rowIndex+4);
        addToSystemMatrix<Scalar>(S, rowInd, nEq, EM_S);
    }
}

template
void Rectangle::addToMatrixAndRHSGhost<real_t>(escript::AbstractSystemMatrix* S, escript::Data& F,
         const std::vector<real_t>& EM_S, const std::vector<real_t>& EM_F,
         bool addS, bool addF, const index_t* rowIndex, int nEq, int nComp) const;
template
void Rectangle::addToMatrixAndRHSGhost<cplx_t>(escript::AbstractSystemMatrix* S, escript::Data& F,
         const std::vector<cplx_t>& EM_S, const std::vector<cplx_t>& EM_F,
         bool addS, bool addF, const index_t* rowIndex, int nEq, int nComp) const;

//protected
template<typename Scalar>
std::vector<Scalar> Rectangle::exchangeGhostCoeff(const escript::Data& coef) const
{
    std::vector<Scalar> out;
    if (!m_ghost || coef.isEmpty())
        return out;

    const long nGhost  = (long) m_ghost->ghosts.elem_count;
    const long nMirror = (long) m_ghost->mirrors.elem_count;
    // In-memory bytes per getSampleDataRO() sample: an expanded Data stores one
    // value per quadrature point, a constant/tagged Data only a single point.
    const size_t sampleSize = (coef.actsExpanded()
                                ? (size_t) coef.getNumDataPointsPerSample() : 1)
                            * (size_t) coef.getDataPointSize();
    if (nGhost == 0 || sampleSize == 0)
        return out;

    const Scalar zero = static_cast<Scalar>(0);
    // mirror_data[m] -> the coefficient sample of the local element that is
    // mirror m (its cumulative local element index is piggy3.local_num).
    std::vector<const void*> mirror_data((size_t) nMirror, nullptr);
    for (long m = 0; m < nMirror; ++m) {
        p4est_quadrant_t* mq = p4est_quadrant_array_index(&m_ghost->mirrors, m);
        const long le = (long) mq->p.piggy3.local_num;
        mirror_data[m] = (const void*) coef.getSampleDataRO(le, zero);
    }
    out.resize((size_t) nGhost * sampleSize);
    p4est_ghost_exchange_custom(p4est, m_ghost, sampleSize * sizeof(Scalar),
                                const_cast<void**>(mirror_data.data()), out.data());
    return out;
}

template std::vector<real_t> Rectangle::exchangeGhostCoeff<real_t>(const escript::Data&) const;
template std::vector<cplx_t> Rectangle::exchangeGhostCoeff<cplx_t>(const escript::Data&) const;

//protected
template<typename Scalar>
std::vector<Scalar> Rectangle::exchangeGhostBoundary(const escript::Data& d,
                        const escript::Data& y, size_t& dSize, size_t& ySize) const
{
    std::vector<Scalar> out;
    // In-memory scalars per getSampleDataRO() sample (expanded: one per boundary
    // quadrature point; constant/tagged: a single point).
    dSize = d.isEmpty()?0:(size_t)(d.actsExpanded()?d.getNumDataPointsPerSample():1)*d.getDataPointSize();
    ySize = y.isEmpty()?0:(size_t)(y.actsExpanded()?y.getNumDataPointsPerSample():1)*y.getDataPointSize();
    if (!m_ghost || (dSize==0 && ySize==0))
        return out;
    const long nGhost  = (long) m_ghost->ghosts.elem_count;
    const long nMirror = (long) m_ghost->mirrors.elem_count;
    if (nGhost == 0)
        return out;

    const Scalar zero = static_cast<Scalar>(0);
    const size_t perSide = 1 + dSize + ySize;
    const size_t perOct  = 4 * perSide;

    // Map octant -> per-side FaceElements sample index, keyed by LOCAL LEAF INDEX
    // (mirrors carry a reliable leaf index in piggy3.local_num; piggy3.which_tree
    // is not dependable for mirrors). Build (treeid,quad->x,quad->y) -> leaf over
    // the local leaves, then re-key the boundary faces by leaf index.
    struct Key { p4est_topidx_t t; p4est_qcoord_t x, y;
                 bool operator==(const Key& o) const { return t==o.t && x==o.x && y==o.y; } };
    struct KeyHash { size_t operator()(const Key& k) const {
        return ((size_t)k.t*73856093u) ^ ((size_t)k.x*19349663u) ^ ((size_t)k.y*83492791u); } };
    std::unordered_map<Key, long, KeyHash> octKey2leaf;
    long leaf = 0;
    for (p4est_topidx_t tt = p4est->first_local_tree; tt <= p4est->last_local_tree; ++tt) {
        p4est_tree_t* tree = p4est_tree_array_index(p4est->trees, tt);
        sc_array_t* quads = &tree->quadrants;
        const long Q = (long) quads->elem_count;
        for (long q = 0; q < Q; ++q, ++leaf) {
            p4est_quadrant_t* qd = p4est_quadrant_array_index(quads, q);
            octKey2leaf[Key{ tt, qd->x, qd->y }] = leaf;
        }
    }

    std::unordered_map<long, std::array<long,4>> fmap;   // leaf -> per-side sample
    const std::vector<borderNodeInfo>* lists[4] =
        { &NodeIDsLeft, &NodeIDsRight, &NodeIDsBottom, &NodeIDsTop };
    for (int s = 0; s < 4; ++s) {
        if (m_faceOffset[s] < 0) continue;
        const std::vector<borderNodeInfo>& Lst = *lists[s];
        for (long k = 0; k < (long) Lst.size(); ++k) {
            auto lit = octKey2leaf.find(Key{ Lst[k].treeid, Lst[k].x, Lst[k].y });
            if (lit == octKey2leaf.end()) continue;
            const long lf = lit->second;
            auto it = fmap.find(lf);
            if (it == fmap.end())
                it = fmap.emplace(lf, std::array<long,4>{{-1,-1,-1,-1}}).first;
            it->second[s] = (long) m_faceOffset[s] + k;
        }
    }

    // Pack each mirror octant's boundary d/y samples by side (keyed by leaf index).
    std::vector<Scalar> mirrorPacked((size_t) nMirror * perOct, zero);
    std::vector<void*> mirror_data((size_t) nMirror, nullptr);
    for (long m = 0; m < nMirror; ++m) {
        p4est_quadrant_t* mq = p4est_quadrant_array_index(&m_ghost->mirrors, m);
        Scalar* base = &mirrorPacked[(size_t) m * perOct];
        mirror_data[m] = (void*) base;
        auto it = fmap.find((long) mq->p.piggy3.local_num);
        if (it == fmap.end()) continue;
        for (int s = 0; s < 4; ++s) {
            const long sample = it->second[s];
            if (sample < 0) continue;
            Scalar* sb = base + (size_t) s * perSide;
            sb[0] = static_cast<Scalar>(1);
            if (dSize) { const Scalar* dp = d.getSampleDataRO(sample, zero);
                         std::copy(dp, dp+dSize, sb+1); }
            if (ySize) { const Scalar* yp = y.getSampleDataRO(sample, zero);
                         std::copy(yp, yp+ySize, sb+1+dSize); }
        }
    }

    out.resize((size_t) nGhost * perOct, zero);
    p4est_ghost_exchange_custom(p4est, m_ghost, perOct * sizeof(Scalar),
                                mirror_data.data(), out.data());
    return out;
}

template std::vector<real_t> Rectangle::exchangeGhostBoundary<real_t>(
        const escript::Data&, const escript::Data&, size_t&, size_t&) const;
template std::vector<cplx_t> Rectangle::exchangeGhostBoundary<cplx_t>(
        const escript::Data&, const escript::Data&, size_t&, size_t&) const;

//protected
void Rectangle::interpolateNodesOnElements(escript::Data& out,
                                           const escript::Data& in,
                                           bool reduced) const
{
    if (in.isComplex()!=out.isComplex())
    {
        throw OxleyException("Programmer Error: in and out parameters do not have the same complexity.");
    }
    if (in.isComplex())
    {
        interpolateNodesOnElementsWorker(out, in, reduced, escript::DataTypes::cplx_t(0));
    }
    else
    {
        interpolateNodesOnElementsWorker(out, in, reduced, escript::DataTypes::real_t(0));      
    }
}

//protected
void Rectangle::interpolateNodesOnFaces(escript::Data& out,
                                           const escript::Data& in,
                                           bool reduced) const
{
    if (in.isComplex()!=out.isComplex())
    {
        throw OxleyException("Programmer Error: in and out parameters do not have the same complexity.");
    }
    if (in.isComplex())
    {
        interpolateNodesOnFacesWorker(out, in, reduced, escript::DataTypes::cplx_t(0));
    }
    else
    {
        interpolateNodesOnFacesWorker(out, in, reduced, escript::DataTypes::real_t(0));      
    }
}


namespace {

/**
   \brief
   Replaces the values sitting in an element's hanging corner slots by the
   values AT those corners.

   lnodes does not number the position at a 2:1 seam, so the slot of a hanging
   corner holds a MASTER instead: the FAR one, the node at the other end of the
   coarse neighbour's edge, a distance 2h away. (That is the same fact
   getGhostRealCornerCoords uses to place it at 2H-A.) The value at the hanging
   position is the mean of the two masters, and BOTH of them are corners of this
   very element - the far master in the hanging slot itself, the near one in the
   slot getHangingNodes reports - so the constraint is local to the element and
   needs no halo, no communication and no node that does not already exist.

   Correcting in place is safe: the near master is p4est's anchor corner, which
   is never itself hanging, so no corrected value is ever read back.

   Without this an element on the fine side of a seam reads a value from 2h away
   as though it sat on its own edge. Nothing crashes and nothing looks wrong -
   grad() is simply wrong by O(1) along every seam, and integrals of anything
   interpolated from nodes are slightly wrong everywhere near one.
*/
template<typename Scalar>
inline void constrainHangingCorners(const int hangingCorner[P4EST_CHILDREN],
                                    Scalar* corner[P4EST_CHILDREN],
                                    dim_t numComp)
{
    for (int n = 0; n < P4EST_CHILDREN; ++n) {
        const int a = hangingCorner[n];
        if (a < 0)
            continue;
        for (dim_t i = 0; i < numComp; ++i)
            corner[n][i] = static_cast<Scalar>(0.5)
                         * (corner[n][i] + corner[a][i]);
    }
}

} // anonymous namespace

// private   
template <typename S> 
void Rectangle::interpolateNodesOnElementsWorker(escript::Data& out,
                                           const escript::Data& in,
                                           bool reduced, S sentinel) const
{
    const dim_t numComp = in.getDataPointSize();
    
    if (reduced) {
        out.requireWrite();
        const S c0 = 0.25;
        std::vector<S> f_00(numComp);
        std::vector<S> f_01(numComp);
        std::vector<S> f_10(numComp);
        std::vector<S> f_11(numComp);

        const int V = nodes->vnodes;
        long e = 0;
        for(p4est_topidx_t treeid = p4est->first_local_tree; treeid <= p4est->last_local_tree; ++treeid) {
            p4est_tree_t * tree = p4est_tree_array_index(p4est->trees, treeid);
            sc_array_t * tquadrants = &tree->quadrants;
            p4est_locidx_t Q = (p4est_locidx_t) tquadrants->elem_count;
            for(int q = 0; q < Q; q++, ++e)
            {
                long ids[4];
                for(int n = 0; n < V; ++n) ids[n] = (long) nodes->element_nodes[(size_t) e * V + n];
                const long quadID = e;

                memcpy(&f_00[0], in.getSampleDataRO(ids[0],sentinel), numComp*sizeof(S));
                memcpy(&f_01[0], in.getSampleDataRO(ids[2],sentinel), numComp*sizeof(S));
                memcpy(&f_10[0], in.getSampleDataRO(ids[1],sentinel), numComp*sizeof(S));
                memcpy(&f_11[0], in.getSampleDataRO(ids[3],sentinel), numComp*sizeof(S));

                // on the fine side of a seam these slots hold masters, not the
                // corner values; see constrainHangingCorners
                int hangingCorner[P4EST_CHILDREN];
                if (getHangingNodes(nodes->face_code[e], hangingCorner)) {
                    S* corner[P4EST_CHILDREN] =
                            { &f_00[0], &f_10[0], &f_01[0], &f_11[0] };
                    constrainHangingCorners(hangingCorner, corner, numComp);
                }

                S* o = out.getSampleDataRW(quadID,sentinel);
                for (index_t i=0; i < numComp; ++i) {
                    o[INDEX2(i,numComp,0)] = c0*(f_00[i] + f_01[i] + f_10[i] + f_11[i]);
                }
            }
        }
    } 
    else 
    {
        out.requireWrite();
        const S c0 = 0.16666666666666666667;
        const S c1 = 0.044658198738520451079;
        const S c2 = 0.62200846792814621559;

        std::vector<S> f_00(numComp);
        std::vector<S> f_01(numComp);
        std::vector<S> f_10(numComp);
        std::vector<S> f_11(numComp);

        const int V = nodes->vnodes;
        long e = 0;
        for(p4est_topidx_t treeid = p4est->first_local_tree; treeid <= p4est->last_local_tree; ++treeid) {
            p4est_tree_t * tree = p4est_tree_array_index(p4est->trees, treeid);
            sc_array_t * tquadrants = &tree->quadrants;
            p4est_locidx_t Q = (p4est_locidx_t) tquadrants->elem_count;
            for(int q = 0; q < Q; q++, ++e)
            {
                long ids[4];
                for(int n = 0; n < V; ++n) ids[n] = (long) nodes->element_nodes[(size_t) e * V + n];
                const long quadId = e;

            #ifdef OXLEY_ENABLE_DEBUG_INTERPOLATE_EXTRA
                std::cout << "interpolateNodesOnElementsWorker quadID: " << quadId << ", node IDs " << 
                                    ids[0] << ", " << ids[2] << ", " << ids[1] << ", " << ids[3] << std::endl;
            #endif

                memcpy(&f_00[0], in.getSampleDataRO(ids[0], sentinel), numComp*sizeof(S));
                memcpy(&f_01[0], in.getSampleDataRO(ids[2], sentinel), numComp*sizeof(S));
                memcpy(&f_10[0], in.getSampleDataRO(ids[1], sentinel), numComp*sizeof(S));
                memcpy(&f_11[0], in.getSampleDataRO(ids[3], sentinel), numComp*sizeof(S));

                // on the fine side of a seam these slots hold masters, not the
                // corner values; see constrainHangingCorners
                int hangingCorner[P4EST_CHILDREN];
                if (getHangingNodes(nodes->face_code[e], hangingCorner)) {
                    S* corner[P4EST_CHILDREN] =
                            { &f_00[0], &f_10[0], &f_01[0], &f_11[0] };
                    constrainHangingCorners(hangingCorner, corner, numComp);
                }

                S* o = out.getSampleDataRW(quadId, sentinel);
                for (index_t i=0; i < numComp; ++i) {
                    o[INDEX2(i,numComp,0)] = c0*(f_01[i] + f_10[i]) + c1*f_11[i] + c2*f_00[i];
                    o[INDEX2(i,numComp,1)] = c0*(f_00[i] + f_11[i]) + c1*f_01[i] + c2*f_10[i];
                    o[INDEX2(i,numComp,2)] = c0*(f_00[i] + f_11[i]) + c1*f_10[i] + c2*f_01[i];
                    o[INDEX2(i,numComp,3)] = c0*(f_01[i] + f_10[i]) + c1*f_00[i] + c2*f_11[i];

            #ifdef OXLEY_ENABLE_DEBUG_INTERPOLATE
                std::cout << quadId << ": ";
                std::cout << c0*(f_01[i] + f_10[i]) + c1*f_11[i] + c2*f_00[i] << ", ";
                std::cout << c0*(f_00[i] + f_11[i]) + c1*f_01[i] + c2*f_10[i] << ", ";
                std::cout << c0*(f_00[i] + f_11[i]) + c1*f_10[i] + c2*f_01[i] << ", ";
                std::cout << c0*(f_01[i] + f_10[i]) + c1*f_00[i] + c2*f_11[i] << std::endl;
            #endif
                }
            }
        }
    }
}

//

void Rectangle::p4est_qcoord_to_vertex_mod (p4est_connectivity_t * connectivity,
                        p4est_topidx_t treeid,
                        p4est_qcoord_t x, p4est_qcoord_t y,
                        double vxyz[3]) const
{
    const double       *vertices = connectivity->vertices;
    double wx[2], wy[2];
    const p4est_topidx_t * vindices = connectivity->tree_to_vertex + P4EST_CHILDREN * treeid;

    vxyz[0] = vxyz[1] = vxyz[2] = 0.;

    wx[1] = (double) x / (double) P4EST_ROOT_LEN;
    wx[0] = 1. - wx[1];

    wy[1] = (double) y / (double) P4EST_ROOT_LEN;
    wy[0] = 1. - wy[1];

    for (int yi = 0; yi < 2; ++yi) {
        double yfactor = wy[yi];
        for (int xi = 0; xi < 2; ++xi) {
            double xfactor = yfactor * wx[xi];
            p4est_topidx_t vindex = *vindices++;
            vxyz[0] += xfactor * vertices[3 * vindex + 0];
            vxyz[1] += xfactor * vertices[3 * vindex + 1];
            vxyz[2] += xfactor * vertices[3 * vindex + 2];
        }
    }
}

bool Rectangle::checkHangingBorderNode(p4est_quadrant_t * quad, p4est_qcoord_t x, p4est_qcoord_t y, 
                                p4est_topidx_t treeid, int n) const
{
    p4est_qcoord_t l = P4EST_QUADRANT_LEN(quad->level);
    double xy[3], xyA[3], xyB[3];
    p4est_qcoord_to_vertex(p4est->connectivity, treeid, x, y, xy);

    bool west  = xy[0] == forestData.m_origin[0];
    bool south = xy[1] == forestData.m_origin[1];
    bool east  = xy[0] == forestData.m_lxy[0];
    bool north = xy[1] == forestData.m_lxy[1];
    int position = p4est_quadrant_child_id(quad);

    if(west) // W boundary
    {
        if(position == 2 && n == 0)
            return true;
        else
            return false;
    }
    else if(east) // E boundary
    {
        if(position == 3 && n == 1)
            return true;
        else
            return false;
    }
    else if(south) // S boundary
    {
        if(position == 1 && n == 0)
            return true;
        else
            return false;
    }
    else if(north) // N boundary
    {
        if(position == 2 && n == 3)
            return true;
        else
            return false;
    }
    else // interior 
    {
        return false;
    }
}

// Note: This function assumes that the node is hanging and is on the border
int Rectangle::getHangingBorderNodeFacecode(p4est_quadrant_t * quad, int8_t level, p4est_qcoord_t x, p4est_qcoord_t y,
                                                p4est_topidx_t treeid, int n , int boundary) const
{
    p4est_qcoord_t l = P4EST_QUADRANT_LEN(level);
    int adj[4][2]={{0,0},{l,0},{0,l},{l,l}};
    double child_xy[4][3] = {{-1}};
    for(int i = 0; i < 4; i++)
        p4est_qcoord_to_vertex(p4est->connectivity, treeid, x+adj[i][0], y+adj[i][1], child_xy[i]);

    bool west  = child_xy[0][0] == forestData.m_origin[0];
    bool south = child_xy[0][1] == forestData.m_origin[1];
    bool east  = child_xy[3][0] == forestData.m_lxy[0];
    bool north = child_xy[3][1] == forestData.m_lxy[1];
    int count = -1;

    // get the parent quadrant
    p4est_quadrant_t * parent;
    p4est_quadrant_t parent_quad;
    parent = &parent_quad;
    p4est_quadrant_parent(quad, parent);
    ESYS_ASSERT(p4est_quadrant_is_valid(parent),"getHangingBorderNodeFacecode: Invalid parent quadrant");

    // get the parent quadrants coordinates
    p4est_qcoord_t l2 = P4EST_QUADRANT_LEN(parent->level);
    int adj2[4][2]={{0,0},{l2,0},{0,l2},{l2,l2}};
    double parent_xy[4][3] = {{-1}};
    for(int i = 0; i < 4; i++)
        p4est_qcoord_to_vertex(p4est->connectivity, treeid, parent->x+adj2[i][0], parent->y+adj2[i][1], parent_xy[i]);

    switch(boundary)
    {
    case 3:
        if((child_xy[2][0] == parent_xy[2][0]) && (child_xy[2][1] == parent_xy[2][1]))
            return 3;
        else if((child_xy[3][0] == parent_xy[3][0]) && (child_xy[3][1] == parent_xy[3][1]))
            return 2;
        else
            throw OxleyException("getHangingBorderNodeFacecode: Unknown error.");
        break;
    case 2:
        if((child_xy[0][0] == parent_xy[0][0]) && (child_xy[0][1] == parent_xy[0][1]))
            return 2;
        else if((child_xy[1][0] == parent_xy[1][0]) && (child_xy[1][1] == parent_xy[1][1]))
            return 3;
        else
            throw OxleyException("getHangingBorderNodeFacecode: Unknown error.");
        break;
    case 1:
        if((child_xy[1][0] == parent_xy[1][0]) && (child_xy[1][1] == parent_xy[1][1]))
            return 2;
        else if((child_xy[3][0] == parent_xy[3][0]) && (child_xy[3][1] == parent_xy[3][1]))
            return 3;
        else
            throw OxleyException("getHangingBorderNodeFacecode: Unknown error.");
        break;
    case 0:
        if((child_xy[0][0] == parent_xy[0][0]) && (child_xy[0][1] == parent_xy[0][1]))
            return 2;
        else if((child_xy[2][0] == parent_xy[2][0]) && (child_xy[2][1] == parent_xy[2][1]))
            return 3;
        else
            throw OxleyException("getHangingBorderNodeFacecode: Unknown error.");
        break;
    default:
        throw OxleyException("getHangingBorderNodeFacecode: Node in interior of mesh.");
    }
}


//private
template <typename S>
void Rectangle::gatherCornersConstrained(const escript::Data& in,
                                         const borderNodeInfo& b, dim_t numComp,
                                         S sentinel, std::vector<S>& f_00,
                                         std::vector<S>& f_10,
                                         std::vector<S>& f_01,
                                         std::vector<S>& f_11) const
{
    memcpy(&f_00[0], in.getSampleDataRO(b.neighbours[0], sentinel), numComp*sizeof(S));
    memcpy(&f_10[0], in.getSampleDataRO(b.neighbours[1], sentinel), numComp*sizeof(S));
    memcpy(&f_01[0], in.getSampleDataRO(b.neighbours[2], sentinel), numComp*sizeof(S));
    memcpy(&f_11[0], in.getSampleDataRO(b.neighbours[3], sentinel), numComp*sizeof(S));

    // A hanging corner never lies ON the domain boundary: it is the midpoint of
    // a face shared by two elements, and the relative interior of a shared face
    // is interior to the domain. So neither value on THIS face is ever hanging.
    // The face gradient, however, combines them with the element's other two
    // corners to get the tangential derivative, and one of those can be - which
    // is why all four are read here and corrected together. A routine that uses
    // only the on-face pair is unaffected either way.
    int hangingCorner[P4EST_CHILDREN];
    if (b.quadIndex >= 0
            && getHangingNodes(nodes->face_code[b.quadIndex], hangingCorner)) {
        S* corner[P4EST_CHILDREN] = { &f_00[0], &f_10[0], &f_01[0], &f_11[0] };
        constrainHangingCorners(hangingCorner, corner, numComp);
    }
}

//private
template <typename S>
void Rectangle::interpolateNodesOnFacesWorker(escript::Data& out,
                                        const escript::Data& in,
                                        bool reduced, S sentinel) const
{
    const dim_t numComp = in.getDataPointSize();

    if (reduced) {
        out.requireWrite();

        std::vector<S> f_00(numComp);
        std::vector<S> f_01(numComp);
        std::vector<S> f_10(numComp);
        std::vector<S> f_11(numComp);

        if (m_faceOffset[0] > -1) {
#pragma omp for nowait
            for (index_t k=0; k<NodeIDsLeft.size(); k++) {
                borderNodeInfo tmp = NodeIDsLeft[k];
                gatherCornersConstrained(in, tmp, numComp, sentinel,
                                         f_00, f_10, f_01, f_11);
                
                S* o = out.getSampleDataRW(m_faceOffset[0]+k, sentinel);
                for (index_t i=0; i < numComp; ++i) {
                    o[INDEX2(i,numComp,0)] = (f_00[i] + f_01[i])/static_cast<S>(2);
                } /* end of component loop i */
            } 
        } /* end of face 0 */
        if (m_faceOffset[1] > -1) {
#pragma omp for nowait
            for (index_t k=0; k<NodeIDsRight.size(); k++) {
                borderNodeInfo tmp = NodeIDsRight[k];
                gatherCornersConstrained(in, tmp, numComp, sentinel,
                                         f_00, f_10, f_01, f_11);
                S* o = out.getSampleDataRW(m_faceOffset[1]+k, sentinel);
                for (index_t i=0; i < numComp; ++i) {
                    o[INDEX2(i,numComp,0)] = (f_10[i] + f_11[i])/static_cast<S>(2);
                } /* end of component loop i */
            } 
        } /* end of face 1 */
        if (m_faceOffset[2] > -1) {
#pragma omp for nowait
            for (index_t k=0; k<NodeIDsBottom.size(); k++) {
                borderNodeInfo tmp = NodeIDsBottom[k];
                gatherCornersConstrained(in, tmp, numComp, sentinel,
                                         f_00, f_10, f_01, f_11);
                S* o = out.getSampleDataRW(m_faceOffset[2]+k, sentinel);
                for (index_t i=0; i < numComp; ++i) {
                    o[INDEX2(i,numComp,0)] = (f_00[i] + f_10[i])/static_cast<S>(2);
                } /* end of component loop i */
            } 
        } /* end of face 2 */
        if (m_faceOffset[3] > -1) {
#pragma omp for nowait
            for (index_t k=0; k<NodeIDsTop.size(); k++) {
                borderNodeInfo tmp = NodeIDsTop[k];
                gatherCornersConstrained(in, tmp, numComp, sentinel,
                                         f_00, f_10, f_01, f_11);
                S* o = out.getSampleDataRW(m_faceOffset[3]+k, sentinel);
                for (index_t i=0; i < numComp; ++i) {
                    o[INDEX2(i,numComp,0)] = (f_01[i] + f_11[i])/static_cast<S>(2);
                } /* end of component loop i */
            } 
        } /* end of face 3 */
    } else {
        out.requireWrite();
        const S c0 = 0.21132486540518711775;
        const S c1 = 0.78867513459481288225;

        std::vector<S> f_00(numComp);
        std::vector<S> f_01(numComp);
        std::vector<S> f_10(numComp);
        std::vector<S> f_11(numComp);
        if (m_faceOffset[0] > -1) {
    #pragma omp for nowait
            for (index_t k=0; k<NodeIDsLeft.size(); k++) {
                borderNodeInfo tmp = NodeIDsLeft[k];
                gatherCornersConstrained(in, tmp, numComp, sentinel,
                                         f_00, f_10, f_01, f_11);
                S* o = out.getSampleDataRW(m_faceOffset[0]+k, sentinel);
                for (index_t i=0; i < numComp; ++i) {
                    o[INDEX2(i,numComp,0)] = c0*f_01[i] + c1*f_00[i];
                    o[INDEX2(i,numComp,1)] = c0*f_00[i] + c1*f_01[i];
                } /* end of component loop i */
            }
        } /* end of face 0 */
        if (m_faceOffset[1] > -1) {
    #pragma omp for nowait
            for (index_t k=0; k<NodeIDsRight.size(); k++) {
                borderNodeInfo tmp = NodeIDsRight[k];
                gatherCornersConstrained(in, tmp, numComp, sentinel,
                                         f_00, f_10, f_01, f_11);
                S* o = out.getSampleDataRW(m_faceOffset[1]+k, sentinel);
                for (index_t i=0; i < numComp; ++i) {
                    o[INDEX2(i,numComp,0)] = c1*f_10[i] + c0*f_11[i];
                    o[INDEX2(i,numComp,1)] = c1*f_11[i] + c0*f_10[i];
                } /* end of component loop i */
            } 
        } /* end of face 1 */
        if (m_faceOffset[2] > -1) {
    #pragma omp for nowait
             for (index_t k=0; k<NodeIDsBottom.size(); k++) {
                borderNodeInfo tmp = NodeIDsBottom[k];
                gatherCornersConstrained(in, tmp, numComp, sentinel,
                                         f_00, f_10, f_01, f_11);
                S* o = out.getSampleDataRW(m_faceOffset[2]+k, sentinel);
                for (index_t i=0; i < numComp; ++i) {
                    o[INDEX2(i,numComp,0)] = c0*f_10[i] + c1*f_00[i];
                    o[INDEX2(i,numComp,1)] = c0*f_00[i] + c1*f_10[i];
                } /* end of component loop i */
            } 
        } /* end of face 2 */
        if (m_faceOffset[3] > -1) {
    #pragma omp for nowait
            for (index_t k=0; k<NodeIDsTop.size(); k++) {
                borderNodeInfo tmp = NodeIDsTop[k];
                gatherCornersConstrained(in, tmp, numComp, sentinel,
                                         f_00, f_10, f_01, f_11);
                S* o = out.getSampleDataRW(m_faceOffset[3]+k, sentinel);
                for (index_t i=0; i < numComp; ++i) {
                    o[INDEX2(i,numComp,0)] = c0*f_11[i] + c1*f_01[i];
                    o[INDEX2(i,numComp,1)] = c0*f_01[i] + c1*f_11[i];
                } /* end of component loop i */
            } 
        } /* end of face 3 */
    }
}

////////////////////////////// inline methods ////////////////////////////////
inline dim_t Rectangle::getDofOfNode(dim_t node) const
{
    // Conforming lnodes numbering: every node is a real DOF and (serially)
    // the DOF id equals the local node id. MPI ownership handled in A6.
    return node;
}

// //protected
// inline dim_t Rectangle::getNumDOFInAxis(unsigned axis) const
// {
//     ESYS_ASSERT(axis < m_numDim, "Invalid axis");
//     return (m_gNE[axis]+1)/m_NX[axis];
//     // return 0; //todo
// }

// protected
// inline index_t Rectangle::getFirstInDim(unsigned axis) const
// {
//     // return m_offset[axis] == 0 ? 0 : 1;
//     return nodes->global_offset == 0 ? 0 : 1;
// }

//protected
inline dim_t Rectangle::getNumNodes() const
{
    // lnodes-based node count (owned + ghost). Replaces the coordinate-hash
    // NodeIDs.size(); for a conforming mesh the two counts agree.
    return nodes ? (dim_t) nodes->num_local_nodes : 0;
}

inline dim_t Rectangle::getNumHangingNodes() const
{
    return num_hanging;
}

//protected
inline dim_t Rectangle::getNumElements() const
{
    long numElements = 0;
    for(p4est_topidx_t treeid = p4est->first_local_tree; treeid <= p4est->last_local_tree; ++treeid) {
        p4est_tree_t * tree = p4est_tree_array_index(p4est->trees, treeid);
        sc_array_t * tquadrants = &tree->quadrants;
        p4est_locidx_t Q = (p4est_locidx_t) tquadrants->elem_count;
        numElements+=Q;
    }
    return numElements;
}

bool Rectangle::isConforming() const
{
    int localHang = 0;
    if (nodes) {
        for (long e = 0; e < (long) nodes->num_local_elements; ++e) {
            if (nodes->face_code[e] != 0) { localHang = 1; break; }
        }
    }
    int anyHang = localHang;
#ifdef ESYS_MPI
    MPI_Allreduce(&localHang, &anyHang, 1, MPI_INT, MPI_MAX, m_mpiInfo->comm);
#endif
    return anyHang == 0;
}

//protected
void Rectangle::farMasterCoords(std::vector<long>& ids,
                                std::vector<double>& xy) const
{
    ids.clear();
    xy.clear();

    // A node referenced by our elements ONLY through hanging slots never gets a
    // coordinate from the ordinary walk, which skips those slots so they cannot
    // overwrite their master's position. That happens under MPI: the far master
    // is a corner of the coarse neighbour and of the adjacent fine element, and
    // the SFC cut can put both on another rank. The node is then left at the
    // origin, and everything read through it - getX, any field built from it,
    // and the average that defines the hanging node - is silently wrong.
    //
    // No halo and no communication are needed to repair it. A hanging corner H
    // is the MIDPOINT of the coarse neighbour's edge; the other end A of that
    // edge is a real corner of this same element (it is the slot getHangingNodes
    // reports), so the far master sits at 2H - A. Both H and A are local.
    const int V = nodes->vnodes;
    long e = 0;
    for (p4est_topidx_t treeid = p4est->first_local_tree;
         treeid <= p4est->last_local_tree; ++treeid) {
        p4est_tree_t* tree = p4est_tree_array_index(p4est->trees, treeid);
        sc_array_t* quads = &tree->quadrants;
        const p4est_locidx_t Q = (p4est_locidx_t) quads->elem_count;
        for (p4est_locidx_t q = 0; q < Q; ++q, ++e) {
            int hangingCorner[P4EST_CHILDREN];
            if (!getHangingNodes(nodes->face_code[e], hangingCorner))
                continue;
            p4est_quadrant_t* quad = p4est_quadrant_array_index(quads, q);
            const p4est_qcoord_t len = P4EST_QUADRANT_LEN(quad->level);
            for (int c = 0; c < V; ++c) {
                const int a = hangingCorner[c];
                if (a < 0)
                    continue;
                double H[3] = {0.,0.,0.}, A[3] = {0.,0.,0.};
                p4est_qcoord_to_vertex(p4est->connectivity, treeid,
                        quad->x + (c & 1)*len, quad->y + ((c >> 1) & 1)*len, H);
                p4est_qcoord_to_vertex(p4est->connectivity, treeid,
                        quad->x + (a & 1)*len, quad->y + ((a >> 1) & 1)*len, A);
                ids.push_back((long) nodes->element_nodes[(size_t) e * V + c]);
                xy.push_back(2.*H[0] - A[0]);
                xy.push_back(2.*H[1] - A[1]);
            }
        }
    }
}

namespace {

/// the two corners of face f, in z-order corner indexing
const int rectFaceCorners[4][2] = { {0,2}, {1,3}, {0,1}, {2,3} };

/// One 2:1 seam, described from the COARSE side. The hanging node sits at the
/// midpoint of that octant's face, and (globalQuad, face) is its identity: a
/// hanging position is the midpoint of exactly one coarse face, and every rank
/// touching the seam derives the same pair, so no agreement protocol is needed.
struct Seam
{
    long   globalQuad;      ///< Q: coarse octant's number in the global z-order
    int    face;            ///< f: which of its faces, p4est numbering
    int    owner;           ///< rank owning the coarse octant, hence the node
    double mid[3];          ///< where the node goes
    long   coarseElem;      ///< local element index of the coarse octant, or -1
    long   fineElem[2];     ///< local element indices of the fine octants, or -1
    int    fineCorner[2];   ///< z-order corner of each fine octant sitting at mid
    int    fineOwner[2];    ///< rank owning each fine octant, or -1
};

struct SeamCtx
{
    const p4est_ghost_t* ghost;
    std::vector<Seam>* out;
};

/// owner rank of ghost quadrant g, from the per-rank ghost offsets
inline int ghostOwner(const p4est_ghost_t* ghost, p4est_locidx_t g, int size)
{
    for (int r = 0; r < size; ++r) {
        if (g >= ghost->proc_offsets[r] && g < ghost->proc_offsets[r+1])
            return r;
    }
    return -1;
}

/// cumulative local index of a local quadrant, matching the order the element
/// walk numbers elements in
inline long localElemIndex(p4est_t* p4est, p4est_topidx_t treeid,
                           p4est_locidx_t quadid)
{
    p4est_tree_t* tree = p4est_tree_array_index(p4est->trees, treeid);
    return (long) tree->quadrants_offset + quadid;
}

/**
   p4est_iterate face callback: records every 2:1 seam.

   p4est_iterate fires a face callback on any face shared by two quadrants,
   hanging or not, and skips only faces whose quadrants are all ghosts. So a
   seam is reported on EVERY rank holding either side of it - which is the
   point, since the coarse side has no local mark of its own (face_code lives
   on the fine side).
*/
void collectSeam(p4est_iter_face_info_t* info, void* user)
{
    SeamCtx* ctx = (SeamCtx*) user;
    if (info->sides.elem_count != 2)
        return;                                     // a domain boundary face

    p4est_iter_face_side_t* s0 = p4est_iter_fside_array_index_int(&info->sides, 0);
    p4est_iter_face_side_t* s1 = p4est_iter_fside_array_index_int(&info->sides, 1);
    if (s0->is_hanging == s1->is_hanging)
        return;                                     // equal levels: no seam
    p4est_iter_face_side_t* coarse = s0->is_hanging ? s1 : s0;
    p4est_iter_face_side_t* fine   = s0->is_hanging ? s0 : s1;

    p4est_quadrant_t* cq = coarse->is.full.quad;
    if (cq == NULL)
        return;                                     // not in the ghost layer

    p4est_t* p4est = info->p4est;
    Seam s;
    s.face = coarse->face;
    s.coarseElem = -1;
    if (!coarse->is.full.is_ghost) {
        s.owner = p4est->mpirank;
        s.coarseElem = localElemIndex(p4est, coarse->treeid,
                                      coarse->is.full.quadid);
        s.globalQuad = (long) p4est->global_first_quadrant[s.owner] + s.coarseElem;
    } else {
        s.owner = ghostOwner(ctx->ghost, coarse->is.full.quadid, p4est->mpisize);
        if (s.owner < 0)
            return;
        s.globalQuad = (long) p4est->global_first_quadrant[s.owner]
                     + (long) cq->p.piggy3.local_num;
    }

    // midpoint of the coarse face, in that octant's tree coordinates
    const p4est_qcoord_t len = P4EST_QUADRANT_LEN(cq->level);
    const p4est_qcoord_t h = len / 2;
    p4est_qcoord_t mx = cq->x, my = cq->y;
    switch (s.face) {
        case 0:  my += h;              break;       // -x
        case 1:  mx += len; my += h;   break;       // +x
        case 2:  mx += h;              break;       // -y
        default: mx += h;   my += len; break;       // +y
    }
    s.mid[0] = s.mid[1] = s.mid[2] = 0.;
    p4est_qcoord_to_vertex(p4est->connectivity, coarse->treeid, mx, my, s.mid);

    // The fine octants are listed in the face's own z-order, so the half that
    // comes first touches the coarse face's LOW corner and its far corner is
    // the midpoint - hence 1-k. Only valid while faces are aligned, which holds
    // for the brick connectivity oxley builds (orientation 0).
    for (int k = 0; k < 2; ++k) {
        s.fineElem[k] = -1;
        s.fineOwner[k] = -1;
        s.fineCorner[k] = rectFaceCorners[fine->face][1-k];
        if (fine->is.hanging.quad[k] == NULL)
            continue;
        if (!fine->is.hanging.is_ghost[k]) {
            s.fineElem[k] = localElemIndex(p4est, fine->treeid,
                                           fine->is.hanging.quadid[k]);
            s.fineOwner[k] = p4est->mpirank;
        } else {
            s.fineOwner[k] = ghostOwner(ctx->ghost, fine->is.hanging.quadid[k],
                                        p4est->mpisize);
        }
    }
    ctx->out->push_back(s);
}

/// Collects every 2:1 seam this rank can see. Builds its own ghost layer if the
/// domain is not already holding one, since the callback needs the far side of
/// a seam whose other half lives on another rank.
void collectSeams(p4est_t* p4est, p4est_ghost_t* keptGhost,
                  std::vector<Seam>& seams)
{
    seams.clear();
    p4est_ghost_t* ghost = keptGhost;
    const bool ownGhost = (ghost == NULL);
    if (ownGhost)
        ghost = p4est_ghost_new(p4est, P4EST_CONNECT_FULL);
    SeamCtx ctx;
    ctx.ghost = ghost;
    ctx.out = &seams;
    p4est_iterate(p4est, ghost, (void*) &ctx, NULL, collectSeam, NULL);
    if (ownGhost)
        p4est_ghost_destroy(ghost);
}

} // anonymous namespace

MeshAccess Rectangle::getMeshAccess(bool materializeHanging) const
{
    MeshAccess m;
    m.numDim = 2;
    m.nodesPerElement = nodes->vnodes;                 // 4 for degree 1
    m.numNodes = nodes->num_local_nodes;
    m.numOwnedNodes = nodes->owned_count;
    m.numElements = nodes->num_local_elements;
    m.globalNodeOffset = (long) nodes->global_offset;
    m.numRealNodes = m.numNodes;
    m.mastersPerConstrainedNode = 2;                   // an edge midpoint

    m.nodeCoords.assign((size_t) m.numNodes * m.numDim, 0.0);
    m.nodeLnodesId.resize(m.numNodes);
    m.elementNodes.resize((size_t) m.numElements * m.nodesPerElement);
    m.elementTags.resize(m.numElements);

    // node tags. m_nodeTags is indexed by local node, the same order used here,
    // and populateSampleIds() has sized it - but a domain that has not been
    // through that yet would leave it short, so copy only what is there.
    m.nodeTags.assign(m.numNodes, 0);
    for (long i = 0; i < m.numNodes && i < (long) m_nodeTags.size(); ++i)
        m.nodeTags[i] = (long) m_nodeTags[i];

    // global node ids: owned nodes are contiguous from global_offset, ghost
    // nodes carry their explicit global id in nonlocal_nodes.
    for (long i = 0; i < m.numOwnedNodes; ++i)
        m.nodeLnodesId[i] = m.globalNodeOffset + i;
    for (long i = m.numOwnedNodes; i < m.numNodes; ++i)
        m.nodeLnodesId[i] = (long) nodes->nonlocal_nodes[i - m.numOwnedNodes];

    // walk the leaves in lnodes element order, filling connectivity, tags and
    // (deduplicated by node index) coordinates.
    //
    // A HANGING corner has no node of its own: its element_nodes slot holds the
    // far master, a node that lies OUTSIDE this element. So the slot's coordinate
    // must not be written - it would move the master to the hanging position -
    // and with materializeHanging the slot is redirected to a node created here.
    const int V = m.nodesPerElement;
    std::vector<bool> haveCoords(m.numNodes, false);
    long e = 0;
    for (p4est_topidx_t treeid = p4est->first_local_tree;
         treeid <= p4est->last_local_tree; ++treeid) {
        p4est_tree_t * tree = p4est_tree_array_index(p4est->trees, treeid);
        sc_array_t * quads = &tree->quadrants;
        const p4est_locidx_t Q = (p4est_locidx_t) quads->elem_count;
        for (p4est_locidx_t q = 0; q < Q; ++q, ++e) {
            p4est_quadrant_t * quad = p4est_quadrant_array_index(quads, q);
            const quadrantData * qd = (const quadrantData *) quad->p.user_data;
            m.elementTags[e] = qd ? qd->quadTag : 0;
            const p4est_qcoord_t len = P4EST_QUADRANT_LEN(quad->level);

            int hangingCorner[P4EST_CHILDREN];
            const bool anyHanging = getHangingNodes(nodes->face_code[e],
                                                    hangingCorner);

            for (int c = 0; c < V; ++c) {
                const long ni = (long) nodes->element_nodes[(size_t) e * V + c];
                m.elementNodes[(size_t) e * V + c] = ni;
                const int cx = c & 1;          // z-order corner: bit0=x, bit1=y
                const int cy = (c >> 1) & 1;
                double xy[3] = {0., 0., 0.};
                p4est_qcoord_to_vertex(p4est->connectivity, treeid,
                                       quad->x + cx * len, quad->y + cy * len, xy);
                if (anyHanging && hangingCorner[c] >= 0)
                    continue;                  // a master, not this corner
                m.nodeCoords[(size_t) ni * m.numDim + 0] = xy[0];
                m.nodeCoords[(size_t) ni * m.numDim + 1] = xy[1];
                haveCoords[ni] = true;
            }

        }
    }

    // Materialise the hanging positions, one per 2:1 seam. Driven by the seam
    // list rather than by face_code, because face_code marks the FINE side only
    // and the coarse side needs the node just as much: the node is a corner of
    // the finer neighbour, so the coarse element does not list it, yet it must
    // become a vertex of that element's simplices.
    std::vector<Seam> seams;
    if (materializeHanging) {
        collectSeams(p4est, m_ghost, seams);
        m.elementFaceHangingNode.assign((size_t) m.numElements * 4, -1);
        for (size_t si = 0; si < seams.size(); ++si) {
            const Seam& s = seams[si];

            // The two masters are the endpoints of the coarse face. Read them
            // off the coarse element when it is local; otherwise off a fine one,
            // where the hanging slot holds the far master and getHangingNodes
            // names the near one. Same two nodes either way.
            long masters[2] = {-1, -1};
            if (s.coarseElem >= 0) {
                for (int k = 0; k < 2; ++k)
                    masters[k] = m.elementNodes[(size_t) s.coarseElem * V
                                              + rectFaceCorners[s.face][k]];
            } else {
                for (int k = 0; k < 2 && masters[0] < 0; ++k) {
                    const long fe = s.fineElem[k];
                    if (fe < 0)
                        continue;
                    int hc[P4EST_CHILDREN];
                    if (!getHangingNodes(nodes->face_code[fe], hc))
                        continue;
                    const int c = s.fineCorner[k];
                    if (hc[c] < 0)
                        continue;
                    masters[0] = (long) nodes->element_nodes[(size_t) fe * V + c];
                    masters[1] = (long) nodes->element_nodes[(size_t) fe * V + hc[c]];
                }
            }

            const long ni = m.numNodes++;
            m.nodeCoords.push_back(s.mid[0]);
            m.nodeCoords.push_back(s.mid[1]);
            m.nodeLnodesId.push_back(-1);   // set by finaliseNodeNumbering
            // no node of the domain tagged this position, so it takes its
            // masters' tag when they agree - see inheritedTag()
            m.nodeTags.push_back(inheritedTag(m, masters, 2));
            m.constrainedNodes.push_back(ni);

            // Who WRITES this node in the output. Not the coarse octant's rank,
            // which owns it for the export numbering: weipa emits the octant
            // mesh, and a hanging node is not a corner of the coarse quad - only
            // of the two fine ones. A coarse-side owner would write a point that
            // appears in none of its own cells and has no value to give it, and
            // the point comes out zero. So the writer is the lowest-numbered
            // rank holding a FINE octant of this seam, which every rank can work
            // out from the ghost layer without communicating.
            int writer = -1;
            for (int k = 0; k < 2; ++k) {
                if (s.fineOwner[k] >= 0 && (writer < 0 || s.fineOwner[k] < writer))
                    writer = s.fineOwner[k];
            }
            m.hangingWriterRank.push_back(writer >= 0 ? writer : s.owner);
            for (int k = 0; k < m.mastersPerConstrainedNode; ++k) {
                m.constraintMasters.push_back(k < 2 ? masters[k] : -1);
                m.constraintWeights.push_back(k < 2 ? 0.5 : 0.);
            }

            // the coarse element reaches it by face, the fine ones by corner
            if (s.coarseElem >= 0)
                m.elementFaceHangingNode[(size_t) s.coarseElem * 4 + s.face] = ni;
            for (int k = 0; k < 2; ++k) {
                if (s.fineElem[k] >= 0)
                    m.elementNodes[(size_t) s.fineElem[k] * V + s.fineCorner[k]] = ni;
            }
        }
    }


    // Any node our elements only ever reference as a hanging slot has no
    // coordinate yet; see farMasterCoords(). (MPI only.)
    {
        std::vector<long> fmIds;
        std::vector<double> fmXY;
        farMasterCoords(fmIds, fmXY);
        for (size_t i = 0; i < fmIds.size(); ++i) {
            const long ni = fmIds[i];
            if (ni < 0 || ni >= m.numRealNodes || haveCoords[ni])
                continue;
            haveCoords[ni] = true;
            m.nodeCoords[(size_t) ni * m.numDim + 0] = fmXY[2*i];
            m.nodeCoords[(size_t) ni * m.numDim + 1] = fmXY[2*i+1];
        }
    }

    // Boundary faces. A quadrant face lies on the domain boundary when the
    // quadrant touches the tree boundary in that direction AND the connectivity
    // sends that tree face back to itself, which is p4est's encoding for "no
    // neighbour". This is topological, unlike updateFaceElementCount() which
    // compares coordinates against the domain extent.
    // Hanging nodes need no treatment here, and cannot in 2D: a hanging node is
    // the midpoint of a face that HAS a finer neighbour, so it lies on an
    // interior face, half an edge away from either end - never on the boundary,
    // and never at an octant corner. A boundary edge is therefore always a plain
    // Line2 on two corners. (3D is not like this: a boundary FACE is still never
    // a seam, but its edges can be, so it can carry hanging edge nodes and become
    // a polygon of up to 8 vertices.)
    {
        // face -> its two corners in z-order indexing, wound so that the domain
        // lies to the LEFT of the directed edge, i.e. the outward normal is the
        // clockwise rotation of the tangent. Matches finley's Mesh_rec4.
        static const int faceCorner[4][2] = {{2,0}, {1,3}, {0,1}, {3,2}};
        static const long faceTag[4] = {1, 2, 10, 20};  // left, right, bottom, top
        const p4est_connectivity_t* conn = p4est->connectivity;
        m.nodesPerFace = 2;
        long le = 0;
        for (p4est_topidx_t treeid = p4est->first_local_tree;
             treeid <= p4est->last_local_tree; ++treeid) {
            p4est_tree_t* tree = p4est_tree_array_index(p4est->trees, treeid);
            sc_array_t* quads = &tree->quadrants;
            const p4est_locidx_t Q = (p4est_locidx_t) quads->elem_count;
            for (p4est_locidx_t q = 0; q < Q; ++q, ++le) {
                p4est_quadrant_t* quad = p4est_quadrant_array_index(quads, q);
                const p4est_qcoord_t len = P4EST_QUADRANT_LEN(quad->level);
                for (int f = 0; f < 4; ++f) {
                    if (conn->tree_to_tree[treeid * 4 + f] != treeid ||
                        conn->tree_to_face[treeid * 4 + f] != f)
                        continue;               // a neighbouring tree is there
                    bool touches;
                    switch (f) {
                        case 0:  touches = (quad->x == 0); break;
                        case 1:  touches = (quad->x + len == P4EST_ROOT_LEN); break;
                        case 2:  touches = (quad->y == 0); break;
                        default: touches = (quad->y + len == P4EST_ROOT_LEN); break;
                    }
                    if (!touches)
                        continue;
                    // from m.elementNodes, not element_nodes: a boundary edge of
                    // a fine element can end at a hanging corner, and the face
                    // must name the node that is actually there
                    for (int c = 0; c < 2; ++c)
                        m.faceNodes.push_back(
                                m.elementNodes[(size_t) le * V + faceCorner[f][c]]);
                    m.faceTags.push_back(faceTag[f]);
                    m.faceElements.push_back(le);
                }
            }
        }
        m.numFaces = (long) m.faceTags.size();
    }

    std::vector<long> ownedPerRank(m_mpiInfo->size);
    for (int i = 0; i < m_mpiInfo->size; ++i)
        ownedPerRank[i] = (long) nodes->global_owned_count[i];

    // The export numbering, derived from replicated data alone: every rank
    // computes the same id for a shared node, with nothing exchanged. Built
    // BEFORE finaliseNodeNumbering, which uses the export id of a materialised
    // node as the key naming it across ranks.
    {
        const int size = m_mpiInfo->size;
        std::vector<long> realOffset(size + 1, 0);
        m.finleyDistribution.assign(size + 1, 0);
        for (int r = 0; r < size; ++r) {
            const long quads = (long) p4est->global_first_quadrant[r+1]
                             - (long) p4est->global_first_quadrant[r];
            realOffset[r+1] = realOffset[r] + ownedPerRank[r];
            m.finleyDistribution[r+1] = m.finleyDistribution[r]
                                      + ownedPerRank[r] + 4 * quads;
        }

        // lnodes nodes keep their position within their owner's block. The
        // owner follows from the lnodes id, since lnodes numbers each rank's
        // owned nodes consecutively - so ghosts need no lookup either.
        m.nodeFinleyId.assign(m.numNodes, -1);
        for (long i = 0; i < m.numRealNodes; ++i) {
            const long g = m.nodeLnodesId[i];
            int r = 0;                              // realOffset is sorted
            while (r + 1 < size && realOffset[r+1] <= g)
                ++r;
            m.nodeFinleyId[i] = m.finleyDistribution[r] + (g - realOffset[r]);
        }

        // hanging nodes sit above their owner's lnodes nodes, at the slot their
        // (octant, face) key names
        for (size_t si = 0; si < seams.size(); ++si) {
            const Seam& s = seams[si];
            const long firstQuad = (long) p4est->global_first_quadrant[s.owner];
            m.nodeFinleyId[m.constrainedNodes[si]] =
                    m.finleyDistribution[s.owner] + ownedPerRank[s.owner]
                  + 4 * (s.globalQuad - firstQuad) + s.face;
        }
    }

    finaliseNodeNumbering(m, ownedPerRank);

    return m;
}

//protected
dim_t Rectangle::getNumFaceElements() const
{
    return m_faceCount[0]+m_faceCount[1]+m_faceCount[2]+m_faceCount[3];
    
// #ifdef ENABLE_OPENMP
//     long numFaceElements[omp_get_num_threads()] = {0};
//     for (p4est_topidx_t t = p4est->first_local_tree; t <= p4est->last_local_tree; t++) // Loop over every tree
//     {
//         p4est_tree_t * currenttree = p4est_tree_array_index(p4est->trees, t);
//         sc_array_t * tquadrants = &currenttree->quadrants;
//         p4est_locidx_t Q = (p4est_locidx_t) tquadrants->elem_count;
//     #pragma omp for
//         for(int q = 0; q < Q; q++) // Loop over every quadrant within the tree
//         {
//             p4est_quadrant_t * quad = p4est_quadrant_array_index(tquadrants, q);
//             int l = quad->level;
//             double xy[3];
//             p4est_qcoord_to_vertex(p4est->connectivity, t, quad->x, quad->y, xy);
//             long e = NodeIDs.find(std::make_pair(xy[0],xy[1]))->second;
//             quadrantData * quaddata = (quadrantData *) quad->p.user_data;
//             numFaceElements[omp_get_thread_num()] += quaddata->m_faceOffset[0] ||
//                                quaddata->m_faceOffset[1] ||
//                                quaddata->m_faceOffset[2] ||
//                                quaddata->m_faceOffset[3];

//         }
//     }
//     long answer = 0;
//     for(int i = 0; i < omp_get_num_threads(); i++)
//         answer+=numFaceElements[i];
//     return answer;
// #else
//     long numFaceElements = 0;
//     for (p4est_topidx_t t = p4est->first_local_tree; t <= p4est->last_local_tree; t++) // Loop over every tree
//     {
//         p4est_tree_t * currenttree = p4est_tree_array_index(p4est->trees, t);
//         sc_array_t * tquadrants = &currenttree->quadrants;
//         p4est_locidx_t Q = (p4est_locidx_t) tquadrants->elem_count;
//         for(int q = 0; q < Q; q++) // Loop over every quadrant within the tree
//         {
//             p4est_quadrant_t * quad = p4est_quadrant_array_index(tquadrants, q);
//             quadrantData * quaddata = (quadrantData *) quad->p.user_data;
//             numFaceElements += quaddata->m_faceOffset[0] ||
//                                quaddata->m_faceOffset[1] ||
//                                quaddata->m_faceOffset[2] ||
//                                quaddata->m_faceOffset[3];

//         }
//     }
//     return numFaceElements;
// #endif
}

dim_t Rectangle::getNumDOF() const
{
    // owned nodes only (each owned node is one real DOF). Ghost/shared nodes
    // are columns, not rows. (MPI: A6.)
    return nodes ? (dim_t) nodes->owned_count : 0;
}

void Rectangle::updateTreeIDs()
{
    treeIDs.clear();
    for(p4est_topidx_t treeid = p4est->first_local_tree; treeid <= p4est->last_local_tree; ++treeid) {
        p4est_tree_t * tree = p4est_tree_array_index(p4est->trees, treeid);
        sc_array_t * tquadrants = &tree->quadrants;
        p4est_locidx_t Q = (p4est_locidx_t) tquadrants->elem_count;
#pragma omp parallel for
        for(int q = 0; q < Q; ++q) { // Loop over the elements attached to the tree
            p4est_quadrant_t * quad = p4est_quadrant_array_index(tquadrants, q);
            treeIDs[std::make_pair(quad->x,quad->y)]=treeid;
        }
    }
}

void Rectangle::updateRowsColumns()
{
    // Conforming lnodes numbering: every local node is a real DOF numbered
    // identically (serially). ownSample() consults m_dofMap, so keep it as the
    // identity map. The legacy hanging-node connectivity build (the old
    // coordinate-hash `indices`/update_RC path) is retired; getConnections now
    // derives the matrix graph directly from lnodes.
    const dim_t n = getNumNodes();
    m_dofMap.assign(n, 0);
    for(dim_t i = 0; i < n; ++i)
        m_dofMap[i] = i;
}

#ifdef ESYS_HAVE_TRILINOS
//protected
esys_trilinos::TrilinosGraph_ptr Rectangle::getTrilinosGraph() const
{   
    // if (m_graph.is_null()) {
    //     m_graph = createTrilinosGraph(myRows, myColumns);
    // }
    // return m_graph;
    m_graph = createTrilinosGraph(myRows, myColumns);
    return m_graph;
}
#endif

#ifdef ESYS_HAVE_PASO
//protected
paso::SystemMatrixPattern_ptr Rectangle::getPasoMatrixPattern(
                                                    bool reducedRowOrder,
                                                    bool reducedColOrder) const
{
    if (m_pattern.get())
        return m_pattern;

    // first call - create pattern, then return
    paso::Connector_ptr conn(getPasoConnector());
    const dim_t numDOF = getNumDOF();
    const dim_t numShared = conn->send->numSharedComponents; //todo
    const dim_t numNeighbours = conn->send->neighbour.size();
    const std::vector<index_t>& offsetInShared(conn->send->offsetInShared);
    const index_t* sendShared = conn->send->shared;

    // these are for the couple blocks
    std::vector<IndexVector> colIndices(numDOF);
    std::vector<IndexVector> rowIndices(numShared);

    for(dim_t i=0; i<numNeighbours; i++) {
        const dim_t start = offsetInShared[i];
        const dim_t end = offsetInShared[i+1];
        for(dim_t j = start; j < end; j++) {
            if (j > start)
                doublyLink(colIndices, rowIndices, sendShared[j-1], j);
            doublyLink(colIndices, rowIndices, sendShared[j], j);
            if (j < end-1)
                doublyLink(colIndices, rowIndices, sendShared[j+1], j);
        }
    }
#pragma omp parallel for
    for(dim_t i = 0; i < numShared; i++) {
        sort(rowIndices[i].begin(), rowIndices[i].end());
    }

    // create main and couple blocks
    paso::Pattern_ptr mainPattern = createPasoPattern(getConnections(), numDOF);
    paso::Pattern_ptr colPattern = createPasoPattern(colIndices, numShared);
    paso::Pattern_ptr rowPattern = createPasoPattern(rowIndices, numDOF);

    // allocate paso distribution
    IndexVector m_nodeDistribution = getNodeDistribution();
    escript::Distribution_ptr distribution(new escript::Distribution(m_mpiInfo, m_nodeDistribution));

    // finally create the system matrix pattern
    m_pattern.reset(new paso::SystemMatrixPattern(MATRIX_FORMAT_DEFAULT,
            distribution, distribution, mainPattern, colPattern, rowPattern,
            conn, conn));
    return m_pattern;
}
#endif // ESYS_HAVE_PASO

//private
void Rectangle::populateSampleIds()
{
    m_nodeTags.assign(getNumNodes(), 0);
    updateTagsInUse(Nodes);

    m_elementTags.assign(getNumElements(), 0);
    updateTagsInUse(Elements);
}

void Rectangle::updateFaceElementCount()
{
    // real face count per side (each face pushed once, deduplicated below). The
    // old code initialised to -1 to cancel a size()-1 consumer loop / a corner
    // duplicate; both are gone now, so count from 0. (A6.)
    for(int i = 0; i < 4; i++)
        m_faceCount[i]=0;

    NodeIDsTop.clear();
    NodeIDsBottom.clear();
    NodeIDsLeft.clear();
    NodeIDsRight.clear();

    const int V = nodes->vnodes;   // 4 corners (z-order matches lxy below)
    long e = 0;                    // running local leaf index (lnodes order)
    for(p4est_topidx_t treeid = p4est->first_local_tree; treeid <= p4est->last_local_tree; ++treeid)
    {
        p4est_tree_t * tree = p4est_tree_array_index(p4est->trees, treeid);
        sc_array_t * tquadrants = &tree->quadrants;
        p4est_locidx_t Q = (p4est_locidx_t) tquadrants->elem_count;
        for(int q = 0; q < Q; ++q, ++e)
        {
            p4est_quadrant_t * quad = p4est_quadrant_array_index(tquadrants, q);
            p4est_qcoord_t l = P4EST_QUADRANT_LEN(quad->level);
            // int k = q - Q + nodeIncrements[treeid - p4est->first_local_tree];
            p4est_qcoord_t lxy[4][2] = {{0,0},{l,0},{0,l},{l,l}};
            double xy[4][3] = {{0}};
            int nodeids[4]={-1};
            bool do_check_yes_no[4]={false};
            for(int n = 0; n < 4; n++)
            {
                p4est_qcoord_to_vertex(p4est->connectivity, treeid, quad->x+lxy[n][0], quad->y+lxy[n][1], xy[n]);
                nodeids[n]=(int) nodes->element_nodes[(size_t) e * V + n];

                if(n==0)
                    do_check_yes_no[n]=true;
                else if(n==1 && xy[n][0]==forestData.m_lxy[0])
                    do_check_yes_no[n]=true;
                else if(n==2 && xy[n][1]==forestData.m_lxy[1])
                    do_check_yes_no[n]=true;
                else if(n==3 && xy[n][0]==forestData.m_lxy[0] && xy[n][1]==forestData.m_lxy[1])
                    do_check_yes_no[n]=true;
                else
                    do_check_yes_no[n]=false;
            }

            for(int n = 0; n < 4; n++)
            {
                if(do_check_yes_no[n] == false)
                    continue;

                borderNodeInfo tmp;
                // tmp.nodeid=NodeIDs.find(std::make_pair(xy[n][0],xy[n][1]))->second;
                tmp.nodeid=nodeids[n];
                tmp.neighbours[0]=nodeids[0];
                tmp.neighbours[1]=nodeids[1];
                tmp.neighbours[2]=nodeids[2];
                tmp.neighbours[3]=nodeids[3];
                tmp.x=quad->x;
                tmp.y=quad->y;
                tmp.level=quad->level;
                tmp.treeid=treeid;
                tmp.quadIndex=e;

                // Push each boundary FACE exactly once, keyed on its canonical
                // corner (SW=0 for left/bottom, SE=1 for right, NW=2 for top).
                // Without this, a corner octant lies on two boundaries and its
                // shared corner would push the same face twice; the old code
                // masked that with a size()-1 loop which drops a REAL face on any
                // rank that owns an edge but not its corner (MPI). (A6.)
                if(n==0 && isLeftBoundaryNode(quad, n, treeid, l))
                {
                    NodeIDsLeft.push_back(tmp);
                    m_faceCount[0]++;
                }

                if(n==1 && isRightBoundaryNode(quad, n, treeid, l))
                {
                    NodeIDsRight.push_back(tmp);
                    m_faceCount[1]++;
                }

                if(n==0 && isBottomBoundaryNode(quad, n, treeid, l))
                {
                    NodeIDsBottom.push_back(tmp);
                    m_faceCount[2]++;
                }

                if(n==2 && isTopBoundaryNode(quad, n, treeid, l))
                {
                    NodeIDsTop.push_back(tmp);
                    m_faceCount[3]++;
                }
            
                #ifdef OXLEY_ENABLE_DEBUG_FACEELEMENTS_POINTS
                    double xyz[3];
                    p4est_qcoord_to_vertex(p4est->connectivity, treeid, quad->x+lxy[n][0], quad->y+lxy[n][1], &xyz[n]);
                    std::cout << nodeids[n] << ": quad (x,y) = " << xyz[0] << ", " << xyz[1] << ") ";
                    if(isLeftBoundaryNode(quad, n, treeid, l))
                        std::cout << "L";
                    if(isRightBoundaryNode(quad, n, treeid, l))
                        std::cout << "R";
                    if(isBottomBoundaryNode(quad, n, treeid, l))
                        std::cout << "B";
                    if(isTopBoundaryNode(quad, n, treeid, l))
                        std::cout << "T";
                    std::cout << std::endl;
                #endif
            }
        }
    }

    // Remove duplicates
    // for(int i = 1; i < NodeIDsLeft.size(); i++)
    //     if((NodeIDsLeft[i].treeid == NodeIDsLeft[i-1].treeid))
    //     {
    //         NodeIDsLeft.erase(NodeIDsLeft.begin()+i);
    //         i--;
    //         m_faceCount[0]--;
    //     }
    // for(int i = 1; i < NodeIDsRight.size(); i++)
    //     if(NodeIDsRight[i].treeid == NodeIDsRight[i-1].treeid)
    //     {
    //         NodeIDsRight.erase(NodeIDsRight.begin()+i);
    //         i--;
    //         m_faceCount[1]--;
    //     }
    // for(int i = 1; i < NodeIDsBottom.size(); i++)
    //     if(NodeIDsBottom[i].treeid == NodeIDsBottom[i-1].treeid)
    //     {
    //         NodeIDsBottom.erase(NodeIDsBottom.begin()+i);
    //         i--;
    //         m_faceCount[2]--;
    //     }
    // for(int i = 1; i < NodeIDsTop.size(); i++)
    //     if(NodeIDsTop[i].treeid == NodeIDsTop[i-1].treeid)
    //     {
    //         NodeIDsTop.erase(NodeIDsTop.begin()+i);
    //         i--;
    //         m_faceCount[3]--;
    //     }


    const index_t LEFT=1, RIGHT=2, BOTTOM=10, TOP=20;
    m_faceTags.clear();
    const index_t faceTag[] = { LEFT, RIGHT, BOTTOM, TOP };
    m_faceOffset.assign(4, -1);
    index_t offset=0;
    for (size_t i=0; i<4; i++) {
        if (m_faceCount[i]>0) {
            m_faceOffset[i]=offset;
            offset+=m_faceCount[i];
            m_faceTags.insert(m_faceTags.end(), m_faceCount[i], faceTag[i]);
        }
    }

#ifdef OXLEY_ENABLE_DEBUG_FACEELEMENTS
    std::cout << "NodeIDsLeft" << std::endl;
    for(int i = 0; i < NodeIDsLeft.size();i++)
        std::cout << NodeIDsLeft[i].nodeid << " ";
    std::cout << std::endl;
    std::cout << "NodeIDsRight" << std::endl;
    for(int i = 0; i < NodeIDsRight.size();i++)
        std::cout << NodeIDsRight[i].nodeid << " ";
    std::cout << std::endl;
    std::cout << "NodeIDsTop" << std::endl;
    for(int i = 0; i < NodeIDsTop.size();i++)
        std::cout << NodeIDsTop[i].nodeid << " ";
    std::cout << std::endl;
    std::cout << "NodeIDsBottom" << std::endl;
    for(int i = 0; i < NodeIDsBottom.size();i++)
        std::cout << NodeIDsBottom[i].nodeid << " ";
    std::cout << std::endl;
    std::cout << "-------------------------------------------------------" << std::endl;
#endif

    // set face tags
    setTagMap("left", LEFT);
    setTagMap("right", RIGHT);
    setTagMap("bottom", BOTTOM);
    setTagMap("top", TOP);
    updateTagsInUse(FaceElements);


    // Update faceElementId
    const dim_t NFE = getNumFaceElements();
    m_faceId.resize(NFE);
    for (dim_t k=0; k<NFE; k++)
        m_faceId[k]=k;
}

// This is a wrapper that converts the p4est node information into an IndexVector
IndexVector Rectangle::getNodeDistribution() const
{
    return m_nodeDistribution;
}

// This is a wrapper that converts the p4est node information into an IndexVector
void Rectangle::updateNodeDistribution() 
{
    m_nodeDistribution.clear();
    m_nodeDistribution.assign(MAXP4ESTNODES,0);

    int counter =0;
    const int V = nodes->vnodes;
    long e = 0;
    for(p4est_topidx_t treeid = p4est->first_local_tree; treeid <= p4est->last_local_tree; ++treeid)
    {
        p4est_tree_t * tree = p4est_tree_array_index(p4est->trees, treeid);
        sc_array_t * tquadrants = &tree->quadrants;
        p4est_locidx_t Q = (p4est_locidx_t) tquadrants->elem_count;
        for(int q = 0; q < Q; ++q, ++e)
        {
            for(int n = 0; n < V; n++)
                m_nodeDistribution[counter++] = (long) nodes->element_nodes[(size_t) e * V + n];
        }
    }
    m_nodeDistribution.shrink_to_fit();
}

// updates m_elementIDs()
void Rectangle::updateElementIds()
{
    // element sample ids are the running local leaf indices [0, numElements)
    const dim_t ne = getNumElements();
    m_elementId.clear();
    m_elementId.resize(ne);
    for(dim_t i = 0; i < ne; ++i)
        m_elementId[i] = i;
    m_elementId.shrink_to_fit();
}

//private
std::vector<IndexVector> Rectangle::getConnections(bool includeShared) const
{
    // returns a vector v of size numDOF where v[i] is a vector with indices
    // of DOFs connected to i (up to 9 in 2D).
    // In other words this method returns the occupied (local) matrix columns
    // for all (local) matrix rows.
    // If includeShared==true then connections to non-owned DOFs are also
    // returned (i.e. indices of the column couplings)

    long numNodes = getNumNodes();
    std::vector< std::vector<escript::DataTypes::index_t> > indices(numNodes);

    // Build the node adjacency graph directly from the lnodes element->node
    // connectivity (no coordinate hashing). Every corner of a leaf is coupled
    // to every other corner of that leaf.
    const int V = nodes->vnodes;   // 4 corners (degree-1 lnodes)
    long e = 0;                    // running local leaf index (lnodes order)
    for(p4est_topidx_t treeid = p4est->first_local_tree; treeid <= p4est->last_local_tree; ++treeid)
    {
        p4est_tree_t * tree = p4est_tree_array_index(p4est->trees, treeid);
        sc_array_t * tquadrants = &tree->quadrants;
        p4est_locidx_t Q = (p4est_locidx_t) tquadrants->elem_count;
        for(int q = 0; q < Q; ++q, ++e) // Loop over all quadrants
        {
            long lni[4];
            for(int i = 0; i < V; i++)
                lni[i] = (long) nodes->element_nodes[(size_t) e * V + i];

            for(int i = 0; i < V; i++)
            {
                for(int j = 0; j < V; j++)
                {
                    bool dup = false;
                    for(int k = 0; k < indices[lni[i]].size(); k++)
                        if(indices[lni[i]][k] == lni[j])
                        {
                            dup = true;
                            break;
                        }
                    if(dup == false)
                        indices[lni[i]].push_back(lni[j]);
                }
            }
        }
    }

    // MPI: add couplings contributed by the ghost element halo. Only OWNED rows
    // matter (they must be complete); ghost elements supply the 2nd-layer column
    // couplings for owned boundary nodes. (A6.)
    const long nDOF = getNumDOF();
    const long nGhost = (long) m_ghostElemNodes.size() / (V ? V : 1);
    for(long g = 0; g < nGhost; ++g)
    {
        const index_t* lni = &m_ghostElemNodes[(size_t) g * V];
        for(int i = 0; i < V; i++)
        {
            const long row = (long) lni[i];
            if(row >= nDOF)          // only owned rows are assembled/kept
                continue;
            for(int j = 0; j < V; j++)
            {
                const index_t col = lni[j];
                bool dup = false;
                for(int k = 0; k < indices[row].size(); k++)
                    if(indices[row][k] == col) { dup = true; break; }
                if(!dup)
                    indices[row].push_back(col);
            }
        }
    }

    // for(int i = 0; i < hanging_face_orientation.size(); i++)
    // {       
    //     // Calculate the node ids
    //     double xy[3]={0};
    //     p4est_qcoord_to_vertex(p4est->connectivity, hanging_face_orientation[i].treeid, 
    //                                                 hanging_face_orientation[i].x, 
    //                                                 hanging_face_orientation[i].y, xy); 
        
    //     bool boundary[4];
    //     boundary[0]  = xy[1] == forestData.m_lxy[1];
    //     boundary[1]  = xy[1] == forestData.m_origin[1];
    //     boundary[2]  = xy[0] == forestData.m_lxy[0];
    //     boundary[3]  = xy[0] == forestData.m_origin[0];
        
    //     bool hangingBorderNode = boundary[0] || boundary[1] || boundary[2] || boundary[3];
    //     if( !hangingBorderNode ) 
    //         continue;

    //     // parameters needed below
    //     int dir;
    //     for(int j = 0; j<4; j++)
    //         if(boundary[j] == true)
    //         {
    //             dir=j;
    //             break;
    //         }

    //     int position=hanging_face_orientation[i].position;

    //     // This node's id
    //     long nodeid = NodeIDs.find(std::make_pair(xy[0],xy[1]))->second;

    //     // then lni0 lni1
    //     p4est_qcoord_t l = P4EST_QUADRANT_LEN(hanging_face_orientation[i].neighbour_l);
    //     // int dx[4][2]=  {{0,0},{0,0},{l,0},{0,l}};
    //     // int dy[4][2]=  {{l,0},{0,l},{l,l},{l,l}};

    //     // p4est_qcoord_to_vertex(p4est->connectivity, parent.tree, 
    //     //                         parent->x+dxy[0][dir], 
    //     //                         parent->y+dxy[1][dir], xy);
    //     // long lni0 = NodeIDs.find(std::make_pair(xy[0],xy[1]))->second;
    //     // p4est_qcoord_to_vertex(p4est->connectivity, parent.tree, 
    //     //                         parent->x+dxy[0][dir], 
    //     //                         parent->y+dxy[1][dir], xy);
    //     // long lni1 = NodeIDs.find(std::make_pair(xy[0],xy[1]))->second;

    //     // get the center
    //     // int dxy2[2][4] 

    //     // 0
    //     // 1,1
    //     // 1
    //     // -1,1
    //     // 2
    //     // 1,-1
    //     // 3
    //     // 0,0
        

    //     // p4est_qcoord_to_vertex(p4est->connectivity, hanging_face_orientation[i].neighbour_tree, 
    //     //                         parent.x+dxy[0][position], 
    //     //                         parent.y+dxy[1][position], xy);
    //     // long lni2 = NodeIDs.find(std::make_pair(xy[0],xy[1]))->second;

    //     // // add info 
    //     // indices[nodeid].push_back(lni0);
    //     // indices[nodeid].push_back(lni1);
    //     // // indices[nodeid].push_back(lni2);

    //     // indices[lni0].push_back(nodeid);
    //     // indices[lni1].push_back(nodeid);
    //     // // indices[lni2].push_back(nodeid);
    // }    

// Sorting
    for(int i = 0; i < numNodes; i++){
        std::sort(indices[i].begin(), indices[i].begin()+indices[i].size());
    }

#ifdef OXLEY_ENABLE_DEBUG_GETCONNECTIONS
    std::cout << "Rectangle::getConnections" << std::endl;
    for(int i = 0; i < numNodes; i++) {
        std::cout << i << ": ";
        for(auto j = 0; j < indices[i].size(); j++)
            std::cout << indices[i][j] << ", ";
        std::cout << std::endl;
    }
#endif

    return indices;
}

bool Rectangle::operator==(const AbstractDomain& other) const
{
    const Rectangle* o=dynamic_cast<const Rectangle*>(&other);
    if (o) {
        return ((p4est_checksum(p4est) == p4est_checksum(o->p4est)));
            // && (forestData == o->forestData)); //TODO
    }
    return false;
}

//protected
void Rectangle::assembleGradient(escript::Data& out,
                                 const escript::Data& in) const
{
    if (out.isComplex() && in.isComplex())
        assembleGradientImpl<cplx_t>(out, in);
    else if (!out.isComplex() && !in.isComplex())
        assembleGradientImpl<real_t>(out, in);
    else
        throw ValueError("Gradient: input & output complexity must match.");
}

//protected
template<typename Scalar>
void Rectangle::assembleGradientImpl(escript::Data& out,
                                     const escript::Data& in) const
{
    const dim_t numComp = in.getDataPointSize();
    
    // Find the maximum level of refinement in the mesh
    int max_level = 0;
    for(p4est_topidx_t tree = p4est->first_local_tree; tree <= p4est->last_local_tree; tree++) {
        p4est_tree_t * tree_t = p4est_tree_array_index(p4est->trees, tree);
        max_level = tree_t->maxlevel > max_level ? tree_t->maxlevel : max_level;
    }
    
    double cx[3][P4EST_MAXLEVEL] = {{0}};
    double cy[3][P4EST_MAXLEVEL] = {{0}};
#pragma omp parallel for
    for(int i = 0; i <= max_level; i++)
    {
        double m_dx[2] = {forestData.m_dx[0][P4EST_MAXLEVEL-i], 
                          forestData.m_dx[1][P4EST_MAXLEVEL-i]};
        cx[0][i] = 0.21132486540518711775/m_dx[0];
        cx[1][i] = 0.78867513459481288225/m_dx[0];
        cx[2][i] = 1./m_dx[0];
        cy[0][i] = 0.21132486540518711775/m_dx[1];
        cy[1][i] = 0.78867513459481288225/m_dx[1];
        cy[2][i] = 1./m_dx[1];
    }
    const Scalar zero = static_cast<Scalar>(0);

    if (out.getFunctionSpace().getTypeCode() == Elements) {
        out.requireWrite();

        std::vector<Scalar> f_00(numComp, zero);
        std::vector<Scalar> f_01(numComp, zero);
        std::vector<Scalar> f_10(numComp, zero);
        std::vector<Scalar> f_11(numComp, zero);

        const int V = nodes->vnodes;
        long e = 0;
        for (p4est_topidx_t t = p4est->first_local_tree; t <= p4est->last_local_tree; t++) // Loop over every tree
        {
            p4est_tree_t * currenttree = p4est_tree_array_index(p4est->trees, t);
            sc_array_t * tquadrants = &currenttree->quadrants;
            p4est_locidx_t Q = (p4est_locidx_t) tquadrants->elem_count;
            for(int q = 0; q < Q; ++q, ++e) // Loop over every quadrant within the tree
            {
                p4est_quadrant_t * quad = p4est_quadrant_array_index(tquadrants, q);
                long l = quad->level;

                long ids[4];
                for(int n = 0; n < V; ++n) ids[n] = (long) nodes->element_nodes[(size_t) e * V + n];

                #ifdef OXLEY_ENABLE_DEBUG_ASSEMBLE_GRADIENT
                    std::cout << "quad id: " << e << std::endl;
                #endif

                memcpy(&f_00[0], in.getSampleDataRO(ids[0], zero), numComp*sizeof(Scalar));
                memcpy(&f_01[0], in.getSampleDataRO(ids[2], zero), numComp*sizeof(Scalar));
                memcpy(&f_10[0], in.getSampleDataRO(ids[1], zero), numComp*sizeof(Scalar));
                memcpy(&f_11[0], in.getSampleDataRO(ids[3], zero), numComp*sizeof(Scalar));

                // on the fine side of a seam these slots hold masters, not the
                // corner values; see constrainHangingCorners
                int hangingCorner[P4EST_CHILDREN];
                if (getHangingNodes(nodes->face_code[e], hangingCorner)) {
                    Scalar* corner[P4EST_CHILDREN] =
                            { &f_00[0], &f_10[0], &f_01[0], &f_11[0] };
                    constrainHangingCorners(hangingCorner, corner, numComp);
                }

                Scalar* o = out.getSampleDataRW(e, zero);
                for(index_t i = 0; i < numComp; ++i) {
                    o[INDEX3(i,0,0,numComp,2)] = (f_10[i]-f_00[i])*cx[1][l] + (f_11[i]-f_01[i])*cx[0][l];
                    o[INDEX3(i,1,0,numComp,2)] = (f_01[i]-f_00[i])*cy[1][l] + (f_11[i]-f_10[i])*cy[0][l];
                    o[INDEX3(i,0,1,numComp,2)] = (f_10[i]-f_00[i])*cx[1][l] + (f_11[i]-f_01[i])*cx[0][l];
                    o[INDEX3(i,1,1,numComp,2)] = (f_01[i]-f_00[i])*cy[0][l] + (f_11[i]-f_10[i])*cy[1][l];
                    o[INDEX3(i,0,2,numComp,2)] = (f_10[i]-f_00[i])*cx[0][l] + (f_11[i]-f_01[i])*cx[1][l];
                    o[INDEX3(i,1,2,numComp,2)] = (f_01[i]-f_00[i])*cy[1][l] + (f_11[i]-f_10[i])*cy[0][l];
                    o[INDEX3(i,0,3,numComp,2)] = (f_10[i]-f_00[i])*cx[0][l] + (f_11[i]-f_01[i])*cx[1][l];
                    o[INDEX3(i,1,3,numComp,2)] = (f_01[i]-f_00[i])*cy[0][l] + (f_11[i]-f_10[i])*cy[1][l];
                }
            }
        }
    } else if (out.getFunctionSpace().getTypeCode() == ReducedElements) {
        out.requireWrite();

        std::vector<Scalar> f_00(numComp, zero);
        std::vector<Scalar> f_01(numComp, zero);
        std::vector<Scalar> f_10(numComp, zero);
        std::vector<Scalar> f_11(numComp, zero);

        const int V = nodes->vnodes;
        long e = 0;
        for (p4est_topidx_t t = p4est->first_local_tree; t <= p4est->last_local_tree; t++) // Loop over every tree
        {
            p4est_tree_t * currenttree = p4est_tree_array_index(p4est->trees, t);
            sc_array_t * tquadrants = &currenttree->quadrants;
            p4est_locidx_t Q = (p4est_locidx_t) tquadrants->elem_count;
            for(int q = 0; q < Q; ++q, ++e) // Loop over every quadrant within the tree
            {
                p4est_quadrant_t * quad = p4est_quadrant_array_index(tquadrants, q);
                long l = quad->level;

                long ids[4];
                for(int n = 0; n < V; ++n) ids[n] = (long) nodes->element_nodes[(size_t) e * V + n];

                memcpy(&f_00[0], in.getSampleDataRO(ids[0], zero), numComp*sizeof(Scalar));
                memcpy(&f_01[0], in.getSampleDataRO(ids[2], zero), numComp*sizeof(Scalar));
                memcpy(&f_10[0], in.getSampleDataRO(ids[1], zero), numComp*sizeof(Scalar));
                memcpy(&f_11[0], in.getSampleDataRO(ids[3], zero), numComp*sizeof(Scalar));

                // on the fine side of a seam these slots hold masters, not the
                // corner values; see constrainHangingCorners
                int hangingCorner[P4EST_CHILDREN];
                if (getHangingNodes(nodes->face_code[e], hangingCorner)) {
                    Scalar* corner[P4EST_CHILDREN] =
                            { &f_00[0], &f_10[0], &f_01[0], &f_11[0] };
                    constrainHangingCorners(hangingCorner, corner, numComp);
                }

                Scalar* o = out.getSampleDataRW(e, zero);

                for(index_t i = 0; i < numComp; ++i) {
                    o[INDEX3(i,0,0,numComp,2)] = (f_10[i] + f_11[i] - f_00[i] - f_01[i])*cx[2][l] * 0.5;
                    o[INDEX3(i,1,0,numComp,2)] = (f_01[i] + f_11[i] - f_00[i] - f_10[i])*cy[2][l] * 0.5;
                } 
            }
        }
    } else if (out.getFunctionSpace().getTypeCode() == FaceElements) {
        out.requireWrite();

        for (p4est_topidx_t t = p4est->first_local_tree; t <= p4est->last_local_tree; t++) // Loop over every tree
        {
            p4est_tree_t * currenttree = p4est_tree_array_index(p4est->trees, t);
            sc_array_t * tquadrants = &currenttree->quadrants;
            p4est_locidx_t Q = (p4est_locidx_t) tquadrants->elem_count;
// #pragma omp parallel for
            for(int q = 0; q < Q; ++q) // Loop over every quadrant within the tree
            {

                p4est_quadrant_t * quad = p4est_quadrant_array_index(tquadrants, q);
                double xy[3];
                p4est_qcoord_to_vertex(p4est->connectivity, t, quad->x, quad->y, xy);
                // long e = getQuadID(NodeIDs.find(std::make_pair(xy[0],xy[1]))->second);
                // quadrantData * quaddata = (quadrantData *) quad->p.user_data;

                std::vector<Scalar> f_00(numComp, zero);
                std::vector<Scalar> f_01(numComp, zero);
                std::vector<Scalar> f_10(numComp, zero);
                std::vector<Scalar> f_11(numComp, zero);

                if (m_faceOffset[0] > -1) {
                    for (index_t k=0; k<NodeIDsLeft.size(); k++) {
                        borderNodeInfo tmp = NodeIDsLeft[k];
                        long l = tmp.level;
                        gatherCornersConstrained(in, tmp, numComp, zero,
                                                 f_00, f_10, f_01, f_11);

                        Scalar* o = out.getSampleDataRW(m_faceOffset[0]+k, zero);
                        for(index_t i = 0; i < numComp; ++i) {
                            o[INDEX3(i,0,0,numComp,2)] = (f_10[i]-f_00[i])*cx[1][l] + (f_11[i]-f_01[i])*cx[0][l];
                            o[INDEX3(i,1,0,numComp,2)] = (f_01[i]-f_00[i])*cy[2][l];
                            o[INDEX3(i,0,1,numComp,2)] = (f_10[i]-f_00[i])*cx[0][l] + (f_11[i]-f_01[i])*cx[1][l];
                            o[INDEX3(i,1,1,numComp,2)] = (f_01[i]-f_00[i])*cy[2][l];
                        } // end of component loop i
                    }
                } // end of face 0
                if (m_faceOffset[1] > -1) {
                    for (index_t k=0; k<NodeIDsRight.size(); k++) {
                        borderNodeInfo tmp = NodeIDsRight[k];
                        long l = tmp.level;
                        gatherCornersConstrained(in, tmp, numComp, zero,
                                                 f_00, f_10, f_01, f_11);

                        Scalar* o = out.getSampleDataRW(m_faceOffset[1]+k, zero);
                        for(index_t i = 0; i < numComp; ++i) {
                            o[INDEX3(i,0,0,numComp,2)] = (f_10[i]-f_00[i])*cx[1][l] + (f_11[i]-f_01[i])*cx[0][l];
                            o[INDEX3(i,1,0,numComp,2)] = (f_11[i]-f_10[i])*cy[2][l];
                            o[INDEX3(i,0,1,numComp,2)] = (f_10[i]-f_00[i])*cx[0][l] + (f_11[i]-f_01[i])*cx[1][l];
                            o[INDEX3(i,1,1,numComp,2)] = (f_11[i]-f_10[i])*cy[2][l];
                        } // end of component loop i
                    }
                } // end of face 1
                if (m_faceOffset[2] > -1) {
                    for (index_t k=0; k<NodeIDsBottom.size(); k++) {
                        borderNodeInfo tmp = NodeIDsBottom[k];
                        long l = tmp.level;
                        gatherCornersConstrained(in, tmp, numComp, zero,
                                                 f_00, f_10, f_01, f_11);

                        Scalar* o = out.getSampleDataRW(m_faceOffset[2]+k, zero);
                        for(index_t i = 0; i < numComp; ++i) {
                            o[INDEX3(i,0,0,numComp,2)] = (f_10[i]-f_00[i])*cx[2][l];
                            o[INDEX3(i,1,0,numComp,2)] = (f_01[i]-f_00[i])*cy[1][l] + (f_11[i]-f_10[i])*cy[0][l];
                            o[INDEX3(i,0,1,numComp,2)] = (f_10[i]-f_00[i])*cx[2][l];
                            o[INDEX3(i,1,1,numComp,2)] = (f_01[i]-f_00[i])*cy[0][l] + (f_11[i]-f_10[i])*cy[1][l];
                        } // end of component loop i
                    }
                } // end of face 2
                if (m_faceOffset[3] > -1) {
                    for (index_t k=0; k<NodeIDsTop.size(); k++) {
                        borderNodeInfo tmp = NodeIDsTop[k];
                        long l = tmp.level;
                        gatherCornersConstrained(in, tmp, numComp, zero,
                                                 f_00, f_10, f_01, f_11);

                        Scalar* o = out.getSampleDataRW(m_faceOffset[3]+k, zero);
                        for(index_t i = 0; i < numComp; ++i) {
                            o[INDEX3(i,0,0,numComp,2)] = (f_11[i]-f_01[i])*cx[2][l];
                            o[INDEX3(i,1,0,numComp,2)] = (f_01[i]-f_00[i])*cy[1][l] + (f_11[i]-f_10[i])*cy[0][l];
                            o[INDEX3(i,0,1,numComp,2)] = (f_11[i]-f_01[i])*cx[2][l];
                            o[INDEX3(i,1,1,numComp,2)] = (f_01[i]-f_00[i])*cy[0][l] + (f_11[i]-f_10[i])*cy[1][l];
                        } // end of component loop i
                    }
                } // end of face 3
            }
        }
    } else if (out.getFunctionSpace().getTypeCode() == ReducedFaceElements) {
        out.requireWrite();

        for (p4est_topidx_t t = p4est->first_local_tree; t <= p4est->last_local_tree; t++) // Loop over every tree
        {
            p4est_tree_t * currenttree = p4est_tree_array_index(p4est->trees, t);
            sc_array_t * tquadrants = &currenttree->quadrants;
            p4est_locidx_t Q = (p4est_locidx_t) tquadrants->elem_count;
// #pragma omp parallel for
            for(int q = 0; q < Q; ++q) // Loop over every quadrant within the tree
            {

                p4est_quadrant_t * quad = p4est_quadrant_array_index(tquadrants, q);
                double xy[3];
                p4est_qcoord_to_vertex(p4est->connectivity, t, quad->x, quad->y, xy);
                // long e = getQuadID(NodeIDs.find(std::make_pair(xy[0],xy[1]))->second);
                // quadrantData * quaddata = (quadrantData *) quad->p.user_data;

                std::vector<Scalar> f_00(numComp, zero);
                std::vector<Scalar> f_01(numComp, zero);
                std::vector<Scalar> f_10(numComp, zero);
                std::vector<Scalar> f_11(numComp, zero);

                if (m_faceOffset[0] > -1) {
                    for (index_t k=0; k<NodeIDsLeft.size(); k++) {
                        borderNodeInfo tmp = NodeIDsLeft[k];
                        long l = tmp.level;
                        gatherCornersConstrained(in, tmp, numComp, zero,
                                                 f_00, f_10, f_01, f_11);

                        Scalar* o = out.getSampleDataRW(m_faceOffset[0]+k, zero);
                        for(index_t i = 0; i < numComp; ++i) {
                            o[INDEX3(i,0,0,numComp,2)] = (f_10[i] + f_11[i] - f_00[i] - f_01[i])*cx[2][l] * 0.5;
                            o[INDEX3(i,1,0,numComp,2)] = (f_01[i]-f_00[i])*cy[2][l];
                        } // end of component loop i
                    }
                } // end of face 0
                if (m_faceOffset[1] > -1) {
                    for (index_t k=0; k<NodeIDsRight.size(); k++) {
                        borderNodeInfo tmp = NodeIDsRight[k];
                        long l = tmp.level;
                        gatherCornersConstrained(in, tmp, numComp, zero,
                                                 f_00, f_10, f_01, f_11);

                        Scalar* o = out.getSampleDataRW(m_faceOffset[1]+k, zero);
                        for(index_t i = 0; i < numComp; ++i) {
                            o[INDEX3(i,0,0,numComp,2)] = (f_10[i] + f_11[i] - f_00[i] - f_01[i])*cx[2][l] * 0.5;
                            o[INDEX3(i,1,0,numComp,2)] = (f_11[i]-f_10[i])*cy[2][l];
                        } // end of component loop i
                    }
                } // end of face 1
                if (m_faceOffset[2] > -1) {
                    for (index_t k=0; k<NodeIDsBottom.size(); k++) {
                        borderNodeInfo tmp = NodeIDsBottom[k];
                        long l = tmp.level;
                        gatherCornersConstrained(in, tmp, numComp, zero,
                                                 f_00, f_10, f_01, f_11);

                        Scalar* o = out.getSampleDataRW(m_faceOffset[2]+k, zero);
                        for(index_t i = 0; i < numComp; ++i) {
                            o[INDEX3(i,0,0,numComp,2)] = (f_10[i]-f_00[i])*cx[2][l];
                            o[INDEX3(i,1,0,numComp,2)] = (f_01[i] + f_11[i] - f_00[i] - f_10[i])*cy[2][l] * 0.5;
                        } // end of component loop i
                    }
                } // end of face 2
                if (m_faceOffset[3] > -1) {
                    for (index_t k=0; k<NodeIDsTop.size(); k++) {
                        borderNodeInfo tmp = NodeIDsTop[k];
                        long l = tmp.level;
                        gatherCornersConstrained(in, tmp, numComp, zero,
                                                 f_00, f_10, f_01, f_11);

                        Scalar* o = out.getSampleDataRW(m_faceOffset[3]+k, zero);
                        for(index_t i = 0; i < numComp; ++i) {
                            o[INDEX3(i,0,0,numComp,2)] = (f_11[i]-f_01[i])*cx[2][l];
                            o[INDEX3(i,1,0,numComp,2)] = (f_01[i] + f_11[i] - f_00[i] - f_10[i])*cy[2][l] * 0.5;
                        } // end of component loop i
                    }
                }
            }
        }
    }
}

//protected
void Rectangle::assembleIntegrate(std::vector<real_t>& integrals,
                                  const escript::Data& arg) const
{
    assembleIntegrateImpl<real_t>(integrals, arg);
}

//protected
void Rectangle::assembleIntegrate(std::vector<cplx_t>& integrals,
                                  const escript::Data& arg) const
{
    assembleIntegrateImpl<cplx_t>(integrals, arg);
}

//private
template<typename Scalar>
void Rectangle::assembleIntegrateImpl(std::vector<Scalar>& integrals,
                                      const escript::Data& arg) const
{
    const dim_t numComp = arg.getDataPointSize();
    const int fs = arg.getFunctionSpace().getTypeCode();
    const Scalar zero = static_cast<Scalar>(0);

    bool HavePointData = arg.getFunctionSpace().getTypeCode() == Points;

#ifdef ESYS_MPI
    if(HavePointData && escript::getMPIRankWorld() == 0) {
#else
    if(HavePointData) {
#endif
        integrals[0] += arg.getNumberOfTaggedValues();
    } else if (fs == Elements && arg.actsExpanded()) {
       
        std::vector<Scalar> int_local(numComp, zero);
        long id = 0;   // running local leaf index (element sample order)
        for(p4est_topidx_t treeid = p4est->first_local_tree; treeid <= p4est->last_local_tree; ++treeid)
        {
            p4est_tree_t * tree = p4est_tree_array_index(p4est->trees, treeid);
            sc_array_t * tquadrants = &tree->quadrants;
            p4est_locidx_t Q = (p4est_locidx_t) tquadrants->elem_count;
            for(int q = 0; q < Q; ++q, ++id)
            {
                p4est_quadrant_t * quad = p4est_quadrant_array_index(tquadrants, q);

                real_t w = forestData.m_dx[0][P4EST_MAXLEVEL-quad->level]
                         * forestData.m_dx[1][P4EST_MAXLEVEL-quad->level]
                         / 4.;

                const Scalar* f = arg.getSampleDataRO(id, zero);
                for (index_t i = 0; i < numComp; ++i) {
                    const Scalar f0 = f[INDEX2(i,0,numComp)];
                    const Scalar f1 = f[INDEX2(i,1,numComp)];
                    const Scalar f2 = f[INDEX2(i,2,numComp)];
                    const Scalar f3 = f[INDEX2(i,3,numComp)];
                    int_local[i] += (f0+f1+f2+f3)*w;
                }
            }
        }
        for (index_t i=0; i<numComp; i++)
        {
            integrals[i] += int_local[i];
        }

    } else if (fs==ReducedElements || (fs==Elements && !arg.actsExpanded())) {
        
        // const
        std::vector<Scalar> int_local(numComp, 0);
        long id = 0;   // running local leaf index (element sample order)
        for(p4est_topidx_t treeid = p4est->first_local_tree; treeid <= p4est->last_local_tree; ++treeid)
        {
            p4est_tree_t * tree = p4est_tree_array_index(p4est->trees, treeid);
            sc_array_t * tquadrants = &tree->quadrants;
            p4est_locidx_t Q = (p4est_locidx_t) tquadrants->elem_count;
            for(int q = 0; q < Q; ++q, ++id)
            {
                p4est_quadrant_t * quad = p4est_quadrant_array_index(tquadrants, q);
                const Scalar* f = arg.getSampleDataRO(id, zero);
                real_t w = forestData.m_dx[0][P4EST_MAXLEVEL-quad->level]
                         * forestData.m_dx[1][P4EST_MAXLEVEL-quad->level];
                for (index_t i = 0; i < numComp; ++i) {
                    int_local[i] += f[i]*w;
                }
            }
        }

        for (index_t i = 0; i < numComp; i++)
            integrals[i] += int_local[i];

    } else if (fs == FaceElements && arg.actsExpanded()) {

#pragma omp parallel
        {
            std::vector<Scalar> int_local(numComp, zero);
            if (m_faceOffset[0] > -1) {
#pragma omp for nowait
                for (index_t k=0; k<NodeIDsLeft.size(); k++) {
                    borderNodeInfo tmp = NodeIDsLeft[k];
                    const real_t w1 = forestData.m_dx[1][P4EST_MAXLEVEL-tmp.level]/2.;
                    const Scalar* f = arg.getSampleDataRO(m_faceOffset[0]+k, zero);
                    for (index_t i=0; i < numComp; ++i) {
                        const Scalar f0 = f[INDEX2(i,0,numComp)];
                        const Scalar f1 = f[INDEX2(i,1,numComp)];
                        int_local[i] += (f0+f1)*w1;
                        #ifdef OXLEY_ENABLE_DEBUG_INTEGRATE
                            std::cout << "w1=" << w1 << ", faceoffset=0 f0=" << f0 << ", f1=" << f1 << std::endl;
                        #endif
                    }  // end of component loop i
                }
            }

            if (m_faceOffset[1] > -1) {
#pragma omp for nowait
                for (index_t k=0; k<NodeIDsRight.size(); k++) {
                    borderNodeInfo tmp = NodeIDsRight[k];
                    const real_t w1 = forestData.m_dx[1][P4EST_MAXLEVEL-tmp.level]/2.;
                    const Scalar* f = arg.getSampleDataRO(m_faceOffset[1]+k, zero);
                    for (index_t i = 0; i < numComp; ++i) {
                        const Scalar f0 = f[INDEX2(i,0,numComp)];
                        const Scalar f1 = f[INDEX2(i,1,numComp)];
                        int_local[i] += (f0+f1)*w1;
                        #ifdef OXLEY_ENABLE_DEBUG_INTEGRATE
                            std::cout << "w1=" << w1 << ", faceoffset=1 f0=" << f0 << ", f1=" << f1 << std::endl;
                        #endif
                    }  // end of component loop i
                }
            }

            if (m_faceOffset[2] > -1) {
#pragma omp for nowait
                for (index_t k=0; k<NodeIDsBottom.size(); k++) {
                    borderNodeInfo tmp = NodeIDsBottom[k];
                    const real_t w0 = forestData.m_dx[0][P4EST_MAXLEVEL-tmp.level]/2.;
                    const Scalar* f = arg.getSampleDataRO(m_faceOffset[2]+k, zero);
                    for (index_t i = 0; i < numComp; ++i) {
                        const Scalar f0 = f[INDEX2(i,0,numComp)];
                        const Scalar f1 = f[INDEX2(i,1,numComp)];
                        int_local[i] += (f0+f1)*w0;
                        #ifdef OXLEY_ENABLE_DEBUG_INTEGRATE
                            std::cout << "w0=" << w0 << ", faceoffset=2 f0=" << f0 << ", f1=" << f1 << std::endl;
                        #endif
                    }  // end of component loop i
                }
            }

            if (m_faceOffset[3] > -1) {
#pragma omp for nowait
                for (index_t k=0; k<NodeIDsTop.size(); k++) {
                    borderNodeInfo tmp = NodeIDsTop[k];
                    const real_t w0 = forestData.m_dx[0][P4EST_MAXLEVEL-tmp.level]/2.;
                    const Scalar* f = arg.getSampleDataRO(m_faceOffset[3]+k, zero);
                    for (index_t i = 0; i < numComp; ++i) {
                        const Scalar f0 = f[INDEX2(i,0,numComp)];
                        const Scalar f1 = f[INDEX2(i,1,numComp)];
                        int_local[i] += (f0+f1)*w0;
                        #ifdef OXLEY_ENABLE_DEBUG_INTEGRATE
                            std::cout << "w0=" << w0 << ", faceoffset=3 f0=" << f0 << ", f1=" << f1 << std::endl;
                        #endif
                    }  // end of component loop i
                }
            }
#pragma omp critical
            for (index_t i = 0; i < numComp; i++)
                integrals[i] += int_local[i];
        } // end of parallel section
    } else if (fs==ReducedFaceElements || (fs==FaceElements && !arg.actsExpanded())) {
        std::vector<Scalar> int_local(numComp, 0);
        if (m_faceOffset[0] > -1) {
#pragma omp for nowait
            for (index_t k=0; k<NodeIDsLeft.size(); k++) {
                borderNodeInfo tmp = NodeIDsLeft[k];
                const Scalar* f = arg.getSampleDataRO(m_faceOffset[0]+k, zero);
                #ifdef OXLEY_ENABLE_DEBUG_INTEGRATE
                    std::cout << "offset=" << m_faceOffset[0]+k << " ";
                #endif
                for (index_t i = 0; i < numComp; ++i) {
                    int_local[i] += f[i]*forestData.m_dx[1][P4EST_MAXLEVEL-tmp.level];
                    #ifdef OXLEY_ENABLE_DEBUG_INTEGRATE
                        std::cout << "faceoffset=0 f[i]=" << f[i] << ", dx=" << forestData.m_dx[1][P4EST_MAXLEVEL-tmp.level] << std::endl;
                    #endif
                }
            }
        }

        if (m_faceOffset[1] > -1) {
#pragma omp for nowait
            for (index_t k=0; k<NodeIDsRight.size(); k++) {
                borderNodeInfo tmp = NodeIDsRight[k];
                const Scalar* f = arg.getSampleDataRO(m_faceOffset[1]+k, zero);
                #ifdef OXLEY_ENABLE_DEBUG_INTEGRATE
                    std::cout << "offset=" << m_faceOffset[1]+k << " ";
                #endif
                for (index_t i = 0; i < numComp; ++i) {
                    int_local[i] += f[i]*forestData.m_dx[1][P4EST_MAXLEVEL-tmp.level];
                    #ifdef OXLEY_ENABLE_DEBUG_INTEGRATE
                        std::cout << "faceoffset=1 f[i]=" << f[i] << ", dx=" << forestData.m_dx[1][P4EST_MAXLEVEL-tmp.level] << std::endl;
                    #endif
                }
            }
        }

        if (m_faceOffset[2] > -1) {
#pragma omp for nowait
            for (index_t k=0; k<NodeIDsBottom.size(); k++) {
                borderNodeInfo tmp = NodeIDsBottom[k];
                const Scalar* f = arg.getSampleDataRO(m_faceOffset[2]+k, zero);
                #ifdef OXLEY_ENABLE_DEBUG_INTEGRATE
                    std::cout << "offset=" << m_faceOffset[2]+k << " ";
                #endif
                for (index_t i = 0; i < numComp; ++i) {
                    int_local[i] += f[i]*forestData.m_dx[0][P4EST_MAXLEVEL-tmp.level];
                    #ifdef OXLEY_ENABLE_DEBUG_INTEGRATE
                        std::cout << "faceoffset=2 f[i]=" << f[i] << ", dx=" << forestData.m_dx[0][P4EST_MAXLEVEL-tmp.level] << std::endl;
                    #endif
                }
            }
        }

        if (m_faceOffset[3] > -1) {
#pragma omp for nowait
            for (index_t k=0; k<NodeIDsTop.size(); k++) {
                borderNodeInfo tmp = NodeIDsTop[k];
                const Scalar* f = arg.getSampleDataRO(m_faceOffset[3]+k, zero);
                #ifdef OXLEY_ENABLE_DEBUG_INTEGRATE
                    std::cout << "offset=" << m_faceOffset[3]+k << " ";
                #endif
                for (index_t i = 0; i < numComp; ++i) {
                    int_local[i] += f[i]*forestData.m_dx[0][P4EST_MAXLEVEL-tmp.level];
                    #ifdef OXLEY_ENABLE_DEBUG_INTEGRATE
                        std::cout << "faceoffset=3 f[i]=" << f[i] << ", dx=" << forestData.m_dx[0][P4EST_MAXLEVEL-tmp.level] << std::endl;
                    #endif
                }
            }
        }

#pragma omp critical
        for (index_t i = 0; i < numComp; i++)
            integrals[i] += int_local[i];
    } // function space selector
}

//protected
void Rectangle::nodesToDOF(escript::Data& out, const escript::Data& in) const
{
    // Nodes -> DegreesOfFreedom: the owned nodes are the DOFs (lnodes orders
    // owned nodes first), so copy the first getNumDOF() node samples. Ghost
    // node values belong to other ranks and are dropped. (MPI: A6.)
    const dim_t numComp = in.getDataPointSize();
    out.requireWrite();
    const dim_t nDOF = getNumDOF();
    const real_t zero = 0;
#pragma omp parallel for
    for (index_t i = 0; i < nDOF; i++) {
        const real_t* src = in.getSampleDataRO(i, zero);
        std::copy(src, src+numComp, out.getSampleDataRW(i, zero));
    }
    return;
    // legacy structured-grid implementation below (dead):

//     const dim_t numComp = in.getDataPointSize();
//     out.requireWrite();

//     const index_t left = getFirstInDim(0);
//     const index_t bottom = getFirstInDim(1);
//     const dim_t nDOF0 = getNumDOFInAxis(0);
//     const dim_t nDOF1 = getNumDOFInAxis(1);
// #pragma omp parallel for
//     for(index_t i=0; i<nDOF1; i++) {
//         for(index_t j=0; j<nDOF0; j++) {
//             const index_t n=j+left+(i+bottom)*m_NN[0];
//             const double* src=in.getSampleDataRO(n);
//             copy(src, src+numComp, out.getSampleDataRW(j+i*nDOF0));
//         }
//     }
}

// Updates m_faceOffset for each quadrant
void Rectangle::updateFaceOffset()
{
    p4est_iterate(p4est, NULL, NULL, update_node_faceoffset, NULL, NULL);
}

void Rectangle::updateMeshInformation()
{
    refineMesh("MARE2DEM");
}

static inline void
brick_linear_to_xyz (p4est_topidx_t ti, const int logx[P4EST_DIM],
                     const int rankx[P4EST_DIM], p4est_topidx_t tx[P4EST_DIM])
{
    int lastlog = 0;

    for(int i = 0; i < P4EST_DIM; i++) {
        tx[i] = 0;
    }

    for(int i = 0; i < P4EST_DIM - 1; i++) {
        p4est_topidx_t tempx[3] = { 0, 0, 0 };
        int logi = logx[rankx[i]] - lastlog;
        int idx[3] = { -1, -1, -1 };
        int c = 0;

        for(int k = 0; k < P4EST_DIM - i; k++) {
            int d = rankx[i + k];
            idx[d] = 0;
        }
        for(int k = 0; k < P4EST_DIM; k++) {
            if (idx[k] == 0) {
                idx[k] = c++;
            }
        }

        for(int j = 0; j < logi; j++) {
            int base = (P4EST_DIM - i) * j;
            int shift = (P4EST_DIM - i - 1) * j;

            for(int k = 0; k < P4EST_DIM; k++) {
                int id = idx[k];

                if (id >= 0) {
                    tempx[k] |= (ti & (1 << (base + id))) >> (shift + id);
                }
            }
        }
        for(int k = 0; k < P4EST_DIM; k++) {
            tx[k] += (tempx[k] << lastlog);
        }
        lastlog += logi;
        ti >>= (P4EST_DIM - i) * logi;
    }
    tx[rankx[P4EST_DIM - 1]] += (ti << lastlog);
}

static inline p4est_topidx_t
brick_xyz_to_linear (const p4est_topidx_t tx[P4EST_DIM],
                     const int logx[P4EST_DIM], const int rankx[P4EST_DIM])
{
    int lastlog = logx[rankx[P4EST_DIM - 2]];
    p4est_topidx_t ti = tx[rankx[P4EST_DIM - 1]] >> lastlog;

    for(int i = P4EST_DIM - 2; i >= 0; i--) {
        p4est_topidx_t tempx[3] = { 0, 0, 0 };
        int logi = (i == 0) ? lastlog : lastlog - logx[rankx[i - 1]];
        int idx[3] = { -1, -1, -1 };
        int c = 0;

        for(int k = 0; k < P4EST_DIM - i; k++) {
            int d = rankx[i + k];

            idx[d] = 0;
        }
        for(int k = 0; k < P4EST_DIM; k++) {
            if (idx[k] == 0) {
                idx[k] = c++;
            }
        }

        ti <<= (P4EST_DIM - i) * logi;
        lastlog -= logi;
        for(int k = 0; k < P4EST_DIM; k++) {
            tempx[k] = tx[k] >> lastlog;
        }
        for(int j = 0; j < logi; j++) {
            int shift = (P4EST_DIM - i - 1) * j;

            for(int k = 0; k < P4EST_DIM; k++) {
                int id = idx[k];

                if (id >= 0) {
                    ti |= (tempx[k] & (1 << j)) << (shift + id);
                }
            }
        }
    }
    return ti;
}

// This is a modified version of p4est_connectivity_new_brick
p4est_connectivity_t * Rectangle::new_rectangle_connectivity(
                            int n0, int n1, int periodic_a, int periodic_b, 
                            double x0, double y0, double x1, double y1)
{
    // Number of nodes 
    const p4est_topidx_t m = (p4est_topidx_t) n0;
    const p4est_topidx_t n = (p4est_topidx_t) n1;    

    ESYS_ASSERT(m > 0 && n > 0, "n0 and n1 must be greater than zero.");
    const p4est_topidx_t num_trees = m * n;
    ESYS_ASSERT(num_trees <= MAXTREES ,"n0*n1 must be less than MAXTREES.");

    // Number of corners in each direction
    P4EST_ASSERT(periodic_a == 0 || periodic_a == 1);
    P4EST_ASSERT(periodic_b == 0 || periodic_b == 1);
    const p4est_topidx_t mc = periodic_a ? m : (m - 1);
    const p4est_topidx_t nc = periodic_b ? n : (n - 1);
    const p4est_topidx_t num_corners = mc * nc;
    const p4est_topidx_t num_vertices = (m + 1) * (n + 1);

    // Corners to Trees        
    const p4est_topidx_t num_ctt = P4EST_CHILDREN * num_corners;

    // Other
    const int periodic[P4EST_DIM] = { periodic_a, periodic_b };
    const p4est_topidx_t max[P4EST_DIM] = { m - 1, n - 1 };

    // Trees to faces, trees to corners    
    p4est_topidx_t tf[P4EST_FACES] = {0};
    p4est_topidx_t tc[P4EST_CHILDREN] = {0};

    // Coordinates
    p4est_topidx_t coord[P4EST_DIM] = {0};
    p4est_topidx_t coord2[P4EST_DIM] = {0};
    p4est_topidx_t ttemp = 0;

    // Counters for the number of vertices
    p4est_topidx_t vcount = 0, vicount = 0;

    // Size of the grid spacing
    double dx = (x1 - x0) / n0;
    double dy = (y1 - y0) / n1;
    
    // Connectivity
    p4est_connectivity_t * conn = p4est_connectivity_new(num_vertices, num_trees, num_corners, num_ctt);
   
    //Corner to Tree offsets
    p4est_topidx_t * ctt_offset = conn->ctt_offset;
#pragma omp parallel for
    for(p4est_topidx_t ti = 0; ti < num_corners + 1; ti++) {
        ctt_offset[ti] = 4 * ti;
    }

    // Tree to Vertices
    p4est_topidx_t * tree_to_vertex = conn->tree_to_vertex;
#pragma omp parallel for
    for(p4est_topidx_t ti = 0; ti < 4 * num_trees; ti++) {
        tree_to_vertex[ti] = -1;
    }

    int logx[P4EST_DIM] = {SC_LOG2_32 (m - 1) + 1, SC_LOG2_32 (n - 1) + 1};
    int c[P4EST_DIM] = {0};
    int rankx[P4EST_DIM] = {0};
    if (logx[0] <= logx[1]) {
        rankx[0] = 0;
        rankx[1] = 1;
    } else {
        rankx[0] = 1;
        rankx[1] = 0;
    }

    p4est_topidx_t n_iter = (1 << logx[0]) * (1 << logx[1]);

    // Allocate memory
    p4est_topidx_t * linear_to_tree =  P4EST_ALLOC(p4est_topidx_t, n_iter);
    p4est_topidx_t * tree_to_corner2 = P4EST_ALLOC(p4est_topidx_t, num_trees);

    p4est_topidx_t tj = 0;
    p4est_topidx_t tk = 0;
    for(p4est_topidx_t ti = 0; ti < n_iter; ti++) {
        brick_linear_to_xyz(ti, logx, rankx, coord);
        p4est_topidx_t tx = coord[0];
        p4est_topidx_t ty = coord[1];
        if (tx < m && ty < n && 1) 
        {
            linear_to_tree[ti] = tj;
            if ((tx < m - 1 || periodic_a) && (ty < n - 1 || periodic_b) && 1)
                tree_to_corner2[tj] = tk++;
            else 
                tree_to_corner2[tj] = -1;
            tj++;
        }
        else 
        {
            linear_to_tree[ti] = -1;
        }
    }
    P4EST_ASSERT(tj == num_trees);
    P4EST_ASSERT(tk == num_corners);

    double * vertices = conn->vertices;
    p4est_topidx_t * tree_to_tree = conn->tree_to_tree;
    int8_t * tree_to_face = conn->tree_to_face;
    p4est_topidx_t * tree_to_corner = conn->tree_to_corner;
    p4est_topidx_t * corner_to_tree = conn->corner_to_tree;
    int8_t * corner_to_corner = conn->corner_to_corner;
#ifdef OXLEY_PRINT_VERTICES
    std::cout << "using the vertices..." << std::endl;
#endif
    for(p4est_topidx_t ti = 0; ti < n_iter; ti++) {
        brick_linear_to_xyz(ti, logx, rankx, coord);
        p4est_topidx_t tx = coord[0];
        p4est_topidx_t ty = coord[1];
        if(tx < m && ty < n) {
            tj = linear_to_tree[ti];
            P4EST_ASSERT(tj >= 0);
            for(int i = 0; i < 2; i++) {
                for(int j = 0; j < 2; j++) {
                    int l = 2 * i + j;
                    coord2[0] = ((tx + ((i == 0) ? (2 * j - 1) : 0)) + m) % m;
                    coord2[1] = ((ty + ((i == 1) ? (2 * j - 1) : 0)) + n) % n;
                    tf[l] = brick_xyz_to_linear (coord2, logx, rankx);
                    P4EST_ASSERT(tf[l] < n_iter);
                    tf[l] = linear_to_tree[tf[l]];
                    P4EST_ASSERT(tf[l] >= 0);
                }
            }
        
            for(int i = 0; i < 4; i++) {
                coord2[0] = ((tx + (((i & 1) == 0) ? -1 : 1)) + m) % m;
                coord2[1] = ((ty + ((((i >> 1) & 1) == 0) ? -1 : 1)) + n) % n;
                tc[i] = brick_xyz_to_linear (coord2, logx, rankx);
                P4EST_ASSERT(tc[i] < n_iter);
                tc[i] = linear_to_tree[tc[i]];
                P4EST_ASSERT(tc[i] >= 0);
            }

            for(int i = 0; i < 2; i++) {
                for(int j = 0; j < 2; j++) {
                    int l = i * 2 + j;
                    if (!periodic[i] &&
                        ((coord[i] == 0 && j == 0) || (coord[i] == max[i] && j == 1))) {
                        tree_to_tree[tj * 4 + l] = tj;
                        tree_to_face[tj * 4 + l] = (int8_t) l;
                    }
                    else {
                        tree_to_tree[tj * 4 + l] = tf[l];
                        tree_to_face[tj * 4 + l] = (int8_t) (i * 2 + (j ^ 1));
                    }
                }
            }

            for(int i = 0; i < 4; i++) {
                if(tree_to_corner != NULL) {
                    c[0] = i & 1;
                    c[1] = (i >> 1) & 1;
                    if ((!periodic[0] &&
                         ((coord[0] == 0 && c[0] == 0) ||
                          (coord[0] == max[0] && c[0] == 1))) ||
                        (!periodic[1] &&
                         ((coord[1] == 0 && c[1] == 0) ||
                          (coord[1] == max[1] && c[1] == 1))) ||
                        0) {
                      tree_to_corner[tj * 4 + i] = -1;
                    }
                    else {
                        switch (i) {
                            case 0:
                                ttemp = tc[0];
                                break;
                            case 1:
                                ttemp = tf[2];
                                break;
                            case 2:
                                ttemp = tf[0];
                                break;
                            case 3:
                                ttemp = tj;
                                break;
                            default:
                            SC_ABORT_NOT_REACHED();
                        }
                        ttemp = tree_to_corner2[ttemp];
                        P4EST_ASSERT(ttemp >= 0);
                        tree_to_corner[tj * 4 + i] = ttemp;
                        corner_to_tree[ttemp * 4 + (4 - 1 - i)] = tj;
                        corner_to_corner[ttemp * 4 + (4 - 1 - i)] = (int8_t) i;
                    }
                }

                if (ty > 0 && ((i >> 1) & 1) == 0) {
                    tree_to_vertex[tj * 4 + i] =
                    tree_to_vertex[tf[2] * 4 + i + 2];
                }
                else if (tx > 0 && (i & 1) == 0) {
                    tree_to_vertex[tj * 4 + i] =
                    tree_to_vertex[tf[0] * 4 + i + 1];
                }
                else {
                    tree_to_vertex[tj * 4 + i] = vcount++;
                    vertices[vicount++] = (double) x0 + (dx * (tx + (i & 1)));
                    vertices[vicount++] = (double) y0 + (dy * (ty + ((i >> 1) & 1)));
                    vertices[vicount++] = 0.;
#ifdef OXLEY_PRINT_VERTICES
                    std::cout << "( " << vertices[vicount-3] << ", " 
                                      << vertices[vicount-2] << ", "
                                      << vertices[vicount-1] << " )" << std::endl;
#endif
                }
            }
        }
    }

    P4EST_ASSERT(vcount == num_vertices);
    P4EST_FREE(linear_to_tree);
    P4EST_FREE(tree_to_corner2);
#ifdef OXLEY_ENABLE_DEBUG
    P4EST_ASSERT(p4est_connectivity_is_valid(conn)); //This is very time consuming
#endif
#ifdef OXLEY_PRINT_VERTICES
    std::cout << "using vertices..." << std::endl;
    for(int i=0;i<vicount;i+=3)
        std::cout << vertices[i] << ", " << vertices[i+1] << std::endl;
#endif
    return conn;
}

void Rectangle::addPoints(const std::vector<double>& coords, const std::vector<int>& tags)
{
    // A Dirac point must be claimed by exactly ONE rank in MPI. Each rank finds
    // its nearest OWNED node (owned nodes are the first getNumDOF() local nodes in
    // lnodes order) and its distance; the globally-nearest rank keeps the point
    // (ties broken by lowest rank). Searching nearest-LOCAL-node on every rank (the
    // old behaviour) placed and assembled every point on every rank -> the source
    // was multiplied by the rank count. (A6.)
    const dim_t nOwned = getNumDOF();
    const MeshAccess m = getMeshAccess();
    const double x0=forestData.m_origin[0], y0=forestData.m_origin[1];
    const double x1=forestData.m_lxy[0],    y1=forestData.m_lxy[1];
    double ext = x1-x0; if(y1-y0>ext) ext=y1-y0;
    const double tol = 1e-8*ext;

    for (int i = 0; i < (int)tags.size(); i++) {
        const double px = coords[i*m_numDim + 0];
        const double py = coords[i*m_numDim + 1];

        double best = std::numeric_limits<double>::max();
        long bestNode = -1;
        // out-of-domain points are claimed by nobody
        if (!(px<x0-tol || px>x1+tol || py<y0-tol || py>y1+tol)) {
            for (long n = 0; n < nOwned; ++n) {
                const double dx = m.nodeCoords[(size_t)n*2 + 0] - px;
                const double dy = m.nodeCoords[(size_t)n*2 + 1] - py;
                const double d2 = dx*dx + dy*dy;
                if (d2 < best) { best = d2; bestNode = n; }
            }
        }

        // globally-nearest distance, then lowest rank achieving it
        double globalBest = best;
        int winner = (bestNode>=0) ? m_mpiInfo->rank : m_mpiInfo->size;
#ifdef ESYS_MPI
        if (m_mpiInfo->size > 1) {
            MPI_Allreduce(&best, &globalBest, 1, MPI_DOUBLE, MPI_MIN, m_mpiInfo->comm);
            int cand = (bestNode>=0 && best==globalBest) ? m_mpiInfo->rank : m_mpiInfo->size;
            MPI_Allreduce(&cand, &winner, 1, MPI_INT, MPI_MIN, m_mpiInfo->comm);
        }
#endif
        if (bestNode >= 0 && m_mpiInfo->rank == winner) {
            m_diracPointNodeIDs.push_back(borrowSampleReferenceIDs(Nodes)[bestNode]);
            DiracPoint dp;
            dp.node = bestNode; //local (owned)
            dp.tag = tags[i];
            m_diracPoints.push_back(dp);
        }
    }
}


// Calculates a Gaussian blur convolution matrix for 2D
// See wiki article on the subject
double* get2DGauss(unsigned radius, double sigma)
{
    double* arr = new double[(radius*2+1)*(radius*2+1)];
    const double common = M_1_PI * 0.5 / (sigma*sigma);
    const int r = static_cast<int>(radius);
    double total = 0;
    for (int y = -r; y <= r; ++y) {
        for (int x = -r; x <= r; ++x) {
            arr[(x+r)+(y+r)*(r*2+1)]=common*exp(-(x*x+y*y)/(2*sigma*sigma));
            total+=arr[(x+r)+(y+r)*(r*2+1)];
        }
    }
    const double invtotal = 1/total;
    for (size_t p=0; p<(radius*2+1)*(radius*2+1); ++p) {
        arr[p] *= invtotal;
    }
    return arr;
}

// applies conv to source to get a point.
// (xp, yp) are the coords in the source matrix not the destination matrix
double Convolve2D(double* conv, double* source, size_t xp, size_t yp,
                  unsigned radius, size_t width)
{
    const size_t bx = xp-radius, by=yp-radius;
    const size_t sbase = bx+by*width;
    double result = 0;
    for (int y=0; y<2*radius+1; ++y) {
        for (int x=0; x<2*radius+1; ++x) {
            result += conv[x+y*(2*radius+1)] * source[sbase + x+y*width];
        }
    }
    return result;
}


/* This is a wrapper for filtered (and non-filtered) randoms
 * For detailed doco see randomFillWorker
 */
escript::Data Rectangle::randomFill(const escript::DataTypes::ShapeType& shape,
                                const escript::FunctionSpace& what, long seed,
                                const bp::tuple& filter) const
{
    int numvals=escript::DataTypes::noValues(shape);
    if (len(filter) > 0 && numvals != 1)
        throw escript::NotImplementedError("Oxley only supports filters for scalar data.");

    escript::Data res = randomFillWorker(shape, seed, filter);
    if (res.getFunctionSpace() != what) {
        escript::Data r(res, what);
        return r;
    }
    return res;
}


/* This routine produces a Data object filled with smoothed random data.
 * The dimensions of the rectangle being filled are internal[0] x internal[1]
 * points. A parameter radius gives the size of the stencil used for the
 * smoothing.  A point on the left hand edge for example, will still require
 * `radius` extra points to the left in order to complete the stencil.
 *
 * All local calculation is done on an array called `src`, with
 * dimensions = ext[0] * ext[1], where ext[i]= internal[i]+2*radius.
 *
 * Now for MPI there is overlap to deal with. We need to share both the
 * overlapping values themselves but also the external region.
 *
 * In a hypothetical 1-D case:
 *
 * 1234567 would be split into two ranks thus:
 * 123(4)  (4)567     [4 being a shared element]
 *
 * If the radius is 2. There will be padding elements on the outside:
 * pp123(4)  (4)567pp
 *
 * To ensure that 4 can be correctly computed on both ranks, values from the
 * other rank need to be known.
 *
 * pp123(4)56   23(4)567pp
 *
 * Now in our case, we set all the values 23456 on the left rank and send them
 * to the right hand rank.
 *
 * So the edges _may_ need to be shared at a distance `inset` from all
 * boundaries.
 *
 * inset=2*radius+1
 * This is to ensure that values at distance `radius` from the
 * shared/overlapped element that Oxley has.
 */
escript::Data Rectangle::randomFillWorker(
                        const escript::DataTypes::ShapeType& shape, long seed,
                        const bp::tuple& filter) const
{
    unsigned int radius=0;  // these are only used by gaussian
    double sigma=0.5;

    unsigned int numvals=escript::DataTypes::noValues(shape);

    if (len(filter) == 0) {
        // nothing special required here yet
    } else if (len(filter) == 3) {
        bp::extract<std::string> ex(filter[0]);
        if (!ex.check() || (ex()!="gaussian")) {
            throw ValueError("Unsupported random filter");
        }
        bp::extract<unsigned int> ex1(filter[1]);
        if (!ex1.check()) {
            throw ValueError("Radius of Gaussian filter must be a positive integer.");
        }
        radius = ex1();
        sigma = 0.5;
        bp::extract<double> ex2(filter[2]);
        if (!ex2.check() || (sigma=ex2()) <= 0) {
            throw ValueError("Sigma must be a positive floating point number.");
        }
    } else {
        throw ValueError("Unsupported random filter for Rectangle.");
    }

    // Unfiltered random data (radius==0, or vector-valued): fill every node of a
    // ContinuousFunction Data directly. The ripley-style m_NN smoothing grid below
    // does not apply to the p4est node layout (m_NN is a stale ripley member for
    // oxley -> the copy loop wrote out of bounds and corrupted the heap under MPI,
    // and returned all zeros in serial). (A6.)
    if (radius == 0 || numvals > 1) {
        escript::FunctionSpace fs(getPtr(), getContinuousFunctionCode());
        escript::Data resdat(0, shape, fs, true);
        escript::DataTypes::RealVectorType& dv = resdat.getExpandedVectorReference();
        escript::randomFillArray(seed, &dv[0], dv.size(), m_mpiInfo);
        return resdat;
    }

    // number of points in the internal region
    // that is, the ones we need smoothed versions of
    const dim_t internal[2] = { m_NN[0], m_NN[1] };
    size_t ext[2];
    ext[0]=(size_t)internal[0]+2*radius; // includes points we need as input
    ext[1]=(size_t)internal[1]+2*radius; // for smoothing

    // now we check to see if the radius is acceptable
    // That is, would not cross multiple ranks in MPI

    if (2*radius >= internal[0]-4) {
        throw ValueError("Radius of gaussian filter is too large for X dimension of a rank");
    }
    if (2*radius >= internal[1]-4) {
        throw ValueError("Radius of gaussian filter is too large for Y dimension of a rank");
    }

    double* src = new double[ext[0]*ext[1]*numvals];
    escript::randomFillArray(seed, src, ext[0]*ext[1]*numvals, m_mpiInfo);

    //TODO
// #ifdef ESYS_MPI
//     if ((internal[0] < 5) || (internal[1] < 5)) {
//         // since the dimensions are equal for all ranks, this exception
//         // will be thrown on all ranks
//         throw OxleyException("Random Data in Oxley requires at least five elements per side per rank.");
//     }
//     dim_t X = m_mpiInfo->rank%m_NX[0];
//     dim_t Y = m_mpiInfo->rank/m_NX[0];
// #endif

// #ifdef ESYS_MPI
//     BlockGrid2 grid(m_NX[0]-1, m_NX[1]-1);
//     // it's +2 not +1 because a whole element is shared (and hence there is
//     // an overlap of two points both of which need to have "radius" points on
//     // either side.
//     size_t inset=2*radius+2;

//     // how wide is the x-dimension between the two insets
//     size_t xmidlen=ext[0]-2*inset;
//     size_t ymidlen=ext[1]-2*inset;

//     Block2 block(ext[0], ext[1], inset, xmidlen, ymidlen, numvals);

//     // a non-tight upper bound on how many we need
//     MPI_Request reqs[40];
//     MPI_Status stats[40];
//     short rused=0;

//     messvec incoms;
//     messvec outcoms;

//     grid.generateInNeighbours(X, Y, incoms);
//     grid.generateOutNeighbours(X, Y, outcoms);

//     block.copyAllToBuffer(src);

//     int comserr = 0;
//     for (size_t i=0; i < incoms.size(); ++i) {
//         message& m = incoms[i];
//         comserr |= MPI_Irecv(block.getInBuffer(m.destbuffid),
//                              block.getBuffSize(m.destbuffid), MPI_DOUBLE,
//                              m.sourceID, m.tag, m_mpiInfo->comm,
//                              reqs+(rused++));
//         block.setUsed(m.destbuffid);
//     }

//     for (size_t i=0; i < outcoms.size(); ++i) {
//         message& m = outcoms[i];
//         comserr |= MPI_Isend(block.getOutBuffer(m.srcbuffid),
//                              block.getBuffSize(m.srcbuffid), MPI_DOUBLE,
//                              m.destID, m.tag, m_mpiInfo->comm, reqs+(rused++));
//     }

//     if (!comserr) {
//         comserr = MPI_Waitall(rused, reqs, stats);
//     }

//     if (comserr) {
//         // Yes this is throwing an exception as a result of an MPI error
//         // and no we don't inform the other ranks that we are doing this.
//         // However, we have no reason to believe coms work at this point anyway
//         throw OxleyException("Error in coms for randomFill");
//     }

//     block.copyUsedFromBuffer(src);
// #endif

    // the truth of either should imply the truth of the other but let's be safe
    if (radius==0 || numvals > 1) {
        escript::FunctionSpace fs(getPtr(), getContinuousFunctionCode());
        escript::Data resdat(0, shape, fs, true);
        // don't need to check for exwrite because we just made it
        escript::DataTypes::RealVectorType& dv = resdat.getExpandedVectorReference();

        // now we need to copy values over
        for (size_t y=0; y < internal[1]; ++y) {
            for (size_t x=0; x < internal[0]; ++x) {
                for (unsigned int i=0; i < numvals; ++i) {
                    dv[i+(x+y*(internal[0]))*numvals]=src[i+(x+y*ext[0])*numvals];
                }
            }
        }
        delete[] src;
        return resdat;
    } else { // filter enabled
        escript::FunctionSpace fs(getPtr(), getContinuousFunctionCode());
        escript::Data resdat(0, escript::DataTypes::scalarShape, fs, true);
        // don't need to check for exwrite because we just made it
        escript::DataTypes::RealVectorType& dv=resdat.getExpandedVectorReference();
        double* convolution=get2DGauss(radius, sigma);
        for (size_t y=0; y < internal[1]; ++y) {
            for (size_t x=0; x < internal[0]; ++x) {
                dv[x+y*(internal[0])] = Convolve2D(convolution, src, x+radius, y+radius, radius, ext[0]);
            }
        }
        delete[] convolution;
        delete[] src;
        return resdat;
    }
}

dim_t Rectangle::findNode(const double *coords) const
{
    // Search the lnodes node coordinates for the closest node (used for Dirac
    // points). Replaces the coordinate-hash lookup.
    // reject points outside the domain bounding box (out-of-range Dirac points)
    const double x0=forestData.m_origin[0], y0=forestData.m_origin[1];
    const double x1=forestData.m_lxy[0],    y1=forestData.m_lxy[1];
    double ext = x1-x0; if(y1-y0>ext) ext=y1-y0;
    const double tol = 1e-8*ext;
    if(coords[0]<x0-tol || coords[0]>x1+tol || coords[1]<y0-tol || coords[1]>y1+tol)
        return -1;
    const MeshAccess m = getMeshAccess();
    long closest = 0;
    double best = std::numeric_limits<double>::max();
    for(long i = 0; i < m.numNodes; ++i)
    {
        const double dx = m.nodeCoords[(size_t) i*2 + 0] - coords[0];
        const double dy = m.nodeCoords[(size_t) i*2 + 1] - coords[1];
        const double d2 = dx*dx + dy*dy;
        if(d2 < best) { best = d2; closest = i; }
    }
    return (dim_t) closest;
}

const long Rectangle::getNodeId(double x, double y)
{
    const double coords[2] = {x, y};
    return (long) findNode(coords);
}

RankVector Rectangle::getOwnerVector(int fsType) const
{
    RankVector owner;
    const int rank = m_mpiInfo->rank;

    // Every local element is uniquely owned by this rank. The p4est SFC
    // partition splits leaves (and hence their face elements) disjointly, so
    // getNumElements()/getNumFaceElements() never include a halo -- the ghost
    // overlap used for assembly lives in m_ghost, not in the sample list. This
    // is the same fact that makes ownSample() return true for element types.
    // (Do NOT copy ripley's m_faceCount/m_NX logic here: ripley's element list
    // DOES carry a halo layer and its m_NX is the rank grid, whereas oxley's
    // m_NX is a physical element width -- that mismatch silently dropped cells
    // from weipa output under MPI.)
    if (fsType == Elements || fsType == ReducedElements) {
        owner.assign(getNumElements(), rank);
    } else if (fsType == FaceElements || fsType == ReducedFaceElements) {
        owner.assign(getNumFaceElements(), rank);
    } else {
        throw ValueError("getOwnerVector: only valid for element types");
    }

    return owner;
}

/**
 * \brief
 * Updates the mesh after refinement
*/
void Rectangle::updateMesh()
{
    oxleytimer.toc("\x1B[31mupdating the Mesh....\x1B[37m");
    // Update the nodes
    p4est_lnodes_destroy(nodes);
    p4est_ghost_t * ghost = p4est_ghost_new(p4est, P4EST_CONNECT_FULL);
    nodes = p4est_lnodes_new(p4est, ghost, 1);
    p4est_ghost_destroy(ghost);
    
    updateNodeIncrements();
    renumberNodes();
    updateRowsColumns();
    updateNodeDistribution();
    updateElementIds();
    updateFaceOffset();
    updateFaceElementCount();
    updateQuadrantIDinformation();
    oxleytimer.toc("\x1B[31mdone\x1B[37m");
}

/**
* \brief
* Toggles automatic mesh updates in the refinement functions
*/
void Rectangle::AutomaticMeshUpdateOnOff(bool new_setting)
{
    autoMeshUpdates = new_setting;    
}

/**
    \brief
    Applies a refinementzone
*/
escript::Domain_ptr Rectangle::apply_refinementzone(RefinementZone R)
{
    oxleytimer.toc("Applying the refinement zone...");

    oxley::Rectangle * newDomain = new Rectangle(*this, m_order);
    int numberOfRefinements = R.getNumberOfOperations();

    newDomain->AutomaticMeshUpdateOnOff(false);

    for(int n = 0; n < numberOfRefinements; n++)
    {
        RefinementType Refinement = R.getRefinement(n);
        //set the refinement level for this refinement
        newDomain->setRefinementLevels(Refinement.levels);
        switch(Refinement.flavour)
        {
            case POINT2D:
            {
                double x=Refinement.x0;
                double y=Refinement.y0;
                newDomain->refinePoint(x,y);
                break;
            }
            case REGION2D:
            {
                double x0=Refinement.x0;
                double y0=Refinement.y0;
                double x1=Refinement.x1;
                double y1=Refinement.y1;
                newDomain->refineRegion(x0,x1,y0,y1);
                break;
            }
            case CIRCLE:
            {
                double x0=Refinement.x0;
                double y0=Refinement.y0;
                double r0=Refinement.r;
                newDomain->refineCircle(x0,y0,r0);
                break;
            }
            case BOUNDARY:
            {
                double dx=Refinement.depth;
                switch(Refinement.b)
                {
                    case NORTH:
                    {
                        newDomain->refineBoundary("TOP",dx);
                        break;
                    }
                    case SOUTH:
                    {
                        newDomain->refineBoundary("BOTTOM",dx);
                        break;
                    }
                    case WEST:
                    {
                        newDomain->refineBoundary("LEFT",dx);
                        break;
                    }
                    case EAST:
                    {
                        newDomain->refineBoundary("RIGHT",dx);
                        break;
                    }
                    case TOP:
                    case BOTTOM:
                    default:
                    {
                        throw OxleyException("Invalid border direction.");
                    }
                }
                break;
            }
            case MASK2D:
            {
                if(n == 0)
                {
                    escript::Data d = *Refinement.data;
                    newDomain->refineMask(d);
                    break;
                }
                else
                {
                    throw OxleyException("Can only apply a mask refinement if it is first in the queue.");
                }
            }
            case MASK3D:
            case SPHERE:
            case POINT3D:
            case REGION3D:
            default:
                throw OxleyException("Unknown refinement algorithm.");
        }
    }

    newDomain->updateMesh();
    newDomain->AutomaticMeshUpdateOnOff(true);

    oxleytimer.toc("done");
    return escript::Domain_ptr(newDomain);
}

// interpolateWorker_Data::interpolateWorker_Data(const escript::Data * s, 
//                                                const escript::Data * t, 
//                                                const oxley::Rectangle * o)
// {
//    source = s;
//    target = t;
//    other = o;
// }

// interpolateWorker_Data::~interpolateWorker_Data()
// {

// }

// instantiate our two supported versions
template
void Rectangle::assembleGradientImpl<real_t>(escript::Data& out,
                                             const escript::Data& in) const;

template
void Rectangle::assembleGradientImpl<cplx_t>(escript::Data& out,
                                             const escript::Data& in) const;


} // end of namespace oxley
