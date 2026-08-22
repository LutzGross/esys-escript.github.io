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
#include <random>
#include <vector>

// #include <oxley/tic.h>

#include <escript/Assert.h>
#include <escript/Data.h>
#include <escript/DataFactory.h>
#include <escript/FunctionSpaceFactory.h>
#include <escript/Random.h>
#include <escript/Utils.h>

#include <oxley/AbstractAssembler.h>
#include <oxley/DefaultAssembler3D.h>
#include <oxley/InitAlgorithms.h>
#include <oxley/Oxley.h>
#include <oxley/OxleyData.h>
#include <oxley/Brick.h>
#include <oxley/MeshIO.h>
#include <oxley/RefinementAlgorithms.h>

// p8est headers will include MPI via sc.h when SC_ENABLE_MPI is defined
#include <p8est.h>
#include <p8est_algorithms.h>
#include <p8est_bits.h>
#include <p8est_connectivity.h>
#include <p8est_extended.h>
#include <p8est_ghost.h>
#include <p8est_iterate.h>
#include <p8est_lnodes.h>
#include <p8est_vtk.h>

#include <array>
#include <sc_containers.h>
#include <sc_mpi.h>

// Include after p8est to get MPI_COMM_WORLD

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
Brick::Brick(escript::JMPI jmpi, int order,
    dim_t n0, dim_t n1, dim_t n2,
    double x0, double y0, double z0,
    double x1, double y1, double z1,
    const std::vector<double>& points,
    const std::vector<int>& tags,
    const TagMap& tagnamestonums,
    const std::vector<int>& refine_level):
    OxleyDomain(3, order, jmpi)
{

    oxleytimer.toc("Creating an oxley::Brick...");

    // MPI communicator passed to base class constructor
    // Caller is responsible for ensuring MPI is initialized

    // n0/n1/n2 are the number of BLOCKS (p8est trees) per axis; refine_level is
    // the subdivision applied to each block. It is either a single value
    // (uniform) or one value per block (row-major over (n0,n1,n2), index
    // (i*n1+j)*n2+k), which lets blocks carry different levels and so introduces
    // hanging nodes at the block seams.
    if(n0 <= 0 || n1 <= 0 || n2 <= 0)
        throw OxleyException("Number of blocks in each spatial dimension must be positive");
    const dim_t num_trees = n0 * n1 * n2;
    if(refine_level.empty())
        throw OxleyException("refine_level must not be empty");
    if(refine_level.size() != 1 && (dim_t) refine_level.size() != num_trees)
        throw OxleyException("refine_level must be a single value or one value per block (n0*n1*n2)");
    for(size_t i = 0; i < refine_level.size(); i++)
        if(refine_level[i] < 0)
            throw OxleyException("refine_level must be non-negative");
    const int min_level = *std::min_element(refine_level.begin(), refine_level.end());
    const int max_level = *std::max_element(refine_level.begin(), refine_level.end());

    // Domain decomposition across MPI ranks is handled by p4est/p8est (see
    // p4est_partition below), not by a Cartesian d0 x d1 x d2 block grid.

    oxleytimer.toc("\t creating connectivity");
    connectivity = new_brick_connectivity(n0, n1, n2, false, false, false, x0, x1, y0, y1, z0, z1);

#ifdef OXLEY_ENABLE_TIMECONSUMING_DEBUG_CHECKS
    std::cout << "In Brick() constructor..." << std::endl;
    std::cout << "Checking connectivity ... " << std::endl;
    if(!p8est_connectivity_is_valid(connectivity))
    // if(!p8est_connectivity_is_valid_fast(connectivity))
        std::cout << "\t\tbroken" << std::endl;
    else
        std::cout << "\t\tOK" << std::endl;
#endif

    // create the forestdata
    forestData = new p8estData;

    // Create the p8est - use the custom communicator.
    // fill_uniform + min_level builds a uniform base mesh where every block is
    // subdivided min_level times. min_quadrants MUST be 0: it is p8est's
    // PER-PROCESSOR minimum, so a positive value forces extra refinement under
    // MPI and makes the mesh depend on the rank count. (A6.)
    p8est_locidx_t min_quadrants = 0;
    int fill_uniform = 1;
    oxleytimer.toc("\t creating p8est...");
    p8est = p8est_new_ext(m_mpiInfo->comm, connectivity, min_quadrants,
            min_level, fill_uniform, sizeof(octantData), &init_brick_data, (void *) forestData);

    // If the blocks do not all share the same level, refine each block up to its
    // own target level and 2:1-balance the base mesh. block_levels is indexed by
    // p8est tree id, so map each tree to its block (i,j,k) via its lower corner
    // connectivity vertex (the trees are Morton-ordered, not row-major).
    if(max_level > min_level) {
        const double dxb = (x1-x0)/n0, dyb = (y1-y0)/n1, dzb = (z1-z0)/n2;
        forestData->block_levels.assign(num_trees, min_level);
        for(p4est_topidx_t t = 0; t < num_trees; t++) {
            const p4est_topidx_t v = connectivity->tree_to_vertex[P8EST_CHILDREN*t + 0];
            const double vx = connectivity->vertices[3*v + 0];
            const double vy = connectivity->vertices[3*v + 1];
            const double vz = connectivity->vertices[3*v + 2];
            int bi = (int) std::lround((vx - x0)/dxb);
            int bj = (int) std::lround((vy - y0)/dyb);
            int bk = (int) std::lround((vz - z0)/dzb);
            if(bi < 0) bi = 0; else if(bi >= n0) bi = n0-1;
            if(bj < 0) bj = 0; else if(bj >= n1) bj = n1-1;
            if(bk < 0) bk = 0; else if(bk >= n2) bk = n2-1;
            forestData->block_levels[t] = refine_level[((size_t) bi*n1 + bj)*n2 + bk];
        }
        const int recursive = 1;
        p8est_refine_ext(p8est, recursive, max_level, refine_to_block_level,
                         init_brick_data, NULL);
        p8est_balance_ext(p8est, P8EST_CONNECT_FULL, init_brick_data, NULL);
    }

#ifdef OXLEY_ENABLE_DEBUG_CHECKS //These checks are turned off by default as they can be very timeconsuming
    std::cout << "Checking p8est ... ";
    if(!p8est_is_valid(p8est))
        std::cout << "broken" << std::endl;
    else
        std::cout << "OK" << std::endl;
#endif

    // Distribute the p8est across the processors. This MUST happen before the
    // lnodes are built: p8est_partition moves octants between ranks, and an lnodes
    // built beforehand still describes the old partition (its num_local_elements
    // and element_nodes refer to octants this rank no longer owns), so every later
    // leaf walk runs off the end of element_nodes. It went unnoticed while every
    // mesh was uniform, because then the partition is already balanced and
    // p8est_partition is a no-op. (Milestone B.)
    oxleytimer.toc("\t partitioning...");
    int allow_coarsening = 0;
    p8est_partition(p8est, allow_coarsening, NULL);

    // Nodes numbering
    oxleytimer.toc("\t creating ghost...");
    ghost = p8est_ghost_new(p8est, P8EST_CONNECT_FULL);
    nodes = p8est_lnodes_new(p8est, ghost, 1);

    oxleytimer.toc("\t saving forestData information...");
    // This information is needed by the assembler
    m_NE[0] = n0;
    m_NE[1] = n1;
    m_NE[2] = n2;
    m_NX[0] = (x1-x0)/n0;
    m_NX[1] = (y1-y0)/n1;
    m_NX[2] = (z1-z0)/n2;

    // Record the physical dimensions of the domain and the location of the origin
    m_blocks[0] = n0;
    m_blocks[1] = n1;
    m_blocks[2] = n2;

    forestData->m_origin[0] = x0;
    forestData->m_origin[1] = y0;
    forestData->m_origin[2] = z0;
    forestData->m_lxyz[0] = x1;
    forestData->m_lxyz[1] = y1;
    forestData->m_lxyz[2] = z1;
    forestData->m_length[0] = x1-x0;
    forestData->m_length[1] = y1-y0;
    forestData->m_length[2] = z1-z0;

    // Periodic boundaries are not implemented: the connectivity is built
    // non-periodic. forestData->periodic stays at its default (false).

    double multiplier =  sqrt(x0*x0+y0*y0+z0*z0);
    if(multiplier == 0.0)
        multiplier = 1.0;

    tuple_tolerance = TOLERANCE * multiplier;

    // Find the grid spacing for each level of refinement in the mesh
#pragma omp parallel for
    for(int i = 0; i <= P8EST_MAXLEVEL; i++){
        double numberOfSubDivisions = (p8est_qcoord_t) (1 << (P8EST_MAXLEVEL - i));
        forestData->m_dx[0][i] = forestData->m_length[0] / numberOfSubDivisions;
        forestData->m_dx[1][i] = forestData->m_length[1] / numberOfSubDivisions;
        forestData->m_dx[2][i] = forestData->m_length[2] / numberOfSubDivisions;
    }

    // max levels of refinement
    forestData->max_levels_refinement = MAXREFINEMENTLEVELS;

    // element order
    m_order = order;

    // Number of dimensions
    m_numDim=3;

    // (p8est_partition happens above, before the lnodes are built)

    // Indices
    indices = new std::vector<IndexVector>;


    // Number the nodes
    updateMesh();

    // Tags
    oxleytimer.toc("\t populating sample ids...");
    // TODO
    // populateSampleIds(); // very slow

    oxleytimer.toc("\t populating tags ("+std::to_string(tagnamestonums.size())+")...");
    for (TagMap::const_iterator i = tagnamestonums.begin(); i != tagnamestonums.end(); i++) {
        setTagMap(i->first, i->second);
    }

    // Dirac points and tags
    oxleytimer.toc("\t adding Dirac points...");
    addPoints(points, tags);
    oxleytimer.toc("\t\tdone");

    // To prevent segmentation faults when using numpy ndarray
#ifdef ESYS_HAVE_BOOST_NUMPY
    oxleytimer.toc("\tpy_initialise");
    Py_Initialize();
    oxleytimer.toc("\tboost numpy initialise");
    boost::python::numpy::initialize(); // Not MPI safe?
#endif

#ifdef ESYS_HAVE_PASO

    /// local array length shared
    // dim_t local_length = p8est->local_num_quadrants;
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

    oxleytimer.toc("Brick initialised");
}

Brick::Brick(oxley::Brick& B, int order, bool update):
    OxleyDomain(3, order)
{
    oxleytimer.toc("In Brick copy constructor");

    oxleytimer.setTime(B.oxleytimer.getTime());

    // This information is needed by the assembler
    m_NE[0] = B.m_NE[0];
    m_NE[1] = B.m_NE[1];
    m_NE[2] = B.m_NE[2];
    m_NX[0] = B.m_NX[0];
    m_NX[1] = B.m_NX[1];
    m_NX[2] = B.m_NX[2];

    // create the forestdata
    forestData = new p8estData;
    m_blocks[0] = B.m_blocks[0];
    m_blocks[1] = B.m_blocks[1];
    m_blocks[2] = B.m_blocks[2];

    forestData->m_origin[0] = B.forestData->m_origin[0];
    forestData->m_origin[1] = B.forestData->m_origin[1];
    forestData->m_origin[2] = B.forestData->m_origin[2];
    forestData->m_lxyz[0] = B.forestData->m_lxyz[0];
    forestData->m_lxyz[1] = B.forestData->m_lxyz[1];
    forestData->m_lxyz[2] = B.forestData->m_lxyz[2];
    forestData->m_length[0] = B.forestData->m_length[0];
    forestData->m_length[1] = B.forestData->m_length[1];
    forestData->m_length[2] = B.forestData->m_length[2];
    forestData->m_NX[0] = B.forestData->m_NX[0];
    forestData->m_NX[1] = B.forestData->m_NX[1];
    forestData->m_NX[2] = B.forestData->m_NX[2];
    forestData->m_gNE[0] = B.forestData->m_gNE[0];
    forestData->m_gNE[1] = B.forestData->m_gNE[1];
    forestData->m_gNE[2] = B.forestData->m_gNE[2];
    forestData->m_NE[0] = B.forestData->m_NE[0];
    forestData->m_NE[1] = B.forestData->m_NE[1];
    forestData->m_NE[2] = B.forestData->m_NE[2];
    forestData->periodic[0] = B.forestData->periodic[0];
    forestData->periodic[1] = B.forestData->periodic[1];
    forestData->periodic[2] = B.forestData->periodic[2];
    forestData->max_levels_refinement = B.forestData->max_levels_refinement;
    forestData->refinement_depth = B.forestData->refinement_depth;
    forestData->refinement_boundaries[0] = B.forestData->refinement_boundaries[0];
    forestData->refinement_boundaries[1] = B.forestData->refinement_boundaries[1];
    forestData->refinement_boundaries[2] = B.forestData->refinement_boundaries[2];
    forestData->refinement_boundaries[3] = B.forestData->refinement_boundaries[3];
    forestData->refinement_boundaries[4] = B.forestData->refinement_boundaries[4];
    forestData->refinement_boundaries[5] = B.forestData->refinement_boundaries[5];

    double multiplier = sqrt(forestData->m_length[0]*forestData->m_length[0]
                           + forestData->m_length[1]*forestData->m_length[1] 
                           + forestData->m_length[2]*forestData->m_length[2]);
    if(multiplier == 0.0)
        multiplier = 1.0;

    tuple_tolerance = TOLERANCE * multiplier;

    m_mpiInfo=B.m_mpiInfo;

    oxleytimer.toc("\t Creating connectivity...");
    // p8est_connectivity_t tempConnectivity(*B.connectivity);
    connectivity = new_brick_connectivity(m_NE[0], m_NE[1], m_NE[2], 
                                          false, false, false, 
                                          forestData->m_origin[0], forestData->m_lxyz[0], 
                                          forestData->m_origin[1], forestData->m_lxyz[1], 
                                          forestData->m_origin[2], forestData->m_lxyz[2]);    

    p8est_locidx_t min_quadrants = 0;  // per-processor min; keep 0 for MPI mesh consistency (A6)
    int min_level = 0;
    int fill_uniform = 1;
    oxleytimer.toc("\t creating p8est...");
    p8est = p8est_new_ext(m_mpiInfo->comm, connectivity, min_quadrants,
        min_level, fill_uniform, sizeof(octantData), &init_brick_data, (void *) forestData);

//These checks are turned off by default as they can be very timeconsuming
#ifdef OXLEY_ENABLE_TIMECONSUMING_DEBUG_CHECKS
    std::cout << "In Brick copy constructor..." << std::endl;
    std::cout << "Checking connectivity ... " << std::endl;;
    // if(!p8est_connectivity_is_valid_fast(connectivity))
    if(!p8est_connectivity_is_valid(connectivity))
        std::cout << "\t broken" << std::endl;
    else
        std::cout << "\t OK" << std::endl;
    std::cout << "Checking p8est ... " << std::endl;
    if(!p8est_is_valid(p8est))
        std::cout << "\t broken" << std::endl;
    else
        std::cout << "\t OK" << std::endl;
#endif
#ifdef OXLEY_COPY_CONSTRUCTOR
    std::cout << "\t\t Checking p8est_is_equal: ";
    if(p8est_is_equal(p8est, B.p8est, false))
        std::cout << "OK" << std::endl;
    else
        std::cout << "broken" << std::endl;
    std::cout << "\t\t Checking p8est_connectivity_is_equal: ";
    if(p8est_connectivity_is_equal(connectivity, B.connectivity))
        std::cout << "OK" << std::endl;
    else
        std::cout << "broken" << std::endl;
#endif

    // Nodes numbering
    oxleytimer.toc("\t creating ghost...");
    ghost = p8est_ghost_new(p8est, P8EST_CONNECT_FULL);
    oxleytimer.toc("\t creating lnodes...");
    nodes = p8est_lnodes_new(p8est, ghost, 1);
    
    oxleytimer.toc("\t populating forestData...");
    // Find the grid spacing for each level of refinement in the mesh
#pragma omp parallel for
    for(int i = 0; i <= P8EST_MAXLEVEL; i++)
    {
        double numberOfSubDivisions = (p8est_qcoord_t) (1 << (P8EST_MAXLEVEL - i));
        forestData->m_dx[0][i] = forestData->m_length[0] / numberOfSubDivisions;
        forestData->m_dx[1][i] = forestData->m_length[1] / numberOfSubDivisions;
        forestData->m_dx[2][i] = forestData->m_length[2] / numberOfSubDivisions;
    }

    // element order
    m_order = order;

    // Number of dimensions
    m_numDim=3;

    // Indices
    oxleytimer.toc("\t creating index vector");
    indices = new std::vector<IndexVector>;

    // Number the nodes, if required
    if(update)
    {
        oxleytimer.toc("\t updating the mesh (inc partitioning) ");
        updateMesh();
    }

    // Tags
    oxleytimer.toc("\t populating sample ids");
    // TODO
    // populateSampleIds();
    
    oxleytimer.toc("\t copying tag map");
    TagMap tempTagMap(B.m_tagMap);
    m_tagMap=tempTagMap;

    // Dirac points and tags
    oxleytimer.toc("\t copying Dirac points");
    std::vector<DiracPoint> diracTemp(B.m_diracPoints);
    m_diracPoints=diracTemp;

#ifdef ESYS_HAVE_BOOST_NUMPY
    Py_Initialize();
    boost::python::numpy::initialize();
#endif


    updateMesh();

    oxleytimer.toc("Brick initialised");
}

/**
   \brief
   Destructor.
*/
Brick::~Brick(){
#ifdef OXLEY_ENABLE_DEBUG_CHECKS
    std::cout << "\033[1;34m[oxley]\033[0m In Brick() destructor" << std::endl;
    std::cout << "\033[1;34m[oxley]\033[0m checking p8est ... " << std::endl;
    if(!p8est_is_valid(p8est))
        std::cout << "\t\tbroken" << std::endl;
    else
        std::cout << "\t\tOK" << std::endl;
    #ifdef OXLEY_ENABLE_TIMECONSUMING_DEBUG_CHECKS
    std::cout << "\033[1;34m[oxley]\033[0m checking connectivity ... " << std::endl;
    if(!p8est_connectivity_is_valid(connectivity))
        std::cout << "\t\tbroken" << std::endl;
    else
        std::cout << "\t\tOK" << std::endl;
    #endif
    std::cout << "\033[1;34m[oxley]\033[0m checking ghost ... " << std::endl;
    if(!p8est_ghost_is_valid(p8est,ghost))
        std::cout << "\t\tbroken" << std::endl;
    else
        std::cout << "\t\tOK" << std::endl;
#endif

    delete forestData;
}

/**
   \brief
   returns a description for this domain
*/
std::string Brick::getDescription() const
{
    return "oxley::brick";
}

/**
   \brief
   writes the current mesh to a file with the given name
   \param filename The name of the file to write to
*/
void Brick::write(const std::string& filename) const
{
    throw OxleyException("write: not supported");
}

void Brick::interpolateAcross(escript::Data& target, const escript::Data& source) const
{
    const Brick *other = dynamic_cast<const Brick *>(target.getDomain().get());
    if (other == NULL)
        throw OxleyException("Invalid interpolation: Domains must both be instances of Brick");
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
    if (other->getRefinementLevels() > getRefinementLevels()) {
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
                }
                break;
            case Elements:
                switch (fsTarget) {
                    case Elements:
                        interpolateElementsToElementsFiner(source, target, *other);
                        return;
                }
                break;
            case ReducedElements:
                switch (fsTarget) {
                    case Elements:
                        interpolateReducedToElementsFiner(source, target, *other);
                        return;
                }
                break;
        }
        msg << " when target is a finer mesh";
    } else {
        switch (fsSource) {
            case Nodes:
                switch (fsTarget) {
                    case Elements:
                        escript::Data elements=escript::Vector(0., escript::function(*this), true);
                        interpolateNodesOnElements(elements, source, false);
                        interpolateElementsToElementsCoarser(elements, target, *other);
                        return;
                }
                break;
            case Elements:
                switch (fsTarget) {
                    case Elements:
                        interpolateElementsToElementsCoarser(source, target, *other);
                        return;
                }
                break;
        }
        msg << " when target is a coarser mesh";
    }
    throw OxleyException(msg.str());
}

void Brick::validateInterpolationAcross(int fsType_source, const escript::AbstractDomain& domain, int fsType_target) const
{
    const Brick *other = dynamic_cast<const Brick *>(&domain);
    if (other == NULL)
        throw OxleyException("Invalid interpolation: domains must both be instances of oxley::Brick");

    // TODO
    // if(!p8est_is_equal(borrow_p8est, other->borrow_p8est(), 0))
    //     throw OxleyException("Invalid interpolation: domains have different p8ests");

    // if(!p8est_connectivity_is_equivalent(borrow_connectivity(),other->borrow_connectivity()))
    //     throw OxleyException("Invalid interpolation: domains have different connectivities");


}

void Brick::interpolateNodesToNodesFiner(const escript::Data& source, escript::Data& target, const Brick& other) const
{

}

void Brick::interpolateNodesToElementsFiner(const escript::Data& source, escript::Data& target, const Brick& other)  const
{

}

void Brick::interpolateElementsToElementsCoarser(const escript::Data& source, escript::Data& target, const Brick& other)  const
{

}

void Brick::interpolateElementsToElementsFiner(const escript::Data& source, escript::Data& target, const Brick& other)  const
{

}

void Brick::interpolateReducedToElementsFiner(const escript::Data& source, escript::Data& target, const Brick& other)  const
{

}

void Brick::interpolateReducedToReducedFiner(const escript::Data& source, escript::Data& target, const Brick& other)  const
{

}

void Brick::setToNormal(escript::Data& out) const
{
    if (out.getFunctionSpace().getTypeCode() == FaceElements) {
        out.requireWrite();
#pragma omp parallel
        {
            if(m_faceOffset[0] > -1) {
#pragma omp for nowait
                for (index_t k=0; k<NodeIDsLeft.size(); k++) {
                    double* o = out.getSampleDataRW(m_faceOffset[0]+k);
                    // set vector at four quadrature points
                    *o++ = -1.; *o++ = 0.; *o++ = 0.;
                    *o++ = -1.; *o++ = 0.; *o++ = 0.;
                    *o++ = -1.; *o++ = 0.; *o++ = 0.;
                    *o++ = -1.; *o++ = 0.; *o = 0.;
                }
            }
            
            if(m_faceOffset[1] > -1) {
#pragma omp for nowait
                for (index_t k=0; k<NodeIDsRight.size(); k++) {
                    double* o = out.getSampleDataRW(m_faceOffset[1]+k);
                    // set vector at four quadrature points
                    *o++ = 1.; *o++ = 0.; *o++ = 0.;
                    *o++ = 1.; *o++ = 0.; *o++ = 0.;
                    *o++ = 1.; *o++ = 0.; *o++ = 0.;
                    *o++ = 1.; *o++ = 0.; *o = 0.;
                }
            }

            if (m_faceOffset[2] > -1) {
#pragma omp for nowait
                for (index_t k=0; k<NodeIDsBottom.size(); k++) {
                    double* o = out.getSampleDataRW(m_faceOffset[2]+k);                    
                    // set vector at four quadrature points
                    *o++ = 0.; *o++ = -1.; *o++ = 0.;
                    *o++ = 0.; *o++ = -1.; *o++ = 0.;
                    *o++ = 0.; *o++ = -1.; *o++ = 0.;
                    *o++ = 0.; *o++ = -1.; *o = 0.;
                }
            }

            if (m_faceOffset[3] > -1) {
#pragma omp for nowait
                for (index_t k=0; k<NodeIDsTop.size(); k++) {
                    double* o = out.getSampleDataRW(m_faceOffset[3]+k);
                    // set vector at four quadrature points
                    *o++ = 0.; *o++ = 1.; *o++ = 0.;
                    *o++ = 0.; *o++ = 1.; *o++ = 0.;
                    *o++ = 0.; *o++ = 1.; *o++ = 0.;
                    *o++ = 0.; *o++ = 1.; *o = 0.;
                }
            }

            if(m_faceOffset[4] > -1) {
#pragma omp for nowait
                for (index_t k=0; k<NodeIDsAbove.size(); k++) {
                    double* o = out.getSampleDataRW(m_faceOffset[4]+k);
                    // set vector at four quadrature points
                    *o++ = 0.; *o++ = 0.; *o++ = -1.;
                    *o++ = 0.; *o++ = 0.; *o++ = -1.;
                    *o++ = 0.; *o++ = 0.; *o++ = -1.;
                    *o++ = 0.; *o++ = 0.; *o = -1.;
                }
            }

            if (m_faceOffset[5] > -1) {
#pragma omp for nowait
                for (index_t k=0; k<NodeIDsBelow.size(); k++) {
                    double* o = out.getSampleDataRW(m_faceOffset[5]+k);
                    // set vector at four quadrature points
                    *o++ = 0.; *o++ = 0.; *o++ = 1.;
                    *o++ = 0.; *o++ = 0.; *o++ = 1.;
                    *o++ = 0.; *o++ = 0.; *o++ = 1.;
                    *o++ = 0.; *o++ = 0.; *o = 1.;
                }
            }
        } // end of parallel
    }
    else if (out.getFunctionSpace().getTypeCode() == ReducedFaceElements) 
    {
        out.requireWrite();
        
#pragma omp parallel
        {
            if (m_faceOffset[0] > -1) {
#pragma omp for nowait
                for (index_t k=0; k<NodeIDsLeft.size(); k++) {
                    double* o = out.getSampleDataRW(m_faceOffset[0]+k);
                    *o++ = -1.;
                    *o++ = 0.;
                    *o = 0.;
                }
            }

            if (m_faceOffset[1] > -1) {
#pragma omp for nowait
                for (index_t k=0; k<NodeIDsRight.size(); k++) {
                    double* o = out.getSampleDataRW(m_faceOffset[1]+k);
                    *o++ = 1.;
                    *o++ = 0.;
                    *o = 0.;
                }
            }

            if (m_faceOffset[2] > -1) {
#pragma omp for nowait
                for (index_t k=0; k<NodeIDsBottom.size(); k++) {
                    double* o = out.getSampleDataRW(m_faceOffset[2]+k);
                    *o++ = 0.;
                    *o++ = -1.;
                    *o = 0.;
                }
            }

            if (m_faceOffset[3] > -1) {
#pragma omp for nowait
                for (index_t k=0; k<NodeIDsTop.size(); k++) {
                    double* o = out.getSampleDataRW(m_faceOffset[3]+k);
                    *o++ = 0.;
                    *o++ = 1.;
                    *o = 0.;
                }
            }

            if (m_faceOffset[4] > -1) {
#pragma omp for nowait
                for (index_t k=0; k<NodeIDsAbove.size(); k++) {
                    double* o = out.getSampleDataRW(m_faceOffset[4]+k);
                    *o++ = 0.;
                    *o++ = 0.;
                    *o = -1.;
                }
            }

            if (m_faceOffset[5] > -1) {
#pragma omp for nowait
                for (index_t k=0; k<NodeIDsBelow.size(); k++) {
                    double* o = out.getSampleDataRW(m_faceOffset[5]+k);
                    *o++ = 0.;
                    *o++ = 0.;
                    *o = 1.;
                }
            }
        }
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

void Brick::setToSize(escript::Data& out) const
{
    if (out.getFunctionSpace().getTypeCode() == Elements
            || out.getFunctionSpace().getTypeCode() == ReducedElements) 
    {
        out.requireWrite();

        // Find the maximum level of refinement in the mesh
        int max_level = 0;
        for(p8est_topidx_t tree = p8est->first_local_tree; tree <= p8est->last_local_tree; tree++) {
            p8est_tree_t * tree_t = p8est_tree_array_index(p8est->trees, tree);
            max_level = tree_t->maxlevel > max_level ? tree_t->maxlevel : max_level;
        }
        // Work out the size at each level
        std::vector<double> size_vect(max_level+1, -1.0);
        for(int i = 0 ; i <= max_level ; i++)
        {
            size_vect[i] = sqrt((forestData->m_dx[0][P8EST_MAXLEVEL-i]*forestData->m_dx[0][P8EST_MAXLEVEL-i]
                                                    +forestData->m_dx[1][P8EST_MAXLEVEL-i]*forestData->m_dx[1][P8EST_MAXLEVEL-i]));
        }

        const dim_t numQuad=out.getNumDataPointsPerSample();
        long id = 0;   // running local leaf index (element sample order)
        for(p8est_topidx_t t = p8est->first_local_tree; t <= p8est->last_local_tree; t++)
        {
            p8est_tree_t * currenttree = p8est_tree_array_index(p8est->trees, t);
            sc_array_t * tquadrants = &currenttree->quadrants;
            p8est_qcoord_t Q = (p8est_qcoord_t) tquadrants->elem_count;
            for (int q = 0; q < Q; ++q, ++id)
            {
                p8est_quadrant_t * quad = p8est_quadrant_array_index(tquadrants, q);
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
        // characteristic face size = sqrt(face area); in-plane sizes are
        // m_NX/2^level (per-block length / 2^level).
        const std::vector<borderNodeInfo>* lists[6] = {
            &NodeIDsLeft, &NodeIDsRight, &NodeIDsBottom, &NodeIDsTop,
            &NodeIDsAbove, &NodeIDsBelow };
        static const int planeAxes[6][2] = {{1,2},{1,2},{0,2},{0,2},{0,1},{0,1}};
        for(int fc = 0; fc < 6; ++fc) {
            if(m_faceOffset[fc] < 0) continue;
            const std::vector<borderNodeInfo>& L = *lists[fc];
            for (index_t k=0; k<L.size(); k++) {
                const double h = (double)(1 << L[k].level);
                const double sz = std::sqrt((m_NX[planeAxes[fc][0]]/h) * (m_NX[planeAxes[fc][1]]/h));
                double* o = out.getSampleDataRW(m_faceOffset[fc]+k);
                std::fill(o, o+numQuad, sz);
            }
        }
    } 
    else 
    {
        std::stringstream msg;
        msg << "setToSize: invalid function space type "
            << out.getFunctionSpace().getTypeCode();
        throw ValueError(msg.str());
    }
}

bool Brick::ownSample(int fsType, index_t id) const
{
    if (getMPISize()==1)
        return true;

    switch (fsType) {
        case Nodes:
        case ReducedNodes: //FIXME: reduced
            return (m_dofMap[id] < getNumDOF());
        case DegreesOfFreedom:
        case ReducedDegreesOfFreedom:
            return true;
        case Elements:
        case ReducedElements:
        case FaceElements:
        case ReducedFaceElements:
        case Points:
            // p8est partitions leaves (and hence boundary faces) uniquely across
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

dim_t Brick::getNumDataPointsGlobal() const
{
    // total number of (owned) nodes across all ranks
    if(!nodes) return 0;
    dim_t total = 0;
    for(int r = 0; r < m_mpiInfo->size; ++r)
        total += (dim_t) nodes->global_owned_count[r];
    return total;
}


/* This is a wrapper for filtered (and non-filtered) randoms
 * For detailed doco see randomFillWorker
*/
escript::Data Brick::randomFill(const escript::DataTypes::ShapeType& shape,
                                    const escript::FunctionSpace& fs,
                                    long seed, const bp::tuple& filter) const
{
    throw OxleyException("randomFill"); //TODO this is temporary
}

void Brick::dump(const std::string& fileName) const
{
    // The bespoke Silo dump was demoted in favour of the standard weipa output
    // path (saveSilo / saveVTK), which consumes the lnodes mesh-access interface
    // and works in parallel. (design A4.8)
    throw OxleyException("Brick::dump: use saveSilo()/saveVTK() (weipa) instead");
}

const dim_t* Brick::borrowSampleReferenceIDs(int fsType) const
{
    switch (fsType) {
        case Nodes:
        case ReducedNodes: //FIXME: reduced
            // node reference ids live in m_nodeId (populated by renumberNodes);
            // myColumns is the Yale-format matrix column array, not node ids.
            return &m_nodeId[0];
        case DegreesOfFreedom:
        case ReducedDegreesOfFreedom: //FIXME: reduced
            // conforming serial: DOF id == node id
            return &m_nodeId[0];
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
            msg << "borrowSampleReferenceIDs: invalid function space type "<<fsType;
            throw ValueError(msg.str());
    }
}

void Brick::writeToVTK(std::string filename, bool writeMesh) const
{
    // Check that the filename is reasonable
    int stringlength=filename.length();
    if(stringlength>4){
        if(filename.compare(stringlength-4,4,".vtk") != 0)
            filename.append(".vtk");
    } else {
        filename.append(".vtk");
    }

    // Write to file
    const char * name = filename.c_str();
    if(writeMesh)
    {
        p8est_vtk_write_file(p8est, NULL, name);
    }
    else
    {
        // Create the context for the VTK file
        p8est_vtk_context_t * context = p8est_vtk_context_new(p8est, name);

        // Continuous point data
        p8est_vtk_context_set_continuous(context, true);

        // Set the scale
        p8est_vtk_context_set_scale(context, 1.0);

        // Write the header
        context = p8est_vtk_write_header(context);

        // Get the point and cell data together
        p8est_locidx_t numquads = p8est->local_num_quadrants;
        sc_array_t * quadTag = sc_array_new_count(sizeof(double), numquads);
        p8est_iterate(p8est, NULL, (void *) quadTag, getQuadTagVector, NULL, NULL, NULL);
        sc_array_t * xcoord = sc_array_new_count(sizeof(double), numquads);
        p8est_iterate(p8est, NULL, (void *) xcoord, getXCoordVector, NULL, NULL, NULL);
        sc_array_t * ycoord = sc_array_new_count(sizeof(double), numquads);
        p8est_iterate(p8est, NULL, (void *) ycoord, getYCoordVector, NULL, NULL, NULL);
        sc_array_t * zcoord = sc_array_new_count(sizeof(double), numquads);
        p8est_iterate(p8est, NULL, (void *) zcoord, getZCoordVector, NULL, NULL, NULL);

        // Cell Data
#ifdef OXLEY_ENABLE_DEBUG
        context = p8est_vtk_write_cell_dataf(context,1,1,0,0,4,0, "tag",quadTag,"x",xcoord,"y",ycoord,"z",zcoord,context);
#else
        context = p8est_vtk_write_cell_dataf(context,0,0,0,0,4,0, "tag",quadTag,"x",xcoord,"y",ycoord,"z",zcoord,context);
#endif
        if(context == NULL)
            throw OxleyException("Error writing cell data");

        // Point Data
        context = p8est_vtk_write_point_dataf(context, 0, 0, context);
        if(context == NULL)
            throw OxleyException("Error writing point data");

        // Write the footer
        if(p8est_vtk_write_footer(context))
                throw OxleyException("Error writing footer.");

        // Cleanup
                // Cleanup
        sc_array_reset(quadTag);
        sc_array_destroy(quadTag);
        sc_array_reset(xcoord);
        sc_array_destroy(xcoord);
        sc_array_reset(ycoord);
        sc_array_destroy(ycoord);
        sc_array_reset(zcoord);
        sc_array_destroy(zcoord);
    }
}

void Brick::saveMesh(std::string filename) 
{
    // see Rectangle::saveMesh: p4est does not record the geometry, so write it
    // alongside and let loadMesh() build the domain itself
    MeshHeader header;
    header.dim = 3;
    header.order = m_order;
    for(int d = 0; d < 3; ++d) {
        header.n[d] = m_blocks[d];
        header.origin[d] = forestData->m_origin[d];
        header.extent[d] = forestData->m_lxyz[d];
    }
    if(m_mpiInfo->rank == 0)
        writeMeshHeader(filename, header);

    std::string fnames=filename+".p8est";
    std::string cnames=filename+".conn";

    const char * fname=fnames.c_str();
    const char * cname=cnames.c_str();

    p8est_deflate_quadrants(p8est, NULL);

// #ifdef ESYS_MPI
//     if(escript::getMPIRankWorld()==0)
//     {
// #endif
        int retval = p8est_connectivity_save(cname, connectivity)==0;
        ESYS_ASSERT(retval!=0,"Failed to save connectivity");
        int save_partition = 0;
        int save_data = 1;
        p8est_save_ext(fname, p8est, save_data, save_partition); // Should abort on file error
// #ifdef ESYS_MPI
//     }
// #endif
}

// updateMesh rebuilds everything derived from the forest - balance, partition,
// ghost, lnodes, then the node/element ids and face counts. It uses only p4est
// and oxley's own bookkeeping: the trilinos guard it sat under was left over
// from when updateRowsColumns built a Tpetra graph, and it made the domain
// unbuildable without trilinos, since the constructors call this.
void Brick::updateMesh()
{
    // Mesh creation
    oxleytimer.toc("Brick::updateMesh ...");
    oxleytimer.toc("\t balancing...");    
    p8est_balance_ext(p8est, P8EST_CONNECT_FULL, init_brick_data, refine_copy_parent_octant);
    #ifdef ESYS_MPI
    MPI_Barrier(m_mpiInfo->comm);
    #endif

    oxleytimer.toc("\t partitioning...");    
    bool partition_for_coarsening = true;
    p8est_partition_ext(p8est, partition_for_coarsening, NULL);
    #ifdef ESYS_MPI
    MPI_Barrier(m_mpiInfo->comm);
    #endif

    oxleytimer.toc("\t destroying old ghost..."); 
    p8est_ghost_destroy(ghost);
    #ifdef ESYS_MPI
    MPI_Barrier(m_mpiInfo->comm);
    #endif

    oxleytimer.toc("\t creating new ghost..."); 
    ghost = p8est_ghost_new(p8est, P8EST_CONNECT_FULL);
    #ifdef ESYS_MPI
    MPI_Barrier(m_mpiInfo->comm);
    #endif

    oxleytimer.toc("\t destroying old lnodes...");   
    p8est_lnodes_destroy(nodes);
    #ifdef ESYS_MPI
    MPI_Barrier(m_mpiInfo->comm);
    #endif

    oxleytimer.toc("\t creating new lnodes..."); 
    nodes = p8est_lnodes_new(p8est, ghost, 1);
    #ifdef ESYS_MPI
    MPI_Barrier(m_mpiInfo->comm);
    #endif


    // addition information
    oxleytimer.toc("\t renumbering nodes");
    renumberNodes();
    updateRowsColumns();
    updateElementIds();
    updateFaceElementCount();
    // oxleytimer.toc("Brick::updateMesh: Finished...");

}

void Brick::updateMeshBackend()
{
    // Mesh creation
    oxleytimer.toc("Brick::updateMeshBackend ...");
    oxleytimer.toc("\t balancing...");    
    p8est_balance_ext(p8est, P8EST_CONNECT_FULL, init_brick_data, refine_copy_parent_octant);
    oxleytimer.toc("\t partitioning...");    
    bool partition_for_coarsening = true;
    p8est_partition_ext(p8est, partition_for_coarsening, NULL);
    oxleytimer.toc("\t destroying old ghost..."); 
    p8est_ghost_destroy(ghost);
    oxleytimer.toc("\t creating new ghost..."); 
    ghost = p8est_ghost_new(p8est, P8EST_CONNECT_FULL);
    oxleytimer.toc("\t destroying old lnodes...");   
    p8est_lnodes_destroy(nodes);
    oxleytimer.toc("\t creating new lnodes..."); 
    nodes = p8est_lnodes_new(p8est, ghost, 1);
}

void Brick::AutomaticMeshUpdateOnOff(bool new_setting)
{
    #ifdef OXLEY_ENABLE_PROFILE_TIMERS_INFORMATIONAL
    oxleytimer.toc("INFO: Disabling automatic mesh update");
    #endif
    autoMeshUpdates = new_setting;
}

// loadMesh uses only p4est; it sat inside the trilinos guard for no
// recorded reason, which left a build without trilinos unable to read a
// mesh at all.
void Brick::loadMesh(std::string filename) 
{
    std::string fnames=filename+".p8est";
    std::string cnames=filename+".conn";

    const char * fname=fnames.c_str();
    const char * cname=cnames.c_str();

    int load_data=true;
    int autopartition=true;
    int broadcasthead=false;

    // Delete the old structure
    p8est_connectivity_destroy(connectivity);
    p8est_destroy(p8est);

    // Load the new information
    // sizeof(octantData), not quadrantData: the forest was created with the 3D
    // payload (p8est_new_ext above), and p4est aborts inside p8est_load_ext if
    // the size it is told does not match the one in the file.
    p8est=p8est_load_ext(fname, m_mpiInfo->comm, sizeof(octantData), load_data, 
                    autopartition, broadcasthead, &forestData, &connectivity);
    ESYS_ASSERT(p8est_is_valid(p8est),"Invalid p8est file");

    // Update the nodes
    p8est_lnodes_destroy(nodes);
    nodes = p8est_lnodes_new(p8est, ghost, 1);
    
    // Update Brick
    if(autoMeshUpdates)
        updateMesh();

    // Need to update these now that the mesh has changed
    z_needs_update=true;
    iz_needs_update=true;
}

// See the note in Rectangle.cpp: this refinement needs only p4est and does
// not belong inside the trilinos guard.

void Brick::refineMesh(std::string algorithmname)
{
    z_needs_update=true;
    iz_needs_update=true;

    if(!algorithmname.compare("uniform"))
    {
        p8est_refine_ext(p8est, true, -1, refine_uniform, init_brick_data, refine_copy_parent_octant);
        p8est_balance_ext(p8est, P8EST_CONNECT_FULL, init_brick_data, refine_copy_parent_octant);
    }
    else if(!algorithmname.compare("MARE2DEM") || !algorithmname.compare("mare2dem"))
    {
        if(adaptive_refinement == true)
        {
            p8est_refine_ext(p8est, true, -1, refine_mare2dem, init_brick_data, refine_copy_parent_octant);
            p8est_balance_ext(p8est, P8EST_CONNECT_FULL, init_brick_data, refine_copy_parent_octant);
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
    if(!p8est_is_valid(p8est))
        throw OxleyException("p8est broke during refinement");
    #ifdef OXLEY_ENABLE_TIMECONSUMING_DEBUG_CHECKS
    if(!p8est_connectivity_is_valid(connectivity))
        throw OxleyException("connectivity broke during refinement");
    #endif
#endif

    // Update
    if(autoMeshUpdates)
        updateMesh();
}

void Brick::refineBoundary(std::string boundaryname, double dx)
{
    z_needs_update=true;
    iz_needs_update=true;

    forestData->refinement_depth = dx;

    if(!boundaryname.compare("north") || !boundaryname.compare("North")
        || !boundaryname.compare("n") || !boundaryname.compare("N")
        || !boundaryname.compare("NORTH"))
    {
        p8est_refine_ext(p8est, true, -1, refine_north, init_brick_data, refine_copy_parent_octant);
        p8est_balance_ext(p8est, P8EST_CONNECT_FULL, init_brick_data, refine_copy_parent_octant);
    } 
    else if(!boundaryname.compare("south") || !boundaryname.compare("South")
        || !boundaryname.compare("s") || !boundaryname.compare("S")
        || !boundaryname.compare("SOUTH"))
    {
        p8est_refine_ext(p8est, true, -1, refine_south, init_brick_data, refine_copy_parent_octant);
        p8est_balance_ext(p8est, P8EST_CONNECT_FULL, init_brick_data, refine_copy_parent_octant);
    }
    else if(!boundaryname.compare("west") || !boundaryname.compare("West")
        || !boundaryname.compare("w") || !boundaryname.compare("W")
        || !boundaryname.compare("WEST"))
    {
        p8est_refine_ext(p8est, true, -1, refine_west, init_brick_data, refine_copy_parent_octant);
        p8est_balance_ext(p8est, P8EST_CONNECT_FULL, init_brick_data, refine_copy_parent_octant);
    }
    else if(!boundaryname.compare("east") || !boundaryname.compare("East")
        || !boundaryname.compare("e") || !boundaryname.compare("E")
        || !boundaryname.compare("EAST"))
    {
        p8est_refine_ext(p8est, true, -1, refine_east, init_brick_data, refine_copy_parent_octant);
        p8est_balance_ext(p8est, P8EST_CONNECT_FULL, init_brick_data, refine_copy_parent_octant);  
    }
    else if(!boundaryname.compare("top") || !boundaryname.compare("Top")
        || !boundaryname.compare("t") || !boundaryname.compare("T")
        || !boundaryname.compare("EAST"))
    {
        p8est_refine_ext(p8est, true, -1, refine_top, init_brick_data, refine_copy_parent_octant);
        p8est_balance_ext(p8est, P8EST_CONNECT_FULL, init_brick_data, refine_copy_parent_octant);  
    }
    else if(!boundaryname.compare("bottom") || !boundaryname.compare("Bottom")
        || !boundaryname.compare("b") || !boundaryname.compare("B")
        || !boundaryname.compare("EAST"))
    {
        p8est_refine_ext(p8est, true, -1, refine_bottom, init_brick_data, refine_copy_parent_octant);
        p8est_balance_ext(p8est, P8EST_CONNECT_FULL, init_brick_data, refine_copy_parent_octant);  
    }
    else {
        throw OxleyException("Unknown boundary name. Please try 'north', 'east', 'south', 'west', 'top' or 'bottom'.");
    }

    // Make sure that nothing went wrong
#ifdef OXLEY_ENABLE_DEBUG
    if(!p8est_is_valid(p8est))
        throw OxleyException("p8est broke during refinement");
    #ifdef OXLEY_ENABLE_TIMECONSUMING_DEBUG_CHECKS
    if(!p8est_connectivity_is_valid(connectivity))
        throw OxleyException("connectivity broke during refinement");
    #endif
#endif
    
    // Update
    if(autoMeshUpdates)
        updateMesh();
}

void Brick::refineRegion(double x0, double x1, double y0, double y1, double z0, double z1)
{
    z_needs_update=true;
    iz_needs_update=true;

    // If the boundaries were not specified by the user, default to the border of the domain
    forestData->refinement_boundaries[0] = x0 == -1 ? forestData->m_origin[0] : x0; 
    forestData->refinement_boundaries[1] = x1 == -1 ? forestData->m_origin[1] : x1;
    forestData->refinement_boundaries[2] = y0 == -1 ? forestData->m_lxyz[0] : y0;
    forestData->refinement_boundaries[3] = y1 == -1 ? forestData->m_lxyz[1] : y1;
    forestData->refinement_boundaries[4] = z0 == -1 ? forestData->m_lxyz[0] : z0;
    forestData->refinement_boundaries[5] = z1 == -1 ? forestData->m_lxyz[1] : z1;

    p8est_refine_ext(p8est, true, -1, refine_region, init_brick_data, refine_copy_parent_octant);
    
    // Make sure that nothing went wrong
#ifdef OXLEY_ENABLE_DEBUG
    if(!p8est_is_valid(p8est))
        throw OxleyException("p8est broke during refinement");
    #ifdef OXLEY_ENABLE_TIMECONSUMING_DEBUG_CHECKS
    if(!p8est_connectivity_is_valid(connectivity))
        throw OxleyException("connectivity broke during refinement");
    #endif
#endif
    
    // Update
    if(autoMeshUpdates)
        updateMesh();
}

void Brick::refinePoint(double x0, double y0, double z0)
{
    z_needs_update=true;
    iz_needs_update=true;

    // Check that the point is inside the domain
    if(    x0 < forestData->m_origin[0] || x0 > forestData->m_lxyz[0] 
        || y0 < forestData->m_origin[1] || y0 > forestData->m_lxyz[1] 
        || z0 < forestData->m_origin[2] || z0 > forestData->m_lxyz[2] )
    {
        std::string message = "INFORMATION: Coordinates lie outside the domain. (" + std::to_string(x0) + ", " + std::to_string(y0) + ", " + std::to_string(z0) + ") skipping point";
        std::cout << message << std::endl;
        throw OxleyException(message);
    }

    // Copy over information required by the backend
    forestData->refinement_boundaries[0] = x0;
    forestData->refinement_boundaries[1] = y0;
    forestData->refinement_boundaries[2] = z0;

    // Make sure that nothing went wrong
#ifdef OXLEY_ENABLE_DEBUG_REFINEPOINT
    if(p8est_is_valid(p8est)!=1)
        throw OxleyException("p8est broke during refinement");
    else
        std::cout << "\033[1;34m[oxley]\033[0m refinePoint: valid p8est" << std::endl;
    #ifdef OXLEY_ENABLE_TIMECONSUMING_DEBUG_CHECKS
    if(p8est_connectivity_is_valid(connectivity)!=1)
        throw OxleyException("connectivity broke during refinement");
    else
        std::cout << "\033[1;34m[oxley]\033[0m refinePoint: valid connectivity" << std::endl;
    #endif
#endif

    p8est_refine_ext(p8est, true, -1, refine_point, init_brick_data, refine_copy_parent_octant);

    // Update
    if(autoMeshUpdates)
        updateMesh();
}

void Brick::refineSphere(double x0, double y0, double z0, double r)
{
    z_needs_update=true;
    iz_needs_update=true;
    
    // Check that the point is inside the domain
    if(x0 < forestData->m_origin[0] || x0 > forestData->m_lxyz[0] 
        || y0 < forestData->m_origin[1] || y0 > forestData->m_lxyz[1] 
        || z0 < forestData->m_origin[2] || z0 > forestData->m_lxyz[2] )
    {
        throw OxleyException("Coordinates lie outside the domain.");
    }

    // If the boundaries were not specified by the user, default to the border of the domain
    forestData->refinement_boundaries[0] = x0;
    forestData->refinement_boundaries[1] = y0;
    forestData->refinement_boundaries[2] = z0;
    forestData->refinement_boundaries[3] = r;
    p8est_refine_ext(p8est, true, -1, refine_sphere, init_brick_data, refine_copy_parent_octant);
    p8est_balance_ext(p8est, P8EST_CONNECT_FULL, init_brick_data, refine_copy_parent_octant);

    // Make sure that nothing went wrong
#ifdef OXLEY_ENABLE_DEBUG
    if(!p8est_is_valid(p8est))
        throw OxleyException("p8est broke during refinement");
    #ifdef OXLEY_ENABLE_TIMECONSUMING_DEBUG_CHECKS
    if(!p8est_connectivity_is_valid(connectivity))
        throw OxleyException("connectivity broke during refinement");
    #endif
#endif

    // Update
    if(autoMeshUpdates)
        updateMesh();
}
void Brick::refineMask(escript::Data mask)
{
    z_needs_update=true;
    iz_needs_update=true;

    // If the boundaries were not specified by the user, default to the border of the domain
    forestData->mask = mask;
    p8est_refine_ext(p8est, true, -1, refine_mask, init_brick_data, refine_copy_parent_octant);
    p8est_balance_ext(p8est, P8EST_CONNECT_FULL, init_brick_data, refine_copy_parent_octant);

    // Make sure that nothing went wrong
#ifdef OXLEY_ENABLE_DEBUG
    if(!p8est_is_valid(p8est))
        throw OxleyException("p8est broke during refinement");
    #ifdef OXLEY_ENABLE_TIMECONSUMING_DEBUG_CHECKS
    if(!p8est_connectivity_is_valid(connectivity))
        throw OxleyException("connectivity broke during refinement");
    #endif
#endif
    
    // Update
    if(autoMeshUpdates)
        updateMesh();
}

escript::Data Brick::getX() const
{
    escript::Data out=escript::Vector(0,escript::continuousFunction(*this),true);
    setToX(out);
    out.setProtection();
    return out;
}

void Brick::print_debug_report(std::string locat)
{
    std::cout << "report for " <<  locat << std::endl;
    std::cout << "p8est = " << &p8est << std::endl;
    if(!p8est_is_valid(p8est))
        std::cout << "WARNING: p8est is invalid" << std::endl;
    std::cout << "forestData = " << &forestData << std::endl;
    std::cout << "connectivity = " << &connectivity << std::endl;
    #ifdef OXLEY_ENABLE_TIMECONSUMING_DEBUG_CHECKS
    if(!p8est_connectivity_is_valid(connectivity))
        std::cout << "WARNING: connectivity is invalid" << std::endl;
    #endif
    std::cout << "temp_data = " << &temp_data << std::endl;
}

Assembler_ptr Brick::createAssembler(std::string type, const DataMap& constants) const
{
#ifndef ESYS_HAVE_TRILINOS
    // The 3D assembler is trilinos-only, and PDEs on an adaptive oxley mesh
    // are solved by exporting to finley anyway - so a build without trilinos
    // simply has no assembler, rather than half of one.
    throw escript::NotImplementedError("oxley was built without trilinos and "
            "has no assembler. Export the mesh with toFinley() and solve on "
            "the finley domain, which is how an adaptive oxley mesh is meant "
            "to be solved on in any case.");
#else
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
            return Assembler_ptr(new DefaultAssembler3D<cplx_t>(shared_from_this()));
        } else {
            return Assembler_ptr(new DefaultAssembler3D<real_t>(shared_from_this()));
        }
    } 
    throw escript::NotImplementedError("oxley::brick does not support the requested assembler");
#endif
}

// return True for a boundary node and False for an internal node
bool Brick::isBoundaryNode(p8est_quadrant_t * quad, int n, p8est_topidx_t treeid, p8est_qcoord_t length) const
{
    //TODO
    double lx = length * ((int) (n % 2) == 1);
    double ly = length * ((int) (n / 2) == 1);
    double lz = length * ((int) (n / 2) == 1); //TODO
    double xy[3];
    p8est_qcoord_to_vertex(p8est->connectivity, treeid, quad->x+lx, quad->y+ly, quad->z+lz, xy);
    return (xy[0] == forestData->m_origin[0]) || (xy[0] == forestData->m_lxyz[0]) 
        || (xy[1] == forestData->m_origin[1]) || (xy[1] == forestData->m_lxyz[1])
        || (xy[2] == forestData->m_origin[2]) || (xy[3] == forestData->m_lxyz[2]);
}

// returns True for a boundary node on the north or east of the domain
bool Brick::isUpperBoundaryNode(p8est_quadrant_t * quad, int n, p8est_topidx_t treeid, p8est_qcoord_t length) const
{
    double lx = length * ((int) (n % 2) == 1);
    double ly = length * ((int) (n / 2) == 1);
    double lz = length * ((int) (n / 2) == 1); //TODO
    double xy[3];
    p8est_qcoord_to_vertex(p8est->connectivity, treeid, quad->x+lx, quad->y+ly, quad->z+lz, xy);
    return (xy[0] == forestData->m_lxyz[0]) || (xy[1] == forestData->m_lxyz[1]) || (xy[2] == forestData->m_lxyz[2]);
}

// returns True for a boundary node on the south or west of the domain
bool Brick::isLowerBoundaryNode(p8est_quadrant_t * quad, int n, p8est_topidx_t treeid, p8est_qcoord_t length) const
{
    double lx = length * ((int) (n % 2) == 1);
    double ly = length * ((int) (n / 2) == 1);
    double lz = length * ((int) (n / 2) == 1); //TODO
    double xy[3];
    p8est_qcoord_to_vertex(p8est->connectivity, treeid, quad->x+lx, quad->y+ly, quad->z+lz, xy);
    return (xy[0] == forestData->m_origin[0]) 
        || (xy[1] == forestData->m_origin[1])
        || (xy[2] == forestData->m_origin[2]);
}

// returns True for a boundary node on the left boundary
bool Brick::isLeftBoundaryNode(p8est_quadrant_t * quad, int n, p8est_topidx_t treeid, p8est_qcoord_t length) const
{
    double lx = length * ((int) (n % 2) == 1);
    double ly = length * ((int) (n / 2) == 1);
    double lz = length * ((int) (n / 2) == 1); //TODO
    double xy[3];
    p8est_qcoord_to_vertex(p8est->connectivity, treeid, quad->x+lx, quad->y+ly, quad->z+lz, xy);
    return (xy[0] == forestData->m_origin[0]);
}

// returns True for a boundary node on the right boundary
bool Brick::isRightBoundaryNode(p8est_quadrant_t * quad, int n, p8est_topidx_t treeid, p8est_qcoord_t length) const
{
    double lx = length * ((int) (n % 2) == 1);
    double ly = length * ((int) (n / 2) == 1);
    double lz = length * ((int) (n / 2) == 1); //TODO
    double xy[3];
    p8est_qcoord_to_vertex(p8est->connectivity, treeid, quad->x+lx, quad->y+ly, quad->z+lz, xy);
    return (xy[0] == forestData->m_lxyz[0]);
}

// returns True for a boundary node on the bottom boundary
bool Brick::isBottomBoundaryNode(p8est_quadrant_t * quad, int n, p8est_topidx_t treeid, p8est_qcoord_t length) const
{
    double lx = length * ((int) (n % 2) == 1);
    double ly = length * ((int) (n / 2) == 1);
    double lz = length * ((int) (n / 2) == 1); //TODO
    double xy[3];
    p8est_qcoord_to_vertex(p8est->connectivity, treeid, quad->x+lx, quad->y+ly, quad->z+lz, xy);
    return (xy[1] == forestData->m_origin[1]);
}

// returns True for a boundary node on the top boundary
bool Brick::isTopBoundaryNode(p8est_quadrant_t * quad, int n, p8est_topidx_t treeid, p8est_qcoord_t length) const
{
    double lx = length * ((int) (n % 2) == 1);
    double ly = length * ((int) (n / 2) == 1);
    double lz = length * ((int) (n / 2) == 1); //TODO
    double xy[3];
    p8est_qcoord_to_vertex(p8est->connectivity, treeid, quad->x+lx, quad->y+ly, quad->z+lz, xy);
    return (xy[1] == forestData->m_lxyz[1]);
}

// returns True for a boundary node on the top boundary
bool Brick::isAboveBoundaryNode(p8est_quadrant_t * quad, int n, p8est_topidx_t treeid, p8est_qcoord_t length) const
{
    double lx = length * ((int) (n % 2) == 1);
    double ly = length * ((int) (n / 2) == 1);
    double lz = length * ((int) (n / 2) == 1); //TODO
    double xy[3];
    p8est_qcoord_to_vertex(p8est->connectivity, treeid, quad->x+lx, quad->y+ly, quad->z+lz, xy);
    return (xy[2] == forestData->m_lxyz[3]); //todo check
}

// returns True for a boundary node on the top boundary
bool Brick::isBelowBoundaryNode(p8est_quadrant_t * quad, int n, p8est_topidx_t treeid, p8est_qcoord_t length) const
{
    double lx = length * ((int) (n % 2) == 1);
    double ly = length * ((int) (n / 2) == 1);
    double lz = length * ((int) (n / 2) == 1); //TODO
    double xy[3];
    p8est_qcoord_to_vertex(p8est->connectivity, treeid, quad->x+lx, quad->y+ly, quad->z+lz, xy);
    return (xy[2] == forestData->m_lxyz[4]); //todo check
}


// return True for a hanging node and False for an non-hanging node
bool Brick::isHangingNode(p8est_lnodes_code_t face_code, int n) const
{
    if(face_code == 0)
    {
        return false;
    }
    else
    {
        // TODO
        int8_t c = face_code & 1;
        int8_t ishanging0 = (face_code >> 2) & 1;
        int8_t ishanging1 = (face_code >> 3) & 1;

        int8_t f0 = p8est_corner_faces[c][0];
        int8_t f1 = p8est_corner_faces[c][1];

        // int d0 = c / 2 == 0;
        // int d1 = c % 2 == 0;
        
        return ((ishanging0 == 1) && (p8est_corner_face_corners[c][f0] == 1)) 
            || ((ishanging1 == 1) && (p8est_corner_face_corners[c][f1] == 1));
    }
}

bool Brick::getHangingInfo(p8est_lnodes_code_t face_code, 
            int hanging_face[P8EST_FACES], int hanging_edge[P8EST_EDGES]) const
{
    //cf. p8est_lnodes.h:162-175

    ESYS_ASSERT(face_code >= 0, "getHangingInfo: Invalid face code.");

    memset(hanging_face, -1, 6 * sizeof (int));
    memset(hanging_edge, -1, 12 * sizeof (int));

    if(face_code) {
        int16_t c = face_code & 0x0007;
        int16_t work = face_code >> 3;

        int16_t cwork = c;
        for(int i = 0; i < 3; ++i) {
            if(work & 0x0001) 
            {
                int f = p8est_corner_faces[c][i];
                hanging_face[f] = p8est_corner_face_corners[c][f];
                for(int j = 0; j < 4; j++) 
                {
                    int e = p8est_face_edges[f][j];
                    hanging_edge[e] = 4;
                }
            }
            work >>= 1;
        }

        for(int i = 0; i < 3; ++i) 
        {
            if (work & 0x0001) 
            {
                int e = p8est_corner_edges[c][i];
                hanging_edge[e] = (hanging_edge[e] == -1) ? 0 : 2;
                hanging_edge[e] += (int) (cwork & 0x0001);
            }
            cwork >>= 1;
            work >>= 1;
        }
        return 1;
    }
    else 
    {
        return 0;
    }
}

/**
   \brief
   Returns true if the node is hanging
*/
bool Brick::hasHanging(p8est_lnodes_code_t face_code) const
{
    ESYS_ASSERT(face_code >= 0, "hasHanging: Invalid face code.");
    return face_code == 0 ? false : true;
}

/**
   \brief
   Returns true if the node is on the boundary
*/
bool Brick::onUpperBoundary(double x, double y, double z) const
{
    // return (forestData->m_origin[0] == x) || (forestData->m_origin[1] == y) || (forestData->m_origin[2] == z) 
    //     || (forestData->m_lxyz[0] == x) || (forestData->m_lxyz[1] == y) || (forestData->m_lxyz[2] == z);
    return (forestData->m_lxyz[0] == x) || (forestData->m_lxyz[1] == y) || (forestData->m_lxyz[2] == z);

}

//protected
void Brick::reset_ghost()
{
    p8est_ghost_destroy(ghost);
    ghost = p8est_ghost_new(p8est, P8EST_CONNECT_FULL);
#ifdef OXLEY_ENABLE_DEBUG
    if(p8est_ghost_is_valid(p8est,ghost)!=1)
        throw OxleyException("\033[1;34m[oxley]\033[0m ghost is broken");
    else
        std::cout << "\033[1;34m[oxley]\033[0m valid ghost" << std::endl;
#endif
}


void Brick::renumberNodes()
{
    oxleytimer.toc("renumberNodes...");

    // The global node numbering now comes from p4est_lnodes (built in the
    // constructor / after refinement); this routine only derives m_nodeId,
    // the global id of each local node, from the lnodes owned/ghost partition.
    // The legacy coordinate-hash containers are retired.
    octantInfo.clear();

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

    // MPI: build the ghost octant halo and extend myColumns with any 2nd-layer
    // ghost nodes (nodes that appear only on ghost octants). (A6.)
    buildParallelOverlap();

    oxleytimer.toc("renumberNodes...Done");
}

//protected
void Brick::buildParallelOverlap()
{
    m_ghostElemNodes.clear();
    if (m_mpiInfo->size <= 1 || !ghost)
        return;

    const int V = nodes->vnodes;                       // 8 corners (degree-1)
    const long nLocal = (long) nodes->num_local_nodes;
    const long nLocalElem = (long) nodes->num_local_elements;

    std::unordered_map<long,long> g2l;
    g2l.reserve((size_t) nLocal * 2);
    for (long i = 0; i < nLocal; ++i)
        g2l[(long) m_nodeId[i]] = i;

    const long nGhost  = (long) ghost->ghosts.elem_count;
    const long nMirror = (long) ghost->mirrors.elem_count;

    // Each local octant's V corner GLOBAL node ids (mirror send source).
    std::vector<p4est_gloidx_t> localElemGN((size_t) nLocalElem * V);
    for (long e = 0; e < nLocalElem; ++e)
        for (int c = 0; c < V; ++c)
            localElemGN[(size_t) e*V + c] =
                (p4est_gloidx_t) m_nodeId[ nodes->element_nodes[(size_t) e*V + c] ];

    std::vector<void*> mirror_data((size_t) nMirror, nullptr);
    for (long m = 0; m < nMirror; ++m) {
        p8est_quadrant_t* mq = p8est_quadrant_array_index(&ghost->mirrors, m);
        const long le = (long) mq->p.piggy3.local_num;
        mirror_data[m] = (void*) &localElemGN[(size_t) le*V];
    }

    std::vector<p4est_gloidx_t> ghostElemGN((size_t) nGhost * V);
    p8est_ghost_exchange_custom(p8est, ghost,
                                (size_t) V * sizeof(p4est_gloidx_t),
                                mirror_data.data(), ghostElemGN.data());

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














// void Brick::renumberNodes()
// {
//     oxleytimer.toc("Starting renumberNodes...");

//     #ifdef OXLEY_ENABLE_DEBUG_RENUMBER_NODES
//     std::cout << "\033[1;34m[oxley]\033[0m Renumbering nodes...." << std::endl;
//     #endif

//     // Clear some variables
//     octantIDs.clear();
//     NodeIDs.clear();
//     // hanging_face_orientation.clear();
//     // hanging_edge_orientation.clear();
//     // hanging_edge_node_connections.clear();
//     // hanging_face_node_connections.clear();
//     // false_node_connections.clear();

//     // //                                                                              
//     // //                       0         1         2         3         4         5
//     // std::vector<std::vector<int>> xface = {{0,0,0,0},{1,1,1,1},{1,0,1,0},{1,0,1,0},{1,0,1,0},{1,0,1,0}};
//     // std::vector<std::vector<int>> yface = {{1,0,1,0},{1,0,1,0},{0,0,0,0},{1,1,1,1},{1,1,0,0},{1,1,0,0}};
//     // std::vector<std::vector<int>> zface = {{1,1,0,0},{1,1,0,0},{1,1,0,0},{1,1,0,0},{0,0,0,0},{1,1,1,1}};
    
//     // // used to calculate up the coordinate of the final node in each octant
//     // //                       0         1         2         3         4         5
//     // std::vector<std::vector<int>> xface0 = {{0,0,0,0},{1,1,1,1},{0,1,0,1},{0,1,0,1},{0,1,0,1},{0,1,0,1}};
//     // std::vector<std::vector<int>> yface0 = {{0,1,0,1},{0,1,0,1},{0,0,0,0},{1,1,1,1},{0,0,1,1},{0,0,1,1}};
//     // std::vector<std::vector<int>> zface0 = {{0,0,1,1},{0,0,1,1},{0,0,1,1},{0,0,1,1},{0,0,0,0},{1,1,1,1}};
    
//     // // used to look up the coordinate of the final node in each octant
//     // //                              0            1            2           3            4          5
//     // std::vector<std::vector<signed int>> xface1 = {{0, 0, 0, 0},{1, 1, 1, 1},{2,-1,0, 1},{2,-1, 0, 1},{0,0,0, 1},{2,-1, 0,0}};
//     // std::vector<std::vector<signed int>> yface1 = {{2,-1, 0, 1},{2,-1, 0, 1},{0, 0,0, 0},{1, 1, 1, 1},{0,0,0,-1},{0, 0,-1,0}};
//     // std::vector<std::vector<signed int>> zface1 = {{0, 0,-1,-1},{0, 0,-1,-1},{0, 0,0,-1},{0, 0,-1, 1},{0,0,0, 0},{1, 1, 1,0}};
//     // std::vector<std::vector<signed int>> xface2 = {{0, 0, 0, 0},{1, 1, 1, 1},{0, 1,2,-1},{0, 0, 2,-1},{0,1,2,-1},{0, 0, 0,0}};
//     // std::vector<std::vector<signed int>> yface2 = {{0, 1, 2,-1},{0, 1, 2,-1},{0, 0,0, 0},{1, 0, 1, 1},{0,2,1, 1},{2, 0, 0,0}};
//     // std::vector<std::vector<signed int>> zface2 = {{2, 2, 1, 1},{2, 2, 1, 1},{0, 2,1, 1},{2, 0, 1, 1},{0,0,0, 0},{1, 0, 0,0}};

//     // // three dimensions, six octant faces, four orientations (relative to parent octant), two edge nodes
//     // //                           0     1     2     3
//     // std::vector<std::vector<std::vector<int>>> xedge =
//     //                                              {{{0,0},{0,0},{0,0},{0,0}}, //0
//     //                                             {{1,1},{1,1},{1,1},{1,1}}, //1
//     //                                             {{1,0},{0,1},{0,1},{1,0}}, //2
//     //                                             {{1,0},{0,1},{0,1},{1,0}}, //3
//     //                                             {{1,0},{0,1},{0,1},{1,0}}, //4
//     //                                             {{1,0},{0,1},{0,1},{1,0}}};//5
//     // //                           0     1     2     3
//     // std::vector<std::vector<std::vector<int>>> yedge =
//     //                                              {{{1,0},{0,1},{0,1},{1,0}}, //0
//     //                                             {{1,0},{0,1},{0,1},{1,0}}, //1
//     //                                             {{0,0},{0,0},{0,0},{0,0}}, //2
//     //                                             {{1,1},{1,1},{1,1},{1,1}}, //3
//     //                                             {{0,1},{0,1},{0,1},{0,1}}, //4
//     //                                             {{0,1},{0,1},{0,1},{0,1}}};//5
//     // //                           0     1     2     3
//     // std::vector<std::vector<std::vector<int>>> zedge =
//     //                                              {{{0,1},{0,1},{0,1},{0,1}}, //0
//     //                                             {{0,1},{0,1},{0,1},{0,1}}, //1
//     //                                             {{0,1},{0,1},{0,1},{0,1}}, //2
//     //                                             {{0,1},{0,1},{0,1},{0,1}}, //3
//     //                                             {{0,0},{0,0},{0,0},{0,0}}, //4
//     //                                             {{1,1},{1,1},{1,1},{1,1}}};//5

//     // // six faces, four orientations (relative to parent octant), two edges
//     // //                                0       1       2       3
//     // std::vector<std::vector<std::vector<int>>> edge_lookup = {{{4,8},  {4,10}, {8,6},  {10,6}}, //0
//     //                                                           {{5,9},  {5,11}, {9,7},  {11,7}}, //1
//     //                                                           {{0,8},  {0,9},  {8,2},  {9,2}},  //2
//     //                                                           {{1,10}, {1,11}, {10,3}, {11,3}}, //3
//     //                                                           {{0,4},  {0,5},  {4,1},  {5,1}},  //4
//     //                                                           {{2,6},  {2,7},  {6,3},  {7,3}}}; //5

//     std::vector<std::vector<float>> lxy_nodes = {{0,0,0},{1,0,0},{0,1,0},{1,1,0},{0,0,1},{1,0,1},{0,1,1},{1,1,1}};

//     // Write in NodeIDs
//     int k = 0;
//     std::vector<DoubleTuple> NormalNodesTmp;
//     NormalNodes.clear();
//     std::vector<long> HangingFaceNodesTmp;
//     HangingFaceNodes.clear();
//     std::vector<long> HangingEdgeNodesTmp;
//     HangingEdgeNodes.clear();


// #ifdef ESYS_MPI 
//     // Do the calculation   
//     oxleytimer.toc("\tDoing the calculation");
//     int guessSize = 0;
//     for(p8est_topidx_t treeid = p8est->first_local_tree; treeid <= p8est->last_local_tree; ++treeid) 
//     {
//         p8est_tree_t * tree = p8est_tree_array_index(p8est->trees, treeid);
//         sc_array_t * tquadrants = &tree->quadrants;
//         p8est_locidx_t Q = (p8est_locidx_t) tquadrants->elem_count;
//         guessSize+=Q;
//     }
//     // guessSize*=8;

//     #ifdef OXLEY_ENABLE_DEBUG_RENUMBERNODES_INFORMATIONAL
//         std::cout << "guessSize = " << guessSize << std::endl;
//     #endif

//     std::vector<double> treevecX;
//     std::vector<double> treevecY;
//     std::vector<double> treevecZ;
//     std::vector<DoubleTuple> localNodes;

//     if(m_mpiInfo->rank == 0)
//     {
//         NormalNodesTmp.resize(guessSize,std::make_tuple(-1.0,-1.0,-1.0));
//     }
//     else
//     {
//         treevecX.resize(guessSize, -1.0);
//         treevecY.resize(guessSize, -1.0);
//         treevecZ.resize(guessSize, -1.0);
//         localNodes.resize(guessSize, std::make_tuple(-1.0,-1.0,-1.0));
//     }
//     long actualSize = 0;

//     #ifdef OXLEY_ENABLE_DEBUG_RENUMBERNODES_INFORMATIONAL
//         std::cout << "starting loop " << std::endl;
//     #endif

//     for(p8est_topidx_t treeid = p8est->first_local_tree; treeid <= p8est->last_local_tree; ++treeid) 
//     {
//         p8est_tree_t * tree = p8est_tree_array_index(p8est->trees, treeid);
//         sc_array_t * tquadrants = &tree->quadrants;
//         p8est_locidx_t Q = (p8est_locidx_t) tquadrants->elem_count;
//         for(int q = 0; q < Q; ++q) { 
//             p8est_quadrant_t * oct = p8est_quadrant_array_index(tquadrants, q);
//             p8est_qcoord_t l = P8EST_QUADRANT_LEN(oct->level);
//             for(int n = 0; n < 8; n++) {
//                 double xyzB[3];
//                 p8est_qcoord_to_vertex(p8est->connectivity, treeid, oct->x+l*lxy_nodes[n][0], oct->y+l*lxy_nodes[n][1], oct->z+l*lxy_nodes[n][2], xyzB);
//                 if(m_mpiInfo->rank == 0)
//                 {
//                     std::tuple<double,double,double> point = std::make_tuple(xyzB[0],xyzB[1],xyzB[2]);
//                     bool got_already = gotAlready(point,NormalNodesTmp);
//                     if(!got_already)
//                     {
//                         if(actualSize < guessSize)
//                         {
//                             NormalNodesTmp[actualSize] = point;
//                         }
//                         else
//                         {
//                             NormalNodesTmp.push_back(point);
//                         }
//                     }
//                     actualSize++;
//                 }
//                 else
//                 {
//                     std::tuple<double,double,double> point = std::make_tuple(xyzB[0],xyzB[1],xyzB[2]);
//                     bool got_already = gotAlready(point,NormalNodesTmp);
//                     if(!got_already)
//                     {
//                         if(actualSize < guessSize)
//                         {
//                             localNodes[actualSize] = point;
//                             treevecX[actualSize] = xyzB[0];
//                             treevecY[actualSize] = xyzB[1];
//                             treevecZ[actualSize] = xyzB[2];
//                         }
//                         else
//                         {
//                             localNodes.push_back(point);
//                             treevecX.push_back(xyzB[0]);
//                             treevecY.push_back(xyzB[1]);
//                             treevecZ.push_back(xyzB[2]);
//                         }
//                     }
//                     actualSize++;
//                 }
//             }
//         }
//     }

//     #ifdef OXLEY_ENABLE_DEBUG_RENUMBERNODES_INFORMATIONAL
//         std::cout << "actualSize = " << actualSize << std::endl;
//         std::cout << "guessSize = " << guessSize << std::endl;
//     #endif

//     if(m_mpiInfo->size > 1)
//     {
//         oxleytimer.toc("\tcommunicating A");
//         // Send info to rank 0
//         if(m_mpiInfo->rank != 0)
//         {
//             // broadcast
//             long numPoints = actualSize;

//             std::cout << "numPoints = " << numPoints << std::endl;

//             MPI_Send(&numPoints, 1, MPI_LONG, 0, 0, m_mpiInfo->comm);
//             MPI_Send(treevecX.data(), numPoints, MPI_DOUBLE, 0, 0, m_mpiInfo->comm);
//             MPI_Send(treevecY.data(), numPoints, MPI_DOUBLE, 0, 0, m_mpiInfo->comm);
//             MPI_Send(treevecZ.data(), numPoints, MPI_DOUBLE, 0, 0, m_mpiInfo->comm);
//         }
//         else if(m_mpiInfo->rank == 0)
//         {
//             for(int i = 1; i < m_mpiInfo->size; i++)
//             {   
//                 std::vector<double> tmpX;
//                 std::vector<double> tmpY;
//                 std::vector<double> tmpZ;
                
//                 long numPoints = 0;

//                 std::cout << "numPoints = " << numPoints << std::endl;

//                 MPI_Recv(&numPoints, 1, MPI_LONG, i, 0, m_mpiInfo->comm, MPI_STATUS_IGNORE);

//                 tmpX.resize(numPoints);
//                 tmpY.resize(numPoints);
//                 tmpZ.resize(numPoints);

//                 MPI_Recv(tmpX.data(), numPoints, MPI_DOUBLE, i, 0, m_mpiInfo->comm, MPI_STATUS_IGNORE);
//                 MPI_Recv(tmpY.data(), numPoints, MPI_DOUBLE, i, 0, m_mpiInfo->comm, MPI_STATUS_IGNORE);
//                 MPI_Recv(tmpZ.data(), numPoints, MPI_DOUBLE, i, 0, m_mpiInfo->comm, MPI_STATUS_IGNORE);
//                 for(int j = 0; j < numPoints; j++)
//                     NormalNodesTmp.push_back(std::make_tuple(tmpX[j],tmpY[j],tmpZ[j]));
//             }
//         }
//         oxleytimer.toc("\t\t\t....done");
//     }
    

//     if(m_mpiInfo->size > 1)
//     {
//         // Send info to the other ranks
//         oxleytimer.toc("\tcommunicating");
//         if(m_mpiInfo->size > 1)
//         {
//             if(m_mpiInfo->rank == 0)
//             {
//                 int numPoints = NormalNodesTmp.size();
//                 for(int i = 1; i < m_mpiInfo->size; i++)
//                 {
//                     MPI_Send(&numPoints, 1, MPI_INT, i, 0, m_mpiInfo->comm);
//                     MPI_Send(NormalNodesTmp.data(), 3*NormalNodesTmp.size(), MPI_DOUBLE, i, 0, m_mpiInfo->comm);
//                 }
//             }
//             else
//             {
//                 int numPoints = 0;
//                 MPI_Recv(&numPoints, 1, MPI_INT, 0, 0, m_mpiInfo->comm, MPI_STATUS_IGNORE);        
//                 std::vector<DoubleTuple> temp(numPoints, std::make_tuple(-1.,-1.,-1.));
//                 MPI_Recv(temp.data(), 3*numPoints, MPI_DOUBLE, 0, 0, m_mpiInfo->comm, MPI_STATUS_IGNORE );
//                 NormalNodesTmp = temp;
//             }
//         }
//         oxleytimer.toc("\t\t\t...done");
 
//         oxleytimer.toc("\tsearching for duplicates");
//         int proportion = NormalNodesTmp.size() / m_mpiInfo->size;
//         int first_point = m_mpiInfo->rank * proportion;
//         int last_point = ((m_mpiInfo->rank + 1) * proportion) - 1;
//         if(m_mpiInfo->rank == m_mpiInfo->size - 1)
//             last_point = NormalNodesTmp.size();

//         std::vector<int> locats(NormalNodesTmp.size(),false);
//         for(int i = first_point; i < last_point; i++)
//         {
//             // skip points already tagged as duplicates
//             if(locats[i]==true)
//                 continue;

//             // get point
//             auto point = NormalNodesTmp[i];
//             bool duplicatePoint = hasDuplicate(point,NormalNodesTmp, true);
//             if(duplicatePoint)
//             {
//                 bool firstOccurance = true;
//                 for(int j = 0; j < NormalNodesTmp.size(); j++) // find location of duplicate points
//                 {
//                     if(     (std::abs(std::get<0>(NormalNodesTmp[j]) - std::get<0>(point)) < tuple_tolerance)
//                          && (std::abs(std::get<1>(NormalNodesTmp[j]) - std::get<1>(point)) < tuple_tolerance) 
//                          && (std::abs(std::get<2>(NormalNodesTmp[j]) - std::get<2>(point)) < tuple_tolerance))
//                     {
//                         if(firstOccurance)
//                         {
//                             firstOccurance = false;
//                             continue;
//                         }
//                         else
//                         {
//                             locats[j] = 1;
//                         }
//                     }
//                 }
//             }
//         }

//         oxleytimer.toc("telling rank 0 which points need to be erased");

//         // tell rank 0 which points need to be erased
//         if(m_mpiInfo->rank == 0)
//         {
//             for(int r = 1; r < m_mpiInfo->size; r++)
//             {
//                 std::vector<int> new_locats(NormalNodesTmp.size(),false);
//                 MPI_Recv(new_locats.data(),NormalNodesTmp.size(),MPI_INT,r,0,m_mpiInfo->comm,MPI_STATUS_IGNORE);

//                 for(int i = 0; i < locats.size(); i++)
//                 {
//                     if(new_locats[i] == 1)
//                         locats[i] = 1;
//                 }
//             }

//             // erase duplicates
//             int num_duplicates = 0;
//             for(int i = 0; i < locats.size(); i++)
//                 if(locats[i] == 1)
//                     num_duplicates++;
//             std::vector<DoubleTuple> NormalNodes_NoDuplicates(NormalNodesTmp.size()-num_duplicates,std::make_tuple(-1.,-1.,-1));
//             int duplicatecounter=0;
//             for(int i = 0; i < NormalNodesTmp.size(); i++)
//             {
//                 if(locats[i] == 0)
//                 {
//                     NormalNodes_NoDuplicates[duplicatecounter]=NormalNodesTmp[i];
//                     duplicatecounter++;
//                 }
//             }
//             NormalNodesTmp=NormalNodes_NoDuplicates;
//         }
//         else
//         {
//             MPI_Send(locats.data(),locats.size(),MPI_INT,0,0,m_mpiInfo->comm);
//         }

//         MPI_Barrier(m_mpiInfo->comm);
//         // send info to the other ranks
//         if(m_mpiInfo->rank == 0)
//         {
//             int num = NormalNodesTmp.size();
//             for(int r = 1; r < m_mpiInfo->size; r++)
//             {
//                 MPI_Send(&num,1,MPI_INT,r,0,m_mpiInfo->comm);
//                 MPI_Send(NormalNodesTmp.data(),3*num,MPI_DOUBLE,r,0,m_mpiInfo->comm);
//             }
//         }
//         else
//         {
//             int num;
//             MPI_Recv(&num,1,MPI_INT,0,0,m_mpiInfo->comm,MPI_STATUS_IGNORE);
//             MPI_Recv(NormalNodesTmp.data(),3*num,MPI_DOUBLE,0,0,m_mpiInfo->comm,MPI_STATUS_IGNORE);
//         }


//         // // copy of local collection of points
//         // oxleytimer.toc("\terasing duplicates");
//         // std::vector<DoubleTuple> points_to_send(last_point - first_point - num_local_duplicate_points,std::make_tuple(-10.,-10.,-10.));
//         // int counter=0;
//         // for(int i = first_point; i < last_point; i++)
//         // {
//         //     if(locats[i] != 1) //if not duplicate
//         //         points_to_send[counter++]=NormalNodesTmp[i];
//         // }

//         // oxleytimer.toc("communicating");
//         // oxleytimer.toc("\ttalking to 0");
//         // if(m_mpiInfo->rank == 0)
//         // {
//         //     NormalNodesTmp=points_to_send;


//         //     for(int r = 1; r < m_mpiInfo->size; r++)
//         //     {
//         //         int num;
//         //         MPI_Recv(&num,1,MPI_INT,r,0,m_mpiInfo->comm,MPI_STATUS_IGNORE);
//         //         std::vector<DoubleTuple> temp;
//         //         temp.assign(num,std::make_tuple(-1.,-1.,-1.));
//         //         MPI_Recv(temp.data(),3*num,MPI_DOUBLE,r,0,m_mpiInfo->comm,MPI_STATUS_IGNORE);
//         //         std::vector<DoubleTuple> temp2;
//         //         temp2.reserve(NormalNodesTmp.size()+temp.size());
//         //         temp2.insert(temp2.end(),NormalNodesTmp.begin(),NormalNodesTmp.end());
//         //         temp2.insert(temp2.end(),temp.begin(),temp.end());
//         //         NormalNodesTmp=temp2;
//         //     }
//         // }
//         // else
//         // {
//         //     int num = points_to_send.size();
//         //     MPI_Send(&num,1,MPI_INT,0,0,m_mpiInfo->comm);
//         //     MPI_Send(points_to_send.data(),3*num,MPI_DOUBLE,0,0,m_mpiInfo->comm);
//         // }
//         // oxleytimer.toc("\t\t\t....done");
//         // oxleytimer.toc("\tlistening to 0");
//         // MPI_Barrier(m_mpiInfo->comm);
//         // if(m_mpiInfo->rank == 0)
//         // {
//         //     int num = NormalNodesTmp.size();
//         //     for(int r = 1; r < m_mpiInfo->size; r++)
//         //     {
//         //         MPI_Send(&num,1,MPI_INT,r,0,m_mpiInfo->comm);
//         //         MPI_Send(NormalNodesTmp.data(),3*num,MPI_DOUBLE,r,0,m_mpiInfo->comm);
//         //     }
//         // }
//         // else
//         // {
//         //     int num;
//         //     MPI_Recv(&num,1,MPI_INT,0,0,m_mpiInfo->comm,MPI_STATUS_IGNORE);
//         //     // ESYS_ASSERT(num==NormalNodesTmp.size(),"A point got lost (or added) somehow");
//         //     MPI_Recv(NormalNodesTmp.data(),3*num,MPI_DOUBLE,0,0,m_mpiInfo->comm,MPI_STATUS_IGNORE);
//         // }
//         // oxleytimer.toc("\t\t\t....done");

//         // // Write information into NodeIDs
//         oxleytimer.toc("Writing NodeIDs");
//         for(int i = 0; i < NormalNodesTmp.size(); i++)
//             NodeIDs[NormalNodesTmp[i]]=i;

//     }
//     else // being run on a single process
//     {
//         oxleytimer.toc("\t\tchecking for duplicates");
//         if(m_mpiInfo->rank == 0) // Redundant for mpi size = 1
//         {
//             int num_duplicates = 0;

//             //check for duplicates
//             std::vector<bool> locats(NormalNodesTmp.size(),false);
//             for(int i = 0; i < NormalNodesTmp.size(); i++)
//             {
//                 // skip points already tagged as duplicates
//                 if(locats[i]==true)
//                     continue;

//                 // get point
//                 auto point = NormalNodesTmp[i];
//                 bool duplicatePoint = hasDuplicate(point,NormalNodesTmp, false);
//                 if(duplicatePoint)
//                 {
//                     num_duplicates++;

//                     bool firstOccurance = true;
//                     for(int j = 0; j < NormalNodesTmp.size(); j++) // find location of duplicate points
//                     {
//                         if((std::abs(std::get<0>(NormalNodesTmp[j]) - std::get<0>(point)) < tuple_tolerance)
//                             && (std::abs(std::get<1>(NormalNodesTmp[j]) - std::get<1>(point)) < tuple_tolerance) 
//                                 && (std::abs(std::get<2>(NormalNodesTmp[j]) - std::get<2>(point)) < tuple_tolerance))
//                         {
//                             if(firstOccurance)
//                             {
//                                 firstOccurance = false;
//                                 continue;
//                             }
//                             else
//                             {
//                                 locats[j] = true;
//                                 break;
//                             }
//                         }
//                     }
//                 }
//             }

//             // erase duplicates
//             std::vector<DoubleTuple> NormalNodes_NoDuplicates(NormalNodesTmp.size()-num_duplicates,std::make_tuple(-1.,-1.,-1));
//             int nodecounter=0;
//             for(int j = 0; j < NormalNodesTmp.size(); j++)
//                 if(locats[j] == false)
//                 {
//                     NormalNodes_NoDuplicates[nodecounter]=NormalNodesTmp[j];
//                     nodecounter++;
//                 }
//             NormalNodesTmp=NormalNodes_NoDuplicates;

//             #ifdef OXLEY_ENABLE_DEBUG_RENUMBERNODES_INFORMATIONAL
//             std::cout << "Size after duplicates erased = " << NormalNodesTmp.size() << std::endl;
//             #endif

//             // Write information into NodeIDs
//             // oxleytimer.toc("Writing in " + std::to_string(NormalNodesTmp.size()) + " NodeIDs");
//             // NodeIDs.clear();
//             // for(int i = 0; i < NormalNodesTmp.size(); i++)
//             // {
//             //     NodeIDs[NormalNodesTmp[i]]=i;
//             // }
//             // oxleytimer.toc("'\tdone'");
//         }
//         oxleytimer.toc("\t\t\t...done");
//     }

//     // Send info to the other ranks
//     oxleytimer.toc("\tcommunicating");
//     MPI_Barrier(m_mpiInfo->comm); // have the other ranks wait until rank 0 arrives
//     if(m_mpiInfo->size > 1)
//     {
//         if(m_mpiInfo->rank == 0)
//         {
//             for(int i = 1; i < m_mpiInfo->size; i++)
//             {
//                 int numPoints = NormalNodesTmp.size();
//                 MPI_Send(&numPoints, 1, MPI_INT, i, 0, m_mpiInfo->comm);
//                 MPI_Send(NormalNodesTmp.data(), 3*NormalNodesTmp.size(), MPI_DOUBLE, i, 0, m_mpiInfo->comm);
//             }
//         }
//         else
//         {
//             int numPoints = 0;
//             MPI_Recv(&numPoints, 1, MPI_INT, 0, 0, m_mpiInfo->comm, MPI_STATUS_IGNORE);
//             NormalNodesTmp.resize(numPoints);
//             MPI_Recv(NormalNodesTmp.data(), 3*numPoints, MPI_DOUBLE, 0, 0, m_mpiInfo->comm, MPI_STATUS_IGNORE );

//             NodeIDs.clear();
//             for(int i = 0; i < NormalNodesTmp.size(); i++)
//                 NodeIDs[NormalNodesTmp[i]]=i;
//         }
//     }
//     oxleytimer.toc("\t\t\t...done");

// #else // no mpi
//     // Do the calculation
//     int numOfTrees=m_NE[0]*m_NE[1]*m_NE[2];
//     int numOfPoints=0;
//     int * increment = new int[numOfTrees+1];
//     increment[0]=0;
//     for(p8est_topidx_t treeid = p8est->first_local_tree; treeid <= p8est->last_local_tree; ++treeid) 
//     {
//         std::vector<std::tuple<double,double,double>> treevec;
//         p8est_tree_t * tree = p8est_tree_array_index(p8est->trees, treeid);
//         sc_array_t * tquadrants = &tree->quadrants;
//         p8est_locidx_t Q = (p8est_locidx_t) tquadrants->elem_count;
//         numOfPoints+=8*Q;
//         increment[treeid+1]=numOfPoints;
//     }

//     double * x = new double [numOfPoints];
//     double * y = new double [numOfPoints];
//     double * z = new double [numOfPoints];

//     oxleytimer.toc("\tdoing the calculation");
//     for(p8est_topidx_t treeid = p8est->first_local_tree; treeid <= p8est->last_local_tree; ++treeid) 
//     {
//         std::vector<std::tuple<double,double,double>> treevec;
//         p8est_tree_t * tree = p8est_tree_array_index(p8est->trees, treeid);
//         sc_array_t * tquadrants = &tree->quadrants;
//         p8est_locidx_t Q = (p8est_locidx_t) tquadrants->elem_count;
//         for(int q = 0; q < Q; ++q) { 
//             p8est_quadrant_t * oct = p8est_quadrant_array_index(tquadrants, q);
//             p8est_qcoord_t l = P8EST_QUADRANT_LEN(oct->level);
//             for(int n = 0; n < 8; n++) {
//                 double xyzB[3];
//                 p8est_qcoord_to_vertex(p8est->connectivity, treeid, oct->x+l*lxy_nodes[n][0], oct->y+l*lxy_nodes[n][1], oct->z+l*lxy_nodes[n][2], xyzB);

//                 x[increment[treeid]+8*q+n]=xyzB[0];
//                 y[increment[treeid]+8*q+n]=xyzB[1];
//                 z[increment[treeid]+8*q+n]=xyzB[2];
//             }
//         }
//     }
//     // Copy the data to NormalNodesTmp
//     oxleytimer.toc("\tcopying data");
//     NormalNodesTmp.assign(numOfPoints,std::make_tuple(-1.,-1.,-1.));
//     for(int i = 0; i < numOfPoints; i++)
//     {
//         NormalNodesTmp[i] = std::make_tuple(x[i],y[i],z[i]);
//     }
//     delete[] x;
//     delete[] y;
//     delete[] z;
//     delete[] increment;
//     oxleytimer.toc("\t\tdone");

//     //check for duplicates
//     oxleytimer.toc("\tchecking for duplicates");
//     bool locats[NormalNodesTmp.size()] = {false};
//     int num_duplicate_points = 0;
//     for(int i = 0; i < NormalNodesTmp.size(); i++)
//     {
//         // skip points already tagged as duplicates
//         if(locats[i] == true)
//             continue;

//         auto point = NormalNodesTmp[i];
//         bool duplicatePoint = hasDuplicate(point,NormalNodesTmp,false);
//         if(duplicatePoint)
//         {
//             num_duplicate_points++;
            
//             // tag duplicates
//             bool firstOccurance = true;
//             for(int j = 0; j < NormalNodesTmp.size(); j++) 
//             {
//                 if(     (std::abs(std::get<0>(NormalNodesTmp[j]) - std::get<0>(point)) < tuple_tolerance)
//                      && (std::abs(std::get<1>(NormalNodesTmp[j]) - std::get<1>(point)) < tuple_tolerance) 
//                      && (std::abs(std::get<2>(NormalNodesTmp[j]) - std::get<2>(point)) < tuple_tolerance))
//                 {
//                     if(firstOccurance)
//                     {
//                         firstOccurance=false;
//                         continue;
//                     }
//                     else
//                     {
//                         locats[j] = true;
//                     }
//                 }
//             }
//         }
//     }
    
//     // remove duplicates
//     std::vector<DoubleTuple> NormalNodes_NoDuplicates(NormalNodesTmp.size()-num_duplicate_points,std::make_tuple(-1.,-1.,-1.))
//     int duplicatecounter=0;
//     for(int i = 0; i < NormalNodesTmp.size(); i++)
//     {
//         if(locats[i] == false)
//         {
//             NormalNodes_NoDuplicates[duplicatecounter]=NormalNodesTmp[i];
//             duplicatecounter++;
//         }
//     }
//     NormalNodesTmp=NormalNodes_NoDuplicates;
    
//     oxleytimer.toc("\t\tdone");
//     // Write information into NodeIDs
//     for(int i = 0; i < NormalNodesTmp.size(); i++)
//         NodeIDs[NormalNodesTmp[i]]=i;
// #endif

// //     oxleytimer.toc("\tinformation loop");

// //     // Now work out the connections
// //     octantIDs.clear();

// #ifdef ESYS_MPI
//     long octantcounter = 0;
//     if(m_mpiInfo->rank == 0)
//     {
//         k=0;
//         // Do this process's nodes first
//         oxleytimer.toc("\tdoing rank 0 nodes");
//             for(p8est_topidx_t treeid = p8est->first_local_tree; treeid <= p8est->last_local_tree; ++treeid) {
//             p8est_tree_t * tree = p8est_tree_array_index(p8est->trees, treeid);
//             sc_array_t * tquadrants = &tree->quadrants;
//             p8est_locidx_t Q = (p8est_locidx_t) tquadrants->elem_count;

//             // Loop over octants
//             for(int q = 0; q < Q; ++q) 
//             { 
//                 p8est_quadrant_t * oct = p8est_quadrant_array_index(tquadrants, q);
//                 p8est_qcoord_t l = P8EST_QUADRANT_LEN(oct->level);

//                 // Save octant information
//                 double xyz[3];
//                 p8est_qcoord_to_vertex(p8est->connectivity, treeid, oct->x, oct->y, oct->z, xyz);
//                 octantIDs.push_back(octantcounter++);
//                 oct_info tmp;
//                 tmp.x=xyz[0];
//                 tmp.y=xyz[1];
//                 tmp.z=xyz[2];
//                 tmp.level=oct->level;
//                 octantInfo.push_back(tmp);

// //                 // Get the hanging edge and face information for this octant
// //                 int hanging_faces[6];
// //                 int hanging_edges[12];
// //                 getHangingInfo(nodes->face_code[k++], hanging_faces, hanging_edges);

// //                 #ifdef OXLEY_ENABLE_DEBUG_RENUMBER_NODES_EXTRA
// //                     std::cout << "Octant  " << k-1 << ":" << std::endl;
// //                     int faceCount=6, edgeCount=12;
// //                     for(int i=0;i<6;i++) faceCount+=hanging_faces[i];
// //                     for(int i=0;i<12;i++) edgeCount+=hanging_edges[i];
// //                     std::cout << "\tHanging faces : ";
// //                     if(faceCount==0)
// //                         std::cout << "[none]";
// //                     else
// //                         for(int i=0; i<6;i++)
// //                             if(hanging_faces[i]!=-1)
// //                                 std::cout << i << "(" << hanging_faces[i] << "), " ;
// //                     std::cout << std::endl;
// //                     std::cout << "\tHanging edges : ";
// //                     if(edgeCount==0)
// //                         std::cout << "[none]";
// //                     else
// //                         for(int i=0; i<12;i++)
// //                             if(hanging_edges[i]!=-1)
// //                                 std::cout << i << "(" << hanging_edges[i] << "), " ;

// //                     std::cout << std::endl;
// //                 #endif
    
// //                 // Record hanging node information
// //                 // Loop over faces
// //                 for(int n = 0; n < 6; n++)
// //                 {
// //                     if(hanging_faces[n] == 0)
// //                         continue;

// //                     // *******************************************************************
// //                     // Corner node
// //                     // *******************************************************************
// //                     double xyzC[3];
// //                     p8est_qcoord_to_vertex(p8est->connectivity, treeid, oct->x, oct->y, oct->z, xyzC);
// //                     auto cornerPoint = std::make_tuple(xyzC[0],xyzC[1],xyzC[2]);
// //                     long nodeidC = NodeIDs.find(cornerPoint)->second;

// //                     // *******************************************************************
// //                     // Face node
// //                     // *******************************************************************
// //                     // shift to the corner of the face
// //                     signed long x_corner = ((int) xface[n][hanging_faces[n]]) * l;
// //                     signed long y_corner = ((int) yface[n][hanging_faces[n]]) * l;
// //                     signed long z_corner = ((int) zface[n][hanging_faces[n]]) * l;

// //                     double xyz[3];
// //                     p8est_qcoord_to_vertex(p8est->connectivity, treeid, 
// //                                             oct->x+x_corner, 
// //                                             oct->y+y_corner, 
// //                                             oct->z+z_corner, xyz);

// //                     // Add to list of nodes
// //                     DoubleTuple point = std::make_tuple(xyz[0],xyz[1],xyz[2]);
// //                     long nodeid = NodeIDs.find(point)->second;
// //                     std::vector<long>::iterator have_node = std::find(HangingFaceNodesTmp.begin(), HangingFaceNodesTmp.end(),nodeid);
// //                     if(have_node == HangingFaceNodesTmp.end())
// //                           HangingFaceNodesTmp.push_back(nodeid);

// //                     // *******************************************************************
// //                     // Get the parent octant
// //                     // *******************************************************************
// //                     p8est_quadrant_t parentOct;
// //                     p8est_quadrant_t * parent = &parentOct;
// //                     p8est_quadrant_parent(oct, parent);
                    
// //                     // *******************************************************************
// //                     // Get the octant neighbouring this face
// //                     // *******************************************************************
// //                     // Get the neighbouring face
// //                     p8est_quadrant_t neighbourOct;
// //                     p8est_quadrant_t * neighbour = &neighbourOct;
// //                     // p8est_quadrant_face_neighbor(parent, n, neighbour);
// //                     int faceNeighbourTree=p8est_quadrant_face_neighbor_extra(parent,treeid,n,neighbour,NULL,connectivity);

// //                     // Save this information
// //                     hangingFaceInfo tmp;
// //                     tmp.x=oct->x;
// //                     tmp.y=oct->y;
// //                     tmp.z=oct->z;
// //                     tmp.level=oct->level;
// //                     tmp.treeid=treeid;
// //                     tmp.face_type=n;
// //                     tmp.neighbour_x=neighbour->x;
// //                     tmp.neighbour_y=neighbour->y;
// //                     tmp.neighbour_z=neighbour->z;
// //                     tmp.neighbour_level=neighbour->level;
// //                     tmp.neighbour_tree=faceNeighbourTree;
// //                     hanging_face_orientation.push_back(tmp);

// //                     // *******************************************************************
// //                     // Calculate the edge nodes
// //                     // *******************************************************************
// //                     signed long x_e_corner[2] = {((int) xedge[n][hanging_faces[n]][0]) * l,((int) xedge[n][hanging_faces[n]][1]) * l};
// //                     signed long y_e_corner[2] = {((int) yedge[n][hanging_faces[n]][0]) * l,((int) yedge[n][hanging_faces[n]][1]) * l};
// //                     signed long z_e_corner[2] = {((int) zedge[n][hanging_faces[n]][0]) * l,((int) zedge[n][hanging_faces[n]][1]) * l};
// //                     double xyzA[2][3];
// //                     p8est_qcoord_to_vertex(p8est->connectivity, treeid, oct->x+x_e_corner[0], oct->y+y_e_corner[0], oct->z+z_e_corner[0], xyzA[0]);
// //                     p8est_qcoord_to_vertex(p8est->connectivity, treeid, oct->x+x_e_corner[1], oct->y+y_e_corner[1], oct->z+z_e_corner[1], xyzA[1]);
// //                     DoubleTuple pointB[2] = {{std::make_tuple(xyzA[0][0],xyzA[0][1],xyzA[0][2])},{std::make_tuple(xyzA[1][0],xyzA[1][1],xyzA[1][2])}};
// //                     long nodeidB[2]= {NodeIDs.find(pointB[0])->second,NodeIDs.find(pointB[1])->second};

// //                     // Add to list of nodes
// //                     std::vector<long>::iterator have_edge;
// //                     have_edge = std::find(HangingEdgeNodesTmp.begin(), HangingEdgeNodesTmp.end(), nodeidB[0]);
// //                     if(have_edge == HangingEdgeNodesTmp.end())
// //                           HangingEdgeNodesTmp.push_back(nodeidB[0]);
// //                     have_edge = std::find(HangingEdgeNodesTmp.begin(), HangingEdgeNodesTmp.end(), nodeidB[1]);
// //                     if(have_edge == HangingEdgeNodesTmp.end())
// //                           HangingEdgeNodesTmp.push_back(nodeidB[1]);

// //                     // *******************************************************************
// //                     // Get the neighbouring edge(s)
// //                     // *******************************************************************
// //                     // Get the neighbours
// //                     p8est_quadrant_t EdgeNeighbourOct;
// //                     p8est_quadrant_t * EdgeNeighbour = &EdgeNeighbourOct;
// //                     int edge_number[2]={edge_lookup[n][hanging_faces[n]][0],edge_lookup[n][hanging_faces[n]][1]};
// //                     sc_array_t quads[2];
// //                     sc_array_t treeids[2];
// //                     sc_array_t nedges[2];
// //                     sc_array_t * pQuads[2] = {&quads[0],&quads[1]};
// //                     sc_array_t * pTreeids[2] = {&treeids[0],&treeids[1]};
// //                     sc_array_t * pNEdges[2] = {&nedges[0],&nedges[1]};
// //                     sc_array_init(pQuads[0], (size_t) sizeof(p4est_quadrant_t));
// //                     sc_array_init(pQuads[1], (size_t) sizeof(p4est_quadrant_t));
// //                     sc_array_init(pTreeids[0], (size_t) sizeof(p4est_topidx_t));
// //                     sc_array_init(pTreeids[1], (size_t) sizeof(p4est_topidx_t));
// //                     sc_array_init(pNEdges[0], (size_t) sizeof(int));
// //                     sc_array_init(pNEdges[1], (size_t) sizeof(int));

// //                     // Check for boundaries            
// //                     bool onBoundary[2];
// //                     onBoundary[0] = xyzA[0][0]<=0 || xyzA[0][0]>=forestData->m_length[0] ||
// //                                     xyzA[0][1]<=0 || xyzA[0][1]>=forestData->m_length[1] ||
// //                                     xyzA[0][2]<=0 || xyzA[0][2]>=forestData->m_length[2];
// //                     onBoundary[1] = xyzA[1][0]<=0 || xyzA[1][0]>=forestData->m_length[0] ||
// //                                     xyzA[1][1]<=0 || xyzA[1][1]>=forestData->m_length[1] ||
// //                                     xyzA[1][2]<=0 || xyzA[1][2]>=forestData->m_length[2];

// //                     hangingEdgeInfo temp;
// //                     if(!onBoundary[0])
// //                     {
// //                         p8est_quadrant_edge_neighbor_extra(parent,treeid,edge_number[0],pQuads[0],pTreeids[0],pNEdges[0],connectivity);
// //                         temp.nodeid=nodeidC;
// //                         temp.x=oct->x;
// //                         temp.y=oct->y;
// //                         temp.z=oct->z;
// //                         temp.level=oct->level;
// //                         temp.treeid=treeid;
// //                         temp.edge_type=n;
// //                         p8est_quadrant_t * tempB = (p8est_quadrant_t *) (*pQuads[0]).array;
// //                         temp.neighbour_x=tempB->x;
// //                         temp.neighbour_y=tempB->y;
// //                         temp.neighbour_z=tempB->z;
// //                         temp.neighbour_level=tempB->level;
// //                         p8est_topidx_t * tempC = (p8est_topidx_t *) (*pTreeids[0]).array;
// //                         temp.neighbour_tree=*tempC;
// //                         hanging_edge_orientation.push_back(temp);
// //                     }
// //                     if(!onBoundary[1])
// //                     {
// //                         p8est_quadrant_edge_neighbor_extra(parent,treeid,edge_number[1],pQuads[1],pTreeids[1],pNEdges[1],connectivity);
// //                         temp.nodeid=nodeidC;
// //                         temp.x=oct->x;
// //                         temp.y=oct->y;
// //                         temp.z=oct->z;
// //                         temp.level=oct->level;
// //                         temp.treeid=treeid;
// //                         temp.edge_type=n;
// //                         p8est_quadrant_t * tempB = (p8est_quadrant_t *) (*pQuads[1]).array;
// //                         temp.neighbour_x=tempB->x;
// //                         temp.neighbour_y=tempB->y;
// //                         temp.neighbour_z=tempB->z;
// //                         temp.neighbour_level=tempB->level;
// //                         p8est_topidx_t * tempC = (p8est_topidx_t *) (*pTreeids[1]).array;
// //                         temp.neighbour_tree=*tempC;
// //                         hanging_edge_orientation.push_back(temp);
// //                     }

// //                     // *******************************************************************
// //                     // Store connection information for updateRowsColumns
// //                     // *******************************************************************

// //                     signed long x_corner0 = xface0[n][hanging_faces[n]] * l;
// //                     signed long y_corner0 = yface0[n][hanging_faces[n]] * l;
// //                     signed long z_corner0 = zface0[n][hanging_faces[n]] * l;
                
// //                     double xyzD[3];
// //                     p8est_qcoord_to_vertex(p8est->connectivity, treeid, oct->x+x_corner0, oct->y+y_corner0, oct->z+z_corner0, xyzD);
// //                     auto otherCornerPoint = std::make_tuple(xyzD[0],xyzD[1],xyzD[2]);
// //                     long nodeA = NodeIDs.find(otherCornerPoint)->second;

// //                     p8est_qcoord_to_vertex(p8est->connectivity, treeid, oct->x+x_corner, oct->y+y_corner, oct->z+z_corner, xyzC);
// //                     auto facecornerPoint = std::make_tuple(xyzC[0],xyzC[1],xyzC[2]);
// //                     long facenodeidB = NodeIDs.find(facecornerPoint)->second;

// //                     long need_edgeA[4] = {std::count(hanging_edge_node_connections.begin(), hanging_edge_node_connections.end(), std::make_pair(nodeA,nodeidB[0])),
// //                                          std::count(hanging_edge_node_connections.begin(), hanging_edge_node_connections.end(), std::make_pair(nodeidB[0],nodeA)),
// //                                          std::count(hanging_edge_node_connections.begin(), hanging_edge_node_connections.end(), std::make_pair(nodeA,nodeidB[1])),
// //                                          std::count(hanging_edge_node_connections.begin(), hanging_edge_node_connections.end(), std::make_pair(nodeidB[1],nodeA))};
// //                     long need_edgeB[4] = {std::count(hanging_face_node_connections.begin(), hanging_face_node_connections.end(), std::make_pair(facenodeidB,nodeidB[0])),
// //                                          std::count(hanging_face_node_connections.begin(), hanging_face_node_connections.end(), std::make_pair(nodeidB[0],facenodeidB)),
// //                                          std::count(hanging_face_node_connections.begin(), hanging_face_node_connections.end(), std::make_pair(facenodeidB,nodeidB[1])),
// //                                          std::count(hanging_face_node_connections.begin(), hanging_face_node_connections.end(), std::make_pair(nodeidB[1],facenodeidB))};

                    
// //                     // std::cout << " index[" << n << "][" << hanging_faces[n] << "]" << std::endl;

// //                     if(need_edgeA[0]+need_edgeA[1] == 0)
// //                     {
// //                         // Calculate fallacious node connection that was generated by p8est
// //                         signed long x_corner1 = xface1[n][hanging_faces[n]] * l;
// //                         signed long y_corner1 = yface1[n][hanging_faces[n]] * l;
// //                         signed long z_corner1 = zface1[n][hanging_faces[n]] * l;
// //                         p8est_qcoord_to_vertex(p8est->connectivity, treeid, oct->x+x_corner1, oct->y+y_corner1, oct->z+z_corner1, xyzD);
// //                         auto otherCornerPoint1 = std::make_tuple(xyzD[0],xyzD[1],xyzD[2]);
// //                         long defunct1 = NodeIDs.find(otherCornerPoint1)->second;

// //                         // std::cout << " defunct connection1 = " << nodeA << ", " << nodeidB[0] << ", " << defunct1 << std::endl;

// //                         hanging_edge_node_connections.push_back(std::make_pair(nodeA,nodeidB[0]));
// //                         false_node_connections.push_back(std::make_pair(nodeA,defunct1));

// //                         hanging_edge_node_connections.push_back(std::make_pair(nodeidB[0],defunct1));
// //                         false_node_connections.push_back(std::make_pair(defunct1,nodeA));
// //                     }
// //                     if(need_edgeA[2]+need_edgeA[3] == 0)
// //                     {
// //                         // Calculate fallacious node connection that was generated by p8est
// //                         signed long x_corner2 = xface2[n][hanging_faces[n]] * l;
// //                         signed long y_corner2 = yface2[n][hanging_faces[n]] * l;
// //                         signed long z_corner2 = zface2[n][hanging_faces[n]] * l;
// //                         p8est_qcoord_to_vertex(p8est->connectivity, treeid, oct->x+x_corner2, oct->y+y_corner2, oct->z+z_corner2, xyzD);
// //                         auto otherCornerPoint2 = std::make_tuple(xyzD[0],xyzD[1],xyzD[2]);
// //                         long defunct2 = NodeIDs.find(otherCornerPoint2)->second;

// //                         // std::cout << " defunct connection2 = " << nodeA << ", " << nodeidB[1] << ", " << defunct2 << std::endl;

// //                         hanging_edge_node_connections.push_back(std::make_pair(nodeA,nodeidB[1]));
// //                         false_node_connections.push_back(std::make_pair(nodeA,defunct2));

// //                         hanging_edge_node_connections.push_back(std::make_pair(nodeidB[1],defunct2));
// //                         false_node_connections.push_back(std::make_pair(defunct2,nodeA));
// //                     }
// //                     if(need_edgeB[0]+need_edgeB[1] == 0)
// //                     {
// //                         hanging_face_node_connections.push_back(std::make_pair(nodeidB[0],facenodeidB));
// //                     }
// //                     if(need_edgeB[2]+need_edgeB[3] == 0)
// //                     {
// //                         hanging_face_node_connections.push_back(std::make_pair(nodeidB[1],facenodeidB));
// //                     }


// //                     // *******************************************************************
// //                     // Debuging output
// //                     // *******************************************************************
// //                     #ifdef OXLEY_ENABLE_DEBUG_RENUMBER_NODES_EXTRA
// //                         std::cout << "\t\toctant = \033[1;36m" << k-1 << "\033[0m, corner node = " << nodeidC << std::endl;
// //                         std::cout << "\t\tface node = \033[1;36m" << facenodeidB << "\033[0m" << std::endl;
// //                         std::cout << "\t\t    index [" << n << "][" << hanging_faces[n] << "]" << std::endl;
// //                         // std::cout << "\t\t    shift was (Dx,Dy,Dz)= (" << x_corner << ", " << y_corner << ", " << z_corner << ")" << std::endl;
// //                         std::cout << "\t\tedge nodes: \033[1;36m" << nodeidB[0] << " & " << nodeidB[1] << "\033[0m" << std::endl;
// //                         std::cout << "\t\t    index[" << (int) n << "][" << (int) hanging_faces[n] << "][-]" << std::endl;
// //                         // std::cout << "\t\t    shift was (Dx,Dy,Dz)= (" << x_e_corner[0] << " ," << y_e_corner[0] << ", " << z_e_corner[0] << ")" << std::endl;
// //                         // std::cout << "\t\t    shift was (Dx,Dy,Dz)= (" << x_e_corner[1] << " ," << y_e_corner[1] << ", " << z_e_corner[1] << ")" << std::endl;
// //                         std::cout << "\t\tcorner nodes: \033[1;36m" << nodeA << "\033[0m" << std::endl;
// //                         std::cout << "\t\t    index [" << n << "][" << hanging_faces[n] << "][-]" << std::endl;
// //                         std::cout << "\t\tFace nodes: \033[1;36m" << facenodeidB << "--" << nodeidB[1] << "--" << nodeA << "--" << nodeidB[0] << "\033[0m" << std::endl;
// //                         std::cout << "\t\tedge numbers: \033[1;36m" << edge_number[0] << " & " << edge_number[1] << "\033[0m" << std::endl;
// //                         std::cout << "\t\t    index [" << n << "][" << hanging_faces[n] << "][-]" << std::endl;
// //                     #endif
// //                 }
//             }
//         }
//     }

//     MPI_Barrier(m_mpiInfo->comm);

//     std::vector<int> treeidVec;
//     std::vector<p8est_quadrant_t> octantVec;
//     int num_of_iterations=0;
//     if(m_mpiInfo->size > 1)
//     {
//         for(p8est_topidx_t treeid = p8est->first_local_tree; treeid <= p8est->last_local_tree; ++treeid) 
//         {
//             p8est_tree_t * tree = p8est_tree_array_index(p8est->trees, treeid);
//             sc_array_t * tquadrants = &tree->quadrants;
//             p8est_locidx_t Q = (p8est_locidx_t) tquadrants->elem_count;
//             num_of_iterations+= ((int) Q);
//         }

//         octantVec.resize(num_of_iterations);
//         treeidVec.resize(num_of_iterations);
//         int count = 0;
//         for(p8est_topidx_t treeid = p8est->first_local_tree; treeid <= p8est->last_local_tree; ++treeid) 
//         {
//             p8est_tree_t * tree = p8est_tree_array_index(p8est->trees, treeid);
//             sc_array_t * tquadrants = &tree->quadrants;
//             p8est_locidx_t Q = (p8est_locidx_t) tquadrants->elem_count;
//             for(int q=0; q<Q; q++)
//             {
//                 p8est_quadrant_t * octant = p8est_quadrant_array_index(tquadrants, q);
//                 octantVec[count] = *octant;
//                 treeidVec[count++]=treeid;
//             }
//         }
//     }

//     if(m_mpiInfo->size > 1)
//         MPI_Barrier(m_mpiInfo->comm);

//     int num_iter[m_mpiInfo->size] = {0};
//     if(m_mpiInfo->size > 1)
//     {
//         if(m_mpiInfo->rank == 0)
//         {
//             for(int r = 1; r < m_mpiInfo->size; r++)
//                 MPI_Recv(&num_iter[r], 1, MPI_INT, r, 0, m_mpiInfo->comm, MPI_STATUS_IGNORE);
//         }
//         else
//         {
//             MPI_Send(&num_of_iterations, 1, MPI_INT, 0, 0, m_mpiInfo->comm);
//         }
//     }

//     if(m_mpiInfo->size > 1)
//         MPI_Barrier(m_mpiInfo->comm);

//     if(m_mpiInfo->size > 1)
//     {
//         if(m_mpiInfo->rank == 0)
//         {
//             for(int r = 1; r < m_mpiInfo->size; r++)
//                 MPI_Send(&num_iter[0], m_mpiInfo->size, MPI_INT, r, 0, m_mpiInfo->comm);
//         }
//         else
//         {
//             MPI_Recv(&num_iter[0], m_mpiInfo->size, MPI_INT, 0, 0, m_mpiInfo->comm, MPI_STATUS_IGNORE);
//         }
//     }

//     if(m_mpiInfo->size > 1)
//         MPI_Barrier(m_mpiInfo->comm);

//     int total_iter = 0;
//     for(int i = 0; i < m_mpiInfo->size ; i++)
//         total_iter += num_iter[i];

//     if(m_mpiInfo->size > 1)
//     {
//         if(m_mpiInfo->rank == 0)
//         {
//             for(int r = 1; r < m_mpiInfo->size; r++)
//                 MPI_Send(&octantcounter,1,MPI_LONG,r,0,m_mpiInfo->comm);
//         }
//         else
//         {
//             MPI_Recv(&octantcounter,1,MPI_LONG,0,0,m_mpiInfo->comm,MPI_STATUS_IGNORE);
//         }


//         oxleytimer.toc("\tDoing other ranks...");
//         for(int r = 1; r < m_mpiInfo->size; r++)
//         {
//             for(int j = 0; j < num_iter[r]; j++) 
//             {
//                 MPI_Barrier(m_mpiInfo->comm);
//                 if(m_mpiInfo->rank == r)
//                 {
//                     int treeid = treeidVec[j];
//                     p8est_quadrant_t oct = octantVec[j];
//                     MPI_Send(&oct.x, 1, MPI_INT, 0, 0, m_mpiInfo->comm);
//                     MPI_Send(&oct.y, 1, MPI_INT, 0, 0, m_mpiInfo->comm);
//                     MPI_Send(&oct.z, 1, MPI_INT, 0, 0, m_mpiInfo->comm);
//                     MPI_Send(&oct.level, 1, MPI_INT8_T, 0, 0, m_mpiInfo->comm);
//                     MPI_Send(&oct.pad8, 1, MPI_INT8_T, 0, 0, m_mpiInfo->comm);
//                     MPI_Send(&oct.pad16, 1, MPI_INT16_T, 0, 0, m_mpiInfo->comm);
//                     MPI_Send(&treeid, 1, MPI_INT, 0, 0, m_mpiInfo->comm);
//                     int hanging_code=nodes->face_code[k++];
//                     MPI_Send(&hanging_code, 1, MPI_INT, 0, 0, m_mpiInfo->comm);
//                 }
//                 else if(m_mpiInfo->rank == 0)
//                 {
//                     p8est_quadrant_t * oct;
//                     p8est_quadrant_t octant;
//                     int x, y, z, treeid;
//                     MPI_Recv(&octant.x, 1, MPI_INT, r, 0, m_mpiInfo->comm, MPI_STATUS_IGNORE);
//                     MPI_Recv(&octant.y, 1, MPI_INT, r, 0, m_mpiInfo->comm, MPI_STATUS_IGNORE);
//                     MPI_Recv(&octant.z, 1, MPI_INT, r, 0, m_mpiInfo->comm, MPI_STATUS_IGNORE);
//                     MPI_Recv(&octant.level, 1, MPI_INT8_T, r, 0, m_mpiInfo->comm, MPI_STATUS_IGNORE);
//                     MPI_Recv(&octant.pad8, 1, MPI_INT8_T, r, 0, m_mpiInfo->comm, MPI_STATUS_IGNORE);
//                     MPI_Recv(&octant.pad16, 1, MPI_INT16_T, r, 0, m_mpiInfo->comm, MPI_STATUS_IGNORE);
//                     MPI_Recv(&treeid, 1, MPI_INT, r, 0, m_mpiInfo->comm, MPI_STATUS_IGNORE);
//                     oct = &octant;
//                     int hanging_code;
//                     MPI_Recv(&hanging_code, 1, MPI_INT, r, 0, m_mpiInfo->comm, MPI_STATUS_IGNORE);

//                     // Save octant information
//                     p8est_qcoord_t l = P8EST_QUADRANT_LEN(oct->level);
//                     double xyz[3];
//                     p8est_qcoord_to_vertex(p8est->connectivity, treeid, oct->x, oct->y, oct->z, xyz);
//                     octantIDs.push_back(octantcounter++);
//                     // oct_info tmp;
//                     // tmp.x=xyz[0];
//                     // tmp.y=xyz[1];
//                     // tmp.z=xyz[2];
//                     // tmp.level=oct->level;
//                     // octantInfo.push_back(tmp);

//                     // // Get the hanging edge and face information for this octant
//                     // int hanging_faces[6];
//                     // int hanging_edges[12];
//                     // // getHangingInfo(nodes->face_code[k++], hanging_faces, hanging_edges);
//                     // getHangingInfo(hanging_code, hanging_faces, hanging_edges);

//                     // #ifdef OXLEY_ENABLE_DEBUG_RENUMBER_NODES_EXTRA
//                     //     std::cout << "Octant  " << k-1 << ":" << std::endl;
//                     //     int faceCount=6, edgeCount=12;
//                     //     for(int i=0;i<6;i++) faceCount+=hanging_faces[i];
//                     //     for(int i=0;i<12;i++) edgeCount+=hanging_edges[i];
//                     //     std::cout << "\tHanging faces : ";
//                     //     if(faceCount==0)
//                     //         std::cout << "[none]";
//                     //     else
//                     //         for(int i=0; i<6;i++)
//                     //             if(hanging_faces[i]!=-1)
//                     //                 std::cout << i << "(" << hanging_faces[i] << "), " ;
//                     //     std::cout << std::endl;
//                     //     std::cout << "\tHanging edges : ";
//                     //     if(edgeCount==0)
//                     //         std::cout << "[none]";
//                     //     else
//                     //         for(int i=0; i<12;i++)
//                     //             if(hanging_edges[i]!=-1)
//                     //                 std::cout << i << "(" << hanging_edges[i] << "), " ;

//                     //     std::cout << std::endl;
//                     // #endif
                
//                     // Record hanging node information
//                     // Loop over faces
//                     // for(int n = 0; n < 6; n++)
//                     // {
//                     //     // If the face does not have any hanging nodes then skip ahead
//                     //     if(hanging_faces[n]==-1)
//                     //         continue;

//                     //     // *******************************************************************
//                     //     // Corner node
//                     //     // *******************************************************************
//                     //     double xyzC[3];
//                     //     p8est_qcoord_to_vertex(p8est->connectivity, treeid, oct->x, oct->y, oct->z, xyzC);
//                     //     auto cornerPoint = std::make_tuple(xyzC[0],xyzC[1],xyzC[2]);
//                     //     long nodeidC = NodeIDs.find(cornerPoint)->second;

//                     //     // *******************************************************************
//                     //     // Face node
//                     //     // *******************************************************************
//                     //     // shift to the corner of the face
//                     //     signed long x_corner = ((int) xface[n][hanging_faces[n]]) * l;
//                     //     signed long y_corner = ((int) yface[n][hanging_faces[n]]) * l;
//                     //     signed long z_corner = ((int) zface[n][hanging_faces[n]]) * l;

//                     //     double xyz[3];
//                     //     p8est_qcoord_to_vertex(p8est->connectivity, treeid, 
//                     //                             oct->x+x_corner, 
//                     //                             oct->y+y_corner, 
//                     //                             oct->z+z_corner, xyz);

//                     //     // Add to list of nodes
//                     //     DoubleTuple point = std::make_tuple(xyz[0],xyz[1],xyz[2]);
//                     //     long nodeid = NodeIDs.find(point)->second;
//                     //     std::vector<long>::iterator have_node = std::find(HangingFaceNodesTmp.begin(), HangingFaceNodesTmp.end(),nodeid);
//                     //     if(have_node == HangingFaceNodesTmp.end())
//                     //           HangingFaceNodesTmp.push_back(nodeid);

//                     //     // *******************************************************************
//                     //     // Get the parent octant
//                     //     // *******************************************************************
//                     //     p8est_quadrant_t parentOct;
//                     //     p8est_quadrant_t * parent = &parentOct;
//                     //     p8est_quadrant_parent(oct, parent);
                        
//                     //     // *******************************************************************
//                     //     // Get the octant neighbouring this face
//                     //     // *******************************************************************
//                     //     // Get the neighbouring face
//                     //     p8est_quadrant_t neighbourOct;
//                     //     p8est_quadrant_t * neighbour = &neighbourOct;
//                     //     // p8est_quadrant_face_neighbor(parent, n, neighbour);
//                     //     int faceNeighbourTree=p8est_quadrant_face_neighbor_extra(parent,treeid,n,neighbour,NULL,connectivity);

//                     //     // Save this information
//                     //     hangingFaceInfo tmp;
//                     //     tmp.x=oct->x;
//                     //     tmp.y=oct->y;
//                     //     tmp.z=oct->z;
//                     //     tmp.level=oct->level;
//                     //     tmp.treeid=treeid;
//                     //     tmp.face_type=n;
//                     //     tmp.neighbour_x=neighbour->x;
//                     //     tmp.neighbour_y=neighbour->y;
//                     //     tmp.neighbour_z=neighbour->z;
//                     //     tmp.neighbour_level=neighbour->level;
//                     //     tmp.neighbour_tree=faceNeighbourTree;
//                     //     hanging_face_orientation.push_back(tmp);

//                     //     // *******************************************************************
//                     //     // Calculate the edge nodes
//                     //     // *******************************************************************
//                     //     signed long x_e_corner[2] = {((int) xedge[n][hanging_faces[n]][0]) * l,((int) xedge[n][hanging_faces[n]][1]) * l};
//                     //     signed long y_e_corner[2] = {((int) yedge[n][hanging_faces[n]][0]) * l,((int) yedge[n][hanging_faces[n]][1]) * l};
//                     //     signed long z_e_corner[2] = {((int) zedge[n][hanging_faces[n]][0]) * l,((int) zedge[n][hanging_faces[n]][1]) * l};
//                     //     double xyzA[2][3];
//                     //     p8est_qcoord_to_vertex(p8est->connectivity, treeid, oct->x+x_e_corner[0], oct->y+y_e_corner[0], oct->z+z_e_corner[0], xyzA[0]);
//                     //     p8est_qcoord_to_vertex(p8est->connectivity, treeid, oct->x+x_e_corner[1], oct->y+y_e_corner[1], oct->z+z_e_corner[1], xyzA[1]);
//                     //     DoubleTuple pointB[2] = {{std::make_tuple(xyzA[0][0],xyzA[0][1],xyzA[0][2])},{std::make_tuple(xyzA[1][0],xyzA[1][1],xyzA[1][2])}};
//                     //     long nodeidB[2]= {NodeIDs.find(pointB[0])->second,NodeIDs.find(pointB[1])->second};

//                     //     // Add to list of nodes
//                     //     std::vector<long>::iterator have_edge;
//                     //     have_edge = std::find(HangingEdgeNodesTmp.begin(), HangingEdgeNodesTmp.end(), nodeidB[0]);
//                     //     if(have_edge == HangingEdgeNodesTmp.end())
//                     //           HangingEdgeNodesTmp.push_back(nodeidB[0]);
//                     //     have_edge = std::find(HangingEdgeNodesTmp.begin(), HangingEdgeNodesTmp.end(), nodeidB[1]);
//                     //     if(have_edge == HangingEdgeNodesTmp.end())
//                     //           HangingEdgeNodesTmp.push_back(nodeidB[1]);

//                     //     // *******************************************************************
//                     //     // Get the neighbouring edge(s)
//                     //     // *******************************************************************
//                     //     // Get the neighbours
//                     //     p8est_quadrant_t EdgeNeighbourOct;
//                     //     p8est_quadrant_t * EdgeNeighbour = &EdgeNeighbourOct;
//                     //     int edge_number[2]={edge_lookup[n][hanging_faces[n]][0],edge_lookup[n][hanging_faces[n]][1]};
//                     //     sc_array_t quads[2];
//                     //     sc_array_t treeids[2];
//                     //     sc_array_t nedges[2];
//                     //     sc_array_t * pQuads[2] = {&quads[0],&quads[1]};
//                     //     sc_array_t * pTreeids[2] = {&treeids[0],&treeids[1]};
//                     //     sc_array_t * pNEdges[2] = {&nedges[0],&nedges[1]};
//                     //     sc_array_init(pQuads[0], (size_t) sizeof(p4est_quadrant_t));
//                     //     sc_array_init(pQuads[1], (size_t) sizeof(p4est_quadrant_t));
//                     //     sc_array_init(pTreeids[0], (size_t) sizeof(p4est_topidx_t));
//                     //     sc_array_init(pTreeids[1], (size_t) sizeof(p4est_topidx_t));
//                     //     sc_array_init(pNEdges[0], (size_t) sizeof(int));
//                     //     sc_array_init(pNEdges[1], (size_t) sizeof(int));

//                     //     // Check for boundaries            
//                     //     bool onBoundary[2];
//                     //     onBoundary[0] = xyzA[0][0]<=0 || xyzA[0][0]>=forestData->m_length[0] ||
//                     //                     xyzA[0][1]<=0 || xyzA[0][1]>=forestData->m_length[1] ||
//                     //                     xyzA[0][2]<=0 || xyzA[0][2]>=forestData->m_length[2];
//                     //     onBoundary[1] = xyzA[1][0]<=0 || xyzA[1][0]>=forestData->m_length[0] ||
//                     //                     xyzA[1][1]<=0 || xyzA[1][1]>=forestData->m_length[1] ||
//                     //                     xyzA[1][2]<=0 || xyzA[1][2]>=forestData->m_length[2];

//                     //     hangingEdgeInfo temp;
//                     //     if(!onBoundary[0])
//                     //     {
//                     //         p8est_quadrant_edge_neighbor_extra(parent,treeid,edge_number[0],pQuads[0],pTreeids[0],pNEdges[0],connectivity);
//                     //         temp.nodeid=nodeidC;
//                     //         temp.x=oct->x;
//                     //         temp.y=oct->y;
//                     //         temp.z=oct->z;
//                     //         temp.level=oct->level;
//                     //         temp.treeid=treeid;
//                     //         temp.edge_type=n;
//                     //         p8est_quadrant_t * tempB = (p8est_quadrant_t *) (*pQuads[0]).array;
//                     //         temp.neighbour_x=tempB->x;
//                     //         temp.neighbour_y=tempB->y;
//                     //         temp.neighbour_z=tempB->z;
//                     //         temp.neighbour_level=tempB->level;
//                     //         p8est_topidx_t * tempC = (p8est_topidx_t *) (*pTreeids[0]).array;
//                     //         temp.neighbour_tree=*tempC;
//                     //         hanging_edge_orientation.push_back(temp);
//                     //     }
//                     //     if(!onBoundary[1])
//                     //     {
//                     //         p8est_quadrant_edge_neighbor_extra(parent,treeid,edge_number[1],pQuads[1],pTreeids[1],pNEdges[1],connectivity);
//                     //         temp.nodeid=nodeidC;
//                     //         temp.x=oct->x;
//                     //         temp.y=oct->y;
//                     //         temp.z=oct->z;
//                     //         temp.level=oct->level;
//                     //         temp.treeid=treeid;
//                     //         temp.edge_type=n;
//                     //         p8est_quadrant_t * tempB = (p8est_quadrant_t *) (*pQuads[1]).array;
//                     //         temp.neighbour_x=tempB->x;
//                     //         temp.neighbour_y=tempB->y;
//                     //         temp.neighbour_z=tempB->z;
//                     //         temp.neighbour_level=tempB->level;
//                     //         p8est_topidx_t * tempC = (p8est_topidx_t *) (*pTreeids[1]).array;
//                     //         temp.neighbour_tree=*tempC;
//                     //         hanging_edge_orientation.push_back(temp);
//                     //     }

//                     //     // *******************************************************************
//                     //     // Store connection information for updateRowsColumns
//                     //     // *******************************************************************

//                     //     signed long x_corner0 = xface0[n][hanging_faces[n]] * l;
//                     //     signed long y_corner0 = yface0[n][hanging_faces[n]] * l;
//                     //     signed long z_corner0 = zface0[n][hanging_faces[n]] * l;
                    
//                     //     double xyzD[3];
//                     //     p8est_qcoord_to_vertex(p8est->connectivity, treeid, oct->x+x_corner0, oct->y+y_corner0, oct->z+z_corner0, xyzD);
//                     //     auto otherCornerPoint = std::make_tuple(xyzD[0],xyzD[1],xyzD[2]);
//                     //     long nodeA = NodeIDs.find(otherCornerPoint)->second;

//                     //     p8est_qcoord_to_vertex(p8est->connectivity, treeid, oct->x+x_corner, oct->y+y_corner, oct->z+z_corner, xyzC);
//                     //     auto facecornerPoint = std::make_tuple(xyzC[0],xyzC[1],xyzC[2]);
//                     //     long facenodeidB = NodeIDs.find(facecornerPoint)->second;

//                     //     long need_edgeA[4] = {std::count(hanging_edge_node_connections.begin(), hanging_edge_node_connections.end(), std::make_pair(nodeA,nodeidB[0])),
//                     //                          std::count(hanging_edge_node_connections.begin(), hanging_edge_node_connections.end(), std::make_pair(nodeidB[0],nodeA)),
//                     //                          std::count(hanging_edge_node_connections.begin(), hanging_edge_node_connections.end(), std::make_pair(nodeA,nodeidB[1])),
//                     //                          std::count(hanging_edge_node_connections.begin(), hanging_edge_node_connections.end(), std::make_pair(nodeidB[1],nodeA))};
//                     //     long need_edgeB[4] = {std::count(hanging_face_node_connections.begin(), hanging_face_node_connections.end(), std::make_pair(facenodeidB,nodeidB[0])),
//                     //                          std::count(hanging_face_node_connections.begin(), hanging_face_node_connections.end(), std::make_pair(nodeidB[0],facenodeidB)),
//                     //                          std::count(hanging_face_node_connections.begin(), hanging_face_node_connections.end(), std::make_pair(facenodeidB,nodeidB[1])),
//                     //                          std::count(hanging_face_node_connections.begin(), hanging_face_node_connections.end(), std::make_pair(nodeidB[1],facenodeidB))};

                        
//                     //     // std::cout << " index[" << n << "][" << hanging_faces[n] << "]" << std::endl;

//                     //     if(need_edgeA[0]+need_edgeA[1] == 0)
//                     //     {
//                     //         // Calculate fallacious node connection that was generated by p8est
//                     //         signed long x_corner1 = xface1[n][hanging_faces[n]] * l;
//                     //         signed long y_corner1 = yface1[n][hanging_faces[n]] * l;
//                     //         signed long z_corner1 = zface1[n][hanging_faces[n]] * l;
//                     //         p8est_qcoord_to_vertex(p8est->connectivity, treeid, oct->x+x_corner1, oct->y+y_corner1, oct->z+z_corner1, xyzD);
//                     //         auto otherCornerPoint1 = std::make_tuple(xyzD[0],xyzD[1],xyzD[2]);
//                     //         long defunct1 = NodeIDs.find(otherCornerPoint1)->second;

//                     //         // std::cout << " defunct connection1 = " << nodeA << ", " << nodeidB[0] << ", " << defunct1 << std::endl;

//                     //         hanging_edge_node_connections.push_back(std::make_pair(nodeA,nodeidB[0]));
//                     //         false_node_connections.push_back(std::make_pair(nodeA,defunct1));

//                     //         hanging_edge_node_connections.push_back(std::make_pair(nodeidB[0],defunct1));
//                     //         false_node_connections.push_back(std::make_pair(defunct1,nodeA));
//                     //     }
//                     //     if(need_edgeA[2]+need_edgeA[3] == 0)
//                     //     {
//                     //         // Calculate fallacious node connection that was generated by p8est
//                     //         signed long x_corner2 = xface2[n][hanging_faces[n]] * l;
//                     //         signed long y_corner2 = yface2[n][hanging_faces[n]] * l;
//                     //         signed long z_corner2 = zface2[n][hanging_faces[n]] * l;
//                     //         p8est_qcoord_to_vertex(p8est->connectivity, treeid, oct->x+x_corner2, oct->y+y_corner2, oct->z+z_corner2, xyzD);
//                     //         auto otherCornerPoint2 = std::make_tuple(xyzD[0],xyzD[1],xyzD[2]);
//                     //         long defunct2 = NodeIDs.find(otherCornerPoint2)->second;

//                     //         // std::cout << " defunct connection2 = " << nodeA << ", " << nodeidB[1] << ", " << defunct2 << std::endl;

//                     //         hanging_edge_node_connections.push_back(std::make_pair(nodeA,nodeidB[1]));
//                     //         false_node_connections.push_back(std::make_pair(nodeA,defunct2));

//                     //         hanging_edge_node_connections.push_back(std::make_pair(nodeidB[1],defunct2));
//                     //         false_node_connections.push_back(std::make_pair(defunct2,nodeA));
//                     //     }
//                     //     if(need_edgeB[0]+need_edgeB[1] == 0)
//                     //     {
//                     //         hanging_face_node_connections.push_back(std::make_pair(nodeidB[0],facenodeidB));
//                     //     }
//                     //     if(need_edgeB[2]+need_edgeB[3] == 0)
//                     //     {
//                     //         hanging_face_node_connections.push_back(std::make_pair(nodeidB[1],facenodeidB));
//                     //     }


//                     //     // *******************************************************************
//                     //     // Debuging output
//                     //     // *******************************************************************
//                     //     #ifdef OXLEY_ENABLE_DEBUG_RENUMBER_NODES_EXTRA
//                     //         std::cout << "\t\toctant = \033[1;36m" << k-1 << "\033[0m, corner node = " << nodeidC << std::endl;
//                     //         std::cout << "\t\tface node = \033[1;36m" << facenodeidB << "\033[0m" << std::endl;
//                     //         std::cout << "\t\t    index [" << n << "][" << hanging_faces[n] << "]" << std::endl;
//                     //         // std::cout << "\t\t    shift was (Dx,Dy,Dz)= (" << x_corner << ", " << y_corner << ", " << z_corner << ")" << std::endl;
//                     //         std::cout << "\t\tedge nodes: \033[1;36m" << nodeidB[0] << " & " << nodeidB[1] << "\033[0m" << std::endl;
//                     //         std::cout << "\t\t    index[" << (int) n << "][" << (int) hanging_faces[n] << "][-]" << std::endl;
//                     //         // std::cout << "\t\t    shift was (Dx,Dy,Dz)= (" << x_e_corner[0] << " ," << y_e_corner[0] << ", " << z_e_corner[0] << ")" << std::endl;
//                     //         // std::cout << "\t\t    shift was (Dx,Dy,Dz)= (" << x_e_corner[1] << " ," << y_e_corner[1] << ", " << z_e_corner[1] << ")" << std::endl;
//                     //         std::cout << "\t\tcorner nodes: \033[1;36m" << nodeA << "\033[0m" << std::endl;
//                     //         std::cout << "\t\t    index [" << n << "][" << hanging_faces[n] << "][-]" << std::endl;
//                     //         std::cout << "\t\tFace nodes: \033[1;36m" << facenodeidB << "--" << nodeidB[1] << "--" << nodeA << "--" << nodeidB[0] << "\033[0m" << std::endl;
//                     //         std::cout << "\t\tedge numbers: \033[1;36m" << edge_number[0] << " & " << edge_number[1] << "\033[0m" << std::endl;
//                     //         std::cout << "\t\t    index [" << n << "][" << hanging_faces[n] << "][-]" << std::endl;
//                     //     #endif
//                     // }
//                 }
//             }
//         }
//     }


//     MPI_Barrier(m_mpiInfo->comm);
//     oxleytimer.toc("\tcommunicating");
//     if(m_mpiInfo->size > 1)
//     {
//         if(m_mpiInfo->rank == 0)
//         {
//             for(int i = 1; i < m_mpiInfo->size; i++)
//             {
//                 int num;
//                 num = octantIDs.size();
//                 MPI_Send(&num, 1, MPI_INT, i, 0, m_mpiInfo->comm);
//                 MPI_Send(octantIDs.data(), num, MPI_LONG, i, 0, m_mpiInfo->comm);

//                 num=hanging_face_orientation.size();
//                 MPI_Send(&num, 1, MPI_INT, i, 0, m_mpiInfo->comm);
//                 for(int j = 0; j < hanging_face_orientation.size(); j++ )
//                 {
//                     hangingFaceInfo t = hanging_face_orientation[j];
//                     MPI_Send(&t.x, 1, MPI_LONG, i, 0, m_mpiInfo->comm);
//                     MPI_Send(&t.y, 1, MPI_LONG, i, 0, m_mpiInfo->comm);
//                     MPI_Send(&t.z, 1, MPI_LONG, i, 0, m_mpiInfo->comm);
//                     MPI_Send(&t.level, 1, MPI_INT8_T, i, 0, m_mpiInfo->comm);
//                     MPI_Send(&t.treeid, 1, MPI_INT, i, 0, m_mpiInfo->comm);
//                     MPI_Send(&t.face_type, 1, MPI_INT8_T, i, 0, m_mpiInfo->comm);
//                     MPI_Send(&t.neighbour_x, 1, MPI_LONG, i, 0, m_mpiInfo->comm);
//                     MPI_Send(&t.neighbour_y, 1, MPI_LONG, i, 0, m_mpiInfo->comm);
//                     MPI_Send(&t.neighbour_z, 1, MPI_LONG, i, 0, m_mpiInfo->comm);
//                     MPI_Send(&t.neighbour_level, 1, MPI_INT, i, 0, m_mpiInfo->comm);
//                     MPI_Send(&t.neighbour_tree, 1, MPI_INT, i, 0, m_mpiInfo->comm);
//                 }

//                 num=hanging_face_node_connections.size();
//                 MPI_Send(&num, 1, MPI_INT, i, 0, m_mpiInfo->comm);
//                 MPI_Send(hanging_face_node_connections.data(), 2*num, MPI_LONG, i, 0, m_mpiInfo->comm);
                
//                 num=hanging_edge_orientation.size();
//                 MPI_Send(&num, 1, MPI_INT, i, 0, m_mpiInfo->comm);
//                 for(int j = 0; j < hanging_edge_orientation.size(); j++ )
//                 {
//                     hangingEdgeInfo t = hanging_edge_orientation[j];
//                     MPI_Send(&t.nodeid, 1, MPI_LONG, i, 0, m_mpiInfo->comm);
//                     MPI_Send(&t.x, 1, MPI_LONG, i, 0, m_mpiInfo->comm);
//                     MPI_Send(&t.y, 1, MPI_LONG, i, 0, m_mpiInfo->comm);
//                     MPI_Send(&t.z, 1, MPI_LONG, i, 0, m_mpiInfo->comm);
//                     MPI_Send(&t.level, 1, MPI_INT8_T, i, 0, m_mpiInfo->comm);
//                     MPI_Send(&t.treeid, 1, MPI_INT, i, 0, m_mpiInfo->comm);
//                     MPI_Send(&t.edge_type, 1, MPI_INT8_T, i, 0, m_mpiInfo->comm);
//                     MPI_Send(&t.neighbour_x, 1, MPI_LONG, i, 0, m_mpiInfo->comm);
//                     MPI_Send(&t.neighbour_y, 1, MPI_LONG, i, 0, m_mpiInfo->comm);
//                     MPI_Send(&t.neighbour_z, 1, MPI_LONG, i, 0, m_mpiInfo->comm);
//                     MPI_Send(&t.neighbour_level, 1, MPI_INT, i, 0, m_mpiInfo->comm);
//                     MPI_Send(&t.neighbour_tree, 1, MPI_INT, i, 0, m_mpiInfo->comm);
//                 }

//                 num=2*hanging_edge_node_connections.size();
//                 MPI_Send(&num, 1, MPI_INT, i, 0, m_mpiInfo->comm);
//                 MPI_Send(hanging_edge_node_connections.data(), num, MPI_LONG, i, 0, m_mpiInfo->comm);

//                 num=2*false_node_connections.size();
//                 MPI_Send(&num, 1, MPI_INT, i, 0, m_mpiInfo->comm);
//                 MPI_Send(false_node_connections.data(), num, MPI_LONG, i, 0, m_mpiInfo->comm);

//                 num=3*NormalNodesTmp.size();
//                 MPI_Send(&num, 1, MPI_INT, i, 0, m_mpiInfo->comm);
//                 MPI_Send(NormalNodesTmp.data(), num, MPI_DOUBLE, i, 0, m_mpiInfo->comm);

//                 num=2*HangingFaceNodesTmp.size();
//                 MPI_Send(&num, 1, MPI_INT, i, 0, m_mpiInfo->comm);
//                 MPI_Send(HangingFaceNodesTmp.data(), num, MPI_LONG, i, 0, m_mpiInfo->comm);

//                 num=HangingEdgeNodesTmp.size();
//                 MPI_Send(&num, 1, MPI_INT, i, 0, m_mpiInfo->comm);
//                 MPI_Send(HangingEdgeNodesTmp.data(), num, MPI_LONG, i, 0, m_mpiInfo->comm);
//             }        
//         }
//         else
//         {
//             int num;
//             MPI_Recv(&num, 1, MPI_INT, 0, 0, m_mpiInfo->comm, MPI_STATUS_IGNORE);

//             octantIDs.resize(num);
//             MPI_Recv(octantIDs.data(), num, MPI_LONG, 0, 0, m_mpiInfo->comm, MPI_STATUS_IGNORE);
//             MPI_Recv(&num, 1, MPI_INT, 0, 0, m_mpiInfo->comm, MPI_STATUS_IGNORE);
//             hanging_face_orientation.clear();

//             for(int j = 0; j < num; j++ )
//             {
//                 hangingFaceInfo t;
//                 MPI_Recv(&t.x, 1, MPI_LONG, 0, 0, m_mpiInfo->comm, MPI_STATUS_IGNORE);
//                 MPI_Recv(&t.y, 1, MPI_LONG, 0, 0, m_mpiInfo->comm, MPI_STATUS_IGNORE);
//                 MPI_Recv(&t.z, 1, MPI_LONG, 0, 0, m_mpiInfo->comm, MPI_STATUS_IGNORE);
//                 MPI_Recv(&t.level, 1, MPI_INT8_T, 0, 0, m_mpiInfo->comm, MPI_STATUS_IGNORE);
//                 MPI_Recv(&t.treeid, 1, MPI_INT, 0, 0, m_mpiInfo->comm, MPI_STATUS_IGNORE);
//                 MPI_Recv(&t.face_type, 1, MPI_INT8_T, 0, 0, m_mpiInfo->comm, MPI_STATUS_IGNORE);
//                 MPI_Recv(&t.neighbour_x, 1, MPI_LONG, 0, 0, m_mpiInfo->comm, MPI_STATUS_IGNORE);
//                 MPI_Recv(&t.neighbour_y, 1, MPI_LONG, 0, 0, m_mpiInfo->comm, MPI_STATUS_IGNORE);
//                 MPI_Recv(&t.neighbour_z, 1, MPI_LONG, 0, 0, m_mpiInfo->comm, MPI_STATUS_IGNORE);
//                 MPI_Recv(&t.neighbour_level, 1, MPI_LONG, 0, 0, m_mpiInfo->comm, MPI_STATUS_IGNORE);
//                 MPI_Recv(&t.neighbour_tree, 1, MPI_INT, 0, 0, m_mpiInfo->comm, MPI_STATUS_IGNORE);
//                 hanging_face_orientation.push_back(t);
//             }

//             MPI_Recv(&num, 1, MPI_INT, 0, 0, m_mpiInfo->comm, MPI_STATUS_IGNORE);
//             std::vector<std::pair<long,long>> temp_face_node_connections;
//             temp_face_node_connections.resize(num);
//             MPI_Recv(temp_face_node_connections.data(), 2*num, MPI_LONG, 0, 0, m_mpiInfo->comm, MPI_STATUS_IGNORE);
//             hanging_face_node_connections=temp_face_node_connections;
//             MPI_Recv(&num, 1, MPI_INT, 0, 0, m_mpiInfo->comm, MPI_STATUS_IGNORE);
//             hanging_edge_orientation.clear();
//             for(int j = 0; j < num; j++ )
//             {
//                 hangingEdgeInfo t;
//                 MPI_Recv(&t.nodeid, 1, MPI_LONG, 0, 0, m_mpiInfo->comm, MPI_STATUS_IGNORE);
//                 MPI_Recv(&t.x, 1, MPI_LONG, 0, 0, m_mpiInfo->comm, MPI_STATUS_IGNORE);
//                 MPI_Recv(&t.y, 1, MPI_LONG, 0, 0, m_mpiInfo->comm, MPI_STATUS_IGNORE);
//                 MPI_Recv(&t.z, 1, MPI_LONG, 0, 0, m_mpiInfo->comm, MPI_STATUS_IGNORE);
//                 MPI_Recv(&t.level, 1, MPI_INT8_T, 0, 0, m_mpiInfo->comm, MPI_STATUS_IGNORE);
//                 MPI_Recv(&t.treeid, 1, MPI_INT, 0, 0, m_mpiInfo->comm, MPI_STATUS_IGNORE);
//                 MPI_Recv(&t.edge_type, 1, MPI_INT8_T, 0, 0, m_mpiInfo->comm, MPI_STATUS_IGNORE);
//                 MPI_Recv(&t.neighbour_x, 1, MPI_LONG, 0, 0, m_mpiInfo->comm, MPI_STATUS_IGNORE);
//                 MPI_Recv(&t.neighbour_y, 1, MPI_LONG, 0, 0, m_mpiInfo->comm, MPI_STATUS_IGNORE);
//                 MPI_Recv(&t.neighbour_z, 1, MPI_LONG, 0, 0, m_mpiInfo->comm, MPI_STATUS_IGNORE);
//                 MPI_Recv(&t.neighbour_level, 1, MPI_LONG, 0, 0, m_mpiInfo->comm, MPI_STATUS_IGNORE);
//                 MPI_Recv(&t.neighbour_tree, 1, MPI_INT, 0, 0, m_mpiInfo->comm, MPI_STATUS_IGNORE);
//                 hanging_edge_orientation.push_back(t);
//             }

//             MPI_Recv(&num, 1, MPI_INT, 0, 0, m_mpiInfo->comm, MPI_STATUS_IGNORE);
//             hanging_edge_node_connections.clear();
//             hanging_edge_node_connections.resize(0.5*num);
//             MPI_Recv(hanging_edge_node_connections.data(), num, MPI_LONG, 0, 0, m_mpiInfo->comm, MPI_STATUS_IGNORE);

//             MPI_Recv(&num, 1, MPI_INT, 0, 0, m_mpiInfo->comm, MPI_STATUS_IGNORE);
//             false_node_connections.clear();
//             false_node_connections.resize(0.5*num);
//             MPI_Recv(false_node_connections.data(), num, MPI_LONG, 0, 0, m_mpiInfo->comm, MPI_STATUS_IGNORE);

//             MPI_Recv(&num, 1, MPI_INT, 0, 0, m_mpiInfo->comm, MPI_STATUS_IGNORE);
//             NormalNodesTmp.clear();
//             NormalNodesTmp.resize(num/3);
//             MPI_Recv(NormalNodesTmp.data(), num, MPI_DOUBLE, 0, 0, m_mpiInfo->comm, MPI_STATUS_IGNORE);

//             MPI_Recv(&num, 1, MPI_INT, 0, 0, m_mpiInfo->comm, MPI_STATUS_IGNORE);
//             HangingFaceNodesTmp.clear();
//             HangingFaceNodesTmp.resize(0.5*num);
//             MPI_Recv(HangingFaceNodesTmp.data(), num, MPI_LONG, 0, 0, m_mpiInfo->comm, MPI_STATUS_IGNORE);

//             MPI_Recv(&num, 1, MPI_INT, 0, 0, m_mpiInfo->comm, MPI_STATUS_IGNORE);
//             HangingEdgeNodesTmp.clear();
//             HangingEdgeNodesTmp.resize(num);
//             MPI_Recv(HangingEdgeNodesTmp.data(), num, MPI_LONG, 0, 0, m_mpiInfo->comm, MPI_STATUS_IGNORE);
//         }
//     }
//     oxleytimer.toc("\t\t\t....done 3");
//     MPI_Barrier(m_mpiInfo->comm);
// #else
//     k=0;
//     for(p8est_topidx_t treeid = p8est->first_local_tree; treeid <= p8est->last_local_tree; ++treeid) {
//         p8est_tree_t * tree = p8est_tree_array_index(p8est->trees, treeid);
//         sc_array_t * tquadrants = &tree->quadrants;
//         p8est_locidx_t Q = (p8est_locidx_t) tquadrants->elem_count;
//         // Loop over octants
//         for(int q = 0; q < Q; ++q) { 
//             p8est_quadrant_t * oct = p8est_quadrant_array_index(tquadrants, q);
//             p8est_qcoord_t l = P8EST_QUADRANT_LEN(oct->level);

//             // Save octant information
//             double xyz[3];
//             p8est_qcoord_to_vertex(p8est->connectivity, treeid, oct->x, oct->y, oct->z, xyz);

//             octantIDs.push_back(NodeIDs.find(std::make_tuple(xyz[0],xyz[1],xyz[2]))->second);
//             oct_info tmp;
//             tmp.x=xyz[0];
//             tmp.y=xyz[1];
//             tmp.z=xyz[2];
//             tmp.level=oct->level;
//             octantInfo.push_back(tmp);

//             // // Get the hanging edge and face information for this octant
//             // int hanging_faces[6];
//             // int hanging_edges[12];
//             // getHangingInfo(nodes->face_code[k++], hanging_faces, hanging_edges);

//             // #ifdef OXLEY_ENABLE_DEBUG_RENUMBER_NODES_EXTRA
//             //     std::cout << "Octant  " << k-1 << ":" << std::endl;
//             //     int faceCount=6, edgeCount=12;
//             //     for(int i=0;i<6;i++) faceCount+=hanging_faces[i];
//             //     for(int i=0;i<12;i++) edgeCount+=hanging_edges[i];
//             //     std::cout << "\tHanging faces : ";
//             //     if(faceCount==0)
//             //         std::cout << "[none]";
//             //     else
//             //         for(int i=0; i<6;i++)
//             //             if(hanging_faces[i]!=-1)
//             //                 std::cout << i << "(" << hanging_faces[i] << "), " ;
//             //     std::cout << std::endl;
//             //     std::cout << "\tHanging edges : ";
//             //     if(edgeCount==0)
//             //         std::cout << "[none]";
//             //     else
//             //         for(int i=0; i<12;i++)
//             //             if(hanging_edges[i]!=-1)
//             //                 std::cout << i << "(" << hanging_edges[i] << "), " ;

//             //     std::cout << std::endl;
//             // #endif
        
//             // // Record hanging node information
//             // // Loop over faces
//             // for(int n = 0; n < 6; n++)
//             // {
//             //     // If the face does not have any hanging nodes then skip ahead
//             //     if(hanging_faces[n]==-1)
//             //         continue;

//             //     // *******************************************************************
//             //     // Corner node
//             //     // *******************************************************************
//             //     double xyzC[3];
//             //     p8est_qcoord_to_vertex(p8est->connectivity, treeid, oct->x, oct->y, oct->z, xyzC);
//             //     auto cornerPoint = std::make_tuple(xyzC[0],xyzC[1],xyzC[2]);
//             //     long nodeidC = NodeIDs.find(cornerPoint)->second;

//             //     // *******************************************************************
//             //     // Face node
//             //     // *******************************************************************
//             //     // shift to the corner of the face
//             //     signed long x_corner = ((int) xface[n][hanging_faces[n]]) * l;
//             //     signed long y_corner = ((int) yface[n][hanging_faces[n]]) * l;
//             //     signed long z_corner = ((int) zface[n][hanging_faces[n]]) * l;

//             //     double xyz[3];
//             //     p8est_qcoord_to_vertex(p8est->connectivity, treeid, 
//             //                             oct->x+x_corner, 
//             //                             oct->y+y_corner, 
//             //                             oct->z+z_corner, xyz);

//             //     // Add to list of nodes
//             //     DoubleTuple point = std::make_tuple(xyz[0],xyz[1],xyz[2]);
//             //     long nodeid = NodeIDs.find(point)->second;
//             //     std::vector<long>::iterator have_node = std::find(HangingFaceNodesTmp.begin(), HangingFaceNodesTmp.end(),nodeid);
//             //     if(have_node == HangingFaceNodesTmp.end())
//             //           HangingFaceNodesTmp.push_back(nodeid);

//             //     // *******************************************************************
//             //     // Get the parent octant
//             //     // *******************************************************************
//             //     p8est_quadrant_t parentOct;
//             //     p8est_quadrant_t * parent = &parentOct;
//             //     p8est_quadrant_parent(oct, parent);
                
//             //     // *******************************************************************
//             //     // Get the octant neighbouring this face
//             //     // *******************************************************************
//             //     // Get the neighbouring face
//             //     p8est_quadrant_t neighbourOct;
//             //     p8est_quadrant_t * neighbour = &neighbourOct;
//             //     // p8est_quadrant_face_neighbor(parent, n, neighbour);
//             //     int faceNeighbourTree=p8est_quadrant_face_neighbor_extra(parent,treeid,n,neighbour,NULL,connectivity);

//             //     // Save this information
//             //     hangingFaceInfo tmp;
//             //     tmp.x=oct->x;
//             //     tmp.y=oct->y;
//             //     tmp.z=oct->z;
//             //     tmp.level=oct->level;
//             //     tmp.treeid=treeid;
//             //     tmp.face_type=n;
//             //     tmp.neighbour_x=neighbour->x;
//             //     tmp.neighbour_y=neighbour->y;
//             //     tmp.neighbour_z=neighbour->z;
//             //     tmp.neighbour_level=neighbour->level;
//             //     tmp.neighbour_tree=faceNeighbourTree;
//             //     hanging_face_orientation.push_back(tmp);

//             //     // *******************************************************************
//             //     // Calculate the edge nodes
//             //     // *******************************************************************
//             //     signed long x_e_corner[2] = {((int) xedge[n][hanging_faces[n]][0]) * l,((int) xedge[n][hanging_faces[n]][1]) * l};
//             //     signed long y_e_corner[2] = {((int) yedge[n][hanging_faces[n]][0]) * l,((int) yedge[n][hanging_faces[n]][1]) * l};
//             //     signed long z_e_corner[2] = {((int) zedge[n][hanging_faces[n]][0]) * l,((int) zedge[n][hanging_faces[n]][1]) * l};
//             //     double xyzA[2][3];
//             //     p8est_qcoord_to_vertex(p8est->connectivity, treeid, oct->x+x_e_corner[0], oct->y+y_e_corner[0], oct->z+z_e_corner[0], xyzA[0]);
//             //     p8est_qcoord_to_vertex(p8est->connectivity, treeid, oct->x+x_e_corner[1], oct->y+y_e_corner[1], oct->z+z_e_corner[1], xyzA[1]);
//             //     DoubleTuple pointB[2] = {{std::make_tuple(xyzA[0][0],xyzA[0][1],xyzA[0][2])},{std::make_tuple(xyzA[1][0],xyzA[1][1],xyzA[1][2])}};
//             //     long nodeidB[2]= {NodeIDs.find(pointB[0])->second,NodeIDs.find(pointB[1])->second};

//             //     // Add to list of nodes
//             //     std::vector<long>::iterator have_edge;
//             //     have_edge = std::find(HangingEdgeNodesTmp.begin(), HangingEdgeNodesTmp.end(), nodeidB[0]);
//             //     if(have_edge == HangingEdgeNodesTmp.end())
//             //           HangingEdgeNodesTmp.push_back(nodeidB[0]);
//             //     have_edge = std::find(HangingEdgeNodesTmp.begin(), HangingEdgeNodesTmp.end(), nodeidB[1]);
//             //     if(have_edge == HangingEdgeNodesTmp.end())
//             //           HangingEdgeNodesTmp.push_back(nodeidB[1]);

//             //     // *******************************************************************
//             //     // Get the neighbouring edge(s)
//             //     // *******************************************************************
//             //     // Get the neighbours
//             //     p8est_quadrant_t EdgeNeighbourOct;
//             //     p8est_quadrant_t * EdgeNeighbour = &EdgeNeighbourOct;
//             //     int edge_number[2]={edge_lookup[n][hanging_faces[n]][0],edge_lookup[n][hanging_faces[n]][1]};
//             //     sc_array_t quads[2];
//             //     sc_array_t treeids[2];
//             //     sc_array_t nedges[2];
//             //     sc_array_t * pQuads[2] = {&quads[0],&quads[1]};
//             //     sc_array_t * pTreeids[2] = {&treeids[0],&treeids[1]};
//             //     sc_array_t * pNEdges[2] = {&nedges[0],&nedges[1]};
//             //     sc_array_init(pQuads[0], (size_t) sizeof(p4est_quadrant_t));
//             //     sc_array_init(pQuads[1], (size_t) sizeof(p4est_quadrant_t));
//             //     sc_array_init(pTreeids[0], (size_t) sizeof(p4est_topidx_t));
//             //     sc_array_init(pTreeids[1], (size_t) sizeof(p4est_topidx_t));
//             //     sc_array_init(pNEdges[0], (size_t) sizeof(int));
//             //     sc_array_init(pNEdges[1], (size_t) sizeof(int));

//             //     // Check for boundaries            
//             //     bool onBoundary[2];
//             //     onBoundary[0] = xyzA[0][0]<=0 || xyzA[0][0]>=forestData->m_length[0] ||
//             //                     xyzA[0][1]<=0 || xyzA[0][1]>=forestData->m_length[1] ||
//             //                     xyzA[0][2]<=0 || xyzA[0][2]>=forestData->m_length[2];
//             //     onBoundary[1] = xyzA[1][0]<=0 || xyzA[1][0]>=forestData->m_length[0] ||
//             //                     xyzA[1][1]<=0 || xyzA[1][1]>=forestData->m_length[1] ||
//             //                     xyzA[1][2]<=0 || xyzA[1][2]>=forestData->m_length[2];

//             //     hangingEdgeInfo temp;
//             //     if(!onBoundary[0])
//             //     {
//             //         p8est_quadrant_edge_neighbor_extra(parent,treeid,edge_number[0],pQuads[0],pTreeids[0],pNEdges[0],connectivity);
//             //         temp.nodeid=nodeidC;
//             //         temp.x=oct->x;
//             //         temp.y=oct->y;
//             //         temp.z=oct->z;
//             //         temp.level=oct->level;
//             //         temp.treeid=treeid;
//             //         temp.edge_type=n;
//             //         p8est_quadrant_t * tempB = (p8est_quadrant_t *) (*pQuads[0]).array;
//             //         temp.neighbour_x=tempB->x;
//             //         temp.neighbour_y=tempB->y;
//             //         temp.neighbour_z=tempB->z;
//             //         temp.neighbour_level=tempB->level;
//             //         p8est_topidx_t * tempC = (p8est_topidx_t *) (*pTreeids[0]).array;
//             //         temp.neighbour_tree=*tempC;
//             //         hanging_edge_orientation.push_back(temp);
//             //     }
//             //     if(!onBoundary[1])
//             //     {
//             //         p8est_quadrant_edge_neighbor_extra(parent,treeid,edge_number[1],pQuads[1],pTreeids[1],pNEdges[1],connectivity);
//             //         temp.nodeid=nodeidC;
//             //         temp.x=oct->x;
//             //         temp.y=oct->y;
//             //         temp.z=oct->z;
//             //         temp.level=oct->level;
//             //         temp.treeid=treeid;
//             //         temp.edge_type=n;
//             //         p8est_quadrant_t * tempB = (p8est_quadrant_t *) (*pQuads[1]).array;
//             //         temp.neighbour_x=tempB->x;
//             //         temp.neighbour_y=tempB->y;
//             //         temp.neighbour_z=tempB->z;
//             //         temp.neighbour_level=tempB->level;
//             //         p8est_topidx_t * tempC = (p8est_topidx_t *) (*pTreeids[1]).array;
//             //         temp.neighbour_tree=*tempC;
//             //         hanging_edge_orientation.push_back(temp);
//             //     }

//             //     // *******************************************************************
//             //     // Store connection information for updateRowsColumns
//             //     // *******************************************************************

//             //     signed long x_corner0 = xface0[n][hanging_faces[n]] * l;
//             //     signed long y_corner0 = yface0[n][hanging_faces[n]] * l;
//             //     signed long z_corner0 = zface0[n][hanging_faces[n]] * l;
            
//             //     double xyzD[3];
//             //     p8est_qcoord_to_vertex(p8est->connectivity, treeid, oct->x+x_corner0, oct->y+y_corner0, oct->z+z_corner0, xyzD);
//             //     auto otherCornerPoint = std::make_tuple(xyzD[0],xyzD[1],xyzD[2]);
//             //     long nodeA = NodeIDs.find(otherCornerPoint)->second;

//             //     p8est_qcoord_to_vertex(p8est->connectivity, treeid, oct->x+x_corner, oct->y+y_corner, oct->z+z_corner, xyzC);
//             //     auto facecornerPoint = std::make_tuple(xyzC[0],xyzC[1],xyzC[2]);
//             //     long facenodeidB = NodeIDs.find(facecornerPoint)->second;

//             //     long need_edgeA[4] = {std::count(hanging_edge_node_connections.begin(), hanging_edge_node_connections.end(), std::make_pair(nodeA,nodeidB[0])),
//             //                          std::count(hanging_edge_node_connections.begin(), hanging_edge_node_connections.end(), std::make_pair(nodeidB[0],nodeA)),
//             //                          std::count(hanging_edge_node_connections.begin(), hanging_edge_node_connections.end(), std::make_pair(nodeA,nodeidB[1])),
//             //                          std::count(hanging_edge_node_connections.begin(), hanging_edge_node_connections.end(), std::make_pair(nodeidB[1],nodeA))};
//             //     long need_edgeB[4] = {std::count(hanging_face_node_connections.begin(), hanging_face_node_connections.end(), std::make_pair(facenodeidB,nodeidB[0])),
//             //                          std::count(hanging_face_node_connections.begin(), hanging_face_node_connections.end(), std::make_pair(nodeidB[0],facenodeidB)),
//             //                          std::count(hanging_face_node_connections.begin(), hanging_face_node_connections.end(), std::make_pair(facenodeidB,nodeidB[1])),
//             //                          std::count(hanging_face_node_connections.begin(), hanging_face_node_connections.end(), std::make_pair(nodeidB[1],facenodeidB))};

                
//             //     // std::cout << " index[" << n << "][" << hanging_faces[n] << "]" << std::endl;

//             //     if(need_edgeA[0]+need_edgeA[1] == 0)
//             //     {
//             //         // Calculate fallacious node connection that was generated by p8est
//             //         signed long x_corner1 = xface1[n][hanging_faces[n]] * l;
//             //         signed long y_corner1 = yface1[n][hanging_faces[n]] * l;
//             //         signed long z_corner1 = zface1[n][hanging_faces[n]] * l;
//             //         p8est_qcoord_to_vertex(p8est->connectivity, treeid, oct->x+x_corner1, oct->y+y_corner1, oct->z+z_corner1, xyzD);
//             //         auto otherCornerPoint1 = std::make_tuple(xyzD[0],xyzD[1],xyzD[2]);
//             //         long defunct1 = NodeIDs.find(otherCornerPoint1)->second;

//             //         // std::cout << " defunct connection1 = " << nodeA << ", " << nodeidB[0] << ", " << defunct1 << std::endl;

//             //         hanging_edge_node_connections.push_back(std::make_pair(nodeA,nodeidB[0]));
//             //         false_node_connections.push_back(std::make_pair(nodeA,defunct1));

//             //         hanging_edge_node_connections.push_back(std::make_pair(nodeidB[0],defunct1));
//             //         false_node_connections.push_back(std::make_pair(defunct1,nodeA));
//             //     }
//             //     if(need_edgeA[2]+need_edgeA[3] == 0)
//             //     {
//             //         // Calculate fallacious node connection that was generated by p8est
//             //         signed long x_corner2 = xface2[n][hanging_faces[n]] * l;
//             //         signed long y_corner2 = yface2[n][hanging_faces[n]] * l;
//             //         signed long z_corner2 = zface2[n][hanging_faces[n]] * l;
//             //         p8est_qcoord_to_vertex(p8est->connectivity, treeid, oct->x+x_corner2, oct->y+y_corner2, oct->z+z_corner2, xyzD);
//             //         auto otherCornerPoint2 = std::make_tuple(xyzD[0],xyzD[1],xyzD[2]);
//             //         long defunct2 = NodeIDs.find(otherCornerPoint2)->second;

//             //         // std::cout << " defunct connection2 = " << nodeA << ", " << nodeidB[1] << ", " << defunct2 << std::endl;

//             //         hanging_edge_node_connections.push_back(std::make_pair(nodeA,nodeidB[1]));
//             //         false_node_connections.push_back(std::make_pair(nodeA,defunct2));

//             //         hanging_edge_node_connections.push_back(std::make_pair(nodeidB[1],defunct2));
//             //         false_node_connections.push_back(std::make_pair(defunct2,nodeA));
//             //     }
//             //     if(need_edgeB[0]+need_edgeB[1] == 0)
//             //     {
//             //         hanging_face_node_connections.push_back(std::make_pair(nodeidB[0],facenodeidB));
//             //     }
//             //     if(need_edgeB[2]+need_edgeB[3] == 0)
//             //     {
//             //         hanging_face_node_connections.push_back(std::make_pair(nodeidB[1],facenodeidB));
//             //     }


//             //     // *******************************************************************
//             //     // Debuging output
//             //     // *******************************************************************
//             //     #ifdef OXLEY_ENABLE_DEBUG_RENUMBER_NODES_EXTRA
//             //         std::cout << "\t\toctant = \033[1;36m" << k-1 << "\033[0m, corner node = " << nodeidC << std::endl;
//             //         std::cout << "\t\tface node = \033[1;36m" << facenodeidB << "\033[0m" << std::endl;
//             //         std::cout << "\t\t    index [" << n << "][" << hanging_faces[n] << "]" << std::endl;
//             //         // std::cout << "\t\t    shift was (Dx,Dy,Dz)= (" << x_corner << ", " << y_corner << ", " << z_corner << ")" << std::endl;
//             //         std::cout << "\t\tedge nodes: \033[1;36m" << nodeidB[0] << " & " << nodeidB[1] << "\033[0m" << std::endl;
//             //         std::cout << "\t\t    index[" << (int) n << "][" << (int) hanging_faces[n] << "][-]" << std::endl;
//             //         // std::cout << "\t\t    shift was (Dx,Dy,Dz)= (" << x_e_corner[0] << " ," << y_e_corner[0] << ", " << z_e_corner[0] << ")" << std::endl;
//             //         // std::cout << "\t\t    shift was (Dx,Dy,Dz)= (" << x_e_corner[1] << " ," << y_e_corner[1] << ", " << z_e_corner[1] << ")" << std::endl;
//             //         std::cout << "\t\tcorner nodes: \033[1;36m" << nodeA << "\033[0m" << std::endl;
//             //         std::cout << "\t\t    index [" << n << "][" << hanging_faces[n] << "][-]" << std::endl;
//             //         std::cout << "\t\tFace nodes: \033[1;36m" << facenodeidB << "--" << nodeidB[1] << "--" << nodeA << "--" << nodeidB[0] << "\033[0m" << std::endl;
//             //         std::cout << "\t\tedge numbers: \033[1;36m" << edge_number[0] << " & " << edge_number[1] << "\033[0m" << std::endl;
//             //         std::cout << "\t\t    index [" << n << "][" << hanging_faces[n] << "][-]" << std::endl;
//             //     #endif
//             }
//         }
//     }
//     oxleytimer.toc("\t\tdone");
// #endif // ESYS_MPI

// //     // Update num_hanging
// //     num_hanging=HangingFaceNodesTmp.size()+HangingEdgeNodesTmp.size();

// #ifdef OXLEY_ENABLE_DEBUG_RENUMBERNODES_INFORMATIONAL
//     std::cout << "Informational: " << std::endl;
//     std::cout << "\tNormal nodes = " << NormalNodesTmp.size() << std::endl;
//     std::cout << "\tHanging face nodes = " << HangingFaceNodesTmp.size() << std::endl;
//     std::cout << "\tHanging edge nodes = " << HangingEdgeNodesTmp.size() << std::endl;
//     std::cout << "\tTotal nodes = " << NormalNodesTmp.size()+HangingFaceNodesTmp.size()+HangingEdgeNodesTmp.size() << std::endl;
// #endif

//     // Populate NodeIDs
//     std::vector<int> is_hanging_tmp;
//     int num_nodes=NormalNodesTmp.size()+HangingFaceNodesTmp.size()+HangingEdgeNodesTmp.size();
//     is_hanging_tmp.resize(num_nodes,0);
//     NodeIDs.clear();

//     for(int i=0;i<NormalNodesTmp.size();i++)
//         NodeIDs[NormalNodesTmp[i]]=i;


//     // for(int i=0;i<HangingFaceNodesTmp.size();i++)
//     //     NodeIDs[HangingFaceNodesTmp[i]]=i;
//     // for(int i=0;i<HangingEdgeNodesTmp.size();i++)
//     //     NodeIDs[HangingEdgeNodesTmp[i]]=i;
//     for(int i=0;i<HangingFaceNodesTmp.size();i++)
//         is_hanging_tmp[HangingFaceNodesTmp[i]]=1;
//     for(int i=0;i<HangingEdgeNodesTmp.size();i++)
//         is_hanging_tmp[HangingEdgeNodesTmp[i]]=1;

//     // Renumber hanging nodes
//     // By custom, hanging numbers are numbered last
//     std::vector<int> new_node_ids(num_nodes,-1);
//     int count1=0;
//     int count2=num_nodes-num_hanging;
//     for(int i = 0; i < num_nodes; i++)
//     {
//         if(is_hanging_tmp[i]==0)
//             new_node_ids[i]=count1++;
//         else
//             new_node_ids[i]=count2++;
//     }
//     #ifdef OXLEY_ENABLE_DEBUG
//     ESYS_ASSERT(count1==num_nodes-num_hanging, "Unknown error (count1)");
//     ESYS_ASSERT(count2==num_nodes, "Unknown error (count2)");
//     #endif

//     // TODO vectorise

//     for(std::pair<DoubleTuple,long> e : NodeIDs)
//     {
//         NodeIDs[e.first]=new_node_ids[e.second];
//     }

// // #ifdef ESYS_MPI
// //     if(m_mpiInfo->size > 0)
// //     {
// //         oxleytimer.toc("communicating... ");
// //         if(m_mpiInfo->rank == 0)
// //         {
// //             for(int r = 1; r < m_mpiInfo->size; r++)
// //             {
// //                 int num = NodeIDs.size();
// //                 MPI_Send(&num, 1, MPI_INT, r, 0, m_mpiInfo->comm);

// //                 for(std::pair<DoubleTuple,long> e : NodeIDs)
// //                 {
// //                     int k = e.second;
// //                     double a=std::get<0>(e.first);
// //                     double b=std::get<1>(e.first);
// //                     double c=std::get<2>(e.first);

// //                     MPI_Send(&k, 1, MPI_LONG, r, 0, m_mpiInfo->comm);
// //                     MPI_Send(&a, 1, MPI_DOUBLE,r,0,m_mpiInfo->comm);
// //                     MPI_Send(&b, 1, MPI_DOUBLE,r,0,m_mpiInfo->comm);
// //                     MPI_Send(&c, 1, MPI_DOUBLE,r,0,m_mpiInfo->comm);                
// //                 }
// //             }
// //         }
// //         else
// //         {
// //             NodeIDs.clear();
// //             int num; 
// //             MPI_Recv(&num, 1, MPI_INT, 0, 0, m_mpiInfo->comm, MPI_STATUS_IGNORE);

// //             long k;
// //             double a, b, c;
// //             for(int i = 0; i < num; i++)
// //             {
// //                 MPI_Recv(&k, 1, MPI_LONG, 0, 0, m_mpiInfo->comm, MPI_STATUS_IGNORE);
// //                 MPI_Recv(&a, 1, MPI_DOUBLE, 0, 0, m_mpiInfo->comm, MPI_STATUS_IGNORE);
// //                 MPI_Recv(&b, 1, MPI_DOUBLE, 0, 0, m_mpiInfo->comm, MPI_STATUS_IGNORE);
// //                 MPI_Recv(&c, 1, MPI_DOUBLE, 0, 0, m_mpiInfo->comm, MPI_STATUS_IGNORE);

// //                 NodeIDs[std::make_tuple(a,b,c)]=k;
// //             }
// //         }
// //         oxleytimer.toc("\t\t\tdone 4");
// //     }
// // #endif

//     NormalNodes.resize(NormalNodesTmp.size(),std::make_tuple(-1.,-1.,-1));
//     for(int i = 0; i<NormalNodesTmp.size(); i++)
//         NormalNodes.push_back(NormalNodesTmp[new_node_ids[i]]);

//     HangingFaceNodes.clear();
//     HangingFaceNodes.resize(HangingFaceNodesTmp.size());
//     for(int i = 0; i<HangingFaceNodesTmp.size();i++)
//         HangingFaceNodes[i]=new_node_ids[HangingFaceNodesTmp[i]]; 

//     HangingEdgeNodes.clear();
//     HangingEdgeNodes.resize(HangingEdgeNodesTmp.size());
//     for(int i = 0; i<HangingEdgeNodesTmp.size();i++)
//         HangingEdgeNodes[i]=new_node_ids[HangingEdgeNodesTmp[i]];

//     // is_hanging.clear();
//     // is_hanging.resize(num_nodes,false);
//     // for(int i = 0; i<num_nodes; i++)
//     //     is_hanging[i]=is_hanging_tmp[new_node_ids[i]];

//     // // hanging_face_orientation

//     // for(int i = 0; i < hanging_edge_orientation.size(); i++)
//     //     hanging_edge_orientation[i].nodeid=new_node_ids[hanging_edge_orientation[i].nodeid];

//     // for(int i = 0; i < hanging_edge_node_connections.size(); i++)
//     // {
//     //     hanging_edge_node_connections[i].first=new_node_ids[hanging_edge_node_connections[i].first];
//     //     hanging_edge_node_connections[i].second=new_node_ids[hanging_edge_node_connections[i].second];
//     // }
    
//     // for(int i = 0; i < hanging_face_node_connections.size(); i++)
//     // {
//     //     hanging_face_node_connections[i].first=new_node_ids[hanging_face_node_connections[i].first];
//     //     hanging_face_node_connections[i].second=new_node_ids[hanging_face_node_connections[i].second];
//     // }
    
//     // for(int i = 0; i < false_node_connections.size(); i++)
//     // {
//     //     false_node_connections[i].first=new_node_ids[false_node_connections[i].first];
//     //     false_node_connections[i].second=new_node_ids[false_node_connections[i].second];
//     // }

//     // Populate m_nodeIDs
//     m_nodeId.clear();
//     m_nodeId.resize(NodeIDs.size());
//     long count=0;
//     for(std::pair<DoubleTuple,long> e : NodeIDs)
//     {
//         m_nodeId[count++]=e.second;
//     }

// #ifdef OXLEY_ENABLE_DEBUG_RENUMBER_NODES_PRINT_OCTANT_INFO
//     std::cout << "There are " << quadrantIDs.size() << " octants" << std::endl;
//     for(int i = 0; i < quadrantInfo.size(); i++)
//     {
//         std::cout << i << ": (" << quadrantInfo[i].x << ", " << quadrantInfo[i].y << "), l =" 
//                 << quadrantInfo[i].level << std::endl;
//     }
// #endif

// #ifdef OXLEY_ENABLE_DEBUG_RENUMBER_NODES_PRINT_NODEIDS
//     std::cout << "Printing NodeIDs " << std::endl;
//     double xyf[NodeIDs.size()][3]={{0}};

//     for(std::pair<DoubleTuple,long> e : NodeIDs)
//     {
//         xyf[e.second][0]=std::get<0>(e.first);
//         xyf[e.second][1]=std::get<1>(e.first);
//         xyf[e.second][2]=std::get<2>(e.first);
//     }

//     for(int i=0; i<NodeIDs.size(); i++)
//         std::cout << i << ": " << xyf[i][0] << ", " << xyf[i][1] << ", " << xyf[i][2] << std::endl;
// #endif

// #ifdef OXLEY_ENABLE_DEBUG_RENUMBER_NODES_PRINT_NODEIDS_HANGING
//     std::cout << "Hanging nodes on edges : ";
//     if(HangingEdgeNodes.size()==0)
//         std::cout << "[none]";
//     else
//         for(int i = 0; i < HangingEdgeNodes.size(); i++)
//             std::cout << HangingEdgeNodes[i] << ", ";
//     std::cout << std::endl;
//     std::cout << "Hanging nodes on faces : ";
//     if(HangingFaceNodes.size()==0)
//         std::cout << "[none]";
//     else
//         for(int i = 0; i < HangingFaceNodes.size(); i++)
//             std::cout << HangingFaceNodes[i] << ", ";
//     std::cout << std::endl;
// #endif

// #ifdef OXLEY_ENABLE_DEBUG_RENUMBER_NODES_PRINT_HANGING_NODE_CONNECTIONS
//     std::cout << "Hanging nodes connections : ";
//     if(hanging_node_connections.size()==0)
//         std::cout << "[none]";
//     else
//         for(int i = 0; i < hanging_node_connections.size(); i++)
//             std::cout << hanging_node_connections[i].first << "---" << hanging_node_connections[i].second << ", ";
//     std::cout << std::endl;
// #endif

// #ifdef ESYS_MPI
//     oxleytimer.toc("waiting at barrier");
//     MPI_Barrier(m_mpiInfo->comm);
// #endif

//     oxleytimer.toc("Nodes renumbered"); 
// }

bool Brick::gotAlready(DoubleTuple point, std::vector<DoubleTuple> vec)
{
    double a = std::get<0>(point);
    double b = std::get<1>(point);
    double c = std::get<2>(point);

    for(int i = 0; i < vec.size(); i++)
    {
        if(    (std::abs(std::get<0>(vec[i]) - a) < tuple_tolerance)
            && (std::abs(std::get<1>(vec[i]) - b) < tuple_tolerance) 
            && (std::abs(std::get<2>(vec[i]) - c) < tuple_tolerance))
                    return true;
    }
    return false;
}

bool Brick::gotPoint(DoubleTuple point, std::vector<DoubleTuple> vec)
{
    double a = std::get<0>(point);
    double b = std::get<1>(point);
    double c = std::get<2>(point);

    for(int i = 0; i < vec.size(); i++)
    {
        if(    (std::abs(std::get<0>(vec[i]) - a) < tuple_tolerance)
            && (std::abs(std::get<1>(vec[i]) - b) < tuple_tolerance) 
            && (std::abs(std::get<2>(vec[i]) - c) < tuple_tolerance))
        {
             return true;
        }
    }
    return false;
}


bool Brick::hasDuplicate(DoubleTuple point, std::vector<DoubleTuple> vec, bool serial)
{
    if(serial)
    {
        int count=0;

        double a = std::get<0>(point);
        double b = std::get<1>(point);
        double c = std::get<2>(point);

        for(int i = 0; i < vec.size(); i++)
        {
            if(    (std::abs(std::get<0>(vec[i]) - a) < tuple_tolerance)
                && (std::abs(std::get<1>(vec[i]) - b) < tuple_tolerance) 
                && (std::abs(std::get<2>(vec[i]) - c) < tuple_tolerance))
            {
                count++;
                if(count > 1)
                    return true;
            }
        }
        return false;
    }

// #ifdef ESYS_MPI

//     int rank;// = m_mpiInfo->rank;
//     MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//     int size; // = m_mpiInfo->size;
//     MPI_Comm_size(MPI_COMM_WORLD, &size);

//     //compiled with MPI but currently being run on a single process
//     if(size==1) 
//     {
//         int thread_num = omp_get_thread_num();
//         bool gotPoint[omp_get_num_threads()]={false};
//         // double tol = 0.0;
//     #pragma omp parallel for 
//         for(int i = 0 ; i < vec.size(); i++)
//             if( (std::abs(std::get<0>(vec[i]) - std::get<0>(point)) < 1e-12)
//                 && (std::abs(std::get<1>(vec[i]) - std::get<1>(point)) < 1e-12) 
//                     && (std::abs(std::get<2>(vec[i]) - std::get<2>(point)) < 1e-12))
//                         gotPoint[thread_num] = true;
//     #pragma omp barrier
//         int count=0;
//         for(int i = 0; i < omp_get_num_threads(); i++ )
//             if(gotPoint[i])
//                 count++;
//         if(count > 1)
//             return true;
//         else
//             return false;
//     }
//     else // more than one mpi process
//     {
//         bool PointInVector=false;
//         int foundPoint = 0;
//         // search our part of the vector for point
//         for(int i = 0 ; i < vec.size(); i++)
//         {
//             if(i % size != rank)
//                 continue;

//             bool test = (std::abs(std::get<0>(vec[i]) - std::get<0>(point)) < 1e-12)
//                                 && (std::abs(std::get<1>(vec[i]) - std::get<1>(point)) < 1e-12) 
//                                     && (std::abs(std::get<2>(vec[i]) - std::get<2>(point)) < 1e-12);

//             if(test==true)
//             {
//                 foundPoint = 1;
//                 break;
//             }
//         }
//         MPI_Barrier(m_mpiInfo->comm);
//         // relay results to 0
//         bool somoneElseFoundPoint=false;
//         if(rank != 0)
//         {
//             MPI_Send(&foundPoint, 1,  MPI_INT, 0, 0, m_mpiInfo->comm);
//         }
//         else
//         {
//             // check to see if someone else found it
//             for(int j = 1; j < size ; j++)
//             {
//                 MPI_Status *status;
//                 int temp=0;
//                 MPI_Recv(&temp, 1, MPI_INT, j, 0, m_mpiInfo->comm, status);
//                 if(temp == 1)
//                     somoneElseFoundPoint=true;
//             }
//         }
//         return somoneElseFoundPoint || foundPoint;
//     }
// #else // NOMPI
    #ifdef OPENMPFLAG
        int thread_num = omp_get_thread_num();
        bool gotPoint[omp_get_num_threads()]={false};
        #pragma omp parallel shared(vec, gotPoint)
        {
        #pragma omp parallel for 
            for(int i = 0 ; i < vec.size(); i++)
                if((std::abs(std::get<0>(vec[i]) - std::get<0>(point)) < tuple_tolerance)
                    && (std::abs(std::get<1>(vec[i]) - std::get<1>(point)) < tuple_tolerance) 
                        && (std::abs(std::get<2>(vec[i]) - std::get<2>(point)) < tuple_tolerance))
                            gotPoint[thread_num] = true;
        }
        #pragma omp barrier
        int count = 0;
        for(int i = 0; i < omp_get_num_threads(); i++ )
            if(gotPoint[i])
                count++;
        if(count > 1)
            return true;
        else
            return false;
    #else // No MPI or OPENMP
        int count=0;
        for(int i = 0 ; i < vec.size(); i++)
            if((std::abs(std::get<0>(vec[i]) - std::get<0>(point)) < tuple_tolerance)
                && (std::abs(std::get<1>(vec[i]) - std::get<1>(point)) < tuple_tolerance) 
                    && (std::abs(std::get<2>(vec[i]) - std::get<2>(point)) < tuple_tolerance))
            {
                count++;
                if(count > 1)
                    return true;
            }
        return false;
    #endif //OPENMPFLAG

// #endif
}

// A faster version of p8est_qcoord_to_vertex
void Brick::p8est_qcoord_to_vertex_fast(p8est_connectivity_t * connectivity,
                        p4est_topidx_t treeid,
                        p4est_qcoord_t x, p4est_qcoord_t y, p4est_qcoord_t z,
                        double vxyz[3])
{
    const double       *vertices = connectivity->vertices;
    const p4est_topidx_t *vindices;
    // int                 xi, yi, zi;
    double              wx[2], wy[2], wz[2];

    vindices = connectivity->tree_to_vertex + 8 * treeid;

    vxyz[0] = vxyz[1] = vxyz[2] = 0.;

    double divisor = ((p4est_qcoord_t) 1 << P8EST_MAXLEVEL);

    wx[1] = (double) x / divisor;
    wx[0] = 1. - wx[1];

    wy[1] = (double) y / divisor;
    wy[0] = 1. - wy[1];

    wz[1] = (double) z / divisor;
    wz[0] = 1. - wz[1];

    // double w[6] = {x/divisor, (1-x)/divisor, y/divisor, (1-y)/divisor, z/divisor, (1-z)/divisor};

    int ii[8][3] ={{0,0,0},{0,0,1},{0,1,0},{0,1,1},{1,0,0},{1,0,1},{1,1,0},{1,1,1}};

    for(int i = 0; i < 8; i++) {
        double xfactor = wz[ii[i][0]] * wy[ii[i][1]] * wx[ii[i][2]];
        p4est_topidx_t vindex = 3 * (*vindices++);
        vxyz[0] += xfactor * vertices[vindex++];
        vxyz[1] += xfactor * vertices[vindex++];
        vxyz[2] += xfactor * vertices[vindex++];
    }
}

//protected
void Brick::assembleCoordinates(escript::Data& arg) const
{    
    if (!arg.isDataPointShapeEqual(1, &m_numDim))
        throw ValueError("assembleCoordinates: Invalid Data object shape");
    if (!arg.numSamplesEqual(1, getNumNodes()))
        throw ValueError("assembleCoordinates: Illegal number of samples in Data object");
    arg.requireWrite();

#ifdef OXLEY_ENABLE_DEBUG_ASSEMBLE_COORDINATES
    std::cout << "assemble coordinates " << std::endl;
    float bounds[6]={0.0};
#endif

    std::vector<bool> duplicates(getNumNodes(),false);

    const int V = nodes->vnodes;   // 8 corners (degree-1 lnodes)
    long e = 0;                    // running local leaf index (lnodes order)
    for(p8est_topidx_t treeid = p8est->first_local_tree; treeid <= p8est->last_local_tree; ++treeid)
    {
        p8est_tree_t * tree = p8est_tree_array_index(p8est->trees, treeid);
        sc_array_t * tquadrants = &tree->quadrants;
        p8est_locidx_t Q = (p8est_locidx_t) tquadrants->elem_count;
        for(int q = 0; q < Q; ++q, ++e)
        {
            p8est_quadrant_t * quad = p8est_quadrant_array_index(tquadrants, q);
            p8est_qcoord_t l = P8EST_QUADRANT_LEN(quad->level);
            int dxyz[8][3] = {{0,0,0},{l,0,0},{0,l,0},{l,l,0},{0,0,l},{l,0,l},{0,l,l},{l,l,l}};

            // A HANGING corner has no node of its own: the slot holds a MASTER,
            // which is a node of the coarse neighbour and lies OUTSIDE this
            // octant. Writing this corner's position into it would move the
            // master to the hanging position, so skip those slots - the master
            // gets its coordinate from an octant where it is a real corner.
            int masterCount[P8EST_CHILDREN], masterSlot[P8EST_CHILDREN][4];
            const bool anyHanging = getHangingNodes(nodes->face_code[e],
                                                    masterCount, masterSlot);

            for(int n = 0; n < 8; ++n)
            {
                if(anyHanging && masterCount[n] > 0)
                    continue;

                double xy[3];
                p8est_qcoord_to_vertex(p8est->connectivity, treeid, quad->x+dxyz[n][0], quad->y+dxyz[n][1], quad->z+dxyz[n][2], xy);
#ifdef OXLEY_ENABLE_DEBUG_ASSEMBLE_COORDINATES
                std::cout << "(" << xy[0] << ", " << xy[1] << ", " << xy[2] << ")" << std::endl;
#endif
                // lnodes local node id for this corner (no coordinate hashing)
                long lni = (long) nodes->element_nodes[(size_t) e * V + n];
                if(duplicates[lni] == true)
                    continue;
                else
                    duplicates[lni] = true;
                double * point = arg.getSampleDataRW(lni);
                point[0] = xy[0];
                point[1] = xy[1];
                point[2] = xy[2];
#ifdef OXLEY_ENABLE_DEBUG_ASSEMBLE_COORDINATES
                std::cout<<"lni="<<lni<<"\tCorner=("<<xy[0]<<","<<xy[1]<<"),"<<std::endl;
                if(xy[0] < bounds[0]) bounds[0]=xy[0];
                if(xy[0] > bounds[1]) bounds[1]=xy[0];
                if(xy[1] < bounds[2]) bounds[2]=xy[1];
                if(xy[1] > bounds[3]) bounds[3]=xy[1];
                if(xy[2] < bounds[4]) bounds[4]=xy[2];
                if(xy[2] > bounds[5]) bounds[5]=xy[2];
#endif
                // }
            }
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
        << " and " << bounds[2] << ", " << bounds[3] << 
        " and " << bounds[4] << ", " << bounds[5] << std::endl;
#endif
}


//private
template<typename Scalar>
void Brick::addToMatrixAndRHS(escript::AbstractSystemMatrix* S, escript::Data& F,
         const std::vector<Scalar>& EM_S, const std::vector<Scalar>& EM_F, 
         bool addS, bool addF, borderNodeInfo quad, int nEq, int nComp) const
{    
    // A boundary face element has 4 nodes (stored in neighbours[0..3]); EM_S is
    // 4x4 and EM_F length 4 (per equation component).
    IndexVector rowIndex(4);
    for(int i = 0; i < 4; i++)
        rowIndex[i] = (index_t) quad.neighbours[i];
    if(addF)
    {
        Scalar* F_p = F.getSampleDataRW(0, static_cast<Scalar>(0));
        for(int i=0; i < 4; i++) {
            if (rowIndex[i]<getNumDOF()) {
                for(int eq=0; eq<nEq; eq++) {
                    F_p[INDEX2(eq, rowIndex[i], nEq)]+=EM_F[INDEX2(eq,i,nEq)];
                }
            }
        }
    }
    if(addS)
        addToSystemMatrix<Scalar>(S, rowIndex, nEq, EM_S);
}


template
void Brick::addToMatrixAndRHS<real_t>(escript::AbstractSystemMatrix* S, escript::Data& F,
         const std::vector<real_t>& EM_S, const std::vector<real_t>& EM_F, 
         bool addS, bool addF, borderNodeInfo firstNode, int nEq, int nComp) const;

template
void Brick::addToMatrixAndRHS<cplx_t>(escript::AbstractSystemMatrix* S, escript::Data& F,
         const std::vector<cplx_t>& EM_S, const std::vector<cplx_t>& EM_F, 
         bool addS, bool addF, borderNodeInfo firstNode, int nEq, int nComp) const;

template
void Brick::addToMatrixAndRHS<real_t>(escript::AbstractSystemMatrix* S, escript::Data& F,
         const std::vector<real_t>& EM_S, const std::vector<real_t>& EM_F, 
         bool addS, bool addF, index_t e, index_t t, int nEq, int nComp) const;

template
void Brick::addToMatrixAndRHS<cplx_t>(escript::AbstractSystemMatrix* S, escript::Data& F,
         const std::vector<cplx_t>& EM_S, const std::vector<cplx_t>& EM_F,
         bool addS, bool addF, index_t e, index_t t, int nEq, int nComp) const;

//protected
template<typename Scalar>
void Brick::addToMatrixAndRHSGhost(escript::AbstractSystemMatrix* S, escript::Data& F,
         const std::vector<Scalar>& EM_S, const std::vector<Scalar>& EM_F,
         bool addS, bool addF, const index_t* rowIndex, int nEq, int nComp) const
{
    // rowIndex are extended-local corner node ids of a ghost (halo) octant.
    // Only OWNED rows (< getNumDOF()) are kept; columns may be 2nd-layer ghosts.
    if(addF)
    {
        Scalar* F_p = F.getSampleDataRW(0, static_cast<Scalar>(0));
        for(int i=0; i<8; i++) {
            if (rowIndex[i]<getNumDOF()) {
                for(int eq=0; eq<nEq; eq++) {
                    F_p[INDEX2(eq, rowIndex[i], nEq)]+=EM_F[INDEX2(eq,i,nEq)];
                }
            }
        }
    }
    if(addS)
    {
        IndexVector rowInd(rowIndex, rowIndex+8);
        addToSystemMatrix<Scalar>(S, rowInd, nEq, EM_S);
    }
}

template
void Brick::addToMatrixAndRHSGhost<real_t>(escript::AbstractSystemMatrix* S, escript::Data& F,
         const std::vector<real_t>& EM_S, const std::vector<real_t>& EM_F,
         bool addS, bool addF, const index_t* rowIndex, int nEq, int nComp) const;
template
void Brick::addToMatrixAndRHSGhost<cplx_t>(escript::AbstractSystemMatrix* S, escript::Data& F,
         const std::vector<cplx_t>& EM_S, const std::vector<cplx_t>& EM_F,
         bool addS, bool addF, const index_t* rowIndex, int nEq, int nComp) const;

//protected
template<typename Scalar>
std::vector<Scalar> Brick::exchangeGhostCoeff(const escript::Data& coef) const
{
    std::vector<Scalar> out;
    if (!ghost || coef.isEmpty())
        return out;
    const long nGhost  = (long) ghost->ghosts.elem_count;
    const long nMirror = (long) ghost->mirrors.elem_count;
    // In-memory bytes per getSampleDataRO() sample: an expanded Data stores one
    // value per quadrature point, a constant/tagged Data only a single point.
    const size_t sampleSize = (coef.actsExpanded()
                                ? (size_t) coef.getNumDataPointsPerSample() : 1)
                            * (size_t) coef.getDataPointSize();
    if (nGhost == 0 || sampleSize == 0)
        return out;
    const Scalar zero = static_cast<Scalar>(0);
    std::vector<const void*> mirror_data((size_t) nMirror, nullptr);
    for (long m = 0; m < nMirror; ++m) {
        p8est_quadrant_t* mq = p8est_quadrant_array_index(&ghost->mirrors, m);
        const long le = (long) mq->p.piggy3.local_num;
        mirror_data[m] = (const void*) coef.getSampleDataRO(le, zero);
    }
    out.resize((size_t) nGhost * sampleSize);
    p8est_ghost_exchange_custom(p8est, ghost, sampleSize * sizeof(Scalar),
                                const_cast<void**>(mirror_data.data()), out.data());
    return out;
}

template std::vector<real_t> Brick::exchangeGhostCoeff<real_t>(const escript::Data&) const;
template std::vector<cplx_t> Brick::exchangeGhostCoeff<cplx_t>(const escript::Data&) const;

//protected
template<typename Scalar>
void Brick::addToMatrixAndRHSGhostFace(escript::AbstractSystemMatrix* S, escript::Data& F,
         const std::vector<Scalar>& EM_S, const std::vector<Scalar>& EM_F,
         bool addS, bool addF, const index_t* rowIndex, int nEq, int nComp) const
{
    // rowIndex are the 4 extended-local face-node ids of a ghost boundary face.
    // Only OWNED rows (< getNumDOF()) are kept.
    if(addF)
    {
        Scalar* F_p = F.getSampleDataRW(0, static_cast<Scalar>(0));
        for(int i=0; i<4; i++) {
            if (rowIndex[i]<getNumDOF()) {
                for(int eq=0; eq<nEq; eq++)
                    F_p[INDEX2(eq, rowIndex[i], nEq)]+=EM_F[INDEX2(eq,i,nEq)];
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
void Brick::addToMatrixAndRHSGhostFace<real_t>(escript::AbstractSystemMatrix* S, escript::Data& F,
         const std::vector<real_t>& EM_S, const std::vector<real_t>& EM_F,
         bool addS, bool addF, const index_t* rowIndex, int nEq, int nComp) const;
template
void Brick::addToMatrixAndRHSGhostFace<cplx_t>(escript::AbstractSystemMatrix* S, escript::Data& F,
         const std::vector<cplx_t>& EM_S, const std::vector<cplx_t>& EM_F,
         bool addS, bool addF, const index_t* rowIndex, int nEq, int nComp) const;

//protected
template<typename Scalar>
std::vector<Scalar> Brick::exchangeGhostBoundary(const escript::Data& d,
                        const escript::Data& y, size_t& dSize, size_t& ySize) const
{
    std::vector<Scalar> out;
    dSize = d.isEmpty()?0:(size_t)(d.actsExpanded()?d.getNumDataPointsPerSample():1)*d.getDataPointSize();
    ySize = y.isEmpty()?0:(size_t)(y.actsExpanded()?y.getNumDataPointsPerSample():1)*y.getDataPointSize();
    if (!ghost || (dSize==0 && ySize==0))
        return out;
    const long nGhost  = (long) ghost->ghosts.elem_count;
    const long nMirror = (long) ghost->mirrors.elem_count;
    if (nGhost == 0)
        return out;

    const Scalar zero = static_cast<Scalar>(0);
    const size_t perSide = 1 + dSize + ySize;
    const size_t perOct  = 6 * perSide;

    // octant coords -> local leaf index (mirror piggy3.local_num is reliable)
    struct Key { p4est_topidx_t t; p4est_qcoord_t x, y, z;
                 bool operator==(const Key& o) const { return t==o.t && x==o.x && y==o.y && z==o.z; } };
    struct KeyHash { size_t operator()(const Key& k) const {
        return ((size_t)k.t*73856093u) ^ ((size_t)k.x*19349663u)
             ^ ((size_t)k.y*83492791u) ^ ((size_t)k.z*49979687u); } };
    std::unordered_map<Key, long, KeyHash> octKey2leaf;
    long leaf = 0;
    for (p8est_topidx_t tt = p8est->first_local_tree; tt <= p8est->last_local_tree; ++tt) {
        p8est_tree_t* tree = p8est_tree_array_index(p8est->trees, tt);
        sc_array_t* quads = &tree->quadrants;
        const long Q = (long) quads->elem_count;
        for (long q = 0; q < Q; ++q, ++leaf) {
            p8est_quadrant_t* qd = p8est_quadrant_array_index(quads, q);
            octKey2leaf[Key{ tt, qd->x, qd->y, qd->z }] = leaf;
        }
    }

    std::unordered_map<long, std::array<long,6>> fmap;   // leaf -> per-side sample
    const std::vector<borderNodeInfo>* lists[6] =
        { &NodeIDsLeft, &NodeIDsRight, &NodeIDsBottom, &NodeIDsTop, &NodeIDsAbove, &NodeIDsBelow };
    for (int s = 0; s < 6; ++s) {
        if (m_faceOffset[s] < 0) continue;
        const std::vector<borderNodeInfo>& Lst = *lists[s];
        for (long k = 0; k < (long) Lst.size(); ++k) {
            auto lit = octKey2leaf.find(Key{ Lst[k].treeid, Lst[k].x, Lst[k].y, Lst[k].z });
            if (lit == octKey2leaf.end()) continue;
            const long lf = lit->second;
            auto it = fmap.find(lf);
            if (it == fmap.end())
                it = fmap.emplace(lf, std::array<long,6>{{-1,-1,-1,-1,-1,-1}}).first;
            it->second[s] = (long) m_faceOffset[s] + k;
        }
    }

    std::vector<Scalar> mirrorPacked((size_t) nMirror * perOct, zero);
    std::vector<void*> mirror_data((size_t) nMirror, nullptr);
    for (long m = 0; m < nMirror; ++m) {
        p8est_quadrant_t* mq = p8est_quadrant_array_index(&ghost->mirrors, m);
        Scalar* base = &mirrorPacked[(size_t) m * perOct];
        mirror_data[m] = (void*) base;
        auto it = fmap.find((long) mq->p.piggy3.local_num);
        if (it == fmap.end()) continue;
        for (int s = 0; s < 6; ++s) {
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
    p8est_ghost_exchange_custom(p8est, ghost, perOct * sizeof(Scalar),
                                mirror_data.data(), out.data());
    return out;
}

template std::vector<real_t> Brick::exchangeGhostBoundary<real_t>(
        const escript::Data&, const escript::Data&, size_t&, size_t&) const;
template std::vector<cplx_t> Brick::exchangeGhostBoundary<cplx_t>(
        const escript::Data&, const escript::Data&, size_t&, size_t&) const;

//protected
void Brick::interpolateNodesOnElements(escript::Data& out,
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
void Brick::interpolateNodesOnFaces(escript::Data& out,
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

/// A face of the octant -> its four corners in z-order indexing, ordered
/// (a0b0, a1b0, a0b1, a1b1) in the face's two in-plane axes (a,b). This is the
/// order borderNodeInfo::neighbours is filled in and the face interpolation
/// and gradient read back. It is NOT the counter-clockwise winding
/// getMeshAccess needs for finley's outward normals - that is a second table.
const int faceCornerSlots[6][4] = {
    {0,2,4,6}, {1,3,5,7}, {0,1,4,5}, {2,3,6,7}, {0,1,2,3}, {4,5,6,7}
};

/**
   \brief
   Replaces the values sitting in an element's hanging corner slots by the
   values AT those corners.

   Same idea as Rectangle's version, but 3D has two kinds of hanging corner and
   that changes the arithmetic in one important way.

   lnodes does not number a hanging position, so its slot holds a MASTER: a node
   of the coarse neighbour lying outside this octant. A corner in the middle of
   a coarse FACE is the average of that face's four corners; one in the middle of
   a coarse EDGE is the average of its two ends. getHangingNodes returns those
   masters as slots of this same element, so the constraint stays element-local -
   no halo, no communication, no node that does not already exist.

   THE 3D DIFFERENCE: the chains overlap. The masters of a face centre are the
   four corners of that coarse face, and two of them are themselves hanging (the
   midpoints of the coarse face's edges). Correcting in place, as Rectangle does,
   would feed a corrected value back in and give
   (A + (A+B)/2 + (A+C)/2 + D)/4 instead of (A+B+C+D)/4. So every value is
   computed from a snapshot of the ORIGINAL slots. In 2D the only master is the
   anchor corner, which never hangs, which is why that version can work in place.

   Without this an element on the fine side of a seam reads values from 2h away
   as though they sat on its own corners. Nothing crashes and nothing looks
   wrong - grad() is simply wrong by O(1) along every seam.
*/
template<typename Scalar>
inline void constrainHangingCorners(const int masterCount[P8EST_CHILDREN],
                                    const int masters[P8EST_CHILDREN][4],
                                    Scalar* corner[P8EST_CHILDREN],
                                    dim_t numComp,
                                    std::vector<Scalar>& scratch)
{
    scratch.resize((size_t) P8EST_CHILDREN * numComp);
    for (int n = 0; n < P8EST_CHILDREN; ++n)
        for (dim_t i = 0; i < numComp; ++i)
            scratch[(size_t) n * numComp + i] = corner[n][i];

    for (int n = 0; n < P8EST_CHILDREN; ++n) {
        const int k = masterCount[n];
        if (k == 0)
            continue;
        for (dim_t i = 0; i < numComp; ++i) {
            Scalar s = static_cast<Scalar>(0);
            for (int j = 0; j < k; ++j)
                s += scratch[(size_t) masters[n][j] * numComp + i];
            corner[n][i] = s / static_cast<Scalar>(k);
        }
    }
}

} // anonymous namespace

//private
template <typename S>
void Brick::gatherCornersConstrained(const escript::Data& in, long quadIndex,
                                     dim_t numComp, S sentinel,
                                     std::vector<S>& corners,
                                     std::vector<S>& scratch) const
{
    const int V = nodes->vnodes;                       // 8
    corners.resize((size_t) P8EST_CHILDREN * numComp);
    for (int n = 0; n < P8EST_CHILDREN; ++n) {
        const index_t id = (index_t) nodes->element_nodes[(size_t) quadIndex * V + n];
        memcpy(&corners[(size_t) n * numComp], in.getSampleDataRO(id, sentinel),
               numComp * sizeof(S));
    }

    int masterCount[P8EST_CHILDREN], masterSlot[P8EST_CHILDREN][4];
    if (getHangingNodes(nodes->face_code[quadIndex], masterCount, masterSlot)) {
        S* corner[P8EST_CHILDREN];
        for (int n = 0; n < P8EST_CHILDREN; ++n)
            corner[n] = &corners[(size_t) n * numComp];
        constrainHangingCorners(masterCount, masterSlot, corner, numComp, scratch);
    }
}

//private
template <typename S>
void Brick::gatherFaceCornersConstrained(const escript::Data& in,
                                         const borderNodeInfo& b, int face,
                                         dim_t numComp, S sentinel,
                                         std::vector<S>& onFace,
                                         std::vector<S>& corners,
                                         std::vector<S>& scratch) const
{
    onFace.resize((size_t) 4 * numComp);

    // A domain-boundary face is never itself a seam - it has no neighbour at
    // all - but unlike 2D its CORNERS can still hang: a coarse octant beside it
    // can have an edge lying in the boundary plane, and the midpoint of that
    // edge is a corner of this face. So the values have to be constrained here
    // too, which needs the element they belong to. A domain built before
    // quadIndex was recorded falls back to reading the slots raw.
    if (b.quadIndex < 0) {
        for (int c = 0; c < 4; ++c)
            memcpy(&onFace[(size_t) c * numComp],
                   in.getSampleDataRO(b.neighbours[c], sentinel),
                   numComp * sizeof(S));
        return;
    }

    gatherCornersConstrained(in, b.quadIndex, numComp, sentinel, corners, scratch);
    for (int c = 0; c < 4; ++c)
        memcpy(&onFace[(size_t) c * numComp],
               &corners[(size_t) faceCornerSlots[face][c] * numComp],
               numComp * sizeof(S));
}

// private
template <typename S>
void Brick::interpolateNodesOnElementsWorker(escript::Data& out,
                                           const escript::Data& in,
                                           bool reduced, S sentinel) const
{
    const dim_t numComp = in.getDataPointSize();
    
    if (reduced) {
        out.requireWrite();

        std::vector<S> f_000(numComp);
        std::vector<S> f_001(numComp);
        std::vector<S> f_010(numComp);
        std::vector<S> f_011(numComp);
        std::vector<S> f_100(numComp);
        std::vector<S> f_101(numComp);
        std::vector<S> f_110(numComp);
        std::vector<S> f_111(numComp);

        std::vector<S> corners, scratch;
        long e = 0;
        for(p8est_topidx_t treeid = p8est->first_local_tree; treeid <= p8est->last_local_tree; ++treeid)
        {
            p8est_tree_t * tree = p8est_tree_array_index(p8est->trees, treeid);
            sc_array_t * tquadrants = &tree->quadrants;
            p8est_locidx_t Q = (p8est_locidx_t) tquadrants->elem_count;
            for(int q = 0; q < Q; q++, ++e)
            {
                const long quadID = e;
                // on the fine side of a seam the slots hold masters, not the
                // corner values; see constrainHangingCorners
                gatherCornersConstrained(in, e, numComp, sentinel, corners, scratch);
                memcpy(&f_000[0], &corners[0*numComp], numComp*sizeof(S));
                memcpy(&f_001[0], &corners[1*numComp], numComp*sizeof(S));
                memcpy(&f_010[0], &corners[2*numComp], numComp*sizeof(S));
                memcpy(&f_011[0], &corners[3*numComp], numComp*sizeof(S));
                memcpy(&f_100[0], &corners[4*numComp], numComp*sizeof(S));
                memcpy(&f_101[0], &corners[5*numComp], numComp*sizeof(S));
                memcpy(&f_110[0], &corners[6*numComp], numComp*sizeof(S));
                memcpy(&f_111[0], &corners[7*numComp], numComp*sizeof(S));
                S* o = out.getSampleDataRW(quadID, sentinel);
                for (index_t i=0; i < numComp; ++i) {
                    o[INDEX2(i,numComp,0)] = (f_000[i] + f_001[i] + f_010[i] + f_011[i] + f_100[i] + f_101[i] + f_110[i] + f_111[i])/static_cast<S>(8);
                } // end of component loop i
            }
        }
    } 
    else 
    {
        out.requireWrite();
        const S c0 = .0094373878376559314545;
        const S c1 = .035220810900864519624;
        const S c2 = .13144585576580214704;
        const S c3 = .49056261216234406855;

        std::vector<S> f_000(numComp);
        std::vector<S> f_001(numComp);
        std::vector<S> f_010(numComp);
        std::vector<S> f_011(numComp);
        std::vector<S> f_100(numComp);
        std::vector<S> f_101(numComp);
        std::vector<S> f_110(numComp);
        std::vector<S> f_111(numComp);

        std::vector<S> corners, scratch;
        long e = 0;
        for(p8est_topidx_t treeid = p8est->first_local_tree; treeid <= p8est->last_local_tree; ++treeid) {
            p8est_tree_t * tree = p8est_tree_array_index(p8est->trees, treeid);
            sc_array_t * tquadrants = &tree->quadrants;
            p8est_locidx_t Q = (p8est_locidx_t) tquadrants->elem_count;
            for(int q = 0; q < Q; q++, ++e)
            {
                const long quadId = e;

                #ifdef OXLEY_ENABLE_DEBUG_INTERPOLATE_QUADIDS
                    std::cout << "interpolateNodesOnElementsWorker quadID: " << quadId << std::endl;
                #endif

                // on the fine side of a seam the slots hold masters, not the
                // corner values; see constrainHangingCorners
                gatherCornersConstrained(in, e, numComp, sentinel, corners, scratch);
                memcpy(&f_000[0], &corners[0*numComp], numComp*sizeof(S));
                memcpy(&f_001[0], &corners[1*numComp], numComp*sizeof(S));
                memcpy(&f_010[0], &corners[2*numComp], numComp*sizeof(S));
                memcpy(&f_011[0], &corners[3*numComp], numComp*sizeof(S));
                memcpy(&f_100[0], &corners[4*numComp], numComp*sizeof(S));
                memcpy(&f_101[0], &corners[5*numComp], numComp*sizeof(S));
                memcpy(&f_110[0], &corners[6*numComp], numComp*sizeof(S));
                memcpy(&f_111[0], &corners[7*numComp], numComp*sizeof(S));
                S* o = out.getSampleDataRW(quadId, sentinel);
            #pragma omp parallel for
                for (index_t i=0; i < numComp; ++i) {
                    o[INDEX2(i,numComp,0)] = f_000[i]*c3 + f_111[i]*c0 + c2*(f_001[i] + f_010[i] + f_100[i]) + c1*(f_011[i] + f_101[i] + f_110[i]);
                    o[INDEX2(i,numComp,1)] = f_011[i]*c0 + f_100[i]*c3 + c2*(f_000[i] + f_101[i] + f_110[i]) + c1*(f_001[i] + f_010[i] + f_111[i]);
                    o[INDEX2(i,numComp,2)] = f_010[i]*c3 + f_101[i]*c0 + c2*(f_000[i] + f_011[i] + f_110[i]) + c1*(f_001[i] + f_100[i] + f_111[i]);
                    o[INDEX2(i,numComp,3)] = f_001[i]*c0 + f_110[i]*c3 + c2*(f_010[i] + f_100[i] + f_111[i]) + c1*(f_000[i] + f_011[i] + f_101[i]);
                    o[INDEX2(i,numComp,4)] = f_001[i]*c3 + f_110[i]*c0 + c2*(f_000[i] + f_011[i] + f_101[i]) + c1*(f_010[i] + f_100[i] + f_111[i]);
                    o[INDEX2(i,numComp,5)] = f_010[i]*c0 + f_101[i]*c3 + c2*(f_001[i] + f_100[i] + f_111[i]) + c1*(f_000[i] + f_011[i] + f_110[i]);
                    o[INDEX2(i,numComp,6)] = f_011[i]*c3 + f_100[i]*c0 + c2*(f_001[i] + f_010[i] + f_111[i]) + c1*(f_000[i] + f_101[i] + f_110[i]);
                    o[INDEX2(i,numComp,7)] = f_000[i]*c0 + f_111[i]*c3 + c2*(f_011[i] + f_101[i] + f_110[i]) + c1*(f_001[i] + f_010[i] + f_100[i]);
                } 
            }
        } 
        
    }
}

//protected


//private
template <typename S>
void Brick::interpolateNodesOnFacesWorker(escript::Data& out, const escript::Data& in,
                                    bool reduced, S sentinel) const
{
    const dim_t numComp = in.getDataPointSize();
    if (reduced) {
        out.requireWrite();

        std::vector<S> f_000(numComp);
        std::vector<S> f_001(numComp);
        std::vector<S> f_010(numComp);
        std::vector<S> f_011(numComp);
        std::vector<S> f_100(numComp);
        std::vector<S> f_101(numComp);
        std::vector<S> f_110(numComp);
        std::vector<S> f_111(numComp);
        std::vector<S> onFace, corners, scratch;

        if (m_faceOffset[0] > -1) {
#pragma omp for nowait
            for (index_t k=0; k<NodeIDsLeft.size(); k++) {
                borderNodeInfo tmp = NodeIDsLeft[k];
                // a corner of a boundary face can hang on a coarse edge in the
                // boundary plane; see gatherFaceCornersConstrained
                gatherFaceCornersConstrained(in, tmp, 0, numComp, sentinel,
                                             onFace, corners, scratch);
                memcpy(&f_000[0], &onFace[0*numComp], numComp*sizeof(S));
                memcpy(&f_001[0], &onFace[2*numComp], numComp*sizeof(S));
                memcpy(&f_010[0], &onFace[1*numComp], numComp*sizeof(S));
                memcpy(&f_011[0], &onFace[3*numComp], numComp*sizeof(S));
                S* o = out.getSampleDataRW(m_faceOffset[0]+k, sentinel);
                for (index_t i=0; i < numComp; ++i) {
                    o[INDEX2(i,numComp,0)] = (f_000[i] + f_001[i] + f_010[i] + f_011[i])/static_cast<S>(4);
                } // end of component loop i
            }
        }

        if (m_faceOffset[1] > -1) {
#pragma omp for nowait
            for (index_t k=0; k<NodeIDsRight.size(); k++) {
                borderNodeInfo tmp = NodeIDsRight[k];
                // a corner of a boundary face can hang on a coarse edge in the
                // boundary plane; see gatherFaceCornersConstrained
                gatherFaceCornersConstrained(in, tmp, 1, numComp, sentinel,
                                             onFace, corners, scratch);
                memcpy(&f_100[0], &onFace[0*numComp], numComp*sizeof(S));
                memcpy(&f_101[0], &onFace[2*numComp], numComp*sizeof(S));
                memcpy(&f_110[0], &onFace[1*numComp], numComp*sizeof(S));
                memcpy(&f_111[0], &onFace[3*numComp], numComp*sizeof(S));
                S* o = out.getSampleDataRW(m_faceOffset[1]+k, sentinel);
                for (index_t i=0; i < numComp; ++i) {
                    o[INDEX2(i,numComp,0)] = (f_100[i] + f_101[i] + f_110[i] + f_111[i])/static_cast<S>(4);
                } // end of component loop i
            }
        }
        if (m_faceOffset[2] > -1) {
#pragma omp for nowait
            for (index_t k=0; k<NodeIDsBottom.size(); k++) {
                borderNodeInfo tmp = NodeIDsBottom[k];
                // a corner of a boundary face can hang on a coarse edge in the
                // boundary plane; see gatherFaceCornersConstrained
                gatherFaceCornersConstrained(in, tmp, 2, numComp, sentinel,
                                             onFace, corners, scratch);
                memcpy(&f_000[0], &onFace[0*numComp], numComp*sizeof(S));
                memcpy(&f_001[0], &onFace[2*numComp], numComp*sizeof(S));
                memcpy(&f_100[0], &onFace[1*numComp], numComp*sizeof(S));
                memcpy(&f_101[0], &onFace[3*numComp], numComp*sizeof(S));
                S* o = out.getSampleDataRW(m_faceOffset[2]+k, sentinel);
                for (index_t i=0; i < numComp; ++i) {
                    o[INDEX2(i,numComp,0)] = (f_000[i] + f_001[i] + f_100[i] + f_101[i])/static_cast<S>(4);
                } // end of component loop i
            }
        }
        if (m_faceOffset[3] > -1) {
#pragma omp for nowait
            for (index_t k=0; k<NodeIDsTop.size(); k++) {
                borderNodeInfo tmp = NodeIDsTop[k];
                // a corner of a boundary face can hang on a coarse edge in the
                // boundary plane; see gatherFaceCornersConstrained
                gatherFaceCornersConstrained(in, tmp, 3, numComp, sentinel,
                                             onFace, corners, scratch);
                memcpy(&f_010[0], &onFace[0*numComp], numComp*sizeof(S));
                memcpy(&f_011[0], &onFace[2*numComp], numComp*sizeof(S));
                memcpy(&f_110[0], &onFace[1*numComp], numComp*sizeof(S));
                memcpy(&f_111[0], &onFace[3*numComp], numComp*sizeof(S));
                S* o = out.getSampleDataRW(m_faceOffset[3]+k, sentinel);
                for (index_t i=0; i < numComp; ++i) {
                    o[INDEX2(i,numComp,0)] = (f_010[i] + f_011[i] + f_110[i] + f_111[i])/static_cast<S>(4);
                } // end of component loop i
            }
        }
        if (m_faceOffset[4] > -1) {
#pragma omp for nowait
            for (index_t k=0; k<NodeIDsAbove.size(); k++) {
                borderNodeInfo tmp = NodeIDsAbove[k];
                // a corner of a boundary face can hang on a coarse edge in the
                // boundary plane; see gatherFaceCornersConstrained
                gatherFaceCornersConstrained(in, tmp, 4, numComp, sentinel,
                                             onFace, corners, scratch);
                memcpy(&f_000[0], &onFace[0*numComp], numComp*sizeof(S));
                memcpy(&f_010[0], &onFace[2*numComp], numComp*sizeof(S));
                memcpy(&f_100[0], &onFace[1*numComp], numComp*sizeof(S));
                memcpy(&f_110[0], &onFace[3*numComp], numComp*sizeof(S));
                S* o = out.getSampleDataRW(m_faceOffset[4]+k, sentinel);
                for (index_t i=0; i < numComp; ++i) {
                    o[INDEX2(i,numComp,0)] = (f_000[i] + f_010[i] + f_100[i] + f_110[i])/static_cast<S>(4);
                } // end of component loop i
            }
        }
        if (m_faceOffset[5] > -1) {
#pragma omp for nowait
            for (index_t k=0; k<NodeIDsBelow.size(); k++) {
                borderNodeInfo tmp = NodeIDsBelow[k];
                // a corner of a boundary face can hang on a coarse edge in the
                // boundary plane; see gatherFaceCornersConstrained
                gatherFaceCornersConstrained(in, tmp, 5, numComp, sentinel,
                                             onFace, corners, scratch);
                memcpy(&f_001[0], &onFace[0*numComp], numComp*sizeof(S));
                memcpy(&f_011[0], &onFace[2*numComp], numComp*sizeof(S));
                memcpy(&f_101[0], &onFace[1*numComp], numComp*sizeof(S));
                memcpy(&f_111[0], &onFace[3*numComp], numComp*sizeof(S));
                S* o = out.getSampleDataRW(m_faceOffset[5]+k, sentinel);
                for (index_t i=0; i < numComp; ++i) {
                    o[INDEX2(i,numComp,0)] = (f_001[i] + f_011[i] + f_101[i] + f_111[i])/static_cast<S>(4);
                } // end of component loop i
            } 
        } 
    } else {
        out.requireWrite();
        const S c0 = 0.044658198738520451079;
        const S c1 = 0.16666666666666666667;
        const S c2 = 0.62200846792814621559;

        std::vector<S> f_000(numComp);
        std::vector<S> f_001(numComp);
        std::vector<S> f_010(numComp);
        std::vector<S> f_011(numComp);
        std::vector<S> f_100(numComp);
        std::vector<S> f_101(numComp);
        std::vector<S> f_110(numComp);
        std::vector<S> f_111(numComp);
        std::vector<S> onFace, corners, scratch;

        //TODO fix tmp.neighbours indices below
        if (m_faceOffset[0] > -1) {
#pragma omp for nowait
            for (index_t k=0; k<NodeIDsLeft.size(); k++) {
                borderNodeInfo tmp = NodeIDsLeft[k];
                // a corner of a boundary face can hang on a coarse edge in the
                // boundary plane; see gatherFaceCornersConstrained
                gatherFaceCornersConstrained(in, tmp, 0, numComp, sentinel,
                                             onFace, corners, scratch);
                memcpy(&f_000[0], &onFace[0*numComp], numComp*sizeof(S));
                memcpy(&f_001[0], &onFace[2*numComp], numComp*sizeof(S));
                memcpy(&f_010[0], &onFace[1*numComp], numComp*sizeof(S));
                memcpy(&f_011[0], &onFace[3*numComp], numComp*sizeof(S));
                S* o = out.getSampleDataRW(m_faceOffset[0]+k, sentinel);
                for (index_t i=0; i < numComp; ++i) {
                    o[INDEX2(i,numComp,0)] = f_000[i]*c2 + f_011[i]*c0 + c1*(f_001[i] + f_010[i]);
                    o[INDEX2(i,numComp,1)] = f_001[i]*c0 + f_010[i]*c2 + c1*(f_000[i] + f_011[i]);
                    o[INDEX2(i,numComp,2)] = f_001[i]*c2 + f_010[i]*c0 + c1*(f_000[i] + f_011[i]);
                    o[INDEX2(i,numComp,3)] = f_000[i]*c0 + f_011[i]*c2 + c1*(f_001[i] + f_010[i]);
                } // end of component loop i
            }
        }
        if (m_faceOffset[1] > -1) {
    #pragma omp for nowait
            for (index_t k=0; k<NodeIDsRight.size(); k++) {
                borderNodeInfo tmp = NodeIDsRight[k];
                // a corner of a boundary face can hang on a coarse edge in the
                // boundary plane; see gatherFaceCornersConstrained
                gatherFaceCornersConstrained(in, tmp, 1, numComp, sentinel,
                                             onFace, corners, scratch);
                memcpy(&f_100[0], &onFace[0*numComp], numComp*sizeof(S));
                memcpy(&f_101[0], &onFace[2*numComp], numComp*sizeof(S));
                memcpy(&f_110[0], &onFace[1*numComp], numComp*sizeof(S));
                memcpy(&f_111[0], &onFace[3*numComp], numComp*sizeof(S));
                S* o = out.getSampleDataRW(m_faceOffset[1]+k, sentinel);
                for (index_t i=0; i < numComp; ++i) {
                    o[INDEX2(i,numComp,0)] = f_100[i]*c2 + f_111[i]*c0 + c1*(f_101[i] + f_110[i]);
                    o[INDEX2(i,numComp,1)] = f_101[i]*c0 + f_110[i]*c2 + c1*(f_100[i] + f_111[i]);
                    o[INDEX2(i,numComp,2)] = f_101[i]*c2 + f_110[i]*c0 + c1*(f_100[i] + f_111[i]);
                    o[INDEX2(i,numComp,3)] = f_100[i]*c0 + f_111[i]*c2 + c1*(f_101[i] + f_110[i]);
                } // end of component loop i
            }
        }
        if (m_faceOffset[2] > -1) {
    #pragma omp for nowait
             for (index_t k=0; k<NodeIDsBottom.size(); k++) {
                borderNodeInfo tmp = NodeIDsBottom[k];
                // a corner of a boundary face can hang on a coarse edge in the
                // boundary plane; see gatherFaceCornersConstrained
                gatherFaceCornersConstrained(in, tmp, 2, numComp, sentinel,
                                             onFace, corners, scratch);
                memcpy(&f_000[0], &onFace[0*numComp], numComp*sizeof(S));
                memcpy(&f_001[0], &onFace[2*numComp], numComp*sizeof(S));
                memcpy(&f_100[0], &onFace[1*numComp], numComp*sizeof(S));
                memcpy(&f_101[0], &onFace[3*numComp], numComp*sizeof(S));
                S* o = out.getSampleDataRW(m_faceOffset[2]+k, sentinel);
                for (index_t i=0; i < numComp; ++i) {
                    o[INDEX2(i,numComp,0)] = f_000[i]*c2 + f_101[i]*c0 + c1*(f_001[i] + f_100[i]);
                    o[INDEX2(i,numComp,1)] = f_001[i]*c0 + f_100[i]*c2 + c1*(f_000[i] + f_101[i]);
                    o[INDEX2(i,numComp,2)] = f_001[i]*c2 + f_100[i]*c0 + c1*(f_000[i] + f_101[i]);
                    o[INDEX2(i,numComp,3)] = f_000[i]*c0 + f_101[i]*c2 + c1*(f_001[i] + f_100[i]);
                } // end of component loop i
            }
        }
        if (m_faceOffset[3] > -1) {
    #pragma omp for nowait
            for (index_t k=0; k<NodeIDsTop.size(); k++) {
            borderNodeInfo tmp = NodeIDsTop[k];
                // a corner of a boundary face can hang on a coarse edge in the
                // boundary plane; see gatherFaceCornersConstrained
                gatherFaceCornersConstrained(in, tmp, 3, numComp, sentinel,
                                             onFace, corners, scratch);
                memcpy(&f_010[0], &onFace[0*numComp], numComp*sizeof(S));
                memcpy(&f_011[0], &onFace[2*numComp], numComp*sizeof(S));
                memcpy(&f_110[0], &onFace[1*numComp], numComp*sizeof(S));
                memcpy(&f_111[0], &onFace[3*numComp], numComp*sizeof(S));
                S* o = out.getSampleDataRW(m_faceOffset[3]+k, sentinel);
                for (index_t i=0; i < numComp; ++i) {
                    o[INDEX2(i,numComp,0)] = f_010[i]*c2 + f_111[i]*c0 + c1*(f_011[i] + f_110[i]);
                    o[INDEX2(i,numComp,1)] = f_011[i]*c0 + f_110[i]*c2 + c1*(f_010[i] + f_111[i]);
                    o[INDEX2(i,numComp,2)] = f_011[i]*c2 + f_110[i]*c0 + c1*(f_010[i] + f_111[i]);
                    o[INDEX2(i,numComp,3)] = f_010[i]*c0 + f_111[i]*c2 + c1*(f_011[i] + f_110[i]);
                } // end of component loop i
            }
        }
        if (m_faceOffset[4] > -1) {
    #pragma omp for nowait
            for (index_t k=0; k<NodeIDsAbove.size(); k++) {
            borderNodeInfo tmp = NodeIDsAbove[k];
                // a corner of a boundary face can hang on a coarse edge in the
                // boundary plane; see gatherFaceCornersConstrained
                gatherFaceCornersConstrained(in, tmp, 4, numComp, sentinel,
                                             onFace, corners, scratch);
                memcpy(&f_000[0], &onFace[0*numComp], numComp*sizeof(S));
                memcpy(&f_010[0], &onFace[2*numComp], numComp*sizeof(S));
                memcpy(&f_100[0], &onFace[1*numComp], numComp*sizeof(S));
                memcpy(&f_110[0], &onFace[3*numComp], numComp*sizeof(S));
                S* o = out.getSampleDataRW(m_faceOffset[4]+k, sentinel);
                for (index_t i=0; i < numComp; ++i) {
                    o[INDEX2(i,numComp,0)] = f_000[i]*c2 + f_110[i]*c0 + c1*(f_010[i] + f_100[i]);
                    o[INDEX2(i,numComp,1)] = f_010[i]*c0 + f_100[i]*c2 + c1*(f_000[i] + f_110[i]);
                    o[INDEX2(i,numComp,2)] = f_010[i]*c2 + f_100[i]*c0 + c1*(f_000[i] + f_110[i]);
                    o[INDEX2(i,numComp,3)] = f_000[i]*c0 + f_110[i]*c2 + c1*(f_010[i] + f_100[i]);
                } // end of component loop i
            }
        }
        if (m_faceOffset[5] > -1) {
    #pragma omp for nowait
            for (index_t k=0; k<NodeIDsBelow.size(); k++) {
            borderNodeInfo tmp = NodeIDsBelow[k];
                // a corner of a boundary face can hang on a coarse edge in the
                // boundary plane; see gatherFaceCornersConstrained
                gatherFaceCornersConstrained(in, tmp, 5, numComp, sentinel,
                                             onFace, corners, scratch);
                memcpy(&f_001[0], &onFace[0*numComp], numComp*sizeof(S));
                memcpy(&f_011[0], &onFace[2*numComp], numComp*sizeof(S));
                memcpy(&f_101[0], &onFace[1*numComp], numComp*sizeof(S));
                memcpy(&f_111[0], &onFace[3*numComp], numComp*sizeof(S));
                S* o = out.getSampleDataRW(m_faceOffset[5]+k, sentinel);
                for (index_t i=0; i < numComp; ++i) {
                    o[INDEX2(i,numComp,0)] = f_001[i]*c2 + f_111[i]*c0 + c1*(f_011[i] + f_101[i]);
                    o[INDEX2(i,numComp,1)] = f_011[i]*c0 + f_101[i]*c2 + c1*(f_001[i] + f_111[i]);
                    o[INDEX2(i,numComp,2)] = f_011[i]*c2 + f_101[i]*c0 + c1*(f_001[i] + f_111[i]);
                    o[INDEX2(i,numComp,3)] = f_001[i]*c0 + f_111[i]*c2 + c1*(f_011[i] + f_101[i]);
                } // end of component loop i
            }
        }
    }
}

//protected
inline dim_t Brick::getNumNodes() const
{
    // lnodes-based node count (owned + ghost). Replaces the coordinate-hash
    // NodeIDs.size(); for a conforming mesh the two counts agree.
    return nodes ? (dim_t) nodes->num_local_nodes : 0;
}

inline dim_t Brick::getNumHangingNodes() const
{
    return num_hanging;
}

//protected
inline dim_t Brick::getNumElements() const
{
    return nodes->num_local_elements;
}

unsigned Brick::forestChecksum() const
{
    // collective; the same value on every rank
    return p8est_checksum(p8est);
}

bool Brick::isConforming() const
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
bool Brick::getHangingNodes(p8est_lnodes_code_t face_code,
                            int masterCount[P8EST_CHILDREN],
                            int masters[P8EST_CHILDREN][4]) const
{
    // Decodes the lnodes face_code, following p8est_lnodes_decode: the low
    // P8EST_DIM bits hold the corner c that touches no hanging face, the next
    // three say which of c's faces hang, the three after that which of its edges
    // do. A corner in the middle of a coarse FACE is the average of that face's
    // four corners, one in the middle of a coarse EDGE the average of its two
    // ends - and those masters are all present among this element's own slots,
    // because a hanging slot already holds the far master. So the masters are
    // slot indices into this element and the weights are all 1/masterCount.
    for (int i = 0; i < P8EST_CHILDREN; ++i)
        masterCount[i] = 0;

    if (!face_code)
        return false;

    static const int ones = P8EST_CHILDREN - 1;        // 7
    const int c = (int) (face_code & ones);
    int work = (int) (face_code >> P8EST_DIM);

    for (int i = 0; i < P8EST_DIM; ++i) {             // hanging faces -> centre
        if (work & 1) {
            const int f = p8est_corner_faces[c][i];
            const int hn = c ^ ones ^ (1 << i);       // opposite c on face f
            masterCount[hn] = 4;
            for (int j = 0; j < 4; ++j)
                masters[hn][j] = p8est_face_corners[f][j];
        }
        work >>= 1;
    }
    for (int i = 0; i < P8EST_DIM; ++i) {             // hanging edges -> midpoint
        if (work & 1) {
            const int ed = p8est_corner_edges[c][i];
            const int hn = c ^ (1 << i);              // far end of edge ed from c
            masterCount[hn] = 2;
            for (int j = 0; j < 2; ++j)
                masters[hn][j] = p8est_edge_corners[ed][j];
        }
        work >>= 1;
    }
    return true;
}

MeshAccess Brick::getMeshAccess(bool materializeHanging) const
{
    MeshAccess m;
    m.numDim = 3;
    m.nodesPerElement = nodes->vnodes;                 // 8 for degree 1
    m.numNodes = nodes->num_local_nodes;
    m.numOwnedNodes = nodes->owned_count;
    m.numElements = nodes->num_local_elements;
    m.globalNodeOffset = (long) nodes->global_offset;
    m.globalElementOffset = (long) p8est->global_first_quadrant[m_mpiInfo->rank];
    m.numRealNodes = m.numNodes;
    m.mastersPerConstrainedNode = 4;                   // a face centre

    m.nodeCoords.assign((size_t) m.numNodes * m.numDim, 0.0);
    m.nodeLnodesId.resize(m.numNodes);
    m.elementNodes.resize((size_t) m.numElements * m.nodesPerElement);
    m.elementTags.resize(m.numElements);

    // node tags, indexed by local node exactly as m_nodeTags is. A domain that
    // has not been through populateSampleIds() leaves it short, so copy only
    // what is there; the nodes materialised later inherit via addHangingNode.
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
    // A HANGING corner has no node of its own: its element_nodes slot holds a
    // master, a node that lies OUTSIDE this element. So the slot's coordinate must
    // not be written - it would move the master to the hanging position - and with
    // materializeHanging the slot is redirected to a node created here.
    const int V = m.nodesPerElement;
    HangingNodeMap hangingOf;           // position -> materialised node index
    long e = 0;
    for (p4est_topidx_t treeid = p8est->first_local_tree;
         treeid <= p8est->last_local_tree; ++treeid) {
        p8est_tree_t * tree = p8est_tree_array_index(p8est->trees, treeid);
        sc_array_t * octs = &tree->quadrants;
        const p4est_locidx_t Q = (p4est_locidx_t) octs->elem_count;
        for (p4est_locidx_t q = 0; q < Q; ++q, ++e) {
            p8est_quadrant_t * oct = p8est_quadrant_array_index(octs, q);
            const octantData * od = (const octantData *) oct->p.user_data;
            m.elementTags[e] = od ? od->octantTag : 0;
            const p4est_qcoord_t len = P8EST_QUADRANT_LEN(oct->level);

            int masterCount[P8EST_CHILDREN], masterSlot[P8EST_CHILDREN][4];
            const bool anyHanging = getHangingNodes(nodes->face_code[e],
                                                   masterCount, masterSlot);

            double cornerPos[P8EST_CHILDREN][3];
            for (int c = 0; c < V; ++c) {
                const long ni = (long) nodes->element_nodes[(size_t) e * V + c];
                m.elementNodes[(size_t) e * V + c] = ni;
                const int cx = c & 1;          // z-order corner: bit0=x,bit1=y,bit2=z
                const int cy = (c >> 1) & 1;
                const int cz = (c >> 2) & 1;
                p8est_qcoord_to_vertex(p8est->connectivity, treeid,
                                       oct->x + cx * len, oct->y + cy * len,
                                       oct->z + cz * len, cornerPos[c]);
                if (anyHanging && masterCount[c] > 0)
                    continue;                  // a master, not this corner
                for (int d = 0; d < m.numDim; ++d)
                    m.nodeCoords[(size_t) ni * m.numDim + d] = cornerPos[c][d];
            }

            if (!materializeHanging || !anyHanging)
                continue;

            for (int c = 0; c < V; ++c) {
                if (masterCount[c] == 0)
                    continue;
                long masters[4];
                double weights[4];
                for (int k = 0; k < masterCount[c]; ++k) {
                    masters[k] = (long) nodes->element_nodes[
                            (size_t) e * V + masterSlot[c][k]];
                    weights[k] = 1. / masterCount[c];
                }
                m.elementNodes[(size_t) e * V + c] =
                        addHangingNode(m, hangingOf, cornerPos[c], masters,
                                       weights, masterCount[c]);
            }
        }
    }

    // Boundary faces. An octant face lies on the domain boundary when the octant
    // touches the tree boundary in that direction AND the connectivity sends
    // that tree face back to itself, which is p8est's encoding for "no
    // neighbour". Topological, unlike updateFaceElementCount().
    // NOTE conforming meshes only: a boundary face of a coarse octant may carry
    // hanging edge nodes, which are not represented here yet.
    {
        // face -> its four corners in z-order indexing, wound counter-clockwise
        // as seen from OUTSIDE, so the right-hand rule gives the outward normal.
        // Matches finley's Mesh_hex8. Faces are -x,+x,-y,+y,-z,+z.
        static const int faceCorner[6][4] = {
            {0,4,6,2}, {1,3,7,5}, {0,1,5,4}, {2,6,7,3}, {0,2,3,1}, {4,5,7,6} };
        static const long faceTag[6] = {1, 2, 10, 20, 100, 200};

        // Faces are collected per direction and concatenated afterwards,
        // because that is the order the domain gives its boundary SAMPLES:
        // m_faceOffset lays them out in blocks -x, +x, -y, +y, -z, +z, each in
        // octant order. Building them octant-major instead - as this did -
        // makes MeshAccess face f a different face from FunctionOnBoundary
        // sample f, which nothing notices until something pairs the two: the
        // boundary transfer read every value off the wrong face. Rectangle
        // carries the same rule; treat it as the invariant.
        std::vector<long> bucketNodes[6], bucketElements[6];
        const p8est_connectivity_t* conn = p8est->connectivity;
        m.nodesPerFace = 4;
        long le = 0;
        for (p4est_topidx_t treeid = p8est->first_local_tree;
             treeid <= p8est->last_local_tree; ++treeid) {
            p8est_tree_t* tree = p8est_tree_array_index(p8est->trees, treeid);
            sc_array_t* octs = &tree->quadrants;
            const p4est_locidx_t Q = (p4est_locidx_t) octs->elem_count;
            for (p4est_locidx_t q = 0; q < Q; ++q, ++le) {
                p8est_quadrant_t* oct = p8est_quadrant_array_index(octs, q);
                const p4est_qcoord_t len = P8EST_QUADRANT_LEN(oct->level);
                for (int f = 0; f < 6; ++f) {
                    if (conn->tree_to_tree[treeid * 6 + f] != treeid ||
                        conn->tree_to_face[treeid * 6 + f] != f)
                        continue;               // a neighbouring tree is there
                    bool touches;
                    switch (f) {
                        case 0:  touches = (oct->x == 0); break;
                        case 1:  touches = (oct->x + len == P8EST_ROOT_LEN); break;
                        case 2:  touches = (oct->y == 0); break;
                        case 3:  touches = (oct->y + len == P8EST_ROOT_LEN); break;
                        case 4:  touches = (oct->z == 0); break;
                        default: touches = (oct->z + len == P8EST_ROOT_LEN); break;
                    }
                    if (!touches)
                        continue;
                    // from m.elementNodes, not element_nodes: a boundary face of
                    // a fine octant can have a hanging corner (its edge midpoint
                    // hangs on the coarse neighbour), and the face must name the
                    // node that is actually there
                    for (int c = 0; c < 4; ++c)
                        bucketNodes[f].push_back(
                                m.elementNodes[(size_t) le * V + faceCorner[f][c]]);
                    bucketElements[f].push_back(le);
                }
            }
        }
        for (int f = 0; f < 6; ++f) {
            m.faceNodes.insert(m.faceNodes.end(), bucketNodes[f].begin(),
                               bucketNodes[f].end());
            m.faceElements.insert(m.faceElements.end(),
                                  bucketElements[f].begin(),
                                  bucketElements[f].end());
            m.faceTags.insert(m.faceTags.end(), bucketElements[f].size(),
                              faceTag[f]);
            m.faceDirections.insert(m.faceDirections.end(),
                                    bucketElements[f].size(), (long) f);
        }
        m.numFaces = (long) m.faceTags.size();
    }

    std::vector<long> ownedPerRank(m_mpiInfo->size);
    for (int i = 0; i < m_mpiInfo->size; ++i)
        ownedPerRank[i] = (long) nodes->global_owned_count[i];
    finaliseNodeNumbering(m, ownedPerRank);

    return m;
}

//protected
dim_t Brick::getNumFaceElements() const
{
    return m_faceCount[0]+m_faceCount[1]
          +m_faceCount[2]+m_faceCount[3]
          +m_faceCount[4]+m_faceCount[5];
}

//protected
inline dim_t Brick::getNumDOF() const
{
    // owned nodes only (each owned node is one real DOF). Ghost/shared nodes
    // are columns, not rows. (MPI: A6.)
    return nodes ? (dim_t) nodes->owned_count : 0;
}

void Brick::updateTreeIDs()
{
    treeIDs.clear();
    for(p8est_topidx_t treeid = p8est->first_local_tree; treeid <= p8est->last_local_tree; ++treeid) {
        p8est_tree_t * tree = p8est_tree_array_index(p8est->trees, treeid);
        sc_array_t * tquadrants = &tree->quadrants;
        p8est_qcoord_t Q = (p8est_qcoord_t) tquadrants->elem_count;
#pragma omp parallel for
        for(int q = 0; q < Q; ++q) { // Loop over the elements attached to the tree
            p8est_quadrant_t * quad = p8est_quadrant_array_index(tquadrants, q);
            treeIDs[std::make_tuple(quad->x,quad->y,quad->z)]=treeid;
        }
    }
}

void Brick::updateRowsColumns()
{
    // Conforming lnodes numbering: every local node is a real DOF numbered
    // identically (serially). ownSample() consults m_dofMap, so keep it the
    // identity map. The legacy hanging-node connectivity build is retired;
    // getConnections now derives the matrix graph directly from lnodes.
    const dim_t n = getNumNodes();
    m_dofMap.assign(n, 0);
    for(dim_t i = 0; i < n; ++i)
        m_dofMap[i] = i;
}

#ifdef ESYS_HAVE_TRILINOS
//protected
esys_trilinos::TrilinosGraph_ptr Brick::getTrilinosGraph() const
{   
    if (m_graph.is_null()) {
        m_graph = createTrilinosGraph(myRows, myColumns);
    }
    return m_graph;
}
#endif

#ifdef ESYS_HAVE_PASO
//protected
paso::SystemMatrixPattern_ptr Brick::getPasoMatrixPattern(
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
void Brick::populateSampleIds()
{
    m_nodeTags.assign(getNumNodes(), 0);
    updateTagsInUse(Nodes);

    m_elementTags.assign(getNumElements(), 0);
    updateTagsInUse(Elements);
}

void Brick::updateFaceElementCount()
{
    oxleytimer.toc("updateFaceElementCount");

    // Build the six boundary-face element lists directly from the lnodes
    // connectivity. A face element stores in neighbours[0..3] the four nodes of
    // the quad face it lies on, ordered (a0b0, a1b0, a0b1, a1b1) in the face's
    // two in-plane axes (a,b) -- the order interpolateNodesOnFaces expects.
    //
    // Face id : boundary        : (a,b) axes : quad z-order corners
    //   0 Left   x = xmin (-x)   : (y,z)      : 0,2,4,6
    //   1 Right  x = xmax (+x)   : (y,z)      : 1,3,5,7
    //   2 Bottom y = ymin (-y)   : (x,z)      : 0,1,4,5
    //   3 Top    y = ymax (+y)   : (x,z)      : 2,3,6,7
    //   4 (z-min, normal -z)     : (x,y)      : 0,1,2,3
    //   5 (z-max, normal +z)     : (x,y)      : 4,5,6,7
    // faceCornerSlots, at the top of this file, is that table.

    for(int i = 0; i < 6; ++i)
        m_faceCount[i] = 0;
    NodeIDsLeft.clear();   NodeIDsRight.clear();
    NodeIDsBottom.clear(); NodeIDsTop.clear();
    NodeIDsAbove.clear();  NodeIDsBelow.clear();
    std::vector<borderNodeInfo>* lists[6] = {
        &NodeIDsLeft, &NodeIDsRight, &NodeIDsBottom, &NodeIDsTop,
        &NodeIDsAbove, &NodeIDsBelow };

    const int V = nodes->vnodes;   // 8
    const double x0 = forestData->m_origin[0], y0 = forestData->m_origin[1], z0 = forestData->m_origin[2];
    const double x1 = forestData->m_lxyz[0],   y1 = forestData->m_lxyz[1],   z1 = forestData->m_lxyz[2];

    long e = 0;
    for(p8est_topidx_t treeid = p8est->first_local_tree; treeid <= p8est->last_local_tree; ++treeid)
    {
        p8est_tree_t * tree = p8est_tree_array_index(p8est->trees, treeid);
        sc_array_t * tquadrants = &tree->quadrants;
        p8est_locidx_t Q = (p8est_locidx_t) tquadrants->elem_count;
        for(int q = 0; q < Q; ++q, ++e)
        {
            p8est_quadrant_t * quad = p8est_quadrant_array_index(tquadrants, q);
            const p8est_qcoord_t len = P8EST_QUADRANT_LEN(quad->level);
            double lo[3], hi[3];
            p8est_qcoord_to_vertex(p8est->connectivity, treeid, quad->x,       quad->y,       quad->z,       lo);
            p8est_qcoord_to_vertex(p8est->connectivity, treeid, quad->x+len,   quad->y+len,   quad->z+len,   hi);

            const bool onFace[6] = {
                lo[0]==x0, hi[0]==x1, lo[1]==y0, hi[1]==y1, lo[2]==z0, hi[2]==z1 };

            for(int fc = 0; fc < 6; ++fc)
            {
                if(!onFace[fc]) continue;
                borderNodeInfo tmp;
                for(int c = 0; c < 4; ++c)
                    tmp.neighbours[c] = (int) nodes->element_nodes[(size_t) e * V + faceCornerSlots[fc][c]];
                tmp.nodeid  = tmp.neighbours[0];
                tmp.x = quad->x; tmp.y = quad->y; tmp.z = quad->z;
                tmp.level = quad->level; tmp.treeid = treeid;
                // the way back to this element's face_code, which says which of
                // the neighbours above are hanging slots holding a master
                tmp.quadIndex = e;
                lists[fc]->push_back(tmp);
                m_faceCount[fc]++;
            }
        }
    }

    // face sample offsets and face tags (x-/x+/y-/y+/z-/z+ -> 1/2/10/20/100/200)
    const index_t faceTag[6] = { 1, 2, 10, 20, 100, 200 };
    m_faceTags.clear();
    m_faceOffset.assign(6, -1);
    index_t offset = 0;
    for(int i = 0; i < 6; ++i) {
        if(m_faceCount[i] > 0) {
            m_faceOffset[i] = offset;
            offset += m_faceCount[i];
            m_faceTags.insert(m_faceTags.end(), m_faceCount[i], faceTag[i]);
        }
    }
    const dim_t NFE = getNumFaceElements();
    m_faceId.resize(NFE);
    for(dim_t k = 0; k < NFE; ++k) m_faceId[k] = k;
    oxleytimer.toc("updateFaceElementCount... done");
}

// This is a wrapper that converts the p8est node information into an IndexVector
IndexVector Brick::getNodeDistribution() const
{
    return m_nodeDistribution;
}

// This is a wrapper that converts the p8est node information into an IndexVector
void Brick::updateNodeDistribution() 
{
    // TODO

    // m_nodeDistribution.clear();
    // m_nodeDistribution.assign(MAXP4ESTNODES,0);

    // int counter =0;
    // for(p8est_topidx_t treeid = p8est->first_local_tree; treeid <= p8est->last_local_tree; ++treeid) 
    // {
    //     p8est_tree_t * tree = p8est_tree_array_index(p8est->trees, treeid);
    //     sc_array_t * tquadrants = &tree->quadrants;
    //     p8est_locidx_t Q = (p8est_locidx_t) tquadrants->elem_count;
    //     for(int q = 0; q < Q; ++q) 
    //     { 
    //         p8est_quadrant_t * quad = p8est_quadrant_array_index(tquadrants, q);
    //         p8est_qcoord_t length = P8EST_QUADRANT_LEN(quad->level);
    //         for(int n = 0; n < 4; n++)
    //         {
    //             double lx = length * ((int) (n % 2) == 1);
    //             double ly = length * ((int) (n / 2) == 1);
    //             double lz = length * ((int) (n / 2) == 1); //TODO
    //             double xy[3];
    //             p8est_qcoord_to_vertex(p8est->connectivity, treeid, quad->x+lx, quad->y+ly, quad->z+lz, xy);
    //             m_nodeDistribution[counter++]=NodeIDs.find(std::make_tuple(xy[0],xy[1],xy[2]))->second;
    //         }
    //     }
    // }
    // m_nodeDistribution.shrink_to_fit();
}

// updates m_elementIDs()
void Brick::updateElementIds()
{
    oxleytimer.toc("updateElementIds");
    // element sample ids are the running local leaf indices [0, numElements)
    const dim_t ne = getNumElements();
    m_elementId.clear();
    m_elementId.resize(ne);
    for(dim_t i = 0; i < ne; ++i)
        m_elementId[i] = i;
    m_elementId.shrink_to_fit();
    oxleytimer.toc("done");
}

std::vector<IndexVector> Brick::getConnections(bool includeShared) const
{
    // returns a vector v of size numNodes where v[i] lists the node indices
    // coupled to node i (the occupied local matrix columns of row i).
    //
    // Built directly from the lnodes element->node connectivity (no coordinate
    // hashing): every corner of a leaf is coupled to every other corner. This
    // replaces the update_RC / coordinate-hash machinery in updateRowsColumns.

    long numNodes = getNumNodes();
    std::vector<IndexVector> indices_out(numNodes);

    const int V = nodes->vnodes;   // 8 corners (degree-1 lnodes)
    long e = 0;                    // running local leaf index (lnodes order)
    for(p8est_topidx_t treeid = p8est->first_local_tree; treeid <= p8est->last_local_tree; ++treeid)
    {
        p8est_tree_t * tree = p8est_tree_array_index(p8est->trees, treeid);
        sc_array_t * tquadrants = &tree->quadrants;
        p8est_locidx_t Q = (p8est_locidx_t) tquadrants->elem_count;
        for(int q = 0; q < Q; ++q, ++e)
        {
            long lni[8];
            for(int i = 0; i < V; i++)
                lni[i] = (long) nodes->element_nodes[(size_t) e * V + i];

            for(int i = 0; i < V; i++)
            {
                for(int j = 0; j < V; j++)
                {
                    bool dup = false;
                    for(int k = 0; k < indices_out[lni[i]].size(); k++)
                        if(indices_out[lni[i]][k] == lni[j])
                        {
                            dup = true;
                            break;
                        }
                    if(dup == false)
                        indices_out[lni[i]].push_back(lni[j]);
                }
            }
        }
    }

    // MPI: add couplings from the ghost octant halo for OWNED rows only, so the
    // Trilinos colMap/graph cover the 2nd-layer columns of owned boundary nodes.
    const long nDOF = getNumDOF();
    const long nGhost = (long) m_ghostElemNodes.size() / (V ? V : 1);
    for(long g = 0; g < nGhost; ++g)
    {
        const index_t* lni = &m_ghostElemNodes[(size_t) g * V];
        for(int i = 0; i < V; i++)
        {
            const long row = (long) lni[i];
            if(row >= nDOF)
                continue;
            for(int j = 0; j < V; j++)
            {
                const index_t col = lni[j];
                bool dup = false;
                for(int k = 0; k < indices_out[row].size(); k++)
                    if(indices_out[row][k] == col) { dup = true; break; }
                if(!dup)
                    indices_out[row].push_back(col);
            }
        }
    }

    for(long i = 0; i < numNodes; i++)
        std::sort(indices_out[i].begin(), indices_out[i].end());

    return indices_out;
}

bool Brick::operator==(const AbstractDomain& other) const
{
    const Brick* o=dynamic_cast<const Brick*>(&other);
    if (o) {
        return ((p8est_checksum(p8est) == p8est_checksum(o->p8est)));
    }
    return false;
}

//protected
void Brick::assembleGradient(escript::Data& out,
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
void Brick::assembleGradientImpl(escript::Data& out,
                                 const escript::Data& in) const
{
    const dim_t numComp = in.getDataPointSize();
    const Scalar zero = static_cast<Scalar>(0);

    // reference shape gradients at the 8 Gauss points (2x2x2) and at the centre.
    const double SQ = 1.0/std::sqrt(3.0);
    const double chi = (1.0 + SQ)/2.0, clo = (1.0 - SQ)/2.0;
    double gradRef[8][8][3];   // node a, gauss g, direction
    double gradCtr[8][3];      // node a, direction (centre)
    for(int a=0;a<8;++a){
        const int pa=a&1, qa=(a>>1)&1, ra=(a>>2)&1;
        const double sx=pa?1.0:-1.0, sy=qa?1.0:-1.0, sz=ra?1.0:-1.0;
        for(int g=0;g<8;++g){
            const int gi=(g>>2)&1, gj=(g>>1)&1, gk=g&1;  // interp gauss order (x slowest)
            const double Xa=(pa==gi)?chi:clo, Ya=(qa==gj)?chi:clo, Za=(ra==gk)?chi:clo;
            gradRef[a][g][0]=sx*Ya*Za; gradRef[a][g][1]=Xa*sy*Za; gradRef[a][g][2]=Xa*Ya*sz;
        }
        gradCtr[a][0]=sx*0.25; gradCtr[a][1]=sy*0.25; gradCtr[a][2]=sz*0.25;
    }

    const int fsType = out.getFunctionSpace().getTypeCode();
    const bool reduced = (fsType == ReducedElements || fsType == ReducedFaceElements);
    const bool onBoundary = (fsType == FaceElements || fsType == ReducedFaceElements);
    if (fsType != Elements && fsType != ReducedElements && !onBoundary)
        throw ValueError("Gradient: unsupported function space.");
    out.requireWrite();

    // On the fine side of a seam a corner slot holds a master rather than the
    // value at that corner, so the values are read through
    // gatherCornersConstrained rather than straight out of the samples.
    std::vector<Scalar> corners, scratch;

    if (!onBoundary) {
        long e = 0;
        for(p8est_topidx_t treeid = p8est->first_local_tree; treeid <= p8est->last_local_tree; ++treeid)
        {
            p8est_tree_t * tree = p8est_tree_array_index(p8est->trees, treeid);
            sc_array_t * tquadrants = &tree->quadrants;
            p8est_locidx_t Q = (p8est_locidx_t) tquadrants->elem_count;
            for(int q = 0; q < Q; ++q, ++e)
            {
                p8est_quadrant_t * quad = p8est_quadrant_array_index(tquadrants, q);
                const double hh = (double)(1 << quad->level);
                const double h[3] = { m_NX[0]/hh, m_NX[1]/hh, m_NX[2]/hh };

                // node values of the input field (indexed by lnode id == dof id serially)
                gatherCornersConstrained(in, e, numComp, zero, corners, scratch);
                const Scalar* f[8];
                for(int a=0;a<8;++a)
                    f[a] = &corners[(size_t) a*numComp];

                Scalar* o = out.getSampleDataRW(e, zero);
                if(reduced) {
                    for(index_t i=0;i<numComp;++i)
                        for(int d=0;d<3;++d) {
                            Scalar s = zero;
                            for(int a=0;a<8;++a) s += f[a][i]*(gradCtr[a][d]/h[d]);
                            o[INDEX2(i,d,numComp)] = s;
                        }
                } else {
                    for(int g=0;g<8;++g)
                        for(index_t i=0;i<numComp;++i)
                            for(int d=0;d<3;++d) {
                                Scalar s = zero;
                                for(int a=0;a<8;++a) s += f[a][i]*(gradRef[a][g][d]/h[d]);
                                o[INDEX3(i,d,g,numComp,3)] = s;
                            }
                }
            }
        }
    } else {
        // The boundary faces. There was no branch for them at all: the loop
        // above wrote one element-shaped sample per octant into a Data holding
        // one face-shaped sample per boundary face, which ran off the end of the
        // buffer and corrupted the heap - reproducibly, on a uniform forest.
        //
        // A face's quadrature points lie ON the face, so the trilinear
        // interpolant is evaluated with the normal coordinate fixed at the face
        // (0 or 1). That is the same convention Rectangle uses in 2D, where it
        // shows up as the tangential derivative being a plain difference of the
        // two values on the face.
        //
        // The point order is the one interpolateNodesOnFaces writes, so that
        // getX() and grad() on FunctionOnBoundary speak about the same points:
        // g = ia + 2*ib over the face's two in-plane axes in increasing order,
        // each running (clo, chi). A reduced face has the single centre point.
        const std::vector<borderNodeInfo>* lists[6] = {
            &NodeIDsLeft, &NodeIDsRight, &NodeIDsBottom, &NodeIDsTop,
            &NodeIDsAbove, &NodeIDsBelow };
        const int numPoints = reduced ? 1 : 4;

        for(int fc = 0; fc < 6; ++fc)
        {
            if(m_faceOffset[fc] < 0)
                continue;
            const int nrm = fc / 2;                     // axis normal to the face
            const double nu = (fc % 2) ? 1.0 : 0.0;     // which end of that axis
            const int axisA = (nrm == 0) ? 1 : 0;       // the other two, in order
            const int axisB = (nrm == 2) ? 1 : 2;
            const std::vector<borderNodeInfo>& L = *lists[fc];

            for(index_t k = 0; k < (index_t) L.size(); ++k)
            {
                const borderNodeInfo& b = L[k];
                if(b.quadIndex < 0)
                    throw OxleyException("Gradient: the boundary face list does "
                            "not say which element a face belongs to.");
                const double hh = (double)(1 << b.level);
                const double h[3] = { m_NX[0]/hh, m_NX[1]/hh, m_NX[2]/hh };

                gatherCornersConstrained(in, b.quadIndex, numComp, zero,
                                         corners, scratch);

                Scalar* o = out.getSampleDataRW(m_faceOffset[fc]+k, zero);
                for(int g = 0; g < numPoints; ++g)
                {
                    double xi[3];
                    xi[nrm] = nu;
                    if(reduced) {
                        xi[axisA] = xi[axisB] = 0.5;
                    } else {
                        xi[axisA] = (g & 1) ? chi : clo;
                        xi[axisB] = (g & 2) ? chi : clo;
                    }
                    double dN[8][3];
                    for(int a = 0; a < 8; ++a) {
                        const int pa=a&1, qa=(a>>1)&1, ra=(a>>2)&1;
                        const double X = pa ? xi[0] : 1.0-xi[0];
                        const double Y = qa ? xi[1] : 1.0-xi[1];
                        const double Z = ra ? xi[2] : 1.0-xi[2];
                        dN[a][0] = (pa?1.0:-1.0)*Y*Z;
                        dN[a][1] = X*(qa?1.0:-1.0)*Z;
                        dN[a][2] = X*Y*(ra?1.0:-1.0);
                    }
                    for(index_t i = 0; i < numComp; ++i)
                        for(int d = 0; d < 3; ++d) {
                            Scalar s = zero;
                            for(int a = 0; a < 8; ++a)
                                s += corners[(size_t) a*numComp+i]*(dN[a][d]/h[d]);
                            if(reduced)
                                o[INDEX2(i,d,numComp)] = s;
                            else
                                o[INDEX3(i,d,g,numComp,3)] = s;
                        }
                }
            }
        }
    }
}



//private
template<typename Scalar>
void Brick::addToMatrixAndRHS(escript::AbstractSystemMatrix* S, escript::Data& F,
         const std::vector<Scalar>& EM_S, const std::vector<Scalar>& EM_F, 
         bool addS, bool addF, index_t e, index_t t, int nEq, int nComp) const
{    
    const int V = nodes->vnodes;   // 8 corners (z-order matches the kernels)
    IndexVector rowIndex(V);
    p8est_tree_t * currenttree = p8est_tree_array_index(p8est->trees, t);
    // global local-leaf index in lnodes order; quadrants_offset is the cumulative
    // number of local octants in the trees before t.
    const long g = (long) currenttree->quadrants_offset + (long) e;
    for(int i = 0; i < V; i++)
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


//protected
void Brick::nodesToDOF(escript::Data& out, const escript::Data& in) const
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

//     const index_t left = (m_offset[0]==0 ? 0 : 1);
//     const index_t bottom = (m_offset[1]==0 ? 0 : 1);
//     const index_t front = (m_offset[2]==0 ? 0 : 1);
//     const dim_t nDOF0 = (m_gNE[0]+1)/m_NX[0];
//     const dim_t nDOF1 = (m_gNE[1]+1)/m_NX[1];
//     const dim_t nDOF2 = (m_gNE[2]+1)/m_NX[2];
// #pragma omp parallel for
//     for (index_t i=0; i<nDOF2; i++) {
//         for (index_t j=0; j<nDOF1; j++) {
//             for (index_t k=0; k<nDOF0; k++) {
//                 const index_t n=k+left+(j+bottom)*m_NN[0]+(i+front)*m_NN[0]*m_NN[1];
//                 const double* src=in.getSampleDataRO(n);
//                 std::copy(src, src+numComp, out.getSampleDataRW(k+j*nDOF0+i*nDOF0*nDOF1));
//             }
//         }
//     }
}

// Updates m_faceOffset for each quadrant
void Brick::updateFaceOffset()
{
    oxleytimer.toc("updateFaceOffset");
    p8est_iterate(p8est, NULL, NULL, update_node_faceoffset, NULL, NULL, NULL);
    oxleytimer.toc("done");
}

void Brick::updateMeshInformation()
{
    refineMesh("MARE2DEM");
}

////////////////////////////// inline methods ////////////////////////////////
inline dim_t Brick::getDofOfNode(dim_t node) const
{
    // Conforming lnodes numbering: every node is a real DOF and (serially)
    // the DOF id equals the local node id. MPI ownership handled in A6.
    return node;
}

dim_t Brick::findNode(const double *coords) const
{
    // Search the lnodes node coordinates for the closest node (used for Dirac
    // points). Replaces the coordinate-hash lookup.
    // reject points outside the domain bounding box (out-of-range Dirac points)
    const double x0=forestData->m_origin[0], y0=forestData->m_origin[1], z0=forestData->m_origin[2];
    const double x1=forestData->m_lxyz[0],   y1=forestData->m_lxyz[1],   z1=forestData->m_lxyz[2];
    double ext = x1-x0; if(y1-y0>ext) ext=y1-y0; if(z1-z0>ext) ext=z1-z0;
    const double tol = 1e-8*ext;
    if(coords[0]<x0-tol || coords[0]>x1+tol || coords[1]<y0-tol || coords[1]>y1+tol
       || coords[2]<z0-tol || coords[2]>z1+tol)
        return -1;
    const MeshAccess m = getMeshAccess();
    long closest = 0;
    double best = std::numeric_limits<double>::max();
    for(long i = 0; i < m.numNodes; ++i)
    {
        const double dx = m.nodeCoords[(size_t) i*3 + 0] - coords[0];
        const double dy = m.nodeCoords[(size_t) i*3 + 1] - coords[1];
        const double dz = m.nodeCoords[(size_t) i*3 + 2] - coords[2];
        const double d2 = dx*dx + dy*dy + dz*dz;
        if(d2 < best) { best = d2; closest = i; }
    }
    return (dim_t) closest;
}

// adds the dirac points and tags 
void Brick::addPoints(const std::vector<double>& coords, const std::vector<int>& tags)
{
    // A Dirac point must be claimed by exactly ONE rank in MPI: find this rank's
    // nearest OWNED node (first getNumDOF() local nodes) and distance, then keep
    // the point only on the globally-nearest rank (ties broken by lowest rank).
    // See Rectangle::addPoints. (A6.)
    const dim_t nOwned = getNumDOF();
    const MeshAccess m = getMeshAccess();
    const double x0=forestData->m_origin[0], y0=forestData->m_origin[1], z0=forestData->m_origin[2];
    const double x1=forestData->m_lxyz[0],   y1=forestData->m_lxyz[1],   z1=forestData->m_lxyz[2];
    double ext = x1-x0; if(y1-y0>ext) ext=y1-y0; if(z1-z0>ext) ext=z1-z0;
    const double tol = 1e-8*ext;

    for (int i = 0; i < (int)tags.size(); i++) {
        const double px = coords[i*m_numDim + 0];
        const double py = coords[i*m_numDim + 1];
        const double pz = coords[i*m_numDim + 2];

        double best = std::numeric_limits<double>::max();
        long bestNode = -1;
        if (!(px<x0-tol || px>x1+tol || py<y0-tol || py>y1+tol || pz<z0-tol || pz>z1+tol)) {
            for (long n = 0; n < nOwned; ++n) {
                const double dx = m.nodeCoords[(size_t)n*3 + 0] - px;
                const double dy = m.nodeCoords[(size_t)n*3 + 1] - py;
                const double dz = m.nodeCoords[(size_t)n*3 + 2] - pz;
                const double d2 = dx*dx + dy*dy + dz*dz;
                if (d2 < best) { best = d2; bestNode = n; }
            }
        }

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

static inline void
brick_linear_to_xyz (p8est_topidx_t ti, const int logx[P8EST_DIM],
                     const int rankx[P8EST_DIM], p8est_topidx_t tx[P8EST_DIM])
{
    int i, j, k;
    int lastlog = 0;

    for (i = 0; i < P8EST_DIM; i++) {
        tx[i] = 0;
    }

    for (i = 0; i < P8EST_DIM - 1; i++) {
        p8est_topidx_t tempx[3] = { 0, 0, 0 };
        int logi = logx[rankx[i]] - lastlog;
        int idx[3] = { -1, -1, -1 };
        int c = 0;

        for (k = 0; k < P8EST_DIM - i; k++) {
            int d = rankx[i + k];
            idx[d] = 0;
        }
    
        for (k = 0; k < P8EST_DIM; k++) {
            if (idx[k] == 0) {
                idx[k] = c++;
            }
        }

        for (j = 0; j < logi; j++) {
            int base = (P8EST_DIM - i) * j;
            int shift = (P8EST_DIM - i - 1) * j;

            for (k = 0; k < P8EST_DIM; k++) {
                int id = idx[k];

                if (id >= 0) {
                    tempx[k] |= (ti & (1 << (base + id))) >> (shift + id);
                }
            }
        }

        for (k = 0; k < P8EST_DIM; k++) {
            tx[k] += (tempx[k] << lastlog);
        }
        lastlog += logi;
        ti >>= (P8EST_DIM - i) * logi;
    }
    tx[rankx[P8EST_DIM - 1]] += (ti << lastlog);
}

static inline p8est_topidx_t
brick_xyz_to_linear (const p8est_topidx_t tx[P8EST_DIM],
                     const int logx[P8EST_DIM], const int rankx[P8EST_DIM])
{
    int i, j, k;
    int lastlog = logx[rankx[P8EST_DIM - 2]];
    p8est_topidx_t ti = tx[rankx[P8EST_DIM - 1]] >> lastlog;

    for (i = P8EST_DIM - 2; i >= 0; i--) {
        p8est_topidx_t tempx[3] = { 0, 0, 0 };
        int logi =  (i == 0) ? lastlog : lastlog - logx[rankx[i - 1]];
        int idx[3] = { -1, -1, -1 };
        int c = 0;

        for (k = 0; k < P8EST_DIM - i; k++) {
            int d = rankx[i + k];
            idx[d] = 0;
        }

        for (k = 0; k < P8EST_DIM; k++) {
            if (idx[k] == 0) {
                idx[k] = c++;
            }
        }

        ti <<= (P8EST_DIM - i) * logi;
        lastlog -= logi;
        for (k = 0; k < P8EST_DIM; k++) {
            tempx[k] = tx[k] >> lastlog;
        }
        for (j = 0; j < logi; j++) {
            int shift = (P8EST_DIM - i - 1) * j;

            for (k = 0; k < P8EST_DIM; k++) {
                int id = idx[k];

                if (id >= 0) {
                    ti |= (tempx[k] & (1 << j)) << (shift + id);
                }
            }
        }
    }
    return ti;
}

p8est_connectivity_t *
Brick::new_brick_connectivity (int n0, int n1, int n2, int periodic_a, int periodic_b, int periodic_c,
                               double x0, double x1, double y0, double y1, double z0, double z1)
{
    const p8est_topidx_t m = (p8est_topidx_t) n0;
    const p8est_topidx_t n = (p8est_topidx_t) n1;
    const p8est_topidx_t p = (p8est_topidx_t) n2;
    ESYS_ASSERT(m > 0 && n > 0 && p > 0, "n0, n1 and n2 must be greater than zero.");

    const p8est_topidx_t mc = periodic_a ? m : (m - 1);
    const p8est_topidx_t nc = periodic_b ? n : (n - 1);
    const p8est_topidx_t pc = periodic_c ? p : (p - 1);
    const p8est_topidx_t num_trees = m * n * p;
    ESYS_ASSERT(num_trees <= MAXTREES ,"n0*n1*n2 must be less than MAXTREES.");

    const p8est_topidx_t num_corners = mc * nc * pc;
    const p8est_topidx_t num_ctt = P8EST_CHILDREN * num_corners;
    const p8est_topidx_t num_edges = m * nc * pc + mc * n * pc + mc * nc * p;
    const p8est_topidx_t num_ett = 4 * num_edges;
    const p8est_topidx_t num_vertices = (m + 1) * (n + 1) * (p + 1);
    const int periodic[P8EST_DIM] = { periodic_a, periodic_b, periodic_c };
    const p8est_topidx_t max[P8EST_DIM] = { m - 1, n - 1, p - 1 };

    p8est_topidx_t  n_iter;
    int logx[P8EST_DIM];
    int rankx[P8EST_DIM];
    int i, j, l;
    p8est_topidx_t  tf[P8EST_FACES], tc[P8EST_CHILDREN];
    p8est_topidx_t  coord[P8EST_DIM], coord2[P8EST_DIM], ttemp;
    p8est_topidx_t *linear_to_tree;
    p8est_topidx_t *tree_to_corner2;
    p8est_topidx_t  vcount = 0, vicount = 0;
    int c[P8EST_DIM];
    p8est_connectivity_t *conn;
    p8est_topidx_t tl;
    p8est_topidx_t tz;
    p8est_topidx_t te[P8EST_EDGES];
    p8est_topidx_t *tree_to_edge2;
    int dir1, dir2;

    // Size of the grid spacing
    double dx = (x1 - x0) / n0;
    double dy = (y1 - y0) / n1;
    double dz = (z1 - z0) / n2;

    conn = p8est_connectivity_new (num_vertices, num_trees, 
                                 num_edges, num_ett,
                                 num_corners, num_ctt);

    double * vertices = conn->vertices;
    p8est_topidx_t * tree_to_vertex = conn->tree_to_vertex;
    p8est_topidx_t * tree_to_tree = conn->tree_to_tree;
    int8_t * tree_to_face = conn->tree_to_face;
    p8est_topidx_t * tree_to_edge = conn->tree_to_edge;
    p8est_topidx_t * ett_offset = conn->ett_offset;
    p8est_topidx_t * edge_to_tree = conn->edge_to_tree;
    int8_t * edge_to_edge = conn->edge_to_edge;
    p8est_topidx_t * tree_to_corner = conn->tree_to_corner;
    p8est_topidx_t * ctt_offset = conn->ctt_offset;
    p8est_topidx_t * corner_to_tree = conn->corner_to_tree;
    int8_t * corner_to_corner = conn->corner_to_corner;

    p8est_topidx_t  ti, tj, tk;
#pragma omp parallel for
    for (ti = 0; ti < num_edges + 1; ti++) {
        ett_offset[ti] = 4 * ti;
    }

#pragma omp parallel for
    for (ti = 0; ti < num_corners + 1; ti++) {
        ctt_offset[ti] = P8EST_CHILDREN * ti;
    }

#pragma omp parallel for
    for (ti = 0; ti < P8EST_CHILDREN * num_trees; ti++) {
        tree_to_vertex[ti] = -1;
    }

    logx[0] = SC_LOG2_32 (m - 1) + 1;
    logx[1] = SC_LOG2_32 (n - 1) + 1;
    n_iter = (1 << logx[0]) * (1 << logx[1]);
    if (logx[0] <= logx[1]) {
        rankx[0] = 0;
        rankx[1] = 1;
    } else {
        rankx[0] = 1;
        rankx[1] = 0;
    }

    logx[2] = SC_LOG2_32 (p - 1) + 1;
    n_iter *= (1 << logx[2]);
    if (logx[2] < logx[rankx[0]]) {
        rankx[2] = rankx[1];
        rankx[1] = rankx[0];
        rankx[0] = 2;
    }
    else if (logx[rankx[1]] <= logx[2]) {
        rankx[2] = 2;
    }
    else {
        rankx[2] = rankx[1];
        rankx[1] = 2;
    }

    linear_to_tree = P4EST_ALLOC(p8est_topidx_t, n_iter);
    tree_to_corner2 = P4EST_ALLOC(p8est_topidx_t, num_trees);
    tree_to_edge2 = P4EST_ALLOC(p8est_topidx_t, 3 * num_trees);

    tj = 0;
    tk = 0;
    tl = 0;

    p8est_topidx_t  tx, ty;

    for (ti = 0; ti < n_iter; ti++) {
        brick_linear_to_xyz (ti, logx, rankx, coord);
        tx = coord[0];
        ty = coord[1];
        tz = coord[2];
        if (tx < m && ty < n && tz < p && 1) 
        {
            linear_to_tree[ti] = tj;
            if ((tx < m - 1 || periodic_a) && (ty < n - 1 || periodic_b) && (tz < p - 1 || periodic_c) && 1) 
            {
                tree_to_corner2[tj] = tk++;
                tree_to_edge2[3 * tj] = tl++;
                tree_to_edge2[3 * tj + 1] = tl++;
                tree_to_edge2[3 * tj + 2] = tl++;
            } 
            else 
            {
                tree_to_corner2[tj] = -1;
                if ((ty < n - 1 || periodic_b) && (tz < p - 1 || periodic_c)) {
                  tree_to_edge2[3 * tj] = tl++;
                }
                else {
                  tree_to_edge2[3 * tj] = -1;
                }
                if ((tx < m - 1 || periodic_a) && (tz < p - 1 || periodic_c)) {
                  tree_to_edge2[3 * tj + 1] = tl++;
                }
                else {
                  tree_to_edge2[3 * tj + 1] = -1;
                }
                if ((tx < m - 1 || periodic_a) && (ty < n - 1 || periodic_b)) {
                  tree_to_edge2[3 * tj + 2] = tl++;
                }
                else {
                  tree_to_edge2[3 * tj + 2] = -1;
                }
            }
            tj++;
        }
        else 
        {
            linear_to_tree[ti] = -1;
        }
    }
    P4EST_ASSERT(tj == num_trees);
    P4EST_ASSERT(tk == num_corners);
    P4EST_ASSERT(tl == num_edges);

    for (ti = 0; ti < n_iter; ti++) {
        brick_linear_to_xyz (ti, logx, rankx, coord);
        tx = coord[0];
        ty = coord[1];
        tz = coord[2];
        if (tx < m && ty < n && tz < p && 1) {
            tj = linear_to_tree[ti];
            P4EST_ASSERT(tj >= 0);
            for (i = 0; i < P8EST_DIM; i++) {
                for (j = 0; j < 2; j++) {
                    l = 2 * i + j;
                    coord2[0] = ((tx + ((i == 0) ? (2 * j - 1) : 0)) + m) % m;
                    coord2[1] = ((ty + ((i == 1) ? (2 * j - 1) : 0)) + n) % n;
                    coord2[2] = ((tz + ((i == 2) ? (2 * j - 1) : 0)) + p) % p;

                    tf[l] = brick_xyz_to_linear (coord2, logx, rankx);
                    P4EST_ASSERT(tf[l] < n_iter);
                    tf[l] = linear_to_tree[tf[l]];
                    P4EST_ASSERT(tf[l] >= 0);
                }
                for (j = 0; j < 4; j++) {
                    l = 4 * i + j;
                    coord2[0] = ((tx + ((i == 0) ? 0 : (2 * (j & 1) - 1))) + m) % m;
                    coord2[1] = ((ty + ((i == 1) ? 0 :
                                      (2 * ((i == 0) ? (j & 1) : (j / 2)) - 1))) +
                               n) % n;
                    coord2[2] = ((tz + ((i == 2) ? 0 : (2 * (j / 2) - 1))) + p) % p;
                    te[l] = brick_xyz_to_linear (coord2, logx, rankx);
                    P4EST_ASSERT(te[l] < n_iter);
                    te[l] = linear_to_tree[te[l]];
                    P4EST_ASSERT(te[l] >= 0);
                }
            }
            for (i = 0; i < P8EST_CHILDREN; i++) {
                coord2[0] = ((tx + (((i & 1) == 0) ? -1 : 1)) + m) % m;
                coord2[1] = ((ty + ((((i >> 1) & 1) == 0) ? -1 : 1)) + n) % n;
                coord2[2] = ((tz + (((i >> 2) == 0) ? -1 : 1)) + p) % p;
                tc[i] = brick_xyz_to_linear (coord2, logx, rankx);
                P4EST_ASSERT(tc[i] < n_iter);
                tc[i] = linear_to_tree[tc[i]];
                P4EST_ASSERT(tc[i] >= 0);
            }
            for (i = 0; i < P8EST_DIM; i++) {
                for (j = 0; j < 2; j++) {
                    l = i * 2 + j;
                    if (!periodic[i] &&
                          ((coord[i] == 0 && j == 0) || (coord[i] == max[i] && j == 1))) {
                        tree_to_tree[tj * P8EST_FACES + l] = tj;
                        tree_to_face[tj * P8EST_FACES + l] = (int8_t) l;
                    }
                    else 
                    {
                        tree_to_tree[tj * P8EST_FACES + l] = tf[l];
                        tree_to_face[tj * P8EST_FACES + l] = (int8_t) (i * 2 + (j ^ 1));
                    }
                }
                if (tree_to_edge != NULL) {
                    /** dir1, dir2 should be in correct z order */
                    dir1 = (i == 0) ? 1 : 0;
                    dir2 = (i == 2) ? 1 : 2;
                    for (j = 0; j < 4; j++) {
                        l = i * 4 + j;
                        if ((!periodic[dir1] &&
                             ((coord[dir1] == 0 && (j & 1) == 0) ||
                              (coord[dir1] == max[dir1] && (j & 1) == 1))) ||
                            (!periodic[dir2] &&
                                ((coord[dir2] == 0 && (j / 2) == 0) ||
                                (coord[dir2] == max[dir2] && (j / 2) == 1)))) 
                        {
                            tree_to_edge[tj * P8EST_EDGES + l] = -1;
                        }
                        else 
                        {
                            switch (j) {
                            case 0:
                              ttemp = tree_to_edge2[te[l] * 3 + i];
                              break;
                            case 1:
                              ttemp = tree_to_edge2[tf[dir2 * 2] * 3 + i];
                              break;
                            case 2:
                              ttemp = tree_to_edge2[tf[dir1 * 2] * 3 + i];
                              break;
                            case 3:
                              ttemp = tree_to_edge2[tj * 3 + i];
                              break;
                            default:
                              SC_ABORT_NOT_REACHED();
                        }
                        P4EST_ASSERT(ttemp >= 0);
                        tree_to_edge[tj * P8EST_EDGES + l] = ttemp;
                        edge_to_tree[4 * ttemp + (3 - j)] = tj;
                        edge_to_edge[4 * ttemp + (3 - j)] = (int8_t) l;
                        }
                    }
                }
            }
            for (i = 0; i < P8EST_CHILDREN; i++) {
                if (tree_to_corner != NULL) {
                    c[0] = i & 1;
                    c[1] = (i >> 1) & 1;
                    c[2] = i >> 2;
                    if ((!periodic[0] &&
                         ((coord[0] == 0 && c[0] == 0) ||
                          (coord[0] == max[0] && c[0] == 1))) ||
                        (!periodic[1] &&
                         ((coord[1] == 0 && c[1] == 0) ||
                          (coord[1] == max[1] && c[1] == 1))) ||
                        (!periodic[2] &&
                         ((coord[2] == 0 && c[2] == 0) ||
                          (coord[2] == max[2] && c[2] == 1))) ||
                        0) 
                    {
                        tree_to_corner[tj * P8EST_CHILDREN + i] = -1;
                    }
                    else 
                    {
                        switch (i) 
                        {
                            case 0:
                                ttemp = tc[0];
                                break;
                            case 1:
                                ttemp = te[0];
                                break;
                            case 2:
                                ttemp = te[4];
                                break;
                            case 3:
                                ttemp = tf[4];
                                break;
                            case 4:
                                ttemp = te[8];
                                break;
                            case 5:
                                ttemp = tf[2];
                                break;
                            case 6:
                                ttemp = tf[0];
                                break;
                            case 7:
                                ttemp = tj;
                                break;
                            default:
                                SC_ABORT_NOT_REACHED ();
                        }
                        ttemp = tree_to_corner2[ttemp];
                        P4EST_ASSERT(ttemp >= 0);
                        tree_to_corner[tj * P8EST_CHILDREN + i] = ttemp;
                        corner_to_tree[ttemp * P8EST_CHILDREN +
                                       (P8EST_CHILDREN - 1 - i)] = tj;
                        corner_to_corner[ttemp * P8EST_CHILDREN +
                                         (P8EST_CHILDREN - 1 - i)] = (int8_t) i;
                    }
                }
                if (tz > 0 && (i >> 2) == 0) {
                  tree_to_vertex[tj * P8EST_CHILDREN + i] =
                    tree_to_vertex[tf[4] * P8EST_CHILDREN + i + 4];
                }
                else
                {
                    if (ty > 0 && ((i >> 1) & 1) == 0) 
                    {
                      tree_to_vertex[tj * P8EST_CHILDREN + i] =
                        tree_to_vertex[tf[2] * P8EST_CHILDREN + i + 2];
                    }
                    else if (tx > 0 && (i & 1) == 0) 
                    {
                        tree_to_vertex[tj * P8EST_CHILDREN + i] =
                        tree_to_vertex[tf[0] * P8EST_CHILDREN + i + 1];
                    }
                    else 
                    {
                        tree_to_vertex[tj * P8EST_CHILDREN + i] = vcount++;
                        vertices[vicount++] = (double) x0 + dx * (tx + (i & 1));
                        vertices[vicount++] = (double) y0 + dy * (ty + ((i >> 1) & 1));
                        vertices[vicount++] = (double) z0 + dz * (tz + (i >> 2));
                        #ifdef OXLEY_PRINT_VERTICES
                        std::cout << "( " << vertices[vicount-3] << ", " 
                                          << vertices[vicount-2] << ", "
                                          << vertices[vicount-1] << " )" << std::endl;
                        #endif
                    }
                }
            }
        }
    }

    P4EST_ASSERT(vcount == num_vertices);
    // P4EST_FREE(linear_to_tree);
    // P4EST_FREE(tree_to_corner2);
    // P4EST_FREE(tree_to_edge2);

#ifdef OXLEY_ENABLE_TIMECONSUMING_DEBUG_CHECKS
    P4EST_ASSERT(p8est_connectivity_is_valid(conn)); //This is very time consuming
#endif

    return conn;
}

// instantiate our two supported versions
template
void Brick::assembleGradientImpl<real_t>(escript::Data& out,
                                         const escript::Data& in) const;

template
void Brick::assembleGradientImpl<cplx_t>(escript::Data& out,
                                         const escript::Data& in) const;

//protected
void Brick::assembleIntegrate(std::vector<real_t>& integrals, const escript::Data& arg) const
{
    assembleIntegrateImpl<real_t>(integrals, arg);
}

//protected
void Brick::assembleIntegrate(std::vector<cplx_t>& integrals, const escript::Data& arg) const
{
    assembleIntegrateImpl<cplx_t>(integrals, arg);
}

//private
template<typename Scalar>
void Brick::assembleIntegrateImpl(std::vector<Scalar>& integrals, const escript::Data& arg) const
{
    const dim_t numComp = arg.getDataPointSize();
    const int fs = arg.getFunctionSpace().getTypeCode();
    const Scalar zero = static_cast<Scalar>(0);

    const bool HavePointData = (fs == Points);

#ifdef ESYS_MPI
    if(HavePointData && escript::getMPIRankWorld() == 0) {
#else
    if(HavePointData) {
#endif
        integrals[0] += arg.getNumberOfTaggedValues();

    } else if (fs == Elements && arg.actsExpanded()) {
        // 2x2x2 Gauss on the trilinear hex: each of the 8 points carries V/8.
        // Element size is m_NX/2^level (per-block length / 2^level); forestData
        // m_dx is total-length based for Brick and must not be used here.
        std::vector<Scalar> int_local(numComp, zero);
        long id = 0;
        for(p8est_topidx_t treeid = p8est->first_local_tree; treeid <= p8est->last_local_tree; ++treeid)
        {
            p8est_tree_t * tree = p8est_tree_array_index(p8est->trees, treeid);
            sc_array_t * tquadrants = &tree->quadrants;
            p8est_locidx_t Q = (p8est_locidx_t) tquadrants->elem_count;
            for(int q = 0; q < Q; ++q, ++id)
            {
                p8est_quadrant_t * quad = p8est_quadrant_array_index(tquadrants, q);
                const double h = (double)(1 << quad->level);
                const real_t w = (m_NX[0]/h) * (m_NX[1]/h) * (m_NX[2]/h) / 8.;
                const Scalar* f = arg.getSampleDataRO(id, zero);
                for (index_t i = 0; i < numComp; ++i) {
                    Scalar s = zero;
                    for (int c = 0; c < 8; ++c)
                        s += f[INDEX2(i,c,numComp)];
                    int_local[i] += s*w;
                }
            }
        }
        for (index_t i = 0; i < numComp; ++i)
            integrals[i] += int_local[i];

    } else if (fs==ReducedElements || (fs==Elements && !arg.actsExpanded())) {
        // single centre point carries the full element volume V.
        std::vector<Scalar> int_local(numComp, zero);
        long id = 0;
        for(p8est_topidx_t treeid = p8est->first_local_tree; treeid <= p8est->last_local_tree; ++treeid)
        {
            p8est_tree_t * tree = p8est_tree_array_index(p8est->trees, treeid);
            sc_array_t * tquadrants = &tree->quadrants;
            p8est_locidx_t Q = (p8est_locidx_t) tquadrants->elem_count;
            for(int q = 0; q < Q; ++q, ++id)
            {
                p8est_quadrant_t * quad = p8est_quadrant_array_index(tquadrants, q);
                const double h = (double)(1 << quad->level);
                const real_t w = (m_NX[0]/h) * (m_NX[1]/h) * (m_NX[2]/h);
                const Scalar* f = arg.getSampleDataRO(id, zero);
                for (index_t i = 0; i < numComp; ++i)
                    int_local[i] += f[i]*w;
            }
        }
        for (index_t i = 0; i < numComp; ++i)
            integrals[i] += int_local[i];

    } else if (fs == FaceElements && arg.actsExpanded()) {
        // Each boundary face element carries area A over 4 Gauss points (A/4 each).
        // Face area uses the in-plane element sizes dx=m_NX/2^level.
        std::vector<Scalar> int_local(numComp, zero);
        const std::vector<borderNodeInfo>* lists[6] = {
            &NodeIDsLeft, &NodeIDsRight, &NodeIDsBottom, &NodeIDsTop,
            &NodeIDsAbove, &NodeIDsBelow };
        static const int planeAxes[6][2] = {{1,2},{1,2},{0,2},{0,2},{0,1},{0,1}};
        for(int fc = 0; fc < 6; ++fc) {
            if(m_faceOffset[fc] < 0) continue;
            const std::vector<borderNodeInfo>& L = *lists[fc];
            for(size_t k = 0; k < L.size(); ++k) {
                const double h = (double)(1 << L[k].level);
                const real_t A = (m_NX[planeAxes[fc][0]]/h) * (m_NX[planeAxes[fc][1]]/h);
                const Scalar* f = arg.getSampleDataRO(m_faceOffset[fc]+k, zero);
                for (index_t i = 0; i < numComp; ++i) {
                    Scalar srow = zero;
                    for (int c = 0; c < 4; ++c)
                        srow += f[INDEX2(i,c,numComp)];
                    int_local[i] += srow * (A/4.);
                }
            }
        }
        for (index_t i = 0; i < numComp; ++i)
            integrals[i] += int_local[i];

    } else if (fs==ReducedFaceElements || (fs==FaceElements && !arg.actsExpanded())) {
        std::vector<Scalar> int_local(numComp, zero);
        const std::vector<borderNodeInfo>* lists[6] = {
            &NodeIDsLeft, &NodeIDsRight, &NodeIDsBottom, &NodeIDsTop,
            &NodeIDsAbove, &NodeIDsBelow };
        static const int planeAxes[6][2] = {{1,2},{1,2},{0,2},{0,2},{0,1},{0,1}};
        for(int fc = 0; fc < 6; ++fc) {
            if(m_faceOffset[fc] < 0) continue;
            const std::vector<borderNodeInfo>& L = *lists[fc];
            for(size_t k = 0; k < L.size(); ++k) {
                const double h = (double)(1 << L[k].level);
                const real_t A = (m_NX[planeAxes[fc][0]]/h) * (m_NX[planeAxes[fc][1]]/h);
                const Scalar* f = arg.getSampleDataRO(m_faceOffset[fc]+k, zero);
                for (index_t i = 0; i < numComp; ++i)
                    int_local[i] += f[i] * A;
            }
        }
        for (index_t i = 0; i < numComp; ++i)
            integrals[i] += int_local[i];

    } else {
        throw OxleyException("assembleIntegrate: unsupported function space");
    }
}
RankVector Brick::getOwnerVector(int fsType) const
{
    RankVector owner;
    const int rank = m_mpiInfo->rank;

    // Serial-correct: every local element is owned by this rank (no ghosts).
    // TODO(MPI, A6): mark ghost elements across the p4est partition boundary.
    if (fsType == Elements || fsType == ReducedElements) {
        owner.assign(getNumElements(), rank);
    } else if (fsType == FaceElements || fsType == ReducedFaceElements) {
        owner.assign(getNumFaceElements(), rank);
    } else {
        throw ValueError("getOwnerVector: only valid for element types");
    }
    return owner;
}

const long Brick::getNodeId(double x, double y, double z)
{
    const double coords[3] = {x, y, z};
    return (long) findNode(coords);
}

// void Brick::getNeighouringNodeIDs(int8_t level, p8est_qcoord_t x, p8est_qcoord_t y, p8est_qcoord_t z, 
//                                                 p8est_topidx_t treeid, long (&ids) [8]) const
// {
//     p8est_qcoord_t l = P8EST_QUADRANT_LEN(level);
//     int adj[8][3] = {{0,0,0},{l,0,0},{0,l,0},{l,l,0},
//                      {0,0,l},{l,0,l},{0,l,l},{l,l,l}};
// #pragma omp parallel for
//     for(int i=0; i<8;i++)
//     {
//         double xy[3];
//         p8est_qcoord_to_vertex(p8est->connectivity, treeid, x+adj[i][0], y+adj[i][1], z+adj[i][1], xy);
//         ids[i]=(long) NodeIDs.find(std::make_tuple(xy[0],xy[1],xy[2]))->second;
//     }
// }

/**
    \brief
    Applies a refinementzone
*/
escript::Domain_ptr Brick::applyRefinement(RefinementQueue& R)
{
    oxleytimer.toc("Applying the refinement zone...");

    bool update=false; // Update after the refinement zones have been applied to save time

    oxleytimer.toc("\tcreating a new Brick...");
    oxley::Brick * newDomain = new Brick(*this, m_order, update);
    oxleytimer.toc("\t\t\t Brick created...");

#ifdef OXLEY_ENABLE_DEBUG_CHECKS 
    std::cout << "In applyRefinement debug checks..." << std::endl;
    std::cout << "Checking connectivity (1) ... ";
    if(!p8est_connectivity_is_valid(connectivity))
        std::cout << "broken" << std::endl;
    else
        std::cout << "OK" << std::endl;
    #ifdef OXLEY_ENABLE_TIMECONSUMING_DEBUG_CHECKS
        std::cout << "Checking connectivity (2) ... ";
        if(!p8est_connectivity_is_equal(connectivity, newDomain->connectivity))
            std::cout << "broken" << std::endl;
        else
            std::cout << "OK" << std::endl;
        std::cout << "Checking connectivity (3) ... ";
        if(!p8est_connectivity_is_equivalent(connectivity, newDomain->connectivity))
            std::cout << "broken" << std::endl;
        else
            std::cout << "OK" << std::endl;
    #endif
    std::cout << "Checking p8est (1) ... ";
    if(!p8est_is_valid(p8est))
        std::cout << "broken" << std::endl;
    else
        std::cout << "OK" << std::endl;
    std::cout << "Checking p8est (2) ... ";
    if(!p8est_is_equal(p8est, newDomain->p8est, false))
        std::cout << "broken" << std::endl;
    else
        std::cout << "OK" << std::endl;
#endif

    newDomain->AutomaticMeshUpdateOnOff(false);

    int numberOfRefinements = R.getNumberOfOperations();

    oxleytimer.toc("\tapplying individual refinements");
    for(int n = 0; n < numberOfRefinements; n++)
    {
        #ifdef OXLEY_ENABLE_PROFILE_TIMERS_INFORMATIONAL
        if(n % 200 == 0)
            oxleytimer.toc("apply_refinementQueue: Updating Mesh (" + std::to_string(n) + " of " + std::to_string(numberOfRefinements) + ")");
        #endif

        RefinementType Refinement = R.getRefinement(n);
        newDomain->setRefinementLevels(Refinement.levels);
        switch(Refinement.flavour)
        {
            case UNIFORM:
            {
                newDomain->refineMesh("uniform");
                break;
            }
            case POINT3D:
            {
                double x=Refinement.x0;
                double y=Refinement.y0;
                double z=Refinement.z0;
                oxleytimer.toc("\trefining point ("+std::to_string(x)+","+std::to_string(y)+","+std::to_string(z)+")");
                newDomain->refinePoint(x,y,z);
                break;
            }
            case REGION3D:
            {
                double x0=Refinement.x0;
                double y0=Refinement.y0;
                double z0=Refinement.z0;
                double x1=Refinement.x1;
                double y1=Refinement.y1;
                double z1=Refinement.z1;
                newDomain->refineRegion(x0, x1, y0, y1, z0, z1);
                break;
            }
            case SPHERE:
            {
                double x0=Refinement.x0;
                double y0=Refinement.y0;
                double z0=Refinement.z0;
                double r0=Refinement.r;
                newDomain->refineSphere(x0,y0,z0,r0);
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
                    {
                        newDomain->refineBoundary("TOP",dx);
                        break;
                    }
                    case BOTTOM:
                    {
                        newDomain->refineBoundary("BOTTOM",dx);
                        break;
                    }
                    default:
                    {
                        throw OxleyException("Invalid border direction.");
                    }
                }
                break;
            }
            case MASK3D:
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
            case MASK2D:
            case CIRCLE:
            case POINT2D:
            case REGION2D:
            default:
                throw OxleyException("Unknown refinement algorithm.");
        }

        newDomain->updateMeshBackend();
    }

    oxleytimer.toc("\tupdating the mesh");
    newDomain->updateMesh();
    newDomain->AutomaticMeshUpdateOnOff(true);

    oxleytimer.toc("done");

    return escript::Domain_ptr(newDomain);
}

int Brick::p8est_connectivity_is_valid_fast(p8est_connectivity_t * conn)
{
    // int                 nvert;
    // int                 face, rface, nface, orientation;
    // int                 errcode, errcount;
    // int                 edge, nedge;
    // int                 flip, nflip, nflip1, nflip2;
    // p8est_topidx_t      nett, aedge, edge_begin, edge_end;
    // int                 corner, ncorner;
    // int                 good, cfound;
    // p8est_topidx_t      vertex, tree, ntree;
    // p8est_topidx_t      acorner, corner_begin, corner_end;
    // p8est_topidx_t      nctt;
    const p8est_topidx_t num_vertices = conn->num_vertices;
    const p8est_topidx_t num_trees = conn->num_trees;
    const p8est_topidx_t *ttv = conn->tree_to_vertex;
    const p8est_topidx_t *ttt = conn->tree_to_tree;
    const int8_t         *ttf = conn->tree_to_face;
    const p8est_topidx_t num_edges = conn->num_edges;
    const p8est_topidx_t *tte = conn->tree_to_edge;
    const p8est_topidx_t *eoff = conn->ett_offset;
    const p8est_topidx_t *ett = conn->edge_to_tree;
    const int8_t         *ete = conn->edge_to_edge;
    const p8est_topidx_t num_ett = eoff[num_edges];
    p8est_edge_info_t ei;
    const p8est_topidx_t num_corners = conn->num_corners;
    const p8est_topidx_t *ttc = conn->tree_to_corner;
    const p8est_topidx_t *coff = conn->ctt_offset;
    const p8est_topidx_t *ctt = conn->corner_to_tree;
    const int8_t *ctc = conn->corner_to_corner;
    const p4est_topidx_t num_ctt = coff[num_corners];
    // p8est_corner_info_t ci;
    // sc_array_t *cta = &ci.corner_transforms;

    std::string message = "";

#ifdef OXLEY_ENABLE_PROFILE_TIMERS
    oxleytimer.toc("checking connectivity ...");
#endif

    // sc_array_t *eta = &ei.edge_transforms;
    // sc_array_init(eta, sizeof(p8est_edge_transform_t));
    
    // sc_array_init(cta, sizeof (p4est_corner_transform_t));

    if(num_vertices == 0 && (conn->vertices != NULL || ttv != NULL)) {
        throw OxleyException("Zero vertices still with arrays\n");
    }
    if(num_vertices > 0  && (conn->vertices == NULL || ttv == NULL)) {
        throw OxleyException("Nonzero vertices missing arrays\n");
    }

#ifdef OXLEY_ENABLE_PROFILE_TIMERS
    oxleytimer.toc("\t checking edges ...");
#endif

#pragma omp parallel for
    for(p4est_topidx_t nett = 0; nett < num_ett; ++nett) {
        if(ett[nett] < 0 || ett[nett] >= num_trees) {
            message += "Edge to tree " + std::to_string(nett) + " out of range\n";
        }
        if(ete[nett] < 0 || ete[nett] >= 24) {
            message += "Edge to edge " + std::to_string(nett) + " out of range\n";
        }
    }

#ifdef OXLEY_ENABLE_PROFILE_TIMERS
    oxleytimer.toc("\t checking corners ...");
#endif
    #pragma omp parallel for
    for(p8est_topidx_t nctt = 0; nctt < num_ctt; ++nctt) {
        if(ctt[nctt] < 0 || ctt[nctt] >= num_trees) {
            std::string message = "Corner to tree " + std::to_string(nctt) + " out of range\n";
            throw OxleyException(message);
        }
        if(ctc[nctt] < 0 || ctc[nctt] >= 8) {
            std::string message = "Corner to corner " + std::to_string(nctt) + " out of range\n";
            throw OxleyException(message);
        }
    }

#ifdef OXLEY_ENABLE_PROFILE_TIMERS
    oxleytimer.toc("\t checking vertices ...");
#endif

    if(num_vertices > 0) {
#pragma omp parallel for
        for(p8est_topidx_t tree = 0; tree < num_trees; ++tree) {
            for(int nvert = 0; nvert < 8; ++nvert) {
                p8est_topidx_t vertex = ttv[tree * 8 + nvert];
                if(vertex < 0 || vertex >= num_vertices) {
                    std::string message = "Tree to vertex out of range " + std::to_string(tree) + " " + std::to_string(nvert);
                    throw OxleyException(message);
                }
            }
        }
    }

    if((conn->tree_to_attr != NULL) != (conn->tree_attr_bytes > 0)) {
        message += "Tree attribute properties inconsistent " + std::to_string(conn->tree_attr_bytes);
    }

#ifdef OXLEY_ENABLE_PROFILE_TIMERS
    oxleytimer.toc("\t checking trees ...");
#endif

#pragma omp parallel for
    for(p4est_topidx_t tree = 0; tree < num_trees; ++tree) {
        for(int face = 0; face < 6; ++face) {
            p8est_topidx_t ntree = ttt[tree * 6 + face];
            if(ntree < 0 || ntree >= num_trees) {
                message += "Tree to tree out of range " + std::to_string(tree) + " " + std::to_string(face) + "\n";
            }
            int rface = (int) ttf[tree * 6 + face];
            if (rface < 0 || rface >= 6 * 4) {
                message += "Tree to face out of range " + std::to_string(tree) + " " + std::to_string(face) + "\n";
            }
            int nface = rface % 6;      /* clamp to a real face index */
            int orientation = rface / 6;        /* 0..P4EST_HALF-1 */
            if(ntree == tree) {
            /* no neighbor across this face or self-periodic */
                if (nface == face && orientation != 0) {
                    message += "Face invalid in " + std::to_string(tree) + " " + std::to_string(face) + "\n";
                }
            }
            if(ntree != tree || nface != face) {
                /* check reciprocity */
                if (ttt[ntree * 6 + nface] != tree) {
                    message += "Tree to tree reciprocity in " + std::to_string(tree) + " " + std::to_string(face) + "\n";
                }
                if ((int) ttf[ntree * 6 + nface] !=
                    face + 6 * orientation) {
                    message += "Tree to face reciprocity in " + std::to_string(tree) + " " + std::to_string(face) + "\n";
                }
            }
        }

        for(p8est_topidx_t aedge = 0; aedge < num_edges; ++aedge) {
            if (eoff[aedge + 1] < eoff[aedge]) {
                message += "Edge offset backwards " + std::to_string(aedge) + "\n";
            }
        }

        if(num_edges > 0) {
            for(int edge = 0; edge < P8EST_EDGES; ++edge) {
                p8est_find_edge_transform(conn, tree, edge, &ei);
                p8est_topidx_t aedge = tte[tree * P8EST_EDGES + edge];
                if(aedge < -1 || aedge >= num_edges) {
                    message += "Tree to edge out of range "  + std::to_string(tree) + " " + std::to_string(edge) + "\n";
                }
                if(aedge == -1) {
                    continue;
                }
                int errcode = 0, errcount = 0;
                int flip = -1, nflip1 = -1, nflip2 = -1;
                p8est_topidx_t edge_begin = eoff[aedge];
                p8est_topidx_t edge_end = eoff[aedge + 1];
                if(edge_begin < 0 || edge_begin >= num_ett || edge_end < 0 || edge_end > num_ett) {
                    message += "Invalid edge range "  + std::to_string(tree) + " " + std::to_string(edge) + "\n";
                }
                for(p8est_topidx_t nett = edge_begin; nett < edge_end; ++nett) {
                    p8est_topidx_t ntree = ett[nett];
                    int nedge = (int) ete[nett] % P8EST_EDGES;
                    if (tte[ntree * P8EST_EDGES + nedge] != aedge) {
                        message += "Invalid edge range "  + std::to_string(tree) + " " + std::to_string(edge)  + " " + std::to_string(nett) + "\n";
                    }
                    int nflip = (int) ete[nett] / P8EST_EDGES;
                    if(ntree == tree && nedge == edge) {
                        if(flip != -1 && nflip == flip) {
                            errcode = 1;
                            break;
                        }
                        flip = nflip;
                        continue;
                    }
                    ++errcount;
                }
                if(errcode > 0) {
                    // std::string message = "Shared edge " + std::to_string(tree) + " " + std::to_string((int) edge) + " " + std::to_string((int) nett) + " " + " inconsistency " +  + std::to_string((int) errcode);
                    // throw OxleyException(message);
                    // throw OxleyException("Shared edge inconsistency");
                    message += "Shared edge inconsistency\n";
                }
                if(flip == -1 || !(nflip1 == -1 || nflip1 == flip) || !(nflip2 == -1 || nflip2 == flip)) {
                    message += "Shared edge " + std::to_string(tree) + " " + std::to_string(edge) + " " + " inconsistent flip " + "\n";
                }
            }
        }

        for(p8est_topidx_t acorner = 0; acorner < num_corners; ++acorner) {
            if (coff[acorner + 1] < coff[acorner]) {
                message += "Corner offset backwards " + std::to_string(acorner) + "\n";
            }
        }

        if(num_corners > 0) {
            for(int corner = 0; corner < 8; ++corner) {
                p8est_corner_info_t ci;
                sc_array_t *cta = &ci.corner_transforms;
                sc_array_init(cta, sizeof(p8est_corner_transform_t));
                p8est_find_corner_transform(conn, tree, corner, &ci);  //here
                p8est_topidx_t acorner = ttc[tree * 8 + corner];
                if (acorner < -1 || acorner >= num_corners) {
                    message += "Tree to corner out of range " + std::to_string(tree) + " " + std::to_string(corner) + "\n";
                }
                if (acorner == -1) {
                    continue;
                }
                int errcode = 0, errcount = 0;
                int cfound = 0;
                p8est_topidx_t corner_begin = coff[acorner];
                p8est_topidx_t corner_end = coff[acorner + 1];
                if(corner_begin < 0 || corner_begin >= num_ctt || corner_end < 0 || corner_end > num_ctt) {
                    message += "Invalid corner range " + std::to_string(tree) + " " + std::to_string(corner) + "\n";
                }
                for(p8est_topidx_t nctt = corner_begin; nctt < corner_end; ++nctt) {
                p8est_topidx_t ntree = ctt[nctt];
                int ncorner = (int) ctc[nctt];
                if(ttc[ntree * P8EST_CHILDREN + ncorner] != acorner) {
                    message += "Invalid corner range " + std::to_string(tree) + " " + std::to_string(corner) + " " + std::to_string(nctt) + "\n";
                }
                if (ntree == tree && ncorner == corner) {
                    if (cfound) {
                        errcode = 1;
                        break;
                    }
                    cfound = 1;
                    continue;
                }
                ++errcount;
                }
                if (errcode > 0) {
                    message += "Invalid corner range " + std::to_string(tree) + " " + std::to_string(corner) + " inconsistency "+ std::to_string(errcode) + "\n";
                }
                if (!cfound) {
                    message += "Shared corner " + std::to_string(tree) + " " + std::to_string(corner) + " inconsistent count B \n";
                }
                sc_array_reset(cta);
            }
        }
    }

    if(message.compare("")!=0)
        throw OxleyException(message);
  
#ifdef OXLEY_ENABLE_PROFILE_TIMERS
    oxleytimer.toc("done ...");
#endif

    // sc_array_reset (eta);
    // sc_array_reset (cta);

    return 1;
}


// Defined here rather than in RefinementQueue.cpp: it needs the concrete
// domain class, and including this header there runs into the OxleyData.h <->
// Brick.h include cycle.
escript::Domain_ptr RefinementQueue3D::apply(escript::Domain_ptr domain)
{
    if(domain.get() == NULL)
        throw OxleyException("RefinementQueue3D::apply: no domain given.");
    Brick * d = dynamic_cast<Brick *>(domain.get());
    if(d == NULL)
        throw OxleyException("RefinementQueue3D::apply: the domain is not a Brick. "
                "Use RefinementQueue2D for a Rectangle.");
    return d->applyRefinement(*this);
}

} // end of namespace oxley
