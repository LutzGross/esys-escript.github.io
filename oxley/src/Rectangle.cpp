/*****************************************************************************
*
* Copyright (c) 2003-2019 by The University of Queensland
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

#include <algorithm>
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

#include <oxley/AbstractAssembler.h>
#include <oxley/DefaultAssembler2D.h>
#include <oxley/domainhelpers.h>
#include <oxley/InitAlgorithms.h>
#include <oxley/Oxley.h>
#include <oxley/OxleyData.h>
#include <oxley/Rectangle.h>
#include <oxley/RefinementAlgorithms.h>

#include <p4est.h>
#include <p4est_algorithms.h>
#include <p4est_bits.h>
#include <p4est_connectivity.h>
#include <p4est_extended.h>
#include <p4est_iterate.h>
#include <p4est_lnodes.h>
#include <p4est_vtk.h>

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
       Constructor
    */
Rectangle::Rectangle(int order,
    dim_t n0, dim_t n1,
    double x0, double y0,
    double x1, double y1,
    int d0, int d1,
    const std::vector<double>& points, 
    const std::vector<int>& tags,
    const TagMap& tagnamestonums,
    int periodic0, int periodic1):
    OxleyDomain(2, order){

    // Possible error: User passes invalid values for the dimensions
    if(n0 <= 0 || n1 <= 0)
        throw OxleyException("Number of elements in each spatial dimension must be positive");

    // Ignore d0 and d1 if we are running in serial
    m_mpiInfo = escript::makeInfo(MPI_COMM_WORLD);
    if(m_mpiInfo->size == 1) {
        d0=1;
        d1=1;
    }

    // If the user did not set the number of divisions manually
    if(d0 == -1 && d1 == -1)
    {
        d0 = m_mpiInfo->size < 3 ? 1 : m_mpiInfo->size / 3;
        d1 = m_mpiInfo->size / d0;

        if(d0*d1 != m_mpiInfo->size)
            throw OxleyException("Could not find values for d0, d1 and d2. Please set them manually.");
    }

    // Create the connectivity
    // const p4est_topidx_t num_vertices = (n0+1)*(n1+1);
    // const p4est_topidx_t num_trees = n0*n1;
    // const p4est_topidx_t num_corners = (n0-1)*(n1-1);
    // const double vertices[P4EST_CHILDREN * 3] = {
    //                                             x0, y0, 0,
    //                                             x1, y0, 0,
    //                                             x0, y1, 0,
    //                                             x1, y1, 0,
    //                                             };
    // const p4est_topidx_t tree_to_vertex[P4EST_CHILDREN] = {0, 1, 2, 3,};
    // const p4est_topidx_t tree_to_tree[P4EST_FACES] = {0, 0, 0, 0,};
    // // const int8_t tree_to_face[P4EST_FACES] = {1, 0, 3, 2,}; //TODO: add in periodic boundary conditions
    // const int8_t tree_to_face[P4EST_FACES] = {0, 1, 2, 3,};
    // const p4est_topidx_t tree_to_corner[P4EST_CHILDREN] = {0, 0, 0, 0,};
    // const p4est_topidx_t ctt_offset[2] = {0, P4EST_CHILDREN,};
    // const p4est_topidx_t corner_to_tree[P4EST_CHILDREN] = {0, 0, 0, 0,};
    // const int8_t corner_to_corner[P4EST_CHILDREN] = {0, 1, 2, 3,};

    // connectivity = p4est_connectivity_new_copy(num_vertices, num_trees,
    //                                   num_corners, vertices, tree_to_vertex,
    //                                   tree_to_tree, tree_to_face,
    //                                   tree_to_corner, ctt_offset,
    //                                   corner_to_tree, corner_to_corner);

    // connectivity = p4est_connectivity_new_brick(n0, n1, false, false); 

    // signed int refinex = n0, refiney = n1, num_refine = 0;
    // long div=1;
    // while(refinex % 2 == 0 && refiney % 2 == 0)
    // {
    //     refinex /= 2;
    //     refiney /= 2;
    //     num_refine++;
    //     div*=2;
    // }
    connectivity = new_rectangle_connectivity(n0, n1, false, false, x0, y0, x1, y1);

#ifdef OXLEY_ENABLE_DEBUG_CHECKS //These checks are turned off by default as they can be very timeconsuming
    std::cout << "In Rectangle() constructor..." << std::endl;
    std::cout << "Checking connectivity ... ";
    if(!p4est_connectivity_is_valid(connectivity))
        std::cout << "broken" << std::endl;
    else
        std::cout << "OK" << std::endl;
#endif

    // Create a p4est
    p4est_locidx_t min_quadrants = n0*n1;
    int min_level = 0;
    int fill_uniform = 1;
    p4est = p4est_new_ext(m_mpiInfo->comm, connectivity, min_quadrants,
            min_level, fill_uniform, sizeof(quadrantData), init_rectangle_data, (void *) &forestData);

#ifdef OXLEY_ENABLE_DEBUG_CHECKS //These checks are turned off by default as they can be very timeconsuming
    std::cout << "Checking p4est ... ";
    if(!p4est_is_valid(p4est))
        std::cout << "broken" << std::endl;
    else
        std::cout << "OK" << std::endl;
#endif

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

    // Whether or not we have periodic boundaries
    forestData.periodic[0] = periodic0;
    forestData.periodic[1] = periodic1;

    // Find the grid spacing for each level of refinement in the mesh
#pragma omp parallel for
    for(int i = 0; i<=P4EST_MAXLEVEL; i++){
        double numberOfSubDivisions = (p4est_qcoord_t) (1 << (P4EST_MAXLEVEL - i));
        forestData.m_dx[0][i] = forestData.m_length[0] / (numberOfSubDivisions);
        forestData.m_dx[1][i] = forestData.m_length[1] / (numberOfSubDivisions);
    }

    // max levels of refinement
    forestData.max_levels_refinement = MAXREFINEMENTLEVELS;
    
    // element order
    m_order = order;

    // initial tag
    // tags[0] = 0;
    // numberOfTags=1;

    // Number of dimensions
    m_numDim=2;

    // Distribute the p4est across the processors
    int allow_coarsening = 0;
    p4est_partition(p4est, allow_coarsening, NULL);

    // Number the nodes
    updateNodeIncrements();
    renumberNodes();
    updateRowsColumns();
    updateNodeDistribution();
    updateElementIds();
    updateFaceOffset();
    updateFaceElementCount();
    // populateDofMap()

    // Tags
    populateSampleIds();
    for (TagMap::const_iterator i = tagnamestonums.begin(); i != tagnamestonums.end(); i++) {
        setTagMap(i->first, i->second);
    }
    
    // Dirac points and tags
    addPoints(points, tags);

    // srand(time(NULL));

    // To prevent segmentation faults when using numpy ndarray
#ifdef ESYS_HAVE_BOOST_NUMPY
    Py_Initialize();
    boost::python::numpy::initialize();
#endif

#ifdef ESYS_HAVE_PASO

    /// local array length shared
    // dim_t local_length = p4est->local_num_quadrants;
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

    // Refine the mesh 
    // refineMesh(num_refine, "uniform");
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
    p4est_connectivity_destroy(connectivity);
    p4est_destroy(p4est);
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

void Rectangle::interpolateNodesToNodesFiner(const escript::Data& source, escript::Data& target, const Rectangle& other) const
{

}

void Rectangle::interpolateNodesToElementsFiner(const escript::Data& source, escript::Data& target, const Rectangle& other)  const
{

}

void Rectangle::interpolateElementsToElementsCoarser(const escript::Data& source, escript::Data& target, const Rectangle& other)  const
{

}

void Rectangle::interpolateElementsToElementsFiner(const escript::Data& source, escript::Data& target, const Rectangle& other)  const
{

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
                for (index_t k=0; k<NodeIDsLeft.size()-1; k++) {
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
                for (index_t k=0; k<NodeIDsBottom.size()-1; k++) {
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
                for (index_t k=0; k<NodeIDsRight.size()-1; k++) {
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
                for (index_t k=0; k<NodeIDsTop.size()-1; k++) {
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
                for (index_t k=0; k<NodeIDsLeft.size()-1; k++) {
                    double* o = out.getSampleDataRW(m_faceOffset[0]+k);
                    *o++ = -1.;
                    *o = 0.;
                }
            }

            if (m_faceOffset[1] > -1) {
#pragma omp for nowait
                for (index_t k=0; k<NodeIDsBottom.size()-1; k++) {
                    double* o = out.getSampleDataRW(m_faceOffset[1]+k);
                    *o++ = 1.;
                    *o = 0.;
                }
            }

            if (m_faceOffset[2] > -1) {
#pragma omp for nowait
                for (index_t k=0; k<NodeIDsRight.size()-1; k++) {
                    double* o = out.getSampleDataRW(m_faceOffset[2]+k);
                    *o++ = 0.;
                    *o = -1.;
                }
            }

            if (m_faceOffset[3] > -1) {
#pragma omp for nowait
                for (index_t k=0; k<NodeIDsTop.size()-1; k++) {
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
}

void Rectangle::setToSize(escript::Data& out) const
{
    if (out.getFunctionSpace().getTypeCode() == Elements
        || out.getFunctionSpace().getTypeCode() == ReducedElements)
    {
        out.requireWrite();
        const dim_t numQuad = out.getNumDataPointsPerSample();
        // const dim_t numElements = getNumElements();

        for(p4est_topidx_t t = p4est->first_local_tree; t <= p4est->last_local_tree; t++) 
        {
            p4est_tree_t * currenttree = p4est_tree_array_index(p4est->trees, t);
            sc_array_t * tquadrants = &currenttree->quadrants;
            p4est_locidx_t Q = (p4est_locidx_t) tquadrants->elem_count;
        #pragma omp parallel for
            for(p4est_locidx_t e = nodes->global_offset; e < Q+nodes->global_offset; e++)
            {
                p4est_quadrant_t * quad = p4est_quadrant_array_index(tquadrants, e);
                // int l = quad->level;
                const double size = sqrt(forestData.m_dx[0][quad->level]*forestData.m_dx[0][quad->level]
                                        +forestData.m_dx[1][quad->level]*forestData.m_dx[1][quad->level]);
                double* o = out.getSampleDataRW(e);
                std::fill(o, o+numQuad, size);
            }
        }
    } 
    else if (out.getFunctionSpace().getTypeCode() == FaceElements
            || out.getFunctionSpace().getTypeCode() == ReducedFaceElements) 
    {
        out.requireWrite();
        const dim_t numQuad=out.getNumDataPointsPerSample();

        for(p4est_topidx_t t = p4est->first_local_tree; t <= p4est->last_local_tree; t++) 
        {
            p4est_tree_t * currenttree = p4est_tree_array_index(p4est->trees, t);
            sc_array_t * tquadrants = &currenttree->quadrants;
            p4est_locidx_t Q = (p4est_locidx_t) tquadrants->elem_count;
        #pragma omp parallel for
            for(p4est_locidx_t e = nodes->global_offset; e < Q+nodes->global_offset; e++)
            {
                // Work out what level this element is on 
                p4est_quadrant_t * quad = p4est_quadrant_array_index(tquadrants, e);
                quadrantData * quaddata = (quadrantData *) quad->p.user_data;

                if (quaddata->m_faceOffset[0]) {
                    double* o = out.getSampleDataRW(e);
                    std::fill(o, o+numQuad, forestData.m_dx[1][quad->level]);
                }
                else if (quaddata->m_faceOffset[1]) {
                    double* o = out.getSampleDataRW(e);
                    std::fill(o, o+numQuad, forestData.m_dx[1][quad->level]);
                }
                else if (quaddata->m_faceOffset[2]) {
                    double* o = out.getSampleDataRW(e);
                    std::fill(o, o+numQuad, forestData.m_dx[0][quad->level]);
                }
                else if (quaddata->m_faceOffset[3]) {
                    double* o = out.getSampleDataRW(e);
                    std::fill(o, o+numQuad, forestData.m_dx[0][quad->level]);
                }
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
            // check ownership of element's bottom left node
            // return (m_dofMap[id%m_NE[0]+m_NN[0]*(id/m_NE[0])] < getNumDOF());
            throw OxleyException("Rectangle::ownSample Currently not implemented.");
            return false;
        case FaceElements:
        case ReducedFaceElements:
            {
                // // determine which face the sample belongs to before
                // // checking ownership of corresponding element's first node
                // dim_t n=0;
                // for (size_t i=0; i<4; i++) {
                //     n+=m_faceCount[i];
                //     if (id<n) {
                //         index_t k;
                //         if (i==1)
                //             k=m_NN[0]-2;
                //         else if (i==3)
                //             k=m_NN[0]*(m_NN[1]-2);
                //         else
                //             k=0;
                //         // determine whether to move right or up
                //         const index_t delta=(i/2==0 ? m_NN[0] : 1);
                //         return (m_dofMap[k+(id-n+m_faceCount[i])*delta] < getNumDOF());
                //     }
                // }
                throw OxleyException("Rectangle::ownSample Currently not implemented.");
                return false;
            }
        default:
            break;
    }

    std::stringstream msg;
    msg << "ownSample: invalid function space type " << fsType;
    throw ValueError(msg.str());
}


dim_t Rectangle::getNumDataPointsGlobal() const
{
    return getNumNodes();
}

void Rectangle::dump(const std::string& fileName) const
{
#ifdef ESYS_HAVE_SILO
    
    // Add the suffix to the filename if required
    std::string fn = fileName;
    if (fileName.length() < 6 || fileName.compare(fileName.length()-5, 5, ".silo") != 0) {
        fn+=".silo";
    }

    int driver=DB_HDF5;

    // Silo file pointer
    DBfile* dbfile = NULL; 

    // The coordinate arrays
    float *pNodex = nullptr;
    float *pNodey = nullptr;
    long int *pNode_ids = nullptr;
    double * pValues = nullptr;

    pNodex = new float[MAXP4ESTNODES];
    pNodey = new float[MAXP4ESTNODES];
    pNode_ids = new long int [MAXP4ESTNODES];
    pValues = new double[MAXP4ESTNODES];

    for(std::pair<DoublePair,long> element : NodeIDs)
    {
        pNodex[element.second]=element.first.first;
        pNodey[element.second]=element.first.second;
        pNode_ids[element.second]=element.second;
    }

    if(current_solution.size() != 0)
        for(std::pair<DoublePair,long> element : NodeIDs)
        {
            pValues[element.second]=current_solution.at(element.second);
        }

    // Array of the coordinate arrays
    float * pCoordinates[2];
    pCoordinates[0]=pNodex;
    pCoordinates[1]=pNodey;

    // Create the file
    dbfile = DBCreate(fn.c_str(), DB_CLOBBER, DB_LOCAL, getDescription().c_str(), driver);
    if (!dbfile)
        throw escript::IOError("dump: Could not create Silo file");

    // create the nodelist
    std::vector<int> nodelist;
    long ids[4]={0};
    for(p4est_topidx_t treeid = p4est->first_local_tree; treeid <= p4est->last_local_tree; ++treeid) {
        p4est_tree_t * tree = p4est_tree_array_index(p4est->trees, treeid);
        sc_array_t * tquadrants = &tree->quadrants;
        p4est_locidx_t Q = (p4est_locidx_t) tquadrants->elem_count;
        for(int q = 0; q < Q; q++)
        {
            p4est_quadrant_t * quad = p4est_quadrant_array_index(tquadrants, q);          
            getNeighouringNodeIDs(quad, treeid, ids);
            nodelist.push_back(ids[0]);
            nodelist.push_back(ids[2]);
            nodelist.push_back(ids[3]);
            nodelist.push_back(ids[1]);
        }
    }

    int* nodelistarray = &nodelist[0];

    // write mesh
    int lnodelist = nodelist.size();
    int shapesize[] = {4};
    int shapecounts[] = {lnodelist/4};
    int nshapetypes = 1;
    int const shapetype[1] = {DB_ZONETYPE_QUAD};

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

    // Node values
    if(current_solution.size() != 0)
        DBPutPointvar1(dbfile, "u", "nodes", pValues, getNumNodes(), DB_DOUBLE, NULL);    

    DBClose(dbfile);

    delete [] pNodex;
    delete [] pNodey;
    delete [] pNode_ids;
    delete [] pValues;

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
            throw OxleyException("Unknown Error.");
        case Elements:
        case ReducedElements:
            return &m_elementId[0];
        case FaceElements:
        case ReducedFaceElements:
            throw OxleyException("borrowSampleReferenceIDs: FaceIDs");
            // return &m_faceId[0];
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

void Rectangle::refineMesh(std::string algorithmname)
{
    forestData.current_solution = &current_solution;    
    forestData.NodeIDs = &NodeIDs;

    if(!algorithmname.compare("uniform"))
    {
        p4est_refine_ext(p4est, true, -1, refine_uniform, init_rectangle_data, refine_copy_parent_quadrant);
        p4est_balance_ext(p4est, P4EST_CONNECT_FULL, init_rectangle_data, refine_copy_parent_quadrant);
    }
    else if(!algorithmname.compare("MARE2DEM") || !algorithmname.compare("mare2dem"))
    {
        if(adaptive_refinement == true)
        {
            p4est_refine_ext(p4est, true, -1, refine_mare2dem, init_rectangle_data, refine_copy_parent_quadrant);
            p4est_balance_ext(p4est, P4EST_CONNECT_FULL, init_rectangle_data, refine_copy_parent_quadrant);
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
    updateNodeIncrements();
    renumberNodes();
    updateRowsColumns();
    updateElementIds();
    updateFaceOffset();
    updateFaceElementCount();
}

void Rectangle::refineBoundary(std::string boundaryname, double dx)
{
    forestData.refinement_depth = dx;

    if(!boundaryname.compare("top") || !boundaryname.compare("Top")
        || !boundaryname.compare("t") || !boundaryname.compare("T")
        || !boundaryname.compare("TOP"))
    {
        p4est_refine_ext(p4est, true, -1, refine_north, init_rectangle_data, refine_copy_parent_quadrant);
        p4est_balance_ext(p4est, P4EST_CONNECT_FULL, init_rectangle_data, refine_copy_parent_quadrant);
    } 
    else if(!boundaryname.compare("bottom") || !boundaryname.compare("Bottom")
        || !boundaryname.compare("b") || !boundaryname.compare("B")
        || !boundaryname.compare("BOTTOM"))
    {
        p4est_refine_ext(p4est, true, -1, refine_south, init_rectangle_data, refine_copy_parent_quadrant);
        p4est_balance_ext(p4est, P4EST_CONNECT_FULL, init_rectangle_data, refine_copy_parent_quadrant);
    }
    else if(!boundaryname.compare("left") || !boundaryname.compare("Left")
        || !boundaryname.compare("l") || !boundaryname.compare("L")
        || !boundaryname.compare("LEFT"))
    {
        p4est_refine_ext(p4est, true, -1, refine_west, init_rectangle_data, refine_copy_parent_quadrant);
        p4est_balance_ext(p4est, P4EST_CONNECT_FULL, init_rectangle_data, refine_copy_parent_quadrant);
    }
    else if(!boundaryname.compare("right") || !boundaryname.compare("Right")
        || !boundaryname.compare("r") || !boundaryname.compare("R")
        || !boundaryname.compare("RIGHT"))
    {
        p4est_refine_ext(p4est, true, -1, refine_east, init_rectangle_data, refine_copy_parent_quadrant);
        p4est_balance_ext(p4est, P4EST_CONNECT_FULL, init_rectangle_data, refine_copy_parent_quadrant);  
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
    updateNodeIncrements();
    renumberNodes();
    updateRowsColumns();
}

void Rectangle::refineRegion(double x0, double x1, double y0, double y1)
{
    // If the boundaries were not specified by the user, default to the border of the domain
    forestData.refinement_boundaries[0] = x0 == -1 ? forestData.m_origin[0] : x0; 
    forestData.refinement_boundaries[1] = x1 == -1 ? forestData.m_origin[1] : x1;
    forestData.refinement_boundaries[2] = y0 == -1 ? forestData.m_lxy[0] : y0;
    forestData.refinement_boundaries[3] = y1 == -1 ? forestData.m_lxy[1] : y1;

    p4est_refine_ext(p4est, true, -1, refine_region, init_rectangle_data, refine_copy_parent_quadrant);
    p4est_balance_ext(p4est, P4EST_CONNECT_FULL, init_rectangle_data, refine_copy_parent_quadrant);

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
    updateNodeIncrements();
    renumberNodes();
    updateRowsColumns();
}

void Rectangle::refinePoint(double x0, double y0)
{
    // Check that the point is inside the domain
    if(x0 < forestData.m_origin[0] || x0 > forestData.m_lxy[0] 
        || y0 < forestData.m_origin[1] || y0 > forestData.m_lxy[1] )
    {
        throw OxleyException("Coordinates lie outside the domain.");
    }

    // If the boundaries were not specified by the user, default to the border of the domain
    forestData.refinement_boundaries[0] = x0;
    forestData.refinement_boundaries[1] = y0;
    p4est_refine_ext(p4est, true, -1, refine_point, init_rectangle_data, refine_copy_parent_quadrant);
    p4est_balance_ext(p4est, P4EST_CONNECT_FULL, init_rectangle_data, refine_copy_parent_quadrant);

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
    updateNodeIncrements();
    renumberNodes();
    updateRowsColumns();
}

void Rectangle::refineCircle(double x0, double y0, double r)
{
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
    p4est_refine_ext(p4est, true, -1, refine_circle, init_rectangle_data, refine_copy_parent_quadrant);
    p4est_balance_ext(p4est, P4EST_CONNECT_FULL, init_rectangle_data, refine_copy_parent_quadrant);

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
    updateNodeIncrements();
    renumberNodes();
    updateRowsColumns();
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
        int8_t c = face_code & 0x03;
        int8_t ishanging0 = (face_code >> 2) & 0x01;
        int8_t ishanging1 = (face_code >> 3) & 0x01;

        int8_t f0 = p4est_corner_faces[c][0];
        int8_t f1 = p4est_corner_faces[c][1];

        // int d0 = c / 2 == 0;
        // int d1 = c % 2 == 0;
        
        return ((ishanging0 == 1) && (p4est_corner_face_corners[c][f0] == 1)) 
            || ((ishanging1 == 1) && (p4est_corner_face_corners[c][f1] == 1));
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
#ifdef OXLEY_ENABLE_DEBUG_NODES
    std::cout << "Renumbering nodes... " << std::endl;
#endif

    NodeIDs.clear();
    quadrantIDs.clear();
// #pragma omp for
    for(p4est_topidx_t treeid = p4est->first_local_tree; treeid <= p4est->last_local_tree; ++treeid) {
        p4est_tree_t * tree = p4est_tree_array_index(p4est->trees, treeid);
        sc_array_t * tquadrants = &tree->quadrants;
        p4est_locidx_t Q = (p4est_locidx_t) tquadrants->elem_count;
        for(int q = 0; q < Q; ++q) { 
            p4est_quadrant_t * quad = p4est_quadrant_array_index(tquadrants, q);
            p4est_qcoord_t l = P4EST_QUADRANT_LEN(quad->level);
            p4est_qcoord_t lxy[4][2] = {{0,0},{l,0},{0,l},{l,l}};
            for(int n = 0; n < 4; n++)
            {
                double xy[3];
                p4est_qcoord_to_vertex(p4est->connectivity, treeid, quad->x+lxy[n][0], quad->y+lxy[n][1], xy);
                if(NodeIDs.count(std::make_pair(xy[0],xy[1]))==0)
                {
#ifdef OXLEY_ENABLE_DEBUG_NODES
                    std::cout << NodeIDs.size() << ": " << xy[0] << ", " << xy[1] << std::endl;
#endif
                    NodeIDs[std::make_pair(xy[0],xy[1])]=NodeIDs.size();
                }

                if(n==0)
                {
                    p4est_qcoord_to_vertex(p4est->connectivity, treeid, quad->x, quad->y, xy);
                    quadrantIDs.push_back(NodeIDs.find(std::make_pair(xy[0],xy[1]))->second);
#ifdef OXLEY_ENABLE_DEBUG_NODES
                    std::cout << "New Quad. ID: " << NodeIDs.find(std::make_pair(xy[0],xy[1]))->second << std::endl;
#endif
                }
            }
        }
    }


    // Populate m_nodeIDs
    m_nodeId.clear();
    m_nodeId.resize(NodeIDs.size());
    int count=0;
    for(std::pair<DoublePair,long> e : NodeIDs)
        m_nodeId[count++]=e.second;

#ifdef OXLEY_PRINT_NODEIDS
    std::cout << "Printing NodeIDs " << std::endl;
    double xyf[MAXP4ESTNODES][2]={{0}};
    for(std::pair<DoublePair,long> e : NodeIDs)
    {
        xyf[e.second][0]=e.first.first;
        xyf[e.second][1]=e.first.second;
    }
    for(int i=0; i<NodeIDs.size(); i++)
        std::cout << i << ": " << xyf[i][0] << ", " << xyf[i][1] << std::endl;
#endif
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

    for(p4est_topidx_t treeid = p4est->first_local_tree; treeid <= p4est->last_local_tree; ++treeid) {
        p4est_tree_t * tree = p4est_tree_array_index(p4est->trees, treeid);
        sc_array_t * tquadrants = &tree->quadrants;
        p4est_locidx_t Q = (p4est_locidx_t) tquadrants->elem_count;

#pragma omp parallel for
        for(int q = 0; q < Q; ++q) { // Loop over the elements attached to the tree
            p4est_quadrant_t * quad = p4est_quadrant_array_index(tquadrants, q);
            p4est_qcoord_t length = P4EST_QUADRANT_LEN(quad->level);

            // Loop over the four corners of the quadrant
            for(int n = 0; n < 4; ++n){
                int k = q - Q + nodeIncrements[treeid - p4est->first_local_tree];
                double lx = length * ((int) (n % 2) == 1);
                double ly = length * ((int) (n / 2) == 1);
                double xy[3];
                p4est_qcoord_to_vertex(p4est->connectivity, treeid, quad->x+lx, quad->y+ly, xy);

                if( (n == 0) 
                  || isHangingNode(nodes->face_code[k], n)
                  || isUpperBoundaryNode(quad, n, treeid, length) 
                )
                {
                    long lni = NodeIDs.find(std::make_pair(xy[0],xy[1]))->second;
                    double * point = arg.getSampleDataRW(lni);
                    point[0] = xy[0];
                    point[1] = xy[1];
#ifdef OXLEY_ENABLE_DEBUG_ASSEMBLE_COORDINATES
                    // float bounds[4]={0.0};
                    if(xy[0] < bounds[0]) bounds[0]=xy[0];
                    if(xy[0] > bounds[1]) bounds[1]=xy[0];
                    if(xy[1] < bounds[2]) bounds[2]=xy[1];
                    if(xy[1] > bounds[3]) bounds[3]=xy[1];
#endif
                }
            }
        }
    }
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
    sc_array_t * tquadrants = &currenttree->quadrants;
    p4est_quadrant_t * quad = p4est_quadrant_array_index(tquadrants, e);
    p4est_qcoord_t l = P4EST_QUADRANT_LEN(quad->level);
    int lxy[4][2] = {{0,0},{l,0},{0,l},{l,l}};

// #pragma omp for
    for(int i = 0; i < 4; i++)
    {
        double xy[3];
        p4est_qcoord_to_vertex(p4est->connectivity, t, quad->x+lxy[i][0], quad->y+lxy[i][1], xy);
        rowIndex[i] = NodeIDs.find(std::make_pair(xy[0],xy[1]))->second;
    }

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
    long rowIndex[4] = {0};
    getNeighouringNodeIDs(quad.quad, quad.treeid, rowIndex);
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

        for(p4est_topidx_t treeid = p4est->first_local_tree; treeid <= p4est->last_local_tree; ++treeid) {
            p4est_tree_t * tree = p4est_tree_array_index(p4est->trees, treeid);
            sc_array_t * tquadrants = &tree->quadrants;
            p4est_locidx_t Q = (p4est_locidx_t) tquadrants->elem_count;
            #pragma omp parallel for
            for(int q = 0; q < Q; q++)
            {
                p4est_quadrant_t * quad = p4est_quadrant_array_index(tquadrants, q);
               
                long ids[4]={0};
                getNeighouringNodeIDs(quad, treeid, ids);
                int quadID=getQuadID(ids[0]);

                memcpy(&f_00[0], in.getSampleDataRO(ids[0],sentinel), numComp*sizeof(S));
                memcpy(&f_01[0], in.getSampleDataRO(ids[1],sentinel), numComp*sizeof(S));
                memcpy(&f_10[0], in.getSampleDataRO(ids[2],sentinel), numComp*sizeof(S));
                memcpy(&f_11[0], in.getSampleDataRO(ids[3],sentinel), numComp*sizeof(S));
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

        for(p4est_topidx_t treeid = p4est->first_local_tree; treeid <= p4est->last_local_tree; ++treeid) {
            p4est_tree_t * tree = p4est_tree_array_index(p4est->trees, treeid);
            sc_array_t * tquadrants = &tree->quadrants;
            p4est_locidx_t Q = (p4est_locidx_t) tquadrants->elem_count;
            #pragma omp parallel for
            for(int q = 0; q < Q; q++)
            {
                p4est_quadrant_t * quad = p4est_quadrant_array_index(tquadrants, q);

                long ids[4]={0};
                getNeighouringNodeIDs(quad, treeid, ids);
                long quadId = getQuadID(ids[0]);

                memcpy(&f_00[0], in.getSampleDataRO(ids[0], sentinel), numComp*sizeof(S));
                memcpy(&f_01[0], in.getSampleDataRO(ids[1], sentinel), numComp*sizeof(S));
                memcpy(&f_10[0], in.getSampleDataRO(ids[2], sentinel), numComp*sizeof(S));
                memcpy(&f_11[0], in.getSampleDataRO(ids[3], sentinel), numComp*sizeof(S));

        #ifdef OXLEY_ENABLE_DEBUG_INTERPOLATE
                std::cout << "interpolateNodesOnElementsWorker quadID: " << quadId << ", node IDs " << ids[0] << ", " << ids[1] << ", " << ids[2] << ", " << ids[3] << std::endl;
        #endif
                
                S* o = out.getSampleDataRW(quadId, sentinel);
                for (index_t i=0; i < numComp; ++i) {
                    o[INDEX2(i,numComp,0)] = c0*(f_01[i] + f_10[i]) + c1*f_11[i] + c2*f_00[i];
                    o[INDEX2(i,numComp,1)] = c0*(f_00[i] + f_11[i]) + c1*f_01[i] + c2*f_10[i];
                    o[INDEX2(i,numComp,2)] = c0*(f_00[i] + f_11[i]) + c1*f_10[i] + c2*f_01[i];
                    o[INDEX2(i,numComp,3)] = c0*(f_01[i] + f_10[i]) + c1*f_00[i] + c2*f_11[i];
                }
            }
        }
    }
}

//protected
void Rectangle::getNeighouringNodeIDs(p4est_quadrant_t * quad, p4est_topidx_t treeid, long (&ids) [4]) const
{
    p4est_qcoord_t l = P4EST_QUADRANT_LEN(quad->level);
    int adj[4][2]={{0,0},{l,0},{0,l},{l,l}};
#pragma omp parallel for
    for(int i=0; i<4;i++)
    {
        double xy[3];
        p4est_qcoord_to_vertex(p4est->connectivity, treeid, quad->x+adj[i][0], quad->y+adj[i][1], xy);
        ids[i]=(long) NodeIDs.find(std::make_pair(xy[0],xy[1]))->second;
    }
}

long Rectangle::getQuadID(long nodeid) const
{
    for(int i = 0; i < quadrantIDs.size(); i++)
        if(quadrantIDs[i]==nodeid)
            return i;
    throw OxleyException("Unknown error");
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
            for (index_t k=0; k<NodeIDsLeft.size()-1; k++) {
                borderNodeInfo tmp = NodeIDsLeft[k];
                memcpy(&f_00[0], in.getSampleDataRO(tmp.neighbours[0], sentinel), numComp*sizeof(S));
                memcpy(&f_01[0], in.getSampleDataRO(tmp.neighbours[2], sentinel), numComp*sizeof(S));
                
                S* o = out.getSampleDataRW(m_faceOffset[0]+k, sentinel);
                for (index_t i=0; i < numComp; ++i) {
                    o[INDEX2(i,numComp,0)] = (f_00[i] + f_01[i])/static_cast<S>(2);
                } /* end of component loop i */
            } 
        } /* end of face 0 */
        if (m_faceOffset[1] > -1) {
#pragma omp for nowait
            for (index_t k=0; k<NodeIDsRight.size()-1; k++) {
                borderNodeInfo tmp = NodeIDsRight[k];
                memcpy(&f_10[0], in.getSampleDataRO(tmp.neighbours[1], sentinel), numComp*sizeof(S));
                memcpy(&f_11[0], in.getSampleDataRO(tmp.neighbours[3], sentinel), numComp*sizeof(S));
                S* o = out.getSampleDataRW(m_faceOffset[1]+k, sentinel);
                for (index_t i=0; i < numComp; ++i) {
                    o[INDEX2(i,numComp,0)] = (f_10[i] + f_11[i])/static_cast<S>(2);
                } /* end of component loop i */
            } 
        } /* end of face 1 */
        if (m_faceOffset[2] > -1) {
#pragma omp for nowait
            for (index_t k=0; k<NodeIDsBottom.size()-1; k++) {
                borderNodeInfo tmp = NodeIDsBottom[k];
                memcpy(&f_00[0], in.getSampleDataRO(tmp.neighbours[0], sentinel), numComp*sizeof(S));
                memcpy(&f_10[0], in.getSampleDataRO(tmp.neighbours[1], sentinel), numComp*sizeof(S));
                S* o = out.getSampleDataRW(m_faceOffset[2]+k, sentinel);
                for (index_t i=0; i < numComp; ++i) {
                    o[INDEX2(i,numComp,0)] = (f_00[i] + f_10[i])/static_cast<S>(2);
                } /* end of component loop i */
            } 
        } /* end of face 2 */
        if (m_faceOffset[3] > -1) {
#pragma omp for nowait
            for (index_t k=0; k<NodeIDsTop.size()-1; k++) {
                borderNodeInfo tmp = NodeIDsTop[k];
                memcpy(&f_01[0], in.getSampleDataRO(tmp.neighbours[2], sentinel), numComp*sizeof(S));
                memcpy(&f_11[0], in.getSampleDataRO(tmp.neighbours[3], sentinel), numComp*sizeof(S));
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
            for (index_t k=0; k<NodeIDsLeft.size()-1; k++) {
                borderNodeInfo tmp = NodeIDsLeft[k];
                memcpy(&f_00[0], in.getSampleDataRO(tmp.neighbours[0], sentinel), numComp*sizeof(S));
                memcpy(&f_01[0], in.getSampleDataRO(tmp.neighbours[2], sentinel), numComp*sizeof(S));
                S* o = out.getSampleDataRW(m_faceOffset[0]+k, sentinel);
                for (index_t i=0; i < numComp; ++i) {
                    o[INDEX2(i,numComp,0)] = c0*f_01[i] + c1*f_00[i];
                    o[INDEX2(i,numComp,1)] = c0*f_00[i] + c1*f_01[i];
                } /* end of component loop i */
            }
        } /* end of face 0 */
        if (m_faceOffset[1] > -1) {
    #pragma omp for nowait
            for (index_t k=0; k<NodeIDsRight.size()-1; k++) {
                borderNodeInfo tmp = NodeIDsRight[k];
                memcpy(&f_10[0], in.getSampleDataRO(tmp.neighbours[1], sentinel), numComp*sizeof(S));
                memcpy(&f_11[0], in.getSampleDataRO(tmp.neighbours[3], sentinel), numComp*sizeof(S));
                S* o = out.getSampleDataRW(m_faceOffset[1]+k, sentinel);
                for (index_t i=0; i < numComp; ++i) {
                    o[INDEX2(i,numComp,0)] = c1*f_10[i] + c0*f_11[i];
                    o[INDEX2(i,numComp,1)] = c1*f_11[i] + c0*f_10[i];
                } /* end of component loop i */
            } 
        } /* end of face 1 */
        if (m_faceOffset[2] > -1) {
    #pragma omp for nowait
             for (index_t k=0; k<NodeIDsBottom.size()-1; k++) {
                borderNodeInfo tmp = NodeIDsTop[k];
                memcpy(&f_00[0], in.getSampleDataRO(tmp.neighbours[0], sentinel), numComp*sizeof(S));
                memcpy(&f_10[0], in.getSampleDataRO(tmp.neighbours[1], sentinel), numComp*sizeof(S));
                S* o = out.getSampleDataRW(m_faceOffset[2]+k, sentinel);
                for (index_t i=0; i < numComp; ++i) {
                    o[INDEX2(i,numComp,0)] = c0*f_10[i] + c1*f_00[i];
                    o[INDEX2(i,numComp,1)] = c0*f_00[i] + c1*f_10[i];
                } /* end of component loop i */
            } 
        } /* end of face 2 */
        if (m_faceOffset[3] > -1) {
    #pragma omp for nowait
            for (index_t k=0; k<NodeIDsTop.size()-1; k++) {
                borderNodeInfo tmp = NodeIDsTop[k];
                memcpy(&f_01[0], in.getSampleDataRO(tmp.neighbours[2], sentinel), numComp*sizeof(S));
                memcpy(&f_11[0], in.getSampleDataRO(tmp.neighbours[3], sentinel), numComp*sizeof(S));
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
    return m_dofMap[node];
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
    return NodeIDs.size();
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
    return getNumNodes();
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
    std::vector<std::vector<long>> * indices;
    indices = new std::vector<std::vector<long>>;
    int initial[] = {0, -1, -1, -1, -1};
    indices->resize(getNumNodes(), std::vector<long>(initial, initial+5));

    update_RC_data * data;
    data = new update_RC_data;
    data->indices = indices;
    data->pNodeIDs = &NodeIDs;
    data->p4est = p4est;
    data->m_origin[0]=forestData.m_origin[0];
    data->m_origin[1]=forestData.m_origin[1];

    // This function loops over all interior faces
    // Note that it does not loop over the nodes on the boundaries
    // x = Lx and y = Ly
    p4est_iterate_ext(p4est, NULL, data, NULL, update_RC, NULL, true);

    // Find the indices of the nodes on the boundaries x = Lx and y = Ly
    for(p4est_topidx_t treeid = p4est->first_local_tree; treeid <= p4est->last_local_tree; ++treeid) {
        p4est_tree_t * tree = p4est_tree_array_index(p4est->trees, treeid);
        sc_array_t * tquadrants = &tree->quadrants;
        p4est_locidx_t Q = (p4est_locidx_t) tquadrants->elem_count;
// #pragma omp parallel for
        for(int q = 0; q < Q; ++q) { // Loop over the elements attached to the tree
            p4est_quadrant_t * quad = p4est_quadrant_array_index(tquadrants, q);
            p4est_qcoord_t length = P4EST_QUADRANT_LEN(quad->level);
            double xy[3];
            p4est_qcoord_to_vertex(p4est->connectivity, treeid, quad->x+length, quad->y, xy);

            // If the node is on the boundary x=Lx or y=Ly
            if(xy[0] == forestData.m_lxy[0]) 
            {
                // Get the node IDs
                long lni0 = NodeIDs.find(std::make_pair(xy[0],xy[1]))->second;
                p4est_qcoord_to_vertex(p4est->connectivity, treeid, quad->x+length, quad->y+length, xy);
                long lni1 = NodeIDs.find(std::make_pair(xy[0],xy[1]))->second;

                std::vector<long> * idx0 = &indices[0][lni0];
                std::vector<long> * idx1 = &indices[0][lni1];

                // Check for duplicates
                bool dup = false;
                for(int i = 1; i < idx0[0][0] + 1; i++)
                    if(idx0[0][i] == lni1)
                        dup = true;

                // Add the new indices
                if(dup == false)
                {
                    idx0[0][0]++;
                    idx1[0][0]++;
                    idx0[0][idx0[0][0]]=lni1;
                    idx1[0][idx1[0][0]]=lni0;
                }
            }

            p4est_qcoord_to_vertex(p4est->connectivity, treeid, quad->x, quad->y+length, xy);
            if(xy[1] == forestData.m_lxy[1])
            {
                // Get the node IDs
                long lni0 = NodeIDs.find(std::make_pair(xy[0],xy[1]))->second;
                p4est_qcoord_to_vertex(p4est->connectivity, treeid, quad->x+length, quad->y+length, xy);
                long lni1 = NodeIDs.find(std::make_pair(xy[0],xy[1]))->second;

                std::vector<long> * idx0 = &indices[0][lni0];
                std::vector<long> * idx1 = &indices[0][lni1];

                // Check for duplicates
                bool dup = false;
                for(int i = 1; i < idx0[0][0] + 1; i++)
                    if(idx0[0][i] == lni1)
                        dup = true;

                // Add the new indices
                if(dup == false)
                {
                    idx0[0][0]++;
                    idx1[0][0]++;
                    idx0[0][idx0[0][0]]=lni1;
                    idx1[0][idx1[0][0]]=lni0;
                }
            }
        }
    }

    // Sorting
#pragma omp for
    for(int i = 0; i < getNumNodes(); i++){
        std::vector<long> * idx0 = &indices[0][i];
        std::sort(indices[0][i].begin()+1, indices[0][i].begin()+idx0[0][0]+1);
    }

#ifdef OXLEY_ENABLE_DEBUG_NODES
    std::cout << "Node connections: " << std::endl;
    // Output for debugging
    for(int i = 0; i < getNumNodes(); i++){
        std::vector<long> * idx0 = &indices[0][i];
        std::cout << i << ": ";
        for(int j = 1; j < idx0[0][0]+1; j++)
            std::cout << idx0[0][j] << ", ";
        std::cout << std::endl;
    }
#endif

    // Convert to CRS format
    myRows.clear();
    myRows.push_back(0);
    myColumns.clear();
    m_dofMap.assign(getNumNodes(), 0);
    long counter = 0;
    for(int i = 0; i < getNumNodes(); i++)
    {
        std::vector<long> * idx0 = &indices[0][i];
        std::vector<long> temp; 
        for(int j = 1; j < idx0[0][0]+1; j++)
        {
            temp.push_back(idx0[0][j]);
            counter++;
        }
        std::sort(temp.begin(),temp.end());
        for(int i = 0; i < temp.size(); i++)
            myColumns.push_back(temp[i]);
        m_dofMap[i] = counter-myRows[i];
        if(i < getNumNodes()-1)
            myRows.push_back(counter);
    }
    myRows.push_back(myColumns.size());
#ifdef OXLEY_ENABLE_DEBUG_NODES
    std::cout << "Converted to Yale format... "<< std::endl;
    std::cout << "COL_INDEX [";
    for(auto i = myColumns.begin(); i < myColumns.end(); i++)
        std::cout << *i << " ";
    std::cout << "]" << std::endl;
    std::cout << "ROW_INDEX [";
    for(auto i = myRows.begin(); i < myRows.end(); i++)
        std::cout << *i << " ";
    std::cout << "]" << std::endl;
    std::cout << "m_dofMap [";
    for(auto i = m_dofMap.begin(); i < m_dofMap.end(); i++)
        std::cout << *i << " ";
    std::cout << "]" << std::endl;
#endif

    delete data;
    delete indices;
}

#ifdef ESYS_HAVE_TRILINOS
//protected
esys_trilinos::const_TrilinosGraph_ptr Rectangle::getTrilinosGraph() const
{   
    if (m_graph.is_null()) {
        m_graph = createTrilinosGraph(myRows, myColumns);
    }
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
    // degrees of freedom are numbered from left to right, bottom to top in
    // each rank, continuing on the next rank (ranks also go left-right,
    // bottom-top).
    // This means rank 0 has id 0...n0-1, rank 1 has id n0...n1-1 etc. which
    // helps when writing out data rank after rank.

    // build node distribution vector first.
    // rank i owns m_nodeDistribution[i+1]-nodeDistribution[i] nodes which is
    // constant for all ranks in this implementation
    
    // m_nodeDistribution.assign(p4est->mpisize+1, 0);
    // for (dim_t k=1; k<m_mpiInfo->size; k++) {
    //     m_nodeDistribution[k]=p4est->mpirank;
    // }
    // m_nodeDistribution[p4est->mpisize]=getNumDataPointsGlobal();
    // try {
        // m_nodeId.resize(getNumNodes());
        // m_dofId.resize(numDOF);
    //     m_elementId.resize(getNumElements());
    // } catch (const std::length_error& le) {
    //     throw OxleyException("The system does not have sufficient memory for a domain of this size.");
    // }

    // // populate face element counts
    // updateFaceElementCount();

    // const dim_t NFE = getNumFaceElements();
    // m_faceId.resize(NFE);

    // const index_t left = getFirstInDim(0);
    // const index_t bottom = getFirstInDim(1);
    // const dim_t nDOF0 = getNumDOFInAxis(0);
    // const dim_t nDOF1 = getNumDOFInAxis(1);
    // const dim_t NE0 = m_NE[0];
    // const dim_t NE1 = m_NE[1];

// #define globalNodeId(x,y) \
//     ((m_offset[0]+x)/nDOF0)*nDOF0*nDOF1+(m_offset[0]+x)%nDOF0 \
//     + ((m_offset[1]+y)/nDOF1)*nDOF0*nDOF1*m_NX[0]+((m_offset[1]+y)%nDOF1)*nDOF0

//     // set corner id's outside the parallel region
    // m_nodeId[0] = globalNodeId(0, 0);
//     m_nodeId[m_NN[0]-1] = globalNodeId(m_NN[0]-1, 0);
//     m_nodeId[m_NN[0]*(m_NN[1]-1)] = globalNodeId(0, m_NN[1]-1);
//     m_nodeId[m_NN[0]*m_NN[1]-1] = globalNodeId(m_NN[0]-1,m_NN[1]-1);
// #undef globalNodeId


        // populate degrees of freedom and own nodes (identical id)
        // for (dim_t i=0; i<nDOF1; i++) {
        //     for (dim_t j=0; j<nDOF0; j++) {
        //         const index_t nodeIdx=j+left+(i+bottom)*m_NN[0];
        //         const index_t dofIdx=j+i*nDOF0;
        //         m_dofId[dofIdx] = m_nodeId[nodeIdx]
        //             = m_nodeDistribution[m_mpiInfo->rank]+dofIdx;
        //     }
        // }

        // populate the rest of the nodes (shared with other ranks)
//         if (m_faceCount[0]==0) { // left column
// #pragma omp for nowait
//             for (dim_t i=0; i<nDOF1; i++) {
//                 const index_t nodeIdx=(i+bottom)*m_NN[0];
//                 const index_t dofId=(i+1)*nDOF0-1;
//                 m_nodeId[nodeIdx]
//                     = m_nodeDistribution[m_mpiInfo->rank-1]+dofId;
//             }
//         }
//         if (m_faceCount[1]==0) { // right column
// #pragma omp for nowait
//             for (dim_t i=0; i<nDOF1; i++) {
//                 const index_t nodeIdx=(i+bottom+1)*m_NN[0]-1;
//                 const index_t dofId=i*nDOF0;
//                 m_nodeId[nodeIdx]
//                     = m_nodeDistribution[m_mpiInfo->rank+1]+dofId;
//             }
//         }
//         if (m_faceCount[2]==0) { // bottom row
// #pragma omp for nowait
//             for (dim_t i=0; i<nDOF0; i++) {
//                 const index_t nodeIdx=i+left;
//                 const index_t dofId=nDOF0*(nDOF1-1)+i;
//                 m_nodeId[nodeIdx]
//                     = m_nodeDistribution[m_mpiInfo->rank-m_NX[0]]+dofId;
//             }
//         }
//         if (m_faceCount[3]==0) { // top row
// #pragma omp for nowait
//             for (dim_t i=0; i<nDOF0; i++) {
//                 const index_t nodeIdx=m_NN[0]*(m_NN[1]-1)+i+left;
//                 const index_t dofId=i;
//                 m_nodeId[nodeIdx]
//                     = m_nodeDistribution[m_mpiInfo->rank+m_NX[0]]+dofId;
//             }
//         }

//         // populate element id's
// #pragma omp for nowait
//         for (dim_t i1=0; i1<NE1; i1++) {
//             for (dim_t i0=0; i0<NE0; i0++) {
//                 m_elementId[i0+i1*NE0]=(m_offset[1]+i1)*m_gNE[0]+m_offset[0]+i0;
//             }
//         }

//         // face elements
// #pragma omp for
//         for (dim_t k=0; k<NFE; k++)
//             m_faceId[k]=k;

    m_nodeTags.assign(getNumNodes(), 0);
    updateTagsInUse(Nodes);

    m_elementTags.assign(getNumElements(), 0);
    updateTagsInUse(Elements);

    // generate face offset vector and set face tags
    // const index_t LEFT=1, RIGHT=2, BOTTOM=10, TOP=20;
    // const index_t faceTag[] = { LEFT, RIGHT, BOTTOM, TOP };
    // m_faceOffset.assign(4, -1);
    // m_faceTags.clear();
    // index_t offset=0;
    // for (size_t i=0; i<4; i++) {
    //     if (m_faceCount[i]>0) {
    //         m_faceOffset[i]=offset;
    //         offset+=m_faceCount[i];
    //         m_faceTags.insert(m_faceTags.end(), m_faceCount[i], faceTag[i]);
    //     }
    // }
    // setTagMap("left", LEFT);
    // setTagMap("right", RIGHT);
    // setTagMap("bottom", BOTTOM);
    // setTagMap("top", TOP);
    // updateTagsInUse(FaceElements);

//     populateDofMap();

}

void Rectangle::updateFaceElementCount()
{
    for(int i = 0; i < 4; i++)
        m_faceCount[i]=-1;

    NodeIDsTop.clear();
    NodeIDsBottom.clear();
    NodeIDsLeft.clear();
    NodeIDsRight.clear();

    for(p4est_topidx_t treeid = p4est->first_local_tree; treeid <= p4est->last_local_tree; ++treeid) 
    {
        p4est_tree_t * tree = p4est_tree_array_index(p4est->trees, treeid);
        sc_array_t * tquadrants = &tree->quadrants;
        p4est_locidx_t Q = (p4est_locidx_t) tquadrants->elem_count;
        for(int q = 0; q < Q; ++q) 
        {
            p4est_quadrant_t * quad = p4est_quadrant_array_index(tquadrants, q);
            p4est_qcoord_t l = P4EST_QUADRANT_LEN(quad->level);
            // int k = q - Q + nodeIncrements[treeid - p4est->first_local_tree];
            p4est_qcoord_t lxy[4][2] = {{0,0},{l,0},{0,l},{l,l}};
            double xy[4][3] = {0};
            int nodeids[4]={-1};
            for(int n = 0; n < 4; n++)
            {
                p4est_qcoord_to_vertex(p4est->connectivity, treeid, quad->x+lxy[n][0], quad->y+lxy[n][1], xy[n]);
                nodeids[n]=NodeIDs.find(std::make_pair(xy[n][0],xy[n][1]))->second;
            }

            for(int n = 0; n < 4; n++)
            {
                borderNodeInfo tmp;
                // tmp.nodeid=NodeIDs.find(std::make_pair(xy[n][0],xy[n][1]))->second;
                tmp.nodeid=nodeids[n];
                tmp.neighbours[0]=nodeids[0];
                tmp.neighbours[1]=nodeids[1];
                tmp.neighbours[2]=nodeids[2];
                tmp.neighbours[3]=nodeids[3];
                tmp.quad=quad;
                tmp.treeid=treeid;

                if(isLeftBoundaryNode(quad, n, treeid, l))
                {
                    NodeIDsLeft.push_back(tmp);
                    m_faceCount[0]++;
                }
                    
                if(isBottomBoundaryNode(quad, n, treeid, l))
                {
                    NodeIDsBottom.push_back(tmp);
                    m_faceCount[2]++;
                }

                if(isRightBoundaryNode(quad, n, treeid, l))
                {
                    NodeIDsRight.push_back(tmp);
                    m_faceCount[1]++;
                }
                    
                if(isTopBoundaryNode(quad, n, treeid, l))
                {
                    NodeIDsTop.push_back(tmp);
                    m_faceCount[3]++;
                }
            
            }
        }
    }

    // Remove duplicates
    for(int i = 1; i < NodeIDsLeft.size()-1; i++)
        if(NodeIDsLeft[i].treeid == NodeIDsLeft[i-1].treeid)
        {
            NodeIDsLeft.erase(NodeIDsLeft.begin()+i);
            i--;
            m_faceCount[0]--;
        }
    for(int i = 1; i < NodeIDsBottom.size()-1; i++)
        if(NodeIDsBottom[i].treeid == NodeIDsBottom[i-1].treeid)
        {
            NodeIDsBottom.erase(NodeIDsBottom.begin()+i);
            i--;
            m_faceCount[2]--;
        }
    for(int i = 1; i < NodeIDsRight.size()-1; i++)
        if(NodeIDsRight[i].treeid == NodeIDsRight[i-1].treeid)
        {
            NodeIDsRight.erase(NodeIDsRight.begin()+i);
            i--;
            m_faceCount[1]--;
        }
    for(int i = 1; i < NodeIDsTop.size()-1; i++)
        if(NodeIDsTop[i].treeid == NodeIDsTop[i-1].treeid)
        {
            NodeIDsTop.erase(NodeIDsTop.begin()+i);
            i--;
            m_faceCount[3]--;
        }


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
    for(int i = 0; i < NodeIDsLeft.size()-1;i++)
        std::cout << NodeIDsLeft[i].nodeid << " ";
    std::cout << std::endl;
    std::cout << "NodeIDsRight" << std::endl;
    for(int i = 0; i < NodeIDsRight.size()-1;i++)
        std::cout << NodeIDsRight[i].nodeid << " ";
    std::cout << std::endl;
    std::cout << "NodeIDsTop" << std::endl;
    for(int i = 0; i < NodeIDsTop.size()-1;i++)
        std::cout << NodeIDsTop[i].nodeid << " ";
    std::cout << std::endl;
    std::cout << "NodeIDsBottom" << std::endl;
    for(int i = 0; i < NodeIDsBottom.size()-1;i++)
        std::cout << NodeIDsBottom[i].nodeid << " ";
    std::cout << std::endl;
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
    for(p4est_topidx_t treeid = p4est->first_local_tree; treeid <= p4est->last_local_tree; ++treeid) 
    {
        p4est_tree_t * tree = p4est_tree_array_index(p4est->trees, treeid);
        sc_array_t * tquadrants = &tree->quadrants;
        p4est_locidx_t Q = (p4est_locidx_t) tquadrants->elem_count;
        for(int q = 0; q < Q; ++q) 
        { 
            p4est_quadrant_t * quad = p4est_quadrant_array_index(tquadrants, q);
            p4est_qcoord_t length = P4EST_QUADRANT_LEN(quad->level);
            for(int n = 0; n < 4; n++)
            {
                double lx = length * ((int) (n % 2) == 1);
                double ly = length * ((int) (n / 2) == 1);
                double xy[3];
                p4est_qcoord_to_vertex(p4est->connectivity, treeid, quad->x+lx, quad->y+ly, xy);
                m_nodeDistribution[counter++]=NodeIDs.find(std::make_pair(xy[0],xy[1]))->second;
            }
        }
    }
    m_nodeDistribution.shrink_to_fit();
}

// updates m_elementIDs()
void Rectangle::updateElementIds()
{
    m_elementId.clear();
    m_elementId.assign(MAXP4ESTNODES,0);
    int count=0;
    for(std::pair<DoublePair,long> e : NodeIDs)
        m_elementId[count++]=e.second;
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

    // Loop over interior quadrants
    // getConnections_data * data;
    // data = new getConnections_data;
    // data->pNodeIDs = &NodeIDs;
    // data->indices = &indices;
    // data->p4est = p4est;
    // p4est_iterate(p4est, NULL, data, update_connections, NULL, NULL);

    // Loop over the quadrants skipped by p4est_iterate  
    for(p4est_topidx_t treeid = p4est->first_local_tree; treeid <= p4est->last_local_tree; ++treeid) 
    {
        p4est_tree_t * tree = p4est_tree_array_index(p4est->trees, treeid);
        sc_array_t * tquadrants = &tree->quadrants;
        p4est_locidx_t Q = (p4est_locidx_t) tquadrants->elem_count;
// #pragma omp parallel for
        for(int q = 0; q < Q; ++q) // Loop over all quadrants
        { 
            p4est_quadrant_t * quad = p4est_quadrant_array_index(tquadrants, q);
            p4est_qcoord_t length = P4EST_QUADRANT_LEN(quad->level);
            for(int n = 0; n < 4; n++)
            {
                // if(!isUpperBoundaryNode(quad, n, treeid, length))
                //     continue;

                double xy[3];
                long lx[4] = {0,length,0,length};
                long ly[4] = {0,0,length,length};
                long lni[4] = {-1};
                for(int i = 0; i < 4; i++)
                {
                    p4est_qcoord_to_vertex(p4est->connectivity, treeid, quad->x+lx[i], quad->y+ly[i], xy);
                    lni[i] = NodeIDs.find(std::make_pair(xy[0],xy[1]))->second;
                }

                for(int i = 0; i < 4; i++)
                {
                    for(int j = 0; j < 4; j++)
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
    }

// Sorting
// #pragma omp parallel for
    for(int i = 0; i < numNodes; i++){
        std::sort(indices[i].begin(), indices[i].begin()+indices[i].size());
    }

// #ifdef OXLEY_ENABLE_DEBUG
//     std::cout << "Rectangle::getConnections" << std::endl;
//     for(int i = 0; i < numNodes; i++) {
//         std::cout << "i:" << i << " ";
//         for(auto j = 0; j < indices[i].size(); j++)
//             std::cout << indices[i][j] << ", ";
//         std::cout << std::endl;
//     }
// #endif

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
        double m_dx[2] = {m_NX[0]*forestData.m_dx[0][P4EST_MAXLEVEL-i], 
                          m_NX[1]*forestData.m_dx[1][P4EST_MAXLEVEL-i]};
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

        for (p4est_topidx_t t = p4est->first_local_tree; t <= p4est->last_local_tree; t++) // Loop over every tree
        {
            p4est_tree_t * currenttree = p4est_tree_array_index(p4est->trees, t);
            sc_array_t * tquadrants = &currenttree->quadrants;
            p4est_locidx_t Q = (p4est_locidx_t) tquadrants->elem_count;
#pragma omp parallel for
            for(int q = 0; q < Q; ++q) // Loop over every quadrant within the tree
            {
                p4est_quadrant_t * quad = p4est_quadrant_array_index(tquadrants, q);
                double xy[3];
                p4est_qcoord_to_vertex(p4est->connectivity, t, quad->x, quad->y, xy);
                long e = getQuadID(NodeIDs.find(std::make_pair(xy[0],xy[1]))->second);

                memcpy(&f_00[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar)); // TODO BUG? 
                memcpy(&f_01[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                memcpy(&f_10[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                memcpy(&f_11[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));

                Scalar* o = out.getSampleDataRW(e, zero);
                for(index_t i = 0; i < numComp; ++i) {
                    o[INDEX3(i,0,0,numComp,2)] = (f_10[i]-f_00[i])*cx[1][i] + (f_11[i]-f_01[i])*cx[0][i];
                    o[INDEX3(i,1,0,numComp,2)] = (f_01[i]-f_00[i])*cy[1][i] + (f_11[i]-f_10[i])*cy[0][i];
                    o[INDEX3(i,0,1,numComp,2)] = (f_10[i]-f_00[i])*cx[1][i] + (f_11[i]-f_01[i])*cx[0][i];
                    o[INDEX3(i,1,1,numComp,2)] = (f_01[i]-f_00[i])*cy[0][i] + (f_11[i]-f_10[i])*cy[1][i];
                    o[INDEX3(i,0,2,numComp,2)] = (f_10[i]-f_00[i])*cx[0][i] + (f_11[i]-f_01[i])*cx[1][i];
                    o[INDEX3(i,1,2,numComp,2)] = (f_01[i]-f_00[i])*cy[1][i] + (f_11[i]-f_10[i])*cy[0][i];
                    o[INDEX3(i,0,3,numComp,2)] = (f_10[i]-f_00[i])*cx[0][i] + (f_11[i]-f_01[i])*cx[1][i];
                    o[INDEX3(i,1,3,numComp,2)] = (f_01[i]-f_00[i])*cy[0][i] + (f_11[i]-f_10[i])*cy[1][i];
                } 
            }
        }
    } else if (out.getFunctionSpace().getTypeCode() == ReducedElements) {
        out.requireWrite();

        std::vector<Scalar> f_00(numComp, zero);
        std::vector<Scalar> f_01(numComp, zero);
        std::vector<Scalar> f_10(numComp, zero);
        std::vector<Scalar> f_11(numComp, zero);

        for (p4est_topidx_t t = p4est->first_local_tree; t <= p4est->last_local_tree; t++) // Loop over every tree
        {
            p4est_tree_t * currenttree = p4est_tree_array_index(p4est->trees, t);
            sc_array_t * tquadrants = &currenttree->quadrants;
            p4est_locidx_t Q = (p4est_locidx_t) tquadrants->elem_count;
#pragma omp parallel for
            for(int q = 0; q < Q; ++q) // Loop over every quadrant within the tree
            {

                p4est_quadrant_t * quad = p4est_quadrant_array_index(tquadrants, q);
                double xy[3];
                p4est_qcoord_to_vertex(p4est->connectivity, t, quad->x, quad->y, xy);
                long e = getQuadID(NodeIDs.find(std::make_pair(xy[0],xy[1]))->second);
                memcpy(&f_00[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                memcpy(&f_01[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                memcpy(&f_10[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                memcpy(&f_11[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));

                Scalar* o = out.getSampleDataRW(e, zero);

                for(index_t i = 0; i < numComp; ++i) {
                    o[INDEX3(i,0,0,numComp,2)] = (f_10[i] + f_11[i] - f_00[i] - f_01[i])*cx[2][i] * 0.5;
                    o[INDEX3(i,1,0,numComp,2)] = (f_01[i] + f_11[i] - f_00[i] - f_10[i])*cy[2][i] * 0.5;
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
#pragma omp parallel for
            for(int q = 0; q < Q; ++q) // Loop over every quadrant within the tree
            {

                p4est_quadrant_t * quad = p4est_quadrant_array_index(tquadrants, q);
                double xy[3];
                p4est_qcoord_to_vertex(p4est->connectivity, t, quad->x, quad->y, xy);
                long e = getQuadID(NodeIDs.find(std::make_pair(xy[0],xy[1]))->second);
                quadrantData * quaddata = (quadrantData *) quad->p.user_data;

                std::vector<Scalar> f_00(numComp, zero);
                std::vector<Scalar> f_01(numComp, zero);
                std::vector<Scalar> f_10(numComp, zero);
                std::vector<Scalar> f_11(numComp, zero);

                if (quaddata->m_faceOffset[0]) {
                    memcpy(&f_00[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_01[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_10[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_11[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    Scalar* o = out.getSampleDataRW(e, zero);
                    for(index_t i = 0; i < numComp; ++i) {
                        o[INDEX3(i,0,0,numComp,2)] = (f_10[i]-f_00[i])*cx[1][i] + (f_11[i]-f_01[i])*cx[0][i];
                        o[INDEX3(i,1,0,numComp,2)] = (f_01[i]-f_00[i])*cy[2][i];
                        o[INDEX3(i,0,1,numComp,2)] = (f_10[i]-f_00[i])*cx[0][i] + (f_11[i]-f_01[i])*cx[1][i];
                        o[INDEX3(i,1,1,numComp,2)] = (f_01[i]-f_00[i])*cy[2][i];
                    } // end of component loop i
                } // end of face 0
                if (quaddata->m_faceOffset[1]) {
                    memcpy(&f_00[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_01[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_10[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_11[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    Scalar* o = out.getSampleDataRW(e, zero);
                    for(index_t i = 0; i < numComp; ++i) {
                        o[INDEX3(i,0,0,numComp,2)] = (f_10[i]-f_00[i])*cx[1][i] + (f_11[i]-f_01[i])*cx[0][i];
                        o[INDEX3(i,1,0,numComp,2)] = (f_11[i]-f_10[i])*cy[2][i];
                        o[INDEX3(i,0,1,numComp,2)] = (f_10[i]-f_00[i])*cx[0][i] + (f_11[i]-f_01[i])*cx[1][i];
                        o[INDEX3(i,1,1,numComp,2)] = (f_11[i]-f_10[i])*cy[2][i];
                    } // end of component loop i
                } // end of face 1
                if (quaddata->m_faceOffset[2]) {
                    memcpy(&f_00[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_01[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_10[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_11[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    Scalar* o = out.getSampleDataRW(e, zero);
                    for(index_t i = 0; i < numComp; ++i) {
                        o[INDEX3(i,0,0,numComp,2)] = (f_10[i]-f_00[i])*cx[2][i];
                        o[INDEX3(i,1,0,numComp,2)] = (f_01[i]-f_00[i])*cy[1][i] + (f_11[i]-f_10[i])*cy[0][i];
                        o[INDEX3(i,0,1,numComp,2)] = (f_10[i]-f_00[i])*cx[2][i];
                        o[INDEX3(i,1,1,numComp,2)] = (f_01[i]-f_00[i])*cy[0][i] + (f_11[i]-f_10[i])*cy[1][i];
                    } // end of component loop i
                } // end of face 2
                if (quaddata->m_faceOffset[3]) {
                    memcpy(&f_00[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_01[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_10[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_11[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    Scalar* o = out.getSampleDataRW(e, zero);
                    for(index_t i = 0; i < numComp; ++i) {
                        o[INDEX3(i,0,0,numComp,2)] = (f_11[i]-f_01[i])*cx[2][i];
                        o[INDEX3(i,1,0,numComp,2)] = (f_01[i]-f_00[i])*cy[1][i] + (f_11[i]-f_10[i])*cy[0][i];
                        o[INDEX3(i,0,1,numComp,2)] = (f_11[i]-f_01[i])*cx[2][i];
                        o[INDEX3(i,1,1,numComp,2)] = (f_01[i]-f_00[i])*cy[0][i] + (f_11[i]-f_10[i])*cy[1][i];
                    } // end of component loop i
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
#pragma omp parallel for
            for(int q = 0; q < Q; ++q) // Loop over every quadrant within the tree
            {

                p4est_quadrant_t * quad = p4est_quadrant_array_index(tquadrants, q);
                double xy[3];
                p4est_qcoord_to_vertex(p4est->connectivity, t, quad->x, quad->y, xy);
                long e = getQuadID(NodeIDs.find(std::make_pair(xy[0],xy[1]))->second);
                quadrantData * quaddata = (quadrantData *) quad->p.user_data;

                std::vector<Scalar> f_00(numComp, zero);
                std::vector<Scalar> f_01(numComp, zero);
                std::vector<Scalar> f_10(numComp, zero);
                std::vector<Scalar> f_11(numComp, zero);

                if (quaddata->m_faceOffset[0]) {
                    memcpy(&f_00[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_01[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_10[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_11[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    Scalar* o = out.getSampleDataRW(e, zero);
                    for(index_t i = 0; i < numComp; ++i) {
                        o[INDEX3(i,0,0,numComp,2)] = (f_10[i] + f_11[i] - f_00[i] - f_01[i])*cx[2][i] * 0.5;
                        o[INDEX3(i,1,0,numComp,2)] = (f_01[i]-f_00[i])*cy[2][i];
                    } // end of component loop i
                } // end of face 0
                if (quaddata->m_faceOffset[1]) {
                    memcpy(&f_00[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_01[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_10[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_11[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    Scalar* o = out.getSampleDataRW(e, zero);
                    for(index_t i = 0; i < numComp; ++i) {
                        o[INDEX3(i,0,0,numComp,2)] = (f_10[i] + f_11[i] - f_00[i] - f_01[i])*cx[2][i] * 0.5;
                        o[INDEX3(i,1,0,numComp,2)] = (f_11[i]-f_10[i])*cy[2][i];
                    } // end of component loop i
                } // end of face 1
                if (quaddata->m_faceOffset[2]) {
                    memcpy(&f_00[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_01[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_10[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_11[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    Scalar* o = out.getSampleDataRW(e, zero);
                    for(index_t i = 0; i < numComp; ++i) {
                        o[INDEX3(i,0,0,numComp,2)] = (f_10[i]-f_00[i])*cx[2][i];
                        o[INDEX3(i,1,0,numComp,2)] = (f_01[i] + f_11[i] - f_00[i] - f_10[i])*cy[2][i] * 0.5;
                    } // end of component loop i
                } // end of face 2
                if (quaddata->m_faceOffset[3]) {
                    memcpy(&f_00[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_01[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_10[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_11[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    Scalar* o = out.getSampleDataRW(e, zero);
                    for(index_t i = 0; i < numComp; ++i) {
                        o[INDEX3(i,0,0,numComp,2)] = (f_11[i]-f_01[i])*cx[2][i];
                        o[INDEX3(i,1,0,numComp,2)] = (f_01[i] + f_11[i] - f_00[i] - f_10[i])*cy[2][i] * 0.5;
                    } // end of component loop i
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
//     const dim_t numComp = arg.getDataPointSize();
//     const index_t left = getFirstInDim(0);
//     const index_t bottom = getFirstInDim(1);
//     const int fs = arg.getFunctionSpace().getTypeCode();
//     const Scalar zero = static_cast<Scalar>(0);

//     bool HavePointData = arg.getFunctionSpace().getTypeCode() == Points;

// #ifdef ESYS_MPI
//     if(HavePointData && escript::getMPIRankWorld() == 0) {
// #else
//     if(HavePointData) {
// #endif
//         integrals[0] += arg.getNumberOfTaggedValues();
//     } else if (fs == Elements && arg.actsExpanded()) {
// #pragma omp parallel
//         {
//             std::vector<Scalar> int_local(numComp, zero);
//             const real_t w = m_dx[0]*m_dx[1]/4.;
// #pragma omp for nowait
//             for (index_t k1 = bottom; k1 < bottom+m_ownNE[1]; ++k1) {
//                 for (index_t k0 = left; k0 < left+m_ownNE[0]; ++k0) {
//                     const Scalar* f = arg.getSampleDataRO(INDEX2(k0, k1, m_NE[0]), zero);
//                     for (index_t i = 0; i < numComp; ++i) {
//                         const Scalar f0 = f[INDEX2(i,0,numComp)];
//                         const Scalar f1 = f[INDEX2(i,1,numComp)];
//                         const Scalar f2 = f[INDEX2(i,2,numComp)];
//                         const Scalar f3 = f[INDEX2(i,3,numComp)];
//                         int_local[i] += (f0+f1+f2+f3)*w;
//                     }  // end of component loop i
//                 } // end of k0 loop
//             } // end of k1 loop
// #pragma omp critical
//             for (index_t i=0; i<numComp; i++)
//                 integrals[i] += int_local[i];
//         } // end of parallel section

//     } else if (fs==ReducedElements || (fs==Elements && !arg.actsExpanded())) {
//         const real_t w = m_dx[0]*m_dx[1];
// #pragma omp parallel
//         {
//             std::vector<Scalar> int_local(numComp, 0);
// #pragma omp for nowait
//             for (index_t k1 = bottom; k1 < bottom+m_ownNE[1]; ++k1) {
//                 for (index_t k0 = left; k0 < left+m_ownNE[0]; ++k0) {
//                     const Scalar* f = arg.getSampleDataRO(INDEX2(k0, k1, m_NE[0]), zero);
//                     for (index_t i = 0; i < numComp; ++i) {
//                         int_local[i] += f[i]*w;
//                     }
//                 }
//             }
// #pragma omp critical
//             for (index_t i = 0; i < numComp; i++)
//                 integrals[i] += int_local[i];
//         } // end of parallel section

//     } else if (fs == FaceElements && arg.actsExpanded()) {
// #pragma omp parallel
//         {
//             std::vector<Scalar> int_local(numComp, zero);
//             const real_t w0 = m_dx[0]/2.;
//             const real_t w1 = m_dx[1]/2.;
//             if (m_faceOffset[0] > -1) {
// #pragma omp for nowait
//                 for (index_t k1 = bottom; k1 < bottom+m_ownNE[1]; ++k1) {
//                     const Scalar* f = arg.getSampleDataRO(m_faceOffset[0]+k1, zero);
//                     for (index_t i=0; i < numComp; ++i) {
//                         const Scalar f0 = f[INDEX2(i,0,numComp)];
//                         const Scalar f1 = f[INDEX2(i,1,numComp)];
//                         int_local[i] += (f0+f1)*w1;
//                     }  // end of component loop i
//                 } // end of k1 loop
//             }

//             if (m_faceOffset[1] > -1) {
// #pragma omp for nowait
//                 for (index_t k1 = bottom; k1 < bottom+m_ownNE[1]; ++k1) {
//                     const Scalar* f = arg.getSampleDataRO(m_faceOffset[1]+k1, zero);
//                     for (index_t i = 0; i < numComp; ++i) {
//                         const Scalar f0 = f[INDEX2(i,0,numComp)];
//                         const Scalar f1 = f[INDEX2(i,1,numComp)];
//                         int_local[i] += (f0+f1)*w1;
//                     }  // end of component loop i
//                 } // end of k1 loop
//             }

//             if (m_faceOffset[2] > -1) {
// #pragma omp for nowait
//                 for (index_t k0 = left; k0 < left+m_ownNE[0]; ++k0) {
//                     const Scalar* f = arg.getSampleDataRO(m_faceOffset[2]+k0, zero);
//                     for (index_t i = 0; i < numComp; ++i) {
//                         const Scalar f0 = f[INDEX2(i,0,numComp)];
//                         const Scalar f1 = f[INDEX2(i,1,numComp)];
//                         int_local[i] += (f0+f1)*w0;
//                     }  // end of component loop i
//                 } // end of k0 loop
//             }

//             if (m_faceOffset[3] > -1) {
// #pragma omp for nowait
//                 for (index_t k0 = left; k0 < left+m_ownNE[0]; ++k0) {
//                     const Scalar* f = arg.getSampleDataRO(m_faceOffset[3]+k0, zero);
//                     for (index_t i = 0; i < numComp; ++i) {
//                         const Scalar f0 = f[INDEX2(i,0,numComp)];
//                         const Scalar f1 = f[INDEX2(i,1,numComp)];
//                         int_local[i] += (f0+f1)*w0;
//                     }  // end of component loop i
//                 } // end of k0 loop
//             }
// #pragma omp critical
//             for (index_t i = 0; i < numComp; i++)
//                 integrals[i] += int_local[i];
//         } // end of parallel section

//     } else if (fs==ReducedFaceElements || (fs==FaceElements && !arg.actsExpanded())) {
// #pragma omp parallel
//         {
//             std::vector<Scalar> int_local(numComp, 0);
//             if (m_faceOffset[0] > -1) {
// #pragma omp for nowait
//                 for (index_t k1 = bottom; k1 < bottom+m_ownNE[1]; ++k1) {
//                     const Scalar* f = arg.getSampleDataRO(m_faceOffset[0]+k1, zero);
//                     for (index_t i = 0; i < numComp; ++i) {
//                         int_local[i] += f[i]*m_dx[1];
//                     }
//                 }
//             }

//             if (m_faceOffset[1] > -1) {
// #pragma omp for nowait
//                 for (index_t k1 = bottom; k1 < bottom+m_ownNE[1]; ++k1) {
//                     const Scalar* f = arg.getSampleDataRO(m_faceOffset[1]+k1, zero);
//                     for (index_t i = 0; i < numComp; ++i) {
//                         int_local[i] += f[i]*m_dx[1];
//                     }
//                 }
//             }

//             if (m_faceOffset[2] > -1) {
// #pragma omp for nowait
//                 for (index_t k0 = left; k0 < left+m_ownNE[0]; ++k0) {
//                     const Scalar* f = arg.getSampleDataRO(m_faceOffset[2]+k0, zero);
//                     for (index_t i = 0; i < numComp; ++i) {
//                         int_local[i] += f[i]*m_dx[0];
//                     }
//                 }
//             }

//             if (m_faceOffset[3] > -1) {
// #pragma omp for nowait
//                 for (index_t k0 = left; k0 < left+m_ownNE[0]; ++k0) {
//                     const Scalar* f = arg.getSampleDataRO(m_faceOffset[3]+k0, zero);
//                     for (index_t i = 0; i < numComp; ++i) {
//                         int_local[i] += f[i]*m_dx[0];
//                     }
//                 }
//             }

// #pragma omp critical
//             for (index_t i = 0; i < numComp; i++)
//                 integrals[i] += int_local[i];
//         } // end of parallel section
//     } // function space selector
}

//protected
void Rectangle::nodesToDOF(escript::Data& out, const escript::Data& in) const
{
    //TODO
    throw OxleyException("nodesToDOF");

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

// Copies the solution information to the mesh
void Rectangle::updateSolutionInformation(escript::Data solution)
{
    for(p4est_topidx_t treeid = p4est->first_local_tree; treeid <= p4est->last_local_tree; ++treeid) 
    {
        p4est_tree_t * tree = p4est_tree_array_index(p4est->trees, treeid);
        sc_array_t * tquadrants = &tree->quadrants;
        p4est_locidx_t Q = (p4est_locidx_t) tquadrants->elem_count;
// #pragma omp parallel for
        for(int q = 0; q < Q; ++q) // Loop over all quadrants
        { 
            p4est_quadrant_t * quad = p4est_quadrant_array_index(tquadrants, q);
            p4est_qcoord_t length = P4EST_QUADRANT_LEN(quad->level);
            for(int n = 0; n < 4; n++)
            {
                double xy[3];
                long lx[4] = {0,length,0,length};
                long ly[4] = {0,0,length,length};
                long lni[4] = {-1};
                for(int i = 0; i < 4; i++)
                {
                    p4est_qcoord_to_vertex(p4est->connectivity, treeid, quad->x+lx[i], quad->y+ly[i], xy);
                    lni[i] = NodeIDs.find(std::make_pair(xy[0],xy[1]))->second;
                    auto temp = solution.getSampleDataRO(i);
                    current_solution[lni[i]]=*temp;
                }                
            }
        }
    }
}

void Rectangle::updateMeshInformation()
{
    refineMesh("MARE2DEM");
}

// Returns the solution information currently stored in Oxley

escript::Data Rectangle::getUpdatedSolution()
{
    return escript::Data(0); //TODO
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
    for (int i = 0; i < tags.size(); i++) {
        dim_t node = findNode(&coords[i * m_numDim]);
        if (node >= 0) {
            m_diracPointNodeIDs.push_back(borrowSampleReferenceIDs(Nodes)[node]);
            DiracPoint dp;
            dp.node = node; //local
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
    escript::randomFillArray(seed, src, ext[0]*ext[1]*numvals);

#ifdef ESYS_MPI
    if ((internal[0] < 5) || (internal[1] < 5)) {
        // since the dimensions are equal for all ranks, this exception
        // will be thrown on all ranks
        throw OxleyException("Random Data in Oxley requires at least five elements per side per rank.");
    }
    dim_t X = m_mpiInfo->rank%m_NX[0];
    dim_t Y = m_mpiInfo->rank/m_NX[0];
#endif

#ifdef ESYS_MPI
    BlockGrid2 grid(m_NX[0]-1, m_NX[1]-1);
    // it's +2 not +1 because a whole element is shared (and hence there is
    // an overlap of two points both of which need to have "radius" points on
    // either side.
    size_t inset=2*radius+2;

    // how wide is the x-dimension between the two insets
    size_t xmidlen=ext[0]-2*inset;
    size_t ymidlen=ext[1]-2*inset;

    Block2 block(ext[0], ext[1], inset, xmidlen, ymidlen, numvals);

    // a non-tight upper bound on how many we need
    MPI_Request reqs[40];
    MPI_Status stats[40];
    short rused=0;

    messvec incoms;
    messvec outcoms;

    grid.generateInNeighbours(X, Y, incoms);
    grid.generateOutNeighbours(X, Y, outcoms);

    block.copyAllToBuffer(src);

    int comserr = 0;
    for (size_t i=0; i < incoms.size(); ++i) {
        message& m = incoms[i];
        comserr |= MPI_Irecv(block.getInBuffer(m.destbuffid),
                             block.getBuffSize(m.destbuffid), MPI_DOUBLE,
                             m.sourceID, m.tag, m_mpiInfo->comm,
                             reqs+(rused++));
        block.setUsed(m.destbuffid);
    }

    for (size_t i=0; i < outcoms.size(); ++i) {
        message& m = outcoms[i];
        comserr |= MPI_Isend(block.getOutBuffer(m.srcbuffid),
                             block.getBuffSize(m.srcbuffid), MPI_DOUBLE,
                             m.destID, m.tag, m_mpiInfo->comm, reqs+(rused++));
    }

    if (!comserr) {
        comserr = MPI_Waitall(rused, reqs, stats);
    }

    if (comserr) {
        // Yes this is throwing an exception as a result of an MPI error
        // and no we don't inform the other ranks that we are doing this.
        // However, we have no reason to believe coms work at this point anyway
        throw RipleyException("Error in coms for randomFill");
    }

    block.copyUsedFromBuffer(src);
#endif

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
    // Check to see if the node is in the map
    if(NodeIDs.count(std::make_pair(coords[0],coords[1]))==1)
        return NodeIDs.find(std::make_pair(coords[0],coords[1]))->second;

    // TODO speed enhancements
    // Otherwise find the nearest element
    double sq_distance = m_NX[0]*m_NX[1];
    long closest = 0;

    for(std::pair<DoublePair,long> e : NodeIDs)
        if(e.first.first*e.first.first+e.first.second*e.first.second < sq_distance)
            closest = e.second;

    return closest;

}

// instantiate our two supported versions
template
void Rectangle::assembleGradientImpl<real_t>(escript::Data& out,
                                             const escript::Data& in) const;

template
void Rectangle::assembleGradientImpl<cplx_t>(escript::Data& out,
                                             const escript::Data& in) const;


} // end of namespace oxley
