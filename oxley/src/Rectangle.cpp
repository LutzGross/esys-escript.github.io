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

#include <ctime>
#include <random>
#include <vector>

#include <escript/Data.h>
#include <escript/DataFactory.h>
#include <escript/FunctionSpaceFactory.h>

#include <oxley/AbstractAssembler.h>
#include <oxley/DefaultAssembler2D.h>
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
    // // const int8_t tree_to_face[P4EST_FACES] = {1, 0, 3, 2,}; //aeae todo: add in periodic boundary conditions
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

    // connectivity = p4est_connectivity_new_brick(n0, n1, false, false); //AEAE this is temporary
    connectivity = new_rectangle_connectivity(n0, n1, false, false, x0, x1, y0, y1);

    if(!p4est_connectivity_is_valid(connectivity))
        throw OxleyException("Could not create a valid connectivity.");

    // Allocate some memory
    forestData = new p4estData;
    forestData = (p4estData *) malloc(sizeof(p4estData));

    // Create a p4est
    p4est_locidx_t min_quadrants = n0*n1;
    int min_level = 0;
    int fill_uniform = 1;
    p4est = p4est_new_ext(m_mpiInfo->comm, connectivity, min_quadrants,
            min_level, fill_uniform, sizeof(p4estData), init_rectangle_data, (void *) &forestData);

    // Nodes numbering
    p4est_ghost_t * ghost = p4est_ghost_new(p4est, P4EST_CONNECT_FULL);
    nodes = p4est_lnodes_new(p4est, ghost, 1);
    p4est_ghost_destroy(ghost);

    // Record the physical dimensions of the domain and the location of the origin
    forestData->m_origin[0] = x0;
    forestData->m_origin[1] = y0;
    forestData->m_lxy[0] = x1;
    forestData->m_lxy[1] = y1;
    forestData->m_length[0] = x1-x0;
    forestData->m_length[1] = y1-y0;

    // Whether or not we have periodic boundaries
    forestData->periodic[0] = periodic0;
    forestData->periodic[1] = periodic1;

    // Find the grid spacing for each level of refinement in the mesh
#pragma omp parallel for
    for(int i = 0; i<=P4EST_MAXLEVEL; i++){
        double numberOfSubDivisions = (p4est_qcoord_t) (1 << (P4EST_MAXLEVEL - i));
        forestData->m_dx[0][i] = forestData->m_length[0] / numberOfSubDivisions;
        forestData->m_dx[1][i] = forestData->m_length[1] / numberOfSubDivisions;
    }

    // max levels of refinement
    forestData->max_levels_refinement = MAXREFINEMENTLEVELS;
    
    // element order
    m_order = order;

    // initial tag
    tags[0] = 0;
    numberOfTags=1;

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
    // populateDofMap();

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
    
}

/**
   \brief
   Destructor.
*/
Rectangle::~Rectangle(){
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
   dumps the mesh to a file with the given name
   \param filename The name of the output file
*/
void Rectangle::dump(const std::string& filename) const
{
    throw OxleyException("dump: not supported");
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

bool Rectangle::probeInterpolationAcross(int fsType_source,
        const escript::AbstractDomain& domain, int fsType_target) const
{
    throw OxleyException("probeInterpolationAcross"); //AE this is temporary
}

void Rectangle::interpolateAcross(escript::Data& target, const escript::Data& source) const
{
    throw OxleyException("interpolateAcross"); //AE this is temporary
}

void Rectangle::setToNormal(escript::Data& out) const
{
    if (out.getFunctionSpace().getTypeCode() == FaceElements) {
        out.requireWrite();
        for (p4est_topidx_t t = p4est->first_local_tree; t <= p4est->last_local_tree; t++) 
        {
            p4est_tree_t * currenttree = p4est_tree_array_index(p4est->trees, t);
            sc_array_t * tquadrants = &currenttree->quadrants;
            p4est_locidx_t Q = (p4est_locidx_t) tquadrants->elem_count;
#pragma omp parallel for
            for (p4est_locidx_t e = nodes->global_offset; e < Q+nodes->global_offset; e++)
            {
                // Work out what level this element is on 
                p4est_quadrant_t * quad = p4est_quadrant_array_index(tquadrants, e);
                quadrantData * quaddata = (quadrantData *) quad->p.user_data;

                if (quaddata->m_faceOffset == 0) {
                    double* o = out.getSampleDataRW(e);
                    // set vector at two quadrature points
                    *o++ = -1.;
                    *o++ = 0.;
                    *o++ = -1.;
                    *o = 0.;
                }

                if (quaddata->m_faceOffset == 1) {
                    double* o = out.getSampleDataRW(e);
                    // set vector at two quadrature points
                    *o++ = 1.;
                    *o++ = 0.;
                    *o++ = 1.;
                    *o = 0.;
                }

                if (quaddata->m_faceOffset == 2) {
                    double* o = out.getSampleDataRW(e);
                    // set vector at two quadrature points
                    *o++ = 0.;
                    *o++ = -1.;
                    *o++ = 0.;
                    *o = -1.;
                }

                if (quaddata->m_faceOffset == 3) {
                    double* o = out.getSampleDataRW(e);
                    // set vector at two quadrature points
                    *o++ = 0.;
                    *o++ = 1.;
                    *o++ = 0.;
                    *o = 1.;
                }
            }
        }
    } else if (out.getFunctionSpace().getTypeCode() == ReducedFaceElements) {
        out.requireWrite();
        for (p4est_topidx_t t = p4est->first_local_tree; t <= p4est->last_local_tree; t++) 
        {
            p4est_tree_t * currenttree = p4est_tree_array_index(p4est->trees, t);
            sc_array_t * tquadrants = &currenttree->quadrants;
            p4est_locidx_t Q = (p4est_locidx_t) tquadrants->elem_count;
#pragma omp parallel for
            for (p4est_locidx_t e = nodes->global_offset; e < Q+nodes->global_offset; e++)
            {
                // Work out what level this element is on 
                p4est_quadrant_t * quad = p4est_quadrant_array_index(tquadrants, e);
                quadrantData * quaddata = (quadrantData *) quad->p.user_data;

                if (quaddata->m_faceOffset == 0) {
                    double* o = out.getSampleDataRW(e);
                    *o++ = -1.;
                    *o = 0.;
                }

                if (quaddata->m_faceOffset == 1) {
                    double* o = out.getSampleDataRW(e);
                    *o++ = 1.;
                    *o = 0.;
                }

                if (quaddata->m_faceOffset == 2) {
                    double* o = out.getSampleDataRW(e);
                    *o++ = 0.;
                    *o = -1.;
                }

                if (quaddata->m_faceOffset == 3) {
                    double* o = out.getSampleDataRW(e);
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
}

void Rectangle::setToSize(escript::Data& out) const
{
    if (out.getFunctionSpace().getTypeCode() == Elements
        || out.getFunctionSpace().getTypeCode() == ReducedElements)
    {
        out.requireWrite();
        const dim_t numQuad = out.getNumDataPointsPerSample();
        const dim_t numElements = getNumElements();

        for (p4est_topidx_t t = p4est->first_local_tree; t <= p4est->last_local_tree; t++) 
        {
            p4est_tree_t * currenttree = p4est_tree_array_index(p4est->trees, t);
            sc_array_t * tquadrants = &currenttree->quadrants;
            p4est_locidx_t Q = (p4est_locidx_t) tquadrants->elem_count;
        #pragma omp parallel for
            for (p4est_locidx_t e = nodes->global_offset; e < Q+nodes->global_offset; e++)
            {
                p4est_quadrant_t * quad = p4est_quadrant_array_index(tquadrants, e);
                // int l = quad->level;
                const double size = sqrt(forestData->m_dx[0][quad->level]*forestData->m_dx[0][quad->level]
                                        +forestData->m_dx[1][quad->level]*forestData->m_dx[1][quad->level]);
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

        for (p4est_topidx_t t = p4est->first_local_tree; t <= p4est->last_local_tree; t++) 
        {
            p4est_tree_t * currenttree = p4est_tree_array_index(p4est->trees, t);
            sc_array_t * tquadrants = &currenttree->quadrants;
            p4est_locidx_t Q = (p4est_locidx_t) tquadrants->elem_count;
        #pragma omp parallel for
            for (p4est_locidx_t e = nodes->global_offset; e < Q+nodes->global_offset; e++)
            {
                // Work out what level this element is on 
                p4est_quadrant_t * quad = p4est_quadrant_array_index(tquadrants, e);
                quadrantData * quaddata = (quadrantData *) quad->p.user_data;

                if (quaddata->m_faceOffset == 0) {
                    double* o = out.getSampleDataRW(e);
                    std::fill(o, o+numQuad, forestData->m_dx[1][quad->level]);
                }
                else if (quaddata->m_faceOffset == 1) {
                    double* o = out.getSampleDataRW(e);
                    std::fill(o, o+numQuad, forestData->m_dx[1][quad->level]);
                }
                else if (quaddata->m_faceOffset == 2) {
                    double* o = out.getSampleDataRW(e);
                    std::fill(o, o+numQuad, forestData->m_dx[0][quad->level]);
                }
                else if (quaddata->m_faceOffset == 3) {
                    double* o = out.getSampleDataRW(e);
                    std::fill(o, o+numQuad, forestData->m_dx[0][quad->level]);
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
    throw OxleyException("ownSample"); //AE this is temporary
}


dim_t Rectangle::getNumDataPointsGlobal() const
{
    return getNumNodes();
}

/* This is a wrapper for filtered (and non-filtered) randoms
 * For detailed doco see randomFillWorker
*/
escript::Data Rectangle::randomFill(const escript::DataTypes::ShapeType& shape,
                                    const escript::FunctionSpace& fs,
                                    long seed, const bp::tuple& filter) const
{
    throw OxleyException("randomFill"); //AE this is temporary
}


const dim_t* Rectangle::borrowSampleReferenceIDs(int fsType) const
{
    throw OxleyException("borrowSampleReferenceIDs"); 
    
    switch (fsType) {
        case Nodes:
        case ReducedNodes:
            return &myColumns[0];
        case DegreesOfFreedom:
        case ReducedDegreesOfFreedom:
            throw OxleyException("Unknown Error.");
        case Elements:
        case ReducedElements:
            throw OxleyException("borrowSampleReferenceIDs: Elements");
            // return &m_elementId[0];
        case FaceElements:
        case ReducedFaceElements:
            throw OxleyException("borrowSampleReferenceIDs: FaceIDs");
            // return &m_faceId[0];
        case Points:
            throw OxleyException("borrowSampleReferenceIDs: m_diracPointNodeIDs");
            // return &m_diracPointNodeIDs[0];
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
#ifdef P4EST_ENABLE_DEBUG
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

void Rectangle::refineMesh(int maxRecursion, std::string algorithmname)
{
    if(!algorithmname.compare("uniform"))
    {
        if(maxRecursion < 0)
            throw OxleyException("Invalid value for maxRecursion");

        if(maxRecursion == 0){
            // This may be unwise as it consumes very large amounts of memory
            // This will perform P4EST_QMAXLEVEL levels of refinement which is defined in p4est.h
            // by default P4EST_QMAXLEVEL = 29
            p4est_refine_ext(p4est, true, -1, refine_uniform, init_rectangle_data, refine_copy_parent_quadrant);
        } else {
            p4est_refine_ext(p4est, true, maxRecursion, refine_uniform, init_rectangle_data, refine_copy_parent_quadrant);
        }
        p4est_balance_ext(p4est, P4EST_CONNECT_FULL, init_rectangle_data, refine_copy_parent_quadrant);
    } 
    else if(!algorithmname.compare("random"))
    {
        if(maxRecursion <= 0)
            throw OxleyException("Invalid value for maxRecursion");
        p4est_refine_ext(p4est, true, maxRecursion, random_refine, init_rectangle_data, refine_copy_parent_quadrant);
        p4est_balance_ext(p4est, P4EST_CONNECT_FULL, init_rectangle_data, refine_copy_parent_quadrant);
    }
    else {
        throw OxleyException("Unknown refinement algorithm name.");
    }

    // Make sure that nothing went wrong
#ifdef P4EST_ENABLE_DEBUG
    if(!p4est_is_valid(p4est))
        throw OxleyException("p4est broke during refinement");
    if(!p4est_connectivity_is_valid(connectivity))
        throw OxleyException("connectivity broke during refinement");
#endif

    bool partition_for_coarsening = true;
    p4est_partition_ext(p4est, partition_for_coarsening, NULL);
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
    return (xy[0] == 0) || (xy[0] == forestData->m_lxy[0]) || (xy[1] == 0) || (xy[1] == forestData->m_lxy[1]);
}

// returns True for a boundary node on the north or east of the domain
bool Rectangle::isUpperBoundaryNode(p4est_quadrant_t * quad, int n, p4est_topidx_t treeid, p4est_qcoord_t length) const
{
    double lx = length * ((int) (n % 2) == 1);
    double ly = length * ((int) (n / 2) == 1);
    double xy[3];
    p4est_qcoord_to_vertex(p4est->connectivity, treeid, quad->x+lx, quad->y+ly, xy);
    return (xy[0] == forestData->m_lxy[0]) || (xy[1] == forestData->m_lxy[1]);
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
        int8_t ishanging1 = face_code >> 3;

        int8_t f0 = p4est_corner_faces[c][0];
        int8_t f1 = p4est_corner_faces[c][1];

        int d0 = c / 2 == 0;
        int d1 = c % 2 == 0;
        
        return ((ishanging0 == 1) && (p4est_face_corners[f0][d0] == n)) 
            || ((ishanging1 == 1) && (p4est_face_corners[f1][d1] == n));
    }
}

//protected
void Rectangle::updateNodeIncrements()
{
    nodeIncrements[0] = 1;
    for(p4est_topidx_t treeid = p4est->first_local_tree+1, k=1; treeid <= p4est->last_local_tree; ++treeid, ++k) 
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
#ifdef P4EST_ENABLE_DEBUG
    std::cout << "Renumbering nodes... " << std::endl;
#endif

    NodeIDs.clear();
#pragma omp for
    for(p4est_topidx_t treeid = p4est->first_local_tree; treeid <= p4est->last_local_tree; ++treeid) {
        p4est_tree_t * tree = p4est_tree_array_index(p4est->trees, treeid);
        sc_array_t * tquadrants = &tree->quadrants;
        p4est_locidx_t Q = (p4est_locidx_t) tquadrants->elem_count;
        for(int q = 0; q < Q; ++q) { 
            p4est_quadrant_t * quad = p4est_quadrant_array_index(tquadrants, q);
            p4est_qcoord_t length = P4EST_QUADRANT_LEN(quad->level);
            int k = q - Q + nodeIncrements[treeid - p4est->first_local_tree];

            // If there are no hanging nodes here, skip to the next quadrant
            if(nodes->face_code[k] != 0)
            {
                // Decode the info
                int8_t face_code = nodes->face_code[k];
                int8_t c = face_code & 0x03;
                int8_t ishanging0 = (face_code >> 2) & 0x01;
                int8_t ishanging1 = face_code >> 3;

                int8_t f0 = p4est_corner_faces[c][0];
                int8_t f1 = p4est_corner_faces[c][1];

                int d0 = c / 2 == 0;
                int d1 = c % 2 == 0;

                int tmp0 = p4est_face_corners[f0][d0];
                int tmp1 = p4est_face_corners[f1][d1];
                
                // Record which nodes are hanging
                if(ishanging0 || ishanging1)
                {
                    p4est_quadrant_t * quad = p4est_quadrant_array_index(tquadrants, q);
                    p4est_qcoord_t length = P4EST_QUADRANT_LEN(quad->level);
                    double xy[3];

                    if(ishanging0)
                    {                   
                        double lx = length * ((int) (tmp0 % 2) == 1);
                        double ly = length * ((int) (tmp0 / 2) == 1);
                        p4est_qcoord_to_vertex(p4est->connectivity, treeid, quad->x+lx, quad->y+ly, xy);
                        if(NodeIDs.count(std::make_pair(xy[0],xy[1]))==0)
                        {
#ifdef P4EST_ENABLE_DEBUG
                            std::cout << NodeIDs.size() << ": " << xy[0] << ", " << xy[1] << std::endl;
#endif
                            NodeIDs[std::make_pair(xy[0],xy[1])]=NodeIDs.size();
                        }
                    }
                    if(ishanging1)
                    {
                        double lx = length * ((int) (tmp1 % 2) == 1);
                        double ly = length * ((int) (tmp1 / 2) == 1);
                        p4est_qcoord_to_vertex(p4est->connectivity, treeid, quad->x+lx, quad->y+ly, xy);
                        if(NodeIDs.count(std::make_pair(xy[0],xy[1]))==0)
                        {
#ifdef P4EST_ENABLE_DEBUG
                            std::cout << NodeIDs.size() << ": " << xy[0] << ", " << xy[1] << std::endl;
#endif
                            NodeIDs[std::make_pair(xy[0],xy[1])]=NodeIDs.size();
                        }
                    }
                }
            }

            for(int n = 0; n < 4; n++)
            {
                double lx = length * ((int) (n % 2) == 1);
                double ly = length * ((int) (n / 2) == 1);
                double xy[3];
                p4est_qcoord_to_vertex(p4est->connectivity, treeid, quad->x+lx, quad->y+ly, xy);
                if( (n == 0) 
                    || ((n == 1) && (xy[0] == forestData->m_lxy[0]))
                    || ((n == 2) && (xy[1] == forestData->m_lxy[1]))
                    || ((n == 3) && (xy[0] == forestData->m_lxy[0]) && (xy[1] == forestData->m_lxy[1]))
                  )
                {
                    if(NodeIDs.count(std::make_pair(xy[0],xy[1]))==0)
                    {
#ifdef P4EST_ENABLE_DEBUG
                        std::cout << NodeIDs.size() << ": " << xy[0] << ", " << xy[1] << std::endl;
#endif
                        NodeIDs[std::make_pair(xy[0],xy[1])]=NodeIDs.size();
                    }
                }
            }
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

    for(p4est_topidx_t treeid = p4est->first_local_tree; treeid <= p4est->last_local_tree; ++treeid) {
        p4est_tree_t * tree = p4est_tree_array_index(p4est->trees, treeid);
        sc_array_t * tquadrants = &tree->quadrants;
        p4est_locidx_t Q = (p4est_locidx_t) tquadrants->elem_count;

#pragma omp parallel for
        for (int q = 0; q < Q; ++q) { // Loop over the elements attached to the tree
            p4est_quadrant_t * quad = p4est_quadrant_array_index(tquadrants, q);
            p4est_qcoord_t length = P4EST_QUADRANT_LEN(quad->level);

            // Loop over the four corners of the quadrant
            for(int n = 0; n < 4; ++n){
                int k = q - Q + nodeIncrements[treeid - p4est->first_local_tree];
                long lidx = nodes->element_nodes[4*k+n];

                if(lidx < nodes->owned_count) // if this process owns the node
                {

                    // if(isHangingNode(nodes->face_code[k], n))
                    // {
                    //     double lx = length * ((int) (n % 2) == 1);
                    //     double ly = length * ((int) (n / 2) == 1);
                    //     double xy[3];
                    //     p4est_qcoord_to_vertex(p4est->connectivity, treeid, quad->x+lx, quad->y+ly, xy);

                    //     std::pair<double,double> coords = std::make_pair(xy[0],xy[1]);
                    //     long lni = hangingNodeIDs.find(std::make_pair(xy[0],xy[1]))->second;
                    //     double * point = arg.getSampleDataRW(lni);
                    //     point[0] = xy[0];
                    //     point[1] = xy[1];

                    //     continue;
                    // }

                    double lx = length * ((int) (n % 2) == 1);
                    double ly = length * ((int) (n / 2) == 1);
                    double xy[3];
                    p4est_qcoord_to_vertex(p4est->connectivity, treeid, quad->x+lx, quad->y+ly, xy);

                    if(     (n == 0) 
                        || ((n == 1) && (xy[0] == forestData->m_lxy[0]))
                        || ((n == 2) && (xy[1] == forestData->m_lxy[1]))
                        || ((n == 3) && (xy[0] == forestData->m_lxy[0]) && (xy[1] == forestData->m_lxy[1]))
                      )
                    {
                        long lni = NodeIDs.find(std::make_pair(xy[0],xy[1]))->second;
                        double * point = arg.getSampleDataRW(lni);
                        point[0] = xy[0];
                        point[1] = xy[1];
                        // std::cout << lni << ": (" << xy[0] << ", " << xy[1] << ")" << std::endl;
                    }
                }
            }
        }
    }
}

// //private
// void Rectangle::populateDofMap()
// {
//     m_dofMap.assign(getNumNodes(), 0);
// #pragma omp parallel for
//     for(index_t i = 0; i < getNumNodes()-1; i++)
//     {
//         m_dofMap[i] = myRows[i+1]-myRows[i];
//     }

// // #ifdef ESYS_HAVE_PASO
// //     createPasoConnector(neighbour, offsetInShared, offsetInShared, sendShared, recvShared);
// // #endif
// }


//private
template<typename Scalar>
void Rectangle::addToMatrixAndRHS(escript::AbstractSystemMatrix* S, escript::Data& F,
         const std::vector<Scalar>& EM_S, const std::vector<Scalar>& EM_F, 
         bool addS, bool addF, index_t firstNode, int nEq, int nComp) const
{    
    IndexVector rowIndex(4);
    // rowIndex[0] = m_dofMap[firstNode];
    // rowIndex[1] = m_dofMap[firstNode+1];
    // rowIndex[2] = m_dofMap[firstNode+m_NN[0]];
    // rowIndex[3] = m_dofMap[firstNode+m_NN[0]+1];

    if(addF)
    {
        Scalar* F_p = F.getSampleDataRW(0, static_cast<Scalar>(0));
        for(p4est_topidx_t treeid = p4est->first_local_tree; treeid <= p4est->last_local_tree; treeid++)
        {
            p4est_tree_t * currenttree = p4est_tree_array_index(p4est->trees, treeid);
            sc_array_t * tquadrants = &currenttree->quadrants;
            p4est_locidx_t Q = (p4est_locidx_t) tquadrants->elem_count;
#pragma omp parallel for
            for (p4est_locidx_t e = nodes->global_offset; e < Q+nodes->global_offset; e++)
                F_p[e]+=EM_F[e];
        }
    }
    if(addS)
    {
        addToSystemMatrix<Scalar>(S, rowIndex, nEq, EM_S);
    }
}

template
void Rectangle::addToMatrixAndRHS<real_t>(escript::AbstractSystemMatrix* S, escript::Data& F,
         const std::vector<real_t>& EM_S, const std::vector<real_t>& EM_F, 
         bool addS, bool addF, index_t firstNode, int nEq, int nComp) const;

template
void Rectangle::addToMatrixAndRHS<cplx_t>(escript::AbstractSystemMatrix* S, escript::Data& F,
         const std::vector<cplx_t>& EM_S, const std::vector<cplx_t>& EM_F, 
         bool addS, bool addF, index_t firstNode, int nEq, int nComp) const;

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
    //todo:
    throw OxleyException("interpolateNodesOnFaces");

    // if (in.isComplex()!=out.isComplex())
    // {
    //     throw OxleyException("Programmer Error: in and out parameters do not have the same complexity.");
    // }
    // if (in.isComplex())
    // {
    //     interpolateNodesOnFacesWorker(out, in, reduced, escript::DataTypes::cplx_t(0));
    // }
    // else
    // {
    //     interpolateNodesOnFacesWorker(out, in, reduced, escript::DataTypes::real_t(0));      
    // }
}


// private   
template <typename S> 
void Rectangle::interpolateNodesOnElementsWorker(escript::Data& out,
                                           const escript::Data& in,
                                           bool reduced, S sentinel) const
{
    const dim_t numComp = in.getDataPointSize();
    const long  numNodes = getNumNodes();
    if (reduced) {
        out.requireWrite();
        const S c0 = 0.25;
        double * fxx = new double[4*numComp*numNodes];

        // This structure is used to store info needed by p4est
        interpolateNodesOnElementsWorker_Data<S> interpolateData;
        interpolateData.fxx = fxx;
        interpolateData.sentinel = sentinel;
        interpolateData.offset = numComp*sizeof(S);

        p4est_iterate(p4est, NULL, &interpolateData, get_interpolateNodesOnElementWorker_data, NULL, NULL);
        for(p4est_topidx_t treeid = p4est->first_local_tree; treeid <= p4est->last_local_tree; treeid++)
        {
            p4est_tree_t * currenttree = p4est_tree_array_index(p4est->trees, treeid);
            sc_array_t * tquadrants = &currenttree->quadrants;
            p4est_locidx_t Q = (p4est_locidx_t) tquadrants->elem_count;
#pragma omp parallel for
            for (p4est_locidx_t e = nodes->global_offset; e < Q+nodes->global_offset; e++)
            {
                S* o = out.getSampleDataRW(e, sentinel);
                for (index_t i=0; i < numComp; ++i) 
                {
                    o[i] = c0*(fxx[INDEX3(0,i,e,numComp,numNodes)] + 
                               fxx[INDEX3(1,i,e,numComp,numNodes)] + 
                               fxx[INDEX3(2,i,e,numComp,numNodes)] + 
                               fxx[INDEX3(3,i,e,numComp,numNodes)]);
                }
            }
        }
        delete[] fxx;
    } else {
        out.requireWrite();
        const S c0 = 0.16666666666666666667;
        const S c1 = 0.044658198738520451079;
        const S c2 = 0.62200846792814621559;
        double * fxx = new double[4*numComp*numNodes];
        // This structure is used to store info needed by p4est
        interpolateNodesOnElementsWorker_Data<S> interpolateData;
        interpolateData.fxx = fxx;
        interpolateData.sentinel = sentinel;
        interpolateData.offset = numComp*sizeof(S);
        p4est_iterate(p4est, NULL, &interpolateData, get_interpolateNodesOnElementWorker_data, NULL, NULL);
        for(p4est_topidx_t treeid = p4est->first_local_tree; treeid <= p4est->last_local_tree; treeid++)
        {
            p4est_tree_t * currenttree = p4est_tree_array_index(p4est->trees, treeid);
            sc_array_t * tquadrants = &currenttree->quadrants;
            p4est_locidx_t Q = (p4est_locidx_t) tquadrants->elem_count;
#pragma omp parallel for
            for (p4est_locidx_t e = nodes->global_offset; e < Q+nodes->global_offset; e++)
            {
                S* o = out.getSampleDataRW(e, sentinel);
                for (index_t i=0; i < numComp; ++i) {
                    o[INDEX2(i,numComp,0)] = c0*(fxx[INDEX3(1,i,e,numComp,numNodes)] +    fxx[INDEX3(2,i,e,numComp,numNodes)]) + 
                                             c1* fxx[INDEX3(3,i,e,numComp,numNodes)] + c2*fxx[INDEX3(0,i,e,numComp,numNodes)];
                    o[INDEX2(i,numComp,1)] = c0*(fxx[INDEX3(2,i,e,numComp,numNodes)] +    fxx[INDEX3(3,i,e,numComp,numNodes)]) + 
                                             c1* fxx[INDEX3(1,i,e,numComp,numNodes)] + c2*fxx[INDEX3(2,i,e,numComp,numNodes)];
                    o[INDEX2(i,numComp,2)] = c0*(fxx[INDEX3(0,i,e,numComp,numNodes)] +    fxx[INDEX3(3,i,e,numComp,numNodes)]) + 
                                             c1* fxx[INDEX3(2,i,e,numComp,numNodes)] + c2*fxx[INDEX3(1,i,e,numComp,numNodes)];
                    o[INDEX2(i,numComp,3)] = c0*(fxx[INDEX3(1,i,e,numComp,numNodes)] +    fxx[INDEX3(2,i,e,numComp,numNodes)]) + 
                                             c1* fxx[INDEX3(0,i,e,numComp,numNodes)] + c2*fxx[INDEX3(3,i,e,numComp,numNodes)];
                }
            }
        }
        delete[] fxx;
    }
}

//private
template <typename S>
void Rectangle::interpolateNodesOnFacesWorker(escript::Data& out,
                                        const escript::Data& in,
                                        bool reduced, S sentinel) const
{
    //todo:
    throw OxleyException("interpolateNodesOnFacesWorker");

//     const dim_t numComp = in.getDataPointSize();
//     const dim_t numNodes = getNumNodes();
//     if (reduced) {
//         out.requireWrite();
//         double * fxx = new double[4*numComp*numNodes];
//         // This structure is used to store info needed by p4est
//         interpolateNodesOnFacesWorker_Data<S> interpolateData;
//         interpolateData.fxx = fxx;
//         interpolateData.sentinel = sentinel;
//         interpolateData.offset = numComp*sizeof(S);
//         for(p4est_topidx_t treeid = p4est->first_local_tree; treeid <= p4est->last_local_tree; treeid++)
//         {
//             p4est_tree_t * currenttree = p4est_tree_array_index(p4est->trees, treeid);
//             sc_array_t * tquadrants = &currenttree->quadrants;
//             p4est_locidx_t Q = (p4est_locidx_t) tquadrants->elem_count;
// #pragma omp parallel for
//             for (p4est_locidx_t e = nodes->global_offset; e < Q+nodes->global_offset; e++)
//             {
//                 p4est_quadrant_t * quad = p4est_quadrant_array_index(tquadrants, e);

//                 if(quad->x == 0)
//                 {
//                     interpolateData.direction=0;
//                     p4est_iterate(p4est, NULL, &interpolateData, get_interpolateNodesOnFacesWorker_data, NULL, NULL);

//                     S* o = out.getSampleDataRW(e, sentinel);
//                     for (index_t i=0; i < numComp; ++i) {
//                         o[e] = (fxx[INDEX3(0,i,e,numComp,numNodes)] + fxx[INDEX3(1,i,e,numComp,numNodes)])/static_cast<S>(2);
//                     }
//                 }

//                 if(quad->x == P4EST_ROOT_LEN)
//                 {
//                     interpolateData.direction=1;
//                     p4est_iterate(p4est, NULL, &interpolateData, get_interpolateNodesOnFacesWorker_data, NULL, NULL);

//                     S* o = out.getSampleDataRW(e, sentinel);
//                     for (index_t i=0; i < numComp; ++i) {
//                         o[e] = (fxx[INDEX3(2,i,e,numComp,numNodes)] + fxx[INDEX3(3,i,e,numComp,numNodes)])/static_cast<S>(2);
//                     }
//                 }

//                 if(quad->y == 0)
//                 {
//                     interpolateData.direction=2;
//                     p4est_iterate(p4est, NULL, &interpolateData, get_interpolateNodesOnFacesWorker_data, NULL, NULL);

//                     S* o = out.getSampleDataRW(e, sentinel);
//                     for (index_t i=0; i < numComp; ++i) {
//                         o[e] = (fxx[INDEX3(0,i,e,numComp,numNodes)] + fxx[INDEX3(2,i,e,numComp,numNodes)])/static_cast<S>(2);
//                     }
//                 }

//                 if(quad->y == P4EST_ROOT_LEN)
//                 {
//                     interpolateData.direction=3;
//                     p4est_iterate(p4est, NULL, &interpolateData, get_interpolateNodesOnFacesWorker_data, NULL, NULL);

//                     S* o = out.getSampleDataRW(e, sentinel);
//                     for (index_t i=0; i < numComp; ++i) {
//                         o[e] = (fxx[INDEX3(1,i,e,numComp,numNodes)] + fxx[INDEX3(3,i,e,numComp,numNodes)])/static_cast<S>(2);
//                     }
//                 }
//             }
//         }
//     } else {
//         out.requireWrite();
//         const S c0 = 0.21132486540518711775;
//         const S c1 = 0.78867513459481288225;

//         double * fxx = new double[4*numComp*numNodes];
//         // This structure is used to store info needed by p4est
//         interpolateNodesOnFacesWorker_Data<S> interpolateData;
//         interpolateData.fxx = fxx;
//         interpolateData.sentinel = sentinel;
//         interpolateData.offset = numComp*sizeof(S);
//         for(p4est_topidx_t treeid = p4est->first_local_tree; treeid <= p4est->last_local_tree; treeid++)
//         {
//             p4est_tree_t * currenttree = p4est_tree_array_index(p4est->trees, treeid);
//             sc_array_t * tquadrants = &currenttree->quadrants;
//             p4est_locidx_t Q = (p4est_locidx_t) tquadrants->elem_count;
// #pragma omp parallel for
//             for (p4est_locidx_t e = nodes->global_offset; e < Q+nodes->global_offset; e++)
//             {
//                 p4est_quadrant_t * quad = p4est_quadrant_array_index(tquadrants, e);

//                 if(quad->x == 0)
//                 {
//                     interpolateData.direction=0;
//                     p4est_iterate(p4est, NULL, &interpolateData, get_interpolateNodesOnFacesWorker_data, NULL, NULL);
//                     S* o = out.getSampleDataRW(e, sentinel);
//                     for (index_t i=0; i < numComp; ++i) {
//                         o[INDEX2(i,numComp,0)] = c0*fxx[INDEX3(1,i,e,numComp,numNodes)] + c1*fxx[INDEX3(0,i,e,numComp,numNodes)];
//                         o[INDEX2(i,numComp,1)] = c0*fxx[INDEX3(0,i,e,numComp,numNodes)] + c1*fxx[INDEX3(1,i,e,numComp,numNodes)];
//                     }
//                 }

//                 if(quad->x == P4EST_ROOT_LEN)
//                 {
//                     interpolateData.direction=1;
//                     p4est_iterate(p4est, NULL, &interpolateData, get_interpolateNodesOnFacesWorker_data, NULL, NULL);
//                     S* o = out.getSampleDataRW(e, sentinel);
//                     for (index_t i=0; i < numComp; ++i) {
//                         o[INDEX2(i,numComp,0)] = c1*fxx[INDEX3(2,i,e,numComp,numNodes)] + c0*fxx[INDEX3(3,i,e,numComp,numNodes)];
//                         o[INDEX2(i,numComp,1)] = c1*fxx[INDEX3(3,i,e,numComp,numNodes)] + c0*fxx[INDEX3(2,i,e,numComp,numNodes)];
//                     }
//                 }

//                 if(quad->y == 0)
//                 {
//                     interpolateData.direction=2;
//                     p4est_iterate(p4est, NULL, &interpolateData, get_interpolateNodesOnFacesWorker_data, NULL, NULL);
//                     S* o = out.getSampleDataRW(e, sentinel);
//                     for (index_t i=0; i < numComp; ++i) {
//                         o[INDEX2(i,numComp,0)] = c0*fxx[INDEX3(2,i,e,numComp,numNodes)] + c1*fxx[INDEX3(0,i,e,numComp,numNodes)];
//                         o[INDEX2(i,numComp,1)] = c0*fxx[INDEX3(0,i,e,numComp,numNodes)] + c1*fxx[INDEX3(2,i,e,numComp,numNodes)];
//                     }
//                 }

//                 if(quad->y == P4EST_ROOT_LEN)
//                 {
//                     interpolateData.direction=3;
//                     p4est_iterate(p4est, NULL, &interpolateData, get_interpolateNodesOnFacesWorker_data, NULL, NULL);
//                     S* o = out.getSampleDataRW(e, sentinel);
//                     for (index_t i=0; i < numComp; ++i) {
//                         o[INDEX2(i,numComp,0)] = c0*fxx[INDEX3(3,i,e,numComp,numNodes)] + c1*fxx[INDEX3(1,i,e,numComp,numNodes)];
//                         o[INDEX2(i,numComp,1)] = c0*fxx[INDEX3(1,i,e,numComp,numNodes)] + c1*fxx[INDEX3(3,i,e,numComp,numNodes)];
//                     }
//                 }
//             }
//         }
//     }
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
    // return NodeIDs.size() + hangingNodeIDs.size();
    return NodeIDs.size();
}

//protected
inline dim_t Rectangle::getNumElements() const
{
    return nodes->num_local_elements;
}

//protected
inline dim_t Rectangle::getNumFaceElements() const
{
    double xy[2];
    double x0 = forestData->m_origin[0];
    double y0 = forestData->m_origin[1];
    double x1 = forestData->m_length[0]+forestData->m_origin[0];
    double y1 = forestData->m_length[1]+forestData->m_origin[1];
    p4est_locidx_t numFaceElements = 0;
    for (p4est_topidx_t t = p4est->first_local_tree; t <= p4est->last_local_tree; t++) 
    {
        p4est_tree_t * currenttree = p4est_tree_array_index(p4est->trees, t);
        sc_array_t * tquadrants = &currenttree->quadrants;
        for (p4est_locidx_t e = 0; e < tquadrants->elem_count; e++)
        {
            p4est_quadrant_t * q = p4est_quadrant_array_index(tquadrants, e);
            p4est_qcoord_t length = P4EST_QUADRANT_LEN(q->level);
            p4est_qcoord_to_vertex(p4est->connectivity, t, q->x,q->y,xy);
            if(xy[0] == x0 || xy[1] == y0){
                numFaceElements++;
                break;
            }
            p4est_qcoord_to_vertex(p4est->connectivity, t, q->x+length,q->y+length,xy);
            if(xy[0] == x1 || xy[1] == y1){
                numFaceElements++;
                break;
            }
        }
    }
    return numFaceElements;
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
        for (int q = 0; q < Q; ++q) { // Loop over the elements attached to the tree
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
    data->phangingNodeIDs = &hangingNodeIDs;
    data->p4est = p4est;

    // This function loops over all interior faces
    // Note that it does not loop over the nodes on the boundaries
    // x = Lx and y = Ly
    p4est_iterate_ext(p4est, NULL, data, NULL, update_RC, NULL, true);

    // Find the indices of the nodes on the boundaries x = Lx and y = Ly
    for(p4est_topidx_t treeid = p4est->first_local_tree; treeid <= p4est->last_local_tree; ++treeid) {
        p4est_tree_t * tree = p4est_tree_array_index(p4est->trees, treeid);
        sc_array_t * tquadrants = &tree->quadrants;
        p4est_locidx_t Q = (p4est_locidx_t) tquadrants->elem_count;
#pragma omp parallel for
        for (int q = 0; q < Q; ++q) { // Loop over the elements attached to the tree
            p4est_quadrant_t * quad = p4est_quadrant_array_index(tquadrants, q);
            p4est_qcoord_t length = P4EST_QUADRANT_LEN(quad->level);
            // Loop over the four corners of the quadrant
            for(int n = 1; n < 3; ++n){
                int k = q - Q + nodeIncrements[treeid - p4est->first_local_tree];
                long lidx = nodes->element_nodes[4*k+n];
                if(lidx < nodes->owned_count) // if this process owns the node
                {
                    double lx = length * ((int) (n % 2) == 1);
                    double ly = length * ((int) (n / 2) == 1);
                    double xy[3];
                    p4est_qcoord_to_vertex(p4est->connectivity, treeid, quad->x+lx, quad->y+ly, xy);

                    // If the node is on the boundary x=Lx or y=Ly
                    if( (n == 1) && (xy[0] == forestData->m_lxy[0]) )
                      
                    {
                        // Get the node IDs
                        long lni0 = NodeIDs.find(std::make_pair(xy[0],xy[1]))->second;
                        p4est_qcoord_to_vertex(p4est->connectivity, treeid, quad->x+lx, quad->y+length, xy);
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
                    else if((n == 2) && (xy[1] == forestData->m_lxy[1]))
                    {
                        // Get the node IDs
                        long lni0 = NodeIDs.find(std::make_pair(xy[0],xy[1]))->second;
                        p4est_qcoord_to_vertex(p4est->connectivity, treeid, quad->x+length, quad->y+ly, xy);
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
        }
    }

    // Sorting
#pragma omp for
    for(int i = 0; i < getNumNodes(); i++){
        std::vector<long> * idx0 = &indices[0][i];
        std::sort(indices[0][i].begin()+1, indices[0][i].begin()+idx0[0][0]+1);
    }

#ifdef P4EST_ENABLE_DEBUG
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
#ifdef P4EST_ENABLE_DEBUG
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

    for (dim_t i=0; i<numNeighbours; i++) {
        const dim_t start = offsetInShared[i];
        const dim_t end = offsetInShared[i+1];
        for (dim_t j = start; j < end; j++) {
            if (j > start)
                doublyLink(colIndices, rowIndices, sendShared[j-1], j);
            doublyLink(colIndices, rowIndices, sendShared[j], j);
            if (j < end-1)
                doublyLink(colIndices, rowIndices, sendShared[j+1], j);
        }
    }
#pragma omp parallel for
    for (dim_t i = 0; i < numShared; i++) {
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
    throw OxleyException("populateSampleIds");

//     // degrees of freedom are numbered from left to right, bottom to top in
//     // each rank, continuing on the next rank (ranks also go left-right,
//     // bottom-top).
//     // This means rank 0 has id 0...n0-1, rank 1 has id n0...n1-1 etc. which
//     // helps when writing out data rank after rank.

//     // build node distribution vector first.
//     // rank i owns m_nodeDistribution[i+1]-nodeDistribution[i] nodes which is
//     // constant for all ranks in this implementation
//     IndexVector m_nodeDistribution = getNodeDistribution();

//     m_nodeDistribution.assign(m_mpiInfo->size+1, 0);
//     const dim_t numDOF=getNumDOF();
//     for (dim_t k=1; k<m_mpiInfo->size; k++) {
//         m_nodeDistribution[k]=k*numDOF;
//     }
//     m_nodeDistribution[m_mpiInfo->size]=getNumDataPointsGlobal();
//     try {
//         m_nodeId.resize(getNumNodes());
//         m_dofId.resize(numDOF);
//         m_elementId.resize(getNumElements());
//     } catch (const std::length_error& le) {
//         throw RipleyException("The system does not have sufficient memory for a domain of this size.");
//     }

//     // populate face element counts
//     //left
//     if (forestData->m_offset[0]==0)
//         forestData->m_faceCount[0]=forestData->m_NE[1];
//     else
//         forestData->m_faceCount[0]=0;
//     //right
//     if (m_mpiInfo->rank%forestData->m_NX[0]==forestData->m_NX[0]-1)
//         forestData->m_faceCount[1]=forestData->m_NE[1];
//     else
//         forestData->m_faceCount[1]=0;
//     //bottom
//     if (forestData->m_offset[1]==0)
//         forestData->m_faceCount[2]=forestData->m_NE[0];
//     else
//         forestData->m_faceCount[2]=0;
//     //top
//     if (m_mpiInfo->rank/forestData->m_NX[0]==forestData->m_NX[1]-1)
//         forestData->m_faceCount[3]=m_NE[0];
//     else
//         forestData->m_faceCount[3]=0;

//     const dim_t NFE = getNumFaceElements();
//     m_faceId.resize(NFE);

//     const index_t left = getFirstInDim(0);
//     const index_t bottom = getFirstInDim(1);
//     const dim_t nDOF0 = getNumDOFInAxis(0);
//     const dim_t nDOF1 = getNumDOFInAxis(1);
//     const dim_t NE0 = m_NE[0];
//     const dim_t NE1 = m_NE[1];

// #define globalNodeId(x,y) \
//     ((m_offset[0]+x)/nDOF0)*nDOF0*nDOF1+(m_offset[0]+x)%nDOF0 \
//     + ((m_offset[1]+y)/nDOF1)*nDOF0*nDOF1*m_NX[0]+((m_offset[1]+y)%nDOF1)*nDOF0

//     // set corner id's outside the parallel region
//     m_nodeId[0] = globalNodeId(0, 0);
//     m_nodeId[m_NN[0]-1] = globalNodeId(m_NN[0]-1, 0);
//     m_nodeId[m_NN[0]*(m_NN[1]-1)] = globalNodeId(0, m_NN[1]-1);
//     m_nodeId[m_NN[0]*m_NN[1]-1] = globalNodeId(m_NN[0]-1,m_NN[1]-1);
// #undef globalNodeId

// #pragma omp parallel
//     {
//         // populate degrees of freedom and own nodes (identical id)
// #pragma omp for nowait
//         for (dim_t i=0; i<nDOF1; i++) {
//             for (dim_t j=0; j<nDOF0; j++) {
//                 const index_t nodeIdx=j+left+(i+bottom)*m_NN[0];
//                 const index_t dofIdx=j+i*nDOF0;
//                 m_dofId[dofIdx] = m_nodeId[nodeIdx]
//                     = m_nodeDistribution[m_mpiInfo->rank]+dofIdx;
//             }
//         }

//         // populate the rest of the nodes (shared with other ranks)
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
//     } // end parallel section

//     m_nodeTags.assign(getNumNodes(), 0);
//     updateTagsInUse(Nodes);

//     m_elementTags.assign(getNumElements(), 0);
//     updateTagsInUse(Elements);

//     // generate face offset vector and set face tags
//     const index_t LEFT=1, RIGHT=2, BOTTOM=10, TOP=20;
//     const index_t faceTag[] = { LEFT, RIGHT, BOTTOM, TOP };
//     m_faceOffset.assign(4, -1);
//     m_faceTags.clear();
//     index_t offset=0;
//     for (size_t i=0; i<4; i++) {
//         if (m_faceCount[i]>0) {
//             m_faceOffset[i]=offset;
//             offset+=m_faceCount[i];
//             m_faceTags.insert(m_faceTags.end(), m_faceCount[i], faceTag[i]);
//         }
//     }

//     setTagMap("left", LEFT);
//     setTagMap("right", RIGHT);
//     setTagMap("bottom", BOTTOM);
//     setTagMap("top", TOP);
//     updateTagsInUse(FaceElements);

//     populateDofMap();
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
    m_nodeDistribution.assign(m_mpiInfo->size+1, 0);

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
                if(p4est->mpirank)
                {
                    double lx = length * ((int) (n % 2) == 1);
                    double ly = length * ((int) (n / 2) == 1);
                    double xy[3];
                    p4est_qcoord_to_vertex(p4est->connectivity, treeid, quad->x+lx, quad->y+ly, xy);
                    m_nodeDistribution[counter++]=NodeIDs.find(std::make_pair(xy[0],xy[1]))->second;
                }
        }
    }
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
    std::vector<IndexVector> indices(numNodes);

#pragma omp parallel for
    for(long r = 0; r <= numNodes; r++)
    {
        for(int c = myRows[r]; c < myRows[r+1]; c++)
            indices[r].push_back(myColumns[c]);
        // sort(indices[r].begin(), indices[r].end());
    }

    return indices;
}

bool Rectangle::operator==(const AbstractDomain& other) const
{
    const Rectangle* o=dynamic_cast<const Rectangle*>(&other);
    if (o) {
        return ((p4est_checksum(p4est) == p4est_checksum(o->p4est))
            && (forestData == o->forestData));
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
    for(p4est_topidx_t tree = p4est->first_local_tree; tree < p4est->last_local_tree; tree++) {
        p4est_tree_t * tree_t = p4est_tree_array_index(p4est->trees, tree);
        max_level = SC_MAX(max_level, tree_t->maxlevel);
    }
    
    double cx[3][P4EST_MAXLEVEL] = {{0}};
    double cy[3][P4EST_MAXLEVEL] = {{0}};
#pragma omp parallel for
    for(int i = 0; i < max_level; i++)
    {
        cx[0][i] = 0.21132486540518711775/forestData->m_dx[0][i];
        cx[1][i] = 0.78867513459481288225/forestData->m_dx[0][i];
        cx[2][i] = 1./forestData->m_dx[0][i];
        cy[0][i] = 0.21132486540518711775/forestData->m_dx[1][i];
        cy[1][i] = 0.78867513459481288225/forestData->m_dx[1][i];
        cy[2][i] = 1./forestData->m_dx[1][i];
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
            for (p4est_locidx_t e = nodes->global_offset; e < Q+nodes->global_offset; e++) // Loop over every quadrant within the tree
            {
                p4est_quadrant_t * quad = p4est_quadrant_array_index(tquadrants, e);
                quadrantData * quaddata = (quadrantData *) quad->p.user_data;

                memcpy(&f_00[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                memcpy(&f_01[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                memcpy(&f_10[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                memcpy(&f_11[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));

                Scalar* o = out.getSampleDataRW(e, zero);
                for (index_t i = 0; i < numComp; ++i) {
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
            for (p4est_locidx_t e = nodes->global_offset; e < Q+nodes->global_offset; e++) // Loop over every quadrant within the tree
            {
                p4est_quadrant_t * quad = p4est_quadrant_array_index(tquadrants, e);
                quadrantData * quaddata = (quadrantData *) quad->p.user_data;

                memcpy(&f_00[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                memcpy(&f_01[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                memcpy(&f_10[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                memcpy(&f_11[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));

                Scalar* o = out.getSampleDataRW(e, zero);

             for (index_t i = 0; i < numComp; ++i) {
                    o[INDEX3(i,0,0,numComp,2)] = (f_10[i] + f_11[i] - f_00[i] - f_01[i])*cx[2][i] * 0.5;
                    o[INDEX3(i,1,0,numComp,2)] = (f_01[i] + f_11[i] - f_00[i] - f_10[i])*cy[2][i] * 0.5;
                } 
            }
        }
    } else if (out.getFunctionSpace().getTypeCode() == FaceElements) {
        out.requireWrite();

        for (p4est_topidx_t t = p4est->first_local_tree; t <= p4est->last_local_tree; t++) 
        {
            p4est_tree_t * currenttree = p4est_tree_array_index(p4est->trees, t);
            sc_array_t * tquadrants = &currenttree->quadrants;
            p4est_locidx_t Q = (p4est_locidx_t) tquadrants->elem_count;
#pragma omp parallel for
            for (p4est_locidx_t e = nodes->global_offset; e < Q+nodes->global_offset; e++)
            {
                // Work out what level this element is on 
                p4est_quadrant_t * quad = p4est_quadrant_array_index(tquadrants, e);
                quadrantData * quaddata = (quadrantData *) quad->p.user_data;

                std::vector<Scalar> f_00(numComp, zero);
                std::vector<Scalar> f_01(numComp, zero);
                std::vector<Scalar> f_10(numComp, zero);
                std::vector<Scalar> f_11(numComp, zero);

                if (quaddata->m_faceOffset == 0) {
                    memcpy(&f_00[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_01[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_10[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_11[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    Scalar* o = out.getSampleDataRW(e, zero);
                    for (index_t i = 0; i < numComp; ++i) {
                        o[INDEX3(i,0,0,numComp,2)] = (f_10[i]-f_00[i])*cx[1][i] + (f_11[i]-f_01[i])*cx[0][i];
                        o[INDEX3(i,1,0,numComp,2)] = (f_01[i]-f_00[i])*cy[2][i];
                        o[INDEX3(i,0,1,numComp,2)] = (f_10[i]-f_00[i])*cx[0][i] + (f_11[i]-f_01[i])*cx[1][i];
                        o[INDEX3(i,1,1,numComp,2)] = (f_01[i]-f_00[i])*cy[2][i];
                    } // end of component loop i
                } // end of face 0
                if (quaddata->m_faceOffset == 1) {
                    memcpy(&f_00[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_01[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_10[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_11[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    Scalar* o = out.getSampleDataRW(e, zero);
                    for (index_t i = 0; i < numComp; ++i) {
                        o[INDEX3(i,0,0,numComp,2)] = (f_10[i]-f_00[i])*cx[1][i] + (f_11[i]-f_01[i])*cx[0][i];
                        o[INDEX3(i,1,0,numComp,2)] = (f_11[i]-f_10[i])*cy[2][i];
                        o[INDEX3(i,0,1,numComp,2)] = (f_10[i]-f_00[i])*cx[0][i] + (f_11[i]-f_01[i])*cx[1][i];
                        o[INDEX3(i,1,1,numComp,2)] = (f_11[i]-f_10[i])*cy[2][i];
                    } // end of component loop i
                } // end of face 1
                if (quaddata->m_faceOffset == 2) {
                    memcpy(&f_00[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_01[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_10[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_11[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    Scalar* o = out.getSampleDataRW(e, zero);
                    for (index_t i = 0; i < numComp; ++i) {
                        o[INDEX3(i,0,0,numComp,2)] = (f_10[i]-f_00[i])*cx[2][i];
                        o[INDEX3(i,1,0,numComp,2)] = (f_01[i]-f_00[i])*cy[1][i] + (f_11[i]-f_10[i])*cy[0][i];
                        o[INDEX3(i,0,1,numComp,2)] = (f_10[i]-f_00[i])*cx[2][i];
                        o[INDEX3(i,1,1,numComp,2)] = (f_01[i]-f_00[i])*cy[0][i] + (f_11[i]-f_10[i])*cy[1][i];
                    } // end of component loop i
                } // end of face 2
                if (quaddata->m_faceOffset == 3) {
                    memcpy(&f_00[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_01[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_10[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_11[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    Scalar* o = out.getSampleDataRW(e, zero);
                    for (index_t i = 0; i < numComp; ++i) {
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

        for (p4est_topidx_t t = p4est->first_local_tree; t <= p4est->last_local_tree; t++) 
        {
            p4est_tree_t * currenttree = p4est_tree_array_index(p4est->trees, t);
            sc_array_t * tquadrants = &currenttree->quadrants;
            p4est_locidx_t Q = (p4est_locidx_t) tquadrants->elem_count;
#pragma omp parallel for
            for (p4est_locidx_t e = nodes->global_offset; e < Q+nodes->global_offset; e++)
            {
                // Work out what level this element is on 
                p4est_quadrant_t * quad = p4est_quadrant_array_index(tquadrants, e);
                quadrantData * quaddata = (quadrantData *) quad->p.user_data;

                std::vector<Scalar> f_00(numComp, zero);
                std::vector<Scalar> f_01(numComp, zero);
                std::vector<Scalar> f_10(numComp, zero);
                std::vector<Scalar> f_11(numComp, zero);

                if (quaddata->m_faceOffset == 0) {
                    memcpy(&f_00[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_01[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_10[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_11[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    Scalar* o = out.getSampleDataRW(e, zero);
                    for (index_t i = 0; i < numComp; ++i) {
                        o[INDEX3(i,0,0,numComp,2)] = (f_10[i] + f_11[i] - f_00[i] - f_01[i])*cx[2][i] * 0.5;
                        o[INDEX3(i,1,0,numComp,2)] = (f_01[i]-f_00[i])*cy[2][i];
                    } // end of component loop i
                } // end of face 0
                if (quaddata->m_faceOffset == 1) {
                    memcpy(&f_00[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_01[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_10[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_11[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    Scalar* o = out.getSampleDataRW(e, zero);
                    for (index_t i = 0; i < numComp; ++i) {
                        o[INDEX3(i,0,0,numComp,2)] = (f_10[i] + f_11[i] - f_00[i] - f_01[i])*cx[2][i] * 0.5;
                        o[INDEX3(i,1,0,numComp,2)] = (f_11[i]-f_10[i])*cy[2][i];
                    } // end of component loop i
                } // end of face 1
                if (quaddata->m_faceOffset == 2) {
                    memcpy(&f_00[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_01[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_10[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_11[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    Scalar* o = out.getSampleDataRW(e, zero);
                    for (index_t i = 0; i < numComp; ++i) {
                        o[INDEX3(i,0,0,numComp,2)] = (f_10[i]-f_00[i])*cx[2][i];
                        o[INDEX3(i,1,0,numComp,2)] = (f_01[i] + f_11[i] - f_00[i] - f_10[i])*cy[2][i] * 0.5;
                    } // end of component loop i
                } // end of face 2
                if (quaddata->m_faceOffset == 3) {
                    memcpy(&f_00[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_01[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_10[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_11[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    Scalar* o = out.getSampleDataRW(e, zero);
                    for (index_t i = 0; i < numComp; ++i) {
                        o[INDEX3(i,0,0,numComp,2)] = (f_11[i]-f_01[i])*cx[2][i];
                        o[INDEX3(i,1,0,numComp,2)] = (f_01[i] + f_11[i] - f_00[i] - f_10[i])*cy[2][i] * 0.5;
                    } // end of component loop i
                }
            }
        }
    }
}

dim_t Rectangle::findNode(const double *coords) const
{
    throw OxleyException("findNode");
    return -1;

    // const dim_t NOT_MINE = -1;
    // //is the found element even owned by this rank
    // // (inside owned or shared elements but will map to an owned element)
    // for (int dim = 0; dim < m_numDim; dim++) {
    //     //allows for point outside mapping onto node
    //     double min = m_origin[dim] + m_offset[dim]* m_dx[dim]
    //             - m_dx[dim]/2. + escript::DataTypes::real_t_eps();
    //     double max = m_origin[dim] + (m_offset[dim] + m_NE[dim])*m_dx[dim]
    //             + m_dx[dim]/2. - escript::DataTypes::real_t_eps();
    //     if (min > coords[dim] || max < coords[dim]) {
    //         return NOT_MINE;
    //     }
    // }
    // // get distance from origin
    // double x = coords[0] - m_origin[0];
    // double y = coords[1] - m_origin[1];

    // //check if the point is even inside the domain
    // if (x < 0 || y < 0 || x > m_length[0] || y > m_length[1])
    //     return NOT_MINE;

    // // distance in elements
    // dim_t ex = (dim_t) floor((x + 0.01*m_dx[0]) / m_dx[0]);
    // dim_t ey = (dim_t) floor((y + 0.01*m_dx[1]) / m_dx[1]);
    // // set the min distance high enough to be outside the element plus a bit
    // dim_t closest = NOT_MINE;
    // double minDist = 1;
    // for (int dim = 0; dim < m_numDim; dim++) {
    //     minDist += m_dx[dim]*m_dx[dim];
    // }
    // //find the closest node
    // for (int dx = 0; dx < 1; dx++) {
    //     double xdist = (x - (ex + dx)*m_dx[0]);
    //     for (int dy = 0; dy < 1; dy++) {
    //         double ydist = (y - (ey + dy)*m_dx[1]);
    //         double total = xdist*xdist + ydist*ydist;
    //         if (total < minDist) {
    //             closest = INDEX2(ex+dx-m_offset[0], ey+dy-m_offset[1], m_NN[0]);
    //             minDist = total;
    //         }
    //     }
    // }
    // //if this happens, we've let a dirac point slip through, which is awful
    // if (closest == NOT_MINE) {
    //     throw RipleyException("Unable to map appropriate dirac point to a node,"
    //             " implementation problem in Rectangle::findNode()");
    // }
    // return closest;
}

//protected
void Rectangle::nodesToDOF(escript::Data& out, const escript::Data& in) const
{
    //AEAE todo
    throw OxleyException("nodesToDOF");

//     const dim_t numComp = in.getDataPointSize();
//     out.requireWrite();

//     const index_t left = getFirstInDim(0);
//     const index_t bottom = getFirstInDim(1);
//     const dim_t nDOF0 = getNumDOFInAxis(0);
//     const dim_t nDOF1 = getNumDOFInAxis(1);
// #pragma omp parallel for
//     for (index_t i=0; i<nDOF1; i++) {
//         for (index_t j=0; j<nDOF0; j++) {
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

static inline void
brick_linear_to_xyz (p4est_topidx_t ti, const int logx[P4EST_DIM],
                     const int rankx[P4EST_DIM], p4est_topidx_t tx[P4EST_DIM])
{
  int                 i, j, k;
  int                 lastlog = 0;

  for (i = 0; i < P4EST_DIM; i++) {
    tx[i] = 0;
  }

  for (i = 0; i < P4EST_DIM - 1; i++) {
    p4est_topidx_t      tempx[3] = { 0, 0, 0 };
    int                 logi = logx[rankx[i]] - lastlog;
    int                 idx[3] = { -1, -1, -1 };
    int                 c = 0;

    for (k = 0; k < P4EST_DIM - i; k++) {
      int                 d = rankx[i + k];

      idx[d] = 0;
    }
    for (k = 0; k < P4EST_DIM; k++) {
      if (idx[k] == 0) {
        idx[k] = c++;
      }
    }

    for (j = 0; j < logi; j++) {
      int                 base = (P4EST_DIM - i) * j;
      int                 shift = (P4EST_DIM - i - 1) * j;

      for (k = 0; k < P4EST_DIM; k++) {
        int                 id = idx[k];

        if (id >= 0) {
          tempx[k] |= (ti & (1 << (base + id))) >> (shift + id);
        }
      }
    }
    for (k = 0; k < P4EST_DIM; k++) {
      tx[k] += (tempx[k] << lastlog);
    }
    lastlog += logi;
    ti >>= (P4EST_DIM - i) * logi;
  }
  tx[rankx[P4EST_DIM - 1]] += (ti << lastlog);
}

static inline       p4est_topidx_t
brick_xyz_to_linear (const p4est_topidx_t tx[P4EST_DIM],
                     const int logx[P4EST_DIM], const int rankx[P4EST_DIM])
{
  int                 i, j, k;
  int                 lastlog = logx[rankx[P4EST_DIM - 2]];
  p4est_topidx_t      ti = tx[rankx[P4EST_DIM - 1]] >> lastlog;

  for (i = P4EST_DIM - 2; i >= 0; i--) {
    p4est_topidx_t      tempx[3] = { 0, 0, 0 };
    int                 logi =
      (i == 0) ? lastlog : lastlog - logx[rankx[i - 1]];
    int                 idx[3] = { -1, -1, -1 };
    int                 c = 0;

    for (k = 0; k < P4EST_DIM - i; k++) {
      int                 d = rankx[i + k];

      idx[d] = 0;
    }
    for (k = 0; k < P4EST_DIM; k++) {
      if (idx[k] == 0) {
        idx[k] = c++;
      }
    }

    ti <<= (P4EST_DIM - i) * logi;
    lastlog -= logi;
    for (k = 0; k < P4EST_DIM; k++) {
      tempx[k] = tx[k] >> lastlog;
    }
    for (j = 0; j < logi; j++) {
      int                 shift = (P4EST_DIM - i - 1) * j;

      for (k = 0; k < P4EST_DIM; k++) {
        int                 id = idx[k];

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
                            int mi, int ni, int periodic_a, int periodic_b, 
                            double x0, double x1, double y0, double y1)
{
    const p4est_topidx_t m = (p4est_topidx_t) mi;
    const p4est_topidx_t n = (p4est_topidx_t) ni;
    const p4est_topidx_t mc = periodic_a ? m : (m - 1);
    const p4est_topidx_t nc = periodic_b ? n : (n - 1);
    const p4est_topidx_t num_trees = m * n;
    const p4est_topidx_t num_corners = mc * nc;
    const p4est_topidx_t num_ctt = P4EST_CHILDREN * num_corners;
    const p4est_topidx_t num_vertices = (m + 1) * (n + 1);
    const int           periodic[P4EST_DIM] = { periodic_a, periodic_b };
    const p4est_topidx_t max[P4EST_DIM] = { m - 1, n - 1 };

    double             *vertices;
    p4est_topidx_t     *tree_to_vertex;
    p4est_topidx_t     *tree_to_tree;
    int8_t             *tree_to_face;
    p4est_topidx_t     *tree_to_corner;
    p4est_topidx_t     *ctt_offset;
    p4est_topidx_t     *corner_to_tree;
    int8_t             *corner_to_corner;
    p4est_topidx_t      n_iter;
    int                 logx[P4EST_DIM];
    int                 rankx[P4EST_DIM];
    int                 i, j, l;
    p4est_topidx_t      ti, tj, tk;
    p4est_topidx_t      tx, ty;
    p4est_topidx_t      tf[P4EST_FACES], tc[P4EST_CHILDREN];
    p4est_topidx_t      coord[P4EST_DIM], coord2[P4EST_DIM], ttemp;
    p4est_topidx_t     *linear_to_tree;
    p4est_topidx_t     *tree_to_corner2;
    p4est_topidx_t      vcount = 0, vicount = 0;
    int                 c[P4EST_DIM];
    p4est_connectivity_t *conn;

    P4EST_ASSERT (m > 0 && n > 0);

    if(num_trees > MAXTREES)
        throw OxleyException("n0*n1 exceeds the maximum number of allowable trees.");
    else
        conn = p4est_connectivity_new(num_vertices, num_trees, num_corners, num_ctt);

    double dx = (x1 - x0) / mi;
    double dy = (y1 - y0) / ni;

    vertices = conn->vertices;
    tree_to_vertex = conn->tree_to_vertex;
    tree_to_tree = conn->tree_to_tree;
    tree_to_face = conn->tree_to_face;
    tree_to_corner = conn->tree_to_corner;
    ctt_offset = conn->ctt_offset;
    corner_to_tree = conn->corner_to_tree;
    corner_to_corner = conn->corner_to_corner;

    for (ti = 0; ti < num_corners + 1; ti++) {
        ctt_offset[ti] = P4EST_CHILDREN * ti;
    }

    for (ti = 0; ti < P4EST_CHILDREN * num_trees; ti++) {
        tree_to_vertex[ti] = -1;
    }

    logx[0] = SC_LOG2_32 (m - 1) + 1;
    logx[1] = SC_LOG2_32 (n - 1) + 1;
    n_iter = (1 << logx[0]) * (1 << logx[1]);
    if (logx[0] <= logx[1]) {
        rankx[0] = 0;
        rankx[1] = 1;
    }
    else {
        rankx[0] = 1;
        rankx[1] = 0;
    }

    linear_to_tree = P4EST_ALLOC (p4est_topidx_t, n_iter);
    tree_to_corner2 = P4EST_ALLOC (p4est_topidx_t, num_trees);

    tj = 0;
    tk = 0;
    for (ti = 0; ti < n_iter; ti++) {
    brick_linear_to_xyz (ti, logx, rankx, coord);
    tx = coord[0];
    ty = coord[1];
    if (tx < m && ty < n &&
        1) {
      linear_to_tree[ti] = tj;
      if ((tx < m - 1 || periodic_a) && (ty < n - 1 || periodic_b) &&
          1) {
        tree_to_corner2[tj] = tk++;
      }
      else {
        tree_to_corner2[tj] = -1;
      }
      tj++;
    }
    else {
      linear_to_tree[ti] = -1;
    }
    }
    P4EST_ASSERT (tj == num_trees);
    P4EST_ASSERT (tk == num_corners);

    for (ti = 0; ti < n_iter; ti++) {
    brick_linear_to_xyz (ti, logx, rankx, coord);
    tx = coord[0];
    ty = coord[1];
    if (tx < m && ty < n &&
        1) {
      tj = linear_to_tree[ti];
      P4EST_ASSERT (tj >= 0);
      for (i = 0; i < P4EST_DIM; i++) {
        for (j = 0; j < 2; j++) {
          l = 2 * i + j;
          coord2[0] = ((tx + ((i == 0) ? (2 * j - 1) : 0)) + m) % m;
          coord2[1] = ((ty + ((i == 1) ? (2 * j - 1) : 0)) + n) % n;
          tf[l] = brick_xyz_to_linear (coord2, logx, rankx);
          P4EST_ASSERT (tf[l] < n_iter);
          tf[l] = linear_to_tree[tf[l]];
          P4EST_ASSERT (tf[l] >= 0);
        }
      }
      for (i = 0; i < P4EST_CHILDREN; i++) {
        coord2[0] = ((tx + (((i & 1) == 0) ? -1 : 1)) + m) % m;
        coord2[1] = ((ty + ((((i >> 1) & 1) == 0) ? -1 : 1)) + n) % n;
        tc[i] = brick_xyz_to_linear (coord2, logx, rankx);
        P4EST_ASSERT (tc[i] < n_iter);
        tc[i] = linear_to_tree[tc[i]];
        P4EST_ASSERT (tc[i] >= 0);
      }
      for (i = 0; i < P4EST_DIM; i++) {
        for (j = 0; j < 2; j++) {
          l = i * 2 + j;
          if (!periodic[i] &&
              ((coord[i] == 0 && j == 0) || (coord[i] == max[i] && j == 1))) {
            tree_to_tree[tj * P4EST_FACES + l] = tj;
            tree_to_face[tj * P4EST_FACES + l] = (int8_t) l;
          }
          else {
            tree_to_tree[tj * P4EST_FACES + l] = tf[l];
            tree_to_face[tj * P4EST_FACES + l] = (int8_t) (i * 2 + (j ^ 1));
          }
        }
      }
      for (i = 0; i < P4EST_CHILDREN; i++) {
        if (tree_to_corner != NULL) {
          c[0] = i & 1;
          c[1] = (i >> 1) & 1;
          if ((!periodic[0] &&
               ((coord[0] == 0 && c[0] == 0) ||
                (coord[0] == max[0] && c[0] == 1))) ||
              (!periodic[1] &&
               ((coord[1] == 0 && c[1] == 0) ||
                (coord[1] == max[1] && c[1] == 1))) ||
              0) {
            tree_to_corner[tj * P4EST_CHILDREN + i] = -1;
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
              SC_ABORT_NOT_REACHED ();
            }
            ttemp = tree_to_corner2[ttemp];
            P4EST_ASSERT (ttemp >= 0);
            tree_to_corner[tj * P4EST_CHILDREN + i] = ttemp;
            corner_to_tree[ttemp * P4EST_CHILDREN + (P4EST_CHILDREN - 1 - i)] = tj;
            corner_to_corner[ttemp * P4EST_CHILDREN + (P4EST_CHILDREN - 1 - i)] = (int8_t) i;
          }
        }
        if (ty > 0 && ((i >> 1) & 1) == 0) {
          tree_to_vertex[tj * P4EST_CHILDREN + i] =
            tree_to_vertex[tf[2] * P4EST_CHILDREN + i + 2];
        }
        else if (tx > 0 && (i & 1) == 0) {
          tree_to_vertex[tj * P4EST_CHILDREN + i] =
            tree_to_vertex[tf[0] * P4EST_CHILDREN + i + 1];
        }
        else {
          tree_to_vertex[tj * P4EST_CHILDREN + i] = vcount++;
          vertices[vicount++] = (double) dx * (tx + (i & 1));
          vertices[vicount++] = (double) dy * (ty + ((i >> 1) & 1));
          vertices[vicount++] = 0.;
        }
      }
    }
    }

    P4EST_ASSERT (vcount == num_vertices);
    P4EST_FREE (linear_to_tree);
    P4EST_FREE (tree_to_corner2);
    P4EST_ASSERT (p4est_connectivity_is_valid (conn));

    return conn;
}

// instantiate our two supported versions
template
void Rectangle::assembleGradientImpl<real_t>(escript::Data& out,
                                             const escript::Data& in) const;

template
void Rectangle::assembleGradientImpl<cplx_t>(escript::Data& out,
                                             const escript::Data& in) const;


} // end of namespace oxley
