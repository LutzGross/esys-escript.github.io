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

#include <oxley/Rectangle.h>
#include <oxley/InitAlgorithms.h>
#include <oxley/OxleyData.h>
#include <oxley/RefinementAlgorithms.h>

#include <p4est_extended.h>
#include <p4est_vtk.h>
#include <p4est_iterate.h>

#include <p4est_algorithms.h>

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

    // ensure number of subdivisions is valid and nodes can be distributed
    // among number of ranks
    if(d0*d1 != m_mpiInfo->size)
        throw OxleyException("Invalid number of spatial subdivisions");

    // These two statements configure the level of verbosity used by p4est
    sc_set_log_defaults(NULL, NULL, LOG_LEVEL);
    p4est_init(NULL, LOG_LEVEL);

    // Create the connectivity
    const p4est_topidx_t num_vertices = 4;
    const p4est_topidx_t num_trees = 1;
    const p4est_topidx_t num_ctt = 0;
    const double vertices[4 * 3] = {    x0, y0, 0,
                                        x1, y0, 0,
                                        x0, y1, 0,
                                        x1, y1, 0,  };
    const p4est_topidx_t tree_to_vertex[4] = {0, 1, 2, 3,};
    const p4est_topidx_t tree_to_tree[4] = {0, 0, 0, 0,};
    const int8_t tree_to_face[4] = {1, 0, 3, 2,};
    const p4est_topidx_t tree_to_corner[4] = {0, 0, 0, 0,};
    const p4est_topidx_t ctt_offset[2] = {0, 4,};
    const p4est_topidx_t corner_to_tree[4] = {0, 0, 0, 0,};
    const int8_t corner_to_corner[4] = {0, 1, 2, 3,};

    connectivity = p4est_connectivity_new_copy (4, 1, 0,
                                        vertices, tree_to_vertex,
                                        tree_to_tree, tree_to_face,
                                        NULL, &num_ctt, NULL, NULL);

    // connectivity = p4est_connectivity_new_brick((int) n0, (int) n1, periodic0, periodic1);

    if(!p4est_connectivity_is_valid(connectivity))
        throw OxleyException("Could not create a valid connectivity.");

    // Allocate some memory
    forestData = new p4estData;
    forestData = (p4estData *) malloc(sizeof(p4estData));

    // p4est = p4est_new(m_mpiInfo->comm, connectivity, datasize, &init_rectangle_data, NULL);
    p4est_locidx_t min_quadrants = (n0*n1) / m_mpiInfo->size;
    int min_level = 0;
    int fill_uniform = 0;
    p4est = p4est_new_ext(m_mpiInfo->comm, connectivity, min_quadrants,
            min_level, fill_uniform, sizeof(p4estData), init_rectangle_data, (void *) &forestData);

    // Record the physical dimensions of the domain and the location of the origin
    forestData->m_origin[0] = x0;
    forestData->m_origin[1] = y0;
    forestData->m_length[0] = x1-x0;
    forestData->m_length[1] = y1-y0;

    // number of elements in each dimension
    forestData->m_gNE[0] = n0;
    forestData->m_gNE[1] = n1;

    // Whether or not we have periodic boundaries
    forestData->periodic[0] = periodic0;
    forestData->periodic[1] = periodic1;

    //number of elements for this rank in each dimension including shared
    forestData->m_NX[0] = d0;
    forestData->m_NX[1] = d1;

    // local number of elements
    forestData->m_NE[0] = forestData->m_gNE[0] / d0;
    forestData->m_NE[1] = forestData->m_gNE[1] / d1;

    // max levels of refinement
    forestData->max_levels_refinement = MAXREFINEMENTLEVELS;

    // element order
    m_order = order;

    // initial tag
    tags[0] = 0;
    numberOfTags=1;

    // Number of dimensions
    m_numDim=2;

    // To prevent segmentation faults from using numpy ndarray
#ifdef ESYS_HAVE_BOOST_NUMPY
    Py_Initialize();
    boost::python::numpy::initialize();
#endif
}

/**
   \brief
   Destructor.
*/
Rectangle::~Rectangle(){

    free(forestData);
    p4est_connectivity_destroy(connectivity);
    p4est_destroy(p4est);
    // sc_finalize();

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
    throw OxleyException("currently: not supported"); //AE this is temporary
}

void Rectangle::interpolateAcross(escript::Data& target, const escript::Data& source) const
{
    throw OxleyException("currently: not supported"); //AE this is temporary
}

void Rectangle::setToNormal(escript::Data& out) const
{
    throw OxleyException("currently: not supported"); //AE this is temporary
}

void Rectangle::setToSize(escript::Data& out) const
{
    throw OxleyException("currently: not supported"); //AE this is temporary
}

bool Rectangle::ownSample(int fsType, index_t id) const
{
    throw OxleyException("currently: not supported"); //AE this is temporary
}


/* This is a wrapper for filtered (and non-filtered) randoms
 * For detailed doco see randomFillWorker
*/
escript::Data Rectangle::randomFill(const escript::DataTypes::ShapeType& shape,
                                    const escript::FunctionSpace& fs,
                                    long seed, const bp::tuple& filter) const
{
    throw OxleyException("currently: not supported"); //AE this is temporary
}


const dim_t* Rectangle::borrowSampleReferenceIDs(int fsType) const
{
    throw OxleyException("currently: not supported"); //AE this is temporary
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

        // quad Tag
        sc_array_t * quadTag = sc_array_new_count(sizeof(double), numquads);
        p4est_iterate(p4est, NULL, (void *) quadTag, getQuadTagVector, NULL, NULL);
        sc_array_t * xcoord = sc_array_new_count(sizeof(double), numquads);
        p4est_iterate(p4est, NULL, (void *) xcoord, getXCoordVector, NULL, NULL);
        sc_array_t * ycoord = sc_array_new_count(sizeof(double), numquads);
        p4est_iterate(p4est, NULL, (void *) ycoord, getYCoordVector, NULL, NULL);

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
    if(!algorithmname.compare("uniform")){
        if(maxRecursion == 0){
            p4est_refine_ext(p4est, true, -1, refine_uniform, NULL, refine_copy_parent_quadrant);
        } else {
            p4est_refine_ext(p4est, true, maxRecursion+2, refine_uniform, NULL, refine_copy_parent_quadrant);
        }
        p4est_balance_ext(p4est, P4EST_CONNECT_FULL, NULL, refine_copy_parent_quadrant);
        int partforcoarsen = 1;
        p4est_partition(p4est, partforcoarsen, NULL);
    } else {
        throw OxleyException("Unknown refinement algorithm name.");
    }
}

escript::Data Rectangle::getX() const
{
    throw OxleyException("Currently not implemented"); //aeae this is temporary
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


} // end of namespace oxley
