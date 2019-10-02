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

#include <oxley/Brick.h>
#include <oxley/OxleyData.h>
#include <oxley/InitAlgorithms.h>
#include <oxley/RefinementAlgorithms.h>

#include <p8est_extended.h>
#include <p8est_vtk.h>
#include <p8est_iterate.h>

namespace bp = boost::python;

namespace oxley {

    /**
       \brief
       Constructor
    */
Brick::Brick(int order,
      dim_t n0, dim_t n1, dim_t n2,
      double x0, double y0, double z0,
      double x1, double y1, double z1,
      int d0, int d1, int d2,
      int periodic0, int periodic1, int periodic2): OxleyDomain(3, order){

    // Possible error: User passes invalid values for the dimensions
    if(n0 <= 0 || n1 <= 0 || n2 <= 0)
        throw OxleyException("Number of elements in each spatial dimension must be positive");

    // Ignore d0 and d1 if we are running in serial
    m_mpiInfo = escript::makeInfo(MPI_COMM_WORLD);
    if(m_mpiInfo->size == 1) {
        d0=1;
        d1=1;
        d2=1;
    }

    // ensure number of subdivisions is valid and nodes can be distributed
    // among number of ranks
    if(d0*d1*d2 != m_mpiInfo->size)
        throw OxleyException("Invalid number of spatial subdivisions");

    // These two statements configure the level of verbosity used by p4est
    sc_init(m_mpiInfo->comm, 1, LOG_BACKTRACE, NULL, LOG_LEVEL);
    p4est_init(NULL, LOG_LEVEL);

    //Create a connectivity
    connectivity = p8est_connectivity_new_brick((int) n0, (int) n1, (int) n2,
        periodic0, periodic1, periodic2);
    if(!p8est_connectivity_is_valid(connectivity))
        throw OxleyException("Could not create a valid connectivity.");

    // Allocate some memory
    forestData = (p8estData *) malloc(sizeof(p8estData));

    // p8est = p8est_new(m_mpiInfo->comm, connectivity, sizeof(forestData), &init_brick_data, NULL);
    p4est_locidx_t min_quadrants = n0*n1*n2 / m_mpiInfo->size;
    int min_level = 0;
    int fill_uniform = 1;
    p8est = p8est_new_ext(m_mpiInfo->comm, connectivity, min_quadrants,
            min_level, fill_uniform, sizeof(p8estData), &init_brick_data, (void *) &forestData);

    // Record the physical dimensions of the domain and the location of the origin
    forestData->m_origin[0] = x0;
    forestData->m_origin[1] = y0;
    forestData->m_origin[2] = z0;
    forestData->m_length[0] = x1-x0;
    forestData->m_length[1] = y1-y0;
    forestData->m_length[2] = z1-z0;

    // number of elements in each dimension
    forestData->m_gNE[0] = n0;
    forestData->m_gNE[1] = n1;
    forestData->m_gNE[2] = n2;

    // Whether or not we have periodic boundaries
    forestData->periodic[0] = periodic0;
    forestData->periodic[1] = periodic1;
    forestData->periodic[2] = periodic2;

    //number of elements for this rank in each dimension including shared
    forestData->m_NX[0] = d0;
    forestData->m_NX[1] = d1;
    forestData->m_NX[2] = d2;

    // local number of elements
    forestData->m_NE[0] = forestData->m_gNE[0] / d0;
    forestData->m_NE[1] = forestData->m_gNE[1] / d1;
    forestData->m_NE[2] = forestData->m_gNE[2] / d2;

    // max levels of refinement
    forestData->max_levels_refinement = MAXREFINEMENTLEVELS;

    // element order
    m_order = order;

    // initial tag
    tags[0] = 0;
    numberOfTags=1;

    // Number of dimensions
    m_numDim=3;

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
Brick::~Brick(){

    free(forestData);
    p8est_connectivity_destroy(connectivity);
    p8est_destroy(p8est);
    sc_finalize();

}

/**
   \brief
   returns a description for this domain
*/
std::string Brick::getDescription() const{

    return "oxley::Brick";

}


/**
   \brief
   dumps the mesh to a file with the given name
   \param filename The name of the output file
*/
void Brick::dump(const std::string& filename) const
{
    throw OxleyException("dump: not supported");
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

bool Brick::probeInterpolationAcross(int fsType_source,
        const escript::AbstractDomain& domain, int fsType_target) const
{
    throw OxleyException("currently: not supported"); //AE this is temporary
}

void Brick::interpolateAcross(escript::Data& target, const escript::Data& source) const
{
    throw OxleyException("currently: not supported"); //AE this is temporary
}

void Brick::setToNormal(escript::Data& out) const
{
    throw OxleyException("currently: not supported"); //AE this is temporary
}

void Brick::setToSize(escript::Data& out) const
{
    throw OxleyException("currently: not supported"); //AE this is temporary
}

bool Brick::ownSample(int fsType, index_t id) const
{
    throw OxleyException("currently: not supported"); //AE this is temporary
}


/* This is a wrapper for filtered (and non-filtered) randoms
 * For detailed doco see randomFillWorker
*/
escript::Data Brick::randomFill(const escript::DataTypes::ShapeType& shape,
                                    const escript::FunctionSpace& fs,
                                    long seed, const bp::tuple& filter) const
{
    throw OxleyException("currently: not supported"); //AE this is temporary
}


const dim_t* Brick::borrowSampleReferenceIDs(int fsType) const
{
    throw OxleyException("currently: not supported"); //AE this is temporary
}

void Brick::writeToVTK(std::string filename, bool writeToVTK) const
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
    if(!writeToVTK)
    {
        p8est_vtk_write_file(p8est, NULL, name);
    }
    else
    {
        // Get the number of quadrants in this thread
        p4est_locidx_t numquads = p8est->local_num_quadrants;

        // This vector will have one value for each quadrant with double precision
        sc_array_t * tagNumber = sc_array_new_size(sizeof(double), numquads);

        // Iterate over the p8est and fill the solution values into the sc_array
        p8est_iterate(p8est, NULL, (void *) tagNumber, getTagVector, NULL, NULL, NULL);

        // Create the context for the VTK file
        p8est_vtk_context_t * context = p8est_vtk_context_new(p8est, name);

        double scale = 1.00;
        p8est_vtk_context_set_scale(context, scale);

        // Start writing the file
        context = p8est_vtk_write_header(context);
        SC_CHECK_ABORT (context != NULL, P8EST_STRING "writeToVTK: Error writing vtk header.");
        context = p8est_vtk_write_cell_dataf(context, 1, 1, 0, 0, 1, 0, "tag", tagNumber, context);
        SC_CHECK_ABORT (context != NULL, P8EST_STRING "writeToVTK: Error writing cell data.");
        context = p8est_vtk_write_point_dataf(context, 0, 0, context);
        SC_CHECK_ABORT (context != NULL, P8EST_STRING "writeToVTK: Error writing point data.");
        int temp = p8est_vtk_write_footer(context);
        SC_CHECK_ABORT(!temp, P8EST_STRING "_vtk: Error writing footer");

        // Cleanup
        // p4est_vtk_context_destroy(context);
        sc_array_destroy(tagNumber);
    }

}

void Brick::refineMesh(int maxRecursion, std::string algorithmname)
{
    if(!algorithmname.compare("uniform")){
        p8est_refine_ext(p8est, true, maxRecursion, refine_uniform, NULL, refine_copy_parent_octant);
        p8est_balance_ext(p8est, P8EST_CONNECT_FULL, NULL, refine_copy_parent_octant);
        int partforcoarsen = 1;
        p8est_partition(p8est, partforcoarsen, NULL);
    } else {
        throw OxleyException("Unknown refinement algorithm name.");
    }
}

} // end of namespace oxley
