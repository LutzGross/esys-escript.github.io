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
#include <random>
#include <vector>

#include <escript/Assert.h>
#include <escript/Data.h>
#include <escript/DataFactory.h>
#include <escript/FunctionSpaceFactory.h>

#include <oxley/AbstractAssembler.h>
#include <oxley/DefaultAssembler3D.h>
#include <oxley/InitAlgorithms.h>
#include <oxley/Oxley.h>
#include <oxley/OxleyData.h>
#include <oxley/Brick.h>
#include <oxley/RefinementAlgorithms.h>

#include <p8est.h>
#include <p8est_algorithms.h>
#include <p8est_bits.h>
#include <p8est_connectivity.h>
#include <p8est_extended.h>
#include <p8est_iterate.h>
#include <p8est_lnodes.h>
#include <p8est_vtk.h>

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
    const std::vector<double>& points, const std::vector<int>& tags,
    int periodic0, int periodic1, int periodic2): 
    OxleyDomain(3, order)
{

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

    // If the user did not set the number of divisions manually
    if(d0 == -1 && d1 == -1 && d2 == -1)
    {
        d0 = m_mpiInfo->size < 3 ? 1 : m_mpiInfo->size / 3;
        d1 = m_mpiInfo->size < 3 ? 1 : m_mpiInfo->size / 3;
        d2 = m_mpiInfo->size / (d0*d1);

        if(d0*d1*d2 != m_mpiInfo->size)
            throw OxleyException("Could not find values for d0, d1 and d2. Please set them manually.");
    }


    // else
    // {
    //     // ensure number of subdivisions chosen by the user is valid and nodes can be distributed
    //     // among number of ranks
    //     if(d0*d1*d2 != m_mpiInfo->size)
    //         throw OxleyException("Invalid number of spatial subdivisions");
    // }

    //Create a connectivity
    // const p4est_topidx_t num_vertices = 8;
    // const p4est_topidx_t num_trees = 1;
    // const p4est_topidx_t num_edges = 3;
    // const p4est_topidx_t num_corners = 1;
    // const double vertices[8 * 3] = {
    //                                 x0, y0, z0,
    //                                 x1, y0, z0,
    //                                 x0, y1, z0,
    //                                 x1, y1, z0,
    //                                 x0, y0, z1,
    //                                 x1, y0, z1,
    //                                 x0, y1, z1,
    //                                 x1, y1, z1,
    //                                 };
    // const p4est_topidx_t tree_to_vertex[8] = {0, 1, 2, 3, 4, 5, 6, 7,};
    // const p4est_topidx_t tree_to_tree[6] = {0, 0, 0, 0, 0, 0,};
    // // const int8_t tree_to_face[6] = {1, 0, 3, 2, 5, 4, };
    // const int8_t tree_to_face[6] = {0, 1, 2, 3, 4, 5 };
    // const p4est_topidx_t tree_to_edge[12] = {0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2,};
    // const p4est_topidx_t ett_offset[4] = {0, 4, 8, 12,};
    // const p4est_topidx_t edge_to_tree[12] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,};
    // const int8_t edge_to_edge[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11,};
    // const p4est_topidx_t tree_to_corner[8] = {0, 0, 0, 0, 0, 0, 0, 0,};
    // const p4est_topidx_t ctt_offset[2] = {0, 8,};
    // const p4est_topidx_t corner_to_tree[8] = {0, 0, 0, 0, 0, 0, 0, 0,};
    // const int8_t corner_to_corner[8] = {0, 1, 2, 3, 4, 5, 6, 7,};
    // connectivity = p8est_connectivity_new_copy(num_vertices, num_trees, num_edges,
    //                                   num_corners, vertices, tree_to_vertex,
    //                                   tree_to_tree, tree_to_face,
    //                                   tree_to_edge, ett_offset,
    //                                   edge_to_tree, edge_to_edge,
    //                                   tree_to_corner, ctt_offset,
    //                                   corner_to_tree, corner_to_corner);

    connectivity = new_brick_connectivity(n0, n1, n2, 0, 0, 0, x0, x1, y0, y1, z0, z1);

    // Create the p8est
    p4est_locidx_t min_quadrants = n0*n1*n2;
    int min_level = 0;
    int fill_uniform = 1;
    p8est = p8est_new_ext(m_mpiInfo->comm, connectivity, min_quadrants,
            min_level, fill_uniform, sizeof(octantData), &init_brick_data, (void *) &forestData);

    // Nodes numbering
    p8est_ghost_t * ghost = p8est_ghost_new(p8est, P8EST_CONNECT_FULL);
    nodes = p8est_lnodes_new(p8est, ghost, 1);
    p8est_ghost_destroy(ghost);

    // This information is needed by the assembler
    m_NE[0] = n0;
    m_NE[1] = n1;
    m_NE[2] = n2;
    m_NX[0] = (x1-x0)/n0;
    m_NX[1] = (y1-y0)/n1;
    m_NX[2] = (z1-z0)/n2;

    // Record the physical dimensions of the domain and the location of the origin
    forestData.m_origin[0] = x0;
    forestData.m_origin[1] = y0;
    forestData.m_origin[2] = z0;
    forestData.m_lxyz[0] = x1;
    forestData.m_lxyz[1] = y1;
    forestData.m_lxyz[2] = z1;
    forestData.m_length[0] = x1-x0;
    forestData.m_length[1] = y1-y0;
    forestData.m_length[2] = z1-z0;
    forestData.m_lxyz[0] = (x1-x0)/n0;
    forestData.m_lxyz[1] = (y1-y0)/n1;
    forestData.m_lxyz[2] = (z1-z0)/n2;

    // Whether or not we have periodic boundaries
    forestData.periodic[0] = periodic0;
    forestData.periodic[1] = periodic1;
    forestData.periodic[2] = periodic2;

    // Find the grid spacing for each level of refinement in the mesh
#pragma omp parallel for
    for(int i = 0; i <= P8EST_MAXLEVEL; i++){
        double numberOfSubDivisions = (p4est_qcoord_t) (1 << (P8EST_MAXLEVEL - i));
        forestData.m_dx[0][i] = forestData.m_length[0] / numberOfSubDivisions;
        forestData.m_dx[1][i] = forestData.m_length[1] / numberOfSubDivisions;
        forestData.m_dx[2][i] = forestData.m_length[2] / numberOfSubDivisions;
    }

    // max levels of refinement
    forestData.max_levels_refinement = MAXREFINEMENTLEVELS;

    // element order
    m_order = order;

    // initial tag
    // tags[0] = 0;
    // numberOfTags=1;

    // Number of dimensions
    m_numDim=3;

    //  // Distribute the p8est across the processors
    int allow_coarsening = 0;
    p8est_partition(p8est, allow_coarsening, NULL);

    // Number the nodes
    updateNodeIncrements();
    renumberNodes();
    updateRowsColumns();
    updateNodeDistribution();

    // Dirac points and tags
    addPoints(points, tags);

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
Brick::~Brick(){

    // free(forestData);
    p8est_connectivity_destroy(connectivity);
    p8est_destroy(p8est);
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

void Brick::interpolateAcross(escript::Data& target, const escript::Data& source) const
{
    throw OxleyException("currently: not supported"); //AE this is temporary
}

void Brick::setToNormal(escript::Data& out) const
{
    if (out.getFunctionSpace().getTypeCode() == FaceElements) 
    {
        out.requireWrite();
        for(p8est_topidx_t t = p8est->first_local_tree; t <= p8est->last_local_tree; t++) 
        {
            p8est_tree_t * currenttree = p8est_tree_array_index(p8est->trees, t);
            sc_array_t * tquadrants = &currenttree->quadrants;
            p4est_qcoord_t Q = (p4est_qcoord_t) tquadrants->elem_count;
#pragma omp parallel for
            for(p4est_locidx_t e = nodes->global_offset; e < Q+nodes->global_offset; e++)
            {
                // Work out what level this element is on 
                p8est_quadrant_t * quad = p8est_quadrant_array_index(tquadrants, e);
                octantData * quaddata = (octantData *) quad->p.user_data;               

                if(quaddata->m_faceOffset[0]) {
                    double* o = out.getSampleDataRW(e);
                    // set vector at four quadrature points
                    *o++ = -1.; *o++ = 0.; *o++ = 0.;
                    *o++ = -1.; *o++ = 0.; *o++ = 0.;
                    *o++ = -1.; *o++ = 0.; *o++ = 0.;
                    *o++ = -1.; *o++ = 0.; *o = 0.;
                }
            
                if(quaddata->m_faceOffset[1]) {
                    double* o = out.getSampleDataRW(e);
                    // set vector at four quadrature points
                    *o++ = 1.; *o++ = 0.; *o++ = 0.;
                    *o++ = 1.; *o++ = 0.; *o++ = 0.;
                    *o++ = 1.; *o++ = 0.; *o++ = 0.;
                    *o++ = 1.; *o++ = 0.; *o = 0.;
                }

                if(quaddata->m_faceOffset[2]) {
                    double* o = out.getSampleDataRW(e);                    
                    // set vector at four quadrature points
                    *o++ = 0.; *o++ = -1.; *o++ = 0.;
                    *o++ = 0.; *o++ = -1.; *o++ = 0.;
                    *o++ = 0.; *o++ = -1.; *o++ = 0.;
                    *o++ = 0.; *o++ = -1.; *o = 0.;
                }

                if(quaddata->m_faceOffset[3]) {
                    double* o = out.getSampleDataRW(e);
                    // set vector at four quadrature points
                    *o++ = 0.; *o++ = 1.; *o++ = 0.;
                    *o++ = 0.; *o++ = 1.; *o++ = 0.;
                    *o++ = 0.; *o++ = 1.; *o++ = 0.;
                    *o++ = 0.; *o++ = 1.; *o = 0.;
                }

                if(quaddata->m_faceOffset[4]) {
                    double* o = out.getSampleDataRW(e);
                    // set vector at four quadrature points
                    *o++ = 0.; *o++ = 0.; *o++ = -1.;
                    *o++ = 0.; *o++ = 0.; *o++ = -1.;
                    *o++ = 0.; *o++ = 0.; *o++ = -1.;
                    *o++ = 0.; *o++ = 0.; *o = -1.;
                }

                if(quaddata->m_faceOffset[5]) {
                    double* o = out.getSampleDataRW(e);
                    // set vector at four quadrature points
                    *o++ = 0.; *o++ = 0.; *o++ = 1.;
                    *o++ = 0.; *o++ = 0.; *o++ = 1.;
                    *o++ = 0.; *o++ = 0.; *o++ = 1.;
                    *o++ = 0.; *o++ = 0.; *o = 1.;
                }
            }
        }
    }
    else if (out.getFunctionSpace().getTypeCode() == ReducedFaceElements) 
    {
        out.requireWrite();
        for(p8est_topidx_t t = p8est->first_local_tree; t <= p8est->last_local_tree; t++) 
        {
            p8est_tree_t * currenttree = p8est_tree_array_index(p8est->trees, t);
            sc_array_t * tquadrants = &currenttree->quadrants;
            p4est_qcoord_t Q = (p4est_qcoord_t) tquadrants->elem_count;
#pragma omp parallel for
            for(p4est_locidx_t e = nodes->global_offset; e < Q+nodes->global_offset; e++)
            {
                // Work out what level this element is on 
                p8est_quadrant_t * quad = p8est_quadrant_array_index(tquadrants, e);
                octantData * quaddata = (octantData *) quad->p.user_data; 

                if(quaddata->m_faceOffset[0]) {
                    double* o = out.getSampleDataRW(e);
                    *o++ = -1.;
                    *o++ = 0.;
                    *o = 0.;
                }

                if(quaddata->m_faceOffset[1]) {
                    double* o = out.getSampleDataRW(e);
                    *o++ = 1.;
                    *o++ = 0.;
                    *o = 0.;
                }

                if(quaddata->m_faceOffset[2]) {
                    double* o = out.getSampleDataRW(e);
                    *o++ = 0.;
                    *o++ = -1.;
                    *o = 0.;
                }

                if(quaddata->m_faceOffset[3]) {
                    double* o = out.getSampleDataRW(e);
                    *o++ = 0.;
                    *o++ = 1.;
                    *o = 0.;
                }

                if(quaddata->m_faceOffset[4]) {
                    double* o = out.getSampleDataRW(e);
                    *o++ = 0.;
                    *o++ = 0.;
                    *o = -1.;
                }

                if(quaddata->m_faceOffset[5]) {
                    double* o = out.getSampleDataRW(e);
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
}

void Brick::setToSize(escript::Data& out) const
{
    if (out.getFunctionSpace().getTypeCode() == Elements
            || out.getFunctionSpace().getTypeCode() == ReducedElements) 
    {
        out.requireWrite();
        const dim_t numQuad=out.getNumDataPointsPerSample();

        for(p8est_topidx_t t = p8est->first_local_tree; t <= p8est->last_local_tree; t++) 
        {
            p4est_tree_t * currenttree = p4est_tree_array_index(p8est->trees, t);
            sc_array_t * tquadrants = &currenttree->quadrants;
            p4est_qcoord_t Q = (p4est_qcoord_t) tquadrants->elem_count;
        #pragma omp parallel for
            for(p4est_qcoord_t e = nodes->global_offset; e < Q+nodes->global_offset; e++)
            {
                // Work out what level this element is on 
                p8est_quadrant_t * quad = p8est_quadrant_array_index(tquadrants, e);
                // octantData * quaddata = (octantData *) quad->p.user_data; 
                const double size = sqrt(forestData.m_dx[0][quad->level]*forestData.m_dx[0][quad->level]
                                        +forestData.m_dx[1][quad->level]*forestData.m_dx[1][quad->level]
                                        +forestData.m_dx[2][quad->level]*forestData.m_dx[2][quad->level]);
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

        for(p8est_topidx_t t = p8est->first_local_tree; t <= p8est->last_local_tree; t++) 
        {
            p4est_tree_t * currenttree = p4est_tree_array_index(p8est->trees, t);
            sc_array_t * tquadrants = &currenttree->quadrants;
            p4est_qcoord_t Q = (p4est_qcoord_t) tquadrants->elem_count;
        #pragma omp parallel for
            for(p4est_qcoord_t e = nodes->global_offset; e < Q+nodes->global_offset; e++)
            {
                // Work out what level this element is on 
                p8est_quadrant_t * quad = p8est_quadrant_array_index(tquadrants, e);
                octantData * quaddata = (octantData *) quad->p.user_data;

                if (quaddata->m_faceOffset[0]) {
                    double* o = out.getSampleDataRW(e);
                    std::fill(o, o+numQuad, forestData.m_dx[1][quad->level]);
                }

                if (quaddata->m_faceOffset[1]) {
                    double* o = out.getSampleDataRW(e);
                    std::fill(o, o+numQuad, forestData.m_dx[1][quad->level]);
                }

                if (quaddata->m_faceOffset[2]) {
                    double* o = out.getSampleDataRW(e);
                    std::fill(o, o+numQuad, forestData.m_dx[0][quad->level]);
                }

                if (quaddata->m_faceOffset[3]) {
                    double* o = out.getSampleDataRW(e);
                    std::fill(o, o+numQuad, forestData.m_dx[0][quad->level]);
                }

                if (quaddata->m_faceOffset[4]) {
                    double* o = out.getSampleDataRW(e);
                    std::fill(o, o+numQuad, forestData.m_dx[2][quad->level]);
                }

                if (quaddata->m_faceOffset[5]) {
                    double* o = out.getSampleDataRW(e);
                    std::fill(o, o+numQuad, forestData.m_dx[2][quad->level]);
                }
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
            {
            // check ownership of element's bottom left node
            // return (m_dofMap[id%m_NE[0]+m_NN[0]*(id/m_NE[0])] < getNumDOF());
            throw OxleyException("Brick::ownSample Currently not implemented.");
            return false;
            }
        case FaceElements:
        case ReducedFaceElements:
            {
                // // check ownership of face element's last node
                // dim_t n=0;
                // for (size_t i=0; i<6; i++) {
                //     n+=m_faceCount[i];
                //     if (id<n) {
                //         const index_t j=id-n+m_faceCount[i];
                //         if (i>=4) { // front or back
                //             const index_t first=(i==4 ? 0 : m_NN[0]*m_NN[1]*(m_NN[2]-1));
                //             return (m_dofMap[first+j%m_NE[0]+1+(j/m_NE[0]+1)*m_NN[0]] < getNumDOF());
                //         } else if (i>=2) { // bottom or top
                //             const index_t first=(i==2 ? 0 : m_NN[0]*(m_NN[1]-1));
                //             return (m_dofMap[first+j%m_NE[0]+1+(j/m_NE[0]+1)*m_NN[0]*m_NN[1]] < getNumDOF());
                //         } else { // left or right
                //             const index_t first=(i==0 ? 0 : m_NN[0]-1);
                //             return (m_dofMap[first+(j%m_NE[1]+1)*m_NN[0]+(j/m_NE[1]+1)*m_NN[0]*m_NN[1]] < getNumDOF());
                //         }
                //     }
                // }
                throw OxleyException("Brick::ownSample Currently not implemented.");
                return false;
            }
        default:
            break;
    }

    std::stringstream msg;
    msg << "ownSample: invalid function space type " << fsType;
    throw ValueError(msg.str());
}

dim_t Brick::getNumDataPointsGlobal() const
{
    return getNumNodes();
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


const dim_t* Brick::borrowSampleReferenceIDs(int fsType) const
{
    switch (fsType) {
        case Nodes:
        case ReducedNodes: //FIXME: reduced
            return &myColumns[0];
        case DegreesOfFreedom:
        case ReducedDegreesOfFreedom: //FIXME: reduced
            throw OxleyException("Unknown Error.");
        case Elements:
        case ReducedElements:
            throw OxleyException("borrowSampleReferenceIDs: Elements");
        case FaceElements:
        case ReducedFaceElements:
            throw OxleyException("borrowSampleReferenceIDs: FaceIDS");
        case Points:
            throw OxleyException("borrowSampleReferenceIDs: Points");
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
        p4est_locidx_t numquads = p8est->local_num_quadrants;
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

void Brick::refineMesh(std::string algorithmname)
{
    double max_levels_refinement = 8; //TODO

    if(!algorithmname.compare("uniform")){
        p8est_refine_ext(p8est, true, max_levels_refinement, refine_uniform, NULL, refine_copy_parent_octant);
        p8est_balance_ext(p8est, P8EST_CONNECT_FULL, NULL, refine_copy_parent_octant);
        int partforcoarsen = 1;
        p8est_partition(p8est, partforcoarsen, NULL);
    } else {
        throw OxleyException("Unknown refinement algorithm name.");
    }

    // Make sure that nothing went wrong
#ifdef OXLEY_ENABLE_DEBUG
    if(!p8est_is_valid(p8est))
        throw OxleyException("p4est broke during refinement");
    if(!p8est_connectivity_is_valid(connectivity))
        throw OxleyException("connectivity broke during refinement");
#endif

    bool partition_for_coarsening = true;
    p8est_partition_ext(p8est, partition_for_coarsening, NULL);

    // Update the nodes
    p8est_lnodes_destroy(nodes);
    p8est_ghost_t * ghost = p8est_ghost_new(p8est, P8EST_CONNECT_FULL);
    nodes = p8est_lnodes_new(p8est, ghost, 1);
    p8est_ghost_destroy(ghost);

    // Update
    updateNodeIncrements();
    renumberNodes();
    updateRowsColumns();
}

void Brick::refineBoundary(std::string boundaryname, double dx)
{
    forestData.refinement_depth = dx;

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
    if(!p8est_connectivity_is_valid(connectivity))
        throw OxleyException("connectivity broke during refinement");
#endif

    bool partition_for_coarsening = true;
    p8est_partition_ext(p8est, partition_for_coarsening, NULL);

    // Update the nodes
    p8est_lnodes_destroy(nodes);
    p8est_ghost_t * ghost = p8est_ghost_new(p8est, P8EST_CONNECT_FULL);
    nodes = p8est_lnodes_new(p8est, ghost, 1);
    p8est_ghost_destroy(ghost);

    // Update
    updateNodeIncrements();
    renumberNodes();
    updateRowsColumns();
}

void Brick::refineRegion(double x0, double x1, double y0, double y1)
{
    
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
    if(!p8est_connectivity_is_valid(connectivity))
        std::cout << "WARNING: connectivity is invalid" << std::endl;
    std::cout << "temp_data = " << &temp_data << std::endl;
}

Assembler_ptr Brick::createAssembler(std::string type, const DataMap& constants) const
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
            return Assembler_ptr(new DefaultAssembler3D<cplx_t>(shared_from_this()));
        } else {
            return Assembler_ptr(new DefaultAssembler3D<real_t>(shared_from_this()));
        }
    } 
    throw escript::NotImplementedError("oxley::rectangle does not support the requested assembler");
}

// return True for a boundary node and False for an internal node
bool Brick::isBoundaryNode(p8est_quadrant_t * quad, int n, p8est_topidx_t treeid, p4est_qcoord_t length) const
{
    // AEAE TODO
    return false;
    // double lx = length * ((int) (n % 2) == 1);
    // double ly = length * ((int) (n / 2) == 1);
    // double xy[3];
    // p4est_qcoord_to_vertex(p4est->connectivity, treeid, quad->x+lx, quad->y+ly, xy);
    // return (xy[0] == 0) || (xy[0] == forestData.m_lxy[0]) || (xy[1] == 0) || (xy[1] == forestData.m_lxy[1]);
}

// returns True for a boundary node on the north or east of the domain
bool Brick::isUpperBoundaryNode(p8est_quadrant_t * quad, int n, p8est_topidx_t treeid, p4est_qcoord_t length) const
{
    // AEAE TODO
    return false;
//     double lx = length * ((int) (n % 2) == 1);
//     double ly = length * ((int) (n / 2) == 1);
//     double xy[3];
//     p4est_qcoord_to_vertex(p4est->connectivity, treeid, quad->x+lx, quad->y+ly, xy);
//     return (xy[0] == forestData.m_lxy[0]) || (xy[1] == forestData.m_lxy[1]);
}

// return True for a hanging node and False for an non-hanging node
bool Brick::isHangingNode(p8est_lnodes_code_t face_code, int n) const
{
    // AEAE TODO
    return false;
    // if(face_code == 0)
    // {
    //     return false;
    // }
    // else
    // {
    //     int8_t c = face_code & 0x03;
    //     int8_t ishanging0 = (face_code >> 2) & 0x01;
    //     int8_t ishanging1 = (face_code >> 3) & 0x01;

    //     int8_t f0 = p4est_corner_faces[c][0];
    //     int8_t f1 = p4est_corner_faces[c][1];

    //     // int d0 = c / 2 == 0;
    //     // int d1 = c % 2 == 0;
        
    //     return ((ishanging0 == 1) && (p4est_corner_face_corners[c][f0] == 1)) 
    //         || ((ishanging1 == 1) && (p4est_corner_face_corners[c][f1] == 1));
    // }
}

//protected
void Brick::updateNodeIncrements()
{
    nodeIncrements[0] = 1;
    for(p8est_topidx_t treeid = p8est->first_local_tree+1, k=1; treeid <= p8est->last_local_tree; ++treeid, ++k) 
    {
        p8est_tree_t * tree = p8est_tree_array_index(p8est->trees, treeid);
        sc_array_t * tquadrants = &tree->quadrants;
        p4est_qcoord_t Q = (p4est_qcoord_t) tquadrants->elem_count;
        nodeIncrements[k] = nodeIncrements[k-1] + Q;
    }
}

void Brick::renumberNodes()
{
    // AEAE TODO

// #ifdef OXLEY_ENABLE_DEBUG_NODES
//     std::cout << "Renumbering nodes... " << std::endl;
// #endif

//     NodeIDs.clear();
// #pragma omp for
//     for(p8est_topidx_t treeid = p8est->first_local_tree; treeid <= p8est->last_local_tree; ++treeid) {
//         p4est_tree_t * tree = p8est_tree_array_index(p8est->trees, treeid);
//         sc_array_t * tquadrants = &tree->quadrants;
//         p4est_qcoord_t Q = (p4est_qcoord_t) tquadrants->elem_count;
//         for(int q = 0; q < Q; ++q) { 
//             p8est_quadrant_t * quad = p8est_quadrant_array_index(tquadrants, q);
//             p4est_qcoord_t length = P8EST_QUADRANT_LEN(quad->level);
//             int k = q - Q + nodeIncrements[treeid - p4est->first_local_tree];
//             for(int n = 0; n < 4; n++)
//             {
//                 // if( (n == 0) 
//                 //     || isHangingNode(nodes->face_code[k], n)
//                 //     || isUpperBoundaryNode(quad, n, treeid, length) 
//                   // )
//                     double lx = length * ((int) (n % 2) == 1);
//                     double ly = length * ((int) (n / 2) == 1);
//                     double xy[3];
//                     p4est_qcoord_to_vertex(p4est->connectivity, treeid, quad->x+lx, quad->y+ly, xy);
//                     if(NodeIDs.count(std::make_pair(xy[0],xy[1]))==0)
//                     {
// #ifdef OXLEY_ENABLE_DEBUG_NODES
//                         std::cout << NodeIDs.size() << ": " << xy[0] << ", " << xy[1] << std::endl;
// #endif
//                         NodeIDs[std::make_pair(xy[0],xy[1])]=NodeIDs.size();
//                     }
//                 // }
//             }
//         }
//     }
}

// // protected
// inline index_t Brick::getFirstInDim(unsigned axis) const
// {
//     // return m_offset[axis] == 0 ? 0 : 1;
//     return nodes->global_offset == 0 ? 0 : 1;
// }

//protected
inline dim_t Brick::getNumNodes() const
{
    return NodeIDs.size();
}

//protected
inline dim_t Brick::getNumElements() const
{
    return nodes->num_local_elements;
}

//protected
inline dim_t Brick::getNumFaceElements() const
{
    double xyz[3];
    double x0 = forestData.m_origin[0];
    double y0 = forestData.m_origin[1];
    double z0 = forestData.m_origin[2];
    double x1 = forestData.m_length[0]+forestData.m_origin[0];
    double y1 = forestData.m_length[1]+forestData.m_origin[1];
    double z1 = forestData.m_length[2]+forestData.m_origin[2];
    p4est_qcoord_t numFaceElements = 0;
    for(p8est_topidx_t t = p8est->first_local_tree; t <= p8est->last_local_tree; t++) 
    {
        p8est_tree_t * currenttree = p8est_tree_array_index(p8est->trees, t);
        sc_array_t * tquadrants = &currenttree->quadrants;
        for(p4est_qcoord_t e = 0; e < tquadrants->elem_count; e++)
        {
            p8est_quadrant_t * q  = p8est_quadrant_array_index(tquadrants, e);
            p4est_qcoord_t length = P8EST_QUADRANT_LEN(q->level);
            p8est_qcoord_to_vertex(p8est->connectivity, t, q->x,q->y,q->z,xyz);
            if(xyz[0] == x0 || xyz[1] == y0 || xyz[2] == z0){
                numFaceElements++;
                break;
            }
            p8est_qcoord_to_vertex(p8est->connectivity, t, q->x+length,q->y+length,q->z+length,xyz);
            if(xyz[0] == x1 || xyz[1] == y1 || xyz[2] == z1){
                numFaceElements++;
                break;
            }
        }
    }
    return numFaceElements;
}

//protected
inline dim_t Brick::getNumDOF() const
{
    return getNumNodes();
}

void Brick::updateTreeIDs()
{
    treeIDs.clear();
    for(p8est_topidx_t treeid = p8est->first_local_tree; treeid <= p8est->last_local_tree; ++treeid) {
        p8est_tree_t * tree = p8est_tree_array_index(p8est->trees, treeid);
        sc_array_t * tquadrants = &tree->quadrants;
        p4est_qcoord_t Q = (p4est_qcoord_t) tquadrants->elem_count;
#pragma omp parallel for
        for(int q = 0; q < Q; ++q) { // Loop over the elements attached to the tree
            p8est_quadrant_t * quad = p8est_quadrant_array_index(tquadrants, q);
            treeIDs[std::make_tuple(quad->x,quad->y,quad->z)]=treeid;
        }
    }
}

void Brick::updateRowsColumns()
{

    // AEAE TODO

//     std::vector<std::vector<long>> * indices;
//     indices = new std::vector<std::vector<long>>;
//     int initial[] = {0, -1, -1, -1, -1, -1, -1};
//     indices->resize(getNumNodes(), std::vector<long>(initial, initial+7));

//     update_RC_data * data;
//     data = new update_RC_data;
//     data->indices = indices;
//     data->pNodeIDs = &NodeIDs;
//     data->p8est = p8est;

//     // This function loops over all interior faces
//     // Note that it does not loop over the nodes on the boundaries
//     // x = Lx and y = Ly
//     p8est_iterate_ext(p4est, NULL, data, NULL, update_RC, NULL, true);

//     // Find the indices of the nodes on the boundaries x = Lx and y = Ly
//     for(p8est_topidx_t treeid = p8est->first_local_tree; treeid <= p8est->last_local_tree; ++treeid) {
//         p4est_tree_t * tree = p8est_tree_array_index(p8est->trees, treeid);
//         sc_array_t * tquadrants = &tree->quadrants;
//         p4est_qcoord_t Q = (p4est_qcoord_t) tquadrants->elem_count;
// // #pragma omp parallel for
//         for(int q = 0; q < Q; ++q) { // Loop over the elements attached to the tree
//             p8est_quadrant_t * quad = p8est_quadrant_array_index(tquadrants, q);
//             p4est_qcoord_t length = P8EST_QUADRANT_LEN(quad->level);
//             // int k = q - Q + nodeIncrements[treeid - p4est->first_local_tree];
//             int n = 1;
//             // long lidx = nodes->element_nodes[4*k+n];
//             double lx = length * ((int) (n % 2) == 1);
//             double ly = length * ((int) (n / 2) == 1);
//             double xy[3];
//             p4est_qcoord_to_vertex(p8est->connectivity, treeid, quad->x+lx, quad->y+ly, xy);

//             // If the node is on the boundary x=Lx or y=Ly
//             if( (xy[0] == forestData.m_lxy[0]) )
//             {
//                 // Get the node IDs
//                 long lni0 = NodeIDs.find(std::make_pair(xy[0],xy[1]))->second;
//                 p4est_qcoord_to_vertex(p4est->connectivity, treeid, quad->x+lx, quad->y+length, xy);
//                 long lni1 = NodeIDs.find(std::make_pair(xy[0],xy[1]))->second;

//                 std::vector<long> * idx0 = &indices[0][lni0];
//                 std::vector<long> * idx1 = &indices[0][lni1];

//                 // Check for duplicates
//                 bool dup = false;
//                 for(int i = 1; i < idx0[0][0] + 1; i++)
//                     if(idx0[0][i] == lni1)
//                         dup = true;

//                 // Add the new indices
//                 if(dup == false)
//                 {
//                     idx0[0][0]++;
//                     idx1[0][0]++;
//                     idx0[0][idx0[0][0]]=lni1;
//                     idx1[0][idx1[0][0]]=lni0;
//                 }
//             }

//             n = 2;
//             // long lidx = nodes->element_nodes[4*k+n];
//             lx = length * ((int) (n % 2) == 1);
//             ly = length * ((int) (n / 2) == 1);
//             p4est_qcoord_to_vertex(p4est->connectivity, treeid, quad->x+lx, quad->y+ly, xy);
//             if((xy[1] == forestData.m_lxy[1]))
//             {
//                 // Get the node IDs
//                 long lni0 = NodeIDs.find(std::make_pair(xy[0],xy[1]))->second;
//                 p4est_qcoord_to_vertex(p4est->connectivity, treeid, quad->x+length, quad->y+ly, xy);
//                 long lni1 = NodeIDs.find(std::make_pair(xy[0],xy[1]))->second;

//                 std::vector<long> * idx0 = &indices[0][lni0];
//                 std::vector<long> * idx1 = &indices[0][lni1];

//                 // Check for duplicates
//                 bool dup = false;
//                 for(int i = 1; i < idx0[0][0] + 1; i++)
//                     if(idx0[0][i] == lni1)
//                         dup = true;

//                 // Add the new indices
//                 if(dup == false)
//                 {
//                     idx0[0][0]++;
//                     idx1[0][0]++;
//                     idx0[0][idx0[0][0]]=lni1;
//                     idx1[0][idx1[0][0]]=lni0;
//                 }
//             }
//         }
//     }

//     // Sorting
// #pragma omp for
//     for(int i = 0; i < getNumNodes(); i++){
//         std::vector<long> * idx0 = &indices[0][i];
//         std::sort(indices[0][i].begin()+1, indices[0][i].begin()+idx0[0][0]+1);
//     }

// #ifdef OXLEY_ENABLE_DEBUG_NODES
//     std::cout << "Node connections: " << std::endl;
//     // Output for debugging
//     for(int i = 0; i < getNumNodes(); i++){
//         std::vector<long> * idx0 = &indices[0][i];
//         std::cout << i << ": ";
//         for(int j = 1; j < idx0[0][0]+1; j++)
//             std::cout << idx0[0][j] << ", ";
//         std::cout << std::endl;
//     }
// #endif

//     // Convert to CRS format
//     myRows.clear();
//     myRows.push_back(0);
//     myColumns.clear();
//     m_dofMap.assign(getNumNodes(), 0);
//     long counter = 0;
//     for(int i = 0; i < getNumNodes(); i++)
//     {
//         std::vector<long> * idx0 = &indices[0][i];
//         std::vector<long> temp; 
//         for(int j = 1; j < idx0[0][0]+1; j++)
//         {
//             temp.push_back(idx0[0][j]);
//             counter++;
//         }
//         std::sort(temp.begin(),temp.end());
//         for(int i = 0; i < temp.size(); i++)
//             myColumns.push_back(temp[i]);
//         m_dofMap[i] = counter-myRows[i];
//         if(i < getNumNodes()-1)
//             myRows.push_back(counter);
//     }
//     myRows.push_back(myColumns.size());
// #ifdef OXLEY_ENABLE_DEBUG_NODES
//     std::cout << "Converted to Yale format... "<< std::endl;
//     std::cout << "COL_INDEX [";
//     for(auto i = myColumns.begin(); i < myColumns.end(); i++)
//         std::cout << *i << " ";
//     std::cout << "]" << std::endl;
//     std::cout << "ROW_INDEX [";
//     for(auto i = myRows.begin(); i < myRows.end(); i++)
//         std::cout << *i << " ";
//     std::cout << "]" << std::endl;
//     std::cout << "m_dofMap [";
//     for(auto i = m_dofMap.begin(); i < m_dofMap.end(); i++)
//         std::cout << *i << " ";
//     std::cout << "]" << std::endl;
// #endif

//     delete data;
//     delete indices;
}

//protected
void Brick::assembleCoordinates(escript::Data& arg) const
{

    // AEAE TODO

    // int numDim = m_numDim;
    // if (!arg.isDataPointShapeEqual(1, &numDim))
    //     throw ValueError("setToX: Invalid Data object shape");
    // if (!arg.numSamplesEqual(1, getNumNodes()))
    //     throw ValueError("setToX: Illegal number of samples in Data object");
    // arg.requireWrite();
    
//     for(p4est_topidx_t treeid = p4est->first_local_tree; treeid <= p4est->last_local_tree; ++treeid) {
//         p4est_tree_t * tree = p4est_tree_array_index(p4est->trees, treeid);
//         sc_array_t * tquadrants = &tree->quadrants;
//         p4est_locidx_t Q = (p4est_locidx_t) tquadrants->elem_count;

// #pragma omp parallel for
//         for(int q = 0; q < Q; ++q) { // Loop over the elements attached to the tree
//             p4est_quadrant_t * quad = p4est_quadrant_array_index(tquadrants, q);
//             p4est_qcoord_t length = P4EST_QUADRANT_LEN(quad->level);

//             // Loop over the four corners of the quadrant
//             for(int n = 0; n < 4; ++n){
//                 int k = q - Q + nodeIncrements[treeid - p4est->first_local_tree];
//                 long lidx = nodes->element_nodes[4*k+n];

//                 double lx = length * ((int) (n % 2) == 1);
//                 double ly = length * ((int) (n / 2) == 1);
//                 double xy[3];
//                 p4est_qcoord_to_vertex(p4est->connectivity, treeid, quad->x+lx, quad->y+ly, xy);

//                 if( (n == 0) 
//                   || isHangingNode(nodes->face_code[k], n)
//                   || isUpperBoundaryNode(quad, n, treeid, length) 
//                 )
//                 {
//                     long lni = NodeIDs.find(std::make_pair(xy[0],xy[1]))->second;
//                     double * point = arg.getSampleDataRW(lni);
//                     point[0] = xy[0];
//                     point[1] = xy[1];
//                 }
//             }
//         }
//     }     
}

#ifdef ESYS_HAVE_TRILINOS
//protected
esys_trilinos::const_TrilinosGraph_ptr Brick::getTrilinosGraph() const
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
    // TODO

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

// bool Brick::operator==(const AbstractDomain& other) const
// {
//     const Brick* o=dynamic_cast<const Brick*>(&other);
//     if (o) {
//         return (p8est_checksum(p8est) == p8est_checksum(o->p8est)
//             && forestData == o->forestData);
//     }
//     return false;
// }

//private
template<typename Scalar>
void Brick::addToMatrixAndRHS(escript::AbstractSystemMatrix* S, escript::Data& F,
         const std::vector<Scalar>& EM_S, const std::vector<Scalar>& EM_F, 
         bool addS, bool addF, index_t e, index_t t, int nEq, int nComp) const
{    
    IndexVector rowIndex(4);
    p8est_tree_t * currenttree = p8est_tree_array_index(p8est->trees, t);
    sc_array_t * tquadrants = &currenttree->quadrants;
    // p4est_locidx_t Q = (p4est_locidx_t) tquadrants->elem_count;
    p8est_quadrant_t * quad = p8est_quadrant_array_index(tquadrants, e);
    p4est_qcoord_t length = P8EST_QUADRANT_LEN(quad->level);
    for(int i = 0; i < 4; i++)
    {
        double lx = length * ((int) (i % 2) == 1);
        double ly = length * ((int) (i / 2) == 1);
        double lz = length * ((int) (i / 2) == 1);
        double xyz[3];
        p8est_qcoord_to_vertex(p8est->connectivity, t, quad->x+lx, quad->y+ly, quad->z+lz, xyz);
        rowIndex[i] = NodeIDs.find(std::make_tuple(xyz[0],xyz[1],xyz[2]))->second;
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


template
void Brick::addToMatrixAndRHS<real_t>(escript::AbstractSystemMatrix* S, escript::Data& F,
         const std::vector<real_t>& EM_S, const std::vector<real_t>& EM_F, 
         bool addS, bool addF, index_t e, index_t t, int nEq, int nComp) const;

template
void Brick::addToMatrixAndRHS<cplx_t>(escript::AbstractSystemMatrix* S, escript::Data& F,
         const std::vector<cplx_t>& EM_S, const std::vector<cplx_t>& EM_F, 
         bool addS, bool addF, index_t e, index_t t, int nEq, int nComp) const;

//protected
void Brick::interpolateNodesOnElements(escript::Data& out,
                                       const escript::Data& in,
                                       bool reduced) const
{
    if (out.isComplex()!=in.isComplex())
    {
        throw OxleyException("Programmer Error: in and out parameters do not have the same complexity.");   
    }
    if (out.isComplex())
    {
        interpolateNodesOnElementsWorker(out, in, reduced, escript::DataTypes::cplx_t(0));    
    }
    else
    {
        interpolateNodesOnElementsWorker(out, in, reduced, escript::DataTypes::real_t(0));    
    }  
}
//protected
void Brick::interpolateNodesOnFaces(escript::Data& out, const escript::Data& in,
                                    bool reduced) const
{
    //todo:
    throw OxleyException("interpolateNodesOnFaces");
    
    // if (out.isComplex()!=in.isComplex())
    // {
    //     throw OxleyException("Programmer Error: in and out parameters do not have the same complexity.");   
    // }
    // if (out.isComplex())
    // {
    //     interpolateNodesOnFacesWorker(out, in, reduced, escript::DataTypes::cplx_t(0));    
    // }
    // else
    // {
    //     interpolateNodesOnFacesWorker(out, in, reduced, escript::DataTypes::real_t(0));    
    // }      
}

//private
template <typename S>
void Brick::interpolateNodesOnFacesWorker(escript::Data& out, const escript::Data& in,
                                    bool reduced, S sentinel) const
{
    //todo:
    throw OxleyException("interpolateNodesOnFacesWorker");
//     const dim_t numComp = in.getDataPointSize();
//     const long  numNodes = getNumNodes();
//     if (reduced) {
//         out.requireWrite();
//         const S c0 = 0.25;
//         double * fxx = new double[4*numComp*numNodes];


//         // This structure is used to store info needed by p4est
//         interpolateNodesOnElementsWorker_Data<S> interpolateData;
//         interpolateData.fxx = fxx;
//         interpolateData.sentinel = sentinel;
//         interpolateData.offset = numComp*sizeof(S);

//         p8est_iterate(p8est, NULL, &interpolateData, get_interpolateNodesOnElementWorker_data, NULL, NULL);
//         for(p8est_topidx_t treeid = p8est->first_local_tree; treeid <= p8est->last_local_tree; treeid++)
//         {
//             p4est_tree_t * currenttree = p8est_tree_array_index(p8est->trees, treeid);
//             sc_array_t * tquadrants = &currenttree->quadrants;
//             p4est_qcoord_t Q = (p4est_qcoord_t) tquadrants->elem_count;
// #pragma omp parallel for
//             for(p4est_qcoord_t e = nodes->global_offset; e < Q+nodes->global_offset; e++)
//             {
//                 S* o = out.getSampleDataRW(e, sentinel);

//             if (m_faceOffset[0] > -1) {
//                 for (index_t k2=0; k2 < m_NE[2]; ++k2) {
//                     for (index_t k1=0; k1 < m_NE[1]; ++k1) {
//                         memcpy(&f_000[0], in.getSampleDataRO(INDEX3(0,k1,k2, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
//                         memcpy(&f_001[0], in.getSampleDataRO(INDEX3(0,k1,k2+1, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
//                         memcpy(&f_010[0], in.getSampleDataRO(INDEX3(0,k1+1,k2, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
//                         memcpy(&f_011[0], in.getSampleDataRO(INDEX3(0,k1+1,k2+1, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
//                         S* o = out.getSampleDataRW(m_faceOffset[0]+INDEX2(k1,k2,m_NE[1]), sentinel);
//                         for (index_t i=0; i < numComp; ++i) {
//                             o[INDEX2(i,numComp,0)] = (f_000[i] + f_001[i] + f_010[i] + f_011[i])/static_cast<S>(4);
//                         } // end of component loop i
//                     } // end of k1 loop
//                 } // end of k2 loop
//             } // end of face 0
//             if (m_faceOffset[1] > -1) {
//                 for (index_t k2=0; k2 < m_NE[2]; ++k2) {
//                     for (index_t k1=0; k1 < m_NE[1]; ++k1) {
//                         memcpy(&f_100[0], in.getSampleDataRO(INDEX3(m_NN[0]-1,k1,k2, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
//                         memcpy(&f_101[0], in.getSampleDataRO(INDEX3(m_NN[0]-1,k1,k2+1, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
//                         memcpy(&f_110[0], in.getSampleDataRO(INDEX3(m_NN[0]-1,k1+1,k2, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
//                         memcpy(&f_111[0], in.getSampleDataRO(INDEX3(m_NN[0]-1,k1+1,k2+1, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
//                         S* o = out.getSampleDataRW(m_faceOffset[1]+INDEX2(k1,k2,m_NE[1]), sentinel);
//                         for (index_t i=0; i < numComp; ++i) {
//                             o[INDEX2(i,numComp,0)] = (f_100[i] + f_101[i] + f_110[i] + f_111[i])/static_cast<S>(4);
//                         } // end of component loop i
//                     } // end of k1 loop
//                 } // end of k2 loop
//             } // end of face 1
//             if (m_faceOffset[2] > -1) {
//                 for (index_t k2=0; k2 < m_NE[2]; ++k2) {
//                     for (index_t k0=0; k0 < m_NE[0]; ++k0) {
//                         memcpy(&f_000[0], in.getSampleDataRO(INDEX3(k0,0,k2, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
//                         memcpy(&f_001[0], in.getSampleDataRO(INDEX3(k0,0,k2+1, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
//                         memcpy(&f_100[0], in.getSampleDataRO(INDEX3(k0+1,0,k2, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
//                         memcpy(&f_101[0], in.getSampleDataRO(INDEX3(k0+1,0,k2+1, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
//                         S* o = out.getSampleDataRW(m_faceOffset[2]+INDEX2(k0,k2,m_NE[0]), sentinel);
//                         for (index_t i=0; i < numComp; ++i) {
//                             o[INDEX2(i,numComp,0)] = (f_000[i] + f_001[i] + f_100[i] + f_101[i])/static_cast<S>(4);
//                         } // end of component loop i
//                     } // end of k0 loop
//                 } // end of k2 loop
//             } // end of face 2
//             if (m_faceOffset[3] > -1) {
//                 for (index_t k2=0; k2 < m_NE[2]; ++k2) {
//                     for (index_t k0=0; k0 < m_NE[0]; ++k0) {
//                         memcpy(&f_010[0], in.getSampleDataRO(INDEX3(k0,m_NN[1]-1,k2, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
//                         memcpy(&f_011[0], in.getSampleDataRO(INDEX3(k0,m_NN[1]-1,k2+1, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
//                         memcpy(&f_110[0], in.getSampleDataRO(INDEX3(k0+1,m_NN[1]-1,k2, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
//                         memcpy(&f_111[0], in.getSampleDataRO(INDEX3(k0+1,m_NN[1]-1,k2+1, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
//                         S* o = out.getSampleDataRW(m_faceOffset[3]+INDEX2(k0,k2,m_NE[0]), sentinel);
//                         for (index_t i=0; i < numComp; ++i) {
//                             o[INDEX2(i,numComp,0)] = (f_010[i] + f_011[i] + f_110[i] + f_111[i])/static_cast<S>(4);
//                         } // end of component loop i
//                     } // end of k0 loop
//                 } // end of k2 loop
//             } // end of face 3
//             if (m_faceOffset[4] > -1) {
//                 for (index_t k1=0; k1 < m_NE[1]; ++k1) {
//                     for (index_t k0=0; k0 < m_NE[0]; ++k0) {
//                         memcpy(&f_000[0], in.getSampleDataRO(INDEX3(k0,k1,0, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
//                         memcpy(&f_010[0], in.getSampleDataRO(INDEX3(k0,k1+1,0, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
//                         memcpy(&f_100[0], in.getSampleDataRO(INDEX3(k0+1,k1,0, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
//                         memcpy(&f_110[0], in.getSampleDataRO(INDEX3(k0+1,k1+1,0, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
//                         S* o = out.getSampleDataRW(m_faceOffset[4]+INDEX2(k0,k1,m_NE[0]), sentinel);
//                         for (index_t i=0; i < numComp; ++i) {
//                             o[INDEX2(i,numComp,0)] = (f_000[i] + f_010[i] + f_100[i] + f_110[i])/static_cast<S>(4);
//                         } // end of component loop i
//                     } // end of k0 loop
//                 } // end of k1 loop
//             } // end of face 4
//             if (m_faceOffset[5] > -1) {
//                 for (index_t k1=0; k1 < m_NE[1]; ++k1) {
//                     for (index_t k0=0; k0 < m_NE[0]; ++k0) {
//                         memcpy(&f_001[0], in.getSampleDataRO(INDEX3(k0,k1,m_NN[2]-1, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
//                         memcpy(&f_011[0], in.getSampleDataRO(INDEX3(k0,k1+1,m_NN[2]-1, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
//                         memcpy(&f_101[0], in.getSampleDataRO(INDEX3(k0+1,k1,m_NN[2]-1, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
//                         memcpy(&f_111[0], in.getSampleDataRO(INDEX3(k0+1,k1+1,m_NN[2]-1, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
//                         S* o = out.getSampleDataRW(m_faceOffset[5]+INDEX2(k0,k1,m_NE[0]), sentinel);
//                         for (index_t i=0; i < numComp; ++i) {
//                             o[INDEX2(i,numComp,0)] = (f_001[i] + f_011[i] + f_101[i] + f_111[i])/static_cast<S>(4);
//                         } // end of component loop i
//                     } // end of k0 loop
//                 } // end of k1 loop
//             } // end of face 5
//     } else {
//         out.requireWrite();
//         const S c0 = 0.044658198738520451079;
//         const S c1 = 0.16666666666666666667;
//         const S c2 = 0.62200846792814621559;
// #pragma omp parallel
//         {
//             vector<S> f_000(numComp);
//             vector<S> f_001(numComp);
//             vector<S> f_010(numComp);
//             vector<S> f_011(numComp);
//             vector<S> f_100(numComp);
//             vector<S> f_101(numComp);
//             vector<S> f_110(numComp);
//             vector<S> f_111(numComp);
//             if (m_faceOffset[0] > -1) {
// #pragma omp for nowait
//                 for (index_t k2=0; k2 < m_NE[2]; ++k2) {
//                     for (index_t k1=0; k1 < m_NE[1]; ++k1) {
//                         memcpy(&f_000[0], in.getSampleDataRO(INDEX3(0,k1,k2, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
//                         memcpy(&f_001[0], in.getSampleDataRO(INDEX3(0,k1,k2+1, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
//                         memcpy(&f_010[0], in.getSampleDataRO(INDEX3(0,k1+1,k2, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
//                         memcpy(&f_011[0], in.getSampleDataRO(INDEX3(0,k1+1,k2+1, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
//                         S* o = out.getSampleDataRW(m_faceOffset[0]+INDEX2(k1,k2,m_NE[1]), sentinel);
//                         for (index_t i=0; i < numComp; ++i) {
//                             o[INDEX2(i,numComp,0)] = f_000[i]*c2 + f_011[i]*c0 + c1*(f_001[i] + f_010[i]);
//                             o[INDEX2(i,numComp,1)] = f_001[i]*c0 + f_010[i]*c2 + c1*(f_000[i] + f_011[i]);
//                             o[INDEX2(i,numComp,2)] = f_001[i]*c2 + f_010[i]*c0 + c1*(f_000[i] + f_011[i]);
//                             o[INDEX2(i,numComp,3)] = f_000[i]*c0 + f_011[i]*c2 + c1*(f_001[i] + f_010[i]);
//                         } // end of component loop i
//                     } // end of k1 loop
//                 } // end of k2 loop
//             } // end of face 0
//             if (m_faceOffset[1] > -1) {
// #pragma omp for nowait
//                 for (index_t k2=0; k2 < m_NE[2]; ++k2) {
//                     for (index_t k1=0; k1 < m_NE[1]; ++k1) {
//                         memcpy(&f_100[0], in.getSampleDataRO(INDEX3(m_NN[0]-1,k1,k2, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
//                         memcpy(&f_101[0], in.getSampleDataRO(INDEX3(m_NN[0]-1,k1,k2+1, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
//                         memcpy(&f_110[0], in.getSampleDataRO(INDEX3(m_NN[0]-1,k1+1,k2, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
//                         memcpy(&f_111[0], in.getSampleDataRO(INDEX3(m_NN[0]-1,k1+1,k2+1, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
//                         S* o = out.getSampleDataRW(m_faceOffset[1]+INDEX2(k1,k2,m_NE[1]), sentinel);
//                         for (index_t i=0; i < numComp; ++i) {
//                             o[INDEX2(i,numComp,0)] = f_100[i]*c2 + f_111[i]*c0 + c1*(f_101[i] + f_110[i]);
//                             o[INDEX2(i,numComp,1)] = f_101[i]*c0 + f_110[i]*c2 + c1*(f_100[i] + f_111[i]);
//                             o[INDEX2(i,numComp,2)] = f_101[i]*c2 + f_110[i]*c0 + c1*(f_100[i] + f_111[i]);
//                             o[INDEX2(i,numComp,3)] = f_100[i]*c0 + f_111[i]*c2 + c1*(f_101[i] + f_110[i]);
//                         } // end of component loop i
//                     } // end of k1 loop
//                 } // end of k2 loop
//             } // end of face 1
//             if (m_faceOffset[2] > -1) {
// #pragma omp for nowait
//                 for (index_t k2=0; k2 < m_NE[2]; ++k2) {
//                     for (index_t k0=0; k0 < m_NE[0]; ++k0) {
//                         memcpy(&f_000[0], in.getSampleDataRO(INDEX3(k0,0,k2, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
//                         memcpy(&f_001[0], in.getSampleDataRO(INDEX3(k0,0,k2+1, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
//                         memcpy(&f_100[0], in.getSampleDataRO(INDEX3(k0+1,0,k2, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
//                         memcpy(&f_101[0], in.getSampleDataRO(INDEX3(k0+1,0,k2+1, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
//                         S* o = out.getSampleDataRW(m_faceOffset[2]+INDEX2(k0,k2,m_NE[0]), sentinel);
//                         for (index_t i=0; i < numComp; ++i) {
//                             o[INDEX2(i,numComp,0)] = f_000[i]*c2 + f_101[i]*c0 + c1*(f_001[i] + f_100[i]);
//                             o[INDEX2(i,numComp,1)] = f_001[i]*c0 + f_100[i]*c2 + c1*(f_000[i] + f_101[i]);
//                             o[INDEX2(i,numComp,2)] = f_001[i]*c2 + f_100[i]*c0 + c1*(f_000[i] + f_101[i]);
//                             o[INDEX2(i,numComp,3)] = f_000[i]*c0 + f_101[i]*c2 + c1*(f_001[i] + f_100[i]);
//                         } // end of component loop i
//                     } // end of k0 loop
//                 } // end of k2 loop
//             } // end of face 2
//             if (m_faceOffset[3] > -1) {
// #pragma omp for nowait
//                 for (index_t k2=0; k2 < m_NE[2]; ++k2) {
//                     for (index_t k0=0; k0 < m_NE[0]; ++k0) {
//                         memcpy(&f_010[0], in.getSampleDataRO(INDEX3(k0,m_NN[1]-1,k2, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
//                         memcpy(&f_011[0], in.getSampleDataRO(INDEX3(k0,m_NN[1]-1,k2+1, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
//                         memcpy(&f_110[0], in.getSampleDataRO(INDEX3(k0+1,m_NN[1]-1,k2, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
//                         memcpy(&f_111[0], in.getSampleDataRO(INDEX3(k0+1,m_NN[1]-1,k2+1, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
//                         S* o = out.getSampleDataRW(m_faceOffset[3]+INDEX2(k0,k2,m_NE[0]), sentinel);
//                         for (index_t i=0; i < numComp; ++i) {
//                             o[INDEX2(i,numComp,0)] = f_010[i]*c2 + f_111[i]*c0 + c1*(f_011[i] + f_110[i]);
//                             o[INDEX2(i,numComp,1)] = f_011[i]*c0 + f_110[i]*c2 + c1*(f_010[i] + f_111[i]);
//                             o[INDEX2(i,numComp,2)] = f_011[i]*c2 + f_110[i]*c0 + c1*(f_010[i] + f_111[i]);
//                             o[INDEX2(i,numComp,3)] = f_010[i]*c0 + f_111[i]*c2 + c1*(f_011[i] + f_110[i]);
//                         } // end of component loop i
//                     } // end of k0 loop
//                 } // end of k2 loop
//             } // end of face 3
//             if (m_faceOffset[4] > -1) {
// #pragma omp for nowait
//                 for (index_t k1=0; k1 < m_NE[1]; ++k1) {
//                     for (index_t k0=0; k0 < m_NE[0]; ++k0) {
//                         memcpy(&f_000[0], in.getSampleDataRO(INDEX3(k0,k1,0, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
//                         memcpy(&f_010[0], in.getSampleDataRO(INDEX3(k0,k1+1,0, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
//                         memcpy(&f_100[0], in.getSampleDataRO(INDEX3(k0+1,k1,0, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
//                         memcpy(&f_110[0], in.getSampleDataRO(INDEX3(k0+1,k1+1,0, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
//                         S* o = out.getSampleDataRW(m_faceOffset[4]+INDEX2(k0,k1,m_NE[0]), sentinel);
//                         for (index_t i=0; i < numComp; ++i) {
//                             o[INDEX2(i,numComp,0)] = f_000[i]*c2 + f_110[i]*c0 + c1*(f_010[i] + f_100[i]);
//                             o[INDEX2(i,numComp,1)] = f_010[i]*c0 + f_100[i]*c2 + c1*(f_000[i] + f_110[i]);
//                             o[INDEX2(i,numComp,2)] = f_010[i]*c2 + f_100[i]*c0 + c1*(f_000[i] + f_110[i]);
//                             o[INDEX2(i,numComp,3)] = f_000[i]*c0 + f_110[i]*c2 + c1*(f_010[i] + f_100[i]);
//                         } // end of component loop i
//                     } // end of k0 loop
//                 } // end of k1 loop
//             } // end of face 4
//             if (m_faceOffset[5] > -1) {
// #pragma omp for nowait
//                 for (index_t k1=0; k1 < m_NE[1]; ++k1) {
//                     for (index_t k0=0; k0 < m_NE[0]; ++k0) {
//                         memcpy(&f_001[0], in.getSampleDataRO(INDEX3(k0,k1,m_NN[2]-1, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
//                         memcpy(&f_011[0], in.getSampleDataRO(INDEX3(k0,k1+1,m_NN[2]-1, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
//                         memcpy(&f_101[0], in.getSampleDataRO(INDEX3(k0+1,k1,m_NN[2]-1, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
//                         memcpy(&f_111[0], in.getSampleDataRO(INDEX3(k0+1,k1+1,m_NN[2]-1, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
//                         S* o = out.getSampleDataRW(m_faceOffset[5]+INDEX2(k0,k1,m_NE[0]), sentinel);
//                         for (index_t i=0; i < numComp; ++i) {
//                             o[INDEX2(i,numComp,0)] = f_001[i]*c2 + f_111[i]*c0 + c1*(f_011[i] + f_101[i]);
//                             o[INDEX2(i,numComp,1)] = f_011[i]*c0 + f_101[i]*c2 + c1*(f_001[i] + f_111[i]);
//                             o[INDEX2(i,numComp,2)] = f_011[i]*c2 + f_101[i]*c0 + c1*(f_001[i] + f_111[i]);
//                             o[INDEX2(i,numComp,3)] = f_001[i]*c0 + f_111[i]*c2 + c1*(f_011[i] + f_101[i]);
//                         } // end of component loop i
//                     } // end of k0 loop
//                 } // end of k1 loop
//             } // end of face 5
//         } // end of parallel section
//     }
}


//private
template <typename S>
void Brick::interpolateNodesOnElementsWorker(escript::Data& out,
                                       const escript::Data& in,
                                       bool reduced, S sentinel) const                                       
{
    const dim_t numComp = in.getDataPointSize();
    const long numNodes = getNumNodes();
    if (reduced) {
        out.requireWrite();
        const S c0 = 0.125;  
        double * fxxx = new double[8*numComp*numNodes];

        // This structure is used to store info needed by p4est
        interpolateNodesOnElementsWorker_Data<S> interpolateData;
        interpolateData.fxx = fxxx;
        interpolateData.sentinel = sentinel;
        interpolateData.offset = numComp*sizeof(S);

        p8est_iterate(p8est, NULL, &interpolateData, get_interpolateNodesOnElementWorker_data, NULL, NULL, NULL);
        for(p8est_topidx_t treeid = p8est->first_local_tree; treeid <= p8est->last_local_tree; treeid++)
        {
            p8est_tree_t * currenttree = p8est_tree_array_index(p8est->trees, treeid);
            sc_array_t * tquadrants = &currenttree->quadrants;
            p4est_qcoord_t Q = (p4est_locidx_t) tquadrants->elem_count;
#pragma omp parallel for
            for(p4est_qcoord_t e = nodes->global_offset; e < Q+nodes->global_offset; e++)
            {
                S* o = out.getSampleDataRW(e, sentinel);
                for (index_t i=0; i < numComp; ++i) {
                    o[i] = c0 * (fxxx[INDEX3(0,i,e,numComp,numNodes)] + 
                                 fxxx[INDEX3(1,i,e,numComp,numNodes)] + 
                                 fxxx[INDEX3(2,i,e,numComp,numNodes)] + 
                                 fxxx[INDEX3(3,i,e,numComp,numNodes)] + 
                                 fxxx[INDEX3(4,i,e,numComp,numNodes)] + 
                                 fxxx[INDEX3(5,i,e,numComp,numNodes)] + 
                                 fxxx[INDEX3(6,i,e,numComp,numNodes)] + 
                                 fxxx[INDEX3(7,i,e,numComp,numNodes)] );
                }
            }
        }
        delete[] fxxx;
    } else {
        out.requireWrite();
        const S c0 = .0094373878376559314545;
        const S c1 = .035220810900864519624;
        const S c2 = .13144585576580214704;
        const S c3 = .49056261216234406855;

        double * fxxx = new double[8*numComp*numNodes];
        // This structure is used to store info needed by p4est
        interpolateNodesOnElementsWorker_Data<S> interpolateData;
        interpolateData.fxx = fxxx;
        interpolateData.sentinel = sentinel;
        interpolateData.offset = numComp*sizeof(S);

        p8est_iterate(p8est, NULL, &interpolateData, get_interpolateNodesOnElementWorker_data, NULL, NULL, NULL);
        for(p8est_topidx_t treeid = p8est->first_local_tree; treeid <= p8est->last_local_tree; treeid++)
        {
            p8est_tree_t * currenttree = p8est_tree_array_index(p8est->trees, treeid);
            sc_array_t * tquadrants = &currenttree->quadrants;
            p4est_qcoord_t Q = (p4est_locidx_t) tquadrants->elem_count;
#pragma omp parallel for
            for(p4est_qcoord_t e = nodes->global_offset; e < Q+nodes->global_offset; e++)
            {
                S* o = out.getSampleDataRW(e, sentinel);
                for (index_t i=0; i < numComp; ++i) {
                    o[INDEX2(i,numComp,0)] = c3*fxxx[INDEX3(0,i,e,numComp,numNodes)] + 
                                             c0*fxxx[INDEX3(7,i,e,numComp,numNodes)] + 
                                             c2*(fxxx[INDEX3(1,i,e,numComp,numNodes)] +
                                                 fxxx[INDEX3(2,i,e,numComp,numNodes)] +  
                                                 fxxx[INDEX3(4,i,e,numComp,numNodes)]) + 
                                             c1*(fxxx[INDEX3(3,i,e,numComp,numNodes)] + 
                                                 fxxx[INDEX3(5,i,e,numComp,numNodes)] +  
                                                 fxxx[INDEX3(6,i,e,numComp,numNodes)]);
                    o[INDEX2(i,numComp,1)] = fxxx[INDEX3(3,i,e,numComp,numNodes)]*c0 + 
                                             fxxx[INDEX3(4,i,e,numComp,numNodes)]*c3 + 
                                             c2*(fxxx[INDEX3(0,i,e,numComp,numNodes)] + 
                                                 fxxx[INDEX3(5,i,e,numComp,numNodes)] + 
                                                 fxxx[INDEX3(6,i,e,numComp,numNodes)]) + 
                                             c1*(fxxx[INDEX3(1,i,e,numComp,numNodes)] + 
                                                 fxxx[INDEX3(2,i,e,numComp,numNodes)] + 
                                                 fxxx[INDEX3(7,i,e,numComp,numNodes)]);
                    o[INDEX2(i,numComp,2)] = fxxx[INDEX3(2,i,e,numComp,numNodes)]*c3 + 
                                             fxxx[INDEX3(5,i,e,numComp,numNodes)]*c0 + 
                                             c2*(fxxx[INDEX3(0,i,e,numComp,numNodes)] + 
                                                 fxxx[INDEX3(3,i,e,numComp,numNodes)] + 
                                                 fxxx[INDEX3(6,i,e,numComp,numNodes)]) + 
                                             c1*(fxxx[INDEX3(1,i,e,numComp,numNodes)] + 
                                                 fxxx[INDEX3(4,i,e,numComp,numNodes)] + 
                                                 fxxx[INDEX3(7,i,e,numComp,numNodes)]);
                    o[INDEX2(i,numComp,3)] = fxxx[INDEX3(1,i,e,numComp,numNodes)]*c0 + 
                                             fxxx[INDEX3(6,i,e,numComp,numNodes)]*c3 + 
                                             c2*(fxxx[INDEX3(2,i,e,numComp,numNodes)] + 
                                                 fxxx[INDEX3(4,i,e,numComp,numNodes)] + 
                                                 fxxx[INDEX3(7,i,e,numComp,numNodes)]) + 
                                             c1*(fxxx[INDEX3(0,i,e,numComp,numNodes)] + 
                                                 fxxx[INDEX3(3,i,e,numComp,numNodes)] + 
                                                 fxxx[INDEX3(5,i,e,numComp,numNodes)]);
                    o[INDEX2(i,numComp,4)] = fxxx[INDEX3(1,i,e,numComp,numNodes)]*c3 + 
                                             fxxx[INDEX3(6,i,e,numComp,numNodes)]*c0 + 
                                             c2*(fxxx[INDEX3(0,i,e,numComp,numNodes)] + 
                                                 fxxx[INDEX3(3,i,e,numComp,numNodes)] + 
                                                 fxxx[INDEX3(5,i,e,numComp,numNodes)]) + 
                                             c1*(fxxx[INDEX3(2,i,e,numComp,numNodes)] + 
                                                 fxxx[INDEX3(4,i,e,numComp,numNodes)] + 
                                                 fxxx[INDEX3(7,i,e,numComp,numNodes)]);
                    o[INDEX2(i,numComp,5)] = fxxx[INDEX3(2,i,e,numComp,numNodes)]*c0 + 
                                             fxxx[INDEX3(5,i,e,numComp,numNodes)]*c3 + 
                                             c2*(fxxx[INDEX3(1,i,e,numComp,numNodes)] + 
                                                 fxxx[INDEX3(4,i,e,numComp,numNodes)] + 
                                                 fxxx[INDEX3(7,i,e,numComp,numNodes)]) + 
                                             c1*(fxxx[INDEX3(0,i,e,numComp,numNodes)] + 
                                                 fxxx[INDEX3(3,i,e,numComp,numNodes)] + 
                                                 fxxx[INDEX3(6,i,e,numComp,numNodes)]);
                    o[INDEX2(i,numComp,6)] = fxxx[INDEX3(3,i,e,numComp,numNodes)]*c3 + 
                                             fxxx[INDEX3(4,i,e,numComp,numNodes)]*c0 + 
                                             c2*(fxxx[INDEX3(1,i,e,numComp,numNodes)] + 
                                                 fxxx[INDEX3(2,i,e,numComp,numNodes)] + 
                                                 fxxx[INDEX3(7,i,e,numComp,numNodes)]) + 
                                             c1*(fxxx[INDEX3(0,i,e,numComp,numNodes)] + 
                                                 fxxx[INDEX3(5,i,e,numComp,numNodes)] + 
                                                 fxxx[INDEX3(6,i,e,numComp,numNodes)]);
                    o[INDEX2(i,numComp,7)] = fxxx[INDEX3(0,i,e,numComp,numNodes)]*c0 + 
                                             fxxx[INDEX3(7,i,e,numComp,numNodes)]*c3 + 
                                             c2*(fxxx[INDEX3(3,i,e,numComp,numNodes)] + 
                                                 fxxx[INDEX3(5,i,e,numComp,numNodes)] + 
                                                 fxxx[INDEX3(6,i,e,numComp,numNodes)]) + 
                                             c1*(fxxx[INDEX3(1,i,e,numComp,numNodes)] + 
                                                 fxxx[INDEX3(2,i,e,numComp,numNodes)] + 
                                                 fxxx[INDEX3(4,i,e,numComp,numNodes)]);
                    
                }
            }
        }
    }
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

    // Find the maximum level of refinement on the mesh
    int max_level = 0;
    for(p8est_topidx_t tree = p8est->first_local_tree; tree < p8est->last_local_tree; tree++) {
        p8est_tree_t * tree_t = p8est_tree_array_index(p8est->trees, tree);
        max_level = SC_MAX(max_level, tree_t->maxlevel);
    }

    // double cx[7][P4EST_MAXLEVEL] = {{0}};
    // double cy[7][P4EST_MAXLEVEL] = {{0}};

    const double C0 = .044658198738520451079;
    const double C1 = .16666666666666666667;
    const double C2 = .21132486540518711775;
    const double C3 = .25;
    const double C4 = .5;
    const double C5 = .62200846792814621559;
    const double C6 = .78867513459481288225;
    // const dim_t NE0 = m_NE[0];
    // const dim_t NE1 = m_NE[1];
    // const dim_t NE2 = m_NE[2];
    const Scalar zero = static_cast<Scalar>(0);

    if (out.getFunctionSpace().getTypeCode() == Elements) 
    {
        out.requireWrite();

        std::vector<Scalar> f_000(numComp, zero);
        std::vector<Scalar> f_001(numComp, zero);
        std::vector<Scalar> f_010(numComp, zero);
        std::vector<Scalar> f_011(numComp, zero);
        std::vector<Scalar> f_100(numComp, zero);
        std::vector<Scalar> f_101(numComp, zero);
        std::vector<Scalar> f_110(numComp, zero);
        std::vector<Scalar> f_111(numComp, zero);

        for(p8est_topidx_t t = p8est->first_local_tree; t <= p8est->last_local_tree; t++) // Loop over every tree
        {
            p4est_tree_t * currenttree = p4est_tree_array_index(p8est->trees, t);
            sc_array_t * tquadrants = &currenttree->quadrants;
            p4est_qcoord_t Q = (p4est_locidx_t) tquadrants->elem_count;
#pragma omp parallel for
            for(p4est_qcoord_t e = nodes->global_offset; e < Q+nodes->global_offset; e++) // Loop over every quadrant within the tree
            {
                // Work out what level this element is on 
                p8est_quadrant_t * quad = p8est_quadrant_array_index(tquadrants, e);
                // octantData * quaddata = (octantData *) quad->p.user_data;
                int l = quad->level;

                memcpy(&f_000[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                memcpy(&f_001[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                memcpy(&f_010[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                memcpy(&f_011[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                memcpy(&f_100[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                memcpy(&f_101[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                memcpy(&f_110[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                memcpy(&f_111[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));

                Scalar* o = out.getSampleDataRW(e, zero);

                for (index_t i=0; i < numComp; ++i) {
                    const Scalar V0 =((f_100[i]-f_000[i])*C5 + (f_111[i]-f_011[i])*C0 + (f_101[i]+f_110[i]-f_001[i]-f_010[i])*C1) / forestData.m_dx[0][l];
                    const Scalar V1 =((f_110[i]-f_010[i])*C5 + (f_101[i]-f_001[i])*C0 + (f_100[i]+f_111[i]-f_000[i]-f_011[i])*C1) / forestData.m_dx[0][l];
                    const Scalar V2 =((f_101[i]-f_001[i])*C5 + (f_110[i]-f_010[i])*C0 + (f_100[i]+f_111[i]-f_000[i]-f_011[i])*C1) / forestData.m_dx[0][l];
                    const Scalar V3 =((f_111[i]-f_011[i])*C5 + (f_100[i]-f_000[i])*C0 + (f_101[i]+f_110[i]-f_001[i]-f_010[i])*C1) / forestData.m_dx[0][l];
                    const Scalar V4 =((f_010[i]-f_000[i])*C5 + (f_111[i]-f_101[i])*C0 + (f_011[i]+f_110[i]-f_001[i]-f_100[i])*C1) / forestData.m_dx[1][l];
                    const Scalar V5 =((f_110[i]-f_100[i])*C5 + (f_011[i]-f_001[i])*C0 + (f_010[i]+f_111[i]-f_000[i]-f_101[i])*C1) / forestData.m_dx[1][l];
                    const Scalar V6 =((f_011[i]-f_001[i])*C5 + (f_110[i]-f_100[i])*C0 + (f_010[i]+f_111[i]-f_000[i]-f_101[i])*C1) / forestData.m_dx[1][l];
                    const Scalar V7 =((f_111[i]-f_101[i])*C5 + (f_010[i]-f_000[i])*C0 + (f_011[i]+f_110[i]-f_001[i]-f_100[i])*C1) / forestData.m_dx[1][l];
                    const Scalar V8 =((f_001[i]-f_000[i])*C5 + (f_111[i]-f_110[i])*C0 + (f_011[i]+f_101[i]-f_010[i]-f_100[i])*C1) / forestData.m_dx[2][l];
                    const Scalar V9 =((f_101[i]-f_100[i])*C5 + (f_011[i]-f_010[i])*C0 + (f_001[i]+f_111[i]-f_000[i]-f_110[i])*C1) / forestData.m_dx[2][l];
                    const Scalar V10=((f_011[i]-f_010[i])*C5 + (f_101[i]-f_100[i])*C0 + (f_001[i]+f_111[i]-f_000[i]-f_110[i])*C1) / forestData.m_dx[2][l];
                    const Scalar V11=((f_111[i]-f_110[i])*C5 + (f_001[i]-f_000[i])*C0 + (f_011[i]+f_101[i]-f_010[i]-f_100[i])*C1) / forestData.m_dx[2][l];

                    o[INDEX3(i,0,0,numComp,3)] = V0;
                    o[INDEX3(i,1,0,numComp,3)] = V4;
                    o[INDEX3(i,2,0,numComp,3)] = V8;
                    o[INDEX3(i,0,1,numComp,3)] = V0;
                    o[INDEX3(i,1,1,numComp,3)] = V5;
                    o[INDEX3(i,2,1,numComp,3)] = V9;
                    o[INDEX3(i,0,2,numComp,3)] = V1;
                    o[INDEX3(i,1,2,numComp,3)] = V4;
                    o[INDEX3(i,2,2,numComp,3)] = V10;
                    o[INDEX3(i,0,3,numComp,3)] = V1;
                    o[INDEX3(i,1,3,numComp,3)] = V5;
                    o[INDEX3(i,2,3,numComp,3)] = V11;
                    o[INDEX3(i,0,4,numComp,3)] = V2;
                    o[INDEX3(i,1,4,numComp,3)] = V6;
                    o[INDEX3(i,2,4,numComp,3)] = V8;
                    o[INDEX3(i,0,5,numComp,3)] = V2;
                    o[INDEX3(i,1,5,numComp,3)] = V7;
                    o[INDEX3(i,2,5,numComp,3)] = V9;
                    o[INDEX3(i,0,6,numComp,3)] = V3;
                    o[INDEX3(i,1,6,numComp,3)] = V6;
                    o[INDEX3(i,2,6,numComp,3)] = V10;
                    o[INDEX3(i,0,7,numComp,3)] = V3;
                    o[INDEX3(i,1,7,numComp,3)] = V7;
                    o[INDEX3(i,2,7,numComp,3)] = V11;
                }
            }
        }
    } 
    else if (out.getFunctionSpace().getTypeCode() == ReducedElements) 
    {
        out.requireWrite();

        std::vector<Scalar> f_000(numComp, zero);
        std::vector<Scalar> f_001(numComp, zero);
        std::vector<Scalar> f_010(numComp, zero);
        std::vector<Scalar> f_011(numComp, zero);
        std::vector<Scalar> f_100(numComp, zero);
        std::vector<Scalar> f_101(numComp, zero);
        std::vector<Scalar> f_110(numComp, zero);
        std::vector<Scalar> f_111(numComp, zero);

        for(p8est_topidx_t t = p8est->first_local_tree; t <= p8est->last_local_tree; t++) // Loop over every tree
        {
            p4est_tree_t * currenttree = p4est_tree_array_index(p8est->trees, t);
            sc_array_t * tquadrants = &currenttree->quadrants;
            p4est_qcoord_t Q = (p4est_locidx_t) tquadrants->elem_count;
#pragma omp parallel for
            for(p4est_qcoord_t e = nodes->global_offset; e < Q+nodes->global_offset; e++) // Loop over every quadrant within the tree
            {
                // Work out what level this element is on 
                p8est_quadrant_t * quad = p8est_quadrant_array_index(tquadrants, e);
                // octantData * quaddata = (octantData *) quad->p.user_data;

                int l = quad->level;

                memcpy(&f_000[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                memcpy(&f_001[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                memcpy(&f_010[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                memcpy(&f_011[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                memcpy(&f_100[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                memcpy(&f_101[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                memcpy(&f_110[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                memcpy(&f_111[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                Scalar* o = out.getSampleDataRW(e, zero);
                for (index_t i=0; i < numComp; ++i) {
                    o[INDEX3(i,0,0,numComp,3)] = (f_100[i]+f_101[i]+f_110[i]+f_111[i]-f_000[i]-f_001[i]-f_010[i]-f_011[i])*C3 / forestData.m_dx[0][l];
                    o[INDEX3(i,1,0,numComp,3)] = (f_010[i]+f_011[i]+f_110[i]+f_111[i]-f_000[i]-f_001[i]-f_100[i]-f_101[i])*C3 / forestData.m_dx[1][l];
                    o[INDEX3(i,2,0,numComp,3)] = (f_001[i]+f_011[i]+f_101[i]+f_111[i]-f_000[i]-f_010[i]-f_100[i]-f_110[i])*C3 / forestData.m_dx[2][l];
                } // end of component loop i
            }
        }
    } 
    else if (out.getFunctionSpace().getTypeCode() == FaceElements) 
    {
        out.requireWrite();

        std::vector<Scalar> f_000(numComp, zero);
        std::vector<Scalar> f_001(numComp, zero);
        std::vector<Scalar> f_010(numComp, zero);
        std::vector<Scalar> f_011(numComp, zero);
        std::vector<Scalar> f_100(numComp, zero);
        std::vector<Scalar> f_101(numComp, zero);
        std::vector<Scalar> f_110(numComp, zero);
        std::vector<Scalar> f_111(numComp, zero);


        for(p8est_topidx_t t = p8est->first_local_tree; t <= p8est->last_local_tree; t++) 
        {
            p8est_tree_t * currenttree = p8est_tree_array_index(p8est->trees, t);
            sc_array_t * tquadrants = &currenttree->quadrants;
            p4est_qcoord_t Q = (p4est_locidx_t) tquadrants->elem_count;
#pragma omp parallel for
            for(p4est_qcoord_t e = nodes->global_offset; e < Q+nodes->global_offset; e++)
            {
                // Work out what level this element is on 
                p8est_quadrant_t * quad = p8est_quadrant_array_index(tquadrants, e);
                octantData * quaddata = (octantData *) quad->p.user_data;

                int l = quad->level;

                if(quaddata->m_faceOffset[0]) 
                {
                    memcpy(&f_000[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_001[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_010[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_011[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_100[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_101[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_110[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_111[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));

                    Scalar* o = out.getSampleDataRW(e, zero);

                    for (index_t i=0; i < numComp; ++i) {
                        const Scalar V0=((f_010[i]-f_000[i])*C6 + (f_011[i]-f_001[i])*C2) / forestData.m_dx[1][l];
                        const Scalar V1=((f_010[i]-f_000[i])*C2 + (f_011[i]-f_001[i])*C6) / forestData.m_dx[1][l];
                        const Scalar V2=((f_001[i]-f_000[i])*C6 + (f_010[i]-f_011[i])*C2) / forestData.m_dx[2][l];
                        const Scalar V3=((f_001[i]-f_000[i])*C2 + (f_011[i]-f_010[i])*C6) / forestData.m_dx[2][l];
                        o[INDEX3(i,0,0,numComp,3)] = ((f_100[i]-f_000[i])*C5 + (f_111[i]-f_011[i])*C0 + (f_101[i]+f_110[i]-f_001[i]-f_010[i])*C1) / forestData.m_dx[0][l];
                        o[INDEX3(i,1,0,numComp,3)] = V0;
                        o[INDEX3(i,2,0,numComp,3)] = V2;
                        o[INDEX3(i,0,1,numComp,3)] = ((f_110[i]-f_010[i])*C5 + (f_101[i]-f_001[i])*C0 + (f_100[i]+f_111[i]-f_000[i]-f_011[i])*C1) / forestData.m_dx[0][l];
                        o[INDEX3(i,1,1,numComp,3)] = V0;
                        o[INDEX3(i,2,1,numComp,3)] = V3;
                        o[INDEX3(i,0,2,numComp,3)] = ((f_101[i]-f_001[i])*C5 + (f_110[i]-f_010[i])*C0 + (f_100[i]+f_111[i]-f_000[i]-f_011[i])*C1) / forestData.m_dx[0][l];
                        o[INDEX3(i,1,2,numComp,3)] = V1;
                        o[INDEX3(i,2,2,numComp,3)] = V2;
                        o[INDEX3(i,0,3,numComp,3)] = ((f_111[i]-f_011[i])*C5 + (f_100[i]-f_000[i])*C0 + (f_101[i]+f_110[i]-f_001[i]-f_010[i])*C1) / forestData.m_dx[0][l];
                        o[INDEX3(i,1,3,numComp,3)] = V1;
                        o[INDEX3(i,2,3,numComp,3)] = V3;
                    } 
                }
            
                if(quaddata->m_faceOffset[1]) 
                {
                    memcpy(&f_000[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_001[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_010[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_011[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_100[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_101[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_110[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_111[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    Scalar* o = out.getSampleDataRW(e, zero);

                    for (index_t i=0; i < numComp; ++i) {
                        const Scalar V0=((f_110[i]-f_100[i])*C6 + (f_111[i]-f_101[i])*C2) / forestData.m_dx[1][l];
                        const Scalar V1=((f_110[i]-f_100[i])*C2 + (f_111[i]-f_101[i])*C6) / forestData.m_dx[1][l];
                        const Scalar V2=((f_101[i]-f_100[i])*C6 + (f_111[i]-f_110[i])*C2) / forestData.m_dx[2][l];
                        const Scalar V3=((f_101[i]-f_100[i])*C2 + (f_111[i]-f_110[i])*C6) / forestData.m_dx[2][l];
                        o[INDEX3(i,0,0,numComp,3)] = ((f_100[i]-f_000[i])*C5 + (f_111[i]-f_011[i])*C0 + (f_101[i]+f_110[i]-f_001[i]-f_010[i])*C1) / forestData.m_dx[0][l];
                        o[INDEX3(i,1,0,numComp,3)] = V0;
                        o[INDEX3(i,2,0,numComp,3)] = V2;
                        o[INDEX3(i,0,1,numComp,3)] = ((f_110[i]-f_010[i])*C5 + (f_101[i]-f_001[i])*C0 + (f_100[i]+f_111[i]-f_000[i]-f_011[i])*C1) / forestData.m_dx[0][l];
                        o[INDEX3(i,1,1,numComp,3)] = V0;
                        o[INDEX3(i,2,1,numComp,3)] = V3;
                        o[INDEX3(i,0,2,numComp,3)] = ((f_101[i]-f_001[i])*C5 + (f_110[i]-f_010[i])*C0 + (f_100[i]+f_111[i]-f_000[i]-f_011[i])*C1) / forestData.m_dx[0][l];
                        o[INDEX3(i,1,2,numComp,3)] = V1;
                        o[INDEX3(i,2,2,numComp,3)] = V2;
                        o[INDEX3(i,0,3,numComp,3)] = ((f_111[i]-f_011[i])*C5 + (f_100[i]-f_000[i])*C0 + (f_101[i]+f_110[i]-f_001[i]-f_010[i])*C1) / forestData.m_dx[0][l];
                        o[INDEX3(i,1,3,numComp,3)] = V1;
                        o[INDEX3(i,2,3,numComp,3)] = V3;
                    } // end of component loop i
                }

                if(quaddata->m_faceOffset[2]) 
                {
                    memcpy(&f_000[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_001[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_010[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_011[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_100[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_101[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_110[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_111[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    Scalar* o = out.getSampleDataRW(e, zero);

                    for (index_t i=0; i < numComp; ++i) {
                        const Scalar V0=((f_100[i]-f_000[i])*C6 + (f_101[i]-f_001[i])*C2) / forestData.m_dx[0][l];
                        const Scalar V1=((f_001[i]-f_000[i])*C6 + (f_101[i]-f_100[i])*C2) / forestData.m_dx[2][l];
                        const Scalar V2=((f_001[i]-f_000[i])*C2 + (f_101[i]-f_100[i])*C6) / forestData.m_dx[2][l];
                        o[INDEX3(i,0,0,numComp,3)] = V0;
                        o[INDEX3(i,1,0,numComp,3)] = ((f_010[i]-f_000[i])*C5 + (f_111[i]-f_101[i])*C0 + (f_011[i]+f_110[i]-f_001[i]-f_100[i])*C1) / forestData.m_dx[1][l];
                        o[INDEX3(i,2,0,numComp,3)] = V1;
                        o[INDEX3(i,0,1,numComp,3)] = V0;
                        o[INDEX3(i,1,1,numComp,3)] = ((f_110[i]-f_100[i])*C5 + (f_011[i]-f_001[i])*C0 + (f_010[i]+f_111[i]-f_000[i]-f_101[i])*C1) / forestData.m_dx[1][l];
                        o[INDEX3(i,2,1,numComp,3)] = V2;
                        o[INDEX3(i,0,2,numComp,3)] = V0;
                        o[INDEX3(i,1,2,numComp,3)] = ((f_011[i]-f_001[i])*C5 + (f_110[i]-f_100[i])*C0 + (f_010[i]+f_111[i]-f_000[i]-f_101[i])*C1) / forestData.m_dx[1][l];
                        o[INDEX3(i,2,2,numComp,3)] = V1;
                        o[INDEX3(i,0,3,numComp,3)] = V0;
                        o[INDEX3(i,1,3,numComp,3)] = ((f_111[i]-f_101[i])*C5 + (f_010[i]-f_000[i])*C0 + (f_011[i]+f_110[i]-f_001[i]-f_100[i])*C1) / forestData.m_dx[1][l];
                        o[INDEX3(i,2,3,numComp,3)] = V2;
                    } // end of component loop i
                } // end of face 2

                if(quaddata->m_faceOffset[3]) 
                {
                    memcpy(&f_000[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_001[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_010[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_011[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_100[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_101[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_110[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_111[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    Scalar* o = out.getSampleDataRW(e, zero);

                    for (index_t i=0; i < numComp; ++i) {
                        const Scalar V0=((f_110[i]-f_010[i])*C6 + (f_111[i]-f_011[i])*C2) / forestData.m_dx[0][l];
                        const Scalar V1=((f_110[i]-f_010[i])*C2 + (f_111[i]-f_011[i])*C6) / forestData.m_dx[0][l];
                        const Scalar V2=((f_011[i]-f_010[i])*C6 + (f_111[i]-f_110[i])*C2) / forestData.m_dx[2][l];
                        const Scalar V3=((f_011[i]-f_010[i])*C2 + (f_111[i]-f_110[i])*C6) / forestData.m_dx[2][l];
                        o[INDEX3(i,0,0,numComp,3)] = V0;
                        o[INDEX3(i,1,0,numComp,3)] = ((f_010[i]-f_000[i])*C5 + (f_111[i]-f_101[i])*C0 + (f_011[i]+f_110[i]-f_001[i]-f_100[i])*C1) / forestData.m_dx[1][l];
                        o[INDEX3(i,2,0,numComp,3)] = V2;
                        o[INDEX3(i,0,1,numComp,3)] = V0;
                        o[INDEX3(i,1,1,numComp,3)] = ((f_110[i]-f_100[i])*C5 + (f_011[i]-f_001[i])*C0 + (f_010[i]+f_111[i]-f_000[i]-f_101[i])*C1) / forestData.m_dx[1][l];
                        o[INDEX3(i,2,1,numComp,3)] = V3;
                        o[INDEX3(i,0,2,numComp,3)] = V1;
                        o[INDEX3(i,1,2,numComp,3)] = ((f_011[i]-f_001[i])*C5 + (f_110[i]-f_100[i])*C0 + (f_010[i]+f_111[i]-f_000[i]-f_101[i])*C1) / forestData.m_dx[1][l];
                        o[INDEX3(i,2,2,numComp,3)] = V2;
                        o[INDEX3(i,0,3,numComp,3)] = V1;
                        o[INDEX3(i,1,3,numComp,3)] = ((f_111[i]-f_101[i])*C5 + (f_010[i]-f_000[i])*C0 + (f_011[i]+f_110[i]-f_001[i]-f_100[i])*C1) / forestData.m_dx[1][l];
                        o[INDEX3(i,2,3,numComp,3)] = V3;
                    } // end of component loop i
                } // end of face 3

                if(quaddata->m_faceOffset[4]) 
                {
                    memcpy(&f_000[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_001[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_010[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_011[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_100[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_101[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_110[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_111[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    Scalar* o = out.getSampleDataRW(e, zero);
                    for (index_t i=0; i < numComp; ++i) {
                        const Scalar V0=((f_100[i]-f_000[i])*C6 + (f_110[i]-f_010[i])*C2) / forestData.m_dx[0][l];
                        const Scalar V1=((f_100[i]-f_000[i])*C2 + (f_110[i]-f_010[i])*C6) / forestData.m_dx[0][l];
                        const Scalar V2=((f_010[i]-f_000[i])*C6 + (f_110[i]-f_100[i])*C2) / forestData.m_dx[1][l];
                        const Scalar V3=((f_010[i]-f_000[i])*C2 + (f_110[i]-f_100[i])*C6) / forestData.m_dx[1][l];
                        o[INDEX3(i,0,0,numComp,3)] = V0;
                        o[INDEX3(i,1,0,numComp,3)] = V2;
                        o[INDEX3(i,2,0,numComp,3)] = ((f_001[i]-f_000[i])*C5 + (f_111[i]-f_110[i])*C0 + (f_011[i]+f_101[i]-f_010[i]-f_100[i])*C1) / forestData.m_dx[2][l];
                        o[INDEX3(i,0,1,numComp,3)] = V0;
                        o[INDEX3(i,1,1,numComp,3)] = V3;
                        o[INDEX3(i,2,1,numComp,3)] = ((f_101[i]-f_100[i])*C5 + (f_011[i]-f_010[i])*C0 + (f_001[i]+f_111[i]-f_000[i]-f_110[i])*C1) / forestData.m_dx[2][l];
                        o[INDEX3(i,0,2,numComp,3)] = V1;
                        o[INDEX3(i,1,2,numComp,3)] = V2;
                        o[INDEX3(i,2,2,numComp,3)] = ((f_011[i]-f_010[i])*C5 + (f_101[i]-f_100[i])*C0 + (f_001[i]+f_111[i]-f_000[i]-f_110[i])*C1) / forestData.m_dx[2][l];
                        o[INDEX3(i,0,3,numComp,3)] = V1;
                        o[INDEX3(i,1,3,numComp,3)] = V3;
                        o[INDEX3(i,2,3,numComp,3)] = ((f_111[i]-f_110[i])*C5 + (f_001[i]-f_000[i])*C0 + (f_011[i]+f_101[i]-f_010[i]-f_100[i])*C1) / forestData.m_dx[2][l];
                    } // end of component loop i
                } // end of face 4

                if(quaddata->m_faceOffset[5]) 
                {
                    memcpy(&f_000[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_001[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_010[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_011[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_100[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_101[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_110[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_111[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    Scalar* o = out.getSampleDataRW(e, zero);
                    for (index_t i=0; i < numComp; ++i) {
                        const Scalar V0=((f_101[i]-f_001[i])*C6 + (f_111[i]-f_011[i])*C2) / forestData.m_dx[0][l];
                        const Scalar V1=((f_101[i]-f_001[i])*C2 + (f_111[i]-f_011[i])*C6) / forestData.m_dx[0][l];
                        const Scalar V2=((f_011[i]-f_001[i])*C6 + (f_111[i]-f_101[i])*C2) / forestData.m_dx[1][l];
                        const Scalar V3=((f_011[i]-f_001[i])*C2 + (f_111[i]-f_101[i])*C6) / forestData.m_dx[1][l];
                        o[INDEX3(i,0,0,numComp,3)] = V0;
                        o[INDEX3(i,1,0,numComp,3)] = V2;
                        o[INDEX3(i,2,0,numComp,3)] = ((f_001[i]-f_000[i])*C5 + (f_111[i]-f_110[i])*C0 + (f_011[i]+f_101[i]-f_010[i]-f_100[i])*C1) / forestData.m_dx[2][l];
                        o[INDEX3(i,0,1,numComp,3)] = V0;
                        o[INDEX3(i,1,1,numComp,3)] = V3;
                        o[INDEX3(i,2,1,numComp,3)] = ((f_011[i]-f_010[i])*C0 + (f_101[i]-f_100[i])*C5 + (f_001[i]+f_111[i]-f_000[i]-f_110[i])*C1) / forestData.m_dx[2][l];
                        o[INDEX3(i,0,2,numComp,3)] = V1;
                        o[INDEX3(i,1,2,numComp,3)] = V2;
                        o[INDEX3(i,2,2,numComp,3)] = ((f_011[i]-f_010[i])*C5 + (f_101[i]-f_100[i])*C0 + (f_001[i]+f_111[i]-f_000[i]-f_110[i])*C1) / forestData.m_dx[2][l];
                        o[INDEX3(i,0,3,numComp,3)] = V1;
                        o[INDEX3(i,1,3,numComp,3)] = V3;
                        o[INDEX3(i,2,3,numComp,3)] = ((f_001[i]-f_000[i])*C0 + (f_111[i]-f_110[i])*C5 + (f_011[i]+f_101[i]-f_010[i]-f_100[i])*C1) / forestData.m_dx[2][l];
                    } // end of component loop i
                } // end of face 5
            }
        }
    } else if (out.getFunctionSpace().getTypeCode() == ReducedFaceElements) {

        out.requireWrite();

        for(p8est_topidx_t t = p8est->first_local_tree; t <= p8est->last_local_tree; t++) 
        {
            p8est_tree_t * currenttree = p8est_tree_array_index(p8est->trees, t);
            sc_array_t * tquadrants = &currenttree->quadrants;
            p4est_qcoord_t Q = (p4est_locidx_t) tquadrants->elem_count;
#pragma omp parallel for
            for(p4est_qcoord_t e = nodes->global_offset; e < Q+nodes->global_offset; e++)
            {
                // Work out what level this element is on 
                p8est_quadrant_t * quad = p8est_quadrant_array_index(tquadrants, e);
                octantData * quaddata = (octantData *) quad->p.user_data;

                int l = quad->level;

                std::vector<Scalar> f_000(numComp, zero);
                std::vector<Scalar> f_001(numComp, zero);
                std::vector<Scalar> f_010(numComp, zero);
                std::vector<Scalar> f_011(numComp, zero);
                std::vector<Scalar> f_100(numComp, zero);
                std::vector<Scalar> f_101(numComp, zero);
                std::vector<Scalar> f_110(numComp, zero);
                std::vector<Scalar> f_111(numComp, zero);



                if(quaddata->m_faceOffset[0]) 
                {
                    memcpy(&f_000[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_001[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_010[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_011[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_100[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_101[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_110[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_111[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    Scalar* o = out.getSampleDataRW(e, zero);
                    for (index_t i=0; i < numComp; ++i) {
                        o[INDEX3(i,0,0,numComp,3)] = (f_100[i]+f_101[i]+f_110[i]+f_111[i]-f_000[i]-f_001[i]-f_010[i]-f_011[i])*C3 / forestData.m_dx[0][l];
                        o[INDEX3(i,1,0,numComp,3)] = (f_010[i]+f_011[i]-f_000[i]-f_001[i])*C4 / forestData.m_dx[1][l];
                        o[INDEX3(i,2,0,numComp,3)] = (f_001[i]+f_011[i]-f_000[i]-f_010[i])*C4 / forestData.m_dx[2][l];
                    } // end of component loop i
                } // end of face 0


                if(quaddata->m_faceOffset[1]) 
                {
                    memcpy(&f_000[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_001[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_010[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_011[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_100[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_101[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_110[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_111[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    Scalar* o = out.getSampleDataRW(e, zero);
                    for (index_t i=0; i < numComp; ++i) {
                        o[INDEX3(i,0,0,numComp,3)] = (f_100[i]+f_101[i]+f_110[i]+f_111[i]-f_000[i]-f_001[i]-f_010[i]-f_011[i])*C3 / forestData.m_dx[0][l];
                        o[INDEX3(i,1,0,numComp,3)] = (f_110[i]+f_111[i]-f_100[i]-f_101[i])*C4 / forestData.m_dx[1][l];
                        o[INDEX3(i,2,0,numComp,3)] = (f_101[i]+f_111[i]-f_100[i]-f_110[i])*C4 / forestData.m_dx[2][l];
                    } // end of component loop i
                } // end of face 1

                if(quaddata->m_faceOffset[2]) 
                {
                    memcpy(&f_000[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_001[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_010[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_011[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_100[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_101[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_110[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_111[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    Scalar* o = out.getSampleDataRW(e, zero);
                    for (index_t i=0; i < numComp; ++i) {
                        o[INDEX3(i,0,0,numComp,3)] = (f_100[i]+f_101[i]-f_000[i]-f_001[i])*C4 / forestData.m_dx[0][l];
                        o[INDEX3(i,1,0,numComp,3)] = (f_010[i]+f_011[i]+f_110[i]+f_111[i]-f_000[i]-f_001[i]-f_100[i]-f_101[i])*C3 / forestData.m_dx[1][l];
                        o[INDEX3(i,2,0,numComp,3)] = (f_001[i]+f_101[i]-f_000[i]-f_100[i])*C4 / forestData.m_dx[2][l];
                    } // end of component loop i
                } // end of face 2


                if(quaddata->m_faceOffset[3]) 
                {
                    memcpy(&f_000[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_001[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_010[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_011[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_100[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_101[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_110[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_111[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    Scalar* o = out.getSampleDataRW(e, zero);
                    for (index_t i=0; i < numComp; ++i) {
                            o[INDEX3(i,0,0,numComp,3)] = (f_110[i]+f_111[i]-f_010[i]-f_011[i])*C4 / forestData.m_dx[0][l];
                            o[INDEX3(i,1,0,numComp,3)] = (f_010[i]+f_011[i]+f_110[i]+f_111[i]-f_000[i]-f_001[i]-f_100[i]-f_101[i])*C3 / forestData.m_dx[1][l];
                            o[INDEX3(i,2,0,numComp,3)] = (f_011[i]+f_111[i]-f_010[i]-f_110[i])*C4 / forestData.m_dx[2][l];
                    } // end of component loop i
                } // end of face 3
    
                if(quaddata->m_faceOffset[4]) 
                {
                    memcpy(&f_000[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_001[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_010[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_011[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_100[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_101[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_110[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_111[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    Scalar* o = out.getSampleDataRW(e, zero);
                    for (index_t i=0; i < numComp; ++i) {
                            o[INDEX3(i,0,0,numComp,3)] = (f_100[i]+f_110[i]-f_000[i]-f_010[i])*C4 / forestData.m_dx[0][l];
                            o[INDEX3(i,1,0,numComp,3)] = (f_010[i]+f_110[i]-f_000[i]-f_100[i])*C4 / forestData.m_dx[1][l];
                            o[INDEX3(i,2,0,numComp,3)] = (f_001[i]+f_011[i]+f_101[i]+f_111[i]-f_000[i]-f_010[i]-f_100[i]-f_110[i])*C4 / forestData.m_dx[2][l];
                    } // end of component loop i
                } // end of face 4

                if(quaddata->m_faceOffset[5]) 
                {
                    memcpy(&f_000[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_001[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_010[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_011[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_100[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_101[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_110[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    memcpy(&f_111[0], in.getSampleDataRO(e, zero), numComp*sizeof(Scalar));
                    Scalar* o = out.getSampleDataRW(e, zero);
                    for (index_t i=0; i < numComp; ++i) {
                            o[INDEX3(i,0,0,numComp,3)] = (f_101[i]+f_111[i]-f_001[i]-f_011[i])*C4 / forestData.m_dx[0][l];
                            o[INDEX3(i,1,0,numComp,3)] = (f_011[i]+f_111[i]-f_001[i]-f_101[i])*C4 / forestData.m_dx[1][l];
                            o[INDEX3(i,2,0,numComp,3)] = (f_001[i]+f_011[i]+f_101[i]+f_111[i]-f_000[i]-f_010[i]-f_100[i]-f_110[i])*C3 / forestData.m_dx[2][l];
                    } // end of component loop i
                } // end of face 5
            }
        }
    }
}

//protected
void Brick::nodesToDOF(escript::Data& out, const escript::Data& in) const
{
    //TODO
    throw OxleyException("nodesToDOF");

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

////////////////////////////// inline methods ////////////////////////////////
inline dim_t Brick::getDofOfNode(dim_t node) const
{
    //TODO
    throw OxleyException("getDofOfNode");
    return -1;
    // return m_dofMap[node];
}

dim_t Brick::findNode(const double *coords) const
{
    //TODO
    throw OxleyException("findNode");
    return -1;

    
}

// This is a wrapper that converts the p4est node information into an IndexVector
IndexVector Brick::getNodeDistribution() const
{
    return m_nodeDistribution;
}

// This is a wrapper that converts the p4est node information into an IndexVector
void Brick::updateNodeDistribution() 
{
    // TODO
    // m_nodeDistribution.clear();
    // m_nodeDistribution.assign(MAXP4ESTNODES,0);

    // int counter =0;
    // for(p4est_topidx_t treeid = p4est->first_local_tree; treeid <= p4est->last_local_tree; ++treeid) 
    // {
    //     p4est_tree_t * tree = p4est_tree_array_index(p4est->trees, treeid);
    //     sc_array_t * tquadrants = &tree->quadrants;
    //     p4est_locidx_t Q = (p4est_locidx_t) tquadrants->elem_count;
    //     for(int q = 0; q < Q; ++q) 
    //     { 
    //         p4est_quadrant_t * quad = p4est_quadrant_array_index(tquadrants, q);
    //         p4est_qcoord_t length = P4EST_QUADRANT_LEN(quad->level);
    //         for(int n = 0; n < 4; n++)
    //         {
    //             double lx = length * ((int) (n % 2) == 1);
    //             double ly = length * ((int) (n / 2) == 1);
    //             double xy[3];
    //             p4est_qcoord_to_vertex(p4est->connectivity, treeid, quad->x+lx, quad->y+ly, xy);
    //             m_nodeDistribution[counter++]=NodeIDs.find(std::make_pair(xy[0],xy[1]))->second;
    //         }
    //     }
    // }
    // m_nodeDistribution.shrink_to_fit();
}

// adds the dirac points and tags 
void Brick::addPoints(const std::vector<double>& coords, const std::vector<int>& tags)
{
    
}

std::vector<IndexVector> Brick::getConnections(bool includeShared) const
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
    for(p8est_topidx_t treeid = p8est->first_local_tree; treeid <= p8est->last_local_tree; ++treeid) 
    {
        p8est_tree_t * tree = p8est_tree_array_index(p8est->trees, treeid);
        sc_array_t * tquadrants = &tree->quadrants;
        long Q = (p4est_qcoord_t) tquadrants->elem_count;
// #pragma omp parallel for
        for(long q = 0; q < Q; ++q) // Loop over all quadrants
        { 
            p8est_quadrant_t * quad = p8est_quadrant_array_index(tquadrants, q);
            long length = P8EST_QUADRANT_LEN(quad->level);
            for(int n = 0; n < 8; n++)
            {
                double xyz[3];
                long lx[8] = {0,length,0,length,0,length,0,length};
                long ly[8] = {0,0,length,length,0,0,length,length};
                long lz[8] = {0,0,0,0,length,length,length,length};
                long lni[4] = {-1};
                for(int i = 0; i < 8; i++)
                {
                    p8est_qcoord_to_vertex(p8est->connectivity, treeid, quad->x+lx[i], quad->y+ly[i], quad->z+lz[i], xyz);
                    lni[i] = NodeIDs.find(std::make_tuple(xyz[0],xyz[1],xyz[2]))->second;
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

static inline void
brick_linear_to_xyz (p4est_topidx_t ti, const int logx[P8EST_DIM],
                     const int rankx[P8EST_DIM], p4est_topidx_t tx[P8EST_DIM])
{
    int i, j, k;
    int lastlog = 0;

    for (i = 0; i < P8EST_DIM; i++) {
        tx[i] = 0;
    }

    for (i = 0; i < P8EST_DIM - 1; i++) {
        p4est_topidx_t tempx[3] = { 0, 0, 0 };
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

static inline p4est_topidx_t
brick_xyz_to_linear (const p4est_topidx_t tx[P8EST_DIM],
                     const int logx[P8EST_DIM], const int rankx[P8EST_DIM])
{
    int i, j, k;
    int lastlog = logx[rankx[P8EST_DIM - 2]];
    p4est_topidx_t ti = tx[rankx[P8EST_DIM - 1]] >> lastlog;

    for (i = P8EST_DIM - 2; i >= 0; i--) {
        p4est_topidx_t tempx[3] = { 0, 0, 0 };
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
    const p8est_topidx_t m = (p4est_topidx_t) n0;
    const p8est_topidx_t n = (p4est_topidx_t) n1;
    const p8est_topidx_t p = (p4est_topidx_t) n2;
    const p8est_topidx_t mc = periodic_a ? m : (m - 1);
    const p8est_topidx_t nc = periodic_b ? n : (n - 1);
    const p8est_topidx_t pc = periodic_c ? p : (p - 1);
    const p8est_topidx_t num_trees = m * n * p;
    const p8est_topidx_t num_corners = mc * nc * pc;
    const p8est_topidx_t num_ctt = P8EST_CHILDREN * num_corners;
    const p8est_topidx_t num_edges = m * nc * pc + mc * n * pc + mc * nc * p;
    const p8est_topidx_t num_ett = 4 * num_edges;
    const p8est_topidx_t num_vertices = (m + 1) * (n + 1) * (p + 1);
    const int periodic[P8EST_DIM] = { periodic_a, periodic_b, periodic_c };

    const p8est_topidx_t max[P8EST_DIM] = { m - 1, n - 1, p - 1 };
    double *vertices;
    p8est_topidx_t *tree_to_vertex;
    p8est_topidx_t *tree_to_tree;
    int8_t *tree_to_face;
    p8est_topidx_t *tree_to_corner;
    p8est_topidx_t *ctt_offset;
    p8est_topidx_t *corner_to_tree;
    int8_t *corner_to_corner;
    p8est_topidx_t  n_iter;
    int logx[P8EST_DIM];
    int rankx[P8EST_DIM];
    int i, j, l;
    p8est_topidx_t  ti, tj, tk;
    p8est_topidx_t  tx, ty;
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
    p8est_topidx_t *tree_to_edge;
    p8est_topidx_t *ett_offset;
    p8est_topidx_t *edge_to_tree;
    int8_t *edge_to_edge;
    p8est_topidx_t *tree_to_edge2;
    int dir1, dir2;

    ESYS_ASSERT(m > 0 && n > 0 && p > 0, "n0, n1 and n2 must be greater than zero.");


    conn = p8est_connectivity_new (num_vertices, num_trees, 
                                 num_edges, num_ett,
                                 num_corners, num_ctt);

    vertices = conn->vertices;
    tree_to_vertex = conn->tree_to_vertex;
    tree_to_tree = conn->tree_to_tree;
    tree_to_face = conn->tree_to_face;
    tree_to_edge = conn->tree_to_edge;
    ett_offset = conn->ett_offset;
    edge_to_tree = conn->edge_to_tree;
    edge_to_edge = conn->edge_to_edge;
    tree_to_corner = conn->tree_to_corner;
    ctt_offset = conn->ctt_offset;
    corner_to_tree = conn->corner_to_tree;
    corner_to_corner = conn->corner_to_corner;

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


    linear_to_tree = P4EST_ALLOC(p4est_topidx_t, n_iter);
    tree_to_corner2 = P4EST_ALLOC(p4est_topidx_t, num_trees);
    tree_to_edge2 = P4EST_ALLOC(p4est_topidx_t, 3 * num_trees);

    tj = 0;
    tk = 0;
    tl = 0;

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
                        vertices[vicount++] = (double) (tx + (i & 1));
                        vertices[vicount++] = (double) (ty + ((i >> 1) & 1));
                        vertices[vicount++] = (double) (tz + (i >> 2));
                    }
                }
            }
        }
    }

    P4EST_ASSERT(vcount == num_vertices);
    P4EST_FREE(linear_to_tree);
    P4EST_FREE(tree_to_corner2);
    P4EST_FREE(tree_to_edge2);

#ifdef OXLEY_ENABLE_DEBUG
    P4EST_ASSERT (p8est_connectivity_is_valid (conn)); //This is very time consuming
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
//     const dim_t numComp = arg.getDataPointSize();
//     const index_t left = (m_offset[0]==0 ? 0 : 1);
//     const index_t bottom = (m_offset[1]==0 ? 0 : 1);
//     const index_t front = (m_offset[2]==0 ? 0 : 1);
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
//         const real_t w_0 = m_dx[0]*m_dx[1]*m_dx[2]/8.;
// #pragma omp parallel
//         {
//             vector<Scalar> int_local(numComp, zero);
// #pragma omp for nowait
//             for (index_t k2 = front; k2 < front+m_ownNE[2]; ++k2) {
//                 for (index_t k1 = bottom; k1 < bottom+m_ownNE[1]; ++k1) {
//                     for (index_t k0 = left; k0 < left+m_ownNE[0]; ++k0) {
//                         const Scalar* f = arg.getSampleDataRO(INDEX3(k0, k1, k2, m_NE[0], m_NE[1]), zero);
//                         for (index_t i = 0; i < numComp; ++i) {
//                             const Scalar f_0 = f[INDEX2(i,0,numComp)];
//                             const Scalar f_1 = f[INDEX2(i,1,numComp)];
//                             const Scalar f_2 = f[INDEX2(i,2,numComp)];
//                             const Scalar f_3 = f[INDEX2(i,3,numComp)];
//                             const Scalar f_4 = f[INDEX2(i,4,numComp)];
//                             const Scalar f_5 = f[INDEX2(i,5,numComp)];
//                             const Scalar f_6 = f[INDEX2(i,6,numComp)];
//                             const Scalar f_7 = f[INDEX2(i,7,numComp)];
//                             int_local[i]+=(f_0+f_1+f_2+f_3+f_4+f_5+f_6+f_7)*w_0;
//                         }  // end of component loop i
//                     } // end of k0 loop
//                 } // end of k1 loop
//             } // end of k2 loop

// #pragma omp critical
//             for (index_t i = 0; i < numComp; i++)
//                 integrals[i] += int_local[i];
//         } // end of parallel section

//     } else if (fs==ReducedElements || (fs==Elements && !arg.actsExpanded())) {
//         const real_t w_0 = m_dx[0]*m_dx[1]*m_dx[2];
// #pragma omp parallel
//         {
//             vector<Scalar> int_local(numComp, zero);
// #pragma omp for nowait
//             for (index_t k2 = front; k2 < front+m_ownNE[2]; ++k2) {
//                 for (index_t k1 = bottom; k1 < bottom+m_ownNE[1]; ++k1) {
//                     for (index_t k0 = left; k0 < left+m_ownNE[0]; ++k0) {
//                         const Scalar* f = arg.getSampleDataRO(INDEX3(k0, k1, k2, m_NE[0], m_NE[1]), zero);
//                         for (index_t i = 0; i < numComp; ++i) {
//                             int_local[i] += f[i]*w_0;
//                         }  // end of component loop i
//                     } // end of k0 loop
//                 } // end of k1 loop
//             } // end of k2 loop

// #pragma omp critical
//             for (index_t i = 0; i < numComp; i++)
//                 integrals[i] += int_local[i];
//         } // end of parallel section

//     } else if (fs == FaceElements && arg.actsExpanded()) {
//         const real_t w_0 = m_dx[1]*m_dx[2]/4.;
//         const real_t w_1 = m_dx[0]*m_dx[2]/4.;
//         const real_t w_2 = m_dx[0]*m_dx[1]/4.;
// #pragma omp parallel
//         {
//             vector<Scalar> int_local(numComp, zero);
//             if (m_faceOffset[0] > -1) {
// #pragma omp for nowait
//                 for (index_t k2 = front; k2 < front+m_ownNE[2]; ++k2) {
//                     for (index_t k1 = bottom; k1 < bottom+m_ownNE[1]; ++k1) {
//                         const Scalar* f = arg.getSampleDataRO(m_faceOffset[0]+INDEX2(k1,k2,m_NE[1]), zero);
//                         for (index_t i = 0; i < numComp; ++i) {
//                             const Scalar f_0 = f[INDEX2(i,0,numComp)];
//                             const Scalar f_1 = f[INDEX2(i,1,numComp)];
//                             const Scalar f_2 = f[INDEX2(i,2,numComp)];
//                             const Scalar f_3 = f[INDEX2(i,3,numComp)];
//                             int_local[i] += (f_0+f_1+f_2+f_3)*w_0;
//                         }  // end of component loop i
//                     } // end of k1 loop
//                 } // end of k2 loop
//             }

//             if (m_faceOffset[1] > -1) {
// #pragma omp for nowait
//                 for (index_t k2 = front; k2 < front+m_ownNE[2]; ++k2) {
//                     for (index_t k1 = bottom; k1 < bottom+m_ownNE[1]; ++k1) {
//                         const Scalar* f = arg.getSampleDataRO(m_faceOffset[1]+INDEX2(k1,k2,m_NE[1]), zero);
//                         for (index_t i = 0; i < numComp; ++i) {
//                             const Scalar f_0 = f[INDEX2(i,0,numComp)];
//                             const Scalar f_1 = f[INDEX2(i,1,numComp)];
//                             const Scalar f_2 = f[INDEX2(i,2,numComp)];
//                             const Scalar f_3 = f[INDEX2(i,3,numComp)];
//                             int_local[i]+=(f_0+f_1+f_2+f_3)*w_0;
//                         }  // end of component loop i
//                     } // end of k1 loop
//                 } // end of k2 loop
//             }

//             if (m_faceOffset[2] > -1) {
// #pragma omp for nowait
//                 for (index_t k2 = front; k2 < front+m_ownNE[2]; ++k2) {
//                     for (index_t k0 = left; k0 < left+m_ownNE[0]; ++k0) {
//                         const Scalar* f = arg.getSampleDataRO(m_faceOffset[2]+INDEX2(k0,k2,m_NE[0]), zero);
//                         for (index_t i = 0; i < numComp; ++i) {
//                             const Scalar f_0 = f[INDEX2(i,0,numComp)];
//                             const Scalar f_1 = f[INDEX2(i,1,numComp)];
//                             const Scalar f_2 = f[INDEX2(i,2,numComp)];
//                             const Scalar f_3 = f[INDEX2(i,3,numComp)];
//                             int_local[i]+=(f_0+f_1+f_2+f_3)*w_1;
//                         }  // end of component loop i
//                     } // end of k1 loop
//                 } // end of k2 loop
//             }

//             if (m_faceOffset[3] > -1) {
// #pragma omp for nowait
//                 for (index_t k2 = front; k2 < front+m_ownNE[2]; ++k2) {
//                     for (index_t k0 = left; k0 < left+m_ownNE[0]; ++k0) {
//                         const Scalar* f = arg.getSampleDataRO(m_faceOffset[3]+INDEX2(k0,k2,m_NE[0]), zero);
//                         for (index_t i = 0; i < numComp; ++i) {
//                             const Scalar f_0 = f[INDEX2(i,0,numComp)];
//                             const Scalar f_1 = f[INDEX2(i,1,numComp)];
//                             const Scalar f_2 = f[INDEX2(i,2,numComp)];
//                             const Scalar f_3 = f[INDEX2(i,3,numComp)];
//                             int_local[i] += (f_0+f_1+f_2+f_3)*w_1;
//                         }  // end of component loop i
//                     } // end of k1 loop
//                 } // end of k2 loop
//             }

//             if (m_faceOffset[4] > -1) {
// #pragma omp for nowait
//                 for (index_t k1 = bottom; k1 < bottom+m_ownNE[1]; ++k1) {
//                     for (index_t k0 = left; k0 < left+m_ownNE[0]; ++k0) {
//                         const Scalar* f = arg.getSampleDataRO(m_faceOffset[4]+INDEX2(k0,k1,m_NE[0]), zero);
//                         for (index_t i = 0; i < numComp; ++i) {
//                             const Scalar f_0 = f[INDEX2(i,0,numComp)];
//                             const Scalar f_1 = f[INDEX2(i,1,numComp)];
//                             const Scalar f_2 = f[INDEX2(i,2,numComp)];
//                             const Scalar f_3 = f[INDEX2(i,3,numComp)];
//                             int_local[i] += (f_0+f_1+f_2+f_3)*w_2;
//                         }  // end of component loop i
//                     } // end of k1 loop
//                 } // end of k2 loop
//             }

//             if (m_faceOffset[5] > -1) {
// #pragma omp for nowait
//                 for (index_t k1 = bottom; k1 < bottom+m_ownNE[1]; ++k1) {
//                     for (index_t k0 = left; k0 < left+m_ownNE[0]; ++k0) {
//                         const Scalar* f = arg.getSampleDataRO(m_faceOffset[5]+INDEX2(k0,k1,m_NE[0]), zero);
//                         for (index_t i = 0; i < numComp; ++i) {
//                             const Scalar f_0 = f[INDEX2(i,0,numComp)];
//                             const Scalar f_1 = f[INDEX2(i,1,numComp)];
//                             const Scalar f_2 = f[INDEX2(i,2,numComp)];
//                             const Scalar f_3 = f[INDEX2(i,3,numComp)];
//                             int_local[i]+=(f_0+f_1+f_2+f_3)*w_2;
//                         }  // end of component loop i
//                     } // end of k1 loop
//                 } // end of k2 loop
//             }

// #pragma omp critical
//             for (index_t i = 0; i < numComp; i++)
//                 integrals[i] += int_local[i];
//         } // end of parallel section

//     } else if (fs==ReducedFaceElements || (fs==FaceElements && !arg.actsExpanded())) {
//         const real_t w_0 = m_dx[1]*m_dx[2];
//         const real_t w_1 = m_dx[0]*m_dx[2];
//         const real_t w_2 = m_dx[0]*m_dx[1];
// #pragma omp parallel
//         {
//             vector<Scalar> int_local(numComp, zero);
//             if (m_faceOffset[0] > -1) {
// #pragma omp for nowait
//                 for (index_t k2 = front; k2 < front+m_ownNE[2]; ++k2) {
//                     for (index_t k1 = bottom; k1 < bottom+m_ownNE[1]; ++k1) {
//                         const Scalar* f = arg.getSampleDataRO(m_faceOffset[0]+INDEX2(k1,k2,m_NE[1]), zero);
//                         for (index_t i = 0; i < numComp; ++i) {
//                             int_local[i] += f[i]*w_0;
//                         }  // end of component loop i
//                     } // end of k1 loop
//                 } // end of k2 loop
//             }

//             if (m_faceOffset[1] > -1) {
// #pragma omp for nowait
//                 for (index_t k2 = front; k2 < front+m_ownNE[2]; ++k2) {
//                     for (index_t k1 = bottom; k1 < bottom+m_ownNE[1]; ++k1) {
//                         const Scalar* f = arg.getSampleDataRO(m_faceOffset[1]+INDEX2(k1,k2,m_NE[1]), zero);
//                         for (index_t i = 0; i < numComp; ++i) {
//                             int_local[i] += f[i]*w_0;
//                         }  // end of component loop i
//                     } // end of k1 loop
//                 } // end of k2 loop
//             }

//             if (m_faceOffset[2] > -1) {
// #pragma omp for nowait
//                 for (index_t k2 = front; k2 < front+m_ownNE[2]; ++k2) {
//                     for (index_t k0 = left; k0 < left+m_ownNE[0]; ++k0) {
//                         const Scalar* f = arg.getSampleDataRO(m_faceOffset[2]+INDEX2(k0,k2,m_NE[0]), zero);
//                         for (index_t i = 0; i < numComp; ++i) {
//                             int_local[i] += f[i]*w_1;
//                         }  // end of component loop i
//                     } // end of k1 loop
//                 } // end of k2 loop
//             }

//             if (m_faceOffset[3] > -1) {
// #pragma omp for nowait
//                 for (index_t k2 = front; k2 < front+m_ownNE[2]; ++k2) {
//                     for (index_t k0 = left; k0 < left+m_ownNE[0]; ++k0) {
//                         const Scalar* f = arg.getSampleDataRO(m_faceOffset[3]+INDEX2(k0,k2,m_NE[0]), zero);
//                         for (index_t i = 0; i < numComp; ++i) {
//                             int_local[i] += f[i]*w_1;
//                         }  // end of component loop i
//                     } // end of k1 loop
//                 } // end of k2 loop
//             }

//             if (m_faceOffset[4] > -1) {
// #pragma omp for nowait
//                 for (index_t k1 = bottom; k1 < bottom+m_ownNE[1]; ++k1) {
//                     for (index_t k0 = left; k0 < left+m_ownNE[0]; ++k0) {
//                         const Scalar* f = arg.getSampleDataRO(m_faceOffset[4]+INDEX2(k0,k1,m_NE[0]), zero);
//                         for (index_t i = 0; i < numComp; ++i) {
//                             int_local[i] += f[i]*w_2;
//                         }  // end of component loop i
//                     } // end of k1 loop
//                 } // end of k2 loop
//             }

//             if (m_faceOffset[5] > -1) {
// #pragma omp for nowait
//                 for (index_t k1 = bottom; k1 < bottom+m_ownNE[1]; ++k1) {
//                     for (index_t k0 = left; k0 < left+m_ownNE[0]; ++k0) {
//                         const Scalar* f = arg.getSampleDataRO(m_faceOffset[5]+INDEX2(k0,k1,m_NE[0]), zero);
//                         for (index_t i = 0; i < numComp; ++i) {
//                             int_local[i] += f[i]*w_2;
//                         }  // end of component loop i
//                     } // end of k1 loop
//                 } // end of k2 loop
//             }

// #pragma omp critical
//             for (index_t i = 0; i < numComp; i++)
//                 integrals[i] += int_local[i];
//         } // end of parallel section
//     } // function space selector
}

} // end of namespace oxley
