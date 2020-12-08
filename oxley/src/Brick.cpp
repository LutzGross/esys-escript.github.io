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
      int periodic0, int periodic1, int periodic2): OxleyDomain(3, order){

    // These two statements configure the level of verbosity used by p4est
    // sc_set_log_defaults(NULL, NULL, LOG_LEVEL);
    // p4est_init(NULL, LOG_LEVEL);

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
    const p4est_topidx_t num_vertices = 8;
    const p4est_topidx_t num_trees = 1;
    const p4est_topidx_t num_edges = 3;
    const p4est_topidx_t num_corners = 1;
    const double vertices[8 * 3] = {
                                    x0, y0, z0,
                                    x1, y0, z0,
                                    x0, y1, z0,
                                    x1, y1, z0,
                                    x0, y0, z1,
                                    x1, y0, z1,
                                    x0, y1, z1,
                                    x1, y1, z1,
                                    };
    const p4est_topidx_t tree_to_vertex[8] = {0, 1, 2, 3, 4, 5, 6, 7,};
    const p4est_topidx_t tree_to_tree[6] = {0, 0, 0, 0, 0, 0,};
    // const int8_t tree_to_face[6] = {1, 0, 3, 2, 5, 4, };
    const int8_t tree_to_face[6] = {0, 1, 2, 3, 4, 5 };
    const p4est_topidx_t tree_to_edge[12] = {0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2,};
    const p4est_topidx_t ett_offset[4] = {0, 4, 8, 12,};
    const p4est_topidx_t edge_to_tree[12] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,};
    const int8_t edge_to_edge[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11,};
    const p4est_topidx_t tree_to_corner[8] = {0, 0, 0, 0, 0, 0, 0, 0,};
    const p4est_topidx_t ctt_offset[2] = {0, 8,};
    const p4est_topidx_t corner_to_tree[8] = {0, 0, 0, 0, 0, 0, 0, 0,};
    const int8_t corner_to_corner[8] = {0, 1, 2, 3, 4, 5, 6, 7,};
    connectivity = p8est_connectivity_new_copy(num_vertices, num_trees, num_edges,
                                      num_corners, vertices, tree_to_vertex,
                                      tree_to_tree, tree_to_face,
                                      tree_to_edge, ett_offset,
                                      edge_to_tree, edge_to_edge,
                                      tree_to_corner, ctt_offset,
                                      corner_to_tree, corner_to_corner);

    if(!p8est_connectivity_is_valid(connectivity))
        throw OxleyException("Could not create a valid connectivity.");

    // Allocate some memory
    forestData = (p8estData *) malloc(sizeof(p8estData));

    // Create the p8est
    p4est_locidx_t min_quadrants = n0*n1*n2 / m_mpiInfo->size;
    int min_level = 0;
    int fill_uniform = 1;
    p8est = p8est_new_ext(m_mpiInfo->comm, connectivity, min_quadrants,
            min_level, fill_uniform, sizeof(octantData), &init_brick_data, (void *) &forestData);

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

    // //number of elements for this rank in each dimension including shared
    // forestData->m_NX[0] = d0;
    // forestData->m_NX[1] = d1;
    // forestData->m_NX[2] = d2;

    // local number of elements
    // forestData->m_NE[0] = forestData->m_gNE[0] / d0;
    // forestData->m_NE[1] = forestData->m_gNE[1] / d1;
    // forestData->m_NE[2] = forestData->m_gNE[2] / d2;
    forestData->m_NE[0] = forestData->m_gNE[0];
    forestData->m_NE[1] = forestData->m_gNE[1];
    forestData->m_NE[2] = forestData->m_gNE[2];

    // max levels of refinement
    forestData->max_levels_refinement = MAXREFINEMENTLEVELS;

    // element order
    m_order = order;

    // initial tag
    tags[0] = 0;
    numberOfTags=1;

    // Number of dimensions
    m_numDim=3;

    // Create the node numbering scheme
    // nodes_ghost = p8est_ghost_new(p8est, P8EST_CONNECT_FULL);
    // nodes = p8est_lnodes_new(p8est, NULL, 3);

    //  // Distribute the p8est across the processors
    int allow_coarsening = 0;
    p8est_partition(p8est, allow_coarsening, NULL);

    // p8est_partition_lnodes(p8est, NULL, 3, allow_coarsening);

    // To prevent segmentation faults when using numpy ndarray
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

    // free(forestData);
    p8est_connectivity_destroy(connectivity);
    p8est_destroy(p8est);
    // p8est_lnodes_destroy(nodes);
    // p8est_ghost_destroy(nodes_ghost);
    // sc_finalize();
}

/**
   \brief
   returns a description for this domain
*/
std::string Brick::getDescription() const{

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
    throw OxleyException("currently: not supported"); //AE this is temporary
}


const dim_t* Brick::borrowSampleReferenceIDs(int fsType) const
{
    throw OxleyException("currently: not supported"); //AE this is temporary
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

        // Write the header
        context = p8est_vtk_write_header(context);

        // Get the point and cell data together
        p4est_locidx_t numquads = p8est->local_num_quadrants;
        sc_array_t * quadTag = sc_array_new_count(sizeof(double), numquads);
        p8est_iterate(p8est, NULL, (void *) quadTag, getQuadTagVector, NULL, NULL, NULL);
        sc_array_t * xCoords = sc_array_new_count(sizeof(double), numquads);
        p8est_iterate(p8est, NULL, (void *) xCoords, getXCoordVector, NULL, NULL, NULL);
        sc_array_t * yCoords = sc_array_new_count(sizeof(double), numquads);
        p8est_iterate(p8est, NULL, (void *) yCoords, getYCoordVector, NULL, NULL, NULL);
        sc_array_t * zCoords = sc_array_new_count(sizeof(double), numquads);
        p8est_iterate(p8est, NULL, (void *) zCoords, getZCoordVector, NULL, NULL, NULL);

        // Cell Data
#ifdef OXLEY_ENABLE_DEBUG
        // context = p8est_vtk_write_cell_dataf(context,1,1,0,0,1,0,"tag",quadTag,context);
        context = p8est_vtk_write_cell_dataf(context,1,1,0,0,4,0,
            "tag",quadTag,"x",xCoords,"y",yCoords,"z",zCoords,context);
#else
        // context = p8est_vtk_write_cell_dataf(context,0,0,0,0,1,0,"tag",quadTag,context);
        context = p8est_vtk_write_cell_dataf(context,0,0,0,0,4,0,
            "tag",quadTag,"x",xCoords,"y",yCoords,"z",zCoords,context);
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
        sc_array_reset(quadTag);
        sc_array_destroy(quadTag);
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

void Brick::refineBoundary(std::string boundary, double dx)
{

}

// protected
inline index_t Brick::getFirstInDim(unsigned axis) const
{
    // return m_offset[axis] == 0 ? 0 : 1;
    return nodes->global_offset == 0 ? 0 : 1;
}

//protected
inline dim_t Brick::getNumNodes() const
{
    // return m_NN[0]*m_NN[1];
    return nodes->num_local_nodes + *nodes->nonlocal_nodes;
}

//protected
inline dim_t Brick::getNumElements() const
{
    // return m_NE[0]*m_NE[1];
}

//protected
inline dim_t Brick::getNumFaceElements() const
{
    // return m_faceCount[0] + m_faceCount[1] + m_faceCount[2] + m_faceCount[3];
}

//protected
inline dim_t Brick::getNumDOF() const
{
    // return (m_gNE[0]+1)/m_NX[0]*(m_gNE[1]+1)/m_NX[1];
    return getNumNodes()+1; //TODO
}

//protected
void Brick::assembleCoordinates(escript::Data& arg) const
{
    int numDim = m_numDim;
    if (!arg.isDataPointShapeEqual(1, &numDim))
        throw ValueError("setToX: Invalid Data object shape");
    if (!arg.numSamplesEqual(1, getNumNodes()))
        throw ValueError("setToX: Illegal number of samples in Data object");
    arg.requireWrite();
    for(int treeid = p8est->first_local_tree; treeid < p8est->last_local_tree; treeid++)
    {
        p8est_tree_t * currenttree = p8est_tree_array_index(p8est->trees, treeid);
        sc_array_t * tquadrants = &currenttree->quadrants;
        p4est_locidx_t Q = (p4est_locidx_t) tquadrants->elem_count;
#pragma omp parallel for
        for (p4est_locidx_t e = nodes->global_offset; e < Q+nodes->global_offset; e++) // Loop over every quadrant within the tree
        {
            p8est_quadrant_t * quad = p8est_quadrant_array_index(tquadrants, e);
            double* point = arg.getSampleDataRW(e);
            p8est_qcoord_to_vertex(connectivity, treeid, quad->x, quad->y, quad->z, &point[0]);
        }
    }        
}

#ifdef ESYS_HAVE_TRILINOS
//protected
esys_trilinos::const_TrilinosGraph_ptr Brick::getTrilinosGraph() const
{
    //TODO

    // if (m_graph.is_null()) {
    //     m_graph = createTrilinosGraph(m_dofId, m_nodeId);
    // }
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

//     if (m_pattern.get())
//         return m_pattern;

//     // first call to this method -> create the pattern, then return it
//     paso::Connector_ptr conn(getPasoConnector());
//     const dim_t nDOF0 = (m_gNE[0]+1)/m_NX[0];
//     const dim_t nDOF1 = (m_gNE[1]+1)/m_NX[1];
//     const dim_t nDOF2 = (m_gNE[2]+1)/m_NX[2];
//     const dim_t numDOF = nDOF0*nDOF1*nDOF2;
//     const dim_t numShared = conn->send->numSharedComponents;
//     const index_t* sendShared = conn->send->shared;
//     const int x = m_mpiInfo->rank%m_NX[0];
//     const int y = m_mpiInfo->rank%(m_NX[0]*m_NX[1])/m_NX[0];
//     const int z = m_mpiInfo->rank/(m_NX[0]*m_NX[1]);

//     // these are for the couple blocks
//     std::vector<IndexVector> colIndices(numDOF);
//     std::vector<IndexVector> rowIndices(numShared);

//     for (dim_t i=0; i < conn->send->neighbour.size(); i++) {
//         const dim_t start = conn->send->offsetInShared[i];
//         const dim_t end = conn->send->offsetInShared[i+1];
//         // location of neighbour rank relative to this rank
//         const int xDiff = conn->send->neighbour[i]%m_NX[0] - x;
//         const int yDiff = conn->send->neighbour[i]%(m_NX[0]*m_NX[1])/m_NX[0] - y;
//         const int zDiff = conn->send->neighbour[i]/(m_NX[0]*m_NX[1]) - z;
        
//         if (xDiff==0 && yDiff==0) {
//             // sharing front or back plane
//             for (dim_t j = start; j < end; j++) {
//                 const dim_t i0 = (j-start)%nDOF0;
//                 const dim_t i1 = (j-start)/nDOF0;
//                 if (i0 > 0) {
//                     if (i1 > 0)
//                         doublyLink(colIndices, rowIndices, sendShared[j-1-nDOF0], j);
//                     doublyLink(colIndices, rowIndices, sendShared[j-1], j);
//                     if (i1 < nDOF1-1)
//                         doublyLink(colIndices, rowIndices, sendShared[j-1+nDOF0], j);
//                 }
//                 if (i1 > 0)
//                     doublyLink(colIndices, rowIndices, sendShared[j-nDOF0], j);
//                 doublyLink(colIndices, rowIndices, sendShared[j], j);
//                 if (i1 < nDOF1-1)
//                     doublyLink(colIndices, rowIndices, sendShared[j+nDOF0], j);
//                 if (i0 < nDOF0-1) {
//                     if (i1 > 0)
//                         doublyLink(colIndices, rowIndices, sendShared[j+1-nDOF0], j);
//                     doublyLink(colIndices, rowIndices, sendShared[j+1], j);
//                     if (i1 < nDOF1-1)
//                         doublyLink(colIndices, rowIndices, sendShared[j+1+nDOF0], j);
//                 }
//             }
//         } else if (xDiff==0 && zDiff==0) {
//             // sharing top or bottom plane
//             for (dim_t j = start; j < end; j++) {
//                 const dim_t i0 = (j-start)%nDOF0;
//                 const dim_t i1 = (j-start)/nDOF0;
//                 if (i0 > 0) {
//                     if (i1 > 0)
//                         doublyLink(colIndices, rowIndices, sendShared[j]-1-nDOF0*nDOF1, j);
//                     doublyLink(colIndices, rowIndices, sendShared[j]-1, j);
//                     if (i1 < nDOF2-1)
//                         doublyLink(colIndices, rowIndices, sendShared[j]-1+nDOF0*nDOF1, j);
//                 }
//                 if (i1 > 0)
//                     doublyLink(colIndices, rowIndices, sendShared[j]-nDOF0*nDOF1, j);
//                 doublyLink(colIndices, rowIndices, sendShared[j], j);
//                 if (i1 < nDOF2-1)
//                     doublyLink(colIndices, rowIndices, sendShared[j]+nDOF0*nDOF1, j);
//                 if (i0 < nDOF0-1) {
//                     if (i1 > 0)
//                         doublyLink(colIndices, rowIndices, sendShared[j]+1-nDOF0*nDOF1, j);
//                     doublyLink(colIndices, rowIndices, sendShared[j]+1, j);
//                     if (i1 < nDOF2-1)
//                         doublyLink(colIndices, rowIndices, sendShared[j]+1+nDOF0*nDOF1, j);
//                 }
//             }
//         } else if (yDiff==0 && zDiff==0) {
//             // sharing left or right plane
//             for (dim_t j = start; j < end; j++) {
//                 const dim_t i0 = (j-start)%nDOF1;
//                 const dim_t i1 = (j-start)/nDOF1;
//                 if (i0 > 0) {
//                     if (i1 > 0)
//                         doublyLink(colIndices, rowIndices, sendShared[j]-nDOF0-nDOF0*nDOF1, j);
//                     doublyLink(colIndices, rowIndices, sendShared[j]-nDOF0, j);
//                     if (i1 < nDOF2-1)
//                         doublyLink(colIndices, rowIndices, sendShared[j]-nDOF0+nDOF0*nDOF1, j);
//                 }
//                 if (i1 > 0)
//                     doublyLink(colIndices, rowIndices, sendShared[j]-nDOF0*nDOF1, j);
//                 doublyLink(colIndices, rowIndices, sendShared[j], j);
//                 if (i1 < nDOF2-1)
//                     doublyLink(colIndices, rowIndices, sendShared[j]+nDOF0*nDOF1, j);
//                 if (i0 < nDOF1-1) {
//                     if (i1 > 0)
//                         doublyLink(colIndices, rowIndices, sendShared[j]+nDOF0-nDOF0*nDOF1, j);
//                     doublyLink(colIndices, rowIndices, sendShared[j]+nDOF0, j);
//                     if (i1 < nDOF2-1)
//                         doublyLink(colIndices, rowIndices, sendShared[j]+nDOF0+nDOF0*nDOF1, j);
//                 }
//             }
//         } else if (xDiff == 0) {
//             // sharing an edge in x direction
//             for (dim_t j = start; j < end; j++) {
//                 if (j > start)
//                     doublyLink(colIndices, rowIndices, sendShared[j]-1, j);
//                 doublyLink(colIndices, rowIndices, sendShared[j], j);
//                 if (j < end-1)
//                     doublyLink(colIndices, rowIndices, sendShared[j]+1, j);
//             }
//         } else if (yDiff == 0) {
//             // sharing an edge in y direction
//             for (dim_t j = start; j < end; j++) {
//                 if (j > start)
//                     doublyLink(colIndices, rowIndices, sendShared[j]-nDOF0, j);
//                 doublyLink(colIndices, rowIndices, sendShared[j], j);
//                 if (j < end-1)
//                     doublyLink(colIndices, rowIndices, sendShared[j]+nDOF0, j);
//             }
//         } else if (zDiff == 0) {
//             // sharing an edge in z direction
//             for (dim_t j = start; j < end; j++) {
//                 if (j > start)
//                     doublyLink(colIndices, rowIndices, sendShared[j]-nDOF0*nDOF1, j);
//                 doublyLink(colIndices, rowIndices, sendShared[j], j);
//                 if (j < end-1)
//                     doublyLink(colIndices, rowIndices, sendShared[j]+nDOF0*nDOF1, j);
//             }
//         } else {
//             // sharing a node
//             doublyLink(colIndices, rowIndices, sendShared[start], start);
//         }
//     }

// #pragma omp parallel for
//     for (dim_t i = 0; i < numShared; i++) {
//         sort(rowIndices[i].begin(), rowIndices[i].end());
//     }

//     // create main and couple blocks
//     paso::Pattern_ptr mainPattern = createPasoPattern(getConnections(), numDOF);
//     paso::Pattern_ptr colPattern = createPasoPattern(colIndices, numShared);
//     paso::Pattern_ptr rowPattern = createPasoPattern(rowIndices, numDOF);

//     // allocate paso distribution
//     escript::Distribution_ptr distribution(new escript::Distribution(
//                                                m_mpiInfo, m_nodeDistribution));

//     // finally create the system matrix pattern
//     m_pattern.reset(new paso::SystemMatrixPattern(MATRIX_FORMAT_DEFAULT,
//             distribution, distribution, mainPattern, colPattern, rowPattern,
//             conn, conn));

//     // useful debug output
//     /*
//     std::cout << "--- colIndices ---" << std::endl;
//     for (size_t i=0; i<colIndices.size(); i++) {
//         std::cout << "colIndices[" << i << "].size()=" << colIndices[i].size() << std::endl;
//     }
//     std::cout << "--- rowIndices ---" << std::endl;
//     for (size_t i=0; i<rowIndices.size(); i++) {
//         std::cout << "rowIndices[" << i << "].size()=" << rowIndices[i].size() << std::endl;
//     }
//     */
//     /*
//     std::cout << "--- main_pattern ---" << std::endl;
//     std::cout << "M=" << mainPattern->numOutput << ", N=" << mainPattern->numInput << std::endl;
//     for (size_t i=0; i<mainPattern->numOutput+1; i++) {
//         std::cout << "ptr[" << i << "]=" << mainPattern->ptr[i] << std::endl;
//     }
//     for (size_t i=0; i<mainPattern->ptr[mainPattern->numOutput]; i++) {
//         std::cout << "index[" << i << "]=" << mainPattern->index[i] << std::endl;
//     }
//     */
//     /*
//     std::cout << "--- colCouple_pattern ---" << std::endl;
//     std::cout << "M=" << colPattern->numOutput << ", N=" << colPattern->numInput << std::endl;
//     for (size_t i=0; i<colPattern->numOutput+1; i++) {
//         std::cout << "ptr[" << i << "]=" << colPattern->ptr[i] << std::endl;
//     }
//     for (size_t i=0; i<colPattern->ptr[colPattern->numOutput]; i++) {
//         std::cout << "index[" << i << "]=" << colPattern->index[i] << std::endl;
//     }
//     */
//     /*
//     std::cout << "--- rowCouple_pattern ---" << std::endl;
//     std::cout << "M=" << rowPattern->numOutput << ", N=" << rowPattern->numInput << std::endl;
//     for (size_t i=0; i<rowPattern->numOutput+1; i++) {
//         std::cout << "ptr[" << i << "]=" << rowPattern->ptr[i] << std::endl;
//     }
//     for (size_t i=0; i<rowPattern->ptr[rowPattern->numOutput]; i++) {
//         std::cout << "index[" << i << "]=" << rowPattern->index[i] << std::endl;
//     }
//     */

//     return m_pattern;
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
    if (out.isComplex()!=in.isComplex())
    {
        throw OxleyException("Programmer Error: in and out parameters do not have the same complexity.");   
    }
    if (out.isComplex())
    {
        interpolateNodesOnFacesWorker(out, in, reduced, escript::DataTypes::cplx_t(0));    
    }
    else
    {
        interpolateNodesOnFacesWorker(out, in, reduced, escript::DataTypes::real_t(0));    
    }      
}

//private
template <typename S>
void Brick::interpolateNodesOnFacesWorker(escript::Data& out, const escript::Data& in,
                                    bool reduced, S sentinel) const
{
    //TODO


//     const dim_t numComp = in.getDataPointSize();
//     if (reduced) {
//         out.requireWrite();
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
//                             o[INDEX2(i,numComp,0)] = (f_000[i] + f_001[i] + f_010[i] + f_011[i])/static_cast<S>(4);
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
//                             o[INDEX2(i,numComp,0)] = (f_100[i] + f_101[i] + f_110[i] + f_111[i])/static_cast<S>(4);
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
//                             o[INDEX2(i,numComp,0)] = (f_000[i] + f_001[i] + f_100[i] + f_101[i])/static_cast<S>(4);
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
//                             o[INDEX2(i,numComp,0)] = (f_010[i] + f_011[i] + f_110[i] + f_111[i])/static_cast<S>(4);
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
//                             o[INDEX2(i,numComp,0)] = (f_000[i] + f_010[i] + f_100[i] + f_110[i])/static_cast<S>(4);
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
//                             o[INDEX2(i,numComp,0)] = (f_001[i] + f_011[i] + f_101[i] + f_111[i])/static_cast<S>(4);
//                         } // end of component loop i
//                     } // end of k0 loop
//                 } // end of k1 loop
//             } // end of face 5
//         } // end of parallel section
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
//     const dim_t numComp = in.getDataPointSize();
//     if (reduced) {
//         out.requireWrite();
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
// #pragma omp for
//             for (index_t k2=0; k2 < m_NE[2]; ++k2) {
//                 for (index_t k1=0; k1 < m_NE[1]; ++k1) {
//                     for (index_t k0=0; k0 < m_NE[0]; ++k0) {
//                         memcpy(&f_000[0], in.getSampleDataRO(INDEX3(k0,k1,k2, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
//                         memcpy(&f_001[0], in.getSampleDataRO(INDEX3(k0,k1,k2+1, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
//                         memcpy(&f_010[0], in.getSampleDataRO(INDEX3(k0,k1+1,k2, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
//                         memcpy(&f_011[0], in.getSampleDataRO(INDEX3(k0,k1+1,k2+1, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
//                         memcpy(&f_100[0], in.getSampleDataRO(INDEX3(k0+1,k1,k2, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
//                         memcpy(&f_101[0], in.getSampleDataRO(INDEX3(k0+1,k1,k2+1, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
//                         memcpy(&f_110[0], in.getSampleDataRO(INDEX3(k0+1,k1+1,k2, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
//                         memcpy(&f_111[0], in.getSampleDataRO(INDEX3(k0+1,k1+1,k2+1, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
//                         S* o = out.getSampleDataRW(INDEX3(k0,k1,k2,m_NE[0],m_NE[1]), sentinel);
//                         for (index_t i=0; i < numComp; ++i) {
//                             o[INDEX2(i,numComp,0)] = (f_000[i] + f_001[i] + f_010[i] + f_011[i] + f_100[i] + f_101[i] + f_110[i] + f_111[i])/static_cast<S>(8);
//                         } // end of component loop i
//                     } // end of k0 loop
//                 } // end of k1 loop
//             } // end of k2 loop
//         } // end of parallel section
//     } else {
//         out.requireWrite();
//         const S c0 = .0094373878376559314545;
//         const S c1 = .035220810900864519624;
//         const S c2 = .13144585576580214704;
//         const S c3 = .49056261216234406855;
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
// #pragma omp for
//             for (index_t k2=0; k2 < m_NE[2]; ++k2) {
//                 for (index_t k1=0; k1 < m_NE[1]; ++k1) {
//                     for (index_t k0=0; k0 < m_NE[0]; ++k0) {
//                         memcpy(&f_000[0], in.getSampleDataRO(INDEX3(k0,k1,k2, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
//                         memcpy(&f_001[0], in.getSampleDataRO(INDEX3(k0,k1,k2+1, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
//                         memcpy(&f_010[0], in.getSampleDataRO(INDEX3(k0,k1+1,k2, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
//                         memcpy(&f_011[0], in.getSampleDataRO(INDEX3(k0,k1+1,k2+1, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
//                         memcpy(&f_100[0], in.getSampleDataRO(INDEX3(k0+1,k1,k2, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
//                         memcpy(&f_101[0], in.getSampleDataRO(INDEX3(k0+1,k1,k2+1, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
//                         memcpy(&f_110[0], in.getSampleDataRO(INDEX3(k0+1,k1+1,k2, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
//                         memcpy(&f_111[0], in.getSampleDataRO(INDEX3(k0+1,k1+1,k2+1, m_NN[0],m_NN[1]), sentinel), numComp*sizeof(S));
//                         S* o = out.getSampleDataRW(INDEX3(k0,k1,k2,m_NE[0],m_NE[1]), sentinel);
//                         for (index_t i=0; i < numComp; ++i) {
//                             o[INDEX2(i,numComp,0)] = f_000[i]*c3 + f_111[i]*c0 + c2*(f_001[i] + f_010[i] + f_100[i]) + c1*(f_011[i] + f_101[i] + f_110[i]);
//                             o[INDEX2(i,numComp,1)] = f_011[i]*c0 + f_100[i]*c3 + c2*(f_000[i] + f_101[i] + f_110[i]) + c1*(f_001[i] + f_010[i] + f_111[i]);
//                             o[INDEX2(i,numComp,2)] = f_010[i]*c3 + f_101[i]*c0 + c2*(f_000[i] + f_011[i] + f_110[i]) + c1*(f_001[i] + f_100[i] + f_111[i]);
//                             o[INDEX2(i,numComp,3)] = f_001[i]*c0 + f_110[i]*c3 + c2*(f_010[i] + f_100[i] + f_111[i]) + c1*(f_000[i] + f_011[i] + f_101[i]);
//                             o[INDEX2(i,numComp,4)] = f_001[i]*c3 + f_110[i]*c0 + c2*(f_000[i] + f_011[i] + f_101[i]) + c1*(f_010[i] + f_100[i] + f_111[i]);
//                             o[INDEX2(i,numComp,5)] = f_010[i]*c0 + f_101[i]*c3 + c2*(f_001[i] + f_100[i] + f_111[i]) + c1*(f_000[i] + f_011[i] + f_110[i]);
//                             o[INDEX2(i,numComp,6)] = f_011[i]*c3 + f_100[i]*c0 + c2*(f_001[i] + f_010[i] + f_111[i]) + c1*(f_000[i] + f_101[i] + f_110[i]);
//                             o[INDEX2(i,numComp,7)] = f_000[i]*c0 + f_111[i]*c3 + c2*(f_011[i] + f_101[i] + f_110[i]) + c1*(f_001[i] + f_010[i] + f_100[i]);
//                         } // end of component loop i
//                     } // end of k0 loop
//                 } // end of k1 loop
//             } // end of k2 loop
//         } // end of parallel section
//     }
}

//protected
void Brick::assembleGradient(escript::Data& out,
                             const escript::Data& in) const
{
    if (out.isComplex() != in.isComplex())
        throw ValueError("Gradient: input & output complexity must match.");
    else if (in.isComplex())
        assembleGradientImpl<cplx_t>(out, in);
    else
        assembleGradientImpl<real_t>(out, in);
}

//protected
template<typename Scalar>
void Brick::assembleGradientImpl(escript::Data& out,
                                 const escript::Data& in) const
{
    //TODO

//     const dim_t numComp = in.getDataPointSize();
//     const double C0 = .044658198738520451079;
//     const double C1 = .16666666666666666667;
//     const double C2 = .21132486540518711775;
//     const double C3 = .25;
//     const double C4 = .5;
//     const double C5 = .62200846792814621559;
//     const double C6 = .78867513459481288225;
//     const dim_t NE0 = m_NE[0];
//     const dim_t NE1 = m_NE[1];
//     const dim_t NE2 = m_NE[2];
//     const Scalar zero = static_cast<Scalar>(0);

//     if (out.getFunctionSpace().getTypeCode() == Elements) {
//         out.requireWrite();
// #pragma omp parallel
//         {
//             vector<Scalar> f_000(numComp, zero);
//             vector<Scalar> f_001(numComp, zero);
//             vector<Scalar> f_010(numComp, zero);
//             vector<Scalar> f_011(numComp, zero);
//             vector<Scalar> f_100(numComp, zero);
//             vector<Scalar> f_101(numComp, zero);
//             vector<Scalar> f_110(numComp, zero);
//             vector<Scalar> f_111(numComp, zero);
// #pragma omp for
//             for (index_t k2=0; k2 < NE2; ++k2) {
//                 for (index_t k1=0; k1 < NE1; ++k1) {
//                     for (index_t k0=0; k0 < NE0; ++k0) {
//                         memcpy(&f_000[0], in.getSampleDataRO(INDEX3(k0,k1,k2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_001[0], in.getSampleDataRO(INDEX3(k0,k1,k2+1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_010[0], in.getSampleDataRO(INDEX3(k0,k1+1,k2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_011[0], in.getSampleDataRO(INDEX3(k0,k1+1,k2+1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_100[0], in.getSampleDataRO(INDEX3(k0+1,k1,k2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_101[0], in.getSampleDataRO(INDEX3(k0+1,k1,k2+1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_110[0], in.getSampleDataRO(INDEX3(k0+1,k1+1,k2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_111[0], in.getSampleDataRO(INDEX3(k0+1,k1+1,k2+1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         Scalar* o = out.getSampleDataRW(INDEX3(k0,k1,k2,NE0,NE1), zero);
//                         for (index_t i=0; i < numComp; ++i) {
//                             const Scalar V0=((f_100[i]-f_000[i])*C5 + (f_111[i]-f_011[i])*C0 + (f_101[i]+f_110[i]-f_001[i]-f_010[i])*C1) / m_dx[0];
//                             const Scalar V1=((f_110[i]-f_010[i])*C5 + (f_101[i]-f_001[i])*C0 + (f_100[i]+f_111[i]-f_000[i]-f_011[i])*C1) / m_dx[0];
//                             const Scalar V2=((f_101[i]-f_001[i])*C5 + (f_110[i]-f_010[i])*C0 + (f_100[i]+f_111[i]-f_000[i]-f_011[i])*C1) / m_dx[0];
//                             const Scalar V3=((f_111[i]-f_011[i])*C5 + (f_100[i]-f_000[i])*C0 + (f_101[i]+f_110[i]-f_001[i]-f_010[i])*C1) / m_dx[0];
//                             const Scalar V4=((f_010[i]-f_000[i])*C5 + (f_111[i]-f_101[i])*C0 + (f_011[i]+f_110[i]-f_001[i]-f_100[i])*C1) / m_dx[1];
//                             const Scalar V5=((f_110[i]-f_100[i])*C5 + (f_011[i]-f_001[i])*C0 + (f_010[i]+f_111[i]-f_000[i]-f_101[i])*C1) / m_dx[1];
//                             const Scalar V6=((f_011[i]-f_001[i])*C5 + (f_110[i]-f_100[i])*C0 + (f_010[i]+f_111[i]-f_000[i]-f_101[i])*C1) / m_dx[1];
//                             const Scalar V7=((f_111[i]-f_101[i])*C5 + (f_010[i]-f_000[i])*C0 + (f_011[i]+f_110[i]-f_001[i]-f_100[i])*C1) / m_dx[1];
//                             const Scalar V8=((f_001[i]-f_000[i])*C5 + (f_111[i]-f_110[i])*C0 + (f_011[i]+f_101[i]-f_010[i]-f_100[i])*C1) / m_dx[2];
//                             const Scalar V9=((f_101[i]-f_100[i])*C5 + (f_011[i]-f_010[i])*C0 + (f_001[i]+f_111[i]-f_000[i]-f_110[i])*C1) / m_dx[2];
//                             const Scalar V10=((f_011[i]-f_010[i])*C5 + (f_101[i]-f_100[i])*C0 + (f_001[i]+f_111[i]-f_000[i]-f_110[i])*C1) / m_dx[2];
//                             const Scalar V11=((f_111[i]-f_110[i])*C5 + (f_001[i]-f_000[i])*C0 + (f_011[i]+f_101[i]-f_010[i]-f_100[i])*C1) / m_dx[2];
//                             o[INDEX3(i,0,0,numComp,3)] = V0;
//                             o[INDEX3(i,1,0,numComp,3)] = V4;
//                             o[INDEX3(i,2,0,numComp,3)] = V8;
//                             o[INDEX3(i,0,1,numComp,3)] = V0;
//                             o[INDEX3(i,1,1,numComp,3)] = V5;
//                             o[INDEX3(i,2,1,numComp,3)] = V9;
//                             o[INDEX3(i,0,2,numComp,3)] = V1;
//                             o[INDEX3(i,1,2,numComp,3)] = V4;
//                             o[INDEX3(i,2,2,numComp,3)] = V10;
//                             o[INDEX3(i,0,3,numComp,3)] = V1;
//                             o[INDEX3(i,1,3,numComp,3)] = V5;
//                             o[INDEX3(i,2,3,numComp,3)] = V11;
//                             o[INDEX3(i,0,4,numComp,3)] = V2;
//                             o[INDEX3(i,1,4,numComp,3)] = V6;
//                             o[INDEX3(i,2,4,numComp,3)] = V8;
//                             o[INDEX3(i,0,5,numComp,3)] = V2;
//                             o[INDEX3(i,1,5,numComp,3)] = V7;
//                             o[INDEX3(i,2,5,numComp,3)] = V9;
//                             o[INDEX3(i,0,6,numComp,3)] = V3;
//                             o[INDEX3(i,1,6,numComp,3)] = V6;
//                             o[INDEX3(i,2,6,numComp,3)] = V10;
//                             o[INDEX3(i,0,7,numComp,3)] = V3;
//                             o[INDEX3(i,1,7,numComp,3)] = V7;
//                             o[INDEX3(i,2,7,numComp,3)] = V11;
//                         } // end of component loop i
//                     } // end of k0 loop
//                 } // end of k1 loop
//             } // end of k2 loop
//         } // end of parallel section
//     } else if (out.getFunctionSpace().getTypeCode() == ReducedElements) {
//         out.requireWrite();
// #pragma omp parallel
//         {
//             vector<Scalar> f_000(numComp, zero);
//             vector<Scalar> f_001(numComp, zero);
//             vector<Scalar> f_010(numComp, zero);
//             vector<Scalar> f_011(numComp, zero);
//             vector<Scalar> f_100(numComp, zero);
//             vector<Scalar> f_101(numComp, zero);
//             vector<Scalar> f_110(numComp, zero);
//             vector<Scalar> f_111(numComp, zero);
// #pragma omp for
//             for (index_t k2=0; k2 < NE2; ++k2) {
//                 for (index_t k1=0; k1 < NE1; ++k1) {
//                     for (index_t k0=0; k0 < NE0; ++k0) {
//                         memcpy(&f_000[0], in.getSampleDataRO(INDEX3(k0,k1,k2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_001[0], in.getSampleDataRO(INDEX3(k0,k1,k2+1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_010[0], in.getSampleDataRO(INDEX3(k0,k1+1,k2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_011[0], in.getSampleDataRO(INDEX3(k0,k1+1,k2+1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_100[0], in.getSampleDataRO(INDEX3(k0+1,k1,k2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_101[0], in.getSampleDataRO(INDEX3(k0+1,k1,k2+1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_110[0], in.getSampleDataRO(INDEX3(k0+1,k1+1,k2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_111[0], in.getSampleDataRO(INDEX3(k0+1,k1+1,k2+1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         Scalar* o = out.getSampleDataRW(INDEX3(k0,k1,k2,NE0,NE1), zero);
//                         for (index_t i=0; i < numComp; ++i) {
//                             o[INDEX3(i,0,0,numComp,3)] = (f_100[i]+f_101[i]+f_110[i]+f_111[i]-f_000[i]-f_001[i]-f_010[i]-f_011[i])*C3 / m_dx[0];
//                             o[INDEX3(i,1,0,numComp,3)] = (f_010[i]+f_011[i]+f_110[i]+f_111[i]-f_000[i]-f_001[i]-f_100[i]-f_101[i])*C3 / m_dx[1];
//                             o[INDEX3(i,2,0,numComp,3)] = (f_001[i]+f_011[i]+f_101[i]+f_111[i]-f_000[i]-f_010[i]-f_100[i]-f_110[i])*C3 / m_dx[2];
//                         } // end of component loop i
//                     } // end of k0 loop
//                 } // end of k1 loop
//             } // end of k2 loop
//         } // end of parallel section
//     } else if (out.getFunctionSpace().getTypeCode() == FaceElements) {
//         out.requireWrite();
// #pragma omp parallel
//         {
//             vector<Scalar> f_000(numComp, zero);
//             vector<Scalar> f_001(numComp, zero);
//             vector<Scalar> f_010(numComp, zero);
//             vector<Scalar> f_011(numComp, zero);
//             vector<Scalar> f_100(numComp, zero);
//             vector<Scalar> f_101(numComp, zero);
//             vector<Scalar> f_110(numComp, zero);
//             vector<Scalar> f_111(numComp, zero);
//             if (m_faceOffset[0] > -1) {
// #pragma omp for nowait
//                 for (index_t k2=0; k2 < NE2; ++k2) {
//                     for (index_t k1=0; k1 < NE1; ++k1) {
//                         memcpy(&f_000[0], in.getSampleDataRO(INDEX3(0,k1,k2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_001[0], in.getSampleDataRO(INDEX3(0,k1,k2+1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_010[0], in.getSampleDataRO(INDEX3(0,k1+1,k2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_011[0], in.getSampleDataRO(INDEX3(0,k1+1,k2+1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_100[0], in.getSampleDataRO(INDEX3(1,k1,k2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_101[0], in.getSampleDataRO(INDEX3(1,k1,k2+1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_110[0], in.getSampleDataRO(INDEX3(1,k1+1,k2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_111[0], in.getSampleDataRO(INDEX3(1,k1+1,k2+1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         Scalar* o = out.getSampleDataRW(m_faceOffset[0]+INDEX2(k1,k2,NE1), zero);
//                         for (index_t i=0; i < numComp; ++i) {
//                             const Scalar V0=((f_010[i]-f_000[i])*C6 + (f_011[i]-f_001[i])*C2) / m_dx[1];
//                             const Scalar V1=((f_010[i]-f_000[i])*C2 + (f_011[i]-f_001[i])*C6) / m_dx[1];
//                             const Scalar V2=((f_001[i]-f_000[i])*C6 + (f_010[i]-f_011[i])*C2) / m_dx[2];
//                             const Scalar V3=((f_001[i]-f_000[i])*C2 + (f_011[i]-f_010[i])*C6) / m_dx[2];
//                             o[INDEX3(i,0,0,numComp,3)] = ((f_100[i]-f_000[i])*C5 + (f_111[i]-f_011[i])*C0 + (f_101[i]+f_110[i]-f_001[i]-f_010[i])*C1) / m_dx[0];
//                             o[INDEX3(i,1,0,numComp,3)] = V0;
//                             o[INDEX3(i,2,0,numComp,3)] = V2;
//                             o[INDEX3(i,0,1,numComp,3)] = ((f_110[i]-f_010[i])*C5 + (f_101[i]-f_001[i])*C0 + (f_100[i]+f_111[i]-f_000[i]-f_011[i])*C1) / m_dx[0];
//                             o[INDEX3(i,1,1,numComp,3)] = V0;
//                             o[INDEX3(i,2,1,numComp,3)] = V3;
//                             o[INDEX3(i,0,2,numComp,3)] = ((f_101[i]-f_001[i])*C5 + (f_110[i]-f_010[i])*C0 + (f_100[i]+f_111[i]-f_000[i]-f_011[i])*C1) / m_dx[0];
//                             o[INDEX3(i,1,2,numComp,3)] = V1;
//                             o[INDEX3(i,2,2,numComp,3)] = V2;
//                             o[INDEX3(i,0,3,numComp,3)] = ((f_111[i]-f_011[i])*C5 + (f_100[i]-f_000[i])*C0 + (f_101[i]+f_110[i]-f_001[i]-f_010[i])*C1) / m_dx[0];
//                             o[INDEX3(i,1,3,numComp,3)] = V1;
//                             o[INDEX3(i,2,3,numComp,3)] = V3;
//                         } // end of component loop i
//                     } // end of k1 loop
//                 } // end of k2 loop
//             } // end of face 0
//             if (m_faceOffset[1] > -1) {
// #pragma omp for nowait
//                 for (index_t k2=0; k2 < NE2; ++k2) {
//                     for (index_t k1=0; k1 < NE1; ++k1) {
//                         memcpy(&f_000[0], in.getSampleDataRO(INDEX3(m_NN[0]-2,k1,k2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_001[0], in.getSampleDataRO(INDEX3(m_NN[0]-2,k1,k2+1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_010[0], in.getSampleDataRO(INDEX3(m_NN[0]-2,k1+1,k2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_011[0], in.getSampleDataRO(INDEX3(m_NN[0]-2,k1+1,k2+1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_100[0], in.getSampleDataRO(INDEX3(m_NN[0]-1,k1,k2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_101[0], in.getSampleDataRO(INDEX3(m_NN[0]-1,k1,k2+1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_110[0], in.getSampleDataRO(INDEX3(m_NN[0]-1,k1+1,k2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_111[0], in.getSampleDataRO(INDEX3(m_NN[0]-1,k1+1,k2+1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         Scalar* o = out.getSampleDataRW(m_faceOffset[1]+INDEX2(k1,k2,NE1), zero);
//                         for (index_t i=0; i < numComp; ++i) {
//                             const Scalar V0=((f_110[i]-f_100[i])*C6 + (f_111[i]-f_101[i])*C2) / m_dx[1];
//                             const Scalar V1=((f_110[i]-f_100[i])*C2 + (f_111[i]-f_101[i])*C6) / m_dx[1];
//                             const Scalar V2=((f_101[i]-f_100[i])*C6 + (f_111[i]-f_110[i])*C2) / m_dx[2];
//                             const Scalar V3=((f_101[i]-f_100[i])*C2 + (f_111[i]-f_110[i])*C6) / m_dx[2];
//                             o[INDEX3(i,0,0,numComp,3)] = ((f_100[i]-f_000[i])*C5 + (f_111[i]-f_011[i])*C0 + (f_101[i]+f_110[i]-f_001[i]-f_010[i])*C1) / m_dx[0];
//                             o[INDEX3(i,1,0,numComp,3)] = V0;
//                             o[INDEX3(i,2,0,numComp,3)] = V2;
//                             o[INDEX3(i,0,1,numComp,3)] = ((f_110[i]-f_010[i])*C5 + (f_101[i]-f_001[i])*C0 + (f_100[i]+f_111[i]-f_000[i]-f_011[i])*C1) / m_dx[0];
//                             o[INDEX3(i,1,1,numComp,3)] = V0;
//                             o[INDEX3(i,2,1,numComp,3)] = V3;
//                             o[INDEX3(i,0,2,numComp,3)] = ((f_101[i]-f_001[i])*C5 + (f_110[i]-f_010[i])*C0 + (f_100[i]+f_111[i]-f_000[i]-f_011[i])*C1) / m_dx[0];
//                             o[INDEX3(i,1,2,numComp,3)] = V1;
//                             o[INDEX3(i,2,2,numComp,3)] = V2;
//                             o[INDEX3(i,0,3,numComp,3)] = ((f_111[i]-f_011[i])*C5 + (f_100[i]-f_000[i])*C0 + (f_101[i]+f_110[i]-f_001[i]-f_010[i])*C1) / m_dx[0];
//                             o[INDEX3(i,1,3,numComp,3)] = V1;
//                             o[INDEX3(i,2,3,numComp,3)] = V3;
//                         } // end of component loop i
//                     } // end of k1 loop
//                 } // end of k2 loop
//             } // end of face 1
//             if (m_faceOffset[2] > -1) {
// #pragma omp for nowait
//                 for (index_t k2=0; k2 < NE2; ++k2) {
//                     for (index_t k0=0; k0 < NE0; ++k0) {
//                         memcpy(&f_000[0], in.getSampleDataRO(INDEX3(k0,0,k2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_001[0], in.getSampleDataRO(INDEX3(k0,0,k2+1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_010[0], in.getSampleDataRO(INDEX3(k0,1,k2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_011[0], in.getSampleDataRO(INDEX3(k0,1,k2+1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_100[0], in.getSampleDataRO(INDEX3(k0+1,0,k2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_101[0], in.getSampleDataRO(INDEX3(k0+1,0,k2+1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_110[0], in.getSampleDataRO(INDEX3(k0+1,1,k2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_111[0], in.getSampleDataRO(INDEX3(k0+1,1,k2+1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         Scalar* o = out.getSampleDataRW(m_faceOffset[2]+INDEX2(k0,k2,NE0), zero);
//                         for (index_t i=0; i < numComp; ++i) {
//                             const Scalar V0=((f_100[i]-f_000[i])*C6 + (f_101[i]-f_001[i])*C2) / m_dx[0];
//                             const Scalar V1=((f_001[i]-f_000[i])*C6 + (f_101[i]-f_100[i])*C2) / m_dx[2];
//                             const Scalar V2=((f_001[i]-f_000[i])*C2 + (f_101[i]-f_100[i])*C6) / m_dx[2];
//                             o[INDEX3(i,0,0,numComp,3)] = V0;
//                             o[INDEX3(i,1,0,numComp,3)] = ((f_010[i]-f_000[i])*C5 + (f_111[i]-f_101[i])*C0 + (f_011[i]+f_110[i]-f_001[i]-f_100[i])*C1) / m_dx[1];
//                             o[INDEX3(i,2,0,numComp,3)] = V1;
//                             o[INDEX3(i,0,1,numComp,3)] = V0;
//                             o[INDEX3(i,1,1,numComp,3)] = ((f_110[i]-f_100[i])*C5 + (f_011[i]-f_001[i])*C0 + (f_010[i]+f_111[i]-f_000[i]-f_101[i])*C1) / m_dx[1];
//                             o[INDEX3(i,2,1,numComp,3)] = V2;
//                             o[INDEX3(i,0,2,numComp,3)] = V0;
//                             o[INDEX3(i,1,2,numComp,3)] = ((f_011[i]-f_001[i])*C5 + (f_110[i]-f_100[i])*C0 + (f_010[i]+f_111[i]-f_000[i]-f_101[i])*C1) / m_dx[1];
//                             o[INDEX3(i,2,2,numComp,3)] = V1;
//                             o[INDEX3(i,0,3,numComp,3)] = V0;
//                             o[INDEX3(i,1,3,numComp,3)] = ((f_111[i]-f_101[i])*C5 + (f_010[i]-f_000[i])*C0 + (f_011[i]+f_110[i]-f_001[i]-f_100[i])*C1) / m_dx[1];
//                             o[INDEX3(i,2,3,numComp,3)] = V2;
//                         } // end of component loop i
//                     } // end of k0 loop
//                 } // end of k2 loop
//             } // end of face 2
//             if (m_faceOffset[3] > -1) {
// #pragma omp for nowait
//                 for (index_t k2=0; k2 < NE2; ++k2) {
//                     for (index_t k0=0; k0 < NE0; ++k0) {
//                         memcpy(&f_000[0], in.getSampleDataRO(INDEX3(k0,m_NN[1]-2,k2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_001[0], in.getSampleDataRO(INDEX3(k0,m_NN[1]-2,k2+1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_010[0], in.getSampleDataRO(INDEX3(k0,m_NN[1]-1,k2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_011[0], in.getSampleDataRO(INDEX3(k0,m_NN[1]-1,k2+1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_100[0], in.getSampleDataRO(INDEX3(k0+1,m_NN[1]-2,k2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_101[0], in.getSampleDataRO(INDEX3(k0+1,m_NN[1]-2,k2+1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_110[0], in.getSampleDataRO(INDEX3(k0+1,m_NN[1]-1,k2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_111[0], in.getSampleDataRO(INDEX3(k0+1,m_NN[1]-1,k2+1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         Scalar* o = out.getSampleDataRW(m_faceOffset[3]+INDEX2(k0,k2,NE0), zero);
//                         for (index_t i=0; i < numComp; ++i) {
//                             const Scalar V0=((f_110[i]-f_010[i])*C6 + (f_111[i]-f_011[i])*C2) / m_dx[0];
//                             const Scalar V1=((f_110[i]-f_010[i])*C2 + (f_111[i]-f_011[i])*C6) / m_dx[0];
//                             const Scalar V2=((f_011[i]-f_010[i])*C6 + (f_111[i]-f_110[i])*C2) / m_dx[2];
//                             const Scalar V3=((f_011[i]-f_010[i])*C2 + (f_111[i]-f_110[i])*C6) / m_dx[2];
//                             o[INDEX3(i,0,0,numComp,3)] = V0;
//                             o[INDEX3(i,1,0,numComp,3)] = ((f_010[i]-f_000[i])*C5 + (f_111[i]-f_101[i])*C0 + (f_011[i]+f_110[i]-f_001[i]-f_100[i])*C1) / m_dx[1];
//                             o[INDEX3(i,2,0,numComp,3)] = V2;
//                             o[INDEX3(i,0,1,numComp,3)] = V0;
//                             o[INDEX3(i,1,1,numComp,3)] = ((f_110[i]-f_100[i])*C5 + (f_011[i]-f_001[i])*C0 + (f_010[i]+f_111[i]-f_000[i]-f_101[i])*C1) / m_dx[1];
//                             o[INDEX3(i,2,1,numComp,3)] = V3;
//                             o[INDEX3(i,0,2,numComp,3)] = V1;
//                             o[INDEX3(i,1,2,numComp,3)] = ((f_011[i]-f_001[i])*C5 + (f_110[i]-f_100[i])*C0 + (f_010[i]+f_111[i]-f_000[i]-f_101[i])*C1) / m_dx[1];
//                             o[INDEX3(i,2,2,numComp,3)] = V2;
//                             o[INDEX3(i,0,3,numComp,3)] = V1;
//                             o[INDEX3(i,1,3,numComp,3)] = ((f_111[i]-f_101[i])*C5 + (f_010[i]-f_000[i])*C0 + (f_011[i]+f_110[i]-f_001[i]-f_100[i])*C1) / m_dx[1];
//                             o[INDEX3(i,2,3,numComp,3)] = V3;
//                         } // end of component loop i
//                     } // end of k0 loop
//                 } // end of k2 loop
//             } // end of face 3
//             if (m_faceOffset[4] > -1) {
// #pragma omp for nowait
//                 for (index_t k1=0; k1 < NE1; ++k1) {
//                     for (index_t k0=0; k0 < NE0; ++k0) {
//                         memcpy(&f_000[0], in.getSampleDataRO(INDEX3(k0,k1,0, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_001[0], in.getSampleDataRO(INDEX3(k0,k1,1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_010[0], in.getSampleDataRO(INDEX3(k0,k1+1,0, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_011[0], in.getSampleDataRO(INDEX3(k0,k1+1,1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_100[0], in.getSampleDataRO(INDEX3(k0+1,k1,0, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_101[0], in.getSampleDataRO(INDEX3(k0+1,k1,1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_110[0], in.getSampleDataRO(INDEX3(k0+1,k1+1,0, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_111[0], in.getSampleDataRO(INDEX3(k0+1,k1+1,1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         Scalar* o = out.getSampleDataRW(m_faceOffset[4]+INDEX2(k0,k1,NE0), zero);
//                         for (index_t i=0; i < numComp; ++i) {
//                             const Scalar V0=((f_100[i]-f_000[i])*C6 + (f_110[i]-f_010[i])*C2) / m_dx[0];
//                             const Scalar V1=((f_100[i]-f_000[i])*C2 + (f_110[i]-f_010[i])*C6) / m_dx[0];
//                             const Scalar V2=((f_010[i]-f_000[i])*C6 + (f_110[i]-f_100[i])*C2) / m_dx[1];
//                             const Scalar V3=((f_010[i]-f_000[i])*C2 + (f_110[i]-f_100[i])*C6) / m_dx[1];
//                             o[INDEX3(i,0,0,numComp,3)] = V0;
//                             o[INDEX3(i,1,0,numComp,3)] = V2;
//                             o[INDEX3(i,2,0,numComp,3)] = ((f_001[i]-f_000[i])*C5 + (f_111[i]-f_110[i])*C0 + (f_011[i]+f_101[i]-f_010[i]-f_100[i])*C1) / m_dx[2];
//                             o[INDEX3(i,0,1,numComp,3)] = V0;
//                             o[INDEX3(i,1,1,numComp,3)] = V3;
//                             o[INDEX3(i,2,1,numComp,3)] = ((f_101[i]-f_100[i])*C5 + (f_011[i]-f_010[i])*C0 + (f_001[i]+f_111[i]-f_000[i]-f_110[i])*C1) / m_dx[2];
//                             o[INDEX3(i,0,2,numComp,3)] = V1;
//                             o[INDEX3(i,1,2,numComp,3)] = V2;
//                             o[INDEX3(i,2,2,numComp,3)] = ((f_011[i]-f_010[i])*C5 + (f_101[i]-f_100[i])*C0 + (f_001[i]+f_111[i]-f_000[i]-f_110[i])*C1) / m_dx[2];
//                             o[INDEX3(i,0,3,numComp,3)] = V1;
//                             o[INDEX3(i,1,3,numComp,3)] = V3;
//                             o[INDEX3(i,2,3,numComp,3)] = ((f_111[i]-f_110[i])*C5 + (f_001[i]-f_000[i])*C0 + (f_011[i]+f_101[i]-f_010[i]-f_100[i])*C1) / m_dx[2];
//                         } // end of component loop i
//                     } // end of k0 loop
//                 } // end of k1 loop
//             } // end of face 4
//             if (m_faceOffset[5] > -1) {
// #pragma omp for nowait
//                 for (index_t k1=0; k1 < NE1; ++k1) {
//                     for (index_t k0=0; k0 < NE0; ++k0) {
//                         memcpy(&f_000[0], in.getSampleDataRO(INDEX3(k0,k1,m_NN[2]-2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_001[0], in.getSampleDataRO(INDEX3(k0,k1,m_NN[2]-1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_010[0], in.getSampleDataRO(INDEX3(k0,k1+1,m_NN[2]-2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_011[0], in.getSampleDataRO(INDEX3(k0,k1+1,m_NN[2]-1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_100[0], in.getSampleDataRO(INDEX3(k0+1,k1,m_NN[2]-2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_101[0], in.getSampleDataRO(INDEX3(k0+1,k1,m_NN[2]-1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_110[0], in.getSampleDataRO(INDEX3(k0+1,k1+1,m_NN[2]-2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_111[0], in.getSampleDataRO(INDEX3(k0+1,k1+1,m_NN[2]-1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         Scalar* o = out.getSampleDataRW(m_faceOffset[5]+INDEX2(k0,k1,NE0), zero);
//                         for (index_t i=0; i < numComp; ++i) {
//                             const Scalar V0=((f_101[i]-f_001[i])*C6 + (f_111[i]-f_011[i])*C2) / m_dx[0];
//                             const Scalar V1=((f_101[i]-f_001[i])*C2 + (f_111[i]-f_011[i])*C6) / m_dx[0];
//                             const Scalar V2=((f_011[i]-f_001[i])*C6 + (f_111[i]-f_101[i])*C2) / m_dx[1];
//                             const Scalar V3=((f_011[i]-f_001[i])*C2 + (f_111[i]-f_101[i])*C6) / m_dx[1];
//                             o[INDEX3(i,0,0,numComp,3)] = V0;
//                             o[INDEX3(i,1,0,numComp,3)] = V2;
//                             o[INDEX3(i,2,0,numComp,3)] = ((f_001[i]-f_000[i])*C5 + (f_111[i]-f_110[i])*C0 + (f_011[i]+f_101[i]-f_010[i]-f_100[i])*C1) / m_dx[2];
//                             o[INDEX3(i,0,1,numComp,3)] = V0;
//                             o[INDEX3(i,1,1,numComp,3)] = V3;
//                             o[INDEX3(i,2,1,numComp,3)] = ((f_011[i]-f_010[i])*C0 + (f_101[i]-f_100[i])*C5 + (f_001[i]+f_111[i]-f_000[i]-f_110[i])*C1) / m_dx[2];
//                             o[INDEX3(i,0,2,numComp,3)] = V1;
//                             o[INDEX3(i,1,2,numComp,3)] = V2;
//                             o[INDEX3(i,2,2,numComp,3)] = ((f_011[i]-f_010[i])*C5 + (f_101[i]-f_100[i])*C0 + (f_001[i]+f_111[i]-f_000[i]-f_110[i])*C1) / m_dx[2];
//                             o[INDEX3(i,0,3,numComp,3)] = V1;
//                             o[INDEX3(i,1,3,numComp,3)] = V3;
//                             o[INDEX3(i,2,3,numComp,3)] = ((f_001[i]-f_000[i])*C0 + (f_111[i]-f_110[i])*C5 + (f_011[i]+f_101[i]-f_010[i]-f_100[i])*C1) / m_dx[2];
//                         } // end of component loop i
//                     } // end of k0 loop
//                 } // end of k1 loop
//             } // end of face 5
//         } // end of parallel section
//     } else if (out.getFunctionSpace().getTypeCode() == ReducedFaceElements) {
//         out.requireWrite();
// #pragma omp parallel
//         {
//             vector<Scalar> f_000(numComp, zero);
//             vector<Scalar> f_001(numComp, zero);
//             vector<Scalar> f_010(numComp, zero);
//             vector<Scalar> f_011(numComp, zero);
//             vector<Scalar> f_100(numComp, zero);
//             vector<Scalar> f_101(numComp, zero);
//             vector<Scalar> f_110(numComp, zero);
//             vector<Scalar> f_111(numComp, zero);
//             if (m_faceOffset[0] > -1) {
// #pragma omp for nowait
//                 for (index_t k2=0; k2 < NE2; ++k2) {
//                     for (index_t k1=0; k1 < NE1; ++k1) {
//                         memcpy(&f_000[0], in.getSampleDataRO(INDEX3(0,k1,k2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_001[0], in.getSampleDataRO(INDEX3(0,k1,k2+1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_010[0], in.getSampleDataRO(INDEX3(0,k1+1,k2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_011[0], in.getSampleDataRO(INDEX3(0,k1+1,k2+1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_100[0], in.getSampleDataRO(INDEX3(1,k1,k2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_101[0], in.getSampleDataRO(INDEX3(1,k1,k2+1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_110[0], in.getSampleDataRO(INDEX3(1,k1+1,k2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_111[0], in.getSampleDataRO(INDEX3(1,k1+1,k2+1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         Scalar* o = out.getSampleDataRW(m_faceOffset[0]+INDEX2(k1,k2,NE1), zero);
//                         for (index_t i=0; i < numComp; ++i) {
//                             o[INDEX3(i,0,0,numComp,3)] = (f_100[i]+f_101[i]+f_110[i]+f_111[i]-f_000[i]-f_001[i]-f_010[i]-f_011[i])*C3 / m_dx[0];
//                             o[INDEX3(i,1,0,numComp,3)] = (f_010[i]+f_011[i]-f_000[i]-f_001[i])*C4 / m_dx[1];
//                             o[INDEX3(i,2,0,numComp,3)] = (f_001[i]+f_011[i]-f_000[i]-f_010[i])*C4 / m_dx[2];
//                         } // end of component loop i
//                     } // end of k1 loop
//                 } // end of k2 loop
//             } // end of face 0
//             if (m_faceOffset[1] > -1) {
// #pragma omp for nowait
//                 for (index_t k2=0; k2 < NE2; ++k2) {
//                     for (index_t k1=0; k1 < NE1; ++k1) {
//                         memcpy(&f_000[0], in.getSampleDataRO(INDEX3(m_NN[0]-2,k1,k2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_001[0], in.getSampleDataRO(INDEX3(m_NN[0]-2,k1,k2+1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_010[0], in.getSampleDataRO(INDEX3(m_NN[0]-2,k1+1,k2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_011[0], in.getSampleDataRO(INDEX3(m_NN[0]-2,k1+1,k2+1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_100[0], in.getSampleDataRO(INDEX3(m_NN[0]-1,k1,k2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_101[0], in.getSampleDataRO(INDEX3(m_NN[0]-1,k1,k2+1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_110[0], in.getSampleDataRO(INDEX3(m_NN[0]-1,k1+1,k2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_111[0], in.getSampleDataRO(INDEX3(m_NN[0]-1,k1+1,k2+1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         Scalar* o = out.getSampleDataRW(m_faceOffset[1]+INDEX2(k1,k2,NE1), zero);
//                         for (index_t i=0; i < numComp; ++i) {
//                             o[INDEX3(i,0,0,numComp,3)] = (f_100[i]+f_101[i]+f_110[i]+f_111[i]-f_000[i]-f_001[i]-f_010[i]-f_011[i])*C3 / m_dx[0];
//                             o[INDEX3(i,1,0,numComp,3)] = (f_110[i]+f_111[i]-f_100[i]-f_101[i])*C4 / m_dx[1];
//                             o[INDEX3(i,2,0,numComp,3)] = (f_101[i]+f_111[i]-f_100[i]-f_110[i])*C4 / m_dx[2];
//                         } // end of component loop i
//                     } // end of k1 loop
//                 } // end of k2 loop
//             } // end of face 1
//             if (m_faceOffset[2] > -1) {
// #pragma omp for nowait
//                 for (index_t k2=0; k2 < NE2; ++k2) {
//                     for (index_t k0=0; k0 < NE0; ++k0) {
//                         memcpy(&f_000[0], in.getSampleDataRO(INDEX3(k0,0,k2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_001[0], in.getSampleDataRO(INDEX3(k0,0,k2+1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_010[0], in.getSampleDataRO(INDEX3(k0,1,k2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_011[0], in.getSampleDataRO(INDEX3(k0,1,k2+1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_100[0], in.getSampleDataRO(INDEX3(k0+1,0,k2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_101[0], in.getSampleDataRO(INDEX3(k0+1,0,k2+1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_110[0], in.getSampleDataRO(INDEX3(k0+1,1,k2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_111[0], in.getSampleDataRO(INDEX3(k0+1,1,k2+1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         Scalar* o = out.getSampleDataRW(m_faceOffset[2]+INDEX2(k0,k2,NE0), zero);
//                         for (index_t i=0; i < numComp; ++i) {
//                             o[INDEX3(i,0,0,numComp,3)] = (f_100[i]+f_101[i]-f_000[i]-f_001[i])*C4 / m_dx[0];
//                             o[INDEX3(i,1,0,numComp,3)] = (f_010[i]+f_011[i]+f_110[i]+f_111[i]-f_000[i]-f_001[i]-f_100[i]-f_101[i])*C3 / m_dx[1];
//                             o[INDEX3(i,2,0,numComp,3)] = (f_001[i]+f_101[i]-f_000[i]-f_100[i])*C4 / m_dx[2];
//                         } // end of component loop i
//                     } // end of k0 loop
//                 } // end of k2 loop
//             } // end of face 2
//             if (m_faceOffset[3] > -1) {
// #pragma omp for nowait
//                 for (index_t k2=0; k2 < NE2; ++k2) {
//                     for (index_t k0=0; k0 < NE0; ++k0) {
//                         memcpy(&f_000[0], in.getSampleDataRO(INDEX3(k0,m_NN[1]-2,k2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_001[0], in.getSampleDataRO(INDEX3(k0,m_NN[1]-2,k2+1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_010[0], in.getSampleDataRO(INDEX3(k0,m_NN[1]-1,k2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_011[0], in.getSampleDataRO(INDEX3(k0,m_NN[1]-1,k2+1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_100[0], in.getSampleDataRO(INDEX3(k0+1,m_NN[1]-2,k2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_101[0], in.getSampleDataRO(INDEX3(k0+1,m_NN[1]-2,k2+1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_110[0], in.getSampleDataRO(INDEX3(k0+1,m_NN[1]-1,k2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_111[0], in.getSampleDataRO(INDEX3(k0+1,m_NN[1]-1,k2+1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         Scalar* o = out.getSampleDataRW(m_faceOffset[3]+INDEX2(k0,k2,NE0), zero);
//                         for (index_t i=0; i < numComp; ++i) {
//                             o[INDEX3(i,0,0,numComp,3)] = (f_110[i]+f_111[i]-f_010[i]-f_011[i])*C4 / m_dx[0];
//                             o[INDEX3(i,1,0,numComp,3)] = (f_010[i]+f_011[i]+f_110[i]+f_111[i]-f_000[i]-f_001[i]-f_100[i]-f_101[i])*C3 / m_dx[1];
//                             o[INDEX3(i,2,0,numComp,3)] = (f_011[i]+f_111[i]-f_010[i]-f_110[i])*C4 / m_dx[2];
//                         } // end of component loop i
//                     } // end of k0 loop
//                 } // end of k2 loop
//             } // end of face 3
//             if (m_faceOffset[4] > -1) {
// #pragma omp for nowait
//                 for (index_t k1=0; k1 < NE1; ++k1) {
//                     for (index_t k0=0; k0 < NE0; ++k0) {
//                         memcpy(&f_000[0], in.getSampleDataRO(INDEX3(k0,k1,0, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_001[0], in.getSampleDataRO(INDEX3(k0,k1,1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_010[0], in.getSampleDataRO(INDEX3(k0,k1+1,0, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_011[0], in.getSampleDataRO(INDEX3(k0,k1+1,1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_100[0], in.getSampleDataRO(INDEX3(k0+1,k1,0, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_101[0], in.getSampleDataRO(INDEX3(k0+1,k1,1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_110[0], in.getSampleDataRO(INDEX3(k0+1,k1+1,0, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_111[0], in.getSampleDataRO(INDEX3(k0+1,k1+1,1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         Scalar* o = out.getSampleDataRW(m_faceOffset[4]+INDEX2(k0,k1,NE0), zero);
//                         for (index_t i=0; i < numComp; ++i) {
//                             o[INDEX3(i,0,0,numComp,3)] = (f_100[i]+f_110[i]-f_000[i]-f_010[i])*C4 / m_dx[0];
//                             o[INDEX3(i,1,0,numComp,3)] = (f_010[i]+f_110[i]-f_000[i]-f_100[i])*C4 / m_dx[1];
//                             o[INDEX3(i,2,0,numComp,3)] = (f_001[i]+f_011[i]+f_101[i]+f_111[i]-f_000[i]-f_010[i]-f_100[i]-f_110[i])*C4 / m_dx[2];
//                         } // end of component loop i
//                     } // end of k0 loop
//                 } // end of k1 loop
//             } // end of face 4
//             if (m_faceOffset[5] > -1) {
// #pragma omp for nowait
//                 for (index_t k1=0; k1 < NE1; ++k1) {
//                     for (index_t k0=0; k0 < NE0; ++k0) {
//                         memcpy(&f_000[0], in.getSampleDataRO(INDEX3(k0,k1,m_NN[2]-2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_001[0], in.getSampleDataRO(INDEX3(k0,k1,m_NN[2]-1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_010[0], in.getSampleDataRO(INDEX3(k0,k1+1,m_NN[2]-2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_011[0], in.getSampleDataRO(INDEX3(k0,k1+1,m_NN[2]-1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_100[0], in.getSampleDataRO(INDEX3(k0+1,k1,m_NN[2]-2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_101[0], in.getSampleDataRO(INDEX3(k0+1,k1,m_NN[2]-1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_110[0], in.getSampleDataRO(INDEX3(k0+1,k1+1,m_NN[2]-2, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         memcpy(&f_111[0], in.getSampleDataRO(INDEX3(k0+1,k1+1,m_NN[2]-1, m_NN[0],m_NN[1]), zero), numComp*sizeof(Scalar));
//                         Scalar* o = out.getSampleDataRW(m_faceOffset[5]+INDEX2(k0,k1,NE0), zero);
//                         for (index_t i=0; i < numComp; ++i) {
//                             o[INDEX3(i,0,0,numComp,3)] = (f_101[i]+f_111[i]-f_001[i]-f_011[i])*C4 / m_dx[0];
//                             o[INDEX3(i,1,0,numComp,3)] = (f_011[i]+f_111[i]-f_001[i]-f_101[i])*C4 / m_dx[1];
//                             o[INDEX3(i,2,0,numComp,3)] = (f_001[i]+f_011[i]+f_101[i]+f_111[i]-f_000[i]-f_010[i]-f_100[i]-f_110[i])*C3 / m_dx[2];
//                         } // end of component loop i
//                     } // end of k0 loop
//                 } // end of k1 loop
//             } // end of face 5
//         } // end of parallel section
//     }
}

//protected
void Brick::nodesToDOF(escript::Data& out, const escript::Data& in) const
{
    //TODO

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
    return -1;
    // return m_dofMap[node];
}

// instantiate our two supported versions
template
void Brick::assembleGradientImpl<real_t>(escript::Data& out,
                                         const escript::Data& in) const;

template
void Brick::assembleGradientImpl<cplx_t>(escript::Data& out,
                                         const escript::Data& in) const;

dim_t Brick::findNode(const double *coords) const
{
    //TODO
    return -1;

    
}

std::vector<IndexVector> Brick::getConnections(bool includeShared) const
{
    //TODO

//     // returns a vector v of size numDOF where v[i] is a vector with indices
//     // of DOFs connected to i (up to 27 in 3D).
//     // In other words this method returns the occupied (local) matrix columns
//     // for all (local) matrix rows.
//     // If includeShared==true then connections to non-owned DOFs are also
//     // returned (i.e. indices of the column couplings)
//     const dim_t nDOF0 = getNumDOFInAxis(0);
//     const dim_t nDOF1 = getNumDOFInAxis(1);
//     const dim_t nDOF2 = getNumDOFInAxis(2);
//     const dim_t numMatrixRows = nDOF0*nDOF1*nDOF2;
//     vector<IndexVector> indices(numMatrixRows);

//     if (includeShared) {
//         const index_t left = getFirstInDim(0);
//         const index_t bottom = getFirstInDim(1);
//         const index_t front = getFirstInDim(2);
//         const dim_t NN0 = m_NN[0];
//         const dim_t NN1 = m_NN[1];
//         const dim_t NN2 = m_NN[2];
// #pragma omp parallel for
//         for (index_t i=0; i < numMatrixRows; i++) {
//             const index_t x = left + i % nDOF0;
//             const index_t y = bottom + i % (nDOF0*nDOF1)/nDOF0;
//             const index_t z = front + i / (nDOF0*nDOF1);
//             // loop through potential neighbours and add to index if positions
//             // are within bounds
//             for (int i2=z-1; i2<z+2; i2++) {
//                 for (int i1=y-1; i1<y+2; i1++) {
//                     for (int i0=x-1; i0<x+2; i0++) {
//                         if (i0>=0 && i1>=0 && i2>=0
//                                 && i0<NN0 && i1<NN1 && i2<NN2) {
//                             indices[i].push_back(m_dofMap[i2*NN0*NN1+i1*NN0+i0]);
//                         }
//                     }
//                 }
//             }
//             sort(indices[i].begin(), indices[i].end());
//         }
//     } else {
// #pragma omp parallel for
//         for (index_t i=0; i < numMatrixRows; i++) {
//             const index_t x = i % nDOF0;
//             const index_t y = i % (nDOF0*nDOF1)/nDOF0;
//             const index_t z = i / (nDOF0*nDOF1);
//             // loop through potential neighbours and add to index if positions
//             // are within bounds
//             for (int i2=z-1; i2<z+2; i2++) {
//                 for (int i1=y-1; i1<y+2; i1++) {
//                     for (int i0=x-1; i0<x+2; i0++) {
//                         if (i0>=0 && i1>=0 && i2>=0
//                                 && i0<nDOF0 && i1<nDOF1 && i2<nDOF2) {
//                             indices[i].push_back(i2*nDOF0*nDOF1+i1*nDOF0+i0);
//                         }
//                     }
//                 }
//             }
//         }
//     }
//     return indices;
}

Assembler_ptr Brick::createAssembler(std::string type, const DataMap& constants) const
{
    // bool isComplex = false;
    // DataMap::const_iterator it;
    // for (it = constants.begin(); it != constants.end(); it++) {
    //     if (!it->second.isEmpty() && it->second.isComplex()) {
    //         isComplex = true;
    //         break;
    //     }
    // }

    // if (type.compare("DefaultAssembler") == 0) {
    //     if (isComplex) {
    //         return Assembler_ptr(new DefaultAssembler3D<cplx_t>(shared_from_this(), m_dx, m_NE, m_NN));
    //     } else {
    //         return Assembler_ptr(new DefaultAssembler3D<real_t>(shared_from_this(), m_dx, m_NE, m_NN));
    //     }
    // } else if (type.compare("WaveAssembler") == 0) {
    //     return Assembler_ptr(new WaveAssembler3D(shared_from_this(), m_dx, m_NE, m_NN, constants));
    // } else if (type.compare("LameAssembler") == 0) {
    //     return Assembler_ptr(new LameAssembler3D(shared_from_this(), m_dx, m_NE, m_NN));
    // } else { //else ifs would go before this for other types
    //     throw RipleyException("Oxley::Brick does not support the requested "
    //                           "assembler");
    // }
}

} // end of namespace oxley
