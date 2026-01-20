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
*
*****************************************************************************/

#include <oxley/OxleyData.h>

#include <p4est.h>
#include <p4est_iterate.h>
#include <p8est.h>
#include <p8est_iterate.h>

#ifdef ESYS_HAVE_BOOST_NUMPY
#include <boost/python/numpy.hpp>
#endif

namespace oxley {

// Debugging info
void print_quad_debug_info(p4est_quadrant_t * quadrant);

// Uniform refinement
int refine_uniform(p4est_t * p4est, p4est_topidx_t tree, p4est_quadrant_t * quadrant);
int refine_uniform(p8est_t * p8est, p4est_topidx_t tree, p8est_quadrant_t * quadrant);

// mare2dem
int refine_mare2dem(p4est_t * p4est, p4est_topidx_t tree, p4est_quadrant_t * quadrant);
int refine_mare2dem(p8est_t * p8est, p8est_topidx_t tree, p8est_quadrant_t * quadrant);

// Boundaries
int refine_north(p4est_t * p4est, p4est_topidx_t tree, p4est_quadrant_t * quadrant);
int refine_south(p4est_t * p4est, p4est_topidx_t tree, p4est_quadrant_t * quadrant);
int refine_west(p4est_t * p4est, p4est_topidx_t tree, p4est_quadrant_t * quadrant);
int refine_east(p4est_t * p4est, p4est_topidx_t tree, p4est_quadrant_t * quadrant);
int refine_north(p8est_t * p8est, p8est_topidx_t tree, p8est_quadrant_t * quadrant);
int refine_south(p8est_t * p8est, p8est_topidx_t tree, p8est_quadrant_t * quadrant);
int refine_west(p8est_t * p8est, p8est_topidx_t tree, p8est_quadrant_t * quadrant);
int refine_east(p8est_t * p8est, p8est_topidx_t tree, p8est_quadrant_t * quadrant);
int refine_top(p8est_t * p8est, p8est_topidx_t tree, p8est_quadrant_t * quadrant);
int refine_bottom(p8est_t * p8est, p8est_topidx_t tree, p8est_quadrant_t * quadrant);
int refine_region(p4est_t * p4est, p4est_topidx_t tree, p4est_quadrant_t * quadrant);
int refine_point(p4est_t * p4est, p4est_topidx_t tree, p4est_quadrant_t * quadrant);
int refine_circle(p4est_t * p4est, p4est_topidx_t tree, p4est_quadrant_t * quadrant);
int refine_region(p8est_t * p8est, p8est_topidx_t tree, p8est_quadrant_t * quadrant);
int refine_point(p8est_t * p8est, p8est_topidx_t tree, p8est_quadrant_t * quadrant);
int refine_sphere(p8est_t * p8est, p8est_topidx_t tree, p8est_quadrant_t * quadrant);
int refine_mask(p4est_t * p4est, p4est_topidx_t tree, p4est_quadrant_t * quadrant);
int refine_mask(p8est_t * p8est, p8est_topidx_t tree, p8est_quadrant_t * quadrant);

/*
 *  \brief
 * Checks that the quadrant is valid (used when debugging)
 */
void gce_first_pass(p4est_iter_volume_info_t * info, void *quad_data);

void gce_first_pass(p8est_iter_volume_info_t * info, void *quad_data);

/*
 *  \brief
 *  gce_second_pass identifies whether a quadrant/octant is in the region being
 *  refined and identifies whether the nodes are above or below the surface
 */
void gce_second_pass(p4est_iter_volume_info_t * info, void *quad_data);

void gce_second_pass(p8est_iter_volume_info_t * info, void *quad_data);

/*
 *  \brief
 *  gce_third_pass updates the tags of quadrants/octants inside the old region
 *  that were not refined in gce_second_pass
 */
void gce_third_pass(p4est_iter_volume_info_t * info, void *quad_data);

void gce_third_pass(p8est_iter_volume_info_t * info, void *quad_data);

/*
 *  \brief
 *  refine_gce returns true if we are refining this quadrant/octant
 */
int refine_gce(p4est_t * p4est, p4est_topidx_t tree, p4est_quadrant_t * quadrant);

int refine_gce(p8est_t * p8est, p4est_topidx_t tree, p8est_quadrant_t * quadrant);

/*
 *  \brief
 *  nodes interpolation
 */
int refine_nodesToNodesFiner(p4est_t * p4est, p4est_topidx_t tree, p4est_quadrant_t * quadrant);

// Copys parent info directly to the children without modification
// void refine_copy_parent_quadrant(p4est_t * p4est, p4est_topidx_t tree,
//                                  int num_outgoing,
//                                  p4est_quadrant_t * outgoing[],
//                                  int num_incoming,
//                                  p4est_quadrant_t * incoming[]);

void refine_copy_parent_octant(p8est_t * p8est, p4est_topidx_t tree,
                                 int num_outgoing,
                                 p8est_quadrant_t * outgoing[],
                                 int num_incoming,
                                 p8est_quadrant_t * incoming[]);

void get_interpolateNodesOnElementWorker_data(p4est_iter_volume_info_t * info, void *fxx);
void get_interpolateNodesOnElementWorker_data(p8est_iter_volume_info_t * info, void *fxx);

void get_interpolateNodesOnFacesWorker_data(p4est_iter_volume_info_t * info, void *fxx);
void get_interpolateNodesOnFacesWorker_data(p8est_iter_volume_info_t * info, void *fxx);

// Updates m_faceOffset for each quadrant / octant
void update_node_faceoffset(p4est_iter_volume_info_t * info, void *fxx);
void update_node_faceoffset(p8est_iter_volume_info_t * info, void *fxx);

// Used when updating myRows and myColumns
void update_RC(p4est_iter_face_info_t *info, void *user_data);
void update_RC(p8est_iter_edge_info_t *info, void *user_data);

// Used by getConnections
void update_connections(p4est_iter_volume_info_t *info, void *user_data);

// interpolation 
void refine_copy_parent_quadrant_data(p4est_t * p4est, p4est_topidx_t tree,
                                 int num_outgoing,
                                 p4est_quadrant_t * outgoing[],
                                 int num_incoming,
                                 p4est_quadrant_t * incoming[]);
void refine_copy_parent_element_data(p4est_t * p4est, p4est_topidx_t tree,
                                 int num_outgoing,
                                 p4est_quadrant_t * outgoing[],
                                 int num_incoming,
                                 p4est_quadrant_t * incoming[]);

} //namespace oxley
