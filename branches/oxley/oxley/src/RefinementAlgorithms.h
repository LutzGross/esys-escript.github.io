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

#include <oxley/OxleyData.h>

#include <p4est.h>
#include <p4est_iterate.h>
#include <p8est.h>
#include <p8est_iterate.h>

#ifdef ESYS_HAVE_BOOST_NUMPY
#include <boost/python/numpy.hpp>
#endif

namespace oxley {

// AEAE
void print_quad_debug_info(p4est_quadrant_t * quadrant);

// Uniform refinement
int refine_uniform(p4est_t * p4est, p4est_topidx_t tree, p4est_quadrant_t * quadrant);

int refine_uniform(p8est_t * p8est, p4est_topidx_t tree, p8est_quadrant_t * quadrant);

// Random refinement
int random_refine(p4est_t * p4est, p4est_topidx_t tree, p4est_quadrant_t * quadrant);

int random_refine(p8est_t * p8est, p4est_topidx_t tree, p8est_quadrant_t * quadrant);

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


// Copys parent info directly to the children without modification
void refine_copy_parent_quadrant(p4est_t * p4est, p4est_topidx_t tree,
                                 int num_outgoing,
                                 p4est_quadrant_t * outgoing[],
                                 int num_incoming,
                                 p4est_quadrant_t * incoming[]);

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

} //namespace oxley
