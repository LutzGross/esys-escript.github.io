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
#include <p4est_extended.h>
#include <p8est.h>
#include <p8est_extended.h>

#ifdef ESYS_HAVE_BOOST_NUMPY
#include <boost/python/numpy.hpp>
#endif

// This file contains various callback functions to initialise the data

namespace oxley {


/*
 *  \brief
 *  Initial data for each Rectangle quadrant
 */
void init_rectangle_data(p4est_t * p4est, p4est_topidx_t which_tree, p4est_quadrant_t * q);

/*
 *  \brief
 *  Initial data for each Brick octant
 */
void init_brick_data(p8est_t * p8est, p4est_topidx_t which_tree, p8est_quadrant_t * q);

/*
 *  \brief
 *  Update quads such that each node is tagged with the right surface environment.
 *  Used by the GCE algorithm
 */
void gce_rectangle_replace(p4est_t *p4est, p4est_topidx_t tree,
    int num_outgoing, p4est_quadrant_t *outgoing[],
    int num_incoming, p4est_quadrant_t *incoming[]);

void gce_init_new_rectangle(p4est_t * p4est, p4est_topidx_t which_tree, p4est_quadrant_t * q);

/*
 *  \brief
 *  Update octants such that each node is tagged with the right surface environment
 *  Used by the GCE algorithm
 */
void gce_brick_replace(p8est_t *p8est, p4est_topidx_t tree,
    int num_outgoing, p8est_quadrant_t *outgoing[],
    int num_incoming, p8est_quadrant_t *incoming[]);

void gce_init_new_brick(p8est_t * p8est, p4est_topidx_t which_tree, p8est_quadrant_t * q);


/*
 *  \brief
 *  used to calculate the neighbouring nodes in the assembler
 */
void getConnections_shared(p4est_iter_face_info_t *info, void *user_data);
void getConnections_notshared(p4est_iter_face_info_t *info, void *user_data);


///////////////////////////////////////////////////////////////////////////

} //namespace oxley
