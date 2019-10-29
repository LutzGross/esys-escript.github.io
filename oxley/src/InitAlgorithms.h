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

/*
 *  \brief
 *  Update octants such that each node is tagged with the right surface environment
 *  Used by the GCE algorithm
 */
void gce_brick_replace(p8est_t *p8est, p4est_topidx_t tree,
    int num_outgoing, p8est_quadrant_t *outgoing[],
    int num_incoming, p8est_quadrant_t *incoming[]);

///////////////////////////////////////////////////////////////////////////

} //namespace oxley
