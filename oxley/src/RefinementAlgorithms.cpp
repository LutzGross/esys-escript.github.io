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

#include <iostream>

#include <oxley/RefinementAlgorithms.h>
#include <oxley/OtherAlgorithms.h>
#include <oxley/OxleyData.h>
#include <oxley/OxleyException.h>

#include <p4est_bits.h>
#include <p8est_bits.h>

// This file contains various callback functions to decide on refinement.
 namespace oxley {

int refine_uniform(p4est_t * p4est, p4est_topidx_t tree, p4est_quadrant_t * quadrant)
{
    return 1;
}

int refine_uniform(p8est_t * p4est, p4est_topidx_t tree, p8est_quadrant_t * quadrant)
{
    return 1;
}

void gce_first_pass(p4est_iter_volume_info_t * info, void *quad_data)
{
    // Get some pointers
    p4est_t * p4est = info->p4est;
    p4estData * forestData = (p4estData *) p4est->user_pointer;
    addSurfaceData * surfaceinfo = (addSurfaceData *) forestData->info;
    p4est_quadrant_t * quad = info->quad;
    quadrantData *quaddata = (quadrantData *) quad->p.user_data;
    p4est_topidx_t tree = info->treeid;

    // Get the coordinates of the four corners
    double xy[4][2];
    p4est_qcoord_to_vertex(p4est->connectivity, tree, quad->x, quad->y, xy[0]);

    // Work out if we are above or below the surface
    long n = surfaceinfo->x.size();
    quaddata->nodeTag = aboveCurve(surfaceinfo->x, surfaceinfo->y,
                                n, xy[0][0], xy[0][1]);

}

void gce_second_pass(p4est_iter_volume_info_t * info, void *quad_data)
{
    // Get some pointers
    p4est_t * p4est = info->p4est;
    p4estData * forestData = (p4estData *) p4est->user_pointer;
    p4est_quadrant_t * quad = info->quad;
    quadrantData *quaddata = (quadrantData *) quad->p.user_data;
    // p4est_topidx_t tree = info->treeid;

    // This variable records whether a node is above or below the curve
    bool ab[4];
    ab[0] = quaddata->nodeTag;

    // Note: For the numbering scheme used by p4est, cf. Burstedde et al. (2011)
    // First face neighbour
    p4est_quadrant_t * neighbourQuadrant = quad;
    quadrantData *neighbourData = (quadrantData *) quad->p.user_data;
    p4est_quadrant_face_neighbor(quad, 3, neighbourQuadrant);
    ab[2] = neighbourData->nodeTag;

    // The corner neighbour
    p4est_quadrant_corner_neighbor(quad, 2, neighbourQuadrant);
    ab[3] = neighbourData->nodeTag;

    // Second face neighbour
    p4est_quadrant_face_neighbor(quad, 1, neighbourQuadrant);
    ab[1] = neighbourData->nodeTag;

    // If at least two of the nodes are on different sides of the curve,
    // record the octant tag for the next step
    bool allNodesAreTheSame = ((ab[0] == ab[2]) && (ab[2] == ab[3])
                            && (ab[3] == ab[1]) && (ab[1] == ab[0]));

    addSurfaceData * surfaceinfo = (addSurfaceData *) forestData->info;
    if(surfaceinfo->oldTag == -1 && !allNodesAreTheSame)
        surfaceinfo->oldTag = quaddata->quadTag;

    // If the quadrant is not going to be refined, and if it is above the curve
    // update the tag information
    if(allNodesAreTheSame && quaddata->nodeTag)
        quaddata->quadTag=surfaceinfo->newTag;
}

int refine_gce(p4est_t * p4est, p4est_topidx_t tree, p4est_quadrant_t * quad)
{
    // Get some pointers
    p4estData * forestData = (p4estData *) p4est->user_pointer;
    quadrantData *quaddata = (quadrantData *) quad->p.user_data;

    // Get the tag info
    addSurfaceData * surfaceinfo = (addSurfaceData *) forestData->info;
    int oldTag = surfaceinfo->oldTag;

    // This variable records whether a node is above or below the curve
    bool ab[4];
    ab[0] = quaddata->nodeTag;

    // First face neighbour
    p4est_quadrant_t * neighbourQuadrant = quad;
    quadrantData *neighbourData = (quadrantData *) quad->p.user_data;
    p4est_quadrant_face_neighbor(quad, 3, neighbourQuadrant);
    ab[2] = neighbourData->nodeTag;

    // The corner neighbour
    p4est_quadrant_corner_neighbor(quad, 2, neighbourQuadrant);
    ab[3] = neighbourData->nodeTag;

    // Second face neighbour
    p4est_quadrant_face_neighbor(quad, 1, neighbourQuadrant);
    ab[1] = neighbourData->nodeTag;

    // If we are in the region being defined
    if(quaddata->quadTag == oldTag)
    {
        bool AllNodesAreTheSame = ((ab[0] == ab[2]) && (ab[2] == ab[3])
                                && (ab[3] == ab[1]) && (ab[1] == ab[0]));
        return !AllNodesAreTheSame;
    }
    else
    {
        return false;
    }
}

void gce_first_pass(p8est_iter_volume_info_t * info, void *quad_data)
{
    // Get some pointers
    p8est_t * p8est = info->p4est;
    p8estData * forestData = (p8estData *) p8est->user_pointer;
    addSurfaceData * surfaceinfo = (addSurfaceData *) forestData->info;
    p8est_quadrant_t * quad = info->quad;
    quadrantData *quaddata = (quadrantData *) quad->p.user_data;
    p4est_topidx_t tree = info->treeid;

    // Get the coordinates of the corner
    double xyz[4][3];
    p8est_qcoord_to_vertex(p8est->connectivity, tree, quad->x, quad->y, quad->z, xyz[0]);

    // Work out if we are above or below the surface
    long nx = surfaceinfo->x.size();
    long ny = surfaceinfo->y.size();
    quaddata->nodeTag = aboveSurface(surfaceinfo->x, surfaceinfo->y, surfaceinfo->z,
                            nx, ny, xyz[0][0], xyz[0][1], xyz[0][2]);
}

void gce_second_pass(p8est_iter_volume_info_t * info, void *quad_data)
{
    // Get some pointers
    p8est_t * p8est = info->p4est;
    p8estData * forestData = (p8estData *) p8est->user_pointer;
    addSurfaceData * surfaceinfo = (addSurfaceData *) forestData->info;
    p8est_quadrant_t * quad = info->quad;
    quadrantData *quaddata = (quadrantData *) quad->p.user_data;
    // p4est_topidx_t tree = info->treeid;

    // This variable records whether a node is above or below the curve
    bool ab[8];
    ab[0] = quaddata->nodeTag;

    // Note: For the numbering scheme used by p8est, cf. Burstedde et al. (2011)
    // The neighbouring octants are on edges 3, 7, 11; faces 1, 3, 5;
    // and corner 7.
    p8est_quadrant_t * neighbourQuadrant = quad;
    quadrantData *neighbourData = (quadrantData *) quad->p.user_data;

    p8est_quadrant_face_neighbor(quad, 1, neighbourQuadrant);
    ab[1] = neighbourData->nodeTag;
    p8est_quadrant_face_neighbor(quad, 3, neighbourQuadrant);
    ab[2] = neighbourData->nodeTag;
    p8est_quadrant_edge_neighbor(quad, 11, neighbourQuadrant);
    ab[3] = neighbourData->nodeTag;
    p8est_quadrant_face_neighbor(quad, 5, neighbourQuadrant);
    ab[4] = neighbourData->nodeTag;
    p8est_quadrant_edge_neighbor(quad, 7, neighbourQuadrant);
    ab[5] = neighbourData->nodeTag;
    p8est_quadrant_edge_neighbor(quad, 3, neighbourQuadrant);
    ab[6] = neighbourData->nodeTag;
    p8est_quadrant_corner_neighbor(quad, 7, neighbourQuadrant);
    ab[7] = neighbourData->nodeTag;

    // If at least two of the nodes are on different sides of the curve,
    // record the octant tag for the next step
    bool allNodesAreTheSame = ((ab[0] == ab[2]) && (ab[2] == ab[3])
                            && (ab[3] == ab[4]) && (ab[4] == ab[5])
                            && (ab[5] == ab[6]) && (ab[6] == ab[7])
                            && (ab[7] == ab[0]));
    if(surfaceinfo->oldTag == -1 && !allNodesAreTheSame)
        surfaceinfo->oldTag = quaddata->quadTag;

    // If the quadrant is not going to be refined, and if it is above the curve
    // update the tag information
    if(allNodesAreTheSame && quaddata->nodeTag)
        quaddata->quadTag=surfaceinfo->newTag;
}

int refine_gce(p8est_t * p8est, p4est_topidx_t tree, p8est_quadrant_t * quad)
{
    // Get some pointers
    p8estData * forestData = (p8estData *) p8est->user_pointer;
    quadrantData *quaddata = (quadrantData *) quad->p.user_data;
    addSurfaceData * surfaceinfo = (addSurfaceData *) forestData->info;

    // Get the tag info
    int oldTag = surfaceinfo->oldTag;

    // This variable records whether a node is above or below the curve
    bool ab[8];
    ab[0] = quaddata->nodeTag;

    // Get the other coordinates
    p8est_quadrant_t * neighbourQuadrant = quad;
    quadrantData *neighbourData = (quadrantData *) quad->p.user_data;
    p8est_quadrant_face_neighbor(quad, 1, neighbourQuadrant);
    ab[1] = neighbourData->nodeTag;
    p8est_quadrant_face_neighbor(quad, 3, neighbourQuadrant);
    ab[2] = neighbourData->nodeTag;
    p8est_quadrant_edge_neighbor(quad, 11, neighbourQuadrant);
    ab[3] = neighbourData->nodeTag;
    p8est_quadrant_face_neighbor(quad, 5, neighbourQuadrant);
    ab[4] = neighbourData->nodeTag;
    p8est_quadrant_edge_neighbor(quad, 7, neighbourQuadrant);
    ab[5] = neighbourData->nodeTag;
    p8est_quadrant_edge_neighbor(quad, 3, neighbourQuadrant);
    ab[6] = neighbourData->nodeTag;
    p8est_quadrant_corner_neighbor(quad, 7, neighbourQuadrant);
    ab[7] = neighbourData->nodeTag;

    // If we are in the region being defined
    if(quaddata->quadTag == oldTag)
    {
        bool allNodesAreTheSame = ((ab[0] == ab[2]) && (ab[2] == ab[3])
                                && (ab[3] == ab[4]) && (ab[4] == ab[5])
                                && (ab[5] == ab[6]) && (ab[6] == ab[7])
                                && (ab[7] == ab[0]));
        return !allNodesAreTheSame;
    }
    else
    {
        return false;
    }
}


} // namespace oxley
