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

// AEAE
void print_quad_debug_info(p4est_iter_volume_info_t * info, p4est_quadrant_t * quadrant)
{
    double xy[2] = {0.0,0.0};
    if(!p4est_quadrant_is_valid(quadrant))
        std::cout << "WARNING! Invalid quadrant: " << info->treeid << "." << info->quadid << std::endl;
    p4est_qcoord_to_vertex(info->p4est->connectivity, info->treeid, quadrant->x, quadrant->y, &xy[0]);
    std::cout << info->treeid << "." << info->quadid <<
        ": qcoords = " << quadrant->x << ", " << quadrant->y <<
        "\tquadlevel = " << quadrant->level << ", pad8=" << quadrant->pad8 << ", pad16=" << quadrant->pad16 <<
        "\t(x,y) = (" << xy[0] << "," << xy[1] << ")";

    quadrantData *quaddata = (quadrantData *) quadrant->p.user_data;
    std::cout << ",\tnodetag = " << quaddata->nodeTag
        << " locat: " << &quadrant
        << std::endl;
}

// AEAE
void print_quad_debug_info(p8est_iter_volume_info_t * info, p8est_quadrant_t * quadrant)
{
    double xyz[2] = {0.0,0.0};
    if(!p8est_quadrant_is_valid(quadrant))
        std::cout << "WARNING! Invalid quadrant" << info->treeid << "." << info->quadid << std::endl;
    p8est_qcoord_to_vertex(info->p4est->connectivity, info->treeid, quadrant->x, quadrant->y, quadrant->z, &xyz[0]);
    std::cout << info->treeid << "." << info->quadid << ": (x,y,z) = (" << xyz[0] << "," << xyz[1] << "," <<xyz[2] << ")" << std::endl;
}

void gce_first_pass(p4est_iter_volume_info_t * info, void *tmp)
{
    // Get some pointers
    // p4est_t * p4est = info->p4est;
    // p4estData * forestData = (p4estData *) info->p4est->user_pointer;
    // addSurfaceData * surfaceinfo = (addSurfaceData *) tmp;

    // p4est_quadrant_t * quad = info->quad;
    // quadrantData *quaddata = (quadrantData *) info->quad->p.user_data;
    // p4est_topidx_t tree = info->treeid;

    // Get the spatial coordinates of the corner node
    double xy[2];
    p4est_qcoord_to_vertex(info->p4est->connectivity, info->treeid, info->quad->x, info->quad->y, xy);

    // Work out if the node is above or below the curve and update nodeTag
    addSurfaceData * surfaceinfo = (addSurfaceData *) tmp;
    long n = surfaceinfo->x.size();
    quadrantData *quaddata = (quadrantData *) info->quad->p.user_data;
    quaddata->nodeTag = aboveCurve(surfaceinfo->x, surfaceinfo->y, n, xy[0], xy[1]);

#ifdef P4EST_ENABLE_DEBUG
    print_quad_debug_info(info, info->quad);
#endif

    if(!p4est_quadrant_is_valid(info->quad))
        throw OxleyException("Invalid quadrant");
}


void gce_second_pass(p4est_iter_volume_info_t * info, void *tmp)
{
    // Note: For the numbering scheme used by p4est, cf. Burstedde et al. (2011)

#ifdef P4EST_ENABLE_DEBUG
    print_quad_debug_info(info, info->quad);
#endif

    // Get some pointers
    quadrantData * quaddata = (quadrantData *) info->quad->p.user_data;

    // a quadrant to temporarily store data
    p4est_quadrant_t * tmpquad;
    quadrantData * tmpquaddata;
    tmpquaddata = (quadrantData *) tmpquad->p.user_data;

    // This variable records whether a node is above, on (true) or below (false) the curve
    bool ab[4] = {false};

    // Work out which of the nodes are above and below the curve
    ab[0] = quaddata->nodeTag;

    // First face neighbour
#ifdef P4EST_ENABLE_DEBUG
    print_quad_debug_info(info, tmpquad);
#endif
    p4est_quadrant_face_neighbor(info->quad, 1, tmpquad);
    ab[1] = tmpquaddata->nodeTag;

    // Second face neighbour
#ifdef P4EST_ENABLE_DEBUG
    print_quad_debug_info(info, tmpquad);
#endif
    p4est_quadrant_face_neighbor(info->quad, 3, tmpquad);
    ab[2] = tmpquaddata->nodeTag;

    // The corner neighbour
#ifdef P4EST_ENABLE_DEBUG
    print_quad_debug_info(info, tmpquad);
#endif
    p4est_quadrant_corner_neighbor(info->quad, 3, tmpquad);
    ab[3] = tmpquaddata->nodeTag;

    // If at least two of the nodes are on different sides of the curve,
    // record the octant tag for the next step
    bool allNodesAreTheSame = ((ab[0] == ab[2]) && (ab[2] == ab[3])
                            && (ab[3] == ab[1]) && (ab[1] == ab[0]));
    addSurfaceData * surfaceinfo = (addSurfaceData *) tmp;
    if(surfaceinfo->oldTag == -1 && !allNodesAreTheSame)
        surfaceinfo->oldTag = quaddata->quadTag;

    // If the quadrant is not going to be refined, and if it is above the curve
    // update the tag information
    if(allNodesAreTheSame && quaddata->nodeTag)
        quaddata->quadTag = surfaceinfo->newTag;
}

int refine_gce(p4est_t * p4est, p4est_topidx_t tree, p4est_quadrant_t * quad)
{
    // Get the tag info
    p4estData * forestData = (p4estData *) p4est->user_pointer;
    addSurfaceData * surfaceinfo = (addSurfaceData *) forestData->info;
    int oldTag = surfaceinfo->oldTag;

    // This variable records whether a node is above or below the curve
    bool ab[4];
    quadrantData *quaddata = (quadrantData *) quad->p.user_data;
    ab[0] = quaddata->nodeTag;

    // First face neighbour
    p4est_quadrant_t * tmpquad = quad;
    quadrantData *neighbourData = (quadrantData *) tmpquad->p.user_data;
    p4est_quadrant_face_neighbor(quad, 1, tmpquad);
    ab[1] = neighbourData->nodeTag;

    // Second face neighbour
    p4est_quadrant_face_neighbor(quad, 3, tmpquad);
    ab[2] = neighbourData->nodeTag;

    // The corner neighbour
    p4est_quadrant_corner_neighbor(quad, 3, tmpquad);
    ab[3] = neighbourData->nodeTag;

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

void refine_copy_parent_quadrant(p4est_t * p4est, p4est_topidx_t tree,
                                 int num_outgoing,
                                 p4est_quadrant_t * outgoing[],
                                 int num_incoming,
                                 p4est_quadrant_t * incoming[])
{
    if(num_outgoing == 1 && num_incoming == 4)
    {
        // Get the parent quadrant
        p4est_quadrant_t parent;
        p4est_quadrant_parent(incoming[0], &parent);

        // parent user data
        quadrantData *parentData = (quadrantData *) parent.p.user_data;

        for(int i = 0; i < P4EST_CHILDREN; i++){
            quadrantData *childData = (quadrantData *) incoming[i]->p.user_data;
            childData->u=parentData->u;
          	childData->quadTag=parentData->quadTag;

            // Update the spatial coordinates
            p4est_qcoord_to_vertex(p4est->connectivity, tree,
                incoming[i]->x, incoming[i]->y, &childData->xy[0]);
        }
    }
    else if(num_outgoing == 4 && num_incoming == 1)
    {
        // Get the parent quadrant
        p4est_quadrant_t * child[4];
        p4est_quadrant_childrenpv(outgoing[0], child);

        // parent user data
        quadrantData *childData = (quadrantData *) child[0]->p.user_data;
        quadrantData *parentData = (quadrantData *) outgoing[0]->p.user_data;

        parentData->u=childData->u;
        parentData->quadTag=childData->quadTag;

        // Update the spatial coordinates
        p4est_qcoord_to_vertex(p4est->connectivity, tree,
            outgoing[0]->x, outgoing[0]->y, &parentData->xy[0]);
    }
    else
    {
        throw OxleyException("refine_copy_parent_quadrant: Unknown error.");
    }
}

void refine_copy_parent_octant(p8est_t * p8est, p4est_topidx_t tree,
                                 int num_outgoing,
                                 p8est_quadrant_t * outgoing[],
                                 int num_incoming,
                                 p8est_quadrant_t * incoming[])
{
    if(num_outgoing == 1 && num_incoming == 8)
    {
        // Get the parent quadrant
        p8est_quadrant_t parent;
        p8est_quadrant_parent(incoming[0], &parent);

        // parent user data
        octantData *parentData = (octantData *) parent.p.user_data;

        for(int i = 0; i < P8EST_CHILDREN; i++)
        {
            octantData *childData = (octantData *) incoming[i]->p.user_data;
            childData->u=parentData->u;
          	childData->octantTag=parentData->octantTag;

            // Update the spatial coordinates
            p8est_qcoord_to_vertex(p8est->connectivity, tree,
                incoming[i]->x, incoming[i]->y, incoming[i]->z, &childData->xyz[0]);
        }
    }
    else if(num_outgoing == 8 && num_incoming == 1)
    {
        // Get the parent quadrant
        p8est_quadrant_t * child[8];
        p8est_quadrant_childrenpv(outgoing[0], child);

        // parent user data
        octantData *childData = (octantData *) child[0]->p.user_data;
        octantData *parentData = (octantData *) incoming[0]->p.user_data;

        parentData->u=childData->u;
        parentData->octantTag=childData->octantTag;

        // Update the spatial coordinates
        p8est_qcoord_to_vertex(p8est->connectivity, tree,
            outgoing[0]->x, outgoing[0]->y, outgoing[0]->z, &parentData->xyz[0]);
    }
    else
    {
        throw OxleyException("refine_copy_parent_quadrant: Unknown error.");
    }
}

} // namespace oxley
