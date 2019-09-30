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

#include <oxley/InitAlgorithms.h>
#include <oxley/OtherAlgorithms.h>
#include <oxley/OxleyException.h>

#include <p4est.h>
#include <p8est.h>
#include <p4est_bits.h>
#include <p8est_bits.h>

namespace oxley {

void init_rectangle_data(p4est_t * p4est, p4est_topidx_t tree, p4est_quadrant_t * q)
{
    // the data associated with each quadrant
    quadrantData *data = (quadrantData *) q->p.user_data;

    data->u=0.0;
    data->quadTag=0;

    // Save the spatial coordinates
    p4est_qcoord_to_vertex(p4est->connectivity, tree, q->x, q->y, &data->xy[0]);
}

void init_brick_data(p8est_t * p8est, p4est_topidx_t tree, p8est_quadrant_t * q)
{
    // the data associated with each quadrant
    octantData *data = (octantData *) q->p.user_data;

    data->u=0.0;
  	data->octantTag=0;

    // Save the spatial coordinates
    p8est_qcoord_to_vertex(p8est->connectivity, tree, q->x, q->y, q->z, &data->xyz[0]);
}

void gce_rectangle_replace(p4est_t *p4est, p4est_topidx_t tree,
    int num_outgoing, p4est_quadrant_t *outgoing[],
    int num_incoming, p4est_quadrant_t *incoming[])
{
    if(num_incoming == P4EST_CHILDREN && num_outgoing == 1)
    {
        // The tag numbers being used
        p4estData * forestData = (p4estData *) p4est->user_pointer;
        addSurfaceData * info = (addSurfaceData *) forestData->info;
        int newTag = info->newTag;
        int oldTag = info->oldTag;

        // This is written in such a way as to minimise the number of calls
        // to "aboveCurve" during refinement.

        // The first quadrant's corner is the same as the parents
        quadrantData *parentData = (quadrantData *) outgoing[0]->p.user_data;
        quadrantData *childData  = (quadrantData *) incoming[0]->p.user_data;
        if(parentData->nodeTag)
            childData->quadTag = oldTag;
        else
            childData->quadTag = newTag;
        childData->xy[0] = parentData->xy[0];
        childData->xy[1] = parentData->xy[1];

        // Loop over the remaining children
        for(int i = 1; i < P4EST_CHILDREN; i++)
        {
            quadrantData *childData = (quadrantData *) incoming[i]->p.user_data;
            childData->xy[0] = incoming[i]->x;
            childData->xy[1] = incoming[i]->y;

            int n = info->y.size();
            childData->nodeTag=aboveCurve(info->x,info->y,
                                n,childData->xy[0],childData->xy[1]);

            if(childData->nodeTag)
                childData->quadTag = oldTag;
            else
                childData->quadTag = newTag;
        }
    }
    else if(num_incoming == 1 && num_outgoing == P4EST_CHILDREN)
    {
        // Coarsening should not occur during the gce algorithm
        throw OxleyException("gce_rectangle_replace: Unexpected attempt to coarsen the mesh.");
    }
    else
    {
        throw OxleyException("gce_rectangle_replace: Unknown error.");
    }
}

void gce_brick_replace(p8est_t *p8est, p4est_topidx_t tree,
    int num_outgoing, p8est_quadrant_t *outgoing[],
    int num_incoming, p8est_quadrant_t *incoming[])
{
    if(num_incoming == P8EST_CHILDREN && num_outgoing == 1)
    {
        // The tag numbers being used
        p8estData * forestData = (p8estData *) p8est->user_pointer;
        addSurfaceData * surfaceinfo = (addSurfaceData *) forestData->info;
        int newTag = surfaceinfo->newTag;
        int oldTag = surfaceinfo->oldTag;

        // This is written in such a way as to minimise the number of calls
        // to "aboveCurve" during refinement.

        // The first quadrant's corner is the same as the parents
        octantData *parentData = (octantData *) outgoing[0]->p.user_data;
        octantData *childData  = (octantData *) incoming[0]->p.user_data;
        if(parentData->nodeTag)
            childData->octantTag = oldTag;
        else
            childData->octantTag = newTag;
        childData->xyz[0] = parentData->xyz[0];
        childData->xyz[1] = parentData->xyz[1];
        childData->xyz[2] = parentData->xyz[2];

        // Loop over the remaining children
        for(int i = 1; i < P8EST_CHILDREN; i++)
        {
            octantData *childData = (octantData *) incoming[i]->p.user_data;
            childData->xyz[0] = incoming[i]->x;
            childData->xyz[1] = incoming[i]->y;
            childData->xyz[2] = incoming[i]->z;

            long nx = surfaceinfo->x.size();
            long ny = surfaceinfo->y.size();
            childData->nodeTag=aboveSurface(surfaceinfo->x,surfaceinfo->y,surfaceinfo->z,
                nx,ny,childData->xyz[0],childData->xyz[1],childData->xyz[2]);

            if(childData->nodeTag)
                childData->octantTag = oldTag;
            else
                childData->octantTag = newTag;
        }
    }
    else if(num_incoming == 1 && num_outgoing == P8EST_CHILDREN)
    {
        // Coarsening should not occur during the gce algorithm
        throw OxleyException("gce_brick_replace: Unexpected attempt to coarsen the mesh.");
    }
    else
    {
        throw OxleyException("gce_brick_replace: Unknown error.");
    }
}

void refine_copy_parent_quadrant(p4est_t * p4est, p4est_topidx_t tree,
    int num_outgoing, p4est_quadrant_t * outgoing[],
    int num_incoming, p4est_quadrant_t * incoming[])
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
    int num_outgoing, p8est_quadrant_t * outgoing[],
    int num_incoming, p8est_quadrant_t * incoming[])
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
