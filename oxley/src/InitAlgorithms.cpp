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
    data->quadTag=0;
    data->owner=p4est->mpirank;
    p4est_qcoord_to_vertex(p4est->connectivity, tree, q->x, q->y, &data->xy[0]);
    data->treeid=tree;
    data->needs_refinement=false;
    data->u_real = nullptr;
    data->u_cplx = nullptr;
}

void init_brick_data(p8est_t * p8est, p4est_topidx_t tree, p8est_quadrant_t * q)
{
    // the data associated with each quadrant
    octantData *data = (octantData *) q->p.user_data; //here
    data->u=0.0;
    data->nodeTag=0;
  	data->octantTag=0;
    p8est_qcoord_to_vertex(p8est->connectivity, tree, q->x, q->y, q->z, &data->xyz[0]);
    data->treeid = -1;
    data->owner = -1;
    data->m_faceOffset[6] = {false};
}

void gce_rectangle_replace(p4est_t *p4est, p4est_topidx_t tree,
    int num_outgoing, p4est_quadrant_t *outgoing[],
    int num_incoming, p4est_quadrant_t *incoming[])
{
    if(num_incoming == 4 && num_outgoing == 1)
    {
        // Get the tag numbers being used
        addSurfaceData * surfaceinfo = (addSurfaceData *)p4est->user_pointer;
        int newTag = surfaceinfo->newTag;
        int oldTag = surfaceinfo->oldTag;

        quadrantData * quaddata = (quadrantData *) outgoing[0]->p.user_data;
        if(aboveCurve(surfaceinfo->x, surfaceinfo->y, p4est->connectivity, tree,
            surfaceinfo->x.size(), incoming[0]->x, incoming[0]->y))
        {
            quaddata->nodeTag = newTag;
        }
        else
        {
            quaddata->nodeTag = oldTag;
        }

        // Update the coordinates
        p4est_qcoord_to_vertex(p4est->connectivity, tree, incoming[0]->x, incoming[0]->y, quaddata->xy);
    }
    else if(num_incoming == 1 && num_outgoing == 4)
    {
        // Get the tag numbers being used
        addSurfaceData * surfaceinfo = (addSurfaceData *)p4est->user_pointer;
        int newTag = surfaceinfo->newTag;
        int oldTag = surfaceinfo->oldTag;

        // Loop over the children
        for(int i = 1; i < 4; i++)
        {
            quadrantData * quaddata = (quadrantData *) incoming[i]->p.user_data;

            // Get coordinates
            p4est_qcoord_to_vertex(p4est->connectivity, tree, incoming[i]->x, incoming[i]->y, quaddata->xy);

            if(aboveCurve(surfaceinfo->x, surfaceinfo->y, p4est->connectivity, tree,
                surfaceinfo->x.size(), incoming[i]->x, incoming[i]->y))
            {
                quaddata->nodeTag = newTag;
            }
            else
            {
                quaddata->nodeTag = oldTag;
            }
        }
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
    if(num_incoming == 8 && num_outgoing == 1)
    {
        // Get the tag numbers being used
        addSurfaceData * surfaceinfo = (addSurfaceData *)p8est->user_pointer;
        int newTag = surfaceinfo->newTag;
        int oldTag = surfaceinfo->oldTag;

        octantData * quaddata = (octantData *) outgoing[0]->p.user_data;
        if(aboveSurface(surfaceinfo->x, surfaceinfo->y, surfaceinfo->z,
            p8est->connectivity, tree, surfaceinfo->x.size(), surfaceinfo->y.size(),
            incoming[0]->x, incoming[0]->y, incoming[0]->z))
        {
            quaddata->nodeTag = newTag;
        }
        else
        {
            quaddata->nodeTag = oldTag;
        }

        // Update the coordinates
        p8est_qcoord_to_vertex(p8est->connectivity, tree,
            incoming[0]->x, incoming[0]->y, incoming[0]->z, quaddata->xyz);
    }
    else if(num_incoming == 1 && num_outgoing == 8)
    {
        // Get the tag numbers being used
        addSurfaceData * surfaceinfo = (addSurfaceData *) p8est->user_pointer;
        int newTag = surfaceinfo->newTag;
        int oldTag = surfaceinfo->oldTag;

        // Get the length of the arrays for later
        long nx = surfaceinfo->x.size();
        long ny = surfaceinfo->y.size();

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
        for(int i = 1; i < 8; i++)
        {
            // Get coordinates
            double xyz[3];
            p8est_qcoord_to_vertex(p8est->connectivity, tree, incoming[i]->x, incoming[i]->y, incoming[i]->z, xyz);

            quadrantData * quaddata = (quadrantData *) incoming[i]->p.user_data;
            if(aboveSurface(surfaceinfo->x, surfaceinfo->y, surfaceinfo->z, p8est->connectivity,
                tree, nx, ny, incoming[i]->x, incoming[i]->y, incoming[i]->z))
            {
                quaddata->nodeTag = newTag;
            }
            else
            {
                quaddata->nodeTag = oldTag;
            }
        }
    }
    else
    {
        throw OxleyException("gce_brick_replace: Unknown error.");
    }
}

void gce_init_new_rectangle(p4est_t * p4est, p4est_topidx_t tree, p4est_quadrant_t * q)
{
    // the data associated with each quadrant
    quadrantData *data = (quadrantData *) q->p.user_data;

    data->quadTag=-1;

    // Save the spatial coordinates
    p4est_qcoord_to_vertex(p4est->connectivity, tree, q->x, q->y, &data->xy[0]);

    data->needs_refinement=false;
}

void gce_init_new_brick(p8est_t * p8est, p4est_topidx_t tree, p8est_quadrant_t * q)
{
    // the data associated with each quadrant
    octantData *data = (octantData *) q->p.user_data;

    data->u=0.0;
  	data->octantTag=-1;

    // Save the spatial coordinates
    p8est_qcoord_to_vertex(p8est->connectivity, tree, q->x, q->y, q->z, &data->xyz[0]);
}


void getConnections_shared(p4est_iter_face_info_t *info, void *indices)
{

}

void getConnections_notshared(p4est_iter_face_info_t *info, void *indices)
{
    // if(info->tree_boundary == True) // On the boundary of the forest, not the interior
    // {
    //     p4est_quadrant_t * quad = inf0->sides[0];

    //     int n0, n1;
        

    //     indices[n0].push_back(n1);
    //     indices[n1].push_back(n0);
    // }
}


} // namespace oxley
