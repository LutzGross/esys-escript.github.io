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
#include <random>
#include <vector>

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

int random_refine(p4est_t * p4est, p4est_topidx_t tree, p4est_quadrant_t * quadrant)
{
    return ((int) rand() % 2) == 0;
}

int random_refine(p8est_t * p4est, p4est_topidx_t tree, p8est_quadrant_t * quadrant)
{
    return ((int) rand() % 2) == 0;
}

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

void print_quad_debug_info(p8est_iter_volume_info_t * info, p8est_quadrant_t * quadrant)
{
    double xyz[3] = {0.0,0.0,0.0};
    if(!p8est_quadrant_is_valid(quadrant))
        std::cout << "WARNING! Invalid quadrant" << info->treeid << "." << info->quadid << std::endl;
    p8est_qcoord_to_vertex(info->p4est->connectivity, info->treeid, quadrant->x, quadrant->y, quadrant->z, &xyz[0]);
    std::cout << info->treeid << "." << info->quadid << ": (x,y,z) = (" << xyz[0] << "," << xyz[1] << "," << xyz[2] << ")" << std::endl;
}

void gce_first_pass(p4est_iter_volume_info_t * info, void *tmp)
{
    if(!p4est_quadrant_is_valid(info->quad))
        throw OxleyException("Invalid quadrant");
}

void gce_second_pass(p4est_iter_volume_info_t * info, void *tmp)
{
    // Note: For the numbering scheme used by p4est, cf. Burstedde et al. (2011)

    // Get some pointers
    quadrantData * quaddata = (quadrantData *) info->quad->p.user_data;
    addSurfaceData * surfaceinfo = (addSurfaceData *) tmp;

    if(surfaceinfo->oldTag == -1)
    {
        // Spatial indices of x, y, z to pass on
        double x = info->quad->x;
        double y = info->quad->y;
        long n = surfaceinfo->x.size();

        // Work out the length of the quadrant
        double l = P4EST_QUADRANT_LEN(info->quad->level);

        // Check that the point is inside the domain defined by the function
        double xy[2][2];
        p4est_qcoord_to_vertex(info->p4est->connectivity, info->treeid, x, y, xy[0]);
        p4est_qcoord_to_vertex(info->p4est->connectivity, info->treeid, x+l, y, xy[1]);
        bool insideDomain = (xy[0][0] >= surfaceinfo->xmin) && (xy[1][0] <= surfaceinfo->xmax);

        if(insideDomain)
        {
            // This variable records whether a node is above, on (true) or below (false) the curve
            signed ab[4] = {0};
            double increment[4][2] = {{0,0},{l,0},{0,l},{l,l}};
#pragma omp parallel for
            for(int i = 0; i < 4; i++)
                ab[i] = aboveCurve(surfaceinfo->x, surfaceinfo->y,
                                    info->p4est->connectivity, info->treeid, n,
                                    x+increment[i][0], y+increment[i][1]);

            signed sum = 0;
            for(int i = 0; i < 4; i++)
                sum += ab[i];

            if(abs(sum) != 4)
                surfaceinfo->oldTag = quaddata->quadTag;
        }
    }
}

void gce_third_pass(p4est_iter_volume_info_t * info, void *tmp)
{
    // Get some pointers
    quadrantData * quaddata = (quadrantData *) info->quad->p.user_data;
    addSurfaceData * surfaceinfo = (addSurfaceData *) tmp;

    double x = info->quad->x;
    double y = info->quad->y;
    long n = surfaceinfo->x.size();

    // Check that the point is inside the domain defined by the function
    double l = P4EST_QUADRANT_LEN(info->quad->level);
    double xy[2][2];
    p4est_qcoord_to_vertex(info->p4est->connectivity, info->treeid, x, y, xy[0]);
    p4est_qcoord_to_vertex(info->p4est->connectivity, info->treeid, x+l, y, xy[1]);
    bool insideDomain = (xy[0][0] >= surfaceinfo->xmin) && (xy[1][0] <= surfaceinfo->xmax);

    bool needToUpdate = (quaddata->quadTag == surfaceinfo->oldTag) || (quaddata->quadTag == -1);

    if(insideDomain && needToUpdate)
    {
        bool above = aboveCurve(surfaceinfo->x, surfaceinfo->y,
            info->p4est->connectivity, info->treeid, n, x, y);
        if(above)
            quaddata->quadTag = surfaceinfo->newTag;
        else
            quaddata->quadTag = surfaceinfo->oldTag;
    }
}

int refine_gce(p4est_t * p4est, p4est_topidx_t tree, p4est_quadrant_t * quad)
{
    // Get some pointers
    addSurfaceData * surfaceinfo = (addSurfaceData *) p4est->user_pointer;
    quadrantData * quaddata = (quadrantData *) quad->p.user_data;

    // If we are in the region being defined
    if(quaddata->quadTag == surfaceinfo->oldTag || quaddata->quadTag == -1)
    {
        // The length of the coordinate vector
        long n = surfaceinfo->x.size();

        // Work out the length of the quadrant
        double l = P4EST_QUADRANT_LEN(quad->level);

        // Spatial indices of x, y, z to pass on
        double x = quad->x;
        double y = quad->y;

        // Check to see if we are inside the domain defined by z[x,y]
        double xy[2][2];
        p4est_qcoord_to_vertex(p4est->connectivity, tree, x, y, xy[0]);
        p4est_qcoord_to_vertex(p4est->connectivity, tree, x+l, y, xy[1]);
        bool insideDomain = (xy[0][0] >= surfaceinfo->xmin) && (xy[1][0] <= surfaceinfo->xmax);
        if(!insideDomain)
            return false;

        // This variable records whether a node is above, on (true) or below (false) the curve
        signed ab[4] = {0};
        double increment[4][2] = {{0,0},{l,0},{0,l},{l,l}};
#pragma omp parallel for
        for(int i = 0; i < 4; i++)
            ab[i] = aboveCurve(surfaceinfo->x, surfaceinfo->y, p4est->connectivity, tree, n, x+increment[i][0], y+increment[i][1]);

        signed sum = 0;
        for(int i = 0; i < 4; i++)
            sum += ab[i];

        return (abs(sum) == 4) ? false : true;
    }
    else
    {
        return false;
    }
}

// This is used when debugging to check to see if there are any invalid quads
void gce_first_pass(p8est_iter_volume_info_t * info, void *quad_data)
{
    if(!p8est_quadrant_is_valid(info->quad))
        throw OxleyException("Invalid quadrant");
}

// This works out the tag of the region we are about to subdivide
void gce_second_pass(p8est_iter_volume_info_t * info, void *tmp)
{
    // Note: For the numbering scheme used by p4est, cf. Burstedde et al. (2011)

    // Get some pointers
    octantData * quaddata = (octantData *) info->quad->p.user_data;
    addSurfaceData * surfaceinfo = (addSurfaceData *) tmp;

    if(surfaceinfo->oldTag == -1)
    {
        // Work out the length of the quadrant
        double l = P8EST_QUADRANT_LEN(info->quad->level);

        // Spatial indices of x, y, z to pass on
        double x = info->quad->x;
        double y = info->quad->y;
        double z = info->quad->z;

        long nx = surfaceinfo->x.size();
        long ny = surfaceinfo->y.size();

        // Check that the point is inside the domain defined by the function
        double xy[3][3];
        p8est_qcoord_to_vertex(info->p4est->connectivity, info->treeid, x, y, z, xy[0]);
        p8est_qcoord_to_vertex(info->p4est->connectivity, info->treeid, x+l, y, z, xy[1]);
        p8est_qcoord_to_vertex(info->p4est->connectivity, info->treeid, x, y+l, z, xy[2]);

        bool insideDomain = (xy[0][0] >= surfaceinfo->xmin) && (xy[1][0] <= surfaceinfo->xmax)
                         && (xy[1][1] >= surfaceinfo->ymin) && (xy[2][1] <= surfaceinfo->ymax);

        if(insideDomain)
        {
            // This variable records whether a node is above, on (true) or below (false) the curve
            signed ab[8] = {0};
            double increment[8][3] = {{0,0,0},{l,0,0},{0,l,0},{l,l,0},{0,0,l},{l,0,l},{0,l,l},{l,l,l}};
#pragma omp parallel for
            for(int i = 0; i < 8; i++)
                ab[i] = aboveSurface(surfaceinfo->x, surfaceinfo->y, surfaceinfo->z,
                    info->p4est->connectivity, info->treeid, nx, ny,
                    x+increment[i][0], y+increment[i][1], z+increment[i][2]);

            signed sum = 0;
            for(int i = 0; i < 8; i++)
                sum += ab[i];

            if(abs(sum) != 8)
                surfaceinfo->oldTag = quaddata->octantTag;
        }
    }
}

// This tags quads that are in the old region but are not being refined
void gce_third_pass(p8est_iter_volume_info_t * info, void *tmp)
{
    // Get some pointers
    octantData * quadData = (octantData *) info->quad->p.user_data;
    addSurfaceData * surfaceinfo = (addSurfaceData *) tmp;

    double x = info->quad->x;
    double y = info->quad->y;
    double z = info->quad->z;
    long nx = surfaceinfo->x.size();
    long ny = surfaceinfo->y.size();
    double l = P8EST_QUADRANT_LEN(info->quad->level);

    // Work out if we are inside the domain defined by z[x,y]
    double xy[3][3];
    p8est_qcoord_to_vertex(info->p4est->connectivity, info->treeid, x, y, z, xy[0]);
    p8est_qcoord_to_vertex(info->p4est->connectivity, info->treeid, x+l, y, z, xy[1]);
    p8est_qcoord_to_vertex(info->p4est->connectivity, info->treeid, x, y+l, z, xy[2]);
    bool insideDomain = (xy[0][0] >= surfaceinfo->xmin) && (xy[1][0] <= surfaceinfo->xmax)
                    && (xy[0][1] >= surfaceinfo->ymin) && (xy[2][1] <= surfaceinfo->ymax);

    bool needToUpdate = (quadData->octantTag == surfaceinfo->oldTag) || quadData->octantTag == -1;
    if(insideDomain && needToUpdate)
    {
        bool above = aboveSurface(surfaceinfo->x, surfaceinfo->y, surfaceinfo->z,
                        info->p4est->connectivity, info->treeid, nx, ny, x, y, z);
        if(above)
           quadData->octantTag = surfaceinfo->newTag;
        else
           quadData->octantTag = surfaceinfo->oldTag;
    }
}

// This decided whether or not to refine a quadrant
int refine_gce(p8est_t * p8est, p4est_topidx_t tree, p8est_quadrant_t * quad)
{
    // Get the tag info
    addSurfaceData * surfaceinfo = (addSurfaceData *) p8est->user_pointer;
    int oldTag = surfaceinfo->oldTag;
    octantData * quaddata = (octantData *) quad->p.user_data;

    if(quaddata->octantTag == oldTag || quaddata->octantTag == -1)
    {
        // Get the length of the quadrant
        double l = P8EST_QUADRANT_LEN(quad->level);

        // Spatial indices of x, y, z to pass on
        p4est_qcoord_t x = quad->x;
        p4est_qcoord_t y = quad->y;
        p4est_qcoord_t z = quad->z;

        long nx = surfaceinfo->x.size();
        long ny = surfaceinfo->y.size();

        // Check to see if we are inside the domain defined by z[x,y]
        double xy[3][3];
        p8est_qcoord_to_vertex(p8est->connectivity, tree, x, y, z, xy[0]);
        p8est_qcoord_to_vertex(p8est->connectivity, tree, x+l, y, z, xy[1]);
        p8est_qcoord_to_vertex(p8est->connectivity, tree, x, y+l, z, xy[2]);
        bool insideDomain = (xy[0][0] >= surfaceinfo->xmin) && (xy[1][0] <= surfaceinfo->xmax)
                        && (xy[0][1] >= surfaceinfo->ymin) && (xy[2][1] <= surfaceinfo->ymax);
        if(!insideDomain)
            return false;

        signed ab[8] = {0};
        double increment[8][3] = {{0,0,0},{l,0,0},{0,l,0},{l,l,0},{0,0,l},{l,0,l},{0,l,l},{l,l,l}};
#pragma omp parallel for
        for(int i = 0; i < 8; i++)
            ab[i] = aboveSurface(surfaceinfo->x, surfaceinfo->y, surfaceinfo->z,
                p8est->connectivity, tree, nx, ny,
                x+increment[i][0], y+increment[i][1], z+increment[i][2]);

        signed sum = 0;
        for(int i = 0; i < 8; i++)
            sum += ab[i];

        return (abs(sum) == 8) ? false : true;
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
    if(num_incoming == 4 && num_outgoing == 1)
    {
        // parent user data
        // quadrantData *childData = (quadrantData *) incoming[0]->p.user_data;
        quadrantData *parentData = (quadrantData *) outgoing[0]->p.user_data;

        // Averaging
        parentData->u = 0.0;
        for(int i = 0; i < 4; i++)
        {
            quadrantData *childData = (quadrantData *) incoming[i]->p.user_data;
            parentData->u+=childData->u;
        }
        parentData->u=0.25*parentData->u;

        // Tags
        quadrantData *childData = (quadrantData *) incoming[0]->p.user_data;
        parentData->quadTag=childData->quadTag;

        // Update the spatial coordinates
        p4est_qcoord_to_vertex(p4est->connectivity, tree, outgoing[0]->x, outgoing[0]->y, &parentData->xy[0]);
    }
    else if(num_incoming == 1 && num_outgoing == 4)
    {
        quadrantData *parentData = (quadrantData *) incoming[0]->p.user_data;

        // Loop over the four children
        for(int i = 0; i < 4; i++){
            quadrantData *childData = (quadrantData *) outgoing[i]->p.user_data;
            childData->u=parentData->u;
          	childData->quadTag=parentData->quadTag;

            // Update the spatial coordinates
            p4est_qcoord_to_vertex(p4est->connectivity, tree, outgoing[i]->x, outgoing[i]->y, &childData->xy[0]);
        }
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
    if(num_incoming == 8 && num_outgoing == 1)
    {
        // Get the parent quadrant
        p8est_quadrant_t * child[8];
        p8est_quadrant_childrenpv(outgoing[0], child);

        // parent user data
        octantData *childData = (octantData *) child[0]->p.user_data;
        octantData *parentData = (octantData *) outgoing[0]->p.user_data;

        parentData->u=childData->u;
        parentData->octantTag=childData->octantTag;

        // Update the spatial coordinates
        p8est_qcoord_to_vertex(p8est->connectivity, tree,
        outgoing[0]->x, outgoing[0]->y, outgoing[0]->z, &parentData->xyz[0]);
    }
    else if(num_incoming == 1 && num_outgoing == 8)
    {
        // Get the parent quadrant
        p8est_quadrant_t parent;
        p8est_quadrant_parent(incoming[0], &parent);

        // parent user data
        octantData *parentData = (octantData *) parent.p.user_data;

        for(int i = 0; i < 8; i++)
        {
            octantData *childData = (octantData *) incoming[i]->p.user_data;
            childData->u=parentData->u;
          	childData->octantTag=parentData->octantTag;

            // Update the spatial coordinates
            p8est_qcoord_to_vertex(p8est->connectivity, tree,
                incoming[i]->x, incoming[i]->y, incoming[i]->z, &childData->xyz[0]);
        }
    }
    else
    {
        throw OxleyException("refine_copy_parent_quadrant: Unknown error.");
    }
}

void get_interpolateNodesOnElementWorker_data(p4est_iter_volume_info_t * info, void *fxx)
{
    //AEAE todo

    // memcpy(&f_00[0], in.getSampleDataRO(INDEX2(k0,k1, m_NN[0]), sentinel), numComp*sizeof(S));
    // memcpy(&f_01[0], in.getSampleDataRO(INDEX2(k0,k1+1, m_NN[0]), sentinel), numComp*sizeof(S));
    // memcpy(&f_10[0], in.getSampleDataRO(INDEX2(k0+1,k1, m_NN[0]), sentinel), numComp*sizeof(S));
    // memcpy(&f_11[0], in.getSampleDataRO(INDEX2(k0+1,k1+1, m_NN[0]), sentinel), numComp*sizeof(S));

}

void get_interpolateNodesOnElementWorker_data(p8est_iter_volume_info_t * info, void *fxx)
{
    //AEAE todo
}


void get_interpolateNodesOnFacesWorker_data(p4est_iter_volume_info_t * info, void *fxx)
{
    //AEAE todo

//     node->x == 0 || node->x == P4EST_ROOT_LEN ||
//           node->y == 0 || node->y == P4EST_ROOT_LEN 

        //direction 0
//                     memcpy(&f_00[0], in.getSampleDataRO(INDEX2(0,k1, m_NN[0]), sentinel), numComp*sizeof(S));
//                     memcpy(&f_01[0], in.getSampleDataRO(INDEX2(0,k1+1, m_NN[0]), sentinel), numComp*sizeof(S));

    // face 1
    // memcpy(&f_10[0], in.getSampleDataRO(INDEX2(m_NN[0]-1,k1, m_NN[0]), sentinel), numComp*sizeof(S));
//                     memcpy(&f_11[0], in.getSampleDataRO(INDEX2(m_NN[0]-1,k1+1, m_NN[0]), sentinel), numComp*sizeof(S));


    //face 2

// memcpy(&f_00[0], in.getSampleDataRO(INDEX2(k0,0, m_NN[0]), sentinel), numComp*sizeof(S));
//                     memcpy(&f_10[0], in.getSampleDataRO(INDEX2(k0+1,0, m_NN[0]), sentinel), numComp*sizeof(S));


    //face 3
// memcpy(&f_01[0], in.getSampleDataRO(INDEX2(k0,m_NN[1]-1, m_NN[0]), sentinel), numComp*sizeof(S));
//                     memcpy(&f_11[0], in.getSampleDataRO(INDEX2(k0+1,m_NN[1]-1, m_NN[0]), sentinel), numComp*sizeof(S));


    /////////////////////////////////////////////////////////
    // shared
    /////////////////////////////////////////////////////////

    //face 0
    // memcpy(&f_00[0], in.getSampleDataRO(INDEX2(0,k1, m_NN[0]), sentinel), numComp*sizeof(S));
    // memcpy(&f_01[0], in.getSampleDataRO(INDEX2(0,k1+1, m_NN[0]), sentinel), numComp*sizeof(S));

    //face 1
    // memcpy(&f_10[0], in.getSampleDataRO(INDEX2(m_NN[0]-1,k1, m_NN[0]), sentinel), numComp*sizeof(S));
    //memcpy(&f_11[0], in.getSampleDataRO(INDEX2(m_NN[0]-1,k1+1, m_NN[0]), sentinel), numComp*sizeof(S));

    //face 2
// memcpy(&f_00[0], in.getSampleDataRO(INDEX2(k0,0, m_NN[0]), sentinel), numComp*sizeof(S));
//                     memcpy(&f_10[0], in.getSampleDataRO(INDEX2(k0+1,0, m_NN[0]), sentinel), numComp*sizeof(S));


    //face 3
// memcpy(&f_01[0], in.getSampleDataRO(INDEX2(k0,m_NN[1]-1, m_NN[0]), sentinel), numComp*sizeof(S));
//                     memcpy(&f_11[0], in.getSampleDataRO(INDEX2(k0+1,m_NN[1]-1, m_NN[0]), sentinel), numComp*sizeof(S));

}

// void get_interpolateNodesOnFacesWorker_data(p8est_iter_volume_info_t * info, void *fxx);
// {
//     //AEAE todo
// }


void update_node_faceoffset(p4est_iter_volume_info_t * info, void *fxx)
{
    quadrantData * quaddata = (quadrantData *) info->quad->p.user_data;

    // faceOffset[i]=-1 if face i is not an external face, otherwise it is
    // the index of that face (where i: 0=left, 1=right, 2=bottom, 3=top)
    // escript::DataTypes::IndexVector m_faceOffset;
    quaddata->m_faceOffset = -1;
    if(info->quad->x == 0)
        quaddata->m_faceOffset = 0;
    if(info->quad->x == P4EST_ROOT_LEN)
        quaddata->m_faceOffset = 1;
    if(info->quad->y == 0)
        quaddata->m_faceOffset = 2;
    if(info->quad->y == P4EST_ROOT_LEN)
        quaddata->m_faceOffset = 3;
}

void update_node_faceoffset(p8est_iter_volume_info_t * info, void *fxx)
{

}

void update_RC(p4est_iter_face_info_t *info, void *user_data)
{
    //Get some pointers
    update_RC_data * data = (update_RC_data *) user_data;
    sc_array_t * sides = &(info->sides);

    p4est_iter_face_side_t * side = p4est_iter_fside_array_index_int(sides, 0);
    
    p4est_quadrant_t * quad;
    if(side->is_hanging)
        quad = side->is.hanging.quad[0];
    else
        quad = side->is.full.quad;

    // First Node
    double xy[3];
    p4est_qcoord_to_vertex(data->p4est->connectivity, side->treeid, quad->x, quad->y, xy);
    long lni0 = data->pNodeIDs->find(std::make_pair(xy[0],xy[1]))->second;

    // Second node
    p4est_qcoord_t length = P4EST_QUADRANT_LEN(quad->level);
    int fn = (int) side->face;
    long lx = ((fn < 2) == 0) * length;
    long ly = ((fn >= 2) == 0) * length;
    p4est_qcoord_to_vertex(data->p4est->connectivity, side->treeid, quad->x+lx, quad->y+ly, xy);
    long lni1 = data->pNodeIDs->find(std::make_pair(xy[0],xy[1]))->second;

    std::vector<long> * idx0 = &data->indices[0][lni0];
    std::vector<long> * idx1 = &data->indices[0][lni1];

    bool dup = false;
    for(int i = 1; i < idx0[0][0] + 1; i++)
        if(idx0[0][i] == lni1)
        {
            dup = true;
            break;
        }

    if(dup == false)
    {
#ifdef DOXLEY_ENABLE_DEBUG
        std::cout << "update_RC " << lni1 << ": (" << xy[0] << ", " << xy[1] << ")" << std::endl; // coordinates
#endif
        idx0[0][0]++;
        idx1[0][0]++;
        idx0[0][idx0[0][0]]=lni1;
        idx1[0][idx1[0][0]]=lni0;
    }
}

void update_connections(p4est_iter_volume_info_t *info, void *user_data)
{
    //Get some pointers
    getConnections_data * data = (getConnections_data *) user_data;
    p4est_quadrant_t * quad = info->quad;
    
    // Coordinates
    p4est_qcoord_t length = P4EST_QUADRANT_LEN(quad->level);
    double xy[3];
    long lx[4] = {0,length,0,length};
    long ly[4] = {0,0,length,length};
    long lni[4] = {-1};
#pragma omp parallel for
    for(int i = 0; i < 4; i++)
    {
        p4est_qcoord_to_vertex(data->p4est->connectivity, info->treeid, quad->x+lx[i], quad->y+ly[i], xy);
        lni[i] = data->pNodeIDs->find(std::make_pair(xy[0],xy[1]))->second;
    }


#pragma omp parallel for
    for(int i = 0; i < 4; i++)
    {
        std::vector<escript::DataTypes::index_t> * temp = &data->indices[0][i];
        for(int j = 0; j < 4; j++)
        {
            bool dup = false;
            for(int k = 0; k < data->indices[0][i].size(); k++)
                if(temp[0][k] == lni[j])
                {
                    dup = true;
                    break;
                }
            if(dup == false)
                temp->push_back(lni[j]);
        }
    }

    // std::cout << "xy = " << xy[0] << ", " << xy[1] << std::endl; // coordinates


}

} // namespace oxley
