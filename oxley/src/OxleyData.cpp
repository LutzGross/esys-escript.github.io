
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


namespace oxley {

void getTagVector(p4est_iter_volume_info_t * info, void *tagVector)
{
    sc_array_t *tagNumber = (sc_array_t *) tagVector;

    // Get the p4est
    p4est_t *p4est = info->p4est;

    // Get the quadrant data
    p4est_quadrant_t *q = info->quad;
    quadrantData *data = (quadrantData *) q->p.user_data;

    // Get a pointer to the tree
    p4est_tree_t * tree = p4est_tree_array_index(p4est->trees, info->treeid);
    p4est_locidx_t local_id = info->quadid;
    local_id += tree->quadrants_offset;

    // Copy over the data
    double * tmp = (double *) sc_array_index(tagNumber, local_id);

    tmp[0] = data->quadTag;
}

void getTagVector(p8est_iter_volume_info_t * info, void *tagVector)
{
    sc_array_t *tagNumber = (sc_array_t *) tagVector;

    // Get the p4est
    p8est_t *p8est = info->p4est;

    // Get the quadrant data
    p8est_quadrant_t *q = info->quad;
    octantData *data = (octantData *) q->p.user_data;

    // Get a pointer to the tree
    p8est_tree_t * tree = p8est_tree_array_index(p8est->trees, info->treeid);
    p4est_locidx_t local_id = info->quadid;
    local_id += tree->quadrants_offset;

    // Copy over the data
    double * tmp = (double *) sc_array_index(tagNumber, local_id);

    tmp[0] = data->octantTag;
}


} //namespace oxley
