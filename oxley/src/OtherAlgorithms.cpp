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
#include <ctime>

#include <oxley/InitAlgorithms.h>
#include <oxley/OtherAlgorithms.h>
#include <oxley/RefinementAlgorithms.h>

#include <p4est.h>
#include <p4est_bits.h>
#include <p4est_extended.h>
#include <p4est_ghost.h>
#include <p4est_iterate.h>
#include <p4est_lnodes.h>

#include <p4est_algorithms.h> //aeae delete this

namespace oxley {

// Forward declarations
long getNewTag(OxleyDomainRect_ptr domain);
long getNewTag(OxleyDomainBrick_ptr domain);

#ifdef ESYS_HAVE_BOOST_NUMPY
void addSurface(OxleyDomainRect_ptr domain)
{
    p4est_t * p4est = domain->borrow_p4est();
    p4estData * forestData = (p4estData *) domain->borrow_forestData();

    // Get a new tag number and name
    addSurfaceData * surfacedata;
    surfacedata = new addSurfaceData;
    surfacedata = (addSurfaceData *) malloc(sizeof(addSurfaceData));
    surfacedata = (addSurfaceData *) domain->borrow_temp_data();
    surfacedata->newTag = getNewTag(domain);

    // Add to the list
    domain->numberOfTags++;
    domain->tags[domain->numberOfTags]=surfacedata->newTag;

    // Work out the boundaries of the function z[x,y]
    surfacedata->xmin = surfacedata->x[0];
    surfacedata->xmax = surfacedata->x[surfacedata->x.size()];
    for(int i = 0; i < surfacedata->x.size(); i++)
    {
        if(surfacedata->x[i] < surfacedata->xmin)
            surfacedata->xmin = surfacedata->x[i];
        if(surfacedata->x[i] > surfacedata->xmax)
            surfacedata->xmax = surfacedata->x[i];
    }

    // Note: the forest must be face balanced for p4est_iterate() to execute
    // a callback function on faces (see p4est_balance()).
    p4est_balance(p4est, P4EST_CONNECT_FACE, gce_init_new_rectangle);
    int partition_for_coarsening = 0;
    p4est_partition_ext(p4est, partition_for_coarsening, NULL);
#ifdef P4EST_ENABLE_DEBUG
    p4est_iterate(p4est, NULL, (void *) surfacedata, gce_first_pass, NULL, NULL);
#endif
    forestData->assign_info(surfacedata);
    p4est_iterate(p4est, NULL, (void *) surfacedata, gce_second_pass, NULL, NULL);
    p4est->user_pointer = surfacedata;
    p4est_refine_ext(p4est, true, forestData->max_levels_refinement,
        refine_gce, gce_init_new_rectangle, gce_rectangle_replace);
    p4est->user_pointer = forestData;

    // Balance and repartition
    p4est_balance(p4est, P4EST_CONNECT_FULL, gce_init_new_rectangle);
    p4est_iterate(p4est, NULL, (void *) surfacedata, gce_third_pass, NULL, NULL);
    p4est_partition_ext(p4est, partition_for_coarsening, NULL);
}
#endif

#ifdef ESYS_HAVE_BOOST_NUMPY
void addSurface(OxleyDomainBrick_ptr domain)
{
    p8est_t * p8est = domain->borrow_p8est();
    p8estData * forestData = (p8estData *) domain->borrow_forestData();

    // Get a new tag number and name
    addSurfaceData * surfacedata;
    surfacedata = new addSurfaceData;
    surfacedata = (addSurfaceData *) malloc(sizeof(addSurfaceData));
    surfacedata = (addSurfaceData *) domain->borrow_temp_data();
    surfacedata->newTag = getNewTag(domain);

    // Add to the list
    domain->numberOfTags++;
    domain->tags[domain->numberOfTags]=surfacedata->newTag;

    // Work out the boundaries of the function z[x,y]
    surfacedata->xmin = surfacedata->x[0];
    surfacedata->xmax = surfacedata->x[surfacedata->x.size()];
    for(int i = 0; i < surfacedata->x.size(); i++)
    {
        if(surfacedata->x[i] < surfacedata->xmin)
            surfacedata->xmin = surfacedata->x[i];
        if(surfacedata->x[i] > surfacedata->xmax)
            surfacedata->xmax = surfacedata->x[i];
    }
    surfacedata->ymin = surfacedata->y[0];
    surfacedata->ymax = surfacedata->y[surfacedata->y.size()];
    for(int i = 0; i < surfacedata->y.size(); i++)
    {
        if(surfacedata->y[i] < surfacedata->ymin)
            surfacedata->ymin = surfacedata->y[i];
        if(surfacedata->y[i] > surfacedata->ymax)
            surfacedata->ymax = surfacedata->y[i];
    }

    // Note: the forest must be face balanced for p4est_iterate() to execute
    // a callback function on faces
    p8est_balance(p8est, P8EST_CONNECT_FACE, gce_init_new_brick);
    int partition_for_coarsening = 0;
    p8est_partition_ext(p8est, partition_for_coarsening, NULL);
#ifdef P4EST_ENABLE_DEBUG
    p8est_iterate(p8est, NULL, (void *) surfacedata, gce_first_pass, NULL, NULL, NULL);
#endif
    forestData->assign_info(surfacedata);
    p8est_iterate(p8est, NULL, (void *) surfacedata, gce_second_pass, NULL, NULL, NULL);
    p8est->user_pointer = surfacedata;
    p8est_refine_ext(p8est, true, forestData->max_levels_refinement,
        refine_gce, gce_init_new_brick, gce_brick_replace);
    p8est->user_pointer = forestData;

    // Balance and repartition
    p8est_balance(p8est, P8EST_CONNECT_FACE, gce_init_new_brick);
    p8est_iterate(p8est, NULL, (void *) surfacedata, gce_third_pass, NULL, NULL, NULL);
    p8est_partition_ext(p8est, partition_for_coarsening, NULL);
}
#endif

// returns a randomly chosen new tag
long getNewTag(OxleyDomainRect_ptr domain)
{
    std::srand(time(0));
    int newTag;
    bool findingTag=true;
    const int * taglist = domain->borrowListOfTagsInUse(1);
    while(findingTag){
        newTag=rand() % MAXTAGS;
        findingTag=false;
        for(int i = 0; i < domain->getNumberOfTagsInUse(1); i++){
            if(taglist[i] == newTag){
                findingTag=true;
                break;
            }
        }
    }
    return newTag;
}

// returns a randomly chosen new tag
long getNewTag(OxleyDomainBrick_ptr domain)
{
    std::srand(time(0));
    int newTag;
    bool findingTag=true;
    const int * taglist = domain->borrowListOfTagsInUse(1);
    while(findingTag){
        newTag=rand() % MAXTAGS;
        findingTag=false;
        for(int i = 0; i < domain->getNumberOfTagsInUse(1); i++){
            if(taglist[i] == newTag){
                findingTag=true;
                break;
            }
        }
    }
    return newTag;
}

bool aboveCurve(double x[], double z[], int nx, double _x, double _z)
{
    // Find the index of the corresponding x nodes
    long ix1 = -1, ix2 = -1;
    for(long i = 0; i < nx; i++){
        if(_x >= x[i]){
            ix1=i;
            ix2=i+1;
            break;
        }
    }

    // Point is outside the domain
    if(ix1 == -1 && ix2 == -1)
        return false;

    // Do the check
    if(x[ix1] == _x) // If the point is on the node
    {
        return _z > z[ix1];
    }
    else // otherwise, interpolate
    {
        double z0=z[ix1], z1=z[ix2];
        double x0=x[ix1], x1=x[ix2];
        double tmp1 = (z0*(x1-_x)+z1*(_x-x0))/(x1-x0);
        return _z > tmp1;
    }
}

bool aboveCurve(std::vector<double> x, std::vector<double> y,
                    p4est_connectivity_t * connectivity, p4est_topidx_t treeid,
                    long n, p4est_qcoord_t qx, p4est_qcoord_t qy)
{
    // Get the spatial coordinates of the point we are considering
    double xy[2] = {-1.0,-1.0};
    p4est_qcoord_to_vertex(connectivity, treeid, qx, qy, xy);
    double _x = xy[0];
    double _y = xy[1];

    // Find the indices of the corresponding x and y nodes
    long ix1 = -1, ix2 = -1;

    for(long i = 0; i < n; i++){
        if(_x > x[i] && _x < x[i+1])
        {
            ix1=i;
            ix2=i+1;
            break;
        } else if (_x > x[i] && i == n)
        {
            ix1=i-1;
            ix2=i;
            break;
        } else if(_x == x[i]) {
            return _y >= y[i] ? true : false;
        }
    }

    // Point is outside the domain
    if(ix1 == -1)
        return false;

    // interpolate
    double x1=x[ix1], x2=x[ix2];
    double y1=y[ix1], y2=y[ix2];
    double tmp1 = _x*(y1-y2)/(x1-x2);
    double tmp2 = (x2*y1-x1*y2)/(x1-x2);
    return _y >= tmp1 - tmp2 ? true : false;
}

double distanceToCurve(double x[], double z[], int nx, double _x, double _z)
{
    // Find the index of the corresponding x nodes
    long ix1, ix2 = 0;

    for(long i = 0; i < nx; i++){
        if(_x >= x[i]){
            ix1=i;
            ix2=i+1;
            break;
        }
    }

    // Do the check
    if(x[ix1] == _x) // If the point is on the node
    {
        return std::abs(_z - z[ix1]);
    }
    else // otherwise, interpolate
    {
        double z0=z[ix1], z1=z[ix2];
        double x0=x[ix1], x1=x[ix2];
        double tmp1 = (z0*(x1-_x)+z1*(_x-x0))/(x1-x0);
        return std::abs(_z - tmp1);
    }
}

bool aboveSurface(double x[], double y[], double z[], int nx, int ny, double _x, double _y, double _z)
{
    // Find the indices of the corresponding x and y nodes
    long ix1, iy1, ix2, iy2 = -1;

    for(long i = 0; i < nx; i++){
        if(_x >= x[i]){
            ix1=i;
            ix2=i+1;
            break;
        }
    }

    for(long j = 0; j < ny; j++){
        if(_y >= y[j]){
            iy1=j;
            iy2=j+1;
            break;
        }
    }

    // If the point is outside the domain of the vectors x and y
    if(ix1 == -1)
        return false;
    if(iy1 == -1)
        return false;

    // Do the check
    if(x[ix1] == _x && y[iy1] == _y) // If the point is on the node
    {
        return _z > z[INDEX2(ix1,iy1,nx)];
    }
    else // otherwise, interpolate
    {
        double q11, q12, q21, q22 = 0;
        q11=z[INDEX2(ix1,ix1,nx)]; q12=z[INDEX2(ix1,iy2,nx)];
        q21=z[INDEX2(ix2,iy1,nx)]; q22=z[INDEX2(ix2,iy2,nx)];
        double x1, x2, y1, y2 = 0;
        x1=x[ix1];y1=y[iy1];x2=x[ix2];y2=y[iy2];

        double tmp1 = ((x2-_x)/(x2-x1))*q11+((_x-x1)/(x2-x1))*q21;
        double tmp2 = ((x2-_x)/(x2-x1))*q12+((_x-x1)/(x2-x1))*q22;
        double tmp3 = ((y2-_y)/(y2-y1))*tmp1+((_y-y1)/(y2-y1))*tmp2;
        return _z > tmp3;
    }
}

bool aboveSurface(std::vector<double> x, std::vector<double> y, std::vector<double> z,
                    p8est_connectivity_t * connectivity, p4est_topidx_t treeid,
                    long nx, long ny, p4est_qcoord_t qx, p4est_qcoord_t qy, p4est_qcoord_t qz)
{
    // Get the spatial coordinates of the point we are considering
    double xyz[3] = {-1.0,-1.0,-1.0};
    p8est_qcoord_to_vertex(connectivity, treeid, qx, qy, qz, xyz);
    double _x = xyz[0];
    double _y = xyz[1];
    double _z = xyz[2];

    // Find the indices of the x and y nodes that box in the point
    long ix1 = -1, iy1 = -1, ix2 = -1, iy2 = -1;
    for(long i = 0; i < nx; i++){
        if(_x > x[i] && _x < x[i+1])
        {
            ix1=i;
            ix2=i+1;
            break;
        }
        else if (_x > x[i] && i == nx)
        {
            ix1=i-1;
            ix2=i;
            break;
        } else if(_x == x[i]) {
            ix1=i;
            break;
        }
    }

    for(long j = 0; j < ny; j++){
        if(_y > y[j] && _y < y[j+1])
        {
            iy1=j;
            iy2=j+1;
            break;
        }
        else if(_y > y[j] && j == ny)
        {
            iy1=j-1;
            iy2=j;
            break;
        } else if(_y == y[j]) {
            iy1=j;
            break;
        }
    }

    // Point is outside the domain
    if(ix1 == -1 && ix2 == -1)
        return false;
    if(iy1 == -1 && iy2 == -1)
        return false;

    // Do the check
    if(x[ix1] == _x && y[iy1] == _y) // If the point is on the node
    {
        return _z - z[INDEX2(ix1,iy1,nx)] < 0 ? true : false;
    }
    else // this is bilinear interpolation
    {
        double q11 = 0.0, q12 = 0.0, q21 = 0.0, q22 = 0.0;
        q11=z[INDEX2(ix1,iy1,nx)]; q12=z[INDEX2(ix1,iy2,nx)];
        q21=z[INDEX2(ix2,iy1,nx)]; q22=z[INDEX2(ix2,iy2,nx)];
        double x1=x[ix1], y1=y[iy1], x2=x[ix2], y2=y[iy2];

        double tmp1 = ((x2-_x)/(x2-x1))*q11+((_x-x1)/(x2-x1))*q21;
        double tmp2 = ((x2-_x)/(x2-x1))*q12+((_x-x1)/(x2-x1))*q22;
        double tmp3 = ((y2-_y)/(y2-y1))*tmp1+((_y-y1)/(y2-y1))*tmp2;

        return _z - tmp3 < 0 ? true : false;
    }
}

double distanceToSurface(double x[], double y[], double z[], int nx, int ny, double _x, double _y, double _z)
{
    // Find the indices of the corresponding x and y nodes
    long ix1, iy1, ix2, iy2 = 0;

    for(long i = 0; i < nx; i++){
        if(_x >= x[i]){
            ix1=i;
            ix2=i+1;
            break;
        }
    }

    for(long j = 0; j < ny; j++){
        if(_y >= y[j]){
            iy1=j;
            iy2=j+1;
            break;
        }
    }

    // Do the check
    if(x[ix1] == _x && y[iy1] == _y) // If the point is on the node
    {
        return abs(_z - z[INDEX2(ix1,iy1,nx)]);
    }
    else // otherwise, interpolate
    {
        double q11 = 0.0, q12 = 0.0, q21 = 0.0, q22 = 0.0;
        q11=z[INDEX2(ix1,ix1,nx)]; q12=z[INDEX2(ix1,iy2,nx)];
        q21=z[INDEX2(ix2,iy1,nx)]; q22=z[INDEX2(ix2,iy2,nx)];
        double x1, x2, y1, y2 = 0;
        x1=x[ix1];y1=y[iy1];x2=x[ix2];y2=y[iy2];

        double tmp1 = ((x2-_x)/(x2-x1))*q11+((_x-x1)/(x2-x1))*q21;
        double tmp2 = ((x2-_x)/(x2-x1))*q12+((_x-x1)/(x2-x1))*q22;
        double tmp3 = ((y2-_y)/(y2-y1))*tmp1+((_y-y1)/(y2-y1))*tmp2;
        return abs(_z - tmp3);
    }
}



} // end namespace oxley
