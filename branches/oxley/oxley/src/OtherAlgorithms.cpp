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

namespace oxley {

// Forward declarations
long getNewTag(OxleyDomainRect_ptr domain);
long getNewTag(OxleyDomainBrick_ptr domain);

#ifdef ESYS_HAVE_BOOST_NUMPY
void addSurface(OxleyDomainRect_ptr domain)
{
    p4estData *forestData = (p4estData *) domain->p4est->user_pointer;

    // Get a new tag number and name
    addSurfaceData * surfacedata = (addSurfaceData *) forestData->info;
    surfacedata->newTag = getNewTag(domain);

    // Add to the list
    domain->numberOfTags++;
    domain->tags[domain->numberOfTags]=surfacedata->newTag;

    // int input_size = input.shape(0);
    // double* input_ptr = reinterpret_cast<double*>(input.get_data());
    // std::vector<double> v(input_size);
    // for (int i = 0; i < input_size; ++i)
    //     v[i] = *(input_ptr + i);

    // In order to minimise the number of function calls, this code has three
    // separate loops
    // The first loop tags the quadrant's corner node, as being above or below
    // the curve and then exits
    // The seconds loop records whether or not a quadrant should be refined
    // based on this info and info from neighbouring quads
    // The third loop does the refinement and updates the new quads
    p4est_iterate(domain->p4est, NULL, NULL, gce_first_pass, NULL, NULL);
#ifdef P4EST_ENABLE_DEBUG
    std::string filename = "first_pass";
    domain->writeToVTK(filename, false);
#endif
    p4est_iterate(domain->p4est, NULL, NULL, gce_second_pass, NULL, NULL);
#ifdef P4EST_ENABLE_DEBUG
    filename = "second_pass";
    domain->writeToVTK(filename, false);
#endif
    p4est_refine_ext(domain->p4est, true, forestData->max_levels_refinement,
        refine_gce, init_rectangle_data, gce_rectangle_replace);
#ifdef P4EST_ENABLE_DEBUG
    filename = "refinement";
    domain->writeToVTK(filename, false);
#endif

    p4est_balance_ext(domain->p4est, P4EST_CONNECT_FULL,
        init_rectangle_data, gce_rectangle_replace);

    int partition_for_coarsening = 0; //Do not allow coarsening while partitioning
    p4est_partition_ext(domain->p4est, partition_for_coarsening, NULL);

    // clean up
    delete surfacedata;
}
#endif

#ifdef ESYS_HAVE_BOOST_NUMPY
void addSurface(OxleyDomainBrick_ptr domain)
{
    p8estData *forestData = (p8estData *) domain->p8est->user_pointer;

    // Get a new tag number and name
    addSurfaceData * surfacedata = (addSurfaceData *) forestData->info;
    surfacedata->newTag = getNewTag(domain);

    // Add to the list
    domain->numberOfTags++;
    domain->tags[domain->numberOfTags]=surfacedata->newTag;

    // In order to minimise the number of function calls, this code has three
    // separate loops
    // The first loop tags the quadrant's corner node, as being above or below
    // the curve and then exits
    // The seconds loop records whether or not a quadrant should be refined
    // based on info from the first loop.
    // The third loop does the refinement and updates the new quads
    p8est_iterate(domain->p8est, NULL, NULL, gce_first_pass, NULL, NULL, NULL);
    p8est_iterate(domain->p8est, NULL, NULL, gce_second_pass, NULL, NULL, NULL);
    p8est_refine_ext(domain->p8est, true, forestData->max_levels_refinement,
        refine_gce, init_brick_data, gce_brick_replace);

    p8est_balance_ext(domain->p8est, P8EST_CONNECT_FULL,
        init_brick_data, gce_brick_replace);

    int partition_for_coarsening = 0; //Do not allow coarsening while partitioning
    p8est_partition_ext(domain->p8est, partition_for_coarsening, NULL);

    // clean up
    delete surfacedata;
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

bool aboveCurve(std::vector<double> x, std::vector<double> z, int nx, double _x, double _z)
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

    if(ix1 == -1 && ix2 == -1)
        return false;

    // Do the check
    if(x[ix1] == _x) // If the point is on the node
    {
        return _z > z[ix1];
    }
    else // otherwise, interpolate
    {
        double z0, z1, x0, x1;
        z0=z[ix1]; z1=z[ix2];
        x0=x[ix1]; x1=x[ix2];
        return _z > (z0*(x1-_x)+z1*(_x-x0))/(x1-x0);
    }
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

    // Point is outside the domain
    if(ix1 == -1 && ix2 == -1)
        return false;
    if(iy1 == -1 && iy2 == -1)
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
                    int nx, int ny, double _x, double _y, double _z)
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

    // Point is outside the domain
    if(ix1 == -1 && ix2 == -1)
        return false;
    if(iy1 == -1 && iy2 == -1)
        return false;

    // Do the check
    if(x[ix1] == _x && y[iy1] == _y) // If the point is on the node
    {
        double temp = _z - z[INDEX2(ix1,iy1,nx)];
        return std::abs(temp);
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
        return std::abs(_z - tmp3);
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
        double q11, q12, q21, q22 = 0;
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
