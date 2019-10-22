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

#include <escript/DataTypes.h>

#include <boost/python/numpy.hpp>

#include <p4est_iterate.h>
#include <p8est_iterate.h>

#ifndef __OXLEY_DATA_H__
#define __OXLEY_DATA_H__

// Macro for array indexing
#define INDEX2(_X1_,_X2_,_N1_) ((_X1_)+(_N1_)*(_X2_))

////////////////////////////////////////////////////////////////////////
// This file contains the data structures used by Rectangle and Brick
////////////////////////////////////////////////////////////////////////

//This structure describes the information that is stored at each
//quadrant / octant in the p4est / p8est
struct quadrantData
{
	double u = 0.0;

	// The quadrant's tag
	long quadTag = 0;

	// Node tag
	double nodeTag = 0;

	// Spatial coordinates of the corner node that defines the quadrant
	double xy[2] = {0.0,0.0};
};

struct octantData
{
	double u = 0.0; // A Scalar variable

	// The octant's tag
	long octantTag = 0;

	// Node tag
	double nodeTag = 0;

	// Spatial coordinates of the corner node that defines the octant
	double xyz[3] = {0.0,0.0,0.0};
};

//This structure describes the information that is stored with the p4est
struct p4estData
{
	// This is here to temporarily store information
	void * info;

	// origin of domain
    double m_origin[2] = {0.0,0.0};

    // side lengths of domain
    double m_length[2] = {0.0,0.0};

    // number of spatial subdivisions
    int m_NX[2] = {0,0};

    // total number of elements in each dimension
    escript::DataTypes::dim_t m_gNE[2] = {0,0};

    // number of elements for this rank in each dimension including shared
    escript::DataTypes::dim_t m_NE[2] = {0,0};

    // periodic boundary conditions
    bool periodic[2] {false, false};

	// maximum levels of recursion to use during refinement
	int max_levels_refinement = 0;
};

struct p8estData
{
	// This is here to temporarily store information
	void * info;

	// origin of domain
    double m_origin[3] = {0.0,0.0,0.0};

    // side lengths of domain
    double m_length[3] = {0.0,0.0,0.0};

    // number of spatial subdivisions
    int m_NX[3] = {0,0,0};

    // total number of elements in each dimension
    escript::DataTypes::dim_t m_gNE[3] = {0,0,0};

    // number of elements for this rank in each dimension including shared
    escript::DataTypes::dim_t m_NE[3] = {0,0,0};

    // periodic boundary conditions
    bool periodic[3] {false, false, false};

	// maximum levels of recursion to use during refinement
	int max_levels_refinement = 0;
};

// This structure temporarily stores information used by the addSurface function
struct addSurfaceData {

	int oldTag = -1;
	int newTag = -1;
	std::vector<double> x;
	std::vector<double> y;
	std::vector<double> z;

};

namespace oxley {

// Call back function that copies quadrant tags onto tagVector
void getQuadTagVector(p4est_iter_volume_info_t * info, void *tagVector);
void getQuadTagVector(p8est_iter_volume_info_t * info, void *tagVector);

// Call back function that copies coordinate info onto tagVector
void getXCoordVector(p4est_iter_volume_info_t * info, void *tagVector);
void getYCoordVector(p4est_iter_volume_info_t * info, void *tagVector);
void getXCoordVector(p8est_iter_volume_info_t * info, void *tagVector);
void getYCoordVector(p8est_iter_volume_info_t * info, void *tagVector);
void getZCoordVector(p8est_iter_volume_info_t * info, void *tagVector);


void getNodeTagVector(p4est_iter_volume_info_t * info, void *tagVector);
void getNodeTagVector(p8est_iter_volume_info_t * info, void *tagVector);

} //namespace oxley

#endif
