/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* https://github.com/LutzGross/esys-escript.github.io
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
*
*****************************************************************************/

#include <oxley/Brick.h>
#include <oxley/Rectangle.h>

#include <boost/python/numpy.hpp>

////////////////////////////////////////////////////////////////////////
// This file contains other algorithms that operate on but are not part of
// Rectangle and Brick
////////////////////////////////////////////////////////////////////////

namespace oxley {

#ifdef ESYS_HAVE_BOOST_NUMPY
void addSurface(OxleyDomainRect_ptr domain);
void addSurface(OxleyDomainBrick_ptr domain);
#endif

// Gets a random new tag number
long getNewTag(OxleyDomainRect_ptr domain);
long getNewTag(OxleyDomainBrick_ptr domain);

// Returns true if point (_x,_z) lies above the surface defined by x[], z[]
// When necessary, this uses linear interpolation
// Array x has nx elements and array z both have nx elements.
signed aboveCurve(std::vector<double> x, std::vector<double> y,
                    p4est_connectivity_t * connectivity, p4est_topidx_t treeid,
                    long n, p4est_qcoord_t _x, p4est_qcoord_t _y);


// Returns the distance between a point and the interpolated curve
double distanceToCurve(double x[], double z[], int nx, double _x, double _z);

// Returns true if point (_x,_y,_z) lies above the surface defined by x[], y[], z[]
// When necessary, this uses bilinear interpolation
// Array x has nx elements and array y has ny elements. Array z has nx*ny elements.
signed aboveSurface(std::vector<double> x, std::vector<double> y, std::vector<double> z,
                    p8est_connectivity_t * connectivity, p4est_topidx_t treeid,
                    long nx, long ny, p4est_qcoord_t _x, p4est_qcoord_t _y, p4est_qcoord_t _z);

// Returns the distance between a point and the interpolated surface
double distanceToSurface(double x[], double y[], double z[],
                    int nx, int ny, double _x, double _y, double _z);

// Returns true if the quadrant/octant is on the boundary
bool onBoundary(p4est_quadrant_t * quadrant);
bool onBoundary(p8est_quadrant_t * octant);

// Returns true if the quadrant/octant has a hanging node
bool isHanging(p4est_quadrant_t * quadrant);
bool isHanging(p8est_quadrant_t * octant);

///////////////////////////////////////////////////////////////////////////

} // end namespace oxley
