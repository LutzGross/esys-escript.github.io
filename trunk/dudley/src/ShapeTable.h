
/*****************************************************************************
*
* Copyright (c) 2003-2017 by The University of Queensland
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

/* Shape Function info
These tables are a much simplified version of content from finley's ShapeFunctions files

This file is not to be included in .h files - only .cpp files should have any use for it
*/

#ifndef __DUDLEY_SHAPETABLE_H__
#define __DUDLEY_SHAPETABLE_H__

#include "Dudley.h"
#include "ElementType.h"

namespace dudley {

// These are constructed from dsdv in ShapeFunction.cpp in finley.
// The first one is just there for functions that want a pointer
static const double DTDV_1D[2][2] = { {-1., 1}, {-1., 1.} };

// The repetition here is a hack to prevent out of bounds access
static const double DTDV_2D[3 * 3][2] = {
    {-1, 1}, {0, -1.}, {0, 1},
    {-1, 1}, {0, -1.}, {0, 1},
    {-1, 1}, {0, -1.}, {0, 1}
};

static const double DTDV_3D[4][3] = {
    {-1, -1, -1},
    { 1,  0,  0},
    { 0,  1,  0},
    { 0,  0,  1}
};

// Index the following by ElementTypeID
// The number of local dimensions (as opposed to dimension of the embedding
// space)
static const int localDims[8] = { 0, 1, 2, 3, 0, 1, 2, 0 };
static const int Dims[8] = { 0, 1, 2, 3, 1, 2, 3, 0 };

// the following lists are only used for face elements defined by
// numNodesOnFace>0

// if the element is allowed as a face element, numNodesOnFace defines the
// number of nodes defining the face
static const int numNodesOnFaceMap[8] = { 1, 2, 3, 4, 1, 2, 4, -1 };

// defines a permutation of the nodes which rotates the nodes on the face
static const int shiftNodesMap[8][4] = { {0}, {1, 0}, {1, 2, 0}, {-1}, {0, 1, 2}, {1, 0, 2}, {1, 2, 0, 3}, {0} };

// reverses the order of the nodes on a face. the permutation has keep 0 fixed.
// shiftNodes={-1} or reverseNodes={-1} are ignored.
static const int reverseNodesMap[8][4] = { {-1}, {-1}, {0, 2, 1}, {-1}, {-1}, {-1}, {0, 2, 1, 3}, {0} };

// [0] is reduced quadrature, [1] is full quadrature
// in order the positions are POINT, LINE, TRI, TET
static const double QuadWeight[4][2] = { {0, 0}, {1., 0.5}, {0.5, 1. / 6}, {1. / 6, 1. / 24} };

// number of quadrature points per element
static const int QuadNums[4][2] = { {0, 0}, {1, 2}, {1, 3}, {1, 4} };

// shape functions at quadrature nodes
bool getQuadShape(dim_t sim, bool reduced, const double **shapearr);

const char* getElementName(ElementTypeId id);

} // namespace dudley

#endif // __DUDLEY_SHAPETABLE_H__

