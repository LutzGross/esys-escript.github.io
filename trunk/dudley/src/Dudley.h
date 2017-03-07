
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

#ifndef __DUDLEY_H__
#define __DUDLEY_H__

/****************************************************************************/

/*    Dudley finite element solver */

/****************************************************************************/

#include <escript/DataTypes.h>

#include <dudley/DudleyException.h>

#include <escript/Data.h>
#include <escript/EsysMPI.h>

namespace dudley {

using escript::DataTypes::index_t;
using escript::DataTypes::dim_t;
using escript::DataTypes::IndexVector;

#define DUDLEY_UNKNOWN -1
#define DUDLEY_DEGREES_OF_FREEDOM 1
#define DUDLEY_NODES 3
#define DUDLEY_ELEMENTS 4
#define DUDLEY_FACE_ELEMENTS 5
#define DUDLEY_POINTS 6
#define DUDLEY_REDUCED_ELEMENTS 10
#define DUDLEY_REDUCED_FACE_ELEMENTS 11

//
// Codes for function space types supported
enum {
    DegreesOfFreedom = DUDLEY_DEGREES_OF_FREEDOM,
    Nodes = DUDLEY_NODES,
    Elements = DUDLEY_ELEMENTS,
    ReducedElements = DUDLEY_REDUCED_ELEMENTS,
    FaceElements = DUDLEY_FACE_ELEMENTS,
    ReducedFaceElements = DUDLEY_REDUCED_FACE_ELEMENTS,
    Points = DUDLEY_POINTS
};

inline bool hasReducedIntegrationOrder(const escript::Data& in)
{
    const int fs = in.getFunctionSpace().getTypeCode();
    return (fs == ReducedElements || fs == ReducedFaceElements);
}

#define DUDLEY_INITIAL_STATUS 0

} // namespace dudley

#endif // __DUDLEY_H__

