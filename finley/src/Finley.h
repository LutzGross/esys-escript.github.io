
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
**
*****************************************************************************/

#ifndef __FINLEY_H__
#define __FINLEY_H__

/****************************************************************************

  Finley finite element solver

*****************************************************************************/

// first include to avoid _POSIX_C_SOURCE redefinition warnings
#include <escript/DataTypes.h>

#include <finley/FinleyException.h>

#include <escript/EsysMPI.h>

#include <vector>

namespace finley {

using escript::DataTypes::dim_t;
using escript::DataTypes::index_t;
using escript::DataTypes::IndexVector;

// real_t clashes with metis real_t !
//using escript::DataTypes::real_t;

//#define Finley_TRACE
#define FINLEY_UNKNOWN -1
#define FINLEY_DEGREES_OF_FREEDOM 1
#define FINLEY_NODES 3
#define FINLEY_ELEMENTS 4
#define FINLEY_FACE_ELEMENTS 5
#define FINLEY_POINTS 6
#define FINLEY_CONTACT_ELEMENTS_1 7
#define FINLEY_CONTACT_ELEMENTS_2 8
#define FINLEY_REDUCED_DEGREES_OF_FREEDOM 2
#define FINLEY_REDUCED_NODES 14
#define FINLEY_REDUCED_ELEMENTS 10
#define FINLEY_REDUCED_FACE_ELEMENTS 11
#define FINLEY_REDUCED_CONTACT_ELEMENTS_1 12
#define FINLEY_REDUCED_CONTACT_ELEMENTS_2 13

enum {
    DegreesOfFreedom = FINLEY_DEGREES_OF_FREEDOM,
    ReducedDegreesOfFreedom = FINLEY_REDUCED_DEGREES_OF_FREEDOM,
    Nodes = FINLEY_NODES,
    ReducedNodes = FINLEY_REDUCED_NODES,
    Elements = FINLEY_ELEMENTS,
    ReducedElements = FINLEY_REDUCED_ELEMENTS,
    FaceElements = FINLEY_FACE_ELEMENTS,
    ReducedFaceElements = FINLEY_REDUCED_FACE_ELEMENTS,
    Points = FINLEY_POINTS,
    ContactElementsZero = FINLEY_CONTACT_ELEMENTS_1,
    ReducedContactElementsZero = FINLEY_REDUCED_CONTACT_ELEMENTS_1,
    ContactElementsOne = FINLEY_CONTACT_ELEMENTS_2,
    ReducedContactElementsOne = FINLEY_REDUCED_CONTACT_ELEMENTS_2
};

#define FINLEY_INITIAL_STATUS 0

} // namespace finley

#endif // __FINLEY_H__

