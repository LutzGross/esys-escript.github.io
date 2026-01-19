
/*****************************************************************************
*
* Copyright (c) 2003-2020 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/

#ifndef __SPECKLEY_SPECKLEY_H__
#define __SPECKLEY_SPECKLEY_H__

/*****************************************************************************
 *  Speckley is a spectral element domain library with regular 
 *  hexagonal/rectangular elements and high order quadrature points
 ****************************************************************************/

#include <speckley/system_dep.h>

#include <escript/EsysMPI.h>

#include <boost/shared_ptr.hpp>
#include <list>
#include <map>
#include <string>
#include <vector>

namespace speckley {

using escript::DataTypes::dim_t;
using escript::DataTypes::index_t;
using escript::DataTypes::real_t;
using escript::DataTypes::cplx_t;

typedef std::pair<index_t,index_t> IndexPair;
typedef std::vector<index_t> IndexVector;
typedef std::vector<real_t> DoubleVector;
typedef std::map<std::string,int> TagMap;

#define SPECKLEY_DEGREES_OF_FREEDOM 1
#define SPECKLEY_NODES 3
#define SPECKLEY_ELEMENTS 4
#define SPECKLEY_FACE_ELEMENTS 5
#define SPECKLEY_POINTS 6
#define SPECKLEY_REDUCED_DEGREES_OF_FREEDOM 2
#define SPECKLEY_REDUCED_NODES 14
#define SPECKLEY_REDUCED_ELEMENTS 10
#define SPECKLEY_REDUCED_FACE_ELEMENTS 11

enum {
    DegreesOfFreedom=SPECKLEY_DEGREES_OF_FREEDOM,
    ReducedDegreesOfFreedom=SPECKLEY_REDUCED_DEGREES_OF_FREEDOM,
    Nodes=SPECKLEY_NODES,
    ReducedNodes=SPECKLEY_REDUCED_NODES,
    Elements=SPECKLEY_ELEMENTS,
    ReducedElements=SPECKLEY_REDUCED_ELEMENTS,
    FaceElements=SPECKLEY_FACE_ELEMENTS,
    ReducedFaceElements=SPECKLEY_REDUCED_FACE_ELEMENTS,
    Points=SPECKLEY_POINTS
};

//quadrature point locations
const double point_locations[][11] = {
        {0.0, 0.5, 1.0},
        {0.0, 0.27639320225, 0.72360679775, 1.0},
        {0.0, 0.172673164646, 0.5, 0.827326835354, 1.0},
        {0.0, 0.117472338035, 0.35738424176, 0.64261575824, 0.882527661965, 1.0},
        {0.0, 0.0848880518607, 0.265575603265, 0.5, 0.734424396735, 0.915111948139, 1.0},
        {0.0, 0.0641299257452, 0.204149909283, 0.395350391049, 0.604649608951, 0.795850090717, 0.935870074255, 1.0},
        {0.0, 0.0501210022943, 0.161406860245, 0.318441268087, 0.5, 0.681558731913, 0.838593139755, 0.949878997706, 1.0},
        {0.0, 0.0402330459168, 0.130613067447, 0.261037525095, 0.417360521167, 0.582639478833, 0.738962474905, 0.869386932553, 0.959766954083, 1.0},
        {0.0, 0.032999284796, 0.107758263168, 0.217382336502, 0.352120932207, 0.5, 0.647879067793, 0.782617663498, 0.892241736832, 0.967000715204, 1.0}};

} // namespace speckley

#endif /* __SPECKLEY_SPECKLEY_H__ */

