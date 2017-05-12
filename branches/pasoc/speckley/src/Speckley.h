
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

enum {
    DegreesOfFreedom=1,
    ReducedDegreesOfFreedom=2,
    Nodes=3,
    ReducedNodes=14,
    Elements=4,
    ReducedElements=10,
    FaceElements=5,
    ReducedFaceElements=11,
    Points=6
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

