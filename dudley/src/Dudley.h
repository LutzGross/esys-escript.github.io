
/*****************************************************************************
*
* Copyright (c) 2003-2016 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
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

#include <escript/EsysMPI.h>
#include <cstring>

namespace dudley {

using escript::DataTypes::index_t;
using escript::DataTypes::dim_t;

#define DUDLEY_UNKNOWN -1
#define DUDLEY_DEGREES_OF_FREEDOM 1
#define DUDLEY_NODES 3
#define DUDLEY_ELEMENTS 4
#define DUDLEY_FACE_ELEMENTS 5
#define DUDLEY_POINTS 6
#define DUDLEY_REDUCED_DEGREES_OF_FREEDOM 2
#define DUDLEY_REDUCED_NODES 14
#define DUDLEY_REDUCED_ELEMENTS 10
#define DUDLEY_REDUCED_FACE_ELEMENTS 11

/* status stuff */
typedef int Dudley_Status_t;
#define Dudley_increaseStatus(self) ((self)->status)++
#define DUDLEY_INITIAL_STATUS 0

} // namespace dudley

#endif // __DUDLEY_H__

