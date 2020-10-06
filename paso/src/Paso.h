
/*****************************************************************************
*
* Copyright (c) 2003-2020 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014-2017 by Centre for Geoscience Computing (GeoComp)
* Development from 2019 by School of Earth and Environmental Sciences
**
*****************************************************************************/


/****************************************************************************/

/*    Paso finite element solver library                                    */

/****************************************************************************/

/*  Copyrights by ACcESS Australia, 2003,2004,2005 */
/*  Author: Lutz Gross, l.gross@uq.edu.au */

/****************************************************************************/

#ifndef __PASO_H__
#define __PASO_H__

#include "system_dep.h"
#include <escript/index.h>
#include <escript/DataTypes.h>
#include <escript/EsysMPI.h>

#include <boost/enable_shared_from_this.hpp>
#include <boost/shared_ptr.hpp>

namespace paso {

// return codes used by the solvers
enum SolverResult {
    NoError = 0,
    MaxIterReached,
    InputError,
    MemoryError,
    Breakdown,
    NegativeNormError,
    Divergence
};

using escript::DataTypes::dim_t;
using escript::DataTypes::index_t;
using escript::DataTypes::real_t;
using escript::DataTypes::cplx_t;

}

#define MATRIX_FORMAT_DEFAULT 1
#define MATRIX_FORMAT_CSC 2
#define MATRIX_FORMAT_BLK1 4
#define MATRIX_FORMAT_OFFSET1 8
#define MATRIX_FORMAT_DIAGONAL_BLOCK 32
#define MATRIX_FORMAT_COMPLEX 64

#define PASO_ONE (double)(1.0)
#define PASO_ZERO (double)(0.0)

#endif // __PASO_H__

