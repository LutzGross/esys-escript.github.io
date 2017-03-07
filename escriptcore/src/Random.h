/*****************************************************************************
*
* Copyright (c) 2013-2017 by The University of Queensland
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

#ifndef __ESCRIPT_RANDOM_H__
#define __ESCRIPT_RANDOM_H__

namespace escript
{
/* \brief put n random doubles (from [0.0, 1.0]) in array (uses OpenMP).
   If using this on Data, then be sure to CHECK_EX_WRITE first
*/
void randomFillArray(long seed, double* array, size_t n);


void patternFillArray2D(size_t x, size_t y, double* array, size_t spacing,
                        size_t basex, size_t basey, size_t numpoints);

// Intended for debugging use only
void patternFillArray(int pattern, size_t x, size_t y, size_t z, double* array,
                      size_t spacing, size_t basex, size_t basey, size_t basez,
                      size_t numpoints);

}

#endif // __ESCRIPT_RANDOM_H__

