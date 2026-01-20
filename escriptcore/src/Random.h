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
**
*****************************************************************************/

#ifndef __ESCRIPT_RANDOM_H__
#define __ESCRIPT_RANDOM_H__

#include "system_dep.h"
#include "EsysMPI.h"

namespace escript
{
/* \brief put n random doubles (from [0.0, 1.0]) in array (uses OpenMP).
   If using this on Data, then be sure to CHECK_EX_WRITE first
*/
ESCRIPT_DLL_API
void randomFillArray(long seed, double* array, size_t n, JMPI mpiInfo);


ESCRIPT_DLL_API
void patternFillArray2D(size_t x, size_t y, double* array, size_t spacing,
                        size_t basex, size_t basey, size_t numpoints);

// Intended for debugging use only
void patternFillArray(int pattern, size_t x, size_t y, size_t z, double* array,
                      size_t spacing, size_t basex, size_t basey, size_t basez,
                      size_t numpoints);

}

#endif // __ESCRIPT_RANDOM_H__

