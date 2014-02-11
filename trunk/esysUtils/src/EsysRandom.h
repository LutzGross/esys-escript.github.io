/*****************************************************************************
*
* Copyright (c) 2013-2014 by University of Queensland
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

#ifndef ESYS_RANDOM_H
#define ESYS_RANDOM_H
namespace esysUtils
{
/* \brief put n random doubles in array (uses OpenMP).
   If using this on Data, then be sure to CHECK_EX_WRITE first
*/
void randomFillArray(long seed, double* array, size_t n);

/* Intended for debugging use only */
void patternFillArray(int pattern, size_t x, size_t y, size_t z, double* array, size_t spacing, size_t basex, size_t basey, size_t basez);

}
#endif
