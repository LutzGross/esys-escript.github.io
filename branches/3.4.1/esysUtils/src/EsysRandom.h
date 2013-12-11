/*****************************************************************************
*
* Copyright (c) 2013 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development since 2012 by School of Earth Sciences
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
}
#endif
