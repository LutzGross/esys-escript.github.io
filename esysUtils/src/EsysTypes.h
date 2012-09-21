
/*****************************************************************************
*
* Copyright (c) 2003-2012 by University of Queensland
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


#ifndef ESYSTYPES_H
#define ESYSTYPES_H
#include "system_dep.h"
/*
 * Use the integer types defined in the 1999 ISO C Standard
 * To specify a suitable Esys integer type
 */
#include <stdint.h>

#if ESYS_INT_BITS==64
typedef int64_t EsysIntType;
#else
typedef int32_t EsysIntType;
#endif

/*
 * A primitive test to ensure the array index type is at least as large
 * as requested. Could put in another test if it is larger.
 * An obscure compile error will result if the array index type isn't large
 * enough 
 */
static char EsysIntType_Too_Small[sizeof(EsysIntType)*8-ESYS_INT_BITS];

#endif
