
/*****************************************************************************
*
* Copyright (c) 2003-2013 by University of Queensland
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


#ifndef INC_PASO_COMMON
#define INC_PASO_COMMON

/************************************************************************************/

/*    Finley finite element solver: common include file       */

/************************************************************************************/

/*   Copyrights by ACcESS Australia, 2003 */

/************************************************************************************/

#include "esysUtils/maths.h"
#include "esysUtils/mem.h"
#include "esysUtils/index.h"
#include "esysUtils/types.h"


#include <float.h>
#include <stdio.h>
#include <limits.h>
#include <stdlib.h>
#include <string.h>

#define LenString_MAX FILENAME_MAX*2
#define LenErrorMsg_MAX LenString_MAX

#ifdef USE_LAPACK

   #ifdef MKL_LAPACK
     #include <mkl_lapack.h>
   #else	/* assuming clapack */
     #include <clapack.h>
   #endif
   
#endif

#endif /* #ifndef INC_PASO_COMMON */
