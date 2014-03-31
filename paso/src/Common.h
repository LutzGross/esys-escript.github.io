
/*****************************************************************************
*
* Copyright (c) 2003-2014 by University of Queensland
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


#ifndef INC_PASO_COMMON
#define INC_PASO_COMMON

/************************************************************************************/

/*    Finley finite element solver: common include file       */

/************************************************************************************/

/*   Copyrights by ACcESS Australia, 2003 */

/************************************************************************************/

#include "esysUtils/error.h"
#include "esysUtils/index.h"
#include "esysUtils/maths.h"
#include "esysUtils/mem.h"
#include "esysUtils/types.h"

#include <boost/enable_shared_from_this.hpp>
#include <boost/shared_ptr.hpp>

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
    extern "C"
    {
     #include <clapack.h>
    }
   #endif
   
#endif


typedef enum 
{
    PASO_AMG_UNDECIDED=-1,
    PASO_AMG_IN_F=0,
    PASO_AMG_IN_C=1
} AMGBlockSelect;

#endif /* #ifndef INC_PASO_COMMON */
