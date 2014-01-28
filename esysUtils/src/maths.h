
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


#ifndef INC_ESYS_MATHS
#define INC_ESYS_MATHS

/************************************************************************************/

/*    Pull in a maths library and define ISNAN      */


/* some system values */
/* FIXME: This is not satisfactory.                                */
/* _ECC, __INTEL_COMPILER, and other                               */
/* intel compiler pre-defines need to be handled                   */
/* (__ICL, __ICC come to mind)                                     */
#if defined(_WIN32) && defined(__INTEL_COMPILER)
#include <mathimf.h>
#else
#include <cmath>
#endif

/*#ifndef NAN
   #define NAN (0.0/0.0)
#endif
*/
/*#define IS_NAN(__VAL__)  ( (__VAL__) == NAN )*/  /* this does not work */
/* #define IS_NAN(__VAL__)  ( ! ( ( (__VAL__) >= 0. ) ||  ( (__VAL__) <= 0. ) ) )  this does not work */

#ifdef isnan
  #define IS_NAN(__VAL__) (isnan(__VAL__))
#elif defined _isnan
  #define IS_NAN(__VAL__) (_isnan(__VAL__))
#else
/* If we do not have isnan then we can't reliably check for NaN - return false */
//  #define IS_NAN(__VAL__) (0)
  // This is not guaranteed to work if the optimiser thinks it can optimise this check away
  #define IS_NAN(__VAL__) ((__VAL__)!=(__VAL__))
#endif


#define EPSILON DBL_EPSILON
#define LARGE_POSITIVE_FLOAT DBL_MAX
#define SMALL_NEGATIVE_FLOAT -DBL_MAX

#endif 
