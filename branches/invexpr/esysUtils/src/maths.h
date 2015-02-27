
/*****************************************************************************
*
* Copyright (c) 2003-2015 by University of Queensland
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


#ifndef INC_ESYS_MATHS
#define INC_ESYS_MATHS

/************************************************************************************/

/*    Pull in a maths library and define ISNAN      */
#include <cmath>

#ifdef isnan
  #define IS_NAN(__VAL__) (isnan(__VAL__))
#elif defined _isnan
  #define IS_NAN(__VAL__) (_isnan(__VAL__))
#else
  // This is not guaranteed to work if the optimiser thinks it can optimise this check away
  #define IS_NAN(__VAL__) (!((__VAL__)==(__VAL__)))
#endif


#define EPSILON DBL_EPSILON
#define LARGE_POSITIVE_FLOAT DBL_MAX
#define SMALL_NEGATIVE_FLOAT -DBL_MAX

#endif 
