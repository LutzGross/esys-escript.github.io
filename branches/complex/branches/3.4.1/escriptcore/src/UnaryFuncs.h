
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


#if !defined escript_UnaryFuncs_20041124_H
#define escript_UnaryFuncs_20041124_H
#include "system_dep.h"

namespace escript {

#ifndef FP_NAN
#define FP_NAN IEEE_NaN()
#endif

#ifndef INFINITY
#define INFINITY IEEE_Infy()
#endif

//======================================================================

inline
double log1p (const double x)
{
  volatile double y;
  y = 1 + x;
  return log(y) - ((y-1)-x)/y ;
}

//======================================================================

inline
float IEEE_NaN()
{
   static unsigned char nan[4]={ 0, 0, 0xc0, 0x7f };
   return *( float *)nan;
}

//======================================================================

inline
double IEEE_Infy()
{
   static unsigned char infy[8]={ 0, 0, 0, 0, 0, 0, 0xf0, 0x7f };
   return *( double *)infy;
}


//======================================================================

#if defined (_WIN32) && !defined(__INTEL_COMPILER)
inline
double
acosh_substitute (const double x)
{
  if (x > 1.0 / SQRT_DBL_EPSILON)
    {
      return log (x) + M_LN2;
    }
  else if (x > 2)
    {
      return log (2 * x - 1 / (sqrt (x * x - 1) + x));
    }
  else if (x > 1)
    {
      double t = x - 1;
      return log1p (t + sqrt (2 * t + t * t));
    }
  else if (x == 1)
    {
      return 0;
    }
  else
    {
      return FP_NAN;
    }
}


//======================================================================

inline
double
asinh_substitute (const double x)
{
  double a = fabs (x);
  double s = (x < 0) ? -1 : 1;

  if (a > 1 / SQRT_DBL_EPSILON)
    {
      return s * (log (a) + M_LN2);
    }
  else if (a > 2)
    {
      return s * log (2 * a + 1 / (a + sqrt (a * a + 1)));
    }
  else if (a > SQRT_DBL_EPSILON)
    {
      double a2 = a * a;
      return s * log1p (a + a2 / (1 + sqrt (1 + a2)));
    }
  else
    {
      return x;
    }
}


//======================================================================

inline
double
atanh_substitute (const double x)
{
  double a = fabs (x);
  double s = (x < 0) ? -1 : 1;

  if (a > 1)
    {
      return FP_NAN;
    }
  else if (a == 1)
    {
      return (x < 0) ? -INFINITY : INFINITY;
    }
  else if (a >= 0.5)
    {
      return s * 0.5 * log1p (2 * a / (1 - a));
    }
  else if (a > DBL_EPSILON)
    {
      return s * 0.5 * log1p (2 * a + 2 * a * a / (1 - a));
    }
  else
    {
      return x;
    }
}
#endif  // windows substitutes for stupid microsoft compiler.


inline
double
fsign(double x)
{
  if (x == 0) {
    return 0;
  } else {
    return x/fabs(x);
  }
}

}

#endif
