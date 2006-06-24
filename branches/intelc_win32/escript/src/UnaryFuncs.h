// $Id$
/* 
 ************************************************************
 *          Copyright 2006 by ACcESS MNRF                   *
 *                                                          *
 *              http://www.access.edu.au                    *
 *       Primary Business: Queensland, Australia            *
 *  Licensed under the Open Software License version 3.0    *
 *     http://www.opensource.org/licenses/osl-3.0.php       *
 *                                                          *
 ************************************************************
*/
                                                                           
#if !defined escript_UnaryFuncs_20041124_H
#define escript_UnaryFuncs_20041124_H
#include "system_dep.h"

namespace escript {

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

} // end of namespace
#endif
