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

#if !defined  escript_Utils_H
#define escript_Utils_H
#include "system_dep.h"

namespace escript {

  /**
     \brief
     some functions

  */

  /**
     \brief
     set the number of threads 
  */
  ESCRIPT_DLL_API void setNumberOfThreads(const int num_threads);

  /**
     \brief
     returns  the number of threads 
  */
  ESCRIPT_DLL_API int getNumberOfThreads();

} // end of namespace
#endif
