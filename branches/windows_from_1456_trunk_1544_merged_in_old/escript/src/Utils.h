
/* $Id$ */

/*******************************************************
 *
 *           Copyright 2003-2007 by ACceSS MNRF
 *       Copyright 2007 by University of Queensland
 *
 *                http://esscc.uq.edu.au
 *        Primary Business: Queensland, Australia
 *  Licensed under the Open Software License version 3.0
 *     http://www.opensource.org/licenses/osl-3.0.php
 *
 *******************************************************/

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
     return the latest SVN version number
  */
  ESCRIPT_DLL_API int getSvnVersion();

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
