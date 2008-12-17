
/*******************************************************
*
* Copyright (c) 2003-2008 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
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
     print a message about how many MPI CPUs and OpenMP threads we're using
  */
  ESCRIPT_DLL_API void printParallelThreadCnt();

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

  /**
     \brief
     returns the total number of available MPI processes for MPI_COMM_WORLD
  */
  ESCRIPT_DLL_API int getMPISizeWorld();

  /**
     \brief
     returns the MPI processor number within MPI_COMM_WORLD
  */
  ESCRIPT_DLL_API int getMPIRankWorld();

} // end of namespace
#endif
