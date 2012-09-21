
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
     return the SVN version number used to build this version.
     \warning Only gives accurate answers for clean checkouts
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
     \warning Use of this method is strongly discouraged. It may be deprecated in future.
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
  /**
     \brief
     returns the maximum value of an integer over all processors within MPI_COMM_WORLD
  */
  ESCRIPT_DLL_API int getMPIWorldMax(const int val);

  /**
     \brief returns sum of an integer over all processors with MPI_COMM_WORLD
  */
  ESCRIPT_DLL_API int getMPIWorldSum(const int val);

  /**
     \brief performs a barrier synchronization across all processors.
  */
  ESCRIPT_DLL_API void MPIBarrierWorld();

  /**
     \brief uses MPI_Comm_spawn to run an external MPI program safely.
  */
  ESCRIPT_DLL_API int runMPIProgram(const boost::python::list args);

 /**
    \brief 
    returns machine precision
 */
 ESCRIPT_DLL_API double getMachinePrecision();
 /*
   \brief
   return largest positive float
 */
 ESCRIPT_DLL_API double getMaxFloat();

 ESCRIPT_DLL_API
 void
 saveDataCSV(const std::string& filename, boost::python::dict arg, const std::string& sep, const std::string& csep,
bool append=false); 


 /**
 Resolve a collection of Data objects together now.
 To get performance improvements, the objects will need to have the same function space and share Dag components.
 \param obj A python list or tuple of Data objects to be resolved together.
 
 */
 ESCRIPT_DLL_API
 void 
 resolveGroup(boost::python::object obj);

} // end of namespace
#endif
