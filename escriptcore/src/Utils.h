
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/

#ifndef __ESCRIPT_UTILS_H__
#define __ESCRIPT_UTILS_H__

#include "system_dep.h"
#include <boost/python/dict.hpp>

#ifdef ESYS_HAVE_BOOST_NUMPY
#include <boost/python/numpy.hpp>
#endif

#include <string>

#ifdef ESYS_HAVE_NETCDF4
// Forward declaration
namespace netCDF {
    class NcFile;
}
#endif

namespace escript {

// Forward declaration
class FunctionSpace;

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
    return the GIT hash used to build this version.
*/
ESCRIPT_DLL_API std::string getGitVersion();

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
    returns the number of threads
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
    returns the machine precision
*/
ESCRIPT_DLL_API double getMachinePrecision();

/*
    \brief
    return largest positive float
*/
ESCRIPT_DLL_API double getMaxFloat();

ESCRIPT_DLL_API void saveDataCSV(const std::string& filename,
                                 boost::python::dict arg,
                                 const std::string& sep,
                                 const std::string& csep,
                                 bool refid=false,
                                 bool append=false);

#ifdef ESYS_HAVE_BOOST_NUMPY
ESCRIPT_DLL_API boost::python::list getNumpy(boost::python::dict arg);
#else
ESCRIPT_DLL_API void getNumpy(boost::python::dict arg);
#endif

#ifdef ESYS_HAVE_BOOST_NUMPY
ESCRIPT_DLL_API boost::python::numpy::ndarray convertToNumpy(escript::Data data);
#else
ESCRIPT_DLL_API void convertToNumpy(escript::Data data);
#endif

// #ifdef ESYS_HAVE_BOOST_NUMPY
// void initBoostNumpy();
// #endif

#ifdef ESYS_HAVE_BOOST_NUMPY
escript::Data numpyToData(boost::python::numpy::ndarray& array, bool isComplex, FunctionSpace& functionspace);
#endif


/**
    \brief
    Resolve a collection of Data objects together now.
    To get performance improvements, the objects will need to have the same
    function space and share Dag components.
    \param obj A python list or tuple of Data objects to be resolved together.
*/
ESCRIPT_DLL_API void resolveGroup(boost::python::object obj);

#ifdef ESYS_HAVE_NETCDF4
/**
    \brief
    Opens a NetCDF file for reading
    \param file NcFile object to be initialized
    \param filename Path to the NetCDF file
    \return true if successful, false otherwise
*/
ESCRIPT_DLL_API bool openNcFile(netCDF::NcFile& file, const std::string& filename);
#endif

} // end of namespace

#endif // __ESCRIPT_UTILS_H__
