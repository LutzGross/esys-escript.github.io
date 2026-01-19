
/*****************************************************************************
*
* Copyright (c) 2003-2020 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/

#include "EscriptParams.h"
#include "EsysException.h"
#include "EsysMPI.h"

#include <cmath> // to test if we know how to check for nan
#include <boost/python/tuple.hpp>

namespace bp = boost::python;

namespace escript {

ESCRIPT_DLL_API EscriptParams escriptParams; // externed in header file


EscriptParams::EscriptParams()
{
    // These #defs are for performance testing only
    // in general, I don't want people tweaking the
    // default value using compiler options
    // I've provided a python interface for that
#ifdef FAUTOLAZYON
    autoLazy = 1;
#else
    autoLazy = 0;
#endif
    lazyStrFmt = 0;
    lazyVerbose = 0;
#ifdef FRESCOLLECTON
    resolveCollective = 1;
#else
    resolveCollective = 0;
#endif
    tooManyLevels = 9;	// this is fairly arbitrary
    tooManyLines = 80;

    // now populate feature set
#ifdef ESYS_HAVE_CUDA
    features.insert("cuda");
#endif
#ifdef ESYS_HAVE_FINLEY
    features.insert("finley");
#endif
#ifdef ESYS_HAVE_LAPACK
    features.insert("lapack");
#endif
#ifdef ESYS_HAVE_MKL
    features.insert("mkl");
#endif
#ifdef ESYS_MPI
    features.insert("mpi");
#endif
        // Given c++11 we assume that this is
        // supported now
    features.insert("NAN_CHECK");
#ifdef ESYS_HAVE_HDF5
    features.insert("hdf5");
#endif
#ifdef _OPENMP
    features.insert("openmp");
#endif
#ifdef ESYS_HAVE_PASO
    features.insert("paso");
#endif
#ifdef ESYS_HAVE_RIPLEY
    features.insert("ripley");
#endif
#ifdef ESYS_HAVE_SILO
    features.insert("silo");
#endif
#ifdef ESYS_HAVE_SPECKLEY
    features.insert("speckley");
#endif
#ifdef ESYS_HAVE_TRILINOS
    features.insert("trilinos");
#endif
#ifdef ESYS_HAVE_UMFPACK
    features.insert("umfpack");
#endif
#ifdef ESYS_HAVE_MUMPS
    features.insert("mumps");
#endif
#ifdef ESYS_HAVE_WEIPA
    features.insert("weipa");
#endif
#ifdef ESYS_HAVE_BOOST_IO
    features.insert("unzip");
#endif
#ifdef ESYS_INDEXTYPE_LONG
    features.insert("longindex");
#endif

#ifdef ESYS_HAVE_BOOST_NUMPY
    features.insert("boostnumpy");
#endif

#ifndef ESYS_NO_SYMPY
    features.insert("sympy");
#endif
}

int EscriptParams::getInt(const std::string& name, int sentinel) const
{
    if (name == "AUTOLAZY")
        return autoLazy;
    else if (name == "LAZY_STR_FMT")
        return lazyStrFmt;
    else if (name == "LAZY_VERBOSE")
        return lazyVerbose;
    else if (name == "RESOLVE_COLLECTIVE")
        return resolveCollective;
    else if (name == "TOO_MANY_LEVELS")
        return tooManyLevels;
    else if (name == "TOO_MANY_LINES")
        return tooManyLines;

    return sentinel;
}

void
EscriptParams::setInt(const std::string& name, int value)
{
    if (name == "AUTOLAZY")
        autoLazy = value;
    else if (name == "LAZY_STR_FMT")
        lazyStrFmt = value;
    else if (name == "LAZY_VERBOSE")
        lazyVerbose = value;
    else if (name == "RESOLVE_COLLECTIVE")
        resolveCollective = value;
    else if (name == "TOO_MANY_LEVELS")
        tooManyLevels = value;
    else if (name == "TOO_MANY_LINES")
        tooManyLines = value;
    else
        throw ValueError("Invalid parameter name - "+name);
}

bp::list EscriptParams::listEscriptParams() const
{
   bp::list l;
   l.append(bp::make_tuple("AUTOLAZY", autoLazy, "{0,1} Operations involving Expanded Data will create lazy results."));
   l.append(bp::make_tuple("LAZY_STR_FMT", lazyStrFmt, "{0,1,2}(TESTING ONLY) change output format for lazy expressions."));
   l.append(bp::make_tuple("LAZY_VERBOSE", lazyVerbose, "{0,1} Print a warning when expressions are resolved because they are too large."));
   l.append(bp::make_tuple("RESOLVE_COLLECTIVE", resolveCollective, "(TESTING ONLY) {0.1} Collective operations will resolve their data."));
   l.append(bp::make_tuple("TOO_MANY_LEVELS", tooManyLevels, "(TESTING ONLY) maximum levels allowed in an expression."));
   l.append(bp::make_tuple("TOO_MANY_LINES", tooManyLines, "Maximum number of lines to output when printing data before printing a summary instead."));
   return l;
}

bool EscriptParams::hasFeature(const std::string& name) const
{
    if (name == "PASO_DIRECT") {
        // This is not in the constructor because escriptparams could be
        // constructed before main (and hence no opportunity to call INIT)
        // This use of MPI_COMM_WORLD is ok.
#ifdef ESYS_MPI
        int size;
        if (MPI_Comm_size(MPI_COMM_WORLD, &size) != MPI_SUCCESS || size > 1)
            return false;
#endif
        return hasFeature("paso") && (hasFeature("umfpack") || hasFeature("mkl") || hasFeature("mumps"));
    }

    return features.count(name) > 0;
}

bp::list EscriptParams::listFeatures() const
{
   bp::list l;
   FeatureSet::const_iterator it;
   for (it = features.begin(); it != features.end(); it++)
       l.append(*it);
   return l;
}

} // end namespace

