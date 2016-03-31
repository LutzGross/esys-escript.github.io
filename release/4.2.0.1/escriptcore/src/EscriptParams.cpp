

/*****************************************************************************
*
* Copyright (c) 2003-2016 by The University of Queensland
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

#define ESNEEDPYTHON
#include "esysUtils/first.h"


#include "EscriptParams.h"
#include <cstring>
#include <boost/python/tuple.hpp>
#include <cmath>                        // to test if we know how to check for nan

#include "esysUtils/EsysException.h"
#include "esysUtils/Esys_MPI.h"

namespace escript
{

EscriptParams escriptParams;                // externed in header file


EscriptParams::EscriptParams()
{
   too_many_lines=80;
   autolazy=0;
   too_many_levels=70;
   lazy_str_fmt=0;
   lazy_verbose=0;
#ifdef USE_NETCDF
   has_netcdf=1;
#else   
   has_netcdf=0;
#endif   
#ifdef USE_LAPACK
   lapack_support=1;
#else
   lapack_support=0;
#endif

    gmsh = gmsh_mpi = 0;
#if defined(GMSH) || defined(GMSH_MPI)
    gmsh = 1;
#endif
    //only mark gmsh as mpi if escript built with mpi, otherwise comm_spawns
    //might just fail terribly
#if defined(GMSH_MPI) && defined(ESYS_MPI)
    gmsh_mpi = 1;
#endif

#ifdef ESYS_MPI
    amg_disabled=true;
#else
    amg_disabled=false;
#endif

    temp_direct_solver=false;   // This variable is to be removed once proper
                                // SolverOptions support is in place
#ifdef MKL
    temp_direct_solver=true;
#endif
#ifdef USE_UMFPACK
    temp_direct_solver=true;
#endif
#ifdef PASTIX
    temp_direct_solver=true;
#endif

                        // These #defs are for performance testing only
                        // in general, I don't want people tweaking the
                        // default value using compiler options
                        // I've provided a python interface for that
#ifdef FAUTOLAZYON
   autolazy=1;
#endif
#ifdef FAUTOLAZYOFF
   autolazy=0;
#endif

#ifdef FRESCOLLECTON
   resolve_collective=1;
#endif
#ifdef FRESCOLLECTOFF
   resolve_collective=0;
#endif
}

int 
EscriptParams::getInt(const char* name, int sentinel) const
{
   if (!strcmp(name,"TOO_MANY_LINES"))
   {
        return too_many_lines;
   }
   if (!strcmp(name,"AUTOLAZY"))
   {
        return autolazy;
   }
   if (!strcmp(name,"TOO_MANY_LEVELS"))
   {
        return too_many_levels;
   }
   if (!strcmp(name,"RESOLVE_COLLECTIVE"))
   {
        return resolve_collective;
   }
   if (!strcmp(name,"LAZY_STR_FMT"))
   {
        return lazy_str_fmt;
   }
   if (!strcmp(name,"LAPACK_SUPPORT"))
   {
        return lapack_support;
   }
   if (!strcmp(name, "NAN_CHECK"))
   {
#ifdef isnan        
        return 1;
#else
        return 0;
#endif
   }
   if (!strcmp(name,"LAZY_VERBOSE"))
   {
        return lazy_verbose;
   }
   if (!strcmp(name, "DISABLE_AMG"))
   {
        return amg_disabled;
   }
   if (!strcmp(name, "MPIBUILD"))
   {
#ifdef ESYS_MPI           
        return 1;
#else
        return 0;
#endif
   }
   if (!strcmp(name, "PASO_DIRECT"))
   {
        // This is not in the constructor because escriptparams could be constructed 
        // before main (and hence no opportunity to call INIT)
        #ifdef ESYS_MPI
            int size;
            if (MPI_Comm_size(MPI_COMM_WORLD, &size)!=MPI_SUCCESS)        // This would break in a subworld
            {
                temp_direct_solver=false;        
            }
            if (size>1)
            {
                temp_direct_solver=false;
            }
        #endif   
        return temp_direct_solver;
   }
    if (!strcmp(name, "NETCDF_BUILD"))
    {
       return has_netcdf; 
    }
    if (!strcmp(name, "GMSH_SUPPORT"))
        return gmsh;
    if (!strcmp(name, "GMSH_MPI"))
        return gmsh_mpi;
   return sentinel;
}
  
void 
EscriptParams::setInt(const char* name, int value)
{
   // Note: there is no way to modify the LAPACK_SUPPORT variable ATM
    if (!strcmp(name,"TOO_MANY_LINES"))
        too_many_lines=value;
    else if (!strcmp(name,"AUTOLAZY"))
        autolazy=!(value==0);        // set to 1 or zero
    else if (!strcmp(name,"TOO_MANY_LEVELS"))
        too_many_levels=value;
    else if (!strcmp(name,"RESOLVE_COLLECTIVE"))
        resolve_collective=value;
    else if (!strcmp(name,"LAZY_STR_FMT"))
        lazy_str_fmt=value;
    else if (!strcmp(name,"LAZY_VERBOSE"))
        lazy_verbose=value;
    else
       throw esysUtils::EsysException("Invalid parameter name");
}

void 
setEscriptParamInt(const char* name, int value)
{
   escriptParams.setInt(name,value);
}


int
getEscriptParamInt(const char* name, int sentinel)
{
   return escriptParams.getInt(name, sentinel);
}

boost::python::list
EscriptParams::listEscriptParams()
{
   using namespace boost::python;
   boost::python::list l;
   l.append(make_tuple("TOO_MANY_LINES", too_many_lines, "Maximum number of lines to output when printing data before printing a summary instead."));
   l.append(make_tuple("AUTOLAZY", autolazy, "{0,1} Operations involving Expanded Data will create lazy results."));
   l.append(make_tuple("RESOLVE_COLLECTIVE",resolve_collective ,"(TESTING ONLY) {0.1} Collective operations will resolve their data."));
   l.append(make_tuple("TOO_MANY_LEVELS", too_many_levels, "(TESTING ONLY) maximum levels allowed in an expression."));
   l.append(make_tuple("LAZY_STR_FMT", lazy_str_fmt, "{0,1,2}(TESTING ONLY) change output format for lazy expressions."));
   l.append(make_tuple("LAZY_VERBOSE", lazy_verbose, "{0,1} Print a warning when expressions are resolved because they are too large."));
   l.append(make_tuple("DISABLE_AMG", amg_disabled, "{0,1} AMG is disabled."));
   l.append(make_tuple("NETCDF_BUILD", has_netcdf, "{0,1} Was this build made with netcdf libraries?"));
   l.append(make_tuple("GMSH_SUPPORT", gmsh, "{0,1} Non-python GMSH support is available."));
   l.append(make_tuple("GMSH_MPI", gmsh_mpi, "{0,1} Both GMSH and escript have MPI capabilities."));
   return l;
}




}        // end namespace
