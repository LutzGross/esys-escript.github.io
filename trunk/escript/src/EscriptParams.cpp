

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

#include "EscriptParams.h"
#include <cstring>
#include <boost/python/tuple.hpp>
#include <cmath>			// to test if we know how to check for nan

namespace escript
{

EscriptParams escriptParams;		// externed in header file


EscriptParams::EscriptParams()
{
   too_many_lines=80;
   autolazy=0;
   too_many_levels=70;
//    too_many_nodes=15000;
   resolve_collective=0;
   lazy_str_fmt=0;
   lazy_verbose=0;

#ifdef USE_LAPACK
   lapack_support=1;
#else
   lapack_support=0;
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


#ifdef ESYS_MPI
	amg_disabled=true;
#else
	amg_disabled=false;
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
//    if (!strcmp(name,"TOO_MANY_NODES"))
//    {
// 	return too_many_nodes;
//    }
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
   return sentinel;
}
  
void 
EscriptParams::setInt(const char* name, int value)
{
   if (!strcmp(name,"TOO_MANY_LINES"))
   {
	too_many_lines=value;
   }
   if (!strcmp(name,"AUTOLAZY"))
   {
	autolazy=!(value==0);	// set to 1 or zero
   }
   if (!strcmp(name,"TOO_MANY_LEVELS"))
   {
	too_many_levels=value;
   }
//    if (!strcmp(name,"TOO_MANY_NODES"))
//    {
// 	too_many_nodes=value;
//    }
   if (!strcmp(name,"RESOLVE_COLLECTIVE"))
   {
	resolve_collective=value;
   }
   if (!strcmp(name,"LAZY_STR_FMT"))
   {
	lazy_str_fmt=value;
   }
   if (!strcmp(name,"LAZY_VERBOSE"))
   {
	lazy_verbose=value;
   }
   // Note: there is no way to modifiy the LAPACK_SUPPORT variable atm
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
//    l.append(make_tuple("TOO_MANY_NODES", too_many_nodes, "(TESTING ONLY) maximum number of nodes in a expression."));
   l.append(make_tuple("LAZY_STR_FMT", lazy_str_fmt, "{0,1,2}(TESTING ONLY) change output format for lazy expressions."));
   l.append(make_tuple("LAZY_VERBOSE", lazy_verbose, "{0,1} Print a warning when expressions are resolved because they are too large."));
   l.append(make_tuple("DISABLE_AMG", amg_disabled, "{0,1} AMG is disabled."));
   return l;
}




}	// end namespace
