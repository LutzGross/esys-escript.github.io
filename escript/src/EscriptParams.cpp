

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

#include "EscriptParams.h"
#include <cstring>
#include <boost/python/tuple.hpp>

namespace escript
{

EscriptParams escriptParams;		// externed in header file


EscriptParams::EscriptParams()
{
   too_many_lines=80;
   autolazy=0;
   too_many_levels=70;
   too_many_nodes=15000;
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
   if (!strcmp(name,"TOO_MANY_NODES"))
   {
	return too_many_nodes;
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
   if (!strcmp(name,"TOO_MANY_NODES"))
   {
	too_many_nodes=value;
   }
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
listEscriptParams()
{
   using namespace boost::python;
   boost::python::list l;
   l.append(make_tuple("TOO_MANY_LINES","Maximum number of lines to output when printing data before printing a summary instead."));
   l.append(make_tuple("AUTOLAZY","{0,1} Operations involving Expanded Data will create lazy results."));
   return l;
}


}	// end namespace
