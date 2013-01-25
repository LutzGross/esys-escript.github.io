

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

#include <cstring>
#include "EscriptParams.h"

namespace escript
{

EscriptParams escriptParams;		// externed in header file


EscriptParams::EscriptParams()
{
   too_many_lines=80;
}

int 
EscriptParams::getInt(const char* name, int sentinel) const
{
   if (!strcmp(name,"TOO_MANY_LINES"))
   {
	return too_many_lines;
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


}	// end namespace
