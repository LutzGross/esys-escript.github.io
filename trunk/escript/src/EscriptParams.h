
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

#ifndef escript_EscriptParams_H
#define escript_EscriptParams_H
#include "system_dep.h"
#include <boost/python/list.hpp>

namespace escript
{

class EscriptParams
{
public:
  ESCRIPT_DLL_API
  EscriptParams();

  ESCRIPT_DLL_API
  int getInt(const char* name, int sentinel=0) const;
  
  ESCRIPT_DLL_API
  void setInt(const char* name, int value);
private:

  // If we get more params we can replace this with a map
	int too_many_lines;
};


extern EscriptParams escriptParams;

/**
  \brief Set the value of a named parameter.
  See listEscriptParams() (showEscriptParams() in python) for available parameters.
*/
ESCRIPT_DLL_API
void setEscriptParamInt(const char* name, int value);

/**
  \brief get the value of a named parameter.
  See listEscriptParams() (showEscriptParams() in python) for available parameters.
*/
ESCRIPT_DLL_API
int getEscriptParamInt(const char* name, int sentinel=0);

/**
  \brief describe available paramters.
  \return a list of tuples (parameter name, description)
*/
ESCRIPT_DLL_API
boost::python::list listEscriptParams();

}
#endif
