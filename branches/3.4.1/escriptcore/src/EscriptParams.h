
/*****************************************************************************
*
* Copyright (c) 2003-2013 by University of Queensland
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

#ifndef escript_EscriptParams_H
#define escript_EscriptParams_H
#include "system_dep.h"
#include <boost/python/list.hpp>
#include "Data.h"    // for the operators

namespace escript
{

class Data;
class DataLazy;

class EscriptParams
{
public:
  ESCRIPT_DLL_API
  EscriptParams();

  ESCRIPT_DLL_API
  int getInt(const char* name, int sentinel=0) const;
  
  ESCRIPT_DLL_API
  void setInt(const char* name, int value);

  boost::python::list
  listEscriptParams();

private:

  // If we get more params we can replace this with a map
	int too_many_lines;
	int autolazy;
	int too_many_levels;
// 	int too_many_nodes;
	int resolve_collective;
	int lazy_str_fmt;
	int lapack_support;
	int lazy_verbose;
	int amg_disabled;

protected: 
  // This is to provide fast access for methods in Data.
  // Its a little bit ugly, needing all those friends but I really want to
  // limit outside access to the char* interface

  int getTOO_MANY_LINES() {return too_many_lines;}
  int getAUTOLAZY() { return autolazy;}
  int getTOO_MANY_LEVELS() {return too_many_levels;}
//   int getTOO_MANY_NODES() {return too_many_nodes;}
  int getRESOLVE_COLLECTIVE() {return resolve_collective;}
  int getLAZY_STR_FMT() {return lazy_str_fmt;}
  int getLAZY_VERBOSE() {return lazy_verbose;}

  friend class escript::Data;
  friend class escript::DataLazy;
  friend Data operator+(const boost::python::api::object&, const escript::Data&);
  friend Data operator-(const boost::python::api::object&, const escript::Data&);
  friend Data operator*(const boost::python::api::object&, const escript::Data&);
  friend Data operator/(const boost::python::api::object&, const escript::Data&);
  friend Data operator+(const escript::Data&, const escript::Data&);
  friend Data operator-(const escript::Data&, const escript::Data&);
  friend Data operator*(const escript::Data&, const escript::Data&);
  friend Data operator/(const escript::Data&, const escript::Data&);
  friend Data operator+(const escript::Data&, const boost::python::api::object&);
  friend Data operator-(const escript::Data&, const boost::python::api::object&);
  friend Data operator*(const escript::Data&, const boost::python::api::object&);
  friend Data operator/(const escript::Data&, const boost::python::api::object&);
  friend Data C_GeneralTensorProduct(escript::Data& arg_0, escript::Data& arg_1,
                     int axis_offset, int transpose);
  friend Data condEval(escript::Data& mask, escript::Data& trueval, escript::Data& falseval);
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
  \return a list of tuples (parameter name, value, description)
*/
ESCRIPT_DLL_API
inline boost::python::list listEscriptParams()
{
   return escriptParams.listEscriptParams();
}



}
#endif
