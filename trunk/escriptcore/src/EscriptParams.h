
/*****************************************************************************
*
* Copyright (c) 2003-2016 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/

#ifndef __ESCRIPT_PARAMS_H__
#define __ESCRIPT_PARAMS_H__

#include "system_dep.h"
#include "Data.h"    // for the operators

#include <boost/python/list.hpp>

namespace escript
{

class Data;
class DataLazy;

class EscriptParams
{
public:
  EscriptParams();

  int getInt(const char* name, int sentinel=0) const;
  
  void setInt(const char* name, int value);

  boost::python::list
  listEscriptParams();

private:

  // If we get more params we can replace this with a map
    int too_many_lines;
    int autolazy;
    int too_many_levels;
    int resolve_collective;
    int lazy_str_fmt;
    int lapack_support;
    int lazy_verbose;
    int amg_disabled;
    int has_netcdf;
    int have_trilinos;
    int have_unzip;
    int gmsh;
    int gmsh_mpi;
    mutable int temp_direct_solver;

protected: 
  // This is to provide fast access for methods in Data.
  // Its a little bit ugly, needing all those friends but I really want to
  // limit outside access to the char* interface

  int getTOO_MANY_LINES() {return too_many_lines;}
  int getAUTOLAZY() { return autolazy;}
  int getTOO_MANY_LEVELS() {return too_many_levels;}
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
void setEscriptParamInt(const char* name, int value);

/**
  \brief get the value of a named parameter.
  See listEscriptParams() (showEscriptParams() in python) for available parameters.
*/
int getEscriptParamInt(const char* name, int sentinel=0);

/**
  \brief describe available paramters.
  \return a list of tuples (parameter name, value, description)
*/
inline boost::python::list listEscriptParams()
{
   return escriptParams.listEscriptParams();
}

} // namespace escript

#endif // __ESCRIPT_PARAMS_H__

