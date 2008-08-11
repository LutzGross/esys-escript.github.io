
/* $Id$ */

/*******************************************************
 *
 *           Copyright 2003-2008 by ACceSS MNRF
 *       Copyright 2008 by University of Queensland
 *
 *                http://esscc.uq.edu.au
 *        Primary Business: Queensland, Australia
 *  Licensed under the Open Software License version 3.0
 *     http://www.opensource.org/licenses/osl-3.0.php
 *
 *******************************************************/

#if !defined escript_DataTypes_20080811_H
#define escript_DataTypes_20080811_H
#include "system_dep.h"
#include "DataVector.h"
#include <vector>

namespace escript {

 namespace DataTypes {
  //
  // Some basic types which define the data values and view shapes.
  typedef DataVector                        ValueType;
  typedef std::vector<int>                  ShapeType;
  typedef std::vector<std::pair<int, int> > RegionType;
  typedef std::vector<std::pair<int, int> > RegionLoopRangeType;
  static const int maxRank=4;

 }   // End DataTypes


} // end of namespace

#endif
