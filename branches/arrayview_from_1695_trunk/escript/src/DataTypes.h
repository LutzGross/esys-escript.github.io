
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

// This file contains static functions moved from DataArrayView
  /**
     \brief
     Calculate the number of values for the given shape.
  */
  ESCRIPT_DLL_API
  int
  noValues(const DataTypes::ShapeType& shape);

  /**
     \brief
     Calculate the number of values for the given region.
  */
  ESCRIPT_DLL_API
  int
  noValues(const DataTypes::RegionLoopRangeType& region);

  /**
     \brief
     Return the given shape as a string.

     \param shape - Input.
  */
  ESCRIPT_DLL_API
  std::string
  shapeToString(const DataTypes::ShapeType& shape);



 }   // End DataTypes


} // end of namespace

#endif
