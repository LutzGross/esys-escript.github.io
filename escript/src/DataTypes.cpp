/* $Id:$ */

/*******************************************************
*
*           Copyright 2003-2007 by ACceSS MNRF
*       Copyright 2007 by University of Queensland
*
*                http://esscc.uq.edu.au
*        Primary Business: Queensland, Australia
*  Licensed under the Open Software License version 3.0
*     http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/

#include <sstream>
#include "DataTypes.h"

namespace escript
{
namespace DataTypes
{

   int
   noValues(const ShapeType& shape) 
   {
      ShapeType::const_iterator i;
      //
      // An empty shape vector means rank 0 which contains 1 value
      int noValues=1;
      for (i=shape.begin();i!=shape.end();i++) {
         noValues*=(*i);
      }
      return noValues;
   }

   int
   noValues(const RegionLoopRangeType& region) 
   {
      //
      // An empty region vector means rank 0 which contains 1 value
      int noValues=1;
      unsigned int i;
      for (i=0;i<region.size();i++) {
         noValues*=region[i].second-region[i].first;
      }
      return noValues;
   }

   std::string
   shapeToString(const DataTypes::ShapeType& shape)
   {
      std::stringstream temp;
      temp << "(";
      unsigned int i;
      for (i=0;i<shape.size();i++) {
         temp << shape[i];
         if (i < shape.size()-1) {
            temp << ",";
         }
      }
      temp << ")";
      return temp.str();
   }

}	// end namespace DataTypes
}	// end namespace escript