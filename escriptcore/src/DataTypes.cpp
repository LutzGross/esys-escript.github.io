
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* https://github.com/LutzGross/esys-escript.github.io
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/

#include "DataTypes.h"
#include "DataException.h"

#include <sstream>
#include <boost/python/extract.hpp>
#include <boost/python/tuple.hpp>

namespace bp = boost::python;

namespace {

using namespace escript;
using namespace escript::DataTypes;

/*
  \brief
  Calculate the slice range from the given python key object
  Used by getSliceRegion - since it is not used anywhere else, I have elected not to declare it
  in the header.
  Returns the python slice object key as a pair of ints where the first 
  member is the start and the second member is the end. the presence of a possible
  step attribute with value different from one will throw an exception

  /param key - Input - key object specifying slice range.
*/
   std::pair<int,int>
   getSliceRange(const bp::object& key, int shape)
   {
      /* default slice range is range of entire shape dimension */
      int s0=0, s1=shape;;
      bp::extract<int> slice_int(key);
      if (slice_int.check()) {
         /* if the key is a single int set start=key and end=key */
         /* in this case, we want to return a rank-1 dimension object from
         this object, taken from a single index value for one of this
         object's dimensions */
         s0=slice_int();
         s1=s0;
      } else {
         /* if key is a pair extract begin and end values */
         bp::extract<int> step(key.attr("step"));
         if (step.check() && step()!=1) {
            throw DataException("Error - Data does not support increments in slicing ");
         } else {
            bp::extract<int> start(key.attr("start"));
            if (start.check()) {
               s0=start();
            }
            bp::extract<int> stop(key.attr("stop"));
            if (stop.check()) {
               s1=stop();
            }
         }
      }
      if (s0 < 0) 
         throw DataException("Error - slice index out of range.");
      if (s0 == s1 && s1 >= shape)
         throw DataException("Error - slice index out of range.");
      if (s0 != s1 &&  s1>shape)
         throw DataException("Error - slice index out of range.");
      if (s0 > s1) 
         throw DataException("Error - lower index must less or equal upper index.");
      return std::pair<int,int>(s0,s1);
   }
} // anonymous namespace


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





   DataTypes::RegionType
   getSliceRegion(const DataTypes::ShapeType& shape, const bp::object& key)
   {
      int slice_rank, i;
      int this_rank=shape.size();
      RegionType out(this_rank);
      /* allow for case where key is singular eg: [1], this implies we
      want to generate a rank-1 dimension object, as opposed to eg: [1,2]
      which implies we want to take a rank dimensional object with one
      dimension of size 1 */
      bp::extract<bp::tuple> key_tuple(key);
      if (key_tuple.check()) {
         slice_rank=bp::extract<int> (key.attr("__len__")());
         /* ensure slice is correctly dimensioned */
         if (slice_rank>this_rank) {
            throw DataException("Error - rank of slices does not match rank of slicee");
         } else {
            /* calculate values for slice range */
            for (i=0;i<slice_rank;i++) {
               out[i]=getSliceRange(key[i],shape[i]);
            }
         }
      } else {
         slice_rank=1;
         if (slice_rank>this_rank) {
            throw DataException("Error - rank of slices does not match rank of slicee");
         } else {
            out[0]=getSliceRange(key,shape[0]);
         }
      }
      for (i=slice_rank;i<this_rank;i++) {
         out[i]=std::pair<int,int>(0,shape[i]);
      }
      return out;
   }

   DataTypes::ShapeType
   getResultSliceShape(const RegionType& region)
   {
      int dimSize;
      ShapeType result;
      RegionType::const_iterator i;
      for (i=region.begin();i!=region.end();i++) {
         dimSize=((i->second)-(i->first));
         if (dimSize!=0) {
            result.push_back(dimSize);
         }
      }
      return result;
   }

   DataTypes::RegionLoopRangeType
   getSliceRegionLoopRange(const DataTypes::RegionType& region) 
   {
      DataTypes::RegionLoopRangeType region_loop_range(region.size());
      unsigned int i;
      for (i=0;i<region.size();i++) {
         if (region[i].first==region[i].second) {
            region_loop_range[i].first=region[i].first;
            region_loop_range[i].second=region[i].second+1;
         } else {
            region_loop_range[i].first=region[i].first;
            region_loop_range[i].second=region[i].second;
         }
      }
      return region_loop_range;
   }


   std::string 
   createShapeErrorMessage(const std::string& messagePrefix,
                                          const DataTypes::ShapeType& other,
					  const DataTypes::ShapeType& thisShape)
   {
      std::stringstream temp;
      temp << messagePrefix
           << " This shape: " << shapeToString(thisShape)
           << " Other shape: " << shapeToString(other);
      return temp.str();
   }

}	// end namespace DataTypes
}	// end namespace escript

