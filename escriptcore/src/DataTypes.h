
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

#ifndef __ESCRIPT_DATATYPES_H__
#define __ESCRIPT_DATATYPES_H__

#include <boost/python/object_fwd.hpp>

#include "system_dep.h"
#include "Assert.h"

#include <complex>
#include <limits>
#include <string>
#include <vector>

namespace escript {

namespace DataTypes {

/**
\namespace escript::DataTypes
\brief Contains the types to represent Shapes, Regions, RegionLoop ranges and
       vectors of data as well as the functions to manipulate them.
\note The contents of the namespace are spread between DataTypes.h and DataVector.h
*/

  //
  // Some basic types which define the data values and view shapes.
  typedef std::vector<int>                  ShapeType;//!< The shape of a single datapoint.
  typedef std::vector<std::pair<int, int> > RegionType;
  typedef std::vector<std::pair<int, int> > RegionLoopRangeType;
  static const int maxRank=4;//!< The maximum number of dimensions a datapoint can have.
  static const ShapeType scalarShape;//!< Use this instead of creating empty shape objects for scalars.
  typedef long vec_size_type;

  /// type of all real-valued scalars in escript
  typedef double real_t;

  /// complex data type
  typedef std::complex<real_t> cplx_t;

  /// type for array/matrix indices used both globally and on each rank
#ifdef ESYS_INDEXTYPE_LONG
  typedef long index_t;
#else
  typedef int index_t;
#endif

  typedef std::vector<index_t> IndexVector;

  typedef index_t dim_t;

  /**
     \brief
     Returns the minimum finite value for the index_t type.
  */
  inline index_t index_t_min()
  {
      return std::numeric_limits<index_t>::min();
  }

  /**
     \brief
     Returns the maximum finite value for the index_t type.
  */
  inline index_t index_t_max()
  {
      return std::numeric_limits<index_t>::max();
  }

  /**
     \brief
     Returns the maximum finite value for the real_t type.
  */
  inline real_t real_t_max()
  {
      return std::numeric_limits<real_t>::max();
  }

  /**
     \brief
     Returns the machine epsilon for the real_t type.
  */
  inline real_t real_t_eps()
  {
      return std::numeric_limits<real_t>::epsilon();
  }

  /**
     \brief
     Calculate the number of values in a datapoint with the given shape.
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

  /**
     \brief
     Determine the shape of the specified slice region.

     \param region - Input - Slice region
  */
  ESCRIPT_DLL_API
  DataTypes::ShapeType
  getResultSliceShape(const DataTypes::RegionType& region);


 /**
     \brief
     Determine the region specified by the given python slice object.

     \param shape - Input - Shape of the object being sliced.
     \param key - Input -
                    python slice object specifying region to be returned.

     The slice object is a tuple of n python slice specifiers, where
     n <= the rank of this Data object. Each slice specifier specifies the
     range of indexes to be sliced from the corresponding dimension. The
     first specifier corresponds to the first dimension, the second to the
     second and so on. Where n < the rank, the remaining dimensions are
     sliced across the full range of their indices.

     Each slice specifier is of the form "a:b", which specifies a slice
     from index a, up to but not including index b. Where index a is omitted
     a is assumed to be 0. Where index b is omitted, b is assumed to be the
     length of this dimension. Where both are omitted (eg: ":") the slice is
     assumed to encompass that entire dimension.

     Where one of the slice specifiers is a single integer, eg: [1], we
     want to generate a rank-1 dimension object, as opposed to eg: [1,2]
     which implies we want to take a rank dimensional object with one
     dimension of size 1.

     The return value is a vector of pairs with length equal to the rank of
     this object. Each pair corresponds to the range of indices from the
     corresponding dimension to be sliced from, as specified in the input
     slice object.

     Examples:

       For a rank 1 object of shape(5):

         getSliceRegion(:)   => < <0,5> >
         getSliceRegion(2:3) => < <2,3> >
         getSliceRegion(:3)  => < <0,3> >
         getSliceRegion(2:)  => < <2,5> >

       For a rank 2 object of shape(4,5):

         getSliceRegion(2:3) => < <2,3> <0,5> >
         getSliceRegion(2)   => < <2,3> <0,5> >
           NB: but return object requested will have rank 1, shape(5), with
               values taken from index 2 of this object's first dimension.

       For a rank 3 object of shape (2,4,6):

         getSliceRegion(0:2,0:4,0:6) => < <0,2> <0,4> <0,6> >
         getSliceRegion(:,:,:)       => < <0,2> <0,4> <0,6> >
         getSliceRegion(0:1)         => < <0,1> <0,4> <0,6> >
         getSliceRegion(:1,0:2)      => < <0,1> <0,2> <0,6> >


     Note: Not unit tested in c++.
  */
   DataTypes::RegionType
   getSliceRegion(const DataTypes::ShapeType& shape, const boost::python::object& key);

  /**
  \brief
   Modify region to copy from in order to
   deal with the case where one range in the region contains identical indexes,
   eg: <<1,1><0,3><0,3>>
   This situation implies we want to copy from an object with rank greater than that of this
   object. eg: we want to copy the values from a two dimensional slice out of a three
   dimensional object into a two dimensional object.
   We do this by taking a slice from the other object where one dimension of
   the slice region is of size 1. So in the above example, we modify the above
   region like so: <<1,2><0,3><0,3>> and take this slice.
  */
  ESCRIPT_DLL_API
  DataTypes::RegionLoopRangeType
  getSliceRegionLoopRange(const DataTypes::RegionType& region);

  /**
   \brief Return the rank (number of dimensions) of the given shape.

   \param shape
   \return the rank.
  */
  inline
  int
  getRank(const DataTypes::ShapeType& shape)
  {
	return shape.size();
  }


  /**
  \brief Compute the offset (in 1D vector) of a given subscript with a shape.

  \param shape - Input - Shape of the datapoint.
  \param i - Input - subscript to locate.
  \return offset relative to the beginning of the datapoint.
  */
  inline
  vec_size_type
  getRelIndex(const DataTypes::ShapeType& shape, vec_size_type i)
  {
  	ESYS_ASSERT(getRank(shape)==1, "Incorrect number of indices for the rank of this object.");
	ESYS_ASSERT(i < DataTypes::noValues(shape), "Invalid index.");
	return i;
  }

  /**
  \brief Compute the offset (in 1D vector) of a given subscript with a shape.

  \param shape - Input - Shape of the datapoint.
  \param i - Input - row
  \param j - Input - column
  \return offset relative to the beginning of the datapoint.
  */
  inline
  vec_size_type
  getRelIndex(const DataTypes::ShapeType& shape, vec_size_type i,
	   vec_size_type j)
  {
	// Warning: This is not C ordering. Do not try to figure out the params by looking at the code
  	ESYS_ASSERT(getRank(shape)==2, "Incorrect number of indices for the rank of this object.");
  	vec_size_type temp=i+j*shape[0];
  	ESYS_ASSERT(temp < DataTypes::noValues(shape), "Invalid index.");
	return temp;
  }

  /**
  \brief Compute the offset (in 1D vector) of a given subscript with a shape.

  \param shape - Input - Shape of the datapoint.
  \param i,j,k - Input - subscripts to locate.
  \return offset relative to the beginning of the datapoint.
  */
  inline
  vec_size_type
  getRelIndex(const DataTypes::ShapeType& shape, vec_size_type i,
	   vec_size_type j, vec_size_type k)
  {
	// Warning: This is not C ordering. Do not try to figure out the params by looking at the code
  	ESYS_ASSERT(getRank(shape)==3, "Incorrect number of indices for the rank of this object.");
  	vec_size_type temp=i+j*shape[0]+k*shape[1]*shape[0];
  	ESYS_ASSERT(temp < DataTypes::noValues(shape), "Invalid index.");
  	return temp;
  }

  /**
  \brief Compute the offset (in 1D vector) of a given subscript with a shape.

  \param shape - Input - Shape of the datapoint.
  \param i,j,k,m - Input - subscripts to locate.
  \return offset relative to the beginning of the datapoint.
  */
  inline
  vec_size_type
  getRelIndex(const DataTypes::ShapeType& shape, vec_size_type i,
	   vec_size_type j, vec_size_type k,
	   vec_size_type m)
  {
	// Warning: This is not C ordering. Do not try to figure out the params by looking at the code
	ESYS_ASSERT(getRank(shape)==4, "Incorrect number of indices for the rank of this object.");
	vec_size_type temp=i+j*shape[0]+k*shape[1]*shape[0]+m*shape[2]*shape[1]*shape[0];
	ESYS_ASSERT(temp < DataTypes::noValues(shape), "Invalid index.");
	return temp;
  }

  /**
     \brief Test if two shapes are equal.
  */
  inline
  bool
  checkShape(const ShapeType& s1, const ShapeType& s2)
  {
	return s1==s2;
  }

  /**
   \brief Produce a string containing two shapes.

   \param messagePrefix - Beginning of the message.
   \param other - displayed in the message as "Other shape"
   \param thisShape - displayed in the message as "This shape"
  */
   ESCRIPT_DLL_API
   std::string
   createShapeErrorMessage(const std::string& messagePrefix,
                                          const DataTypes::ShapeType& other,
					  const DataTypes::ShapeType& thisShape);

   inline
   bool
   checkOffset(vec_size_type offset, int size, int noval)
   {
      return (size >= (offset+noval));
   }

 }   // End of namespace DataTypes

} // End of namespace escript

#endif // __ESCRIPT_DATATYPES_H__

