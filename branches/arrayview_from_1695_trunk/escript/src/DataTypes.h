
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
#include <boost/python/object.hpp>

namespace escript {

 namespace DataTypes {
  //
  // Some basic types which define the data values and view shapes.
  typedef DataVector                        ValueType;
  typedef std::vector<int>                  ShapeType;
  typedef std::vector<std::pair<int, int> > RegionType;
  typedef std::vector<std::pair<int, int> > RegionLoopRangeType;
  static const int maxRank=4;
  static const ShapeType scalarShape;

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

  /**
     \brief
     Determine the shape of the specified slice region.

     \param region - Input -
                       Slice region.
  */
  ESCRIPT_DLL_API
  DataTypes::ShapeType
  getResultSliceShape(const DataTypes::RegionType& region);


 /**
     \brief
     Determine the region specified by the given python slice object.

     \param key - Input -
                    python slice object specifying region to be returned.

     The slice object is a tuple of n python slice specifiers, where
     n <= the rank of this Data object. Each slice specifier specifies the
     range of indexes to be sliced from the corresponding dimension. The
     first specifier corresponds to the first dimension, the second to the
     second and so on. Where n < the rank, the remaining dimensions are
     sliced across the full range of their indicies.

     Each slice specifier is of the form "a:b", which specifies a slice
     from index a, up to but not including index b. Where index a is ommitted
     a is assumed to be 0. Where index b is ommitted, b is assumed to be the
     length of this dimension. Where both are ommitted (eg: ":") the slice is
     assumed to encompass that entire dimension.

     Where one of the slice specifiers is a single integer, eg: [1], we
     want to generate a rank-1 dimension object, as opposed to eg: [1,2]
     which implies we want to take a rank dimensional object with one
     dimension of size 1.

     The return value is a vector of pairs with length equal to the rank of
     this object. Each pair corresponds to the range of indicies from the
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

  */
   ESCRIPT_DLL_API
   DataTypes::RegionType
   getSliceRegion(const DataTypes::ShapeType& shape, const boost::python::object& key);



  ESCRIPT_DLL_API
  inline
  int
  getRank(const DataTypes::ShapeType& shape)
  {
	return shape.size();
  }

  ESCRIPT_DLL_API
  inline
  DataTypes::ValueType::size_type
  getRelIndex(const DataTypes::ShapeType& shape, DataTypes::ValueType::size_type i)
  {
  	EsysAssert((getRank(shape)==1),"Incorrect number of indices for the rank of this object.");
	EsysAssert((i < DataTypes::noValues(shape)), "Error - Invalid index.");
	return i;
  }

  ESCRIPT_DLL_API
  inline
  DataTypes::ValueType::size_type
  getRelIndex(const DataTypes::ShapeType& shape, DataTypes::ValueType::size_type i,
	   DataTypes::ValueType::size_type j)
  {
  	EsysAssert((getRank()==2),"Incorrect number of indices for the rank of this object.");
  	DataTypes::ValueType::size_type temp=i+j*shape[0];
  	EsysAssert((temp < DataTypes::noValues(shape)), "Error - Invalid index.");
	return temp;
  }

  ESCRIPT_DLL_API
  inline
  DataTypes::ValueType::size_type
  getRelIndex(const DataTypes::ShapeType& shape, DataTypes::ValueType::size_type i,
	   DataTypes::ValueType::size_type j, DataTypes::ValueType::size_type k)
  {
  	EsysAssert((getRank()==3),"Incorrect number of indices for the rank of this object.");
  	DataTypes::ValueType::size_type temp=i+j*shape[0]+k*shape[1]*shape[0];
  	EsysAssert((temp < DataTypes::noValues(shape)), "Error - Invalid index.");
  	return temp;
  }

  ESCRIPT_DLL_API
  inline
  DataTypes::ValueType::size_type
  getRelIndex(const DataTypes::ShapeType& shape, DataTypes::ValueType::size_type i,
	   DataTypes::ValueType::size_type j, DataTypes::ValueType::size_type k,
	   DataTypes::ValueType::size_type m)
  {
	EsysAssert((getRank()==4),"Incorrect number of indices for the rank of this object.");
	DataTypes::ValueType::size_type temp=i+j*shape[0]+k*shape[1]*shape[0]+m*shape[2]*shape[1]*shape[0];
	EsysAssert((temp < DataTypes::noValues(m_shape)), "Error - Invalid index.");
	return temp;
  }


 }   // End namespace DataTypes


} // end of namespace escipt

#endif
