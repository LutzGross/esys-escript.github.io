// $Id$
/* 
 ******************************************************************************
 *                                                                            *
 *       COPYRIGHT  ACcESS 2004 -  All Rights Reserved                        *
 *                                                                            *
 * This software is the property of ACcESS. No part of this code              *
 * may be copied in any form or by any means without the expressed written    *
 * consent of ACcESS.  Copying, use or modification of this software          *
 * by any unauthorised person is illegal unless that person has a software    *
 * license agreement with ACcESS.                                             *
 *                                                                            *
 ******************************************************************************
*/
                                                                           
#if !defined escript_DataArrayView_20040323_H
#define escript_DataArrayView_20040323_H

#include "esysUtils/EsysAssert.h"

#include <vector>
#include <boost/python/numeric.hpp>
#include <boost/python/object.hpp>
#include <boost/shared_ptr.hpp>

#include <iostream>

namespace escript {

/**
   \brief
   DataArrayView provides a view of data allocated externally. 

   Description:
   DataArrayView provides a view of data allocated externally. The 
   external data should provide sufficient values so that the dimensions
   specified for the view will be satisfied. The viewer can update
   values within the external data but cannot resize the external data.
*/

class DataArrayView {

  friend bool operator==(const DataArrayView& left, const DataArrayView& right);
  friend bool operator!=(const DataArrayView& left, const DataArrayView& right);

 public:

  typedef std::vector<double> ValueType;
  typedef std::vector<int> ShapeType;
  typedef std::vector<std::pair<int, int> > RegionType;

  /**
     \brief
     Default constructor for DataArrayView.

     Description:
     Default constructor for DataArrayView. Creates an
     empty view with no associated data.
  */
  DataArrayView();

  /**
     \brief
     Constructor for DataArrayView.

     Description:
     Constructor for DataArrayView.

     \param data - Input - Container holding data to be viewed. This must
                           be large enough for the specified shape.
     \param viewShape - Input - The shape of the view.
     \param offset - Input - Offset into the data at which view should start.
  */
  DataArrayView(ValueType& data,
                const ShapeType& viewShape,
                int offset=0);

  /**
     \brief
     Copy constructor for DataArrayView.

     Description:
     Copy constructor for DataArrayView.

     \param other - Input - DataArrayView to copy.

     NOTE: The copy references the same data container.
  */
  DataArrayView(const DataArrayView& other);

  /**
     \brief
     Check if view is empty.
  */
  bool
  isEmpty() const;

  /**
     \brief
     Copy from a numarray into the data managed by DataArrayView.
  */
  void
  copy(const boost::python::numeric::array& value);

  /**
     \brief
     Copy from another DataArrayViewinto the data managed by this
     DataArrayView at the default offset.
  */
  void
  copy(const DataArrayView& other);

  /**
     \brief
     Copy from another DataArrayView into this view at the 
     given offset.
  */
  void
  copy(ValueType::size_type offset,
       const DataArrayView& other);

  /**
     \brief
     Copy from another DataArrayView into this view at the 
     given offsets.

     \param thisOffset - Input - Offset into this view's data array to copy to.
     \param other - Input - View for the data to copy.
     \param otherOffset - Input - Offset into the other view to copy data from.
  */
  void
  copy(ValueType::size_type thisOffset,
       const DataArrayView& other,
       ValueType::size_type otherOffset);

  /**
     \brief
     Copy the given single value over the entire view.

     \param offset - Input - Offset into this view's data array to copy to.
     \param value - Input - Value to copy.
  */
  void
  copy(ValueType::size_type offset,
       ValueType::value_type value);

  /**
    \brief
    Return this view's offset to the start of data.
  */
  ValueType::size_type
  getOffset() const;

  /**
    \brief
     Return the rank of the data.
  */
  int
  getRank() const;

  /**
     \brief
     Return the shape of the data.
  */
  const ShapeType&
  getShape() const;

  /**
     \brief
     Calculate the number of values for the given shape.
  */
  static int
  noValues(const ShapeType& shape);

  /**
     \brief
     Return the number of values for the current view.
  */
  int
  noValues() const;

  /**
     \brief
     Return true if the given shape is the same as this object's shape.
  */
  bool
  checkShape(const DataArrayView::ShapeType& other) const;

  /**
     \brief
     Return a reference to the underlying data.
  */
  ValueType&
  getData() const;

  /**
     \brief
     Return the value with the given index for the view. This takes into account
     the offset. Effectively this returns each value of the view in sequential order.
  */
  ValueType::reference
  getData(ValueType::size_type i) const;

  /**
     \brief
     Create a shape error message. Normally used when there is a shape 
     mismatch.

     \param messagePrefix - Input - First part of the error message. 
     \param other - Input - The other shape. 
  */
  std::string
  createShapeErrorMessage(const std::string& messagePrefix,
                          const DataArrayView::ShapeType& other) const;

  /**
     \brief
     Set the offset into the array.
  */
  void
  setOffset(ValueType::size_type offset);

  /**
     \brief
     Return the shape of a data view following the slice specified.

     \param region - Input - Slice region.
  */
  static DataArrayView::ShapeType
  getResultSliceShape(const RegionType& region);

  /**
     \brief
     Copy a slice specified by the given region from the given view
     into this view.
     The given region must be the same shape as this view.

     \param other - Input - view to copy from.
     \param region - Input - region in other to copy.
  */
  void
  copySlice(const DataArrayView& other,
            const RegionType& region);

  /**
     \brief
     Copy a slice specified by the given region from the given view into this view.
     The given region must be the same shape as this view.

     \param thisOffset - Input - use this offset into this object instead of the default.
     \param other - Input - view to copy from.
     \param otherOffset - Input - use this offset into the given object instead of the default.
     \param region - Input - region in other to copy.
  */
  void
  copySlice(ValueType::size_type thisOffset,
            const DataArrayView& other,
            ValueType::size_type otherOffset,
            const RegionType& region);

 /**
     \brief
     Copy into a slice from the given view.
     This view must have the same rank as the slice region.

     \param other - Input - Data to copy from.
     \param region - Input - Slice region.
  */
  void
  copySliceFrom(const DataArrayView& other,
                const RegionType& region);

  /**
     \brief
     Copy into a slice from the given value.
     This view must have the same rank as the slice region.

     \param thisOffset - Input - use this view offset instead of the default.
     \param other - Input - Data to copy from.
     \param otherOffset - Input - use this slice offset instead of the default.
     \param region - Input - Slice region.
  */
  void
  copySliceFrom(ValueType::size_type thisOffset,
                const DataArrayView& other,
                ValueType::size_type otherOffset,
                const RegionType& region);

  /**
     \brief
     Determine the shape of the result array for a matrix multiplication.
  */
  static ShapeType
  determineResultShape(const DataArrayView& left,
                       const DataArrayView& right);

  /**
     \brief
     Determine the region specified by the given python slice object.

     \param key - Input - python slice object specifying region to be returned.

     \description

     The slice object is a tuple of n python slice specifiers, where
     n <= the rank of this Data object. Each slice specifier specifies the
     range of indexes to be sliced from the corresponding dimension. The
     first specifier corresponds to the first dimension, the second to the
     second and so on. Where n < the rank, the remaining dimensions are
     sliced across the full range of their indicies.

     Each slice specifier is of the form "a:b", which specifies a slice
     from index a, up to but not including index b. Where index a is ommitted
     a is assumed to be 0. Where index b is ommitted, b is assumed to be the
     length of this dimension.

     The return value is a vector of pairs with length equal to the rank of
     this object. Each pair corresponds to the range of indexes from the
     corresponding dimension to be sliced from, as specified in the input
     slice object.

     Examples:

       For a rank 1 object of shape(5):

         getSliceRegion(:)   => < <0,5> >
         getSliceRegion(2:3) => < <2,3> >
         getSliceRegion(:3)  => < <0,3> >
         getSliceRegion(2:)  => < <2,5> >

       For a rank 3 object of shape (2,4,6):

         getSliceRegion(0:2,0:4,0:6) => < <0,2> <0,4> <0,6> >
         getSliceRegion(:,:,:)       => < <0,2> <0,4> <0,6> >
         getSliceRegion(0:1)         => < <0,1> <0,4> <0,6> >
         getSliceRegion(:1,0:2)      => < <0,1> <0,2> <0,6> >
  */
  DataArrayView::RegionType
  getSliceRegion(const boost::python::object& key) const;

  // *******************************************************************
  // NOTE: The following relIndex functions are a hack. The indexing
  // mechanism should be split out from DataArrayView to get the flexability
  // needed.

  /**
    \brief
    Return the 1 dimensional index of the element at position i,j,k,m
    ignoring the offset.
  */
  ValueType::size_type
  relIndex(ValueType::size_type i,
           ValueType::size_type j,
           ValueType::size_type k,
           ValueType::size_type m) const;

  /**
    \brief
    Return the 1 dimensional index of the element at position i,j,k
    ignoring the offset.
  */
  ValueType::size_type
  relIndex(ValueType::size_type i,
           ValueType::size_type j,
           ValueType::size_type k) const;

  /**
    \brief
    Return the 1 dimensional index of the element at position i,j
    ignoring the offset.
  */
  ValueType::size_type
  relIndex(ValueType::size_type i,
           ValueType::size_type j) const;

  // ********************************************************************

  /**
    \brief
    Return the 1 dimensional index of the element at position i,j,k,m.
  */
  ValueType::size_type
  index(ValueType::size_type i,
        ValueType::size_type j,
        ValueType::size_type k,
        ValueType::size_type m) const;

  /**
    \brief
    Return the 1 dimensional index of the element at position i,j,k.
  */
  ValueType::size_type
  index(ValueType::size_type i,
        ValueType::size_type j,
        ValueType::size_type k) const;

  /**
    \brief
    Return the 1 dimensional index of the element at position i,j.
  */
  ValueType::size_type
  index(ValueType::size_type i,
        ValueType::size_type j) const;

  /**
    \brief
    Return the 1 dimensional index of the element at position i.
  */
  ValueType::size_type
  index(ValueType::size_type i) const;

  /**
    \brief
    Return a reference to the element at position i.
  */
  ValueType::reference
  operator()(ValueType::size_type i);

  ValueType::const_reference
  operator()(ValueType::size_type i) const;

  /**
    \brief
    Return a reference to the element at position i,j.
  */
  ValueType::reference
  operator()(ValueType::size_type i,
             ValueType::size_type j);

  ValueType::const_reference
  operator()(ValueType::size_type i,
             ValueType::size_type j) const;

  /**
    \brief
    Return a reference to the element at position i,j,k.
  */
  ValueType::reference
  operator()(ValueType::size_type i,
             ValueType::size_type j,
             ValueType::size_type k);

  ValueType::const_reference
  operator()(ValueType::size_type i,
             ValueType::size_type j,
             ValueType::size_type k) const;

 /**
    \brief
    Return a reference to the element at position i,j,k,m.
  */
  ValueType::reference
  operator()(ValueType::size_type i,
             ValueType::size_type j,
             ValueType::size_type k,
             ValueType::size_type m);

  ValueType::const_reference
  operator()(ValueType::size_type i,
             ValueType::size_type j,
             ValueType::size_type k,
             ValueType::size_type m) const;

  /**
    \brief
    Return a reference for the only element, assuming rank 0.
  */
  ValueType::reference
  operator()();

  ValueType::const_reference
  operator()() const;

  /**
     \brief
     Perform the unary operation using the given offset
     instead of the offset defined within the view.
  */
  template <class UnaryFunction>
  void
  unaryOp(ValueType::size_type leftOffset,
          UnaryFunction operation);

  /**
     \brief
     Perform the unary operation using the view's offset.
  */
  template <class UnaryFunction>
  void
  unaryOp(UnaryFunction operation);

  /**
     \brief
     Perform the given data point reduction operation on the data point
     specified by the given offset into the view. Reduces all elements of
     the data point using the given operation, returning the result as a 
     scalar.

     Called by escript::dp_algorithm.
  */
  template <class UnaryFunction>
  double
  dp_algorithm(ValueType::size_type leftOffset,
               UnaryFunction operation);

  /**
     \brief
     Perform the given data point reduction operation on the data point
     specified by the default offset into the view. Reduces all elements of
     the data point using the given operation, returning the result as a 
     scalar.

     Called by escript::dp_algorithm.
  */
  template <class UnaryFunction>
  double
  dp_algorithm(UnaryFunction operation);

  /**
     \brief
     Perform the given operation and return a result.
  */
  template <class UnaryFunction>
  double
  algorithm(ValueType::size_type leftOffset,
            UnaryFunction operation) const;

  /**
     \brief
     Perform the given operation and return a result. Use the default offset.
  */
  template <class UnaryFunction>
  double
  algorithm(UnaryFunction operation) const;

  /**
     \brief
     Perform the binary operation. Version which applies a double value 
     to all values within the view. The given offset is used instead of 
     the default offset specified within the view.
  */
  template <class BinaryFunction>
  void
  binaryOp(ValueType::size_type leftOffset,
           double right,
           BinaryFunction operation);

  /**
     \brief
     Perform the binary operation. Version which applies a double value 
     to all values within the view. 
  */
  template <class BinaryFunction>
  void
  binaryOp(double right,
           BinaryFunction operation);

  /**
     \brief
     Perform the binary operation. The given offsets override the default 
     offsets specified within the views.
  */
  template <class BinaryFunction>
  void
  binaryOp(ValueType::size_type leftOffset,
           const DataArrayView& right,
           ValueType::size_type rightOffset,
           BinaryFunction operation);

  /**
     \brief
     Perform the binary operation.
  */
  template <class BinaryFunction>
  void
  binaryOp(const DataArrayView& right,
           BinaryFunction operation);

  /**
     \brief
     Return the data as a string. Not recommended for very large objects.
     \param suffix - Input - Suffix appended to index display.
  */
  std::string
  toString(const std::string& suffix=std::string("")) const;

  /**
     \brief
     Return the given shape as a string.

     \param shape - Input.
  */
  static std::string
  shapeToString(const DataArrayView::ShapeType& shape);

  /**
     \brief
     Perform matrix multiply.

     Description:
     Perform matrix multiply.

     \param left - Input - The left hand side.
     \param right - Input - The right hand side.
     \param result - Output - The result of the operation.
  */
  static void
  matMult(const DataArrayView& left,
          const DataArrayView& right,
          DataArrayView& result);

 protected:

 private:

  //
  static const int m_maxRank=4;

  //
  // The data values for the view. NOTE: This points to data external to the view.
  ValueType* m_data;

  //
  // The offset into the data array used by different views.
  ValueType::size_type m_offset;

  //
  // The shape of the data.
  ShapeType m_shape;

  //
  // The number of values needed for the array.
  int m_noValues;

};

/**
  \brief
  Calculate the slice range from the given python key object
  Used by DataArrayView::getSliceRegion.
  Returns the python slice object key as a pair of ints where the first 
  member is the start and the second member is the end. the presence of a possible
  step attribute with value different from one will throw an exception

  /param key - Input - key object specifying slice range.
*/
std::pair<int,int>
getSliceRange(const boost::python::object& key,
              const int shape);

inline
DataArrayView::ValueType::size_type
DataArrayView::getOffset() const
{
  return m_offset;
}

inline
DataArrayView::ValueType&
DataArrayView::getData() const 
{
  EsysAssert(!isEmpty(),"Error - View is empty");
  return *m_data;
}

inline
DataArrayView::ValueType::reference
DataArrayView::getData(ValueType::size_type i) const 
{
  //
  // don't do any checking to allow one past the end of the vector providing 
  // the equivalent of end()
  return (*m_data)[i+m_offset];
}

template <class UnaryFunction>
inline
void
DataArrayView::unaryOp(ValueType::size_type leftOffset,
                       UnaryFunction operation)
{
  for (ValueType::size_type i=0;i<(noValues(m_shape));i++) {
    (*m_data)[i+leftOffset]=operation((*m_data)[i+leftOffset]);
  }
}

template <class UnaryFunction>
inline
void
DataArrayView::unaryOp(UnaryFunction operation)
{
  unaryOp(m_offset,operation);
}

template <class UnaryFunction>
inline
double
DataArrayView::dp_algorithm(ValueType::size_type leftOffset,
                            UnaryFunction operation)
{
  operation.resetResult();
  for (ValueType::size_type i=0;i<(noValues(m_shape));i++) {
    operation((*m_data)[i+leftOffset]);
  }
  return operation.getResult();
}

template <class UnaryFunction>
inline
double
DataArrayView::dp_algorithm(UnaryFunction operation)
{
  return dp_algorithm(m_offset,operation);
}

template <class UnaryFunction>
inline
double
DataArrayView::algorithm(ValueType::size_type leftOffset,
                         UnaryFunction operation) const
{
  for (ValueType::size_type i=0;i<(noValues(m_shape));++i) {
    operation((*m_data)[i+leftOffset]);
  }
  return operation.getResult();
}

template <class UnaryFunction>
inline
double
DataArrayView::algorithm(UnaryFunction operation) const
{
  return algorithm(m_offset,operation);
}

template <class BinaryFunction>
inline
void
DataArrayView::binaryOp(ValueType::size_type leftOffset,
                        const DataArrayView& right, 
                        ValueType::size_type rightOffset,
                        BinaryFunction operation)
{
  EsysAssert(getShape()==right.getShape(),
	     "Error - Right hand shape: " << shapeToString(right.getShape()) << " doesn't match the left: " << shapeToString(getShape()));
  for (ValueType::size_type i=0;i<noValues();++i) {
    (*m_data)[i+leftOffset]=operation((*m_data)[i+leftOffset],(*right.m_data)[i+rightOffset]);
  }
}

template <class BinaryFunction>
inline
void
DataArrayView::binaryOp(const DataArrayView& right,
                        BinaryFunction operation)
{
  //
  // version using the offsets specified within the view
  binaryOp(m_offset,right,right.getOffset(),operation);
}

template <class BinaryFunction>
inline
void
DataArrayView::binaryOp(ValueType::size_type leftOffset,
                        double right,
                        BinaryFunction operation)
{
  //
  // If a scalar is to be applied to the entire array force the caller to
  // explicitly specify a single value
  for (ValueType::size_type i=0;i<(noValues(m_shape));++i) {
    (*m_data)[i+leftOffset]=operation((*m_data)[i+leftOffset],right);
  }
}

template <class BinaryFunction>
inline
void
DataArrayView::binaryOp(double right,
                        BinaryFunction operation)
{
  //
  // use the default offset
  binaryOp(m_offset,right,operation);
}

inline
DataArrayView::ValueType::size_type 
DataArrayView::index(ValueType::size_type i) const 
{
    //EsysAssert((i >= 0) && (i < noValues(m_shape)), "Invalid index.");
    EsysAssert((i < noValues(m_shape)), "Invalid index.");
    return (i+m_offset);
}

inline
DataArrayView::ValueType::size_type 
DataArrayView::relIndex(ValueType::size_type i,
                        ValueType::size_type j) const
{
  ValueType::size_type temp=i+j*m_shape[0];
  //EsysAssert((temp >= 0 || temp < noValues(m_shape)), "Invalid index.");
  EsysAssert((temp < noValues(m_shape)), "Invalid index.");
  //
  // no offset
  return temp;
}

inline
DataArrayView::ValueType::size_type
DataArrayView::index(ValueType::size_type i,
		     ValueType::size_type j) const
{
  ValueType::size_type temp=i+j*m_shape[0];
  //EsysAssert((temp >= 0 || temp < noValues(m_shape)), "Invalid index.");
  EsysAssert((temp < noValues(m_shape)), "Invalid index.");
  return (temp+m_offset);
}

inline
DataArrayView::ValueType::size_type 
DataArrayView::relIndex(ValueType::size_type i,
			ValueType::size_type j,
			ValueType::size_type k) const 
{
    ValueType::size_type temp=i+j*m_shape[0]+k*m_shape[1]*m_shape[0];
    //EsysAssert((temp >= 0 || temp < noValues(m_shape)), "Invalid index.");
    EsysAssert((temp < noValues(m_shape)), "Invalid index.");
    //
    // no offset
    return temp;
}

inline
DataArrayView::ValueType::size_type 
DataArrayView::index(ValueType::size_type i,
		     ValueType::size_type j,
		     ValueType::size_type k) const 
{
    ValueType::size_type temp=i+j*m_shape[0]+k*m_shape[1]*m_shape[0];
    //EsysAssert((temp >= 0 || temp < noValues(m_shape)), "Invalid index.");
    EsysAssert((temp < noValues(m_shape)), "Invalid index.");
    return (temp+m_offset);
}

inline
DataArrayView::ValueType::size_type 
DataArrayView::relIndex(ValueType::size_type i,
                        ValueType::size_type j,
                        ValueType::size_type k,
                        ValueType::size_type m) const
{
  ValueType::size_type temp=i+j*m_shape[0]+k*m_shape[1]*m_shape[0]+m*m_shape[2]*m_shape[1]*m_shape[0];
  //EsysAssert((temp >= 0 || temp < noValues(m_shape)), "Invalid index.");
  EsysAssert((temp < noValues(m_shape)), "Invalid index.");
  //
  // no offset
  return temp;
}

inline
DataArrayView::ValueType::size_type
DataArrayView::index(ValueType::size_type i,
		     ValueType::size_type j,
		     ValueType::size_type k,
		     ValueType::size_type m) const
{
  ValueType::size_type temp=i+j*m_shape[0]+k*m_shape[1]*m_shape[0]+m*m_shape[2]*m_shape[1]*m_shape[0];
  //EsysAssert((temp >= 0 || temp < noValues(m_shape)), "Invalid index.");
  EsysAssert((temp < noValues(m_shape)), "Invalid index.");
  return (temp+m_offset);
}

inline
DataArrayView::ValueType::reference
DataArrayView::operator()(ValueType::size_type i,
                          ValueType::size_type j,
                          ValueType::size_type k,
                          ValueType::size_type m)
{
    EsysAssert((getRank()==4),
	       "Incorrect number of indices for the rank of this object.");
    return (*m_data)[index(i,j,k,m)];
}

inline
DataArrayView::ValueType::const_reference
DataArrayView::operator()(ValueType::size_type i,
                          ValueType::size_type j,
                          ValueType::size_type k,
                          ValueType::size_type m) const
{
    EsysAssert((getRank()==4),
	       "Incorrect number of indices for the rank of this object.");
    return (*m_data)[index(i,j,k,m)];
}

inline
DataArrayView::ValueType::reference
DataArrayView::operator()(ValueType::size_type i,
                          ValueType::size_type j,
                          ValueType::size_type k)
{
    EsysAssert((getRank()==3),
	       "Incorrect number of indices for the rank of this object.");
    return (*m_data)[index(i,j,k)];
}

inline
DataArrayView::ValueType::const_reference
DataArrayView::operator()(ValueType::size_type i,
                          ValueType::size_type j,
                          ValueType::size_type k) const
{
    EsysAssert((getRank()==3),
	       "Incorrect number of indices for the rank of this object.");
    return (*m_data)[index(i,j,k)];
}

inline
DataArrayView::ValueType::reference
DataArrayView::operator()(ValueType::size_type i,
                          ValueType::size_type j)
{
    EsysAssert((getRank()==2),
	       "Incorrect number of indices for the rank of this object.");
    return (*m_data)[index(i,j)];
}

inline
DataArrayView::ValueType::const_reference 
DataArrayView::operator()(ValueType::size_type i,
                          ValueType::size_type j) const
{
    EsysAssert((getRank()==2),
	       "Incorrect number of indices for the rank of this object.");
    return (*m_data)[index(i,j)];
}

inline
DataArrayView::ValueType::reference
DataArrayView::operator()(ValueType::size_type i)
{
    EsysAssert((getRank()==1),
	       "Incorrect number of indices for the rank of this object.");
    return (*m_data)[index(i)];
}

inline
DataArrayView::ValueType::const_reference
DataArrayView::operator()(ValueType::size_type i) const
{
    EsysAssert((getRank()==1),
	       "Incorrect number of indices for the rank of this object.");
    return (*m_data)[index(i)];
}

inline
DataArrayView::ValueType::reference
DataArrayView::operator()()
{
    EsysAssert((getRank()==0),
	       "Incorrect number of indices for the rank of this object.");
    return (*m_data)[m_offset];
}

inline
DataArrayView::ValueType::const_reference
DataArrayView::operator()() const
{
    EsysAssert((getRank()==0),
	       "Incorrect number of indices for the rank of this object.");
    return (*m_data)[m_offset];
}

bool operator==(const DataArrayView& left, const DataArrayView& right);
bool operator!=(const DataArrayView& left, const DataArrayView& right);

} // end of namespace

#endif
