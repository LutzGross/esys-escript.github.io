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

#include "EsysAssert.h"

#include "DataVector.h"
#include "LocalOps.h"

#include <boost/python/numeric.hpp>
#include <boost/python/object.hpp>

#include <vector>

namespace escript {

/**
   \brief
   DataArrayView provides a view of external data, configured according
   to the given shape information and offset. 

   Description:
   DataArrayView provides a view of data allocated externally. The 
   external data should provide sufficient values so that the dimensions
   specified for the view shape will be satisfied. The viewer can update
   values within the external data but cannot resize the external data.

   The view provided represents a single n-dimensional data-point
   comprised of values taken from the given data array, starting at the
   specified offset and extending for as many values as are necessary to
   satisfy the given shape. The default offset can be changed, or different
   offsets specified, in order to provide views of other data-points in
   the underlying data array.
*/

class DataArrayView {

  friend bool operator==(const DataArrayView& left, const DataArrayView& right);
  friend bool operator!=(const DataArrayView& left, const DataArrayView& right);

 public:

  //
  // Some basic types which define the data values and view shapes.
  typedef DataVector                        ValueType;
  //typedef std::vector<double>               ValueType;
  typedef std::vector<int>                  ShapeType;
  typedef std::vector<std::pair<int, int> > RegionType;
  typedef std::vector<std::pair<int, int> > RegionLoopRangeType;

  /**
     \brief
     Default constructor for DataArrayView.

     Description:
     Default constructor for DataArrayView.
     Creates an empty view with no associated data array.

     This is essentially useless but here for completeness.
  */
  DataArrayView();

  /**
     \brief
     Constructor for DataArrayView.

     Description:
     Constructor for DataArrayView.

     \param data - Input -
                Array holding data to be viewed. This must be large enough
                to supply sufficient values for the specified shape starting at
                the given offset.
     \param viewShape - Input -
                The shape of the view to create.
     \param offset - Input -
                Offset into the data at which view should start.
  */
  DataArrayView(ValueType& data,
                const ShapeType& viewShape,
                int offset=0);

  /**
     \brief
     Copy constructor for DataArrayView.

     Description:
     Copy constructor for DataArrayView.

     \param other - Input -
                DataArrayView to copy.

     NOTE: The copy references the same data array.
  */
  DataArrayView(const DataArrayView& other);

  /**
     \brief
     Copy from a numarray into the data array viewed by this.
     This must have same shape as the given value - will throw if not.
  */
  void
  copy(const boost::python::numeric::array& value);

  /**
     \brief
     Copy data from another DataArrayView into the data array viewed by this
     starting at the default offset in both views.
     The shapes of both views must be the same - will throw if not.
     NB: views may or may not reference same underlying data array!
  */
  void
  copy(const DataArrayView& other);

  /**
     \brief
     Copy data from another DataArrayView into this view starting at the 
     given offset in this view and the default offset in the other view.
     The shapes of both views must be the same - will throw if not.
     NB: views may or may not reference same underlying data array!
  */
  void
  copy(ValueType::size_type offset,
       const DataArrayView& other);

  /**
     \brief
     Copy data from another DataArrayView into this view starting at the 
     given offsets in each view.
     The shapes of both views must be compatible - will throw if not.
     NB: views may or may not reference same underlying data array!

     \param thisOffset - Input -
                   Offset into this view's data array to copy to.
     \param other - Input -
                   View to copy data from.
     \param otherOffset - Input -
                   Offset into the other view's data array to copy from.
  */
  void
  copy(ValueType::size_type thisOffset,
       const DataArrayView& other,
       ValueType::size_type otherOffset);

  /**
     \brief
     Copy the given single value over the entire view starting at the given
     offset in the view's data array.

     \param offset - Input -
                   Offset into this view's data array to copy to.
     \param value - Input -
                   Value to copy.
  */
  void
  copy(ValueType::size_type offset,
       ValueType::value_type value);

  /**
     \brief
     Check if view is empty. ie: does not point to any actual data.
  */
  bool
  isEmpty() const;

  /**
    \brief
    Return this view's offset into the viewed data array.
  */
  ValueType::size_type
  getOffset() const;

  /**
     \brief
     Set the offset into the data array at which this view is to start.
     This could be used to step through the underlying data array by incrementing
     the offset by noValues successively. Thus this view would provide a moving
     window on the underlying data with the given shape.
  */
  void
  setOffset(ValueType::size_type offset);

  /**
     \brief
     Increment the offset by the number of values in the shape, if possible. Thus
     moving the view onto the next data point of the given shape in the underlying
     data array.
  */
  void
  incrOffset();

  /**
     \brief
     Check the (given) offset will not result in two few values being available in
     the underlying data array for this view's shape.
  */
  bool
  checkOffset() const;

  bool
  checkOffset(ValueType::size_type offset) const;

  /**
    \brief
     Return the rank of the shape of this view.
  */
  int
  getRank() const;

  /**
     \brief
     Return the number of values for the shape of this view.
  */
  int
  noValues() const;

  /**
     \brief
     Calculate the number of values for the given shape or region.
     This is purely a utility method and has no bearing on this view.
  */
  static
  int
  noValues(const ShapeType& shape);

  static
  int
  noValues(const RegionLoopRangeType& region);

  /**
     \brief
     Return a reference to the underlying data array.
  */
  ValueType&
  getData() const;

  /**
     \brief
     Return a reference to the data value with the given
     index in this view. This takes into account the offset.
  */
  ValueType::reference
  getData(ValueType::size_type i) const;

  /**
     \brief
     Return the shape of this view.
  */
  const
  ShapeType&
  getShape() const;

  /**
     \brief
     Return true if the given shape is the same as this view's shape.
  */
  bool
  checkShape(const ShapeType& other) const;

  /**
     \brief
     Create a shape error message. Normally used when there is a shape
     mismatch between this shape and the other shape.

     \param messagePrefix - Input -
                       First part of the error message.
     \param other - Input -
                       The other shape.
  */
  std::string
  createShapeErrorMessage(const std::string& messagePrefix,
                          const ShapeType& other) const;

  /**
     \brief
     Return the viewed data as a formatted string.
     Not recommended for very large objects!

     \param suffix - Input -
                       Suffix appended to index display.
  */
  std::string
  toString(const std::string& suffix=std::string("")) const;

  /**
     \brief
     Return the given shape as a string.
     This is purely a utility method and has no bearing on this view.

     \param shape - Input.
  */
  static
  std::string
  shapeToString(const ShapeType& shape);

  /**
    \brief
    Return the 1 dimensional index into the data array of the only element
    in the view, *ignoring the offset*.
    Assumes a rank 0 view.
  */
  ValueType::size_type
  relIndex() const;

  /**
    \brief
    Return the 1 dimensional index into the data array of the element at
    position i in the view, *ignoring the offset*.
    Assumes a rank 1 view.
  */
  ValueType::size_type
  relIndex(ValueType::size_type i) const;

  /**
    \brief
    Return the 1 dimensional index into the data array of the element at
    position i,j in the view, *ignoring the offset*.
    Assumes a rank 2 view.
  */
  ValueType::size_type
  relIndex(ValueType::size_type i,
           ValueType::size_type j) const;

  /**
    \brief
    Return the 1 dimensional index into the data array of the element at
    position i,j,k in the view, *ignoring the offset*.
    Assumes a rank 3 view.
  */
  ValueType::size_type
  relIndex(ValueType::size_type i,
           ValueType::size_type j,
           ValueType::size_type k) const;

  /**
    \brief
    Return the 1 dimensional index into the data array of the element at
    position i,j,k,m in the view, *ignoring the offset*.
    Assumes a rank 4 view.
  */
  ValueType::size_type
  relIndex(ValueType::size_type i,
           ValueType::size_type j,
           ValueType::size_type k,
           ValueType::size_type m) const;

  /**
    \brief
    Return the 1 dimensional index into the data array of the only element 
    in the view.
    Assumes a rank 0 view.
  */
  ValueType::size_type
  index() const;

  /**
    \brief
    Return the 1 dimensional index into the data array of the element at
    position i in the view.
    Assumes a rank 1 view.
  */
  ValueType::size_type
  index(ValueType::size_type i) const;

  /**
    \brief
    Return the 1 dimensional index into the data array of the element at
    position i,j in the view.
    Assumes a rank 2 view.
  */
  ValueType::size_type
  index(ValueType::size_type i,
        ValueType::size_type j) const;

  /**
    \brief
    Return the 1 dimensional index into the data array of the element at
    position i,j,k in the view.
    Assumes a rank 3 view.
  */
  ValueType::size_type
  index(ValueType::size_type i,
        ValueType::size_type j,
        ValueType::size_type k) const;

  /**
    \brief
    Return the 1 dimensional index into the data array of the element at
    position i,j,k,m in the view.
    Assumes a rank 4 view.
  */
  ValueType::size_type
  index(ValueType::size_type i,
        ValueType::size_type j,
        ValueType::size_type k,
        ValueType::size_type m) const;

  /**
    \brief
    Return a reference for the only element in the view.
    Assumes a rank 0 view.
  */
  ValueType::reference
  operator()();

  ValueType::const_reference
  operator()() const;

  /**
    \brief
    Return a reference to the element at position i in the view.
    Assumes a rank 1 view.
  */
  ValueType::reference
  operator()(ValueType::size_type i);

  ValueType::const_reference
  operator()(ValueType::size_type i) const;

  /**
    \brief
    Return a reference to the element at position i,j in the view.
    Assumes a rank 2 view.
  */
  ValueType::reference
  operator()(ValueType::size_type i,
             ValueType::size_type j);

  ValueType::const_reference
  operator()(ValueType::size_type i,
             ValueType::size_type j) const;

  /**
    \brief
    Return a reference to the element at position i,j,k in the view.
    Assumes a rank 3 view.
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
    Return a reference to the element at position i,j,k,m in the view.
    Assumes a rank 4 view.
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
     Determine the shape of the specified slice region.
     This is purely a utility method and has no bearing on this view.

     \param region - Input -
                       Slice region.
  */
  static
  ShapeType
  getResultSliceShape(const RegionType& region);

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
  RegionType
  getSliceRegion(const boost::python::object& key) const;

  /**
     \brief
     Copy a data slice specified by the given region from the given view
     into this view, using the default offsets in both views.

     \param other - Input -
                      View to copy data from.
     \param region - Input -
                      Region in other view to copy data from.
  */
  void
  copySlice(const DataArrayView& other,
            const RegionLoopRangeType& region);

  /**
     \brief
     Copy a data slice specified by the given region and offset from the
     given view into this view at the given offset.

     \param thisOffset - Input -
                      Copy the slice to this offset in this view.
     \param other - Input -
                      View to copy data from.
     \param otherOffset - Input -
                      Copy the slice from this offset in the given view.
     \param region - Input -
                      Region in other view to copy data from.
  */
  void
  copySlice(ValueType::size_type thisOffset,
            const DataArrayView& other,
            ValueType::size_type otherOffset,
            const RegionLoopRangeType& region);

  /**
     \brief
     Copy data into a slice specified by the given region in this view from
     the given view, using the default offsets in both views.

     \param other - Input -
                  View to copy data from.
     \param region - Input -
                  Region in this view to copy data to.
  */
  void
  copySliceFrom(const DataArrayView& other,
                const RegionLoopRangeType& region);

  /**
     \brief
     Copy data into a slice specified by the given region and offset in
     this view from the given view at the given offset.

     \param thisOffset - Input -
                    Copy the slice to this offset in this view.
     \param other - Input -
                    View to copy data from.
     \param otherOffset - Input -
                    Copy the slice from this offset in the given view.
     \param region - Input -
                    Region in this view to copy data to.
  */
  void
  copySliceFrom(ValueType::size_type thisOffset,
                const DataArrayView& other,
                ValueType::size_type otherOffset,
                const RegionLoopRangeType& region);

  /**
     \brief
     Perform the unary operation on the data point specified by the view's
     default offset. Applies the specified operation to each value in the data
     point. 

     Called by escript::unaryOp.

     \param operation - Input -
                  Operation to apply. Operation must be a pointer to a function.
  */
  template <class UnaryFunction>
  void
  unaryOp(UnaryFunction operation);

  /**
     \brief
     Perform the unary operation on the data point specified by the given
     offset. Applies the specified operation to each value in the data
     point. Operation must be a pointer to a function.

     Called by escript::unaryOp.

     \param offset - Input -
                  Apply the operation to data starting at this offset in this view.
     \param operation - Input -
                  Operation to apply. Must be a pointer to a function.
  */
  template <class UnaryFunction>
  void
  unaryOp(ValueType::size_type offset,
          UnaryFunction operation);

  /**
     \brief
     Perform the binary operation on the data points specified by the default
     offsets in this view and in view "right". Applies the specified operation
     to corresponding values in both data points. Operation must be a pointer
     to a function.

     Called by escript::binaryOp.

     \param right - Input -
                  View to act as RHS in given binary operation.
     \param operation - Input -
                  Operation to apply. Must be a pointer to a function.
  */
  template <class BinaryFunction>
  void
  binaryOp(const DataArrayView& right,
           BinaryFunction operation);

  /**
     \brief
     Perform the binary operation on the data points specified by the given
     offsets in this view and in view "right". Applies the specified operation
     to corresponding values in both data points. Operation must be a pointer
     to a function.

     Called by escript::binaryOp.

     \param leftOffset - Input -
                  Apply the operation to data starting at this offset in this view.
     \param right - Input -
                  View to act as RHS in given binary operation.
     \param rightOffset - Input -
                  Apply the operation to data starting at this offset in right.
     \param operation - Input -
                  Operation to apply. Must be a pointer to a function.
  */
  template <class BinaryFunction>
  void
  binaryOp(ValueType::size_type leftOffset,
           const DataArrayView& right,
           ValueType::size_type rightOffset,
           BinaryFunction operation);

  /**
     \brief
     Perform the binary operation on the data point specified by the default
     offset in this view using the scalar value "right". Applies the specified
     operation to values in the data point. Operation must be a pointer to
     a function.

     Called by escript::binaryOp.

     \param right - Input -
                  Value to act as RHS in given binary operation.
     \param operation - Input -
                  Operation to apply. Must be a pointer to a function.
  */
  template <class BinaryFunction>
  void
  binaryOp(double right,
           BinaryFunction operation);

  /**
     \brief
     Perform the binary operation on the data point specified by the given
     offset in this view using the scalar value "right". Applies the specified
     operation to values in the data point. Operation must be a pointer
     to a function.

     Called by escript::binaryOp.

     \param offset - Input -
                  Apply the operation to data starting at this offset in this view.
     \param right - Input -
                  Value to act as RHS in given binary operation.
     \param operation - Input -
                  Operation to apply. Must be a pointer to a function.
  */
  template <class BinaryFunction>
  void
  binaryOp(ValueType::size_type offset,
           double right,
           BinaryFunction operation);

  /**
     \brief
     Perform the given data point reduction operation on the data point
     specified by the default offset into the view. Reduces all elements of
     the data point using the given operation, returning the result as a 
     scalar. Operation must be a pointer to a function.

     Called by escript::algorithm.

     \param operation - Input -
                  Operation to apply. Must be a pointer to a function.
  */
  template <class BinaryFunction>
  double
  reductionOp(BinaryFunction operation,
              double initial_value) const;

  /**
     \brief
     Perform the given data point reduction operation on the data point
     specified by the given offset into the view. Reduces all elements of
     the data point using the given operation, returning the result as a 
     scalar. Operation must be a pointer to a function.

     Called by escript::algorithm.

     \param offset - Input -
                  Apply the operation to data starting at this offset in this view.
     \param operation - Input -
                  Operation to apply. Must be a pointer to a function.
  */
  template <class BinaryFunction>
  double
  reductionOp(ValueType::size_type offset,
              BinaryFunction operation,
              double initial_value) const;

  /**
     \brief
     solves a local eigenvalue problem 

     \param in - Input - matrix 
     \param inOffset - Input - offset into in
     \param ev - Output - The eigenvalues
     \param inOffset - Input - offset into ev
  */
  static
  inline
  void
  DataArrayView::eigenvalues(DataArrayView& in,
                             ValueType::size_type inOffset,
                             DataArrayView& ev,
                             ValueType::size_type evOffset)
 {
   double in00,in10,in20,in01,in11,in21,in02,in12,in22;
   double ev0,ev1,ev2;
   int s=in.getShape()[0];
   if (s==1) {
      in00=(*(in.m_data))[inOffset+in.index(0,0)];
      eigenvalues1(in00,&ev0);
      (*(ev.m_data))[evOffset+ev.index(0)]=ev0;

   } else  if (s==2) {
      in00=(*(in.m_data))[inOffset+in.index(0,0)];
      in10=(*(in.m_data))[inOffset+in.index(1,0)];
      in01=(*(in.m_data))[inOffset+in.index(0,1)];
      in11=(*(in.m_data))[inOffset+in.index(1,1)];
      eigenvalues2(in00,(in01+in10)/2.,in11,&ev0,&ev1);
      (*(ev.m_data))[evOffset+ev.index(0)]=ev0;
      (*(ev.m_data))[evOffset+ev.index(1)]=ev1;

   } else  if (s==3) {
      in00=(*(in.m_data))[inOffset+in.index(0,0)];
      in10=(*(in.m_data))[inOffset+in.index(1,0)];
      in20=(*(in.m_data))[inOffset+in.index(2,0)];
      in01=(*(in.m_data))[inOffset+in.index(0,1)];
      in11=(*(in.m_data))[inOffset+in.index(1,1)];
      in21=(*(in.m_data))[inOffset+in.index(2,1)];
      in02=(*(in.m_data))[inOffset+in.index(0,2)];
      in12=(*(in.m_data))[inOffset+in.index(1,2)];
      in22=(*(in.m_data))[inOffset+in.index(2,2)];
      eigenvalues3(in00,(in01+in10)/2.,(in02+in20)/2.,in11,(in21+in12)/2.,in22,
                 &ev0,&ev1,&ev2);
      (*(ev.m_data))[evOffset+ev.index(0)]=ev0;
      (*(ev.m_data))[evOffset+ev.index(1)]=ev1;
      (*(ev.m_data))[evOffset+ev.index(2)]=ev2;

   }
 }

  /**
     \brief
     solves a local eigenvalue problem 

     \param in - Input - matrix 
     \param inOffset - Input - offset into in
     \param ev - Output - The eigenvalues
     \param evOffset - Input - offset into ev
     \param V - Output - The eigenvectors
     \param VOffset - Input - offset into V
     \param tol - Input - eigenvalues with relative difference tol are treated as equal
  */
  static
  inline
  void
  DataArrayView::eigenvalues_and_eigenvectors(DataArrayView& in,
                                              ValueType::size_type inOffset,
                                              DataArrayView& ev, 
                                              ValueType::size_type evOffset,
                                              DataArrayView& V, 
                                              ValueType::size_type VOffset,
                                              const double tol=1.e-13)
  {
   double ev0,ev1,ev2;
   int s=in.getShape()[0];
   if (s==1) {
      eigenvalues1(in(0,0),&ev0);
      ev(0)=ev0;

   } else  if (s==2) {
      eigenvalues2(in(0,0),(in(0,1)+in(1,0))/2.,in(1,1),
                   &ev0,&ev1);
      ev(0)=ev0;
      ev(1)=ev1;

   } else  if (s==3) {
      eigenvalues3(in(0,0),(in(0,1)+in(1,0))/2.,(in(0,2)+in(2,0))/2.,in(1,1),(in(2,1)+in(1,2))/2.,in(2,2),
                 &ev0,&ev1,&ev2);
      ev(0)=ev0;
      ev(1)=ev1;
      ev(2)=ev2;

   }
 }
 /**
     \brief
     Perform a matrix multiply of the given views.
     This is purely a utility method and has no bearing on this view.

     NB: Only multiplies together the two given views, ie: two data-points,
     would need to call this over all data-points to multiply the entire
     Data objects involved.

     \param left - Input - The left hand side.
     \param right - Input - The right hand side.
     \param result - Output - The result of the operation.
  */
  static
  void
  matMult(const DataArrayView& left,
          const DataArrayView& right,
          DataArrayView& result);

  /**
     \brief
     Determine the shape of the result array for a matrix multiplication
     of the given views.
     This is purely a utility method and has no bearing on this view.
  */
  static
  ShapeType
  determineResultShape(const DataArrayView& left,
                       const DataArrayView& right);

 protected:

 private:

  //
  // The maximum rank allowed for the shape of any view.
  static const int m_maxRank=4;

  //
  // The data values for the view.
  // NOTE: This points to data external to the view.
  // This is just a pointer to an array of ValueType.
  ValueType* m_data;

  //
  // The offset into the data array used by different views.
  // This is simply an integer specifying a position in the data array
  // pointed to by m_data.
  ValueType::size_type m_offset;

  //
  // The shape of the data.
  // This is simply an STL vector specifying the lengths of each dimension
  // of the shape as ints.
  ShapeType m_shape;

  //
  // The number of values needed for the array.
  // This can be derived from m_shape by multiplying the size of each dimension, but
  // is stored here for convenience.
  int m_noValues;

};

bool operator==(const DataArrayView& left, const DataArrayView& right);
bool operator!=(const DataArrayView& left, const DataArrayView& right);

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
DataArrayView::RegionLoopRangeType
getSliceRegionLoopRange(const DataArrayView::RegionType& region);

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

/**
   Inline function definitions.
*/

template <class UnaryFunction>
inline
void
DataArrayView::unaryOp(UnaryFunction operation)
{
  unaryOp(m_offset,operation);
}

template <class UnaryFunction>
inline
void
DataArrayView::unaryOp(ValueType::size_type offset,
                       UnaryFunction operation)
{
  EsysAssert((!isEmpty()&&checkOffset(offset)),
               "Error - Couldn't perform unaryOp due to insufficient storage.");
  for (ValueType::size_type i=0;i<noValues();i++) {
    (*m_data)[offset+i]=operation((*m_data)[offset+i]);
  }
}

template <class BinaryFunction>
inline
void
DataArrayView::binaryOp(const DataArrayView& right,
                        BinaryFunction operation)
{
  binaryOp(m_offset,right,right.getOffset(),operation);
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
	     "Error - Couldn't perform binaryOp due to shape mismatch,");
  EsysAssert((!isEmpty()&&checkOffset(leftOffset)),
             "Error - Couldn't perform binaryOp due to insufficient storage in left object.");
  EsysAssert((!right.isEmpty()&&right.checkOffset(rightOffset)),
             "Error - Couldn't perform binaryOp due to insufficient storage in right object.");
  for (ValueType::size_type i=0;i<noValues();i++) {
    (*m_data)[leftOffset+i]=operation((*m_data)[leftOffset+i],(*right.m_data)[rightOffset+i]);
  }
}

template <class BinaryFunction>
inline
void
DataArrayView::binaryOp(double right,
                        BinaryFunction operation)
{
  binaryOp(m_offset,right,operation);
}

template <class BinaryFunction>
inline
void
DataArrayView::binaryOp(ValueType::size_type offset,
                        double right,
                        BinaryFunction operation)
{
  EsysAssert((!isEmpty()&&checkOffset(offset)),
             "Error - Couldn't perform binaryOp due to insufficient storage in left object.");
  for (ValueType::size_type i=0;i<noValues();i++) {
    (*m_data)[offset+i]=operation((*m_data)[offset+i],right);
  }
}

template <class BinaryFunction>
inline
double
DataArrayView::reductionOp(BinaryFunction operation,
                           double initial_value) const
{
  return reductionOp(m_offset,operation,initial_value);
}

template <class BinaryFunction>
inline
double
DataArrayView::reductionOp(ValueType::size_type offset,
                           BinaryFunction operation,
                           double initial_value) const
{
  EsysAssert((!isEmpty()&&checkOffset(offset)),
               "Error - Couldn't perform reductionOp due to insufficient storage.");
  double current_value=initial_value;
  for (ValueType::size_type i=0;i<noValues();i++) {
    current_value=operation(current_value,(*m_data)[offset+i]);
  }
  return current_value;
}

inline
DataArrayView::ValueType::size_type 
DataArrayView::relIndex() const 
{
  EsysAssert((getRank()==0),"Incorrect number of indices for the rank of this object.");
  return 0;
}

inline
DataArrayView::ValueType::size_type 
DataArrayView::index() const 
{
  EsysAssert((getRank()==0),"Incorrect number of indices for the rank of this object.");
  return (m_offset);
}

inline
DataArrayView::ValueType::size_type 
DataArrayView::relIndex(ValueType::size_type i) const 
{
  EsysAssert((getRank()==1),"Incorrect number of indices for the rank of this object.");
  EsysAssert((i < noValues(m_shape)), "Error - Invalid index.");
  return i;
}

inline
DataArrayView::ValueType::size_type 
DataArrayView::index(ValueType::size_type i) const 
{
  EsysAssert((getRank()==1),"Incorrect number of indices for the rank of this object.");
  EsysAssert((i < noValues(m_shape)), "Error - Invalid index.");
  return (m_offset+i);
}

inline
DataArrayView::ValueType::size_type 
DataArrayView::relIndex(ValueType::size_type i,
                        ValueType::size_type j) const
{
  EsysAssert((getRank()==2),"Incorrect number of indices for the rank of this object.");
  ValueType::size_type temp=i+j*m_shape[0];
  EsysAssert((temp < noValues(m_shape)), "Error - Invalid index.");
  return temp;
}

inline
DataArrayView::ValueType::size_type
DataArrayView::index(ValueType::size_type i,
		     ValueType::size_type j) const
{
  EsysAssert((getRank()==2),"Incorrect number of indices for the rank of this object.");
  ValueType::size_type temp=i+j*m_shape[0];
  EsysAssert((temp < noValues(m_shape)), "Error - Invalid index.");
  return (m_offset+temp);
}

inline
DataArrayView::ValueType::size_type 
DataArrayView::relIndex(ValueType::size_type i,
			ValueType::size_type j,
			ValueType::size_type k) const 
{
  EsysAssert((getRank()==3),"Incorrect number of indices for the rank of this object.");
  ValueType::size_type temp=i+j*m_shape[0]+k*m_shape[1]*m_shape[0];
  EsysAssert((temp < noValues(m_shape)), "Error - Invalid index.");
  return temp;
}

inline
DataArrayView::ValueType::size_type 
DataArrayView::index(ValueType::size_type i,
		     ValueType::size_type j,
		     ValueType::size_type k) const 
{
  EsysAssert((getRank()==3),"Incorrect number of indices for the rank of this object.");
  ValueType::size_type temp=i+j*m_shape[0]+k*m_shape[1]*m_shape[0];
  EsysAssert((temp < noValues(m_shape)), "Error - Invalid index.");
  return (m_offset+temp);
}

inline
DataArrayView::ValueType::size_type 
DataArrayView::relIndex(ValueType::size_type i,
                        ValueType::size_type j,
                        ValueType::size_type k,
                        ValueType::size_type m) const
{
  EsysAssert((getRank()==4),"Incorrect number of indices for the rank of this object.");
  ValueType::size_type temp=i+j*m_shape[0]+k*m_shape[1]*m_shape[0]+m*m_shape[2]*m_shape[1]*m_shape[0];
  EsysAssert((temp < noValues(m_shape)), "Error - Invalid index.");
  return temp;
}

inline
DataArrayView::ValueType::size_type
DataArrayView::index(ValueType::size_type i,
		     ValueType::size_type j,
		     ValueType::size_type k,
		     ValueType::size_type m) const
{
  EsysAssert((getRank()==4),"Incorrect number of indices for the rank of this object.");
  ValueType::size_type temp=i+j*m_shape[0]+k*m_shape[1]*m_shape[0]+m*m_shape[2]*m_shape[1]*m_shape[0];
  EsysAssert((temp < noValues(m_shape)), "Error - Invalid index.");
  return (m_offset+temp);
}

inline
DataArrayView::ValueType::reference
DataArrayView::operator()()
{
  return (*m_data)[index()];
}

inline
DataArrayView::ValueType::const_reference
DataArrayView::operator()() const
{
  return (*m_data)[index()];
}

inline
DataArrayView::ValueType::reference
DataArrayView::operator()(ValueType::size_type i)
{
  return (*m_data)[index(i)];
}

inline
DataArrayView::ValueType::const_reference
DataArrayView::operator()(ValueType::size_type i) const
{
  return (*m_data)[index(i)];
}

inline
DataArrayView::ValueType::reference
DataArrayView::operator()(ValueType::size_type i,
                          ValueType::size_type j)
{
  return (*m_data)[index(i,j)];
}

inline
DataArrayView::ValueType::const_reference 
DataArrayView::operator()(ValueType::size_type i,
                          ValueType::size_type j) const
{
  return (*m_data)[index(i,j)];
}

inline
DataArrayView::ValueType::reference
DataArrayView::operator()(ValueType::size_type i,
                          ValueType::size_type j,
                          ValueType::size_type k)
{
  return (*m_data)[index(i,j,k)];
}

inline
DataArrayView::ValueType::const_reference
DataArrayView::operator()(ValueType::size_type i,
                          ValueType::size_type j,
                          ValueType::size_type k) const
{
  return (*m_data)[index(i,j,k)];
}

inline
DataArrayView::ValueType::reference
DataArrayView::operator()(ValueType::size_type i,
                          ValueType::size_type j,
                          ValueType::size_type k,
                          ValueType::size_type m)
{
  return (*m_data)[index(i,j,k,m)];
}

inline
DataArrayView::ValueType::const_reference
DataArrayView::operator()(ValueType::size_type i,
                          ValueType::size_type j,
                          ValueType::size_type k,
                          ValueType::size_type m) const
{
  return (*m_data)[index(i,j,k,m)];
}

} // end of namespace

#endif
