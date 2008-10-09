
/*******************************************************
*
* Copyright (c) 2003-2008 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


#if !defined escript_DataConstant_20040323_H
#define escript_DataConstant_20040323_H
#include "system_dep.h"

#include "DataReady.h"

#include <boost/python/numeric.hpp>

namespace escript {

/**
   \brief
   DataConstant stores a single data point which represents the entire
   function space.

   Description:
   DataConstant stores a single data point which represents the entire
   function space.
*/
class DataConstant : public DataReady  {
typedef DataReady parent;
 public:

  /**
     \brief
     Constructor for DataConstant objects.

     Description:
     Constructor for DataConstant objects.

     \param value - Input - Data value for a single point.
     \param what - Input - A description of what this data object represents.
  */
  ESCRIPT_DLL_API
  DataConstant(const boost::python::numeric::array& value,
               const FunctionSpace& what);

  /**
     \brief
     Copy constructor. Performs a deep copy.
  */
  ESCRIPT_DLL_API
  DataConstant(const DataConstant& other);

//  /**
/*     \brief
     Alternative constructor for DataConstant objects.

     Description:
     Alternative Constructor for DataConstant objects.
     \param value - Input - Data value for a single point.
     \param what - Input - A description of what this data object represents.*/
//  */
//   ESCRIPT_DLL_API
//   DataConstant(const DataArrayView& value,
//                const FunctionSpace& what);

  /**
     \brief
     Alternative constructor for DataConstant objects.

     Description:
     Alternative Constructor for DataConstant objects.
     \param other - Input - Data object to copy from.
     \param region - Input - region to copy.
  */
  ESCRIPT_DLL_API
  DataConstant(const DataConstant& other,
               const DataTypes::RegionType& region);

  /**
     \brief
     Alternative constructor for DataConstant objects.

     Description:
     Alternative Constructor for DataConstant objects.
     \param what - Input - A description of what this data object represents.
     \param shape - Input - the shape of each data-point.
     \param data - the data values for each data-point.
  */
  ESCRIPT_DLL_API
  DataConstant(const FunctionSpace& what,
               const DataTypes::ShapeType &shape,
               const DataTypes::ValueType &data);

  /**
     \brief
     Write the data as a string.
  */
  ESCRIPT_DLL_API
  std::string
  toString() const;

  /**
     \brief Return a deep copy of the current object.
  */
  ESCRIPT_DLL_API
  virtual
  DataAbstract*
  deepCopy();


 /**
     \brief
     dumps the object into a netCDF file
  */
  ESCRIPT_DLL_API
  virtual
  void
  dump(const std::string fileName) const;

 /**
     \brief
    sets all values to zero
  */
  ESCRIPT_DLL_API
  virtual
  void
  setToZero();

  /**
     \brief
     Return the offset for the given sample. This is a somewhat artificial notion
     but returns the offset in bytes for the given point into the container
     holding the point data. Only really necessary to avoid many DataArrayView
     objects.
     \param sampleNo - Input - sample number.
     \param dataPointNo - Input - data point number for the sample.
   */
  ESCRIPT_DLL_API
  virtual
  DataTypes::ValueType::size_type
  getPointOffset(int sampleNo,
                 int dataPointNo) const;

//  /**
/*     \brief
     Return a view into the data for the data point specified.
     \param sampleNo - Input - sample number.
     \param dataPointNo - Input - data point number for the sample.*/
//  */
//   ESCRIPT_DLL_API
//   virtual
//   DataArrayView
//   getDataPoint(int sampleNo,
//                int dataPointNo);

  /**
     \brief
     Return the number of doubles stored for the Data object.
  */
  ESCRIPT_DLL_API
  virtual
  DataTypes::ValueType::size_type
  getLength() const;

  /**
     \brief
     Factory method that returns a newly created DataConstant object
     sliced from the specified region of this object.
     The caller is reponsible for managing the object created.
     \param region - Input - region to slice from this object.
  */
  ESCRIPT_DLL_API
  virtual
  DataAbstract*
  getSlice(const DataTypes::RegionType& region) const;

  /**
     \brief
     Copy the specified region from the given value.
     \param value - Input - Data object to copy from.
     \param region - Input - Region to copy.
  */
  ESCRIPT_DLL_API
  virtual
  void
  setSlice(const DataAbstract* value,
           const DataTypes::RegionType& region);


  /**
     \brief
     Computes a symmetric matrix (A + AT) / 2

     \param ev - Output - symmetric matrix

  */
  ESCRIPT_DLL_API
  virtual void
  symmetric(DataAbstract* ev);

  /**
     \brief
     Computes a nonsymmetric matrix (A - AT) / 2

     \param ev - Output - nonsymmetric matrix

  */
  ESCRIPT_DLL_API
  virtual void
  nonsymmetric(DataAbstract* ev);

  /**
     \brief
     Computes the trace of a matrix

     \param ev - Output - trace of matrix

  */
  ESCRIPT_DLL_API
  virtual void
  trace(DataAbstract* ev, int axis_offset);

  /**
     \brief
     Transpose each data point of this Data object around the given axis.

     \param ev - Output - transpose of matrix

  */
  ESCRIPT_DLL_API
  virtual void
  transpose(DataAbstract* ev, int axis_offset);

  /**
     \brief
     swaps components axis0 and axis1

     \param ev - Output - swapped components

  */
  ESCRIPT_DLL_API
  virtual void
  swapaxes(DataAbstract* ev, int axis0, int axis1);


  /**
     \brief
     solves the eigenvalue problem this*V=ev*V for the eigenvalues ev

     \param ev - Output - eigenvalues in increasing order at each data point

  */
  ESCRIPT_DLL_API
  virtual void
  eigenvalues(DataAbstract* ev);

  /**
     \brief
     solves the eigenvalue problem this*V=ev*V for the eigenvalues ev and eigenvectors V

     \param ev - Output - eigenvalues in increasing order at each data point
     \param V - Output - corresponding eigenvectors. They are normalized such that their length is one
                         and the first nonzero component is positive.
     \param tol - Input - eigenvalue with relative distance tol are treated as equal.

  */

  ESCRIPT_DLL_API
  virtual void
  eigenvalues_and_eigenvectors(DataAbstract* ev,DataAbstract* V,const double tol=1.e-13);


  /**
     \brief
     Return a reference to the underlying DataVector.
  */
  ESCRIPT_DLL_API
  DataTypes::ValueType&
  getVector();

  ESCRIPT_DLL_API
  const DataTypes::ValueType&
  getVector() const;


 protected:

 private:
  //
  // the actual data
  DataTypes::ValueType m_data;

};

inline
DataTypes::ValueType&
DataConstant::getVector()
{
  return m_data;
}

inline
const DataTypes::ValueType&
DataConstant::getVector() const
{
  return m_data;
}


} // end of namespace
#endif
