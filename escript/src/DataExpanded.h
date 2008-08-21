
/* $Id$ */

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

#if !defined escript_DataExpanded_20040323_H
#define escript_DataExpanded_20040323_H
#include "system_dep.h"

#include "DataAbstract.h"
#include "DataBlocks2D.h"
#include "DataArrayView.h"

#include <boost/python/numeric.hpp>

namespace escript {

//
// Forward declarations of other Data types.
class DataConstant;
class DataTagged;

/**
   \brief
   Give a short description of what DataExpanded does.

   Description:
   Give a detailed description of DataExpanded.

   Template Parameters:
   For templates describe any conditions that the parameters used in the
   template must satisfy.
*/

class DataExpanded : public DataAbstract {

 public:

  /**
     \brief
     Constructor for DataExpanded.

     Description:
     Constructor for DataExpanded.

     The given single data value is copied to all the data points in 
     this data object, where the number of data points is defined by
     the given function space.

     \param value - Input - A single data value.
     \param what - Input - A description of what this data represents.
  */
  ESCRIPT_DLL_API
  DataExpanded(const boost::python::numeric::array& value,
               const FunctionSpace& what);

  /**
     \brief
     Alternative constructor for DataExpanded.

     Description:
     Alternative Constructor for DataExpanded.

     The given single data value is copied to all the data points in 
     this data object, where the number of data points is defined by
     the given function space.

     \param value - Input - A single data value.
     \param what - Input - A description of what this data represents.
  */
  ESCRIPT_DLL_API
  DataExpanded(const DataArrayView& value,
               const FunctionSpace& what);

  /**
     \brief
     Alternative constructor for DataExpanded that copies a slice from
     another DataExpanded.

     \param other - Input - DataExpanded object to slice from.
     \param region - Input - region to copy.
  */
  ESCRIPT_DLL_API
  DataExpanded(const DataExpanded& other,
               const DataTypes::RegionType& region);

  /**
     \brief
     Alternative constructor for DataExpanded objects.

     Description:
     Alternative Constructor for DataExpanded objects.
     \param what - Input - A description of what this data object represents.
     \param shape - Input - the shape of each data-point.
     \param data - the array of data values for the data-points.
  */
  ESCRIPT_DLL_API
  DataExpanded(const FunctionSpace& what,
               const DataTypes::ShapeType &shape,
               const DataTypes::ValueType &data);

  /**
     \brief
     Copy constructor for DataExpanded.
     Performs a deep copy from another DataExpanded.
  */
  ESCRIPT_DLL_API
  DataExpanded(const DataExpanded& other);

  /**
     \brief
     Copy constructor for DataExpanded.
     Construct a DataExpanded from a DataConstant.
  */
  ESCRIPT_DLL_API
  DataExpanded(const DataConstant& other);

  /**
     \brief
     Copy constructor for DataExpanded.
     Construct a DataExpanded from a DataTagged.
  */
  ESCRIPT_DLL_API
  DataExpanded(const DataTagged& other);

  /**
     \brief
     Default destructor for DataExpanded.
  */
  ESCRIPT_DLL_API
  virtual
  ~DataExpanded();

  /**
     \brief
     Return a textual representation of the data.
  */
  ESCRIPT_DLL_API
  virtual
  std::string
  toString() const;
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
     Return the offset for the given given data point. This returns
     the offset in bytes for the given point into the container
     holding the point data.

     \param sampleNo - Input - sample number.
     \param dataPointNo - Input - data point number.
  */
  ESCRIPT_DLL_API
  virtual
  DataTypes::ValueType::size_type
  getPointOffset(int sampleNo,
                 int dataPointNo) const;

  /**
     \brief
     Return a view into the data array for the data point specified.

     NOTE: Construction of the DataArrayView is a relatively expensive 
     operation.

     \param sampleNo - Input - sample number.
     \param dataPointNo - Input - data point number.
     \return DataArrayView for the data point.
  */
  ESCRIPT_DLL_API
  DataArrayView
  getDataPoint(int sampleNo,
               int dataPointNo);

  /**
     \brief
     Return the number of doubles stored for the Data.
  */
  ESCRIPT_DLL_API
  virtual
  ValueType::size_type
  getLength() const;

  /**
     \brief
     Factory method that returns a newly created DataExpanded.
     The caller is reponsible for managing the object created.

     \param region - Input - Region to copy.
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
    Archive the underlying data values to the file referenced
    by ofstream. A count of the number of values expected to be written
    is provided as a cross-check.

    The return value indicates success (0) or otherwise (1).
  */
  ESCRIPT_DLL_API
  int
  archiveData(std::ofstream& archiveFile,
              const DataTypes::ValueType::size_type noValues) const;

  /**
    \brief
    Extract the number of values specified by noValues from the file
    referenced by ifstream to the underlying data structure.

    The return value indicates success (0) or otherwise (1).
  */
  ESCRIPT_DLL_API
  int
  extractData(std::ifstream& archiveFile,
              const DataTypes::ValueType::size_type noValues);

  /**
     \brief
     setTaggedValue

     Description:
     uses tag to set a new value

     \param tagKey - Input - Integer key.
     \param value - Input - Single DataArrayView value to be assigned to the tag.
  */
  ESCRIPT_DLL_API
  virtual
  void
  setTaggedValue(int tagKey,
                 const DataArrayView& value);



  /**
     \brief
     setTaggedValue

     Description:
     uses tag to set a new value

     \param tagKey - Input - Integer key.
     \param pointshape - Input - The shape of the value parameter
     \param value - Input - .
  */
  void  
  setTaggedValue(int tagKey,
 	         const DataTypes::ShapeType& pointshape,
                 const DataTypes::ValueType& value);



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

     \param ev - Output - trace of your matrix

  */
  ESCRIPT_DLL_API
  virtual void
  trace(DataAbstract* ev, int axis_offset);

  /**
     \brief
     Transpose each data point of this Data object around the given axis.

     \param ev - Output - transpose of your matrix

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
 *      \brief
 *           reorders data sample ordered by reference_ids to the ordering of the functions space
 *
 *                \param reference_ids - Input - reference_ids used for current ordering
 *                  */
  ESCRIPT_DLL_API
  virtual void
  reorderByReferenceIDs(int *reference_ids);



 protected:

 private:

  /**
     \brief
     Common initialisation called from constructors.

     Description:
     Common initialisation called from constructors.

     Resizes the underlying data array to provide sufficient storage for the
     given shape and number of data points, and creates the corresponding
     DataArrayView of this data.

     \param shape - Input - The shape of the point data.
     \param noSamples - Input - number of samples.
     \param noDataPointsPerSample - Input - number of data points per sample.
  */
  void
  initialise(const DataTypes::ShapeType& shape,
             int noSamples,
             int noDataPointsPerSample);

  /**
     \brief
     Copy the given data point value to all data points in this object.

     Description:
     Copy the given data point to all data points in this object.

     \param value Input - A single data point value.
  */
  void
  copy(const DataArrayView& value);

  /**
     \brief
     Copy the given data point value given a numarray object to all data points in this object.

     Description:
     Copy the given data point value given a numarray object to all data points in this object.

     \param value Input - A single data point value.
  */
  void
  copy(const boost::python::numeric::array& value);

  /**
     \brief
     Copy the numarray object to the data points in this object.

     Description:
     Copy the numarray object to the data points in this object.

     \param value Input - new values for the data points
  */
  void
  copyAll(const boost::python::numeric::array& value);

  /**
     \brief
     Copy a double value to the data point dataPointNo of sample sampleNo in this object.

     Description:
     Copy a double value to the data point dataPointNo of sample sampleNo in this object.

     \param sampleNo Input - sample number
     \param dataPointNo Input - data point of the sample
     \param value Input - new values for the data point
  */
  ESCRIPT_DLL_API
  virtual void
  copyToDataPoint(const int sampleNo, const int dataPointNo, const double value);


  /**
     \brief
     Copy the numarray object to the data point dataPointNo of sample sampleNo in this object.

     Description:
     Copy the numarray object to the data point dataPointNo of sample sampleNo in this object.

     \param sampleNo Input - sample number
     \param dataPointNo Input - data point of the sample
     \param value Input - new values for the data point
  */
  void
  copyToDataPoint(const int sampleNo, const int dataPointNo, const boost::python::numeric::array& value);

  //
  // The main data storage array, a 2D array of data blocks.
  // noSamples * noDataPointsPerSample
  DataBlocks2D m_data;

};

} // end of namespace

#endif
