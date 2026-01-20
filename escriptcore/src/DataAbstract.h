
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/

#ifndef __ESCRIPT_DATAABSTRACT_H__
#define __ESCRIPT_DATAABSTRACT_H__

#include "system_dep.h"
#include "DataTypes.h"
#include "DataVector.h"
#include "FunctionSpace.h"

#include <boost/scoped_ptr.hpp>

#include "DataException.h"

#include <string>
#include <fstream>
#include <vector>

#ifdef ESYS_HAVE_HDF5
  #include <H5Cpp.h>
#endif

#include "Pointers.h"

namespace escript {

/**
   \brief
   DataAbstract provides an abstract interface for the class of containers
   which hold ESyS data.

   Description:
   DataAbstract provides an abstract interface for the class of containers
   which hold ESyS data. The container may be thought of as a 2 dimensional
   array of data points where one dimension corresponds to the number of samples
   and the other to the number of data points per sample as defined by the function
   space associated with each Data object. The data points themselves are arrays of
   reals or complexes of rank 0-4.
*/

class DataAbstract;

typedef POINTER_WRAPPER_CLASS(DataAbstract) DataAbstract_ptr;
typedef POINTER_WRAPPER_CLASS(const DataAbstract) const_DataAbstract_ptr;

class DataReady;

typedef POINTER_WRAPPER_CLASS(DataReady) DataReady_ptr;
typedef POINTER_WRAPPER_CLASS(const DataReady) const_DataReady_ptr;

class ESCRIPT_DLL_API DataAbstract : public REFCOUNT_BASE_CLASS(DataAbstract)
{

 public:

  typedef DataTypes::ShapeType ShapeType;

   /**
   \brief Return shared pointer managing this object.

   If there is not already a shared pointer managing this object then create one.
   Once a shared pointer is created for an object, the deallocation of the object
   must be handled by shared_ptr.

   \warning So, do not call this on an automatic object.
   Do not call this in a method where you do not pass the shared_pointer out and
   you need the object to outlast the method.

   Note: This is _not_ equivalent to weak_ptr::lock.

   */
   DataAbstract_ptr getPtr();
   const_DataAbstract_ptr getPtr() const;



  /**
     \brief
     Constructor for DataAbstract.

     \param what - Input - The functionspace to use.
     \param shape - Input - Shape of each data value.
     \param isDataEmpty - Input - Is this an instance of DataEmpty (for internal use only)
  */
  DataAbstract(const FunctionSpace& what, const ShapeType& shape, bool isDataEmpty=false,bool isCplx=false);

  /**
    \brief
    Destructor for DataAbstract.
  */
  virtual
  ~DataAbstract();

  /**
     \brief
     Write the data as a string.
  */
  virtual
  std::string
  toString() const = 0;

  /**
     \brief Return a deep copy of the current object.
  */
  virtual
  DataAbstract*
  deepCopy() const =0 ;
  
  /**
     \brief Return an object with the same type, domain (and tags if appropriate)
     as this, but all values are zeroed.
  */
  virtual
  DataAbstract*
  zeroedCopy() const =0 ;
  
  

  /**
     \brief Return a data object with all points resolved.
  */
  virtual
  DataReady_ptr
  resolve()=0;

 /**
     \brief
     dumps the object into a HDF5 file
  */
#ifdef ESYS_HAVE_HDF5
  virtual
  void
  dump_hdf5(const H5::Group h5_grp) const;
#endif


  /**
     \brief
     Return the number of data points per sample.
  */
  int
  getNumDPPSample() const;

  /**
     \brief
     Return the number of samples.
  */
  int
  getNumSamples() const;
  
  bool
  hasNoSamples() const
  {
      return getNumSamples()==0;
  }

  /**
     \brief
     Return the shape information for the point data.

     The omission of a non-constant form is deliberate.
  */
  const DataTypes::ShapeType&
  getShape() const;

  /**
     \brief
     Return the rank information for the point data.
  */
  unsigned int
  getRank() const;



 /**
    \brief
    Return the offset for the given sample. This returns the offset for the given
    point into the container holding the point data.

    \param sampleNo - Input - sample number.
    \param dataPointNo - Input - data point number.
  */
  virtual
  DataTypes::RealVectorType::size_type
  getPointOffset(int sampleNo,
                 int dataPointNo) const = 0;



  /**
     \brief
     Return the number of doubles stored for this Data object.
  */
  virtual
  DataTypes::RealVectorType::size_type
  getLength() const = 0;

  /**
     \brief
     Return the real sample data for the given tag key.
     NB: If the data isn't tagged an exception will be thrown.
  */
  virtual
  DataTypes::real_t*
  getSampleDataByTag(int tag, DataTypes::real_t dummy=0);

  /**
     \brief
     Return the complex sample data for the given tag key.
     NB: If the data isn't tagged an exception will be thrown.
  */
  virtual
  DataTypes::cplx_t*
  getSampleDataByTag(int tag, DataTypes::cplx_t dummy);

  /**
     \brief Return number of tagged values stored in the data object
     \warning results are only meaningful for DataTagged. All other types return 0.
     This functionality is only currently used by reducers and should
     not be exposed to Python without making it more generally applicable
  */
  virtual
  size_t
  getTagCount() const;

  /**
     \brief
     Check this and the given RHS operands are compatible. Throws
     an exception if they aren't.

     \param right - Input - The right hand side.
  */
  void
  operandCheck(const DataAbstract& right) const;

  /**
     \brief
     Return true if a valid sample point number.
  */
  bool
  validSamplePointNo(int samplePointNo) const;

  /**
     \brief
     Return true if a valid sample number.
  */
  bool
  validSampleNo(int sampleNo) const;


  /**
     \brief
     Return the function space associated with this Data object.
  */
  const
  FunctionSpace&
  getFunctionSpace() const;

  /**
     \brief
     Return the given slice from this object.

     NB: The caller is responsible for managing the object created.
  */
  virtual
  DataAbstract*
  getSlice(const DataTypes::RegionType& region) const = 0;



  /**
     \brief
     setTaggedValue

     Description:
     Assign the given value to the given tag.

     NB: If the data isn't tagged an exception will be thrown.

     \param tagKey - Input - Integer key.
     \param pointshape - Input - the shape of the value parameter.
     \param value - Input - vector to copy data value from
     \param dataOffset - Input - Offset within value to begin copying from

     The final parameter is to allow for the case whete the vector contains
     multiple data values.
  */
  virtual
  void
  setTaggedValue(int tagKey,
                 const DataTypes::ShapeType& pointshape,
                 const DataTypes::RealVectorType& value,
                 int dataOffset=0);

  virtual
  void
  setTaggedValue(int tagKey,
                 const DataTypes::ShapeType& pointshape,
                 const DataTypes::CplxVectorType& value,
                 int dataOffset=0);

  /**
     \brief
     Copy a double value to the data point dataPointNo of sample sampleNo in this object.

     Description:
     Copy a double value to the data point dataPointNo of sample sampleNo in this object.

     \param sampleNo Input - sample number
     \param dataPointNo Input - data point of the sample
     \param value Input - new values for the data point
  */
  virtual void
  copyToDataPoint(const int sampleNo, const int dataPointNo, const DataTypes::real_t value);

  virtual void
  copyToDataPoint(const int sampleNo, const int dataPointNo, const DataTypes::cplx_t value);

  /**
     \brief
     Copy the array object to the data point dataPointNo of sample sampleNo in this object.

     \param sampleNo Input - sample number
     \param dataPointNo Input - data point of the sample
     \param value Input - new values for the data point
  */
  virtual void
  copyToDataPoint(const int sampleNo, const int dataPointNo, const WrappedArray& value);


  /**
     \brief
     Return the tag number associated with the given data-point number.

     If the object cannot be referenced by tag numbers, an exception
     will be thrown.
  */
  virtual
  int
  getTagNumber(int dpno);

  /**
     \brief
     Computes a symmetric matrix (A + AT) / 2

     \param ev - Output - a symmetric matrix

  */
  virtual void
  symmetric(DataAbstract* ev);

  /**
     \brief
     Computes a antisymmetric matrix (A - AT) / 2

     \param ev - Output - a nonsymmetric matrix

  */
  virtual void
  antisymmetric(DataAbstract* ev);

  /**
     \brief
     Computes a symmetric matrix (A + A*) / 2

     \param ev - Output - an hermitian matrix

  */
  virtual void
  hermitian(DataAbstract* ev);

  /**
     \brief
     Computes a antisymmetric matrix (A - A*) / 2

     \param ev - Output - an antihermitian matrix

  */
  virtual void
  antihermitian(DataAbstract* ev);

  /**
     \brief
     Computes the trace of a matrix

     \param ev - Output - the trace of a matrix
     \param axis_offset
  */
  virtual void
  trace(DataAbstract* ev, int axis_offset);

  /**
     \brief
     Transpose each data point of this Data object around the given axis.

     \param ev - Output - the transpose of a matrix
     \param axis_offset
  */
  virtual void
  transpose(DataAbstract* ev, int axis_offset);

  /**
     \brief
     swaps components axis0 and axis1

     \param ev - Output - swapped components
     \param axis0
     \param axis1
  */
  virtual void
  swapaxes(DataAbstract* ev, int axis0, int axis1);
  /**
     \brief
     solves the eigenvalue problem this*V=ev*V for the eigenvalues ev

     \param ev - Output - eigenvalues in increasing order at each data point

  */
  virtual void
  eigenvalues(DataAbstract* ev);

  /**
    \brief invert square matricies
    \param out - Where to store the results
    \return errorcode (0 indicates success)
  */
  virtual int
  matrixInverse(DataAbstract* out) const;

  /**
     \brief
     sets values to zero

  */
  virtual void
  setToZero();

  /**
     \brief
     solves the eigenvalue problem this*V=ev*V for the eigenvalues ev and eigenvectors V

     \param ev - Output - eigenvalues in increasing order at each data point
     \param V - Output - corresponding eigenvectors. They are normalized such that their length is one
                         and the first nonzero component is positive.
     \param tol - Input - eigenvalue with relative distance tol are treated as equal.

  */

  virtual void
  eigenvalues_and_eigenvectors(DataAbstract* ev,DataAbstract* V,const double tol=1.e-13);

  /**
     \brief
     reorders data sample ordered by reference_ids to the ordering of the functions space

     \param reference_ids - Input - reference_ids used for current ordering
  */
  virtual void
  reorderByReferenceIDs(DataTypes::dim_t *reference_ids);



  /**
        \brief
        Return the number of values in the shape for this object.
  */
  unsigned int
  getNoValues() const;


  bool isLazy() const;  // a test to determine if this object is an instance of DataLazy

  virtual
  bool
  isConstant() const {return false;}

  virtual
  bool
  isExpanded() const {return false;}


  /**
     \brief
     Return true if this Data is expanded or resolves to expanded.
     That is, if it has a separate value for each datapoint in the sample.
  */
  virtual
  bool
  actsExpanded() const {return false;}

  virtual
  bool
  isTagged() const {return false;}

  bool isEmpty() const; // a fast test to determine if this object is an instance of DataEmpty

  /**
   \brief true if the components of datapoints are complex
  */
  bool isComplex() const;

#ifdef SLOWSHARECHECK   
  
  // For this to be threadsafe, we need to be sure that this is the
  // only way shared-ness is tested.
  /**
        \brief Is this object owned by more than one Data object
  */
  bool
  isShared() const
  {
    bool shared=false;
    #pragma omp critical        // because two treads could try
    {                   // this check at the same time
      try               // and shared_from_this increments count
      {
        shared=shared_from_this().use_count()>2;
      }
      catch (...)
      {
      }
    }
    return shared;
  }
#endif    

#ifdef EXWRITECHK
  bool exclusivewritecalled;    // used to check for some potential programming faults
                                // involving shared data.
                                // This flag only asserts that exclusive write has been called
                                // on this object, it does not definitively guarantee that
                                // sharing has not occurred since that call
                                // This flag is for internal use only may be removed without warning
#endif

/*
 * Make the object complex
*/
 virtual void complicate();

protected:
    friend class DataLazy;

  //
  // The number of samples in this Data object.
  // This is derived directly from the FunctionSpace.
  int m_noSamples;

  //
  // The number of data points per sample in this Data object.
  // This is derived directly from the FunctionSpace.
  int m_noDataPointsPerSample;

  //
  // is the data made of complex components
  bool m_iscompl;
private:

  //
  // A FunctionSpace which provides a description of the data associated
  // with this Data object.
  FunctionSpace m_functionSpace;

  //
  // The shape of the points stored in this view
  DataTypes::ShapeType m_shape;

  //
  // The number of values in each point
  unsigned int m_novalues;

  //
  // The rank of the points stored in this view
  unsigned int m_rank;

  //
  // Is this an instance of DataEmpty?
  bool m_isempty;
};

inline
bool
DataAbstract::isEmpty() const
{
        return m_isempty;
}

inline
bool
DataAbstract::validSamplePointNo(int samplePointNo) const
{
  return ((0 <= samplePointNo) && (samplePointNo < m_noDataPointsPerSample));
}

inline
bool
DataAbstract::validSampleNo(int sampleNo) const
{
  return ((0 <= sampleNo) && (sampleNo < m_noSamples));
}

inline
int
DataAbstract::getNumDPPSample() const
{
  if (isEmpty())
  {
        throw DataException("Error - Operations (getNumDPPSample) not permitted on instances of DataEmpty.");
  }
  return m_noDataPointsPerSample;
}

inline
int
DataAbstract::getNumSamples() const
{
  if (isEmpty())
  {
        throw DataException("Error - Operations (getNumSamples) not permitted on instances of DataEmpty.");
  }
  return m_noSamples;
}

inline
const
FunctionSpace&
DataAbstract::getFunctionSpace() const
{
  return m_functionSpace;
}

inline
const DataTypes::ShapeType&
DataAbstract::getShape() const
{
        if (isEmpty())
        {
                throw DataException("Error - Operations (getShape) not permitted on instances of DataEmpty.");
        }
        return m_shape;
}

inline
unsigned int
DataAbstract::getRank() const
{
        if (isEmpty())
        {
                throw DataException("Error - Operations (getRank) not permitted on instances of DataEmpty.");
        }
        return m_rank;
}

inline
unsigned int
DataAbstract::getNoValues() const
{
        if (isEmpty())
        {
                throw DataException("Error - Operations (getNoValues) not permitted on instances of DataEmpty.");
        }
        return m_novalues;
}

} // end of namespace

#endif // __ESCRIPT_DATAABSTRACT_H__

