
/*******************************************************
*
* Copyright (c) 2003-2009 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


#if !defined escript_DataAbstract_20040315_H
#define escript_DataAbstract_20040315_H
#include "system_dep.h"

#include "DataTypes.h"
#include "FunctionSpace.h"

#include <boost/scoped_ptr.hpp>

#include "DataException.h"

#include <string>
#include <fstream>
#include <vector>

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
   doubles of rank 0-4.
*/

class DataAbstract;

typedef POINTER_WRAPPER_CLASS(DataAbstract) DataAbstract_ptr;
typedef POINTER_WRAPPER_CLASS(const DataAbstract) const_DataAbstract_ptr;

class DataReady;

typedef POINTER_WRAPPER_CLASS(DataReady) DataReady_ptr;
typedef POINTER_WRAPPER_CLASS(const DataReady) const_DataReady_ptr;

class DataAbstract : public REFCOUNT_BASE_CLASS(DataAbstract)
{

 public:

  typedef DataTypes::ValueType ValueType;
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
   ESCRIPT_DLL_API
   DataAbstract_ptr getPtr();
   ESCRIPT_DLL_API
   const_DataAbstract_ptr getPtr() const; 



  /**
     \brief
     Constructor for DataAbstract.

     \param what - Input - The functionspace to use.
     \param shape - Input - Shape of each data value.
     \param isDataEmpty - Input - Is this an instance of DataEmpty (for internal use only)
  */
  ESCRIPT_DLL_API
  DataAbstract(const FunctionSpace& what, const ShapeType& shape, bool isDataEmpty=false);

  /**
    \brief
    Destructor for DataAbstract.
  */
  ESCRIPT_DLL_API
  virtual
  ~DataAbstract();

  /**
     \brief
     Write the data as a string.
  */
  ESCRIPT_DLL_API
  virtual
  std::string
  toString() const = 0;

  /**
     \brief Return a deep copy of the current object.
  */
  ESCRIPT_DLL_API
  virtual
  DataAbstract*
  deepCopy()=0;

  /**
     \brief Return a data object with all points resolved.
  */
  ESCRIPT_DLL_API
  virtual
  DataReady_ptr
  resolve()=0;

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
     Return the number of data points per sample.
  */
  ESCRIPT_DLL_API
  int
  getNumDPPSample() const;

  /**
     \brief
     Return the number of samples.
  */
  ESCRIPT_DLL_API
  int
  getNumSamples() const;

  /**
     \brief 
     Return the shape information for the point data.

     The omission of a non-constant form is deliberate.
  */
  ESCRIPT_DLL_API
  const DataTypes::ShapeType& 
  getShape() const;

  /**
     \brief 
     Return the rank information for the point data.
  */
  ESCRIPT_DLL_API
  unsigned int 
  getRank() const;



 /**
    \brief
    Return the offset for the given sample. This returns the offset for the given
    point into the container holding the point data. 

    \param sampleNo - Input - sample number.
    \param dataPointNo - Input - data point number.
  */
  ESCRIPT_DLL_API
  virtual
  ValueType::size_type
  getPointOffset(int sampleNo,
                 int dataPointNo) const = 0;

  ESCRIPT_DLL_API
  virtual
  ValueType::size_type
  getPointOffset(int sampleNo,
                 int dataPointNo) = 0;


  /**
     \brief
     Return the number of doubles stored for this Data object.
  */
  ESCRIPT_DLL_API
  virtual
  ValueType::size_type
  getLength() const = 0;

  /**
     \brief
     Return the sample data for the given tag key.
     NB: If the data isn't tagged an exception will be thrown.
  */
  ESCRIPT_DLL_API
  virtual
  double*
  getSampleDataByTag(int tag);

  /**
    This method is used primarily for LazyData.
    \return the size of the buffer required to evaulate a sample for this object.
  */
  ESCRIPT_DLL_API
  virtual size_t
  getSampleBufferSize() const=0;



  /**
     \brief
     Check this and the given RHS operands are compatible. Throws
     an exception if they aren't.

     \param right - Input - The right hand side.
  */
  ESCRIPT_DLL_API
  void
  operandCheck(const DataAbstract& right) const;

  /**
     \brief
     Return true if a valid sample point number.
  */
  ESCRIPT_DLL_API
  bool
  validSamplePointNo(int samplePointNo) const;

  /**
     \brief
     Return true if a valid sample number.
  */
  ESCRIPT_DLL_API
  bool
  validSampleNo(int sampleNo) const;


  /**
     \brief
     Return the function space associated with this Data object.
  */
  ESCRIPT_DLL_API
  const
  FunctionSpace&
  getFunctionSpace() const;

  /**
     \brief
     Return the given slice from this object.

     NB: The caller is responsible for managing the object created.
  */
  ESCRIPT_DLL_API
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
  ESCRIPT_DLL_API
  virtual
  void
  setTaggedValue(int tagKey,
		 const DataTypes::ShapeType& pointshape,
                 const DataTypes::ValueType& value,
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
  ESCRIPT_DLL_API
  virtual void
  copyToDataPoint(const int sampleNo, const int dataPointNo, const double value);

  /**
     \brief
     Copy the array object to the data point dataPointNo of sample sampleNo in this object.

     \param sampleNo Input - sample number
     \param dataPointNo Input - data point of the sample
     \param value Input - new values for the data point
  */
  ESCRIPT_DLL_API
  virtual void
  copyToDataPoint(const int sampleNo, const int dataPointNo, const WrappedArray& value);


  /**
     \brief
     Return the tag number associated with the given data-point number.

     If the object cannot be referenced by tag numbers, an exception
     will be thrown.
  */
  ESCRIPT_DLL_API
  virtual
  int
  getTagNumber(int dpno);

  /**
     \brief
     Computes a symmetric matrix (A + AT) / 2

     \param ev - Output - a symmetric matrix

  */
  ESCRIPT_DLL_API
  virtual void
  symmetric(DataAbstract* ev);

  /**
     \brief
     Computes a nonsymmetric matrix (A - AT) / 2

     \param ev - Output - a nonsymmetric matrix

  */
  ESCRIPT_DLL_API
  virtual void
  nonsymmetric(DataAbstract* ev);

  /**
     \brief
     Computes the trace of a matrix

     \param ev - Output - the trace of a matrix
     \param axis_offset
  */
  ESCRIPT_DLL_API
  virtual void
  trace(DataAbstract* ev, int axis_offset);

  /**
     \brief
     Transpose each data point of this Data object around the given axis.

     \param ev - Output - the transpose of a matrix
     \param axis_offset
  */
  ESCRIPT_DLL_API
  virtual void
  transpose(DataAbstract* ev, int axis_offset);

  /**
     \brief
     swaps components axis0 and axis1

     \param ev - Output - swapped components
     \param axis0
     \param axis1
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
     sets values to zero

  */
  ESCRIPT_DLL_API
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

  ESCRIPT_DLL_API
  virtual void
  eigenvalues_and_eigenvectors(DataAbstract* ev,DataAbstract* V,const double tol=1.e-13);

  /**
     \brief
     reorders data sample ordered by reference_ids to the ordering of the functions space

     \param reference_ids - Input - reference_ids used for current ordering
  */
  ESCRIPT_DLL_API
  virtual void
  reorderByReferenceIDs(int *reference_ids);



  /**
	\brief
	Return the number of values in the shape for this object.
  */
  ESCRIPT_DLL_API
  unsigned int
  getNoValues() const;


  ESCRIPT_DLL_API
  bool isLazy() const;	// a test to determine if this object is an instance of DataLazy

  ESCRIPT_DLL_API
  virtual
  bool
  isConstant() const {return false;}

  ESCRIPT_DLL_API
  virtual
  bool
  isExpanded() const {return false;}


  /**
     \brief
     Return true if this Data is expanded or resolves to expanded.
     That is, if it has a separate value for each datapoint in the sample.
  */
  ESCRIPT_DLL_API
  virtual
  bool
  actsExpanded() const {return false;}

  ESCRIPT_DLL_API
  virtual
  bool
  isTagged() const {return false;}

  ESCRIPT_DLL_API
  bool isEmpty() const;	// a fast test to determine if this object is an instance of DataEmpty


  /**
  	\warning should only be used in single threaded code (or inside a single/critical section)
  */
  void
  addOwner(Data*);

  /**
  	\warning should only be used in single threaded code (or inside a single/critical section)
  */
  void
  removeOwner(Data*);

  /**
	\brief Is this object owned by more than one Data object
  */
  ESCRIPT_DLL_API
  bool
  isShared() const
  {
	return m_lazyshared || (m_owners.size()>1);
  }

 protected:

   /**
   \brief Returns true if this object is not shared.
   For internal use only. - It may not be particularly fast
   */
   ESCRIPT_DLL_API
   bool checkNoSharing() const;

   /**
   \brief Marks this DataAbstract shared as LazyData
   For internal use only.
   */
   void
   makeLazyShared();	

   friend class DataLazy;

 private:

  //
  // The number of samples in this Data object.
  // This is derived directly from the FunctionSpace.
  int m_noSamples;

  //
  // The number of data points per sample in this Data object.
  // This is derived directly from the FunctionSpace.
  int m_noDataPointsPerSample;

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

public:			// these should be private once I have finished debugging
  std::vector<Data*> m_owners;
  bool m_lazyshared;
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
     	throw DataException("Error - Operations not permitted on instances of DataEmpty.");
  }
  return m_noDataPointsPerSample;
}

inline
int
DataAbstract::getNumSamples() const
{
  if (isEmpty())
  {
     	throw DataException("Error - Operations not permitted on instances of DataEmpty.");
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
		throw DataException("Error - Operations not permitted on instances of DataEmpty.");
	}
	return m_shape;
}

inline 
unsigned int
DataAbstract::getRank() const
{
	if (isEmpty())
	{
		throw DataException("Error - Operations not permitted on instances of DataEmpty.");
	}
	return m_rank;
}

inline
unsigned int
DataAbstract::getNoValues() const
{	
	if (isEmpty())
	{
		throw DataException("Error - Operations not permitted on instances of DataEmpty.");
	}
	return m_novalues;
}


} // end of namespace

#endif
