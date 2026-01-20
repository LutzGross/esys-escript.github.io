
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


#if !defined escript_DataExpanded_20040323_H
#define escript_DataExpanded_20040323_H
#include "system_dep.h"

#include "DataReady.h"

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

class DataExpanded : public DataReady {

typedef DataReady parent;

 public:

  /**
     \brief
     Constructor for DataExpanded.

     Description:
     Constructor for DataExpanded.

     The given single data value is copied to all the data points in 
     this data object, where the number of data points is defined by
     the given function space.

     \param value - Input - The value of a single data point.
     \param what - Input - A description of what this data represents.
  */
  ESCRIPT_DLL_API
  explicit DataExpanded(const WrappedArray& value,
               const FunctionSpace& what);

  /**
     \brief
     Alternative constructor for DataExpanded that copies a slice from
     another DataExpanded.

     \param other - Input - DataExpanded object to slice from.
     \param region - Input - region to copy.
  */
  ESCRIPT_DLL_API
  explicit DataExpanded(const DataExpanded& other,
               const DataTypes::RegionType& region);

  /**
     \brief
     Alternative constructor for DataExpanded objects.

     Description:
     Alternative Constructor for DataExpanded objects.
     \param what - Input - A description of what this data object represents.
     \param shape - Input - the shape of each data-point.
     \param data - the array of data values for the data-points.

TODO Note that this constructor will also copy data to all points if it only contains enough elements to hold a single point.  ie this is the merge of two separate constructors.
  */
  ESCRIPT_DLL_API
  explicit DataExpanded(const FunctionSpace& what,
               const DataTypes::ShapeType &shape,
               const DataTypes::RealVectorType &data);
  
  
  ESCRIPT_DLL_API
  explicit DataExpanded(const FunctionSpace& what,
               const DataTypes::ShapeType &shape,
               const DataTypes::CplxVectorType &data);
  

	       
  ESCRIPT_DLL_API
  explicit DataExpanded(const FunctionSpace& what,
               const DataTypes::ShapeType &shape,
               const DataTypes::real_t data);	       
  
  ESCRIPT_DLL_API
  explicit DataExpanded(const FunctionSpace& what,
               const DataTypes::ShapeType &shape,
               const DataTypes::cplx_t data);	       
  
	       
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
  explicit DataExpanded(const DataConstant& other);

  /**
     \brief
     Copy constructor for DataExpanded.
     Construct a DataExpanded from a DataTagged.
  */
  ESCRIPT_DLL_API
  explicit DataExpanded(const DataTagged& other);

  /**
     \brief
     Default destructor for DataExpanded.
  */
  ESCRIPT_DLL_API
  virtual
  ~DataExpanded();

  ESCRIPT_DLL_API
  bool
  isExpanded() const 
  {
    return true;
  };

  ESCRIPT_DLL_API
  bool
  actsExpanded() const
  {
    return true;
  }

  /**
  \brief Return true if any value in the data contains a NaN. 
  */
  ESCRIPT_DLL_API
  bool
  hasNaN() const;

  /**
  \brief replaces all NaN values with value 
  */
  ESCRIPT_DLL_API
  void
  replaceNaN(DataTypes::real_t value);

  ESCRIPT_DLL_API
  void
  replaceNaN(DataTypes::cplx_t value);
  
  /**
   \brief Return true if data contains Inf or -Inf 
  */
  ESCRIPT_DLL_API
  virtual bool
  hasInf() const;

  /**
  \brief replaces all (+/-)Inf values with value 
  */
  ESCRIPT_DLL_API
  virtual void
  replaceInf(DataTypes::real_t value);
  
  /**
  \brief replaces all (+/-)Inf values with value 
  */
  ESCRIPT_DLL_API
  virtual void
  replaceInf(DataTypes::cplx_t value);     
  
    
  /**
     \brief
     Return a textual representation of the data.
  */
  ESCRIPT_DLL_API
  virtual
  std::string
  toString() const;

  /**
     \brief Return a deep copy of the current object.
  */
  ESCRIPT_DLL_API
  virtual
  DataAbstract*
  deepCopy() const;

  /**
     \brief Return an object with the same type, domain (and tags if appropriate)
     as this, but all values are zeroed.
  */  
  ESCRIPT_DLL_API
  virtual
  DataAbstract*
  zeroedCopy() const;    
  
 /**
     \brief
     dumps the object into a HDF5 file
  */
 #ifdef ESYS_HAVE_HDF5
  ESCRIPT_DLL_API
  virtual
  void
  dump_hdf5(const H5::Group h5_grp) const;
#endif

  /**
    \brief invert square matricies
    \param out - Where to store the results
    \return errorcode (0 indicates success)
  */
  ESCRIPT_DLL_API
  virtual int
  matrixInverse(DataAbstract* out) const;

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
  DataTypes::RealVectorType::size_type
  getPointOffset(int sampleNo,
                 int dataPointNo) const;

//   ESCRIPT_DLL_API
//   virtual
//   DataTypes::RealVectorType::size_type
//   getPointOffset(int sampleNo,
//                  int dataPointNo);

  /**
     \brief
     Return a a reference to the underlying DataVector.
  */

  ESCRIPT_DLL_API
  DataTypes::RealVectorType&
  getVectorRW();

  ESCRIPT_DLL_API
  const DataTypes::RealVectorType&
  getVectorRO() const;

  ESCRIPT_DLL_API
  DataTypes::CplxVectorType&
  getVectorRWC();

  ESCRIPT_DLL_API
  const DataTypes::CplxVectorType&
  getVectorROC() const;
  
  virtual DataTypes::RealVectorType&
  getTypedVectorRW(DataTypes::real_t dummy);  
  
  virtual const DataTypes::RealVectorType&
  getTypedVectorRO(DataTypes::real_t dummy) const;

  virtual DataTypes::CplxVectorType&
  getTypedVectorRW(DataTypes::cplx_t dummy);
  
  virtual const DataTypes::CplxVectorType&
  getTypedVectorRO(DataTypes::cplx_t dummy) const;    

  /**
     \brief
     Return the number of doubles stored for the Data.
  */
  ESCRIPT_DLL_API
  virtual
  DataTypes::RealVectorType::size_type
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
     setTaggedValue

     Description:
     uses tag to set a new value

     \param tagKey - Input - Integer key.
     \param pointshape - Input - The shape of the value parameter
     \param value - Input - 
     \param dataOffset - Input - where in the value parameter to start reading the data point value.
  */
  void  
  setTaggedValue(int tagKey,
 	         const DataTypes::ShapeType& pointshape,
                 const DataTypes::RealVectorType& value,
		 int dataOffset=0);

  void  
  setTaggedValue(int tagKey,
 	         const DataTypes::ShapeType& pointshape,
                 const DataTypes::CplxVectorType& value,
		 int dataOffset=0);

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
     Computes a antisymmetric matrix (A - AT) / 2

     \param ev - Output - nonsymmetric matrix

  */
  ESCRIPT_DLL_API
  virtual void
  antisymmetric(DataAbstract* ev);

  /**
     \brief
     Computes an hermitian matrix (A + A*) / 2

     \param ev - Output - hermitian matrix

  */
  ESCRIPT_DLL_API
  virtual void
  hermitian(DataAbstract* ev);

  /**
     \brief
     Computes an antihermitian matrix (A - A*) / 2

     \param ev - Output - antihermitian matrix

  */
  ESCRIPT_DLL_API
  virtual void
  antihermitian(DataAbstract* ev);



  /**
     \brief
     Computes the trace of a matrix

     \param ev - Output - trace of your matrix
     \param axis_offset - 

  */
  ESCRIPT_DLL_API
  virtual void
  trace(DataAbstract* ev, int axis_offset);

  /**
     \brief
     Transpose each data point of this Data object around the given axis.

     \param ev - Output - transpose of your matrix
     \param axis_offset - 
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
  reorderByReferenceIDs(DataTypes::dim_t *reference_ids);

  ESCRIPT_DLL_API
  void
  complicate();
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

     \param noSamples - Input - number of samples.
     \param noDataPointsPerSample - Input - number of data points per sample.
     \param cplx - Input - is this data complex?
  */
  void
  initialise(int noSamples,
             int noDataPointsPerSample,
	     bool cplx
	    );

  /**
     \brief
     Copy the given data point value to all data points in this object.

     Description:
     Copy the given data point to all data points in this object.

     \param value Input - A single data point value.
  */
  void
  copy(const DataConstant& value);



  /**
     \brief
     Copy the given data point value to all data points in this object.

     \param value Input - A single data point value.
  */

  void
  copy(const WrappedArray& value);


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
  copyToDataPoint(const int sampleNo, const int dataPointNo, const DataTypes::real_t value);

  ESCRIPT_DLL_API
  virtual void
  copyToDataPoint(const int sampleNo, const int dataPointNo, const DataTypes::cplx_t value);  
  

  /**
     \brief
     Copy the value to the data point dataPointNo of sample sampleNo in this object.

     \param sampleNo Input - sample number
     \param dataPointNo Input - data point of the sample
     \param value Input - new values for the data point
  */
  ESCRIPT_DLL_API
  virtual void
  copyToDataPoint(const int sampleNo, const int dataPointNo, const WrappedArray& value);

  //
  // The main data storage array, a 2D array of data blocks.
  // noSamples * noDataPointsPerSample
  DataTypes::RealVectorType m_data_r;
  DataTypes::CplxVectorType m_data_c;
};

} // end of namespace

#endif
