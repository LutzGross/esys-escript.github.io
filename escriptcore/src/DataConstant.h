
/*****************************************************************************
*
* Copyright (c) 2003-2020 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/


#if !defined escript_DataConstant_20040323_H
#define escript_DataConstant_20040323_H
#include "system_dep.h"

#include "DataReady.h"
#include "WrappedArray.h"
#ifdef ESYS_HAVE_HDF5
  #include <H5Cpp.h>
#endif


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
  explicit DataConstant(const WrappedArray& value,
               const FunctionSpace& what);


  /**
     \brief
     Copy constructor. Performs a deep copy.
  */
  ESCRIPT_DLL_API
  DataConstant(const DataConstant& other);


  /**
     \brief
     Alternative constructor for DataConstant objects.

     Description:
     Alternative Constructor for DataConstant objects.
     \param other - Input - Data object to copy from.
     \param region - Input - region to copy.
  */
  ESCRIPT_DLL_API
  explicit DataConstant(const DataConstant& other,
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
  explicit DataConstant(const FunctionSpace& what,
               const DataTypes::ShapeType &shape,
               const DataTypes::RealVectorType &data);
  
  explicit DataConstant(const FunctionSpace& what,
               const DataTypes::ShapeType &shape,
               const DataTypes::CplxVectorType &data);    

  ESCRIPT_DLL_API
  explicit DataConstant(const FunctionSpace& what,
                           const DataTypes::ShapeType &shape,
                           const DataTypes::real_t v);
               
  explicit DataConstant(const FunctionSpace& what,
                           const DataTypes::ShapeType &shape,
                           const DataTypes::cplx_t v);
               
  ESCRIPT_DLL_API
  bool
  isConstant() const 
  {
    return true;
  };

  /**
  \brief Return true if the value contains a NaN. 
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
     dumps the object into an HDF5 file
  */
 #ifdef ESYS_HAVE_HDF5
  ESCRIPT_DLL_API
  virtual
  void
  dump_hdf5(const H5::Group h5_grp) const;
#endif
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
  DataTypes::RealVectorType::size_type
  getPointOffset(int sampleNo,
                 int dataPointNo) const;

  ESCRIPT_DLL_API
  virtual
  DataTypes::RealVectorType::size_type
  getPointOffset(int sampleNo,
                 int dataPointNo);

  /**
     \brief
     Return the number of doubles stored for the Data object.
  */
  ESCRIPT_DLL_API
  virtual
  DataTypes::RealVectorType::size_type
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
     Computes an anti-hermitian matrix (A - A*) / 2

     \param ev - Output - antihermitian matrix

  */
  ESCRIPT_DLL_API
  virtual void
  antihermitian(DataAbstract* ev);

  /**
     \brief
     Computes the trace of a matrix

     \param ev - Output - trace of matrix
     \param axis_offset

  */
  ESCRIPT_DLL_API
  virtual void
  trace(DataAbstract* ev, int axis_offset);

  /**
     \brief
     Transpose each data point of this Data object around the given axis.

     \param ev - Output - transpose of matrix
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
    \brief invert square matricies
    \param out - Where to store the results
    \return errorcode (0 indicates success)
  */
  ESCRIPT_DLL_API
  virtual int
  matrixInverse(DataAbstract* out) const;

  /**
     \brief
     Return a reference to the underlying DataVector.
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



  ESCRIPT_DLL_API
  virtual DataTypes::RealVectorType&
  getTypedVectorRW(DataTypes::real_t dummy);  
  
  ESCRIPT_DLL_API
  virtual const DataTypes::RealVectorType&
  getTypedVectorRO(DataTypes::real_t dummy) const;

  ESCRIPT_DLL_API
  virtual DataTypes::CplxVectorType&
  getTypedVectorRW(DataTypes::cplx_t dummy);
  
  ESCRIPT_DLL_API
  virtual const DataTypes::CplxVectorType&
  getTypedVectorRO(DataTypes::cplx_t dummy) const;  



  
  /**
   * \brief Convert from real data to complex data.
  */ 
  ESCRIPT_DLL_API
  void complicate();

 protected:

 private:
  //
  // the actual data
  DataTypes::RealVectorType m_data_r;
  DataTypes::CplxVectorType m_data_c;

};

} // end of namespace
#endif
