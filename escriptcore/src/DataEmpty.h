
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

#ifndef __ESCRIPT_DATAEMPTY_H__
#define __ESCRIPT_DATAEMPTY_H__

#include "system_dep.h"
#include "DataReady.h"

namespace escript {

/**
   \brief
   Implements the DataAbstract interface for an empty Data object.

   Description:
   Implements the DataAbstract interface for an empty Data object.
*/

class DataEmpty : public DataReady {
typedef DataReady parent;
 public:

  /**
     \brief
     Default constructor for DataEmpty.

     Description:
     Default constructor for DataEmpty.

  */
  ESCRIPT_DLL_API
  DataEmpty();

  /**
     \brief
     Destructor for DataEmpty.
  */
  ESCRIPT_DLL_API
  virtual
  ~DataEmpty();

  /**
     \brief
     Return a textual representation of the Data object.
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
     Return the offset for the given sample.
     NB: This will throw an exception as obviously an empty Data object contains no
     samples. An implementation is required by parent DataAbstract class.
     \param sampleNo - Input - Sample number.
     \param dataPointNo - Input - data-point number.
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
     As this is an empty Data object, this method will always return 0.
  */
  ESCRIPT_DLL_API
  virtual
  DataTypes::RealVectorType::size_type
  getLength() const;

  /**
     \brief
     Factory method that returns a newly created DataEmpty sliced from the
     current Data object according to the specified region.
     NB: This will throw an exception as obviously an empty Data object contains no
     data to slice from. An implementation is required by parent DataAbstract class.
  */
  ESCRIPT_DLL_API
  virtual
  DataAbstract*
  getSlice(const DataTypes::RegionType& region) const;

  /**
     \brief
     Set the current Data object according to the specified slice from the
     given input value.
     NB: This will throw an exception as obviously an empty Data object contains no
     data to slice to. An implementation is required by parent DataAbstract class.
     \param value Input - Data to copy from
     \param region Input - Region to copy.
  */
  ESCRIPT_DLL_API
  virtual
  void
  setSlice(const DataAbstract* value,
           const DataTypes::RegionType& region);

  /**
    \brief invert square matricies
    \param out - Where to store the results
    \return errorcode (0 indicates success)
  */
  ESCRIPT_DLL_API
  int
  matrixInverse(DataAbstract* out) const;

  void
  dump(const std::string fileName) const;

  ESCRIPT_DLL_API
  bool
  hasNaN() const
  {
	return false;
  }

  ESCRIPT_DLL_API
  void
  replaceNaN(DataTypes::real_t value)
  {
  
  }
  
  ESCRIPT_DLL_API
  void
  replaceNaN(DataTypes::cplx_t value)
  {
  
  }  
  
  /**
   \brief Return true if data contains Inf or -Inf 
  */
  ESCRIPT_DLL_API
  bool
  hasInf() const
  {
      return false;
  }

  /**
  \brief replaces all (+/-)Inf values with value 
  */
  ESCRIPT_DLL_API
  void
  replaceInf(DataTypes::real_t value)
  {
  }
  
  /**
  \brief replaces all (+/-)Inf values with value 
  */
  ESCRIPT_DLL_API
  void
  replaceInf(DataTypes::cplx_t value) 
  {
  }
  

 protected:

 /**
	\brief Provide access to underlying storage. Internal use only!
  */
  ESCRIPT_DLL_API
  virtual DataTypes::RealVectorType&
  getVectorRW();


  ESCRIPT_DLL_API
  virtual const DataTypes::RealVectorType&
  getVectorRO() const;

  ESCRIPT_DLL_API
  virtual DataTypes::CplxVectorType&
  getVectorRWC();


  ESCRIPT_DLL_API
  virtual const DataTypes::CplxVectorType&
  getVectorROC() const;
  
  virtual DataTypes::RealVectorType&
  getTypedVectorRW(DataTypes::real_t dummy);  
  
  virtual const DataTypes::RealVectorType&
  getTypedVectorRO(DataTypes::real_t dummy) const;

  virtual DataTypes::CplxVectorType&
  getTypedVectorRW(DataTypes::cplx_t dummy);
  
  virtual const DataTypes::CplxVectorType&
  getTypedVectorRO(DataTypes::cplx_t dummy) const;      
  
  
 private:

//  /**
/*     \brief
     Throw a standard exception. This function is called if an attempt
     is made to use functions of DataEmpty that are not valid.*/
//  */
//   void
//   throwStandardException(const std::string& functionName) const;

};

} // end of namespace

#endif // __ESCRIPT_DATAEMPTY_H__

