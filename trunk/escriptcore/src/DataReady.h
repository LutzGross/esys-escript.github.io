
/*****************************************************************************
*
* Copyright (c) 2003-2016 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/


#if !defined escript_DataReady_20081008_H
#define escript_DataReady_20081008_H
#include "system_dep.h"

#include "DataAbstract.h"

namespace escript {

// Anything which requires getVector should be moved down here



/**
  \class DataReady 
  Base class for Data which do not need to be resolved.
  Internally this means they have getVectorRO and getVectorRW methods.
*/
class DataReady : public DataAbstract
{
typedef DataAbstract parent; 
public:
   DataReady(const FunctionSpace& what, const ShapeType& shape, bool isDataEmpty=false);
   ~DataReady(){};



  /**
     \brief
     Return the sample data for the given sample number.
  */
  ESCRIPT_DLL_API
  double*
  getSampleDataRW(DataTypes::RealVectorType::size_type sampleNo);

  ESCRIPT_DLL_API
  const double*
  getSampleDataRO(DataTypes::RealVectorType::size_type sampleNo) const;

  /**
	\brief Provide access to underlying storage. Internal use only!
  */

  ESCRIPT_DLL_API
  virtual DataTypes::RealVectorType&
  getVectorRW()=0;


  ESCRIPT_DLL_API
  virtual const DataTypes::RealVectorType&
  getVectorRO() const=0;

  ESCRIPT_DLL_API
  virtual DataTypes::CplxVectorType&
  getVectorRWC()=0;


  ESCRIPT_DLL_API
  virtual const DataTypes::CplxVectorType&
  getVectorROC() const=0;
  
  /**
     \brief These versions use the type system rather than method name to determine return type
  */
  ESCRIPT_DLL_API
  virtual DataTypes::RealVectorType&
  getTypedVectorRW(DataTypes::real_t dummy)=0;  
  
  ESCRIPT_DLL_API
  virtual const DataTypes::RealVectorType&
  getTypedVectorRO(DataTypes::real_t dummy) const=0;

  ESCRIPT_DLL_API
  virtual DataTypes::CplxVectorType&
  getTypedVectorRW(DataTypes::cplx_t dummy)=0;
  
  ESCRIPT_DLL_API
  virtual const DataTypes::CplxVectorType&
  getTypedVectorRO(DataTypes::cplx_t dummy) const=0;  
  

  
  
  
  /**
  \brief return true if data contains NaN.
  \warning This is dependent on the ability to reliably detect NaNs on your compiler.
   See the nancheck function in LocalOps for details.
  */
  ESCRIPT_DLL_API
  virtual bool
  hasNaN() const=0;

  /**
  \brief replaces all NaN values with value 
  */
  ESCRIPT_DLL_API
  virtual void
  replaceNaN(DataTypes::real_t value) = 0;
  
  /**
  \brief replaces all NaN values with value 
  */
  ESCRIPT_DLL_API
  virtual void
  replaceNaN(DataTypes::cplx_t value) = 0;  
  
  /**
     \brief
     Copy the specified region from the given object.

     \param value - Input - Data to copy from
     \param region - Input - Region to copy.
  */
  ESCRIPT_DLL_API
  virtual
  void
  setSlice(const DataAbstract* value,
           const DataTypes::RegionType& region) = 0;


 /**
     \brief get a reference to the beginning of a data point
 */
  ESCRIPT_DLL_API
  DataTypes::RealVectorType::const_reference
  getDataAtOffsetRO(DataTypes::RealVectorType::size_type i) const;


  ESCRIPT_DLL_API
  DataTypes::RealVectorType::reference
  getDataAtOffsetRW(DataTypes::RealVectorType::size_type i);
  
  ESCRIPT_DLL_API
  DataTypes::CplxVectorType::const_reference
  getDataAtOffsetROC(DataTypes::CplxVectorType::size_type i) const;


  ESCRIPT_DLL_API
  DataTypes::CplxVectorType::reference
  getDataAtOffsetRWC(DataTypes::CplxVectorType::size_type i);  
  
  
  
  
  
  ESCRIPT_DLL_API
  DataReady_ptr 
  resolve();

};


inline
DataTypes::real_t*
DataReady::getSampleDataRW(DataTypes::RealVectorType::size_type sampleNo)
{
  return &(getVectorRW()[getPointOffset(sampleNo,0)]);		// exclusive write checks will be done in getVectorRW()
}

inline const DataTypes::real_t*
DataReady::getSampleDataRO(DataTypes::RealVectorType::size_type sampleNo) const
{
  return &(getVectorRO()[getPointOffset(sampleNo,0)]);		
}


inline
DataTypes::RealVectorType::const_reference
DataReady::getDataAtOffsetRO(DataTypes::RealVectorType::size_type i) const
{
   return getVectorRO()[i];
}

inline
DataTypes::RealVectorType::reference
DataReady::getDataAtOffsetRW(DataTypes::RealVectorType::size_type i)	// exclusive write checks will be done in getVectorRW()
{
   return getVectorRW()[i];
}


inline
DataTypes::CplxVectorType::const_reference
DataReady::getDataAtOffsetROC(DataTypes::CplxVectorType::size_type i) const
{
   return getVectorROC()[i];
}

inline
DataTypes::CplxVectorType::reference
DataReady::getDataAtOffsetRWC(DataTypes::CplxVectorType::size_type i)	// exclusive write checks will be done in getVectorRW()
{
   return getVectorRWC()[i];
}



}

#endif
