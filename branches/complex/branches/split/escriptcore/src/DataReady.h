
/*****************************************************************************
*
* Copyright (c) 2003-2014 by University of Queensland
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
  getSampleDataRW(ValueType::size_type sampleNo);

  ESCRIPT_DLL_API
  const double*
  getSampleDataRO(ValueType::size_type sampleNo) const;

  /**
	\brief Provide access to underlying storage. Internal use only!
  */

  ESCRIPT_DLL_API
  virtual DataTypes::ValueType&
  getVectorRW()=0;


  ESCRIPT_DLL_API
  virtual const DataTypes::ValueType&
  getVectorRO() const=0;

  /**
  \brief return true if data contains NaN.
  \warning This is dependent on the ability to reliably detect NaNs on your compiler.
   See the nancheck function in LocalOps for details.
  */
  ESCRIPT_DLL_API
  virtual bool
  hasNaN() const=0;

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
  DataTypes::ValueType::const_reference
  getDataAtOffsetRO(DataTypes::ValueType::size_type i) const;


  ESCRIPT_DLL_API
  DataTypes::ValueType::reference
  getDataAtOffsetRW(DataTypes::ValueType::size_type i);

  ESCRIPT_DLL_API
  DataReady_ptr 
  resolve();

};


inline
DataAbstract::ValueType::value_type*
DataReady::getSampleDataRW(ValueType::size_type sampleNo)
{
  return &(getVectorRW()[getPointOffset(sampleNo,0)]);		// exclusive write checks will be done in getVectorRW()
}

inline const double*
DataReady::getSampleDataRO(ValueType::size_type sampleNo) const
{
  return &(getVectorRO()[getPointOffset(sampleNo,0)]);		
}


inline
DataTypes::ValueType::const_reference
DataReady::getDataAtOffsetRO(DataTypes::ValueType::size_type i) const
{
   return getVectorRO()[i];
}

inline
DataTypes::ValueType::reference
DataReady::getDataAtOffsetRW(DataTypes::ValueType::size_type i)	// exclusive write checks will be done in getVectorRW()
{
   return getVectorRW()[i];
}



}

#endif
