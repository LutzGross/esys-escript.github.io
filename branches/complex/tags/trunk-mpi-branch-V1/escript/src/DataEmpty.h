//$Id$
/* 
 ************************************************************
 *          Copyright 2006 by ACcESS MNRF                   *
 *                                                          *
 *              http://www.access.edu.au                    *
 *       Primary Business: Queensland, Australia            *
 *  Licensed under the Open Software License version 3.0    *
 *     http://www.opensource.org/licenses/osl-3.0.php       *
 *                                                          *
 ************************************************************
*/

#if !defined escript_DataEmpty_20040726_H
#define escript_DataEmpty_20040726_H
#include "system_dep.h"

#include "DataAbstract.h"

namespace escript {

/**
   \brief
   Implements the DataAbstract interface for an empty Data object.

   Description:
   Implements the DataAbstract interface for an empty Data object.
*/

class DataEmpty : public DataAbstract {

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
     \brief
     Return the offset for the given sample.
     NB: This will throw an exception as obviously an empty Data object contains no
     samples. An implementation is required by parent DataAbstract class.
     \param sampleNo - Input - Sample number.
     \param dataPointNo - Input - data-point number.
   */
  ESCRIPT_DLL_API
  virtual
  DataArrayView::ValueType::size_type
  getPointOffset(int sampleNo,
                 int dataPointNo) const;

  /**
     \brief
     Return a view into the data for the data point specified.
     NB: This will throw an exception as obviously an empty Data object contains no
     data points. An implementation is required by parent DataAbstract class.
     \param sampleNo - Input - Sample number.
     \param dataPointNo - Input - data-point number.
  */
  ESCRIPT_DLL_API
  virtual
  DataArrayView
  getDataPoint(int sampleNo,
               int dataPointNo);

  /**
     \brief
     Return the number of doubles stored for the Data object.
     As this is an empty Data object, this method will always return 0.
  */
  ESCRIPT_DLL_API
  virtual
  ValueType::size_type
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
  getSlice(const DataArrayView::RegionType& region) const;

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
           const DataArrayView::RegionType& region);

 protected:

 private:

  /**
     \brief
     Throw a standard exception. This function is called if an attempt
     is made to use functions of DataEmpty that are not valid.
  */
  void
  throwStandardException(const std::string& functionName) const;

};

} // end of namespace

#endif
