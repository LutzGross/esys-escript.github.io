
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


#if !defined escript_DataVector_H
#define escript_DataVector_H
#include "system_dep.h"

#include "esysUtils/EsysAssert.h"

#include <vector>
#include <iostream>
#include <fstream>

#include "DataTypes.h"

#include "DataVectorTaipan.h"
#include "DataVectorAlt.h"

// ensure that nobody else tries to instantiate the complex version
extern template class escript::DataTypes::DataVectorAlt<escript::DataTypes::cplx_t>;


namespace escript {

// Functions in DataTypes:: which manipulate DataVectors
namespace DataTypes
{
  
  // This is the main version we had
  //typedef DataVectorTaipan DataVector;
  typedef escript::DataTypes::DataVectorAlt<real_t> FloatVectorType;//!< Vector to store underlying data.
  typedef escript::DataTypes::DataVectorAlt<cplx_t> CplxVectorType;

  /**
     \brief
     Copy a data slice specified by the given region and offset from the
     "other" view into the "left" view at the given offset.
     
     \param left - vector to copy into
     \param leftShape - shape of datapoints for the left vector
     \param leftOffset - location within left to start copying to
     \param other - vector to copy from
     \param otherShape - shape of datapoints for the other vector
     \param otherOffset - location within other vector to start copying from
     \param region - Input -
                      Region in other view to copy data from.
  */
   ESCRIPT_DLL_API
   void
   copySlice(FloatVectorType& left,
			    const ShapeType& leftShape,
			    vec_size_type leftOffset,
                            const FloatVectorType& other,
			    const ShapeType& otherShape,
                            vec_size_type otherOffset,
                            const RegionLoopRangeType& region);

  /**
     \brief
     Copy data into a slice specified by the given region and offset in
     the left vector from the other vector at the given offset.

     \param left - vector to copy into
     \param leftShape - shape of datapoints for the left vector
     \param leftOffset - location within left to start copying to
     \param other - vector to copy from
     \param otherShape - shape of datapoints for the other vector
     \param otherOffset - location within other vector to start copying from
     \param region - Input -
                      Region in the left vector to copy data to.
  */
   ESCRIPT_DLL_API
   void
   copySliceFrom(FloatVectorType& left,
				const ShapeType& leftShape,
				vec_size_type leftOffset,
                                const FloatVectorType& other,
				const ShapeType& otherShape,
                                vec_size_type otherOffset,
                                const RegionLoopRangeType& region);


   /**
      \brief Display a single value (with the specified shape) from the data.

     Despite its similar name this function behaves differently to pointToString.
     There are no prefixes or (i,j,k) identifiers on each field. each datapoint is printed without
     new lines.
     It also works with double* rather than vectors so be careful what you pass it.

     \param os - stream to write to
     \param data - vector containing the datapoint
     \param shape - shape of the datapoint
     \param offset - start of the datapoint within data
     \param needsep - Does this output need to start with a separator
     \param sep - separator string to print between components
   */
   void
   pointToStream(std::ostream& os, const FloatVectorType::ElementType* data,const ShapeType& shape, int offset, bool needsep=true, const std::string& sep=",");

   /**
      \brief Display a single value (with the specified shape) from the data.

     \param data - vector containing the datapoint
     \param shape - shape of the datapoint
     \param offset - start of the datapoint within data
     \param prefix - string to prepend to the output
   */
   std::string
   pointToString(const FloatVectorType& data,const ShapeType& shape, int offset, const std::string& prefix);


   std::string
   pointToString(const CplxVectorType& data,const ShapeType& shape, int offset, const std::string& prefix);

   /**
      \brief  Copy a point from one vector to another. Note: This version does not check to see if shapes are the same.

   \param dest - vector to copy to
   \param doffset - beginning of the target datapoint in dest
   \param nvals - the number of values comprising the datapoint
   \param src - vector to copy from
   \param soffset - beginning of the datapoint in src
   */
   void copyPoint(FloatVectorType& dest, vec_size_type doffset, vec_size_type nvals, const FloatVectorType& src, vec_size_type soffset);  
  
  
  
  
}

 
} // end of namespace

#endif
