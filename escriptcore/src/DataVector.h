
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

#ifndef __ESCRIPT_DATAVECTOR_H__
#define __ESCRIPT_DATAVECTOR_H__

#include "system_dep.h"
#include "DataTypes.h"
#include "Assert.h"
#include "DataVectorAlt.h"
#include "DataVectorTaipan.h"

#ifdef ESYS_HAVE_BOOST_NUMPY
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#endif

// ensure that nobody else tries to instantiate the complex version
extern template class escript::DataTypes::DataVectorAlt<escript::DataTypes::cplx_t>;


namespace escript {

// Functions in DataTypes:: which manipulate DataVectors
namespace DataTypes
{
  
  // This is the main version we had
  //typedef DataVectorTaipan DataVector;
  typedef escript::DataTypes::DataVectorAlt<real_t> RealVectorType;//!< Vector to store underlying data.
  typedef escript::DataTypes::DataVectorAlt<cplx_t> CplxVectorType;

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
   pointToStream(std::ostream& os, const RealVectorType::ElementType* data,const ShapeType& shape, int offset, bool needsep=true, const std::string& sep=",");

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
   pointToStream(std::ostream& os, const CplxVectorType::ElementType* data,const ShapeType& shape, int offset, bool needsep=true, const std::string& sep=",");


#ifdef ESYS_HAVE_BOOST_NUMPY

   void
   pointToNumpyArrayOld(boost::python::numpy::ndarray& dataArray, const RealVectorType::ElementType* data, const ShapeType& shape, long offset, long numsamples, long dpps, long numdata);


   void
   pointToNumpyArrayOld(boost::python::numpy::ndarray& dataArray, const CplxVectorType::ElementType* data, const ShapeType& shape, long offset, long numsamples, long dpps, long numdata);


   /**
      \brief Display a single value (with the specified shape) from the data.

     \param dataArray - numpy ndarray that data is being written to
     \param data - vector containing the datapoint
     \param shape - shape of the datapoint
     \param offset - start of the datapoint within data
     \param numsamples - index looping over the number of samples
     \param dpps - index looping over the data points in the sample
     \param numdata - total amount of data being copied
   */
   void
   pointToNumpyArray(boost::python::numpy::ndarray& dataArray, const RealVectorType::ElementType* data, const ShapeType& shape, int offset, int d, int index);


   void
   pointToNumpyArray(boost::python::numpy::ndarray& dataArray, const CplxVectorType::ElementType* data, const ShapeType& shape, int offset, int d, int index);
#endif

   /**
      \brief Display a single value (with the specified shape) from the data.

     \param data - vector containing the datapoint
     \param shape - shape of the datapoint
     \param offset - start of the datapoint within data
     \param prefix - string to prepend to the output
   */
   std::string
   pointToString(const RealVectorType& data,const ShapeType& shape, int offset, const std::string& prefix);


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
   void copyPoint(RealVectorType& dest, vec_size_type doffset, vec_size_type nvals, const RealVectorType& src, vec_size_type soffset);

   /**
      \brief  Copy a point from one vector to another. Note: This version does not check to see if shapes are the same.

   \param dest - vector to copy to
   \param doffset - beginning of the target datapoint in dest
   \param nvals - the number of values comprising the datapoint
   \param src - vector to copy from
   \param soffset - beginning of the datapoint in src
   */
   void copyPoint(CplxVectorType& dest, vec_size_type doffset, vec_size_type nvals, const CplxVectorType& src, vec_size_type soffset);

   /**
    * \brief copy data from a real vector to a complex vector
    * The complex vector will be resized as needed and any previous
    * values will be replaced.
   */
   void fillComplexFromReal(const RealVectorType& r, CplxVectorType& c);

  /**
     \brief
     Copy a data slice specified by the given region and offset from the
     "other" vector into the "left" vector at the given offset.

     \param left - vector to copy into
     \param leftShape - shape of datapoints for the left vector
     \param leftOffset - location within left to start copying to
     \param other - vector to copy from
     \param otherShape - shape of datapoints for the other vector
     \param otherOffset - location within other vector to start copying from
     \param region - Input -
                      Region in other vector to copy data from.
  */
   template <class VEC>
   // ESCRIPT_DLL_API
   void
   copySlice(VEC& left,
             const ShapeType& leftShape,
             typename VEC::size_type leftOffset,
             const VEC& other,
             const ShapeType& otherShape,
             typename VEC::size_type otherOffset,
             const RegionLoopRangeType& region)
   {
      //
      // Make sure vectors are not empty

      ESYS_ASSERT(!left.size()==!1, "left data is empty."); //Note: !1=0, but clang returns an error if the rhs is 0 here
      ESYS_ASSERT(!other.size()==!1, "other data is empty.");

      //
      // Check the vector to be sliced from is compatible with the region to be sliced,
      // and that the region to be sliced is compatible with this vector:
      ESYS_ASSERT(checkOffset(leftOffset,left.size(),noValues(leftShape)),
                 "offset incompatible with this vector.");
      ESYS_ASSERT(otherOffset+noValues(leftShape)<=other.size(),
                 "offset incompatible with other vector.");

      ESYS_ASSERT(getRank(otherShape)==region.size(),
                 "slice not same rank as vector to be sliced from.");

      ESYS_ASSERT(noValues(leftShape)==noValues(getResultSliceShape(region)),
                 "slice shape not compatible shape for this vector.");

      //
      // copy the values in the specified region of the other vector into this vector

      // the following loops cannot be parallelised due to the numCopy counter
      int numCopy=0;

      switch (region.size()) {
      case 0:
         /* this case should never be encountered,
         as python will never pass us an empty region.
         here for completeness only, allows slicing of a scalar */
//          (*m_data)[leftOffset+numCopy]=(*other.m_data)[otherOffset+other.relIndex()];

         left[leftOffset+numCopy]=other[otherOffset];
         numCopy++;
         break;
      case 1:
         for (int i=region[0].first;i<region[0].second;i++) {
            left[leftOffset+numCopy]=other[otherOffset+getRelIndex(otherShape,i)];
            numCopy++;
         }
         break;
      case 2:
         for (int j=region[1].first;j<region[1].second;j++) {
            for (int i=region[0].first;i<region[0].second;i++) {
/*               (*m_data)[leftOffset+numCopy]=(*other.m_data)[otherOffset+other.relIndex(i,j)];*/
               left[leftOffset+numCopy]=other[otherOffset+getRelIndex(otherShape,i,j)];
               numCopy++;
            }
         }
         break;
      case 3:
         for (int k=region[2].first;k<region[2].second;k++) {
            for (int j=region[1].first;j<region[1].second;j++) {
               for (int i=region[0].first;i<region[0].second;i++) {
//                  (*m_data)[leftOffset+numCopy]=(*other.m_data)[otherOffset+other.relIndex(i,j,k)];
                  left[leftOffset+numCopy]=other[otherOffset+getRelIndex(otherShape,i,j,k)];
                  numCopy++;
               }
            }
         }
         break;
      case 4:
         for (int l=region[3].first;l<region[3].second;l++) {
            for (int k=region[2].first;k<region[2].second;k++) {
               for (int j=region[1].first;j<region[1].second;j++) {
                  for (int i=region[0].first;i<region[0].second;i++) {
/*                     (*m_data)[leftOffset+numCopy]=(*other.m_data)[otherOffset+other.relIndex(i,j,k,l)];*/
                     left[leftOffset+numCopy]=other[otherOffset+getRelIndex(otherShape,i,j,k,l)];
                     numCopy++;
                  }
               }
            }
         }
         break;
      default:
         std::stringstream mess;
         mess << "Error - (copySlice) Invalid slice region rank: " << region.size();
         throw DataException(mess.str());
      }
   }

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
   template<typename VEC>
   // ESCRIPT_DLL_API
   void
   copySliceFrom(VEC& left,
                 const ShapeType& leftShape,
                 typename VEC::size_type leftOffset,
                 const VEC& other,
                 const ShapeType& otherShape,
                 typename VEC::size_type otherOffset,
                 const RegionLoopRangeType& region)
   {
      //
      // Make sure vectors are not empty

      ESYS_ASSERT(left.size()!=0, "this vector is empty.");
      ESYS_ASSERT(other.size()!=0, "other vector is empty.");

      //
      // Check this vector is compatible with the region to be sliced,
      // and that the region to be sliced is compatible with the other vector:

      ESYS_ASSERT(checkOffset(otherOffset,other.size(),noValues(otherShape)),
                 "offset incompatible with other vector.");
      ESYS_ASSERT(leftOffset+noValues(otherShape)<=left.size(),
                 "offset incompatible with this vector.");

      ESYS_ASSERT(getRank(leftShape)==region.size(),
                 "slice not same rank as this vector.");

      ESYS_ASSERT(getRank(otherShape)==0 || noValues(otherShape)==noValues(getResultSliceShape(region)),
                 "slice shape not compatible shape for other vector.");

      //
      // copy the values in the other vector into the specified region of this vector

      // allow for case where other vector is a scalar
      if (getRank(otherShape)==0) {

         // the following loops cannot be parallelised due to the numCopy counter
         int numCopy=0;

         switch (region.size()) {
         case 0:
            /* this case should never be encountered,
            as python will never pass us an empty region.
            here for completeness only, allows slicing of a scalar */
            //(*m_data)[leftOffset+relIndex()]=(*other.m_data)[otherOffset];
            left[leftOffset]=other[otherOffset];
            numCopy++;
            break;
         case 1:
            for (int i=region[0].first;i<region[0].second;i++) {
               left[leftOffset+getRelIndex(leftShape,i)]=other[otherOffset];
               numCopy++;
            }
            break;
         case 2:
            for (int j=region[1].first;j<region[1].second;j++) {
               for (int i=region[0].first;i<region[0].second;i++) {
                  left[leftOffset+getRelIndex(leftShape,i,j)]=other[otherOffset];
                  numCopy++;
               }
            }
            break;
         case 3:
            for (int k=region[2].first;k<region[2].second;k++) {
               for (int j=region[1].first;j<region[1].second;j++) {
                  for (int i=region[0].first;i<region[0].second;i++) {
                     left[leftOffset+getRelIndex(leftShape,i,j,k)]=other[otherOffset];
                     numCopy++;
                  }
               }
            }
            break;
         case 4:
            for (int l=region[3].first;l<region[3].second;l++) {
               for (int k=region[2].first;k<region[2].second;k++) {
                  for (int j=region[1].first;j<region[1].second;j++) {
                     for (int i=region[0].first;i<region[0].second;i++) {
                        left[leftOffset+getRelIndex(leftShape,i,j,k,l)]=other[otherOffset];
                        numCopy++;
                     }
                  }
               }
            }
            break;
         default:
            numCopy = 0;
            std::stringstream mess;
            mess << "Error - (copySliceFrom) Invalid slice region rank: " << region.size();
            throw DataException(mess.str());
         }

      } else {

         // the following loops cannot be parallelised due to the numCopy counter
         int numCopy=0;

         switch (region.size()) {
         case 0:
            /* this case should never be encountered,
            as python will never pass us an empty region.
            here for completeness only, allows slicing of a scalar */
            //(*m_data)[leftOffset+relIndex()]=(*other.m_data)[otherOffset+numCopy];
            left[leftOffset]=other[otherOffset+numCopy];
            numCopy++;
            break;
         case 1:
            for (int i=region[0].first;i<region[0].second;i++) {
               left[leftOffset+getRelIndex(leftShape,i)]=other[otherOffset+numCopy];
               numCopy++;
            }
            break;
         case 2:
            for (int j=region[1].first;j<region[1].second;j++) {
               for (int i=region[0].first;i<region[0].second;i++) {
                  left[leftOffset+getRelIndex(leftShape,i,j)]=other[otherOffset+numCopy];
                  numCopy++;
               }
            }
            break;
         case 3:
            for (int k=region[2].first;k<region[2].second;k++) {
               for (int j=region[1].first;j<region[1].second;j++) {
                  for (int i=region[0].first;i<region[0].second;i++) {
                     left[leftOffset+getRelIndex(leftShape,i,j,k)]=other[otherOffset+numCopy];
                     numCopy++;
                  }
               }
            }
            break;
         case 4:
            for (int l=region[3].first;l<region[3].second;l++) {
               for (int k=region[2].first;k<region[2].second;k++) {
                  for (int j=region[1].first;j<region[1].second;j++) {
                     for (int i=region[0].first;i<region[0].second;i++) {
                        left[leftOffset+getRelIndex(leftShape,i,j,k,l)]=other[otherOffset+numCopy];
                        numCopy++;
                     }
                  }
               }
            }
            break;
         default:
            numCopy = 0;
            std::stringstream mess;
            mess << "Error - (copySliceFrom) Invalid slice region rank: " << region.size();
            throw DataException(mess.str());
         }

      }

   }
}

} // end of namespace

#endif // __ESCRIPT_DATAVECTOR_H__

