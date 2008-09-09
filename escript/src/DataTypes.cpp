/* $Id:$ */

/*******************************************************
*
*           Copyright 2003-2007 by ACceSS MNRF
*       Copyright 2007 by University of Queensland
*
*                http://esscc.uq.edu.au
*        Primary Business: Queensland, Australia
*  Licensed under the Open Software License version 3.0
*     http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/

#include <sstream>
#include <boost/python/extract.hpp>
#include <boost/python/tuple.hpp>
#include "DataException.h"
#include "DataTypes.h"

using namespace boost::python;

namespace escript
{
namespace DataTypes
{

   int
   noValues(const ShapeType& shape) 
   {
      ShapeType::const_iterator i;
      //
      // An empty shape vector means rank 0 which contains 1 value
      int noValues=1;
      for (i=shape.begin();i!=shape.end();i++) {
         noValues*=(*i);
      }
      return noValues;
   }

   int
   noValues(const RegionLoopRangeType& region) 
   {
      //
      // An empty region vector means rank 0 which contains 1 value
      int noValues=1;
      unsigned int i;
      for (i=0;i<region.size();i++) {
         noValues*=region[i].second-region[i].first;
      }
      return noValues;
   }

   std::string
   shapeToString(const DataTypes::ShapeType& shape)
   {
      std::stringstream temp;
      temp << "(";
      unsigned int i;
      for (i=0;i<shape.size();i++) {
         temp << shape[i];
         if (i < shape.size()-1) {
            temp << ",";
         }
      }
      temp << ")";
      return temp.str();
   }


/*
  \brief
  Calculate the slice range from the given python key object
  Used by getSliceRegion - since it is not used anywhere else, I have elected not to declare it
  in the header.
  Returns the python slice object key as a pair of ints where the first 
  member is the start and the second member is the end. the presence of a possible
  step attribute with value different from one will throw an exception

  /param key - Input - key object specifying slice range.
*/
   std::pair<int,int>
   getSliceRange(const boost::python::object& key,
              const int shape)
   {
      /* default slice range is range of entire shape dimension */
      int s0=0, s1=shape;;
      extract<int> slice_int(key);
      if (slice_int.check()) {
         /* if the key is a single int set start=key and end=key */
         /* in this case, we want to return a rank-1 dimension object from
         this object, taken from a single index value for one of this
         object's dimensions */
         s0=slice_int();
         s1=s0;
      } else {
         /* if key is a pair extract begin and end values */
         extract<int> step(key.attr("step"));
         if (step.check() && step()!=1) {
            throw DataException("Error - Data does not support increments in slicing ");
         } else {
            extract<int> start(key.attr("start"));
            if (start.check()) {
               s0=start();
            }
            extract<int> stop(key.attr("stop"));
            if (stop.check()) {
               s1=stop();
            }
         }
      }
      if (s0 < 0) 
         throw DataException("Error - slice index out of range.");
      if (s0 == s1 && s1 >= shape)
         throw DataException("Error - slice index out of range.");
      if (s0 != s1 &&  s1>shape)
         throw DataException("Error - slice index out of range.");
      if (s0 > s1) 
         throw DataException("Error - lower index must less or equal upper index.");
      return std::pair<int,int>(s0,s1);
   }



   DataTypes::RegionType
   getSliceRegion(const DataTypes::ShapeType& shape, const boost::python::object& key)
   {
      int slice_rank, i;
      int this_rank=shape.size();
      RegionType out(this_rank);
      /* allow for case where key is singular eg: [1], this implies we
      want to generate a rank-1 dimension object, as opposed to eg: [1,2]
      which implies we want to take a rank dimensional object with one
      dimension of size 1 */
      extract<tuple> key_tuple(key);
      if (key_tuple.check()) {
         slice_rank=extract<int> (key.attr("__len__")());
         /* ensure slice is correctly dimensioned */
         if (slice_rank>this_rank) {
            throw DataException("Error - rank of slices does not match rank of slicee");
         } else {
            /* calculate values for slice range */
            for (i=0;i<slice_rank;i++) {
               out[i]=getSliceRange(key[i],shape[i]);
            }
         }
      } else {
         slice_rank=1;
         if (slice_rank>this_rank) {
            throw DataException("Error - rank of slices does not match rank of slicee");
         } else {
            out[0]=getSliceRange(key,shape[0]);
         }
      }
      for (i=slice_rank;i<this_rank;i++) {
         out[i]=std::pair<int,int>(0,shape[i]);
      }
      return out;
   }

   DataTypes::ShapeType
   getResultSliceShape(const RegionType& region)
   {
      int dimSize;
      ShapeType result;
      RegionType::const_iterator i;
      for (i=region.begin();i!=region.end();i++) {
         dimSize=((i->second)-(i->first));
         if (dimSize!=0) {
            result.push_back(dimSize);
         }
      }
      return result;
   }

   DataTypes::RegionLoopRangeType
   getSliceRegionLoopRange(const DataTypes::RegionType& region) 
   {
      DataTypes::RegionLoopRangeType region_loop_range(region.size());
      unsigned int i;
      for (i=0;i<region.size();i++) {
         if (region[i].first==region[i].second) {
            region_loop_range[i].first=region[i].first;
            region_loop_range[i].second=region[i].second+1;
         } else {
            region_loop_range[i].first=region[i].first;
            region_loop_range[i].second=region[i].second;
         }
      }
      return region_loop_range;
   }


   std::string 
   createShapeErrorMessage(const std::string& messagePrefix,
                                          const DataTypes::ShapeType& other,
					  const DataTypes::ShapeType& thisShape)
   {
      std::stringstream temp;
      temp << messagePrefix
           << " This shape: " << shapeToString(thisShape)
           << " Other shape: " << shapeToString(other);
      return temp.str();
   }


// Additional slice operations

   inline
   bool
   checkOffset(ValueType::size_type offset, int size, int noval)
   {
      return (size >= (offset+noval));
   }


   void
   copySlice(ValueType& left,
			    const ShapeType& leftShape,
			    ValueType::size_type thisOffset,
                            const ValueType& other,
			    const ShapeType& otherShape,
                            ValueType::size_type otherOffset,
                            const RegionLoopRangeType& region)
   {
      //
      // Make sure views are not empty

      EsysAssert(!left.size()==0,
                 "Error - left data is empty.");
      EsysAssert(!other.size()==0,
                 "Error - other data is empty.");

      //
      // Check the view to be sliced from is compatible with the region to be sliced,
      // and that the region to be sliced is compatible with this view:
      EsysAssert(checkOffset(thisOffset,left.size(),noValues(leftShape)),
                 "Error - offset incompatible with this view.");
      EsysAssert(otherOffset+noValues(leftShape)<=other.size(),
                 "Error - offset incompatible with other view.");

      EsysAssert(getRank(otherShape)==region.size(),
                 "Error - slice not same rank as view to be sliced from.");

      EsysAssert(noValues(leftShape)==noValues(getResultSliceShape(region)),
                 "Error - slice shape not compatible shape for this view.");

      //
      // copy the values in the specified region of the other view into this view

      // the following loops cannot be parallelised due to the numCopy counter
      int numCopy=0;

      switch (region.size()) {
      case 0:
         /* this case should never be encountered, 
         as python will never pass us an empty region.
         here for completeness only, allows slicing of a scalar */
//          (*m_data)[thisOffset+numCopy]=(*other.m_data)[otherOffset+other.relIndex()];

         left[thisOffset+numCopy]=other[otherOffset];
         numCopy++;
         break;
      case 1:
         for (int i=region[0].first;i<region[0].second;i++) {
            left[thisOffset+numCopy]=other[otherOffset+getRelIndex(otherShape,i)];
            numCopy++;
         }
         break;
      case 2:
         for (int j=region[1].first;j<region[1].second;j++) {
            for (int i=region[0].first;i<region[0].second;i++) {
/*               (*m_data)[thisOffset+numCopy]=(*other.m_data)[otherOffset+other.relIndex(i,j)];*/
               left[thisOffset+numCopy]=other[otherOffset+getRelIndex(otherShape,i,j)];
               numCopy++;
            }
         }
         break;
      case 3:
         for (int k=region[2].first;k<region[2].second;k++) {
            for (int j=region[1].first;j<region[1].second;j++) {
               for (int i=region[0].first;i<region[0].second;i++) {
//                  (*m_data)[thisOffset+numCopy]=(*other.m_data)[otherOffset+other.relIndex(i,j,k)];
                  left[thisOffset+numCopy]=other[otherOffset+getRelIndex(otherShape,i,j,k)];
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
/*                     (*m_data)[thisOffset+numCopy]=(*other.m_data)[otherOffset+other.relIndex(i,j,k,l)];*/
                     left[thisOffset+numCopy]=other[otherOffset+getRelIndex(otherShape,i,j,k,l)];
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


   void
   copySliceFrom(ValueType& left,
				const ShapeType& leftShape,
				ValueType::size_type thisOffset,
                                const ValueType& other,
				const ShapeType& otherShape,
                                ValueType::size_type otherOffset,
                                const RegionLoopRangeType& region)
   {
      //
      // Make sure views are not empty

      EsysAssert(left.size()!=0,
                 "Error - this view is empty.");
      EsysAssert(other.size()!=0,
                 "Error - other view is empty.");

      //
      // Check this view is compatible with the region to be sliced,
      // and that the region to be sliced is compatible with the other view:

      EsysAssert(checkOffset(otherOffset,other.size(),noValues(otherShape)),
                 "Error - offset incompatible with other view.");
      EsysAssert(thisOffset+noValues(otherShape)<=left.size(),
                 "Error - offset incompatible with this view.");

      EsysAssert(getRank(leftShape)==region.size(),
                 "Error - slice not same rank as this view.");

      EsysAssert(getRank(otherShape)==0 || noValues(otherShape)==noValues(getResultSliceShape(region)),
                 "Error - slice shape not compatible shape for other view.");

      //
      // copy the values in the other view into the specified region of this view

      // allow for case where other view is a scalar
      if (getRank(otherShape)==0) {

         // the following loops cannot be parallelised due to the numCopy counter
         int numCopy=0;

         switch (region.size()) {
         case 0:
            /* this case should never be encountered, 
            as python will never pass us an empty region.
            here for completeness only, allows slicing of a scalar */
            //(*m_data)[thisOffset+relIndex()]=(*other.m_data)[otherOffset];
	    left[thisOffset]=other[otherOffset];
            numCopy++;
            break;
         case 1:
            for (int i=region[0].first;i<region[0].second;i++) {
               left[thisOffset+getRelIndex(leftShape,i)]=other[otherOffset];
               numCopy++;
            }
            break;
         case 2:
            for (int j=region[1].first;j<region[1].second;j++) {
               for (int i=region[0].first;i<region[0].second;i++) {
                  left[thisOffset+getRelIndex(leftShape,i,j)]=other[otherOffset];
                  numCopy++;
               }
            }
            break;
         case 3:
            for (int k=region[2].first;k<region[2].second;k++) {
               for (int j=region[1].first;j<region[1].second;j++) {
                  for (int i=region[0].first;i<region[0].second;i++) {
                     left[thisOffset+getRelIndex(leftShape,i,j,k)]=other[otherOffset];
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
                        left[thisOffset+getRelIndex(leftShape,i,j,k,l)]=other[otherOffset];
                        numCopy++;
                     }
                  }
               }
            }
            break;
         default:
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
            //(*m_data)[thisOffset+relIndex()]=(*other.m_data)[otherOffset+numCopy];
	    left[thisOffset]=other[otherOffset+numCopy];
            numCopy++;
            break;
         case 1:
            for (int i=region[0].first;i<region[0].second;i++) {
               left[thisOffset+getRelIndex(leftShape,i)]=other[otherOffset+numCopy];
               numCopy++;
            }
            break;
         case 2:
            for (int j=region[1].first;j<region[1].second;j++) {
               for (int i=region[0].first;i<region[0].second;i++) {
                  left[thisOffset+getRelIndex(leftShape,i,j)]=other[otherOffset+numCopy];
                  numCopy++;
               }
            }
            break;
         case 3:
            for (int k=region[2].first;k<region[2].second;k++) {
               for (int j=region[1].first;j<region[1].second;j++) {
                  for (int i=region[0].first;i<region[0].second;i++) {
                     left[thisOffset+getRelIndex(leftShape,i,j,k)]=other[otherOffset+numCopy];
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
                        left[thisOffset+getRelIndex(leftShape,i,j,k,l)]=other[otherOffset+numCopy];
                        numCopy++;
                     }
                  }
               }
            }
            break;
         default:
            std::stringstream mess;
            mess << "Error - (copySliceFrom) Invalid slice region rank: " << region.size();
            throw DataException(mess.str());
         }

      }

   }

   std::string
   pointToString(const ValueType& data,const ShapeType& shape, int offset, const std::string& suffix)
   {
      using namespace std;
      EsysAssert(data.size()>0,"Error - Data object is empty.");
      stringstream temp;
      string finalSuffix=suffix;
      if (suffix.length() > 0) {
         finalSuffix+=" ";
      }
      switch (getRank(shape)) {
      case 0:
         temp << finalSuffix << data[0];
         break;
      case 1:
         for (int i=0;i<shape[0];i++) {
            temp << finalSuffix << "(" << i <<  ") " << data[i+offset];
            if (i!=(shape[0]-1)) {
               temp << endl;
            }
         }
         break;
      case 2:
         for (int i=0;i<shape[0];i++) {
            for (int j=0;j<shape[1];j++) {
               temp << finalSuffix << "(" << i << "," << j << ") " << data[offset+getRelIndex(shape,i,j)];
               if (!(i==(shape[0]-1) && j==(shape[1]-1))) {
                  temp << endl;
               }
            }
         }
         break;
      case 3:
         for (int i=0;i<shape[0];i++) {
            for (int j=0;j<shape[1];j++) {
               for (int k=0;k<shape[2];k++) {
                  temp << finalSuffix << "(" << i << "," << j << "," << k << ") " << data[offset+getRelIndex(shape,i,j,k)];
                  if (!(i==(shape[0]-1) && j==(shape[1]-1) && k==(shape[2]-1))) {
                     temp << endl;
                  }
               }
            }
         }
         break;
      case 4:
         for (int i=0;i<shape[0];i++) {
            for (int j=0;j<shape[1];j++) {
               for (int k=0;k<shape[2];k++) {
                  for (int l=0;l<shape[3];l++) {
                     temp << finalSuffix << "(" << i << "," << j << "," << k << "," << l << ") " << data[offset+getRelIndex(shape,i,j,k,l)];
                     if (!(i==(shape[0]-1) && j==(shape[1]-1) && k==(shape[2]-1) && l==(shape[3]-1))) {
                        temp << endl;
                     }
                  }
               }
            }
         }
         break;
      default:
         stringstream mess;
         mess << "Error - (toString) Invalid rank: " << getRank(shape);
         throw DataException(mess.str());
      }
      return temp.str();
   }


   void copyPoint(ValueType& dest, ValueType::size_type doffset, ValueType::size_type nvals, const ValueType& src, ValueType::size_type soffset)
   {
      EsysAssert((dest.size()>0&&src.size()>0&&checkOffset(doffset,dest.size(),nvals)),
                 "Error - Couldn't copy due to insufficient storage.");
//       EsysAssert((checkShape(other.getShape())),
//                  createShapeErrorMessage("Error - Couldn't copy due to shape mismatch.",other.getShape(),m_shape));
      if (checkOffset(doffset,dest.size(),nvals) && checkOffset(soffset,src.size(),nvals)) {
         memcpy(&dest[doffset],&src[soffset],sizeof(double)*nvals);
      } else {
         throw DataException("Error - invalid offset specified.");
      }



   } 

}	// end namespace DataTypes
}	// end namespace escript