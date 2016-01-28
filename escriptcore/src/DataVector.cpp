
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

#define ESNEEDPYTHON
#include "esysUtils/first.h"


#include "DataVector.h"

#include "Taipan.h"
#include "DataException.h"
#include <boost/python/extract.hpp>
#include "DataTypes.h"
#include "WrappedArray.h"

#include <cassert>

using namespace std;
using namespace escript;
using namespace boost::python;
using namespace DataTypes;

namespace escript {

// Additional slice operations

   void
   DataTypes::copySlice(DataTypes::FloatVectorType& left,
			    const DataTypes::ShapeType& leftShape,
			    DataTypes::vec_size_type thisOffset,
                            const DataTypes::FloatVectorType& other,
			    const DataTypes::ShapeType& otherShape,
                            DataTypes::vec_size_type otherOffset,
                            const DataTypes::RegionLoopRangeType& region)
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
   DataTypes::copySliceFrom(DataTypes::FloatVectorType& left,
				const ShapeType& leftShape,
				DataTypes::vec_size_type thisOffset,
                                const FloatVectorType& other,
				const DataTypes::ShapeType& otherShape,
                                DataTypes::vec_size_type otherOffset,
                                const DataTypes::RegionLoopRangeType& region)
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


   void
   DataTypes::pointToStream(std::ostream& os, const FloatVectorType::ElementType* data,const ShapeType& shape, int offset, bool needsep, const std::string& sep)
   {
      using namespace std;
      EsysAssert(data!=0, "Error - data is null");
//      EsysAssert(data.size()>0,"Error - Data object is empty.");
      switch (getRank(shape)) {
      case 0:
	 if (needsep)
	 {
		os << sep;
	 }
	 else
	 {
		needsep=true;
	 }
         os << data[offset];
         break;
      case 1:
         for (int i=0;i<shape[0];i++) {
	    if (needsep)
	    {
		os << sep;
	    }
	    else
	    {
		needsep=true;
	    }
	    os << data[i+offset];
         }
         break;
      case 2:
         for (int i=0;i<shape[0];i++) {
            for (int j=0;j<shape[1];j++) {
		if (needsep)
		{
			os << sep;
		}
		else
		{
			needsep=true;
		}
                os << data[offset+getRelIndex(shape,i,j)];
            }
         }
         break;
      case 3:
         for (int i=0;i<shape[0];i++) {
            for (int j=0;j<shape[1];j++) {
               for (int k=0;k<shape[2];k++) {
		   if (needsep)
		   {
			os << sep;
		   }
		   else
		   {
			needsep=true;
		   }
                   os << data[offset+getRelIndex(shape,i,j,k)];
               }
            }
         }
         break;
      case 4:
         for (int i=0;i<shape[0];i++) {
            for (int j=0;j<shape[1];j++) {
               for (int k=0;k<shape[2];k++) {
                  for (int l=0;l<shape[3];l++) {
			if (needsep)
			{
				os << sep;
			}
			else
			{
				needsep=true;
			}
			os << data[offset+getRelIndex(shape,i,j,k,l)];
                  }
               }
            }
         }
         break;
      default:
         stringstream mess;
         mess << "Error - (pointToStream) Invalid rank: " << getRank(shape);
         throw DataException(mess.str());
      }
   }


   std::string
   DataTypes::pointToString(const CplxVectorType& data,const ShapeType& shape, int offset, const std::string& prefix)
   {
      using namespace std;
      EsysAssert(data.size()>0,"Error - Data object is empty.");
      stringstream temp;
      string finalPrefix=prefix;
      if (prefix.length() > 0) {
         finalPrefix+=" ";
      }
      switch (getRank(shape)) {
      case 0:
         temp << finalPrefix << data[offset];
         break;
      case 1:
         for (int i=0;i<shape[0];i++) {
            temp << finalPrefix << "(" << i <<  ") " << data[i+offset];
            if (i!=(shape[0]-1)) {
               temp << endl;
            }
         }
         break;
      case 2:
         for (int i=0;i<shape[0];i++) {
            for (int j=0;j<shape[1];j++) {
               temp << finalPrefix << "(" << i << "," << j << ") " << data[offset+getRelIndex(shape,i,j)];
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
                  temp << finalPrefix << "(" << i << "," << j << "," << k << ") " << data[offset+getRelIndex(shape,i,j,k)];
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
                     temp << finalPrefix << "(" << i << "," << j << "," << k << "," << l << ") " << data[offset+getRelIndex(shape,i,j,k,l)];
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

   std::string
   DataTypes::pointToString(const FloatVectorType& data,const ShapeType& shape, int offset, const std::string& prefix)
   {
      using namespace std;
      EsysAssert(data.size()>0,"Error - Data object is empty.");
      stringstream temp;
      string finalPrefix=prefix;
      if (prefix.length() > 0) {
         finalPrefix+=" ";
      }
      switch (getRank(shape)) {
      case 0:
         temp << finalPrefix << data[offset];
         break;
      case 1:
         for (int i=0;i<shape[0];i++) {
            temp << finalPrefix << "(" << i <<  ") " << data[i+offset];
            if (i!=(shape[0]-1)) {
               temp << endl;
            }
         }
         break;
      case 2:
         for (int i=0;i<shape[0];i++) {
            for (int j=0;j<shape[1];j++) {
               temp << finalPrefix << "(" << i << "," << j << ") " << data[offset+getRelIndex(shape,i,j)];
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
                  temp << finalPrefix << "(" << i << "," << j << "," << k << ") " << data[offset+getRelIndex(shape,i,j,k)];
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
                     temp << finalPrefix << "(" << i << "," << j << "," << k << "," << l << ") " << data[offset+getRelIndex(shape,i,j,k,l)];
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


   void DataTypes::copyPoint(FloatVectorType& dest, FloatVectorType::size_type doffset, FloatVectorType::size_type nvals, const FloatVectorType& src, FloatVectorType::size_type soffset)
   {
      EsysAssert((dest.size()>0&&src.size()>0&&checkOffset(doffset,dest.size(),nvals)),
                 "Error - Couldn't copy due to insufficient storage.");
      if (checkOffset(doffset,dest.size(),nvals) && checkOffset(soffset,src.size(),nvals)) {
         memcpy(&dest[doffset],&src[soffset],sizeof(double)*nvals);
      } else {
         throw DataException("Error - invalid offset specified.");
      }
   } 
   
   /**
    * \brief copy data from a real vector to a complex vector
    * The complex vector will be resized as needed and any previous
    * values will be replaced.
   */
   void DataTypes::fillComplexFromReal(FloatVectorType& r, CplxVectorType& c)
   {
       if (c.size()!=r.size())
       {
	   c.resize(r.size(), 0, 1);
       }
       size_t limit=r.size();
       #pragma omp parallel for schedule(static)
       for (size_t i=0;i<limit;++i)
       {
	   c[i]=r[i];
       }
   }

} // end of namespace
