
/*****************************************************************************
*
* Copyright (c) 2003-2018 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/

#include "DataVector.h"

#include "DataException.h"
#include "DataTypes.h"
#include "Taipan.h"
#include "WrappedArray.h"

#include <boost/python/extract.hpp>

using namespace std;
using namespace escript;
using namespace boost::python;
using namespace DataTypes;

namespace {

inline
void cplxout(std::ostream& os, const DataTypes::cplx_t& c)
{
    os << c.real();
    if (c.imag()>=0)
    {
        os << '+';
    }
    os << c.imag() << 'j';
}

}


namespace escript {

   void
   DataTypes::pointToStream(std::ostream& os, const CplxVectorType::ElementType* data,const ShapeType& shape, int offset, bool needsep, const std::string& sep)
   {
      using namespace std;
      ESYS_ASSERT(data!=0, "Error - data is null");
//      ESYS_ASSERT(data.size()>0,"Error - Data object is empty.");
      switch (getRank(shape)) {
      case 0:
         if (needsep){
            os << sep;
         } else {
            needsep = true;
         }
         cplxout(os,data[offset]);
         break;
      case 1:
         for (int i=0;i<shape[0];i++) {
            if (needsep) {
               os << sep;
            } else {
               needsep=true;
            }
            cplxout(os,data[i+offset]);
         }
         break;
      case 2:
         for (int i=0;i<shape[0];i++) {
            for (int j=0;j<shape[1];j++) {
               if (needsep){
                  os << sep;
               } else {
                  needsep=true;
               }
               cplxout(os,data[offset+getRelIndex(shape,i,j)]);
            }
         }
         break;
      case 3:
         for (int i=0;i<shape[0];i++) {
            for (int j=0;j<shape[1];j++) {
               for (int k=0;k<shape[2];k++) {
                  if (needsep) {
                     os << sep;
                  } else { 
                     needsep=true;
                  }
                  cplxout(os,data[offset+getRelIndex(shape,i,j,k)]);
               }
            }
         }
         break;
      case 4:
         for (int i=0;i<shape[0];i++) {
            for (int j=0;j<shape[1];j++) {
               for (int k=0;k<shape[2];k++) {
                  for (int l=0;l<shape[3];l++) {
                     if (needsep) {
                        os << sep;
                     } else {
                        needsep=true;
                     }
                     cplxout(os,data[offset+getRelIndex(shape,i,j,k,l)]);
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


   void
   DataTypes::pointToStream(std::ostream& os, const RealVectorType::ElementType* data,const ShapeType& shape, int offset, bool needsep, const std::string& sep)
   {
      using namespace std;
      ESYS_ASSERT(data!=0, "Error - data is null");
//      ESYS_ASSERT(data.size()>0,"Error - Data object is empty.");
      switch (getRank(shape)) {
      case 0:
         if (needsep) {
            os << sep;
         } else {
            needsep=true;
         }
         os << data[offset];
         break;
      case 1:
         for (int i=0;i<shape[0];i++) {
            if (needsep) {
               os << sep;
            } else {
               needsep=true;
            }
            os << data[i+offset];
         }
         break;
      case 2:
         for (int i=0;i<shape[0];i++) {
            for (int j=0;j<shape[1];j++) {
               if (needsep) {
                  os << sep;
               } else {
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
                  if (needsep) {
                     os << sep;
                  } else {
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
                     if (needsep) {
                        os << sep;
                     } else {
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
      ESYS_ASSERT(data.size()>0,"Error - Data object is empty.");
      stringstream temp;
      string finalPrefix=prefix;
      if (prefix.length() > 0) {
         finalPrefix+=" ";
      }
      switch (getRank(shape)) {
      case 0:
         temp << finalPrefix;
         cplxout(temp,data[offset]);
         break;
      case 1:
         for (int i=0;i<shape[0];i++) {
            temp << finalPrefix << "(" << i <<  ") ";
            cplxout(temp,data[i+offset]);
            if (i!=(shape[0]-1)) {
               temp << endl;
            }
         }
         break;
      case 2:
         for (int i=0;i<shape[0];i++) {
            for (int j=0;j<shape[1];j++) {
               temp << finalPrefix << "(" << i << "," << j << ") ";
               cplxout(temp,data[offset+getRelIndex(shape,i,j)]);
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
                  temp << finalPrefix << "(" << i << "," << j << "," << k << ") ";
                  cplxout(temp,data[offset+getRelIndex(shape,i,j,k)]);
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
                     temp << finalPrefix << "(" << i << "," << j << "," << k << "," << l << ") ";
                     cplxout(temp,data[offset+getRelIndex(shape,i,j,k,l)]);
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
   DataTypes::pointToString(const RealVectorType& data,const ShapeType& shape, int offset, const std::string& prefix)
   {
      using namespace std;
      ESYS_ASSERT(data.size()>0,"Error - Data object is empty.");
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


   void DataTypes::copyPoint(RealVectorType& dest, RealVectorType::size_type doffset, RealVectorType::size_type nvals, const RealVectorType& src, RealVectorType::size_type soffset)
   {
      ESYS_ASSERT((dest.size()>0&&src.size()>0&&checkOffset(doffset,dest.size(),nvals)),
                 "Error - Couldn't copy due to insufficient storage.");
      if (checkOffset(doffset,dest.size(),nvals) && checkOffset(soffset,src.size(),nvals)) {
         memcpy(&dest[doffset],&src[soffset],sizeof(real_t)*nvals);
      } else {
         throw DataException("Error - invalid offset specified.");
      }
   }

   void DataTypes::copyPoint(CplxVectorType& dest, CplxVectorType::size_type doffset, CplxVectorType::size_type nvals, const CplxVectorType& src, CplxVectorType::size_type soffset)
   {
      ESYS_ASSERT((dest.size()>0&&src.size()>0&&checkOffset(doffset,dest.size(),nvals)),
                 "Error - Couldn't copy due to insufficient storage.");
      if (checkOffset(doffset,dest.size(),nvals) && checkOffset(soffset,src.size(),nvals)) {
         memcpy(&dest[doffset],&src[soffset],sizeof(cplx_t)*nvals);
      } else {
         throw DataException("Error - invalid offset specified.");
      }
   }

   /**
    * \brief copy data from a real vector to a complex vector
    * The complex vector will be resized as needed and any previous
    * values will be replaced.
   */
   void DataTypes::fillComplexFromReal(const RealVectorType& r, CplxVectorType& c)
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

