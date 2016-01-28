#include "DataVectorAlt.h"

/* This file exists to provide a custom implementation of complex methods for DataVectorAlt
   It also explicitly instantiates the complex version of the template to ensure linkage
*/   

namespace escript
{
  
namespace DataTypes
{

// Please make sure that any implementation changes here are reflected in the generic version in the .h file
template<>
void 
DataVectorAlt<DataTypes::cplx_t>::copyFromArrayToOffset(const WrappedArray& value, size_type offset, size_type copies)
{
  const DataTypes::ShapeType& tempShape=value.getShape();
  size_type len=DataTypes::noValues(tempShape);
  if (offset+len*copies>size())
  {
     std::ostringstream ss;
     ss << "Error - not enough room for that DataPoint at that offset. (";
     ss << "offset=" << offset << " + " << " len=" << len << " >= " << size();
     throw DataException(ss.str());
  }
  size_type si=0,sj=0,sk=0,sl=0;
  switch (value.getRank())
  {
  case 0:	
	for (size_type z=0;z<copies;++z)
	{
	   m_array_data[offset+z]=value.getEltC();
	}
	break;
  case 1:
	for (size_type z=0;z<copies;++z)
	{
	   for (size_t i=0;i<tempShape[0];++i)
	   {
	      m_array_data[offset+i]=value.getEltC(i);
	   }
	   offset+=len;
	}
	break;
  case 2:
	si=tempShape[0];
	sj=tempShape[1];
	for (size_type z=0;z<copies;++z)
	{
           for (size_type i=0;i<si;i++)
	   {
              for (size_type j=0;j<sj;j++)
	      {
                 m_array_data[offset+DataTypes::getRelIndex(tempShape,i,j)]=value.getEltC(i,j);
              }
           }
	   offset+=len;
	}
	break;
  case 3:
	si=tempShape[0];
	sj=tempShape[1];
	sk=tempShape[2];
	for (size_type z=0;z<copies;++z) 
	{
          for (size_type i=0;i<si;i++)
	  {
            for (size_type j=0;j<sj;j++)
	    {
              for (size_type k=0;k<sk;k++)
	      {
                 m_array_data[offset+DataTypes::getRelIndex(tempShape,i,j,k)]=value.getEltC(i,j,k);
              }
            }
          }
	  offset+=len;
	}
	break;
  case 4:
	si=tempShape[0];
	sj=tempShape[1];
	sk=tempShape[2];
	sl=tempShape[3];
	for (size_type z=0;z<copies;++z)
	{
          for (size_type i=0;i<si;i++)
	  {
            for (size_type j=0;j<sj;j++)
	    {
              for (size_type k=0;k<sk;k++)
	      {
                 for (size_type l=0;l<sl;l++)
		 {
                    m_array_data[offset+DataTypes::getRelIndex(tempShape,i,j,k,l)]=value.getEltC(i,j,k,l);
                 }
              }
            }
          }
	  offset+=len;
	}
	break;
  default:
	std::ostringstream oss;
	oss << "Error - unknown rank. Rank=" << value.getRank();
	throw DataException(oss.str());
  }
}

template class DataVectorAlt<DataTypes::cplx_t>;

}	// end namespace
}	// end namespace
