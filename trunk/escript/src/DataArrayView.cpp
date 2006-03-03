// $Id$

/*
 ******************************************************************************
 *                                                                            *
 *       COPYRIGHT  ACcESS 2004 -  All Rights Reserved                        *
 *                                                                            *
 * This software is the property of ACcESS. No part of this code              *
 * may be copied in any form or by any means without the expressed written    *
 * consent of ACcESS.  Copying, use or modification of this software          *
 * by any unauthorised person is illegal unless that person has a software    *
 * license agreement with ACcESS.                                             *
 *                                                                            *
 ******************************************************************************
*/

#include "DataArrayView.h"
#include "DataException.h"

#include <sstream>

#include <boost/python/extract.hpp>

using namespace std;
using namespace boost::python;

namespace escript {

DataArrayView::DataArrayView():
    m_shape(),
    m_offset(0),
    m_noValues(0)
{
  m_data=0;
}

DataArrayView::DataArrayView(ValueType& data,
                             const ShapeType& viewShape,
                             int offset):
    m_data(&data),
    m_shape(viewShape),
    m_offset(offset),
    m_noValues(noValues(viewShape))
{
    //
    // check the shape rank and size and throw an exception if a
    // shape incompatible with the given data array is specified
    if (viewShape.size() > m_maxRank) {
      stringstream temp;
      temp << "Error- Couldn't construct DataArrayView, invalid rank."
	   << "Input data rank: " << viewShape.size() << " "
	   << "Maximum rank allowed for DataArrayView: " << m_maxRank;
      throw DataException(temp.str());
    }
    if (m_data->size() < (m_offset+noValues())) {
      stringstream temp;
      temp << "Error- Couldn't construct DataArrayView, insufficient data values. "
	   << "Shape requires: " << noValues()+m_offset << " values."
	   << "Data values available: " << m_data->size();
      throw DataException(temp.str());
    }
}

DataArrayView::DataArrayView(const DataArrayView& other):
    m_data(other.m_data),
    m_offset(other.m_offset),
    m_shape(other.m_shape),
    m_noValues(other.m_noValues)
{
}

bool
DataArrayView::isEmpty() const
{
  return (m_data==0);
}

void
DataArrayView::copy(const boost::python::numeric::array& value)
{
    ShapeType tempShape;    
    for (int i=0; i<value.getrank(); i++) {
      tempShape.push_back(extract<int>(value.getshape()[i]));
    }

    EsysAssert((!isEmpty()&&checkShape(tempShape)),
	       createShapeErrorMessage("Error - Couldn't copy due to shape mismatch.",tempShape));

    if (value.getrank()==0) {
      (*this)()=extract<double>(value[value.getshape()]);
    } else if (value.getrank()==1) {
      for (ValueType::size_type i=0;i<tempShape[0];i++) {
	(*this)(i)=extract<double>(value[i]);
      }
    } else if (value.getrank()==2) {
      for (ValueType::size_type i=0;i<tempShape[0];i++) {
	for (ValueType::size_type j=0;j<tempShape[1];j++) {
	  (*this)(i,j)=extract<double>(value[i][j]);
	}
      }
    } else if (value.getrank()==3) {
      for (ValueType::size_type i=0;i<tempShape[0];i++) {
	for (ValueType::size_type j=0;j<tempShape[1];j++) {
	  for (ValueType::size_type k=0;k<tempShape[2];k++) {
	    (*this)(i,j,k)=extract<double>(value[i][j][k]);
	  }
	}
      }
    } else if (value.getrank()==4) {
      for (ValueType::size_type i=0;i<tempShape[0];i++) {
	for (ValueType::size_type j=0;j<tempShape[1];j++) {
	  for (ValueType::size_type k=0;k<tempShape[2];k++) {
	    for (ValueType::size_type l=0;l<tempShape[3];l++) {
	      (*this)(i,j,k,l)=extract<double>(value[i][j][k][l]);
	    }
	  }
	}
      }
    }
}
 
void
DataArrayView::copy(const DataArrayView& other)
{
  copy(m_offset,other);
}

void
DataArrayView::copy(ValueType::size_type offset,
                    const DataArrayView& other)
{
    EsysAssert((!isEmpty()&&!other.isEmpty()&&checkOffset(offset)),
	       "Error - Couldn't copy due to insufficient storage.");
    EsysAssert((checkShape(other.getShape())),
	       createShapeErrorMessage("Error - Couldn't copy due to shape mismatch.",other.getShape()));
    if (checkOffset(offset)) {
      memcpy(&(*m_data)[offset],&(*other.m_data)[other.m_offset],sizeof(double)*noValues());
    } else {
      throw DataException("Error - invalid offset specified.");
    }
}

void
DataArrayView::copy(ValueType::size_type thisOffset, 
                    const DataArrayView& other, 
                    ValueType::size_type otherOffset)
{
    EsysAssert((!isEmpty()&&!other.isEmpty()&&checkOffset(thisOffset)&&other.checkOffset(otherOffset)),
	       "Error - Couldn't copy due to insufficient storage.");
    EsysAssert((checkShape(other.getShape())),
	       createShapeErrorMessage("Error - Couldn't copy due to shape mismatch.",other.getShape()));
    if (checkOffset(thisOffset)&&other.checkOffset(otherOffset)) {
      memcpy(&(*m_data)[thisOffset],&(*other.m_data)[otherOffset],sizeof(double)*noValues());
    } else {
      throw DataException("Error - invalid offset specified.");
    }
}

void
DataArrayView::copy(ValueType::size_type offset,
                    ValueType::value_type value)
{
    EsysAssert((!isEmpty()&&checkOffset(offset)),
	       "Error - Couldn't copy due to insufficient storage.");
    if (checkOffset(offset)) {
      vector<double> temp(noValues(),value);
      memcpy(&(*m_data)[offset],&temp[0],sizeof(double)*noValues());
    } else {
      throw DataException("Error - invalid offset specified.");
    }
}

int
DataArrayView::getRank() const
{
  return m_shape.size();
}

const
DataArrayView::ShapeType&
DataArrayView::getShape() const
{
  return m_shape;
}

int
DataArrayView::noValues(const ShapeType& shape) 
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
DataArrayView::noValues(const RegionLoopRangeType& region) 
{
  //
  // An empty region vector means rank 0 which contains 1 value
  int noValues=1;
  for (int i=0;i<region.size();i++) {
    noValues*=region[i].second-region[i].first;
  }
  return noValues;
}
 
int
DataArrayView::noValues() const
{
  return m_noValues;
}

bool
DataArrayView::checkShape(const DataArrayView::ShapeType& other) const
{
  return (m_shape==other);
}

string 
DataArrayView::createShapeErrorMessage(const string& messagePrefix,
                                       const DataArrayView::ShapeType& other) const
{
  stringstream temp;
  temp << messagePrefix
       << " This shape: " << shapeToString(m_shape)
       << " Other shape: " << shapeToString(other);
  return temp.str();
}

DataArrayView::ValueType::size_type
DataArrayView::getOffset() const
{
  return m_offset;
}

void
DataArrayView::setOffset(ValueType::size_type offset)
{
  EsysAssert((checkOffset(offset)), "Error - Invalid offset.");
  if (checkOffset(offset)) {
    m_offset=offset;
  } else {
    throw DataException("Error - invalid offset specified.");
  }
}

void
DataArrayView::incrOffset()
{
  EsysAssert((checkOffset(m_offset+noValues())), "Error - Cannot increment offset.");
  if (checkOffset(m_offset+noValues())) {
    m_offset=m_offset+noValues();
  } else {
    throw DataException("Error - Cannot increment offset.");
  }
}

bool
DataArrayView::checkOffset() const
{
  return checkOffset(m_offset);
}

bool
DataArrayView::checkOffset(ValueType::size_type offset) const
{
  return (m_data->size() >= (offset+noValues()));
}

DataArrayView::ValueType&
DataArrayView::getData() const
{
  EsysAssert(!isEmpty(),"Error - View is empty");
  return *m_data;
}

DataArrayView::ValueType::reference
DataArrayView::getData(ValueType::size_type i) const
{
  //
  // don't do any index checking so allow one past the end of the
  // vector providing the equivalent of end()
  return (*m_data)[m_offset+i];
}

DataArrayView::ShapeType
DataArrayView::getResultSliceShape(const RegionType& region)
{
  int dimSize;
  RegionType::const_iterator i;
  ShapeType result;
  for (i=region.begin();i!=region.end();i++) {
    dimSize=((i->second)-(i->first));
    if (dimSize!=0) {
      result.push_back(dimSize);
    }
  }
  return result;
}

DataArrayView::RegionType
DataArrayView::getSliceRegion(const boost::python::object& key) const
{
  int slice_rank, i;
  int this_rank=getRank();
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
           out[i]=getSliceRange(key[i],getShape()[i]);
         }
      }
  } else {
      slice_rank=1;
      if (slice_rank>this_rank) {
           throw DataException("Error - rank of slices does not match rank of slicee");
      } else {
           out[0]=getSliceRange(key,getShape()[0]);
      }
  }
  for (i=slice_rank;i<this_rank;i++) {
    out[i]=std::pair<int,int>(0,getShape()[i]);
  }
  return out;
}

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
  if (s0 == s1 and s1 >= shape)
     throw DataException("Error - slice index out of range.");
  if (s0 != s1 and  s1>shape)
     throw DataException("Error - slice index out of range.");
  if (s0 > s1) 
     throw DataException("Error - lower index must less or equal upper index.");
  return std::pair<int,int>(s0,s1);
}

DataArrayView::RegionLoopRangeType
getSliceRegionLoopRange(const DataArrayView::RegionType& region) 
{
    DataArrayView::RegionLoopRangeType region_loop_range(region.size());
    for (int i=0;i<region.size();i++) {
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

void
DataArrayView::copySlice(const DataArrayView& other, 
                         const RegionLoopRangeType& region)
{
  //
  // version of copySlice that uses the default offsets
  copySlice(m_offset,other,other.m_offset,region);
}

void
DataArrayView::copySlice(ValueType::size_type thisOffset,
                         const DataArrayView& other,
                         ValueType::size_type otherOffset,
                         const RegionLoopRangeType& region)
{

    //
    // Make sure views are not empty

    EsysAssert(!isEmpty(),
               "Error - this view is empty.");
    EsysAssert(!other.isEmpty(),
               "Error - other view is empty.");

    //
    // Check the view to be sliced from is compatible with the region to be sliced,
    // and that the region to be sliced is compatible with this view:

    EsysAssert(checkOffset(thisOffset),
               "Error - offset incompatible with this view.");
    EsysAssert(otherOffset+noValues()<=other.getData().size(),
               "Error - offset incompatible with other view.");

    EsysAssert(other.getRank()==region.size(),
               "Error - slice not same rank as view to be sliced from.");

    EsysAssert(noValues()==noValues(getResultSliceShape(region)),
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
      (*m_data)[thisOffset+numCopy]=(*other.m_data)[otherOffset+other.relIndex()];
      numCopy++;
      break;
    case 1:
      for (int i=region[0].first;i<region[0].second;i++) {
	(*m_data)[thisOffset+numCopy]=(*other.m_data)[otherOffset+other.relIndex(i)];
        numCopy++;
      }
      break;
    case 2:
      for (int j=region[1].first;j<region[1].second;j++) {
	for (int i=region[0].first;i<region[0].second;i++) {
	  (*m_data)[thisOffset+numCopy]=(*other.m_data)[otherOffset+other.relIndex(i,j)];
	  numCopy++;
	}
      }
      break;
    case 3:
      for (int k=region[2].first;k<region[2].second;k++) {
	for (int j=region[1].first;j<region[1].second;j++) {
	  for (int i=region[0].first;i<region[0].second;i++) {
	    (*m_data)[thisOffset+numCopy]=(*other.m_data)[otherOffset+other.relIndex(i,j,k)];
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
	      (*m_data)[thisOffset+numCopy]=(*other.m_data)[otherOffset+other.relIndex(i,j,k,l)];
	      numCopy++;
	    }
	  }
	}
      }
      break;
    default:
      stringstream mess;
      mess << "Error - (copySlice) Invalid slice region rank: " << region.size();
      throw DataException(mess.str());
    }
}
 
void
DataArrayView::copySliceFrom(const DataArrayView& other,
                             const RegionLoopRangeType& region_loop_range)
{
    //
    // version of copySliceFrom that uses the default offsets
    copySliceFrom(m_offset,other,other.m_offset,region_loop_range);
}

void
DataArrayView::copySliceFrom(ValueType::size_type thisOffset,
                             const DataArrayView& other,
                             ValueType::size_type otherOffset,
                             const RegionLoopRangeType& region)
{

    //
    // Make sure views are not empty

    EsysAssert(!isEmpty(),
               "Error - this view is empty.");
    EsysAssert(!other.isEmpty(),
               "Error - other view is empty.");

    //
    // Check this view is compatible with the region to be sliced,
    // and that the region to be sliced is compatible with the other view:

    EsysAssert(other.checkOffset(otherOffset),
               "Error - offset incompatible with other view.");
    EsysAssert(thisOffset+other.noValues()<=getData().size(),
               "Error - offset incompatible with this view.");

    EsysAssert(getRank()==region.size(),
               "Error - slice not same rank as this view.");

    EsysAssert(other.getRank()==0 || other.noValues()==noValues(getResultSliceShape(region)),
               "Error - slice shape not compatible shape for other view.");

    //
    // copy the values in the other view into the specified region of this view

    // allow for case where other view is a scalar
    if (other.getRank()==0) {

        // the following loops cannot be parallelised due to the numCopy counter
        int numCopy=0;

        switch (region.size()) {
        case 0:
        /* this case should never be encountered, 
           as python will never pass us an empty region.
           here for completeness only, allows slicing of a scalar */
          (*m_data)[thisOffset+relIndex()]=(*other.m_data)[otherOffset];
	  numCopy++;
          break;
        case 1:
          for (int i=region[0].first;i<region[0].second;i++) {
	    (*m_data)[thisOffset+relIndex(i)]=(*other.m_data)[otherOffset];
	    numCopy++;
          }
          break;
        case 2:
          for (int j=region[1].first;j<region[1].second;j++) {
	    for (int i=region[0].first;i<region[0].second;i++) {
	      (*m_data)[thisOffset+relIndex(i,j)]=(*other.m_data)[otherOffset];
	      numCopy++;
	    }
          }
          break;
        case 3:
          for (int k=region[2].first;k<region[2].second;k++) {
	    for (int j=region[1].first;j<region[1].second;j++) {
	      for (int i=region[0].first;i<region[0].second;i++) {
	        (*m_data)[thisOffset+relIndex(i,j,k)]=(*other.m_data)[otherOffset];
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
	          (*m_data)[thisOffset+relIndex(i,j,k,l)]=(*other.m_data)[otherOffset];
	          numCopy++;
	        }
	      }
	    }
          }
          break;
        default:
          stringstream mess;
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
          (*m_data)[thisOffset+relIndex()]=(*other.m_data)[otherOffset+numCopy];
	  numCopy++;
          break;
        case 1:
          for (int i=region[0].first;i<region[0].second;i++) {
	    (*m_data)[thisOffset+relIndex(i)]=(*other.m_data)[otherOffset+numCopy];
	    numCopy++;
          }
          break;
        case 2:
          for (int j=region[1].first;j<region[1].second;j++) {
	    for (int i=region[0].first;i<region[0].second;i++) {
	      (*m_data)[thisOffset+relIndex(i,j)]=(*other.m_data)[otherOffset+numCopy];
	      numCopy++;
	    }
          }
          break;
        case 3:
          for (int k=region[2].first;k<region[2].second;k++) {
	    for (int j=region[1].first;j<region[1].second;j++) {
	      for (int i=region[0].first;i<region[0].second;i++) {
	        (*m_data)[thisOffset+relIndex(i,j,k)]=(*other.m_data)[otherOffset+numCopy];
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
	          (*m_data)[thisOffset+relIndex(i,j,k,l)]=(*other.m_data)[otherOffset+numCopy];
	          numCopy++;
	        }
	      }
	    }
          }
          break;
        default:
          stringstream mess;
          mess << "Error - (copySliceFrom) Invalid slice region rank: " << region.size();
          throw DataException(mess.str());
        }

    }

}
 
string
DataArrayView::toString(const string& suffix) const
{
    EsysAssert(!isEmpty(),"Error - View is empty.");
    stringstream temp;
    string finalSuffix=suffix;
    if (suffix.length() > 0) {
      finalSuffix+=" ";
    }
    switch (getRank()) {
    case 0:
      temp << finalSuffix << (*this)();
      break;
    case 1:
      for (int i=0;i<getShape()[0];i++) {
	temp << "(" << i << ") " << finalSuffix << (*this)(i);
	if (i!=(getShape()[0]-1)) {
	  temp << endl;
	}
      }
      break;
    case 2:
      for (int i=0;i<getShape()[0];i++) {
	for (int j=0;j<getShape()[1];j++) {
	  temp << "(" << i << "," << j << ") " << finalSuffix << (*this)(i,j);
	  if (!(i==(getShape()[0]-1) && j==(getShape()[1]-1))) {
	    temp << endl;
	  }
	}
      }
      break;
    case 3:
      for (int i=0;i<getShape()[0];i++) {
	for (int j=0;j<getShape()[1];j++) {
	  for (int k=0;k<getShape()[2];k++) {
	    temp << "(" << i << "," << j << "," << k << ") " << finalSuffix << (*this)(i,j,k);
	    if (!(i==(getShape()[0]-1) && j==(getShape()[1]-1) && k==(getShape()[2]-1))) {
	      temp << endl;
	    }
	  }
	}
      }
      break;
    case 4:
      for (int i=0;i<getShape()[0];i++) {
	for (int j=0;j<getShape()[1];j++) {
	  for (int k=0;k<getShape()[2];k++) {
	    for (int l=0;l<getShape()[3];l++) {
	      temp << "(" << i << "," << j << "," << k << "," << l << ") " << finalSuffix << (*this)(i,j,k,l);
	      if (!(i==(getShape()[0]-1) && j==(getShape()[1]-1) && k==(getShape()[2]-1) && l==(getShape()[3]-1))) {
	        temp << endl;
	      }
            }
	  }
	}
      }
      break;
    default:
      stringstream mess;
      mess << "Error - (toString) Invalid rank: " << getRank();
      throw DataException(mess.str());
    }
    return temp.str();
}

string
DataArrayView::shapeToString(const DataArrayView::ShapeType& shape)
{
    stringstream temp;
    temp << "(";
    for (int i=0;i<shape.size();i++) {
      temp << shape[i];
      if (i < shape.size()-1) {
	temp << ",";
      }
    }
    temp << ")";
    return temp.str();
}

void 
DataArrayView::eigenvalues(const DataArrayView& in,
                           DataArrayView& ev)
{
   EsysAssert(!in.isEmpty(),
               "Error - in view is empty.");
   EsysAssert(!ev.isEmpty(),
               "Error - ev view is empty.");
   EsysAssert(in.getRank()==2,
               "Error - in has not rank 2.");
   EsysAssert(in.getShape()[0]==in.getShape()[1],
               "Error - in is not square");
   EsysAssert(ev.getRank()==1,
               "Error - ev has not rank 1.");
   EsysAssert(ev.getShape()[0]==in.getShape()[0],
               "Error - ev is too short");
   EsysAssert(0<ev.getShape()[0] && ev.getShape()[0]<4,
               "Error - dimension is restricted to 3.");
   double ev0,ev1,ev2;
   int s=in.getShape()[0];
   if (s==1) {
      eigenvalues1(in(0,0),&ev0);
      ev(0)=ev0;

   } else  if (s==2) {
      eigenvalues2(in(0,0),(in(0,1)+in(1,0))/2.,in(1,1),
                   &ev0,&ev1);
      ev(0)=ev0;
      ev(1)=ev1;

   } else  if (s==3) {
      eigenvalues3(in(0,0),(in(0,1)+in(1,0))/2.,(in(0,2)+in(2,0))/2.,in(1,1),(in(2,1)+in(1,2))/2.,in(2,2),
                 &ev0,&ev1,&ev2);
      ev(0)=ev0;
      ev(1)=ev1;
      ev(2)=ev2;

   }
}

void
DataArrayView::matMult(const DataArrayView& left, 
                       const DataArrayView& right, 
                       DataArrayView& result)
{

    if (left.getRank()==0 || right.getRank()==0) {
      stringstream temp;
      temp << "Error - (matMult) Invalid for rank 0 objects.";
      throw DataException(temp.str());
    }

    if (left.getShape()[left.getRank()-1] != right.getShape()[0]) {
      stringstream temp;
      temp << "Error - (matMult) Dimension: " << left.getRank() 
	   << ", size: " << left.getShape()[left.getRank()-1] 
	   << " of LHS and dimension: 1, size: " << right.getShape()[0]
	   << " of RHS don't match.";
      throw DataException(temp.str());
    }

    int outputRank = left.getRank()+right.getRank()-2;

    if (outputRank < 0) {
      stringstream temp;
      temp << "Error - (matMult) LHS and RHS cannot be multiplied "
	   << "as they have incompatable rank.";
      throw DataException(temp.str());
    }

    if (outputRank != result.getRank()) {
      stringstream temp;
      temp << "Error - (matMult) Rank of result array is: " 
	   << result.getRank() 
	   << " it must be: " << outputRank;
      throw DataException(temp.str());
    }

    for (int i=0; i<(left.getRank()-1); i++) {
      if (left.getShape()[i] != result.getShape()[i]) {
	stringstream temp;
	temp << "Error - (matMult) Dimension: " << i 
	     << " of LHS and result array don't match.";
	throw DataException(temp.str());
      }
    }

    for (int i=1; i<right.getRank(); i++) {
      if (right.getShape()[i] != result.getShape()[i+left.getRank()-2]) {
	stringstream temp;
	temp << "Error - (matMult) Dimension: " << i
	     << ", size: " << right.getShape()[i]
	     << " of RHS and dimension: " << i+left.getRank()-1 
	     << ", size: " << result.getShape()[i+left.getRank()-1]
	     << " of result array don't match.";
	throw DataException(temp.str());
      }
    }

    switch (left.getRank()) {

      case 1:
        switch (right.getRank()) {
          case 1:
            result()=0;
            for (int i=0;i<left.getShape()[0];i++) {
              result()+=left(i)*right(i);
            }
            break;
          case 2:
            for (int i=0;i<result.getShape()[0];i++) {
              result(i)=0;
              for (int j=0;j<right.getShape()[0];j++) {
                result(i)+=left(j)*right(j,i);
              }
            }
            break;
          default:
            stringstream temp; temp << "Error - (matMult) Invalid rank. Programming error.";
            throw DataException(temp.str());
            break;
        }
        break;

      case 2:
        switch (right.getRank()) {
          case 1:
            result()=0;
              for (int i=0;i<left.getShape()[0];i++) {
               result(i)=0;
               for (int j=0;j<left.getShape()[1];j++) {
                 result(i)+=left(i,j)*right(i);
               }
            }
	    break;
          case 2:
            for (int i=0;i<result.getShape()[0];i++) {
              for (int j=0;j<result.getShape()[1];j++) {
                result(i,j)=0;
                for (int jR=0;jR<right.getShape()[0];jR++) {
                  result(i,j)+=left(i,jR)*right(jR,j);
                }
              }
            }
            break;
          default:
            stringstream temp; temp << "Error - (matMult) Invalid rank. Programming error.";
            throw DataException(temp.str());
            break;
        }
        break;

      default:
        stringstream temp; temp << "Error - (matMult) Not supported for rank: " << left.getRank();
        throw DataException(temp.str());
        break;
    }

}

DataArrayView::ShapeType
DataArrayView::determineResultShape(const DataArrayView& left,
                                    const DataArrayView& right)
{
    ShapeType temp;
    for (int i=0; i<(left.getRank()-1); i++) {
      temp.push_back(left.getShape()[i]);
    }
    for (int i=1; i<right.getRank(); i++) {
      temp.push_back(right.getShape()[i]);
    }
    return temp;
}

bool operator==(const DataArrayView& left, const DataArrayView& right)
{
    //
    // The views are considered equal if they have the same shape and each value
    // of the view is the same.
    bool result=(left.m_shape==right.m_shape);
    if (result) {
      for (int i=0;i<left.noValues();i++) {
	result=result && (left.getData(i)==right.getData(i));
      }
    }
    return result;
}

bool operator!=(const DataArrayView& left, const DataArrayView& right)
{
   return (!operator==(left,right));
}

}  // end of namespace
