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

#include "escript/Data/DataArrayView.h" 
#include "escript/Data/DataException.h"

#include <sstream>
#include <vector>
#include <iostream>
#include <boost/python/extract.hpp>
#include <boost/python/object.hpp>

using namespace std;
using namespace boost::python;

namespace escript {

DataArrayView::DataArrayView():
    m_shape(),
    m_offset(),
    m_noValues()
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
    // check the shape and throw an exception if an illegal shape is specified
    if (viewShape.size() > m_maxRank) {
      stringstream temp;
      temp << "Error- Couldn't construct DataArrayView, invalid rank." 
	   << "Input data rank: " << viewShape.size() 
	   << " Maximum rank allowed for DataArrayView: " 
	   << m_maxRank;
      throw DataException(temp.str());
    }
    if (m_data->size() < (noValues()+m_offset)) {
      stringstream temp;
      temp << "Error- Couldn't construct DataArrayView, insufficient storage." 
	   << " Shape requires: " << noValues()+m_offset << " values." 
	   << " Data values available: " << m_data->size();
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
    return(m_data==0);
}

void
DataArrayView::copy(const boost::python::numeric::array& value)
{
    DataArrayView::ShapeType tempShape;    
    for (int i=0; i<value.getrank(); ++i) {
      tempShape.push_back(extract<int>(value.getshape()[i]));
    }
    EsysAssert((!isEmpty()&&checkShape(tempShape)),
	       createShapeErrorMessage("Error - Couldn't copy due to shape mismatch.",tempShape));

    if (value.getrank()==0) {
      (*this)()=extract<double>(value[value.getshape()]);
    } else if (value.getrank()==1) {
      for (ValueType::size_type i=0;i<tempShape[0];++i) {
	(*this)(i)=extract<double>(value[i]);
      }
    } else if (value.getrank()==2) {
      for (ValueType::size_type i=0;i<tempShape[0];++i) {
	for (ValueType::size_type j=0;j<tempShape[1];++j) {
	  (*this)(i,j)=extract<double>(value[i][j]);
	}
      }
    } else if (value.getrank()==3) {
      for (ValueType::size_type i=0;i<tempShape[0];++i) {
	for (ValueType::size_type j=0;j<tempShape[1];++j) {
	  for (ValueType::size_type k=0;k<tempShape[2];++k) {
	    (*this)(i,j,k)=extract<double>(value[i][j][k]);
	  }
	}
      }
    } else if (value.getrank()==4) {
      for (ValueType::size_type i=0;i<tempShape[0];++i) {
	for (ValueType::size_type j=0;j<tempShape[1];++j) {
	  for (ValueType::size_type k=0;k<tempShape[2];++k) {
	    for (ValueType::size_type m=0;m<tempShape[3];++m) {
	      (*this)(i,j,k,m)=extract<double>(value[i][j][k][m]);
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
    EsysAssert((!isEmpty()&&checkShape(other.getShape())),
	       createShapeErrorMessage("Error - Couldn't copy due to shape mismatch.",other.getShape()));
    memcpy(&(*m_data)[offset],&(*other.m_data)[other.m_offset],sizeof(double)*noValues());
}

void
DataArrayView::copy(ValueType::size_type thisOffset, 
                    const DataArrayView& other, 
                    ValueType::size_type otherOffset)
{
    EsysAssert((!isEmpty()&&checkShape(other.getShape())),
	       createShapeErrorMessage("Error - Couldn't copy due to shape mismatch.",other.getShape()));
    memcpy(&(*m_data)[thisOffset],&(*other.m_data)[otherOffset],sizeof(double)*noValues());
}

void
DataArrayView::copy(ValueType::size_type offset,
                    ValueType::value_type value)
{
    //
    // fill the entire view with the single value
    EsysAssert(!isEmpty(),"Error - View is empty");
    ValueType temp(noValues(),value);
    memcpy(&(*m_data)[offset],&temp[0],sizeof(double)*noValues());
}

int
DataArrayView::getRank() const
{
    return m_shape.size();
}

const DataArrayView::ShapeType&
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
    for (i=shape.begin();i!=shape.end();++i) {
      noValues*=(*i);
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
 
void
DataArrayView::setOffset(ValueType::size_type offset)
{
    EsysAssert((m_data->size() >= (noValues()+offset)), "Invalid offset");
    m_offset=offset;
}

DataArrayView::ShapeType
DataArrayView::getResultSliceShape(const RegionType& region)
{
    int dimSize;
    RegionType::const_iterator i;
    DataArrayView::ShapeType result;
    for (i=region.begin();i!=region.end();i++) {
      dimSize=((i->second) - (i->first));
      if (dimSize!=0) {
	result.push_back(dimSize);
      }
    }
    return result;
}

DataArrayView::RegionType
DataArrayView::getSliceRegion(const boost::python::object& key) const
{
  int slice_rank, out_rank, i;
  int this_rank=getRank();
  /* allow for case where key is singular eg: [1], this implies we
     want to generate a rank-1 dimension object, as opposed to eg: [1,2]
     which implies we want to take a rank dimensional object with one
     dimension of size 1 */
  extract<tuple> key_tuple(key);
  if (key_tuple.check()) {
      slice_rank=extract<int> (key.attr("__len__")());
  } else {
      slice_rank=1;
  }
  /* ensure slice is correctly dimensioned */
  if (slice_rank>=this_rank) {
       out_rank=slice_rank;
  } else {
       out_rank=this_rank;
  }
  DataArrayView::RegionType out(out_rank);
  /* calculate values for slice range */
  if (key_tuple.check()) {
    for (i=0;i<slice_rank;i++) {
      out[i]=getSliceRange(key[i],getShape()[i]);
    }
  } else {
    out[0]=getSliceRange(key,getShape()[0]);
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
  return std::pair<int,int>(s0,s1);
}

void
DataArrayView::copySlice(const DataArrayView& other, 
                         const RegionType& region)
{
    //
    // version of copySlice that uses the default offsets
    copySlice(m_offset,other,other.m_offset,region);
}

void
DataArrayView::copySlice(ValueType::size_type thisOffset,
                         const DataArrayView& other,
                         ValueType::size_type otherOffset,
                         const RegionType& region)
{

    //
    // Make sure this object is not empty
    EsysAssert(!isEmpty(),"Error - this data object is empty.");

    //
    // Modify region to copy from in order to
    // deal with the case where one range in the region contains identical indexes,
    // eg: <<1,1><0,3><0,3>>
    // This situation implies we want to copy from an object with rank greater than that of this
    // object. eg: we want to copy the values from a two dimensional slice out of a three
    // dimensional object into a two dimensional object.
    // We do this by taking a slice from the other object where one dimension of
    // the slice region is of size 1. So in the above example, we modify the above 
    // region like so: <<1,2><0,3><0,3>> and take this slice.

    RegionType slice_region=region;
    for (int i=0;i<slice_region.size();i++) {
      if (slice_region[i].first==slice_region[i].second) {
        (slice_region[i].second)++;
      }
    }

    //
    // Check the object to be sliced from is compatible with the region to be sliced,
    // and that the region to be sliced is compatible with this object:
    // 1: the region to be sliced should be the same rank as the other object.
    // 2: the other object should have at least as many datapoints at the region requires.
    // 3: the region should have the same number of datapoints as this object

    DataArrayView::ShapeType slice_shape=getResultSliceShape(slice_region);
    EsysAssert(other.getRank()==slice_region.size(),
               "Error - slice rank incompatible with object to be sliced from.");
    EsysAssert((noValues(slice_shape)<=other.noValues()),
               "Error - slice shape incompatible with object to be sliced from.");
    EsysAssert((noValues(slice_shape)==noValues()),
               "Error - slice shape incompatible with this object.");

    //
    // copy the values in the specified slice_region of the other object into this object

    int numCopy=0;
    switch (slice_region.size()) {
    case 0:
      /* this case should never be encountered, as python will never pass us an empty slice */
      (*m_data)[thisOffset]=(*other.m_data)[otherOffset];
      break;
    case 1:
      for (int i=slice_region[0].first;i<slice_region[0].second;i++) {
	(*m_data)[numCopy+thisOffset]=(*other.m_data)[i+otherOffset];
	numCopy++;
      }
      break;
    case 2:
      for (int j=slice_region[1].first;j<slice_region[1].second;j++) {
	for (int i=slice_region[0].first;i<slice_region[0].second;i++) {
	  (*m_data)[numCopy+thisOffset]=(*other.m_data)[other.relIndex(i,j)+otherOffset];
	  numCopy++;
	}
      }
      break;
    case 3:
      for (int k=slice_region[2].first;k<slice_region[2].second;k++) {
	for (int j=slice_region[1].first;j<slice_region[1].second;j++) {
	  for (int i=slice_region[0].first;i<slice_region[0].second;i++) {
	    (*m_data)[numCopy+thisOffset]=(*other.m_data)[other.relIndex(i,j,k)+otherOffset];
	    numCopy++;
	  }
	}
      }
      break;
    case 4:
      for (int l=slice_region[3].first;l<slice_region[3].second;l++) {
	for (int k=slice_region[2].first;k<slice_region[2].second;k++) {
	  for (int j=slice_region[1].first;j<slice_region[1].second;j++) {
	    for (int i=slice_region[0].first;i<slice_region[0].second;i++) {
	      (*m_data)[numCopy+thisOffset]=(*other.m_data)[other.relIndex(i,j,k,l)+otherOffset];
	      numCopy++;
	    }
	  }
	}
      }
      break;
    default:
      stringstream mess;
      mess << "Error - (copySlice) Invalid rank: " << slice_region.size() << " for region.";
      throw DataException(mess.str());
    }
}
 
void
DataArrayView::copySliceFrom(const DataArrayView& other,
                             const RegionType& region)
{
    //
    // version of slice that uses the default offset
    copySliceFrom(m_offset,other,other.m_offset,region);
}

void
DataArrayView::copySliceFrom(ValueType::size_type thisOffset,
                             const DataArrayView& other,
                             ValueType::size_type otherOffset,
                             const RegionType& region)
{

    //
    // Make sure this object is not empty
    EsysAssert(!isEmpty(),"Error - this data object is empty.");

    //
    // Modify region to copy in order to
    // deal with the case where one range in the region contains identical idexes,
    // eg: <<1,1><0,3><0,3>>
    // This situation implies we want to copy from an object with rank less than that of this
    // object. eg: we want to copy the values from a two dimensional object into a two dimensional
    // slice of a three dimensional object.
    // We do this by copying to a slice from this object where one dimension of
    // the slice region is of size 1. So in the above example, we modify the above 
    // region like so: <<1,2><0,3><0,3>> and copy into this slice.

    RegionType slice_region=region;
    for (int i=0;i<slice_region.size();i++) {
      if (slice_region[i].first==slice_region[i].second) {
        (slice_region[i].second)++;
      }
    }

    //
    // Check the object to be coped from is compatible with the slice region to be copied
    // into, and that the region to be sliced into is compatible with this object:
    // 1: the region to be sliced into should be the same rank as this object.
    // 2: the other object should have exactly the same number of data points as the
    //    region to be sliced into requires.
    // 3: the region can be a higher rank than the other object, provided 2 holds.
    // 4: this object must have at least as many data points as the region requires.

    //
    // Check region size and shape
    DataArrayView::ShapeType slice_shape=getResultSliceShape(slice_region);
    EsysAssert((getRank()==slice_region.size()),"Error - Invalid slice region.");
    EsysAssert((noValues(slice_shape)==other.noValues()),
               "Error - slice shape incompatible with object to be sliced from.");
    EsysAssert(other.getRank()<=slice_region.size(),
               "Error - slice rank incompatible with object to be sliced from.");
    EsysAssert((noValues(slice_shape)<=noValues()),
               "Error - slice shape incompatible with this object.");

    //
    // copy the values in the other object into the specified slice_region of this object

    int numCopy=0;
    switch (slice_region.size()) {
    case 0:
      /* this case should never be encountered, as python will never pass us an empty slice */
      (*m_data)[thisOffset]=(*other.m_data)[otherOffset];
      break;
    case 1:
      for (int i=slice_region[0].first;i<slice_region[0].second;i++) {
	(*m_data)[i+thisOffset]=(*other.m_data)[numCopy+otherOffset];
	numCopy++;
      }
      break;
    case 2:
      for (int j=slice_region[1].first;j<slice_region[1].second;j++) {
	for (int i=slice_region[0].first;i<slice_region[0].second;i++) {
	  (*m_data)[relIndex(i,j)+thisOffset]=(*other.m_data)[numCopy+otherOffset];
	  numCopy++;
	}
      }
      break;
    case 3:
      for (int k=slice_region[2].first;k<slice_region[2].second;k++) {
	for (int j=slice_region[1].first;j<slice_region[1].second;j++) {
	  for (int i=slice_region[0].first;i<slice_region[0].second;i++) {
	    (*m_data)[relIndex(i,j,k)+thisOffset]=(*other.m_data)[numCopy+otherOffset];
	    numCopy++;
	  }
	}
      }
      break;
    case 4:
      for (int l=slice_region[3].first;l<slice_region[3].second;l++) {
	for (int k=slice_region[2].first;k<slice_region[2].second;k++) {
	  for (int j=slice_region[1].first;j<slice_region[1].second;j++) {
	    for (int i=slice_region[0].first;i<slice_region[0].second;i++) {
	      (*m_data)[relIndex(i,j,k,l)+thisOffset]=(*other.m_data)[numCopy+otherOffset];
	      numCopy++;
	    }
	  }
	}
      }
      break;
    default:
      stringstream mess;
      mess << "Error - (copySliceFrom) Invalid rank: " << slice_region.size() << " for region.";
      throw DataException(mess.str());
    }
}
 
DataArrayView::ShapeType
DataArrayView::determineResultShape(const DataArrayView& left,
                                    const DataArrayView& right)
{
    DataArrayView::ShapeType temp;
    for (int i=0; i<(left.getRank()-1); i++) {
      temp.push_back(left.getShape()[i]);
    }
    for (int i=1; i<right.getRank(); i++) {
      temp.push_back(right.getShape()[i]);
    }
    return temp;
}

string
DataArrayView::toString(const string& suffix) const
{
    EsysAssert(!isEmpty(),"Error - View is empty");
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
DataArrayView::matMult(const DataArrayView& left, 
                       const DataArrayView& right, 
                       DataArrayView& result)
{
    if (left.getRank()==0 || right.getRank()==0) {
      stringstream temp;
      temp << "Error - (matMult) Invalid for rank 0 objects.";
      throw DataException(temp.str());
    }

    if (left.getShape()[left.getRank()-1]!=right.getShape()[0]) {
      stringstream temp;
      temp << "Error - (matMult) Dimension: " << left.getRank() 
	   << ", size: " << left.getShape()[left.getRank()-1] 
	   << " of LHS and dimension: 1, size: " << right.getShape()[0]
	   << " of RHS don't match.";
      throw DataException(temp.str());
    }
    int outputRank=left.getRank()+right.getRank()-2;
    
    if (outputRank < 0) {
      stringstream temp;
      temp << "Error - (matMult) Left and right operands cannot be multiplied "
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

    for (int i=0; i<(left.getRank()-1); ++i) {
      if (left.getShape()[i]!=result.getShape()[i]) {
	stringstream temp;
	temp << "Error - (matMult) Dimension: " << i 
	     << " of LHS and output array don't match.";
	throw DataException(temp.str());
      }
    }

    for (int i=1; i<right.getRank(); ++i) {
      if (right.getShape()[i] != result.getShape()[i+left.getRank()-2]) {
	stringstream temp;
	temp << "Error - (matMult) Dimension: " << i
	     << ", size: " << right.getShape()[i]
	     << " of RHS and dimension: " << i+left.getRank()-1 
	     << ", size: " << result.getShape()[i+left.getRank()-1]
	     << " of output array don't match.";
	throw DataException(temp.str());
      }
    }

    switch (left.getRank()) {
    case 1:
      switch (right.getRank()) {
      case 1:
	result()=0;
	for (int i=0;i<left.getShape()[0];++i) {
	  result()+=left(i)*right(i);
	}
	break;
      case 2:
        for (int i=0;i<result.getShape()[0];++i) {
	  result(i)=0;
	  for (int j=0;j<right.getShape()[0];++j) {
	    result(i)+=left(j)*right(j,i);
	  }
	}
        break;
      default:
	stringstream temp;
	temp << "Error - (matMult) Invalid rank. Programming error.";
        throw DataException(temp.str());
        break;
      }
      break;
    case 2:
      switch (right.getRank()) {
      case 1:
	result()=0;
	for (int i=0;i<left.getShape()[0];++i) {
	  result(i)=0;
	  for (int j=0;j<left.getShape()[1];++j) {
	    result(i)+=left(i,j)*right(i);
	  }
	}
	break;
      case 2:
        for (int i=0;i<result.getShape()[0];++i) {
	  for (int j=0;j<result.getShape()[1];++j) {
	    result(i,j)=0;
	    for (int jR=0;jR<right.getShape()[0];++jR) {
	      result(i,j)+=left(i,jR)*right(jR,j);
	    }
	  }
	}
        break;
      default:
	stringstream temp;
	temp << "Error - (matMult) Invalid rank. Programming error.";
        throw DataException(temp.str());
        break;
      }
      break;
    default:
      stringstream temp;
      temp << "Error - (matMult) Not supported for rank: " << left.getRank();
      throw DataException(temp.str());
      break;
    }
}

bool operator==(const DataArrayView& left, const DataArrayView& right)
{
    //
    // The views are considered equal if they have the same shape and each value
    // of the view is the same.
    bool result=(left.m_shape==right.m_shape);
    if (result) {
      for (int i=0;i<left.noValues();++i) {
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
