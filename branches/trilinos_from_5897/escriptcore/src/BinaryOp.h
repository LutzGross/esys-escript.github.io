
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

#ifndef __ESCRIPT_BINARYOP_H__
#define __ESCRIPT_BINARYOP_H__

#include "system_dep.h"
#include "DataTypes.h"
#include "DataConstant.h"
#include "DataExpanded.h"
#include "DataMaths.h"
#include "DataTagged.h"

/**
\file BinaryOp.h 
\brief Describes binary operations performed on instances of DataAbstract.

For operations on DataVector see DataMaths.h.
For operations on double* see LocalOps.h.
*/

namespace escript {

template <class LVEC, class RVEC>
inline void binaryOpDataReadyHelperTC(DataTagged& left, const DataConstant& right, 
		     escript::ESFunction operation)
{
  //
  // perform the operation on each tagged value
  const DataTagged::DataMapType& lookup=left.getTagLookup();
  DataTagged::DataMapType::const_iterator i;
  DataTagged::DataMapType::const_iterator lookupEnd=lookup.end();
  
  typename LVEC::ElementType dummyL=0;
  
  LVEC& leftVec=left.getTypedVectorRW(dummyL);
  const DataTypes::ShapeType& leftShape=left.getShape();
  const DataTypes::ShapeType& rightShape=right.getShape();
  typename RVEC::ElementType rvalue=0;
  rvalue=right.getTypedVectorRO(rvalue)[0];		// for rank==0
  const RVEC& rightVec=right.getTypedVectorRO(rvalue);   // for rank>0
  if (right.getRank()==0) {
    for (i=lookup.begin();i!=lookupEnd;i++) {
      DataMaths::binaryOpVector(leftVec,leftShape,i->second,rvalue,operation);
    }
  } else {
    for (i=lookup.begin();i!=lookupEnd;i++) {
      DataMaths::binaryOpVector(leftVec, leftShape, i->second,rightVec,rightShape,0,operation);
    }
  }
  //
  // finally perform the operation on the default value
  if (right.getRank()==0) {
    DataMaths::binaryOpVector(leftVec,leftShape,left.getDefaultOffset(),rvalue,operation);
  } else {
    DataMaths::binaryOpVector(leftVec,leftShape,left.getDefaultOffset(),rightVec,rightShape,0,operation);
  }
}

  
/**
   \brief
   Perform the given binary operation.
   \param left Input/Output - The left hand side.
   \param right Input - The right hand side.
   \param operation Input - The operation to perform.
*/
inline void binaryOpDataReady(DataTagged& left, const DataConstant& right, 
		     escript::ESFunction operation)
{
  if (left.isComplex())
  {
      if (right.isComplex())
      {
	  binaryOpDataReadyHelperTC<DataTypes::CplxVectorType, DataTypes::CplxVectorType>(left, right, operation);
      }
      else
      {
	  binaryOpDataReadyHelperTC<DataTypes::CplxVectorType, DataTypes::RealVectorType>(left, right, operation);	
      }    
  }
  else	// left is real
  {
      if (right.isComplex())
      {
	  throw DataException("Programming error: binaryOpDataReady - LHS is real but RHS is complex");
      }
      else	// right is real
      {
	  binaryOpDataReadyHelperTC<DataTypes::RealVectorType, DataTypes::RealVectorType>(left, right, operation);	
      }        
  }
}


/**
   \brief apply the binary op to each value in left and the single value right.

   The value in right will be assumed to begin at offset 0
*/
template <class LVEC, class RVEC>
inline void binaryOpDataReadyHelperTV(DataTagged& left, const RVEC& right, 
		     const DataTypes::ShapeType& shape,
		     escript::ESFunction operation)
{
  //
  // perform the operation on each tagged value
  const DataTagged::DataMapType& lookup=left.getTagLookup();
  DataTagged::DataMapType::const_iterator i;
  DataTagged::DataMapType::const_iterator lookupEnd=lookup.end();
  typename LVEC::ElementType dummy=0;
  LVEC& lvec=left.getTypedVectorRW(dummy);
  const DataTypes::ShapeType& lshape=left.getShape();
  if (DataTypes::getRank(shape)==0) {
    for (i=lookup.begin();i!=lookupEnd;i++) {
      DataMaths::binaryOpVector(lvec, lshape,i->second,right[0],operation);
    }
  } else {
    for (i=lookup.begin();i!=lookupEnd;i++) {
      DataMaths::binaryOpVector(lvec, lshape, i->second,right,shape,0,operation);
    }
  }
  //
  // finally perform the operation on the default value
  if (DataTypes::getRank(shape)==0) {
    DataMaths::binaryOpVector(lvec,lshape,left.getDefaultOffset(),right[0],operation);
  } else {
    DataMaths::binaryOpVector(lvec,lshape,left.getDefaultOffset(),right, shape,0,operation);
  }
}




/**
   \brief apply the binary op to each value in left and the single value right.

   The value in right will be assumed to begin at offset 0
*/
inline void binaryOpDataReady(DataTagged& left, const DataTypes::RealVectorType& right, 
		     const DataTypes::ShapeType& shape,
		     escript::ESFunction operation)
{
  if (left.isComplex())
  {
      binaryOpDataReadyHelperTV<DataTypes::CplxVectorType, DataTypes::RealVectorType>(left, right, shape, operation);	
  }
  else	// left is real
  {
      binaryOpDataReadyHelperTV<DataTypes::RealVectorType, DataTypes::RealVectorType>(left, right, shape,  operation);	        
  }
}


inline void binaryOpDataReady(DataTagged& left, const DataTypes::CplxVectorType& right, 
		     const DataTypes::ShapeType& shape,
		     escript::ESFunction operation)
{
  if (left.isComplex())
  {
      binaryOpDataReadyHelperTV<DataTypes::CplxVectorType, DataTypes::CplxVectorType>(left, right, shape, operation);
  }
  else	// left is real
  {
      throw DataException("Programming error: binaryOpDataReady - LHS is real but RHS is complex");
  }
}



template <class LVEC, class RVEC>
inline void binaryOpDataReadyHelperTT(DataTagged& left, const DataTagged& right, 
		     escript::ESFunction operation)
{
  using namespace DataMaths;

  int right_rank=right.getRank();
  //
  // Add the right hand tag keys which can't currently be found on the left
  const DataTagged::DataMapType& rightLookup=right.getTagLookup();
  DataTagged::DataMapType::const_iterator i;
  DataTagged::DataMapType::const_iterator rightLookupEnd=rightLookup.end();
  for (i=rightLookup.begin();i!=rightLookupEnd;i++) {
    //
    // If the left does not already have a value assigned to this tag,
    // add the right hand tag to the left hand tag list and assign
    // the left's default value.
    if (!left.isCurrentTag(i->first)) {
      left.addTag(i->first);
    }
  }
  typename LVEC::ElementType dummyl=0;
  typename RVEC::ElementType dummyr=0;
  LVEC& leftVec=left.getTypedVectorRW(dummyl);
  const DataTypes::ShapeType& leftShape=left.getShape();
  //
  // Perform the operation.
  const DataTagged::DataMapType& leftLookup=left.getTagLookup();
  DataTagged::DataMapType::const_iterator leftLookupEnd=leftLookup.end();
  for (i=leftLookup.begin();i!=leftLookupEnd;i++) {
    if (right_rank==0) {
       binaryOpVector(leftVec,leftShape,i->second, right.getDataByTagRO(i->first,0, dummyr),operation);

    } else {	// rank>0
       binaryOpVector(leftVec,leftShape,left.getOffsetForTag(i->first),right.getTypedVectorRO(dummyr), right.getShape(), right.getOffsetForTag(i->first), operation);
    }
  }
  //
  // finally perform the operation on the default value
  if (right_rank==0) {
     binaryOpVector(leftVec,leftShape, left.getDefaultOffset(), right.getTypedVectorRO(dummyr)[0],operation);
  } else {
     binaryOpVector(leftVec,leftShape, left.getDefaultOffset(), right.getTypedVectorRO(dummyr), right.getShape(), right.getDefaultOffset(), operation);
  }
}


inline void binaryOpDataReady(DataTagged& left, const DataTagged& right, 
		     escript::ESFunction operation)
{
  if (left.isComplex())
  {
      if (right.isComplex())
      {
	  binaryOpDataReadyHelperTT<DataTypes::CplxVectorType, DataTypes::CplxVectorType>(left, right, operation);
      }
      else
      {
	  binaryOpDataReadyHelperTT<DataTypes::CplxVectorType, DataTypes::RealVectorType>(left, right, operation);	
      }    
  }
  else	// left is real
  {
      if (right.isComplex())
      {
	  throw DataException("Programming error: binaryOpDataReady - LHS is real but RHS is complex");
      }
      else	// right is real
      {
	  binaryOpDataReadyHelperTT<DataTypes::RealVectorType, DataTypes::RealVectorType>(left, right, operation);	
      }        
  }
}

template <class LSCALAR, class RSCALAR>
inline void binaryOpDataReadyHelperCC(DataConstant& left, const DataConstant& right, 
		     escript::ESFunction operation)
{
  LSCALAR dummyl=0;
  RSCALAR dummyr=0;
  if (right.getRank()==0) {
    RSCALAR r=right.getTypedVectorRO(dummyr)[0];
    DataMaths::binaryOpVector(left.getTypedVectorRW(dummyl), left.getShape(),0, r,operation);
  } else {
    DataMaths::binaryOpVector(left.getTypedVectorRW(dummyl), left.getShape(),0, right.getTypedVectorRO(dummyr),right.getShape(),0,operation);
  }
}

void binaryOpDataCCC(DataConstant& result, const DataConstant& left, const DataConstant& right, 
		     escript::ESFunction operation);
void binaryOpDataTCT(DataTagged& result, const DataConstant& left, const DataTagged& right, 
		     escript::ESFunction operation);
void binaryOpDataTTC(DataTagged& result, const DataTagged& left, const DataConstant& right, 
		     escript::ESFunction operation);
void binaryOpDataTTT(DataTagged& result, const DataTagged& left, const DataTagged& right, 
		     escript::ESFunction operation);
void binaryOpDataEEC(DataExpanded& result, const DataExpanded& left, const DataConstant& right, 
		     escript::ESFunction operation);
void binaryOpDataECE(DataExpanded& result, const DataConstant& left, const DataExpanded& right, 
		     escript::ESFunction operation);
void binaryOpDataEEE(DataExpanded& result, const DataExpanded& left, const DataExpanded& right, 
		     escript::ESFunction operation);
void binaryOpDataETE(DataExpanded& result, const DataTagged& left, const DataExpanded& right, 
		     escript::ESFunction operation);
void binaryOpDataEET(DataExpanded& result, const DataExpanded& left, const DataTagged& right, 
 		     escript::ESFunction operation);


inline void binaryOpDataReady(DataConstant& left, const DataConstant& right, 
		     escript::ESFunction operation)
{
  if (left.isComplex())
  {
      if (right.isComplex())
      {
	  binaryOpDataReadyHelperCC<DataTypes::cplx_t, DataTypes::cplx_t>(left, right, operation);
      }
      else
      {
	  binaryOpDataReadyHelperCC<DataTypes::cplx_t, DataTypes::real_t>(left, right, operation);	
      }    
  }
  else	// left is real
  {
      if (right.isComplex())
      {
	  throw DataException("Programming error: binaryOpDataReady - LHS is real but RHS is complex");
      }
      else	// right is real
      {
	  binaryOpDataReadyHelperCC<DataTypes::real_t, DataTypes::real_t>(left, right, operation);	
      }        
  }  
}


template <typename LVEC, typename RVEC>
inline void binaryOpDataReadyHelperER(DataExpanded& left, const DataReady& right, 
		     escript::ESFunction operation)
{
  typename LVEC::ElementType dummyl=0;
  typename RVEC::ElementType dummyr=0;
  int i,j;
  typename LVEC::size_type numDPPSample=left.getNumDPPSample();
  typename LVEC::size_type numSamples=left.getNumSamples();
  if (right.getRank()==0) {

    const DataTypes::ShapeType& leftShape=left.getShape();
    LVEC& leftVec=left.getTypedVectorRW(dummyl);
    //
    // This will call the double version of binaryOp
    #pragma omp parallel for private(i,j) schedule(static)
    for (i=0;i<numSamples;i++) {
      for (j=0;j<numDPPSample;j++) {
	DataMaths::binaryOpVector(leftVec,leftShape,left.getPointOffset(i,j), right.getTypedVectorRO(dummyr)[right.getPointOffset(i,j)]  ,operation);
      }
    }
  } else {
    #pragma omp parallel for private(i,j) schedule(static)
    for (i=0;i<numSamples;i++) {
      for (j=0;j<numDPPSample;j++) {
	DataMaths::binaryOpVector(left.getTypedVectorRW(dummyl),left.getShape(),left.getPointOffset(i,j), right.getTypedVectorRO(dummyr), right.getShape(),right.getPointOffset(i,j), operation);
      }
    }
  }
}   

inline void binaryOpDataReady(DataExpanded& left, const DataReady& right, 
		     escript::ESFunction operation)
{
  if (left.isComplex())
  {
      if (right.isComplex())
      {
	  binaryOpDataReadyHelperER<DataTypes::CplxVectorType, DataTypes::CplxVectorType>(left, right, operation);
      }
      else
      {
	  binaryOpDataReadyHelperER<DataTypes::CplxVectorType, DataTypes::RealVectorType>(left, right, operation);	
      }    
  }
  else	// left is real
  {
      if (right.isComplex())
      {
	  throw DataException("Programming error: binaryOpDataReady - LHS is real but RHS is complex");
      }
      else	// right is real
      {
	  binaryOpDataReadyHelperER<DataTypes::RealVectorType, DataTypes::RealVectorType>(left, right, operation);	
      }        
  }    
}  



// -----------------------------------

/**
   \brief
   Perform the given binary operation.
   \param left Input/Output - The left hand side.
   \param right Input - The right hand side.
   \param operation Input - The operation to perform.
*/
template <class BinaryFunction>
inline void binaryOp(DataTagged& left, const DataConstant& right, 
		     BinaryFunction operation)
{
//  binaryOp(left,right.getPointDataView(),operation);
  //
  // perform the operation on each tagged value
  const DataTagged::DataMapType& lookup=left.getTagLookup();
  DataTagged::DataMapType::const_iterator i;
  DataTagged::DataMapType::const_iterator lookupEnd=lookup.end();
  DataTypes::RealVectorType& leftVec=left.getVectorRW();
  const DataTypes::ShapeType& leftShape=left.getShape();
  const DataTypes::ShapeType& rightShape=right.getShape();
  double rvalue=right.getVectorRO()[0];		// for rank==0
  const DataTypes::RealVectorType& rightVec=right.getVectorRO();   // for rank>0
  if (right.getRank()==0) {
    for (i=lookup.begin();i!=lookupEnd;i++) {
      DataMaths::binaryOp(leftVec,leftShape,i->second,rvalue,operation);
    }
  } else {
    for (i=lookup.begin();i!=lookupEnd;i++) {
      DataMaths::binaryOp(leftVec, leftShape, i->second,rightVec,rightShape,0,operation);
    }
  }
  //
  // finally perform the operation on the default value
  if (right.getRank()==0) {
    DataMaths::binaryOp(leftVec,leftShape,left.getDefaultOffset(),rvalue,operation);
  } else {
    DataMaths::binaryOp(leftVec,leftShape,left.getDefaultOffset(),rightVec,rightShape,0,operation);
  }
}

/**
   \brief apply the binary op to each value in left and the single value right.

   The value in right will be assumed to begin at offset 0
*/
template <class BinaryFunction>
inline void binaryOp(DataTagged& left, const DataTypes::RealVectorType& right, 
		     const DataTypes::ShapeType& shape,
		     BinaryFunction operation)
{
  //
  // perform the operation on each tagged value
  const DataTagged::DataMapType& lookup=left.getTagLookup();
  DataTagged::DataMapType::const_iterator i;
  DataTagged::DataMapType::const_iterator lookupEnd=lookup.end();
  DataTypes::RealVectorType& lvec=left.getVectorRW();
  const DataTypes::ShapeType& lshape=left.getShape();
  if (DataTypes::getRank(shape)==0) {
    for (i=lookup.begin();i!=lookupEnd;i++) {
      DataMaths::binaryOp(lvec, lshape,i->second,right[0],operation);
    }
  } else {
    for (i=lookup.begin();i!=lookupEnd;i++) {
      DataMaths::binaryOp(lvec, lshape, i->second,right,shape,0,operation);
    }
  }
  //
  // finally perform the operation on the default value
  if (DataTypes::getRank(shape)==0) {
    DataMaths::binaryOp(lvec,lshape,left.getDefaultOffset(),right[0],operation);
  } else {
    DataMaths::binaryOp(lvec,lshape,left.getDefaultOffset(),right, shape,0,operation);
  }
}




template <class BinaryFunction>
inline void binaryOp(DataTagged& left, const DataTagged& right, 
		     BinaryFunction operation)
{
  using namespace DataMaths;

  int right_rank=right.getRank();
  //
  // Add the right hand tag keys which can't currently be found on the left
  const DataTagged::DataMapType& rightLookup=right.getTagLookup();
  DataTagged::DataMapType::const_iterator i;
  DataTagged::DataMapType::const_iterator rightLookupEnd=rightLookup.end();
  for (i=rightLookup.begin();i!=rightLookupEnd;i++) {
    //
    // If the left does not already have a value assigned to this tag,
    // add the right hand tag to the left hand tag list and assign
    // the left's default value.
    if (!left.isCurrentTag(i->first)) {
      left.addTag(i->first);
    }
  }
  DataTypes::RealVectorType& leftVec=left.getVectorRW();
  const DataTypes::ShapeType& leftShape=left.getShape();
  //
  // Perform the operation.
  const DataTagged::DataMapType& leftLookup=left.getTagLookup();
  DataTagged::DataMapType::const_iterator leftLookupEnd=leftLookup.end();
  for (i=leftLookup.begin();i!=leftLookupEnd;i++) {
    if (right_rank==0) {
       binaryOp(leftVec,leftShape,i->second, right.getDataByTagRO(i->first,0),operation);

    } else {	// rank>0
       binaryOp(leftVec,leftShape,left.getOffsetForTag(i->first),right.getVectorRO(), right.getShape(), right.getOffsetForTag(i->first), operation);
    }
  }
  //
  // finally perform the operation on the default value
  if (right_rank==0) {
     binaryOp(leftVec,leftShape, left.getDefaultOffset(), right.getVectorRO()[0],operation);
  } else {
     binaryOp(leftVec,leftShape, left.getDefaultOffset(), right.getVectorRO(), right.getShape(), right.getDefaultOffset(), operation);
  }
}

template <class BinaryFunction>
inline void binaryOp(DataConstant& left, const DataConstant& right, 
		     BinaryFunction operation)
{
  if (right.getRank()==0) {
    double r=right.getVectorRO()[0];
    DataMaths::binaryOp(left.getVectorRW(), left.getShape(),0, r,operation);
  } else {
    DataMaths::binaryOp(left.getVectorRW(), left.getShape(),0, right.getVectorRO(),right.getShape(),0,operation);
  }

}



template <class BinaryFunction>
inline void binaryOp(DataExpanded& left, const DataReady& right, 
		     BinaryFunction operation)
{
  int i,j;
  DataTypes::RealVectorType::size_type numDPPSample=left.getNumDPPSample();
  DataTypes::RealVectorType::size_type numSamples=left.getNumSamples();
  if (right.getRank()==0) {

    const DataTypes::ShapeType& leftShape=left.getShape();
    DataTypes::RealVectorType& leftVec=left.getVectorRW();
    //
    // This will call the double version of binaryOp
    #pragma omp parallel for private(i,j) schedule(static)
    for (i=0;i<numSamples;i++) {
      for (j=0;j<numDPPSample;j++) {
	DataMaths::binaryOp(leftVec,leftShape,left.getPointOffset(i,j), right.getVectorRO()[right.getPointOffset(i,j)]  ,operation);
      }
    }
  } else {
    #pragma omp parallel for private(i,j) schedule(static)
    for (i=0;i<numSamples;i++) {
      for (j=0;j<numDPPSample;j++) {
	DataMaths::binaryOp(left.getVectorRW(),left.getShape(),left.getPointOffset(i,j), right.getVectorRO(), right.getShape(),right.getPointOffset(i,j), operation);
      }
    }
  }
}


} // end of namespace

#endif // __ESCRIPT_BINARYOP_H__

