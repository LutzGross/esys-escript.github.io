
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




#include "BinaryDataReadyOps.h"
#include "DataTagged.h"

#include <sstream>
using namespace escript;
using namespace std;

namespace escript
{

template <class ResSCALAR, class LSCALAR, class RSCALAR>
inline void binaryOpDataReadyHelperCCC(DataConstant& res, const DataConstant& left, const DataConstant& right, 
		     escript::ES_optype operation)
{
  ResSCALAR resdummy=0;
  LSCALAR dummyl=0;
  RSCALAR dummyr=0;
  DataTypes::RealVectorType::size_type valcount=1*DataTypes::noValues(res.getShape());	// since we are constant, by definition
	  // there is only one datapoint
  
  
  if (right.getRank()==0) 
  {		
    escript::binaryOpVectorRightScalar(res.getTypedVectorRW(resdummy), 0, 1, valcount, 
			      left.getTypedVectorRO(dummyl), 0,
			      &right.getTypedVectorRO(dummyr)[0], true,
			      operation, true);
  }
  else if (left.getRank()==0)
  {
    escript::binaryOpVectorLeftScalar(res.getTypedVectorRW(resdummy), 0, 1, valcount, 
			      &left.getTypedVectorRO(dummyl)[0], true,		// left is const so it only has one sample of one data point (and from the if we know that sample is rank0)
			      right.getTypedVectorRO(dummyr), 0,
			      operation, true);
  }
  else
  {
    escript::binaryOpVector(res.getTypedVectorRW(resdummy), 0, 1, valcount, 
			      left.getTypedVectorRO(dummyl), 0, false,
			      right.getTypedVectorRO(dummyr), 0, false,
			      operation);
  }
}


void binaryOpDataCCC(DataConstant& result, const DataConstant& left, const DataConstant& right, 
		     escript::ES_optype operation)
{
  bool cplxresult=left.isComplex() || right.isComplex();
  if (result.isComplex()!=cplxresult)
  {
      ostringstream oss;
      oss << "Programming error: result has unexpected complexity ";
      oss << result.isComplex() << "==" << left.isComplex() << "||";
      oss << right.isComplex();
      throw DataException(oss.str());
  }
  
  if (left.isComplex())
  {
      if (right.isComplex())
      {
	  binaryOpDataReadyHelperCCC<DataTypes::cplx_t, DataTypes::cplx_t, DataTypes::cplx_t>(result, left, right, operation);
      }
      else
      {
	  binaryOpDataReadyHelperCCC<DataTypes::cplx_t, DataTypes::cplx_t, DataTypes::real_t>(result, left, right, operation);	
      }    
  }
  else	// left is real
  {
      if (right.isComplex())
      {
	  binaryOpDataReadyHelperCCC<DataTypes::cplx_t, DataTypes::real_t, DataTypes::cplx_t>(result, left, right, operation);	
      }
      else	// right is real
      {
	  binaryOpDataReadyHelperCCC<DataTypes::real_t, DataTypes::real_t, DataTypes::real_t>(result, left, right, operation);	
      }        
  }    
}

template <class ResSCALAR, class LSCALAR, class RSCALAR>
inline void binaryOpDataReadyHelperTCT(DataTagged& res, const DataConstant& left, const DataTagged& right, 
		     escript::ES_optype operation)
{
  ResSCALAR resdummy=0;
  LSCALAR dummyl=0;
  RSCALAR dummyr=0;
  DataTypes::RealVectorType::size_type valcount=DataTypes::noValues(res.getShape());	// there is only one datapoint per sample

  // self update is not a possibility here because res and left are different types
  if (res.getTagCount()!=0)		// no tags
  {
      throw DataException("Programming error: result must have no tags for binaryOpDataReadyTCT");
  }	

  if (res.getTagCount()==0)
  {
      const DataTagged::DataMapType& lookup_1=right.getTagLookup();
      DataTagged::DataMapType::const_iterator i; // i->first is a tag, i->second is an offset into memory
      for (i=lookup_1.begin();i!=lookup_1.end();i++)
      {
	res.addTag(i->first);
      }
  }
  // so now we know that both tagged objects have the same tags (perhaps not in the same order though)
  if (right.getRank()==0) 	// scalar op on the right
  {		
      // This will process the default value (which we know is stored in location 0)
      escript::binaryOpVectorRightScalar(res.getTypedVectorRW(resdummy), 0, 
				1, valcount,
			      left.getTypedVectorRO(dummyl), 0,
			      &right.getTypedVectorRO(dummyr)[0], 0,
			      operation, false);
      const DataTagged::DataMapType& lookup_re=res.getTagLookup();
      DataTagged::DataMapType::const_iterator i; // i->first is a tag, i->second is an offset into memory
      for (i=lookup_re.begin();i!=lookup_re.end();i++)
      {
	  DataTypes::RealVectorType::size_type resoffset=i->second;
	  DataTypes::RealVectorType::size_type rightoffset=right.getOffsetForTag(i->first);
	  escript::binaryOpVectorRightScalar(res.getTypedVectorRW(resdummy), resoffset, 
				    1, valcount,
				  left.getTypedVectorRO(dummyl), 0,
				  &right.getTypedVectorRO(dummyr)[rightoffset], 0,
				  operation, false);	  
      }        
  }
  else if (left.getRank()==0)	// scalar op on the left
  {
      escript::binaryOpVectorLeftScalar(res.getTypedVectorRW(resdummy), 0, 
				1, valcount,
			      &left.getTypedVectorRO(dummyl)[0], 0,
			      right.getTypedVectorRO(dummyr), 0,
			      operation, false);
      const DataTagged::DataMapType& lookup_re=res.getTagLookup();
      DataTagged::DataMapType::const_iterator i; // i->first is a tag, i->second is an offset into memory
      for (i=lookup_re.begin();i!=lookup_re.end();i++)
      {
	  DataTypes::RealVectorType::size_type resoffset=i->second;
	  DataTypes::RealVectorType::size_type rightoffset=right.getOffsetForTag(i->first);	  
	  escript::binaryOpVectorLeftScalar(res.getTypedVectorRW(resdummy), resoffset, 
				    1, valcount, 
				  &left.getTypedVectorRO(dummyl)[0], 0,
				  right.getTypedVectorRO(dummyr), rightoffset,
				  operation, false);	  
      }        
  }
  else
  {
      // This will process the default value (which we know is stored in location 0)
      escript::binaryOpVector(res.getTypedVectorRW(resdummy), 0, 1, valcount, 
			      left.getTypedVectorRO(dummyl), 0, true,
			      right.getTypedVectorRO(dummyr), 0, false,
			      operation);
      const DataTagged::DataMapType& lookup_1=right.getTagLookup();
      DataTagged::DataMapType::const_iterator i; // i->first is a tag, i->second is an offset into memory
      for (i=lookup_1.begin();i!=lookup_1.end();i++)
      {
	  DataTypes::RealVectorType::size_type resoffset=right.getOffsetForTag(i->first);
	  escript::binaryOpVector(res.getTypedVectorRW(resdummy), resoffset, 1, valcount, 
				  left.getTypedVectorRO(dummyl), 0, true,
				  right.getTypedVectorRO(dummyr), i->second, false,
				  operation);	  
      }      
      
  }
}


void binaryOpDataTCT(DataTagged& result, const DataConstant& left, const DataTagged& right, 
		     escript::ES_optype operation)
{
  bool cplxresult=left.isComplex() || right.isComplex();
  if (result.isComplex()!=cplxresult)
  {
      ostringstream oss;
      oss << "Programming error: result has unexpected complexity ";
      oss << result.isComplex() << "==" << left.isComplex() << "||";
      oss << right.isComplex();
      throw DataException(oss.str());
  }
  
  if (left.isComplex())
  {
      if (right.isComplex())
      {
	  binaryOpDataReadyHelperTCT<DataTypes::cplx_t, DataTypes::cplx_t, DataTypes::cplx_t>(result, left, right, operation);
      }
      else
      {
	  binaryOpDataReadyHelperTCT<DataTypes::cplx_t, DataTypes::cplx_t, DataTypes::real_t>(result, left, right, operation);	
      }    
  }
  else	// left is real
  {
      if (right.isComplex())
      {
	  binaryOpDataReadyHelperTCT<DataTypes::cplx_t, DataTypes::real_t, DataTypes::cplx_t>(result, left, right, operation);
      }
      else	// right is real
      {
	  binaryOpDataReadyHelperTCT<DataTypes::real_t, DataTypes::real_t, DataTypes::real_t>(result, left, right, operation);	
      }        
  }    
}


template <class ResSCALAR, class LSCALAR, class RSCALAR>
inline void binaryOpDataReadyHelperECE(DataExpanded& res, const DataConstant& left, const DataExpanded& right, 
		     escript::ES_optype operation)
{
  ResSCALAR resdummy=0;
  LSCALAR dummyl=0;
  RSCALAR dummyr=0;
  DataTypes::RealVectorType::size_type valcount=res.getNumDPPSample()*DataTypes::noValues(res.getShape());
  
  // if both sides are rank0, then that should be handled normally rather than with a special case
  // hence we check for that possibility first

  if (right.getRank()==left.getRank())		// both zero or both equal and non-zero 
  {
    escript::binaryOpVector(res.getTypedVectorRW(resdummy), 0, 
			      res.getNumSamples()*res.getNumDPPSample(),DataTypes::noValues(res.getShape()) , 
			      left.getTypedVectorRO(dummyl), 0, true,
			      right.getTypedVectorRO(dummyr), 0, false,
			      operation);
  } else if (right.getRank()==0) 
  {
      // This is a tricky one. There are lots of individual scalars on the RHS, each of which need to be 
      // multiplied by a number of values which make up a single const
      // fiddling by pretending samples are smaller but there are more of them
    escript::binaryOpVectorRightScalar(res.getTypedVectorRW(resdummy), 0, 
					 res.getNumSamples()*res.getNumDPPSample(), DataTypes::noValues(res.getShape()), 
			      left.getTypedVectorRO(dummyl), 0,
			      &right.getTypedVectorRO(dummyr)[0], false,
			      operation,
			      true);
  }
  else  // if (left.getRank()==0)
  {
    escript::binaryOpVectorLeftScalar(res.getTypedVectorRW(resdummy), 0, right.getNumSamples(), valcount, 
			      &left.getTypedVectorRO(dummyl)[0], true,
			      right.getTypedVectorRO(dummyr), false,
			      operation,
			      false);
  }
}


void binaryOpDataECE(DataExpanded& result, const DataConstant& left, const DataExpanded& right, 
		     escript::ES_optype operation)
{
  bool cplxresult=left.isComplex() || right.isComplex();
  if (result.isComplex()!=cplxresult)
  {
      ostringstream oss;
      oss << "Programming error: result has unexpected complexity ";
      oss << result.isComplex() << "==" << left.isComplex() << "||";
      oss << right.isComplex();
      throw DataException(oss.str());
  }
  
  if (left.isComplex())
  {
      if (right.isComplex())
      {
	  binaryOpDataReadyHelperECE<DataTypes::cplx_t, DataTypes::cplx_t, DataTypes::cplx_t>(result, left, right, operation);
      }
      else
      {
	  binaryOpDataReadyHelperECE<DataTypes::cplx_t, DataTypes::cplx_t, DataTypes::real_t>(result, left, right, operation);	
      }    
  }
  else	// left is real
  {
      if (right.isComplex())
      {
	  binaryOpDataReadyHelperECE<DataTypes::cplx_t, DataTypes::real_t, DataTypes::cplx_t>(result, left, right, operation);
      }
      else	// right is real
      {
	  binaryOpDataReadyHelperECE<DataTypes::real_t, DataTypes::real_t, DataTypes::real_t>(result, left, right, operation);	
      }        
  }    
}

template <class ResSCALAR, class LSCALAR, class RSCALAR>
inline void binaryOpDataReadyHelperEET(DataExpanded& res, const DataExpanded& left, const DataTagged& right, 
		     escript::ES_optype operation)
{
  ResSCALAR resdummy=0;
  LSCALAR dummyl=0;
  RSCALAR dummyr=0;

    escript::binaryOpVectorTagged(res.getTypedVectorRW(resdummy), 
			      res.getNumSamples(),res.getNumDPPSample(), DataTypes::noValues(res.getShape()), 
			      left.getTypedVectorRO(dummyl), left.getRank()==0,
			      right.getTypedVectorRO(dummyr), right.getRank()==0,
			      false,	// right object is the tagged one
			      right,	// source of tags
			      operation);  
  
}


void binaryOpDataEET(DataExpanded& result, const DataExpanded& left, const DataTagged& right, 
		     escript::ES_optype operation)
{
  bool cplxresult=left.isComplex() || right.isComplex();
  if (result.isComplex()!=cplxresult)
  {
      ostringstream oss;
      oss << "Programming error: result has unexpected complexity ";
      oss << result.isComplex() << "==" << left.isComplex() << "||";
      oss << right.isComplex();
      throw DataException(oss.str());
  }
  
  if (left.isComplex())
  {
      if (right.isComplex())
      {
	  binaryOpDataReadyHelperEET<DataTypes::cplx_t, DataTypes::cplx_t, DataTypes::cplx_t>(result, left, right, operation);
      }
      else
      {
	  binaryOpDataReadyHelperEET<DataTypes::cplx_t, DataTypes::cplx_t, DataTypes::real_t>(result, left, right, operation);	
      }    
  }
  else	// left is real
  {
      if (right.isComplex())
      {
	  binaryOpDataReadyHelperEET<DataTypes::cplx_t, DataTypes::real_t, DataTypes::cplx_t>(result, left, right, operation);	
      }
      else	// right is real
      {
	  binaryOpDataReadyHelperEET<DataTypes::real_t, DataTypes::real_t, DataTypes::real_t>(result, left, right, operation);	
      }        
  }    
}


template <class ResSCALAR, class LSCALAR, class RSCALAR>
inline void binaryOpDataReadyHelperETE(DataExpanded& res, const DataTagged& left, const DataExpanded& right, 
		     escript::ES_optype operation)
{
  ResSCALAR resdummy=0;
  LSCALAR dummyl=0;
  RSCALAR dummyr=0;
  escript::binaryOpVectorTagged(res.getTypedVectorRW(resdummy), 
			      res.getNumSamples(),res.getNumDPPSample(), DataTypes::noValues(res.getShape()), 
			      left.getTypedVectorRO(dummyl), left.getRank()==0,
			      right.getTypedVectorRO(dummyr), right.getRank()==0,
			      true,	// left object is the tagged one
			      left,	// source of tags
			      operation);  
  
}


void binaryOpDataETE(DataExpanded& result, const DataTagged& left, const DataExpanded& right, 
		     escript::ES_optype operation)
{
  bool cplxresult=left.isComplex() || right.isComplex();
  if (result.isComplex()!=cplxresult)
  {
      ostringstream oss;
      oss << "Programming error: result has unexpected complexity ";
      oss << result.isComplex() << "==" << left.isComplex() << "||";
      oss << right.isComplex();
      throw DataException(oss.str());
  }
  
  if (left.isComplex())
  {
      if (right.isComplex())
      {
	  binaryOpDataReadyHelperETE<DataTypes::cplx_t, DataTypes::cplx_t, DataTypes::cplx_t>(result, left, right, operation);
      }
      else
      {
	  binaryOpDataReadyHelperETE<DataTypes::cplx_t, DataTypes::cplx_t, DataTypes::real_t>(result, left, right, operation);	
      }    
  }
  else	// left is real
  {
      if (right.isComplex())
      {
	  binaryOpDataReadyHelperETE<DataTypes::cplx_t, DataTypes::real_t, DataTypes::cplx_t>(result, left, right, operation);	
      }
      else	// right is real
      {
	  binaryOpDataReadyHelperETE<DataTypes::real_t, DataTypes::real_t, DataTypes::real_t>(result, left, right, operation);	
      }        
  }    
}




template <class ResSCALAR, class LSCALAR, class RSCALAR>
inline void binaryOpDataReadyHelperTTC(DataTagged& res, const DataTagged& left, const DataConstant& right, 
		     escript::ES_optype operation)
{
  ResSCALAR resdummy=0;
  LSCALAR dummyl=0;
  RSCALAR dummyr=0;
  DataTypes::RealVectorType::size_type valcount=DataTypes::noValues(res.getShape());	// there is only one datapoint per sample

  // We need to consider two possibilities:
  //   1) we are dealing with a new result object (which won't have tags)
  //   2) we are storing the result back into the same object (eg +=)
  // for case 1, we need to add tags in the correct order and then calculate on a per tag basis
  // for case 2, we just need to calculate tags

  
  // first let's exclude anything but our two cases
  if ((&res!=&left) &&		// self update 
    (res.getTagCount()!=0))		// no tags
  {
      throw DataException("binaryOpDataReadyTTC expects a=(a op b) or c=(a op b)");
  }	

  if (res.getTagCount()==0)
  {
      const DataTagged::DataMapType& lookup_1=left.getTagLookup();
      DataTagged::DataMapType::const_iterator i; // i->first is a tag, i->second is an offset into memory
      for (i=lookup_1.begin();i!=lookup_1.end();i++)
      {
	res.addTag(i->first);
      }
  }
  // so now we know that both tagged objects have the same tags (perhaps not in the same order though)
  if (right.getRank()==0) 	// scalar op on the right
  {
    
        // This will process the default value (which we know is stored in location 0)
      escript::binaryOpVectorRightScalar(res.getTypedVectorRW(resdummy), 0, 
				1, valcount,
			      left.getTypedVectorRO(dummyl), 0,
			      &right.getTypedVectorRO(dummyr)[0], 0,
			      operation, false);
      const DataTagged::DataMapType& lookup_re=res.getTagLookup();
      DataTagged::DataMapType::const_iterator i; // i->first is a tag, i->second is an offset into memory
      for (i=lookup_re.begin();i!=lookup_re.end();i++)
      {
	  DataTypes::RealVectorType::size_type resoffset=i->second;
	  DataTypes::RealVectorType::size_type leftoffset=left.getOffsetForTag(i->first);
	  escript::binaryOpVectorRightScalar(res.getTypedVectorRW(resdummy), resoffset, 
				    1, valcount,
				  left.getTypedVectorRO(dummyl), leftoffset,
				  &right.getTypedVectorRO(dummyr)[0], 0,
				  operation, false);     
      }  
  }
  else if (left.getRank()==0)	// scalar op on the left
  {
      escript::binaryOpVectorLeftScalar(res.getTypedVectorRW(resdummy), 0, 
				1, valcount,
			      &left.getTypedVectorRO(dummyl)[0], 0,
			      right.getTypedVectorRO(dummyr), 0,
			      operation, false);
      const DataTagged::DataMapType& lookup_re=res.getTagLookup();
      DataTagged::DataMapType::const_iterator i; // i->first is a tag, i->second is an offset into memory
      for (i=lookup_re.begin();i!=lookup_re.end();i++)
      {
	  DataTypes::RealVectorType::size_type resoffset=i->second;
	  DataTypes::RealVectorType::size_type leftoffset=left.getOffsetForTag(i->first);	  
	  escript::binaryOpVectorLeftScalar(res.getTypedVectorRW(resdummy), resoffset, 
				    1, valcount, 
				  &left.getTypedVectorRO(dummyl)[leftoffset], 0,
				  right.getTypedVectorRO(dummyr), 0,
				  operation, false);	  
      }       
  }
  else
  {
      // This will process the default value (which we know is stored in location 0)
      escript::binaryOpVector(res.getTypedVectorRW(resdummy), 0, 1, valcount, 
			      left.getTypedVectorRO(dummyl), 0, true,
			      right.getTypedVectorRO(dummyr), 0, false,
			      operation);
      const DataTagged::DataMapType& lookup_re=res.getTagLookup();
      DataTagged::DataMapType::const_iterator i; // i->first is a tag, i->second is an offset into memory
      for (i=lookup_re.begin();i!=lookup_re.end();i++)
      {
	  DataTypes::RealVectorType::size_type resoffset=i->second;
	  DataTypes::RealVectorType::size_type leftoffset=left.getOffsetForTag(i->first);	  
	  escript::binaryOpVector(res.getTypedVectorRW(resdummy), resoffset, 1, valcount, 
				  left.getTypedVectorRO(dummyl), leftoffset, true,
				  right.getTypedVectorRO(dummyr), 0, false,
				  operation);	  
      }      
  }  
}


void binaryOpDataTTC(DataTagged& result, const DataTagged& left, const DataConstant& right, 
		     escript::ES_optype operation)
{
  bool cplxresult=left.isComplex() || right.isComplex();
  if (result.isComplex()!=cplxresult)
  {
      ostringstream oss;
      oss << "Programming error: result has unexpected complexity ";
      oss << result.isComplex() << "==" << left.isComplex() << "||";
      oss << right.isComplex();
      throw DataException(oss.str());
  }
  
  if (left.isComplex())
  {
      if (right.isComplex())
      {
	  binaryOpDataReadyHelperTTC<DataTypes::cplx_t, DataTypes::cplx_t, DataTypes::cplx_t>(result, left, right, operation);
      }
      else
      {
	  binaryOpDataReadyHelperTTC<DataTypes::cplx_t, DataTypes::cplx_t, DataTypes::real_t>(result, left, right, operation);	
      }    
  }
  else	// left is real
  {
      if (right.isComplex())
      {
	  binaryOpDataReadyHelperTTC<DataTypes::cplx_t, DataTypes::real_t, DataTypes::cplx_t>(result, left, right, operation);	
      }
      else	// right is real
      {
	  binaryOpDataReadyHelperTTC<DataTypes::real_t, DataTypes::real_t, DataTypes::real_t>(result, left, right, operation);	
      }        
  }    
}


template <class ResSCALAR, class LSCALAR, class RSCALAR>
inline void binaryOpDataReadyHelperTTT(DataTagged& res, const DataTagged& left, const DataTagged& right, 
		     escript::ES_optype operation)
{
  ResSCALAR resdummy=0;
  LSCALAR dummyl=0;
  RSCALAR dummyr=0;
  DataTypes::RealVectorType::size_type valcount=DataTypes::noValues(res.getShape());	// there is only one datapoint per sample

  // We need to consider two possibilities:
  //   1) we are dealing with a new result object (which won't have tags)
  //   2) we are storing the result back into the same object (eg +=)
  // for case 1, we need to add tags in the correct order and then calculate on a per tag basis
  // for case 2, we just need to calculate tags

  
  // first let's exclude anything but our two cases
  if ((&res!=&left) &&		// self update 
    (res.getTagCount()!=0))		// no tags
  {
      throw DataException("binaryOpDataReadyTTT expects a=(a op b) or c=(a op b)");
  }	

  // add tags from both sides
  if (res.getTagCount()==0)
  {
      const DataTagged::DataMapType& lookup_1=left.getTagLookup();
      DataTagged::DataMapType::const_iterator i; // i->first is a tag, i->second is an offset into memory
      for (i=lookup_1.begin();i!=lookup_1.end();i++)
      {
	res.addTag(i->first);
      }

      const DataTagged::DataMapType& lookup_r=right.getTagLookup();
      for (i=lookup_r.begin();i!=lookup_r.end();i++)
      {
	res.addTag(i->first);
      }
  }
  else	// result already has tags in it
  {	// add tags from right, any duplicates are silently ignored by addTag
      const DataTagged::DataMapType& lookup_r=right.getTagLookup();
      for (auto i=lookup_r.begin();i!=lookup_r.end();i++)
      {
	res.addTag(i->first);
      }      
  }
  // so now we know that both tagged objects have the same tags (perhaps not in the same order though)
  if (right.getRank()==0) 	// scalar op on the right
  {		// we'll reuse this code by pretending samples are 1 value long
      // This will process the default value (which we know is stored in location 0)
      escript::binaryOpVector(res.getTypedVectorRW(resdummy), 0, 
				valcount, 1,	// (arguments reversed from normal) 
			      left.getTypedVectorRO(dummyl), 0, false,
			      right.getTypedVectorRO(dummyr), 0, true,
			      operation);
      const DataTagged::DataMapType& lookup_re=res.getTagLookup();
      DataTagged::DataMapType::const_iterator i; // i->first is a tag, i->second is an offset into memory
      for (i=lookup_re.begin();i!=lookup_re.end();i++)
      {
	  DataTypes::RealVectorType::size_type resoffset=i->second;
	  DataTypes::RealVectorType::size_type leftoffset=left.getOffsetForTag(i->first);
	  DataTypes::RealVectorType::size_type rightoffset=right.getOffsetForTag(i->first);	
	  escript::binaryOpVector(res.getTypedVectorRW(resdummy), resoffset, 
				    valcount, 1, // (arguments reversed from normal) 
				  left.getTypedVectorRO(dummyl), leftoffset, false,
				  right.getTypedVectorRO(dummyr), rightoffset, true,
				  operation);	  
      }        
  }
  else if (left.getRank()==0)	// scalar op on the left
  {
      escript::binaryOpVector(res.getTypedVectorRW(resdummy), 0, 
				valcount, 1,	// (arguments reversed from normal) 
			      left.getTypedVectorRO(dummyl), 0, true,
			      right.getTypedVectorRO(dummyr), 0, false,
			      operation);
      const DataTagged::DataMapType& lookup_re=res.getTagLookup();
      DataTagged::DataMapType::const_iterator i; // i->first is a tag, i->second is an offset into memory
      for (i=lookup_re.begin();i!=lookup_re.end();i++)
      {
	  DataTypes::RealVectorType::size_type resoffset=i->second;
	  DataTypes::RealVectorType::size_type leftoffset=left.getOffsetForTag(i->first);
	  DataTypes::RealVectorType::size_type rightoffset=right.getOffsetForTag(i->first);	  
	  escript::binaryOpVector(res.getTypedVectorRW(resdummy), resoffset, 
				    valcount, 1, // (arguments reversed from normal) 
				  left.getTypedVectorRO(dummyl), leftoffset, true,
				  right.getTypedVectorRO(dummyr), rightoffset, false,
				  operation);	  
      }        
  }
  else
  {
      // This will process the default value (which we know is stored in location 0)
      escript::binaryOpVector(res.getTypedVectorRW(resdummy), 0, 1, valcount, 
			      left.getTypedVectorRO(dummyl), 0, false,
			      right.getTypedVectorRO(dummyr), 0, false,
			      operation);
      const DataTagged::DataMapType& lookup_re=res.getTagLookup();
      DataTagged::DataMapType::const_iterator i; // i->first is a tag, i->second is an offset into memory
      for (i=lookup_re.begin();i!=lookup_re.end();i++)
      {
	  DataTypes::RealVectorType::size_type resoffset=res.getOffsetForTag(i->first);
	  DataTypes::RealVectorType::size_type leftoffset=left.getOffsetForTag(i->first);
	  DataTypes::RealVectorType::size_type rightoffset=right.getOffsetForTag(i->first);
	  escript::binaryOpVector(res.getTypedVectorRW(resdummy), resoffset, 1, valcount, 
				  left.getTypedVectorRO(dummyl), leftoffset, false,
				  right.getTypedVectorRO(dummyr), rightoffset, false,
				  operation);	  
      }      
  }    
}


void binaryOpDataTTT(DataTagged& result, const DataTagged& left, const DataTagged& right, 
		     escript::ES_optype operation)
{
  bool cplxresult=left.isComplex() || right.isComplex();
  if (result.isComplex()!=cplxresult)
  {
      ostringstream oss;
      oss << "Programming error: result has unexpected complexity ";
      oss << result.isComplex() << "==" << left.isComplex() << "||";
      oss << right.isComplex();
      throw DataException(oss.str());
  }
  
  if (left.isComplex())
  {
      if (right.isComplex())
      {
	  binaryOpDataReadyHelperTTT<DataTypes::cplx_t, DataTypes::cplx_t, DataTypes::cplx_t>(result, left, right, operation);
      }
      else
      {
	  binaryOpDataReadyHelperTTT<DataTypes::cplx_t, DataTypes::cplx_t, DataTypes::real_t>(result, left, right, operation);	
      }    
  }
  else	// left is real
  {
      if (right.isComplex())
      {
	  binaryOpDataReadyHelperTTT<DataTypes::cplx_t, DataTypes::real_t, DataTypes::cplx_t>(result, left, right, operation);		
      }
      else	// right is real
      {
	  binaryOpDataReadyHelperTTT<DataTypes::real_t, DataTypes::real_t, DataTypes::real_t>(result, left, right, operation);	
      }        
  }    
}

template <class ResSCALAR, class LSCALAR, class RSCALAR>
inline void binaryOpDataReadyHelperEEC(DataExpanded& res, const DataExpanded& left, const DataConstant& right, 
		     escript::ES_optype operation)
{
  ResSCALAR resdummy=0;
  LSCALAR dummyl=0;
  RSCALAR dummyr=0;
  DataTypes::RealVectorType::size_type valcount=res.getNumDPPSample()*DataTypes::noValues(res.getShape());

  if (left.hasNoSamples() || right.hasNoSamples())
  {
      return;
  }
  
  if (right.getRank()==0) 
  {	
    escript::binaryOpVectorRightScalar(res.getTypedVectorRW(resdummy), 0, res.getNumSamples(), valcount,
			      left.getTypedVectorRO(dummyl), 0,
			      &right.getTypedVectorRO(dummyr)[0], true,
			      operation,
			      false);
  }
  else if (left.getRank()==0)
  {
    escript::binaryOpVectorLeftScalar(res.getTypedVectorRW(resdummy), 0, 	// "shrink" the samples to make this work
					res.getNumSamples()*res.getNumDPPSample(),DataTypes::noValues(res.getShape()) , 
			      &left.getTypedVectorRO(dummyl)[0], 0,
			      right.getTypedVectorRO(dummyr), false,
			      operation,
			      true);
  } 
  else //(right.getRank()==left.getRank())
  {
    escript::binaryOpVector(res.getTypedVectorRW(resdummy), 0, 
			      res.getNumSamples()*res.getNumDPPSample(),DataTypes::noValues(res.getShape()) ,
			      left.getTypedVectorRO(dummyl), 0, false,
			      right.getTypedVectorRO(dummyr), 0, true,
			      operation);
  }
  
}


void binaryOpDataEEC(DataExpanded& result, const DataExpanded& left, const DataConstant& right, 
		     escript::ES_optype operation)
{
  bool cplxresult=left.isComplex() || right.isComplex();
  if (result.isComplex()!=cplxresult)
  {
      ostringstream oss;
      oss << "Programming error: result has unexpected complexity ";
      oss << result.isComplex() << "==" << left.isComplex() << "||";
      oss << right.isComplex();
      throw DataException(oss.str());
  }
  
  if (left.isComplex())
  {
      if (right.isComplex())
      {
	  binaryOpDataReadyHelperEEC<DataTypes::cplx_t, DataTypes::cplx_t, DataTypes::cplx_t>(result, left, right, operation);
      }
      else
      {
	  binaryOpDataReadyHelperEEC<DataTypes::cplx_t, DataTypes::cplx_t, DataTypes::real_t>(result, left, right, operation);	
      }    
  }
  else	// left is real
  {
      if (right.isComplex())
      {
	  binaryOpDataReadyHelperEEC<DataTypes::cplx_t, DataTypes::real_t, DataTypes::cplx_t>(result, left, right, operation);	
      }
      else	// right is real
      {
	  binaryOpDataReadyHelperEEC<DataTypes::real_t, DataTypes::real_t, DataTypes::real_t>(result, left, right, operation);	
      }        
  }    
}


template <class ResSCALAR, class LSCALAR, class RSCALAR>
inline void binaryOpDataReadyHelperEEE(DataExpanded& res, const DataExpanded& left, const DataExpanded& right, 
		     escript::ES_optype operation)
{
  ResSCALAR resdummy=0;
  LSCALAR dummyl=0;
  RSCALAR dummyr=0;
  DataTypes::RealVectorType::size_type valcount=res.getNumDPPSample()*DataTypes::noValues(res.getShape());
  
  if (left.hasNoSamples() || right.hasNoSamples())
  {
      return;
  }
  
  if (left.getRank()==right.getRank())
  {
    escript::binaryOpVector(res.getTypedVectorRW(resdummy), 0, res.getNumSamples(), valcount,
			      left.getTypedVectorRO(dummyl), 0, false,
			      right.getTypedVectorRO(dummyr), 0, false,
			      operation);
  }
  else if (right.getRank()==0) 
  {
      // in this case, we need to deal with the fact that the loops in binaryOpVector are in terms of samples
      // but for this case, we need to change arguments once per datapoint.
      // To get around this, we'll pretend that the samples are smaller (and only contain one datapoint)
      // this should hopefully lead to the same openmp thread division
    escript::binaryOpVectorRightScalar(res.getTypedVectorRW(resdummy), 0, 
					 res.getNumSamples()*res.getNumDPPSample(), DataTypes::noValues(res.getShape()), 
			      left.getTypedVectorRO(dummyl), 0, 
			      &right.getTypedVectorRO(dummyr)[0], false,
			      operation, false);
  }
  else // if (left.getRank()==0)
  {
    escript::binaryOpVectorLeftScalar(res.getTypedVectorRW(resdummy), 0, 
					res.getNumSamples()*res.getNumDPPSample(), DataTypes::noValues(res.getShape()),
			      &left.getTypedVectorRO(dummyl)[0], false,
			      right.getTypedVectorRO(dummyr), 0,
			      operation, false);
  }
}


void binaryOpDataEEE(DataExpanded& result, const DataExpanded& left, const DataExpanded& right, 
		     escript::ES_optype operation)
{
  bool cplxresult=left.isComplex() || right.isComplex();
  if (result.isComplex()!=cplxresult)
  {
      ostringstream oss;
      oss << "Programming error: result has unexpected complexity ";
      oss << result.isComplex() << "==" << left.isComplex() << "||";
      oss << right.isComplex();
      throw DataException(oss.str());
  }
  
  if (left.isComplex())
  {
      if (right.isComplex())
      {
	  binaryOpDataReadyHelperEEE<DataTypes::cplx_t, DataTypes::cplx_t, DataTypes::cplx_t>(result, left, right, operation);
      }
      else
      {
	  binaryOpDataReadyHelperEEE<DataTypes::cplx_t, DataTypes::cplx_t, DataTypes::real_t>(result, left, right, operation);	
      }    
  }
  else	// left is real
  {
      if (right.isComplex())
      {
	  binaryOpDataReadyHelperEEE<DataTypes::cplx_t, DataTypes::real_t, DataTypes::cplx_t>(result, left, right, operation);	
      }
      else	// right is real
      {
	  binaryOpDataReadyHelperEEE<DataTypes::real_t, DataTypes::real_t, DataTypes::real_t>(result, left, right, operation);	
      }        
  }    
}

}
