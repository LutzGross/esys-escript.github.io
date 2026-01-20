
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* https://github.com/LutzGross/esys-escript.github.io
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/

#include "DataVectorOps.h"
#include "DataTypes.h"

#include <sstream>

namespace
{
const int SUCCESS=0;
const int BADRANK=1;
const int NOTSQUARE=2;
const int SHAPEMISMATCH=3;
const int NOINVERSE=4;
const int NEEDLAPACK=5;
const int ERRFACTORISE=6;
const int ERRINVERT=7;
}

namespace escript
{

  void
  matMult(const DataTypes::RealVectorType& left, 
	  const DataTypes::ShapeType& leftShape,
	  DataTypes::RealVectorType::size_type leftOffset,
          const DataTypes::RealVectorType& right,
   	  const DataTypes::ShapeType& rightShape,
	  DataTypes::RealVectorType::size_type rightOffset,
          DataTypes::RealVectorType& result,
	  const DataTypes::ShapeType& resultShape)
   {
      using namespace escript::DataTypes;
      using namespace std; 

      int leftRank=getRank(leftShape);
      int rightRank=getRank(rightShape);
      int resultRank=getRank(resultShape);
      if (leftRank==0 || rightRank==0) {
         stringstream temp;
         temp << "Error - (matMult) Invalid for rank 0 objects.";
         throw DataException(temp.str());
      }

      if (leftShape[leftRank-1] != rightShape[0]) {
         stringstream temp;
         temp << "Error - (matMult) Dimension: " << leftRank 
              << ", size: " << leftShape[leftRank-1] 
              << " of LHS and dimension: 1, size: " << rightShape[0]
              << " of RHS don't match.";
         throw DataException(temp.str());
      }

      int outputRank = leftRank+rightRank-2;

      if (outputRank < 0) {
         stringstream temp;
         temp << "Error - (matMult) LHS and RHS cannot be multiplied "
              << "as they have incompatible rank.";
         throw DataException(temp.str());
      }

      if (outputRank != resultRank) {
         stringstream temp;
         temp << "Error - (matMult) Rank of result array is: " 
              << resultRank 
              << " it must be: " << outputRank;
         throw DataException(temp.str());
      }

      for (int i=0; i<(leftRank-1); i++) {
         if (leftShape[i] != resultShape[i]) {
            stringstream temp;
            temp << "Error - (matMult) Dimension: " << i 
                 << " of LHS and result array don't match.";
            throw DataException(temp.str());
         }
      }

      for (int i=1; i<rightRank; i++) {
         if (rightShape[i] != resultShape[i+leftRank-2]) {
            stringstream temp;
            temp << "Error - (matMult) Dimension: " << i
                 << ", size: " << rightShape[i]
                 << " of RHS and dimension: " << i+leftRank-1 
                 << ", size: " << resultShape[i+leftRank-1]
                 << " of result array don't match.";
            throw DataException(temp.str());
         }
      }

      switch (leftRank) {

      case 1:
         switch (rightRank) {
         case 1:
            result[0]=0;
            for (int i=0;i<leftShape[0];i++) {
               result[0]+=left[i+leftOffset]*right[i+rightOffset];
            }
            break;
         case 2:
            for (int i=0;i<resultShape[0];i++) {
               result[i]=0;
               for (int j=0;j<rightShape[0];j++) {
                  result[i]+=left[j+leftOffset]*right[getRelIndex(rightShape,j,i)+rightOffset];
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
         switch (rightRank) {
         case 1:
            result[0]=0;
            for (int i=0;i<leftShape[0];i++) {
               result[i]=0;
               for (int j=0;j<leftShape[1];j++) {
                  result[i]+=left[leftOffset+getRelIndex(leftShape,i,j)]*right[i+rightOffset];
               }
            }
	    break;
         case 2:
            for (int i=0;i<resultShape[0];i++) {
               for (int j=0;j<resultShape[1];j++) {
                  result[getRelIndex(resultShape,i,j)]=0;
                  for (int jR=0;jR<rightShape[0];jR++) {
                     result[getRelIndex(resultShape,i,j)]+=left[leftOffset+getRelIndex(leftShape,i,jR)]*right[rightOffset+getRelIndex(rightShape,jR,j)];
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
         stringstream temp; temp << "Error - (matMult) Not supported for rank: " << leftRank;
         throw DataException(temp.str());
         break;
      }

   }


   DataTypes::ShapeType
   determineResultShape(const DataTypes::ShapeType& left,
                       const DataTypes::ShapeType& right)
   {
      DataTypes::ShapeType result;
      for (int i=0; i<(DataTypes::getRank(left)-1); i++) {
         result.push_back(left[i]);
      }
      for (int i=1; i<DataTypes::getRank(right); i++) {
         result.push_back(right[i]);
      }
      return result;
   }




void matrixInverseError(int err)
{
    switch (err)
    {
    case 0: break;	// not an error
    case BADRANK: throw DataException("matrix_inverse: input and output must be rank 2.");
    case NOTSQUARE: throw DataException("matrix_inverse: matrix must be square.");
    case SHAPEMISMATCH: throw DataException("matrix_inverse: programmer error input and output must be the same shape.");
    case NOINVERSE: throw DataException("matrix_inverse: argument not invertible.");
    case NEEDLAPACK:throw DataException("matrix_inverse: matrices larger than 3x3 require lapack support."); 
    case ERRFACTORISE: throw DataException("matrix_inverse: argument not invertible (factorise stage).");
    case ERRINVERT: throw DataException("matrix_inverse: argument not invertible (inverse stage).");
    default:
	throw DataException("matrix_inverse: unknown error.");
    }
}



// Copied from the python version in util.py
int
matrix_inverse(const DataTypes::RealVectorType& in, 
	    const DataTypes::ShapeType& inShape,
            DataTypes::RealVectorType::size_type inOffset,
            DataTypes::RealVectorType& out,
	    const DataTypes::ShapeType& outShape,
            DataTypes::RealVectorType::size_type outOffset,
	    int count,
	    LapackInverseHelper& helper)
{
    using namespace DataTypes;
    using namespace std;
    int inRank=getRank(inShape);
    int outRank=getRank(outShape);
    int size=DataTypes::noValues(inShape);
    if ((inRank!=2) || (outRank!=2))
    {
	return BADRANK;		
    }
    if (inShape[0]!=inShape[1])
    {
	return NOTSQUARE; 		
    }
    if (inShape!=outShape)
    {
	return SHAPEMISMATCH;	
    }
    if (inShape[0]==1)
    {
	for (int i=0;i<count;++i)
	{
	    if (in[inOffset+i]!=0)
	    {
	    	out[outOffset+i]=1/in[inOffset+i];
	    }
	    else
	    {
		return NOINVERSE;
	    }
	}
    }
    else if (inShape[0]==2)
    {
	int step=0;
	for (int i=0;i<count;++i)
	{	
          double A11=in[inOffset+step+getRelIndex(inShape,0,0)];
          double A12=in[inOffset+step+getRelIndex(inShape,0,1)];
          double A21=in[inOffset+step+getRelIndex(inShape,1,0)];
          double A22=in[inOffset+step+getRelIndex(inShape,1,1)];
          double D = A11*A22-A12*A21;
	  if (D!=0)
	  {
          	D=1/D;
		out[outOffset+step+getRelIndex(inShape,0,0)]= A22*D;
         	out[outOffset+step+getRelIndex(inShape,1,0)]=-A21*D;
          	out[outOffset+step+getRelIndex(inShape,0,1)]=-A12*D;
          	out[outOffset+step+getRelIndex(inShape,1,1)]= A11*D;
	  }
	  else
	  {
		return NOINVERSE;
	  }
	  step+=size;
	}
    }
    else if (inShape[0]==3)
    {
	int step=0;
	for (int i=0;i<count;++i)
	{	
          double A11=in[inOffset+step+getRelIndex(inShape,0,0)];
          double A21=in[inOffset+step+getRelIndex(inShape,1,0)];
          double A31=in[inOffset+step+getRelIndex(inShape,2,0)];
          double A12=in[inOffset+step+getRelIndex(inShape,0,1)];
          double A22=in[inOffset+step+getRelIndex(inShape,1,1)];
          double A32=in[inOffset+step+getRelIndex(inShape,2,1)];
          double A13=in[inOffset+step+getRelIndex(inShape,0,2)];
          double A23=in[inOffset+step+getRelIndex(inShape,1,2)];
          double A33=in[inOffset+step+getRelIndex(inShape,2,2)];
          double D = A11*(A22*A33-A23*A32)+ A12*(A31*A23-A21*A33)+A13*(A21*A32-A31*A22);
	  if (D!=0)
	  {
		D=1/D;
          	out[outOffset+step+getRelIndex(inShape,0,0)]=(A22*A33-A23*A32)*D;
          	out[outOffset+step+getRelIndex(inShape,1,0)]=(A31*A23-A21*A33)*D;
          	out[outOffset+step+getRelIndex(inShape,2,0)]=(A21*A32-A31*A22)*D;
          	out[outOffset+step+getRelIndex(inShape,0,1)]=(A13*A32-A12*A33)*D;
          	out[outOffset+step+getRelIndex(inShape,1,1)]=(A11*A33-A31*A13)*D;
          	out[outOffset+step+getRelIndex(inShape,2,1)]=(A12*A31-A11*A32)*D;
          	out[outOffset+step+getRelIndex(inShape,0,2)]=(A12*A23-A13*A22)*D;
          	out[outOffset+step+getRelIndex(inShape,1,2)]=(A13*A21-A11*A23)*D;
          	out[outOffset+step+getRelIndex(inShape,2,2)]=(A11*A22-A12*A21)*D;
          }
	  else
	  {
		return NOINVERSE;
	  }
	  step+=size;
	}
    }
    else	// inShape[0] >3  (or negative but that can hopefully never happen)
    {
#ifndef ESYS_HAVE_LAPACK
	return NEEDLAPACK;
#else
	int step=0;
	
	
	for (int i=0;i<count;++i)
	{
		// need to make a copy since blas overwrites its input
		for (int j=0;j<size;++j)
		{
		    out[outOffset+step+j]=in[inOffset+step+j];
		}
		double* arr=&(out[outOffset+step]);
		int res=helper.invert(arr);
		if (res!=0)
		{
		    return res;
		}
		step+=size;
	}
#endif
    }
    return SUCCESS;
}


// --------------------------------------------------------

template <>
void
binaryOpVectorTagged(DataTypes::RealVectorType& res,				// where result is to be stored
	  const typename DataTypes::RealVectorType::size_type samplesToProcess,	// number of samples to be updated in the result
	  const typename DataTypes::RealVectorType::size_type DPPSample,	// number of datapoints per sample
	  const typename DataTypes::RealVectorType::size_type DPSize,		// datapoint size
		
	  const DataTypes::RealVectorType& left, 				// LHS of calculation
	  const bool leftscalar,
	  const DataTypes::RealVectorType& right, 				// RHS of the calculation
	  const bool rightscalar,		
	  const bool lefttagged,			// true if left object is the tagged one
	  const DataTagged& tagsource,			// where to get tag offsets from
	  escript::ES_optype operation)		// operation to perform	  
{
  typename DataTypes::RealVectorType::size_type lstep=leftscalar?1:DPSize;
  typename DataTypes::RealVectorType::size_type rstep=rightscalar?1:DPSize;
  typename DataTypes::RealVectorType::size_type limit=samplesToProcess*DPPSample;
  switch (operation)
  {
    case ADD:
      #pragma omp parallel for
      for (typename DataTypes::RealVectorType::size_type i=0;i<limit;++i)
      {
	  typename DataTypes::RealVectorType::size_type leftbase=(lefttagged?tagsource.getPointOffset(i/DPPSample,0):i*lstep);	// only one of these
	  typename DataTypes::RealVectorType::size_type rightbase=(lefttagged?i*rstep:tagsource.getPointOffset(i/DPPSample,0));	// will apply
	  
	  
	  for (typename DataTypes::RealVectorType::size_type j=0;j<DPSize;++j)
	  {
	      res[i*DPSize+j]=left[leftbase+j*(!leftscalar)]+right[rightbase+j*(!rightscalar)];
	  }
	
      }
      break;
    case POW:
      #pragma omp parallel for
      for (typename DataTypes::RealVectorType::size_type i=0;i<limit;++i)
      {
	  typename DataTypes::RealVectorType::size_type leftbase=(lefttagged?tagsource.getPointOffset(i/DPPSample,0):i*lstep);	// only one of these
	  typename DataTypes::RealVectorType::size_type rightbase=(lefttagged?i*rstep:tagsource.getPointOffset(i/DPPSample,0));	// will apply
	  
	  
	  for (typename DataTypes::RealVectorType::size_type j=0;j<DPSize;++j)
	  {
	      res[i*DPSize+j]=pow(left[leftbase+j*(!leftscalar)],right[rightbase+j*(!rightscalar)]);
	  }
	
      }
      break;      
    case SUB:
      #pragma omp parallel for
      for (typename DataTypes::RealVectorType::size_type i=0;i<limit;++i)
      {
	  typename DataTypes::RealVectorType::size_type leftbase=(lefttagged?tagsource.getPointOffset(i/DPPSample,0):i*lstep);	// only one of these
	  typename DataTypes::RealVectorType::size_type rightbase=(lefttagged?i*rstep:tagsource.getPointOffset(i/DPPSample,0));	// will apply
	  
	  
	  for (typename DataTypes::RealVectorType::size_type j=0;j<DPSize;++j)
	  {
	      res[i*DPSize+j]=left[leftbase+j*(!leftscalar)]-right[rightbase+j*(!rightscalar)];
	  }
	
      }
      break;      
    case MUL:
      #pragma omp parallel for
      for (typename DataTypes::RealVectorType::size_type i=0;i<limit;++i)
      {
	  typename DataTypes::RealVectorType::size_type leftbase=(lefttagged?tagsource.getPointOffset(i/DPPSample,0):i*lstep);	// only one of these
	  typename DataTypes::RealVectorType::size_type rightbase=(lefttagged?i*rstep:tagsource.getPointOffset(i/DPPSample,0));	// will apply
	  
	  
	  for (typename DataTypes::RealVectorType::size_type j=0;j<DPSize;++j)
	  {
	      res[i*DPSize+j]=left[leftbase+j*(!leftscalar)]*right[rightbase+j*(!rightscalar)];
	  }
	
      }
      break;      
    case DIV:
      #pragma omp parallel for
      for (typename DataTypes::RealVectorType::size_type i=0;i<limit;++i)
      {
	  typename DataTypes::RealVectorType::size_type leftbase=(lefttagged?tagsource.getPointOffset(i/DPPSample,0):i*lstep);	// only one of these
	  typename DataTypes::RealVectorType::size_type rightbase=(lefttagged?i*rstep:tagsource.getPointOffset(i/DPPSample,0));	// will apply
	  
	  
	  for (typename DataTypes::RealVectorType::size_type j=0;j<DPSize;++j)
	  {
	      res[i*DPSize+j]=left[leftbase+j*(!leftscalar)]/right[rightbase+j*(!rightscalar)];
	  }
	
      }
      break;      
    case LESS:
      #pragma omp parallel for
      for (typename DataTypes::RealVectorType::size_type i=0;i<limit;++i)
      {
	  typename DataTypes::RealVectorType::size_type leftbase=(lefttagged?tagsource.getPointOffset(i/DPPSample,0):i*lstep);	// only one of these
	  typename DataTypes::RealVectorType::size_type rightbase=(lefttagged?i*rstep:tagsource.getPointOffset(i/DPPSample,0));	// will apply
	  
	  
	  for (typename DataTypes::RealVectorType::size_type j=0;j<DPSize;++j)
	  {
	      res[i*DPSize+j]=left[leftbase+j*(!leftscalar)]<right[rightbase+j*(!rightscalar)];
	  }
	
      }
      break;      
    case GREATER:
      #pragma omp parallel for
      for (typename DataTypes::RealVectorType::size_type i=0;i<limit;++i)
      {
	  typename DataTypes::RealVectorType::size_type leftbase=(lefttagged?tagsource.getPointOffset(i/DPPSample,0):i*lstep);	// only one of these
	  typename DataTypes::RealVectorType::size_type rightbase=(lefttagged?i*rstep:tagsource.getPointOffset(i/DPPSample,0));	// will apply
	  
	  
	  for (typename DataTypes::RealVectorType::size_type j=0;j<DPSize;++j)
	  {
	      res[i*DPSize+j]=left[leftbase+j*(!leftscalar)]>right[rightbase+j*(!rightscalar)];
	  }
	
      }
      break;      
    case GREATER_EQUAL:
      #pragma omp parallel for
      for (typename DataTypes::RealVectorType::size_type i=0;i<limit;++i)
      {
	  typename DataTypes::RealVectorType::size_type leftbase=(lefttagged?tagsource.getPointOffset(i/DPPSample,0):i*lstep);	// only one of these
	  typename DataTypes::RealVectorType::size_type rightbase=(lefttagged?i*rstep:tagsource.getPointOffset(i/DPPSample,0));	// will apply
	  
	  
	  for (typename DataTypes::RealVectorType::size_type j=0;j<DPSize;++j)
	  {
	      res[i*DPSize+j]=left[leftbase+j*(!leftscalar)]>=right[rightbase+j*(!rightscalar)];
	  }
	
      }
      break;      
    case LESS_EQUAL:
      #pragma omp parallel for
      for (typename DataTypes::RealVectorType::size_type i=0;i<limit;++i)
      {
	  typename DataTypes::RealVectorType::size_type leftbase=(lefttagged?tagsource.getPointOffset(i/DPPSample,0):i*lstep);	// only one of these
	  typename DataTypes::RealVectorType::size_type rightbase=(lefttagged?i*rstep:tagsource.getPointOffset(i/DPPSample,0));	// will apply
	  
	  
	  for (typename DataTypes::RealVectorType::size_type j=0;j<DPSize;++j)
	  {
	      res[i*DPSize+j]=left[leftbase+j*(!leftscalar)]<=right[rightbase+j*(!rightscalar)];
	  }
	
      }
      break;      
    default:
      throw DataException("Unsupported binary operation");    
  }  
}

template <>
void
binaryOpVectorRightScalar(DataTypes::RealVectorType& res,				// where result is to be stored
	  typename DataTypes::RealVectorType::size_type resOffset,		// offset in the result vector to start storing results
	  const typename DataTypes::RealVectorType::size_type samplesToProcess,	// number of samples to be updated in the result
	  const typename DataTypes::RealVectorType::size_type sampleSize,		// number of values in each sample
	  const DataTypes::RealVectorType& left, 				// LHS of calculation
	  typename DataTypes::RealVectorType::size_type leftOffset,		// where to start reading LHS values
	  const DataTypes::real_t* right, 			// RHS of the calculation
	  const bool rightreset,			// true if RHS is providing a single sample of 1 value only
	  escript::ES_optype operation,		// operation to perform
	  bool singleleftsample)			// set to false for normal operation
{
  size_t substep=(rightreset?0:1);  
  switch (operation)
  {
    case ADD:
      #pragma omp parallel for
      for (typename DataTypes::RealVectorType::size_type i=0;i<samplesToProcess;++i)
      {
	  typename DataTypes::RealVectorType::size_type leftbase=leftOffset+(singleleftsample?0:i*sampleSize);
	  const DataTypes::real_t* rpos=right+(rightreset?0:i*substep);	
	  
	  for (typename DataTypes::RealVectorType::size_type j=0;j<sampleSize;++j)
	  {
	      res[i*sampleSize+resOffset+j]=left[leftbase+j]+*rpos;
	  }
      }
      break;
    case POW:
      #pragma omp parallel for
      for (typename DataTypes::RealVectorType::size_type i=0;i<samplesToProcess;++i)
      {
	  typename DataTypes::RealVectorType::size_type leftbase=leftOffset+(singleleftsample?0:i*sampleSize);
	  const DataTypes::real_t* rpos=right+(rightreset?0:i*substep);	
	  
	  for (typename DataTypes::RealVectorType::size_type j=0;j<sampleSize;++j)
	  {
	      res[i*sampleSize+resOffset+j]=pow(left[leftbase+j],*rpos);
	  }
      }
      break;      
    case SUB:
      #pragma omp parallel for
      for (typename DataTypes::RealVectorType::size_type i=0;i<samplesToProcess;++i)
      {
	  typename DataTypes::RealVectorType::size_type leftbase=leftOffset+(singleleftsample?0:i*sampleSize);
	  const DataTypes::real_t* rpos=right+(rightreset?0:i*substep);	
	  
	  for (typename DataTypes::RealVectorType::size_type j=0;j<sampleSize;++j)
	  {
	      res[i*sampleSize+resOffset+j]=left[leftbase+j]-*rpos;
	  }
      }
      break;      
    case MUL:
      #pragma omp parallel for
      for (typename DataTypes::RealVectorType::size_type i=0;i<samplesToProcess;++i)
      {
	  typename DataTypes::RealVectorType::size_type leftbase=leftOffset+(singleleftsample?0:i*sampleSize);
	  const DataTypes::real_t* rpos=right+(rightreset?0:i*substep);	
	  
	  for (typename DataTypes::RealVectorType::size_type j=0;j<sampleSize;++j)
	  {
	      res[i*sampleSize+resOffset+j]=left[leftbase+j] * *rpos;
	  }
      }
      break;      
    case DIV:
      #pragma omp parallel for
      for (typename DataTypes::RealVectorType::size_type i=0;i<samplesToProcess;++i)
      {
	  typename DataTypes::RealVectorType::size_type leftbase=leftOffset+(singleleftsample?0:i*sampleSize);
	  const DataTypes::real_t* rpos=right+(rightreset?0:i*substep);	
	  
	  for (typename DataTypes::RealVectorType::size_type j=0;j<sampleSize;++j)
	  {
	      res[i*sampleSize+resOffset+j]=left[leftbase+j]/ *rpos;
	  }
      }
      break;      
    case LESS:
      #pragma omp parallel for
      for (typename DataTypes::RealVectorType::size_type i=0;i<samplesToProcess;++i)
      {
	  typename DataTypes::RealVectorType::size_type leftbase=leftOffset+(singleleftsample?0:i*sampleSize);
	  const DataTypes::real_t* rpos=right+(rightreset?0:i*substep);	
	  
	  for (typename DataTypes::RealVectorType::size_type j=0;j<sampleSize;++j)
	  {
	      res[i*sampleSize+resOffset+j]=left[leftbase+j]<*rpos;
	  }
      }
      break;      
    case GREATER:
      #pragma omp parallel for
      for (typename DataTypes::RealVectorType::size_type i=0;i<samplesToProcess;++i)
      {
	  typename DataTypes::RealVectorType::size_type leftbase=leftOffset+(singleleftsample?0:i*sampleSize);
	  const DataTypes::real_t* rpos=right+(rightreset?0:i*substep);	
	  
	  for (typename DataTypes::RealVectorType::size_type j=0;j<sampleSize;++j)
	  {
	      res[i*sampleSize+resOffset+j]=left[leftbase+j]>*rpos;
	  }
      }
      break;      
    case GREATER_EQUAL:
      #pragma omp parallel for
      for (typename DataTypes::RealVectorType::size_type i=0;i<samplesToProcess;++i)
      {
	  typename DataTypes::RealVectorType::size_type leftbase=leftOffset+(singleleftsample?0:i*sampleSize);
	  const DataTypes::real_t* rpos=right+(rightreset?0:i*substep);	
	  
	  for (typename DataTypes::RealVectorType::size_type j=0;j<sampleSize;++j)
	  {
	      res[i*sampleSize+resOffset+j]=left[leftbase+j]>=*rpos;
	  }
      }
      break;      
    case LESS_EQUAL:
      #pragma omp parallel for
      for (typename DataTypes::RealVectorType::size_type i=0;i<samplesToProcess;++i)
      {
	  typename DataTypes::RealVectorType::size_type leftbase=leftOffset+(singleleftsample?0:i*sampleSize);
	  const DataTypes::real_t* rpos=right+(rightreset?0:i*substep);	
	  
	  for (typename DataTypes::RealVectorType::size_type j=0;j<sampleSize;++j)
	  {
	      res[i*sampleSize+resOffset+j]=left[leftbase+j]<=*rpos;
	  }
      }
      break;      
    default:
      throw DataException("Unsupported binary operation");    
  }  
}


template <>
void
binaryOpVectorLeftScalar(DataTypes::RealVectorType& res,				// where result is to be stored
	  typename DataTypes::RealVectorType::size_type resOffset,		// offset in the result vector to start storing results
	  const typename DataTypes::RealVectorType::size_type samplesToProcess,	// number of samples to be updated in the result
	  const typename DataTypes::RealVectorType::size_type sampleSize,		// number of values in each sample
	  const DataTypes::real_t* left, 				// LHS of calculation
          const bool leftreset,				// true if LHS is providing a single sample of 1 value only
	  const DataTypes::RealVectorType& right, 				// RHS of the calculation
	  typename DataTypes::RealVectorType::size_type rightOffset,		// where to start reading RHS values
	  escript::ES_optype operation,		// operation to perform
	  bool singlerightsample)			// right consists of a single sample
{
  size_t substep=(leftreset?0:1);
  switch (operation)
  {
    case ADD:
      #pragma omp parallel for
      for (typename DataTypes::RealVectorType::size_type i=0;i<samplesToProcess;++i)
      {
	  typename DataTypes::RealVectorType::size_type rightbase=rightOffset+(singlerightsample?0:i*sampleSize);
	  const DataTypes::real_t* lpos=left+(leftreset?0:i*substep);
	  for (typename DataTypes::RealVectorType::size_type j=0;j<sampleSize;++j)
	  {
	      res[i*sampleSize+resOffset+j]=*lpos+right[rightbase+j];
	  }	
      }
      break;
    case POW:
      #pragma omp parallel for
      for (typename DataTypes::RealVectorType::size_type i=0;i<samplesToProcess;++i)
      {
	  typename DataTypes::RealVectorType::size_type rightbase=rightOffset+(singlerightsample?0:i*sampleSize);
	  const DataTypes::real_t* lpos=left+(leftreset?0:i*substep);
	  for (typename DataTypes::RealVectorType::size_type j=0;j<sampleSize;++j)
	  {
	      res[i*sampleSize+resOffset+j]=pow(*lpos,right[rightbase+j]);
	  }	
      }
      break;      
    case SUB:
      #pragma omp parallel for
      for (typename DataTypes::RealVectorType::size_type i=0;i<samplesToProcess;++i)
      {
	  typename DataTypes::RealVectorType::size_type rightbase=rightOffset+(singlerightsample?0:i*sampleSize);
	  const DataTypes::real_t* lpos=left+(leftreset?0:i*substep);
	  for (typename DataTypes::RealVectorType::size_type j=0;j<sampleSize;++j)
	  {
	      res[i*sampleSize+resOffset+j]=*lpos-right[rightbase+j];
	  }	
      }
      break;      
    case MUL:
      #pragma omp parallel for
      for (typename DataTypes::RealVectorType::size_type i=0;i<samplesToProcess;++i)
      {
	  typename DataTypes::RealVectorType::size_type rightbase=rightOffset+(singlerightsample?0:i*sampleSize);
	  const DataTypes::real_t* lpos=left+(leftreset?0:i*substep);
	  for (typename DataTypes::RealVectorType::size_type j=0;j<sampleSize;++j)
	  {
	      res[i*sampleSize+resOffset+j]=*lpos*right[rightbase+j];
	  }	
      }
      break;      
    case DIV:
      #pragma omp parallel for
      for (typename DataTypes::RealVectorType::size_type i=0;i<samplesToProcess;++i)
      {
	  typename DataTypes::RealVectorType::size_type rightbase=rightOffset+(singlerightsample?0:i*sampleSize);
	  const DataTypes::real_t* lpos=left+(leftreset?0:i*substep);
	  for (typename DataTypes::RealVectorType::size_type j=0;j<sampleSize;++j)
	  {
	      res[i*sampleSize+resOffset+j]=*lpos/right[rightbase+j];
	  }	
      }
      break;      
    case LESS:
      #pragma omp parallel for
      for (typename DataTypes::RealVectorType::size_type i=0;i<samplesToProcess;++i)
      {
	  typename DataTypes::RealVectorType::size_type rightbase=rightOffset+(singlerightsample?0:i*sampleSize);
	  const DataTypes::real_t* lpos=left+(leftreset?0:i*substep);
	  for (typename DataTypes::RealVectorType::size_type j=0;j<sampleSize;++j)
	  {
	      res[i*sampleSize+resOffset+j]=*lpos<right[rightbase+j];
	  }	
      }
      break;      
    case GREATER:
      #pragma omp parallel for
      for (typename DataTypes::RealVectorType::size_type i=0;i<samplesToProcess;++i)
      {
	  typename DataTypes::RealVectorType::size_type rightbase=rightOffset+(singlerightsample?0:i*sampleSize);
	  const DataTypes::real_t* lpos=left+(leftreset?0:i*substep);
	  for (typename DataTypes::RealVectorType::size_type j=0;j<sampleSize;++j)
	  {
	      res[i*sampleSize+resOffset+j]=*lpos>right[rightbase+j];
	  }	
      }
      break;      
    case GREATER_EQUAL:
      #pragma omp parallel for
      for (typename DataTypes::RealVectorType::size_type i=0;i<samplesToProcess;++i)
      {
	  typename DataTypes::RealVectorType::size_type rightbase=rightOffset+(singlerightsample?0:i*sampleSize);
	  const DataTypes::real_t* lpos=left+(leftreset?0:i*substep);
	  for (typename DataTypes::RealVectorType::size_type j=0;j<sampleSize;++j)
	  {
	      res[i*sampleSize+resOffset+j]=*lpos>=right[rightbase+j];
	  }	
      }
      break;      
    case LESS_EQUAL:
      #pragma omp parallel for
      for (typename DataTypes::RealVectorType::size_type i=0;i<samplesToProcess;++i)
      {
	  typename DataTypes::RealVectorType::size_type rightbase=rightOffset+(singlerightsample?0:i*sampleSize);
	  const DataTypes::real_t* lpos=left+(leftreset?0:i*substep);
	  for (typename DataTypes::RealVectorType::size_type j=0;j<sampleSize;++j)
	  {
	      res[i*sampleSize+resOffset+j]=*lpos<=right[rightbase+j];
	  }	
      }
      break;      
    default:
      throw DataException("Unsupported binary operation");    
  }  
}

template <>
void
binaryOpVector(DataTypes::RealVectorType& res,				// where result is to be stored
	  typename DataTypes::RealVectorType::size_type resOffset,		// offset in the result vector to start storing results
	  const typename DataTypes::RealVectorType::size_type samplesToProcess,	// number of samples to be updated in the result
	  const typename DataTypes::RealVectorType::size_type sampleSize,		// number of values in each sample
	  const DataTypes::RealVectorType& left, 				// LHS of calculation
	  typename DataTypes::RealVectorType::size_type leftOffset,		// where to start reading LHS values
	  const bool leftreset,				// Is LHS only supplying a single sample instead of a bunch of them
	  const DataTypes::RealVectorType& right, 				// RHS of the calculation
	  typename DataTypes::RealVectorType::size_type rightOffset,		// where to start reading RHS values
	  const bool rightreset,			// Is RHS only supplying a single sample instead of a bunch of them
	  escript::ES_optype operation)		// operation to perform
{
  switch (operation)
  {
    case ADD:
      #pragma omp parallel for
      for (typename DataTypes::RealVectorType::size_type i=0;i<samplesToProcess;++i)
      {
	  typename DataTypes::RealVectorType::size_type leftbase=leftOffset+(leftreset?0:i*sampleSize);
	  typename DataTypes::RealVectorType::size_type rightbase=rightOffset+(rightreset?0:i*sampleSize);
	  for (typename DataTypes::RealVectorType::size_type j=0;j<sampleSize;++j)
	  {
	      res[i*sampleSize+resOffset+j]=left[leftbase+j]+right[rightbase+j];
	  }
	
      }
      break;
    case POW:
      #pragma omp parallel for
      for (typename DataTypes::RealVectorType::size_type i=0;i<samplesToProcess;++i)
      {
	  typename DataTypes::RealVectorType::size_type leftbase=leftOffset+(leftreset?0:i*sampleSize);
	  typename DataTypes::RealVectorType::size_type rightbase=rightOffset+(rightreset?0:i*sampleSize);
	  for (typename DataTypes::RealVectorType::size_type j=0;j<sampleSize;++j)
	  {
	      res[i*sampleSize+resOffset+j]=pow(left[leftbase+j],right[rightbase+j]);
	  }
	
      }
      break;      
    case SUB:
      #pragma omp parallel for
      for (typename DataTypes::RealVectorType::size_type i=0;i<samplesToProcess;++i)
      {
	  typename DataTypes::RealVectorType::size_type leftbase=leftOffset+(leftreset?0:i*sampleSize);
	  typename DataTypes::RealVectorType::size_type rightbase=rightOffset+(rightreset?0:i*sampleSize);
	  for (typename DataTypes::RealVectorType::size_type j=0;j<sampleSize;++j)
	  {
	      res[i*sampleSize+resOffset+j]=left[leftbase+j]-right[rightbase+j];
	  }
	
      }
      break;      
    case MUL:
      #pragma omp parallel for
      for (typename DataTypes::RealVectorType::size_type i=0;i<samplesToProcess;++i)
      {
	  typename DataTypes::RealVectorType::size_type leftbase=leftOffset+(leftreset?0:i*sampleSize);
	  typename DataTypes::RealVectorType::size_type rightbase=rightOffset+(rightreset?0:i*sampleSize);
	  for (typename DataTypes::RealVectorType::size_type j=0;j<sampleSize;++j)
	  {
	      res[i*sampleSize+resOffset+j]=left[leftbase+j]*right[rightbase+j];
	  }
	
      }
      break;      
    case DIV:
      #pragma omp parallel for
      for (typename DataTypes::RealVectorType::size_type i=0;i<samplesToProcess;++i)
      {
	  typename DataTypes::RealVectorType::size_type leftbase=leftOffset+(leftreset?0:i*sampleSize);
	  typename DataTypes::RealVectorType::size_type rightbase=rightOffset+(rightreset?0:i*sampleSize);
	  for (typename DataTypes::RealVectorType::size_type j=0;j<sampleSize;++j)
	  {
	      res[i*sampleSize+resOffset+j]=left[leftbase+j]/right[rightbase+j];
	  }
	
      }
      break;      
    case LESS:
      #pragma omp parallel for
      for (typename DataTypes::RealVectorType::size_type i=0;i<samplesToProcess;++i)
      {
	  typename DataTypes::RealVectorType::size_type leftbase=leftOffset+(leftreset?0:i*sampleSize);
	  typename DataTypes::RealVectorType::size_type rightbase=rightOffset+(rightreset?0:i*sampleSize);
	  for (typename DataTypes::RealVectorType::size_type j=0;j<sampleSize;++j)
	  {
	      res[i*sampleSize+resOffset+j]=left[leftbase+j]<right[rightbase+j];
	  }
	
      }
      break;      
    case GREATER:
      #pragma omp parallel for
      for (typename DataTypes::RealVectorType::size_type i=0;i<samplesToProcess;++i)
      {
	  typename DataTypes::RealVectorType::size_type leftbase=leftOffset+(leftreset?0:i*sampleSize);
	  typename DataTypes::RealVectorType::size_type rightbase=rightOffset+(rightreset?0:i*sampleSize);
	  for (typename DataTypes::RealVectorType::size_type j=0;j<sampleSize;++j)
	  {
	      res[i*sampleSize+resOffset+j]=left[leftbase+j]>right[rightbase+j];
	  }
	
      }
      break;      
    case GREATER_EQUAL:
      #pragma omp parallel for
      for (typename DataTypes::RealVectorType::size_type i=0;i<samplesToProcess;++i)
      {
	  typename DataTypes::RealVectorType::size_type leftbase=leftOffset+(leftreset?0:i*sampleSize);
	  typename DataTypes::RealVectorType::size_type rightbase=rightOffset+(rightreset?0:i*sampleSize);
	  for (typename DataTypes::RealVectorType::size_type j=0;j<sampleSize;++j)
	  {
	      res[i*sampleSize+resOffset+j]=left[leftbase+j]>=right[rightbase+j];
	  }
	
      }
      break;      
    case LESS_EQUAL:
      #pragma omp parallel for
      for (typename DataTypes::RealVectorType::size_type i=0;i<samplesToProcess;++i)
      {
	  typename DataTypes::RealVectorType::size_type leftbase=leftOffset+(leftreset?0:i*sampleSize);
	  typename DataTypes::RealVectorType::size_type rightbase=rightOffset+(rightreset?0:i*sampleSize);
	  for (typename DataTypes::RealVectorType::size_type j=0;j<sampleSize;++j)
	  {
	      res[i*sampleSize+resOffset+j]=left[leftbase+j]<=right[rightbase+j];
	  }
	
      }
      break;      
    default:
      throw DataException("Unsupported binary operation");    
  }  
}

  /**
     \brief
     computes an hermitian matrix from your square matrix A: (A + adjoint(A)) / 2

     \param in - vector containing the matrix A
     \param inShape - shape of the matrix A
     \param inOffset - the beginning of A within the vector in
     \param ev - vector to store the output matrix
     \param evShape - expected shape of the output matrix
     \param evOffset - starting location for storing ev in vector ev
  */
  void
   hermitian(const DataTypes::CplxVectorType& in, 
	    const DataTypes::ShapeType& inShape,
            DataTypes::CplxVectorType::size_type inOffset,
            DataTypes::CplxVectorType& ev, 
	    const DataTypes::ShapeType& evShape,
            DataTypes::CplxVectorType::size_type evOffset)
  {
   if (DataTypes::getRank(inShape) == 2) {
     int i0, i1;
     int s0=inShape[0];
     int s1=inShape[1];
     for (i0=0; i0<s0; i0++) {
       for (i1=0; i1<s1; i1++) {
         ev[evOffset+DataTypes::getRelIndex(evShape,i0,i1)] = (in[inOffset+DataTypes::getRelIndex(inShape,i0,i1)] + std::conj(in[inOffset+DataTypes::getRelIndex(inShape,i1,i0)])) / 2.0;
       }
     }
    }
    else if (DataTypes::getRank(inShape) == 4) {
      int i0, i1, i2, i3;
      int s0=inShape[0];
      int s1=inShape[1];
      int s2=inShape[2];
      int s3=inShape[3];
      for (i0=0; i0<s0; i0++) {
        for (i1=0; i1<s1; i1++) {
          for (i2=0; i2<s2; i2++) {
            for (i3=0; i3<s3; i3++) {
              ev[evOffset+DataTypes::getRelIndex(evShape,i0,i1,i2,i3)] = (in[inOffset+DataTypes::getRelIndex(inShape,i0,i1,i2,i3)] + std::conj(in[inOffset+DataTypes::getRelIndex(inShape,i2,i3,i0,i1)])) / 2.0;
            }
          }
        }
      }
    }
   }

  /**
     \brief
     computes a antihermitian matrix from your square matrix A: (A - adjoint(A)) / 2

     \param in - vector containing the matrix A
     \param inShape - shape of the matrix A
     \param inOffset - the beginning of A within the vector in
     \param ev - vector to store the output matrix
     \param evShape - expected shape of the output matrix
     \param evOffset - starting location for storing ev in vector ev
  */
   void
   antihermitian(const DataTypes::CplxVectorType& in, 
 	    const DataTypes::ShapeType& inShape,
             typename DataTypes::CplxVectorType::size_type inOffset,
             DataTypes::CplxVectorType& ev, 
 	    const DataTypes::ShapeType& evShape,
             typename DataTypes::CplxVectorType::size_type evOffset)  
   {
    if (DataTypes::getRank(inShape) == 2) {
      int i0, i1;
      int s0=inShape[0];
      int s1=inShape[1];
      for (i0=0; i0<s0; i0++) {
        for (i1=0; i1<s1; i1++) {
          ev[evOffset+DataTypes::getRelIndex(evShape,i0,i1)] = (in[inOffset+DataTypes::getRelIndex(inShape,i0,i1)] - std::conj(in[inOffset+DataTypes::getRelIndex(inShape,i1,i0)])) / 2.0;
        }
      }
    }
   else if (DataTypes::getRank(inShape) == 4) {
     int i0, i1, i2, i3;
     int s0=inShape[0];
     int s1=inShape[1];
     int s2=inShape[2];
     int s3=inShape[3];
     for (i0=0; i0<s0; i0++) {
       for (i1=0; i1<s1; i1++) {
         for (i2=0; i2<s2; i2++) {
           for (i3=0; i3<s3; i3++) {
             ev[evOffset+DataTypes::getRelIndex(evShape,i0,i1,i2,i3)] = (in[inOffset+DataTypes::getRelIndex(inShape,i0,i1,i2,i3)] - std::conj(in[inOffset+DataTypes::getRelIndex(inShape,i2,i3,i0,i1)])) / 2.0;
           }
         }
       }
     }
   }
  }


}    // end namespace

