
/*******************************************************
*
* Copyright (c) 2003-2011 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


#include "DataTypes.h"
#include "DataMaths.h"
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
namespace DataMaths
{

  void
  matMult(const DataTypes::ValueType& left, 
	  const DataTypes::ShapeType& leftShape,
	  DataTypes::ValueType::size_type leftOffset,
          const DataTypes::ValueType& right,
   	  const DataTypes::ShapeType& rightShape,
	  DataTypes::ValueType::size_type rightOffset,
          DataTypes::ValueType& result,
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
              << "as they have incompatable rank.";
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
matrix_inverse(const DataTypes::ValueType& in, 
	    const DataTypes::ShapeType& inShape,
            DataTypes::ValueType::size_type inOffset,
            DataTypes::ValueType& out,
	    const DataTypes::ShapeType& outShape,
            DataTypes::ValueType::size_type outOffset,
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
#ifndef USE_LAPACK
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

}    // end namespace
}    // end namespace

