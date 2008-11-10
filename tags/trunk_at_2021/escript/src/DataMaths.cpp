
/*******************************************************
*
* Copyright (c) 2003-2008 by University of Queensland
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


}    // end namespace
}    // end namespace

