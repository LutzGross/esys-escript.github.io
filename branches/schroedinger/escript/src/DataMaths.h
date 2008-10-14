
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


#ifndef escript_DataMaths_20080822_H
#define escript_DataMaths_20080822_H
#include "DataAbstract.h"
#include "DataException.h"
#include "LocalOps.h"

namespace escript
{
namespace DataMaths
{

/**
\namespace escript::DataMaths
\brief Contains maths operations performed on data vectors.

In order to properly identify the datapoints, in most cases, the vector, shape and offset of the point must all be supplied.
Note that vector in this context refers to a data vector storing datapoints not a mathematical vector. (However, datapoints within the data vector could represent scalars, vectors, matricies, ...).
*/

//  /**
//     \brief
//     Perform the unary operation on the data point specified by the view's
//     default offset. Applies the specified operation to each value in the data
//     point. 
//
//     Called by escript::unaryOp.
//
//     \param operation - Input -
//                  Operation to apply. Operation must be a pointer to a function.
//  */
//   template <class UnaryFunction>
//   void
//   unaryOp(DataAbstract& data, UnaryFunction operation);

// I'm going to try not to have the above function. It relies on the value of the offset already 
// being known.   I don't want that, offsets need to be explicit.



  /**
     \brief
     Perform the unary operation on the data point specified by the given
     offset. Applies the specified operation to each value in the data
     point. Operation must be a pointer to a function.

     Called by escript::unaryOp.

     \param data - vector containing the datapoint
     \param shape - shape of the point
     \param offset - offset of the point within data
     \param operation - Input -
                  Operation to apply. Must be a pointer to a function.
  */
  template <class UnaryFunction>
  void
  unaryOp(DataTypes::ValueType& data, const DataTypes::ShapeType& shape,
          DataTypes::ValueType::size_type offset,
          UnaryFunction operation);

//  /**
//     \brief
//     Perform the binary operation on the data points specified by the default
//     offsets in this view and in view "right". Applies the specified operation
//     to corresponding values in both data points. Operation must be a pointer
//     to a function.
//
//     Called by escript::binaryOp.
//
//     \param right - Input -
//                  View to act as RHS in given binary operation.
//     \param operation - Input -
//                  Operation to apply. Must be a pointer to a function.
//  */
//   template <class BinaryFunction>
//   void
//   binaryOp(DataAbstract& left, const DataAbstract& right,
//            BinaryFunction operation);

// trying to avoid having this one as well. Again, implicit offset

  /**
     \brief
     Perform the binary operation on the data points specified by the given
     offsets in the "left" and "right" vectors. Applies the specified operation
     to corresponding values in both data points. Operation must be a pointer
     to a function.

     Called by escript::binaryOp.
     \param left,right - vectors containing the datapoints
     \param leftShape,rightShape - shapes of datapoints in the vectors
     \param leftOffset,rightOffset - beginnings of datapoints in the vectors
     \param operation - Input -
                  Operation to apply. Must be a pointer to a function.
  */
  template <class BinaryFunction>
  void
  binaryOp(DataTypes::ValueType& left, 
	   const DataTypes::ShapeType& leftShape, 
           DataTypes::ValueType::size_type leftOffset,
           const DataTypes::ValueType& right, 
           const DataTypes::ShapeType& rightShape,
           DataTypes::ValueType::size_type rightOffset,
           BinaryFunction operation);

  /**
     \brief
     Perform the binary operation on the data point specified by the given
     offset in the vector using the scalar value "right". Applies the specified
     operation to values in the data point. Operation must be a pointer
     to a function.

     Called by escript::binaryOp.

     \param left - vector containing the datapoints
     \param shape - shape of datapoint in the vector
     \param offset - beginning of datapoint in the vector
     \param right - scalar value for the right hand side of the operation
     \param operation - Input -
                  Operation to apply. Must be a pointer to a function.
  */
  template <class BinaryFunction>
  void
  binaryOp(DataTypes::ValueType& left, 
           const DataTypes::ShapeType& shape,
 	   DataTypes::ValueType::size_type offset,
           double right,
           BinaryFunction operation);

//  /**
//     \brief
//     Perform the given data point reduction operation on the data point
//     specified by the default offset into the view. Reduces all elements of
//     the data point using the given operation, returning the result as a 
//     scalar. Operation must be a pointer to a function.
//
//     Called by escript::algorithm.
//
//     \param operation - Input -
//                  Operation to apply. Must be a pointer to a function.
//  */
//   template <class BinaryFunction>
//   double
//   reductionOp(const DataAbstract& left, BinaryFunction operation,
//               double initial_value);

// implicit offset

  /**
     \brief
     Perform the given data point reduction operation on the data point
     specified by the given offset into the view. Reduces all elements of
     the data point using the given operation, returning the result as a 
     scalar. Operation must be a pointer to a function.

     Called by escript::algorithm.

     \param left - vector containing the datapoint
     \param shape - shape of datapoints in the vector
     \param offset - beginning of datapoint in the vector
     \param operation - Input -
                  Operation to apply. Must be a pointer to a function.
     \param initial_value 
  */
  template <class BinaryFunction>
  double
  reductionOp(const DataTypes::ValueType& left, 
	      const DataTypes::ShapeType& shape,
 	      DataTypes::ValueType::size_type offset,
              BinaryFunction operation,
              double initial_value);

 /**
     \brief
     Perform a matrix multiply of the given views.

     NB: Only multiplies together the two given datapoints,
     would need to call this over all data-points to multiply the entire
     Data objects involved.

     \param left,right - vectors containing the datapoints
     \param leftShape,rightShape - shapes of datapoints in the vectors
     \param leftOffset,rightOffset - beginnings of datapoints in the vectors
     \param result - Vector to store the resulting datapoint in
     \param resultShape - expected shape of the resulting datapoint
  */
  ESCRIPT_DLL_API
  void
  matMult(const DataTypes::ValueType& left, 
	  const DataTypes::ShapeType& leftShape,
	  DataTypes::ValueType::size_type leftOffset,
          const DataTypes::ValueType& right,
	  const DataTypes::ShapeType& rightShape,
	  DataTypes::ValueType::size_type rightOffset,
          DataTypes::ValueType& result,
	  const DataTypes::ShapeType& resultShape);
// Hmmmm why is there no offset for the result??




  /**
     \brief
     Determine the shape of the result array for a matrix multiplication
     of the given views.

     \param left,right - shapes of the left and right matricies
     \return the shape of the matrix which would result from multiplying left and right
  */
  ESCRIPT_DLL_API
  DataTypes::ShapeType
  determineResultShape(const DataTypes::ShapeType& left,
                       const DataTypes::ShapeType& right);

  /**
     \brief
     computes a symmetric matrix from your square matrix A: (A + transpose(A)) / 2

     \param in - vector containing the matrix A
     \param inShape - shape of the matrix A
     \param inOffset - the beginning of A within the vector in
     \param ev - vector to store the output matrix
     \param evShape - expected shape of the output matrix
     \param evOffset - starting location for storing ev in vector ev
  */
  ESCRIPT_DLL_API
  inline
  void
  symmetric(const DataTypes::ValueType& in, 
	    const DataTypes::ShapeType& inShape,
            DataTypes::ValueType::size_type inOffset,
            DataTypes::ValueType& ev, 
	    const DataTypes::ShapeType& evShape,
            DataTypes::ValueType::size_type evOffset)
  {
   if (DataTypes::getRank(inShape) == 2) {
     int i0, i1;
     int s0=inShape[0];
     int s1=inShape[1];
     for (i0=0; i0<s0; i0++) {
       for (i1=0; i1<s1; i1++) {
         ev[evOffset+DataTypes::getRelIndex(evShape,i0,i1)] = (in[inOffset+DataTypes::getRelIndex(inShape,i0,i1)] + in[inOffset+DataTypes::getRelIndex(inShape,i1,i0)]) / 2.0;
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
             ev[evOffset+DataTypes::getRelIndex(evShape,i0,i1,i2,i3)] = (in[inOffset+DataTypes::getRelIndex(inShape,i0,i1,i2,i3)] + in[inOffset+DataTypes::getRelIndex(inShape,i2,i3,i0,i1)]) / 2.0;
           }
         }
       }
     }
   }
  }

  /**
     \brief
     computes a nonsymmetric matrix from your square matrix A: (A - transpose(A)) / 2

     \param in - vector containing the matrix A
     \param inShape - shape of the matrix A
     \param inOffset - the beginning of A within the vector in
     \param ev - vector to store the output matrix
     \param evShape - expected shape of the output matrix
     \param evOffset - starting location for storing ev in vector ev
  */
  ESCRIPT_DLL_API
  inline
  void
  nonsymmetric(const DataTypes::ValueType& in, 
	       const DataTypes::ShapeType& inShape,
               DataTypes::ValueType::size_type inOffset,
               DataTypes::ValueType& ev, 
	       const DataTypes::ShapeType& evShape,
               DataTypes::ValueType::size_type evOffset)
  {
   if (DataTypes::getRank(inShape) == 2) {
     int i0, i1;
     int s0=inShape[0];
     int s1=inShape[1];
     for (i0=0; i0<s0; i0++) {
       for (i1=0; i1<s1; i1++) {
         ev[evOffset+DataTypes::getRelIndex(evShape,i0,i1)] = (in[inOffset+DataTypes::getRelIndex(inShape,i0,i1)] - in[inOffset+DataTypes::getRelIndex(inShape,i1,i0)]) / 2.0;
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
             ev[evOffset+DataTypes::getRelIndex(evShape,i0,i1,i2,i3)] = (in[inOffset+DataTypes::getRelIndex(inShape,i0,i1,i2,i3)] - in[inOffset+DataTypes::getRelIndex(inShape,i2,i3,i0,i1)]) / 2.0;
           }
         }
       }
     }
   }
  }

  /**
     \brief
     computes the trace of a matrix

     \param in - vector containing the input matrix
     \param inShape - shape of the input matrix
     \param inOffset - the beginning of the input matrix within the vector "in"
     \param ev - vector to store the output matrix
     \param evShape - expected shape of the output matrix
     \param evOffset - starting location for storing the output matrix in vector ev
     \param axis_offset
  */
  inline
  void
  trace(const DataTypes::ValueType& in, 
	    const DataTypes::ShapeType& inShape,
            DataTypes::ValueType::size_type inOffset,
            DataTypes::ValueType& ev,
	    const DataTypes::ShapeType& evShape,
            DataTypes::ValueType::size_type evOffset,
	    int axis_offset)
  {
   if (DataTypes::getRank(inShape) == 2) {
     int s0=inShape[0]; // Python wrapper limits to square matrix
     int i;
     for (i=0; i<s0; i++) {
       ev[evOffset/*+DataTypes::getRelIndex(evShape)*/] += in[inOffset+DataTypes::getRelIndex(inShape,i,i)];
     }
   }
   else if (DataTypes::getRank(inShape) == 3) {
     if (axis_offset==0) {
       int s0=inShape[0];
       int s2=inShape[2];
       int i0, i2;
       for (i0=0; i0<s0; i0++) {
         for (i2=0; i2<s2; i2++) {
           ev[evOffset+DataTypes::getRelIndex(evShape,i2)] += in[inOffset+DataTypes::getRelIndex(inShape,i0,i0,i2)];
         }
       }
     }
     else if (axis_offset==1) {
       int s0=inShape[0];
       int s1=inShape[1];
       int i0, i1;
       for (i0=0; i0<s0; i0++) {
         for (i1=0; i1<s1; i1++) {
           ev[evOffset+DataTypes::getRelIndex(evShape,i0)] += in[inOffset+DataTypes::getRelIndex(inShape,i0,i1,i1)];
         }
       }
     }
   }
   else if (DataTypes::getRank(inShape) == 4) {
     if (axis_offset==0) {
       int s0=inShape[0];
       int s2=inShape[2];
       int s3=inShape[3];
       int i0, i2, i3;
       for (i0=0; i0<s0; i0++) {
         for (i2=0; i2<s2; i2++) {
           for (i3=0; i3<s3; i3++) {
             ev[evOffset+DataTypes::getRelIndex(evShape,i2,i3)] += in[inOffset+DataTypes::getRelIndex(inShape,i0,i0,i2,i3)];
           }
         }
       }
     }
     else if (axis_offset==1) {
       int s0=inShape[0];
       int s1=inShape[1];
       int s3=inShape[3];
       int i0, i1, i3;
       for (i0=0; i0<s0; i0++) {
         for (i1=0; i1<s1; i1++) {
           for (i3=0; i3<s3; i3++) {
             ev[evOffset+DataTypes::getRelIndex(evShape,i0,i3)] += in[inOffset+DataTypes::getRelIndex(inShape,i0,i1,i1,i3)];
           }
         }
       }
     }
     else if (axis_offset==2) {
       int s0=inShape[0];
       int s1=inShape[1];
       int s2=inShape[2];
       int i0, i1, i2;
       for (i0=0; i0<s0; i0++) {
         for (i1=0; i1<s1; i1++) {
           for (i2=0; i2<s2; i2++) {
             ev[evOffset+DataTypes::getRelIndex(evShape,i0,i1)] += in[inOffset+DataTypes::getRelIndex(inShape,i0,i1,i2,i2)];
           }
         }
       }
     }
   }
  }

  /**
     \brief
     Transpose each data point of this Data object around the given axis.

     \param in - vector containing the input matrix
     \param inShape - shape of the input matrix
     \param inOffset - the beginning of the input matrix within the vector "in"
     \param ev - vector to store the output matrix
     \param evShape - expected shape of the output matrix
     \param evOffset - starting location for storing the output matrix in vector ev
     \param axis_offset
  */
  ESCRIPT_DLL_API
  inline
  void
  transpose(const DataTypes::ValueType& in, 
	    const DataTypes::ShapeType& inShape,
            DataTypes::ValueType::size_type inOffset,
            DataTypes::ValueType& ev,
            const DataTypes::ShapeType& evShape,
            DataTypes::ValueType::size_type evOffset,
	    int axis_offset)
  {
   int inRank=DataTypes::getRank(inShape);
   if ( inRank== 4) {
     int s0=evShape[0];
     int s1=evShape[1];
     int s2=evShape[2];
     int s3=evShape[3];
     int i0, i1, i2, i3;
     if (axis_offset==1) {
       for (i0=0; i0<s0; i0++) {
         for (i1=0; i1<s1; i1++) {
           for (i2=0; i2<s2; i2++) {
             for (i3=0; i3<s3; i3++) {
               ev[evOffset+DataTypes::getRelIndex(evShape,i0,i1,i2,i3)] = in[inOffset+DataTypes::getRelIndex(inShape,i3,i0,i1,i2)];
             }
           }
         }
       }
     }
     else if (axis_offset==2) {
       for (i0=0; i0<s0; i0++) {
         for (i1=0; i1<s1; i1++) {
           for (i2=0; i2<s2; i2++) {
             for (i3=0; i3<s3; i3++) {
               ev[evOffset+DataTypes::getRelIndex(evShape,i0,i1,i2,i3)] = in[inOffset+DataTypes::getRelIndex(inShape,i2,i3,i0,i1)];
             }
           }
         }
       }
     }
     else if (axis_offset==3) {
       for (i0=0; i0<s0; i0++) {
         for (i1=0; i1<s1; i1++) {
           for (i2=0; i2<s2; i2++) {
             for (i3=0; i3<s3; i3++) {
               ev[evOffset+DataTypes::getRelIndex(evShape,i0,i1,i2,i3)] = in[inOffset+DataTypes::getRelIndex(inShape,i1,i2,i3,i0)];
             }
           }
         }
       }
     }
     else {
       for (i0=0; i0<s0; i0++) {
         for (i1=0; i1<s1; i1++) {
           for (i2=0; i2<s2; i2++) {
             for (i3=0; i3<s3; i3++) {
               ev[evOffset+DataTypes::getRelIndex(evShape,i0,i1,i2,i3)] = in[inOffset+DataTypes::getRelIndex(inShape,i0,i1,i2,i3)];
             }
           }
         }
       }
     }
   }
   else if (inRank == 3) {
     int s0=evShape[0];
     int s1=evShape[1];
     int s2=evShape[2];
     int i0, i1, i2;
     if (axis_offset==1) {
       for (i0=0; i0<s0; i0++) {
         for (i1=0; i1<s1; i1++) {
           for (i2=0; i2<s2; i2++) {
             ev[evOffset+DataTypes::getRelIndex(evShape,i0,i1,i2)] = in[inOffset+DataTypes::getRelIndex(inShape,i2,i0,i1)];
           }
         }
       }
     }
     else if (axis_offset==2) {
       for (i0=0; i0<s0; i0++) {
         for (i1=0; i1<s1; i1++) {
           for (i2=0; i2<s2; i2++) {
             ev[evOffset+DataTypes::getRelIndex(evShape,i0,i1,i2)] = in[inOffset+DataTypes::getRelIndex(inShape,i1,i2,i0)];
           }
         }
       }
     }
     else {
       // Copy the matrix unchanged
       for (i0=0; i0<s0; i0++) {
         for (i1=0; i1<s1; i1++) {
           for (i2=0; i2<s2; i2++) {
             ev[evOffset+DataTypes::getRelIndex(evShape,i0,i1,i2)] = in[inOffset+DataTypes::getRelIndex(inShape,i0,i1,i2)];
           }
         }
       }
     }
   }
   else if (inRank == 2) {
     int s0=evShape[0];
     int s1=evShape[1];
     int i0, i1;
     if (axis_offset==1) {
       for (i0=0; i0<s0; i0++) {
         for (i1=0; i1<s1; i1++) {
           ev[evOffset+DataTypes::getRelIndex(evShape,i0,i1)] = in[inOffset+DataTypes::getRelIndex(inShape,i1,i0)];
         }
       }
     }
     else {
       for (i0=0; i0<s0; i0++) {
         for (i1=0; i1<s1; i1++) {
           ev[evOffset+DataTypes::getRelIndex(evShape,i0,i1)] = in[inOffset+DataTypes::getRelIndex(inShape,i0,i1)];
         }
       }
     }
   }
   else if (inRank == 1) {
     int s0=evShape[0];
     int i0;
     for (i0=0; i0<s0; i0++) {
       ev[evOffset+DataTypes::getRelIndex(evShape,i0)] = in[inOffset+DataTypes::getRelIndex(inShape,i0)];
     }
   }
   else if (inRank == 0) {
     ev[evOffset/*+DataTypes::getRelIndex(evShape,)*/] = in[inOffset/*+DataTypes::getRelIndex(inShape,)*/];
   }
   else {
      throw DataException("Error - DataArrayView::transpose can only be calculated for rank 0, 1, 2, 3 or 4 objects.");
   }
  }

  /**
     \brief
     swaps the components axis0 and axis1.

     \param in - vector containing the input matrix
     \param inShape - shape of the input matrix
     \param inOffset - the beginning of the input matrix within the vector "in"
     \param ev - vector to store the output matrix
     \param evShape - expected shape of the output matrix
     \param evOffset - starting location for storing the output matrix in vector ev
     \param axis0 - axis index
     \param axis1 - axis index
  */
  ESCRIPT_DLL_API
  inline
  void
  swapaxes(const DataTypes::ValueType& in, 
	   const DataTypes::ShapeType& inShape,
           DataTypes::ValueType::size_type inOffset,
           DataTypes::ValueType& ev,
	   const DataTypes::ShapeType& evShape,
           DataTypes::ValueType::size_type evOffset,
           int axis0, 
	   int axis1)
  {
     int inRank=DataTypes::getRank(inShape);
     if (inRank == 4) {
     int s0=evShape[0];
     int s1=evShape[1];
     int s2=evShape[2];
     int s3=evShape[3];
     int i0, i1, i2, i3;
     if (axis0==0) {
        if (axis1==1) {
            for (i0=0; i0<s0; i0++) {
              for (i1=0; i1<s1; i1++) {
                for (i2=0; i2<s2; i2++) {
                  for (i3=0; i3<s3; i3++) {
                    ev[evOffset+DataTypes::getRelIndex(evShape,i0,i1,i2,i3)] = in[inOffset+DataTypes::getRelIndex(inShape,i1,i0,i2,i3)];
                  }
                }
              }
            }
        } else if (axis1==2) {
            for (i0=0; i0<s0; i0++) {
              for (i1=0; i1<s1; i1++) {
                for (i2=0; i2<s2; i2++) {
                  for (i3=0; i3<s3; i3++) {
                    ev[evOffset+DataTypes::getRelIndex(evShape,i0,i1,i2,i3)] = in[inOffset+DataTypes::getRelIndex(inShape,i2,i1,i0,i3)];
                  }
                }
              }
            }

        } else if (axis1==3) {
            for (i0=0; i0<s0; i0++) {
              for (i1=0; i1<s1; i1++) {
                for (i2=0; i2<s2; i2++) {
                  for (i3=0; i3<s3; i3++) {
                    ev[evOffset+DataTypes::getRelIndex(evShape,i0,i1,i2,i3)] = in[inOffset+DataTypes::getRelIndex(inShape,i3,i1,i2,i0)];
                  }
                }
              }
            }
        }
     } else if (axis0==1) {
        if (axis1==2) {
            for (i0=0; i0<s0; i0++) {
              for (i1=0; i1<s1; i1++) {
                for (i2=0; i2<s2; i2++) {
                  for (i3=0; i3<s3; i3++) {
                    ev[evOffset+DataTypes::getRelIndex(evShape,i0,i1,i2,i3)] = in[inOffset+DataTypes::getRelIndex(inShape,i0,i2,i1,i3)];
                  }
                }
              }
            }
        } else if (axis1==3) {
            for (i0=0; i0<s0; i0++) {
              for (i1=0; i1<s1; i1++) {
                for (i2=0; i2<s2; i2++) {
                  for (i3=0; i3<s3; i3++) {
                    ev[evOffset+DataTypes::getRelIndex(evShape,i0,i1,i2,i3)] = in[inOffset+DataTypes::getRelIndex(inShape,i0,i3,i2,i1)];
                  }
                }
              }
            }
        }
     } else if (axis0==2) {
        if (axis1==3) {
            for (i0=0; i0<s0; i0++) {
              for (i1=0; i1<s1; i1++) {
                for (i2=0; i2<s2; i2++) {
                  for (i3=0; i3<s3; i3++) {
                    ev[evOffset+DataTypes::getRelIndex(evShape,i0,i1,i2,i3)] = in[inOffset+DataTypes::getRelIndex(inShape,i0,i1,i3,i2)];
                  }
                }
              }
            }
        }
     }

   } else if ( inRank == 3) {
     int s0=evShape[0];
     int s1=evShape[1];
     int s2=evShape[2];
     int i0, i1, i2;
     if (axis0==0) {
        if (axis1==1) {
           for (i0=0; i0<s0; i0++) {
             for (i1=0; i1<s1; i1++) {
               for (i2=0; i2<s2; i2++) {
                 ev[evOffset+DataTypes::getRelIndex(evShape,i0,i1,i2)] = in[inOffset+DataTypes::getRelIndex(inShape,i1,i0,i2)];
               }
             }
           }
        } else if (axis1==2) {
           for (i0=0; i0<s0; i0++) {
             for (i1=0; i1<s1; i1++) {
               for (i2=0; i2<s2; i2++) {
                 ev[evOffset+DataTypes::getRelIndex(evShape,i0,i1,i2)] = in[inOffset+DataTypes::getRelIndex(inShape,i2,i1,i0)];
               }
             }
           }
       }
     } else if (axis0==1) {
        if (axis1==2) {
           for (i0=0; i0<s0; i0++) {
             for (i1=0; i1<s1; i1++) {
               for (i2=0; i2<s2; i2++) {
                 ev[evOffset+DataTypes::getRelIndex(evShape,i0,i1,i2)] = in[inOffset+DataTypes::getRelIndex(inShape,i0,i2,i1)];
               }
             }
           }
        }
     }
   } else if ( inRank == 2) {
     int s0=evShape[0];
     int s1=evShape[1];
     int i0, i1;
     if (axis0==0) {
        if (axis1==1) {
           for (i0=0; i0<s0; i0++) {
             for (i1=0; i1<s1; i1++) {
                 ev[evOffset+DataTypes::getRelIndex(evShape,i0,i1)] = in[inOffset+DataTypes::getRelIndex(inShape,i1,i0)];
             }
           }
        }
    }
  } else {
      throw DataException("Error - DataArrayView::swapaxes can only be calculated for rank 2, 3 or 4 objects.");
  }
 }

  /**
     \brief
     solves a local eigenvalue problem 

     \param in - vector containing the input matrix
     \param inShape - shape of the input matrix
     \param inOffset - the beginning of the input matrix within the vector "in"
     \param ev - vector to store the eigenvalues
     \param evShape - expected shape of the eigenvalues
     \param evOffset - starting location for storing the eigenvalues in vector ev
  */
  ESCRIPT_DLL_API
  inline
  void
  eigenvalues(DataTypes::ValueType& in, 
	      const DataTypes::ShapeType& inShape,
              DataTypes::ValueType::size_type inOffset,
              DataTypes::ValueType& ev,
	      const DataTypes::ShapeType& evShape,
              DataTypes::ValueType::size_type evOffset)
  {
   double in00,in10,in20,in01,in11,in21,in02,in12,in22;
   double ev0,ev1,ev2;
   int s=inShape[0];
   if (s==1) {
      in00=in[inOffset+DataTypes::getRelIndex(inShape,0,0)];
      eigenvalues1(in00,&ev0);
      ev[evOffset+DataTypes::getRelIndex(evShape,0)]=ev0;

   } else  if (s==2) {
      in00=in[inOffset+DataTypes::getRelIndex(inShape,0,0)];
      in10=in[inOffset+DataTypes::getRelIndex(inShape,1,0)];
      in01=in[inOffset+DataTypes::getRelIndex(inShape,0,1)];
      in11=in[inOffset+DataTypes::getRelIndex(inShape,1,1)];
      eigenvalues2(in00,(in01+in10)/2.,in11,&ev0,&ev1);
      ev[evOffset+DataTypes::getRelIndex(evShape,0)]=ev0;
      ev[evOffset+DataTypes::getRelIndex(evShape,1)]=ev1;

   } else  if (s==3) {
      in00=in[inOffset+DataTypes::getRelIndex(inShape,0,0)];
      in10=in[inOffset+DataTypes::getRelIndex(inShape,1,0)];
      in20=in[inOffset+DataTypes::getRelIndex(inShape,2,0)];
      in01=in[inOffset+DataTypes::getRelIndex(inShape,0,1)];
      in11=in[inOffset+DataTypes::getRelIndex(inShape,1,1)];
      in21=in[inOffset+DataTypes::getRelIndex(inShape,2,1)];
      in02=in[inOffset+DataTypes::getRelIndex(inShape,0,2)];
      in12=in[inOffset+DataTypes::getRelIndex(inShape,1,2)];
      in22=in[inOffset+DataTypes::getRelIndex(inShape,2,2)];
      eigenvalues3(in00,(in01+in10)/2.,(in02+in20)/2.,in11,(in21+in12)/2.,in22,
                 &ev0,&ev1,&ev2);
      ev[evOffset+DataTypes::getRelIndex(evShape,0)]=ev0;
      ev[evOffset+DataTypes::getRelIndex(evShape,1)]=ev1;
      ev[evOffset+DataTypes::getRelIndex(evShape,2)]=ev2;

   }
  }

  /**
     \brief
     solves a local eigenvalue problem 

     \param in - vector containing the input matrix
     \param inShape - shape of the input matrix
     \param inOffset - the beginning of the input matrix within the vector "in"
     \param ev - vector to store the eigenvalues
     \param evShape - expected shape of the eigenvalues
     \param evOffset - starting location for storing the eigenvalues in ev
     \param V - vector to store the eigenvectors
     \param VShape - expected shape of the eigenvectors
     \param VOffset - starting location for storing the eigenvectors in V
     \param tol - Input - eigenvalues with relative difference tol are treated as equal
  */
  ESCRIPT_DLL_API
  inline
  void
  eigenvalues_and_eigenvectors(const DataTypes::ValueType& in, const DataTypes::ShapeType& inShape,
                               DataTypes::ValueType::size_type inOffset,
                               DataTypes::ValueType& ev, const DataTypes::ShapeType& evShape, 
                               DataTypes::ValueType::size_type evOffset,
                               DataTypes::ValueType& V, const DataTypes::ShapeType& VShape,
                               DataTypes::ValueType::size_type VOffset,
                               const double tol=1.e-13)
  {
   double in00,in10,in20,in01,in11,in21,in02,in12,in22;
   double V00,V10,V20,V01,V11,V21,V02,V12,V22;
   double ev0,ev1,ev2;
   int s=inShape[0];
   if (s==1) {
      in00=in[inOffset+DataTypes::getRelIndex(inShape,0,0)];
      eigenvalues_and_eigenvectors1(in00,&ev0,&V00,tol);
      ev[evOffset+DataTypes::getRelIndex(evShape,0)]=ev0;
      V[inOffset+DataTypes::getRelIndex(VShape,0,0)]=V00;
   } else  if (s==2) {
      in00=in[inOffset+DataTypes::getRelIndex(inShape,0,0)];
      in10=in[inOffset+DataTypes::getRelIndex(inShape,1,0)];
      in01=in[inOffset+DataTypes::getRelIndex(inShape,0,1)];
      in11=in[inOffset+DataTypes::getRelIndex(inShape,1,1)];
      eigenvalues_and_eigenvectors2(in00,(in01+in10)/2.,in11,
                   &ev0,&ev1,&V00,&V10,&V01,&V11,tol);
      ev[evOffset+DataTypes::getRelIndex(evShape,0)]=ev0;
      ev[evOffset+DataTypes::getRelIndex(evShape,1)]=ev1;
      V[inOffset+DataTypes::getRelIndex(VShape,0,0)]=V00;
      V[inOffset+DataTypes::getRelIndex(VShape,1,0)]=V10;
      V[inOffset+DataTypes::getRelIndex(VShape,0,1)]=V01;
      V[inOffset+DataTypes::getRelIndex(VShape,1,1)]=V11;
   } else  if (s==3) {
      in00=in[inOffset+DataTypes::getRelIndex(inShape,0,0)];
      in10=in[inOffset+DataTypes::getRelIndex(inShape,1,0)];
      in20=in[inOffset+DataTypes::getRelIndex(inShape,2,0)];
      in01=in[inOffset+DataTypes::getRelIndex(inShape,0,1)];
      in11=in[inOffset+DataTypes::getRelIndex(inShape,1,1)];
      in21=in[inOffset+DataTypes::getRelIndex(inShape,2,1)];
      in02=in[inOffset+DataTypes::getRelIndex(inShape,0,2)];
      in12=in[inOffset+DataTypes::getRelIndex(inShape,1,2)];
      in22=in[inOffset+DataTypes::getRelIndex(inShape,2,2)];
      eigenvalues_and_eigenvectors3(in00,(in01+in10)/2.,(in02+in20)/2.,in11,(in21+in12)/2.,in22,
                 &ev0,&ev1,&ev2,
                 &V00,&V10,&V20,&V01,&V11,&V21,&V02,&V12,&V22,tol);
      ev[evOffset+DataTypes::getRelIndex(evShape,0)]=ev0;
      ev[evOffset+DataTypes::getRelIndex(evShape,1)]=ev1;
      ev[evOffset+DataTypes::getRelIndex(evShape,2)]=ev2;
      V[inOffset+DataTypes::getRelIndex(VShape,0,0)]=V00;
      V[inOffset+DataTypes::getRelIndex(VShape,1,0)]=V10;
      V[inOffset+DataTypes::getRelIndex(VShape,2,0)]=V20;
      V[inOffset+DataTypes::getRelIndex(VShape,0,1)]=V01;
      V[inOffset+DataTypes::getRelIndex(VShape,1,1)]=V11;
      V[inOffset+DataTypes::getRelIndex(VShape,2,1)]=V21;
      V[inOffset+DataTypes::getRelIndex(VShape,0,2)]=V02;
      V[inOffset+DataTypes::getRelIndex(VShape,1,2)]=V12;
      V[inOffset+DataTypes::getRelIndex(VShape,2,2)]=V22;

   }
 }


/**
   Inline function definitions.
*/

// template <class UnaryFunction>
// inline
// void
// DataArrayView::unaryOp(UnaryFunction operation)
// {
//   unaryOp(m_offset,operation);
// }


inline
bool
checkOffset(const DataTypes::ValueType& data,
	    const DataTypes::ShapeType& shape,
	    DataTypes::ValueType::size_type offset)
{
	return (data.size() >= (offset+DataTypes::noValues(shape))); 
}

template <class UnaryFunction>
inline
void
unaryOp(DataTypes::ValueType& data, const DataTypes::ShapeType& shape,
          DataTypes::ValueType::size_type offset,
          UnaryFunction operation)
{
  EsysAssert((data.size()>0)&&checkOffset(data,shape,offset),
               "Error - Couldn't perform unaryOp due to insufficient storage.");
  DataTypes::ValueType::size_type nVals=DataTypes::noValues(shape);
  for (DataTypes::ValueType::size_type i=0;i<nVals;i++) {
    data[offset+i]=operation(data[offset+i]);
  }
}


template <class BinaryFunction>
inline
void
binaryOp(DataTypes::ValueType& left, 
			const DataTypes::ShapeType& leftShape,
			DataTypes::ValueType::size_type leftOffset,
                        const DataTypes::ValueType& right,
			const DataTypes::ShapeType& rightShape,
                        DataTypes::ValueType::size_type rightOffset,
                        BinaryFunction operation)
{
  EsysAssert(leftShape==rightShape,
	     "Error - Couldn't perform binaryOp due to shape mismatch,");
  EsysAssert(((left.size()>0)&&checkOffset(left,leftShape, leftOffset)),
             "Error - Couldn't perform binaryOp due to insufficient storage in left object.");
  EsysAssert(((right.size()>0)&&checkOffset(right,rightShape,rightOffset)),
             "Error - Couldn't perform binaryOp due to insufficient storage in right object.");
  for (DataTypes::ValueType::size_type i=0;i<DataTypes::noValues(leftShape);i++) {
    left[leftOffset+i]=operation(left[leftOffset+i],right[rightOffset+i]);
  }
}

template <class BinaryFunction>
inline
void
binaryOp(DataTypes::ValueType& left, 
			const DataTypes::ShapeType& leftShape,
			DataTypes::ValueType::size_type offset,
                        double right,
                        BinaryFunction operation)
{
  EsysAssert(((left.size()>0)&&checkOffset(left,leftShape,offset)),
             "Error - Couldn't perform binaryOp due to insufficient storage in left object.");
  for (DataTypes::ValueType::size_type i=0;i<DataTypes::noValues(leftShape);i++) {
    left[offset+i]=operation(left[offset+i],right);
  }
}

template <class BinaryFunction>
inline
double
reductionOp(const DataTypes::ValueType& left, 
			   const DataTypes::ShapeType& leftShape,
			   DataTypes::ValueType::size_type offset,
                           BinaryFunction operation,
                           double initial_value)
{
  EsysAssert(((left.size()>0)&&checkOffset(left,leftShape,offset)),
               "Error - Couldn't perform reductionOp due to insufficient storage.");
  double current_value=initial_value;
  for (DataTypes::ValueType::size_type i=0;i<DataTypes::noValues(leftShape);i++) {
    current_value=operation(current_value,left[offset+i]);
  }
  return current_value;
}



}  // end namespace DataMath
}  // end namespace escript
#endif

