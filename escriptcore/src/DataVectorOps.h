
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

#ifndef __ESCRIPT_DATAMATHS_H__
#define __ESCRIPT_DATAMATHS_H__

#include "DataAbstract.h"
#include "DataException.h"
#include "ArrayOps.h"
#include "LapackInverseHelper.h"
#include "DataTagged.h"
#include <complex>
/**
\file DataVectorOps.h
\brief Describes binary operations performed on DataVector.


For operations on DataReady see BinaryDataReadyOp.h.
For operations on double* see ArrayOps.h.
*/


namespace escript
{

/**
In order to properly identify the datapoints, in most cases, the vector, shape and offset of the point must all be supplied.
Note that vector in this context refers to a data vector storing datapoints not a mathematical vector. (However, datapoints within the data vector could represent scalars, vectors, matricies, ...).
*/


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
  matMult(const DataTypes::RealVectorType& left,
          const DataTypes::ShapeType& leftShape,
          DataTypes::RealVectorType::size_type leftOffset,
          const DataTypes::RealVectorType& right,
          const DataTypes::ShapeType& rightShape,
          DataTypes::RealVectorType::size_type rightOffset,
          DataTypes::RealVectorType& result,
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
  template<typename VEC>
  inline
  void
  symmetric(const VEC& in,
            const DataTypes::ShapeType& inShape,
            typename VEC::size_type inOffset,
            VEC& ev,
            const DataTypes::ShapeType& evShape,
            typename VEC::size_type evOffset)
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
      computes a antisymmetric matrix from your square matrix A: (A - transpose(A)) / 2

      \param in - vector containing the matrix A
      \param inShape - shape of the matrix A
      \param inOffset - the beginning of A within the vector in
      \param ev - vector to store the output matrix
      \param evShape - expected shape of the output matrix
      \param evOffset - starting location for storing ev in vector ev
   */
   template<typename VEC>
   inline
   void
   antisymmetric(const VEC& in,
            const DataTypes::ShapeType& inShape,
             typename VEC::size_type inOffset,
             VEC& ev,
            const DataTypes::ShapeType& evShape,
             typename VEC::size_type evOffset)
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
            DataTypes::CplxVectorType::size_type evOffset);

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
             typename DataTypes::CplxVectorType::size_type evOffset);

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
  template <class VEC>
  inline
  void
  trace(const VEC& in,
            const DataTypes::ShapeType& inShape,
            typename VEC::size_type inOffset,
            VEC& ev,
            const DataTypes::ShapeType& evShape,
            typename VEC::size_type evOffset,
            int axis_offset)
  {
   for (int j=0;j<DataTypes::noValues(evShape);++j)
   {
      ev[evOffset+j]=0;
   }
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
  // ESCRIPT_DLL_API can only export template specialisation!
  template <class VEC>
  inline
  void
  transpose(const VEC& in,
            const DataTypes::ShapeType& inShape,
            typename VEC::size_type inOffset,
            VEC& ev,
            const DataTypes::ShapeType& evShape,
            typename VEC::size_type evOffset,
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
  // ESCRIPT_DLL_API can only export template specialisation!
  template <class VEC>
  inline
  void
  swapaxes(const VEC& in,
           const DataTypes::ShapeType& inShape,
           typename VEC::size_type inOffset,
           VEC& ev,
           const DataTypes::ShapeType& evShape,
           typename VEC::size_type evOffset,
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
  eigenvalues(const DataTypes::RealVectorType& in,
              const DataTypes::ShapeType& inShape,
              typename DataTypes::RealVectorType::size_type inOffset,
              DataTypes::RealVectorType& ev,
              const DataTypes::ShapeType& evShape,
              typename DataTypes::RealVectorType::size_type evOffset)
  {
   typename DataTypes::RealVectorType::ElementType in00,in10,in20,in01,in11,in21,in02,in12,in22;
   typename DataTypes::RealVectorType::ElementType ev0,ev1,ev2;
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

  inline
  void
  eigenvalues(const DataTypes::CplxVectorType& in,
              const DataTypes::ShapeType& inShape,
              typename DataTypes::CplxVectorType::size_type inOffset,
              DataTypes::CplxVectorType& ev,
              const DataTypes::ShapeType& evShape,
              typename DataTypes::CplxVectorType::size_type evOffset)
  {
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-variable"
   typename DataTypes::CplxVectorType::ElementType in00,in10,in20,in01,in11,in21,in02,in12,in22;
   typename DataTypes::CplxVectorType::ElementType ev0,ev1,ev2;
#pragma clang diagnostic pop
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
     // this doesn't work yet
//       in00=in[inOffset+DataTypes::getRelIndex(inShape,0,0)];
//       in10=in[inOffset+DataTypes::getRelIndex(inShape,1,0)];
//       in20=in[inOffset+DataTypes::getRelIndex(inShape,2,0)];
//       in01=in[inOffset+DataTypes::getRelIndex(inShape,0,1)];
//       in11=in[inOffset+DataTypes::getRelIndex(inShape,1,1)];
//       in21=in[inOffset+DataTypes::getRelIndex(inShape,2,1)];
//       in02=in[inOffset+DataTypes::getRelIndex(inShape,0,2)];
//       in12=in[inOffset+DataTypes::getRelIndex(inShape,1,2)];
//       in22=in[inOffset+DataTypes::getRelIndex(inShape,2,2)];
//       eigenvalues3(in00,(in01+in10)/2.,(in02+in20)/2.,in11,(in21+in12)/2.,in22,
//                  &ev0,&ev1,&ev2);
//       ev[evOffset+DataTypes::getRelIndex(evShape,0)]=ev0;
//       ev[evOffset+DataTypes::getRelIndex(evShape,1)]=ev1;
//       ev[evOffset+DataTypes::getRelIndex(evShape,2)]=ev2;

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
  eigenvalues_and_eigenvectors(const DataTypes::RealVectorType& in, const DataTypes::ShapeType& inShape,
                               DataTypes::RealVectorType::size_type inOffset,
                               DataTypes::RealVectorType& ev, const DataTypes::ShapeType& evShape,
                               DataTypes::RealVectorType::size_type evOffset,
                               DataTypes::RealVectorType& V, const DataTypes::ShapeType& VShape,
                               DataTypes::RealVectorType::size_type VOffset,
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

template <class VEC>
inline
bool
checkOffset(const VEC& data,
            const DataTypes::ShapeType& shape,
            typename VEC::size_type offset)
{
        return (data.size() >= (offset+DataTypes::noValues(shape)));
}

/**
 * This assumes that all data involved have the same points per sample and same shape
*/
template <class ResVEC, class LVEC, class RSCALAR>
void
binaryOpVectorRightScalar(ResVEC& res,                          // where result is to be stored
          typename ResVEC::size_type resOffset,         // offset in the result vector to start storing results
          const typename ResVEC::size_type samplesToProcess,    // number of samples to be updated in the result
          const typename ResVEC::size_type sampleSize,          // number of values in each sample
          const LVEC& left,                             // LHS of calculation
          typename LVEC::size_type leftOffset,          // where to start reading LHS values
          const RSCALAR* right,                         // RHS of the calculation
          const bool rightreset,                        // true if RHS is providing a single sample of 1 value only
          escript::ES_optype operation,         // operation to perform
          bool singleleftsample)                        // set to false for normal operation
{
    size_t substep=(rightreset?0:1);
    switch (operation)
    {
        case ADD:
        {
#pragma omp parallel for
          for (typename ResVEC::size_type i=0;i<samplesToProcess;++i)
          {
              typename LVEC::size_type leftbase=leftOffset+(singleleftsample?0:i*sampleSize);
              const RSCALAR* rpos=right+(rightreset?0:i*substep);

              for (typename ResVEC::size_type j=0;j<sampleSize;++j)
              {
                  res[i*sampleSize+resOffset+j]=left[leftbase+j]+*rpos;
              }
          }
        }
        break;
        case POW:
        {
#pragma omp parallel for
          for (typename ResVEC::size_type i=0;i<samplesToProcess;++i)
          {
              typename LVEC::size_type leftbase=leftOffset+(singleleftsample?0:i*sampleSize);
              const RSCALAR* rpos=right+(rightreset?0:i*substep);

              for (typename ResVEC::size_type j=0;j<sampleSize;++j)
              {
                  res[i*sampleSize+resOffset+j]=pow(left[leftbase+j],*rpos);
              }
          }
        }
        break;
        case SUB:
        {
#pragma omp parallel for
          for (typename ResVEC::size_type i=0;i<samplesToProcess;++i)
          {
              typename LVEC::size_type leftbase=leftOffset+(singleleftsample?0:i*sampleSize);
              const RSCALAR* rpos=right+(rightreset?0:i*substep);

              for (typename ResVEC::size_type j=0;j<sampleSize;++j)
              {
                  res[i*sampleSize+resOffset+j]=left[leftbase+j]-*rpos;
              }
          }
        }
        break;
        case MUL:
        {
#pragma omp parallel for
          for (typename ResVEC::size_type i=0;i<samplesToProcess;++i)
          {
              typename LVEC::size_type leftbase=leftOffset+(singleleftsample?0:i*sampleSize);
              const RSCALAR* rpos=right+(rightreset?0:i*substep);

              for (typename ResVEC::size_type j=0;j<sampleSize;++j)
              {
                  res[i*sampleSize+resOffset+j]=left[leftbase+j] * *rpos;
              }
          }
        }
        break;
        case DIV:
        {
#pragma omp parallel for
          for (typename ResVEC::size_type i=0;i<samplesToProcess;++i)
          {
              typename LVEC::size_type leftbase=leftOffset+(singleleftsample?0:i*sampleSize);
              const RSCALAR* rpos=right+(rightreset?0:i*substep);

              for (typename ResVEC::size_type j=0;j<sampleSize;++j)
              {
                  res[i*sampleSize+resOffset+j]=left[leftbase+j]/ *rpos;
              }
          }
        }
        break;
        default:
            throw DataException("Unsupported binary operation");
    }
}

template<>
void
binaryOpVectorRightScalar(DataTypes::RealVectorType& res,                               // where result is to be stored
          typename DataTypes::RealVectorType::size_type resOffset,              // offset in the result vector to start storing results
          const typename DataTypes::RealVectorType::size_type samplesToProcess, // number of samples to be updated in the result
          const typename DataTypes::RealVectorType::size_type sampleSize,               // number of values in each sample
          const DataTypes::RealVectorType& left,                                // LHS of calculation
          typename DataTypes::RealVectorType::size_type leftOffset,             // where to start reading LHS values
          const DataTypes::real_t* right,                       // RHS of the calculation
          const bool rightreset,                        // true if RHS is providing a single sample of 1 value only
          escript::ES_optype operation,         // operation to perform
          bool singleleftsample);

/**
 * This assumes that all data involved have the same points per sample and same shape
*/
template <class ResVEC, class LSCALAR, class RVEC>
void
binaryOpVectorLeftScalar(ResVEC& res,                           // where result is to be stored
          typename ResVEC::size_type resOffset,         // offset in the result vector to start storing results
          const typename ResVEC::size_type samplesToProcess,    // number of samples to be updated in the result
          const typename ResVEC::size_type sampleSize,          // number of values in each sample
          const LSCALAR* left,                          // LHS of calculation
          const bool leftreset,                         // true if LHS is providing a single sample of 1 value only
          const RVEC& right,                            // RHS of the calculation
          typename RVEC::size_type rightOffset,         // where to start reading RHS values
          escript::ES_optype operation,         // operation to perform
          bool singlerightsample)                       // right consists of a single sample
{
    size_t substep=(leftreset?0:1);
    switch (operation)
    {
        case ADD:
        {
#pragma omp parallel for
          for (typename ResVEC::size_type i=0;i<samplesToProcess;++i)
          {
              typename RVEC::size_type rightbase=rightOffset+(singlerightsample?0:i*sampleSize);
              const LSCALAR* lpos=left+(leftreset?0:i*substep);
              for (typename ResVEC::size_type j=0;j<sampleSize;++j)
              {
                  res[i*sampleSize+resOffset+j]=*lpos+right[rightbase+j];
              }
          }
        }
        break;
        case POW:
        {
#pragma omp parallel for
          for (typename ResVEC::size_type i=0;i<samplesToProcess;++i)
          {
              typename RVEC::size_type rightbase=rightOffset+(singlerightsample?0:i*sampleSize);
              const LSCALAR* lpos=left+(leftreset?0:i*substep);
              for (typename ResVEC::size_type j=0;j<sampleSize;++j)
              {
                  res[i*sampleSize+resOffset+j]=pow(*lpos,right[rightbase+j]);
              }
          }
        }
        break;
        case SUB:
        {
#pragma omp parallel for
          for (typename ResVEC::size_type i=0;i<samplesToProcess;++i)
          {
              typename RVEC::size_type rightbase=rightOffset+(singlerightsample?0:i*sampleSize);
              const LSCALAR* lpos=left+(leftreset?0:i*substep);
              for (typename ResVEC::size_type j=0;j<sampleSize;++j)
              {
                  res[i*sampleSize+resOffset+j]=*lpos-right[rightbase+j];
              }
          }
        }
        break;
        case MUL:
        {
#pragma omp parallel for
          for (typename ResVEC::size_type i=0;i<samplesToProcess;++i)
          {
              typename RVEC::size_type rightbase=rightOffset+(singlerightsample?0:i*sampleSize);
              const LSCALAR* lpos=left+(leftreset?0:i*substep);
              for (typename ResVEC::size_type j=0;j<sampleSize;++j)
              {
                  res[i*sampleSize+resOffset+j]=*lpos*right[rightbase+j];
              }
          }
        }
        break;
        case DIV:
        {
#pragma omp parallel for
          for (typename ResVEC::size_type i=0;i<samplesToProcess;++i)
          {
              typename RVEC::size_type rightbase=rightOffset+(singlerightsample?0:i*sampleSize);
              const LSCALAR* lpos=left+(leftreset?0:i*substep);
              for (typename ResVEC::size_type j=0;j<sampleSize;++j)
              {
                  res[i*sampleSize+resOffset+j]=*lpos/right[rightbase+j];
              }
          }
        }
        break;
        default:
            throw DataException("Unsupported binary operation");
    }
}

template <>
void
binaryOpVectorLeftScalar(DataTypes::RealVectorType& res,                                // where result is to be stored
          typename DataTypes::RealVectorType::size_type resOffset,              // offset in the result vector to start storing results
          const typename DataTypes::RealVectorType::size_type samplesToProcess, // number of samples to be updated in the result
          const typename DataTypes::RealVectorType::size_type sampleSize,               // number of values in each sample
          const DataTypes::real_t* left,                                // LHS of calculation
          const bool leftreset,                         // true if LHS is providing a single sample of 1 value only
          const DataTypes::RealVectorType& right,                               // RHS of the calculation
          typename DataTypes::RealVectorType::size_type rightOffset,            // where to start reading RHS values
          escript::ES_optype operation,         // operation to perform
          bool singlerightsample);                      // right consists of a single sample

/**
 * This assumes that all data involved have the same points per sample and same shape
*/
template <class ResVEC, class LVEC, class RVEC>
void
binaryOpVector(ResVEC& res,                             // where result is to be stored
          typename ResVEC::size_type resOffset,         // offset in the result vector to start storing results
          const typename ResVEC::size_type samplesToProcess,    // number of samples to be updated in the result
          const typename ResVEC::size_type sampleSize,          // number of values in each sample
          const LVEC& left,                             // LHS of calculation
          typename LVEC::size_type leftOffset,          // where to start reading LHS values
          const bool leftreset,                         // Is LHS only supplying a single sample instead of a bunch of them
          const RVEC& right,                            // RHS of the calculation
          typename RVEC::size_type rightOffset,         // where to start reading RHS values
          const bool rightreset,                        // Is RHS only supplying a single sample instead of a bunch of them
          escript::ES_optype operation)         // operation to perform
{
    switch (operation)
    {
        case ADD:
        {
#pragma omp parallel for
          for (typename ResVEC::size_type i=0;i<samplesToProcess;++i)
          {
              typename LVEC::size_type leftbase=leftOffset+(leftreset?0:i*sampleSize);
              typename RVEC::size_type rightbase=rightOffset+(rightreset?0:i*sampleSize);
              for (typename ResVEC::size_type j=0;j<sampleSize;++j)
              {
                  res[i*sampleSize+resOffset+j]=left[leftbase+j]+right[rightbase+j];
              }
          }
        }
        break;
        case POW:
        {
#pragma omp parallel for
          for (typename ResVEC::size_type i=0;i<samplesToProcess;++i)
          {
              typename LVEC::size_type leftbase=leftOffset+(leftreset?0:i*sampleSize);
              typename RVEC::size_type rightbase=rightOffset+(rightreset?0:i*sampleSize);
              for (typename ResVEC::size_type j=0;j<sampleSize;++j)
              {
                  res[i*sampleSize+resOffset+j]=pow(left[leftbase+j],right[rightbase+j]);
              }
          }
        }
        break;
        case SUB:
        {
#pragma omp parallel for
          for (typename ResVEC::size_type i=0;i<samplesToProcess;++i)
          {
              typename LVEC::size_type leftbase=leftOffset+(leftreset?0:i*sampleSize);
              typename RVEC::size_type rightbase=rightOffset+(rightreset?0:i*sampleSize);
              for (typename ResVEC::size_type j=0;j<sampleSize;++j)
              {
                  res[i*sampleSize+resOffset+j]=left[leftbase+j]-right[rightbase+j];
              }
          }
        }
        break;
        case MUL:
        {
#pragma omp parallel for
          for (typename ResVEC::size_type i=0;i<samplesToProcess;++i)
          {
              typename LVEC::size_type leftbase=leftOffset+(leftreset?0:i*sampleSize);
              typename RVEC::size_type rightbase=rightOffset+(rightreset?0:i*sampleSize);
              for (typename ResVEC::size_type j=0;j<sampleSize;++j)
              {
                  res[i*sampleSize+resOffset+j]=left[leftbase+j]*right[rightbase+j];
              }
          }
        }
        break;
        case DIV:
        {
#pragma omp parallel for
          for (typename ResVEC::size_type i=0;i<samplesToProcess;++i)
          {
              typename LVEC::size_type leftbase=leftOffset+(leftreset?0:i*sampleSize);
              typename RVEC::size_type rightbase=rightOffset+(rightreset?0:i*sampleSize);
              for (typename ResVEC::size_type j=0;j<sampleSize;++j)
              {
                  res[i*sampleSize+resOffset+j]=left[leftbase+j]/right[rightbase+j];
              }
          }
        }
        break;
        default:
            throw DataException("Unsupported binary operation");
    }
}

template <>
void
binaryOpVector(DataTypes::RealVectorType& res,                          // where result is to be stored
          typename DataTypes::RealVectorType::size_type resOffset,              // offset in the result vector to start storing results
          const typename DataTypes::RealVectorType::size_type samplesToProcess, // number of samples to be updated in the result
          const typename DataTypes::RealVectorType::size_type sampleSize,               // number of values in each sample
          const DataTypes::RealVectorType& left,                                // LHS of calculation
          typename DataTypes::RealVectorType::size_type leftOffset,             // where to start reading LHS values
          const bool leftreset,                         // Is LHS only supplying a single sample instead of a bunch of them
          const DataTypes::RealVectorType& right,                               // RHS of the calculation
          typename DataTypes::RealVectorType::size_type rightOffset,            // where to start reading RHS values
          const bool rightreset,                        // Is RHS only supplying a single sample instead of a bunch of them
          escript::ES_optype operation);                // operation to perform

#define OPVECLAZYBODY(X)     \
    for (size_t j=0;j<onumsteps;++j)\
    {\
      for (size_t i=0;i<numsteps;++i,res+=resultStep) \
      { \
          for (size_t s=0; s<chunksize; ++s)\
          {\
              res[s] = X;\
          }\
/*        tensor_binary_operation< TYPE >(chunksize, &((*left)[lroffset]), &((*right)[rroffset]), resultp, X);*/ \
          lroffset+=leftstep; \
          rroffset+=rightstep; \
      }\
      lroffset+=oleftstep;\
      rroffset+=orightstep;\
    }

/**
 * This assumes that all data involved have the same points per sample and same shape
 * This version is to be called from within DataLazy.
 * It does not have openmp around loops because it will be evaluating individual samples
 * (Which will be done within an enclosing openmp region.
*/
template <class ResELT, class LELT, class RELT>
void
binaryOpVectorLazyArithmeticHelper(ResELT* res,
                         const LELT* left,
                         const RELT* right,
                         const size_t chunksize,
                         const size_t onumsteps,
                         const size_t numsteps,
                         const size_t resultStep,
                         const size_t leftstep,
                         const size_t rightstep,
                         const size_t oleftstep,
                         const size_t orightstep,
                         size_t lroffset,
                         size_t rroffset,
                         escript::ES_optype operation)          // operation to perform
{
    switch (operation)
    {
        case ADD:
            OPVECLAZYBODY((left[lroffset+s]+right[rroffset+s]));
            break;
        case POW:
            OPVECLAZYBODY(pow(left[lroffset+s],right[rroffset+s]))
            break;
        case SUB:
            OPVECLAZYBODY(left[lroffset+s]-right[rroffset+s])
            break;
        case MUL:
            OPVECLAZYBODY(left[lroffset+s]*right[rroffset+s])
            break;
        case DIV:
            OPVECLAZYBODY(left[lroffset+s]/right[rroffset+s])
            break;
        default:
            ESYS_ASSERT(false, "Invalid operation. This should never happen!");
            // I can't throw here because this will be called inside a parallel section
    }
}


/**
 * This assumes that all data involved have the same points per sample and same shape
 * This version is to be called from within DataLazy.
 * It does not have openmp around loops because it will be evaluating individual samples
 * (Which will be done within an enclosing openmp region.
*/
template <class ResELT, class LELT, class RELT>
void
binaryOpVectorLazyRelationalHelper(ResELT* res,
                         const LELT* left,
                         const RELT* right,
                         const size_t chunksize,
                         const size_t onumsteps,
                         const size_t numsteps,
                         const size_t resultStep,
                         const size_t leftstep,
                         const size_t rightstep,
                         const size_t oleftstep,
                         const size_t orightstep,
                         size_t lroffset,
                         size_t rroffset,
                         escript::ES_optype operation)          // operation to perform
{
    switch (operation)
    {
        case LESS:
            OPVECLAZYBODY(left[lroffset+s]<right[rroffset+s])
            break;
        case GREATER:
            OPVECLAZYBODY(left[lroffset+s]>right[rroffset+s])
            break;
        case GREATER_EQUAL:
            OPVECLAZYBODY(left[lroffset+s]>=right[rroffset+s])
            break;
        case LESS_EQUAL:
            OPVECLAZYBODY(left[lroffset+s]<=right[rroffset+s])
            break;
        default:
            ESYS_ASSERT(false, "Invalid operation. This should never happen!");
            // I can't throw here because this will be called inside a parallel section
    }
}

/**
 * This assumes that all data involved have the same points per sample and same shape
*/
/* trying to make a single version for all Tagged+Expanded interactions */
template <class ResVEC, class LVEC, class RVEC>
void
binaryOpVectorTagged(ResVEC& res,                               // where result is to be stored
          const typename ResVEC::size_type samplesToProcess,    // number of samples to be updated in the result
          const typename ResVEC::size_type DPPSample, // number of datapoints per sample
          const typename ResVEC::size_type DPSize,    // datapoint size
          const LVEC& left,                           // LHS of calculation
          bool leftscalar,
          const RVEC& right,                          // RHS of the calculation
          bool rightscalar,
          bool lefttagged,                            // true if left object is the tagged one
          const DataTagged& tagsource,                // where to get tag offsets from
          escript::ES_optype operation)               // operation to perform
{
    typename ResVEC::size_type lstep=leftscalar?1:DPSize;
    typename ResVEC::size_type rstep=rightscalar?1:DPSize;
    typename ResVEC::size_type limit=samplesToProcess*DPPSample;
    switch (operation)
    {
        case ADD:
        {
#pragma omp parallel for
          for (typename ResVEC::size_type i=0;i<limit;++i)
          {
              typename LVEC::size_type leftbase=(lefttagged?tagsource.getPointOffset(i/DPPSample,0):i*lstep);       // only one of these
              typename RVEC::size_type rightbase=(lefttagged?i*rstep:tagsource.getPointOffset(i/DPPSample,0));      // will apply

              for (typename ResVEC::size_type j=0;j<DPSize;++j)
              {
                  res[i*DPSize+j]=left[leftbase+j*(!leftscalar)]+right[rightbase+j*(!rightscalar)];
              }
          }
        }
        break;
        case POW:
        {
#pragma omp parallel for
          for (typename ResVEC::size_type i=0;i<limit;++i)
          {
              typename LVEC::size_type leftbase=(lefttagged?tagsource.getPointOffset(i/DPPSample,0):i*lstep);       // only one of these
              typename RVEC::size_type rightbase=(lefttagged?i*rstep:tagsource.getPointOffset(i/DPPSample,0));      // will apply

              for (typename ResVEC::size_type j=0;j<DPSize;++j)
              {
                  res[i*DPSize+j]=pow(left[leftbase+j*(!leftscalar)],right[rightbase+j*(!rightscalar)]);
              }
          }
        }
        break;
        case SUB:
        {
#pragma omp parallel for
          for (typename ResVEC::size_type i=0;i<limit;++i)
          {
              typename LVEC::size_type leftbase=(lefttagged?tagsource.getPointOffset(i/DPPSample,0):i*lstep);       // only one of these
              typename RVEC::size_type rightbase=(lefttagged?i*rstep:tagsource.getPointOffset(i/DPPSample,0));      // will apply

              for (typename ResVEC::size_type j=0;j<DPSize;++j)
              {
                  res[i*DPSize+j]=left[leftbase+j*(!leftscalar)]-right[rightbase+j*(!rightscalar)];
              }
          }
        }
        break;
        case MUL:
        {
#pragma omp parallel for
          for (typename ResVEC::size_type i=0;i<limit;++i)
          {
              typename LVEC::size_type leftbase=(lefttagged?tagsource.getPointOffset(i/DPPSample,0):i*lstep);       // only one of these
              typename RVEC::size_type rightbase=(lefttagged?i*rstep:tagsource.getPointOffset(i/DPPSample,0));      // will apply

              for (typename ResVEC::size_type j=0;j<DPSize;++j)
              {
                  res[i*DPSize+j]=left[leftbase+j*(!leftscalar)]*right[rightbase+j*(!rightscalar)];
              }
          }
        }
        break;
        case DIV:
        {
#pragma omp parallel for
          for (typename ResVEC::size_type i=0;i<limit;++i)
          {
              typename LVEC::size_type leftbase=(lefttagged?tagsource.getPointOffset(i/DPPSample,0):i*lstep);       // only one of these
              typename RVEC::size_type rightbase=(lefttagged?i*rstep:tagsource.getPointOffset(i/DPPSample,0));      // will apply

              for (typename ResVEC::size_type j=0;j<DPSize;++j)
              {
                  res[i*DPSize+j]=left[leftbase+j*(!leftscalar)]/right[rightbase+j*(!rightscalar)];
              }
          }
        }
        break;
        default:
            throw DataException("Unsupported binary operation");
    }
}

template<>
void
binaryOpVectorTagged(DataTypes::RealVectorType& res,                            // where result is to be stored
          const typename DataTypes::RealVectorType::size_type samplesToProcess, // number of samples to be updated in the result
          const typename DataTypes::RealVectorType::size_type DPPSample,        // number of datapoints per sample
          const typename DataTypes::RealVectorType::size_type DPSize,           // datapoint size
          const DataTypes::RealVectorType& left,                                // LHS of calculation
          const bool leftscalar,
          const DataTypes::RealVectorType& right,                               // RHS of the calculation
          const bool rightscalar,
          const bool lefttagged,                        // true if left object is the tagged one
          const DataTagged& tagsource,                  // where to get tag offsets from
          escript::ES_optype operation);




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
inline
DataTypes::real_t
reductionOpVector(const DataTypes::RealVectorType& left,
                           const DataTypes::ShapeType& leftShape,
                           DataTypes::RealVectorType::size_type offset,
                           BinaryFunction operation,
                           DataTypes::real_t initial_value)
{
  ESYS_ASSERT((left.size()>0)&&checkOffset(left,leftShape,offset),
         "Couldn't perform reductionOp due to insufficient storage.");
  DataTypes::real_t current_value=initial_value;
  for (DataTypes::RealVectorType::size_type i=0;i<DataTypes::noValues(leftShape);i++) {
    current_value=operation(current_value,left[offset+i]);
  }
  return current_value;
}

template <class BinaryFunction>
inline
DataTypes::real_t
reductionOpVector(const DataTypes::CplxVectorType& left,
                           const DataTypes::ShapeType& leftShape,
                           DataTypes::CplxVectorType::size_type offset,
                           BinaryFunction operation,
                           DataTypes::real_t initial_value)
{
  ESYS_ASSERT((left.size()>0)&&checkOffset(left,leftShape,offset),
         "Couldn't perform reductionOp due to insufficient storage.");
  DataTypes::real_t current_value=initial_value;
  for (DataTypes::RealVectorType::size_type i=0;i<DataTypes::noValues(leftShape);i++) {
    current_value=operation(current_value,left[offset+i]);
  }
  return current_value;
}


/**
     \brief
     computes the inverses of square (up to 3x3) matricies

     \param in - vector containing the input matricies
     \param inShape - shape of the input matricies
     \param inOffset - the beginning of the input matricies within the vector "in"
     \param out - vector to store the inverses
     \param outShape - expected shape of the inverses
     \param outOffset - starting location for storing the inverses in out
     \param count - number of matricies to invert
     \param helper - associated working storage

     \exception DataException if input and output are not the correct shape or if any of the matricies are not invertible.
     \return 0 on success, on failure the return value should be passed to matrixInverseError(int err).
*/
int
matrix_inverse(const DataTypes::RealVectorType& in,
            const DataTypes::ShapeType& inShape,
            DataTypes::RealVectorType::size_type inOffset,
            DataTypes::RealVectorType& out,
            const DataTypes::ShapeType& outShape,
            DataTypes::RealVectorType::size_type outOffset,
            int count,
            LapackInverseHelper& helper);

/**
   \brief
   throws an appropriate exception based on failure of matrix_inverse.

   \param err - error code returned from matrix_inverse
   \warning do not call in a parallel region since it throws.
*/
void
matrixInverseError(int err);

/**
   \brief returns true if the vector contains NaN

*/
inline
bool
vectorHasNaN(const DataTypes::RealVectorType& in, DataTypes::RealVectorType::size_type inOffset, size_t count)
{
        for (size_t z=inOffset;z<inOffset+count;++z)
        {
            if (nancheck(in[z]))
            {
                return true;
            }
        }
        return false;
}

inline
bool
vectorHasNaN(const DataTypes::CplxVectorType& in, DataTypes::CplxVectorType::size_type inOffset, size_t count)
{
        for (size_t z=inOffset;z<inOffset+count;++z)
        {
            if (nancheck(in[z]))
            {
                return true;
            }
        }
        return false;
}


}  // end namespace escript

#endif // __ESCRIPT_DATAMATHS_H__

