/*
// @HEADER
// ************************************************************************
//             FEI: Finite Element Interface to Linear Solvers
//                  Copyright (2005) Sandia Corporation.
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Alan Williams (william@sandia.gov) 
//
// ************************************************************************
// @HEADER
*/


#ifndef _fei_MatrixTraits_Aztec_hpp_
#define _fei_MatrixTraits_Aztec_hpp_

#include <fei_trilinos_macros.hpp>

#ifdef HAVE_FEI_AZTECOO

//
//IMPORTANT NOTE: Make sure that wherever this file is included from, it
//appears BEFORE any include of fei_Vector_Impl.hpp or fei_Matrix_Impl.hpp !!!
//
#include <fei_MatrixTraits.hpp>
#include <snl_fei_BlockMatrixTraits.hpp>
#include <fei_VectorTraits_Aztec.hpp>
#include <fei_Include_Trilinos.hpp>
#include <fei_Vector_Impl.hpp>

namespace fei {
  /** Declare an AztecDMSR_Matrix specialization of the
      fei::MatrixTraits struct.

      This allows AztecDMSR_Matrix to be used as the template parameter
      of the fei::Matrix class.
  */
  template<>
  struct MatrixTraits<AztecDMSR_Matrix> {
    static const char* typeName()
      { return("fei::AztecDMSR_Matrix"); }

    static int setValues(AztecDMSR_Matrix* mat, double scalar)
      {
        return( mat->put(scalar) );
      }

    static int getNumLocalRows(AztecDMSR_Matrix* mat, int& numRows)
    {
      numRows = mat->getAztec_Map().localSize();
      return(0);
    }

    static int getRowLength(AztecDMSR_Matrix* mat, int row, int& length)
      {
        length = mat->rowLength(row);
        if (length < 0) return(-1);
        return( 0 );
      }

    static int copyOutRow(AztecDMSR_Matrix* mat,
                      int row, int len, double* coefs, int* indices)
      {
        int dummy;
        mat->getRow(row, dummy, coefs, indices);
        return(0);
      }

    static int putValuesIn(AztecDMSR_Matrix* mat,
                     int numRows, const int* rows,
                     int numCols, const int* cols,
                     const double* const* values,
                           bool sum_into)
      {
        int err = 0;
        if (sum_into) {
          err = mat->sumIntoRow(numRows, rows, numCols, cols, coefs);
        }
        else {
          for(int i=0; i<numRows; ++i) {
            err = mat->putRow(rows[i], numCols, values[i], cols);
            if (err != 0) {
              return(err);
            }
          }
        }
        return(err);
      }

    static int globalAssemble(AztecDMSR_Matrix* mat)
    {
      if (!mat->isFilled()) {
        int err = mat->fillComplete();
        if (err != 0) {
          fei::console_out() << "MatrixTraits<AztecDMSR_Matrix>::globalAssemble"
                   << " ERROR in mat->fillComplete" << FEI_ENDL;
          return(-1);
        }
      }

      return( 0 );
    }

    static int matvec(AztecDMSR_Matrix* mat,
                      fei::Vector* x,
                      fei::Vector* y)
    {
      fei::Vector_Impl<Aztec_LSVector>* avx =
        dynamic_cast<fei::Vector_Impl<Aztec_LSVector>* >(x);
      fei::Vector_Impl<Aztec_LSVector>* avy =
        dynamic_cast<fei::Vector_Impl<Aztec_LSVector>* >(y);

      if (avx == NULL || avy == NULL) {
        return(-1);
      }

      Aztec_LSVector* ax = avx->getUnderlyingVector();
      Aztec_LSVector* ay = avy->getUnderlyingVector();

      return( mat->matvec(*ax, *ay) );
    }

  };//struct MatrixTraits<AztecDMSR_Matrix>
}//namespace fei

namespace snl_fei {
  /** Declare an fei::AztecDVBR_Matrix specialization of the
      snl_fei::BlockMatrixTraits struct.

      This allows fei::AztecDVBR_Matrix to be used as the template parameter
      for the fei::Matrix class.
  */
  template<>
  struct BlockMatrixTraits<AztecDVBR_Matrix> {
    static const char* typeName()
      { return("fei::AztecDVBR_Matrix"); }

    static int putScalar(AztecDVBR_Matrix* mat, double scalar)
      {
        return( mat->put(scalar) );
      }

    static int getRowLength(AztecDVBR_Matrix* mat, int row, int& length)
      {
        return(-1);
      }

    static int getPointRowLength(AztecDVBR_Matrix* mat, int row, int& length)
    {
      return(-1);
    }

    static int copyOutRow(AztecDVBR_Matrix* mat,
                          int row, int numBlkCols,
                          int rowDim,
                          int* blkCols,
                          int* colDims,
                          double* coefs,
                          int coefsLen,
                          int& blkRowLength)
      {
        return(-1);
      }

    static int copyOutPointRow(AztecDVBR_Matrix* mat,
                               int firstLocalOffset,
                               int row,
                               int len,
                               double* coefs,
                               int* indices,
                               int& rowLength)
      {
        return(-1);
      }

    static int sumIn(AztecDVBR_Matrix* mat,
                     int blockRow,
                     int rowDim,
                     int numBlockCols,
                     const int* blockCols,
                     const int* colDims,
                     int LDA,
                     const double* values)
    {
      return(-1);
    }

    static int copyIn(AztecDVBR_Matrix* mat,
                      int blockRow,
                      int rowDim,
                      int numBlockCols,
                      const int* blockCols,
                      const int* colDims,
                      int LDA,
                      const double* values)
    {
      return(-1);
    }

    static int sumIn(AztecDVBR_Matrix* mat,
                     int row,
                     int rowDim,
                     int numCols,
                     const int* cols,
                     const int* LDAs,
                     const int* colDims,
                     const double* const* values)
      {
        return(-1);
      }

    static int copyIn(AztecDVBR_Matrix* mat,
                      int row,
                      int rowDim,
                      int numCols,
                      const int* cols,
                      const int* LDAs,
                      const int* colDims,
                      const double* const* values)
      {
        return(-1);
      }

    static int globalAssemble(AztecDVBR_Matrix* mat)
    {
      return( mat->loadComplete() );
    }
  };//struct BlockMatrixTraits<AztecDVBR_Matrix>
}//namespace snl_fei
#endif //HAVE_FEI_AZTECOO
#endif // _fei_MatrixTraits_Aztec_hpp_
