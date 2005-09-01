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
                                                                           
#if !defined  finley_SystemMatrixAdapter_20040610_H
#define finley_SystemMatrixAdapter_20040610_H

#include "finley/CPPAdapter/SystemMatrixAdapter.h"
#include "escript/Data/AbstractSystemMatrix.h"
#include "escript/Data/Data.h"
extern "C" {
#include "finley/finleyC/System.h"
}
#include <boost/python/dict.hpp>
#include <boost/shared_ptr.hpp>

namespace finley {

class SystemMatrixAdapter:public escript::AbstractSystemMatrix {

/**
   \brief
   Wrapper for Finley_SystemMatrix. 

   Description:
   Wrapper for Finley_SystemMatrix.
*/

 public:

  /**
     /brief
     Default Constructor for SystemMatrixAdapter.
     NB: Only throws an exception.
  */
  SystemMatrixAdapter();

  /**
     /brief
     Constructor for SystemMatrixAdapter.
  */
  SystemMatrixAdapter(Finley_SystemMatrix* system_matrix,
                      const int row_blocksize,
                      const escript::FunctionSpace& row_functionspace,
                      const int column_blocksize,
                      const escript::FunctionSpace& colum_functionspace);


  /**
     \brief
     Destructor for SystemMatrixAdapter. As specified in the constructor
     this deallocates the pointer given to the constructor.
  */
  ~SystemMatrixAdapter();

  /**
     \brief
     Returns the pointer to the system matrix.
  */
  Finley_SystemMatrix* getFinley_SystemMatrix() const;

  /**
     \brief
     Returns the system matrix as a const AbstractSystemMatrix&.
  */
  inline const escript::AbstractSystemMatrix& asAbstractSystemMatrix() const
  {
     return dynamic_cast<const escript::AbstractSystemMatrix&>(*this);
  }

  /**
     \brief
     Returns a system matrix as a const SystemMatrixAdapter&.
  */
  inline static const SystemMatrixAdapter& asSystemMatrixAdapter(const AbstractSystemMatrix& systemmatrix)
  {
     return dynamic_cast<const SystemMatrixAdapter&>(systemmatrix);
  }

  /**
    \brief
    nullifyRowsAndCols - calls Finley_SystemMatrix_nullifyRowsAndCols.
  */
  void nullifyRowsAndCols(const escript::Data& row_q, const escript::Data& col_q, const double mdv) const;

  /**
     \brief writes the matrix to a file using the Matrix Market file format
  */
  virtual void saveMM(const std::string& fileName) const;

  /**
     \brief writes the matrix to a file using the Harwell-Boeing file format
  */
  virtual void saveHB(const std::string& fileName) const;

  /**
     \brief sets the matrix entries to zero
  */
  virtual void resetValues() const;

 protected:

 private:

   /**
      \brief
      solves the linear system this*out=in
   */
   virtual void setToSolution(escript::Data& out, const escript::Data& in, const boost::python::dict& options) const;

   /**
       \brief
       performs y+=this*x
   */
   virtual void ypAx(escript::Data& y, const escript::Data& x) const;

   //
   // pointer to the externally created finley mesh - system_matrix.
   //
   boost::shared_ptr<Finley_SystemMatrix> m_system_matrix;

};

} // end of namespace
#endif
