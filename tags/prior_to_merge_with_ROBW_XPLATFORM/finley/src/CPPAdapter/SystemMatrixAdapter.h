// $Id$
/* 
 ************************************************************
 *          Copyright 2006 by ACcESS MNRF                   *
 *                                                          *
 *              http://www.access.edu.au                    *
 *       Primary Business: Queensland, Australia            *
 *  Licensed under the Open Software License version 3.0    *
 *     http://www.opensource.org/licenses/osl-3.0.php       *
 *                                                          *
 ************************************************************
*/

#if !defined  finley_SystemMatrixAdapter_20040610_H
#define finley_SystemMatrixAdapter_20040610_H

extern "C" {
#include "SystemMatrix.h"
#include "Options.h"
}

#include "FinleyAdapterException.h"
#include "FinleyError.h"

#include "AbstractSystemMatrix.h"
#include "Data.h"
#include "UtilC.h"

#include <boost/python/dict.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/python/extract.hpp>

namespace finley {

class SystemMatrixAdapter:public escript::AbstractSystemMatrix {

/**
   \brief
   Wrapper for Paso_SystemMatrix. 

   Description:
   Wrapper for Paso_SystemMatrix.
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
  SystemMatrixAdapter(Paso_SystemMatrix* system_matrix,
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
  Paso_SystemMatrix* getPaso_SystemMatrix() const;

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
    nullifyRowsAndCols - calls Paso_SystemMatrix_nullifyRowsAndCols.
  */
  void nullifyRowsAndCols(escript::Data& row_q, escript::Data& col_q, const double mdv) const;

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

  /**
     \brief maps escript options onto Paso options:
  */
  static int mapOptionToPaso(const int option);

 protected:

 private:

   /**
      \brief
      solves the linear system this*out=in
   */
   virtual void setToSolution(escript::Data& out, escript::Data& in, const boost::python::dict& options) const;

   /**
       \brief
       performs y+=this*x
   */
   virtual void ypAx(escript::Data& y, escript::Data& x) const;

   //
   // pointer to the externally created finley mesh - system_matrix.
   //
   boost::shared_ptr<Paso_SystemMatrix> m_system_matrix;

};

} // end of namespace
#endif