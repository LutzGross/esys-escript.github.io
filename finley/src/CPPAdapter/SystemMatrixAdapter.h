
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


#if !defined  finley_SystemMatrixAdapter_20040610_H
#define finley_SystemMatrixAdapter_20040610_H
#include "system_dep.h"

extern "C" {
#include "paso/SystemMatrix.h"
#include "paso/Options.h"
}

#include "FinleyAdapterException.h"
#include "FinleyError.h"

#include "escript/AbstractSystemMatrix.h"
#include "escript/Data.h"
#include "escript/UtilC.h"

#include <boost/python/object.hpp>
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
  FINLEY_DLL_API
  SystemMatrixAdapter();

  /**
     /brief
     Constructor for SystemMatrixAdapter.
  */
  FINLEY_DLL_API
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
  FINLEY_DLL_API
  ~SystemMatrixAdapter();

  /**
     \brief
     Returns the pointer to the system matrix.
  */
  FINLEY_DLL_API
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
  FINLEY_DLL_API
  void nullifyRowsAndCols(escript::Data& row_q, escript::Data& col_q, const double mdv) const;

  /**
     \brief writes the matrix to a file using the Matrix Market file format
  */
  FINLEY_DLL_API
  virtual void saveMM(const std::string& fileName) const;

  /**
     \brief writes the matrix to a file using the Harwell-Boeing file format
  */
  FINLEY_DLL_API
  virtual void saveHB(const std::string& fileName) const;

  /**
     \brief sets the matrix entries to zero
  */
  FINLEY_DLL_API
  virtual void resetValues() const;

  /**
     \brief maps escript options onto Paso options:
  */
  FINLEY_DLL_API
  static int mapOptionToPaso(const int option);

  /**
     \brief extract paso options from SolutionOptions class
  */
 
  FINLEY_DLL_API
  static void escriptToPasoOptions(Paso_Options* paso_options, const boost::python::object& options);

  /**
     \brief copied diagonistic data back to the solver option.
  */
 
  FINLEY_DLL_API
  static void pasoToEscriptOptions(const Paso_Options* paso_options,boost::python::object& options);
 
  /**
     \brief prints information about a system matrix
  */
  FINLEY_DLL_API
  void Print_Matrix_Info(const bool) const;

 protected:

 private:

   /**
      \brief
      solves the linear system this*out=in
   */
   FINLEY_DLL_API
   virtual void setToSolution(escript::Data& out, escript::Data& in, boost::python::object& options) const;

   /**
       \brief
       performs y+=this*x
   */
   FINLEY_DLL_API
   virtual void ypAx(escript::Data& y, escript::Data& x) const;

   //
   // pointer to the externally created finley mesh - system_matrix.
   //
   boost::shared_ptr<Paso_SystemMatrix> m_system_matrix;

};

} // end of namespace
#endif
