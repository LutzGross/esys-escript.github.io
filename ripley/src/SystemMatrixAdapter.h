
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

#ifndef __RIPLEY_SYSTEMMATRIXADAPTER_H__
#define __RIPLEY_SYSTEMMATRIXADAPTER_H__

#include <ripley/system_dep.h>

extern "C" {
#include <paso/SystemMatrix.h>
#include <paso/Options.h>
}

#include <escript/AbstractSystemMatrix.h>
#include <escript/Data.h>
#include <escript/UtilC.h>

#include <boost/python/object.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/python/extract.hpp>

namespace ripley {

class SystemMatrixAdapter: public escript::AbstractSystemMatrix {

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
  */
  RIPLEY_DLL_API
  SystemMatrixAdapter();

  /**
     /brief
     Constructor for SystemMatrixAdapter.
  */
  RIPLEY_DLL_API
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
  RIPLEY_DLL_API
  ~SystemMatrixAdapter();

  /**
     \brief
     Returns the pointer to the system matrix.
  */
  RIPLEY_DLL_API
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
  RIPLEY_DLL_API
  void nullifyRowsAndCols(escript::Data& row_q, escript::Data& col_q, const double mdv) const;

  /**
     \brief writes the matrix to a file using the Matrix Market file format
  */
  RIPLEY_DLL_API
  virtual void saveMM(const std::string& fileName) const;

  /**
     \brief writes the matrix to a file using the Harwell-Boeing file format
  */
  RIPLEY_DLL_API
  virtual void saveHB(const std::string& fileName) const;

  /**
     \brief sets the matrix entries to zero
  */
  RIPLEY_DLL_API
  virtual void resetValues() const;

  /**
     \brief maps escript options onto Paso options:
  */
  RIPLEY_DLL_API
  static int mapOptionToPaso(const int option);

  /**
     \brief extract paso options from SolutionOptions class
  */
 
  RIPLEY_DLL_API
  static void escriptToPasoOptions(Paso_Options* paso_options, const boost::python::object& options);

  /**
     \brief copied diagonistic data back to the solver option.
  */
 
  RIPLEY_DLL_API
  static void pasoToEscriptOptions(const Paso_Options* paso_options,boost::python::object& options);
 
  /**
     \brief prints information about a system matrix
  */
  RIPLEY_DLL_API
  void Print_Matrix_Info(const bool) const;

 private:

   /**
      \brief
      solves the linear system this*out=in
   */
   RIPLEY_DLL_API
   virtual void setToSolution(escript::Data& out, escript::Data& in, boost::python::object& options) const;

   /**
       \brief
       performs y+=this*x
   */
   RIPLEY_DLL_API
   virtual void ypAx(escript::Data& y, escript::Data& x) const;

   //
   // pointer to the externally created ripley mesh - system_matrix.
   //
   boost::shared_ptr<Paso_SystemMatrix> m_system_matrix;
};

} // end of namespace ripley

#endif // __RIPLEY_SYSTEMMATRIXADAPTER_H__

