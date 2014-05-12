
/*****************************************************************************
*
* Copyright (c) 2003-2014 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/


#if !defined  SystemMatrixAdapter_20040610_H
#define SystemMatrixAdapter_20040610_H
#include "system_dep.h"

#include "paso/SystemMatrix.h"
#include "paso/Options.h"

#include "PasoException.h"

#include "escript/AbstractSystemMatrix.h"
#include "escript/Data.h"
#include "escript/UtilC.h"   

#include <boost/python/object.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/python/extract.hpp>

namespace paso {

PASOWRAP_DLL_API
class SystemMatrixAdapter: public escript::AbstractSystemMatrix
{

/**
   \brief
   Wrapper for paso::SystemMatrix. 

   Description:
   Wrapper for paso::SystemMatrix.
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
  SystemMatrixAdapter(SystemMatrix_ptr system_matrix,
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
  SystemMatrix_ptr getPaso_SystemMatrix() const;

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
    nullifyRowsAndCols - calls SystemMatrix::nullifyRowsAndCols.
  */
  void nullifyRowsAndCols(escript::Data& row_q, escript::Data& col_q, const double mdv) const;

  /**
     \brief writes the matrix to a file using the Matrix Market file format
  */
  virtual void saveMM(const std::string& filename) const;

  /**
     \brief writes the matrix to a file using the Harwell-Boeing file format
  */
  virtual void saveHB(const std::string& filename) const;

  /**
     \brief sets the matrix entries to zero
  */
  virtual void resetValues() const;

  /**
     \brief maps escript options onto Paso options
  */
  static int mapOptionToPaso(const int option);

  /**
     \brief returns the identifier of the matrix type to be used for the global
     stiffness matrix when the given solver, preconditioner and package is
     used
  */
  static int getSystemMatrixTypeId(const int solver, const int preconditioner, const int package, const bool symmetry, const esysUtils::JMPI& mpiInfo);

  /**
     \brief extract paso options from SolutionOptions class
  */
 
  static void escriptToPasoOptions(Options* paso_options, const boost::python::object& options);

  /**
     \brief copies diagnostic data back to the solver options
  */
 
  static void pasoToEscriptOptions(const Options* paso_options, boost::python::object& options);
 
  /**
     \brief prints information about a system matrix
  */
  void Print_Matrix_Info(const bool) const;

 protected:

 private:

   /**
      \brief
      solves the linear system this*out=in
   */
   virtual void setToSolution(escript::Data& out, escript::Data& in, boost::python::object& options) const;

   /**
       \brief
       performs y+=this*x
   */
   virtual void ypAx(escript::Data& y, escript::Data& x) const;

   //
   // pointer to the externally created system_matrix.
   //
   SystemMatrix_ptr m_system_matrix;

};

} // end of namespace
#endif
