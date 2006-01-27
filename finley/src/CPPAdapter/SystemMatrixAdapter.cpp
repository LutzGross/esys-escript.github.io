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
 $Id$
*/
extern "C" {
#include "SystemMatrix.h"
#include "Options.h"
}
#include "escript/Data/Data.h"
#include "escript/Data/UtilC.h"
#include "SystemMatrixAdapter.h" 
#include "FinleyAdapterException.h" 
#include "FinleyError.h"
#include <boost/python/extract.hpp>

using namespace std;

namespace finley {

struct null_deleter
{
  void operator()(void const *ptr) const
  {
  }
};


SystemMatrixAdapter::SystemMatrixAdapter()
{
   throw FinleyAdapterException("Error - Illegal to generate default SystemMatrixAdapter.");
}

SystemMatrixAdapter::SystemMatrixAdapter(Paso_SystemMatrix* system_matrix,
                                         const int row_blocksize,
                                         const escript::FunctionSpace& row_functionspace,
                                         const int column_blocksize,
                                         const escript::FunctionSpace& column_functionspace):
AbstractSystemMatrix(row_blocksize,row_functionspace,column_blocksize,column_functionspace)
{
    m_system_matrix.reset(system_matrix,null_deleter());
}

SystemMatrixAdapter::~SystemMatrixAdapter()
{ 
    if (m_system_matrix.unique()) {
        Paso_SystemMatrix* mat=m_system_matrix.get();
        Paso_SystemMatrix_dealloc(mat);
    }
}

Paso_SystemMatrix* SystemMatrixAdapter::getPaso_SystemMatrix() const 
{
   return m_system_matrix.get();
}

void SystemMatrixAdapter::ypAx(escript::Data& y,escript::Data& x) const 
{
   Paso_SystemMatrix* mat=getPaso_SystemMatrix();

  if ( x.getDataPointSize()  != getColumnBlockSize()) {
   throw FinleyAdapterException("matrix vector product : column block size does not match the number of components in input.");
  } else if (y.getDataPointSize() != getRowBlockSize()) {
   throw FinleyAdapterException("matrix vector product : row block size does not match the number of components in output.");
  } else if ( x.getFunctionSpace()  != getColumnFunctionSpace()) {
   throw FinleyAdapterException("matrix vector product : column function space and function space of input don't match.");
  } else if (y.getFunctionSpace() != getRowFunctionSpace()) {
   throw FinleyAdapterException("matrix vector product : row function space and function space of output don't match.");
  }
  x.expand();
  y.expand();
  double* x_dp=x.getSampleData(0);
  double* y_dp=y.getSampleData(0);
  Paso_SystemMatrix_MatrixVector(1., mat,x_dp, 1.,y_dp);
  checkPasoError();
}

int SystemMatrixAdapter::mapOptionToPaso(const int option)  {
   switch (option) {
       case  ESCRIPT_DEFAULT:
          return PASO_DEFAULT;
       case  ESCRIPT_DIRECT:
          return PASO_DIRECT;
       case  ESCRIPT_CHOLEVSKY:
          return PASO_CHOLEVSKY;
       case  ESCRIPT_PCG:
          return PASO_PCG;
       case  ESCRIPT_CR:
          return PASO_CR;
       case  ESCRIPT_CGS:
          return PASO_CGS;
       case  ESCRIPT_BICGSTAB:
          return PASO_BICGSTAB;
       case  ESCRIPT_SSOR:
          return PASO_SSOR;
       case  ESCRIPT_ILU0:
          return PASO_ILU0;
       case  ESCRIPT_ILUT:
          return PASO_ILUT;
       case  ESCRIPT_JACOBI:
          return PASO_JACOBI;
       case  ESCRIPT_GMRES:
          return PASO_GMRES;
       case  ESCRIPT_PRES20:
          return PASO_PRES20;
       case  ESCRIPT_NO_REORDERING:
          return PASO_NO_REORDERING;
       case  ESCRIPT_MINIMUM_FILL_IN:
          return PASO_MINIMUM_FILL_IN;
       case  ESCRIPT_NESTED_DISSECTION:
          return PASO_NESTED_DISSECTION;
       case  ESCRIPT_SCSL:
          return PASO_SCSL;
       case  ESCRIPT_MKL:
          return PASO_MKL;
       case  ESCRIPT_UMFPACK:
          return PASO_UMFPACK;
       case  ESCRIPT_ITERATIVE:
          return PASO_ITERATIVE;
       case  ESCRIPT_PASO:
          return PASO_PASO;
       case  ESCRIPT_LUMPING:
          return PASO_LUMPING;
       case  ESCRIPT_AMG:
          return PASO_AMG;
       case  ESCRIPT_RILU:
          return PASO_RILU;
       default:
           stringstream temp;
           temp << "Error - Cannot map option value "<< option << " onto Paso";
           throw FinleyAdapterException(temp.str());
    }
}

void SystemMatrixAdapter::setToSolution(escript::Data& out,escript::Data& in, const boost::python::dict& options) const
{
    Paso_SystemMatrix* mat=getPaso_SystemMatrix();
    Paso_Options paso_options;
    Paso_Options_setDefaults(&paso_options);
    // extract options 
    #define EXTRACT(__key__,__val__,__type__) if ( options.has_key(__key__)) paso_options.__val__=boost::python::extract<__type__>(options.get(__key__))
    #define EXTRACT_OPTION(__key__,__val__,__type__) if ( options.has_key(__key__)) paso_options.__val__=mapOptionToPaso(boost::python::extract<__type__>(options.get(__key__)))
    EXTRACT("verbose",verbose,int);
    EXTRACT_OPTION("reordering",reordering,int);
    EXTRACT(ESCRIPT_TOLERANCE_KEY,tolerance,double);
    EXTRACT_OPTION(ESCRIPT_METHOD_KEY,method,int);
    EXTRACT(ESCRIPT_SYMMETRY_KEY,symmetric,int);
    EXTRACT_OPTION(ESCRIPT_PACKAGE_KEY,package,int);
    EXTRACT_OPTION("preconditioner",preconditioner,int);
    EXTRACT("iter_max",iter_max,int);
    EXTRACT("drop_tolerance",drop_tolerance,double);
    EXTRACT("drop_storage",drop_storage,double);
    EXTRACT("truncation",truncation,double);
    EXTRACT("restart",restart,double);
    #undef EXTRACT
    #undef EXTRACT_OPTION
    if ( out.getDataPointSize()  != getColumnBlockSize()) {
     throw FinleyAdapterException("solve : column block size does not match the number of components of solution.");
    } else if ( in.getDataPointSize() != getRowBlockSize()) {
     throw FinleyAdapterException("solve : row block size does not match the number of components of  right hand side.");
    } else if ( out.getFunctionSpace()  != getColumnFunctionSpace()) {
     throw FinleyAdapterException("solve : column function space and function space of solution don't match.");
    } else if (in.getFunctionSpace() != getRowFunctionSpace()) {
     throw FinleyAdapterException("solve : row function space and function space of right hand side don't match.");
    }
    out.expand();
    in.expand();
    double* out_dp=out.getSampleData(0);
    double* in_dp=in.getSampleData(0);
    Paso_solve(mat,out_dp,in_dp,&paso_options);
    checkPasoError();
}

void SystemMatrixAdapter::nullifyRowsAndCols(escript::Data& row_q,escript::Data& col_q, const double mdv) const
{
    Paso_SystemMatrix* mat = getPaso_SystemMatrix();
    if ( col_q.getDataPointSize()  != getColumnBlockSize()) {
     throw FinleyAdapterException("nullifyRowsAndCols : column block size does not match the number of components of column mask.");
    } else if ( row_q.getDataPointSize() != getRowBlockSize()) {
     throw FinleyAdapterException("nullifyRowsAndCols : row block size does not match the number of components of row mask.");
    } else if ( col_q.getFunctionSpace()  != getColumnFunctionSpace()) {
     throw FinleyAdapterException("nullifyRowsAndCols : column function space and function space of column mask don't match.");
    } else if (row_q.getFunctionSpace() != getRowFunctionSpace()) {
     throw FinleyAdapterException("nullifyRowsAndCols : row function space and function space of row mask don't match.");
    }
    row_q.expand();
    col_q.expand();
    double* row_q_dp=row_q.getSampleData(0);
    double* col_q_dp=col_q.getSampleData(0);
    Paso_SystemMatrix_nullifyRowsAndCols(mat,row_q_dp,col_q_dp, mdv);
    checkPasoError();
}

void SystemMatrixAdapter::saveMM(const std::string& fileName) const
{
    char fName[fileName.size()+1];
    strcpy(fName,fileName.c_str());
    Paso_SystemMatrix* mat = getPaso_SystemMatrix();
    Paso_SystemMatrix_saveMM(mat,fName);
    checkPasoError();
}

void SystemMatrixAdapter::saveHB(const std::string& fileName) const
{
    char fName[fileName.size()+1];
    strcpy(fName,fileName.c_str());
    Paso_SystemMatrix* mat = getPaso_SystemMatrix();
    Paso_SystemMatrix_saveHB(mat,fName);
    checkPasoError();
}

void SystemMatrixAdapter::resetValues() const
{
   Paso_SystemMatrix* mat = getPaso_SystemMatrix();
   Paso_SystemMatrix_setValues(mat,0.);
   Paso_solve_free(mat);
   checkPasoError();
}

}  // end of namespace
