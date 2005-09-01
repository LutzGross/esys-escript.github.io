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
#include "finley/finleyC/System.h"
}
#include "escript/Data/Data.h"
#include "finley/CPPAdapter/SystemMatrixAdapter.h" 
#include "finley/CPPAdapter/FinleyAdapterException.h" 
#include "finley/CPPAdapter/FinleyError.h"
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

SystemMatrixAdapter::SystemMatrixAdapter(Finley_SystemMatrix* system_matrix,
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
        Finley_SystemMatrix* mat=m_system_matrix.get();
        Finley_SystemMatrix_dealloc(mat);
    }
}

Finley_SystemMatrix* SystemMatrixAdapter::getFinley_SystemMatrix() const 
{
   return m_system_matrix.get();
}

void SystemMatrixAdapter::ypAx(escript::Data& y, const escript::Data& x) const 
{
   Finley_SystemMatrix* mat=getFinley_SystemMatrix();
   Finley_SystemMatrixVector(&(y.getDataC()),mat,&(x.getDataC()));
   checkFinleyError();
}

void SystemMatrixAdapter::setToSolution(escript::Data& out, const escript::Data& in, const boost::python::dict& options) const
{
    Finley_SolverOptions finley_options;
    Finley_SystemMatrix_setDefaults(&finley_options);
    // extract options 
    #define EXTRACT(__key__,__val__,__type__) if ( options.has_key(__key__)) finley_options.__val__=boost::python::extract<__type__>(options.get(__key__))
    EXTRACT("verbose",verbose,int);
    EXTRACT("reordering",reordering,int);
    EXTRACT(ESCRIPT_TOLERANCE_KEY,tolerance,double);
    EXTRACT(ESCRIPT_METHOD_KEY,method,int);
    EXTRACT(ESCRIPT_SYMMETRY_KEY,symmetric,int);
    EXTRACT("preconditioner",preconditioner,int);
    EXTRACT("iter_max",iter_max,int);
    EXTRACT("drop_tolerance",drop_tolerance,double);
    EXTRACT("drop_storage",drop_storage,double);
    EXTRACT("truncation",truncation,double);
    EXTRACT("restart",restart,double);
    #undef EXTRACT
    Finley_SystemMatrix_solve(getFinley_SystemMatrix(),&(out.getDataC()),&(in.getDataC()),&finley_options);
    checkFinleyError();
}

void SystemMatrixAdapter::nullifyRowsAndCols(const escript::Data& row_q, const escript::Data& col_q, const double mdv) const
{
    Finley_SystemMatrix* system_matrix_ptr = getFinley_SystemMatrix();
    escriptDataC row_qC = row_q.getDataC();
    escriptDataC col_qC = col_q.getDataC();

    Finley_SystemMatrix_nullifyRowsAndCols(system_matrix_ptr, &row_qC, &col_qC, mdv);
}

void SystemMatrixAdapter::saveMM(const std::string& fileName) const
{
    char fName[fileName.size()+1];
    strcpy(fName,fileName.c_str());
    Finley_SystemMatrix* system_matrix_ptr = getFinley_SystemMatrix();
    Finley_SystemMatrix_saveMM(system_matrix_ptr,fName);
    checkFinleyError();
}

void SystemMatrixAdapter::saveHB(const std::string& fileName) const
{
    char fName[fileName.size()+1];
    strcpy(fName,fileName.c_str());
    Finley_SystemMatrix* system_matrix_ptr = getFinley_SystemMatrix();
    Finley_SystemMatrix_saveHB(system_matrix_ptr,fName);
    checkFinleyError();
}

void SystemMatrixAdapter::resetValues() const
{
   Finley_SystemMatrix* system_matrix_ptr = getFinley_SystemMatrix();
   Finley_SystemMatrix_setValues(system_matrix_ptr,0.);
   Finley_SystemMatrix_solve_free(system_matrix_ptr);
   checkFinleyError();
}

}  // end of namespace
