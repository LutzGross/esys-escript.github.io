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

SystemMatrixAdapter::SystemMatrixAdapter()
{
   throw FinleyAdapterException("Error - Illegal to generate default SystemMatrixAdapter.");
}

SystemMatrixAdapter::SystemMatrixAdapter(const Finley_SystemMatrix* system_matrix,
                                         const int row_blocksize,
                                         const escript::FunctionSpace& row_functionspace,
                                         const int column_blocksize,
                                         const escript::FunctionSpace& column_functionspace):
AbstractSystemMatrix(row_blocksize,row_functionspace,column_blocksize,column_functionspace),
m_system_matrix(system_matrix)
{
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
    EXTRACT("tolerance",tolerance,double);
    EXTRACT("iterative_method",iterative_method,int);
    EXTRACT("preconditioner",preconditioner,int);
    EXTRACT("iter_max",iter_max,int);
    EXTRACT("drop_tolerance",drop_tolerance,double);
    EXTRACT("drop_storage",drop_storage,double);
    EXTRACT("iterative",iterative,int);
    #undef EXTRACT
    if (finley_options.iterative) {
       Finley_SystemMatrix_iterative(getFinley_SystemMatrix(),&(out.getDataC()),&(in.getDataC()),&finley_options);
    } else {
       Finley_SystemMatrix_solve(getFinley_SystemMatrix(),&(out.getDataC()),&(in.getDataC()),&finley_options);
    }
    checkFinleyError();
}

void SystemMatrixAdapter::nullifyRowsAndCols(const escript::Data& row_q, const escript::Data& col_q, const double mdv) const
{
    Finley_SystemMatrix* system_matrix_ptr = getFinley_SystemMatrix();
    escriptDataC row_qC = row_q.getDataC();
    escriptDataC col_qC = col_q.getDataC();

    Finley_SystemMatrix_nullifyRowsAndCols(system_matrix_ptr, &row_qC, &col_qC, mdv);
}

}  // end of namespace
