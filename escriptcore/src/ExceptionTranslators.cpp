
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

#include "ExceptionTranslators.h" 

namespace escript {

void AssertionErrorTranslator(const EsysException& e)
{
    PyErr_SetString(PyExc_AssertionError, e.what());
}

void IOErrorTranslator(const EsysException& e) 
{
    PyErr_SetString(PyExc_IOError, e.what());
}

void NotImplementedErrorTranslator(const EsysException& e) 
{
    PyErr_SetString(PyExc_NotImplementedError, e.what());
}

void RuntimeErrorTranslator(const EsysException& e) 
{
    PyErr_SetString(PyExc_RuntimeError, e.what());
}

void ValueErrorTranslator(const EsysException& e) 
{
    PyErr_SetString(PyExc_ValueError, e.what());
}

}  // end of namespace

