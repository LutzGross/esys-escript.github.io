
/*****************************************************************************
*
* Copyright (c) 2003-2016 by The University of Queensland
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

