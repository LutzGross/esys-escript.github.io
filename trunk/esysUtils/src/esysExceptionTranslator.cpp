
/*****************************************************************************
*
* Copyright (c) 2003-2012 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development since 2012 by School of Earth Sciences
*
*****************************************************************************/


#include "system_dep.h"
#include "esysExceptionTranslator.h" 
#include <iostream>

using namespace std;

namespace esysUtils {

void esysExceptionTranslator(EsysException const& e) 
  {
    PyErr_SetString(PyExc_RuntimeError,e.what());
  }

}  // end of namespace
