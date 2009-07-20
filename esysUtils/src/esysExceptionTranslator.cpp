
/*******************************************************
*
* Copyright (c) 2003-2009 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


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
