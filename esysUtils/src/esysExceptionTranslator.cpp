
/* $Id$ */

/*******************************************************
 *
 *           Copyright 2003-2007 by ACceSS MNRF
 *       Copyright 2007 by University of Queensland
 *
 *                http://esscc.uq.edu.au
 *        Primary Business: Queensland, Australia
 *  Licensed under the Open Software License version 3.0
 *     http://www.opensource.org/licenses/osl-3.0.php
 *
 *******************************************************/

#include "system_dep.h"
#include "esysExceptionTranslator.h" 
#include <iostream>

using namespace std;

namespace esysUtils {

void esysExceptionTranslator(const EsysException & e) 
  {
	  cout << "translated!\n";
    PyErr_SetString(PyExc_RuntimeError,e.what());
  }

}  // end of namespace
