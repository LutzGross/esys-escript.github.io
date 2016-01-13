
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


#if !defined  esysUtils_esysExceptionTranslator_20040419_H
#define esysUtils_esysExceptionTranslator_20040419_H
#include "system_dep.h"

#include "EsysException.h"

#include "boost/python/errors.hpp"

#include <iostream>

namespace esysUtils {
  /**
     \brief
     Function which translates an EsysException into a python exception
  */
  ESYSUTILS_DLL_API
  void esysExceptionTranslator(EsysException const& e);
} // end of namespace

#endif
