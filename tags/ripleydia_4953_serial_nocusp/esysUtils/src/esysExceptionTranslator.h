
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


#if !defined  esysUtils_esysExceptionTranslator_20040419_H
#define esysUtils_esysExceptionTranslator_20040419_H
#include "system_dep.h"

#include "EsysException.h"

#include "boost/python/errors.hpp"

#include <iostream>

namespace esysUtils {
  /**
     \brief
     Function which translates an EsysException into a python RuntimeError
  */
  ESYSUTILS_DLL_API
  void RuntimeErrorTranslator(EsysException const& e);

  /**
     \brief
     Function which translates an EsysException into a python ValueError
  */
  ESYSUTILS_DLL_API
  void ValueErrorTranslator(EsysException const& e);
} // end of namespace

#endif
