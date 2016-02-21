
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

#ifndef __ESCRIPT_EXCEPTIONTRANSLATORS_H__
#define __ESCRIPT_EXCEPTIONTRANSLATORS_H__

#include "esysUtils/EsysException.h"

namespace escript {
  /**
     \brief
     Function which translates an EsysException into a python RuntimeError
  */
  void RuntimeErrorTranslator(const esysUtils::EsysException& e);

  /**
     \brief
     Function which translates an EsysException into a python ValueError
  */
  void ValueErrorTranslator(const esysUtils::EsysException& e);
} // end of namespace

#endif // __ESCRIPT_EXCEPTIONTRANSLATORS_H__

