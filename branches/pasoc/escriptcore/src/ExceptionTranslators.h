
/*****************************************************************************
*
* Copyright (c) 2003-2017 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/

#ifndef __ESCRIPT_EXCEPTIONTRANSLATORS_H__
#define __ESCRIPT_EXCEPTIONTRANSLATORS_H__

#include "DataTypes.h"
#include "EsysException.h"

// put this within all boost module definitions so exceptions are translated
// properly
#define REGISTER_ESCRIPT_EXCEPTION_TRANSLATORS \
    register_exception_translator<escript::AssertException>(&escript::AssertionErrorTranslator);\
    register_exception_translator<escript::IOError>(&escript::IOErrorTranslator);\
    register_exception_translator<escript::NotImplementedError>(&escript::NotImplementedErrorTranslator);\
    register_exception_translator<escript::ValueError>(&escript::ValueErrorTranslator)

namespace escript {

  /**
     \brief
     Function which translates an EsysException into a python AssertionError
  */
  void AssertionErrorTranslator(const EsysException& e);

  /**
     \brief
     Function which translates an EsysException into a python IOError
  */
  void IOErrorTranslator(const EsysException& e);

  /**
     \brief
     Function which translates an EsysException into a python NotImplementedError
  */
  void NotImplementedErrorTranslator(const EsysException& e);

  /**
     \brief
     Function which translates an EsysException into a python RuntimeError
  */
  void RuntimeErrorTranslator(const EsysException& e);

  /**
     \brief
     Function which translates an EsysException into a python ValueError
  */
  void ValueErrorTranslator(const EsysException& e);

} // end of namespace

#endif // __ESCRIPT_EXCEPTIONTRANSLATORS_H__

