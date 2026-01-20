
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/

#ifndef __ESCRIPT_EXCEPTIONTRANSLATORS_H__
#define __ESCRIPT_EXCEPTIONTRANSLATORS_H__

#include "system_dep.h"
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
  ESCRIPT_DLL_API
  void AssertionErrorTranslator(const EsysException& e);

  /**
     \brief
     Function which translates an EsysException into a python IOError
  */
  ESCRIPT_DLL_API
  void IOErrorTranslator(const EsysException& e);

  /**
     \brief
     Function which translates an EsysException into a python NotImplementedError
  */
  ESCRIPT_DLL_API
  void NotImplementedErrorTranslator(const EsysException& e);

  /**
     \brief
     Function which translates an EsysException into a python RuntimeError
  */
  ESCRIPT_DLL_API
  void RuntimeErrorTranslator(const EsysException& e);

  /**
     \brief
     Function which translates an EsysException into a python ValueError
  */
  ESCRIPT_DLL_API
  void ValueErrorTranslator(const EsysException& e);

} // end of namespace

#endif // __ESCRIPT_EXCEPTIONTRANSLATORS_H__

