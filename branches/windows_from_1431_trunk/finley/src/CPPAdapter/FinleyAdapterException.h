
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

#if !defined  finley_FinleyAdapterException_20040526_H
#define finley_FinleyAdapterException_20040526_H
#include "system_dep.h"

#include "esysUtils/EsysException.h"

namespace finley
{

  /**
  \brief
  FinleyAdapterException exception class.

  Description:
  FinleyAdapterException exception class.
  The class provides a public function returning the exception name
  */
  class FinleyAdapterException : public esysUtils::EsysException
  {

  public:
    /**
    \brief
    Default constructor for the exception.
    */
    FINLEY_DLL_API
    FinleyAdapterException() : EsysException() {}
    /**
    \brief
    Constructor for the exception.
    */
    FINLEY_DLL_API
    FinleyAdapterException(const char *cstr) : EsysException(cstr) {}
    /**
    \brief
    Constructor for the exception.
    */
    FINLEY_DLL_API
    FinleyAdapterException(const std::string &str) : EsysException(str) {}

    /// Destructor
    FINLEY_DLL_API
    virtual ~FinleyAdapterException() throw() {}
    /**
    \brief
    Returns the name of the exception.
    */
    FINLEY_DLL_API
    virtual const std::string & exceptionName() const;

  private:

    //
    // the exception name is immutable and class-wide.
    // Inheritor note; you need one of these too.
    // and an overloaded exceptionName() in your .cpp implementation file. 
    static const std::string exceptionNameValue;
  };

} // end of namespace
#endif
