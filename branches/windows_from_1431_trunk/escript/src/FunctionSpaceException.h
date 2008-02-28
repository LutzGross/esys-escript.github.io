
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

#if !defined  escript_FunctionSpaceException_20040602_H
#define escript_FunctionSpaceException_20040602_H
#include "system_dep.h"

#include "esysUtils/EsysException.h"

namespace escript
{

  /**
  \brief
  FunctionSpaceException exception class.

  Description:
  FunctionSpaceException exception class.
  The class provides a public function returning the exception name
  */
  class FunctionSpaceException : public esysUtils::EsysException
  {

  protected:

     typedef EsysException Parent;

  public:
    /**
    \brief
    Default constructor for the exception.
    */
    ESCRIPT_DLL_API
    FunctionSpaceException() : Parent() {}
    /**
    \brief
    Constructor for the exception.
    */
    ESCRIPT_DLL_API
    FunctionSpaceException(const char *cstr) : Parent(cstr) {}
    /**
    \brief
    Constructor for the exception.
    */
    ESCRIPT_DLL_API
    FunctionSpaceException(const std::string &str) : Parent(str) {}

    /// Destructor
    ESCRIPT_DLL_API
    virtual ~FunctionSpaceException() THROW() {}
    /**
    \brief
    Returns the name of the exception.
    */
    ESCRIPT_DLL_API
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
