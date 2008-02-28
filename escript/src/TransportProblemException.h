
/* $Id$ */

/*******************************************************
*
*       Copyright 2007 by University of Queensland
*
*                http://esscc.uq.edu.au
*        Primary Business: Queensland, Australia
*  Licensed under the Open Software License version 3.0
*     http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/

#if !defined  escript_TransportProblemException_20040608_H
#define escript_TransportProblemException_20040608_H

#include "system_dep.h"
#include "esysUtils/EsysException.h"

namespace escript
{

  /**
  \brief
  TransportProblemException exception class.

  Description:
  TransportProblemException exception class.
  The class provides a public function returning the exception name
  */
  class TransportProblemException : public esysUtils::EsysException {

  protected:

     typedef EsysException Parent;

  public:
    /**
    \brief
    Default constructor for the exception.
    */
    ESCRIPT_DLL_API
    TransportProblemException() : Parent() {}
    /**
    \brief
    Constructor for the exception.
    */
    ESCRIPT_DLL_API
    TransportProblemException(const char *cstr) : Parent(cstr) {}
    /**
    \brief
    Constructor for the exception.
    */
    ESCRIPT_DLL_API
    TransportProblemException(const std::string &str) : Parent(str) {}

    /// Destructor
    ESCRIPT_DLL_API
    virtual ~TransportProblemException() THROW() {}
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
