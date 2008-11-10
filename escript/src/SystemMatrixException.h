
/*******************************************************
*
* Copyright (c) 2003-2008 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


#if !defined  escript_SystemMatrixException_20040608_H
#define escript_SystemMatrixException_20040608_H

#include "system_dep.h"
#include "esysUtils/EsysException.h"

namespace escript
{

  /**
  \brief
  SystemMatrixException exception class.

  Description:
  SystemMatrixException exception class.
  The class provides a public function returning the exception name
  */
  class SystemMatrixException : public esysUtils::EsysException
  {

  protected:

     typedef EsysException Parent;

  public:
    /**
    \brief
    Default constructor for the exception.
    */
    ESCRIPT_DLL_API
    SystemMatrixException() : Parent() { updateMessage();}
    /**
    \brief
    Constructor for the exception.
    */
    ESCRIPT_DLL_API
    SystemMatrixException(const char *cstr) : Parent(cstr) { updateMessage();}
    /**
    \brief
    Constructor for the exception.
    */
    ESCRIPT_DLL_API
    SystemMatrixException(const std::string &str) :
    Parent(str) { updateMessage();}
    /**
    \brief
    Copy Constructor for the exception.
    */
    ESCRIPT_DLL_API
    SystemMatrixException(const SystemMatrixException &other) : Parent(other)
      {
        updateMessage();
      }

    ESCRIPT_DLL_API
    inline SystemMatrixException &
    operator=(const SystemMatrixException &other ) THROW(NO_ARG)
       {
         Parent::operator=(other);
         updateMessage();
         return *this;
       }

    /// Destructor
    ESCRIPT_DLL_API
    virtual ~SystemMatrixException() THROW(NO_ARG) {}
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
