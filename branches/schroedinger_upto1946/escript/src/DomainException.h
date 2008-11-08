
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


#if !defined  escript_DomainException_20040608_H
#define escript_DomainException_20040608_H
#include "system_dep.h"

#include "esysUtils/EsysException.h"

namespace escript
{

  /**
  \brief
  DomainException exception class.

  Description:
  DomainException exception class.
  The class provides a public function returning the exception name
  */
  class DomainException : public esysUtils::EsysException
  {

  protected:

    typedef EsysException Parent;

  public:
    /**
    \brief
    Default constructor for the exception.
    */
    ESCRIPT_DLL_API
    DomainException() : Parent() { updateMessage();}
    /**
    \brief
    Constructor for the exception.
    */
    ESCRIPT_DLL_API
    DomainException(const char *cstr) : Parent(cstr) { updateMessage();}
    /**
    \brief
    Constructor for the exception.
    */
    ESCRIPT_DLL_API
    DomainException(const std::string &str) : Parent(str) { updateMessage();}
    /**
    \brief
    Copy Constructor for the exception.
    */
    ESCRIPT_DLL_API
    DomainException(const DomainException &other) : Parent(other)
      {
        updateMessage();
      }

    ESCRIPT_DLL_API
    inline virtual DomainException &
    operator=(const DomainException &other ) THROW(/**/)
       {
         Parent::operator=(other);
         updateMessage();
         return *this;
       }

    /// Destructor
    ESCRIPT_DLL_API
    virtual ~DomainException() THROW(/**/) {}
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
