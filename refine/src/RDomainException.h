
/*******************************************************
*
* Copyright (c) 2011 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


#if !defined  rdomainException_H
#define rdomainException_H
#include "system_dep.h"

#include "esysUtils/EsysException.h"

namespace refine
{

  /**
  \brief
  RDomainException exception class.

  Description:
  RDomainException exception class.
  The class provides a public function returning the exception name
  */
  class RDomainException : public esysUtils::EsysException
  {

  protected:

     typedef EsysException Parent;

  public:
    /**
    \brief
    Default constructor for the exception.
    */
    REFINE_DLL_API
    RDomainException() : Parent() { updateMessage();}
    /**
    \brief
    Constructor for the exception.
    */
    REFINE_DLL_API
    RDomainException(const char *cstr) : Parent(cstr) { updateMessage();}
    /**
    \brief
    Constructor for the exception.
    */
    REFINE_DLL_API
    RDomainException(const std::string &str) :
    Parent(str) { updateMessage();}
    /**
    \brief
    Copy Constructor for the exception.
    */
    REFINE_DLL_API
    RDomainException(const RDomainException &other) : Parent(other)
      {
        updateMessage();
      }

    /// Destructor
    REFINE_DLL_API
    virtual ~RDomainException() THROW(NO_ARG) {}

    /**
    \brief
    Assignment operator.
    */
    REFINE_DLL_API
    inline RDomainException &
    operator=(const RDomainException &other ) THROW(NO_ARG)
       {
         Parent::operator=(other);
         updateMessage();
         return *this;
       }

    /**
    \brief
    Returns the name of the exception.
    */
    REFINE_DLL_API
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
