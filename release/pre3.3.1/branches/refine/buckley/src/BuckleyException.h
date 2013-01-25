
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

namespace buckley
{

  /**
  \brief
  BuckleyException exception class.

  Description:
  BuckleyException exception class.
  The class provides a public function returning the exception name
  */
  class BuckleyException : public esysUtils::EsysException
  {

  protected:

     typedef EsysException Parent;

  public:
    /**
    \brief
    Default constructor for the exception.
    */
    BUCKLEY_DLL_API
    BuckleyException() : Parent() { updateMessage();}
    /**
    \brief
    Constructor for the exception.
    */
    BUCKLEY_DLL_API
    BuckleyException(const char *cstr) : Parent(cstr) { updateMessage();}
    /**
    \brief
    Constructor for the exception.
    */
    BUCKLEY_DLL_API
    BuckleyException(const std::string &str) :
    Parent(str) { updateMessage();}
    /**
    \brief
    Copy Constructor for the exception.
    */
    BUCKLEY_DLL_API
    BuckleyException(const BuckleyException &other) : Parent(other)
      {
        updateMessage();
      }

    /// Destructor
    BUCKLEY_DLL_API
    virtual ~BuckleyException() THROW(NO_ARG) {}

    /**
    \brief
    Assignment operator.
    */
    BUCKLEY_DLL_API
    inline BuckleyException &
    operator=(const BuckleyException &other ) THROW(NO_ARG)
       {
         Parent::operator=(other);
         updateMessage();
         return *this;
       }

    /**
    \brief
    Returns the name of the exception.
    */
    BUCKLEY_DLL_API
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
