
/*******************************************************
*
* Copyright (c) 2003-2009 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


#if !defined escript_EsysAssertException_20040330_H
#define escript_EsysAssertException_20040330_H
#include "system_dep.h"

#include "EsysException.h"


namespace esysUtils {

  /**
  \brief
  EsysAssertException exception class.

  Description:
  EsysAssertException exception class.
  The class provides a public function returning the exception name.
  */
  class EsysAssertException : public EsysException
  {

  protected:

     typedef EsysException Parent;

  public:
    /**
    \brief
    Default constructor for the exception.
    */
    ESYSUTILS_DLL_API
    inline
    EsysAssertException() : Parent() { updateMessage();}
    /**
    \brief
    Constructor for the exception.
    */
    ESYSUTILS_DLL_API
    inline
    EsysAssertException(const char *cstr) : Parent(cstr) { updateMessage();}
    /**
    \brief
    Constructor for the exception.
    */
    ESYSUTILS_DLL_API
    inline
    EsysAssertException(const std::string &str) :
    Parent(str) { updateMessage();}
    /**
    \brief
    Copy Constructor for the exception.
    */
    ESYSUTILS_DLL_API
    inline
    EsysAssertException(const EsysAssertException &other) : Parent(other)
      {
        updateMessage();
      }

    ESYSUTILS_DLL_API
    inline 
    EsysAssertException &
    operator=(const EsysAssertException &other) THROW(NO_ARG)
        {
           Parent::operator=(other);
           updateMessage();   
           return *this;
        }


    /// Destructor
    ESYSUTILS_DLL_API
    virtual ~EsysAssertException() THROW(NO_ARG) {}

    /**
    \brief
    Returns the name of the exception.
    */
    ESYSUTILS_DLL_API
    virtual const std::string & exceptionName() const;

    /**
    \brief
    Builds a formatted message and throws an EsysAssertException.
    */
    ESYSUTILS_DLL_API
    static void assertFailure (const std::string& assertion,
                               const std::string& date,
                               const std::string& file,
                               int line, const std::string& errDesc);
  private:

    //
    // the exception name is immutable and class-wide.
    // Inheritor note; you need one of these too.
    // and an overloaded exceptionName() in your .cpp implementation file. 
    static const std::string exceptionNameValue;
  };

} // end of namespace
 
#endif
