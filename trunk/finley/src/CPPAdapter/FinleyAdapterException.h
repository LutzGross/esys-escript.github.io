
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

  protected:

     typedef EsysException Parent;

  public:
    /**
    \brief
    Default constructor for the exception.
    */
    FINLEY_DLL_API
    FinleyAdapterException() : Parent() { updateMessage();}
    /**
    \brief
    Constructor for the exception.
    */
    FINLEY_DLL_API
    FinleyAdapterException(const char *cstr) : Parent(cstr) { updateMessage();}
    /**
    \brief
    Constructor for the exception.
    */
    FINLEY_DLL_API
    FinleyAdapterException(const std::string &str) :
    Parent(str) { updateMessage();}
    /**
    \brief
    Copy Constructor for the exception.
    */
    FINLEY_DLL_API
    FinleyAdapterException(const FinleyAdapterException &other) : Parent(other)
      {
        updateMessage();
      }

    /// Destructor
    FINLEY_DLL_API
    virtual ~FinleyAdapterException() THROW(NO_ARG) {}

    /**
    \brief
    Assignment operator.
    */
    FINLEY_DLL_API
    inline FinleyAdapterException &
    operator=(const FinleyAdapterException &other ) THROW(NO_ARG)
       {
         Parent::operator=(other);
         updateMessage();
         return *this;
       }

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
