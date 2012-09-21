
/*****************************************************************************
*
* Copyright (c) 2003-2012 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development since 2012 by School of Earth Sciences
*
*****************************************************************************/


#if !defined  dudley_DudleyAdapterException_20040526_H
#define dudley_DudleyAdapterException_20040526_H
#include "system_dep.h"

#include "esysUtils/EsysException.h"

namespace dudley
{

  /**
  \brief
  DudleyAdapterException exception class.

  Description:
  DudleyAdapterException exception class.
  The class provides a public function returning the exception name
  */
  class DudleyAdapterException : public esysUtils::EsysException
  {

  protected:

     typedef EsysException Parent;

  public:
    /**
    \brief
    Default constructor for the exception.
    */
    DUDLEY_DLL_API
    DudleyAdapterException() : Parent() { updateMessage();}
    /**
    \brief
    Constructor for the exception.
    */
    DUDLEY_DLL_API
    DudleyAdapterException(const char *cstr) : Parent(cstr) { updateMessage();}
    /**
    \brief
    Constructor for the exception.
    */
    DUDLEY_DLL_API
    DudleyAdapterException(const std::string &str) :
    Parent(str) { updateMessage();}
    /**
    \brief
    Copy Constructor for the exception.
    */
    DUDLEY_DLL_API
    DudleyAdapterException(const DudleyAdapterException &other) : Parent(other)
      {
        updateMessage();
      }

    /// Destructor
    DUDLEY_DLL_API
    virtual ~DudleyAdapterException() THROW(NO_ARG) {}

    /**
    \brief
    Assignment operator.
    */
    DUDLEY_DLL_API
    inline DudleyAdapterException &
    operator=(const DudleyAdapterException &other ) THROW(NO_ARG)
       {
         Parent::operator=(other);
         updateMessage();
         return *this;
       }

    /**
    \brief
    Returns the name of the exception.
    */
    DUDLEY_DLL_API
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
