
/*******************************************************
*
* Copyright (c) 2003-2010 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/

/* File extracted from finley and modified */

#if !defined  PasoException_20040526_H
#define PasoException_20040526_H
#include "system_dep.h"

#include "esysUtils/EsysException.h"

namespace paso
{

  /**
  \brief
  PasoException exception class.

  Description:
  PasoException exception class.
  The class provides a public function returning the exception name
  */
  class PasoException : public esysUtils::EsysException
  {

  protected:

     typedef EsysException Parent;

  public:
    /**
    \brief
    Default constructor for the exception.
    */
    PASOWRAP_DLL_API
    PasoException() : Parent() { updateMessage();}
    /**
    \brief
    Constructor for the exception.
    */
    PASOWRAP_DLL_API
    PasoException(const char *cstr) : Parent(cstr) { updateMessage();}
    /**
    \brief
    Constructor for the exception.
    */
    PASOWRAP_DLL_API
    PasoException(const std::string &str) :
    Parent(str) { updateMessage();}
    /**
    \brief
    Copy Constructor for the exception.
    */
    PASOWRAP_DLL_API
    PasoException(const PasoException &other) : Parent(other)
      {
        updateMessage();
      }

    /// Destructor
    PASOWRAP_DLL_API
    virtual ~PasoException() THROW(NO_ARG) {}

    /**
    \brief
    Assignment operator.
    */
    PASOWRAP_DLL_API
    inline PasoException &
    operator=(const PasoException &other ) THROW(NO_ARG)
       {
         Parent::operator=(other);
         updateMessage();
         return *this;
       }

    /**
    \brief
    Returns the name of the exception.
    */
    PASOWRAP_DLL_API
    virtual const std::string & exceptionName() const;

  private:

    //
    // the exception name is immutable and class-wide.
    // Inheritor note; you need one of these too.
    // and an overloaded exceptionName() in your .cpp implementation file. 
    static const std::string exceptionNameValue;
  };

  PASOWRAP_DLL_API
  void checkPasoError(); 
  
  
} // end of namespace
#endif
