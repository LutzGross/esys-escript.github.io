
/*****************************************************************************
*
* Copyright (c) 2014 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/

// Adapted from FunctionSpaceException.h

#ifndef  __ESCRIPT_SPLITWORLDEXCEPTION_H__
#define  __ESCRIPT_SPLITWORLDEXCEPTION_H__

#include "system_dep.h"
#include "esysUtils/EsysException.h"

namespace escript
{

  /**
  \brief
  SplitWorldException exception class.

  Description:
  SplitWorldException exception class.
  The class provides a public function returning the exception name
  */
  class SplitWorldException : public esysUtils::EsysException
  {

  protected:

     typedef EsysException Parent;

  public:
    /**
    \brief
    Default constructor for the exception.
    */
    ESCRIPT_DLL_API
    SplitWorldException() : Parent() { updateMessage();}
    /**
    \brief
    Constructor for the exception.
    */
    ESCRIPT_DLL_API
    SplitWorldException(const char *cstr) : Parent(cstr) { updateMessage();}
    /**
    \brief
    Constructor for the exception.
    */
    ESCRIPT_DLL_API
    SplitWorldException(const std::string &str) :
    Parent(str) { updateMessage();}
    /**
    \brief
    Copy Constructor for the exception.
    */
    ESCRIPT_DLL_API
    SplitWorldException(const SplitWorldException &other) : Parent(other)
      {
        updateMessage();
      }

    ESCRIPT_DLL_API
    inline SplitWorldException &
    operator=(const SplitWorldException &other ) THROW(NO_ARG)
       {
         Parent::operator=(other);
         updateMessage();
         return *this;
       }


    /// Destructor
    ESCRIPT_DLL_API
    virtual ~SplitWorldException() THROW(NO_ARG) {}
    /**
    \brief
    Returns the name of the exception.
    */
    ESCRIPT_DLL_API
    virtual const std::string& exceptionName() const;

  private:

    //
    // the exception name is immutable and class-wide.
    // Inheritor note; you need one of these too.
    // and an overloaded exceptionName() in your .cpp implementation file. 
    static const std::string exceptionNameValue;
  };

} // end of namespace

#endif //  __ESCRIPT_SPLITWORLDEXCEPTION_H__

