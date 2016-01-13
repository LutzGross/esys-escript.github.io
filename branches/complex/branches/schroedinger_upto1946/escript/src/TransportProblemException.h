
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
    TransportProblemException() : Parent() { updateMessage();}
    /**
    \brief
    Constructor for the exception.
    */
    ESCRIPT_DLL_API
    TransportProblemException(const char *cstr) :
    Parent(cstr) { updateMessage();}
    /**
    \brief
    Constructor for the exception.
    */
    ESCRIPT_DLL_API
    TransportProblemException(const std::string &str) :
    Parent(str) { updateMessage();}
    /**
    \brief
    Copy Constructor for the exception.
    */
    ESCRIPT_DLL_API
    TransportProblemException(const TransportProblemException &other) :
    Parent(other)
      {
        updateMessage();
      }

    ESCRIPT_DLL_API
    inline TransportProblemException &
    operator=(const TransportProblemException &other ) THROW_ANY
       {
         Parent::operator=(other);
         updateMessage();
         return *this;
       }

    /// Destructor
    ESCRIPT_DLL_API
    virtual ~TransportProblemException() THROW_ANY {}
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
