
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


#ifndef ESYSEXCEPTION_H
#define ESYSEXCEPTION_H
#include "system_dep.h"

#include <string>
#include <exception>
#include <iostream>

namespace esysUtils
{
  /**
  \page esys_exception Esys Exceptions
  A base class for exception classes used within Esys system.

  \version 1.0.0 

  \section class_desc Class Description:
  A base class for exception classes used within Esys system.

  \section class_limits Class Limitations:
  None

  \section class_conds Class Conditions of Use:
  None

  \section throws Throws:
  None

  */
  class EsysException : public std::exception
  {

  protected:

     typedef std::exception Parent;


  public:
    /**
    \brief
    Default Constructor. Creates an exception with no message.
    */
    ESYSUTILS_DLL_API
    EsysException();

    /**
    * \brief
    Constructor which creates a EsysException with the given message
     
    @param exceptionReason Input - Exception message.
    */
    ESYSUTILS_DLL_API
    EsysException(const std::string &exceptionReason);

    /**
    * \brief
    Constructor which creates a EsysException with the given message

    @param cStr - Exception message.
    */
    ESYSUTILS_DLL_API
    EsysException( const char *cStr );

    /**
    * \brief
    Copy constructor   

    @param other Input - EsysException
    */
    ESYSUTILS_DLL_API
    EsysException(const EsysException &other);

    /// Destructor
    ESYSUTILS_DLL_API
    virtual ~EsysException() THROW_ANY;

    /**
    \brief
    Assignment needed to override any automatic assignment
    of std::exception, which can potentially copy around char *'s,
    causeing trouble in some implementations of STL.
    It will only copy the reason string, and update the message.

    @return re-assigned exception.
    */
    ESYSUTILS_DLL_API
    EsysException &
    operator=(const EsysException &other) THROW_ANY;

    /**
    \brief
    Return the exception message in the form
    &lt;Exception Name&gt;: &lt;Exception Message&gt;

    @return the exception message.
    */
    inline
    const std::string & toString() const;

    /**
    \brief
    Return the name of the exception. This is expected to be overloaded
    in derived classes with the derived class name.

    @return the name of the exception.
    */
    ESYSUTILS_DLL_API
    virtual const std::string & exceptionName() const;

    /**
    \brief
    Return a reference to the string that contains the exception reason.
     
    @return the string for the exception reason.
    */
    inline
    const std::string& reason() const;

    /**
    \brief
    set the string for the reason for the exception.
    This allows ousiders to modify m_reason, but the practice is discouraged.
    If string insertions are required, use string methods.
    */
    inline
    void setReason(const std::string &new_reason);

    /**
    \brief
    Return a description of the exception in the same format as the toString
    method.

    @return a description of the exception.
    */
    ESYSUTILS_DLL_API
    inline
    virtual const char* what() const THROW_ANY;


    /**
    \brief
    update m_exceptionMessage after a reason update.
    **/
    inline
    void updateMessage();


  private:
    //
    // the exception reason
    std::string m_reason;

    //
    // the full exception message 
    std::string m_exceptionMessage;

    //
    // the exception name is immutable and class-wide.
    // Inheritor note; you need one of these too.
    // and an overloaded exceptionName() in your .cpp implementation file. 
    static const std::string exceptionNameValue;

  };

  /**
  \brief
  Stream insertion (print) operator for EsysExceptions

  @param output Input - Output stream.
  @param inException Input - The exception to be inserted into the output 
  stream.
  */ 
  ESYSUTILS_DLL_API
  std::ostream &operator<<(std::ostream &output, EsysException &inException);


  ////////////////////////////////////////////////////////////////////

  const std::string & EsysException::reason() const
  {
    return m_reason;
  }
  
  // return the message as a std::string
  const std::string & EsysException::toString() const
  {
    return m_exceptionMessage;
  }

  void EsysException::setReason(const std::string &new_reason)
  {
    m_reason = new_reason;
    updateMessage();
  }

  const char*  EsysException::what() const THROW_ANY
  {
    return m_exceptionMessage.c_str();
  }

  void EsysException::updateMessage()
  {
    m_exceptionMessage = exceptionName() + ": " + m_reason;
  }

}

#endif
