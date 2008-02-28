
/* $Id$ */

/*******************************************************
 *
 *           Copyright 2003-2007 by ACceSS MNRF
 *       Copyright 2007 by University of Queensland
 *
 *                http://esscc.uq.edu.au
 *        Primary Business: Queensland, Australia
 *  Licensed under the Open Software License version 3.0
 *     http://www.opensource.org/licenses/osl-3.0.php
 *
 *******************************************************/

#ifndef ESYSEXCEPTION_H
#define ESYSEXCEPTION_H
#include "system_dep.h"

#include <string>
#include <exception>
#include <algorithm>

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
  class EsysException:public std::exception {

  public:
    /**
    \brief
    Default Constructor. Creates an exception with no message.
    */
    EsysException();

    /**
    * \brief
    Constructor which creates a EsysException with the given message
     
    @param exceptionReason Input - Exception message.
    */
    EsysException(const std::string &exceptionReason);

    /**
    * \brief
    Constructor which creates a EsysException with the given message

    @param cStr - Exception message.
    */
    EsysException( const char *cStr );

    /**
    * \brief
    Copy constructor   

    @param inException Input - EsysException
    */
    EsysException(const EsysException &inException);

    /// Destructor
    virtual ~EsysException();

    /**
    \brief
    Assignment operator.

    @param inException Input - Exception to be copied.
    */  
    EsysException &operator=(const EsysException &inException);

    /**
    \brief
    Return the exception message in the form
    &lt;Exception Name&gt;: &lt;Exception Message&gt;

    @return the exception message.
    */
    const std::string & toString() const;

    /**
    \brief
    Return the name of the exception. This is expected to be overloaded
    in derived classes with the derived class name.

    @return the name of the exception.
    */
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
    void setReason(const std::string &new_reason);

    /**
    \brief
    Return a description of the exception in the same format as the toString
    method.

    @return a description of the exception.
    */
    virtual const char* what() const throw();


  private:
    /**
    \brief
    update m_exceptionMessage after a reason update.
    **/
    inline void updateMessage();

    //
    // the exception reason
    std::string m_reason;

    //
    // the full exception message 
    std::string m_exceptionMessage;

    //
    // the exception name is immutable and class-wide.
    // Inheritor note; you need one of these too.
    // and an overloaded exceptionName in your .cpp implementation file. 
    static const std::string exceptionNameValue;

  };

  /**
  \brief
  Stream insertion (print) operator for EsysExceptions

  @param output Input - Output stream.
  @param inException Input - The exception to be inserted into the output 
  stream.
  */ 
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

  const char*  EsysException::what() const
  {
    return m_exceptionMessage.c_str();
  }



  void EsysException::updateMessage()
  {
    m_exceptionMessage = exceptionName() + ": " + m_reason;
  }

}

#endif
