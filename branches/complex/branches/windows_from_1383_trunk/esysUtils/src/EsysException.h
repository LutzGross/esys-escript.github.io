
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
#include <sstream>
#include <exception>
#include <algorithm>

namespace esysUtils {
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
  std::string toString() const;
  /**
     \brief
     Return the name of the exception. This is expected to be overloaded
     in derived classes with the derived class name.

     @return the name of the exception.
  */
  virtual std::string exceptionName() const;
  /**
     \brief
     Return the ostrstream that contains the exception message.
     This is useful for entering or adding to the message as the stream
     insertion operators can then be used.
     
     @return the ostrstream for the exception reason.
  */
  std::ostringstream& reason();
  /**
     \brief
     Return a description of the exception in the same format as the toString
     method.

     @return a description of the exception.
  */
  virtual const char* what() const throw();
 private:
  //
  // the exception message
  std::ostringstream m_reason;
  //
  // the full exception message kept in permenant storage
  mutable std::string m_exceptionMessage;
};
/**
   \brief
   Stream insertion operator for EsysExceptions

   @param output Input - Output stream.
   @param inException Input - The exception to be inserted into the output 
                            stream.
*/ 
std::ostream &operator<<(std::ostream &output, EsysException &inException);

}

#endif
