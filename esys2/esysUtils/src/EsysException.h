/* 
 ******************************************************************************
 *                                                                            *
 *       COPYRIGHT  ACcESS 2004 -  All Rights Reserved                        *
 *                                                                            *
 * This software is the property of ACcESS. No part of this code              *
 * may be copied in any form or by any means without the expressed written    *
 * consent of ACcESS.  Copying, use or modification of this software          *
 * by any unauthorised person is illegal unless that person has a software    *
 * license agreement with ACcESS.                                             *
 *                                                                            *
 ******************************************************************************
*/
       
#ifndef ESYSEXCEPTION_H
#define ESYSEXCEPTION_H

#include <string>
#include <sstream>
#include <exception>

namespace esysUtils {
/**
   @memo
   A base class for exception classes used within Esys system.

   @version 1.0.0 

   @doc

   Class Description:
   A base class for exception classes used within Esys system.

   Class Limitations:
   None

   Class Conditions of Use:
   None

   Throws:
   None

*/
class EsysException:public std::exception {
  public:
  /**
     @memo
     Default Constructor. Creates an exception with no message.
  */
  EsysException();
  /**
     @memo
     Constructor which creates a EsysException with the given message
     
     @param exceptionReason Input - Exception message.
  */
  EsysException(const std::string &exceptionReason);
  /**
     @memo
     Constructor which creates a EsysException with the given message
     @param cStr - Exception message.
  */
  EsysException( const char *cStr );
  /**
     @memo
     Copy constructor   
     @param inException Input - EsysException
  */
  EsysException(const EsysException &inException);
  /// Destructor
  virtual ~EsysException() throw();
  /**
     @memo
     Assignment operator.
     @param inException Input - Exception to be copied.
  */  
  EsysException &operator=(const EsysException &inException);
  /**
     @memo
     Return the exception message in the form
     <Exception Name>: <Exception Message>
     @return the exception message.
  */
  std::string toString() const;
  /**
     @memo
     Return the name of the exception. This is expected to be overloaded
     in derived classes with the derived class name.
     @return the name of the exception.
  */
  virtual std::string exceptionName() const;
  /**
     @memo
     Return the ostrstream that contains the exception message.
     This is useful for entering or adding to the message as the stream
     insertion operators can then be used.
     @return the ostrstream for the exception reason.
  */
  std::ostringstream& reason();
  /**
     @memo
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
   @memo
   Stream insertion operator for EsysExceptions
   @param output Input - Output stream.
   @param inException Input - The exception to be inserted into the output 
                            stream.
*/ 
std::ostream &operator<<(std::ostream &output, EsysException &inException);

}

#endif






