//
// Permission to reproduce and create derivative works from the Software ("Software Derivative Works") 
// is hereby granted to you under the copyright of Michael Feathers. Michael Feathers also grants you 
// the right to distribute the Software and Software Derivative Works. 
//
// Michael Feathers licenses the Software to you on an "AS IS" basis, without warranty of any kind.  
// Michael Feathers HEREBY EXPRESSLY DISCLAIMS ALL WARRANTIES OR CONDITIONS, EITHER EXPRESS OR IMPLIED, 
// INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OR CONDITIONS OF MERCHANTABILITY, NON INFRINGEMENT 
// AND FITNESS FOR A PARTICULAR PURPOSE.&nbsp; You are solely responsible for determining the appropriateness 
// of using the Software and assume all risks associated with the use and distribution of this Software, 
// including but not limited to the risks of program errors, damage to or loss of data, programs or 
// equipment, and unavailability or interruption of operations.&nbsp; MICHAEL FEATHERS WILL NOT BE 
// LIABLE FOR ANY DIRECT DAMAGES OR FOR ANY SPECIAL, INCIDENTAL, OR INDIRECT DAMAGES OR FOR ANY ECONOMIC 
// CONSEQUENTIAL DAMAGES (INCLUDING LOST PROFITS OR SAVINGS), EVEN IF MICHAEL FEATHERS HAD BEEN ADVISED 
// OF THE POSSIBILITY OF SUCH DAMAGE.&nbsp; Michael Feathers will not be liable for the loss of, or damage 
// to, your records or data, or any damages claimed by you based on a third party claim. 
//
// You agree to distribute the Software and any Software Derivatives under a license agreement that: 
//
//  1) is sufficient to notify all licensees of the Software and Software Derivatives that Michael 
//     Feathers assumes no liability for any claim that may arise regarding the Software or 
//     Software Derivatives, and 
//  2) that disclaims all warranties, both express and implied, from Michael Feathers regarding the 
//     Software and Software Derivatives.&nbsp; (If you include this Agreement with any distribution 
//     of the Software and Software Derivatives you will have meet this requirement) You agree that 
//     you will not delete any copyright notices in the Software. 
//
// This Agreement is the exclusive statement of your rights in the Software as provided by Michael 
// Feathers. Except for the licenses granted to you in the second paragraph above, no other licenses 
// are granted hereunder, by estoppel, implication or otherwise. 
//
#ifndef CPPUNIT_CPPUNITEXCEPTION_H
#define CPPUNIT_CPPUNITEXCEPTION_H

/* 
 * CppUnitException is an exception that serves
 * descriptive strings through its what () method
 *
 */

#include <exception>
#include <string>
#include "CppUnitTestNamespace.h"

BEGIN_NAMESPACE_CPPUNITTEST

#define CPPUNIT_UNKNOWNFILENAME        "<unknown>"
#define CPPUNIT_UNKNOWNLINENUMBER      (-1)

//#if  defined (__sgi)  &&  ( _COMPILER_VERSION < 730 )
//class exception {
// public:
//  virtual ~exception() {}
//  virtual const char* what() const throw() { return ""; }
//};
//#endif

class CppUnitException : public std::exception
{
public:
                        CppUnitException (std::string  message    = "", 
                                          long         lineNumber = CPPUNIT_UNKNOWNLINENUMBER, 
                                          std::string  fileName   = CPPUNIT_UNKNOWNFILENAME);
                        CppUnitException (const CppUnitException& other);

    virtual             ~CppUnitException () throw() ;

    CppUnitException&   operator= (const CppUnitException& other);

    const char          *what() const throw ();

    long                lineNumber ();
    std::string         fileName ();

private:
    std::string         m_message;
    long                m_lineNumber;
    std::string         m_fileName;

};


// Construct the exception
inline CppUnitException::CppUnitException (const CppUnitException& other)
: std::exception (other)
{ 
    m_message       = other.m_message; 
    m_lineNumber    = other.m_lineNumber;
    m_fileName      = other.m_fileName;
} 

inline CppUnitException::CppUnitException (std::string message, long lineNumber, std::string fileName)
: m_message (message), m_lineNumber (lineNumber), m_fileName (fileName)
{}


// Destruct the exception
inline CppUnitException::~CppUnitException () throw()
{}


// Perform an assignment
inline CppUnitException& CppUnitException::operator= (const CppUnitException& other) 
{ 
#if defined MSWIN
   //
   // ms visual c++ compiler doesn't like the std::exception::operator=
   // approach and the SG 7.2.1 compiler doesn't accept the this->exception
   // style
   //
	this->exception::operator= (other);
#else
	std::exception::operator= (other);
#endif

    if (&other != this) {
        m_message       = other.m_message; 
        m_lineNumber    = other.m_lineNumber;
        m_fileName      = other.m_fileName;
    }

    return *this; 
}


// Return descriptive message
inline const char *CppUnitException::what() const throw ()
{ return m_message.c_str (); }

// The line on which the error occurred
inline long CppUnitException::lineNumber ()
{ return m_lineNumber; }


// The file in which the error occurred
inline std::string CppUnitException::fileName ()
{ return m_fileName; }

END_NAMESPACE_CPPUNITTEST

#endif





