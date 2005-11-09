
#ifndef CPPUNIT_CPPUNITEXCEPTION_H
#define CPPUNIT_CPPUNITEXCEPTION_H

/* 
 * CppUnitException is an exception that serves
 * descriptive strings through its what () method
 *
 */

#include <exception>
#include <string>

#include "CppUnitTest/CppUnitTestNamespace.h"
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





