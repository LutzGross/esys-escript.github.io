#ifndef CPPUNIT_TESTFAILURE_H
#define CPPUNIT_TESTFAILURE_H

#include "CppUnitTest/Guards.h"
#include "CppUnitTest/CppUnitException.h"
#include "CppUnitTest/Test.h"

#include "CppUnitTest/CppUnitTestNamespace.h"
BEGIN_NAMESPACE_CPPUNITTEST

/*
 * A TestFailure collects a failed test together with
 * the caught exception.
 *
 * TestFailure assumes lifetime control for any exception
 * passed to it.  The lifetime of tests is handled by
 * their TestSuite (if they have been added to one) or
 * whomever creates them.
 * 
 * see TestResult
 * see TestSuite
 *
 */

class TestFailure 
{
    REFERENCEOBJECT (TestFailure)

public:
                        TestFailure (Test *failedTest, CppUnitException *thrownException);
                        ~TestFailure ();

    Test                *failedTest ();
    CppUnitException    *thrownException ();
    std::string         toString ();

protected:
    Test                *m_failedTest;
    CppUnitException    *m_thrownException;

};


// Constructs a TestFailure with the given test and exception.
inline TestFailure::TestFailure (Test *failedTest, CppUnitException *thrownException)
 : m_failedTest (failedTest), m_thrownException (thrownException) 
{}


// Deletes the owned exception.
inline TestFailure::~TestFailure ()
{ delete m_thrownException; }


// Gets the failed test.
inline Test *TestFailure::failedTest ()
{ return m_failedTest; }


// Gets the thrown exception.
inline CppUnitException *TestFailure::thrownException ()
{ return m_thrownException; }

END_NAMESPACE_CPPUNITTEST

#endif


