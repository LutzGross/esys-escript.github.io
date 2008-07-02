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
#ifndef CPPUNIT_TESTFAILURE_H
#define CPPUNIT_TESTFAILURE_H

#include "Guards.h"
#include "CppUnitException.h"
#include "Test.h"
#include "CppUnitTestNamespace.h"

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


