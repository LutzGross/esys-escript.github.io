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
#include <iostream>
#include <sstream>
#include <stdexcept>
#ifdef _WIN32 && __INTEL_COMPILER
#include <mathimf.h>
#else
#include <math.h>
#endif

#include "TestCase.h"
#include "TestResult.h"
#include "estring.h"

USING_NAMESPACE_CPPUNITTEST

// Create a default TestResult
TestResult *TestCase::defaultResult ()
{ return new TestResult; } 


// Check for a failed general assertion 
void TestCase::assertImplementation (bool          condition,
                                     std::string   conditionExpression,
                                     long          lineNumber,
                                     std::string   fileName)
{ 
   if (!condition) 
      throw CppUnitException (conditionExpression, lineNumber, fileName); 
}

void TestCase::setCheckErrors(const std::vector<std::string>& checkErrors) {
   m_errorLog=checkErrors;
}

const std::vector<std::string>& TestCase::getCheckErrors() {
   return(m_errorLog);
}

// Check for a failed check 
void TestCase::checkImplementation (bool          condition,
                                    std::string   conditionExpression,
                                    long          lineNumber,
                                    std::string   fileName)
{ 
   if (!condition) {
      std::ostringstream errorLog;
      errorLog << "Check Failed - " << std::endl
               << "File : " << fileName << std::endl
               << "Line Number : " << lineNumber << std::endl
               << "Condition : " << conditionExpression << std::endl
               << std::endl;
      m_errorLog.push_back(errorLog.str());
   }
}


// Check for a failed equality assertion 
void TestCase::assertEquals (long        expected, 
                             long        actual,
                             long        lineNumber,
                             std::string fileName)
{ 
   if (expected != actual) 
      assertImplementation (false, notEqualsMessage(expected, actual), lineNumber, fileName); 
}


// Check for a failed equality assertion
void TestCase::assertEquals (double        expected, 
                             double        actual, 
                             double        delta,
                             long          lineNumber,
                             std::string   fileName)
{ 
   if (fabs (expected - actual) > delta) 
      assertImplementation (false, notEqualsMessage(expected, actual), lineNumber, fileName); 

}


// Run the test and catch any exceptions that are triggered by it 
void TestCase::run (TestResult *result)
{
   //
   // flag, if true the test has been a success
   bool success=false;
   result->startTest (this);
   try {
      m_errorLog.clear();
      setUp ();
      runTest ();
      success=true;
   }

   catch (CppUnitException& e) {
      CppUnitException *copy = new CppUnitException (e);
      result->addFailure (this, copy);
    
   }
   catch (std::exception& e) {
     std::cout << "TestCase Caught exception: " << e.what() << std::endl;
     result->addError (this, new CppUnitException (e.what ()));
    
   }
    
   catch (...) {
      CppUnitException *e = new CppUnitException ("unknown exception");
      result->addError (this, e);
   }
  
   if(m_errorLog.size()>0) {
      std::ostringstream errorLog;
      errorLog << std::endl;
      for(int i=0; i<m_errorLog.size(); ++i) {
         errorLog << m_errorLog[i];
      }
      errorLog << std::endl;
      result->addFailure (this, new CppUnitException (errorLog.str()));
   }

   try {
      tearDown ();
   }

   catch (CppUnitException& e) {
      std::ostringstream temp;
      temp << "Exception from tearDown: " << e.what () << std::flush;
      CppUnitException *copy = new CppUnitException (temp.str().c_str());
      result->addFailure (this, copy);
      success=false;
   }
   catch (std::exception& e) {
      std::ostringstream temp;
      temp << "Exception from tearDown: " << e.what () << std::flush;
      result->addError (this, new CppUnitException (temp.str().c_str()));
      success=false;
   }
    
   catch (...) {
      std::ostringstream temp;
      temp << "Exception from tearDown: Unknown exception." << std::flush;
      CppUnitException *e = new CppUnitException (temp.str().c_str());
      result->addError (this, e);
      success=false;
   }
  
   result->endTest (this);
   if (success) {
      std::cerr << " OK" << std::endl;
   }
   else {
      setFailure();
   }
}


// A default run method 
TestResult *TestCase::run ()
{
   TestResult *result = defaultResult ();

   run (result);
   return result;

}


// All the work for runTest is deferred to subclasses 
void TestCase::runTest ()
{
}


// Build a message about a failed equality check 
std::string TestCase::notEqualsMessage (long expected, long actual)
{ 
   return "expected: " + estring (expected) + " but was: " + estring (actual); 
}


// Build a message about a failed equality check 
std::string TestCase::notEqualsMessage (double expected, double actual)
{ 
   return "expected: " + estring (expected) + " but was: " + estring (actual); 
}

