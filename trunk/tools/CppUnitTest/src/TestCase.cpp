#include <iostream>
#include <sstream>
#include <stdexcept>
#include <math.h>

#include "CppUnitTest/TestCase.h"
#include "CppUnitTest/TestResult.h"
#include "CppUnitTest/estring.h"

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

