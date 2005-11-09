
#ifndef CPPUNIT_TEST_H
#define CPPUNIT_TEST_H

#include <string>
#include <vector>

#include "CppUnitTest/CppUnitTestNamespace.h"
BEGIN_NAMESPACE_CPPUNITTEST

class TestResult;

/*
 * A Test can be run and collect its results.
 * See TestResult.
 * 
 */

class Test
{
 public:
   Test(std::string name):testName(name) {testStatus=true;};
   Test():testName("") {testStatus=true;};
   virtual                ~Test () ;
  
   virtual void           run (TestResult *result)    = 0;
   virtual int            countTestCases ()           = 0;
   virtual std::string    toString()                  = 0;
   virtual std::string    name();
   virtual void           setArgs(const std::vector<std::string>& args);
   virtual void           setFailure() {testStatus=false;}
   virtual bool           hasFailure() {return (!testStatus);}
   virtual const std::vector<std::string>& getArgs();

 private:
   //
   // the name of the test
   std::string testName;
   std::vector<std::string> testArgs;  
  
   // indicates if this test or any of its contained tests have failed
   // if testStatus is true all tests are good, otherwise some test has failed
   bool testStatus;
};

inline Test::~Test ()
{}



// Runs a test and collects its result in a TestResult instance.
inline void Test::run (TestResult *result)
{}


// Counts the number of test cases that will be run by this test.
inline int Test::countTestCases ()
{ return 0; }


// Returns the name of the test instance. 
inline std::string Test::toString ()
{ return ""; }

//
// Return the name of the test
inline std::string Test::name() {
   return testName;
}

//
// set the argument list
inline void Test::setArgs(const std::vector<std::string>& args) {
   //
   // copy the arguments
   //
   testArgs = args;
}

//
// get the argument list
inline const std::vector<std::string>& Test::getArgs() {
   return(testArgs);
}


END_NAMESPACE_CPPUNITTEST

#endif

