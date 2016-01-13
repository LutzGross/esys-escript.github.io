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
#ifndef CPPUNIT_TEST_H
#define CPPUNIT_TEST_H

#include <string>
#include <vector>
#include "CppUnitTestNamespace.h"

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

