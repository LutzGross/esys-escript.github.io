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
/*
 * A command line based tool to run tests.
 * TestRunner expects as its only argument the name of a TestCase class.
 * TestRunner prints out a trace as the tests are executed followed by a
 * summary at the end.
 *
 * You can add to the tests that the TestRunner knows about by 
 * making additional calls to "addTest (...)" in main.
 *
 * Here is the synopsis:
 *
 * TestRunner [-wait] ExampleTestCase
 *
 */


#include <vector>
#include <cstdlib>

#include "TestRunner.h"
#include "TextTestResult.h"
#include "Test.h"

USING_NAMESPACE_CPPUNITTEST

void TestRunner::run (int ac, char **av)
{
    std::string  testCase;
    int     numberOfTests = 0;
    std::vector<std::string> tCases;
    std::vector<std::string> argList;
    bool failuresOccurred=false;
    //
    // process the arguments
    //
    // suites to run are specified by suite=aaa,bbb,ccc
    // more than one suite argument can be given
    // no attempt is made to remove duplicates
    //
    std::string suiteKey = "suite=";
    std::string suiteSep = ",";
    for (int i = 0; i < ac; ++i) {
       std::string testArg = av[i];
       std::string::size_type loc = testArg.find(suiteKey);
       if(i>0 && loc==0) {
	  testArg.erase(0,suiteKey.size());
	  while(testArg.size()>0) {
	     loc = testArg.find(suiteSep);
	     std::string suiteName;
	     if(loc==std::string::npos) {
		suiteName = testArg;
		testArg.erase();
	     }
	     else if(loc>0) {
		suiteName = std::string(testArg,0,loc);
		testArg.erase(0,loc+suiteSep.size());
	     }
	     else {
		testArg.erase(0,suiteSep.size());
	     }
	     if(suiteName.size()>0) {
		tCases.push_back(suiteName);
	     }
	  }
       }
       else {
	  argList.push_back(testArg);
       }
    }
    mappings::iterator it;
    Test *testToRun;
    if (tCases.size()==0) {
      //
      // no particular cases have been selected so run all cases
      for (it = m_mappings.begin (); it!=m_mappings.end(); ++it) {
	std::cout << "Suite: " << (*it).first << std::endl;
	testToRun = (*it).second;
	testToRun->setArgs(argList);
	run (testToRun);
	if(testToRun->hasFailure()) {
	   failuresOccurred=true;
	}
	numberOfTests++;
      }
    } else {
      std::vector<std::string>::const_iterator tci;
      for (tci=tCases.begin(); tci!=tCases.end(); ++tci) {
        testToRun=NULL;
        for (it = m_mappings.begin(); 
	     it != m_mappings.end() && testToRun==NULL;
	     ++it) {
	  testCase=*tci;
	  if ((*it).first == (*tci)) {
            std::cout << "Suite: " << (*it).first << std::endl;
	    testToRun = (*it).second;
	    testToRun->setArgs(argList);
	    run (testToRun);
	    if(testToRun->hasFailure()) {
	       failuresOccurred=true;
	    }
	    numberOfTests++;
	  }
	}
	if (!testToRun) {
	   std::cout << "Test " << testCase << " not found." << std::endl;
	   failuresOccurred=true;
	} 
      }
    }
    if(failuresOccurred) {
       exit(1);
    }
}


TestRunner::~TestRunner ()
{
    for (mappings::iterator it = m_mappings.begin ();
             it != m_mappings.end ();
             ++it)
        delete it->second;

}


void TestRunner::run (Test *test)
{
    TextTestResult  result;

    test->run (&result);

    std::cout << result << std::endl;
}

