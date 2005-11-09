

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

#include "CppUnitTest/TestRunner.h"
#include "CppUnitTest/TextTestResult.h"
#include "CppUnitTest/Test.h"

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

