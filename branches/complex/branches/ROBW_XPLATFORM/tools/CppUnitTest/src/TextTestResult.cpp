
#include <iostream>
#include <iomanip>

#include "CppUnitTest/TextTestResult.h"
#include "CppUnitTest/CppUnitException.h"
#include "CppUnitTest/estring.h"

USING_NAMESPACE_CPPUNITTEST

void TextTestResult::addError (Test *test, CppUnitException *e)
{
    TestResult::addError (test, e);
    std::cout << " Error: ";
    std::cout << e->what() << std::endl; 
}

void TextTestResult::addFailure (Test *test, CppUnitException *e)
{
    TestResult::addFailure (test, e);
    std::cout << " Failure" << std::endl;

}

void TextTestResult::startTest (Test *test)
{
    TestResult::startTest (test);
    std::cout << "  Test: " << test->name();
}




void TextTestResult::printErrors (std::ostream& stream)
{
    if (testErrors () != 0) {

        if (testErrors () == 1)
            stream << "There was " << testErrors () << " error: " << std::endl;
        else
            stream << "There were " << testErrors () << " errors: " << std::endl;

        int i = 1;

        for (std::vector<TestFailure *>::iterator it = errors ().begin (); it != errors ().end (); ++it) {
            TestFailure             *failure    = *it;
            CppUnitException        *e          = failure->thrownException ();

            stream << i 
                   << ") "
                   << (*failure).failedTest()->name() << " "
                   << "line: " << (e ? estring (e->lineNumber ()) : std::string("")) << " "
                   << (e ? e->fileName () : std::string("")) << " "
                   << "\"" << failure->thrownException ()->what () << "\""
                   << std::endl;
            i++;
        }
    }

}

void TextTestResult::printFailures (std::ostream& stream) 
{
    if (testFailures () != 0) {
        if (testFailures () == 1)
            stream << "There was " << testFailures () << " failure: " << std::endl;
        else
            stream << "There were " << testFailures () << " failures: " << std::endl;

        int i = 1;

        for (std::vector<TestFailure *>::iterator it = failures ().begin (); it != failures ().end (); ++it) {
            TestFailure             *failure    = *it;
            CppUnitException        *e          = failure->thrownException ();

            stream << i 
                   << ") "
		   << (*failure).failedTest()->name() << " "
                   << "line: " << (e ? estring (e->lineNumber ()) : std::string("")) << " "
                   << (e ? e->fileName () : std::string("")) << " "
                   << "\"" << failure->thrownException ()->what () << "\""
                   << std::endl;
            i++;
        }
    }

}


void TextTestResult::print (std::ostream& stream) 
{
    printHeader (stream);
    printErrors (stream);
    printFailures (stream);

}


void TextTestResult::printHeader (std::ostream& stream)
{
    if (wasSuccessful ())
        std::cout << std::endl << "OK (" << runTests () << " tests)" << std::endl;
    else
        std::cout << std::endl
             << "!!!FAILURES!!!" << std::endl
             << "Test Results:" << std::endl
             << "Run:  "
             << runTests ()
             << "   Failures: "
             << testFailures ()
             << "   Errors: "
             << testErrors ()
             << std::endl;

}
