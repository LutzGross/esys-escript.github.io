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
#include <iomanip>

#include "TextTestResult.h"
#include "CppUnitException.h"
#include "estring.h"

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
