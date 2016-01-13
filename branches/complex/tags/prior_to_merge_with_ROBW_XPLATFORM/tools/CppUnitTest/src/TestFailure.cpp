
#include "CppUnitTest/TestFailure.h"
#include "CppUnitTest/Test.h"

USING_NAMESPACE_CPPUNITTEST

// Returns a short description of the failure.
std::string TestFailure::toString () 
{ 
    return m_failedTest->toString () + ": " + m_thrownException->what ();
}
