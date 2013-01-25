
#ifndef CPPUNIT_TESTCASE_H
#define CPPUNIT_TESTCASE_H

#include <string>
#include <vector>
//#include <typeinfo>

#include "Guards.h"

#include "Test.h"

#include "CppUnitException.h"

#include "CppUnitTest/CppUnitTestNamespace.h"
BEGIN_NAMESPACE_CPPUNITTEST

class TestResult;

/*
 * A test case defines the fixture to run multiple tests. To define a test case
 * 1) implement a subclass of TestCase
 * 2) define instance variables that store the state of the fixture
 * 3) initialize the fixture state by overriding setUp
 * 4) clean-up after a test by overriding tearDown.
 *
 * Each test runs in its own fixture so there
 * can be no side effects among test runs.
 * Here is an example:
 * 
 * class MathTest : public TestCase {
 *     protected: int m_value1;
 *     protected: int m_value2;
 *
 *     public: MathTest (std::string name)
 *                 : TestCase (name) {
 *     }
 *
 *     protected: void setUp () {
 *         m_value1 = 2;
 *         m_value2 = 3;
 *     }
 * }
 * 
 *
 * For each test implement a method which interacts
 * with the fixture. Verify the expected results with assertions specified
 * by calling assert on the expression you want to test:
 * 
 *    protected: void testAdd () {
 *        int result = value1 + value2;
 *        assert (result == 5);
 *    }
 * 
 * Once the methods are defined you can run them. To do this, use
 * a TestCaller.
 *
 * Test *test = new TestCaller<MathTest>("testAdd", MathTest::testAdd);
 * test->run ();
 *
 *
 * The tests to be run can be collected into a TestSuite. CppUnit provides
 * different test runners which can run a test suite and collect the results.
 * The test runners expect a static method suite as the entry
 * point to get a test to run.
 * 
 * public: static MathTest::suite () {
 *      TestSuite *suiteOfTests = new TestSuite;
 *      suiteOfTests->addTest(new TestCaller<MathTest>("testAdd", testAdd));
 *      suiteOfTests->addTest(new TestCaller<MathTest>("testDivideByZero", testDivideByZero));
 *      return suiteOfTests;
 *  }
 * 
 * Note that the caller of suite assumes lifetime control
 * for the returned suite.
 *
 * see TestResult, TestSuite and TestCaller
 *
 */


class TestCase : public Test 
{
   REFERENCEOBJECT (TestCase)

      public:
   TestCase         (std::string Name);
   ~TestCase        ();

   virtual void        run              (TestResult *result);
   virtual TestResult  *run             ();
   virtual int         countTestCases   ();
   std::string         toString         ();

   virtual void        setUp            ();
   virtual void        tearDown         ();

   void setCheckErrors(const std::vector<std::string>& checkErrors);

   const std::vector<std::string>& getCheckErrors();


 protected:
   virtual void        runTest          ();

   TestResult          *defaultResult   ();
   void                checkImplementation 
      (bool         condition, 
       std::string  conditionExpression = "",
       long         lineNumber = CPPUNIT_UNKNOWNLINENUMBER,
       std::string  fileName = CPPUNIT_UNKNOWNFILENAME);

   void                assertImplementation 
      (bool         condition, 
       std::string  conditionExpression = "",
       long         lineNumber = CPPUNIT_UNKNOWNLINENUMBER,
       std::string  fileName = CPPUNIT_UNKNOWNFILENAME);

   void                assertEquals     (long         expected, 
                                         long         actual,
                                         long         lineNumber = CPPUNIT_UNKNOWNLINENUMBER,
                                         std::string  fileName = CPPUNIT_UNKNOWNFILENAME);

   void                assertEquals     (double       expected, 
                                         double       actual, 
                                         double       delta, 
                                         long         lineNumber = CPPUNIT_UNKNOWNLINENUMBER,
                                         std::string  fileName = CPPUNIT_UNKNOWNFILENAME);

   std::string         notEqualsMessage (long         expected, 
                                         long         actual);

   std::string         notEqualsMessage (double       expected, 
                                         double       actual);

    
 private:

   std::vector<std::string> m_errorLog;

};


// Constructs a test case
inline TestCase::TestCase (std::string name) : Test(name)  
{}


// Destructs a test case
inline TestCase::~TestCase ()
{}


// Returns a count of all the tests executed
inline int TestCase::countTestCases ()
{ return 1; }

// A hook for fixture set up
inline void TestCase::setUp ()
{}


// A hook for fixture tear down
inline void TestCase::tearDown ()
{}


// Returns the name of the test case instance
inline std::string TestCase::toString () {return name();}
//{ const std::type_info& thisClass = typeid (*this); return std::string (thisClass.name ()) + "." + name (); }



// A set of macros which allow us to get the line number
// and file name at the point of an error.
// Just goes to show that preprocessors do have some
// redeeming qualities.

#define CPPUNIT_SOURCEANNOTATION

#ifdef CPPUNIT_SOURCEANNOTATION

#undef assert
#define assert(condition)\
    (this->assertImplementation ((condition),(#condition),\
        __LINE__, __FILE__))

#undef CHECK
#define CHECK(condition)\
    (this->checkImplementation ((condition),(#condition),\
        __LINE__, __FILE__))

#else

#undef assert
#define assert(condition)\
    (this->assertImplementation ((condition),"",\
        __LINE__, __FILE__))

#undef CHECK
#define CHECK(condition)\
    (this->checkImplementation ((condition),"",\
        __LINE__, __FILE__))

#endif


// Macros for primitive value comparisons
#define assertDoublesEqual(expected,actual,delta)\
(this->assertEquals ((expected),\
        (actual),(delta),__LINE__,__FILE__))

#define assertLongsEqual(expected,actual)\
(this->assertEquals ((expected),\
        (actual),__LINE__,__FILE__))

END_NAMESPACE_CPPUNITTEST

#endif













