#ifndef CPPUNIT_TESTCALLER_H
#define CPPUNIT_TESTCALLER_H

#include "CppUnitTest/CppUnitTestNamespace.h"
BEGIN_NAMESPACE_CPPUNITTEST

#include "Guards.h"
#include "TestCase.h"

/* #include <stl.h> */

/* 
 * A test caller provides access to a test case method 
 * on a test case class.  Test callers are useful when 
 * you want to run an individual test or add it to a 
 * suite.
 * 
 * Here is an example:
 * 
 * class MathTest : public TestCase {
 *         ...
 *     public:
 *         void         setUp ();
 *         void         tearDown ();
 *
 *         void         testAdd ();
 *         void         testSubtract ();
 * };
 *
 * Test *MathTest::suite () {
 *     TestSuite *suite = new TestSuite;
 *
 *     suite->addTest (new TestCaller<MathTest> ("testAdd", testAdd));
 *     return suite;
 * }
 *
 * You can use a TestCaller to bind any test method on a TestCase
 * class, as long as it returns accepts void and returns void.
 * 
 * See TestCase
 */


template <class Fixture> class TestCaller : public TestCase
{ 
 private:
   //
   // explicitly don't allow these
   TestCaller (const TestCaller<Fixture>& other); 
   TestCaller<Fixture>& operator= (const TestCaller<Fixture>& other); 

   typedef void             (Fixture::*TestMethod)();
    
 public:
   TestCaller (std::string name, TestMethod test)
      : TestCase (name), m_fixture (new Fixture (name)), m_test (test)
      {};
   ~TestCaller () {delete m_fixture;};

 protected:
   void                    runTest () 
      { 
         try {
            (m_fixture->*m_test)();
            setCheckErrors(m_fixture->getCheckErrors());
         }
         catch(CppUnitException& cpp) {
            setCheckErrors(m_fixture->getCheckErrors());
            throw CppUnitException (cpp.what(), cpp.lineNumber(),
                                    cpp.fileName());
         }
         catch (std::exception& e) {
            setCheckErrors(m_fixture->getCheckErrors());
	    //
	    // rethrow the exception
            throw;
         }

      }  

   void                    setUp ()
      { m_fixture->setArgs(getArgs());
      m_fixture->setUp ();}

   void                    tearDown ()
      { m_fixture->tearDown (); }

 private:
   Fixture*                 m_fixture;
   TestMethod               m_test;

};

END_NAMESPACE_CPPUNITTEST

#endif


