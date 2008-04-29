
/* $Id$ */

/*******************************************************
 *
 *           Copyright 2003-2007 by ACceSS MNRF
 *       Copyright 2007 by University of Queensland
 *
 *                http://esscc.uq.edu.au
 *        Primary Business: Queensland, Australia
 *  Licensed under the Open Software License version 3.0
 *     http://www.opensource.org/licenses/osl-3.0.php
 *
 *******************************************************/

#include <iostream>

#include "EsysExceptionTestCase.h"
#include "esysUtils/EsysException.h"

using namespace std;
using namespace CppUnitTest;
using namespace esysUtils;

class DerivedEx : public EsysException {

   typedef EsysException Parent;

public:

   /// Default Constructor for Exception
   DerivedEx() : Parent() { updateMessage(); }

   /// Constructor for Exception
   DerivedEx(const char *cstr) : Parent(cstr) { updateMessage(); }

   /// Constructor for Exception
   DerivedEx(const string &str) : Parent(str) { updateMessage(); }

   // Copy Constructor.
   DerivedEx(const DerivedEx &other): Parent(other) { updateMessage(); } 

   virtual DerivedEx &
   operator=(const DerivedEx &other) THROW_ANY
      {
         Parent::operator=(other);
         return *this;
      }

   /// Return the exception name
   virtual const string & exceptionName() const
      {
         return rhubarb;
      }
        
   static const string rhubarb;
};

const string DerivedEx::rhubarb("DerivedException");

void EsysExceptionTestCase::testCase0() {

	assert(true);

}

void EsysExceptionTestCase::testCase1() {

	EsysException defEx;

	//
	// default exception text should have contents of exceptionName near start
	//
	string defString = defEx.toString();
	assert(defString.find(defEx.exceptionName()) != string::npos);
	assert(defString.find(defEx.exceptionName()) < 4);

	//
	// default exception text shouldn't be much longer than contents of exception name
	//
	assert(defString.size() > defEx.exceptionName().size());
	assert((defString.size() - defEx.exceptionName().size()) < 10);

	string ex1Text("My first funny exception message.");
	EsysException ex1(ex1Text);

	//
	// exception text should have contents of exceptionName near start
	//
	string ex1String = ex1.toString();
	assert(ex1String.find(ex1.exceptionName()) != string::npos);
	assert(defString.find(ex1.exceptionName()) < 4);

	//
	// exception text should contain entered exception message
	//
	assert(ex1String.find(ex1Text) != string::npos);

	//
	// copy constructed exception should match original
	//
	EsysException copyEx(ex1);
	string copyString = copyEx.toString();
	assert(ex1String == copyString);

	//
	// copy assigned exception should match original
	//
	EsysException assEx;
	assEx = ex1;
	string assString = assEx.toString();
	assert(ex1String == assString);

	//
	// check throw/catch mechanism
	//
	string ex2Text("My second funny exception message.");
	try {

		EsysException ex2(ex2Text);
		throw(ex2);

	}

	catch(EsysException& e) {

		//
		// exception text should have contents of exceptionName near start
		//
		string eString = e.toString();
		assert(eString.find(e.exceptionName()) != string::npos);
		assert(defString.find(e.exceptionName()) < 4);

		//
		// exception text should contain entered exception message
		//
		assert(eString.find(ex2Text) != string::npos);

	}

}

//
// test derived EsysException
//
void EsysExceptionTestCase::testCase2() {

	DerivedEx defEx;
	//
	// default exception text should have contents of exceptionName near start
	//
	string defString = defEx.toString();
	assert(defString.find(defEx.exceptionName()) != string::npos);
	assert(defString.find(defEx.exceptionName()) < 4);

	//
	// default exception text shouldn't be much longer than contents of exception name
	//
	assert(defString.size() > defEx.exceptionName().size());
	assert((defString.size() - defEx.exceptionName().size()) < 10);

	string ex1Text("asdjhieurncidhfjsnfkjefkjndfjkhsdrdfjksdhfweh");
	DerivedEx ex1(ex1Text);
	//
	// exception text should have contents of exceptionName near start
	//
	string ex1String = ex1.toString();
	assert(ex1String.find(ex1.exceptionName()) != string::npos);
	assert(defString.find(ex1.exceptionName()) < 4);

	//
	// exception text should contain entered exception message
	//
	assert(ex1String.find(ex1Text) != string::npos);

	//
	// copy constructed exception should match original
	//
	DerivedEx copyEx(ex1);
	string copyString = copyEx.toString();
	assert(ex1String == copyString);

	//
	// copy assigned exception should match original
	//
	DerivedEx assEx;
	assEx = ex1;
	string assString = assEx.toString();
	assert(ex1String == assString);

	//
	// check throw/catch mechanism
	//
	string ex2Text("pjkkjhdfbnkjerbkjsduflfkjahalkgjlklhjhj");
	try {

		DerivedEx ex2(ex2Text);
		throw(ex2);

	}

	catch(DerivedEx& e) {

		//
		// exception text should have contents of exceptionName near start
		//
		string eString = e.toString();
		assert(eString.find(e.exceptionName()) != string::npos);
		assert(defString.find(e.exceptionName()) < 4);

		//
		// exception text should contain entered exception message
		//
		assert(eString.find(ex2Text) != string::npos);
	}

	//
	// check throw/catch mechanism
	//
	string ex3Text("irfjvniouf;iarhglAKSDIghlAKSDghladg");
	try {

		DerivedEx ex3(ex3Text);
		throw(ex3);

	}
	catch(EsysException& e) {

		//
		// exception text should have contents of exceptionName near start
		//
		DerivedEx ex4;
		std::string eString = e.toString();
		assert(eString.find(ex4.exceptionName()) != string::npos);
		assert(defString.find(ex4.exceptionName()) < 4);

		//
		// exception text should contain entered exception message
		//
		assert(eString.find(ex3Text) != string::npos);

	}

	//
	// test to see if exception name gets lost on rethrow
	//
        try {
	  try {
	    DerivedEx ex4("D ex4 text.");
	    throw ex4;
	  }
	  catch (EsysException& e) {
	    cout << endl << e.toString() << endl;
	    throw;
	  }
        }
        catch (EsysException& e) {
	  cout << e.toString() << endl;
	}

	cout << "Test EsysException may be caught as a std::exception" << endl;
       	try {
	  DerivedEx ex4("Exception caught as std::exception");
	  throw ex4;
       	}
 	catch (exception& e) {
          // cout << e.what() << endl;
          assert(e.what() == string("DerivedException: Exception caught"
                                    " as std::exception")
                 );
  	}
	catch (...) {
           //
           // if the exception is caught here there is a problem
	   assert(false);
	}
}

TestSuite* EsysExceptionTestCase::suite ()
{
  //
  // create the suite of tests to perform.
  TestSuite *testSuite = new TestSuite ("EsysExceptionTestCase");

  testSuite->addTest (new TestCaller<EsysExceptionTestCase>("testCase0",&EsysExceptionTestCase::testCase0));
  testSuite->addTest (new TestCaller<EsysExceptionTestCase>("testCase1",&EsysExceptionTestCase::testCase1));
  testSuite->addTest (new TestCaller<EsysExceptionTestCase>("testCase2",&EsysExceptionTestCase::testCase2));
  return testSuite;
}
