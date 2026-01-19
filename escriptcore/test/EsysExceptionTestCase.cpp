
/*****************************************************************************
*
* Copyright (c) 2003-2020 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/


#include "EsysExceptionTestCase.h"
#include "escript/EsysException.h"

#include <cppunit/TestCaller.h>
#include <iostream>

using namespace std;
using namespace CppUnit;
using namespace escript;

class DerivedEx : public EsysException
{
    typedef EsysException Parent;
public:
    DerivedEx(const string& str) : Parent(str) {}
};

void EsysExceptionTestCase::testCase1()
{
    string ex1Text("My first funny exception message.");
    EsysException ex1(ex1Text);

    string ex1String = ex1.what();

    //
    // exception text should contain entered exception message
    //
    CPPUNIT_ASSERT(ex1String.find(ex1Text) != string::npos);

    //
    // copy constructed exception should match original
    //
    EsysException copyEx(ex1);
    string copyString = copyEx.what();
    CPPUNIT_ASSERT(ex1String == copyString);

    //
    // check throw/catch mechanism
    //
    string ex2Text("My second funny exception message.");
    try {
        EsysException ex2(ex2Text);
        throw(ex2);
    } catch (EsysException& e) {
        //
        // exception text should contain entered exception message
        //
        string eString = e.what();
        CPPUNIT_ASSERT(eString.find(ex2Text) != string::npos);
    }
}

//
// test derived EsysException
//
void EsysExceptionTestCase::testCase2()
{
    string ex1Text("asdjhieurncidhfjsnfkjefkjndfjkhsdrdfjksdhfweh");
    DerivedEx ex1(ex1Text);

    //
    // exception text should contain entered exception message
    //
    string ex1String = ex1.what();
    CPPUNIT_ASSERT(ex1String.find(ex1Text) != string::npos);

    //
    // copy constructed exception should match original
    //
    DerivedEx copyEx(ex1);
    string copyString = copyEx.what();
    CPPUNIT_ASSERT(ex1String == copyString);

    //
    // check throw/catch mechanism
    //
    string ex2Text("pjkkjhdfbnkjerbkjsduflfkjahalkgjlklhjhj");
    try {

        DerivedEx ex2(ex2Text);
        throw(ex2);
    } catch (DerivedEx& e) {
        //
        // exception text should contain entered exception message
        //
        string eString = e.what();
        CPPUNIT_ASSERT(eString.find(ex2Text) != string::npos);
    }

    //
    // check throw/catch mechanism
    //
    string ex3Text("irfjvniouf;iarhglAKSDIghlAKSDghladg");
    try {

        DerivedEx ex3(ex3Text);
        throw(ex3);
    } catch (EsysException& e) {
        //
        // exception text should contain entered exception message
        //
        std::string eString = e.what();
        CPPUNIT_ASSERT(eString.find(ex3Text) != string::npos);
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
            cout << endl << e.what() << endl;
            throw;
        }
    } catch (EsysException& e) {
        cout << e.what() << endl;
    }

    cout << "Test EsysException may be caught as a std::exception" << endl;
    try {
        DerivedEx ex4("Exception caught as std::exception");
        throw ex4;
    } catch (exception& e) {
        // cout << e.what() << endl;
        CPPUNIT_ASSERT(e.what() == string("Exception caught as std::exception"));
    } catch (...) {
        //
        // if the exception is caught here there is a problem
        CPPUNIT_ASSERT(false);
    }
}

TestSuite* EsysExceptionTestCase::suite()
{
    //
    // create the suite of tests to perform.
    TestSuite *testSuite = new TestSuite("EsysExceptionTestCase");

    testSuite->addTest(new TestCaller<EsysExceptionTestCase>(
                "testCase1",&EsysExceptionTestCase::testCase1));
    testSuite->addTest(new TestCaller<EsysExceptionTestCase>(
                "testCase2",&EsysExceptionTestCase::testCase2));
    return testSuite;
}

