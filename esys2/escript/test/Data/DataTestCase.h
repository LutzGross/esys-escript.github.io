// $Id$
/* 
 *****************************************************************************
 *                                                                           *
 *       COPYRIGHT  ACcESS  -  All Rights Reserved                           *
 *                                                                           *
 * This software is the property of ACcESS. No part of this code             *
 * may be copied in any form or by any means without the expressed written   *
 * consent of ACcESS.  Copying, use or modification of this software         *
 * by any unauthorised person is illegal unless that person has a software   *
 * license agreement with ACcESS.                                            *
 *                                                                           *
 *****************************************************************************
*/
#if !defined  DataTestCase_20040624_H
#define  DataTestCase_20040624_H

#include "CppUnitTest/TestCase.h"
#include "CppUnitTest/TestSuite.h"
#include "CppUnitTest/TestCaller.h"

class DataTestCase : public CppUnitTest::TestCase
{
 public:

  //
  // setUp is called before each test method to set up test state
  void setUp();

  //
  // tearDown is called after each test method is called.
  void tearDown(); 

  //
  // A test method must return void and have no arguments

  // Data class
  void testAll();

  void testMore();

  //
  // Test the constructors
  void testConstructors();

  //
  // Tests for DataConstant
  void testDataConstant();

  //
  // Tests for DataTagged
  void testDataTagged();

  //
  // Tests for DataTagged exceptions
  void testDataTaggedExceptions();

  //
  // Tests for slicing
  void testSlicing();
  void testOperations();

  DataTestCase (std::string name) : TestCase (name) {}
  ~DataTestCase() {}

  //
  //
  // return the suite of tests to perform
  //

  static CppUnitTest::TestSuite* suite ();

};

#endif
