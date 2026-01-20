
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/


#if !defined  DataTestCase_20040624_H
#define  DataTestCase_20040624_H

#include <cppunit/TestFixture.h>
#include <cppunit/TestSuite.h>

#define REL_TOL ((double)1.e-10)

class DataTestCase : public CppUnit::TestFixture
{
public:

  void testSome();
  void testConstructors();
  void testDataConstant();
  void testDataTagged();
  void testDataTaggedExceptions();
  void testSlicing();
  void testOperations();
  void testMoreOperations();
  void testMemAlloc();
  void testCopying();
  void testResolveType();
  void testBinary();
  void testComplexSamples();
  static CppUnit::TestSuite* suite();

private:
  void testCopyingWorker(bool delayed);
  void testSlicingWorker(bool delayed);
  void testSomeDriver(bool autolazy);
};

#endif

