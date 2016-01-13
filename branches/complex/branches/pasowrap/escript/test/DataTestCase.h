
/*******************************************************
*
* Copyright (c) 2003-2010 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


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

  static CppUnit::TestSuite* suite();

private:
  void testCopyingWorker(bool delayed);
  void testSlicingWorker(bool delayed);
  void testSomeDriver(bool autolazy);
};

#endif

