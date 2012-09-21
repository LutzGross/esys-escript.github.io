
/*****************************************************************************
*
* Copyright (c) 2003-2012 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development since 2012 by School of Earth Sciences
*
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

  static CppUnit::TestSuite* suite();

private:
  void testCopyingWorker(bool delayed);
  void testSlicingWorker(bool delayed);
  void testSomeDriver(bool autolazy);
};

#endif

