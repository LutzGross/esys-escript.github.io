
/*****************************************************************************
*
* Copyright (c) 2003-2017 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
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
  void testComplexSamples();
  static CppUnit::TestSuite* suite();

private:
  void testCopyingWorker(bool delayed);
  void testSlicingWorker(bool delayed);
  void testSomeDriver(bool autolazy);
};

#endif

