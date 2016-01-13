
/*******************************************************
*
* Copyright (c) 2003-2012 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


#if !defined  DataTaggedTestCase_20040616_H
#define  DataTaggedTestCase_20040616_H

#include <cppunit/TestFixture.h>
#include <cppunit/TestSuite.h>


class DataTaggedTestCase : public CppUnit::TestFixture
{
public:
  void testAll();
  void testAddTaggedValues();
  void testSetTaggedValue();
  void testCopyConstructors();
  void testOperations();
  void testGetSlice();
  void testSetSlice();
//   void testFunctionSpaces();

  static CppUnit::TestSuite* suite();
};

#endif

