
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


#if !defined MeshAdapterTestCase_20040705_H
#define MeshAdapterTestCase_20040705_H

#include <cppunit/TestFixture.h>
#include <cppunit/TestSuite.h>

class MeshAdapterTestCase : public CppUnit::TestFixture
{
public:
  void testAll();

  static CppUnit::TestSuite* suite();
};

#endif

