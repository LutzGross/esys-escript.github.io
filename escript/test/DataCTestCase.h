
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


#if !defined DataCTestCase_20040611_H
#define DataCTestCase_20040611_H

#include <cppunit/TestFixture.h>
#include <cppunit/TestSuite.h>

class DataCTestCase : public CppUnit::TestFixture
{
public:

  void testAll();

  static CppUnit::TestSuite* suite();
};

#endif

