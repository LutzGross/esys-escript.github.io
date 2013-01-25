
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


#if !defined multi_arrayTestCase_20040319_H
#define multi_arrayTestCase_20040319_H

#include <cppunit/TestFixture.h>
#include <cppunit/TestSuite.h>

class multi_arrayTestCase : public CppUnit::TestFixture
{
public:
  void testAll();

  static CppUnit::TestSuite* suite();
};

#endif

