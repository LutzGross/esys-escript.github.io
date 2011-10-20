
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


#if !defined ESYSEXCEPTIONTESTCASE_H
#define ESYSEXCEPTIONTESTCASE_H

#include <cppunit/TestFixture.h>
#include <cppunit/TestSuite.h>

class EsysExceptionTestCase : public CppUnit::TestFixture
{
public:
   void testCase0();
   void testCase1();
   void testCase2();

   static CppUnit::TestSuite* suite();
};

#endif

