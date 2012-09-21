
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

