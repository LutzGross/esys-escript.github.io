
/*****************************************************************************
*
* Copyright (c) 2003-2020 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/


#if !defined ESYSEXCEPTIONTESTCASE_H
#define ESYSEXCEPTIONTESTCASE_H

#include <cppunit/TestFixture.h>
#include <cppunit/TestSuite.h>

class EsysExceptionTestCase : public CppUnit::TestFixture
{
public:
   void testCase1();
   void testCase2();

   static CppUnit::TestSuite* suite();
};

#endif

