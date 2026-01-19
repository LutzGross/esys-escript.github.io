
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


#ifndef __FINLEYDOMAIN_TESTCASE_H__
#define __FINLEYDOMAIN_TESTCASE_H__

#include <cppunit/TestFixture.h>
#include <cppunit/TestSuite.h>

class FinleyDomainTestCase : public CppUnit::TestFixture
{
public:
    void testAll();
    
    static CppUnit::TestSuite* suite();
};

#endif // __FINLEYDOMAIN_TESTCASE_H__

