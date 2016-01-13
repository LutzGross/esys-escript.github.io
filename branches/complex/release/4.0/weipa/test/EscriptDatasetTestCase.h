
/*****************************************************************************
*
* Copyright (c) 2003-2014 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/


#ifndef _ESCRIPTDATASETTESTCASE_H_
#define _ESCRIPTDATASETTESTCASE_H_

#include <escript/AbstractDomain.h>

#include <cppunit/TestFixture.h>
#include <cppunit/TestSuite.h>

class EscriptDatasetTestCase : public CppUnit::TestFixture
{
public:
    void testBase();
#if USE_DUDLEY
    void testDudley();
#endif
#if USE_FINLEY
    void testFinley();
#endif
#if USE_RIPLEY
    void testRipley();
#endif
#if USE_SPECKLEY
    void testSpeckley();
#endif

    static CppUnit::TestSuite* suite();

private:
    void runDomainTests(escript::Domain_ptr dom);
    void checkVTKfile(std::string filename);
    int getDataArrayLength(std::istream& is);
};

#endif

