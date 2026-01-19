
/*****************************************************************************
*
* Copyright (c) 2014-2018 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/


#ifndef _RIPLEY_SYSTEMMATRIXTESTCASE_H_
#define _RIPLEY_SYSTEMMATRIXTESTCASE_H_

#include <escript/AbstractDomain.h>
#include <escript/AbstractSystemMatrix.h>

#include <cppunit/TestFixture.h>
#include <cppunit/TestSuite.h>

class SystemMatrixTestCase : public CppUnit::TestFixture
{
public:
    virtual void setUp();

    void testSpMV_CPU_blocksize1_nonsymmetric();
    void testSpMV_CPU_blocksize2_nonsymmetric();
    void testSpMV_CPU_blocksize3_nonsymmetric();
    void testSpMV_CPU_blocksize4_nonsymmetric();
    void testSpMV_CPU_blocksize1_symmetric();
    void testSpMV_CPU_blocksize2_symmetric();
    void testSpMV_CPU_blocksize3_symmetric();
    void testSpMV_CPU_blocksize4_symmetric();

    static CppUnit::TestSuite* suite();

private:
    escript::ASM_ptr createMatrix(int blocksize, bool symmetric);
    escript::Data createInputVector(int blocksize);

    escript::JMPI mpiInfo;
    escript::Domain_ptr domain;
};

#endif // _RIPLEY_SYSTEMMATRIXTESTCASE_H_

