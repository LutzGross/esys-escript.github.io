
/*****************************************************************************
*
* Copyright (c) 2016 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/


#if !defined  DataCombinationsTestCase_20040624_H
#define  DataCombinationsTestCase_20040624_H

#include <cppunit/TestFixture.h>
#include <cppunit/TestSuite.h>

#define REL_TOL ((double)1.e-10)

class DataCombinationsTestCase : public CppUnit::TestFixture
{
public:

  void testNonUpdate();
  void testUpdate();  

  static CppUnit::TestSuite* suite();

};

#endif
