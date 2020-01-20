
/*****************************************************************************
*
* Copyright (c) 2003-2020 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014-2017 by Centre for Geoscience Computing (GeoComp)
* Development from 2019 by School of Earth and Environmental Sciences
**
*****************************************************************************/


#if !defined  DataVectorTestCase_20050324_H
#define  DataVectorTestCase_20050324_H

#include <cppunit/TestFixture.h>
#include <cppunit/TestSuite.h>

class DataVectorTestCase : public CppUnit::TestFixture
{
public:
  void testAll();

  static CppUnit::TestSuite* suite();
};

#endif

