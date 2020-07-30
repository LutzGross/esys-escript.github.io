
/*****************************************************************************
*
* Copyright (c) 2003-2018 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/


#if !defined  TaipanTestCase_20050427_H
#define  TaipanTestCase_20050427_H

#include <cppunit/TestFixture.h>
#include <cppunit/TestSuite.h>

class TaipanTestCase : public CppUnit::TestFixture
{
public:
  void testAll();
  void testN1();
  void testN0();

  static CppUnit::TestSuite* suite();
};

#endif

