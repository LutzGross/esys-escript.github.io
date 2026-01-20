
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* https://github.com/LutzGross/esys-escript.github.io
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/

#include <escript/DataTypes.h>

#include "multi_arrayTestCase.h"

#include <boost/multi_array.hpp>
#include <cppunit/TestCaller.h>
#include <iostream>

using namespace CppUnit;
using namespace std;


void multi_arrayTestCase::testAll()
{
  //
  // Test boost multi_array
  typedef boost::multi_array<double, 3> ArrayType;
  typedef ArrayType::index_range range;
  ArrayType testArr(boost::extents[1][2][3]);
  cout << endl;
  for (ArrayType::index i=0;i!=1;++i) {
    for (ArrayType::index j=0;j!=2;++j) {
      for (ArrayType::index k=0;k!=3;++k) {
	testArr[i][j][k]=k+j*3+i*(1*2);
        cout << "(" << i << "," << j << "," << k << ") " 
	     << testArr[i][j][k] << endl;
      }
    }
  }
  ArrayType::index_gen indices;
  ArrayType::array_view<2>::type OneView=testArr[ indices[range()][range()][1] ];
  for (ArrayType::index i=0;i!=1;++i) {
    for (ArrayType::index j=0;j!=2;++j) {
      cout << "(" << i << "," << j << ",1) " << OneView[i][j] << endl;
    }
  }

}

TestSuite* multi_arrayTestCase::suite ()
{
  TestSuite *testSuite = new TestSuite("multi_arrayTestCase");

  testSuite->addTest(new TestCaller<multi_arrayTestCase>(
              "testAll",&multi_arrayTestCase::testAll));
  return testSuite;
}

