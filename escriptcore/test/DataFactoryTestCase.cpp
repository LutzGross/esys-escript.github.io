
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

#include <escript/Data.h>
#include "DataFactoryTestCase.h"

#include <escript/DataFactory.h>

#include <cppunit/TestCaller.h>
#include <iostream>

using namespace CppUnit;
using namespace escript;
using namespace std;

void DataFactoryTestCase::testAll()
{
  cout << endl;
  cout << "\tUse DataFactory functions to create some Data objects:" << endl;

  {
    cout << "\tCreate Data (DataConstant) object with Scalar data points." << endl;
    Data scalar=Scalar(1.3);
    //cout << scalar.toString() << endl;
    CPPUNIT_ASSERT(scalar.isConstant());
    CPPUNIT_ASSERT(scalar.getDataPointRank()==0);
    CPPUNIT_ASSERT(scalar.getDataPointShape().empty());
  }

  {
    cout << "\tCreate DataExpanded object with Scalar data points." << endl;
    Data scalar=Scalar(1.5,FunctionSpace(),true);
    //cout << scalar.toString() << endl;
    CPPUNIT_ASSERT(scalar.isExpanded());
    CPPUNIT_ASSERT(scalar.getDataPointRank()==0);
    CPPUNIT_ASSERT(scalar.getDataPointShape().empty());
  }

  {
    cout << "\tCreate Data (DataConstant) object with Vector data points." << endl;
    Data vector=Vector(1.3);
    //cout << vector.toString() << endl;
    CPPUNIT_ASSERT(vector.isConstant());
    CPPUNIT_ASSERT(vector.getDataPointRank()==1);
    CPPUNIT_ASSERT(vector.getDataPointShape()[0]==1);;
  }

  {
    cout << "\tCreate Data Expanded object with Vector data points." << endl;
    Data vector=Vector(1.5,FunctionSpace(),true);
    //cout << vector.toString() << endl;
    CPPUNIT_ASSERT(vector.isExpanded());
    CPPUNIT_ASSERT(vector.getDataPointRank()==1);
    CPPUNIT_ASSERT(vector.getDataPointShape()[0]==1);;
  }

  {
    cout << "\tCreate Data (DataConstant) object with Tensor data points." << endl;
    Data tensor=Tensor(1.3);
    //cout << tensor.toString() << endl;
    CPPUNIT_ASSERT(tensor.isConstant());
    CPPUNIT_ASSERT(tensor.getDataPointRank()==2);
    CPPUNIT_ASSERT(tensor.getDataPointShape()[0]==1);;
    CPPUNIT_ASSERT(tensor.getDataPointShape()[1]==1);;
  }

  {
    cout << "\tCreate Data Expanded object with Tensor data points." << endl;
    Data tensor=Tensor(1.5,FunctionSpace(),true);
    //cout << tensor.toString() << endl;
    CPPUNIT_ASSERT(tensor.isExpanded());
    CPPUNIT_ASSERT(tensor.getDataPointRank()==2);
    CPPUNIT_ASSERT(tensor.getDataPointShape()[0]==1);;
    CPPUNIT_ASSERT(tensor.getDataPointShape()[1]==1);;
  }

  {
    cout << "\tCreate Data (DataConstant) object with Tensor3 data points." << endl;
    Data tensor3=Tensor3(1.3);
    //cout << tensor3.toString() << endl;
    CPPUNIT_ASSERT(tensor3.isConstant());
    CPPUNIT_ASSERT(tensor3.getDataPointRank()==3);
    CPPUNIT_ASSERT(tensor3.getDataPointShape()[0]==1);;
    CPPUNIT_ASSERT(tensor3.getDataPointShape()[1]==1);;
    CPPUNIT_ASSERT(tensor3.getDataPointShape()[2]==1);;
  }

  {
    cout << "\tCreate Data Expanded object with Tensor3 data points." << endl;
    Data tensor3=Tensor3(1.5,FunctionSpace(),true);
    //cout << tensor3.toString() << endl;
    CPPUNIT_ASSERT(tensor3.isExpanded());
    CPPUNIT_ASSERT(tensor3.getDataPointRank()==3);
    CPPUNIT_ASSERT(tensor3.getDataPointShape()[0]==1);;
    CPPUNIT_ASSERT(tensor3.getDataPointShape()[1]==1);;
    CPPUNIT_ASSERT(tensor3.getDataPointShape()[2]==1);;
  }

  {
    cout << "\tCreate Data (DataConstant) object with Tensor4 data points." << endl;
    Data tensor4=Tensor4(1.3);
    //cout << tensor4.toString() << endl;
    CPPUNIT_ASSERT(tensor4.isConstant());
    CPPUNIT_ASSERT(tensor4.getDataPointRank()==4);
    CPPUNIT_ASSERT(tensor4.getDataPointShape()[0]==1);;
    CPPUNIT_ASSERT(tensor4.getDataPointShape()[1]==1);;
    CPPUNIT_ASSERT(tensor4.getDataPointShape()[2]==1);;
    CPPUNIT_ASSERT(tensor4.getDataPointShape()[3]==1);;
  }

  {
    cout << "\tCreate Data Expanded object with Tensor4 data points." << endl;
    Data tensor4=Tensor4(1.5,FunctionSpace(),true);
    //cout << tensor4.toString() << endl;
    CPPUNIT_ASSERT(tensor4.isExpanded());
    CPPUNIT_ASSERT(tensor4.getDataPointRank()==4);
    CPPUNIT_ASSERT(tensor4.getDataPointShape()[0]==1);;
    CPPUNIT_ASSERT(tensor4.getDataPointShape()[1]==1);;
    CPPUNIT_ASSERT(tensor4.getDataPointShape()[2]==1);;
    CPPUNIT_ASSERT(tensor4.getDataPointShape()[3]==1);;
  }

}

TestSuite* DataFactoryTestCase::suite()
{
  TestSuite *testSuite = new TestSuite("DataFactoryTestCase");

  testSuite->addTest(new TestCaller<DataFactoryTestCase>(
              "testAll",&DataFactoryTestCase::testAll));
  return testSuite;
}

