// $Id$
/* 
 *****************************************************************************
 *                                                                           *
 *       COPYRIGHT  ACcESS  -  All Rights Reserved                           *
 *                                                                           *
 * This software is the property of ACcESS. No part of this code             *
 * may be copied in any form or by any means without the expressed written   *
 * consent of ACcESS.  Copying, use or modification of this software         *
 * by any unauthorised person is illegal unless that person has a software   *
 * license agreement with ACcESS.                                            *
 *                                                                           *
 *****************************************************************************
*/
#include "escript/DataFactory.h"
#include "escript/Data.h"

#include "DataFactoryTestCase.h"

#include <iostream>

using namespace CppUnitTest;
using namespace escript;
using namespace std;

void DataFactoryTestCase::setUp() {
  //
  // This is called before each test is run
 
}

void DataFactoryTestCase::tearDown() {
  //
  // This is called after each test has been run
 
}

void DataFactoryTestCase::testAll() {

  cout << endl;

  cout << "\tUse DataFactory functions to create some Data objects:" << endl;

  {
    cout << "\tCreate Data (DataConstant) object with Scalar data points." << endl;
    Data scalar=Scalar(1.3);
    //cout << scalar.toString() << endl;
    assert(scalar.isConstant());
    assert(scalar.getDataPointRank()==0);
    assert(scalar.getDataPointShape().empty());
  }

  {
    cout << "\tCreate DataExpanded object with Scalar data points." << endl;
    Data scalar=Scalar(1.5,FunctionSpace(),true);
    //cout << scalar.toString() << endl;
    assert(scalar.isExpanded());
    assert(scalar.getDataPointRank()==0);
    assert(scalar.getDataPointShape().empty());
  }

  {
    cout << "\tCreate Data (DataConstant) object with Vector data points." << endl;
    Data vector=Vector(1.3);
    //cout << vector.toString() << endl;
    assert(vector.isConstant());
    assert(vector.getDataPointRank()==1);
    assert(vector.getDataPointShape()[0]==1);;
  }

  {
    cout << "\tCreate Data Expanded object with Vector data points." << endl;
    Data vector=Vector(1.5,FunctionSpace(),true);
    //cout << vector.toString() << endl;
    assert(vector.isExpanded());
    assert(vector.getDataPointRank()==1);
    assert(vector.getDataPointShape()[0]==1);;
  }

  {
    cout << "\tCreate Data (DataConstant) object with Tensor data points." << endl;
    Data tensor=Tensor(1.3);
    //cout << tensor.toString() << endl;
    assert(tensor.isConstant());
    assert(tensor.getDataPointRank()==2);
    assert(tensor.getDataPointShape()[0]==1);;
    assert(tensor.getDataPointShape()[1]==1);;
  }

  {
    cout << "\tCreate Data Expanded object with Tensor data points." << endl;
    Data tensor=Tensor(1.5,FunctionSpace(),true);
    //cout << tensor.toString() << endl;
    assert(tensor.isExpanded());
    assert(tensor.getDataPointRank()==2);
    assert(tensor.getDataPointShape()[0]==1);;
    assert(tensor.getDataPointShape()[1]==1);;
  }

  {
    cout << "\tCreate Data (DataConstant) object with Tensor3 data points." << endl;
    Data tensor3=Tensor3(1.3);
    //cout << tensor3.toString() << endl;
    assert(tensor3.isConstant());
    assert(tensor3.getDataPointRank()==3);
    assert(tensor3.getDataPointShape()[0]==1);;
    assert(tensor3.getDataPointShape()[1]==1);;
    assert(tensor3.getDataPointShape()[2]==1);;
  }

  {
    cout << "\tCreate Data Expanded object with Tensor3 data points." << endl;
    Data tensor3=Tensor3(1.5,FunctionSpace(),true);
    //cout << tensor3.toString() << endl;
    assert(tensor3.isExpanded());
    assert(tensor3.getDataPointRank()==3);
    assert(tensor3.getDataPointShape()[0]==1);;
    assert(tensor3.getDataPointShape()[1]==1);;
    assert(tensor3.getDataPointShape()[2]==1);;
  }

  {
    cout << "\tCreate Data (DataConstant) object with Tensor4 data points." << endl;
    Data tensor4=Tensor4(1.3);
    //cout << tensor4.toString() << endl;
    assert(tensor4.isConstant());
    assert(tensor4.getDataPointRank()==4);
    assert(tensor4.getDataPointShape()[0]==1);;
    assert(tensor4.getDataPointShape()[1]==1);;
    assert(tensor4.getDataPointShape()[2]==1);;
    assert(tensor4.getDataPointShape()[3]==1);;
  }

  {
    cout << "\tCreate Data Expanded object with Tensor4 data points." << endl;
    Data tensor4=Tensor4(1.5,FunctionSpace(),true);
    //cout << tensor4.toString() << endl;
    assert(tensor4.isExpanded());
    assert(tensor4.getDataPointRank()==4);
    assert(tensor4.getDataPointShape()[0]==1);;
    assert(tensor4.getDataPointShape()[1]==1);;
    assert(tensor4.getDataPointShape()[2]==1);;
    assert(tensor4.getDataPointShape()[3]==1);;
  }

}

TestSuite* DataFactoryTestCase::suite ()
{
  //
  // create the suite of tests to perform.
  TestSuite *testSuite = new TestSuite ("DataFactoryTestCase");

  testSuite->addTest (new TestCaller< DataFactoryTestCase>("testAll",&DataFactoryTestCase::testAll));
  return testSuite;
}

