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

#include "escriptcpp/DataExpanded.h"
#include "escriptcpp/DataArrayView.h"
#include "escriptcpp/DataAlgorithm.h"
#include "DataAlgorithmAdapterTestCase.h"

#include <iostream>
#include <algorithm>
#include <math.h>
#include <limits>

using namespace CppUnitTest;
using namespace std;
using namespace escript;

void DataAlgorithmAdapterTestCase::setUp() {
  //
  // This is called before each test is run
 
}

void DataAlgorithmAdapterTestCase::tearDown() {
  //
  // This is called after each test has been run
 
}

void DataAlgorithmAdapterTestCase::testAll() {

  cout << endl;

  cout << "\tTesting FMax." << endl;

  FMax fmax;
  assert(fmax(5,6)==6);
  assert(fmax(5,-6)==5);
  assert(fmax(0,0)==0);
  assert(fmax(15,-96)==15);

  DataAlgorithmAdapter<FMax> sup(numeric_limits<double>::max()*-1);
  sup.resetResult();
  sup(-1);
  sup(-2);
  sup(-14);
  sup(3);
  assert(sup.getResult()==3);

  cout << "\tTesting AbsMax." << endl;

  AbsMax absmax;
  assert(absmax(5,6)==6);
  assert(absmax(5,-6)==6);
  assert(absmax(0,0)==0);
  assert(absmax(15,-96)==96);

  DataAlgorithmAdapter<AbsMax> Lsup(0);
  Lsup.resetResult();
  Lsup(-2);
  Lsup(2);
  Lsup(5);
  Lsup(-10);
  assert(Lsup.getResult()==10);

  cout << "\tTesting FMin." << endl;

  FMin fmin;
  assert(fmin(5,6)==5);
  assert(fmin(5,-6)==-6);
  assert(fmin(0,0)==0);
  assert(fmin(15,-96)==-96);

  DataAlgorithmAdapter<FMin> inf(numeric_limits<double>::max());
  inf.resetResult();
  inf(1);
  inf(12);
  inf(2);
  inf(99);
  assert(inf.getResult()==1);

  cout << "\tTesting Length." << endl;

  Length lngth;
  assert(lngth(5,6)==std::sqrt(61.0));
  assert(lngth(5,-6)==std::sqrt(61.0));
  assert(lngth(0,0)==std::sqrt(0.0));
  assert(lngth(15,-96)==std::sqrt(9441.0));

  DataAlgorithmAdapter<Length> length(0);
  length.resetResult();
  length(2);
  length(4);
  length(6);
  length(8);
  assert(length.getResult()==std::sqrt(120.0));
  length.resetResult();
  length(1.5);
  length(2.5);
  length(3.5);
  length(4.5);
  assert(length.getResult()==std::sqrt(41.0));

  cout << "\tTesting Trace." << endl;

  Trace trce;
  assert(trce(5,6)==11);
  assert(trce(5,-6)==-1);
  assert(trce(0,0)==0);
  assert(trce(15,-96)==-81);

  DataAlgorithmAdapter<Trace> trace(0);
  trace.resetResult();
  trace(1);
  trace(2);
  trace(3);
  trace(4);
  trace(5);
  assert(trace.getResult()==15);
  trace.resetResult();
  trace(1.5);
  trace(2.5);
  trace(3.5);
  trace(4.5);
  trace(5.5);
  assert(trace.getResult()==17.5);

}

void DataAlgorithmAdapterTestCase::testAlgorithm() {

  cout << endl;

  {

    cout << "\tTest algorithm on Data objects with a single rank 2 data-point." << endl;

    // define the shape for the DataArrayView
    DataArrayView::ShapeType shape;
    shape.push_back(2);
    shape.push_back(3);

    // allocate the data for the DataArrayView
    DataArrayView::ValueType dataArray(DataArrayView::noValues(shape),0);

    // construct DataArrayView
    DataArrayView dataView(dataArray,shape);

    // assign values to the data point
    for (int i=0;i<shape[0];i++) {
      for (int j=0;j<shape[1];j++) {
        dataView(i,j)=dataView.index(i,j);
      }
    }

    // create a few Data objects from the created DataArrayView
    DataExpanded dataExp(dataView,FunctionSpace());
    DataConstant dataCon(dataView,FunctionSpace());
    DataTagged   dataTag(dataCon);

    // test algorithm on DataExpanded
    FMin fmin_func;
    assert(escript::algorithm(dataExp,fmin_func,numeric_limits<double>::max())==0);
    FMax fmax_func;
    assert(escript::algorithm(dataExp,fmax_func,numeric_limits<double>::max()*-1)==5);

    // test algorithm on DataTagged
    assert(escript::algorithm(dataTag,fmin_func,numeric_limits<double>::max())==0);
    assert(escript::algorithm(dataTag,fmax_func,numeric_limits<double>::max()*-1)==5);

    // test algorithm on DataConstant
    assert(escript::algorithm(dataCon,fmin_func,numeric_limits<double>::max())==0);
    assert(escript::algorithm(dataCon,fmax_func,numeric_limits<double>::max()*-1)==5);

  }

}

void DataAlgorithmAdapterTestCase::testDpAlgorithm() {

  cout << endl;

  {

    cout << "\tTest dp_algorithm on Data objects with a single rank 2 data-point." << endl;

    // define the shapes for the DataArrayViews
    DataArrayView::ShapeType shape;
    shape.push_back(2);
    shape.push_back(3);
    DataArrayView::ShapeType shape2;

    // allocate the data for the DataArrayViews
    DataArrayView::ValueType dataArray(DataArrayView::noValues(shape),0);
    DataArrayView::ValueType dataArray2(DataArrayView::noValues(shape2),0);

    // construct DataArrayViews
    DataArrayView dataView(dataArray,shape);
    DataArrayView dataView2(dataArray2,shape2);

    // assign values to the data point
    for (int i=0;i<shape[0];i++) {
      for (int j=0;j<shape[1];j++) {
        dataView(i,j)=dataView.index(i,j);
      }
    }

    // create a few Data objects from the created DataArrayViews
    DataExpanded dataExp(dataView,FunctionSpace());
    DataConstant dataCon(dataView,FunctionSpace());
    DataTagged   dataTag(dataCon);

    // and create Data objects to receive the results of the dp_algorithm calls
    DataExpanded dataExp2(dataView2,FunctionSpace());
    DataConstant dataCon2(dataView2,FunctionSpace());
    DataTagged   dataTag2(dataCon2);

    // test dp_algorithm on DataExpanded
    FMin fmin_func;
    escript::dp_algorithm(dataExp,dataExp2,fmin_func,numeric_limits<double>::max());
    assert(dataExp2.getDataPoint(0,0)()==0);
    FMax fmax_func;
    escript::dp_algorithm(dataExp,dataExp2,fmax_func,numeric_limits<double>::max()*-1);
    assert(dataExp2.getDataPoint(0,0)()==5);

    // test dp_algorithm on DataTagged
    escript::dp_algorithm(dataTag,dataTag2,fmin_func,numeric_limits<double>::max());
    assert(dataTag2.getDataPoint(0,0)()==0);
    escript::dp_algorithm(dataTag,dataTag2,fmax_func,numeric_limits<double>::max()*-1);
    assert(dataTag2.getDataPoint(0,0)()==5);

    // test dp_algorithm on DataConstant
    escript::dp_algorithm(dataCon,dataCon2,fmin_func,numeric_limits<double>::max());
    assert(dataCon2.getDataPoint(0,0)()==0);
    escript::dp_algorithm(dataCon,dataCon2,fmax_func,numeric_limits<double>::max()*-1);
    assert(dataCon2.getDataPoint(0,0)()==5);

  }

}

TestSuite* DataAlgorithmAdapterTestCase::suite ()
{
  //
  // create the suite of tests to perform.
  TestSuite *testSuite = new TestSuite ("DataAlgorithmAdapterTestCase");

  testSuite->addTest (new TestCaller< DataAlgorithmAdapterTestCase>("testAll",&DataAlgorithmAdapterTestCase::testAll));
  testSuite->addTest (new TestCaller<DataAlgorithmAdapterTestCase>("testAlgorithm",&DataAlgorithmAdapterTestCase::testAlgorithm));
  testSuite->addTest (new TestCaller<DataAlgorithmAdapterTestCase>("testDpAlgorithm",&DataAlgorithmAdapterTestCase::testDpAlgorithm));
  return testSuite;
}
