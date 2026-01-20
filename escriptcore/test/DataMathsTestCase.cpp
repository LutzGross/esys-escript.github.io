
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
#include <escript/DataVectorOps.h>
#include "DataMathsTestCase.h"
#include <escript/DataTypes.h>
#include <escript/DataVector.h>

#include <cppunit/TestCaller.h>
#include <iostream>

using namespace CppUnit;
using namespace escript;
using namespace std;
using namespace escript::DataTypes;


void DataMathsTestCase::testMatMult()
{
    cout << endl;
    cout << "\tTest result shape." << endl;

    DataTypes::ShapeType leftShape;
    leftShape.push_back(1);
    leftShape.push_back(3);
    DataTypes::RealVectorType leftData(DataTypes::noValues(leftShape),0);
//     DataArrayView leftDataView(leftData,leftShape);

    DataTypes::ShapeType rightShape;
    rightShape.push_back(3);
    rightShape.push_back(2);
    DataTypes::RealVectorType rightData(DataTypes::noValues(rightShape),0);
//     DataArrayView rightDataView(rightData,rightShape);

    DataTypes::ShapeType resultShape=determineResultShape(leftShape,rightShape);

    CPPUNIT_ASSERT(resultShape.size()==2);
    CPPUNIT_ASSERT(resultShape[0]==1);
    CPPUNIT_ASSERT(resultShape[1]==2);

    DataTypes::RealVectorType resultData(DataTypes::noValues(resultShape),0);

    cout << "\tTest matrix multiplication.";
    double aValue=0.0;
    for (int i=0;i<leftShape[0];i++) {
      for (int j=0;j<leftShape[1];j++) {
	leftData[getRelIndex(leftShape,i,j)]=++aValue;
      }
    }
    aValue=0.0;
    for (int i=0;i<rightShape[0];i++) {
      for (int j=0;j<rightShape[1];j++) {
	rightData[getRelIndex(rightShape,i,j)]=++aValue;
      }
    }

    matMult(leftData,leftShape,0,rightData,rightShape,0,resultData, resultShape);
    CPPUNIT_ASSERT((resultData[0]==22) && (resultData[1]==28));

    cout << endl;
}

void DataMathsTestCase::testReductionOp()
{

  {
    cout << endl;
    cout << "\tTest reductionOp on scalar Data.";

    // define the shape for the DataArrayView
    DataTypes::ShapeType shape;

    // allocate the data for the DataArrayView
    int npoints=4;
    DataTypes::RealVectorType data(DataTypes::noValues(shape)*npoints,0);

    int offset=0;
    // step the view along each data point in the underlying data
    for (int p=0;p<npoints;p++) {

      // assign values to the data point
      data[offset]=p;

      // apply a reduction operation to this data point and check the results
      FMax fmax_func;
//       CPPUNIT_ASSERT(std::abs(dataView.reductionOpVector(fmax_func,numeric_limits<double>::max()*-1)-p)<=REL_TOL*p);
      CPPUNIT_ASSERT(std::abs(reductionOpVector(data,shape,offset,fmax_func,numeric_limits<double>::max()*-1)-p)<=REL_TOL*p);



      if (p<npoints-1) {
        ++offset;
      }

    }

  }

  {
    cout << endl;
    cout << "\tTest reductionOp on shape (2,3) Data.";

    // define the shape for the DataArrayView
    DataTypes::ShapeType shape;
    shape.push_back(2);
    shape.push_back(3);

    // allocate the data for the DataArrayView
    int npoints=4;
    DataTypes::RealVectorType data(DataTypes::noValues(shape)*npoints,0);

    int offset=0;
    // step the view along each data point in the underlying data
    for (int p=0;p<npoints;p++) {

      // assign values to the data point
      for (int i=0;i<shape[0];i++) {
        for (int j=0;j<shape[1];j++) {
          data[offset+getRelIndex(shape,i,j)]=offset+getRelIndex(shape,i,j);
        }
      }

      // apply a reduction operation to this data point and check the results
      FMin fmin_func;
      CPPUNIT_ASSERT(std::abs(reductionOpVector(data,shape,offset,fmin_func,numeric_limits<double>::max())-offset)<=REL_TOL*std::abs(offset));

      if (p<npoints-1) {
        offset+=noValues(shape);
      }

    }

  }

  {
    cout << endl;
    cout << "\tTest reductionOp on shape (9,8,5,11) Data.";

    // define the shape for the DataArrayView
    DataTypes::ShapeType shape;
    shape.push_back(9);
    shape.push_back(8);
    shape.push_back(5);
    shape.push_back(11);

    // allocate the data for the DataArrayView
    int npoints=4;
    DataTypes::RealVectorType data(DataTypes::noValues(shape)*npoints,0);

    int offset=0;
    // step the view along each data point in the underlying data
    for (int p=0;p<npoints;p++) {

      // assign values to the data point
      for (int i=0;i<shape[0];i++) {
        for (int j=0;j<shape[1];j++) {
          for (int k=0;k<shape[2];k++) {
            for (int l=0;l<shape[3];l++) {
              data[offset+getRelIndex(shape,i,j,k,l)]=offset+getRelIndex(shape,i,j,k,l);
            }
          }
        }
      }

      // apply a reduction operation to this data point and check the results
      AbsMax<DataTypes::real_t> absmax_func;
      CPPUNIT_ASSERT(reductionOpVector(data,shape,offset,absmax_func,0)==offset+getRelIndex(shape,8,7,4,10));

      if (p<npoints-1) {
        offset+=noValues(shape);
      }

    }

  }

  cout << endl;

}

TestSuite* DataMathsTestCase::suite()
{
  // create the suite of tests to perform.
  TestSuite *testSuite = new TestSuite("DataMathsTestCase");
  testSuite->addTest(new TestCaller<DataMathsTestCase>(
              "testReductionOp",&DataMathsTestCase::testReductionOp));
  testSuite->addTest(new TestCaller<DataMathsTestCase>(
              "testMatMult",&DataMathsTestCase::testMatMult));
  return testSuite;
}

