
/*******************************************************
*
* Copyright (c) 2003-2009 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


#include "escript/DataAlgorithm.h"
#include "escript/DataVector.h"
#include "esysUtils/EsysException.h"

#include "DataMathsTestCase.h"
#include "escript/DataTypes.h"

#include <iostream>

using namespace CppUnitTest;
using namespace esysUtils;
using namespace escript;
using namespace std;
using namespace escript::DataTypes;
using namespace escript::DataMaths;

void DataMathsTestCase::setUp() {
  //
  // This is called before each test is run

}

void DataMathsTestCase::tearDown() {
  //
  // This is called after each test has been run

}


void DataMathsTestCase::testMatMult()
{

  {
    cout << endl;
    cout << "\tTest result shape." << endl;

    DataTypes::ShapeType leftShape;
    leftShape.push_back(1);
    leftShape.push_back(3);
    DataTypes::ValueType leftData(DataTypes::noValues(leftShape),0);
//     DataArrayView leftDataView(leftData,leftShape);

    DataTypes::ShapeType rightShape;
    rightShape.push_back(3);
    rightShape.push_back(2);
    DataTypes::ValueType rightData(DataTypes::noValues(rightShape),0);
//     DataArrayView rightDataView(rightData,rightShape);

    DataTypes::ShapeType resultShape=DataMaths::determineResultShape(leftShape,rightShape);

    assert(resultShape.size()==2);
    assert(resultShape[0]==1);
    assert(resultShape[1]==2);

    DataTypes::ValueType resultData(DataTypes::noValues(resultShape),0);

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

    DataMaths::matMult(leftData,leftShape,0,rightData,rightShape,0,resultData, resultShape);
    assert((resultData[0]==22) && (resultData[1]==28));
  }

  cout << endl;

}

void DataMathsTestCase::testUnaryOp()
{

  // This typedef allows function names to be cast to pointers
  // to unary functions of the appropriate type.
  typedef double (*UnaryDFunPtr)(double);

  {
    cout << endl;
    cout << "\tTest unaryOp on scalar Data.";

    // define the shape for the Data
    DataTypes::ShapeType shape;

    // allocate the data for the Data
    int npoints=4;
    DataTypes::ValueType data(DataTypes::noValues(shape)*npoints,0);

    double tmp;
    int offset=0;
    // step the view along each data point in the underlying data
    for (int p=0;p<npoints;p++) {

      // assign values to the data point
      data[offset]=p;

      // apply a unary operation to this data point
      unaryOp(data,shape,offset,(UnaryDFunPtr)std::sin);

      // check the results
      tmp = std::sin((double)p);
      assert(std::abs(data[offset]-tmp)<=REL_TOL*std::abs(tmp));

      if (p<npoints-1) {
        offset++;
      }

    }

  }

  {
    cout << endl;
    cout << "\tTest unaryOp on shape (2,3) Data.";

    // define the shape for the DataArrayView
    DataTypes::ShapeType shape;
    shape.push_back(2);
    shape.push_back(3);

    // allocate the data for the DataArrayView
    int npoints=4;
    DataTypes::ValueType data(DataTypes::noValues(shape)*npoints,0);


    int offset=0;
    // step the view along each data point in the underlying data
    for (int p=0;p<npoints;p++) {

      // assign values to the data point
      for (int i=0;i<shape[0];i++) {
        for (int j=0;j<shape[1];j++) {
          data[offset+getRelIndex(shape,i,j)]=offset+getRelIndex(shape,i,j);
        }
      }

      // apply a unary operation to this data point
//      dataView.unaryOp((UnaryDFunPtr)std::sqrt);
      unaryOp(data,shape,offset,(UnaryDFunPtr)std::sqrt);

      // check the results
      for (int i=0;i<shape[0];i++) {
        for (int j=0;j<shape[1];j++) {
//           assert(std::abs(dataView(i,j)-std::sqrt((double)dataView.index(i,j)))<=REL_TOL*std::sqrt((double)dataView.index(i,j)));
          assert(std::abs(data[offset+getRelIndex(shape,i,j)]-std::sqrt((double)offset+getRelIndex(shape,i,j)))<=REL_TOL*std::sqrt((double)offset+getRelIndex(shape,i,j)));
        }
      }

      if (p<npoints-1) {
        offset++;
      }

    }

  }

  {
    cout << endl;
    cout << "\tTest unaryOp on shape (9,8,5,11) Data.";

    // define the shape for the DataArrayView
    DataTypes::ShapeType shape;
    shape.push_back(9);
    shape.push_back(8);
    shape.push_back(5);
    shape.push_back(11);

    // allocate the data for the DataArrayView
    int npoints=4;
    DataTypes::ValueType data(DataTypes::noValues(shape)*npoints,0);

    int offset=0;
    // step the view along each data point in the underlying data
    for (int p=0;p<npoints;p++) {

      // assign values to the data point
      for (int i=0;i<shape[0];i++) {
        for (int j=0;j<shape[1];j++) {
          for (int k=0;k<shape[2];k++) {
            for (int l=0;l<shape[3];l++) {
              data[offset+getRelIndex(shape,i,j,k,l)]=data[offset+getRelIndex(shape,i,j,k,l)]+1;
            }
          }
        }
      }

      // apply a unary operation to this data point
//       dataView.unaryOp((UnaryDFunPtr)std::log);
      unaryOp(data,shape,offset,(UnaryDFunPtr)std::log);

      // check the results
      for (int i=0;i<shape[0];i++) {
        for (int j=0;j<shape[1];j++) {
          for (int k=0;k<shape[2];k++) {
            for (int l=0;l<shape[3];l++) {
              assert(std::abs(data[offset+getRelIndex(shape,i,j,k,l)]-std::log(1+(double)data[offset+getRelIndex(shape,i,j,k,l)]))<=REL_TOL*std::abs(std::log(1+(double)data[offset+getRelIndex(shape,i,j,k,l)])));
            }
          }
        }
      }

      if (p<npoints-1) {
        offset++;
      }

    }

  }

  cout << endl;

}

void DataMathsTestCase::testBinaryOp()
{

  // This typedef allows function names to be cast to pointers
  // to binary functions of the appropriate type.
  typedef double (*BinaryDFunPtr)(double,double);

  {
    cout << endl;
    cout << "\tTest binaryOp on scalar Data.";

    // define the shape for the DataArrayViews
    DataTypes::ShapeType shape;

    // allocate the data for the DataArrayViews
    int npoints=4;
    DataTypes::ValueType data1(DataTypes::noValues(shape)*npoints,0);
    DataTypes::ValueType data2(DataTypes::noValues(shape)*npoints,0);

    int offset=0;
    // step the views along each data point in the underlying data
    for (int p=0;p<npoints;p++) {

      // assign values to the data points
      data1[offset]=p;
      data2[offset]=p;

      // apply a binary operation to these data points
      binaryOp(data1,scalarShape,offset,data2,scalarShape,offset, plus<double>());

      // check the results
      assert(data1[offset]==p+p);

      if (p<npoints-1) {
	++offset;
      }

    }

  }

  {
    cout << endl;
    cout << "\tTest binaryOp on shape (2,3) Data.";

    // define the shape for the DataArrayViews
    DataTypes::ShapeType shape;
    shape.push_back(2);
    shape.push_back(3);

    // allocate the data for the DataArrayViews
    int npoints=4;
    DataTypes::ValueType data1(DataTypes::noValues(shape)*npoints,0);
    DataTypes::ValueType data2(DataTypes::noValues(shape)*npoints,0);

    int offset=0;
    // step the views along each data point in the underlying data
    for (int p=0;p<npoints;p++) {

      // assign values to the data points
      for (int i=0;i<shape[0];i++) {
        for (int j=0;j<shape[1];j++) {
          data1[offset+getRelIndex(shape,i,j)]=offset+getRelIndex(shape,i,j);
          data2[offset+getRelIndex(shape,i,j)]=offset+getRelIndex(shape,i,j);
        }
      }

      // apply a binary operation to these data points
/*      dataView1.binaryOp(dataView2,multiplies<double>());*/
      binaryOp(data1,shape,offset,data2,shape,offset,multiplies<double>());

      // check the results
      for (int i=0;i<shape[0];i++) {
        for (int j=0;j<shape[1];j++) {
          assert(data1[offset+getRelIndex(shape,i,j)]==(offset+getRelIndex(shape,i,j))*(offset+getRelIndex(shape,i,j)));
        }
      }

      if (p<npoints-1) {
	offset+=noValues(shape);
      }

    }

  }

  {
    cout << endl;
    cout << "\tTest binaryOp on shape (9,8,5,11) Data.";

    // define the shape for the DataArrayViews
    DataTypes::ShapeType shape;
    shape.push_back(9);
    shape.push_back(8);
    shape.push_back(5);
    shape.push_back(11);

    // allocate the data for the DataArrayViews
    int npoints=4;
    DataTypes::ValueType data1(DataTypes::noValues(shape)*npoints,0);
    DataTypes::ValueType data2(DataTypes::noValues(shape)*npoints,0);

    int offset=0;
    // step the views along each data point in the underlying data
    for (int p=0;p<npoints;p++) {

      // assign values to the data points
      for (int i=0;i<shape[0];i++) {
        for (int j=0;j<shape[1];j++) {
          for (int k=0;k<shape[2];k++) {
            for (int l=0;l<shape[3];l++) {
              data1[offset+getRelIndex(shape,i,j,k,l)]=offset+getRelIndex(shape,i,j,k,l);
              data2[offset+getRelIndex(shape,i,j,k,l)]=offset+getRelIndex(shape,i,j,k,l);
            }
          }
        }
      }

      // apply a binary operation to these data points
//      dataView1.binaryOp(dataView2,multiplies<double>());
      binaryOp(data1,shape,offset,data2,shape,offset,multiplies<double>());

      // check the results
      for (int i=0;i<shape[0];i++) {
        for (int j=0;j<shape[1];j++) {
          for (int k=0;k<shape[2];k++) {
            for (int l=0;l<shape[3];l++) {
              assert(data1[offset+getRelIndex(shape,i,j,k,l)]==(offset+getRelIndex(shape,i,j,k,l))*(offset+getRelIndex(shape,i,j,k,l)));
            }
          }
        }
      }

      if (p<npoints-1) {
	offset+=noValues(shape);
      }

    }

  }

  {
    cout << endl;
    cout << "\tTest binaryOp on scalar Data and single value.";

    // define the shape for the DataArrayView
    DataTypes::ShapeType shape;

    // allocate the data for the DataArrayView
    int npoints=4;
    DataTypes::ValueType data(DataTypes::noValues(shape)*npoints,0);

    int offset=0;
    // step the view along each data point in the underlying data
    for (int p=0;p<npoints;p++) {

      // assign values to the data point
      data[offset]=p;

      // apply a binary operation to this data point
//       dataView.binaryOp(4.9,plus<double>());
      binaryOp(data,shape,offset,4.9,plus<double>());

      // check the results
      assert(data[offset]==4.9+p);

      if (p<npoints-1) {
        ++offset;
      }

    }

  }

  {
    cout << endl;
    cout << "\tTest binaryOp on shape (2,3) Data and single value.";

    // define the shape for the DataArrayView
    DataTypes::ShapeType shape;
    shape.push_back(2);
    shape.push_back(3);

    // allocate the data for the DataArrayView
    int npoints=4;
    DataTypes::ValueType data(DataTypes::noValues(shape)*npoints,0);

    int offset=0;
    // step the view along each data point in the underlying data
    for (int p=0;p<npoints;p++) {

      // assign values to the data point
      for (int i=0;i<shape[0];i++) {
        for (int j=0;j<shape[1];j++) {
          data[offset+getRelIndex(shape,i,j)]=offset+getRelIndex(shape,i,j);
        }
      }

      // apply a binary operation to the data point
//       dataView.binaryOp(5.8,multiplies<double>());
      binaryOp(data,shape,offset,5.8,multiplies<double>());

      double tmp;
      // check the results
      for (int i=0;i<shape[0];i++) {
        for (int j=0;j<shape[1];j++) {
          tmp=5.8*(offset+getRelIndex(shape,i,j));
          assert(std::abs(data[offset+getRelIndex(shape,i,j)]-tmp)<=REL_TOL*std::abs(tmp));
        }
      }

      if (p<npoints-1) {
        offset+=noValues(shape);
      }

    }

  }

  {
    cout << endl;
    cout << "\tTest binaryOp on shape (9,8,5,11) Data and single value.";

    // define the shape for the DataArrayView
    DataTypes::ShapeType shape;
    shape.push_back(9);
    shape.push_back(8);
    shape.push_back(5);
    shape.push_back(11);

    // allocate the data for the DataArrayView
    int npoints=4;
    DataTypes::ValueType data(DataTypes::noValues(shape)*npoints,0);

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

      // apply a binary operation to the data point
//       dataView.binaryOp(5.4,multiplies<double>());
      binaryOp(data,shape,offset,5.4,multiplies<double>());

      double tmp;
      // check the results
      for (int i=0;i<shape[0];i++) {
        for (int j=0;j<shape[1];j++) {
          for (int k=0;k<shape[2];k++) {
            for (int l=0;l<shape[3];l++) {
              tmp=5.4*(offset+getRelIndex(shape,i,j,k,l));
              assert(std::abs(data[offset+getRelIndex(shape,i,j,k,l)]-tmp)<=REL_TOL*std::abs(tmp));
            }
          }
        }
      }

      if (p<npoints-1) {
        offset+=noValues(shape);
      }

    }

  }

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
    DataTypes::ValueType data(DataTypes::noValues(shape)*npoints,0);

    int offset=0;
    // step the view along each data point in the underlying data
    for (int p=0;p<npoints;p++) {

      // assign values to the data point
      data[offset]=p;

      // apply a reduction operation to this data point and check the results
      FMax fmax_func;
//       assert(std::abs(dataView.reductionOp(fmax_func,numeric_limits<double>::max()*-1)-p)<=REL_TOL*p);
      assert(std::abs(reductionOp(data,shape,offset,fmax_func,numeric_limits<double>::max()*-1)-p)<=REL_TOL*p);



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
    DataTypes::ValueType data(DataTypes::noValues(shape)*npoints,0);

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
      assert(std::abs(reductionOp(data,shape,offset,fmin_func,numeric_limits<double>::max())-offset)<=REL_TOL*std::abs(offset));

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
    DataTypes::ValueType data(DataTypes::noValues(shape)*npoints,0);

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
      AbsMax absmax_func;
      assert(reductionOp(data,shape,offset,absmax_func,0)==offset+getRelIndex(shape,8,7,4,10));

      if (p<npoints-1) {
        offset+=noValues(shape);
      }

    }

  }

  cout << endl;

}

TestSuite* DataMathsTestCase::suite ()
{
  //
  // create the suite of tests to perform.
  TestSuite *testSuite = new TestSuite ("DataMathsTestCase");

  testSuite->addTest (new TestCaller< DataMathsTestCase>("testUnaryOp",&DataMathsTestCase::testUnaryOp));
  testSuite->addTest (new TestCaller< DataMathsTestCase>("testBinaryOp",&DataMathsTestCase::testBinaryOp));
  testSuite->addTest (new TestCaller< DataMathsTestCase>("testReductionOp",&DataMathsTestCase::testReductionOp));
  testSuite->addTest (new TestCaller< DataMathsTestCase>("testMatMult",&DataMathsTestCase::testMatMult));
  return testSuite;
}
