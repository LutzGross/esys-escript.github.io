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
#include "escript/Data/DataArrayView.h"
#include "escript/Data/DataArray.h"
#include "DataArrayViewTestCase.h"
#include "esysUtils/EsysException.h"

#include <iostream>

using namespace CppUnitTest;
using namespace esysUtils;
using namespace escript;
using namespace std;

void DataArrayViewTestCase::setUp() {
  //
  // This is called before each test is run
 
}

void DataArrayViewTestCase::tearDown() {
  //
  // This is called after each test has been run
 
}

void DataArrayViewTestCase::testSlicing() {

  // Test some arbitrary slicing

  cout << endl;

  {
    DataArrayView::RegionType region;
    region.push_back(DataArrayView::RegionType::value_type(1,5));
    DataArrayView::ShapeType resultShape;
    resultShape.push_back(4);
    assert(DataArrayView::getResultSliceShape(region)==resultShape);
    region.push_back(DataArrayView::RegionType::value_type(2,5));
    resultShape.push_back(3);
    assert(DataArrayView::getResultSliceShape(region)==resultShape);
  }

  {
    // Say an empty region can select a scalar
    DataArrayView::RegionType region;
    DataArrayView::ShapeType sourceShape;
    DataArray source(sourceShape,2.0);
    DataArray target;
    target.getView().copySlice(source.getView(),region);
    assert(source.getView()==target.getView());
  }

  {
    // try a slice from a 1 dimensional array
    DataArrayView::RegionType region;
    region.push_back(DataArrayView::RegionType::value_type(2,4));
    DataArrayView::ShapeType sourceShape;
    sourceShape.push_back(6);
    DataArray source(sourceShape);
    for (int i=0;i<sourceShape[0];++i) {
      source.getView()(i)=i;
    }
    DataArray target(DataArrayView::getResultSliceShape(region));
    target.getView().copySlice(source.getView(),region);
    for (int i=region[0].first;i<region[0].second;++i) {
      assert(source.getView()(i)==target.getView()(i-2));
    }
  }

  {
    // try a slice from a 2 dimensional array
    DataArrayView::RegionType region;
    region.push_back(DataArrayView::RegionType::value_type(2,4));
    region.push_back(DataArrayView::RegionType::value_type(0,2));
    DataArrayView::ShapeType sourceShape;
    sourceShape.push_back(6);
    sourceShape.push_back(3);
    DataArray source(sourceShape);
    int val=0;
    for (int i=0;i<sourceShape[0];++i) {
      for (int j=0;j<sourceShape[1];++j) {
	source.getView()(i,j)=val++;
      }
    }
    DataArray target(DataArrayView::getResultSliceShape(region));
    target.getView().copySlice(source.getView(),region);
    for (int i=region[0].first;i<region[0].second;++i) {
      for (int j=region[1].first;j<region[1].second;++j) {
	assert(source.getView()(i,j)==target.getView()(i-2,j));
      }
    }
  }

}

void DataArrayViewTestCase::testShapeToString() {
  cout << endl;
  DataArrayView::ShapeType shape;
  assert(DataArrayView::shapeToString(shape)=="()");
  shape.push_back(5);
  assert(DataArrayView::shapeToString(shape)=="(5)");
  shape.push_back(2);
  assert(DataArrayView::shapeToString(shape)=="(5,2)");
}

void DataArrayViewTestCase::testScalarView() {
  //
  // Create a vector containing data for three scalars
  // and check three scalar views return the appropriate data
  DataArrayView::ShapeType vShape;
  DataArrayView::ValueType vData;
  vData.push_back(0);
  vData.push_back(1);
  vData.push_back(2);
  int sZero=0;
  int sOne=1;
  int sTwo=2;
  DataArrayView zView(vData,vShape,sZero);
  DataArrayView oView(vData,vShape,sOne);
  DataArrayView tView(vData,vShape,sTwo);
  vShape.push_back(3);
  DataArrayView oneVView(vData,vShape,0);
  assert(zView()==0);
  assert(oView()==1);
  assert(tView()==2);
  assert(tView.noValues()==1);
  assert(oView.getRank()==0);
  //
  // copy the one view to the zero view
  zView.copy(oView);
  assert(zView==oView);
  zView.checkShape(oView.getShape());
  cout << endl << "\tTest shape mismatch functions." << endl;
  if (!zView.checkShape(oneVView.getShape())) {
    //cout << zView.createShapeErrorMessage("Error - Shape mismatch.", oneVView.getShape()) << endl;
    assert(true);
  } else {
    assert(false);
  }
  zView.unaryOp(negate<double>());
  assert(zView()==-1);
  zView.binaryOp(oView,plus<double>());
  assert(zView()==0);
  //
  // test view.operator!=
  assert(zView!=oView);
}

void DataArrayViewTestCase::testAll()
{

  {
    cout << endl;

    cout << "\tTest empty DataArrayView." << endl;

    DataArrayView defView;

    assert(defView.getOffset()==0);
    assert(defView.getShape().empty());
    assert(defView.noValues()==0);
    assert(defView.getRank()==0);
  }

  {
    cout << endl;

    // define the shape for the data array view
    cout << "\tTest shape (5):" << endl;
    DataArrayView::ShapeType shape;
    shape.push_back(5);
    int offset=5;

    // allocate the space external to the DataArrayView
    DataArrayView::ValueType data(DataArrayView::noValues(shape)+offset,0);

    // test constructor
    cout << "\tTest DataArrayView constructor." << endl;
    DataArrayView dataView(data,shape,offset);

    // test shape set correctly
    cout << "\tTest getShape." << endl;
    assert(dataView.getShape()==shape);

    // Assign values to the data
    cout << "\tAssign values via () operator." << endl;
    for (int i=0;i<shape[0];++i) {
      //cout << i << ":" << dataView.index(i) << endl;
      dataView(i)=dataView.index(i);
      assert(dataView(i)==dataView.index(i));
    }

    // test operator==
    cout << "\tTest == operator." << endl;
    assert(dataView==dataView);
  }

  {
    cout << endl;

    // define the shape for the data array view
    cout << "\tTest shape (2,3):" << endl;
    DataArrayView::ShapeType shape;
    shape.push_back(2);
    shape.push_back(3);

    // allocate the space external to the DataArrayView
    DataArrayView::ValueType data(DataArrayView::noValues(shape),0);

    // test constructor
    cout << "\tTest DataArrayView constructor." << endl;
    DataArrayView dataView(data,shape);

    // test shape set correctly
    cout << "\tTest getShape." << endl;
    assert(dataView.getShape()==shape);

    // Assign values to the data
    cout << "\tAssign values via () operator." << endl;
    for (int i=0;i<shape[0];++i) {
      for (int j=0;j<shape[1];++j) {
        //cout << i << "," << j << ":" << dataView.index(i,j) << endl;
	dataView(i,j)=dataView.index(i,j);
	assert(dataView(i,j)==dataView.index(i,j));
      }
    }

    // test operator==
    cout << "\tTest == operator." << endl;
    assert(dataView==dataView);
  }

  {
    cout << endl;

    // define the shape for the data array view
    cout << "\tTest shape (2,3,4):" << endl;
    DataArrayView::ShapeType shape;
    shape.push_back(2);
    shape.push_back(3);
    shape.push_back(4);

    // allocate the space external to the DataArrayView
    DataArrayView::ValueType data(DataArrayView::noValues(shape),0);

    // test constructor
    cout << "\tTest DataArrayView constructor." << endl;
    DataArrayView dataView(data,shape);

    // test shape set correctly
    cout << "\tTest getShape." << endl;
    assert(dataView.getShape()==shape);

    // Assign values to the data
    cout << "\tAssign values via () operator." << endl;
    for (int i=0;i<shape[0];++i) {
      for (int j=0;j<shape[1];++j) {
	for (int k=0;k<shape[2];++k) {
          //cout << i << "," << j << "," << k << ":" << dataView.index(i,j,k) << endl;
	  dataView(i,j,k)=dataView.index(i,j,k);
          assert(dataView(i,j,k)==dataView.index(i,j,k));
	}
      }
    }

    // test operator==
    cout << "\tTest == operator." << endl;
    assert(dataView==dataView);
  }

#if defined DOASSERT 
  {
    cout << endl;

    cout << "\tTest too many indices for rank of array exception." << endl;

    DataArrayView::ShapeType shape;
    DataArrayView::ValueType data(DataArrayView::noValues(shape),0);
    DataArrayView dataView(data,shape);

    // Should be a scalar
    dataView()=1;

    try {
      dataView(1)=1;
      assert(false);
    }
    catch (EsysException& e) {
      //cout << "\t" << e.what() << endl;
      assert(true);
    }

    try {
      dataView(1,1)=1;
      assert(false);
    }
    catch (EsysException& e) {
      //cout << "\t" << e.what() << endl;
      assert(true);
    }

    try {
      dataView(1,1,1)=1;
      assert(false);
    }
    catch (EsysException& e) {
      //cout << "\t" << e.what() << endl;
      assert(true);
    }
  }

  {
    cout << endl;

    cout << "\tTest invalid index exception." << endl;

    DataArrayView::ShapeType shape;
    shape.push_back(4);
    DataArrayView::ValueType data(DataArrayView::noValues(shape),0);
    DataArrayView dataView(data,shape);

    try {
      dataView(4000)=1;
      assert(false);
    }
    catch (EsysException& e) {
      //cout << "\t" << e.what() << endl;
      assert(true);
    }
  }

#endif

  {
    cout << endl;

    cout << "\tTest insufficient data exception." << endl;

    DataArrayView::ShapeType shape;
    DataArrayView::ValueType data;

    try {
      DataArrayView dataView(data,shape);
      assert(false);
    }
    catch (EsysException& e) {
      //cout << "\t" << e.what() << endl;
      assert(true);
    }
  }

  {
    cout << endl;

    cout << "\tTest matrix multiplication." << endl;

    DataArrayView::ShapeType leftShape;
    leftShape.push_back(1);
    leftShape.push_back(3);
    DataArrayView::ValueType leftData(DataArrayView::noValues(leftShape),0);
    DataArrayView leftDataView(leftData,leftShape);

    DataArrayView::ShapeType rightShape;
    rightShape.push_back(3);
    rightShape.push_back(2);
    DataArrayView::ValueType rightData(DataArrayView::noValues(rightShape),0);
    DataArrayView rightDataView(rightData,rightShape);

    DataArrayView::ShapeType resultShape=DataArrayView::determineResultShape(leftDataView,rightDataView);

    for (int i=0;i<resultShape.size();++i) {
      cout << "\tDimension: " << i << " size: " << resultShape[i] << endl;
    }

    DataArrayView::ValueType resultData(DataArrayView::noValues(resultShape),0);
    DataArrayView resultDataView(resultData,resultShape);
  
    cout << "\tTest result shape." << endl;

    assert(resultShape.size()==2);
    assert(resultShape[0]==1);
    assert(resultShape[1]==2);

    // Assign some values
    cout << "\tAssign values." << endl;
    double aValue=0.0;
    for (int i=0;i<leftShape[0];++i) {
      for (int j=0;j<leftShape[1];++j) {
	leftDataView(i,j)=++aValue;
      }
    }
    aValue=0.0;
    for (int i=0;i<rightShape[0];++i) {
      for (int j=0;j<rightShape[1];++j) {
	rightDataView(i,j)=++aValue;
      }
    }

    cout << "\tDo matrix multiplication." << endl;
    DataArrayView::matMult(leftDataView,rightDataView,resultDataView);

    cout << "\tCheck result." << endl;
    //need to hand build result matrix and compare with generated result here.

    /*
    cout << "Create a vector." << endl;
    DataArray v;
    DataArray::ShapeType vShape;
    vShape.push_back(3);
    v.setShape(vShape);
    double aValue=0.0;
    for (int i=0;i<vShape[0];++i) {
      v(i)=++aValue;
    }

    cout << "Create a matrix." << endl;
    DataArray mat;
    DataArray::ShapeType mShape;
    mShape.push_back(3);
    mShape.push_back(2);
    mat.setShape(mShape);
    aValue=0.0;
    for (int i=0;i<mShape[0];++i) {
      for (int j=0;j<mShape[1];++j) {
	mat(i,j)=++aValue;
      }
    }

    // [1,2,3] x |1, 2| = [22,28]
    //           |3, 4|
    //           |5, 6|

    cout << "Test multiplication of matrix and vector." << endl;
    DataArray result=DataArray::matMult(v,mat);
    assert(fabs(result(0) - 22) < 0.1);
    assert(fabs(result(1) - 28) < 0.1);
    */
  }  

}

TestSuite* DataArrayViewTestCase::suite ()
{
  //
  // create the suite of tests to perform.
  TestSuite *testSuite = new TestSuite ("DataArrayViewTestCase");

  testSuite->addTest (new TestCaller< DataArrayViewTestCase>("testAll",&DataArrayViewTestCase::testAll));
  testSuite->addTest (new TestCaller< DataArrayViewTestCase>("testShapeToString",&DataArrayViewTestCase::testShapeToString));
  testSuite->addTest (new TestCaller< DataArrayViewTestCase>("testScalarView",&DataArrayViewTestCase::testScalarView));
  testSuite->addTest (new TestCaller< DataArrayViewTestCase>("testSlicing",&DataArrayViewTestCase::testSlicing));
  return testSuite;
}

