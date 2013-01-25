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

#include "DataArrayView.h"
#include "DataArray.h"
#include "EsysException.h"

#include "DataArrayTestCase.h"

#include <iostream>
#include <vector>

using namespace CppUnitTest;
using namespace esysUtils;
using namespace escript;
using namespace std;

void DataArrayTestCase::setUp() {
  //
  // This is called before each test is run
 
}

void DataArrayTestCase::tearDown() {
  //
  // This is called after each test has been run
 
}

void DataArrayTestCase::testAll() {
  //
  // The test code may be entered here
  // There is nothing special about the function name, it may be renamed to
  // something more suitable. 
  // As many test methods as desired may be added.

  cout << endl;

  {
    cout << endl;
    cout << "\tTest default DataArray constructor - default value." << endl;

    cout << "\tConstruct Default DataArray." << endl;
    DataArrayView::ShapeType newShape;
    DataArray newArray;
    DataArrayView newView = newArray.getView();

    cout << "\tCheck DataArray view shape." << endl;
    assert(newView.getShape()==newShape);
    assert(newView.getRank()==0);
    assert(newView.noValues()==1);

    cout << "\tCheck DataArray data vector values." << endl;
    DataArrayView::ValueType arrayData = newArray.getData();
    DataArrayView::ValueType viewData = newView.getData();
    assert(viewData==arrayData);
    int noValues = newView.noValues();
    for (int i=0; i<noValues; i++) {
      assert(arrayData[i]==0.0);
    }
    assert(newView()==0.0);
  }

  {
    cout << endl;
    cout << "\tTest default DataArray constructor - value 5.0." << endl;

    cout << "\tConstruct Default DataArray." << endl;
    DataArrayView::ShapeType newShape;
    DataArray newArray(5.0);
    DataArrayView newView = newArray.getView();

    cout << "\tCheck DataArray view shape." << endl;
    assert(newView.getShape()==newShape);
    assert(newView.getRank()==0);
    assert(newView.noValues()==1);

    cout << "\tCheck DataArray data vector values." << endl;
    DataArrayView::ValueType arrayData = newArray.getData();
    DataArrayView::ValueType viewData = newView.getData();
    assert(viewData==arrayData);
    int noValues = newView.noValues();
    for (int i=0; i<noValues; i++) {
      assert(arrayData[i]==5.0);
    }
    assert(newView()==5.0);
  }

  {
    cout << endl;
    cout << "\tTest DataArray shape constructor, shape() - default value." << endl;

    cout << "\tConstruct new DataArray with rank 0 view." << endl;
    DataArrayView::ShapeType newShape;
    DataArray newArray(newShape);
    DataArrayView newView = newArray.getView();

    cout << "\tCheck DataArray view shape." << endl;
    assert(newView.getShape()==newShape);
    assert(newView.getRank()==0);
    assert(newView.noValues()==1);

    cout << "\tCheck DataArray data vector values." << endl;
    DataArrayView::ValueType arrayData = newArray.getData();
    DataArrayView::ValueType viewData = newView.getData();
    assert(viewData==arrayData);
    int noValues = newView.noValues();
    for (int i=0; i<noValues; i++) {
      assert(arrayData[i]==0.0);
    }
    assert(newView()==0.0);
  }

  {
    cout << endl;
    cout << "\tTest DataArray shape constructor, shape() - value 5.0." << endl;

    cout << "\tConstruct new DataArray with rank 0 view." << endl;
    DataArrayView::ShapeType newShape;
    DataArray newArray(newShape, 5.0);
    DataArrayView newView = newArray.getView();

    cout << "\tCheck DataArray view shape." << endl;
    assert(newView.getShape()==newShape);
    assert(newView.getRank()==0);
    assert(newView.noValues()==1);

    cout << "\tCheck DataArray data vector values." << endl;
    DataArrayView::ValueType arrayData = newArray.getData();
    DataArrayView::ValueType viewData = newView.getData();
    assert(viewData==arrayData);
    int noValues = newView.noValues();
    for (int i=0; i<noValues; i++) {
      assert(arrayData[i]==5.0);
    }
    assert(newView()==5.0);
  }

  {
    cout << endl;
    cout << "\tTest DataArray shape constructor, shape(5)." << endl;

    cout << "\tConstruct new DataArray with rank 1 view." << endl;
    DataArrayView::ShapeType newShape;
    newShape.push_back(5);
    DataArray newArray(newShape);
    DataArrayView newView = newArray.getView();

    cout << "\tCheck DataArray view shape." << endl;
    assert(newView.getShape()==newShape);
    assert(newView.getRank()==1);
    assert(newView.noValues()==5);

    cout << "\tCheck DataArray data vector values." << endl;
    DataArrayView::ValueType::size_type index;
    DataArrayView::ValueType arrayData = newArray.getData();
    DataArrayView::ValueType viewData = newView.getData();
    assert(viewData==arrayData);
    for (int i=0; i<newShape[0]; i++) {
      index=newView.index(i);
      assert(arrayData[index]==0.0);
      assert(newView(i)==0.0);
    }
  }
  
  {
    cout << endl;
    cout << "\tTest DataArray shape constructor, shape(2,3)." << endl;

    cout << "\tConstruct new DataArray with rank 2 view." << endl;
    DataArrayView::ShapeType newShape;
    newShape.push_back(2);
    newShape.push_back(3);
    DataArray newArray(newShape);
    DataArrayView newView = newArray.getView();

    cout << "\tCheck DataArray view shape." << endl;
    assert(newView.getShape()==newShape);
    assert(newView.getRank()==2);
    assert(newView.noValues()==6);

    cout << "\tCheck DataArray data vector values." << endl;
    DataArrayView::ValueType::size_type index;
    DataArrayView::ValueType arrayData = newArray.getData();
    DataArrayView::ValueType viewData = newView.getData();
    assert(viewData==arrayData);
    for (int i=0; i<newShape[0]; i++) {
      for (int j=0; j<newShape[1]; j++) {
        index=newView.index(i,j);
        assert(arrayData[index]==0.0);
        assert(newView(i,j)==0.0);
      }
    }
  }
  
  {
    cout << endl;
    cout << "\tTest DataArray shape constructor, shape(2,3,4)." << endl;

    cout << "\tConstruct new DataArray with rank 3 view." << endl;
    DataArrayView::ShapeType newShape;
    newShape.push_back(2);
    newShape.push_back(3);
    newShape.push_back(4);
    DataArray newArray(newShape);
    DataArrayView newView = newArray.getView();

    cout << "\tCheck DataArray view shape." << endl;
    assert(newView.getShape()==newShape);
    assert(newView.getRank()==3);
    assert(newView.noValues()==24);

    cout << "\tCheck DataArray data vector values." << endl;
    DataArrayView::ValueType::size_type index;
    DataArrayView::ValueType viewData = newView.getData();
    DataArrayView::ValueType arrayData = newArray.getData();
    assert(viewData==arrayData);
    for (int i=0; i<newShape[0]; i++) {
      for (int j=0; j<newShape[1]; j++) {
        for (int k=0; k<newShape[2]; k++) {
          index=newView.index(i,j,k);
          assert(arrayData[index]==0.0);
          assert(newView(i,j,k)==0.0);
        }
      }
    }
  }
  
  {
    cout << endl;
    cout << "\tTest DataArray shape constructor, shape(2,3,4,7)." << endl;

    cout << "\tConstruct new DataArray with rank 4 view." << endl;
    DataArrayView::ShapeType newShape;
    newShape.push_back(2);
    newShape.push_back(3);
    newShape.push_back(4);
    newShape.push_back(7);
    DataArray newArray(newShape);
    DataArrayView newView = newArray.getView();

    cout << "\tCheck DataArray view shape." << endl;
    assert(newView.getShape()==newShape);
    assert(newView.getRank()==4);
    assert(newView.noValues()==168);

    cout << "\tCheck DataArray data vector values." << endl;
    DataArrayView::ValueType::size_type index;
    DataArrayView::ValueType viewData = newView.getData();
    DataArrayView::ValueType arrayData = newArray.getData();
    assert(viewData==arrayData);
    for (int i=0; i<newShape[0]; i++) {
      for (int j=0; j<newShape[1]; j++) {
        for (int k=0; k<newShape[2]; k++) {
          for (int l=0; l<newShape[3]; l++) {
            index=newView.index(i,j,k,l);
            assert(arrayData[index]==0.0);
            assert(newView(i,j,k,l)==0.0);
          }
        }
      }
    }
  }

  {
    cout << endl;
    cout << "\tTest DataArray copy constructor, default DataArray." << endl;

    cout << "\tConstruct Default DataArray." << endl;
    DataArray oldArray;

    cout << "\tCopy Default DataArray." << endl;
    DataArrayView::ShapeType newShape;
    DataArray newArray(oldArray);
    DataArrayView newView = newArray.getView();

    cout << "\tCheck DataArray view shape." << endl;
    assert(newView.getShape()==newShape);
    assert(newView.getRank()==0);
    assert(newView.noValues()==1);

    cout << "\tCheck DataArray data vector values." << endl;
    DataArrayView::ValueType arrayData = newArray.getData();
    DataArrayView::ValueType viewData = newView.getData();
    assert(viewData==arrayData);
    assert(arrayData[0]==0.0);
    assert(newView()==0.0);
  }
 
  {
    cout << endl;
    cout << "\tTest DataArray copy constructor, shape(2,3)." << endl;

    cout << "\tConstruct new DataArray with rank 2 view." << endl;
    DataArrayView::ShapeType oldShape;
    oldShape.push_back(2);
    oldShape.push_back(3);
    DataArray oldArray(oldShape);

    cout << "\tCopy rank 2 DataArray." << endl;
    DataArrayView::ShapeType newShape = oldShape;
    DataArray newArray(oldArray);
    DataArrayView newView = newArray.getView();

    cout << "\tCheck DataArray view shape." << endl;
    assert(newView.getShape()==newShape);
    assert(newView.getRank()==2);
    assert(newView.noValues()==6);

    cout << "\tCheck DataArray data vector values." << endl;
    DataArrayView::ValueType::size_type index;
    DataArrayView::ValueType viewData = newView.getData();
    DataArrayView::ValueType arrayData = newArray.getData();
    assert(viewData==arrayData);
    for (int i=0; i<newShape[0]; i++) {
      for (int j=0; j<newShape[1]; j++) {
        index=newView.index(i,j);
        assert(arrayData[index]==0.0);
        assert(newView(i,j)==0.0);
      }
    }
  }

  {
    cout << endl;
    cout << "\tTest DataArray copy constructor, shape(2,3,4,7)." << endl;

    cout << "\tConstruct new DataArray with rank 4 view." << endl;
    DataArrayView::ShapeType oldShape;
    oldShape.push_back(2);
    oldShape.push_back(3);
    oldShape.push_back(4);
    oldShape.push_back(7);
    DataArray oldArray(oldShape);

    cout << "\tCopy rank 4 DataArray." << endl;
    DataArrayView::ShapeType newShape = oldShape;
    DataArray newArray(oldArray);
    DataArrayView newView = newArray.getView();

    cout << "\tCheck DataArray view shape." << endl;
    assert(newView.getShape()==newShape);
    assert(newView.getRank()==4);
    assert(newView.noValues()==168);

    cout << "\tCheck DataArray data vector values." << endl;
    DataArrayView::ValueType::size_type index;
    DataArrayView::ValueType viewData = newView.getData();
    DataArrayView::ValueType arrayData = newArray.getData();
    assert(viewData==arrayData);
    for (int i=0; i<newShape[0]; i++) {
      for (int j=0; j<newShape[1]; j++) {
        for (int k=0; k<newShape[2]; k++) {
          for (int l=0; l<newShape[3]; l++) {
            index=newView.index(i,j,k,l);
            assert(arrayData[index]==0.0);
            assert(newView(i,j,k,l)==0.0);
          }
        }
      }
    }
  }

  {
    cout << endl;
    cout << "\tTest DataArray DataArrayView constructor, default DataArray." << endl;

    cout << "\tConstruct Default DataArray." << endl;
    DataArray oldArray;

    cout << "\tExtract Default DataArray DataArrayView." << endl;
    DataArrayView oldView = oldArray.getView();

    cout << "\tCopy Default DataArray from DataArrayView." << endl;
    DataArrayView::ShapeType newShape;
    DataArray newArray(oldView);
    DataArrayView newView = newArray.getView();

    cout << "\tCheck DataArray view shape." << endl;
    assert(newView.getShape()==newShape);
    assert(newView.getRank()==0);
    assert(newView.noValues()==1);

    cout << "\tCheck DataArray data vector values." << endl;
    DataArrayView::ValueType arrayData = newArray.getData();
    DataArrayView::ValueType viewData = newView.getData();
    assert(viewData==arrayData);
    assert(arrayData[0]==0.0);
    assert(newView()==0.0);
  }
 
  {
    cout << endl;
    cout << "\tTest DataArray DataArrayView constructor, shape(2,3)." << endl;

    cout << "\tConstruct new DataArray with rank 2 view." << endl;
    DataArrayView::ShapeType oldShape;
    oldShape.push_back(2);
    oldShape.push_back(3);
    DataArray oldArray(oldShape);

    cout << "\tExtract rank 2 DataArray DataArrayView." << endl;
    DataArrayView oldView = oldArray.getView();

    cout << "\tCopy rank 2 DataArray from DataArrayView." << endl;
    DataArrayView::ShapeType newShape = oldShape;
    DataArray newArray(oldView);
    DataArrayView newView = newArray.getView();

    cout << "\tCheck DataArray view shape." << endl;
    assert(newView.getShape()==newShape);
    assert(newView.getRank()==2);
    assert(newView.noValues()==6);

    cout << "\tCheck DataArray data vector values." << endl;
    DataArrayView::ValueType::size_type index;
    DataArrayView::ValueType viewData = newView.getData();
    DataArrayView::ValueType arrayData = newArray.getData();
    assert(viewData==arrayData);
    for (int i=0; i<newShape[0]; i++) {
      for (int j=0; j<newShape[1]; j++) {
        index=newView.index(i,j);
        assert(arrayData[index]==0.0);
        assert(newView(i,j)==0.0);
      }
    }
  }

  {
    cout << endl;
    cout << "\tTest DataArray DataArrayView constructor, shape(2,3,4,7)." << endl;

    cout << "\tConstruct new DataArray with rank 4 view." << endl;
    DataArrayView::ShapeType oldShape;
    oldShape.push_back(2);
    oldShape.push_back(3);
    oldShape.push_back(4);
    oldShape.push_back(7);
    DataArray oldArray(oldShape);

    cout << "\tExtract rank 4 DataArray DataArrayView." << endl;
    DataArrayView oldView = oldArray.getView();

    cout << "\tCopy rank 4 DataArray from DataArrayView." << endl;
    DataArrayView::ShapeType newShape = oldShape;
    DataArray newArray(oldView);
    DataArrayView newView = newArray.getView();

    cout << "\tCheck DataArray view shape." << endl;
    assert(newView.getShape()==newShape);
    assert(newView.getRank()==4);
    assert(newView.noValues()==168);

    cout << "\tCheck DataArray data vector values." << endl;
    DataArrayView::ValueType::size_type index;
    DataArrayView::ValueType viewData = newView.getData();
    DataArrayView::ValueType arrayData = newArray.getData();
    assert(viewData==arrayData);
    for (int i=0; i<newShape[0]; i++) {
      for (int j=0; j<newShape[1]; j++) {
        for (int k=0; k<newShape[2]; k++) {
          for (int l=0; l<newShape[3]; l++) {
            index=newView.index(i,j,k,l);
            assert(arrayData[index]==0.0);
            assert(newView(i,j,k,l)==0.0);
          }
        }
      }
    }
  }

}

TestSuite* DataArrayTestCase::suite ()
{
  //
  // create the suite of tests to perform.
  TestSuite *testSuite = new TestSuite ("DataArrayTestCase");

  testSuite->addTest (new TestCaller< DataArrayTestCase>("testAll",&DataArrayTestCase::testAll));
  return testSuite;
}
