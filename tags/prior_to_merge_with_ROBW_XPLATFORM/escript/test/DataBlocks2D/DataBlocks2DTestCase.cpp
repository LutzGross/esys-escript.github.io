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

#include "DataBlocks2D.h"
#include "EsysException.h"

#include "DataBlocks2DTestCase.h"

#include <iostream>

using namespace std;
using namespace CppUnitTest;
using namespace escript;
using namespace esysUtils;

void DataBlocks2DTestCase::setUp() {
  //
  // This is called before each test is run
 
}

void DataBlocks2DTestCase::tearDown() {
  //
  // This is called after each test has been run
 
}

void DataBlocks2DTestCase::testAll() {

  cout << endl;

  cout << "\tTest DataBlocks2D constructor for various dimension values:" << endl;

  {
    cout << "\t\tnumRows = 1, numCols = 1, blockSize = 20." << endl;
    int numRows=1;
    int numCols=1;
    int blockSize=20;
    DataBlocks2D myData(numRows,numCols,blockSize);
    int i = numRows-1;
    int j = numCols-1;
    assert(myData.index(i,j) == (i*numCols+j)*blockSize);
    assert(myData.size() == numRows*numCols*blockSize);
  }

  {
    cout << "\t\tnumRows = 3, numCols = 5, blockSize = 20." << endl;
    int numRows=3;
    int numCols=5;
    int blockSize=20;
    DataBlocks2D myData(numRows,numCols,blockSize);
    int i = numRows-1;
    int j = numCols-1;
    assert(myData.index(i,j) == (i*numCols+j)*blockSize);
    assert(myData.size() == numRows*numCols*blockSize);
  }

  {
    cout << "\t\tnumRows = 3, numCols = 5, blockSize = 1." << endl;
    int numRows=3;
    int numCols=5;
    int blockSize=1;
    DataBlocks2D myData(numRows,numCols,blockSize);
    int i = numRows-1;
    int j = numCols-1;
    assert(myData.index(i,j) == (i*numCols+j)*blockSize);
    assert(myData.size() == numRows*numCols*blockSize);
  }

  {
    cout << "\t\tnumRows = 1, numCols = 1, blockSize = 1." << endl;
    int numRows=1;
    int numCols=1;
    int blockSize=1;
    DataBlocks2D myData(numRows,numCols,blockSize);
    int i = numRows-1;
    int j = numCols-1;
    assert(myData.index(i,j) == (i*numCols+j)*blockSize);
    assert(myData.size() == numRows*numCols*blockSize);
  }

  {
    cout << "\tTest DataBlocks2D.index and DataBlocks2D operator[] for blockSize = 3." << endl;
    int numRows=10;
    int numCols=8;
    int blockSize=3;
    DataBlocks2D myData(numRows,numCols,blockSize);
    int val=0;
    for (int i=0; i<numRows; i++) {
      for (int j=0; j<numCols; j++) {
        for (int k=0; k<blockSize; k++) {
	  myData[myData.index(i,j)+k] = val;
          val++;
        }
      }
    }
    val=0;
    for (int i=0; i<numRows; i++) {
      for (int j=0; j<numCols; j++) {
        for (int k=0; k<blockSize; k++) {
	  assert(myData[myData.index(i,j)+k] == val);
          val++;
        }
      }
    }
  }

  {
    cout << "\tTest DataBlocks2D exception for numRows = 0." << endl;
    int numRows=0;
    int numCols=8;
    int blockSize=10;
    try {
        DataBlocks2D myData(numRows,numCols,blockSize);
        assert(false);
    }
    catch(EsysException& e) {
        assert(true);
    }
  }

  {
    cout << "\tTest DataBlocks2D exception for numCols = 0." << endl;
    int numRows=10;
    int numCols=0;
    int blockSize=10;
    try {
        DataBlocks2D myData(numRows,numCols,blockSize);
        assert(false);
    }
    catch(EsysException& e) {
        assert(true);
    }
  }

  {
    cout << "\tTest DataBlocks2D exception for blockSize = 0." << endl;
    int numRows=10;
    int numCols=8;
    int blockSize=0;
    try {
        DataBlocks2D myData(numRows,numCols,blockSize);
        assert(false);
    }
    catch(EsysException& e) {
        assert(true);
    }
  }

  {
    cout << "\tTest getNumRows, getNumCols and getBlockSize." << endl;
    int numRows=1;
    int numCols=1;
    int blockSize=1;
    DataBlocks2D myData(numRows,numCols,blockSize);
    assert(myData.getNumRows() == numRows);
    assert(myData.getNumCols() == numCols);
    assert(myData.getBlockSize() == blockSize);
  }

  {
    cout << "\tTest resize." << endl;
    int numRows=1;
    int numCols=1;
    int blockSize=1;
    DataBlocks2D myData;
    myData.resize(numRows,numCols,blockSize);
    assert(myData.getNumRows() == numRows);
    assert(myData.getNumCols() == numCols);
    assert(myData.getBlockSize() == blockSize);
  }

  {
    cout << "\tTest = operator, swap, and copy constructor." << endl;
    DataBlocks2D myData1;
    DataBlocks2D myData2(1, 1, 1);
    int val=0;
    for (int i=0; i<1; i++) {
      for (int j=0; j<1; j++) {
        for (int k=0; k<1; k++) {
	  myData2[myData2.index(i,j)+k] = val;
          val++;
        }
      }
    }
    myData1 = myData2;
    for (int i=0; i<myData1.getNumRows(); i++) {
      for (int j=0; j<myData1.getNumCols(); j++) {
	assert(myData1(i,j) == myData2(i,j));
      }
    }
  }

  {
    cout << "\tTest DOASSERT exception." << endl;
    DataBlocks2D myData;
    try {
      myData.index(1,2);
      assert(false);
    }
    catch (EsysException& e) {
      assert(true);
    }
  }

}

TestSuite* DataBlocks2DTestCase::suite ()
{
  //
  // create the suite of tests to perform.
  TestSuite *testSuite = new TestSuite ("DataBlocks2DTestCase");

  testSuite->addTest (new TestCaller< DataBlocks2DTestCase>("testAll",&DataBlocks2DTestCase::testAll));
  return testSuite;
}
