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
#include "escript/Data/DataBlocks2D.h"
#include "esysUtils/EsysException.h"

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

  {
    cout << "\tTest DOASSERT exception." << endl;
    DataBlocks2D myData;
    try {
      myData.index(1,2);
      assert(false);
    }
    catch (EsysException& e) {
      //std::cout << e.what() << std::endl;
      assert(true);
    }
  }

  cout << "\tTest DataBlocks2D constructor for various size arrays:" << endl;

  {
    cout << "\t\tnumRows = 1, numCols = 1, blockSize = 20." << endl;
    int numRows=1;
    int numCols=1;
    int blockSize=20;
    DataBlocks2D myData(numRows,numCols,blockSize);
    assert(myData.index(numRows-1,numCols-1) == ((numRows-1)*(numCols-1)*blockSize));
    assert(myData.getData().size() == numRows*numCols*blockSize);
  }

  {
    cout << "\t\tnumRows = 3, numCols = 5, blockSize = 20." << endl;
    int numRows=3;
    int numCols=5;
    int blockSize=20;
    DataBlocks2D myData(numRows,numCols,blockSize);
    assert(myData.index(numRows-1,numCols-1) == ((numRows*numCols-1)*blockSize));
    assert(myData.getData().size() == numRows*numCols*blockSize);
  }

  {
    cout << "\t\tnumRows = 3, numCols = 5, blockSize = 1." << endl;
    int numRows=3;
    int numCols=5;
    int blockSize=1;
    DataBlocks2D myData(numRows,numCols,blockSize);
    assert(myData.index(numRows-1,numCols-1) == ((numRows*numCols-1)*blockSize));
    assert(myData.getData().size() == numRows*numCols*blockSize);
  }

  {
    cout << "\t\tnumRows = 1, numCols = 1, blockSize = 1." << endl;
    int numRows=1;
    int numCols=1;
    int blockSize=1;
    DataBlocks2D myData(numRows,numCols,blockSize);
    assert(myData.index(numRows-1,numCols-1) == ((numRows*numCols-1)*blockSize));
    assert(myData.getData().size() == numRows*numCols*blockSize);
  }

  {
    cout << "\tTest DataBlocks2D.index and DataBlocks2D.getData for blockSize = 3." << endl;
    int numRows=10;
    int numCols=8;
    int blockSize=3;
    DataBlocks2D myData(numRows,numCols,blockSize);
    for (int iR=0;iR<numRows;++iR) {
      for (int iC=0;iC<numCols;++iC) {
        //cout << "i=" << iR << " j=" << iC << " index=" << myData.index(iR,iC) << endl;
	assert(myData.getData()[myData.index(iR,iC)] == 0);
      }
    }
  }

  {
    cout << "\tTest DataBlocks2D.getData for blockSize = 0 - will generate exceptions." << endl;
    //
    // every attempted access for a 0 size block should cause an exception
    int numRows=10;
    int numCols=8;
    int blockSize=0;
    int exceptionCount=0;
    DataBlocks2D myData(numRows,numCols,blockSize);
    for (int iR=0;iR<numRows;++iR) {
      for (int iC=0;iC<numCols;++iC) {
	try {
	  myData.getData()[myData.index(iR,iC)];
          assert(false);
	}
	catch(EsysException& e) {
          //cout << e.toString() << endl;
	  ++exceptionCount;
          assert(true);
	}
      }
    }
    assert(exceptionCount == numRows*numCols);
  }

  {
    cout << "\tTest getNumRows and getNumCols." << endl;
    int numRows=1;
    int numCols=1;
    int blockSize=1;
    DataBlocks2D myData(numRows,numCols,blockSize);
    assert(myData.getNumRows() == numRows);
    assert(myData.getNumCols() == numCols);
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
  }

  {
    cout << "\tTest = operator, swap, and copy constructor." << endl;
    DataBlocks2D myData1;
    DataBlocks2D myData2(1, 1, 1);
    myData1 = myData2;
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
