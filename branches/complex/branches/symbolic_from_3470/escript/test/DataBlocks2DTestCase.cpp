
/*******************************************************
*
* Copyright (c) 2003-2012 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


#include "DataBlocks2DTestCase.h"
#include "escript/DataBlocks2D.h"
#include "esysUtils/EsysException.h"

#include <cppunit/TestCaller.h>
#include <iostream>

using namespace std;
using namespace CppUnit;
using namespace escript;
using namespace esysUtils;

void DataBlocks2DTestCase::testAll()
{
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
    CPPUNIT_ASSERT(myData.index(i,j) == (i*numCols+j)*blockSize);
    CPPUNIT_ASSERT(myData.size() == numRows*numCols*blockSize);
  }

  {
    cout << "\t\tnumRows = 3, numCols = 5, blockSize = 20." << endl;
    int numRows=3;
    int numCols=5;
    int blockSize=20;
    DataBlocks2D myData(numRows,numCols,blockSize);
    int i = numRows-1;
    int j = numCols-1;
    CPPUNIT_ASSERT(myData.index(i,j) == (i*numCols+j)*blockSize);
    CPPUNIT_ASSERT(myData.size() == numRows*numCols*blockSize);
  }

  {
    cout << "\t\tnumRows = 3, numCols = 5, blockSize = 1." << endl;
    int numRows=3;
    int numCols=5;
    int blockSize=1;
    DataBlocks2D myData(numRows,numCols,blockSize);
    int i = numRows-1;
    int j = numCols-1;
    CPPUNIT_ASSERT(myData.index(i,j) == (i*numCols+j)*blockSize);
    CPPUNIT_ASSERT(myData.size() == numRows*numCols*blockSize);
  }

  {
    cout << "\t\tnumRows = 1, numCols = 1, blockSize = 1." << endl;
    int numRows=1;
    int numCols=1;
    int blockSize=1;
    DataBlocks2D myData(numRows,numCols,blockSize);
    int i = numRows-1;
    int j = numCols-1;
    CPPUNIT_ASSERT(myData.index(i,j) == (i*numCols+j)*blockSize);
    CPPUNIT_ASSERT(myData.size() == numRows*numCols*blockSize);
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
	  CPPUNIT_ASSERT(myData[myData.index(i,j)+k] == val);
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
        CPPUNIT_FAIL("Exception not thrown");
    }
    catch(EsysException&) {
        CPPUNIT_ASSERT(true);
    }
  }

  {
    cout << "\tTest DataBlocks2D exception for numCols = 0." << endl;
    int numRows=10;
    int numCols=0;
    int blockSize=10;
    try {
        DataBlocks2D myData(numRows,numCols,blockSize);
        CPPUNIT_FAIL("Exception not thrown");
    }
    catch(EsysException&) {
        CPPUNIT_ASSERT(true);
    }
  }

  {
    cout << "\tTest DataBlocks2D exception for blockSize = 0." << endl;
    int numRows=10;
    int numCols=8;
    int blockSize=0;
    try {
        DataBlocks2D myData(numRows,numCols,blockSize);
        CPPUNIT_FAIL("Exception not thrown");
    }
    catch(EsysException&) {
        CPPUNIT_ASSERT(true);
    }
  }

  {
    cout << "\tTest getNumRows, getNumCols and getBlockSize." << endl;
    int numRows=1;
    int numCols=1;
    int blockSize=1;
    DataBlocks2D myData(numRows,numCols,blockSize);
    CPPUNIT_ASSERT(myData.getNumRows() == numRows);
    CPPUNIT_ASSERT(myData.getNumCols() == numCols);
    CPPUNIT_ASSERT(myData.getBlockSize() == blockSize);
  }

  {
    cout << "\tTest resize." << endl;
    int numRows=1;
    int numCols=1;
    int blockSize=1;
    DataBlocks2D myData;
    myData.resize(numRows,numCols,blockSize);
    CPPUNIT_ASSERT(myData.getNumRows() == numRows);
    CPPUNIT_ASSERT(myData.getNumCols() == numCols);
    CPPUNIT_ASSERT(myData.getBlockSize() == blockSize);
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
	CPPUNIT_ASSERT(myData1(i,j) == myData2(i,j));
      }
    }
  }

#if defined DOASSERT
  {
    cout << "\tTest DOASSERT exception." << endl;
    DataBlocks2D myData;
    CPPUNIT_ASSERT_THROW(myData.index(1,2), EsysException);
  }
#endif
}

TestSuite* DataBlocks2DTestCase::suite()
{
  TestSuite *testSuite = new TestSuite("DataBlocks2DTestCase");

  testSuite->addTest(new TestCaller<DataBlocks2DTestCase>(
              "testAll",&DataBlocks2DTestCase::testAll));
  return testSuite;
}

