
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

#include "DataTypesTestCase.h"

#include <escript/DataVector.h>
#include <escript/EsysException.h>

#include <cppunit/TestCaller.h>
#include <iostream>

using namespace CppUnit;
using namespace escript;
using namespace escript::DataTypes;
using namespace std;


void DataTypesTestCase::testShapeFns()
{
  cout << "\n\tTest Shape functions." << endl;
  ShapeType shape;
  shape.push_back(2);
  ShapeType s1=shape;
  CPPUNIT_ASSERT(checkShape(s1,shape));
  shape.push_back(3);
  ShapeType s2=shape;
  CPPUNIT_ASSERT(checkShape(s2,shape));
  shape.push_back(4);
  ShapeType s3=shape;
  CPPUNIT_ASSERT(checkShape(s3,shape));
  ShapeType s4=shape;
  s4.push_back(5);

  cout << "\t\tNumber of Values." << endl;
  CPPUNIT_ASSERT(noValues(shape)==24);
  CPPUNIT_ASSERT(noValues(scalarShape)==1);
  cout << "\t\tGet Rank." << endl;
  CPPUNIT_ASSERT(getRank(scalarShape)==0);
  CPPUNIT_ASSERT(getRank(shape)==3);


  cout << "\t\tNumber of Values for RegionLoopRangeType." << endl;
  RegionLoopRangeType rlr;
  CPPUNIT_ASSERT(noValues(rlr)==1);
  rlr.push_back(std::pair<int,int>(0,2));
  CPPUNIT_ASSERT(noValues(rlr)==2);
  rlr.push_back(std::pair<int,int>(1,4));
  CPPUNIT_ASSERT(noValues(rlr)==6);
  rlr.push_back(std::pair<int,int>(2,7));
  CPPUNIT_ASSERT(noValues(rlr)==30);
  cout << "\t\tRelative index methods.\n";
  CPPUNIT_ASSERT(getRelIndex(s1,1)==1);
  CPPUNIT_ASSERT(getRelIndex(s2,1,2)==5);
  CPPUNIT_ASSERT(getRelIndex(shape,1,1,1)==9);
  CPPUNIT_ASSERT(getRelIndex(s4,2,1,1,1)==34);
  cout << "\t\tCreateShapeErrorMessage.\n";
  CPPUNIT_ASSERT(createShapeErrorMessage("prefix",s1,s2)==string("prefix This shape: (2,3) Other shape: (2)"));

  cout << "\t\tgetSliceRegionLoopRange." << endl;
  RegionType r;
  r.push_back(std::pair<int,int>(1,1));
  r.push_back(std::pair<int,int>(0,3));
  r.push_back(std::pair<int,int>(0,3));
  RegionLoopRangeType rl;
  rl.push_back(std::pair<int,int>(1,2));
  rl.push_back(std::pair<int,int>(0,3));
  rl.push_back(std::pair<int,int>(0,3));
  RegionLoopRangeType rt=getSliceRegionLoopRange(r);
  CPPUNIT_ASSERT(rt==rl);


/*  
  
#ifdef DOASSERT
// The errors we are testing for are triggered by ESysAssert which is only defined when DOASSERT is.

  cout << "\t\tInvalid index.(too many)" << endl;
  // test too many indices
  CPPUNIT_ASSERT_THROW(getRelIndex(s1,1,1), EsysException);
  CPPUNIT_ASSERT_THROW(getRelIndex(s2,1,1,1), EsysException);
  CPPUNIT_ASSERT_THROW(getRelIndex(s3,1,1,1,1), EsysException);
  // too few indices
  cout << "\t\tInvalid index.(too few)" << endl;
  CPPUNIT_ASSERT_THROW(getRelIndex(s2,1), EsysException);
  CPPUNIT_ASSERT_THROW(getRelIndex(s3,1,1), EsysException);
  // indices too large
  cout << "\t\tInvalid index.(too large for shape)" << endl;
  CPPUNIT_ASSERT_THROW(getRelIndex(s1,10), EsysException);
  CPPUNIT_ASSERT_THROW(getRelIndex(s3,2,4,4), EsysException);
#endif
*/

}

void DataTypesTestCase::testResultSliceShape() {

  cout << endl;
  cout << "\tTest getResultSliceShape method." << endl;

  DataTypes::RegionType region;
  DataTypes::ShapeType resultShape;

  region.push_back(DataTypes::RegionType::value_type(1,5));
  resultShape.push_back(4);
  CPPUNIT_ASSERT(DataTypes::getResultSliceShape(region)==resultShape);

  region.push_back(DataTypes::RegionType::value_type(2,5));
  resultShape.push_back(3);
  CPPUNIT_ASSERT(DataTypes::getResultSliceShape(region)==resultShape);

  region.push_back(DataTypes::RegionType::value_type(3,9));
  resultShape.push_back(6);
  CPPUNIT_ASSERT(DataTypes::getResultSliceShape(region)==resultShape);

  region.push_back(DataTypes::RegionType::value_type(1,7));
  resultShape.push_back(6);
  CPPUNIT_ASSERT(DataTypes::getResultSliceShape(region)==resultShape);

}

void DataTypesTestCase::testSlicing() {

  using namespace DataTypes;
  {

    cout << endl;
    cout << "\tSlice a scalar to a scalar.";

    // Define slice region.
    DataTypes::RegionType region;

    // Define shape of views.
    DataTypes::ShapeType sourceShape;

    // Create source and target views.
    int len = DataTypes::noValues(sourceShape);
    RealVectorType sourceData(len, 2.0, len);
//    DataArrayView sourceView(sourceData, sourceShape);
    RealVectorType targetData(1, 2.0, 1);
//    DataArrayView targetView(targetData, DataTypes::ShapeType());

    // Copy source view to target view.
//    targetView.copySlice(sourceView,region);


   DataTypes::copySlice(targetData, DataTypes::scalarShape, 0, sourceData, sourceShape, 0,region); 

    // Check results of copy.
//     CPPUNIT_ASSERT(sourceView==targetView);
    CPPUNIT_ASSERT(targetData==sourceData);
  }

  {
    cout << endl;
    cout << "\tSlice a scalar from a scalar.";

    // Define slice region.
    DataTypes::RegionType region;

    // Define shape of views.
    DataTypes::ShapeType sourceShape;

    // Create source and target views.
    int len = DataTypes::noValues(sourceShape);
    RealVectorType sourceData(len, 2.0, len);
//     DataArrayView sourceView(sourceData, sourceShape);
    RealVectorType targetData(1, 2.0, 1);
//     DataArrayView targetView(targetData, DataTypes::ShapeType());

    // Copy source view to target view.
    DataTypes::copySliceFrom(targetData, DataTypes::scalarShape, 0, sourceData, sourceShape, 0,region);

    // Check results of copy.
//     CPPUNIT_ASSERT(sourceView==targetView);
    CPPUNIT_ASSERT(sourceData==targetData);
  }

  {
    cout << endl;
    cout << "\tSlice a 1 dimensional slice to a 1 dimensional array.";

    // Define slice region.
    DataTypes::RegionType region;
    region.push_back(DataTypes::RegionType::value_type(2,4));

    // Define shapes of views.
    DataTypes::ShapeType sourceShape;
    sourceShape.push_back(6);
    DataTypes::ShapeType targetShape = DataTypes::getResultSliceShape(region);

    // Create source and target views.
    int len = DataTypes::noValues(sourceShape);
    RealVectorType sourceData(len, 2.0, len);
    //DataArrayView sourceView(sourceData, sourceShape);
    for (int i=0;i<sourceShape[0];i++) {
      sourceData[i]=i;
    }

    len = DataTypes::noValues(targetShape);
    RealVectorType targetData(len, 2.0, len);
    //DataArrayView targetView(targetData, targetShape);

    // Copy source view to target view.
//     targetView.copySlice(sourceView,region);
    DataTypes::copySlice(targetData, targetShape,0,sourceData, sourceShape,0,region);

    // Check results of copy.
    for (int i=region[0].first;i<region[0].second;i++) {
//       CPPUNIT_ASSERT(sourceView(i)==
//              targetView(i-region[0].first));
	CPPUNIT_ASSERT(sourceData[i]==targetData[i-region[0].first]);
    }
  }

  {
    cout << endl;
    cout << "\tSlice a 1 dimensional slice to a scalar.";

    // Define slice region.
    DataTypes::RegionType region;
    region.push_back(DataTypes::RegionType::value_type(2,3));

    // Define shapes of views.
    DataTypes::ShapeType sourceShape;
    sourceShape.push_back(6);
    DataTypes::ShapeType targetShape;

    // Create source and target views.
    int len = DataTypes::noValues(sourceShape);
    RealVectorType sourceData(len, 2.0, len);
//     DataArrayView sourceView(sourceData, sourceShape);

    for (int i=0;i<sourceShape[0];i++) {
      sourceData[i]=i;
    }

    len = DataTypes::noValues(targetShape);
    RealVectorType targetData(len, 2.0, len);
//     DataArrayView targetView(targetData, targetShape);

    // Copy source view to target view.
//     targetView.copySlice(sourceView,region);
    DataTypes::copySlice(targetData, targetShape, 0, sourceData, sourceShape,0,region);

    // Check results of copy.
    for (int i=region[0].first;i<region[0].second;i++) {
/*      CPPUNIT_ASSERT(sourceView(i)==
             targetView());*/
      CPPUNIT_ASSERT(sourceData[i]==
             targetData[0]);

    }
  }

  {
    cout << endl;
    cout << "\tSlice a 1 dimensional slice from a 1 dimensional array.";

    // Define slice region.
    DataTypes::RegionType region;
    region.push_back(DataTypes::RegionType::value_type(2,4));

    // Define shapes of views.
    DataTypes::ShapeType sourceShape = DataTypes::getResultSliceShape(region);
    DataTypes::ShapeType targetShape;
    targetShape.push_back(6);

    // Create source and target views.
    int len = DataTypes::noValues(sourceShape);
    RealVectorType sourceData(len, 2.0, len);
//     DataArrayView sourceView(sourceData, sourceShape);
    for (int i=0;i<sourceShape[0];i++) {
      sourceData[i]=i;
    }

    len = DataTypes::noValues(targetShape);
    RealVectorType targetData(len, 2.0, len);
//     DataArrayView targetView(targetData, targetShape);

    // Copy source view to target view.
//     targetView.copySliceFrom(sourceView,region);
    DataTypes::copySliceFrom(targetData,targetShape,0,sourceData, sourceShape, 0,region);

    // Check results of copy.
    for (int i=region[0].first;i<region[0].second;i++) {
/*      CPPUNIT_ASSERT(sourceView(i-region[0].first)==
             targetView(i));*/
      CPPUNIT_ASSERT(sourceData[i-region[0].first]==
             targetData[i]);
    }
  }

  {
    cout << endl;
    cout << "\tSlice a 1 dimensional slice from a scalar.";

    // Define slice region.
    DataTypes::RegionType region;
    region.push_back(DataTypes::RegionType::value_type(2,4));

    // Define shapes of views.
    DataTypes::ShapeType sourceShape;
    DataTypes::ShapeType targetShape;
    targetShape.push_back(6);

    // Create source and target views.
    int len = DataTypes::noValues(sourceShape);
    RealVectorType sourceData(len, 2.0, len);
//     DataArrayView sourceView(sourceData, sourceShape);
    sourceData[0]=5;

    len = DataTypes::noValues(targetShape);
    RealVectorType targetData(len, 2.0, len);
//     DataArrayView targetView(targetData, targetShape);

    // Copy source view to target view.
//     targetView.copySliceFrom(sourceView,region);
    DataTypes::copySliceFrom(targetData,targetShape,0,sourceData, sourceShape,0,region);	

    // Check results of copy.
    for (int i=region[0].first;i<region[0].second;i++) {
      CPPUNIT_ASSERT(sourceData[0]==
             targetData[i]);
    }
  }

  {
    cout << endl;
    cout << "\tSlice a 2 dimensional slice to a 2 dimensional array.";

    // Define slice region.
    DataTypes::RegionType region;
    region.push_back(DataTypes::RegionType::value_type(2,4));
    region.push_back(DataTypes::RegionType::value_type(0,2));

    // Define shapes of views.
    DataTypes::ShapeType sourceShape;
    sourceShape.push_back(6);
    sourceShape.push_back(3);
    DataTypes::ShapeType targetShape = DataTypes::getResultSliceShape(region);

    // Create source and target views.
    int len = DataTypes::noValues(sourceShape);
    RealVectorType sourceData(len, 2.0, len);
//     DataArrayView sourceView(sourceData, sourceShape);
    int val=0;
    for (int i=0;i<sourceShape[0];i++) {
      for (int j=0;j<sourceShape[1];j++) {
	sourceData[getRelIndex(sourceShape,i,j)]=val++;
      }
    }

    len = DataTypes::noValues(targetShape);
    RealVectorType targetData(len, 2.0, len);
//     DataArrayView targetView(targetData, targetShape);

    // Copy source view to target view.
    copySlice(targetData, targetShape,0,sourceData, sourceShape,0,region);

    // Check results of copy.
    for (int i=region[0].first;i<region[0].second;i++) {
      for (int j=region[1].first;j<region[1].second;j++) {
	CPPUNIT_ASSERT(sourceData[getRelIndex(sourceShape,i,j)]==
               targetData[getRelIndex(targetShape,i-region[0].first,j-region[1].first)]);
      }
    }
  }

  {
    cout << endl;
    cout << "\tSlice a 2 dimensional slice to a 1 dimensional array.";

    // Define slice region.
    DataTypes::RegionType region;
    region.push_back(DataTypes::RegionType::value_type(0,3));
    region.push_back(DataTypes::RegionType::value_type(1,2));

    // Define shapes of views.
    DataTypes::ShapeType sourceShape;
    sourceShape.push_back(3);
    sourceShape.push_back(6);
    DataTypes::ShapeType targetShape;
    targetShape.push_back(3);

    // Create source and target views.
    int len = DataTypes::noValues(sourceShape);
    RealVectorType sourceData(len, 2.0, len);
//     DataArrayView sourceView(sourceData, sourceShape);
    int val=0;
    for (int i=0;i<sourceShape[0];i++) {
      for (int j=0;j<sourceShape[1];j++) {
	sourceData[getRelIndex(sourceShape,i,j)]=val++;
      }
    }

    len = DataTypes::noValues(targetShape);
    RealVectorType targetData(len, 2.0, len);
//     DataArrayView targetView(targetData, targetShape);

    // Copy source view to target view.
//     targetView.copySlice(sourceView,region);
    copySlice(targetData, targetShape,0,sourceData,sourceShape,0,region);

    // Check results of copy.
    for (int i=region[0].first;i<region[0].second;i++) {
      for (int j=region[1].first;j<region[1].second;j++) {
	CPPUNIT_ASSERT(sourceData[getRelIndex(sourceShape,i,j)]==
               targetData[getRelIndex(targetShape,i-region[0].first)]);
      }
    }
  }

  {
    cout << endl;
    cout << "\tSlice a 2 dimensional slice to a scalar.";

    // Define slice region.
    DataTypes::RegionType region;
    region.push_back(DataTypes::RegionType::value_type(2,3));
    region.push_back(DataTypes::RegionType::value_type(1,2));

    // Define shapes of views.
    DataTypes::ShapeType sourceShape;
    sourceShape.push_back(3);
    sourceShape.push_back(6);
    DataTypes::ShapeType targetShape;

    // Create source and target views.
    int len = DataTypes::noValues(sourceShape);
    RealVectorType sourceData(len, 2.0, len);
//     DataArrayView sourceView(sourceData, sourceShape);
    int val=0;
    for (int i=0;i<sourceShape[0];i++) {
      for (int j=0;j<sourceShape[1];j++) {
	sourceData[getRelIndex(sourceShape,i,j)]=val++;
      }
    }

    len = DataTypes::noValues(targetShape);
    RealVectorType targetData(len, 2.0, len);
//     DataArrayView targetView(targetData, targetShape);

    // Copy source view to target view.
//     targetView.copySlice(sourceView,region);
    copySlice(targetData, targetShape,0,sourceData, sourceShape,0,region);


    // Check results of copy.
    for (int i=region[0].first;i<region[0].second;i++) {
      for (int j=region[1].first;j<region[1].second;j++) {
/*	CPPUNIT_ASSERT(sourceView(i,j)==
               targetView());*/
	CPPUNIT_ASSERT(sourceData[getRelIndex(sourceShape,i,j)]==
               targetData[0]);
      }
    }
  }

  {
    cout << endl;
    cout << "\tSlice a 2 dimensional slice from a 2 dimensional array.";

    // Define slice region.
    DataTypes::RegionType region;
    region.push_back(DataTypes::RegionType::value_type(2,4));
    region.push_back(DataTypes::RegionType::value_type(0,2));

    // Define shapes of views.
    DataTypes::ShapeType sourceShape = DataTypes::getResultSliceShape(region);
    DataTypes::ShapeType targetShape;
    targetShape.push_back(6);
    targetShape.push_back(3);

    // Create source and target views.
    int len = DataTypes::noValues(sourceShape);
    RealVectorType sourceData(len, 2.0, len);
//     DataArrayView sourceView(sourceData, sourceShape);
    int val=0;
    for (int i=0;i<sourceShape[0];i++) {
      for (int j=0;j<sourceShape[1];j++) {
	sourceData[getRelIndex(sourceShape,i,j)]=val++;
      }
    }

    len = DataTypes::noValues(targetShape);
    RealVectorType targetData(len, 2.0, len);
//     DataArrayView targetView(targetData, targetShape);

    // Copy source view to target view.
//     targetView.copySliceFrom(sourceView,region);
    copySliceFrom(targetData,targetShape,0,sourceData, sourceShape,0,region);

    // Check results of copy.
    for (int i=region[0].first;i<region[0].second;i++) {
      for (int j=region[1].first;j<region[1].second;j++) {
// 	CPPUNIT_ASSERT(sourceView(i-region[0].first,j-region[1].first)==
//                targetView(i,j));
	CPPUNIT_ASSERT(sourceData[getRelIndex(sourceShape,i-region[0].first,j-region[1].first)]==
               targetData[getRelIndex(targetShape,i,j)]);
      }
    }
  }

  {
    cout << endl;
    cout << "\tSlice a 2 dimensional slice from a scalar.";

    // Define slice region.
    DataTypes::RegionType region;
    region.push_back(DataTypes::RegionType::value_type(2,4));
    region.push_back(DataTypes::RegionType::value_type(0,2));

    // Define shapes of views.
    DataTypes::ShapeType sourceShape;
    DataTypes::ShapeType targetShape;
    targetShape.push_back(6);
    targetShape.push_back(3);

    // Create source and target views.
    int len = DataTypes::noValues(sourceShape);
    RealVectorType sourceData(len, 2.0, len);
//     DataArrayView sourceView(sourceData, sourceShape);
    sourceData[0]=5;

    len = DataTypes::noValues(targetShape);
    RealVectorType targetData(len, 2.0, len);
//     DataArrayView targetView(targetData, targetShape);

    // Copy source view to target view.
//     targetView.copySliceFrom(sourceView,region);
    copySliceFrom(targetData,targetShape,0,sourceData,sourceShape,0,region);

    // Check results of copy.
    for (int i=region[0].first;i<region[0].second;i++) {
      for (int j=region[1].first;j<region[1].second;j++) {
	CPPUNIT_ASSERT(sourceData[0]==
               targetData[getRelIndex(targetShape,i,j)]);
      }
    }
  }

  {
    cout << endl;
    cout << "\tSlice a 3 dimensional slice to a 3 dimensional array.";

    // Define slice region.
    DataTypes::RegionType region;
    region.push_back(DataTypes::RegionType::value_type(2,4));
    region.push_back(DataTypes::RegionType::value_type(0,2));
    region.push_back(DataTypes::RegionType::value_type(5,9));

    // Define shapes of views.
    DataTypes::ShapeType sourceShape;
    sourceShape.push_back(6);
    sourceShape.push_back(3);
    sourceShape.push_back(13);
    DataTypes::ShapeType targetShape = DataTypes::getResultSliceShape(region);

    // Create source and target views.
    int len = DataTypes::noValues(sourceShape);
    RealVectorType sourceData(len, 2.0, len);
//     DataArrayView sourceView(sourceData, sourceShape);
    int val=0;
    for (int i=0;i<sourceShape[0];i++) {
      for (int j=0;j<sourceShape[1];j++) {
        for (int k=0;k<sourceShape[2];k++) {
	  sourceData[getRelIndex(sourceShape,i,j,k)]=val++;
        }
      }
    }

    len = DataTypes::noValues(targetShape);
    RealVectorType targetData(len, 2.0, len);
//     DataArrayView targetView(targetData, targetShape);

    // Copy source view to target view.
//     targetView.copySlice(sourceView,region);
    copySlice(targetData,targetShape,0,sourceData,sourceShape,0,region);

    // Check results of copy.
    for (int i=region[0].first;i<region[0].second;i++) {
      for (int j=region[1].first;j<region[1].second;j++) {
        for (int k=region[2].first;k<region[2].second;k++) {
// 	  CPPUNIT_ASSERT(sourceView(i,j,k)==
//                  targetView(i-region[0].first,j-region[1].first,k-region[2].first));
	  CPPUNIT_ASSERT(sourceData[getRelIndex(sourceShape,i,j,k)]==
                 targetData[getRelIndex(targetShape,i-region[0].first,j-region[1].first,k-region[2].first)]);
        }
      }
    }
  }

  {
    cout << endl;
    cout << "\tSlice a 3 dimensional slice to a 2 dimensional array.";

    // Define slice region.
    DataTypes::RegionType region;
    region.push_back(DataTypes::RegionType::value_type(2,4));
    region.push_back(DataTypes::RegionType::value_type(0,1));
    region.push_back(DataTypes::RegionType::value_type(5,9));

    // Define shapes of views.
    DataTypes::ShapeType sourceShape;
    sourceShape.push_back(6);
    sourceShape.push_back(3);
    sourceShape.push_back(13);
    DataTypes::ShapeType targetShape;
    targetShape.push_back(2);
    targetShape.push_back(4);

    // Create source and target views.
    int len = DataTypes::noValues(sourceShape);
    RealVectorType sourceData(len, 2.0, len);
//     DataArrayView sourceView(sourceData, sourceShape);
    int val=0;
    for (int i=0;i<sourceShape[0];i++) {
      for (int j=0;j<sourceShape[1];j++) {
        for (int k=0;k<sourceShape[2];k++) {
	  sourceData[getRelIndex(sourceShape,i,j,k)]=val++;
        }
      }
    }

    len = DataTypes::noValues(targetShape);
    RealVectorType targetData(len, 2.0, len);
//     DataArrayView targetView(targetData, targetShape);

    // Copy source view to target view.
//     targetView.copySlice(sourceView,region);
    copySlice(targetData, targetShape,0,sourceData,sourceShape,0,region);

    // Check results of copy.
    for (int i=region[0].first;i<region[0].second;i++) {
      for (int j=region[1].first;j<region[1].second;j++) {
        for (int k=region[2].first;k<region[2].second;k++) {
/*	  CPPUNIT_ASSERT(sourceView(i,j,k)==
                 targetView(i-region[0].first,k-region[2].first));*/	 
          CPPUNIT_ASSERT(sourceData[getRelIndex(sourceShape,i,j,k)]==
                 targetData[getRelIndex(targetShape,i-region[0].first,k-region[2].first)]);
        }
      }
    }
  }

  {
    cout << endl;
    cout << "\tSlice a 3 dimensional slice to a 1 dimensional array.";

    // Define slice region.
    DataTypes::RegionType region;
    region.push_back(DataTypes::RegionType::value_type(3,4));
    region.push_back(DataTypes::RegionType::value_type(0,1));
    region.push_back(DataTypes::RegionType::value_type(5,9));

    // Define shapes of views.
    DataTypes::ShapeType sourceShape;
    sourceShape.push_back(6);
    sourceShape.push_back(3);
    sourceShape.push_back(13);
    DataTypes::ShapeType targetShape;
    targetShape.push_back(4);

    // Create source and target views.
    int len = DataTypes::noValues(sourceShape);
    RealVectorType sourceData(len, 2.0, len);
//     DataArrayView sourceView(sourceData, sourceShape);
    int val=0;
    for (int i=0;i<sourceShape[0];i++) {
      for (int j=0;j<sourceShape[1];j++) {
        for (int k=0;k<sourceShape[2];k++) {
	  sourceData[getRelIndex(sourceShape,i,j,k)]=val++;
        }
      }
    }

    len = DataTypes::noValues(targetShape);
    RealVectorType targetData(len, 2.0, len);
//     DataArrayView targetView(targetData, targetShape);

    // Copy source view to target view.
//     targetView.copySlice(sourceView,region);
    copySlice(targetData, targetShape,0,sourceData, sourceShape, 0, region);

    // Check results of copy.
    for (int i=region[0].first;i<region[0].second;i++) {
      for (int j=region[1].first;j<region[1].second;j++) {
        for (int k=region[2].first;k<region[2].second;k++) {
	  CPPUNIT_ASSERT(sourceData[getRelIndex(sourceShape,i,j,k)]==
                 targetData[getRelIndex(targetShape,k-region[2].first)]);
        }
      }
    }
  }

  {
    cout << endl;
    cout << "\tSlice a 3 dimensional slice to a scalar.";

    // Define slice region.
    DataTypes::RegionType region;
    region.push_back(DataTypes::RegionType::value_type(3,4));
    region.push_back(DataTypes::RegionType::value_type(0,1));
    region.push_back(DataTypes::RegionType::value_type(5,6));

    // Define shapes of views.
    DataTypes::ShapeType sourceShape;
    sourceShape.push_back(6);
    sourceShape.push_back(3);
    sourceShape.push_back(13);
    DataTypes::ShapeType targetShape;

    // Create source and target views.
    int len = DataTypes::noValues(sourceShape);
    RealVectorType sourceData(len, 2.0, len);
//     DataArrayView sourceView(sourceData, sourceShape);
    int val=0;
    for (int i=0;i<sourceShape[0];i++) {
      for (int j=0;j<sourceShape[1];j++) {
        for (int k=0;k<sourceShape[2];k++) {
	  sourceData[getRelIndex(sourceShape,i,j,k)]=val++;
        }
      }
    }

    len = DataTypes::noValues(targetShape);
    RealVectorType targetData(len, 2.0, len);
//     DataArrayView targetView(targetData, targetShape);

    // Copy source view to target view.
//     targetView.copySlice(sourceView,region);
    copySlice(targetData, targetShape,0,sourceData, sourceShape, 0, region);

    // Check results of copy.
    for (int i=region[0].first;i<region[0].second;i++) {
      for (int j=region[1].first;j<region[1].second;j++) {
        for (int k=region[2].first;k<region[2].second;k++) {
	  CPPUNIT_ASSERT(sourceData[getRelIndex(sourceShape, i,j,k)]==
                 targetData[0]);
        }
      }
    }
  }

  {
    cout << endl;
    cout << "\tSlice a 3 dimensional slice from a 3 dimensional array.";

    // Define slice region.
    DataTypes::RegionType region;
    region.push_back(DataTypes::RegionType::value_type(4,7));
    region.push_back(DataTypes::RegionType::value_type(2,5));
    region.push_back(DataTypes::RegionType::value_type(6,7));

    // Define shapes of views.
    DataTypes::ShapeType sourceShape = DataTypes::getResultSliceShape(region);
    DataTypes::ShapeType targetShape;
    targetShape.push_back(11);
    targetShape.push_back(8);
    targetShape.push_back(9);

    // Create source and target views.
    int len = DataTypes::noValues(sourceShape);
    RealVectorType sourceData(len, 2.0, len);
//     DataArrayView sourceView(sourceData, sourceShape);
    int val=0;
    for (int i=0;i<sourceShape[0];i++) {
      for (int j=0;j<sourceShape[1];j++) {
        for (int k=0;k<sourceShape[2];k++) {
	  sourceData[getRelIndex(sourceShape,i,j,k)]=val++;
        }
      }
    }

    len = DataTypes::noValues(targetShape);
    RealVectorType targetData(len, 2.0, len);
//     DataArrayView targetView(targetData, targetShape);

    // Copy source view to target view.
//     targetView.copySliceFrom(sourceView,region);
    copySliceFrom(targetData, targetShape,0,sourceData,sourceShape,0, region);

    // Check results of copy.
    for (int i=region[0].first;i<region[0].second;i++) {
      for (int j=region[1].first;j<region[1].second;j++) {
        for (int k=region[2].first;k<region[2].second;k++) {
	  CPPUNIT_ASSERT(sourceData[getRelIndex(sourceShape,i-region[0].first,j-region[1].first,k-region[2].first)]==
                 targetData[getRelIndex(targetShape,i,j,k)]);
        }
      }
    }
  }

  {
    cout << endl;
    cout << "\tSlice a 3 dimensional slice from a scalar.";

    // Define slice region.
    DataTypes::RegionType region;
    region.push_back(DataTypes::RegionType::value_type(4,7));
    region.push_back(DataTypes::RegionType::value_type(2,5));
    region.push_back(DataTypes::RegionType::value_type(6,7));

    // Define shapes of views.
    DataTypes::ShapeType sourceShape;
    DataTypes::ShapeType targetShape;
    targetShape.push_back(11);
    targetShape.push_back(8);
    targetShape.push_back(9);

    // Create source and target views.
    int len = DataTypes::noValues(sourceShape);
    RealVectorType sourceData(len, 2.0, len);
//     DataArrayView sourceView(sourceData, sourceShape);
    sourceData[0]=5;

    len = DataTypes::noValues(targetShape);
    RealVectorType targetData(len, 2.0, len);

    // Copy source view to target view.
//     targetView.copySliceFrom(sourceView,region);
    copySliceFrom(targetData, targetShape,0,sourceData, sourceShape, 0, region);

    // Check results of copy.
    for (int i=region[0].first;i<region[0].second;i++) {
      for (int j=region[1].first;j<region[1].second;j++) {
        for (int k=region[2].first;k<region[2].second;k++) {
	  CPPUNIT_ASSERT(sourceData[0]==
                 targetData[getRelIndex(targetShape,i,j,k)]);
        }
      }
    }
  }

  {
    cout << endl;
    cout << "\tSlice a 4 dimensional slice to a 4 dimensional array.";

    // Define slice region.
    DataTypes::RegionType region;
    region.push_back(DataTypes::RegionType::value_type(2,4));
    region.push_back(DataTypes::RegionType::value_type(0,2));
    region.push_back(DataTypes::RegionType::value_type(5,9));
    region.push_back(DataTypes::RegionType::value_type(3,5));

    // Define shapes of views.
    DataTypes::ShapeType sourceShape;
    sourceShape.push_back(6);
    sourceShape.push_back(3);
    sourceShape.push_back(13);
    sourceShape.push_back(9);
    DataTypes::ShapeType targetShape = DataTypes::getResultSliceShape(region);

    // Create source and target views.
    int len = DataTypes::noValues(sourceShape);
    RealVectorType sourceData(len, 2.0, len);
//     DataArrayView sourceView(sourceData, sourceShape);
    int val=0;
    for (int i=0;i<sourceShape[0];i++) {
      for (int j=0;j<sourceShape[1];j++) {
        for (int k=0;k<sourceShape[2];k++) {
          for (int l=0;l<sourceShape[3];l++) {
	    sourceData[getRelIndex(sourceShape,i,j,k,l)]=val++;
          }
        }
      }
    }

    len = DataTypes::noValues(targetShape);
    RealVectorType targetData(len, 2.0, len);
//     DataArrayView targetView(targetData, targetShape);

    // Copy source view to target view.
//     targetView.copySlice(sourceView,region);
    copySlice(targetData, targetShape, 0, sourceData, sourceShape,0, region);

    // Check results of copy.
    for (int i=region[0].first;i<region[0].second;i++) {
      for (int j=region[1].first;j<region[1].second;j++) {
        for (int k=region[2].first;k<region[2].second;k++) {
          for (int l=region[3].first;l<region[3].second;l++) {
// 	    CPPUNIT_ASSERT(sourceView(i,j,k,l)==
//                    targetView(i-region[0].first,j-region[1].first,k-region[2].first,l-region[3].first));
	    CPPUNIT_ASSERT(sourceData[getRelIndex(sourceShape,i,j,k,l)]==
                   targetData[getRelIndex(targetShape,i-region[0].first,j-region[1].first,k-region[2].first,l-region[3].first)]);
          }
        }
      }
    }
  }

  {
    cout << endl;
    cout << "\tSlice a 4 dimensional slice to a 3 dimensional array.";

    // Define slice region.
    DataTypes::RegionType region;
    region.push_back(DataTypes::RegionType::value_type(2,4));
    region.push_back(DataTypes::RegionType::value_type(0,2));
    region.push_back(DataTypes::RegionType::value_type(5,6));
    region.push_back(DataTypes::RegionType::value_type(3,5));

    // Define shapes of views.
    DataTypes::ShapeType sourceShape;
    sourceShape.push_back(6);
    sourceShape.push_back(3);
    sourceShape.push_back(13);
    sourceShape.push_back(9);
    DataTypes::ShapeType targetShape;
    targetShape.push_back(2);
    targetShape.push_back(2);
    targetShape.push_back(2);

    // Create source and target views.
    int len = DataTypes::noValues(sourceShape);
    RealVectorType sourceData(len, 2.0, len);
//     DataArrayView sourceView(sourceData, sourceShape);
    int val=0;
    for (int i=0;i<sourceShape[0];i++) {
      for (int j=0;j<sourceShape[1];j++) {
        for (int k=0;k<sourceShape[2];k++) {
          for (int l=0;l<sourceShape[3];l++) {
	    sourceData[getRelIndex(sourceShape,i,j,k,l)]=val++;
          }
        }
      }
    }

    len = DataTypes::noValues(targetShape);
    RealVectorType targetData(len, 2.0, len);
//     DataArrayView targetView(targetData, targetShape);

    // Copy source view to target view.
//     targetView.copySlice(sourceView,region);
    copySlice(targetData, targetShape, 0, sourceData, sourceShape, 0, region);

    // Check results of copy.
    for (int i=region[0].first;i<region[0].second;i++) {
      for (int j=region[1].first;j<region[1].second;j++) {
        for (int k=region[2].first;k<region[2].second;k++) {
          for (int l=region[3].first;l<region[3].second;l++) {
	    CPPUNIT_ASSERT(sourceData[getRelIndex(sourceShape,i,j,k,l)]==
                   targetData[getRelIndex(targetShape,i-region[0].first,j-region[1].first,l-region[3].first)]);
          }
        }
      }
    }
  }

  {
    cout << endl;
    cout << "\tSlice a 4 dimensional slice to a 2 dimensional array.";

    // Define slice region.
    DataTypes::RegionType region;
    region.push_back(DataTypes::RegionType::value_type(2,4));
    region.push_back(DataTypes::RegionType::value_type(0,2));
    region.push_back(DataTypes::RegionType::value_type(5,6));
    region.push_back(DataTypes::RegionType::value_type(4,5));

    // Define shapes of views.
    DataTypes::ShapeType sourceShape;
    sourceShape.push_back(6);
    sourceShape.push_back(3);
    sourceShape.push_back(13);
    sourceShape.push_back(9);
    DataTypes::ShapeType targetShape;
    targetShape.push_back(2);
    targetShape.push_back(2);

    // Create source and target views.
    int len = DataTypes::noValues(sourceShape);
    RealVectorType sourceData(len, 2.0, len);
//     DataArrayView sourceView(sourceData, sourceShape);
    int val=0;
    for (int i=0;i<sourceShape[0];i++) {
      for (int j=0;j<sourceShape[1];j++) {
        for (int k=0;k<sourceShape[2];k++) {
          for (int l=0;l<sourceShape[3];l++) {
	    sourceData[getRelIndex(sourceShape,i,j,k,l)]=val++;
          }
        }
      }
    }

    len = DataTypes::noValues(targetShape);
    RealVectorType targetData(len, 2.0, len);
//     DataArrayView targetView(targetData, targetShape);

    // Copy source view to target view.
//     targetView.copySlice(sourceView,region);
    copySlice(targetData, targetShape, 0, sourceData, sourceShape, 0, region);

    // Check results of copy.
    for (int i=region[0].first;i<region[0].second;i++) {
      for (int j=region[1].first;j<region[1].second;j++) {
        for (int k=region[2].first;k<region[2].second;k++) {
          for (int l=region[3].first;l<region[3].second;l++) {
	    CPPUNIT_ASSERT(sourceData[getRelIndex(sourceShape,i,j,k,l)]==
                   targetData[getRelIndex(targetShape,i-region[0].first,j-region[1].first)]);
          }
        }
      }
    }
  }

  {
    cout << endl;
    cout << "\tSlice a 4 dimensional slice to a 1 dimensional array.";

    // Define slice region.
    DataTypes::RegionType region;
    region.push_back(DataTypes::RegionType::value_type(3,4));
    region.push_back(DataTypes::RegionType::value_type(0,2));
    region.push_back(DataTypes::RegionType::value_type(5,6));
    region.push_back(DataTypes::RegionType::value_type(4,5));

    // Define shapes of views.
    DataTypes::ShapeType sourceShape;
    sourceShape.push_back(6);
    sourceShape.push_back(3);
    sourceShape.push_back(13);
    sourceShape.push_back(9);
    DataTypes::ShapeType targetShape;
    targetShape.push_back(2);

    // Create source and target views.
    int len = DataTypes::noValues(sourceShape);
    RealVectorType sourceData(len, 2.0, len);
//     DataArrayView sourceView(sourceData, sourceShape);
    int val=0;
    for (int i=0;i<sourceShape[0];i++) {
      for (int j=0;j<sourceShape[1];j++) {
        for (int k=0;k<sourceShape[2];k++) {
          for (int l=0;l<sourceShape[3];l++) {
	    sourceData[getRelIndex(sourceShape,i,j,k,l)]=val++;
          }
        }
      }
    }

    len = DataTypes::noValues(targetShape);
    RealVectorType targetData(len, 2.0, len);
//     DataArrayView targetView(targetData, targetShape);

    // Copy source view to target view.
//     targetView.copySlice(sourceView,region);
    copySlice(targetData,targetShape,0, sourceData, sourceShape,0,region);

    // Check results of copy.
    for (int i=region[0].first;i<region[0].second;i++) {
      for (int j=region[1].first;j<region[1].second;j++) {
        for (int k=region[2].first;k<region[2].second;k++) {
          for (int l=region[3].first;l<region[3].second;l++) {
	    CPPUNIT_ASSERT(sourceData[getRelIndex(sourceShape,i,j,k,l)]==
                   targetData[getRelIndex(targetShape,j-region[1].first)]);
          }
        }
      }
    }
  }

  {
    cout << endl;
    cout << "\tSlice a 4 dimensional slice to a scalar.";

    // Define slice region.
    DataTypes::RegionType region;
    region.push_back(DataTypes::RegionType::value_type(3,4));
    region.push_back(DataTypes::RegionType::value_type(1,2));
    region.push_back(DataTypes::RegionType::value_type(5,6));
    region.push_back(DataTypes::RegionType::value_type(4,5));

    // Define shapes of views.
    DataTypes::ShapeType sourceShape;
    sourceShape.push_back(6);
    sourceShape.push_back(3);
    sourceShape.push_back(13);
    sourceShape.push_back(9);
    DataTypes::ShapeType targetShape;

    // Create source and target views.
    int len = DataTypes::noValues(sourceShape);
    RealVectorType sourceData(len, 2.0, len);
    //DataArrayView sourceView(sourceData, sourceShape);
    int val=0;
    for (int i=0;i<sourceShape[0];i++) {
      for (int j=0;j<sourceShape[1];j++) {
        for (int k=0;k<sourceShape[2];k++) {
          for (int l=0;l<sourceShape[3];l++) {
	    sourceData[getRelIndex(sourceShape,i,j,k,l)]=val++;
          }
        }
      }
    }

    len = DataTypes::noValues(targetShape);
    RealVectorType targetData(len, 2.0, len);
//     DataArrayView targetView(targetData, targetShape);

    // Copy source view to target view.
//     targetView.copySlice(sourceView,region);
    copySlice(targetData, targetShape, 0, sourceData, sourceShape,0,region);

    // Check results of copy.
    for (int i=region[0].first;i<region[0].second;i++) {
      for (int j=region[1].first;j<region[1].second;j++) {
        for (int k=region[2].first;k<region[2].second;k++) {
          for (int l=region[3].first;l<region[3].second;l++) {
	    CPPUNIT_ASSERT(sourceData[getRelIndex(sourceShape,i,j,k,l)]==
                   targetData[0]);
          }
        }
      }
    }
  }

  {
    cout << endl;
    cout << "\tSlice a 4 dimensional slice from a 4 dimensional array.";

    // Define slice region.
    DataTypes::RegionType region;
    region.push_back(DataTypes::RegionType::value_type(14,37));
    region.push_back(DataTypes::RegionType::value_type(22,57));
    region.push_back(DataTypes::RegionType::value_type(63,71));
    region.push_back(DataTypes::RegionType::value_type(23,51));

    // Define shapes of views.
    DataTypes::ShapeType sourceShape = DataTypes::getResultSliceShape(region);
    DataTypes::ShapeType targetShape;
    targetShape.push_back(50);
    targetShape.push_back(65);
    targetShape.push_back(80);
    targetShape.push_back(90);

    // Create source and target views.
    int len = DataTypes::noValues(sourceShape);
    RealVectorType sourceData(len, 2.0, len);
//     DataArrayView sourceView(sourceData, sourceShape);
    int val=0;
    for (int i=0;i<sourceShape[0];i++) {
      for (int j=0;j<sourceShape[1];j++) {
        for (int k=0;k<sourceShape[2];k++) {
          for (int l=0;l<sourceShape[3];l++) {
	    sourceData[getRelIndex(sourceShape,i,j,k,l)]=val++;
          }
        }
      }
    }

    len = DataTypes::noValues(targetShape);
    RealVectorType targetData(len, 2.0, len);
//     DataArrayView targetView(targetData, targetShape);

    // Copy source view to target view.
//     targetView.copySliceFrom(sourceView,region);
    copySliceFrom(targetData, targetShape, 0, sourceData, sourceShape,0,region);

    // Check results of copy.
    for (int i=region[0].first;i<region[0].second;i++) {
      for (int j=region[1].first;j<region[1].second;j++) {
        for (int k=region[2].first;k<region[2].second;k++) {
          for (int l=region[3].first;l<region[3].second;l++) {
	    CPPUNIT_ASSERT(sourceData[getRelIndex(sourceShape,i-region[0].first,j-region[1].first,k-region[2].first,l-region[3].first)]==
                   targetData[getRelIndex(targetShape,i,j,k,l)]);
          }
        }
      }
    }
  }

  {
    cout << endl;
    cout << "\tSlice a 4 dimensional slice from a scalar.";

    // Define slice region.
    DataTypes::RegionType region;
    region.push_back(DataTypes::RegionType::value_type(14,37));
    region.push_back(DataTypes::RegionType::value_type(22,57));
    region.push_back(DataTypes::RegionType::value_type(63,71));
    region.push_back(DataTypes::RegionType::value_type(23,51));

    // Define shapes of views.
    DataTypes::ShapeType sourceShape;
    DataTypes::ShapeType targetShape;
    targetShape.push_back(50);
    targetShape.push_back(65);
    targetShape.push_back(80);
    targetShape.push_back(90);

    // Create source and target views.
    int len = DataTypes::noValues(sourceShape);
    RealVectorType sourceData(len, 2.0, len);
//     DataArrayView sourceView(sourceData, sourceShape);
    sourceData[0]=5;

    len = DataTypes::noValues(targetShape);
    RealVectorType targetData(len, 2.0, len);
//     DataArrayView targetView(targetData, targetShape);

    // Copy source view to target view.
//     targetView.copySliceFrom(sourceView,region);
    copySliceFrom(targetData, targetShape, 0, sourceData, sourceShape, 0, region);

    // Check results of copy.
    for (int i=region[0].first;i<region[0].second;i++) {
      for (int j=region[1].first;j<region[1].second;j++) {
        for (int k=region[2].first;k<region[2].second;k++) {
          for (int l=region[3].first;l<region[3].second;l++) {
	    CPPUNIT_ASSERT(sourceData[0]==
                   targetData[getRelIndex(targetShape,i,j,k,l)]);
          }
        }
      }
    }
  }

  cout << endl;

}

void DataTypesTestCase::testShapeToString() {

  cout << endl;
  cout << "\tTest shapeToString for a variety of shapes." << endl;

  DataTypes::ShapeType shape;
  CPPUNIT_ASSERT(DataTypes::shapeToString(shape)=="()");
  shape.push_back(5);
  CPPUNIT_ASSERT(DataTypes::shapeToString(shape)=="(5)");
  shape.push_back(2);
  CPPUNIT_ASSERT(DataTypes::shapeToString(shape)=="(5,2)");
  shape.push_back(9);
  CPPUNIT_ASSERT(DataTypes::shapeToString(shape)=="(5,2,9)");
  shape.push_back(4);
  CPPUNIT_ASSERT(DataTypes::shapeToString(shape)=="(5,2,9,4)");

}


TestSuite* DataTypesTestCase::suite()
{
  //
  // create the suite of tests to perform.
  TestSuite *testSuite = new TestSuite("DataTypesTestCase");
  testSuite->addTest(new TestCaller<DataTypesTestCase>(
              "testShapeToString",&DataTypesTestCase::testShapeToString));
  testSuite->addTest(new TestCaller<DataTypesTestCase>(
              "testResultSliceShape",&DataTypesTestCase::testResultSliceShape));
  testSuite->addTest(new TestCaller<DataTypesTestCase>(
              "testSlicing",&DataTypesTestCase::testSlicing));
  testSuite->addTest(new TestCaller<DataTypesTestCase>(
              "testShapeFunctions",&DataTypesTestCase::testShapeFns));
  return testSuite;
}

