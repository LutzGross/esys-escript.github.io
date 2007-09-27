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
#include "escript/DataArray.h"
#include "escript/DataArrayView.h"
#include "escript/DataAlgorithm.h"
#include "esysUtils/EsysException.h"

#include "DataArrayViewTestCase.h"

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

void DataArrayViewTestCase::testResultSliceShape() {

  cout << endl;
  cout << "\tTest getResultSliceShape method." << endl;

  DataArrayView::RegionType region;
  DataArrayView::ShapeType resultShape;

  region.push_back(DataArrayView::RegionType::value_type(1,5));
  resultShape.push_back(4);
  assert(DataArrayView::getResultSliceShape(region)==resultShape);

  region.push_back(DataArrayView::RegionType::value_type(2,5));
  resultShape.push_back(3);
  assert(DataArrayView::getResultSliceShape(region)==resultShape);

  region.push_back(DataArrayView::RegionType::value_type(3,9));
  resultShape.push_back(6);
  assert(DataArrayView::getResultSliceShape(region)==resultShape);

  region.push_back(DataArrayView::RegionType::value_type(1,7));
  resultShape.push_back(6);
  assert(DataArrayView::getResultSliceShape(region)==resultShape);

}

void DataArrayViewTestCase::testSlicing() {

  {
    cout << endl;
    cout << "\tSlice a scalar to a scalar.";

    // Define slice region.
    DataArrayView::RegionType region;

    // Define shape of views.
    DataArrayView::ShapeType sourceShape;

    // Create source and target views.
    DataArray source(sourceShape,2.0);
    DataArrayView sourceView=source.getView();

    DataArray target;
    DataArrayView targetView=target.getView();

    // Copy source view to target view.
    targetView.copySlice(sourceView,region);

    // Check results of copy.
    assert(sourceView==targetView);
  }

  {
    cout << endl;
    cout << "\tSlice a scalar from a scalar.";

    // Define slice region.
    DataArrayView::RegionType region;

    // Define shape of views.
    DataArrayView::ShapeType sourceShape;

    // Create source and target views.
    DataArray source(sourceShape,2.0);
    DataArrayView sourceView=source.getView();

    DataArray target;
    DataArrayView targetView=target.getView();

    // Copy source view to target view.
    targetView.copySliceFrom(sourceView,region);

    // Check results of copy.
    assert(sourceView==targetView);
  }

  {
    cout << endl;
    cout << "\tSlice a 1 dimensional slice to a 1 dimensional array.";

    // Define slice region.
    DataArrayView::RegionType region;
    region.push_back(DataArrayView::RegionType::value_type(2,4));

    // Define shapes of views.
    DataArrayView::ShapeType sourceShape;
    sourceShape.push_back(6);
    DataArrayView::ShapeType targetShape = DataArrayView::getResultSliceShape(region);

    // Create source and target views.
    DataArray source(sourceShape);
    DataArrayView sourceView=source.getView();
    for (int i=0;i<sourceShape[0];i++) {
      sourceView(i)=i;
    }

    DataArray target(targetShape);
    DataArrayView targetView=target.getView();

    // Copy source view to target view.
    targetView.copySlice(sourceView,region);

    // Check results of copy.
    for (int i=region[0].first;i<region[0].second;i++) {
      assert(sourceView(i)==
             targetView(i-region[0].first));
    }
  }

  {
    cout << endl;
    cout << "\tSlice a 1 dimensional slice to a scalar.";

    // Define slice region.
    DataArrayView::RegionType region;
    region.push_back(DataArrayView::RegionType::value_type(2,3));

    // Define shapes of views.
    DataArrayView::ShapeType sourceShape;
    sourceShape.push_back(6);
    DataArrayView::ShapeType targetShape;

    // Create source and target views.
    DataArray source(sourceShape);
    DataArrayView sourceView=source.getView();
    for (int i=0;i<sourceShape[0];i++) {
      sourceView(i)=i;
    }

    DataArray target(targetShape);
    DataArrayView targetView=target.getView();

    // Copy source view to target view.
    targetView.copySlice(sourceView,region);

    // Check results of copy.
    for (int i=region[0].first;i<region[0].second;i++) {
      assert(sourceView(i)==
             targetView());
    }
  }

  {
    cout << endl;
    cout << "\tSlice a 1 dimensional slice from a 1 dimensional array.";

    // Define slice region.
    DataArrayView::RegionType region;
    region.push_back(DataArrayView::RegionType::value_type(2,4));

    // Define shapes of views.
    DataArrayView::ShapeType sourceShape = DataArrayView::getResultSliceShape(region);
    DataArrayView::ShapeType targetShape;
    targetShape.push_back(6);

    // Create source and target views.
    DataArray source(sourceShape);
    DataArrayView sourceView=source.getView();
    for (int i=0;i<sourceShape[0];i++) {
      sourceView(i)=i;
    }

    DataArray target(targetShape);
    DataArrayView targetView=target.getView();

    // Copy source view to target view.
    targetView.copySliceFrom(sourceView,region);

    // Check results of copy.
    for (int i=region[0].first;i<region[0].second;i++) {
      assert(sourceView(i-region[0].first)==
             targetView(i));
    }
  }

  {
    cout << endl;
    cout << "\tSlice a 1 dimensional slice from a scalar.";

    // Define slice region.
    DataArrayView::RegionType region;
    region.push_back(DataArrayView::RegionType::value_type(2,4));

    // Define shapes of views.
    DataArrayView::ShapeType sourceShape;
    DataArrayView::ShapeType targetShape;
    targetShape.push_back(6);

    // Create source and target views.
    DataArray source(sourceShape);
    DataArrayView sourceView=source.getView();
    sourceView()=5;

    DataArray target(targetShape);
    DataArrayView targetView=target.getView();

    // Copy source view to target view.
    targetView.copySliceFrom(sourceView,region);

    // Check results of copy.
    for (int i=region[0].first;i<region[0].second;i++) {
      assert(sourceView()==
             targetView(i));
    }
  }

  {
    cout << endl;
    cout << "\tSlice a 2 dimensional slice to a 2 dimensional array.";

    // Define slice region.
    DataArrayView::RegionType region;
    region.push_back(DataArrayView::RegionType::value_type(2,4));
    region.push_back(DataArrayView::RegionType::value_type(0,2));

    // Define shapes of views.
    DataArrayView::ShapeType sourceShape;
    sourceShape.push_back(6);
    sourceShape.push_back(3);
    DataArrayView::ShapeType targetShape = DataArrayView::getResultSliceShape(region);

    // Create source and target views.
    DataArray source(sourceShape);
    DataArrayView sourceView=source.getView();
    int val=0;
    for (int i=0;i<sourceShape[0];i++) {
      for (int j=0;j<sourceShape[1];j++) {
	sourceView(i,j)=val++;
      }
    }

    DataArray target(targetShape);
    DataArrayView targetView=target.getView();

    // Copy source view to target view.
    targetView.copySlice(sourceView,region);

    // Check results of copy.
    for (int i=region[0].first;i<region[0].second;i++) {
      for (int j=region[1].first;j<region[1].second;j++) {
	assert(sourceView(i,j)==
               targetView(i-region[0].first,j-region[1].first));
      }
    }
  }

  {
    cout << endl;
    cout << "\tSlice a 2 dimensional slice to a 1 dimensional array.";

    // Define slice region.
    DataArrayView::RegionType region;
    region.push_back(DataArrayView::RegionType::value_type(0,3));
    region.push_back(DataArrayView::RegionType::value_type(1,2));

    // Define shapes of views.
    DataArrayView::ShapeType sourceShape;
    sourceShape.push_back(3);
    sourceShape.push_back(6);
    DataArrayView::ShapeType targetShape;
    targetShape.push_back(3);

    // Create source and target views.
    DataArray source(sourceShape);
    DataArrayView sourceView=source.getView();
    int val=0;
    for (int i=0;i<sourceShape[0];i++) {
      for (int j=0;j<sourceShape[1];j++) {
	sourceView(i,j)=val++;
      }
    }

    DataArray target(targetShape);
    DataArrayView targetView=target.getView();

    // Copy source view to target view.
    targetView.copySlice(sourceView,region);

    // Check results of copy.
    for (int i=region[0].first;i<region[0].second;i++) {
      for (int j=region[1].first;j<region[1].second;j++) {
	assert(sourceView(i,j)==
               targetView(i-region[0].first));
      }
    }
  }

  {
    cout << endl;
    cout << "\tSlice a 2 dimensional slice to a scalar.";

    // Define slice region.
    DataArrayView::RegionType region;
    region.push_back(DataArrayView::RegionType::value_type(2,3));
    region.push_back(DataArrayView::RegionType::value_type(1,2));

    // Define shapes of views.
    DataArrayView::ShapeType sourceShape;
    sourceShape.push_back(3);
    sourceShape.push_back(6);
    DataArrayView::ShapeType targetShape;

    // Create source and target views.
    DataArray source(sourceShape);
    DataArrayView sourceView=source.getView();
    int val=0;
    for (int i=0;i<sourceShape[0];i++) {
      for (int j=0;j<sourceShape[1];j++) {
	sourceView(i,j)=val++;
      }
    }

    DataArray target(targetShape);
    DataArrayView targetView=target.getView();

    // Copy source view to target view.
    targetView.copySlice(sourceView,region);
    // Check results of copy.
    for (int i=region[0].first;i<region[0].second;i++) {
      for (int j=region[1].first;j<region[1].second;j++) {
	assert(sourceView(i,j)==
               targetView());
      }
    }
  }

  {
    cout << endl;
    cout << "\tSlice a 2 dimensional slice from a 2 dimensional array.";

    // Define slice region.
    DataArrayView::RegionType region;
    region.push_back(DataArrayView::RegionType::value_type(2,4));
    region.push_back(DataArrayView::RegionType::value_type(0,2));

    // Define shapes of views.
    DataArrayView::ShapeType sourceShape = DataArrayView::getResultSliceShape(region);
    DataArrayView::ShapeType targetShape;
    targetShape.push_back(6);
    targetShape.push_back(3);

    // Create source and target views.
    DataArray source(sourceShape);
    DataArrayView sourceView=source.getView();
    int val=0;
    for (int i=0;i<sourceShape[0];i++) {
      for (int j=0;j<sourceShape[1];j++) {
	sourceView(i,j)=val++;
      }
    }

    DataArray target(targetShape);
    DataArrayView targetView=target.getView();

    // Copy source view to target view.
    targetView.copySliceFrom(sourceView,region);

    // Check results of copy.
    for (int i=region[0].first;i<region[0].second;i++) {
      for (int j=region[1].first;j<region[1].second;j++) {
	assert(sourceView(i-region[0].first,j-region[1].first)==
               targetView(i,j));
      }
    }
  }

  {
    cout << endl;
    cout << "\tSlice a 2 dimensional slice from a scalar.";

    // Define slice region.
    DataArrayView::RegionType region;
    region.push_back(DataArrayView::RegionType::value_type(2,4));
    region.push_back(DataArrayView::RegionType::value_type(0,2));

    // Define shapes of views.
    DataArrayView::ShapeType sourceShape;
    DataArrayView::ShapeType targetShape;
    targetShape.push_back(6);
    targetShape.push_back(3);

    // Create source and target views.
    DataArray source(sourceShape);
    DataArrayView sourceView=source.getView();
    sourceView()=5;

    DataArray target(targetShape);
    DataArrayView targetView=target.getView();

    // Copy source view to target view.
    targetView.copySliceFrom(sourceView,region);

    // Check results of copy.
    for (int i=region[0].first;i<region[0].second;i++) {
      for (int j=region[1].first;j<region[1].second;j++) {
	assert(sourceView()==
               targetView(i,j));
      }
    }
  }

  {
    cout << endl;
    cout << "\tSlice a 3 dimensional slice to a 3 dimensional array.";

    // Define slice region.
    DataArrayView::RegionType region;
    region.push_back(DataArrayView::RegionType::value_type(2,4));
    region.push_back(DataArrayView::RegionType::value_type(0,2));
    region.push_back(DataArrayView::RegionType::value_type(5,9));

    // Define shapes of views.
    DataArrayView::ShapeType sourceShape;
    sourceShape.push_back(6);
    sourceShape.push_back(3);
    sourceShape.push_back(13);
    DataArrayView::ShapeType targetShape = DataArrayView::getResultSliceShape(region);

    // Create source and target views.
    DataArray source(sourceShape);
    DataArrayView sourceView=source.getView();
    int val=0;
    for (int i=0;i<sourceShape[0];i++) {
      for (int j=0;j<sourceShape[1];j++) {
        for (int k=0;k<sourceShape[2];k++) {
	  sourceView(i,j,k)=val++;
        }
      }
    }

    DataArray target(targetShape);
    DataArrayView targetView=target.getView();

    // Copy source view to target view.
    targetView.copySlice(sourceView,region);

    // Check results of copy.
    for (int i=region[0].first;i<region[0].second;i++) {
      for (int j=region[1].first;j<region[1].second;j++) {
        for (int k=region[2].first;k<region[2].second;k++) {
	  assert(sourceView(i,j,k)==
                 targetView(i-region[0].first,j-region[1].first,k-region[2].first));
        }
      }
    }
  }

  {
    cout << endl;
    cout << "\tSlice a 3 dimensional slice to a 2 dimensional array.";

    // Define slice region.
    DataArrayView::RegionType region;
    region.push_back(DataArrayView::RegionType::value_type(2,4));
    region.push_back(DataArrayView::RegionType::value_type(0,1));
    region.push_back(DataArrayView::RegionType::value_type(5,9));

    // Define shapes of views.
    DataArrayView::ShapeType sourceShape;
    sourceShape.push_back(6);
    sourceShape.push_back(3);
    sourceShape.push_back(13);
    DataArrayView::ShapeType targetShape;
    targetShape.push_back(2);
    targetShape.push_back(4);

    // Create source and target views.
    DataArray source(sourceShape);
    DataArrayView sourceView=source.getView();
    int val=0;
    for (int i=0;i<sourceShape[0];i++) {
      for (int j=0;j<sourceShape[1];j++) {
        for (int k=0;k<sourceShape[2];k++) {
	  sourceView(i,j,k)=val++;
        }
      }
    }

    DataArray target(targetShape);
    DataArrayView targetView=target.getView();

    // Copy source view to target view.
    targetView.copySlice(sourceView,region);

    // Check results of copy.
    for (int i=region[0].first;i<region[0].second;i++) {
      for (int j=region[1].first;j<region[1].second;j++) {
        for (int k=region[2].first;k<region[2].second;k++) {
	  assert(sourceView(i,j,k)==
                 targetView(i-region[0].first,k-region[2].first));
        }
      }
    }
  }

  {
    cout << endl;
    cout << "\tSlice a 3 dimensional slice to a 1 dimensional array.";

    // Define slice region.
    DataArrayView::RegionType region;
    region.push_back(DataArrayView::RegionType::value_type(3,4));
    region.push_back(DataArrayView::RegionType::value_type(0,1));
    region.push_back(DataArrayView::RegionType::value_type(5,9));

    // Define shapes of views.
    DataArrayView::ShapeType sourceShape;
    sourceShape.push_back(6);
    sourceShape.push_back(3);
    sourceShape.push_back(13);
    DataArrayView::ShapeType targetShape;
    targetShape.push_back(4);

    // Create source and target views.
    DataArray source(sourceShape);
    DataArrayView sourceView=source.getView();
    int val=0;
    for (int i=0;i<sourceShape[0];i++) {
      for (int j=0;j<sourceShape[1];j++) {
        for (int k=0;k<sourceShape[2];k++) {
	  sourceView(i,j,k)=val++;
        }
      }
    }

    DataArray target(targetShape);
    DataArrayView targetView=target.getView();

    // Copy source view to target view.
    targetView.copySlice(sourceView,region);

    // Check results of copy.
    for (int i=region[0].first;i<region[0].second;i++) {
      for (int j=region[1].first;j<region[1].second;j++) {
        for (int k=region[2].first;k<region[2].second;k++) {
	  assert(sourceView(i,j,k)==
                 targetView(k-region[2].first));
        }
      }
    }
  }

  {
    cout << endl;
    cout << "\tSlice a 3 dimensional slice to a scalar.";

    // Define slice region.
    DataArrayView::RegionType region;
    region.push_back(DataArrayView::RegionType::value_type(3,4));
    region.push_back(DataArrayView::RegionType::value_type(0,1));
    region.push_back(DataArrayView::RegionType::value_type(5,6));

    // Define shapes of views.
    DataArrayView::ShapeType sourceShape;
    sourceShape.push_back(6);
    sourceShape.push_back(3);
    sourceShape.push_back(13);
    DataArrayView::ShapeType targetShape;

    // Create source and target views.
    DataArray source(sourceShape);
    DataArrayView sourceView=source.getView();
    int val=0;
    for (int i=0;i<sourceShape[0];i++) {
      for (int j=0;j<sourceShape[1];j++) {
        for (int k=0;k<sourceShape[2];k++) {
	  sourceView(i,j,k)=val++;
        }
      }
    }

    DataArray target(targetShape);
    DataArrayView targetView=target.getView();

    // Copy source view to target view.
    targetView.copySlice(sourceView,region);

    // Check results of copy.
    for (int i=region[0].first;i<region[0].second;i++) {
      for (int j=region[1].first;j<region[1].second;j++) {
        for (int k=region[2].first;k<region[2].second;k++) {
	  assert(sourceView(i,j,k)==
                 targetView());
        }
      }
    }
  }

  {
    cout << endl;
    cout << "\tSlice a 3 dimensional slice from a 3 dimensional array.";

    // Define slice region.
    DataArrayView::RegionType region;
    region.push_back(DataArrayView::RegionType::value_type(4,7));
    region.push_back(DataArrayView::RegionType::value_type(2,5));
    region.push_back(DataArrayView::RegionType::value_type(6,7));

    // Define shapes of views.
    DataArrayView::ShapeType sourceShape = DataArrayView::getResultSliceShape(region);
    DataArrayView::ShapeType targetShape;
    targetShape.push_back(11);
    targetShape.push_back(8);
    targetShape.push_back(9);

    // Create source and target views.
    DataArray source(sourceShape);
    DataArrayView sourceView=source.getView();
    int val=0;
    for (int i=0;i<sourceShape[0];i++) {
      for (int j=0;j<sourceShape[1];j++) {
        for (int k=0;k<sourceShape[2];k++) {
	  sourceView(i,j,k)=val++;
        }
      }
    }

    DataArray target(targetShape);
    DataArrayView targetView=target.getView();

    // Copy source view to target view.
    targetView.copySliceFrom(sourceView,region);

    // Check results of copy.
    for (int i=region[0].first;i<region[0].second;i++) {
      for (int j=region[1].first;j<region[1].second;j++) {
        for (int k=region[2].first;k<region[2].second;k++) {
	  assert(sourceView(i-region[0].first,j-region[1].first,k-region[2].first)==
                 targetView(i,j,k));
        }
      }
    }
  }

  {
    cout << endl;
    cout << "\tSlice a 3 dimensional slice from a scalar.";

    // Define slice region.
    DataArrayView::RegionType region;
    region.push_back(DataArrayView::RegionType::value_type(4,7));
    region.push_back(DataArrayView::RegionType::value_type(2,5));
    region.push_back(DataArrayView::RegionType::value_type(6,7));

    // Define shapes of views.
    DataArrayView::ShapeType sourceShape;
    DataArrayView::ShapeType targetShape;
    targetShape.push_back(11);
    targetShape.push_back(8);
    targetShape.push_back(9);

    // Create source and target views.
    DataArray source(sourceShape);
    DataArrayView sourceView=source.getView();
    sourceView()=5;

    DataArray target(targetShape);
    DataArrayView targetView=target.getView();

    // Copy source view to target view.
    targetView.copySliceFrom(sourceView,region);

    // Check results of copy.
    for (int i=region[0].first;i<region[0].second;i++) {
      for (int j=region[1].first;j<region[1].second;j++) {
        for (int k=region[2].first;k<region[2].second;k++) {
	  assert(sourceView()==
                 targetView(i,j,k));
        }
      }
    }
  }

  {
    cout << endl;
    cout << "\tSlice a 4 dimensional slice to a 4 dimensional array.";

    // Define slice region.
    DataArrayView::RegionType region;
    region.push_back(DataArrayView::RegionType::value_type(2,4));
    region.push_back(DataArrayView::RegionType::value_type(0,2));
    region.push_back(DataArrayView::RegionType::value_type(5,9));
    region.push_back(DataArrayView::RegionType::value_type(3,5));

    // Define shapes of views.
    DataArrayView::ShapeType sourceShape;
    sourceShape.push_back(6);
    sourceShape.push_back(3);
    sourceShape.push_back(13);
    sourceShape.push_back(9);
    DataArrayView::ShapeType targetShape = DataArrayView::getResultSliceShape(region);

    // Create source and target views.
    DataArray source(sourceShape);
    DataArrayView sourceView=source.getView();
    int val=0;
    for (int i=0;i<sourceShape[0];i++) {
      for (int j=0;j<sourceShape[1];j++) {
        for (int k=0;k<sourceShape[2];k++) {
          for (int l=0;l<sourceShape[3];l++) {
	    sourceView(i,j,k,l)=val++;
          }
        }
      }
    }

    DataArray target(targetShape);
    DataArrayView targetView=target.getView();

    // Copy source view to target view.
    targetView.copySlice(sourceView,region);

    // Check results of copy.
    for (int i=region[0].first;i<region[0].second;i++) {
      for (int j=region[1].first;j<region[1].second;j++) {
        for (int k=region[2].first;k<region[2].second;k++) {
          for (int l=region[3].first;l<region[3].second;l++) {
	    assert(sourceView(i,j,k,l)==
                   targetView(i-region[0].first,j-region[1].first,k-region[2].first,l-region[3].first));
          }
        }
      }
    }
  }

  {
    cout << endl;
    cout << "\tSlice a 4 dimensional slice to a 3 dimensional array.";

    // Define slice region.
    DataArrayView::RegionType region;
    region.push_back(DataArrayView::RegionType::value_type(2,4));
    region.push_back(DataArrayView::RegionType::value_type(0,2));
    region.push_back(DataArrayView::RegionType::value_type(5,6));
    region.push_back(DataArrayView::RegionType::value_type(3,5));

    // Define shapes of views.
    DataArrayView::ShapeType sourceShape;
    sourceShape.push_back(6);
    sourceShape.push_back(3);
    sourceShape.push_back(13);
    sourceShape.push_back(9);
    DataArrayView::ShapeType targetShape;
    targetShape.push_back(2);
    targetShape.push_back(2);
    targetShape.push_back(2);

    // Create source and target views.
    DataArray source(sourceShape);
    DataArrayView sourceView=source.getView();
    int val=0;
    for (int i=0;i<sourceShape[0];i++) {
      for (int j=0;j<sourceShape[1];j++) {
        for (int k=0;k<sourceShape[2];k++) {
          for (int l=0;l<sourceShape[3];l++) {
	    sourceView(i,j,k,l)=val++;
          }
        }
      }
    }

    DataArray target(targetShape);
    DataArrayView targetView=target.getView();

    // Copy source view to target view.
    targetView.copySlice(sourceView,region);

    // Check results of copy.
    for (int i=region[0].first;i<region[0].second;i++) {
      for (int j=region[1].first;j<region[1].second;j++) {
        for (int k=region[2].first;k<region[2].second;k++) {
          for (int l=region[3].first;l<region[3].second;l++) {
	    assert(sourceView(i,j,k,l)==
                   targetView(i-region[0].first,j-region[1].first,l-region[3].first));
          }
        }
      }
    }
  }

  {
    cout << endl;
    cout << "\tSlice a 4 dimensional slice to a 2 dimensional array.";

    // Define slice region.
    DataArrayView::RegionType region;
    region.push_back(DataArrayView::RegionType::value_type(2,4));
    region.push_back(DataArrayView::RegionType::value_type(0,2));
    region.push_back(DataArrayView::RegionType::value_type(5,6));
    region.push_back(DataArrayView::RegionType::value_type(4,5));

    // Define shapes of views.
    DataArrayView::ShapeType sourceShape;
    sourceShape.push_back(6);
    sourceShape.push_back(3);
    sourceShape.push_back(13);
    sourceShape.push_back(9);
    DataArrayView::ShapeType targetShape;
    targetShape.push_back(2);
    targetShape.push_back(2);

    // Create source and target views.
    DataArray source(sourceShape);
    DataArrayView sourceView=source.getView();
    int val=0;
    for (int i=0;i<sourceShape[0];i++) {
      for (int j=0;j<sourceShape[1];j++) {
        for (int k=0;k<sourceShape[2];k++) {
          for (int l=0;l<sourceShape[3];l++) {
	    sourceView(i,j,k,l)=val++;
          }
        }
      }
    }

    DataArray target(targetShape);
    DataArrayView targetView=target.getView();

    // Copy source view to target view.
    targetView.copySlice(sourceView,region);

    // Check results of copy.
    for (int i=region[0].first;i<region[0].second;i++) {
      for (int j=region[1].first;j<region[1].second;j++) {
        for (int k=region[2].first;k<region[2].second;k++) {
          for (int l=region[3].first;l<region[3].second;l++) {
	    assert(sourceView(i,j,k,l)==
                   targetView(i-region[0].first,j-region[1].first));
          }
        }
      }
    }
  }

  {
    cout << endl;
    cout << "\tSlice a 4 dimensional slice to a 1 dimensional array.";

    // Define slice region.
    DataArrayView::RegionType region;
    region.push_back(DataArrayView::RegionType::value_type(3,4));
    region.push_back(DataArrayView::RegionType::value_type(0,2));
    region.push_back(DataArrayView::RegionType::value_type(5,6));
    region.push_back(DataArrayView::RegionType::value_type(4,5));

    // Define shapes of views.
    DataArrayView::ShapeType sourceShape;
    sourceShape.push_back(6);
    sourceShape.push_back(3);
    sourceShape.push_back(13);
    sourceShape.push_back(9);
    DataArrayView::ShapeType targetShape;
    targetShape.push_back(2);

    // Create source and target views.
    DataArray source(sourceShape);
    DataArrayView sourceView=source.getView();
    int val=0;
    for (int i=0;i<sourceShape[0];i++) {
      for (int j=0;j<sourceShape[1];j++) {
        for (int k=0;k<sourceShape[2];k++) {
          for (int l=0;l<sourceShape[3];l++) {
	    sourceView(i,j,k,l)=val++;
          }
        }
      }
    }

    DataArray target(targetShape);
    DataArrayView targetView=target.getView();

    // Copy source view to target view.
    targetView.copySlice(sourceView,region);

    // Check results of copy.
    for (int i=region[0].first;i<region[0].second;i++) {
      for (int j=region[1].first;j<region[1].second;j++) {
        for (int k=region[2].first;k<region[2].second;k++) {
          for (int l=region[3].first;l<region[3].second;l++) {
	    assert(sourceView(i,j,k,l)==
                   targetView(j-region[1].first));
          }
        }
      }
    }
  }

  {
    cout << endl;
    cout << "\tSlice a 4 dimensional slice to a scalar.";

    // Define slice region.
    DataArrayView::RegionType region;
    region.push_back(DataArrayView::RegionType::value_type(3,4));
    region.push_back(DataArrayView::RegionType::value_type(1,2));
    region.push_back(DataArrayView::RegionType::value_type(5,6));
    region.push_back(DataArrayView::RegionType::value_type(4,5));

    // Define shapes of views.
    DataArrayView::ShapeType sourceShape;
    sourceShape.push_back(6);
    sourceShape.push_back(3);
    sourceShape.push_back(13);
    sourceShape.push_back(9);
    DataArrayView::ShapeType targetShape;

    // Create source and target views.
    DataArray source(sourceShape);
    DataArrayView sourceView=source.getView();
    int val=0;
    for (int i=0;i<sourceShape[0];i++) {
      for (int j=0;j<sourceShape[1];j++) {
        for (int k=0;k<sourceShape[2];k++) {
          for (int l=0;l<sourceShape[3];l++) {
	    sourceView(i,j,k,l)=val++;
          }
        }
      }
    }

    DataArray target(targetShape);
    DataArrayView targetView=target.getView();

    // Copy source view to target view.
    targetView.copySlice(sourceView,region);

    // Check results of copy.
    for (int i=region[0].first;i<region[0].second;i++) {
      for (int j=region[1].first;j<region[1].second;j++) {
        for (int k=region[2].first;k<region[2].second;k++) {
          for (int l=region[3].first;l<region[3].second;l++) {
	    assert(sourceView(i,j,k,l)==
                   targetView());
          }
        }
      }
    }
  }

  {
    cout << endl;
    cout << "\tSlice a 4 dimensional slice from a 4 dimensional array.";

    // Define slice region.
    DataArrayView::RegionType region;
    region.push_back(DataArrayView::RegionType::value_type(14,37));
    region.push_back(DataArrayView::RegionType::value_type(22,57));
    region.push_back(DataArrayView::RegionType::value_type(63,71));
    region.push_back(DataArrayView::RegionType::value_type(23,51));

    // Define shapes of views.
    DataArrayView::ShapeType sourceShape = DataArrayView::getResultSliceShape(region);
    DataArrayView::ShapeType targetShape;
    targetShape.push_back(50);
    targetShape.push_back(65);
    targetShape.push_back(80);
    targetShape.push_back(90);

    // Create source and target views.
    DataArray source(sourceShape);
    DataArrayView sourceView=source.getView();
    int val=0;
    for (int i=0;i<sourceShape[0];i++) {
      for (int j=0;j<sourceShape[1];j++) {
        for (int k=0;k<sourceShape[2];k++) {
          for (int l=0;l<sourceShape[3];l++) {
	    sourceView(i,j,k,l)=val++;
          }
        }
      }
    }

    DataArray target(targetShape);
    DataArrayView targetView=target.getView();

    // Copy source view to target view.
    targetView.copySliceFrom(sourceView,region);

    // Check results of copy.
    for (int i=region[0].first;i<region[0].second;i++) {
      for (int j=region[1].first;j<region[1].second;j++) {
        for (int k=region[2].first;k<region[2].second;k++) {
          for (int l=region[3].first;l<region[3].second;l++) {
	    assert(sourceView(i-region[0].first,j-region[1].first,k-region[2].first,l-region[3].first)==
                   targetView(i,j,k,l));
          }
        }
      }
    }
  }

  {
    cout << endl;
    cout << "\tSlice a 4 dimensional slice from a scalar.";

    // Define slice region.
    DataArrayView::RegionType region;
    region.push_back(DataArrayView::RegionType::value_type(14,37));
    region.push_back(DataArrayView::RegionType::value_type(22,57));
    region.push_back(DataArrayView::RegionType::value_type(63,71));
    region.push_back(DataArrayView::RegionType::value_type(23,51));

    // Define shapes of views.
    DataArrayView::ShapeType sourceShape;
    DataArrayView::ShapeType targetShape;
    targetShape.push_back(50);
    targetShape.push_back(65);
    targetShape.push_back(80);
    targetShape.push_back(90);

    // Create source and target views.
    DataArray source(sourceShape);
    DataArrayView sourceView=source.getView();
    sourceView()=5;

    DataArray target(targetShape);
    DataArrayView targetView=target.getView();

    // Copy source view to target view.
    targetView.copySliceFrom(sourceView,region);

    // Check results of copy.
    for (int i=region[0].first;i<region[0].second;i++) {
      for (int j=region[1].first;j<region[1].second;j++) {
        for (int k=region[2].first;k<region[2].second;k++) {
          for (int l=region[3].first;l<region[3].second;l++) {
	    assert(sourceView()==
                   targetView(i,j,k,l));
          }
        }
      }
    }
  }

  cout << endl;

}

void DataArrayViewTestCase::testShapeToString() {

  cout << endl;
  cout << "\tTest shapeToString for a variety of shapes." << endl;

  DataArrayView::ShapeType shape;
  assert(DataArrayView::shapeToString(shape)=="()");
  shape.push_back(5);
  assert(DataArrayView::shapeToString(shape)=="(5)");
  shape.push_back(2);
  assert(DataArrayView::shapeToString(shape)=="(5,2)");
  shape.push_back(9);
  assert(DataArrayView::shapeToString(shape)=="(5,2,9)");
  shape.push_back(4);
  assert(DataArrayView::shapeToString(shape)=="(5,2,9,4)");

}

void DataArrayViewTestCase::testScalarView() {

  cout << endl;
  cout << "\tTest functionality of scalar views." << endl;

  // Create a vector containing enough data for three scalars
  // and check three scalar views return the appropriate data
  DataArrayView::ShapeType vShape;
  DataArrayView::ValueType vData(3);
  vData[0]=0;
  vData[1]=1;
  vData[2]=2;

  // create the three scalar views
  DataArrayView zView(vData,vShape,0);
  DataArrayView oView(vData,vShape,1);
  DataArrayView tView(vData,vShape,2);

  // check attributes of the three scalar views
  assert(zView()==0);
  assert(oView()==1);
  assert(tView()==2);
  assert(zView.noValues()==1);
  assert(oView.noValues()==1);
  assert(tView.noValues()==1);
  assert(zView.getRank()==0);
  assert(oView.getRank()==0);
  assert(tView.getRank()==0);

  // copy the one view to the zero view
  zView.copy(oView);
  assert(zView==oView);
  zView.checkShape(oView.getShape());

  // create a single vector view of all the data
  vShape.push_back(3);
  DataArrayView oneVView(vData,vShape,0);

  // test shape mismatch functions
  if (!zView.checkShape(oneVView.getShape())) {
    assert(true);
  } else {
    assert(false);
  }

  // test some unary ops
  zView.unaryOp(negate<double>());
  assert(zView()==-1);
  zView.binaryOp(oView,plus<double>());
  assert(zView()==0);

  // test operator !=
  assert(zView!=oView);

}

void DataArrayViewTestCase::testAll()
{

  {
    cout << endl;
    cout << "\tTest empty DataArrayView.";

    // default constructor
    DataArrayView defView;

    // check all attributes
    assert(defView.isEmpty());
    assert(defView.getOffset()==0);
    assert(defView.getRank()==0);
    assert(defView.noValues()==0);
    assert(defView.getShape().empty());
    assert(defView.checkShape(DataArrayView::ShapeType()));
  }

  {
    cout << endl;
    cout << "\tTest DataArrayView - shape (5).";

    // define the shape for the DataArrayView
    DataArrayView::ShapeType shape;
    shape.push_back(5);

    // allocate the data for the DataArrayView
    DataArrayView::ValueType data(DataArrayView::noValues(shape),0);

    // constructor
    int offset=0;
    DataArrayView dataView(data,shape,offset);

    // check all attributes
    assert(!dataView.isEmpty());
    assert(dataView.getOffset()==0);
    assert(dataView.getRank()==1);
    assert(dataView.noValues()==5);
    assert(dataView.getShape()==shape);
    assert(dataView.checkShape(shape));

    // assign values to the data
    for (int i=0;i<shape[0];i++) {
      dataView(i)=dataView.index(i);
      assert(dataView(i)==i);
    }

  }

  {
    cout << endl;
    cout << "\tTest DataArrayView - shape (2,3).";

    // define the shape for the DataArrayView
    DataArrayView::ShapeType shape;
    shape.push_back(2);
    shape.push_back(3);

    // allocate the data for the DataArrayView
    int npoints=4;
    DataArrayView::ValueType data(DataArrayView::noValues(shape)*npoints,0);

    // constructor
    DataArrayView dataView(data,shape);

    // check all attributes
    assert(!dataView.isEmpty());
    assert(dataView.getOffset()==0);
    assert(dataView.getRank()==2);
    assert(dataView.noValues()==6);
    assert(dataView.getShape()==shape);
    assert(dataView.checkShape(shape));

    // step the view along each block of this shape in the underlying data
    for (int p=0;p<npoints;p++) {

      assert(dataView.getOffset()==p*dataView.noValues());

      // assign values to the data
      for (int i=0;i<shape[0];i++) {
        for (int j=0;j<shape[1];j++) {
          dataView(i,j)=dataView.index(i,j);
          assert(dataView(i,j)==dataView.index(i,j));
        }
      }

      if (p<npoints-1) {
        dataView.incrOffset();
      }

    }

  }

  {
    cout << endl;
    cout << "\tTest DataArrayView - shape (4,7,9).";

    // define the shape for the DataArrayView
    DataArrayView::ShapeType shape;
    shape.push_back(4);
    shape.push_back(7);
    shape.push_back(9);

    // allocate the data for the DataArrayView
    int npoints=10;
    DataArrayView::ValueType data(DataArrayView::noValues(shape)*npoints,0);

    // constructor
    DataArrayView dataView(data,shape);

    // check all attributes
    assert(!dataView.isEmpty());
    assert(dataView.getOffset()==0);
    assert(dataView.getRank()==3);
    assert(dataView.noValues()==252);
    assert(dataView.getShape()==shape);
    assert(dataView.checkShape(shape));

    // step the view along each block of this shape in the underlying data
    for (int p=0;p<npoints;p++) {

      assert(dataView.getOffset()==p*dataView.noValues());

      // assign values to the data
      for (int i=0;i<shape[0];i++) {
        for (int j=0;j<shape[1];j++) {
          for (int k=0;k<shape[2];k++) {
            dataView(i,j,k)=dataView.index(i,j,k);
            assert(dataView(i,j,k)==dataView.index(i,j,k));
          }
        }
      }

      if (p<npoints-1) {
        dataView.incrOffset();
      }

    }

  }

  {
    cout << endl;
    cout << "\tTest DataArrayView - shape (12,4,5,14).";

    // define the shape for the DataArrayView
    DataArrayView::ShapeType shape;
    shape.push_back(12);
    shape.push_back(4);
    shape.push_back(5);
    shape.push_back(14);

    // allocate the data for the DataArrayView
    int npoints=100;
    DataArrayView::ValueType data(DataArrayView::noValues(shape)*npoints,0);

    // constructor
    DataArrayView dataView(data,shape);

    // check all attributes
    assert(!dataView.isEmpty());
    assert(dataView.getOffset()==0);
    assert(dataView.getRank()==4);
    assert(dataView.noValues()==3360);
    assert(dataView.getShape()==shape);
    assert(dataView.checkShape(shape));

    // step the view along each block of this shape in the underlying data
    for (int p=0;p<npoints;p++) {

      assert(dataView.getOffset()==p*dataView.noValues());

      // assign values to the data
      for (int i=0;i<shape[0];i++) {
        for (int j=0;j<shape[1];j++) {
          for (int k=0;k<shape[2];k++) {
            for (int l=0;l<shape[3];l++) {
              dataView(i,j,k,l)=dataView.index(i,j,k,l);
              assert(dataView(i,j,k,l)==dataView.index(i,j,k,l));
            }
          }
        }
      }

      if (p<npoints-1) {
        dataView.incrOffset();
      }

    }

  }

  {
    cout << endl;
    cout << "\tTest DataArrayView copy constructor - shape (5).";

    // define the shape for the DataArrayView
    DataArrayView::ShapeType shape;
    shape.push_back(5);

    // allocate the data for the DataArrayView
    DataArrayView::ValueType data(DataArrayView::noValues(shape),0);

    // constructor
    DataArrayView dataView(data,shape);

    // copy constructor
    DataArrayView dataViewCopy(dataView);

    // check all attributes
    assert(!dataViewCopy.isEmpty());
    assert(dataViewCopy.getOffset()==0);
    assert(dataViewCopy.getRank()==1);
    assert(dataViewCopy.noValues()==5);
    assert(dataViewCopy.getShape()==shape);
    assert(dataViewCopy.checkShape(shape));

    // check data
    assert(dataView.getData()==dataViewCopy.getData());
    for (int i=0;i<dataView.getData().size();i++) {
      assert(dataView.getData(i)==dataViewCopy.getData(i));
    }

  }

  {
    cout << endl;
    cout << "\tTest DataArrayView copy constructor - shape (5,6,7).";

    // define the shape for the DataArrayView
    DataArrayView::ShapeType shape;
    shape.push_back(5);
    shape.push_back(6);
    shape.push_back(7);

    // allocate the data for the DataArrayView
    DataArrayView::ValueType data(DataArrayView::noValues(shape),0);

    // constructor
    DataArrayView dataView(data,shape);

    // copy constructor
    DataArrayView dataViewCopy(dataView);

    // check all attributes
    assert(!dataViewCopy.isEmpty());
    assert(dataViewCopy.getOffset()==0);
    assert(dataViewCopy.getRank()==3);
    assert(dataViewCopy.noValues()==210);
    assert(dataViewCopy.getShape()==shape);
    assert(dataViewCopy.checkShape(shape));

    // check data
    assert(dataView.getData()==dataViewCopy.getData());
    for (int i=0;i<dataView.getData().size();i++) {
      assert(dataView.getData(i)==dataViewCopy.getData(i));
    }

  }

  {
    cout << endl;
    cout << "\tTest DataArrayView copy method - shape (2,3).";

    // define the shape for the DataArrayViews
    DataArrayView::ShapeType shape;
    shape.push_back(2);
    shape.push_back(3);

    // allocate the data for the DataArrayViews
    int npoints=4;
    DataArrayView::ValueType data1(DataArrayView::noValues(shape)*npoints,0);
    DataArrayView::ValueType data2(DataArrayView::noValues(shape)*npoints,0);

    // construct two views
    DataArrayView dataView1(data1,shape);
    DataArrayView dataView2(data2,shape);

    // step the views along each block of this shape in the underlying data arrays
    for (int p=0;p<npoints;p++) {

      // assign values to the data underlying the first view
      for (int i=0;i<shape[0];i++) {
        for (int j=0;j<shape[1];j++) {
          dataView1(i,j)=dataView1.index(i,j);
        }
      }

      // copy data from the first view to the second
      dataView2.copy(dataView1);

      // check the data underlying the second view
      for (int i=0;i<shape[0];i++) {
        for (int j=0;j<shape[1];j++) {
          assert(dataView2(i,j)==dataView1.index(i,j));
        }
      }

      if (p<npoints-1) {
        dataView1.incrOffset();
        dataView2.incrOffset();
      }

    }

  }

  {
    cout << endl;
    cout << "\tTest DataArrayView copy method - shape (2,3,8,9).";

    // define the shape for the DataArrayViews
    DataArrayView::ShapeType shape;
    shape.push_back(2);
    shape.push_back(3);
    shape.push_back(8);
    shape.push_back(9);

    // allocate the data for the DataArrayViews
    int npoints=10;
    DataArrayView::ValueType data1(DataArrayView::noValues(shape)*npoints,0);
    DataArrayView::ValueType data2(DataArrayView::noValues(shape)*npoints,0);

    // construct two views
    DataArrayView dataView1(data1,shape);
    DataArrayView dataView2(data2,shape);

    // step the views along each block of this shape in the underlying data arrays
    for (int p=0;p<npoints;p++) {

      // assign values to the data underlying the first view
      for (int i=0;i<shape[0];i++) {
        for (int j=0;j<shape[1];j++) {
          for (int k=0;k<shape[2];k++) {
            for (int l=0;l<shape[3];l++) {
              dataView1(i,j,k,l)=dataView1.index(i,j,k,l);
            }
          }
        }
      }

      // copy data from the first view to the second
      dataView2.copy(dataView1);

      // check the data underlying the second view
      for (int i=0;i<shape[0];i++) {
        for (int j=0;j<shape[1];j++) {
          for (int k=0;k<shape[2];k++) {
            for (int l=0;l<shape[3];l++) {
              assert(dataView2(i,j,k,l)==dataView1.index(i,j,k,l));
            }
          }
        }
      }

      if (p<npoints-1) {
        dataView1.incrOffset();
        dataView2.incrOffset();
      }

    }

  }

  {
    cout << endl;
    cout << "\tTest DataArrayView copy with offset method - shape (2,3).";

    // define the shape for the DataArrayViews
    DataArrayView::ShapeType shape;
    shape.push_back(2);
    shape.push_back(3);

    // allocate the data for the DataArrayViews
    int npoints=4;
    DataArrayView::ValueType data1(DataArrayView::noValues(shape)*npoints,0);
    DataArrayView::ValueType data2(DataArrayView::noValues(shape)*npoints,0);

    // construct two views
    DataArrayView dataView1(data1,shape);
    DataArrayView dataView2(data2,shape);

    // assign values to the data underlying the first view
    for (int i=0;i<shape[0];i++) {
      for (int j=0;j<shape[1];j++) {
        dataView1(i,j)=dataView1.index(i,j);
      }
    }

    // copy data from the first view to the second at an offset
    dataView2.copy(12,dataView1,0);

    // check the data underlying the second view
    dataView2.setOffset(12);
    for (int i=0;i<shape[0];i++) {
      for (int j=0;j<shape[1];j++) {
        assert(dataView2(i,j)==dataView1.index(i,j));
      }
    }

  }

  {
    cout << endl;
    cout << "\tTest DataArrayView copy with value method - shape (5,8).";

    // define the shape for the DataArrayView
    DataArrayView::ShapeType shape;
    shape.push_back(5);
    shape.push_back(8);

    // allocate the data for the DataArrayView
    int npoints=4;
    DataArrayView::ValueType data(DataArrayView::noValues(shape)*npoints,0);

    // construct view
    DataArrayView dataView(data,shape);

    // copy a value to the view at an offset
    dataView.copy(80,5);

    // check the data underlying the second view
    dataView.setOffset(80);
    for (int i=0;i<shape[0];i++) {
      for (int j=0;j<shape[1];j++) {
        assert(dataView(i,j)==5);
      }
    }

  }

  #if defined DOASSERT
  {
    cout << endl;
    cout << "\tTest too many indices for shape exception.";

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
      assert(true);
    }

    try {
      dataView(1,1)=1;
      assert(false);
    }
    catch (EsysException& e) {
      assert(true);
    }

    try {
      dataView(1,1,1)=1;
      assert(false);
    }
    catch (EsysException& e) {
      assert(true);
    }

    try {
      dataView(1,1,1,1)=1;
      assert(false);
    }
    catch (EsysException& e) {
      assert(true);
    }

  }
  #endif

  #if defined DOASSERT
  {
    cout << endl;
    cout << "\tTest invalid index exception.";

    DataArrayView::ShapeType shape;
    shape.push_back(4);
    DataArrayView::ValueType data(DataArrayView::noValues(shape),0);
    DataArrayView dataView(data,shape);

    try {
      dataView(4000)=1;
      assert(false);
    }
    catch (EsysException& e) {
      assert(true);
    }
  }
  #endif

  {
    cout << endl;
    cout << "\tTest insufficient data exception.";

    DataArrayView::ShapeType shape;
    DataArrayView::ValueType data;

    try {
      DataArrayView dataView(data,shape);
      assert(false);
    }
    catch (EsysException& e) {
      assert(true);
    }
  }

  cout << endl;

}

void DataArrayViewTestCase::testMatMult()
{

  {
    cout << endl;
    cout << "\tTest result shape." << endl;

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

    assert(resultShape.size()==2);
    assert(resultShape[0]==1);
    assert(resultShape[1]==2);

    DataArrayView::ValueType resultData(DataArrayView::noValues(resultShape),0);
    DataArrayView resultDataView(resultData,resultShape);
  
    cout << "\tTest matrix multiplication.";
    double aValue=0.0;
    for (int i=0;i<leftShape[0];i++) {
      for (int j=0;j<leftShape[1];j++) {
	leftDataView(i,j)=++aValue;
      }
    }
    aValue=0.0;
    for (int i=0;i<rightShape[0];i++) {
      for (int j=0;j<rightShape[1];j++) {
	rightDataView(i,j)=++aValue;
      }
    }

    DataArrayView::matMult(leftDataView,rightDataView,resultDataView);
  }

  cout << endl;

}

void DataArrayViewTestCase::testUnaryOp()
{

  // This typedef allows function names to be cast to pointers
  // to unary functions of the appropriate type.
  typedef double (*UnaryDFunPtr)(double);

  {
    cout << endl;
    cout << "\tTest unaryOp on scalar DataArrayView.";

    // define the shape for the DataArrayView
    DataArrayView::ShapeType shape;

    // allocate the data for the DataArrayView
    int npoints=4;
    DataArrayView::ValueType data(DataArrayView::noValues(shape)*npoints,0);

    // constructor
    DataArrayView dataView(data,shape);

    double tmp;
    // step the view along each data point in the underlying data
    for (int p=0;p<npoints;p++) {

      // assign values to the data point
      dataView()=p;

      // apply a unary operation to this data point
      dataView.unaryOp((UnaryDFunPtr)std::sin);

      // check the results
      tmp = std::sin((double)p);
      assert(std::abs(dataView()-tmp)<=REL_TOL*std::abs(tmp));

      if (p<npoints-1) {
        dataView.incrOffset();
      }

    }

  }

  {
    cout << endl;
    cout << "\tTest unaryOp on shape (2,3) DataArrayView.";

    // define the shape for the DataArrayView
    DataArrayView::ShapeType shape;
    shape.push_back(2);
    shape.push_back(3);

    // allocate the data for the DataArrayView
    int npoints=4;
    DataArrayView::ValueType data(DataArrayView::noValues(shape)*npoints,0);

    // constructor
    DataArrayView dataView(data,shape);

    // step the view along each data point in the underlying data
    for (int p=0;p<npoints;p++) {

      // assign values to the data point
      for (int i=0;i<shape[0];i++) {
        for (int j=0;j<shape[1];j++) {
          dataView(i,j)=dataView.index(i,j);
        }
      }

      // apply a unary operation to this data point
      dataView.unaryOp((UnaryDFunPtr)std::sqrt);

      // check the results
      for (int i=0;i<shape[0];i++) {
        for (int j=0;j<shape[1];j++) {
          assert(std::abs(dataView(i,j)-std::sqrt((double)dataView.index(i,j)))<=REL_TOL*std::sqrt((double)dataView.index(i,j)));
        }
      }

      if (p<npoints-1) {
        dataView.incrOffset();
      }

    }

  }

  {
    cout << endl;
    cout << "\tTest unaryOp on shape (9,8,5,11) DataArrayView.";

    // define the shape for the DataArrayView
    DataArrayView::ShapeType shape;
    shape.push_back(9);
    shape.push_back(8);
    shape.push_back(5);
    shape.push_back(11);

    // allocate the data for the DataArrayView
    int npoints=4;
    DataArrayView::ValueType data(DataArrayView::noValues(shape)*npoints,0);

    // constructor
    DataArrayView dataView(data,shape);

    // step the view along each data point in the underlying data
    for (int p=0;p<npoints;p++) {

      // assign values to the data point
      for (int i=0;i<shape[0];i++) {
        for (int j=0;j<shape[1];j++) {
          for (int k=0;k<shape[2];k++) {
            for (int l=0;l<shape[3];l++) {
              dataView(i,j,k,l)=dataView.index(i,j,k,l)+1;
            }
          }
        }
      }

      // apply a unary operation to this data point
      dataView.unaryOp((UnaryDFunPtr)std::log);

      // check the results
      for (int i=0;i<shape[0];i++) {
        for (int j=0;j<shape[1];j++) {
          for (int k=0;k<shape[2];k++) {
            for (int l=0;l<shape[3];l++) {
              assert(std::abs(dataView(i,j,k,l)-std::log(1+(double)dataView.index(i,j,k,l)))<=REL_TOL*std::abs(std::log(1+(double)dataView.index(i,j,k,l))));
            }
          }
        }
      }

      if (p<npoints-1) {
        dataView.incrOffset();
      }

    }

  }

  cout << endl;

}

void DataArrayViewTestCase::testBinaryOp()
{

  // This typedef allows function names to be cast to pointers
  // to binary functions of the appropriate type.
  typedef double (*BinaryDFunPtr)(double,double);

  {
    cout << endl;
    cout << "\tTest binaryOp on scalar DataArrayViews.";

    // define the shape for the DataArrayViews
    DataArrayView::ShapeType shape;

    // allocate the data for the DataArrayViews
    int npoints=4;
    DataArrayView::ValueType data1(DataArrayView::noValues(shape)*npoints,0);
    DataArrayView::ValueType data2(DataArrayView::noValues(shape)*npoints,0);

    // constructor
    DataArrayView dataView1(data1,shape);
    DataArrayView dataView2(data2,shape);

    // step the views along each data point in the underlying data
    for (int p=0;p<npoints;p++) {

      // assign values to the data points
      dataView1()=p;
      dataView2()=p;

      // apply a binary operation to these data points
      dataView1.binaryOp(dataView2,plus<double>());

      // check the results
      assert(dataView1()==p+p);

      if (p<npoints-1) {
        dataView1.incrOffset();
        dataView2.incrOffset();
      }

    }

  }

  {
    cout << endl;
    cout << "\tTest binaryOp on shape (2,3) DataArrayViews.";

    // define the shape for the DataArrayViews
    DataArrayView::ShapeType shape;
    shape.push_back(2);
    shape.push_back(3);

    // allocate the data for the DataArrayViews
    int npoints=4;
    DataArrayView::ValueType data1(DataArrayView::noValues(shape)*npoints,0);
    DataArrayView::ValueType data2(DataArrayView::noValues(shape)*npoints,0);

    // constructor
    DataArrayView dataView1(data1,shape);
    DataArrayView dataView2(data2,shape);

    // step the views along each data point in the underlying data
    for (int p=0;p<npoints;p++) {

      // assign values to the data points
      for (int i=0;i<shape[0];i++) {
        for (int j=0;j<shape[1];j++) {
          dataView1(i,j)=dataView1.index(i,j);
          dataView2(i,j)=dataView2.index(i,j);
        }
      }

      // apply a binary operation to these data points
      dataView1.binaryOp(dataView2,multiplies<double>());

      // check the results
      for (int i=0;i<shape[0];i++) {
        for (int j=0;j<shape[1];j++) {
          assert(dataView1(i,j)==dataView1.index(i,j)*dataView2.index(i,j));
        }
      }

      if (p<npoints-1) {
        dataView1.incrOffset();
        dataView2.incrOffset();
      }

    }

  }

  {
    cout << endl;
    cout << "\tTest binaryOp on shape (9,8,5,11) DataArrayViews.";

    // define the shape for the DataArrayViews
    DataArrayView::ShapeType shape;
    shape.push_back(9);
    shape.push_back(8);
    shape.push_back(5);
    shape.push_back(11);

    // allocate the data for the DataArrayViews
    int npoints=4;
    DataArrayView::ValueType data1(DataArrayView::noValues(shape)*npoints,0);
    DataArrayView::ValueType data2(DataArrayView::noValues(shape)*npoints,0);

    // constructor
    DataArrayView dataView1(data1,shape);
    DataArrayView dataView2(data2,shape);

    // step the views along each data point in the underlying data
    for (int p=0;p<npoints;p++) {

      // assign values to the data points
      for (int i=0;i<shape[0];i++) {
        for (int j=0;j<shape[1];j++) {
          for (int k=0;k<shape[2];k++) {
            for (int l=0;l<shape[3];l++) {
              dataView1(i,j,k,l)=dataView1.index(i,j,k,l);
              dataView2(i,j,k,l)=dataView2.index(i,j,k,l);
            }
          }
        }
      }

      // apply a binary operation to these data points
      dataView1.binaryOp(dataView2,multiplies<double>());

      // check the results
      for (int i=0;i<shape[0];i++) {
        for (int j=0;j<shape[1];j++) {
          for (int k=0;k<shape[2];k++) {
            for (int l=0;l<shape[3];l++) {
              assert(dataView1(i,j,k,l)==dataView1.index(i,j,k,l)*dataView2.index(i,j,k,l));
            }
          }
        }
      }

      if (p<npoints-1) {
        dataView1.incrOffset();
        dataView2.incrOffset();
      }

    }

  }

  {
    cout << endl;
    cout << "\tTest binaryOp on scalar DataArrayView and single value.";

    // define the shape for the DataArrayView
    DataArrayView::ShapeType shape;

    // allocate the data for the DataArrayView
    int npoints=4;
    DataArrayView::ValueType data(DataArrayView::noValues(shape)*npoints,0);

    // constructor
    DataArrayView dataView(data,shape);

    // step the view along each data point in the underlying data
    for (int p=0;p<npoints;p++) {

      // assign values to the data point
      dataView()=p;

      // apply a binary operation to this data point
      dataView.binaryOp(4.9,plus<double>());

      // check the results
      assert(dataView()==4.9+p);

      if (p<npoints-1) {
        dataView.incrOffset();
      }

    }

  }

  {
    cout << endl;
    cout << "\tTest binaryOp on shape (2,3) DataArrayView and single value.";

    // define the shape for the DataArrayView
    DataArrayView::ShapeType shape;
    shape.push_back(2);
    shape.push_back(3);

    // allocate the data for the DataArrayView
    int npoints=4;
    DataArrayView::ValueType data(DataArrayView::noValues(shape)*npoints,0);

    // constructor
    DataArrayView dataView(data,shape);

    // step the view along each data point in the underlying data
    for (int p=0;p<npoints;p++) {

      // assign values to the data point
      for (int i=0;i<shape[0];i++) {
        for (int j=0;j<shape[1];j++) {
          dataView(i,j)=dataView.index(i,j);
        }
      }

      // apply a binary operation to the data point
      dataView.binaryOp(5.8,multiplies<double>());

      double tmp;
      // check the results
      for (int i=0;i<shape[0];i++) {
        for (int j=0;j<shape[1];j++) {
          tmp=5.8*dataView.index(i,j);
          assert(std::abs(dataView(i,j)-tmp)<=REL_TOL*std::abs(tmp));
        }
      }

      if (p<npoints-1) {
        dataView.incrOffset();
      }

    }

  }

  {
    cout << endl;
    cout << "\tTest binaryOp on shape (9,8,5,11) DataArrayView and single value.";

    // define the shape for the DataArrayView
    DataArrayView::ShapeType shape;
    shape.push_back(9);
    shape.push_back(8);
    shape.push_back(5);
    shape.push_back(11);

    // allocate the data for the DataArrayView
    int npoints=4;
    DataArrayView::ValueType data(DataArrayView::noValues(shape)*npoints,0);

    // constructor
    DataArrayView dataView(data,shape);

    // step the view along each data point in the underlying data
    for (int p=0;p<npoints;p++) {

      // assign values to the data point
      for (int i=0;i<shape[0];i++) {
        for (int j=0;j<shape[1];j++) {
          for (int k=0;k<shape[2];k++) {
            for (int l=0;l<shape[3];l++) {
              dataView(i,j,k,l)=dataView.index(i,j,k,l);
            }
          }
        }
      }

      // apply a binary operation to the data point
      dataView.binaryOp(5.4,multiplies<double>());

      double tmp;
      // check the results
      for (int i=0;i<shape[0];i++) {
        for (int j=0;j<shape[1];j++) {
          for (int k=0;k<shape[2];k++) {
            for (int l=0;l<shape[3];l++) {
              tmp=5.4*dataView.index(i,j,k,l);
              assert(std::abs(dataView(i,j,k,l)-tmp)<=REL_TOL*std::abs(tmp));
            }
          }
        }
      }

      if (p<npoints-1) {
        dataView.incrOffset();
      }

    }

  }

  cout << endl;

}

void DataArrayViewTestCase::testReductionOp()
{

  {
    cout << endl;
    cout << "\tTest reductionOp on scalar DataArrayView.";

    // define the shape for the DataArrayView
    DataArrayView::ShapeType shape;

    // allocate the data for the DataArrayView
    int npoints=4;
    DataArrayView::ValueType data(DataArrayView::noValues(shape)*npoints,0);

    // constructor
    DataArrayView dataView(data,shape);

    // step the view along each data point in the underlying data
    for (int p=0;p<npoints;p++) {

      // assign values to the data point
      dataView()=p;

      // apply a reduction operation to this data point and check the results
      FMax fmax_func;
      assert(std::abs(dataView.reductionOp(fmax_func,numeric_limits<double>::max()*-1)-p)<=REL_TOL*p);

      if (p<npoints-1) {
        dataView.incrOffset();
      }

    }

  }

  {
    cout << endl;
    cout << "\tTest reductionOp on shape (2,3) DataArrayView.";

    // define the shape for the DataArrayView
    DataArrayView::ShapeType shape;
    shape.push_back(2);
    shape.push_back(3);

    // allocate the data for the DataArrayView
    int npoints=4;
    DataArrayView::ValueType data(DataArrayView::noValues(shape)*npoints,0);

    // constructor
    DataArrayView dataView(data,shape);

    // step the view along each data point in the underlying data
    for (int p=0;p<npoints;p++) {

      // assign values to the data point
      for (int i=0;i<shape[0];i++) {
        for (int j=0;j<shape[1];j++) {
          dataView(i,j)=dataView.index(i,j);
        }
      }

      // apply a reduction operation to this data point and check the results
      FMin fmin_func;
      assert(std::abs(dataView.reductionOp(fmin_func,numeric_limits<double>::max())-dataView.index(0,0))<=REL_TOL*std::abs(dataView.index(0,0)));

      if (p<npoints-1) {
        dataView.incrOffset();
      }

    }

  }

  {
    cout << endl;
    cout << "\tTest reductionOp on shape (9,8,5,11) DataArrayView.";

    // define the shape for the DataArrayView
    DataArrayView::ShapeType shape;
    shape.push_back(9);
    shape.push_back(8);
    shape.push_back(5);
    shape.push_back(11);

    // allocate the data for the DataArrayView
    int npoints=4;
    DataArrayView::ValueType data(DataArrayView::noValues(shape)*npoints,0);

    // constructor
    DataArrayView dataView(data,shape);

    // step the view along each data point in the underlying data
    for (int p=0;p<npoints;p++) {

      // assign values to the data point
      for (int i=0;i<shape[0];i++) {
        for (int j=0;j<shape[1];j++) {
          for (int k=0;k<shape[2];k++) {
            for (int l=0;l<shape[3];l++) {
              dataView(i,j,k,l)=dataView.index(i,j,k,l);
            }
          }
        }
      }

      // apply a reduction operation to this data point and check the results
      AbsMax absmax_func;
      assert(dataView.reductionOp(absmax_func,0)==dataView.index(8,7,4,10));

      if (p<npoints-1) {
        dataView.incrOffset();
      }

    }

  }

  cout << endl;

}

TestSuite* DataArrayViewTestCase::suite ()
{
  //
  // create the suite of tests to perform.
  TestSuite *testSuite = new TestSuite ("DataArrayViewTestCase");

  testSuite->addTest (new TestCaller< DataArrayViewTestCase>("testAll",&DataArrayViewTestCase::testAll));
  testSuite->addTest (new TestCaller< DataArrayViewTestCase>("testShapeToString",&DataArrayViewTestCase::testShapeToString));
  testSuite->addTest (new TestCaller< DataArrayViewTestCase>("testScalarView",&DataArrayViewTestCase::testScalarView));
  testSuite->addTest (new TestCaller< DataArrayViewTestCase>("testResultSliceShape",&DataArrayViewTestCase::testResultSliceShape));
  testSuite->addTest (new TestCaller< DataArrayViewTestCase>("testSlicing",&DataArrayViewTestCase::testSlicing));
  testSuite->addTest (new TestCaller< DataArrayViewTestCase>("testUnaryOp",&DataArrayViewTestCase::testUnaryOp));
  testSuite->addTest (new TestCaller< DataArrayViewTestCase>("testBinaryOp",&DataArrayViewTestCase::testBinaryOp));
  testSuite->addTest (new TestCaller< DataArrayViewTestCase>("testReductionOp",&DataArrayViewTestCase::testReductionOp));
  testSuite->addTest (new TestCaller< DataArrayViewTestCase>("testMatMult",&DataArrayViewTestCase::testMatMult));
  return testSuite;
}
