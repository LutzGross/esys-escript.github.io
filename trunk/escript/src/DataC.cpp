// $Id$
/*
 ******************************************************************************
 *                                                                            *
 *       COPYRIGHT  ACcESS 2004 -  All Rights Reserved                        *
 *                                                                            *
 * This software is the property of ACcESS. No part of this code              *
 * may be copied in any form or by any means without the expressed written    *
 * consent of ACcESS.  Copying, use or modification of this software          *
 * by any unauthorised person is illegal unless that person has a software    *
 * license agreement with ACcESS.                                             *
 *                                                                            *
 ******************************************************************************
*/

extern "C" {
#include "DataC.h"
}

#include "Data.h"
#include "DataArrayView.h"

int getFunctionSpaceType(struct escriptDataC* data) 
{
  escript::Data* temp=(escript::Data*)(data->m_dataPtr);
  return temp->getFunctionSpace().getTypeCode();
}


int isDataPointShapeEqual(struct escriptDataC* data, int rank, int* dimensions)
{
  if (data == (struct escriptDataC*)0) {
       return true;
  } else {
     escript::Data* temp=(escript::Data*)(data->m_dataPtr);
     if (temp->isEmpty()) {
        return true;
     } else {
          escript::DataArrayView::ShapeType givenShape(&dimensions[0],&dimensions[rank]);
          return (temp->getPointDataView().getShape()==givenShape);
     }
  }
}

int numSamplesEqual(struct escriptDataC* data, int numDataPointsPerSample,
		    int numSamples)
{
  if (data == (struct escriptDataC*)0) {
       return true;
  } else {
     escript::Data* temp=(escript::Data*)(data->m_dataPtr);
     if (temp->isEmpty()) {
        return true;
     } else {
        int result=(numDataPointsPerSample==temp->getNumDataPointsPerSample());
        result=result && (numSamples==temp->getNumSamples());
        return result;
     }
  }
}

int getDataPointRank(struct escriptDataC* data)
{
  if (data == (struct escriptDataC*)0) {
       return 0;
  } else {
       escript::Data* temp=(escript::Data*)(data->m_dataPtr);
       return temp->getDataPointRank();
  }
}

int getDataPointShape(struct escriptDataC* data,int i)
{
  if (data == (struct escriptDataC*)0) {
       return 0;
  } else {
     escript::Data* temp=(escript::Data*)(data->m_dataPtr);
     int rank = temp->getDataPointRank();
     if (i<0 || i>=rank) {
        return 1;
     } else {
        const escript::DataArrayView::ShapeType view=temp->getDataPointShape();
        return view[i];
     }
  }
}

int getDataPointSize(struct escriptDataC* data)
{
  escript::Data* temp=(escript::Data*)(data->m_dataPtr);
  return temp->getDataPointSize();
}

int getLength(struct escriptDataC* data)
{
  escript::Data* temp=(escript::Data*)(data->m_dataPtr);
  return temp->getLength();
}

int isExpanded(struct escriptDataC* data)
{
  if (data == (struct escriptDataC*)0) {
       return false;
  } else {
     escript::Data* temp=(escript::Data*)(data->m_dataPtr);
     if (temp->isEmpty()) {
        return false;
     } else {
        return temp->isExpanded();
     }
  }
}

int isEmpty(escriptDataC* data) 
{
  if (data == (struct escriptDataC*)0) {
       return true;
  } else {
      escript::Data* temp=(escript::Data*)(data->m_dataPtr);
      return temp->isEmpty();
  }
}


double* getSampleData(struct escriptDataC* data, int sampleNo)
{
  if (data == (struct escriptDataC*)0) {
       return NULL;
  } else {
      escript::Data* temp=(escript::Data*)(data->m_dataPtr);
     if (temp->isEmpty()) {
        return NULL;
     } else {
        return temp->getSampleData(sampleNo);
     }
  }
}
