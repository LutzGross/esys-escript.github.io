// $Id$
/*
 ************************************************************
 *          Copyright 2006 by ACcESS MNRF                   *
 *                                                          *
 *              http://www.access.edu.au                    *
 *       Primary Business: Queensland, Australia            *
 *  Licensed under the Open Software License version 3.0    *
 *     http://www.opensource.org/licenses/osl-3.0.php       *
 *                                                          *
 ************************************************************
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
	printf("ksteube NumDataPointsPerSample=%d NumSamples=%d\n", temp->getNumDataPointsPerSample(), temp->getNumSamples());
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
