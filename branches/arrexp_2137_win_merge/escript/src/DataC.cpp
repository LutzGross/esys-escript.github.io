
/*******************************************************
*
* Copyright (c) 2003-2008 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


extern "C" {
#include "DataC.h"
}

#include "Data.h"
#include "DataTypes.h"

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
          escript::DataTypes::ShapeType givenShape(&dimensions[0],&dimensions[rank]);
          return (temp->getDataPointShape()==givenShape);
     }
  }
}

int  getNumDataPointsPerSample(struct escriptDataC* data) 
{
  if (data == (struct escriptDataC*)0) {
       return 0;
  } else {
     escript::Data* temp=(escript::Data*)(data->m_dataPtr);
     if (temp->isEmpty()) {
        return 0;
     } else {
          return (temp->getNumDataPointsPerSample());
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
        const escript::DataTypes::ShapeType view=temp->getDataPointShape();
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
        return temp->actsExpanded();
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

// The unusual (for me) ordering of __const here is because I'm not sure
// whether gcc would try to interpret __const as a function attribute rather than
// a modifier on the return value. Putting it here should remove any ambiguity
// I have used const rather than __const in the cpp because only c++ will be reading the cpp.
double const* getSampleDataRO(struct escriptDataC* data, int sampleNo, void* buffer)
{
  if (data == (struct escriptDataC*)0) {
       return NULL;
  } else {
      escript::Data* temp=(escript::Data*)(data->m_dataPtr);
     if (temp->isEmpty()) {
        return NULL;
     } else {
        return temp->getSampleDataRO(sampleNo,reinterpret_cast<escript::DataTypes::ValueType*>(buffer));
     }
  }
}

double* getSampleDataRW(struct escriptDataC* data, int sampleNo)
{
  if (data == (struct escriptDataC*)0) {
       return NULL;
  } else {
      escript::Data* temp=(escript::Data*)(data->m_dataPtr);
     if (temp->isEmpty()) {
        return NULL;
     } else {
        return temp->getSampleDataRW(sampleNo);
     }
  }
}

const double* getSampleDataROFast(struct escriptDataC* data, int sampleNo, void* buffer)
{
  escript::Data* temp=(escript::Data*)(data->m_dataPtr);
  return temp->getSampleDataRO(sampleNo, reinterpret_cast<escript::DataTypes::ValueType*>(buffer));
}

double* getSampleDataRWFast(struct escriptDataC* data, int sampleNo)
{
  escript::Data* temp=(escript::Data*)(data->m_dataPtr);
  return temp->getSampleDataRW(sampleNo);
}

void* allocSampleBuffer(escriptDataC* data)
{
  if (data == (struct escriptDataC*)0) {
     return NULL;
  } else {
     escript::Data* temp=(escript::Data*)(data->m_dataPtr);
     return temp->allocSampleBuffer();
  }
}

// Not going to the c++ member for this because I don't need an instance to do this
void freeSampleBuffer(void* buffer)
{
  if (buffer!=NULL)
  {
    delete (reinterpret_cast<escript::DataTypes::ValueType*>(buffer));
  }
}

void requireWrite(escriptDataC* data)
{
  if (data == (struct escriptDataC*)0) {
       return;
  } else {
      ((escript::Data*)(data->m_dataPtr))->requireWrite();
  }
}
