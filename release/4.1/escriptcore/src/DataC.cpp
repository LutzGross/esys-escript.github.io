
/*****************************************************************************
*
* Copyright (c) 2003-2015 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/

#define ESNEEDPYTHON
#include "esysUtils/first.h"

#include "DataC.h"

#include "Data.h"
#include "DataTypes.h"

int getFunctionSpaceType(const escript::Data* data) 
{
  return data->getFunctionSpace().getTypeCode();
}


int isDataPointShapeEqual(const escript::Data* data, int rank, const int* dimensions)
{
  if (data == 0) {
       return 1;
  } else {
     return data->isDataPointShapeEqual(rank, dimensions);
  }
}

int getNumDataPointsPerSample(const escript::Data* data) 
{
  if (data == 0) {
       return 0;
  } else {
     if (data->isEmpty()) {
        return 0;
     } else {
          return (data->getNumDataPointsPerSample());
     }
  }
}

int numSamplesEqual(const escript::Data* data, int numDataPointsPerSample,
                    dim_t numSamples)
{
  if (data == 0) {
     return 1;
  } else {
     return data->numSamplesEqual(numDataPointsPerSample, numSamples);
  }
}

int getDataPointRank(const escript::Data* data)
{
  if (data == (const escript::Data*)0) {
       return 0;
  } else {
       return data->getDataPointRank();
  }
}

int getDataPointShape(const escript::Data* data,int i)
{
  if (data == 0) {
       return 0;
  } else {
     int rank = data->getDataPointRank();
     if (i<0 || i>=rank) {
        return 1;
     } else {
        const escript::DataTypes::ShapeType& view=data->getDataPointShape();
        return view[i];
     }
  }
}

int getDataPointSize(const escript::Data* data)
{
  return data->getDataPointSize();
}

int getLength(const escript::Data* data)
{
  return data->getLength();
}

int isExpanded(const escript::Data* data)
{
  if (data == 0) {
       return false;
  } else {
     if (data->isEmpty()) {
        return false;
     } else {
        return data->actsExpanded();
     }
  }
}

int isEmpty(const escript::Data* data) 
{
  if (data == 0) {
       return true;
  } else {
      return data->isEmpty();
  }
}

double const* getSampleDataRO(const escript::Data* data, int sampleNo)
{
  if (data == 0) {
       return NULL;
  } else {
     if (data->isEmpty()) {
        return NULL;
     } else {
        return data->getSampleDataRO(sampleNo);
     }
  }
}

double* getSampleDataRW(escript::Data* data, int sampleNo)
{
  if (data == 0) {
       return NULL;
  } else {
     if (data->isEmpty()) {
        return NULL;
     } else {
        return data->getSampleDataRW(sampleNo);
     }
  }
}

const double* getSampleDataROFast(const escript::Data* data, int sampleNo)
{
  return data->getSampleDataRO(sampleNo);
}

double* getSampleDataRWFast(escript::Data* data, int sampleNo)
{
  return data->getSampleDataRW(sampleNo);
}

double* getDataRW(escript::Data* data)
{
  
  if (data->getNumSamples()>0)
  {
     requireWrite(data);
     return getSampleDataRWFast(data,0);
  }
  return 0;
}


void requireWrite(escript::Data* data)
{
  if (data == 0) {
       return;
  } else {
      data->requireWrite();
  }
}
