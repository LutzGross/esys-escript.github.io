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
#include "CompareFuncs.h"

/*
 * Return true if all of the calls match the expected values.
 * Return false if results don't match.
 */

int compareTypeCode(struct escriptDataC* data, int typeResult) 
{
  int result=0;
  result+=(getFunctionSpaceType(data)==typeResult);
  return result;
}

int compareNumSamples(struct escriptDataC* data, int numDataPointsPerSample, int numSamples)
{
  return numSamplesEqual(data,numDataPointsPerSample,numSamples);
}

int compareIsExpanded(struct escriptDataC* data, int expanded)
{
  return (isExpanded(data)==expanded);
}

int compareIsEmpty(struct escriptDataC* data, int empty)
{
  return (isEmpty(data)==empty);
}

/*
int comparePointShape(struct escriptDataC* data, int rank, int* dimensions)
{
  return pointShapeEqual(data,rank,dimensions);
}

int compareSampleDataWrite(struct escriptDataC* data, int sampleNo, double* sampleData)
{
  return (getSampleDataWrite(data,sampleNo)==sampleData);
}

int compareSampleDataRead(struct escriptDataC* data, int tag, double* sampleData)
{
  int sampleNo=0;
  return (getSampleDataRead(data,tag,sampleNo)==sampleData);
}
*/
