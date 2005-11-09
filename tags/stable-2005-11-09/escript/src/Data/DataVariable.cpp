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

#include "escript/Data/DataVariable.h" 

using namespace std;
using namespace escript;

namespace escript {

DataVariable::DataVariable() :
  op(nullop),
  leftArg(0),
  rightArg(0),
  opBuffer(0)
{
  // cout << "DataVariable default constructor." << endl;
}

DataVariable::DataVariable(Data* data) :
  op(idop),
  leftArg(data),
  rightArg(0),
  opBuffer(0)
{

  int numSamples = leftArg->getNumSamples();
  int numValsPerSample = leftArg->getNumDataPointsPerSample() * leftArg->getDataPointSize();
  int bufSize = numSamples * numValsPerSample;

  // cout << "DataVariable constructor. " << bufSize << endl;
  opBuffer = new double[bufSize];
  for (int i=0; i<bufSize; i++) {
    opBuffer[i] = 0.0;
  }
}

DataVariable::~DataVariable()
{
  // cout << "DataVariable destructor";

  if (opBuffer!=0) {
    // cout << ": freeing buffer";
    delete[] opBuffer;
  }

  // cout << "." << endl;
}

Data
DataVariable::evaluate()
{
  // cout << "DataVariable evaluate." << endl;

  if (op==idop) {

    return *leftArg;

  } else {

    Data expTemp(0,leftArg->getDataPointShape(),leftArg->getFunctionSpace(),true);

    int numSamples = leftArg->getNumSamples();
    int numValsPerSample = leftArg->getNumDataPointsPerSample() * leftArg->getDataPointSize();

    for (int sampleNo=0; sampleNo<numSamples; sampleNo++) {

      double* sampleBuf = evaluate_samp(sampleNo);
      double* sampleLoc = expTemp.getSampleData(sampleNo);

      for (int i=0; i<numValsPerSample; i++) {
        sampleLoc[i] = sampleBuf[i];
      }

    }

    return expTemp;

  }

  throw DataException("DataVariable: nothing to evaluate");
  return *leftArg;

}

double*
DataVariable::evaluate_samp(int sampleNo)
{

  // cout << "DataVariable(" << sampleNo << ") evaluate." << endl;

  if (op==idop) {

    return leftArg->getSampleData(sampleNo);

  } else {

    int numValsPerSample = leftArg->getNumDataPointsPerSample() * leftArg->getDataPointSize();
    int offset = sampleNo * numValsPerSample;

    double* leftBuff = leftArg->getSampleData(sampleNo);
    double* rightBuff = rightArg->evaluate_samp(sampleNo);

    if (op == sumop) {
      for (int i=0; i<numValsPerSample; i++) {
        opBuffer[offset+i] = leftBuff[i] + rightBuff[i];
      }
    }

    if (op == diffop) {
      for (int i=0; i<numValsPerSample; i++) {
        opBuffer[offset+i] = leftBuff[i] - rightBuff[i];
      }
    }

    return opBuffer+offset;

  }

  throw DataException("DataVariable: nothing to evaluate(i)");
  return 0;

}

void
DataVariable::sum(DataVariable* right)
{
  // cout << "DataVariable sum." << endl;
  if (leftArg != 0) {
    if (leftArg->getNumDataPointsPerSample() == right->leftArg->getNumDataPointsPerSample() &&
        leftArg->getDataPointSize() == right->leftArg->getDataPointSize() &&
        leftArg->getNumSamples() == right->leftArg->getNumSamples()) {
      op=sumop;
      rightArg=right;
    } else {
      throw DataException("DataVariable: incompatible operands");
    }
  } else {
    throw DataException("DataVariable: cannot apply op to empty object.");
  }
}

void
DataVariable::diff(DataVariable* right)
{
  cout << "DataVariable diff." << endl;
  if (leftArg != 0) {
    if (leftArg->getNumDataPointsPerSample() == right->leftArg->getNumDataPointsPerSample() &&
        leftArg->getDataPointSize() == right->leftArg->getDataPointSize() &&
        leftArg->getNumSamples() == right->leftArg->getNumSamples()) {
      op=diffop;
      rightArg=right;
    } else {
      throw DataException("DataVariable: incompatible operands");
    }
  } else {
    throw DataException("DataVariable: cannot apply op to empty object.");
  }
}

}  // end of namespace
