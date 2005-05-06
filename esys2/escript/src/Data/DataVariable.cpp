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

DataVariable::DataVariable() {
}

DataVariable::DataVariable(Data& data) {
}

DataVariable::~DataVariable() {
}

void
DataVariable::evaluate() {
}

DataVariable operator+(const DataVariable& left, const DataVariable& right) {
  DataVariable d;
  return d;
}

}  // end of namespace
