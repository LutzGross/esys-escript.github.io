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

#include <iostream>

#include "escript/Data/DataVector.h" 

using namespace std;

namespace escript {

DataVector::DataVector() {
}

DataVector::DataVector(const DataVector& other) {
  m_data=other.m_data;
}

DataVector::DataVector(ValueType::size_type size, ValueType::value_type val) {
  resize(size, val);
}

DataVector::~DataVector() {
}

} // end of namespace
