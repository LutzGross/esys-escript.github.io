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

#include "esysUtils/esysExceptionTranslator.h" 

#include <iostream>

using namespace std;

namespace esysUtils {

  void esysExceptionTranslator(EsysException const& e) 
  {
    PyErr_SetString(PyExc_StandardError,e.what());
  }

}  // end of namespace
