/*=============================================================================
 
 ******************************************************************************
 *                                                                            *
 *       COPYRIGHT  ???????  -  All Rights Reserved                           *
 *                                                                            *
 * This software is the property of ??????????????.  No part of this code     *
 * may be copied in any form or by any means without the expressed written    *
 * consent of ???????????.  Copying, use or modification of this software     *
 * by any unauthorised person is illegal unless that                          *
 * person has a software license agreement with ?????????????.                *
 *                                                                            *
 ******************************************************************************
 
*********************************************************************************/

#ifndef DATAGENERAL_H
#define DATAGENERAL_H

#include <string>
#include <Python.h>
#include <boost/python/object.hpp>
#include <boost/python/list.hpp>
#include <boost/python/numeric.hpp>
#include "esys/DataAbstract.h"

/**
   @memo
   Data

   @version 1.0.0 

   @doc

   Class Description:
   The general case for the Data class abstraction.

   Class Limitations:
   None

   Class Conditions of Use:
   None

   Throws:
   None

*/
class DataGeneral {
  public:
  /**
     @memo
     Constructor which copies data from a python numarray.
     
     @param exceptionReason Input - Exception message.
  */
  DataGeneral(const boost::python::numeric::array& data);

 private:

};

#endif






