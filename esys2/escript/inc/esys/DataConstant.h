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

#ifndef DATACONSTANT_H
#define DATACONSTANT_H

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
   Data

   Class Limitations:
   None

   Class Conditions of Use:
   None

   Throws:
   None

*/
class DataConstant:public DataAbstract {
  public:
  /**
     @memo
     Constructor which creates a Data object with only a constant value
     
     @param constValue Input - The constant value
  */
  DataConstant(double constValue);

 private:

};

#endif






