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

#ifndef DATAFACTORY_H
#define DATAFACTORY_H

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
   DataFactory is a class which can not be instantiated. It is simply
   a holder for various static factory methods that create Data objects.

   Class Limitations:
   None

   Class Conditions of Use:
   None

   Throws:
   None

*/
class DataFactory {
  public:
  /**
     @memo
     Factory method which creates a Data object with only a constant value
     
     @param constValue Input - The constant value
  */
  static DataAbstract* createData(double constValue);
  /**
     @memo
     Factory method which copies data from a python numarray.
     
     @param data Input - The input data.
  */
  static DataAbstract* createData(const boost::python::numeric::array& data);

 private:
  /**
     @memo
     Don't allow construction.
     
     @param exceptionReason Input - Exception message.
  */
  DataFactory();

};

#endif






