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

#ifndef SERIALACCESSOR_H
#define SERIALACCESSOR_H

#include "esys/AbstractAccessor.h"

#include <string>
#include <Python.h>
#include <boost/python/object.hpp>
#include <boost/python/list.hpp>
#include <boost/python/numeric.hpp>

/**
   @memo
   Data

   @version 1.0.0 

   @doc

   Class Description:
   Base class for classes that provide a consistant interface to 
   an underlying vector of data.

   Class Limitations:
   None

   Class Conditions of Use:
   None

   Throws:
   None

*/
class SerialAccessor:public AbstractAccessor {
  public:
  /**
     @memo
     Constructor which sets the shape of the point data.

  */
  SerialAccessor(const std::vector<int>& pointDataShape, int blockSize, int numBlocks);
  /**
     @memo
     Method to copy from pointData into values

  */
  virtual void copy(const boost::python::numeric::array& pointData,ValueContainerType& values) const;
  /**
     @memo
     Method to perform the given operation on each value
  */
  virtual void evaluate(ValueContainerType& values, void (*operation)(double)) const;
 private:

};

#endif
