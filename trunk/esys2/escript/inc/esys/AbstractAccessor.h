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

#ifndef ABSTRACTACCESSOR_H
#define ABSTRACTACCESSOR_H

#include <string>
#include <vector>
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
class AbstractAccessor {
  public:

  typedef std::vector<double> ValueContainerType;

  /**
     @memo
     Constructor which sets the shape of the point data.

  */
  AbstractAccessor(const std::vector<int>& pointDataShape, int blockSize, int numBlocks);
  /**
     @memo
     Destructor
  */
  virtual ~AbstractAccessor();
  /**
     @memo
     Method to copy from pointData into values

  */
  virtual void copy(const boost::python::numeric::array& pointData,ValueContainerType& values) const=0;
  /**
     @memo
     Method to perform the given operation on each value
  */
  virtual void evaluate(ValueContainerType& values, void (*operation)(double)) const=0;
 private:

};

#endif
