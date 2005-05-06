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

#if !defined escript_DataVariable_20050422_H
#define escript_DataVariable_20050422_H

#include "escript/Data/FunctionSpace.h"
#include "escript/Data/DataException.h"
#include "escript/Data/DataArrayView.h"
#include "escript/Data/DataVector.h"
#include "escript/Data/Data.h"

namespace escript {

/**
   \brief
   Give a short description of what DataVariable does.

   Description:
   Give a detailed description of DataVariable

   Template Parameters:
   For templates describe any conditions that the parameters used in the
   template must satisfy
*/

class DataVariable {

 public:

  typedef enum {id, add} OpCode;

  /**
     \brief
     Default constructor for DataVariable

     Description:
     Default constructor for DataVariable

     Preconditions:
     Describe any preconditions

     Throws:
     Describe any exceptions thrown
  */
  DataVariable();

  /**
     \brief
     Constructor for DataVariable

     Description:
     Constructor for DataVariable

     Preconditions:
     Describe any preconditions

     Throws:
     Describe any exceptions thrown
  */
  DataVariable(Data& data);

  /**
     \brief
     Destructor for DataVariable

     Description:
     Destructor for DataVariable

     Preconditions:
     Describe any preconditions

     Throws:
     Describe any exceptions thrown
  */
  ~DataVariable();

  /**
     \brief
     Evaluator for DataVariable

     Description:
     Evaluator for DataVariable

     Preconditions:
     Describe any preconditions

     Throws:
     Describe any exceptions thrown
  */
  void evaluate();

 protected:

 private:

  FunctionSpace* functionSpace;

  DataArrayView::ShapeType shape;

  DataArrayView::ValueType left;

  DataVariable* right;

  OpCode op;

};

DataVariable operator+(const DataVariable& left, const DataVariable& right);

} // end of namespace

#endif
