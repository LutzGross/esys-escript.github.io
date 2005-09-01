// $Id$
/* 
 ******************************************************************************
 *                                                                            *
 *       COPYRIGHT  ACcESS 2005 -  All Rights Reserved                        *
 *                                                                            *
 * This software is the property of ACcESS. No part of this code              *
 * may be copied in any form or by any means without the expressed written    *
 * consent of ACcESS.  Copying, use or modification of this software          *
 * by any unauthorised person is illegal unless that person has a software    *
 * license agreement with ACcESS.                                             *
 *                                                                            *
 ******************************************************************************
*/
                                                                           
#if !defined bruce_Bruce_20050829_H
#define bruce_Bruce_20050829_H

#include "escript/Data/AbstractContinuousDomain.h"
#include "escript/Data/FunctionSpaceFactory.h"
#include "escript/Data/FunctionSpace.h"
#include "escript/Data/Data.h"

#include <string>

namespace bruce {

/**
   \brief
   Bruce implements the AbstractContinuousDomain
   interface for the Bruce library.

   Description:
   Bruce implements the AbstractContinuousDomain
   interface for the Bruce library.
*/

class Bruce : public escript::AbstractContinuousDomain {

 public:

  //
  // Codes for function space types supported
  static const int Nodes;
  static const int Elements;

  /**
     \brief
     Constructor for Bruce.

     Description:
     Constructor for Bruce.
  */
  Bruce();

  /**
     \brief
     Copy constructor.
  */
  Bruce(const Bruce& other);

  /**
     \brief
     Destructor for Bruce.
  */
  ~Bruce();

  /**
     \brief
     Return this as an AbstractContinuousDomain.
  */
  inline
  const AbstractContinuousDomain&
  asAbstractContinuousDomain() const 
  {
     return *(static_cast<const AbstractContinuousDomain*>(this));
  }

  /**
     \brief
     Return this as an AbstractDomain.
  */
  inline
  const AbstractDomain&
  asAbstractDomain() const 
  {
     return *(static_cast<const AbstractDomain*>(this));
  }

  /**
     \brief
     Returns true if the given integer is a valid function space type
     for this domain.
  */
  virtual
  bool
  isValidFunctionSpaceType(int functionSpaceType) const;

  /**
     \brief
     Return the spatial dimension of the mesh.
  */
  virtual
  int
  getDim() const;

  /**
     \brief
     Return the number of data points per sample, and the number of samples
     needed to represent data on parts of the mesh.
  */
  virtual
  std::pair<int,int>
  getDataShape(int functionSpaceCode) const;

  /**
     \brief
     Returns the locations in the domain of the FEM nodes.
  */
  virtual
  escript::Data
  getX() const;

  /**
     \brief
     Returns the element size.
  */
  virtual
  escript::Data
  getSize() const;

  /**
     \brief
     Comparison operators.
  */
  virtual bool operator==(const AbstractDomain& other) const;
  virtual bool operator!=(const AbstractDomain& other) const;

 protected:

 private:

};

} // end of namespace

#endif
