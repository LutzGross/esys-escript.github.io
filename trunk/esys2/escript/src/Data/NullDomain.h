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
                                                                           
#if !defined  escript_NullDomain_20040604_H
#define escript_NullDomain_20040604_H

#include "escript/Data/AbstractDomain.h"

#include <string>

namespace escript {

/**
   \brief
   NullDomain provides a null value for domain. Needed for the construction
   of a default FunctionSpace.

   Description:
   NullDomain provides a null value for domain. Needed for the construction
   of a default FunctionSpace.
*/
class NullDomain : public AbstractDomain {

 public:

  /**
     \brief
     Default constructor for NullDomain

     Description:
     Default constructor for NullDomain

  */
  NullDomain();
  /**
     \brief
     Returns true if the given integer is a valid function space type
     for this domain.
  */
  virtual bool isValidFunctionSpaceType(int functionSpaceType) const;
  /**
     \brief
     Return a description for this domain
  */
  virtual std::string getDescription() const;
  /**
     \brief
     Return a continuous FunctionSpace
  */
  virtual int getContinuousFunctionCode() const;
  /**
     \brief
     Return a functon FunctionSpace
  */
  virtual int getFunctionCode() const;
  /**
     \brief
     Return a function on boundary FunctionSpace
  */
  virtual int getFunctionOnBoundaryCode() const;
  /**
     \brief
     Return a FunctionSpace
  */
  virtual int getFunctionOnContactZeroCode() const;
  /**
     \brief
     Return a FunctionSpace
  */
  virtual int getFunctionOnContactOneCode() const;
  /**
     \brief
     Return a FunctionSpace
  */
  virtual int getSolutionCode() const;
  /**
     \brief
     Return a FunctionSpace
  */
  virtual int getReducedSolutionCode() const;
  /**
     \brief
     Return a FunctionSpace
  */
  virtual int getDiracDeltaFunctionCode() const;
 /**
     \brief
     Return the number of data points per sample, and the number of samples as a pair.
     \param functionSpaceCode Input - Code for the function space type.
     \return pair, first - number of data points per sample, second - number of samples
  */
  virtual std::pair<int,int> getDataShape(int functionSpaceCode) const;
  /**
     \brief
     Return the tag key for the given sample number.
     \param functionSpaceType Input - The function space type.
     \param sampleNo Input - The sample number.
  */
  virtual int getTagFromSampleNo(int functionSpaceType, int sampleNo) const;
  /**
     \brief
  */
  virtual int getDim() const;
 protected:

 private:
};

} // end of namespace
#endif
