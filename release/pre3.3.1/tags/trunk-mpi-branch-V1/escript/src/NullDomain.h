/* 
 ************************************************************
 *          Copyright 2006 by ACcESS MNRF                   *
 *                                                          *
 *              http://www.access.edu.au                    *
 *       Primary Business: Queensland, Australia            *
 *  Licensed under the Open Software License version 3.0    *
 *     http://www.opensource.org/licenses/osl-3.0.php       *
 *                                                          *
 ************************************************************
*/

#if !defined escript_NullDomain_20040604_H
#define escript_NullDomain_20040604_H
#include "system_dep.h"

#include "AbstractDomain.h"

#include <string>

namespace escript {

/**
   \brief
   NullDomain provides a null value for domain. Needed for the construction
   of a default FunctionSpace.

   Description:
   NullDomain provides a null value for domain. Needed for the construction
   of a default FunctionSpace. Inherits from AbstractDomain and overrides its
   methods.
*/

class NullDomain : public AbstractDomain {

 public:

  /**
     \brief
     Default constructor for NullDomain.

     Description:
     Default constructor for NullDomain.

  */
  ESCRIPT_DLL_API
  NullDomain();

  /**
     \brief
     Returns true if the given integer is a valid function space type
     for this domain.
  */
  ESCRIPT_DLL_API
  virtual bool isValidFunctionSpaceType(int functionSpaceType) const;

  /**
     \brief
     Return a description for this domain.
  */
  ESCRIPT_DLL_API
  virtual std::string getDescription() const;

  /**
     \brief
     Return a continuous FunctionSpace.
  */
  ESCRIPT_DLL_API
  virtual int getContinuousFunctionCode() const;

  /**
     \brief
     Return a function FunctionSpace.
  */
  ESCRIPT_DLL_API
  virtual int getFunctionCode() const;

  /**
     \brief
     Return a function on boundary FunctionSpace.
  */
  ESCRIPT_DLL_API
  virtual int getFunctionOnBoundaryCode() const;

  /**
     \brief
     Return a FunctionSpace.
  */
  ESCRIPT_DLL_API
  virtual int getFunctionOnContactZeroCode() const;

  /**
     \brief
     Return a FunctionSpace.
  */
  ESCRIPT_DLL_API
  virtual int getFunctionOnContactOneCode() const;

  /**
     \brief
     Return a FunctionSpace.
  */
  ESCRIPT_DLL_API
  virtual int getSolutionCode() const;

  /**
     \brief
     Return a FunctionSpace.
  */
  ESCRIPT_DLL_API
  virtual int getReducedSolutionCode() const;

  /**
     \brief
     Return a FunctionSpace.
  */
  ESCRIPT_DLL_API
  virtual int getDiracDeltaFunctionCode() const;

  /**
     \brief
     Return the number of data points per sample, and the number of samples as a pair.
     \param functionSpaceCode Input - Code for the function space type.
     \return pair, first - number of data points per sample, second - number of samples
  */
  ESCRIPT_DLL_API
  virtual std::pair<int,int> getDataShape(int functionSpaceCode) const;

  /**
     \brief
     Return the tag key for the given sample number.
     \param functionSpaceType Input - The function space type.
     \param sampleNo Input - The sample number.
  */
  ESCRIPT_DLL_API
  virtual int getTagFromSampleNo(int functionSpaceType, int sampleNo) const;

  /**
     \brief
     Return a borrowed pointer to the sample reference number id list
     \param functionSpaceType Input - The function space type.
  */
  ESCRIPT_DLL_API
  virtual int* borrowSampleReferenceIDs(int functionSpaceType) const;

  /**
     \brief
  */
  ESCRIPT_DLL_API
  virtual int getDim() const;

  /**
     \brief
     Return true if given domains are equal.
  */
  ESCRIPT_DLL_API
  virtual bool operator==(const AbstractDomain& other) const;
  ESCRIPT_DLL_API
  virtual bool operator!=(const AbstractDomain& other) const;

 protected:

 private:

};

} // end of namespace

#endif
