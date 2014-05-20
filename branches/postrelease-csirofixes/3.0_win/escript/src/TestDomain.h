
/*******************************************************
*
* Copyright (c) 2003-2009 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


#if !defined escript_TestDomain_20090618_H
#define escript_TestDomain_20090618_H
#include "system_dep.h"

#include "AbstractDomain.h"
#include "FunctionSpace.h"

#include <string>

namespace escript {

/**
   \brief
   (Testing use only) Provides a domain to wrap a collection of values.

   This domain provides more functionality than NullDomain in that it can 
   store varying numbers of samples and points per sample.
   
   It currently supports a single function space which does not support tagging.
   No effort has been made to make this work with MPI

   \warning This class exists to support testing and should not be used
   as a general domain without ensuring that it works the way you expect.

*/

class TestDomain : public AbstractDomain {

 public:

  /**
     \brief
     Default constructor for TestDomain.

     Description:
     Default constructor for TestDomain.

  */
  ESCRIPT_DLL_API
  TestDomain(int pointspersample, int numsamples);

  ESCRIPT_DLL_API
  ~TestDomain();

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
     Return a description for the given function space type code.
  */
  ESCRIPT_DLL_API
  virtual std::string functionSpaceTypeAsString(int functionSpaceType) const;

  /**
     \brief
     Interpolates data given on source onto target where source and target have to be given on the same domain.
     TestDomain only has one FunctionSpace so this makes target a shallow copy of source.
  */
  ESCRIPT_DLL_API
  virtual void interpolateOnDomain(escript::Data& target,const escript::Data& source) const;
  ESCRIPT_DLL_API
  virtual bool probeInterpolationOnDomain(int functionSpaceType_source,int functionSpaceType_target) const;

  /**
     \brief
     Interpolates data given on source onto target where source and target are given on different domains.
     We do not permit interpolation into the TestDomain so this method always throws.
  */
  ESCRIPT_DLL_API
  virtual void interpolateACross(escript::Data& target, const escript::Data& source) const;
  ESCRIPT_DLL_API
  virtual bool probeInterpolationACross(int functionSpaceType_source,const AbstractDomain& targetDomain, int functionSpaceType_target) const;


  /**
     \brief
     Return a continuous FunctionSpace.
  */
  ESCRIPT_DLL_API
  virtual int getDefaultCode() const;


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
  virtual const int* borrowSampleReferenceIDs(int functionSpaceType) const;

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

  /**
     \brief Checks if this domain allows tags for the specified functionSpaceCode.
  */
  ESCRIPT_DLL_API
  virtual
  bool canTag(int functionSpaceCode) const;

  /**
      \brief
          return the number of tags in use.
      For this class the answer is always 1(the default tag).
  */
  ESCRIPT_DLL_API
  virtual int getNumberOfTagsInUse(int functionSpaceCode) const;

  /**
     \brief  returns a pointer to an array with the tags used.
     For this class the answer will always be {0} 
  */
  ESCRIPT_DLL_API 
  virtual const int* borrowListOfTagsInUse(int functionSpaceCode) const;

 protected:

 private:
  int m_samples;	// number of samples
  int m_dpps;		// data points per sample
  int* m_samplerefids;	// sample reference ids
};

ESCRIPT_DLL_API
FunctionSpace
getTestDomainFunctionSpace(int dpps, int samples);

} // end of namespace

#endif
