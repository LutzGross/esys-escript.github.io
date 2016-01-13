
/*****************************************************************************
*
* Copyright (c) 2003-2014 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/


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
   This domain supports a single type of FunctionSpace for which canTag is true.
   This compromise is needed to allow the default contructor of DataTagged to 
   have a FunctionSpace which supports tagging.
   See notes on the borrowListOfTagsInUse() method.
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
  \brief get the communicator for this domain.
  Returns 0 on non-MPI builds
  Routine must be implemented by the DomainAdapter
  */
  ESCRIPT_DLL_API
  virtual
#ifdef ESYS_MPI
  MPI_Comm
#else
  unsigned int
#endif
  getMPIComm() const
  {
#ifdef ESYS_MPI
    return MPI_COMM_WORLD;
#else
    return -1;
#endif    
  }
  
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
     NullDomain only has one FunctionSpace so this makes target a shallow copy of source.
  */
  ESCRIPT_DLL_API
  virtual void interpolateOnDomain(escript::Data& target,const escript::Data& source) const;
  ESCRIPT_DLL_API
  virtual bool probeInterpolationOnDomain(int functionSpaceType_source,int functionSpaceType_target) const;

  /**
     \brief
     Interpolates data given on source onto target where source and target are given on different domains.
     We do not permit interpolation into the NullDomain so this method always throws.
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
  virtual int getDiracDeltaFunctionsCode() const;

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

  ESCRIPT_DLL_API
  bool supportsContactElements() const;
  
  ESCRIPT_DLL_API
  virtual escript::Data randomFill(const DataTypes::ShapeType& shape,
       const FunctionSpace& what, long seed, const boost::python::tuple& filter) const;     
  
 protected:

 private:
 
};

} // end of namespace

#endif
