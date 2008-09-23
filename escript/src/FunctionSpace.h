
/* $Id$ */

/*******************************************************
 *
 *           Copyright 2003-2007 by ACceSS MNRF
 *       Copyright 2007 by University of Queensland
 *
 *                http://esscc.uq.edu.au
 *        Primary Business: Queensland, Australia
 *  Licensed under the Open Software License version 3.0
 *     http://www.opensource.org/licenses/osl-3.0.php
 *
 *******************************************************/

#if !defined escript_FunctionSpace_20040323_H
#define escript_FunctionSpace_20040323_H
#include "system_dep.h"

#include "AbstractDomain.h"
#include "NullDomain.h"

#include <string>
#include <list>

// This stops the friend declaration below from looking like a
// declaration of escript::FunctionSpaceTestCase.

class FunctionSpaceTestCase;

namespace escript {

//
// Forward declaration for class Data.
class Data;

/**
   \brief
   Give a short description of what FunctionSpace does.

   Description:
   Give a detailed description of FunctionSpace.

   Template Parameters:
   For templates describe any conditions that the parameters used in the
   template must satisfy.
*/

class FunctionSpace {

  // These are using operator=()
  friend class AbstractSystemMatrix;
  friend class AbstractTransportProblem;
  friend class ::FunctionSpaceTestCase;

 public:
  /**
     \brief
     Default constructor for FunctionSpace.

     Description:
     Default constructor for FunctionSpace
     Generates a function space with a null domain.

     Preconditions:
     Describe any preconditions.

     Throws:
     Describe any exceptions thrown.
  */
  ESCRIPT_DLL_API
  FunctionSpace();

  /**
     \brief
     Constructor for FunctionSpace.

     Description:
     Constructor for FunctionSpace.

     NOTE: The FunctionSpace class relies on the domain existing
     for the lifetime of the FunctionSpace object. ie: domain must
     be an externally managed object (!).
  */
  ESCRIPT_DLL_API
  FunctionSpace(const AbstractDomain& domain,
                int functionSpaceType);

  /**
    \brief
    Return the function space type code.

    Note: The meaning of the code depends on the domain object the FunctionSpace is built on.
  */
  ESCRIPT_DLL_API
  int
  getTypeCode() const;

  /**
   \brief
   Return the function space domain.
  */
  ESCRIPT_DLL_API
  const
  AbstractDomain&
  getDomain() const;

 /**
 \brief Checks if this functionspace support tags
 */
  ESCRIPT_DLL_API
  bool
  canTag() const;


  /**
     \brief assigns new tag newTag to all samples with a positive
     value of mask for any its sample point.

  */
  ESCRIPT_DLL_API
  void setTags(const int newTag, const escript::Data& mask) const;


  /**
   \brief
   Return the shape of the data needed to represent the function space.
  */
  ESCRIPT_DLL_API
  std::pair<int,int>
  getDataShape() const;

  /**
   \brief
   Comparison operator.
   Return true if function spaces are equal.
   ie: Same domain and same function space type.
  */
  ESCRIPT_DLL_API
  bool
  operator==(const FunctionSpace& other) const;

  ESCRIPT_DLL_API
  bool
  operator!=(const FunctionSpace& other) const;

  /**
   \brief
   Return a text description of the function space.
   NOTE: In the python interface, the reference return value policy is
   copy const reference. This mens that if the python string value
   is assigned to a python variable, it's value (python) value will not
   change, even if operator=() changes this Function space object.
   Consequently, it is possible for python to keep a wrong value for this
   string!
   Worse still is that operator=() is exposed to python in escriptcpp.cpp.
  */
  ESCRIPT_DLL_API
  std::string
  toString() const;

   //#define DEBUG_PY_STRINGS

#ifdef DEBUG_PY_STRINGS
  /**
   \brief
   Return a text description of the function space
   as a python string.
   NOTE: This code was used to debug a conversion of
   std::string to python string problem on windows. 
   An alternative approach was sought.
  */
  ESCRIPT_DLL_API
  PyObject *
  toPyString() const;
#endif

  /**
   \brief
   Return the tag associated with the given sample number.
  */
  ESCRIPT_DLL_API
  int
  getTagFromSampleNo(int sampleNo) const;

  /**
   \brief
   Return the tag associated with the given data-point number.
  */
  ESCRIPT_DLL_API
  int
  getTagFromDataPointNo(int dataPointNo) const;

  /**
   \brief
   Return the reference number associated with the given data-point number.
  */
  ESCRIPT_DLL_API
  int getReferenceIDFromDataPointNo(int dataPointNo) const;

  /**
   \brief
   Return the reference number associated with the given sample number.
   This function is not efficient. It is better to first call 
   borrowSampleReferenceIDs and then when iterating over sampleNo to use sampleNo as an offset.
  */
  ESCRIPT_DLL_API
  inline
  int
  getReferenceIDOfSample(int sampleNo) const
  {
      return borrowSampleReferenceIDs()[sampleNo];
  }
  /**
   \brief
   Return a borrowed reference to the list of sample reference IDs
  */
  ESCRIPT_DLL_API
  int*
  borrowSampleReferenceIDs() const;

  /**
   \brief
   Return the spatial locations of the data points.
  */
  ESCRIPT_DLL_API
  escript::Data
  getX() const;
 
  /**
   \brief
   Return the surface normal field.
  */
  ESCRIPT_DLL_API
  escript::Data
  getNormal() const;

  /**
   \brief
   Return the sample size (e.g. the diameter of elements, radius of particles).
  */
  ESCRIPT_DLL_API
  escript::Data
  getSize() const;

  /**
   \brief
   Return the number of samples.
  */
  ESCRIPT_DLL_API
  inline
  int
  getNumSamples() const {
     return getDataShape().second;
  }

  /**
   \brief
   Return the number of data points per sample.
  */
  ESCRIPT_DLL_API
  inline
  int
  getNumDPPSample() const {
     return getNumDataPointsPerSample();
  }

  ESCRIPT_DLL_API
  inline
  int
  getNumDataPointsPerSample() const {
     return getDataShape().first;
  }

  /**
   \brief
   Return the spatial dimension of the underlying domain.
  */
  ESCRIPT_DLL_API
  inline
  int
  getDim() const {
      return getDomain().getDim();
  }
  /**
   \brief
   Returns a list of the tags used in this function space
  */
  ESCRIPT_DLL_API
  boost::python::list
  getListOfTags() const;
  /**
   \brief
   Returns an stl::list of the tags used in this function space
  */
  ESCRIPT_DLL_API
  std::list<int>
  getListOfTagsSTL() const;

  /**
     \brief
        return the number of tags in use and a pointer to an array with the number of tags in use
  */
  ESCRIPT_DLL_API
  virtual int getNumberOfTagsInUse() const;

  ESCRIPT_DLL_API
  virtual int* borrowListOfTagsInUse() const;



 protected:

 private:
  /**
   \brief
   Assignment operator. 
   NOTE: Assignment copies the domain object pointer
   as this object is managed externally to this class.
   Also, breaks the non-mutability of FunctionSpace
   assumed by toString().
  */
  ESCRIPT_DLL_API
  FunctionSpace&
  operator=(const FunctionSpace& other);

  //
  // static null domain value
  static const NullDomain nullDomainValue;

  //
  // function space domain
  const AbstractDomain*  m_domain;

  //
  // function space type code.
  int m_functionSpaceType;

};

} // end of namespace

#endif
