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
#include <vector>

namespace bruce {

/**
   \brief
   Bruce implements a structured AbstractContinuousDomain.

   Description:
   Bruce implements a structured AbstractContinuousDomain.
*/

class Bruce : public escript::AbstractContinuousDomain {

 public:

  //
  // Codes for function space types supported
  static const int ContinuousFunction;
  static const int Function;

  //
  // Type of FunctionSpaceNamesMap
  typedef std::map<int, std::string> FunctionSpaceNamesMapType;

  //
  // Types for the dimension vectors
  typedef std::vector<double> DimVec;

  /**
     \brief
     Default constructor for Bruce.

     Description:
     Default constructor for Bruce.
     Creates a null Bruce object.
  */
  Bruce();

  /**
     \brief
     Constructor for Bruce.

     Description:
     Constructor for Bruce.

     The point "origin" specifies the location of the origin
     of the domain specified by this object. The dimensionality of this
     point determines the dimensionality of the space the domain occupies.

     The vectors v0,v1,v2 specify the axis in
     of the domain of this Bruce object. If v2 is an empty vector, this
     object is a two dimensional domain. If v1 is also an empty vector,
     this object is a one dimensional domain. If v0 is also an empty
     vector, this is a point domain.

     The integers n0,n1,n2 specify the dumber of data-points along each
     axis in the domain.
  */
  Bruce(DimVec v0, DimVec v1, DimVec v2,
        int n0, int n1, int n2,
        DimVec origin);

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
     Return a description for this domain.
  */
  virtual
  inline
  std::string
  getDescription() const
  {
    return "Bruce";
  }

  /**
     \brief
     Returns true if the given integer is a valid function space type
     for this domain.
  */
  virtual
  bool
  isValidFunctionSpaceType(int functionSpaceCode) const;

  /**
     \brief
     Return a description for the given function space type code.
  */
  virtual
  std::string
  functionSpaceTypeAsString(int functionSpaceCode) const;

  /**
     \brief
     Return a continuous FunctionSpace code.
  */
  virtual
  inline
  int
  getContinuousFunctionCode() const
  {
    return ContinuousFunction;
  }

  /**
     \brief
     Return a function FunctionSpace code.
  */
  virtual
  inline
  int
  getFunctionCode() const
  {
    return Function;
  }

  /**
     \brief
     Return the spatial dimension of the mesh.
  */
  virtual
  inline
  int
  getDim() const
  {
    return m_origin.size();
  }

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
     Copies the location of data points on the domain into out.
  */
  virtual
  void
  setToX(escript::Data& out) const;

  /**
     \brief
     Returns the element size.
  */
  virtual
  escript::Data
  getSize() const;

  /**
     \brief
     Copies the size of samples into out.
  */
  virtual
  void
  setToSize(escript::Data& out) const;

  /**
     \brief
     Comparison operators.
  */
  virtual bool operator==(const AbstractDomain& other) const;
  virtual bool operator!=(const AbstractDomain& other) const;

  /*
     \brief
     Return the tag key for the given sample number.
  */
  virtual
  int
  getTagFromSampleNo(int functionSpaceCode, int sampleNo) const;

  /**
     \brief
     Return the reference number of the given sample number.
  */
  virtual
  int
  getReferenceNoFromSampleNo(int functionSpaceCode, int sampleNo) const;


 protected:

  /**
     \brief
     Build the table of function space type names.
  */
  void
  setFunctionSpaceTypeNames();

  /**
     \brief
     Ensure the parameters supplied to the constructor are valid.
  */
  bool
  checkParameters();

  /**
     \brief
     Check if all components of vector are zero.
  */
  static
  bool
  isZero(DimVec vec);

 private:

  //
  // vectors describing axis of the domain
  DimVec m_v0, m_v1, m_v2;

  //
  // number of data points in each axial direction of the domain
  int m_n0, m_n1, m_n2;

  //
  // the coordinates of the origin of the domain
  DimVec m_origin;

  //
  // map from FunctionSpace codes to names
  static FunctionSpaceNamesMapType m_functionSpaceTypeNames;

};

} // end of namespace

#endif
