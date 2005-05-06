// $Id$
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
                                                                           
#if !defined  finley_MeshAdapter_20040526_H
#define finley_MeshAdapter_20040526_H

#include "escript/Data/AbstractContinuousDomain.h"
#include "escript/Data/Data.h"
#include "escript/Data/FunctionSpace.h"
extern "C" {
#include "finley/finleyC/Mesh.h"
}
#include "finley/CPPAdapter/SystemMatrixAdapter.h"
#include <boost/shared_ptr.hpp>
#include <boost/python/object.hpp>
#include <map>
#include <vector>
#include <string>

namespace finley {

/**
   \brief
   MeshAdapter implements the AbstractContinuousDomain
   interface for the Finley library.

   Description:
   MeshAdapter implements the AbstractContinuousDomain
   interface for the Finley library.
*/

class MeshAdapter:public escript::AbstractContinuousDomain {

 public:

  //
  // Codes for function space types supported
  static const int DegreesOfFreedom;
  static const int ReducedDegreesOfFreedom;
  static const int Nodes;
  static const int Elements;
  static const int FaceElements;
  static const int Points;
  static const int ContactElementsZero;
  static const int ContactElementsOne;

  /**
     \brief
     Constructor for MeshAdapter

     Description:
     Constructor for MeshAdapter. The pointer passed to MeshAdapter
     is deleted using a call to Finley_Mesh_deallocate in the
     MeshAdapter destructor.

     Throws:
     May throw an exception derived from EsysException

     \param finleyMesh Input - A pointer to the externally constructed 
                               finley mesh.The pointer passed to MeshAdapter
                               is deleted using a call to 
                               Finley_Mesh_deallocate in the MeshAdapter 
                               destructor.
  */
  MeshAdapter(Finley_Mesh* finleyMesh=0);
  /**
     \brief
     Copy constructor.
  */
  MeshAdapter(const MeshAdapter& in);
  /**
     \brief
     Destructor for MeshAdapter. As specified in the constructor
     this calls  Finley_Mesh_deallocate for the pointer given to the 
     constructor.
  */
  ~MeshAdapter();
  /**
     \brief
     return this as an AbstractContinuousDomain.
  */
  inline const AbstractContinuousDomain& asAbstractContinuousDomain() const 
  {
     return *(static_cast<const AbstractContinuousDomain*>(this));
  }

  /**
     \brief
     return this as an AbstractDomain.
  */
  inline const AbstractDomain& asAbstractDomain() const 
  {
     return *(static_cast<const AbstractDomain*>(this));
  }
  /**
     \brief
     Write the current mesh to a file with the given name.
     \param fileName Input - The name of the file to write to.
  */
  void write(const std::string& fileName) const;
  /**
     \brief
     return the pointer of the underlying finley mesh structure
  */
  Finley_Mesh* getFinley_Mesh() const;
  /**
     \brief
     Return the tag list indexed by sampleNo. 
     \param functionSpaceType Input
     \param tagList Output
     \param numTags Output
  */

  void getTagList(int functionSpaceType, int** tagList, 
			       int* numTags) const;
  /**
     \brief
     Return the reference number list indexed by sampleNo. 
     \param functionSpaceType Input
     \param referenceNoList Output
     \param numReferenceNo Output
  */
  void getReferenceNoList(int functionSpaceType, int** referenceNoList, 
			       int* numReferenceNo) const;

   /**
     \brief
     Return the tag key for the given sample number.
     \param functionSpaceType Input - The function space type.
     \param sampleNo Input - The sample number.
  */
  int getTagFromSampleNo(int functionSpaceType, int sampleNo) const;
  /**
     \brief
     Return the reference number of  the given sample number.
     \param functionSpaceType Input - The function space type.
     \param sampleNo Input - The sample number.
  */
  int getReferenceNoFromSampleNo(int functionSpaceType, int sampleNo) const;

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
     Return a description for the given function space type code
  */
  virtual std::string functionSpaceTypeAsString(int functionSpaceType) const;

  /**
     \brief
     Build the table of function space type names
  */
  void setFunctionSpaceTypeNames();
  /**
     \brief
     Return a continuous FunctionSpace code
  */
  virtual int getContinuousFunctionCode() const;
  /**
     \brief
     Return a functon FunctionSpace code
  */
  virtual int getFunctionCode() const;
  /**
     \brief
     Return a function on boundary FunctionSpace code
  */
  virtual int getFunctionOnBoundaryCode() const;
  /**
     \brief
     Return a FunctionOnContactZero code
  */
  virtual int getFunctionOnContactZeroCode() const;
  /**
     \brief
     Return a FunctionOnContactOne code
  */
  virtual int getFunctionOnContactOneCode() const;
  /**
     \brief
     Return a Solution code
  */
  virtual int getSolutionCode() const;
  /**
     \brief
     Return a ReducedSolution code
  */
  virtual int getReducedSolutionCode() const;
  /**
     \brief
     Return a DiracDeltaFunction code
  */
  virtual int getDiracDeltaFunctionCode() const;
  //
  //
  typedef std::map<int, std::string> FunctionSpaceNamesMapType;
  /**
     \brief
  */
  virtual int getDim() const;
 /**
     \brief
     Return the number of data points per sample, and the number of samples as a pair.
     \param functionSpaceCode Input -
  */
  virtual std::pair<int,int> getDataShape(int functionSpaceCode) const;

  /**
     \brief
     copies the location of data points into arg. The domain of arg has to match this.
     has to be implemented by the actual Domain adapter.
  */
  virtual void setToX(escript::Data& arg) const;
  /**
     \brief
     assigns new location to the domain
  */
  virtual void setNewX(const escript::Data& arg);
  /**
     \brief
     interpolates data given on source onto target where source and target have to be given on the same domain.
  */
  virtual void interpolateOnDomain(escript::Data& target,const escript::Data& source) const;
  virtual bool probeInterpolationOnDomain(int functionSpaceType_source,int functionSpaceType_target) const;
  /**
     \brief
     interpolates data given on source onto target where source and target are given on different domains.
     has to be implemented by the actual Domain adapter.
  */
  virtual void interpolateACross(escript::Data& target, const escript::Data& source) const;
  virtual bool probeInterpolationACross(int functionSpaceType_source,const AbstractDomain& targetDomain, int functionSpaceType_target) const;
  /**
     \brief
     copies the surface normals at data points into out.  The actual function space to be considered
     is defined by out. out has to be defined on this.
  */
  virtual void setToNormal(escript::Data& out) const;
  /**
     \brief
     copies the size of samples into out. The actual function space to be considered
     is defined by out. out has to be defined on this.
  */
  virtual void setToSize(escript::Data& out) const;

  /**
     \brief
     copies the gradient of arg into grad. The actual function space to be considered
     for the gradient is defined by grad. arg and grad have to be defined on this.
  */
  virtual void setToGradient(escript::Data& grad,const escript::Data& arg) const;

  /**
     \brief
     copies the integrals of the function defined by arg into integrals.
     arg has to be defined on this.
  */
  virtual void setToIntegrals(std::vector<double>& integrals,const escript::Data& arg) const;

 /**
     \brief
     return the identifier of the matrix type to be used for the global stiffness matrix when a particular solver, preconditioner
     and symmetric matrix is used.
     \param solver 
     \param symmetry 
  */
  virtual int getSystemMatrixTypeId(const int solver, const bool symmetry) const;

  /**
     \brief
     returns true if data on this domain and a function space of type functionSpaceCode has to 
     considered as cell centered data.
  */
  virtual bool isCellOriented(int functionSpaceCode) const;
  /**
     \brief
     saves data arg to an OpenDX input file.
     considered as cell centered data.
  */
  virtual void saveDX(const std::string& filename,const escript::Data& arg) const;
  /**
     \brief
     saves data arg to a VTK input file.
     considered as cell centered data.
  */
  virtual void saveVTK(const std::string& filename,const escript::Data& arg) const;
  /**
     \brief
     returns the function space representation of the type functionSpaceCode on this domain
     as a vtkObject.
  */
  // vtkObject createVtkObject(int functionSpaceCode) const;
  /**
     \brief
     adds a PDE onto the stiffness matrix mat and a rhs 
  */
  virtual void addPDEToSystem(
                     SystemMatrixAdapter& mat, escript::Data& rhs,
                     const escript::Data& A, const escript::Data& B, const escript::Data& C, 
                     const escript::Data& D, const escript::Data& X, const escript::Data& Y,
                     const escript::Data& d, const escript::Data& y,
                     const escript::Data& d_contact, const escript::Data& y_contact) const;
  /**
     \brief
     adds a PDE onto the stiffness matrix mat and a rhs 
  */
  virtual void addPDEToRHS(escript::Data& rhs,
                     const escript::Data& X, const escript::Data& Y,
                     const escript::Data& y, const escript::Data& y_contact) const;
  /**
     \brief
    creates a SystemMatrixAdapter stiffness matrix an initializes it with zeros:
  */
  SystemMatrixAdapter newSystemMatrix(
                      const int row_blocksize,
                      const escript::FunctionSpace& row_functionspace,
                      const int column_blocksize,
                      const escript::FunctionSpace& column_functionspace,
                      const int type) const;

  /**
     \brief returns locations in the FEM nodes
  */
  virtual escript::Data getX() const;
  /**
     \brief return boundary normals at the quadrature point on the face elements
  */
  virtual escript::Data getNormal() const;
  /**
     \brief returns the element size
  */
  virtual escript::Data getSize() const;

  virtual bool operator==(const AbstractDomain& other) const;
  virtual bool operator!=(const AbstractDomain& other) const;

 protected:

 private:
  //
  // pointer to the externally created finley mesh
  boost::shared_ptr<Finley_Mesh> m_finleyMesh;
 
  static FunctionSpaceNamesMapType m_functionSpaceTypeNames;

};

} // end of namespace
#endif
