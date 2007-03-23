// $Id$
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
                                                                           
#if !defined finley_MeshAdapter_20040526_H
#define finley_MeshAdapter_20040526_H
#include "system_dep.h"

extern "C" {
#include "../Mesh.h"
#include "../Finley.h"
#include "../Assemble.h"
#include "paso/SystemMatrix.h"
}

#include "FinleyError.h"
#include "FinleyAdapterException.h"

#include "SystemMatrixAdapter.h"
#include "escript/AbstractContinuousDomain.h"
#include "escript/FunctionSpace.h"
#include "escript/FunctionSpaceFactory.h"

#include <boost/shared_ptr.hpp>
#include <boost/python/dict.hpp>
#include <boost/python/extract.hpp>

#include <map>
#include <vector>
#include <string>
#include <sstream>

//
// forward declarations
class Data;

//using namespace escript;

namespace finley {

struct null_deleter
{
  void operator()(void const *ptr) const
  {
  }
};


/**
   \brief
   MeshAdapter implements the AbstractContinuousDomain
   interface for the Finley library.

   Description:
   MeshAdapter implements the AbstractContinuousDomain
   interface for the Finley library.
*/

class MeshAdapter : public escript::AbstractContinuousDomain {

 public:

  //
  // Codes for function space types supported
  static const int DegreesOfFreedom;
  static const int ReducedDegreesOfFreedom;
  static const int Nodes;
  static const int Elements;
  static const int ReducedElements;
  static const int FaceElements;
  static const int ReducedFaceElements;
  static const int Points;
  static const int ContactElementsZero;
  static const int ReducedContactElementsZero;
  static const int ContactElementsOne;
  static const int ReducedContactElementsOne;

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
  FINLEY_DLL_API
  MeshAdapter(Finley_Mesh* finleyMesh=0);

  /**
     \brief
     Copy constructor.
  */
  FINLEY_DLL_API
  MeshAdapter(const MeshAdapter& in);

  /**
     \brief
     Destructor for MeshAdapter. As specified in the constructor
     this calls Finley_Mesh_deallocate for the pointer given to the 
     constructor.
  */
  FINLEY_DLL_API
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
  FINLEY_DLL_API
  void write(const std::string& fileName) const;

  /**
     \brief
     return the pointer to the underlying finley mesh structure
  */
  FINLEY_DLL_API
  Finley_Mesh* getFinley_Mesh() const;

   /**
     \brief
     Return the tag key for the given sample number.
     \param functionSpaceType Input - The function space type.
     \param sampleNo Input - The sample number.
  */
  FINLEY_DLL_API
  int getTagFromSampleNo(int functionSpaceType, int sampleNo) const;

  /**
     \brief
     Return the reference number of  the given sample number.
     \param functionSpaceType Input - The function space type.
  */
  FINLEY_DLL_API
  int* borrowSampleReferenceIDs(int functionSpaceType) const;

  /**
     \brief
     Returns true if the given integer is a valid function space type
     for this domain.
  */
  FINLEY_DLL_API
  virtual bool isValidFunctionSpaceType(int functionSpaceType) const;

  /**
     \brief
     Return a description for this domain
  */
  FINLEY_DLL_API
  virtual std::string getDescription() const;

  /**
     \brief
     Return a description for the given function space type code
  */
  FINLEY_DLL_API
  virtual std::string functionSpaceTypeAsString(int functionSpaceType) const;

  /**
     \brief
     Build the table of function space type names
  */
  FINLEY_DLL_API
  void setFunctionSpaceTypeNames();

  /**
     \brief
     Return a continuous FunctionSpace code
  */
  FINLEY_DLL_API
  virtual int getContinuousFunctionCode() const;

  /**
     \brief
     Return a function FunctionSpace code
  */
  FINLEY_DLL_API
  virtual int getFunctionCode() const;

  /**
     \brief
     Return a function with reduced integration order FunctionSpace code
  */
  FINLEY_DLL_API
  virtual int getReducedFunctionCode() const;

  /**
     \brief
     Return a function on boundary FunctionSpace code
  */
  FINLEY_DLL_API
  virtual int getFunctionOnBoundaryCode() const;

  /**
     \brief
     Return a function on boundary with reduced integration order FunctionSpace code
  */
  FINLEY_DLL_API
  virtual int getReducedFunctionOnBoundaryCode() const;

  /**
     \brief
     Return a FunctionOnContactZero code
  */
  FINLEY_DLL_API
  virtual int getFunctionOnContactZeroCode() const;

  /**
     \brief
     Return a FunctionOnContactZero code  with reduced integration order
  */
  FINLEY_DLL_API
  virtual int getReducedFunctionOnContactZeroCode() const;

  /**
     \brief
     Return a FunctionOnContactOne code
  */
  FINLEY_DLL_API
  virtual int getFunctionOnContactOneCode() const;

  /**
     \brief
     Return a FunctionOnContactOne code  with reduced integration order
  */
  FINLEY_DLL_API
  virtual int getReducedFunctionOnContactOneCode() const;

  /**
     \brief
     Return a Solution code
  */
  FINLEY_DLL_API
  virtual int getSolutionCode() const;

  /**
     \brief
     Return a ReducedSolution code
  */
  FINLEY_DLL_API
  virtual int getReducedSolutionCode() const;

  /**
     \brief
     Return a DiracDeltaFunction code
  */
  FINLEY_DLL_API
  virtual int getDiracDeltaFunctionCode() const;

  /**
		 5B
     \brief
  */
  typedef std::map<int, std::string> FunctionSpaceNamesMapType;

  /**
     \brief
  */
  FINLEY_DLL_API
  virtual int getDim() const;

  /**
     \brief
     Return the number of data points per sample, and the number of samples as a pair.
     \param functionSpaceCode Input -
  */
  FINLEY_DLL_API
  virtual std::pair<int,int> getDataShape(int functionSpaceCode) const;

  /**
     \brief
     copies the location of data points into arg. The domain of arg has to match this.
     has to be implemented by the actual Domain adapter.
  */
  FINLEY_DLL_API
  virtual void setToX(escript::Data& arg) const;

  /**
     \brief
     sets a map from a clear tag name to a tag key
     \param name Input - tag name.
     \param tag Input - tag key.
  */
  FINLEY_DLL_API
  virtual void setTagMap(const std::string& name,  int tag);

  /**
     \brief
     Return the tag key for tag name.
     \param name Input - tag name
  */
  FINLEY_DLL_API
  virtual int getTag(const std::string& name) const;

  /**
     \brief
     Returns true if name is a defined tage name. 
     \param name Input - tag name to be checked.
  */
  FINLEY_DLL_API
  virtual bool isValidTagName(const std::string& name) const;

  /**
     \brief
     Returns all tag names in a single string sperated by commas
  */
  FINLEY_DLL_API
  virtual std::string showTagNames() const;

  /**
     \brief
     assigns new location to the domain
  */
  FINLEY_DLL_API
  virtual void setNewX(const escript::Data& arg);

  /**
     \brief
     interpolates data given on source onto target where source and target have to be given on the same domain.
  */
  FINLEY_DLL_API
  virtual void interpolateOnDomain(escript::Data& target,const escript::Data& source) const;
  FINLEY_DLL_API
  virtual bool probeInterpolationOnDomain(int functionSpaceType_source,int functionSpaceType_target) const;

  /**
     \brief
     interpolates data given on source onto target where source and target are given on different domains.
     has to be implemented by the actual Domain adapter.
  */
  FINLEY_DLL_API
  virtual void interpolateACross(escript::Data& target, const escript::Data& source) const;
  FINLEY_DLL_API
  virtual bool probeInterpolationACross(int functionSpaceType_source,const AbstractDomain& targetDomain, int functionSpaceType_target) const;

  /**
     \brief
     copies the surface normals at data points into out. The actual function space to be considered
     is defined by out. out has to be defined on this.
  */
  FINLEY_DLL_API
  virtual void setToNormal(escript::Data& out) const;

  /**
     \brief
     copies the size of samples into out. The actual function space to be considered
     is defined by out. out has to be defined on this.
  */
  FINLEY_DLL_API
  virtual void setToSize(escript::Data& out) const;

  /**
     \brief
     copies the gradient of arg into grad. The actual function space to be considered
     for the gradient is defined by grad. arg and grad have to be defined on this.
  */
  FINLEY_DLL_API
  virtual void setToGradient(escript::Data& grad,const escript::Data& arg) const;

  /**
     \brief
     copies the integrals of the function defined by arg into integrals.
     arg has to be defined on this.
  */
  FINLEY_DLL_API
  virtual void setToIntegrals(std::vector<double>& integrals,const escript::Data& arg) const;

  /**
     \brief
     return the identifier of the matrix type to be used for the global stiffness matrix when a particular solver, package
     and symmetric matrix is used.
     \param solver 
     \param symmetry 
  */
  FINLEY_DLL_API
  virtual int getSystemMatrixTypeId(const int solver, const int package, const bool symmetry) const;

  /**
     \brief
     returns true if data on this domain and a function space of type functionSpaceCode has to 
     considered as cell centered data.
  */
  FINLEY_DLL_API
  virtual bool isCellOriented(int functionSpaceCode) const;

  /**
     \brief
     Saves a dictonary of Data objects to an OpenDX input file. The keywords are used as identifier
                                                                                                                                                                        
     This has to be implemented by the actual Domain adapter.
  */
  FINLEY_DLL_API
  virtual void saveDX(const std::string& filename,const boost::python::dict& arg) const;


  /**
     \brief
     Saves a dictonary of Data objects to an VTK XML input file. The keywords are used as identifier
                                                                                                                                                                        
     This has to be implemented by the actual Domain adapter.
  */
  FINLEY_DLL_API
  virtual void saveVTK(const std::string& filename,const boost::python::dict& arg) const;

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
  FINLEY_DLL_API
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
  FINLEY_DLL_API
  virtual void addPDEToRHS(escript::Data& rhs,
                     const escript::Data& X, const escript::Data& Y,
                     const escript::Data& y, const escript::Data& y_contact) const;

  /**
     \brief
    creates a SystemMatrixAdapter stiffness matrix an initializes it with zeros:
  */
  FINLEY_DLL_API
  SystemMatrixAdapter newSystemMatrix(
                      const int row_blocksize,
                      const escript::FunctionSpace& row_functionspace,
                      const int column_blocksize,
                      const escript::FunctionSpace& column_functionspace,
                      const int type) const;

  /**
     \brief returns locations in the FEM nodes
  */
  FINLEY_DLL_API
  virtual escript::Data getX() const;

  /**
     \brief return boundary normals at the quadrature point on the face elements
  */
  FINLEY_DLL_API
  virtual escript::Data getNormal() const;

  /**
     \brief returns the element size
  */
  FINLEY_DLL_API
  virtual escript::Data getSize() const;

  /**
     \brief comparison operators
  */
  FINLEY_DLL_API
  virtual bool operator==(const AbstractDomain& other) const;
  FINLEY_DLL_API
  virtual bool operator!=(const AbstractDomain& other) const;

  /**
     \brief assigns new tag newTag to all samples of functionspace with a positive
     value of mask for any its sample point.

  */
  FINLEY_DLL_API
  virtual void setTags(const int functionSpaceType, const int newTag, const escript::Data& mask) const;

 protected:

 private:

  //
  // pointer to the externally created finley mesh
  boost::shared_ptr<Finley_Mesh> m_finleyMesh;
 
  static FunctionSpaceNamesMapType m_functionSpaceTypeNames;

};

} // end of namespace

#endif
