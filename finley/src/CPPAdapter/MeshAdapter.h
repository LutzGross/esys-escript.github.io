
/*****************************************************************************
*
* Copyright (c) 2003-2012 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development since 2012 by School of Earth Sciences
*
*****************************************************************************/


#if !defined finley_MeshAdapter_20040526_H
#define finley_MeshAdapter_20040526_H
#include "system_dep.h"

extern "C" {
#include "finley/Mesh.h"
#include "finley/Finley.h"
#include "finley/Assemble.h"
#include "esysUtils/Esys_MPI.h"
}

#include "FinleyError.h"
#include "FinleyAdapterException.h"

#include <pasowrap/SystemMatrixAdapter.h>
#include <pasowrap/TransportProblemAdapter.h>
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

namespace finley {
  
// These are friends implemented in MeshAdapterFactory.cpp  
// They are only fwd declared here so that vis.studio will accept the freind decls
FINLEY_DLL_API
escript::Domain_ptr brick(int n0,int n1,int n2,int order,
		    double l0,double l1,double l2,
		    int periodic0,int periodic1,
		    int periodic2,
		    int integrationOrder,
		    int reducedIntegrationOrder,
		    int useElementsOnFace,
                    int useFullElementOrder,
                     int optimize, 
		    const std::vector<double>& points,
		    const std::vector<int>& tags,
		    const std::map<std::string, int>& tagnamestonums			   
		    );
		    
FINLEY_DLL_API		    
escript::Domain_ptr  rectangle(int n0,int n1,int order,
			double l0, double l1,
			int periodic0,int periodic1,
			int integrationOrder,
                        int reducedIntegrationOrder,
			int useElementsOnFace,
		        int useFullElementOrder,
                        int optimize,
			const std::vector<double>& points,
			const std::vector<int>& tags,
			const std::map<std::string, int>& tagnamestonums);	  
  
}



namespace finley {

struct null_deleter
{
  void operator()(void const *ptr) const
  {
  }
};


/**
   \brief implements the AbstractContinuousDomain interface for the Finley
          library.
*/

class MeshAdapter : public escript::AbstractContinuousDomain {

 public:

  //
  // Codes for function space types supported
  static const int DegreesOfFreedom;
  static const int ReducedDegreesOfFreedom;
  static const int Nodes;
  static const int ReducedNodes;
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
     is deleted using a call to Finley_Mesh_free in the
     MeshAdapter destructor.

     Throws:
     May throw an exception derived from EsysException

     \param finleyMesh Input - A pointer to the externally constructed 
                               finley mesh.The pointer passed to MeshAdapter
                               is deleted using a call to 
                               Finley_Mesh_free in the MeshAdapter 
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
     this calls Finley_Mesh_free for the pointer given to the 
     constructor.
  */
  FINLEY_DLL_API
  ~MeshAdapter();

  /**
     \brief
     return the number of processors used for this domain
  */
  FINLEY_DLL_API
  virtual int getMPISize() const;
  /**
     \brief
     return the number MPI rank of this processor
  */

  FINLEY_DLL_API
  virtual int getMPIRank() const;

  /**
     \brief
     If compiled for MPI then execute an MPI_Barrier, else do nothing
  */

  FINLEY_DLL_API
  virtual void MPIBarrier() const;

  /**
     \brief
     Return true if on MPI processor 0, else false
  */

  FINLEY_DLL_API
  virtual bool onMasterProcessor() const;

  FINLEY_DLL_API
#ifdef ESYS_MPI
  MPI_Comm
#else
  unsigned int
#endif
  getMPIComm() const;

  /**
     \brief
     Write the current mesh to a file with the given name.
     \param fileName Input - The name of the file to write to.
  */
  FINLEY_DLL_API
  void write(const std::string& fileName) const;

  /**
     \brief
     \param full
  */
  FINLEY_DLL_API
  void Print_Mesh_Info(const bool full=false) const;

  /**
     \brief
     dumps the mesh to a file with the given name.
     \param fileName Input - The name of the file
  */
  FINLEY_DLL_API
  void dump(const std::string& fileName) const;

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
  const int* borrowSampleReferenceIDs(int functionSpaceType) const;

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
     Return a continuous on reduced order nodes FunctionSpace code
  */
  FINLEY_DLL_API
  virtual int getReducedContinuousFunctionCode() const;

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
     Return a DiracDeltaFunctions code
  */
  FINLEY_DLL_API
  virtual int getDiracDeltaFunctionsCode() const;

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
      Returns a status indicator of the domain. The status identifier should be unique over 
      the live time if the object but may be updated if changes to the domain happen, e.g. 
      modifications to its geometry. 

     This has to be implemented by the actual Domain adapter.
  */
  FINLEY_DLL_API
  virtual StatusType getStatus() const;


  /**
     \brief
     Return the number of data points summed across all MPI processes
  */
  FINLEY_DLL_API
  virtual int getNumDataPointsGlobal() const;

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
    \brief given a vector of FunctionSpace typecodes, pass back a code which then can all be interpolated to.
    \return true is result is valid, false if not
  */
  FINLEY_DLL_API
  bool
  commonFunctionSpace(const std::vector<int>& fs, int& resultcode) const;

  /**
     \brief
     interpolates data given on source onto target where source and target are given on different domains.
     has to be implemented by the actual Domain adapter.
  */
  FINLEY_DLL_API
  virtual void interpolateACross(escript::Data& target, const escript::Data& source) const;

  /**
  \brief determines whether interpolation from source to target is possible.
  Must be implemented by the actual Domain adapter
  */
  FINLEY_DLL_API
  virtual bool probeInterpolationACross(int functionSpaceType_source,const escript::AbstractDomain& targetDomain, int functionSpaceType_target) const;

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
     return the identifier of the matrix type to be used for the global stiffness matrix when a particular solver, package, perconditioner,
     and symmetric matrix is used.
     \param solver 
     \param preconditioner
     \param package
     \param symmetry 
  */
  FINLEY_DLL_API
  virtual int getSystemMatrixTypeId(const int solver, const int preconditioner, const int package, const bool symmetry) const;

  /**
     \brief
     return the identifier of the transport problem type to be used when a particular solver, perconditioner, package
     and symmetric matrix is used.
     \param solver 
     \param preconditioner
     \param package
     \param symmetry 
  */
  FINLEY_DLL_API
  virtual int getTransportTypeId(const int solver, const int preconditioner, const int package, const bool symmetry) const;

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
  virtual void saveVTK(const std::string& filename,const boost::python::dict& arg,  const std::string& metadata, const std::string& metadata_schema) const;

  FINLEY_DLL_API
  virtual bool ownSample(int fs_code, index_t id) const;

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
                     escript::AbstractSystemMatrix& mat, escript::Data& rhs,
                     const escript::Data& A, const escript::Data& B, const escript::Data& C, 
                     const escript::Data& D, const escript::Data& X, const escript::Data& Y,
                     const escript::Data& d, const escript::Data& y,
                     const escript::Data& d_contact, const escript::Data& y_contact,
                     const escript::Data& d_dirac, const escript::Data& y_dirac) const;
  /**
     \brief
     adds a PDE onto the lumped stiffness matrix matrix
  */
  FINLEY_DLL_API
  virtual void addPDEToLumpedSystem(
                     escript::Data& mat,
                     const escript::Data& D, 
                     const escript::Data& d,
                     const escript::Data& d_dirac,
		     const bool useHRZ) const;

  /**
     \brief
     adds a PDE onto the stiffness matrix mat and a rhs 
  */
  FINLEY_DLL_API
  virtual void addPDEToRHS(escript::Data& rhs,
                     const escript::Data& X, const escript::Data& Y,
                     const escript::Data& y, const escript::Data& y_contact, const escript::Data& y_dirac) const;
  /**
     \brief
     adds a PDE onto a transport problem
  */

  FINLEY_DLL_API
  virtual void addPDEToTransportProblem(
                     escript::AbstractTransportProblem& tp, escript::Data& source, 
                     const escript::Data& M,
                     const escript::Data& A, const escript::Data& B, const escript::Data& C,const  escript::Data& D,
                     const  escript::Data& X,const  escript::Data& Y,
                     const escript::Data& d, const escript::Data& y,
                     const escript::Data& d_contact,const escript::Data& y_contact, const escript::Data& d_dirac,const escript::Data& y_dirac) const;


  /**
     \brief
    creates a SystemMatrixAdapter stiffness matrix and initializes it with zeros:
  */
  FINLEY_DLL_API
  escript::ASM_ptr newSystemMatrix(
                      const int row_blocksize,
                      const escript::FunctionSpace& row_functionspace,
                      const int column_blocksize,
                      const escript::FunctionSpace& column_functionspace,
                      const int type) const;
  /**
   \brief 
    creates a TransportProblemAdapter 

  */

  FINLEY_DLL_API
  escript::ATP_ptr newTransportProblem(
                      const int blocksize,
                      const escript::FunctionSpace& functionspace,
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
  virtual bool operator==(const escript::AbstractDomain& other) const;
  FINLEY_DLL_API
  virtual bool operator!=(const escript::AbstractDomain& other) const;

  /**
     \brief assigns new tag newTag to all samples of functionspace with a positive
     value of mask for any its sample point.

  */
  FINLEY_DLL_API
  virtual void setTags(const int functionSpaceType, const int newTag, const escript::Data& mask) const;

  /**
      \brief
          return the number of tags in use and a pointer to an array with the number of tags in use
  */
  FINLEY_DLL_API
  virtual int getNumberOfTagsInUse(int functionSpaceCode) const;

  FINLEY_DLL_API 
  virtual const int* borrowListOfTagsInUse(int functionSpaceCode) const;


  /**
     \brief Checks if this domain allows tags for the specified functionSpaceCode.
  */
  FINLEY_DLL_API
  virtual
  bool canTag(int functionSpaceCode) const;

   /**
   \brief returns the approximation order used for a function space functionSpaceCode
   */

  FINLEY_DLL_API
  virtual 
  int getApproximationOrder(const int functionSpaceCode) const;

  FINLEY_DLL_API
  bool supportsContactElements() const;
  
  
  private:
  
  /**
   \brief  adds points to support more Dirac delta function.
   
   Do NOT call these at any time other than construction!
   Using them later creates consistancy problems
   */
  FINLEY_DLL_API
  void addDiracPoints( const std::vector<double>& points, const std::vector<int>& tags) const;
//  FINLEY_DLL_API
//  void addDiracPoint( const boost::python::list& points, const int tag=-1) const;
//   FINLEY_DLL_API
//   void addDiracPointWithTagName( const boost::python::list& points, const std::string& tag) const;

 protected:

 private:
  void extractArgsFromDict(const boost::python::dict& arg, int& numData,
                             char**& names, escriptDataC*& data,
                             escriptDataC**& dataPtr) const;

  //
  // pointer to the externally created finley mesh
  boost::shared_ptr<Finley_Mesh> m_finleyMesh;
  
  // This is only provided so that the friends below can add tags during construction
  // do not use for any other purpose
  boost::shared_ptr<Finley_Mesh> getMesh()
  {
      return m_finleyMesh;
  }
 
  static FunctionSpaceNamesMapType m_functionSpaceTypeNames;

  friend escript::Domain_ptr finley::brick(int n0,int n1,int n2,int order,
		    double l0,double l1,double l2,
		    int periodic0,int periodic1,
		    int periodic2,
		    int integrationOrder,
		    int reducedIntegrationOrder,
		    int useElementsOnFace,
                    int useFullElementOrder,
                     int optimize, 
		    const std::vector<double>& points,
		    const std::vector<int>& tags,
		    const std::map<std::string, int>& tagnamestonums			   
		    );
		    
		    
friend   escript::Domain_ptr  finley::rectangle(int n0,int n1,int order,
			double l0, double l1,
			int periodic0,int periodic1,
			int integrationOrder,
                        int reducedIntegrationOrder,
			int useElementsOnFace,
		        int useFullElementOrder,
                        int optimize,
			const std::vector<double>& points,
			const std::vector<int>& tags,
			const std::map<std::string, int>& tagnamestonums);						
};

// Do not use this class. It is a convenience wrapper for the dataexporter.
class FINLEY_DLL_API ReferenceElementSetWrapper {
 public:
  ReferenceElementSetWrapper(Finley_ElementTypeId id, index_t order,
                             index_t reducedOrder);
  ~ReferenceElementSetWrapper();

  Finley_ReferenceElementSet* getElementSet() const { return m_refSet; }

 private:
  Finley_ReferenceElementSet* m_refSet;
};


} // end of namespace

#endif
