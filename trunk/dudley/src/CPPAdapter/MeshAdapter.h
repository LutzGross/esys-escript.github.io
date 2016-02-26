
/*****************************************************************************
*
* Copyright (c) 2003-2016 by The University of Queensland
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

#if !defined dudley_MeshAdapter_20040526_H
#define dudley_MeshAdapter_20040526_H
#include "system_dep.h"

#include "dudley/Mesh.h"
#include "dudley/Dudley.h"
#include "dudley/Assemble.h"

#include "esysUtils/Esys_MPI.h"
#include "escript/AbstractContinuousDomain.h"
#include "escript/FunctionSpace.h"
#include "escript/FunctionSpaceFactory.h"

#include <boost/shared_ptr.hpp>
#include <boost/python/dict.hpp>

#include <map>
#include <vector>
#include <string>

namespace dudley {

struct null_deleter
{
  void operator()(void const *ptr) const
  {
  }
};


/**
   \brief
   MeshAdapter implements the AbstractContinuousDomain
   interface for the Dudley library.

   Description:
   MeshAdapter implements the AbstractContinuousDomain
   interface for the Dudley library.
*/
class DUDLEY_DLL_API MeshAdapter : public escript::AbstractContinuousDomain
{
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
     is deleted using a call to Dudley_Mesh_free in the
     MeshAdapter destructor.

     Throws:
     May throw an exception derived from EsysException

     \param dudleyMesh Input - A pointer to the externally constructed 
                               dudley mesh.The pointer passed to MeshAdapter
                               is deleted using a call to 
                               Dudley_Mesh_free in the MeshAdapter 
                               destructor.
    */
    MeshAdapter(Dudley_Mesh* dudleyMesh=0);

    /**
     \brief
     Copy constructor.
    */
    MeshAdapter(const MeshAdapter& in);

    /**
     \brief
     Destructor for MeshAdapter. As specified in the constructor
     this calls Dudley_Mesh_free for the pointer given to the 
     constructor.
    */
    ~MeshAdapter();

    /**
     \brief
     returns a reference to the MPI information wrapper for this domain
    */
    virtual esysUtils::JMPI getMPI() const;

    /**
     \brief
     return the number of processors used for this domain
    */
    virtual int getMPISize() const;

    /**
     \brief
     return the number MPI rank of this processor
    */

    virtual int getMPIRank() const;

    /**
     \brief
     If compiled for MPI then execute an MPI_Barrier, else do nothing
    */
    virtual void MPIBarrier() const;

    /**
     \brief
     Return true if on MPI processor 0, else false
    */
    virtual bool onMasterProcessor() const;

    MPI_Comm getMPIComm() const;

    /**
     \brief
     Write the current mesh to a file with the given name.
     \param fileName Input - The name of the file to write to.
    */
    void write(const std::string& fileName) const;

    /**
     \brief
     \param full
    */
    void Print_Mesh_Info(const bool full=false) const;

    /**
     \brief
     dumps the mesh to a file with the given name.
     \param fileName Input - The name of the file
    */
    void dump(const std::string& fileName) const;

    /**
     \brief
     return the pointer to the underlying dudley mesh structure
    */
    Dudley_Mesh* getDudley_Mesh() const;

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
    */
    const int* borrowSampleReferenceIDs(int functionSpaceType) const;

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
     Return a continuous on reduced order nodes FunctionSpace code
    */
    virtual int getReducedContinuousFunctionCode() const;

    /**
     \brief
     Return a function FunctionSpace code
    */
    virtual int getFunctionCode() const;

    /**
     \brief
     Return a function with reduced integration order FunctionSpace code
    */
    virtual int getReducedFunctionCode() const;

    /**
     \brief
     Return a function on boundary FunctionSpace code
    */
    virtual int getFunctionOnBoundaryCode() const;

    /**
     \brief
     Return a function on boundary with reduced integration order FunctionSpace code
    */
    virtual int getReducedFunctionOnBoundaryCode() const;

    /**
     \brief
     Return a FunctionOnContactZero code
    */
    virtual int getFunctionOnContactZeroCode() const;

    /**
     \brief
     Return a FunctionOnContactZero code  with reduced integration order
    */
    virtual int getReducedFunctionOnContactZeroCode() const;

    /**
     \brief
     Return a FunctionOnContactOne code
    */
    virtual int getFunctionOnContactOneCode() const;

    /**
     \brief
     Return a FunctionOnContactOne code  with reduced integration order
    */
    virtual int getReducedFunctionOnContactOneCode() const;

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
     Return a DiracDeltaFunctions code
    */
    virtual int getDiracDeltaFunctionsCode() const;

    /**
     \brief
    */
    typedef std::map<int, std::string> FunctionSpaceNamesMapType;

    /**
     \brief
    */
    virtual int getDim() const;

    /**
     \brief
      Returns a status indicator of the domain. The status identifier should be unique over 
      the live time if the object but may be updated if changes to the domain happen, e.g. 
      modifications to its geometry. 

     This has to be implemented by the actual Domain adapter.
    */
    virtual StatusType getStatus() const;

    /**
     \brief
     Return the number of data points summed across all MPI processes
    */
    virtual int getNumDataPointsGlobal() const;

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
     sets a map from a clear tag name to a tag key
     \param name Input - tag name.
     \param tag Input - tag key.
    */
    virtual void setTagMap(const std::string& name, int tag);

    /**
     \brief
     Return the tag key for tag name.
     \param name Input - tag name
    */
    virtual int getTag(const std::string& name) const;

    /**
     \brief
     Returns true if name is a defined tage name. 
     \param name Input - tag name to be checked.
    */
    virtual bool isValidTagName(const std::string& name) const;

    /**
     \brief
     Returns all tag names in a single string sperated by commas
    */
    virtual std::string showTagNames() const;

    /**
     \brief
     assigns new location to the domain
    */
    virtual void setNewX(const escript::Data& arg);

    /**
     \brief
     interpolates data given on source onto target where source and target have to be given on the same domain.
    */
    virtual void interpolateOnDomain(escript::Data& target, const escript::Data& source) const;

    virtual bool probeInterpolationOnDomain(int functionSpaceType_source,int functionSpaceType_target) const;

    virtual signed char preferredInterpolationOnDomain(int functionSpaceType_source,int functionSpaceType_target) const;

    /**
    \brief given a vector of FunctionSpace typecodes, pass back a code which then can all be interpolated to.
    \return true is result is valid, false if not
    */
    bool commonFunctionSpace(const std::vector<int>& fs, int& resultcode) const;

    /**
     \brief
     interpolates data given on source onto target where source and target are given on different domains.
    */
    virtual void interpolateAcross(escript::Data& target, const escript::Data& source) const;

    /**
     \brief determines whether interpolation from source to target is possible.
  */
  virtual bool probeInterpolationAcross(int functionSpaceType_source,const escript::AbstractDomain& targetDomain, int functionSpaceType_target) const;

    /**
     \brief
     copies the surface normals at data points into out. The actual function space to be considered
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
     return the identifier of the matrix type to be used for the global
     stiffness matrix when a particular solver, package, preconditioner,
     and symmetric matrix is used.
     
     \param options a SolverBuddy instance with the desired options set
    */
    virtual int getSystemMatrixTypeId(const boost::python::object& options) const;

    /**
     \brief
     return the identifier of the transport problem type to be used when a particular solver, perconditioner, package
     and symmetric matrix is used.
     \param solver 
     \param preconditioner
     \param package
     \param symmetry 
    */
    virtual int getTransportTypeId(int solver, int preconditioner, int package,
                                   bool symmetry) const;

    /**
     \brief
     returns true if data on this domain and a function space of type functionSpaceCode has to 
     considered as cell centered data.
    */
    virtual bool isCellOriented(int functionSpaceCode) const;

    virtual bool ownSample(int fs_code, index_t id) const;

    /**
     \brief
     adds a PDE onto the stiffness matrix mat and a rhs 
    */
    virtual void addPDEToSystem(
                     escript::AbstractSystemMatrix& mat, escript::Data& rhs,
                     const escript::Data& A, const escript::Data& B,
                     const escript::Data& C, const escript::Data& D,
                     const escript::Data& X, const escript::Data& Y,
                     const escript::Data& d, const escript::Data& y,
                     const escript::Data& d_contact,
                     const escript::Data& y_contact,
                     const escript::Data& d_dirac,
                     const escript::Data& y_dirac) const;

    /**
     \brief
     adds a PDE onto the lumped stiffness matrix matrix
    */
    virtual void addPDEToLumpedSystem(escript::Data& mat,
                                      const escript::Data& D,
                                      const escript::Data& d,
                                      const escript::Data& d_dirac,
                                      const bool useHRZ) const;

    /**
     \brief
     adds a PDE onto the stiffness matrix mat and a rhs 
    */
    virtual void addPDEToRHS(escript::Data& rhs, const escript::Data& X,
                             const escript::Data& Y, const escript::Data& y,
                             const escript::Data& y_contact,
                             const escript::Data& y_dirac) const;

    /**
     \brief
     adds a PDE onto a transport problem
    */
    virtual void addPDEToTransportProblem(
                     escript::AbstractTransportProblem& tp, escript::Data& source, 
                     const escript::Data& M,
                     const escript::Data& A, const escript::Data& B, const escript::Data& C,const  escript::Data& D,
                     const  escript::Data& X,const  escript::Data& Y,
                     const escript::Data& d, const escript::Data& y,
                     const escript::Data& d_contact,const escript::Data& y_contact,
                     const escript::Data& d_dirac,const escript::Data& y_dirac) const;

    /**
     \brief
     creates a stiffness matrix and initializes it with zeros
    */
    escript::ASM_ptr newSystemMatrix(
                      int row_blocksize,
                      const escript::FunctionSpace& row_functionspace,
                      int column_blocksize,
                      const escript::FunctionSpace& column_functionspace,
                      int type) const;

    /**
     \brief 
      creates a TransportProblemAdapter 
    */
    escript::ATP_ptr newTransportProblem(
                      const int blocksize,
                      const escript::FunctionSpace& functionspace,
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

    /**
     \brief comparison operators
    */
    virtual bool operator==(const escript::AbstractDomain& other) const;
    virtual bool operator!=(const escript::AbstractDomain& other) const;

    /**
     \brief assigns new tag newTag to all samples of functionspace with a positive
     value of mask for any its sample point.
    */
    virtual void setTags(int functionSpaceType, int newTag, const escript::Data& mask) const;

    /**
      \brief
          return the number of tags in use and a pointer to an array with the number of tags in use
    */
    virtual int getNumberOfTagsInUse(int functionSpaceCode) const;

    virtual const int* borrowListOfTagsInUse(int functionSpaceCode) const;

    /**
     \brief Checks if this domain allows tags for the specified functionSpaceCode.
    */
    virtual bool canTag(int functionSpaceCode) const;

    /**
     \brief returns the approximation order used for a function space functionSpaceCode
    */
    virtual int getApproximationOrder(int functionSpaceCode) const;

    bool supportsContactElements() const;

    virtual escript::Data randomFill(const escript::DataTypes::ShapeType& shape,
       const escript::FunctionSpace& what, long seed, const boost::python::tuple& filter) const;         
  
private:
    void extractArgsFromDict(const boost::python::dict& arg, int& numData,
                             char**& names, escript::Data*& data,
                             escript::Data**& dataPtr) const;

    //
    // pointer to the externally created dudley mesh
    boost::shared_ptr<Dudley_Mesh> m_dudleyMesh;
 
    static FunctionSpaceNamesMapType m_functionSpaceTypeNames;

};

} // end of namespace

#endif

