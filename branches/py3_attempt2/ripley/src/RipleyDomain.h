
/*******************************************************
*
* Copyright (c) 2003-2012 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/

#ifndef __RIPLEY_DOMAIN_H__
#define __RIPLEY_DOMAIN_H__

#include <ripley/Ripley.h>
#include <ripley/RipleyException.h>
#include <escript/AbstractContinuousDomain.h>
#include <escript/Data.h>
#include <escript/FunctionSpace.h>

struct Paso_Pattern;
struct Paso_SystemMatrixPattern;
struct Paso_SystemMatrix;

namespace ripley {

/**
   \brief
   RipleyDomain extends the AbstractContinuousDomain interface
   for the Ripley library and is the base class for Rectangle and Brick.
*/

class RipleyDomain : public escript::AbstractContinuousDomain
{
public:

    /**
       \brief recovers mesh from a file created with the dump() method
       \param filename the name of the file
    */
    RIPLEY_DLL_API
    static escript::Domain_ptr loadMesh(const std::string& filename);

    /**
       \brief reads a mesh from a file created with the write() method
       \param filename the name of the file
    */
    RIPLEY_DLL_API
    static escript::Domain_ptr readMesh(const std::string& filename);

    /**
       \brief
       Constructor with number of dimensions. Allocates MPI info structure.
    */
    RIPLEY_DLL_API
    RipleyDomain(dim_t dim);

    /**
       \brief
       Destructor for RipleyDomain
    */
    RIPLEY_DLL_API
    ~RipleyDomain();

    /**
       \brief
       returns a description for this domain
    */
    RIPLEY_DLL_API
    virtual std::string getDescription() const;

    /**
       \brief
       returns the number of processors used for this domain
    */
    RIPLEY_DLL_API
    virtual int getMPISize() const { return m_mpiInfo->size; }

    /**
       \brief
       returns the MPI rank of this processor
    */
    RIPLEY_DLL_API
    virtual int getMPIRank() const { return m_mpiInfo->rank; }

    /**
       \brief
       if compiled for MPI then executes an MPI_Barrier, else does nothing
    */
    RIPLEY_DLL_API
    virtual void MPIBarrier() const {
#ifdef ESYS_MPI
        MPI_Barrier(m_mpiInfo->comm);
#endif
    }

    /**
       \brief
       returns true if on MPI processor 0, else false
    */
    RIPLEY_DLL_API
    virtual bool onMasterProcessor() const { return getMPIRank()==0; }

    /**
       \brief
       returns the MPI communicator
    */
    RIPLEY_DLL_API
#ifdef ESYS_MPI
    MPI_Comm
#else
    unsigned int
#endif
    getMPIComm() const {
#ifdef ESYS_MPI
        return m_mpiInfo->comm;
#else
        return 0;
#endif
    }

    /**
       \brief
       returns true if the argument is a valid function space type for this
       domain
    */
    RIPLEY_DLL_API
    virtual bool isValidFunctionSpaceType(int fsType) const;

    /**
       \brief
       returns a description for the given function space type code
    */
    RIPLEY_DLL_API
    virtual std::string functionSpaceTypeAsString(int fsType) const;

    /**
       \brief
       returns the number of spatial dimensions of the domain
    */
    RIPLEY_DLL_API
    virtual int getDim() const { return m_numDim; }

    /**
       \brief equality operator
    */
    RIPLEY_DLL_API
    virtual bool operator==(const escript::AbstractDomain& other) const;

    /**
       \brief inequality operator
    */
    RIPLEY_DLL_API
    virtual bool operator!=(const escript::AbstractDomain& other) const {
        return !(operator==(other));
    }

    /**
       \brief
       writes the current mesh to a file with the given name
       \param filename The name of the file to write to
    */
    RIPLEY_DLL_API
    void write(const std::string& filename) const;

    /**
       \brief
       dumps the mesh to a file with the given name
       \param filename The name of the output file
    */
    RIPLEY_DLL_API
    void dump(const std::string& filename) const;

    /**
       \brief
       returns the number of data points per sample, and the number of samples
       as a pair.
       \param fsType The function space type
    */
    RIPLEY_DLL_API
    virtual std::pair<int,int> getDataShape(int fsType) const;

    /**
       \brief
       returns the tag key for the given sample number
       \param fsType The function space type
       \param sampleNo The sample number
    */
    RIPLEY_DLL_API
    int getTagFromSampleNo(int fsType, int sampleNo) const;

    /**
       \brief
       sets a map from a clear tag name to a tag key
       \param name tag name
       \param tag tag key
    */
    RIPLEY_DLL_API
    virtual void setTagMap(const std::string& name, int tag) {
        m_tagMap[name] = tag;
    }

    /**
       \brief
       returns the tag key for tag name
       \param name tag name
    */
    RIPLEY_DLL_API
    virtual int getTag(const std::string& name) const {
        if (m_tagMap.find(name) != m_tagMap.end()) {
            return m_tagMap.find(name)->second;
        } else {
            throw RipleyException("getTag: invalid tag name");
        }
    }

    /**
       \brief
       returns true if name is a defined tag name
       \param name tag name to be checked
    */
    RIPLEY_DLL_API
    virtual bool isValidTagName(const std::string& name) const {
        return (m_tagMap.find(name)!=m_tagMap.end());
    }

    /**
       \brief
       returns all tag names in a single string separated by commas
    */
    RIPLEY_DLL_API
    virtual std::string showTagNames() const;

    /**
       \brief
       returns the array of reference numbers for a function space type
       \param fsType The function space type
    */
    RIPLEY_DLL_API
    const int* borrowSampleReferenceIDs(int fsType) const;

    /**
       \brief
       assigns new location to the domain.
       \note This is not supported in Ripley
    */
    RIPLEY_DLL_API
    virtual void setNewX(const escript::Data& arg);

    /**
       \brief
       interpolates data given on source onto target where source and target
       have to be given on the same domain
    */
    RIPLEY_DLL_API
    virtual void interpolateOnDomain(escript::Data& target, const escript::Data& source) const;

    /**
       \brief
       returns true if data on fsType_source can be interpolated onto
       fsType_target, false otherwise
    */
    RIPLEY_DLL_API
    virtual bool probeInterpolationOnDomain(int fsType_source, int fsType_target) const;

    /**
       \brief
       given a vector of FunctionSpace type codes, passes back a code which all
       can be interpolated to
       \return true if result is valid, false if not
    */
    RIPLEY_DLL_API
    bool
    commonFunctionSpace(const std::vector<int>& fs, int& resultcode) const;

    /**
       \brief
       interpolates data given on source onto target where source and target
       are given on different domains
    */
    RIPLEY_DLL_API
    virtual void interpolateACross(escript::Data& target, const escript::Data& source) const;

    /**
       \brief
       determines whether interpolation from source to target is possible
    */
    RIPLEY_DLL_API
    virtual bool probeInterpolationACross(int, const escript::AbstractDomain&, int) const;

    /**
       \brief
       returns locations in the FEM nodes
    */
    RIPLEY_DLL_API
    virtual escript::Data getX() const;

    /**
       \brief
       returns boundary normals at the quadrature point on the face elements
    */
    RIPLEY_DLL_API
    virtual escript::Data getNormal() const;

    /**
       \brief returns the element size
    */
    RIPLEY_DLL_API
    virtual escript::Data getSize() const;

    /**
       \brief
       copies the location of data points into arg. The domain of arg has to
       match this domain.
    */
    RIPLEY_DLL_API
    virtual void setToX(escript::Data& arg) const;

    /**
       \brief
       copies the surface normals at data points into out. The actual function
       space to be considered is defined by out. out has to be defined on this
       domain.
    */
    RIPLEY_DLL_API
    virtual void setToNormal(escript::Data& out) const;

    /**
       \brief
       copies the size of samples into out. The actual function space to be
       considered is defined by out. out has to be defined on this domain.
    */
    RIPLEY_DLL_API
    virtual void setToSize(escript::Data& out) const;

    /**
       \brief
       copies the gradient of 'in' into 'out'. The actual function space to be
       considered for the gradient is defined by 'in'. Both arguments have to
       be defined on this domain.
    */
    RIPLEY_DLL_API
    virtual void setToGradient(escript::Data& out, const escript::Data& in) const;

    /**
       \brief
       returns true if this rank owns the sample id on given function space
    */
    RIPLEY_DLL_API
    virtual bool ownSample(int fsType, index_t id) const;

    /**
       \brief
       returns the number of data points summed across all MPI processes
    */
    RIPLEY_DLL_API
    virtual int getNumDataPointsGlobal() const;

    /**
       \brief
       assigns new tag newTag to all samples of given function space with a
       positive value of mask for any of its sample points
    */
    RIPLEY_DLL_API
    virtual void setTags(const int fsType, const int newTag, const escript::Data& mask) const;

    /**
       \brief
       returns true if data on this domain and given function space type has
       to be considered as cell centered data
    */
    RIPLEY_DLL_API
    virtual bool isCellOriented(int fsType) const;

    /**
       \brief
       returns a status indicator of the domain. The status identifier should
       be unique over the lifetime of the object but may be updated if changes
       to the domain happen, e.g. modifications to its geometry.
    */
    RIPLEY_DLL_API
    virtual StatusType getStatus() const { return m_status; }

    /**
       \brief
       returns the number of tags in use for a function space type
    */
    RIPLEY_DLL_API
    virtual int getNumberOfTagsInUse(int fsType) const;

    /**
       \brief
       returns a pointer to the list of tags in use for a function space type
    */
    RIPLEY_DLL_API
    virtual const int* borrowListOfTagsInUse(int fsType) const;

    /**
       \brief
       checks if this domain allows tags for the specified function space type
    */
    RIPLEY_DLL_API
    virtual bool canTag(int fsType) const;

    /**
       \brief
       returns the approximation order used for a function space
    */
    RIPLEY_DLL_API
    virtual int getApproximationOrder(const int fsType) const { return 1; }

    /**
       \brief
       returns true if this domain supports contact elements, false otherwise
    */
    RIPLEY_DLL_API
    virtual bool supportsContactElements() const { return false; }

    /**
       \brief
       returns a continuous FunctionSpace code
    */
    RIPLEY_DLL_API
    virtual int getContinuousFunctionCode() const { return Nodes; }

    /**
       \brief
       returns a continuous on reduced order nodes FunctionSpace code
    */
    RIPLEY_DLL_API
    virtual int getReducedContinuousFunctionCode() const { return ReducedNodes; }

    /**
       \brief
       returns a function FunctionSpace code
    */
    RIPLEY_DLL_API
    virtual int getFunctionCode() const { return Elements; }

    /**
       \brief
       returns a function with reduced integration order FunctionSpace code
    */
    RIPLEY_DLL_API
    virtual int getReducedFunctionCode() const { return ReducedElements; }

    /**
       \brief
       returns a function on boundary FunctionSpace code
    */
    RIPLEY_DLL_API
    virtual int getFunctionOnBoundaryCode() const { return FaceElements; }

    /**
       \brief
       returns a function on boundary with reduced integration order
       FunctionSpace code
    */
    RIPLEY_DLL_API
    virtual int getReducedFunctionOnBoundaryCode() const { return ReducedFaceElements; }

    /**
       \brief
       return a FunctionOnContactZero code
    */
    RIPLEY_DLL_API
    virtual int getFunctionOnContactZeroCode() const {
        throw RipleyException("Ripley does not support contact elements");
    }

    /**
       \brief
       returns a FunctionOnContactZero code with reduced integration order
    */
    RIPLEY_DLL_API
    virtual int getReducedFunctionOnContactZeroCode() const {
        throw RipleyException("Ripley does not support contact elements");
    }

    /**
       \brief
       returns a FunctionOnContactOne code
    */
    RIPLEY_DLL_API
    virtual int getFunctionOnContactOneCode() const {
        throw RipleyException("Ripley does not support contact elements");
    }

    /**
       \brief
       returns a FunctionOnContactOne code with reduced integration order
    */
    RIPLEY_DLL_API
    virtual int getReducedFunctionOnContactOneCode() const {
        throw RipleyException("Ripley does not support contact elements");
    }

    /**
       \brief
       returns a Solution FunctionSpace code
    */
    RIPLEY_DLL_API
    virtual int getSolutionCode() const { return DegreesOfFreedom; }

    /**
       \brief
       returns a ReducedSolution FunctionSpace code
    */
    RIPLEY_DLL_API
    virtual int getReducedSolutionCode() const { return ReducedDegreesOfFreedom; }

    /**
       \brief
       returns a DiracDeltaFunctions FunctionSpace code
    */
    RIPLEY_DLL_API
    virtual int getDiracDeltaFunctionsCode() const { return Points; }

    /**
       \brief
       returns the identifier of the matrix type to be used for the global
       stiffness matrix when a particular solver, package, preconditioner,
       and symmetric matrix is used
       \param solver
       \param preconditioner
       \param package
       \param symmetry
    */
    RIPLEY_DLL_API
    virtual int getSystemMatrixTypeId(const int solver, const int preconditioner, const int package, const bool symmetry) const;

    /**
       \brief
       returns the identifier of the transport problem type to be used when a
       particular solver, preconditioner, package and symmetric matrix is used
       \param solver
       \param preconditioner
       \param package
       \param symmetry
    */
    RIPLEY_DLL_API
    virtual int getTransportTypeId(const int solver, const int preconditioner, const int package, const bool symmetry) const;

    /**
       \brief
       copies the integrals of the function defined by arg into integrals.
       arg has to be defined on this domain.
    */
    RIPLEY_DLL_API
    virtual void setToIntegrals(std::vector<double>& integrals, const escript::Data& arg) const;

    /**
       \brief
       adds a PDE onto the stiffness matrix mat and rhs
    */
    RIPLEY_DLL_API
    virtual void addPDEToSystem(escript::AbstractSystemMatrix& mat,
            escript::Data& rhs, const escript::Data& A, const escript::Data& B,
            const escript::Data& C, const escript::Data& D,
            const escript::Data& X, const escript::Data& Y,
            const escript::Data& d, const escript::Data& y,
            const escript::Data& d_contact, const escript::Data& y_contact,
            const escript::Data& d_dirac, const escript::Data& y_dirac) const;

    /**
       \brief
       adds a PDE onto rhs
    */
    RIPLEY_DLL_API
    virtual void addPDEToRHS(escript::Data& rhs, const escript::Data& X,
            const escript::Data& Y, const escript::Data& y,
            const escript::Data& y_contact, const escript::Data& y_dirac) const;

    /**
       \brief
       adds a PDE onto a transport problem
    */
    RIPLEY_DLL_API
    virtual void addPDEToTransportProblem(escript::AbstractTransportProblem& tp,
            escript::Data& source, const escript::Data& M,
            const escript::Data& A, const escript::Data& B,
            const escript::Data& C, const escript::Data& D,
            const escript::Data& X, const escript::Data& Y,
            const escript::Data& d, const escript::Data& y,
            const escript::Data& d_contact, const escript::Data& y_contact,
            const escript::Data& d_dirac, const escript::Data& y_dirac) const;


    /**
       \brief
       creates a stiffness matrix and initializes it with zeros
    */
    RIPLEY_DLL_API
    virtual escript::ASM_ptr newSystemMatrix(const int row_blocksize,
            const escript::FunctionSpace& row_functionspace,
            const int column_blocksize,
            const escript::FunctionSpace& column_functionspace, const int type) const;

    /**
     \brief
      creates a transport problem
    */
    RIPLEY_DLL_API
    virtual escript::ATP_ptr newTransportProblem(
            const int blocksize, const escript::FunctionSpace& functionspace,
            const int type) const;

    /**
       \brief
       writes information about the mesh to standard output
       \param full whether to print additional data
    */
    RIPLEY_DLL_API
    virtual void Print_Mesh_Info(const bool full=false) const;

    /**
       \brief
       returns the number of nodes per MPI rank in each dimension
    */
    RIPLEY_DLL_API
    virtual IndexVector getNumNodesPerDim() const;

    /**
       \brief
       returns the number of elements per MPI rank in each dimension
    */
    RIPLEY_DLL_API
    virtual IndexVector getNumElementsPerDim() const;

    /**
       \brief
       returns the number of face elements in the order
       (left,right,bottom,top,[front,back]) on current MPI rank
    */
    RIPLEY_DLL_API
    virtual IndexVector getNumFacesPerBoundary() const;

    /**
       \brief
       returns the node distribution vector
    */
    RIPLEY_DLL_API
    virtual IndexVector getNodeDistribution() const;

    /**
       \brief
       returns the number of spatial subdivisions in each dimension
    */
    RIPLEY_DLL_API
    virtual IndexVector getNumSubdivisionsPerDim() const;

    /**
       \brief
       returns the first coordinate value and the node spacing along given
       dimension as a pair
    */
    RIPLEY_DLL_API
    virtual std::pair<double,double> getFirstCoordAndSpacing(dim_t dim) const;

protected:
    /// returns the number of nodes per MPI rank
    virtual dim_t getNumNodes() const;

    /// returns the number of elements per MPI rank
    virtual dim_t getNumElements() const;

    /// returns the number of degrees of freedom per MPI rank
    virtual dim_t getNumDOF() const;

    /// returns the number of face elements on current MPI rank
    virtual dim_t getNumFaceElements() const;

    /// inserts the nodes that share an element with 'node' into 'index' and
    /// returns the number of these neighbours
    virtual dim_t insertNeighbourNodes(IndexVector& index, index_t node) const;

    /// populates the data object 'arg' with the node coordinates
    virtual void assembleCoordinates(escript::Data& arg) const;

    /// computes the gradient of 'in' and puts the result in 'out'
    virtual void assembleGradient(escript::Data& out, escript::Data& in) const;

    /// copies the integrals of the function defined by 'arg' into 'integrals'
    virtual void assembleIntegrate(std::vector<double>& integrals, escript::Data& arg) const;

    /// assembles a single PDE into the system matrix 'mat' and the right hand
    /// side 'rhs'
    virtual void assemblePDESingle(Paso_SystemMatrix* mat, escript::Data& rhs,
            const escript::Data& A, const escript::Data& B,
            const escript::Data& C, const escript::Data& D,
            const escript::Data& X, const escript::Data& Y) const;

    /// assembles boundary conditions of a single PDE into the system matrix
    /// 'mat' and the right hand side 'rhs'
    virtual void assemblePDEBoundarySingle(Paso_SystemMatrix* mat,
            escript::Data& rhs, const escript::Data& d,
            const escript::Data& y) const;

    /// assembles a single PDE with reduced order into the system matrix 'mat'
    /// and the right hand side 'rhs'
    virtual void assemblePDESingleReduced(Paso_SystemMatrix* mat,
            escript::Data& rhs, const escript::Data& A, const escript::Data& B,
            const escript::Data& C, const escript::Data& D,
            const escript::Data& X, const escript::Data& Y) const;

    /// assembles boundary conditions of a single PDE with reduced order into
    /// the system matrix 'mat' and the right hand side 'rhs'
    virtual void assemblePDEBoundarySingleReduced(Paso_SystemMatrix* mat,
            escript::Data& rhs, const escript::Data& d,
            const escript::Data& y) const;

    /// assembles a system of PDEs into the system matrix 'mat' and the right
    /// hand side 'rhs'
    virtual void assemblePDESystem(Paso_SystemMatrix* mat, escript::Data& rhs,
            const escript::Data& A, const escript::Data& B,
            const escript::Data& C, const escript::Data& D,
            const escript::Data& X, const escript::Data& Y) const;

    /// assembles boundary conditions of a system of PDEs into the system
    /// matrix 'mat' and the right hand side 'rhs'
    virtual void assemblePDEBoundarySystem(Paso_SystemMatrix* mat,
            escript::Data& rhs, const escript::Data& d,
            const escript::Data& y) const;

    /// assembles a system of PDEs with reduced order into the system matrix
    /// 'mat' and the right hand side 'rhs'
    virtual void assemblePDESystemReduced(Paso_SystemMatrix* mat,
            escript::Data& rhs, const escript::Data& A, const escript::Data& B,
            const escript::Data& C, const escript::Data& D,
            const escript::Data& X, const escript::Data& Y) const;

    /// assembles boundary conditions of a system of PDEs with reduced order
    /// into the system matrix 'mat' and the right hand side 'rhs'
    virtual void assemblePDEBoundarySystemReduced(Paso_SystemMatrix* mat,
            escript::Data& rhs, const escript::Data& d,
            const escript::Data& y) const;

    /// returns the Paso system matrix pattern
    virtual Paso_SystemMatrixPattern* getPattern(bool reducedRowOrder,
            bool reducedColOrder) const;

    /// interpolates data on nodes in 'in' onto (reduced) elements in 'out'
    virtual void interpolateNodesOnElements(escript::Data& out,
                                       escript::Data& in, bool reduced) const;

    /// interpolates data on nodes in 'in' onto (reduced) face elements in 'out'
    virtual void interpolateNodesOnFaces(escript::Data& out, escript::Data& in,
                                         bool reduced) const;

    /// converts data on nodes in 'in' to degrees of freedom in 'out'
    virtual void nodesToDOF(escript::Data& out, escript::Data& in) const;

    /// converts data on degrees of freedom in 'in' to nodes in 'out'
    virtual void dofToNodes(escript::Data& out, escript::Data& in) const;

    /// copies data in 'in' to 'out' (both must be on same function space)
    void copyData(escript::Data& out, escript::Data& in) const;

    /// averages data in 'in' to 'out' (from non-reduced to reduced fs)
    void averageData(escript::Data& out, escript::Data& in) const;

    // this is const because setTags is const
    void updateTagsInUse(int fsType) const;

    /// allocates and returns a Paso pattern structure
    Paso_Pattern* createPasoPattern(const IndexVector& ptr,
            const IndexVector& index, const dim_t M, const dim_t N) const;

    /// creates the pattern for the main block of the system matrix
    Paso_Pattern* createMainPattern() const;

    /// creates the pattern for the column and row couple blocks of the system
    /// matrix. colIndices[i] contains all IDs of DOFs that are connected with
    /// DOF i but remote and 'N' is the total number of remote components
    void createCouplePatterns(const std::vector<IndexVector>& colIndices,
                              const dim_t N, Paso_Pattern** colPattern,
                              Paso_Pattern** rowPattern) const;

    void addToSystemMatrix(Paso_SystemMatrix* in, const IndexVector& nodes_Eq,
            dim_t num_Eq, const IndexVector& nodes_Sol, dim_t num_Sol,
            const std::vector<double>& array) const;

    dim_t m_numDim;
    StatusType m_status;
    Esys_MPIInfo *m_mpiInfo;
    TagMap m_tagMap;
    mutable IndexVector m_nodeTags, m_nodeTagsInUse;
    mutable IndexVector m_elementTags, m_elementTagsInUse;
    mutable IndexVector m_faceTags, m_faceTagsInUse;

private:
    /// calls the right PDE assembly routines after performing input checks
    virtual void assemblePDE(Paso_SystemMatrix* mat, escript::Data& rhs,
            const escript::Data& A, const escript::Data& B,
            const escript::Data& C, const escript::Data& D,
            const escript::Data& X, const escript::Data& Y) const;

    /// calls the right PDE boundary assembly routines after performing input
    /// checks
    virtual void assemblePDEBoundary(Paso_SystemMatrix* mat,
            escript::Data& rhs, const escript::Data& d,
            const escript::Data& y) const;
};

} // end of namespace ripley

#endif // __RIPLEY_DOMAIN_H__
