
/*****************************************************************************
*
* Copyright (c) 2003-2018 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/

#ifndef __DUDLEY_DOMAIN_H__
#define __DUDLEY_DOMAIN_H__

/****************************************************************************

   Dudley: Domain

   A mesh is built from nodes and elements which describe the domain, surface,
   and point sources (the latter are needed to establish links with other
   codes, in particular particle codes). The nodes are stored in a NodeFile
   and elements in ElementFiles. Dudley domains have three ElementFiles
   containing the elements, surface and point sources, respectively.
   Notice that the surface elements do not necessarily cover the entire
   surface of the domain.

   The element type is either Tri3 or Tet4 depending on dimensionality
   and also determines the type of surface elements to be used.

   The numbering of the nodes starts with 0.

   Important: it is assumed that every node appears in at least one element or
   surface element and that any node used in an element, surface element or as
   a point is specified in the NodeFile, see also resolveNodeIds.

   All nodes and elements are tagged. The tag allows to group nodes and
   elements. A typical application is to mark surface elements on a
   certain portion of the domain with the same tag. All these surface
   elements can then be assigned the same value e.g. for the pressure.

*****************************************************************************/

#include <dudley/Dudley.h>
#include <dudley/ElementFile.h>
#include <dudley/NodeFile.h>
#include <dudley/Util.h>

#include <escript/AbstractContinuousDomain.h>
#include <escript/FunctionSpace.h>
#include <escript/FunctionSpaceFactory.h>

#ifdef ESYS_HAVE_PASO
#include <paso/SystemMatrixPattern.h>
#endif
#ifdef ESYS_HAVE_TRILINOS
#include <trilinoswrap/types.h>
#endif

#include <map>
#include <string>
#include <vector>

namespace dudley {

typedef std::map<std::string, int> TagMap;

enum SystemMatrixType {
    SMT_PASO = 1<<8,
    SMT_TRILINOS = 1<<10,
    SMT_COMPLEX = 1<<16,
    SMT_UNROLL = 1<<17
};

/**
    \brief
    DudleyDomain implements the AbstractContinuousDomain interface for the
    Dudley library.
*/
class DudleyDomain : public escript::AbstractContinuousDomain
{
public:
    /**
     \brief
     recovers domain from a dump file
     \param filename the name of the file
    */
    static escript::Domain_ptr load(const std::string& filename);

    /**
     \brief
     reads a mesh from a fly file. For MPI parallel runs fans out the mesh
     to multiple processes.
     \param mpiInfo the MPI information structure
     \param fileName the name of the file
     \param optimize whether to optimize the node labels
    */
    static escript::Domain_ptr read(escript::JMPI mpiInfo,
                                    const std::string& filename, bool optimize);

    /**
     \brief
     reads a gmsh mesh file.
     \param mpiInfo the MPI information structure
     \param filename the name of the gmsh file
     \param numDim spatial dimensionality
     \param optimize whether to optimize the node labels 
    */
    static escript::Domain_ptr readGmsh(escript::JMPI mpiInfo,
                                        const std::string& filename, int numDim,
                                        bool optimize);

    /**
     \brief
     Creates a 2-dimensional rectangular domain.

     \param NE0 Input - number of elements in first dimension
     \param NE1 Input - number of elements in second dimension
     \param l0 Input - length of domain in first dimension (width)
     \param l1 Input - length of domain in second dimension (height)
     \param optimize Input - whether to optimize node/DOF labelling
     \param jmpi Input - Shared pointer to MPI Information to be used
    */
    static escript::Domain_ptr create2D(dim_t NE0, dim_t NE1, double l0,
                                        double l1, bool optimize,
                                        escript::JMPI jmpi);

    /**
     \brief
     Creates a 3-dimensional rectangular domain.

     \param NE0 Input - number of elements in first dimension
     \param NE1 Input - number of elements in second dimension
     \param NE2 Input - number of elements in third dimension
     \param l0 Input - length of domain in first dimension (width)
     \param l1 Input - length of domain in second dimension (height)
     \param l2 Input - length of domain in third dimension (depth)
     \param optimize Input - whether to optimize node/DOF labelling
     \param jmpi Input - Shared pointer to MPI Information to be used
    */
    static escript::Domain_ptr create3D(dim_t NE0, dim_t NE1, dim_t NE2,
                                        double l0, double l1, double l2,
                                        bool optimize, escript::JMPI jmpi);

    /**
     \brief
     Constructor for DudleyDomain

     \param name a descriptive name for the domain
     \param numDim dimensionality of the domain (2 or 3)
     \param jmpi shared pointer to MPI Information to be used
    */
    DudleyDomain(const std::string& name, int numDim, escript::JMPI jmpi);

    /**
     \brief
     Copy constructor.
    */
    DudleyDomain(const DudleyDomain& in);

    /**
     \brief
     Destructor for DudleyDomain
    */
    ~DudleyDomain();

    /**
     \brief
     returns a pointer to this domain's node file
    */
    NodeFile* getNodes() const { return m_nodes; }

    /**
     \brief
     replaces the element file by `elements`
    */
    void setElements(ElementFile* elements);

    /**
     \brief
     returns a pointer to this domain's element file
    */
    ElementFile* getElements() const { return m_elements; }

    /**
     \brief
     replaces the face element file by `elements`
    */
    void setFaceElements(ElementFile* elements);

    /**
     \brief
     returns a pointer to this domain's face element file
    */
    ElementFile* getFaceElements() const { return m_faceElements; }

    /**
     \brief
     replaces the point element file by `elements`
    */
    void setPoints(ElementFile* elements);

    /**
     \brief
     returns a pointer to this domain's point (nodal) element file
    */
    ElementFile* getPoints() const { return m_points; }

    /**
     \brief
     returns a reference to the MPI information wrapper for this domain
    */
    virtual escript::JMPI getMPI() const { return m_mpiInfo; }

    /**
     \brief
     returns the number of processors used for this domain
    */
    virtual int getMPISize() const { return m_mpiInfo->size; }

    /**
     \brief
     returns the number MPI rank of this processor
    */
    virtual int getMPIRank() const { return m_mpiInfo->rank; }

    /**
     \brief
     If compiled for MPI then execute an MPI_Barrier, else do nothing
    */
    virtual void MPIBarrier() const;

    /**
     \brief
     returns true if on MPI processor 0, else false
    */
    virtual bool onMasterProcessor() const { return getMPIRank() == 0; }

    MPI_Comm getMPIComm() const { return m_mpiInfo->comm; }

    /**
     \brief
     writes the current mesh to a file with the given name in the fly file
     format.
     \param fileName Input - The name of the file to write to.
    */
    void write(const std::string& fileName) const;

    /**
     \brief
     \param full whether to include coordinate values and id's
    */
    void Print_Mesh_Info(bool full=false) const;

    /**
     \brief
     dumps the mesh to a file with the given name.
     \param fileName Input - The name of the file
    */
    void dump(const std::string& fileName) const;

    /**
     \brief
     Return the tag key for the given sample number.
     \param functionSpaceType Input - The function space type.
     \param sampleNo Input - The sample number.
    */
    int getTagFromSampleNo(int functionSpaceType, index_t sampleNo) const;

    /**
     \brief
     Return the reference number of  the given sample number.
     \param functionSpaceType Input - The function space type.
    */
    const index_t* borrowSampleReferenceIDs(int functionSpaceType) const;

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
    virtual int getDim() const { return m_nodes->numDim; }

    /**
     \brief
      Returns a status indicator of the domain. The status identifier should be unique over
      the live time if the object but may be updated if changes to the domain happen, e.g.
      modifications to its geometry.
    */
    virtual StatusType getStatus() const;

    /**
     \brief
     Return the number of data points summed across all MPI processes
    */
    virtual dim_t getNumDataPointsGlobal() const;

    /**
     \brief
     Return the number of data points per sample, and the number of samples as a pair.
     \param functionSpaceCode Input -
    */
    virtual std::pair<int,dim_t> getDataShape(int functionSpaceCode) const;

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
     Returns true if name is a defined tag name.
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
    virtual void interpolateOnDomain(escript::Data& target,
                                     const escript::Data& source) const;

    virtual bool probeInterpolationOnDomain(int functionSpaceType_source,
                                           int functionSpaceType_target) const;

    virtual signed char preferredInterpolationOnDomain(int functionSpaceType_source, int functionSpaceType_target) const;

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
    virtual bool probeInterpolationAcross(int functionSpaceType_source,
                                  const escript::AbstractDomain& targetDomain,
                                  int functionSpaceType_target) const;

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
    virtual void setToGradient(escript::Data& grad, const escript::Data& arg) const;

    /**
     \brief
     copies the integrals of the function defined by arg into integrals.
     arg has to be defined on this.
    */
    virtual void setToIntegrals(std::vector<escript::DataTypes::real_t>& integrals,
                                const escript::Data& arg) const;
    virtual void setToIntegrals(std::vector<escript::DataTypes::cplx_t>& integrals,
                                const escript::Data& arg) const;

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

    virtual bool ownSample(int fsCode, index_t id) const;

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
                                      bool useHRZ) const;

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
                     escript::AbstractTransportProblem& tp,
                     escript::Data& source, const escript::Data& M,
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
      creates a TransportProblem
    */
    escript::ATP_ptr newTransportProblem(int blocksize,
                                   const escript::FunctionSpace& functionspace,
                                   int type) const;

    /**
     \brief returns locations in the FEM nodes
    */
    virtual escript::Data getX() const;

    /**
     \brief returns boundary normals at the quadrature point on the face
            elements
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
     \brief assigns new tag newTag to all samples of functionspace with a
            positive value of mask for any its sample point.
    */
    virtual void setTags(int functionSpaceType, int newTag,
                         const escript::Data& mask) const;

    /**
      \brief
       returns the number of tags in use and a pointer to an array with the
       number of tags in use
    */
    virtual int getNumberOfTagsInUse(int functionSpaceCode) const;

    virtual const int* borrowListOfTagsInUse(int functionSpaceCode) const;

    /**
     \brief Checks if this domain allows tags for the specified
            functionSpace code.
    */
    virtual bool canTag(int functionSpaceCode) const;

    /**
     \brief returns the approximation order used for a function space functionSpaceCode
    */
    virtual int getApproximationOrder(int functionSpaceCode) const;

    virtual bool supportsContactElements() const { return false; }

    virtual escript::Data randomFill(const escript::DataTypes::ShapeType& shape,
                                const escript::FunctionSpace& what, long seed,
                                const boost::python::tuple& filter) const;

    void createMappings(const std::vector<index_t>& dofDistribution,
                        const std::vector<index_t>& nodeDistribution);

    /// assigns new node reference numbers to all element files.
    /// If k is the old node, the new node is newNode[k-offset].
    void relabelElementNodes(const index_t* newNode, index_t offset);

    template<typename Scalar>
    void setToIntegralsWorker(std::vector<Scalar>& integrals,
                              const escript::Data& arg) const;

#ifdef ESYS_HAVE_PASO
    /// returns a reference to the paso matrix pattern
    paso::SystemMatrixPattern_ptr getPasoPattern() const;
#endif

#ifdef ESYS_HAVE_TRILINOS
    /// returns a Trilinos CRS graph suitable to build a sparse matrix.
    esys_trilinos::const_TrilinosGraph_ptr getTrilinosGraph() const {
        return m_nodes->getTrilinosGraph();
    }
#endif

private:
    void prepare(bool optimize);

    /// Initially the element nodes refer to the numbering defined by the
    /// global id assigned to the nodes in the NodeFile. It is also not ensured
    /// that all nodes referred by an element are actually available on the
    /// process. At the output, a local node labeling is used and all nodes are
    /// available. In particular the numbering of the element nodes is between
    /// 0 and Nodes->numNodes.
    /// The function does not create a distribution of the degrees of freedom.
    void resolveNodeIds();

#ifdef ESYS_HAVE_PASO
    paso::SystemMatrixPattern_ptr makePasoPattern() const;
#endif

    void createColoring(const index_t* dofMap);
    void distributeByRankOfDOF(const IndexVector& distribution);
    void markNodes(std::vector<short>& mask, index_t offset) const;
    void optimizeDOFDistribution(IndexVector& distribution);
    void optimizeDOFLabeling(const IndexVector& distribution);
    void optimizeElementOrdering();
    void updateTagList();
    void printElementInfo(const ElementFile* e, const std::string& title,
                          const std::string& defaultType, bool full) const;

    void writeElementInfo(std::ostream& stream, const ElementFile* e,
                          const std::string& defaultType) const;

    /// MPI information
    escript::JMPI m_mpiInfo;
    /// domain description
    std::string m_name;
    /// the table of the nodes
    NodeFile* m_nodes;
    /// the table of the elements
    ElementFile* m_elements;
    /// the table of face elements
    ElementFile* m_faceElements;
    /// the table of points (treated as elements of dimension 0)
    ElementFile* m_points;
    /// the tag map mapping names to tag keys
    TagMap m_tagMap;
#ifdef ESYS_HAVE_PASO
    // pointer to the sparse matrix pattern
    mutable paso::SystemMatrixPattern_ptr pasoPattern;
#endif

    static FunctionSpaceNamesMapType m_functionSpaceTypeNames;
};

} // end of namespace

#endif // __DUDLEY_DOMAIN_H__

