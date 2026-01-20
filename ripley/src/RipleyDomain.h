
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* https://github.com/LutzGross/esys-escript.github.io
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/

#ifndef __RIPLEY_DOMAIN_H__
#define __RIPLEY_DOMAIN_H__

#include <ripley/Ripley.h>
#include <ripley/AbstractAssembler.h>
#include <ripley/RipleyException.h>

#include <escript/AbstractContinuousDomain.h>
#include <escript/Data.h>
#include <escript/FunctionSpace.h>

#ifdef ESYS_HAVE_PASO
#include <paso/Coupler.h>
#include <paso/SystemMatrix.h>
#endif

#ifdef ESYS_HAVE_TRILINOS
#include <trilinoswrap/types.h>
#endif

#include <boost/python/list.hpp>
#include <boost/python/tuple.hpp>

namespace ripley {

enum assembler_t {
    DEFAULT_ASSEMBLER,
    WAVE_ASSEMBLER,
    LAME_ASSEMBLER
};

enum SystemMatrixType {
    SMT_PASO = 1<<8,
    SMT_CUSP = 1<<9,
    SMT_TRILINOS = 1<<10,
    SMT_SYMMETRIC = 1<<15,
    SMT_COMPLEX = 1<<16,
    SMT_UNROLL = 1<<17
};

enum DecompositionPolicy {
    DECOMP_ADD_ELEMENTS,
    DECOMP_EXPAND,
    DECOMP_STRICT
};

/**
   \brief
   Structure that wraps parameters for the grid reading routines.
*/
struct ReaderParameters
{
    /// the (global) offset into the data object to start writing into
    std::vector<dim_t> first;
    /// the number of values to read from file
    std::vector<dim_t> numValues;
    /// how often to write each value from the file into the data object
    /// (e.g. to supersample)
    std::vector<int> multiplier;
    /// if non-zero, values are written from last index to first index
    std::vector<int> reverse;
    /// byte order in the file (used by binary reader only)
    int byteOrder;
    /// data type in the file (used by binary reader only)
    int dataType;
};

/**
   \brief
   A struct to contain a dirac point's information.
*/
struct DiracPoint
{
    dim_t node;
    int tag;
};

/**
   \brief
   RipleyDomain extends the AbstractContinuousDomain interface
   for the Ripley library and is the base class for Rectangle and Brick.
*/

class RIPLEY_DLL_API RipleyDomain : public escript::AbstractContinuousDomain
{
public:
    /**
       \brief
       Constructor with number of dimensions. Allocates MPI info structure.
    */
    RipleyDomain(dim_t dim);

    /**
       \brief
       Constructor with number of dimensions and custom MPI communicator.
    */
    RipleyDomain(dim_t dim, escript::JMPI jmpi);

    /**
       \brief
       Destructor
    */
    ~RipleyDomain();

    static void setDecompositionPolicy(DecompositionPolicy value);
    static DecompositionPolicy getDecompositionPolicy();

    /**
       \brief
       returns true if the argument is a valid function space type for this
       domain
    */
    virtual bool isValidFunctionSpaceType(int fsType) const;

    /**
       \brief
       returns a description for the given function space type code
    */
    virtual std::string functionSpaceTypeAsString(int fsType) const;

    /**
       \brief
       returns the number of spatial dimensions of the domain
    */
    virtual int getDim() const { return m_numDim; }

    /**
       \brief equality operator
    */
    virtual bool operator==(const escript::AbstractDomain& other) const;

    /**
       \brief inequality operator
    */
    virtual bool operator!=(const escript::AbstractDomain& other) const {
        return !(operator==(other));
    }

    /**
       \brief
       returns the number of data points per sample, and the number of samples
       as a pair.
       \param fsType The function space type
    */
    virtual std::pair<int,dim_t> getDataShape(int fsType) const;

    /**
       \brief
       returns the tag key for the given sample number
       \param fsType The function space type
       \param sampleNo The sample number
    */
    int getTagFromSampleNo(int fsType, dim_t sampleNo) const;

    /**
       \brief
       sets a map from a clear tag name to a tag key
       \param name tag name
       \param tag tag key
    */
    virtual void setTagMap(const std::string& name, int tag) {
        m_tagMap[name] = tag;
    }

    /**
       \brief
       returns the tag key for tag name
       \param name tag name
    */
    virtual int getTag(const std::string& name) const {
        if (m_tagMap.find(name) != m_tagMap.end()) {
            return m_tagMap.find(name)->second;
        } else {
            throw escript::ValueError("getTag: invalid tag name");
        }
    }

    /**
       \brief
       returns true if name is a defined tag name
       \param name tag name to be checked
    */
    virtual bool isValidTagName(const std::string& name) const {
        return (m_tagMap.find(name)!=m_tagMap.end());
    }

    /**
       \brief
       returns all tag names in a single string separated by commas
    */
    virtual std::string showTagNames() const;

    /**
       \brief
       assigns new location to the domain.
       \note This is not supported in Ripley
    */
    virtual void setNewX(const escript::Data& arg);

    /**
       \brief
       interpolates data given on source onto target where source and target
       have to be given on the same domain
    */
    virtual void interpolateOnDomain(escript::Data& target, const escript::Data& source) const;

    /**
       \brief
       returns true if data on fsType_source can be interpolated onto
       fsType_target, false otherwise
    */
    virtual bool probeInterpolationOnDomain(int fsType_source, int fsType_target) const;

    /**
       \brief Preferred direction of interpolation.
       If you really need to test for a particular direction, then use
       probeInterpolation.
       \return 0 for not possible, 1 for possible and preferred, -1 other
               direction preferred (does not mean this direction is possible)
    */
    virtual signed char preferredInterpolationOnDomain(int fsType_source,
                                                       int fsType_target) const;

    /**
       \brief
       given a vector of FunctionSpace type codes, passes back a code which all
       can be interpolated to
       \return true if result is valid, false if not
    */
    bool
    commonFunctionSpace(const std::vector<int>& fs, int& resultcode) const;

    /**
       \brief
       interpolates data given on source onto target where source and target
       are given on different domains
    */
    virtual void interpolateAcross(escript::Data& target,
                                   const escript::Data& source) const;

    /**
       \brief
       determines whether interpolation from source to target is possible
    */
    virtual bool probeInterpolationAcross(int, const escript::AbstractDomain&, int) const;

    /**
       \brief
       returns locations in the FEM nodes
    */
    virtual escript::Data getX() const;

#ifdef ESYS_HAVE_BOOST_NUMPY
    /**
       \brief
       returns locations in the FEM nodes as a numpy ndarray
    */
    virtual boost::python::numpy::ndarray getNumpyX() const;

    /**
     \brief returns connectivity information as a numpy ndarray
    */
    virtual boost::python::numpy::ndarray getConnectivityInfo() const;
#endif

    /**
       \brief
       returns boundary normals at the quadrature point on the face elements
    */
    virtual escript::Data getNormal() const;

    /**
       \brief returns the element size
    */
    virtual escript::Data getSize() const;

    /**
       \brief
       copies the location of data points into arg. The domain of arg has to
       match this domain.
    */
    virtual void setToX(escript::Data& arg) const;

    /**
       \brief
       copies the gradient of 'in' into 'out'. The actual function space to be
       considered for the gradient is defined by 'in'. Both arguments have to
       be defined on this domain.
    */
    virtual void setToGradient(escript::Data& out, const escript::Data& in) const;

    /**
       \brief
       assigns new tag newTag to all samples of given function space with a
       positive value of mask for any of its sample points
    */
    virtual void setTags(int fsType, int newTag, const escript::Data& mask) const;

    /**
       \brief
       returns true if data on this domain and given function space type has
       to be considered as cell centered data
    */
    virtual bool isCellOriented(int fsType) const;

    /**
       \brief
       returns a status indicator of the domain. The status identifier should
       be unique over the lifetime of the object but may be updated if changes
       to the domain happen, e.g. modifications to its geometry.
    */
    virtual StatusType getStatus() const { return m_status; }

    /**
       \brief
       returns the number of tags in use for a function space type
    */
    virtual int getNumberOfTagsInUse(int fsType) const;

    /**
       \brief
       returns a pointer to the list of tags in use for a function space type
    */
    virtual const int* borrowListOfTagsInUse(int fsType) const;

    /**
       \brief
       checks if this domain allows tags for the specified function space type
    */
    virtual bool canTag(int fsType) const;

    /**
       \brief
       returns the approximation order used for a function space
    */
    virtual int getApproximationOrder(int fsType) const { return 1; }

    /**
       \brief
       returns true if this domain supports contact elements, false otherwise
    */
    virtual bool supportsContactElements() const { return false; }

    /**
       \brief
       returns a continuous FunctionSpace code
    */
    virtual int getContinuousFunctionCode() const { return Nodes; }

    /**
       \brief
       returns a continuous on reduced order nodes FunctionSpace code
    */
    virtual int getReducedContinuousFunctionCode() const { return ReducedNodes; }

    /**
       \brief
       returns a function FunctionSpace code
    */
    virtual int getFunctionCode() const { return Elements; }

    /**
       \brief
       returns a function with reduced integration order FunctionSpace code
    */
    virtual int getReducedFunctionCode() const { return ReducedElements; }

    /**
       \brief
       returns a function on boundary FunctionSpace code
    */
    virtual int getFunctionOnBoundaryCode() const { return FaceElements; }

    /**
       \brief
       returns a function on boundary with reduced integration order
       FunctionSpace code
    */
    virtual int getReducedFunctionOnBoundaryCode() const { return ReducedFaceElements; }

    /**
       \brief
       return a FunctionOnContactZero code
    */
    virtual int getFunctionOnContactZeroCode() const {
        throw escript::NotImplementedError("Ripley does not support contact elements");
    }

    /**
       \brief
       returns a FunctionOnContactZero code with reduced integration order
    */
    virtual int getReducedFunctionOnContactZeroCode() const {
        throw escript::NotImplementedError("Ripley does not support contact elements");
    }

    /**
       \brief
       returns a FunctionOnContactOne code
    */
    virtual int getFunctionOnContactOneCode() const {
        throw escript::NotImplementedError("Ripley does not support contact elements");
    }

    /**
       \brief
       returns a FunctionOnContactOne code with reduced integration order
    */
    virtual int getReducedFunctionOnContactOneCode() const {
        throw escript::NotImplementedError("Ripley does not support contact elements");
    }

    /**
       \brief
       returns a Solution FunctionSpace code
    */
    virtual int getSolutionCode() const { return DegreesOfFreedom; }

    /**
       \brief
       returns a ReducedSolution FunctionSpace code
    */
    virtual int getReducedSolutionCode() const { return ReducedDegreesOfFreedom; }

    /**
       \brief
       returns a DiracDeltaFunctions FunctionSpace code
    */
    virtual int getDiracDeltaFunctionsCode() const { return Points; }

    /**
       \brief
       returns the identifier of the matrix type to be used for the global
       stiffness matrix when a particular solver, package, preconditioner,
       and symmetric matrix is used
       \param options a python object containing the solver, package,
                preconditioner and symmetry
    */
    virtual int getSystemMatrixTypeId(const boost::python::object& options) const;

    /**
       \brief
       returns the identifier of the transport problem type to be used when a
       particular solver, preconditioner, package and symmetric matrix is used
       \param solver
       \param preconditioner
       \param package
       \param symmetry
    */
    virtual int getTransportTypeId(int solver, int preconditioner, int package,
                                   bool symmetry) const;

    /**
       \brief
       copies the integrals of the function defined by arg into integrals.
       arg has to be defined on this domain.
    */
    virtual void setToIntegrals(std::vector<real_t>& integrals,
                                const escript::Data& arg) const;
    virtual void setToIntegrals(std::vector<cplx_t>& integrals,
                                const escript::Data& arg) const;


    /**
       \brief
       adds a PDE onto the stiffness matrix mat and rhs, used for custom
       solvers with varying arguments counts and so on
    */
    virtual void addToSystem(escript::AbstractSystemMatrix& mat,
                             escript::Data& rhs, const DataMap& data,
                             Assembler_ptr assembler) const;

    /**
       \brief
       a wrapper for addToSystem that allows calling from Python
    */
    virtual void addToSystemFromPython(escript::AbstractSystemMatrix& mat,
            escript::Data& rhs, const boost::python::list& data,
            Assembler_ptr assembler) const;

    /**
       \brief
       adds a PDE onto rhs, used for custom
       solvers with varying arguments counts and so on
    */
    virtual void addToRHS(escript::Data& rhs, const DataMap& data,
                          Assembler_ptr assembler) const;

    /**
       \brief
       a wrapper for addToRHS that allows calling from Python
    */
    virtual void addToRHSFromPython(escript::Data& rhs,
                                    const boost::python::list& data,
                                    Assembler_ptr assembler) const;

    /**
       \brief
       adds a PDE onto a transport problem
    */
    virtual void addPDEToTransportProblem(escript::AbstractTransportProblem& tp,
                                    escript::Data& source, const DataMap& data,
                                    Assembler_ptr assembler) const;
    /**
      \brief Do not use. Throws exception.
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
       adds a PDE onto a transport problem
    */
    void addPDEToTransportProblemFromPython(
                        escript::AbstractTransportProblem& tp,
                        escript::Data& source, const boost::python::list& data,
                        Assembler_ptr assembler) const;

    /**
       \brief
       creates a stiffness matrix and initializes it with zeros
    */
    virtual escript::ASM_ptr newSystemMatrix(int row_blocksize,
            const escript::FunctionSpace& row_functionspace,
            int column_blocksize,
            const escript::FunctionSpace& column_functionspace, int type) const;

    /**
     \brief
      creates a transport problem
    */
    virtual escript::ATP_ptr newTransportProblem(int blocksize,
                                  const escript::FunctionSpace& functionspace,
                                  int type) const;

    /**
       \brief
       writes information about the mesh to standard output
       \param full whether to print additional data
    */
    virtual void Print_Mesh_Info(bool full=false) const;


    /************************************************************************/

    /**
       \brief
       writes the current mesh to a file with the given name
       \param filename The name of the file to write to
    */
    virtual void write(const std::string& filename) const = 0;

    /**
       \brief
       returns a description for this domain
    */
    virtual std::string getDescription() const = 0;

    /**
       \brief
       dumps the mesh to a file with the given name
       \param filename The name of the output file
    */
    void dump(const std::string& filename) const = 0;

    /**
       \brief
       returns the array of reference numbers for a function space type
       \param fsType The function space type
    */
    const dim_t* borrowSampleReferenceIDs(int fsType) const = 0;

    /**
       \brief
       copies the surface normals at data points into out. The actual function
       space to be considered is defined by out. out has to be defined on this
       domain.
    */
    virtual void setToNormal(escript::Data& out) const = 0;

    /**
       \brief
       copies the size of samples into out. The actual function space to be
       considered is defined by out. out has to be defined on this domain.
    */
    virtual void setToSize(escript::Data& out) const = 0;

    /**
       \brief
       reads grid data from a netCDF file into a Data object
    */
    virtual void readNcGrid(escript::Data& out, std::string filename,
            std::string varname, const ReaderParameters& params) const = 0;

    /**
       \brief
       reads grid data from a raw binary file into a Data object
    */
    virtual void readBinaryGrid(escript::Data& out, std::string filename,
                                const ReaderParameters& params) const = 0;

    /**
       \brief
       reads grid data from a compressed raw binary file into a Data object
    */
    virtual void readBinaryGridFromZipped(escript::Data& out,
               std::string filename, const ReaderParameters& params) const = 0;

    /**
       \brief
       writes a Data object to a file in raw binary format
    */
    virtual void writeBinaryGrid(const escript::Data& in, std::string filename,
                                 int byteOrder, int dataType) const = 0;

    /**
       \brief
       returns true if this rank owns the sample id on given function space
    */
    virtual bool ownSample(int fsType, index_t id) const = 0;

    /**
       \brief
       returns the number of data points summed across all MPI processes
    */
    virtual dim_t getNumDataPointsGlobal() const = 0;

    /**
       \brief
       returns the number of nodes per MPI rank in each dimension
    */
    virtual const dim_t* getNumNodesPerDim() const = 0;

    /**
       \brief
       returns the number of elements per MPI rank in each dimension
    */
    virtual const dim_t* getNumElementsPerDim() const = 0;

    /**
       \brief
       returns the number of face elements in the order
       (left,right,bottom,top,[front,back]) on current MPI rank
    */
    virtual const dim_t* getNumFacesPerBoundary() const = 0;

    /**
       \brief
       returns the node distribution vector
    */
    virtual IndexVector getNodeDistribution() const = 0;

    /**
       \brief
       returns the number of spatial subdivisions in each dimension
    */
    virtual const int* getNumSubdivisionsPerDim() const = 0;

    /**
       \brief
       returns the index'th coordinate value in given dimension for this rank
    */
    virtual double getLocalCoordinate(index_t index, int dim) const = 0;

    /**
       \brief
       returns the tuple (origin, spacing, number_of_elements)
    */
    virtual boost::python::tuple getGridParameters() const = 0;

    /**
       \brief
       returns the lengths of the domain
    */
    virtual const double *getLength() const = 0;

    /**
       \brief
       returns the lengths of an element
    */
    virtual const double *getElementLength() const = 0;

    /**
       \brief
       returns a vector of rank numbers where vec[i]=n means that rank n
       'owns' element/face element i.
    */
    virtual RankVector getOwnerVector(int fsType) const = 0;

    /**
       \brief
       returns true if this domain can handle the specified tuple of filter
       options.
    */
    virtual bool supportsFilter(const boost::python::tuple& t) const;

    /**
       \brief
    */
    virtual Assembler_ptr createAssembler(std::string type,
                                          const DataMap& options) const {
        throw escript::NotImplementedError("Domain does not support custom assemblers");
    }

    /**
       \brief
    */
    Assembler_ptr createAssemblerFromPython(std::string type,
                                     const boost::python::list& options) const;

protected:
    int m_numDim;
    StatusType m_status;
    TagMap m_tagMap;
    mutable std::vector<int> m_nodeTags, m_nodeTagsInUse;
    mutable std::vector<int> m_elementTags, m_elementTagsInUse;
    mutable std::vector<int> m_faceTags, m_faceTagsInUse;
    // mutable std::vector<int> m_pointsTagsInUse;
    std::vector<DiracPoint> m_diracPoints;
    IndexVector m_diracPointNodeIDs; //for borrowSampleID
    assembler_t assembler_type;

    /// copies data in 'in' to 'out' (both must be on same function space)
    template<typename Scalar>
    void copyData(escript::Data& out, const escript::Data& in) const;

    /// averages data in 'in' to 'out' (from non-reduced to reduced fs)
    template<typename Scalar>
    void averageData(escript::Data& out, const escript::Data& in) const;

    /// copies data in 'in' to 'out' (from reduced to non-reduced fs)
    template<typename Scalar>
    void multiplyData(escript::Data& out, const escript::Data& in) const;

    // this is const because setTags is const
    void updateTagsInUse(int fsType) const;

#ifdef ESYS_HAVE_PASO
    /// creates a Paso connector
    void createPasoConnector(const RankVector& neighbour,
                             const IndexVector& offsetInSharedSend,
                             const IndexVector& offsetInSharedRecv,
                             const IndexVector& sendShared,
                             const IndexVector& recvShared);

    /// returns a Paso connector required for data transfer and distributed
    /// system matrices
    paso::Connector_ptr getPasoConnector() const { return m_connector; }

    /// allocates and returns a Paso pattern structure
    paso::Pattern_ptr createPasoPattern(const std::vector<IndexVector>& indices,
                                        dim_t N) const;
#endif

#ifdef ESYS_HAVE_TRILINOS
    /// creates and returns a Trilinos CRS graph suitable to build a sparse
    /// matrix
    esys_trilinos::const_TrilinosGraph_ptr createTrilinosGraph(
            const IndexVector& myRows, const IndexVector& myColumns) const;
#endif

    template<typename Scalar>
    void addToSystemMatrix(escript::AbstractSystemMatrix* mat,
                           const IndexVector& nodes, dim_t numEq,
                           const std::vector<Scalar>& array) const;

    void addPoints(const std::vector<double>& coords,
                   const std::vector<int>& tags);

    /***********************************************************************/

    /// returns the number of nodes per MPI rank
    virtual dim_t getNumNodes() const = 0;

    /// returns the number of elements per MPI rank
    virtual dim_t getNumElements() const = 0;

    /// returns the number of degrees of freedom per MPI rank
    virtual dim_t getNumDOF() const = 0;

    /// returns the number of face elements on current MPI rank
    virtual dim_t getNumFaceElements() const = 0;

    /// returns the indices of the occupied matrix diagonals
    virtual IndexVector getDiagonalIndices(bool upperOnly) const = 0;

    /// populates the data object 'arg' with the node coordinates
    virtual void assembleCoordinates(escript::Data& arg) const = 0;

    /// computes the gradient of 'in' and puts the result in 'out'
    virtual void assembleGradient(escript::Data& out, const escript::Data& in) const = 0;

    /// copies the integrals of the function defined by 'arg' into 'integrals'
    virtual void assembleIntegrate(std::vector<real_t>& integrals, const escript::Data& arg) const = 0;
    virtual void assembleIntegrate(std::vector<cplx_t>& integrals, const escript::Data& arg) const = 0;

#ifdef ESYS_HAVE_TRILINOS
    /// returns the Trilinos matrix graph
    virtual esys_trilinos::const_TrilinosGraph_ptr getTrilinosGraph() const = 0;
#endif

    /// returns occupied matrix column indices for all matrix rows
    virtual std::vector<IndexVector> getConnections(bool includeShared) const = 0;

#ifdef ESYS_HAVE_PASO
    /// returns the Paso system matrix pattern
    virtual paso::SystemMatrixPattern_ptr getPasoMatrixPattern(
                                              bool reducedRowOrder,
                                              bool reducedColOrder) const = 0;
#endif

    /// interpolates data on nodes in 'in' onto (reduced) elements in 'out'
    virtual void interpolateNodesOnElements(escript::Data& out,
                                            const escript::Data& in,
                                            bool reduced) const = 0;

    /// interpolates data on nodes in 'in' onto (reduced) face elements in 'out'
    virtual void interpolateNodesOnFaces(escript::Data& out,
                                         const escript::Data& in,
                                         bool reduced) const = 0;

    /// converts data on nodes in 'in' to degrees of freedom in 'out'
    virtual void nodesToDOF(escript::Data& out, const escript::Data& in) const = 0;

    /// converts data on degrees of freedom in 'in' to nodes in 'out'
    template<typename Scalar>
    void dofToNodes(escript::Data& out, const escript::Data& in) const;

    virtual dim_t getDofOfNode(dim_t node) const = 0;

private:
    static DecompositionPolicy m_decompPolicy;

#ifdef ESYS_HAVE_PASO
    // Paso connector used by the system matrix and to interpolate DOF to
    // nodes
    paso::Connector_ptr m_connector;

    /// paso version of adding element matrices to System Matrix
    template <typename T>
    void addToPasoMatrix(paso::SystemMatrix<T>* in, const IndexVector& nodes,
                         dim_t numEq, const std::vector<T>& array) const;
#endif

    /// calls the right PDE assembly routines after performing input checks
    void assemblePDE(escript::AbstractSystemMatrix* mat, escript::Data& rhs,
            const DataMap& coefs, Assembler_ptr assembler) const;

    /// calls the right PDE boundary assembly routines after performing input
    /// checks
    void assemblePDEBoundary(escript::AbstractSystemMatrix* mat,
                             escript::Data& rhs, const DataMap& coefs,
                             Assembler_ptr assembler) const;

    void assemblePDEDirac(escript::AbstractSystemMatrix* mat,
                          escript::Data& rhs, const DataMap& coefs,
                          Assembler_ptr assembler) const;

    template<typename Scalar>
    void setToIntegralsWorker(std::vector<Scalar>& integrals,
                              const escript::Data& arg) const;

    /// finds the node that the given point coordinates belong to
    virtual dim_t findNode(const double *coords) const = 0;
};

} // end of namespace ripley

#endif // __RIPLEY_DOMAIN_H__

