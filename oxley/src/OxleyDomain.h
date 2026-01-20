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

#ifndef __OXLEY_DOMAIN_H__
#define __OXLEY_DOMAIN_H__

#ifndef __OXLEY_EXCEPTION_H__

#include <oxley/Oxley.h>
#include <oxley/OxleyException.h>
#include <oxley/AbstractAssembler.h>
#include <oxley/domainhelpers.h>
#include <oxley/tictoc.h>

#include <escript/EsysMPI.h>
#include <escript/AbstractContinuousDomain.h>

#ifdef ESYS_HAVE_PASO
#include <paso/SystemMatrix.h>
#endif

#ifdef ESYS_HAVE_BOOST_NUMPY
#include <boost/python/numpy.hpp>
#endif

#include <boost/python/tuple.hpp>
#include <boost/python/to_python_converter.hpp>

#ifdef ESYS_HAVE_TRILINOS
#include <trilinoswrap/types.h>
#endif

#include <p4est/p4est.h>

namespace oxley {

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
  A struct used by boost to convert Teuchos RCP arrays to numpy arrays
*/
// struct convert_Teuchos_RCP_to_Python_tuple
// {
//     typedef Tpetra::CrsMatrix<cplx_t,esys_trilinos::LO,esys_trilinos::GO,esys_trilinos::NT> crs_matrix_type;
//     static PyObject* convert(Teuchos::RCP<crs_matrix_type> A)
//     {
//         // Set up the view
//         Teuchos::ArrayView<const esys_trilinos::GO> value_view;
//         Teuchos::ArrayView<const cplx_t> column_view; 
//         A->getGlobalRowView(0,value_view,column_view);
//         // auto value_view_1d  = Kokkos::subview(value_view, Kokkos::ALL(), 0);
//         // auto column_view_1d = Kokkos::subview(column_view, Kokkos::ALL(), 0);

//         // Write in the values
//         boost::python::list new_array;
//         for(int row = 0; row < value_view.getNumRows(); row++)
//         {
//             boost::python::tuple data = boost::python::make_tuple(row,column_view[row],value_view[row]);
//             new_array.append(data);
//         }

//         // return boost::python::make_tuple(true, boost::python::handle<>(new_array));
//         return boost::python::make_tuple(true, boost::python::handle<>(new_array));
//     }
// };

/*
This class is the parent of Oxley Rectangle and Brick
*/
class OxleyDomain : public escript::AbstractContinuousDomain
{
public:
    /**
       \brief
       Constructor
    */
    OxleyDomain(dim_t dim, int order, escript::JMPI jmpi = escript::JMPI());

    /**
       \brief
       Destructor
    */
    ~OxleyDomain();

    static void setDecompositionPolicy(DecompositionPolicy value);
    static DecompositionPolicy getDecompositionPolicy();

    /**
       \brief
       returns a description for this domain
    */
    virtual std::string getDescription() const;

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
    virtual int getTagFromSampleNo(int fsType, dim_t sampleNo) const;

    /**
       \brief
       sets a map from a clear tag name to a tag key
       \param name tag name
       \param tag tag key
    */
    virtual void setTagMap(const std::string& name, int tag)
    {
        m_tagMap[name] = tag;
    }

    /**
       \brief
       writes the current mesh to a file with the given name
       \param filename The name of the file to write to
    */
    virtual void write(const std::string& filename) const = 0;

    /**
       \brief
       dumps the mesh to a file with the given name
       \param filename The name of the output file
    */
    virtual void dump(const std::string& filename) const = 0;

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
       assigns new tag newTag to all samples of given function space with a
       positive value of mask for any of its sample points
    */
    virtual void setTags(int fsType, int newTag, const escript::Data& mask) const;

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
       returns true if name is a defined tag name
       \param name tag name to be checked
    */
    virtual bool isValidTagName(const std::string& name) const
    {
        return (m_tagMap.find(name)!=m_tagMap.end());
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
       returns all tag names in a single string separated by commas
    */
    virtual std::string showTagNames() const;


    // this is const because setTags is const
    virtual void updateTagsInUse(int fsType) const;

    /**
       \brief
       returns the array of reference numbers for a function space type
       \param fsType The function space type
    */
    const index_t* borrowSampleReferenceIDs(int fsType) const = 0;

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
                                   const escript::Data& source) const = 0;

    /**
       \brief
       determines whether interpolation from source to target is possible
    */
    virtual bool probeInterpolationAcross(int, const escript::AbstractDomain&, int) const;

    /**
       \brief
       returns locations of the nodes
    */
    virtual escript::Data getX() const;

    /**
       \brief
       returns boundary normals at the quadrature point on the face elements
    */
    virtual escript::Data getNormal()
    const;

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
    virtual void setToGradient(escript::Data& out,
            const escript::Data& in) const;

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
       returns the approximation order used for a function space
    */
    virtual int getApproximationOrder(int fsType) const { return 1; }

    /**
       \brief
       writes the mesh to a vtk file
    */
    virtual void writeToVTK(std::string filename, bool writeTagInfo) const;

    /**
       \brief
       writes the mesh to file
    */
    virtual void saveMesh(std::string filename) = 0;

    /**
       \brief
       writes the mesh to file
    */
    virtual void loadMesh(std::string filename) = 0;

    /**
       \brief
       sets the number of levels of refinement
    */
    virtual void setRefinementLevels(int refinementlevels) = 0;

    /**
       \brief
       refines the mesh using enum RefinementAlgorithm
    */
    virtual void refineMesh(std::string algorithm);

    /**
       \brief
       refines the mesh near a boundary
       \param maxRecursion Max levels of recursion
       \param algorithmname The algorithm to use
    */
    virtual void refineBoundary(std::string boundary, double dx);

    /**
       \brief
       refines the mesh within the interior of a region bound by 
       x0, x1, y0, y1
       \param x0 boundary of the region
       \param x1 boundary of the region
       \param y0 boundary of the region
       \param y1 boundary of the region
    */
    virtual void refineRegion(double x0, double x1, double y0, double y1);

    /**
       \brief
       refines the mesh within the interior of a region bound by 
       x0, x1, y0, y1
       \param x0 boundary of the region
       \param x1 boundary of the region
       \param y0 boundary of the region
       \param y1 boundary of the region
       \param z0 boundary of the region
       \param z1 boundary of the region
    */
    virtual void refineRegion(double x0, double x1, double y0, double y1, double z0, double z1);

    /**
       \brief
       refines the mesh around the point
       x0, y1
       \param x0 
       \param y1 
    */
    virtual void refinePoint(double x0, double y0);

    /**
       \brief
       refines the mesh around the point
       x0, y0
       \param x0 
       \param y0
       \param z0
    */
    virtual void refinePoint(double x0, double y0, double z0);

    /**
       \brief
       refines a circle on the mesh centered at (x0, y0) with radius r
       \param x0 
       \param y1 
       \param r
    */
    virtual void refineCircle(double x0, double y0, double r);

    /**
       \brief
       refines a sphere on the mesh centered at (x0, y0, z0) with radius r
       \param x0 
       \param y1 
       \param r
    */
    virtual void refineSphere(double x0, double y0, double z0, double r);

    /**
       \brief
       returns a Data object containing the coordinate information
    */
    virtual int getNumVertices() const = 0;

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
       copies the integrals of the function defined by arg into integrals.
       arg has to be defined on this domain.
    */
    virtual void setToIntegrals(std::vector<real_t>& integrals,
                                const escript::Data& arg) const;
    virtual void setToIntegrals(std::vector<cplx_t>& integrals,
                                const escript::Data& arg) const;

    /// copies the integrals of the function defined by 'arg' into 'integrals'
    virtual void assembleIntegrate(std::vector<real_t>& integrals, const escript::Data& arg) const = 0;
    virtual void assembleIntegrate(std::vector<cplx_t>& integrals, const escript::Data& arg) const = 0;


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
      finalises the matrix system
   */
   #ifdef ESYS_HAVE_TRILINOS
   void makeZ(bool complex);
   template<typename S> void makeZworker(const S half,Teuchos::RCP<Tpetra::CrsMatrix<S,esys_trilinos::LO,esys_trilinos::GO,esys_trilinos::NT>>& Z,
               Teuchos::RCP<const Tpetra::Map<>>,Teuchos::RCP<const Tpetra::Map<>>);
   void makeIZ(bool complex);
   template<typename S> void makeIZworker(Teuchos::RCP<Tpetra::CrsMatrix<S,esys_trilinos::LO,esys_trilinos::GO,esys_trilinos::NT>> &Z,
               Teuchos::RCP<const Tpetra::Map<>>,Teuchos::RCP<const Tpetra::Map<>>);

   bool z_needs_update=false;
   bool iz_needs_update=false;
   // Teuchos::RCP<Tpetra::CrsMatrix<cplx_t,esys_trilinos::LO,esys_trilinos::GO,esys_trilinos::NT>> * getZ(bool complex);
   // Teuchos::RCP<Tpetra::CrsMatrix<cplx_t,esys_trilinos::LO,esys_trilinos::GO,esys_trilinos::NT>> * getIZ(bool complex);
   
   void finaliseA(escript::AbstractSystemMatrix& mat, bool isComplex);
   template<typename S>
   void finaliseAworker(escript::AbstractSystemMatrix& mat, 
      Teuchos::RCP<Tpetra::CrsMatrix<S,esys_trilinos::LO,esys_trilinos::GO,esys_trilinos::NT>>& Z);
   
   
   escript::Data finaliseRhs(escript::Data& rhs);
   template<typename S> 
   void finaliseRhsworker(escript::Data& rhs, 
      Teuchos::RCP<Tpetra::CrsMatrix<S,esys_trilinos::LO,esys_trilinos::GO,esys_trilinos::NT>>& Z);

   void resetRhs(escript::Data& rhs) const;
   #endif //ESYS_HAVE_TRILINOS
   
   void saveFsType(escript::Data &rhs);
   int getOrigFsType() { return origFsTypecode;};
   int origFsTypecode=-1;
   


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
       return a FunctionOnContactZero code
    */
    virtual int getFunctionOnContactZeroCode() const {
        throw escript::NotImplementedError("Oxley does not support contact elements");
    }

    /**
       \brief
       returns a FunctionOnContactZero code with reduced integration order
    */
    virtual int getReducedFunctionOnContactZeroCode() const {
        throw escript::NotImplementedError("Oxley does not support contact elements");
    }

    /**
       \brief
       returns a FunctionOnContactOne code
    */
    virtual int getFunctionOnContactOneCode() const {
        throw escript::NotImplementedError("Oxley does not support contact elements");
    }

    /**
       \brief
       returns a FunctionOnContactOne code with reduced integration order
    */
    virtual int getReducedFunctionOnContactOneCode() const {
        throw escript::NotImplementedError("Oxley does not support contact elements");
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
       creates a stiffness matrix and initializes it with zeros
    */
    virtual escript::ASM_ptr newSystemMatrix(int row_blocksize,
            const escript::FunctionSpace& row_functionspace,
            int column_blocksize,
            const escript::FunctionSpace& column_functionspace, int type) const;

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
       returns the identifier of the matrix type to be used for the global
       stiffness matrix when a particular solver, package, preconditioner,
       and symmetric matrix is used
       \param options a python object containing the solver, package,
                preconditioner and symmetry
    */
    virtual int getSystemMatrixTypeId(const boost::python::object& options) const;

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

    /**
       Copies the solution information onto the mesh
    */
    virtual void updateSolutionInformation(escript::Data solution);

    /**
       Refines / coarsens the mesh based on the solution information
    */
    virtual void updateMeshInformation();

    /**
       Sets adaptive refinement on or off
    */
    virtual void setAdaptiveRefinement(bool) = 0;

    //List of tags currently in use
    int tags[MAXTAGS] = {-1};

    // Current number of tags
    int numberOfTags = 0;

    // Flag that determines if adaptive refinement is enabled or not
    bool adaptive_refinement = false;

    // Converter used by Boost
    // Converts the Teuchos CRS matrix to a boost::numpy array
    // typedef Tpetra::CrsMatrix<cplx_t,esys_trilinos::LO,esys_trilinos::GO,esys_trilinos::NT> crs_matrix_type;
    #ifdef ESYS_HAVE_TRILINOS
    typedef Tpetra::CrsMatrix<real_t,esys_trilinos::LO,esys_trilinos::GO,esys_trilinos::NT> real_matrix_type;
    typedef Tpetra::CrsMatrix<cplx_t,esys_trilinos::LO,esys_trilinos::GO,esys_trilinos::NT> cplx_matrix_type;
    Teuchos::RCP<Teuchos::ParameterList> params = Teuchos::parameterList();

    void initZ(bool complex);
    void initIZ(bool complex);
    // void updateZ();
    // void updateIZ();

    IndexVector zYaleRows;
    IndexVector zYaleCols;
    std::vector<IndexVector> zconnections;
    esys_trilinos::TrilinosGraph_ptr zgraph;
    IndexVector izYaleRows;
    IndexVector izYaleCols;
    std::vector<IndexVector> izconnections;
    esys_trilinos::TrilinosGraph_ptr izgraph;
    Teuchos::RCP<const Tpetra::Map<>> zccolMap;
    Teuchos::RCP<const Tpetra::Map<>> zcrowMap;
    Teuchos::RCP<const Tpetra::Map<>> zrcolMap;
    Teuchos::RCP<const Tpetra::Map<>> zrrowMap;
    Teuchos::RCP<const Tpetra::Map<>> izccolMap;
    Teuchos::RCP<const Tpetra::Map<>> izcrowMap;
    Teuchos::RCP<const Tpetra::Map<>> izrcolMap;
    Teuchos::RCP<const Tpetra::Map<>> izrrowMap;
    Teuchos::RCP<const Tpetra::Map<>> zdomainMap;
    Teuchos::RCP<const Tpetra::Map<>> izdomainMap;
    Teuchos::RCP<const Tpetra::Map<>> zrangeMap;
    Teuchos::RCP<const Tpetra::Map<>> izrangeMap;

    Teuchos::RCP<Tpetra::Map<>> f_map;
    Teuchos::RCP<Tpetra::Map<>> g_map;
    Tpetra::MultiVector<cplx_t,esys_trilinos::LO,esys_trilinos::GO,esys_trilinos::NT> fc;
    Tpetra::MultiVector<cplx_t,esys_trilinos::LO,esys_trilinos::GO,esys_trilinos::NT> gc;
    Tpetra::MultiVector<real_t,esys_trilinos::LO,esys_trilinos::GO,esys_trilinos::NT> fr;
    Tpetra::MultiVector<real_t,esys_trilinos::LO,esys_trilinos::GO,esys_trilinos::NT> gr;
    
    Teuchos::RCP<esys_trilinos::VectorType<real_t> > rlclData;
    Teuchos::RCP<esys_trilinos::VectorType<cplx_t> > clclData;
    Teuchos::RCP<esys_trilinos::VectorType<real_t> > rgblData;
    Teuchos::RCP<esys_trilinos::VectorType<cplx_t> > cgblData;
    
    // #ifdef ESYS_MPI
    // const Teuchos::RCP<const Teuchos::Comm<int>> tril_comm = esys_trilinos::TeuchosCommFromEsysComm(m_mpiInfo->comm);
    // const Teuchos::RCP<const Teuchos::Comm<int>> tril_comm = Teuchos::RCP<const Teuchos::SerialComm<int>>();
    // #else
    // const Teuchos::RCP<const Teuchos::Comm<int>> tril_comm = Teuchos::RCP<const Teuchos::SerialComm<int>>();
    // #endif

    Teuchos::RCP<real_matrix_type> rZ;
    Teuchos::RCP<real_matrix_type> rIZ;
    Teuchos::RCP<cplx_matrix_type> cZ;
    Teuchos::RCP<cplx_matrix_type> cIZ;
    // Teuchos::RCP<Tpetra::CrsMatrix<real_t,esys_trilinos::LO,esys_trilinos::GO,esys_trilinos::NT>> * pZ;
    // Teuchos::RCP<Tpetra::CrsMatrix<real_t,esys_trilinos::LO,esys_trilinos::GO,esys_trilinos::NT>> * pIZ;
    // Teuchos::RCP<Tpetra::CrsMatrix<cplx_t,esys_trilinos::LO,esys_trilinos::GO,esys_trilinos::NT>> * cpZ;
    // Teuchos::RCP<Tpetra::CrsMatrix<cplx_t,esys_trilinos::LO,esys_trilinos::GO,esys_trilinos::NT>> * cpIZ;
    #endif //ESYS_HAVE_TRILINOS

    /// stores the hanging node information
    std::vector<DoublePair> hanging_faces; 

    std::vector<std::pair<long,long>> hanging_edge_node_connections;
    std::vector<std::pair<long,long>> hanging_face_node_connections;
    std::vector<std::pair<long,long>> false_node_connections;

    virtual dim_t getNumNodes() const;

    /**
       \brief
       returns a vector of rank numbers where vec[i]=n means that rank n
       'owns' element/face element i.
    */
    virtual RankVector getOwnerVector(int fsType) const = 0;

    /**
       \brief
       returns the number of face elements in the order
       (left,right,bottom,top) on current MPI rank
    */
    virtual const dim_t* getNumFacesPerBoundary() const = 0;

    // stopwatch used when profiling 
    TicTocClock oxleytimer;

protected:

    // element order
    int m_order;

    // number of dimensions
    int m_numDim;

    // Status
    StatusType m_status;

    //max levels of refinement
    int m_refinement_levels;

    // A second MPI handle for p4est.
    sc_MPI_Comm m_p4est_mpiInfo;

    // Dirac point Node IDs
    IndexVector m_diracPointNodeIDs; //for borrowSampleID

    /// returns the number of nodes per MPI rank
    virtual dim_t getNumHangingNodes() const;

    // /// returns the number of hanging nodes per MPI rank
    // virtual int getNumHangingNodes() const = 0;

    /// returns the number of elements per MPI rank
    virtual dim_t getNumElements() const = 0;

    /// returns the number of degrees of freedom per MPI rank
    virtual dim_t getNumDOF() const = 0;

    /// returns the number of face elements on current MPI rank
    virtual dim_t getNumFaceElements() const = 0;

    #ifdef ESYS_HAVE_BOOST_NUMPY    
      virtual boost::python::numpy::ndarray getNumpyX() const;
    #endif

    // Tagmap
    TagMap m_tagMap;
    mutable std::vector<int> m_nodeTags, m_nodeTagsInUse;
    mutable std::vector<int> m_elementTags, m_elementTagsInUse;
    mutable std::vector<int> m_faceTags, m_faceTagsInUse;
    std::vector<DiracPoint> m_diracPoints;
    IndexVector m_faceOffset;

    // Function sused by the assembler
    template<typename Scalar>
    void addToSystemMatrix(escript::AbstractSystemMatrix* mat,
                           const IndexVector& nodes, dim_t numEq,
                           const std::vector<Scalar>& array) const;

    /// computes the gradient of 'in' and puts the result in 'out'
    virtual void assembleGradient(escript::Data& out, const escript::Data& in) const = 0;

    /// populates the data object 'arg' with the node coordinates
    virtual void assembleCoordinates(escript::Data& arg) const = 0;

    // adds the dirac points and tags 
    virtual void addPoints(const std::vector<double>& coords, const std::vector<int>& tags) = 0;

#ifdef ESYS_HAVE_TRILINOS
    /// returns the Trilinos matrix graph
    virtual esys_trilinos::TrilinosGraph_ptr getTrilinosGraph() const = 0;

    /// creates and returns a Trilinos CRS graph suitable to build a sparse
    /// matrix
    esys_trilinos::TrilinosGraph_ptr createTrilinosGraph(
            const IndexVector& myRows,  const IndexVector& myColumns) const;

    esys_trilinos::TrilinosGraph_ptr createTrilinosGraph(
            const IndexVector& myRows,  const IndexVector& myColumns, dim_t dof, dim_t numRows, std::vector<IndexVector> connections) const;
#endif

    /// returns occupied matrix column indices for all matrix rows
    virtual std::vector<IndexVector> getConnections(bool includeShared) const = 0;

//protected
#ifdef ESYS_HAVE_PASO
    /// returns the Paso system matrix pattern
    virtual paso::SystemMatrixPattern_ptr getPasoMatrixPattern(bool reducedRowOrder, bool reducedColOrder) const = 0;

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

    /// copies data in 'in' to 'out' (both must be on same function space)
    template<typename Scalar>
    void copyData(escript::Data& out, const escript::Data& in) const;

    /// averages data in 'in' to 'out' (from non-reduced to reduced fs)
    template<typename Scalar>
    void averageData(escript::Data& out, const escript::Data& in) const;

    /// copies data in 'in' to 'out' (from reduced to non-reduced fs)
    template<typename Scalar>
    void multiplyData(escript::Data& out, const escript::Data& in) const;

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

    #ifdef ESYS_HAVE_TRILINOS
    /// calls the right PDE assembly routine after performing input checks
    void assemblePDEHanging(Teuchos::RCP<Tpetra::CrsMatrix<double,esys_trilinos::LO,esys_trilinos::GO,esys_trilinos::NT>>* mat,
                          Assembler_ptr assembler) const;
    #endif //ESYS_HAVE_TRILINOS

    template<typename Scalar>
    void setToIntegralsWorker(std::vector<Scalar>& integrals,
                              const escript::Data& arg) const;

    /// finds the node that the given point coordinates belong to
    virtual dim_t findNode(const double *coords) const = 0;
};

#define POINTER_WRAPPER_CLASS(x) boost::shared_ptr<x>
typedef POINTER_WRAPPER_CLASS(OxleyDomain) OxleyDomain_ptr;
class Rectangle;
typedef POINTER_WRAPPER_CLASS(Rectangle) OxleyDomainRect_ptr;
class Brick;
typedef POINTER_WRAPPER_CLASS(Brick) OxleyDomainBrick_ptr;

} // end of namespace oxley


#endif //__OXLEY_EXCEPTION_H__


#endif //__OXLEY_DOMAIN_H__
