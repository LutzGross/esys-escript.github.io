
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

#ifndef __RIPLEY_RECTANGLE_H__
#define __RIPLEY_RECTANGLE_H__

#ifdef _WIN32 // for M_1_PI
#define _USE_MATH_DEFINES
#include <math.h>
#endif

#ifdef ESYS_MPI
#include <escript/EsysMPI.h>
#endif

#ifdef ESYS_HAVE_MPI4PY
#include <mpi4py/mpi4py.h>
#endif

#include <ripley/RipleyDomain.h>

namespace ripley {

/**
   \brief
   Rectangle is the 2-dimensional implementation of a RipleyDomain.
*/
class RIPLEY_DLL_API Rectangle: public RipleyDomain
{
    template<class Scalar> friend class DefaultAssembler2D;
    friend class WaveAssembler2D;
    friend class LameAssembler2D;
public:

    /**
       \brief creates a rectangular mesh with n0 x n1 elements over the
              rectangle [x0,x1] x [y0,y1].
       \param n0,n1 number of elements in each dimension
       \param x0,y0,x1,y1 coordinates of bottom-left and top-right corners
       \param d0,d1 number of subdivisions in each dimension
    */
    Rectangle(dim_t n0, dim_t n1, double x0, double y0, double x1, double y1,
              int d0=-1, int d1=-1,
              const std::vector<double>& points = std::vector<double>(),
              const std::vector<int>& tags = std::vector<int>(),
              const TagMap& tagnamestonums = TagMap()
 	    );

    /**
       \brief creates a rectangular mesh with custom MPI communicator
       \param jmpi MPI communicator info
    */
    Rectangle(escript::JMPI jmpi, dim_t n0, dim_t n1, double x0, double y0, double x1, double y1,
              int d0=-1, int d1=-1,
              const std::vector<double>& points = std::vector<double>(),
              const std::vector<int>& tags = std::vector<int>(),
              const TagMap& tagnamestonums = TagMap()
	    );

    /**
       \brief
       Destructor.
    */
    ~Rectangle();

    /**
       \brief
       returns a description for this domain
    */
    virtual std::string getDescription() const;

    /**
       \brief equality operator
    */
    virtual bool operator==(const escript::AbstractDomain& other) const;

    /**
       \brief
       writes the current mesh to a file with the given name
       \param filename The name of the file to write to
    */
    virtual void write(const std::string& filename) const;

    /**
       \brief
       dumps the mesh to a file with the given name
       \param filename The name of the output file
    */
    void dump(const std::string& filename) const;

    /**
    */
    virtual void readNcGrid(escript::Data& out, std::string filename,
            std::string varname, const ReaderParameters& params) const;

    /**
    */
    virtual void readBinaryGrid(escript::Data& out, std::string filename,
                                const ReaderParameters& params) const;

    /**
    */
    virtual void readBinaryGridFromZipped(escript::Data& out, std::string filename,
                                const ReaderParameters& params) const;

    /**
    */
    virtual void writeBinaryGrid(const escript::Data& in,
                                 std::string filename,
                                 int byteOrder, int dataType) const;

    /**
       \brief
       returns the array of reference numbers for a function space type
       \param fsType The function space type
    */
    const dim_t* borrowSampleReferenceIDs(int fsType) const;

    /**
       \brief
       returns true if this rank owns the sample id.
    */
    virtual bool ownSample(int fsType, index_t id) const;

    /**
       \brief
       copies the surface normals at data points into out. The actual function
       space to be considered is defined by out. out has to be defined on this
       domain.
    */
    virtual void setToNormal(escript::Data& out) const;

    /**
       \brief
       copies the size of samples into out. The actual function space to be
       considered is defined by out. out has to be defined on this domain.
    */
    virtual void setToSize(escript::Data& out) const;

    /**
       \brief
       returns the number of data points summed across all MPI processes
    */
    virtual dim_t getNumDataPointsGlobal() const;

    /**
       \brief
       writes information about the mesh to standard output
       \param full whether to print additional data
    */
    virtual void Print_Mesh_Info(const bool full=false) const;

    /**
       \brief
       returns the number of nodes per MPI rank in each dimension
    */
    virtual const dim_t* getNumNodesPerDim() const { return m_NN; }

    /**
       \brief
       returns the number of elements per MPI rank in each dimension
    */
    virtual const dim_t* getNumElementsPerDim() const { return m_NE; }

    /**
       \brief
       returns the number of face elements in the order
       (left,right,bottom,top) on current MPI rank
    */
    virtual const dim_t* getNumFacesPerBoundary() const { return m_faceCount; }

    /**
       \brief
       returns the node distribution vector
    */
    virtual IndexVector getNodeDistribution() const { return m_nodeDistribution; }

    /**
       \brief
       returns the number of spatial subdivisions in each dimension
    */
    virtual const int* getNumSubdivisionsPerDim() const { return m_NX; }

    /**
       \brief
       returns the index'th coordinate value in given dimension for this rank
    */
    virtual double getLocalCoordinate(index_t index, int dim) const;

    /**
       \brief
       returns the tuple (origin, spacing, number_of_elements)
    */
    virtual boost::python::tuple getGridParameters() const;

    /**
     * \brief 
       returns a Data object filled with random data passed through filter.
    */ 
    virtual escript::Data randomFill(const escript::DataTypes::ShapeType& shape,
                                 const escript::FunctionSpace& what, long seed,
                                 const boost::python::tuple& filter) const;
    
    /**
       \brief
       creates and returns an assembler of the requested type.
    */
    virtual Assembler_ptr createAssembler(std::string type,
                                          const DataMap& options) const;

    /**
       \brief
       returns the lengths of the domain
    */
    const double *getLength() const { return m_length; }

    /**
       \brief
       returns the lengths of an element
    */
    const double *getElementLength() const { return m_dx; }

    /**
       \brief
       returns a vector of rank numbers where vec[i]=n means that rank n
       'owns' element/face element i.
    */
    virtual RankVector getOwnerVector(int fsType) const;

protected:
    virtual dim_t getNumNodes() const;
    virtual dim_t getNumElements() const;
    virtual dim_t getNumFaceElements() const;
    virtual dim_t getNumDOF() const;
    virtual dim_t getNumDOFInAxis(unsigned axis) const;
    virtual index_t getFirstInDim(unsigned axis) const;

    virtual IndexVector getDiagonalIndices(bool upperOnly) const;
    virtual void assembleCoordinates(escript::Data& arg) const;
    virtual void assembleGradient(escript::Data& out,
                                  const escript::Data& in) const;
    virtual void assembleIntegrate(std::vector<real_t>& integrals,
                                   const escript::Data& arg) const;
    virtual void assembleIntegrate(std::vector<cplx_t>& integrals,
                                   const escript::Data& arg) const;
    virtual std::vector<IndexVector> getConnections(bool includeShared=false) const;

#ifdef ESYS_HAVE_TRILINOS
    virtual esys_trilinos::const_TrilinosGraph_ptr getTrilinosGraph() const;
#endif

#ifdef ESYS_HAVE_PASO
    virtual paso::SystemMatrixPattern_ptr getPasoMatrixPattern(
                             bool reducedRowOrder, bool reducedColOrder) const;
#endif
    virtual void interpolateNodesOnElements(escript::Data& out,
                                  const escript::Data& in, bool reduced) const;	  
    virtual void interpolateNodesOnFaces(escript::Data& out,
                                         const escript::Data& in,
                                         bool reduced) const;

    virtual void nodesToDOF(escript::Data& out, const escript::Data& in) const;
    virtual dim_t getDofOfNode(dim_t node) const;

    virtual void populateSampleIds();
    virtual void populateDofMap();

    template<typename Scalar>
    void assembleGradientImpl(escript::Data& out,
                              const escript::Data& in) const;

    template<typename Scalar>
    void addToMatrixAndRHS(escript::AbstractSystemMatrix* S, escript::Data& F,
           const std::vector<Scalar>& EM_S, const std::vector<Scalar>& EM_F,
           bool addS, bool addF, index_t firstNode, int nEq=1, int nComp=1) const;

    template<typename ValueType>
    void readBinaryGridImpl(escript::Data& out, const std::string& filename,
                            const ReaderParameters& params) const;

#ifdef ESYS_HAVE_BOOST_IO
    template<typename ValueType>
    void readBinaryGridZippedImpl(escript::Data& out, 
            const std::string& filename, const ReaderParameters& params) const;
#endif

    template<typename ValueType>
    void writeBinaryGridImpl(const escript::Data& in,
                             const std::string& filename, int byteOrder) const;

    virtual dim_t findNode(const double *coords) const;

    
    escript::Data randomFillWorker(const escript::DataTypes::ShapeType& shape,
       long seed, const boost::python::tuple& filter) const;    
    
    /// total number of elements in each dimension
    dim_t m_gNE[2];

    /// origin of domain
    double m_origin[2];

    /// side lengths of domain
    double m_length[2];

    /// grid spacings / cell sizes of domain
    double m_dx[2];

    /// number of spatial subdivisions
    int m_NX[2];

    /// number of elements for this rank in each dimension including shared
    dim_t m_NE[2];

    /// number of own elements for this rank in each dimension
    dim_t m_ownNE[2];

    /// number of nodes for this rank in each dimension
    dim_t m_NN[2];

    /// first node on this rank is at (offset0,offset1) in global mesh
    dim_t m_offset[2];

    /// number of face elements per edge (left, right, bottom, top)
    dim_t m_faceCount[4];

    /// faceOffset[i]=-1 if face i is not an external face, otherwise it is
    /// the index of that face (where i: 0=left, 1=right, 2=bottom, 3=top)
    IndexVector m_faceOffset;

    /// vector of sample reference identifiers
    IndexVector m_dofId;
    IndexVector m_nodeId;
    IndexVector m_elementId;
    IndexVector m_faceId;

    // vector with first node id on each rank
    IndexVector m_nodeDistribution;

    // vector that maps each node to a DOF index (used for the coupler)
    IndexVector m_dofMap;

#ifdef ESYS_HAVE_PASO
    // the Paso System Matrix pattern
    mutable paso::SystemMatrixPattern_ptr m_pattern;
#endif

#ifdef ESYS_HAVE_TRILINOS
    /// Trilinos graph structure, cached for efficiency
    mutable esys_trilinos::const_TrilinosGraph_ptr m_graph;
#endif
private:
    template<typename Scalar>
    void assembleIntegrateImpl(std::vector<Scalar>& integrals, const escript::Data& arg) const;

    template <typename S>
    void interpolateNodesOnElementsWorker(escript::Data& out,
                                  const escript::Data& in, bool reduced, S sentinel) const;	  
    template <typename S>
    void interpolateNodesOnFacesWorker(escript::Data& out,
                                         const escript::Data& in,
                                         bool reduced, S sentinel) const;  
};

////////////////////////////// inline methods ////////////////////////////////
inline dim_t Rectangle::getDofOfNode(dim_t node) const
{
    return m_dofMap[node];
}

inline dim_t Rectangle::getNumDataPointsGlobal() const
{
    return (m_gNE[0]+1)*(m_gNE[1]+1);
}

inline double Rectangle::getLocalCoordinate(index_t index, int dim) const
{
    ESYS_ASSERT(dim>=0 && dim<2, "'dim' out of bounds");
    ESYS_ASSERT(index>=0 && index<m_NN[dim], "'index' out of bounds");
    return m_origin[dim]+m_dx[dim]*(m_offset[dim]+index);
}

inline boost::python::tuple Rectangle::getGridParameters() const
{
    return boost::python::make_tuple(
            boost::python::make_tuple(m_origin[0], m_origin[1]),
            boost::python::make_tuple(m_dx[0], m_dx[1]),
            boost::python::make_tuple(m_gNE[0], m_gNE[1]));
}

//protected
inline dim_t Rectangle::getNumDOF() const
{
    return (m_gNE[0]+1)/m_NX[0]*(m_gNE[1]+1)/m_NX[1];
}

//protected
inline dim_t Rectangle::getNumDOFInAxis(unsigned axis) const
{
    ESYS_ASSERT(axis < m_numDim, "Invalid axis");
    return (m_gNE[axis]+1)/m_NX[axis];
}

//protected
inline dim_t Rectangle::getNumNodes() const
{
    return m_NN[0]*m_NN[1];
}

//protected
inline dim_t Rectangle::getNumElements() const
{
    return m_NE[0]*m_NE[1];
}

//protected
inline dim_t Rectangle::getNumFaceElements() const
{
    return m_faceCount[0] + m_faceCount[1] + m_faceCount[2] + m_faceCount[3];
}

//protected
inline index_t Rectangle::getFirstInDim(unsigned axis) const
{
    return m_offset[axis] == 0 ? 0 : 1;
}

} // end of namespace ripley

#endif // __RIPLEY_RECTANGLE_H__

