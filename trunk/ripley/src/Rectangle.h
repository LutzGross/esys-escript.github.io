
/*****************************************************************************
*
* Copyright (c) 2003-2013 by University of Queensland
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

#ifndef __RIPLEY_RECTANGLE_H__
#define __RIPLEY_RECTANGLE_H__

#include <ripley/RipleyDomain.h>

struct Paso_Connector;

namespace ripley {

/**
   \brief
   Rectangle is the 2-dimensional implementation of a RipleyDomain.
*/
class RIPLEY_DLL_API Rectangle: public RipleyDomain
{
public:

    /**
       \brief creates a rectangular mesh with n0 x n1 elements over the
              rectangle [x0,x1] x [y0,y1].
       \param n0,n1 number of elements in each dimension
       \param x0,y0,x1,y1 coordinates of bottom-left and top-right corners
       \param d0,d1 number of subdivisions in each dimension
    */
    Rectangle(int n0, int n1, double x0, double y0, double x1, double y1,
              int d0=-1, int d1=-1);

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
    virtual void writeBinaryGrid(const escript::Data& in,
                                 std::string filename,
                                 int byteOrder, int dataType) const;

    /**
       \brief
       returns the array of reference numbers for a function space type
       \param fsType The function space type
    */
    const int* borrowSampleReferenceIDs(int fsType) const;

    /**
       \brief
       returns true if this rank owns the sample id.
    */
    virtual bool ownSample(int fs_code, index_t id) const;

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
    virtual int getNumDataPointsGlobal() const;

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
    virtual const int* getNumNodesPerDim() const { return m_NN; }

    /**
       \brief
       returns the number of elements per MPI rank in each dimension
    */
    virtual const int* getNumElementsPerDim() const { return m_NE; }

    /**
       \brief
       returns the number of face elements in the order
       (left,right,bottom,top) on current MPI rank
    */
    virtual const int* getNumFacesPerBoundary() const { return m_faceCount; }

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
    virtual double getLocalCoordinate(int index, int dim) const;

    /**
       \brief
       returns the tuple (origin, spacing, number_of_elements)
    */
    virtual boost::python::tuple getGridParameters() const;
    
    /**
     * \brief 
       Returns a Data object filled with random data passed through filter.
    */ 
    escript::Data randomFill(long seed, const boost::python::tuple& filter) const;
    

protected:
    virtual dim_t getNumNodes() const;
    virtual dim_t getNumElements() const;
    virtual dim_t getNumFaceElements() const;
    virtual dim_t getNumDOF() const;
    virtual dim_t insertNeighbourNodes(IndexVector& index, index_t node) const;
    virtual void assembleCoordinates(escript::Data& arg) const;
    virtual void assembleGradient(escript::Data& out, escript::Data& in) const;
    virtual void assembleIntegrate(DoubleVector& integrals, escript::Data& arg) const;
    virtual void assemblePDESingle(Paso_SystemMatrix* mat, escript::Data& rhs,
            const escript::Data& A, const escript::Data& B,
            const escript::Data& C, const escript::Data& D,
            const escript::Data& X, const escript::Data& Y) const;
    virtual void assemblePDEBoundarySingle(Paso_SystemMatrix* mat,
            escript::Data& rhs, const escript::Data& d,
            const escript::Data& y) const;
    virtual void assemblePDESingleReduced(Paso_SystemMatrix* mat,
            escript::Data& rhs, const escript::Data& A, const escript::Data& B,
            const escript::Data& C, const escript::Data& D,
            const escript::Data& X, const escript::Data& Y) const;
    virtual void assemblePDEBoundarySingleReduced(Paso_SystemMatrix* mat,
            escript::Data& rhs, const escript::Data& d,
            const escript::Data& y) const;
    virtual void assemblePDESystem(Paso_SystemMatrix* mat, escript::Data& rhs,
            const escript::Data& A, const escript::Data& B,
            const escript::Data& C, const escript::Data& D,
            const escript::Data& X, const escript::Data& Y) const;
    virtual void assemblePDEBoundarySystem(Paso_SystemMatrix* mat,
            escript::Data& rhs, const escript::Data& d,
            const escript::Data& y) const;
    virtual void assemblePDESystemReduced(Paso_SystemMatrix* mat,
            escript::Data& rhs, const escript::Data& A, const escript::Data& B,
            const escript::Data& C, const escript::Data& D,
            const escript::Data& X, const escript::Data& Y) const;
    virtual void assemblePDEBoundarySystemReduced(Paso_SystemMatrix* mat,
            escript::Data& rhs, const escript::Data& d,
            const escript::Data& y) const;
    virtual Paso_SystemMatrixPattern* getPattern(bool reducedRowOrder, bool reducedColOrder) const;
    virtual void interpolateNodesOnElements(escript::Data& out,
                                       escript::Data& in, bool reduced) const;
    virtual void interpolateNodesOnFaces(escript::Data& out, escript::Data& in,
                                         bool reduced) const;
    virtual void nodesToDOF(escript::Data& out, escript::Data& in) const;
    virtual void dofToNodes(escript::Data& out, escript::Data& in) const;

private:
    void populateSampleIds();
    void createPattern();
    void addToMatrixAndRHS(Paso_SystemMatrix* S, escript::Data& F,
           const DoubleVector& EM_S, const DoubleVector& EM_F,
           bool addS, bool addF, int firstNode, int nEq=1, int nComp=1) const;

    template<typename ValueType>
    void readBinaryGridImpl(escript::Data& out, const std::string& filename,
                            const ReaderParameters& params) const;

    template<typename ValueType>
    void writeBinaryGridImpl(const escript::Data& in,
                             const std::string& filename, int byteOrder) const;

    /// total number of elements in each dimension
    dim_t m_gNE[2];

    /// origin of domain
    double m_origin[2];

    /// side lengths of domain
    double m_length[2];

    /// grid spacings / cell sizes of domain
    double m_dx[2];

    /// number of spatial subdivisions
    dim_t m_NX[2];

    /// number of elements for this rank in each dimension including shared
    dim_t m_NE[2];

    /// number of own elements for this rank in each dimension
    dim_t m_ownNE[2];

    /// number of nodes for this rank in each dimension
    dim_t m_NN[2];

    /// first node on this rank is at (offset0,offset1) in global mesh
    dim_t m_offset[2];

    /// number of face elements per edge (left, right, bottom, top)
    int m_faceCount[4];

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

    // Paso connector used by the system matrix and to interpolate DOF to
    // nodes
    Paso_Connector* m_connector;

    // the Paso System Matrix pattern
    Paso_SystemMatrixPattern* m_pattern;
};

////////////////////////////// inline methods ////////////////////////////////

inline int Rectangle::getNumDataPointsGlobal() const
{
    return (m_gNE[0]+1)*(m_gNE[1]+1);
}

inline double Rectangle::getLocalCoordinate(int index, int dim) const
{
    EsysAssert((dim>=0 && dim<2), "'dim' out of bounds");
    EsysAssert((index>=0 && index<m_NN[dim]), "'index' out of bounds");
    return m_origin[dim]+m_dx[dim]*(m_offset[dim]+index);
}

inline boost::python::tuple Rectangle::getGridParameters() const
{
    return boost::python::make_tuple(
            boost::python::make_tuple(m_origin[0], m_origin[1]),
            boost::python::make_tuple(m_dx[0], m_dx[1]),
            boost::python::make_tuple(m_gNE[0], m_gNE[1]));
}

inline Paso_SystemMatrixPattern* Rectangle::getPattern(bool reducedRowOrder,
                                                   bool reducedColOrder) const
{
    // TODO: reduced
    return m_pattern;
}


//protected
inline dim_t Rectangle::getNumDOF() const
{
    return (m_gNE[0]+1)/m_NX[0]*(m_gNE[1]+1)/m_NX[1];
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


} // end of namespace ripley

#endif // __RIPLEY_RECTANGLE_H__

