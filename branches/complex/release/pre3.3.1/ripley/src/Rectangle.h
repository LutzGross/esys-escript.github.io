
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
    virtual void readBinaryGrid(escript::Data& out, std::string filename,
                                const std::vector<int>& first,
                                const std::vector<int>& numValues) const;

    /**
    */
    virtual void readNcGrid(escript::Data& out, std::string filename,
            std::string varname, const std::vector<int>& first,
            const std::vector<int>& numValues) const;

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
    virtual int getNumDataPointsGlobal() const { return (m_gNE0+1)*(m_gNE1+1); }

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
    virtual IndexVector getNumNodesPerDim() const;

    /**
       \brief
       returns the number of elements per MPI rank in each dimension
    */
    virtual IndexVector getNumElementsPerDim() const;

    /**
       \brief
       returns the number of face elements in the order
       (left,right,bottom,top,[front,back]) on current MPI rank
    */
    virtual IndexVector getNumFacesPerBoundary() const;

    /**
       \brief
       returns the node distribution vector
    */
    virtual IndexVector getNodeDistribution() const { return m_nodeDistribution; }

    /**
       \brief
       returns the number of spatial subdivisions in each dimension
    */
    virtual IndexVector getNumSubdivisionsPerDim() const;

    /**
       \brief
       returns the first coordinate value and the node spacing along given
       dimension as a pair
    */
    virtual std::pair<double,double> getFirstCoordAndSpacing(dim_t dim) const;

protected:
    virtual dim_t getNumNodes() const { return m_N0*m_N1; }
    virtual dim_t getNumElements() const { return m_NE0*m_NE1; }
    virtual dim_t getNumFaceElements() const;
    virtual dim_t getNumDOF() const;
    virtual dim_t insertNeighbourNodes(IndexVector& index, index_t node) const;
    virtual void assembleCoordinates(escript::Data& arg) const;
    virtual void assembleGradient(escript::Data& out, escript::Data& in) const;
    virtual void assembleIntegrate(std::vector<double>& integrals, escript::Data& arg) const;
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
           const std::vector<double>& EM_S, const std::vector<double>& EM_F,
           bool addS, bool addF, index_t firstNode, dim_t nEq=1, dim_t nComp=1) const;

    /// total number of elements in each dimension
    dim_t m_gNE0, m_gNE1;

    /// location of domain
    double m_x0, m_y0;

    /// side lengths of domain
    double m_l0, m_l1;

    /// number of spatial subdivisions
    int m_NX, m_NY;

    /// number of elements for this rank in each dimension including shared
    dim_t m_NE0, m_NE1;

    /// number of own elements for this rank in each dimension
    dim_t m_ownNE0, m_ownNE1;

    /// number of nodes for this rank in each dimension
    dim_t m_N0, m_N1;

    /// first node on this rank is at (offset0,offset1) in global mesh
    dim_t m_offset0, m_offset1;

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

} // end of namespace ripley

#endif // __RIPLEY_RECTANGLE_H__

