
/*******************************************************
*
* Copyright (c) 2003-2011 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/

#ifndef __RIPLEY_BRICK_H__
#define __RIPLEY_BRICK_H__

#include <ripley/RipleyDomain.h>

namespace ripley {

/**
   \brief
   Brick is the 3-dimensional implementation of a RipleyDomain.
*/
class Brick: public RipleyDomain
{
public:

    /**
       \brief creates a hexagonal mesh with n0 x n1 x n2 elements over the
              brick [0,l0] x [0,l1] x [0,l2].
       \param n0,n1,n2 number of elements in each dimension
       \param l0,l1,l2 length of each side of brick
       \param d0,d1,d2 number of subdivisions in each dimension
    */
    RIPLEY_DLL_API
    Brick(int n0, int n1, int n2, double l0, double l1, double l2, int d0,
          int d1, int d2);

    /**
       \brief
       Destructor.
    */
    RIPLEY_DLL_API
    ~Brick();

    /**
       \brief
       returns a description for this domain
    */
    RIPLEY_DLL_API
    virtual std::string getDescription() const;

    /**
       \brief equality operator
    */
    RIPLEY_DLL_API
    virtual bool operator==(const escript::AbstractDomain& other) const;

    /**
       \brief
       dumps the mesh to a file with the given name
       \param filename The name of the output file
    */
    RIPLEY_DLL_API
    void dump(const std::string& filename) const;

    /**
       \brief
       returns the reference number of the given sample number
       \param fsType The function space type
    */
    RIPLEY_DLL_API
    const int* borrowSampleReferenceIDs(int fsType) const;

    /**
       \brief
       returns true if this rank owns the sample id.
    */
    RIPLEY_DLL_API
    virtual bool ownSample(int fs_code, index_t id) const;

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
       copies the integrals of the function defined by arg into integrals.
       arg has to be defined on this domain.
    */
    RIPLEY_DLL_API
    virtual void setToIntegrals(std::vector<double>& integrals, const escript::Data& arg) const;

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
       returns the number of data points summed across all MPI processes
    */
    RIPLEY_DLL_API
    virtual int getNumDataPointsGlobal() const { return (m_gNE0+1)*(m_gNE1+1)*(m_gNE2+1); }

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
    virtual IndexVector getNodeDistribution() const { return m_nodeDistribution; }

    /**
       \brief
       returns the first coordinate value and the node spacing along given
       dimension as a pair
    */
    RIPLEY_DLL_API
    virtual std::pair<double,double> getFirstCoordAndSpacing(dim_t dim) const;

protected:
    virtual dim_t getNumNodes() const { return m_N0*m_N1*m_N2; }
    virtual dim_t getNumElements() const { return m_NE0*m_NE1*m_NE2; }
    virtual dim_t getNumFaceElements() const;
    virtual dim_t getNumDOF() const;
    virtual void assembleCoordinates(escript::Data& arg) const;
    virtual Paso_SystemMatrixPattern* getPattern(bool reducedRowOrder, bool reducedColOrder) const;
    virtual void interpolateNodesOnElements(escript::Data& out,
                                       escript::Data& in, bool reduced) const;
    virtual void interpolateNodesOnFaces(escript::Data& out, escript::Data& in,
                                         bool reduced) const;

private:
    void populateSampleIds();

    /// total number of elements in each dimension
    dim_t m_gNE0, m_gNE1, m_gNE2;

    /// side lengths of domain
    double m_l0, m_l1, m_l2;

    /// number of spatial subdivisions
    int m_NX, m_NY, m_NZ;

    /// number of elements for this rank in each dimension
    dim_t m_NE0, m_NE1, m_NE2;

    /// number of nodes for this rank in each dimension
    dim_t m_N0, m_N1, m_N2;

    /// first node on this rank is at (offset0,offset1,offset2) in global mesh
    dim_t m_offset0, m_offset1, m_offset2;

    /// faceOffset[i]=-1 if face i is not an external face, otherwise it is
    /// the index of that face (where i: 0=left, 1=right, 2=bottom, 3=top,
    /// 4=front, 5=back)
    IndexVector m_faceOffset;

    /// vector of sample reference identifiers
    IndexVector m_nodeId;
    IndexVector m_elementId;
    IndexVector m_faceId;

    // vector with first node id on each rank
    IndexVector m_nodeDistribution;
};

} // end of namespace ripley

#endif // __RIPLEY_BRICK_H__

