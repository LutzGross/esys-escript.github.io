
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

#ifndef __RIPLEY_MULTIBRICK_H__
#define __RIPLEY_MULTIBRICK_H__

#include <ripley/Brick.h>

namespace ripley {

/**
   \brief
   Brick is the 3-dimensional implementation of a RipleyDomain.
*/
class RIPLEY_DLL_API MultiBrick: public Brick
{
    template<class Scalar> friend class DefaultAssembler3D;
    friend class WaveAssembler3D;
    friend class LameAssembler3D;
public:

    /**
       \brief creates a hexagonal mesh with n0 x n1 x n2 elements over the
              brick [x0,x1] x [y0,y1] x [z0,z1].
       \param n0,n1,n2 number of elements in each dimension
       \param x0,y0,z0,x1,y1,z1 coordinates of corner nodes of the brick
       \param d0,d1,d2 number of subdivisions in each dimension
    */
    MultiBrick(dim_t n0, dim_t n1, dim_t n2, double x0, double y0, double z0,
          double x1, double y1, double z1, int d0=-1, int d1=-1, int d2=-1,
          const std::vector<double>& points = std::vector<double>(),
          const std::vector<int>& tags = std::vector<int>(),
          const TagMap& tagnamestonums = TagMap(),
              unsigned int subdivisions = 1
 	    );

    /**
       \brief creates a multibrick mesh with custom MPI communicator
       \param jmpi MPI communicator info
    */
    MultiBrick(escript::JMPI jmpi, dim_t n0, dim_t n1, dim_t n2, double x0, double y0, double z0,
          double x1, double y1, double z1, int d0=-1, int d1=-1, int d2=-1,
          const std::vector<double>& points = std::vector<double>(),
          const std::vector<int>& tags = std::vector<int>(),
          const TagMap& tagnamestonums = TagMap(),
              unsigned int subdivisions = 1
 	    );

    /**
       \brief
       Destructor.
    */
    ~MultiBrick();

    /**
       \brief
       interpolates data given on source onto target where source and target
       are given on different domains
    */
    virtual void interpolateAcross(escript::Data& target,
                                   const escript::Data& source) const;

    /**
       \brief
       Checks that the given interpolation is possible, throw RipleyExceptions
       if not
    */
    void validateInterpolationAcross(int fsType_source,
              const escript::AbstractDomain& domain, int fsType_target) const;
    
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
       writes information about the mesh to standard output
       \param full whether to print additional data
    */
    virtual void Print_Mesh_Info(const bool full=false) const;

    /**
       \brief
       returns the number of times each root element has been subdivided
    */
    virtual unsigned int getNumSubdivisionsPerElement() const { return m_subdivisions; }

    /**
       \brief
       returns the number of face elements in the order
       (left,right,bottom,top,front,back) on current MPI rank
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
       returns a vector of rank numbers where vec[i]=n means that rank n
       'owns' element/face element i.
    */
    virtual RankVector getOwnerVector(int fsType) const;

protected:
    virtual IndexVector getDiagonalIndices(bool upperOnly) const;
    virtual void interpolateNodesToNodesFiner(const escript::Data& source, escript::Data& target, const MultiBrick& other) const;
    virtual void interpolateNodesToElementsFiner(const escript::Data& source, escript::Data& target, const MultiBrick& other) const;

    virtual void interpolateElementsToElementsCoarser(const escript::Data& source, escript::Data& target, const MultiBrick& other) const;
    virtual void interpolateElementsToElementsFiner(const escript::Data& source, escript::Data& target, const MultiBrick& other) const;

    virtual void interpolateReducedToElementsFiner(const escript::Data& source, escript::Data& target, const MultiBrick& other) const;
    virtual void interpolateReducedToReducedFiner(const escript::Data& source, escript::Data& target, const MultiBrick& other) const;

    virtual void populateSampleIds();
    virtual std::vector<IndexVector> getConnections(bool includeShared) const;

    dim_t findNode(const double *coords) const;

    const unsigned int m_subdivisions;
private:
    template <typename S>
    void interpolateNodesToNodesFinerWorker(const escript::Data& source, escript::Data& target, const MultiBrick& other, S sentinel) const;
    template <typename S>    
    void interpolateNodesToElementsFinerWorker(const escript::Data& source, escript::Data& target, const MultiBrick& other, S sentinel) const;

    template <typename S>    
    void interpolateElementsToElementsCoarserWorker(const escript::Data& source, escript::Data& target, const MultiBrick& other, S sentinel) const;
    template <typename S>
    void interpolateElementsToElementsFinerWorker(const escript::Data& source, escript::Data& target, const MultiBrick& other, S sentinel) const;

    template <typename S>
    void interpolateReducedToElementsFinerWorker(const escript::Data& source, escript::Data& target, const MultiBrick& other, S sentinel) const;
    template <typename S>
    void interpolateReducedToReducedFinerWorker(const escript::Data& source, escript::Data& target, const MultiBrick& other, S sentinel) const;
  
    
};




} // end of namespace ripley

#endif // __RIPLEY_MULTIBRICK_H__

