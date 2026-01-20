
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/

#ifndef __RIPLEY_MULTIRECTANGLE_H__
#define __RIPLEY_MULTIRECTANGLE_H__

#include <ripley/Rectangle.h>

#ifdef ESYS_MPI
#include <escript/EsysMPI.h>
#endif

namespace ripley {

/**
   \brief
   Rectangle is the 2-dimensional implementation of a RipleyDomain.
*/
class RIPLEY_DLL_API MultiRectangle: public Rectangle
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
    MultiRectangle(dim_t n0, dim_t n1, double x0, double y0, double x1, double y1,
              int d0=-1, int d1=-1,
              const std::vector<double>& points = std::vector<double>(),
              const std::vector<int>& tags = std::vector<int>(),
              const TagMap& tagnamestonums = TagMap(),
              unsigned int subdivisions = 1
 	    );

    /**
       \brief creates a multirectangle mesh with custom MPI communicator
       \param jmpi MPI communicator info
    */
    MultiRectangle(escript::JMPI jmpi, dim_t n0, dim_t n1, double x0, double y0, double x1, double y1,
              int d0=-1, int d1=-1,
              const std::vector<double>& points = std::vector<double>(),
              const std::vector<int>& tags = std::vector<int>(),
              const TagMap& tagnamestonums = TagMap(),
              unsigned int subdivisions = 1
 	    );

    /**
       \brief
       Destructor.
    */
    ~MultiRectangle();

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
       returns the number of times each root element has been subdivided
    */
    virtual unsigned int getNumSubdivisionsPerElement() const { return m_subdivisions; }

    /**
       \brief
       returns a vector of rank numbers where vec[i]=n means that rank n
       'owns' element/face element i.
    */
    virtual RankVector getOwnerVector(int fsType) const;

protected:
    virtual void interpolateNodesToNodesFiner(const escript::Data& source, escript::Data& target, const MultiRectangle& other) const;
    virtual void interpolateNodesToElementsFiner(const escript::Data& source, escript::Data& target, const MultiRectangle& other) const;

    virtual void interpolateElementsToElementsCoarser(const escript::Data& source, escript::Data& target, const MultiRectangle& other) const;
    virtual void interpolateElementsToElementsFiner(const escript::Data& source, escript::Data& target, const MultiRectangle& other) const;

    virtual void interpolateReducedToElementsFiner(const escript::Data& source, escript::Data& target, const MultiRectangle& other) const;
    virtual void interpolateReducedToReducedFiner(const escript::Data& source, escript::Data& target, const MultiRectangle& other) const;
#ifdef ESYS_HAVE_PASO
    virtual paso::SystemMatrixPattern_ptr getPasoMatrixPattern(
                                                    bool reducedRowOrder,
                                                    bool reducedColOrder) const;
#endif
    virtual index_t getFirstInDim(unsigned axis) const;
    virtual void populateSampleIds();
    virtual dim_t getNumDOFInAxis(unsigned axis) const;
    virtual dim_t getNumDOF() const;
    virtual void populateDofMap();

    virtual dim_t findNode(const double *coords) const;

private:
    mutable std::vector<IndexVector> m_colIndices;
    mutable std::vector<IndexVector> m_rowIndices;
    unsigned int m_subdivisions;
    
    template <typename S>
    void interpolateNodesToNodesFinerWorker(const escript::Data& source, escript::Data& target, const MultiRectangle& other, S sentinel) const;
    template <typename S>
    void interpolateNodesToElementsFinerWorker(const escript::Data& source, escript::Data& target, const MultiRectangle& other, S sentinel) const;

    template <typename S>
    void interpolateElementsToElementsCoarserWorker(const escript::Data& source, escript::Data& target, const MultiRectangle& other, S sentinel) const;
    template <typename S>
    void interpolateElementsToElementsFinerWorker(const escript::Data& source, escript::Data& target, const MultiRectangle& other, S sentinel) const;

    template <typename S>
    void interpolateReducedToElementsFinerWorker(const escript::Data& source, escript::Data& target, const MultiRectangle& other, S sentinel) const;
    template <typename S>
    void interpolateReducedToReducedFinerWorker(const escript::Data& source, escript::Data& target, const MultiRectangle& other, S sentinel) const;
    
    
    
};

//protected
inline dim_t MultiRectangle::getNumDOF() const
{
    return getNumDOFInAxis(0)*getNumDOFInAxis(1);
}

//protected
inline dim_t MultiRectangle::getNumDOFInAxis(unsigned axis) const
{
    ESYS_ASSERT(axis < m_numDim, "Invalid axis");
    dim_t res = m_ownNE[axis] + 1;
    if (m_offset[axis] + m_NE[axis] < m_gNE[axis]) {
        res--;
    }
    return res;
}

//protected
index_t MultiRectangle::getFirstInDim(unsigned axis) const
{
    return m_offset[axis] == 0 ? 0 : m_subdivisions;
}

} // end of namespace ripley

#endif // __RIPLEY_MULTIRECTANGLE_H__

