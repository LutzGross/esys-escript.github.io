
/*****************************************************************************
*
* Copyright (c) 2003-2016 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/

#ifndef __RIPLEY_MULTIRECTANGLE_H__
#define __RIPLEY_MULTIRECTANGLE_H__

#include <ripley/Rectangle.h>

namespace ripley {

/**
   \brief
   Rectangle is the 2-dimensional implementation of a RipleyDomain.
*/
class RIPLEY_DLL_API MultiRectangle: public Rectangle
{
    friend class DefaultAssembler2D;
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
              escript::SubWorld_ptr w=escript::SubWorld_ptr(),
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

#ifdef USE_BOOSTIO
    virtual void readBinaryGridFromZipped(escript::Data& out, std::string filename,
                                const ReaderParameters& params) const;
#endif

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
    virtual paso::SystemMatrixPattern_ptr getPasoMatrixPattern(
                                                    bool reducedRowOrder,
                                                    bool reducedColOrder) const;
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
};

//protected
inline dim_t MultiRectangle::getNumDOF() const
{
    return getNumDOFInAxis(0)*getNumDOFInAxis(1);
}

//protected
inline dim_t MultiRectangle::getNumDOFInAxis(unsigned axis) const
{
    EsysAssert((axis < m_numDim), "Invalid axis");
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

