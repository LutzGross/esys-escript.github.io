
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

#include <ripley/MultiRectangle.h>
#include <ripley/blocktools.h>
#include <ripley/domainhelpers.h>

#include <escript/DataFactory.h>
#include <escript/FunctionSpaceFactory.h>

#define FIRST_QUAD 0.21132486540518711775
#define SECOND_QUAD 0.78867513459481288225

#include <algorithm>
#include <iomanip>
#include <limits>

using std::vector;
using std::string;

namespace ripley {

MultiRectangle::MultiRectangle(dim_t n0, dim_t n1, double x0, double y0,
                     double x1, double y1, int d0, int d1,
                     const vector<double>& points,
                     const vector<int>& tags,
                     const TagMap& tagnamestonums,
                     escript::SubWorld_ptr w, unsigned int subdivisions)
     : Rectangle(n0,n1, x0,y0, x1,y1, d0,d1, points, tags, tagnamestonums, w),
       m_subdivisions(subdivisions)
{
    if (subdivisions == 0 || (subdivisions & (subdivisions - 1)) != 0)
        throw RipleyException("Element subdivisions must be a power of two");

    dim_t oldNN[2] = {0};

    if (d0 == 0 || d1 == 0)
        throw RipleyException("Domain subdivisions must be positive");

    for (int i = 0; i < 2; i++) {
        m_NE[i] *= subdivisions;
        oldNN[i] = m_NN[i];
        m_NN[i] = m_NE[i] + 1;
        m_gNE[i] *= subdivisions;
        m_ownNE[i] *= subdivisions;
        m_dx[i] /= subdivisions;
        m_faceCount[i] *= subdivisions;
        m_faceCount[2+i] *= subdivisions;
        m_offset[i] *= subdivisions;
    }
    populateSampleIds();
    
    const dim_t nDirac = m_diracPoints.size();
#pragma omp parallel for
    for (int i = 0; i < nDirac; i++) {
        const dim_t node = m_diracPoints[i].node;
        const dim_t x = node % oldNN[0];
        const dim_t y = node / oldNN[0];
        m_diracPoints[i].node = INDEX2(x*subdivisions, y*subdivisions, m_NN[0]);
        m_diracPointNodeIDs[i] = m_diracPoints[i].node;
    }
}

MultiRectangle::~MultiRectangle()
{
}

void MultiRectangle::validateInterpolationAcross(int fsType_source,
        const escript::AbstractDomain& domain, int fsType_target) const
{
    const MultiRectangle *other = dynamic_cast<const MultiRectangle *>(&domain);
    if (other == NULL)
        throw RipleyException("Invalid interpolation: domains must both be instances of MultiRectangle");

    const double *len = other->getLength();
    const int *subdivs = other->getNumSubdivisionsPerDim();
    const dim_t *elements = other->getNumElementsPerDim();
    const unsigned int level = other->getNumSubdivisionsPerElement();
    const unsigned int factor = m_subdivisions > level ? m_subdivisions/level : level/m_subdivisions;
    if ((factor & (factor - 1)) != 0) //factor == 2**x
        throw RipleyException("Invalid interpolation: elemental subdivisions of each domain must be powers of two");

    if (other->getMPIComm() != m_mpiInfo->comm)
        throw RipleyException("Invalid interpolation: Domains are on different communicators");
    for (int i = 0; i < m_numDim; i++) {
        if (m_length[i] != len[i]) {
            throw RipleyException("Invalid interpolation: domain length mismatch");
        }
        if (m_NX[i] != subdivs[i]) {
            throw RipleyException("Invalid interpolation: domain process subdivision mismatch");
        }
        if (m_subdivisions > level) {
            if (m_NE[i]/elements[i] != factor) {
                std::cerr << "m_ownNE[i]/elements[i] = "
                    << m_ownNE[i]/elements[i] << " != " << factor << std::endl;
                throw RipleyException("Invalid interpolation: element factor mismatch");
            }
        } else {
            if (elements[i]/m_NE[i] != factor) {
                throw RipleyException("Invalid interpolation: element factor mismatch");
            }
        }
    }
}

void MultiRectangle::interpolateNodesToNodesFiner(const escript::Data& source,
        escript::Data& target, const MultiRectangle& other) const
{
    const int scaling = other.getNumSubdivisionsPerElement()/m_subdivisions;
    const dim_t NN0 = m_NN[0], NN1 = m_NN[1], otherNN0 = other.getNumNodesPerDim()[0];
    const dim_t numComp = source.getDataPointSize();
    target.requireWrite();
#pragma omp parallel for
    for (dim_t ny = 0; ny < NN1 - 1; ny++) { //source nodes
        for (dim_t nx = 0; nx < NN0 - 1; nx++) {
            const double *x0y0 = source.getSampleDataRO(ny*NN0 + nx);
            const double *x0y1 = source.getSampleDataRO((ny+1)*NN0 + nx);
            const double *x1y0 = source.getSampleDataRO(ny*NN0 + nx + 1);
            const double *x1y1 = source.getSampleDataRO((ny+1)*NN0 + nx + 1);
            const double origin[2] = {getLocalCoordinate(nx, 0), getLocalCoordinate(ny, 1)};
            for (int sy = 0; sy < scaling + 1; sy++) { //target nodes
                const double y = (other.getLocalCoordinate(ny*scaling+sy, 1) - origin[1]) / m_dx[1];
                for (int sx = 0; sx < scaling + 1; sx++) {
                    const double x = (other.getLocalCoordinate(nx*scaling+sx, 0) - origin[0]) / m_dx[0];
                    double *out = target.getSampleDataRW(nx*scaling+sx + (ny*scaling+sy)*otherNN0);
                    for (int comp = 0; comp < numComp; comp++) {
                        out[comp] = x0y0[comp]*(1-x)*(1-y) + x1y0[comp]*x*(1-y) + x0y1[comp]*(1-x)*y + x1y1[comp]*x*y;
                    }
                }
            }
        }
    }
}

void MultiRectangle::interpolateReducedToElementsFiner(const escript::Data& source,
        escript::Data& target, const MultiRectangle& other) const
{
    const int scaling = other.getNumSubdivisionsPerElement()/m_subdivisions;
    const dim_t numComp = source.getDataPointSize();
    target.requireWrite();
    //for each of ours
#pragma omp parallel for
    for (dim_t ey = 0; ey < m_NE[1]; ey++) {
        for (dim_t ex = 0; ex < m_NE[0]; ex++) {
            const double *in = source.getSampleDataRO(ex + ey*m_NE[0]);
            //for each subelement
            for (dim_t sy = 0; sy < scaling; sy++) {
                const dim_t ty = ey*scaling + sy;
                for (dim_t sx = 0; sx < scaling; sx++) {
                    const dim_t tx = ex*scaling + sx;
                    double *out = target.getSampleDataRW(tx + ty*m_NE[0]*scaling);
                    for (dim_t comp = 0; comp < numComp; comp++) {
                        const double quadvalue = in[comp];
                        out[INDEX3(comp, 0, 0, numComp, 2)] = quadvalue;
                        out[INDEX3(comp, 0, 1, numComp, 2)] = quadvalue;
                        out[INDEX3(comp, 1, 0, numComp, 2)] = quadvalue;
                        out[INDEX3(comp, 1, 1, numComp, 2)] = quadvalue;
                    }
                }
            }
        }
    }
}

void MultiRectangle::interpolateReducedToReducedFiner(const escript::Data& source,
        escript::Data& target, const MultiRectangle& other) const
{
    const int scaling = other.getNumSubdivisionsPerElement()/m_subdivisions;
    const dim_t numComp = source.getDataPointSize();
    target.requireWrite();
    //for each of ours
#pragma omp parallel for
    for (dim_t ey = 0; ey < m_NE[1]; ey++) {
        for (dim_t ex = 0; ex < m_NE[0]; ex++) {
            const double *in = source.getSampleDataRO(ex + ey*m_NE[0]);
            //for each subelement
            for (dim_t sy = 0; sy < scaling; sy++) {
                const dim_t ty = ey*scaling + sy;
                for (dim_t sx = 0; sx < scaling; sx++) {
                    const dim_t tx = ex*scaling + sx;
                    double *out = target.getSampleDataRW(tx + ty*m_NE[0]*scaling);
                    for (dim_t comp = 0; comp < numComp; comp++) {
                        out[comp] = in[comp];
                    }
                }
            }
        }
    }
}

void MultiRectangle::interpolateNodesToElementsFiner(const escript::Data& source,
        escript::Data& target, const MultiRectangle& other) const
{
    const int scaling = other.getNumSubdivisionsPerElement()/m_subdivisions;
    const dim_t NE0 = m_NE[0], NE1 = m_NE[1];
    const dim_t numComp = source.getDataPointSize();
    target.requireWrite();
#pragma omp parallel for
    for (dim_t ey = 0; ey < NE1; ey++) { //source nodes
        for (dim_t ex = 0; ex < NE0; ex++) {
            const double *x0y0 = source.getSampleDataRO(ey*(NE0+1) + ex);
            const double *x0y1 = source.getSampleDataRO((ey+1)*(NE0+1) + ex);
            const double *x1y0 = source.getSampleDataRO(ey*(NE0+1) + ex + 1);
            const double *x1y1 = source.getSampleDataRO((ey+1)*(NE0+1) + ex + 1);
            const double origin[2] = {getLocalCoordinate(ex, 0), getLocalCoordinate(ey, 1)};
            for (int sy = 0; sy < scaling; sy++) { //target elements
                for (int sx = 0; sx < scaling; sx++) {
                    const double x1 = (other.getLocalCoordinate(ex*scaling+sx, 0) - origin[0]) / m_dx[0] + FIRST_QUAD/scaling;
                    const double x2 = x1 + (SECOND_QUAD - FIRST_QUAD)/scaling;
                    const double y1 = (other.getLocalCoordinate(ey*scaling+sy, 1) - origin[1]) / m_dx[1] + FIRST_QUAD/scaling;
                    const double y2 = y1 + (SECOND_QUAD - FIRST_QUAD)/scaling;
                    double *out = target.getSampleDataRW(ex*scaling+sx + (ey*scaling+sy)*NE0*scaling);
                    for (int comp = 0; comp < numComp; comp++) {
                        out[INDEX3(comp, 0, 0, numComp, 2)] = x0y0[comp]*(1-x1)*(1-y1) + x1y0[comp]*x1*(1-y1) + x0y1[comp]*(1-x1)*y1 + x1y1[comp]*x1*y1;
                        out[INDEX3(comp, 0, 1, numComp, 2)] = x0y0[comp]*(1-x1)*(1-y2) + x1y0[comp]*x1*(1-y2) + x0y1[comp]*(1-x1)*y2 + x1y1[comp]*x1*y2;
                        out[INDEX3(comp, 1, 0, numComp, 2)] = x0y0[comp]*(1-x2)*(1-y1) + x1y0[comp]*x2*(1-y1) + x0y1[comp]*(1-x2)*y1 + x1y1[comp]*x2*y1;
                        out[INDEX3(comp, 1, 1, numComp, 2)] = x0y0[comp]*(1-x2)*(1-y2) + x1y0[comp]*x2*(1-y2) + x0y1[comp]*(1-x2)*y2 + x1y1[comp]*x2*y2;
                    }
                }
            }
        }
    }
}

void MultiRectangle::interpolateElementsToElementsCoarser(const escript::Data& source,
        escript::Data& target, const MultiRectangle& other) const
{
    const int scaling = m_subdivisions/other.getNumSubdivisionsPerElement();
    const double scaling_volume = (1./scaling)*(1./scaling);
    const dim_t *theirNE = other.getNumElementsPerDim();
    const dim_t numComp = source.getDataPointSize();

    vector<double> points(scaling*2, 0);
    vector<double> first_lagrange(scaling*2, 1);
    vector<double> second_lagrange(scaling*2, 1);
    
    for (int i = 0; i < scaling*2; i+=2) {
        points[i] = (i/2 + FIRST_QUAD)/scaling;
        points[i+1] = (i/2 + SECOND_QUAD)/scaling;
    }
    
    for (int i = 0; i < scaling*2; i++) {
        first_lagrange[i] = (points[i] - SECOND_QUAD) / (FIRST_QUAD - SECOND_QUAD);
        second_lagrange[i] = (points[i] - FIRST_QUAD) / (SECOND_QUAD - FIRST_QUAD);
    }
    target.requireWrite();
    //for each of theirs
#pragma omp parallel for
    for (dim_t ty = 0; ty < theirNE[1]; ty++) {
        for (dim_t tx = 0; tx < theirNE[0]; tx++) {
            double *out = target.getSampleDataRW(tx + ty*theirNE[0]);
            //for each subelement
            for (dim_t sy = 0; sy < scaling; sy++) {
                const dim_t ey = ty*scaling + sy;
                for (dim_t sx = 0; sx < scaling; sx++) {
                    const dim_t ex = tx*scaling + sx;
                    const double *in = source.getSampleDataRO(ex + ey*m_NE[0]);
                    for (int quad = 0; quad < 4; quad++) {
                        int lx = sx*2 + quad%2;
                        int ly = sy*2 + quad/2;
                        for (dim_t comp = 0; comp < numComp; comp++) {
                            const double quadvalue = scaling_volume * in[comp + quad*numComp];
                            out[INDEX3(comp, 0, 0, numComp, 2)] += quadvalue * first_lagrange[lx] * first_lagrange[ly];
                            out[INDEX3(comp, 0, 1, numComp, 2)] += quadvalue * first_lagrange[lx] * second_lagrange[ly];
                            out[INDEX3(comp, 1, 0, numComp, 2)] += quadvalue * second_lagrange[lx] * first_lagrange[ly];
                            out[INDEX3(comp, 1, 1, numComp, 2)] += quadvalue * second_lagrange[lx] * second_lagrange[ly];
                        }
                    }
                }
            }
        }
    }
}


void MultiRectangle::interpolateElementsToElementsFiner(const escript::Data& source,
        escript::Data& target, const MultiRectangle& other) const
{
    const int scaling = other.getNumSubdivisionsPerElement()/m_subdivisions;
    const dim_t numComp = source.getDataPointSize();

    vector<double> points(scaling*2, 0);
    vector<double> lagranges(scaling*4, 1);

    for (int i = 0; i < scaling*2; i+=2) {
        points[i] = (i/2 + FIRST_QUAD)/scaling;
        points[i+1] = (i/2 + SECOND_QUAD)/scaling;
    }
    for (int i = 0; i < scaling*2; i++) {
        lagranges[i] = (points[i] - SECOND_QUAD) / (FIRST_QUAD - SECOND_QUAD);
        lagranges[i + 2*scaling] = (points[i] - FIRST_QUAD) / (SECOND_QUAD - FIRST_QUAD);
    }
    target.requireWrite();
    //for each of ours
#pragma omp parallel for
    for (dim_t ey = 0; ey < m_NE[1]; ey++) {
        for (dim_t ex = 0; ex < m_NE[0]; ex++) {
            const double *in = source.getSampleDataRO(ex + ey*m_NE[0]);
            //for each subelement
            for (dim_t sy = 0; sy < scaling; sy++) {
                const dim_t ty = ey*scaling + sy;
                for (dim_t sx = 0; sx < scaling; sx++) {
                    const dim_t tx = ex*scaling + sx;
                    double *out = target.getSampleDataRW(tx + ty*m_NE[0]*scaling);
                    for (int quad = 0; quad < 4; quad++) {
                        const int lx = scaling*2*(quad%2) + sx*2;
                        const int ly = scaling*2*(quad/2) + sy*2;
                        for (dim_t comp = 0; comp < numComp; comp++) {
                            const double quadvalue = in[comp + quad*numComp];
                            out[INDEX3(comp, 0, 0, numComp, 2)] += quadvalue * lagranges[lx] * lagranges[ly];
                            out[INDEX3(comp, 0, 1, numComp, 2)] += quadvalue * lagranges[lx] * lagranges[ly+1];
                            out[INDEX3(comp, 1, 0, numComp, 2)] += quadvalue * lagranges[lx+1] * lagranges[ly];
                            out[INDEX3(comp, 1, 1, numComp, 2)] += quadvalue * lagranges[lx+1] * lagranges[ly+1];
                        }
                    }
                }
            }
        }
    }
}

void MultiRectangle::interpolateAcross(escript::Data& target,
                                     const escript::Data& source) const
{
    const MultiRectangle *other =
                dynamic_cast<const MultiRectangle *>(target.getDomain().get());
    if (other == NULL)
        throw RipleyException("Invalid interpolation: Domains must both be instances of MultiRectangle");
    //shouldn't ever happen, but I want to know if it does
    if (other == this)
        throw RipleyException("interpolateAcross: this domain is the target");
        
    validateInterpolationAcross(source.getFunctionSpace().getTypeCode(),
            *(target.getDomain().get()), target.getFunctionSpace().getTypeCode());
    int fsSource = source.getFunctionSpace().getTypeCode();
    int fsTarget = target.getFunctionSpace().getTypeCode();

    std::stringstream msg;
    msg << "Invalid interpolation: interpolation not implemented for function space "
        << functionSpaceTypeAsString(fsSource)
        << " -> "
        << functionSpaceTypeAsString(fsTarget);
    if (other->getNumSubdivisionsPerElement() > getNumSubdivisionsPerElement()) {
        switch (fsSource) {
            case Nodes:
                switch (fsTarget) {
                    case Nodes:
                    case ReducedNodes:
                    case DegreesOfFreedom:
                    case ReducedDegreesOfFreedom:
                        interpolateNodesToNodesFiner(source, target, *other);
                        return;
                    case Elements:
                        interpolateNodesToElementsFiner(source, target, *other);
                        return;
                }
                break;
            case Elements:
                switch (fsTarget) {
                    case Elements:
                        interpolateElementsToElementsFiner(source, target, *other);
                        return;
                }
                break;
            case ReducedElements:
                switch (fsTarget) {
                    case Elements:
                        interpolateReducedToElementsFiner(source, target, *other);
                        return;
                }
                break;
        }
        msg << " when target is a finer mesh";
    } else {
        switch (fsSource) {
            case Nodes:
                switch (fsTarget) {
                    case Elements:
                        escript::Data elements=escript::Vector(0., escript::function(*this), true);
                        interpolateNodesOnElements(elements, source, false);
                        interpolateElementsToElementsCoarser(elements, target, *other);
                        return;
                }
                break;
            case Elements:
                switch (fsTarget) {
                    case Elements:
                        interpolateElementsToElementsCoarser(source, target, *other);
                        return;
                }
                break;
        }
        msg << " when target is a coarser mesh";
    }
    throw RipleyException(msg.str());
}

string MultiRectangle::getDescription() const
{
    return "ripley::MultiRectangle";
}

bool MultiRectangle::operator==(const AbstractDomain& other) const
{
    const MultiRectangle* o=dynamic_cast<const MultiRectangle*>(&other);
    if (o) {
        return (RipleyDomain::operator==(other) &&
                m_gNE[0]==o->m_gNE[0] && m_gNE[1]==o->m_gNE[1]
                && m_origin[0]==o->m_origin[0] && m_origin[1]==o->m_origin[1]
                && m_length[0]==o->m_length[0] && m_length[1]==o->m_length[1]
                && m_NX[0]==o->m_NX[0] && m_NX[1]==o->m_NX[1]
                && m_subdivisions==o->m_subdivisions);
    }

    return false;
}

void MultiRectangle::readNcGrid(escript::Data& out, string filename, string varname,
            const ReaderParameters& params) const
{
    if (m_subdivisions != 1)
        throw RipleyException("Non-parent MultiRectangles cannot read datafiles");
    Rectangle::readNcGrid(out, filename, varname, params);
}

void MultiRectangle::readBinaryGrid(escript::Data& out, string filename,
                               const ReaderParameters& params) const
{
    if (m_subdivisions != 1)
        throw RipleyException("Non-parent MultiRectangles cannot read datafiles");
    Rectangle::readBinaryGrid(out, filename, params);
}

#ifdef USE_BOOSTIO
void MultiRectangle::readBinaryGridFromZipped(escript::Data& out, string filename,
                               const ReaderParameters& params) const
{
    if (m_subdivisions != 1)
        throw RipleyException("Non-parent MultiRectangles cannot read datafiles");
    Rectangle::readBinaryGridFromZipped(out, filename, params);
}
#endif

void MultiRectangle::writeBinaryGrid(const escript::Data& in, string filename,
                                int byteOrder, int dataType) const
{
    if (m_subdivisions != 1)
        throw RipleyException("Non-parent MultiRectangles cannot read datafiles");
    Rectangle::writeBinaryGrid(in, filename, byteOrder, dataType);
}

void MultiRectangle::dump(const string& fileName) const
{
    if (m_subdivisions != 1)
        throw RipleyException("Non-parent MultiRectangles dump not implemented");
    Rectangle::dump(fileName);
}

void MultiRectangle::populateDofMap()
{
    const index_t left = getFirstInDim(0);
    const index_t bottom = getFirstInDim(1);
    const dim_t nDOF0 = getNumDOFInAxis(0);
    const dim_t nDOF1 = getNumDOFInAxis(1);
    // populate node->DOF mapping with own degrees of freedom.
    // The rest is assigned in the loop further down
    m_dofMap.assign(getNumNodes(), -7);
#pragma omp parallel for
    for (index_t i=bottom; i<bottom+nDOF1; i++) {
        for (index_t j=left; j<left+nDOF0; j++) {
            m_dofMap[i*m_NN[0]+j]=(i-bottom)*nDOF0+j-left;
        }
    }

    // build list of shared components and neighbours by looping through
    // all potential neighbouring ranks and checking if positions are
    // within bounds
    const dim_t numDOF=nDOF0*nDOF1;
    m_colIndices.clear();
    m_rowIndices.clear();
    m_colIndices.resize(numDOF);
    m_rowIndices.resize(getNumNodes() - numDOF);

    RankVector neighbour;
    IndexVector offsetInSharedSend(1,0);
    IndexVector offsetInSharedRecv(1,0);
    IndexVector sendShared, recvShared;
    const int x=m_mpiInfo->rank%m_NX[0];
    const int y=m_mpiInfo->rank/m_NX[0];
    // numShared will contain the number of shared DOFs after the following
    // blocks
    dim_t numShared=0;
    // sharing bottom edge
    if (y > 0) {
        neighbour.push_back((y-1)*m_NX[0] + x);
        //joining edge, send and recv
        offsetInSharedSend.push_back(offsetInSharedSend.back()+nDOF0);
        offsetInSharedRecv.push_back(offsetInSharedRecv.back()+nDOF0*m_subdivisions);
        for (dim_t i=0; i < nDOF0; i++, numShared++) {
            sendShared.push_back(i);
            recvShared.push_back(numDOF+numShared);
            m_dofMap[i+left]=numDOF + numShared;
            const dim_t ind = i;
            if (i > 0)
                doublyLink(m_colIndices, m_rowIndices, ind - 1, numShared);
            doublyLink(m_colIndices, m_rowIndices, ind, numShared);
            if (i < nDOF0 - 1)
                doublyLink(m_colIndices, m_rowIndices, ind + 1, numShared);
        }
    
        for (unsigned sy = 1; sy < m_subdivisions; sy++) {
            for (dim_t i=0; i < nDOF0; i++, numShared++) {
                recvShared.push_back(numDOF+numShared);
                m_dofMap[left + i + sy*m_NN[0]] = numDOF + numShared;
                const dim_t ind = i;
                if (i > 0)
                    doublyLink(m_colIndices, m_rowIndices, ind - 1, numShared);
                doublyLink(m_colIndices, m_rowIndices, ind, numShared);
                if (i < nDOF0 - 1)
                    doublyLink(m_colIndices, m_rowIndices, ind + 1, numShared);
            }
        }
    }
    // sharing top edge
    if (y < m_NX[1] - 1) {
        neighbour.push_back((y+1)*m_NX[0] + x);
        offsetInSharedSend.push_back(offsetInSharedSend.back()+nDOF0*m_subdivisions);
        offsetInSharedRecv.push_back(offsetInSharedRecv.back()+nDOF0);
        // add to send only
        for (unsigned sy = 0; sy < m_subdivisions; sy++) {
            for (dim_t i=0; i < nDOF0; i++) {
                sendShared.push_back(numDOF-nDOF0*(m_subdivisions - sy) + i);
            }
        }
        for (dim_t i=0; i < nDOF0; i++, numShared++) {
            recvShared.push_back(numDOF+numShared);
            m_dofMap[m_NN[0]*(m_NN[1]-1)+left+i]=numDOF+numShared;
            const dim_t ind = numDOF-nDOF0+i;
            if (i > 0)
                doublyLink(m_colIndices, m_rowIndices, ind - 1, numShared);
            doublyLink(m_colIndices, m_rowIndices, ind, numShared);
            if (i < nDOF0 - 1)
                doublyLink(m_colIndices, m_rowIndices, ind + 1, numShared);
        }
    }
    // sharing left edge
    if (x > 0) {
        neighbour.push_back(y*m_NX[0] + x-1);
        offsetInSharedSend.push_back(offsetInSharedSend.back()+nDOF1);
        offsetInSharedRecv.push_back(offsetInSharedRecv.back()+nDOF1*m_subdivisions);
        for (dim_t i=0; i < nDOF1; i++, numShared++) {
            for (unsigned sx = 0; sx < m_subdivisions - 1; sx++, numShared++) {
                recvShared.push_back(numDOF+numShared);
                m_dofMap[(bottom+i)*m_NN[0] + sx] = numDOF + numShared;
                const dim_t ind = i*nDOF0;
                if (i > 0)
                    doublyLink(m_colIndices, m_rowIndices, ind - nDOF0, numShared);
                doublyLink(m_colIndices, m_rowIndices, ind, numShared);
                if (i < nDOF1 - 1)
                    doublyLink(m_colIndices, m_rowIndices, ind + nDOF0, numShared);
            }
            sendShared.push_back(i*nDOF0);
            recvShared.push_back(numDOF + numShared);
            m_dofMap[(bottom+i)*m_NN[0] + m_subdivisions - 1]=numDOF + numShared;
            const dim_t ind = i*nDOF0;
            if (i > 0)
                doublyLink(m_colIndices, m_rowIndices, ind - nDOF0, numShared);
            doublyLink(m_colIndices, m_rowIndices, ind, numShared);
            if (i < nDOF1 - 1)
                doublyLink(m_colIndices, m_rowIndices, ind + nDOF0, numShared);
        }
    }
    // sharing right edge
    if (x < m_NX[0] - 1) {
        neighbour.push_back(y*m_NX[0] + x+1);
        offsetInSharedSend.push_back(offsetInSharedSend.back()+nDOF1*m_subdivisions);
        offsetInSharedRecv.push_back(offsetInSharedRecv.back()+nDOF1);
        for (dim_t i=0; i < nDOF1; i++, numShared++) {
            for (unsigned sx = 0; sx < m_subdivisions - 1; sx++) {
                sendShared.push_back((i+1)*nDOF0-(m_subdivisions - sx));
            }
            sendShared.push_back((i+1)*nDOF0-1);
            recvShared.push_back(numDOF+numShared);
            m_dofMap[(bottom+1+i)*m_NN[0]- 1]=numDOF+numShared;
            const dim_t ind = (i+1)*nDOF0 - 1;
            if (i > 0)
                doublyLink(m_colIndices, m_rowIndices, ind - nDOF0, numShared);
            doublyLink(m_colIndices, m_rowIndices, ind, numShared);
            if (i < nDOF1 - 1)
                doublyLink(m_colIndices, m_rowIndices, ind + nDOF0, numShared);
        }
    }
    // sharing bottom-left node
    if (x > 0 && y > 0) {
        neighbour.push_back((y-1)*m_NX[0] + x-1);
        // sharing a node
        offsetInSharedSend.push_back(offsetInSharedSend.back()+1);
        offsetInSharedRecv.push_back(offsetInSharedRecv.back()+m_subdivisions*m_subdivisions);
        for (unsigned sy = 0; sy < m_subdivisions; sy++) {
            for (unsigned sx = 0; sx < m_subdivisions; sx++, numShared++) {
                m_dofMap[sx + sy*m_NN[0]] = numDOF + numShared;
                recvShared.push_back(numDOF+numShared);
                doublyLink(m_colIndices, m_rowIndices, 0, numShared);
            }
        }
        sendShared.push_back(0);
    }
    // sharing top-left node
    if (x > 0 && y < m_NX[1]-1) {
        neighbour.push_back((y+1)*m_NX[0] + x-1);
        offsetInSharedSend.push_back(offsetInSharedSend.back()+m_subdivisions);
        offsetInSharedRecv.push_back(offsetInSharedRecv.back()+m_subdivisions);
        for (int s = 0; s < m_subdivisions; s++, numShared++) {
            sendShared.push_back(numDOF - (m_subdivisions - s)*nDOF0);
            recvShared.push_back(numDOF + numShared);
            m_dofMap[m_NN[0]*(m_NN[1]-1) + s] = numDOF + numShared;
            if (s > 0)
                doublyLink(m_colIndices, m_rowIndices, numDOF - (m_subdivisions - s + 1)*nDOF0, numShared);
            doublyLink(m_colIndices, m_rowIndices, numDOF - (m_subdivisions - s)*nDOF0, numShared);
            if (s < m_subdivisions - 1)
                doublyLink(m_colIndices, m_rowIndices, numDOF - (m_subdivisions - s - 1)*nDOF0, numShared);            
        }
    }
    // sharing bottom-right node
    if (x < m_NX[0]-1 && y > 0) {
        neighbour.push_back((y-1)*m_NX[0] + x+1);
        offsetInSharedSend.push_back(offsetInSharedSend.back()+m_subdivisions);
        offsetInSharedRecv.push_back(offsetInSharedRecv.back()+m_subdivisions);
        for (int s = 0; s < m_subdivisions; s++, numShared++) {
            recvShared.push_back(numDOF+numShared);
            m_dofMap[(s+1)*m_NN[0] - 1] = numDOF + numShared;
            sendShared.push_back(nDOF0-(m_subdivisions-s));
            const dim_t ind = nDOF0 - (m_subdivisions - s);
            if (s > 0)
                doublyLink(m_colIndices, m_rowIndices, ind - 1, numShared);
            doublyLink(m_colIndices, m_rowIndices, ind, numShared);
            if (s < m_subdivisions - 1)
                doublyLink(m_colIndices, m_rowIndices, ind + 1, numShared);
        }
    }
    // sharing top-right node
    if (x < m_NX[0]-1 && y < m_NX[1]-1) {
        neighbour.push_back((y+1)*m_NX[0] + x+1);
        offsetInSharedSend.push_back(offsetInSharedSend.back()+m_subdivisions*m_subdivisions);
        offsetInSharedRecv.push_back(offsetInSharedRecv.back()+1);
        for (unsigned sy = 0; sy < m_subdivisions; sy++) {
            for (unsigned sx = 0; sx < m_subdivisions; sx++) {
                sendShared.push_back(numDOF-(m_subdivisions - sy - 1)*nDOF0 - (m_subdivisions - sx));
            }
        }
        recvShared.push_back(numDOF+numShared);
        m_dofMap[m_NN[0]*m_NN[1]-1]=numDOF+numShared;
        doublyLink(m_colIndices, m_rowIndices, numDOF-1, numShared);
        ++numShared;
    }

    // TODO: paso::SharedComponents should take vectors to avoid this
    int* neighPtr = NULL;
    index_t* sendPtr = NULL;
    index_t* recvPtr = NULL;
    if (neighbour.size() > 0) {
        neighPtr = &neighbour[0];
        sendPtr = &sendShared[0];
        recvPtr = &recvShared[0];
    }
    // create connector
    paso::SharedComponents_ptr snd_shcomp(new paso::SharedComponents(
            numDOF, neighbour.size(), neighPtr, sendPtr,
            &offsetInSharedSend[0], 1, 0, m_mpiInfo));
    paso::SharedComponents_ptr rcv_shcomp(new paso::SharedComponents(
            numDOF, neighbour.size(), neighPtr, recvPtr,
            &offsetInSharedRecv[0], 1, 0, m_mpiInfo));
    m_connector.reset(new paso::Connector(snd_shcomp, rcv_shcomp));
}

void MultiRectangle::populateSampleIds()
{
    // degrees of freedom are numbered from left to right, bottom to top in
    // each rank, continuing on the next rank (ranks also go left-right,
    // bottom-top).
    // This means rank 0 has id 0...n0-1, rank 1 has id n0...n1-1 etc. which
    // helps when writing out data rank after rank.

    // build node distribution vector first.
    // rank i owns m_nodeDistribution[i+1]-nodeDistribution[i] nodes.
    // Unlike regular ripley domains this is NOT constant for all ranks so
    // we do an Allgather (we could have also computed per rank but it's a bit
    // involved)
    m_nodeDistribution.assign(m_mpiInfo->size+1, 0);
    dim_t numDOF=getNumDOF();
    if (m_mpiInfo->size > 1) {
#if ESYS_MPI
        MPI_Allgather(&numDOF, 1, MPI_DIM_T, &m_nodeDistribution[0], 1,
                      MPI_DIM_T, m_mpiInfo->comm);

        // accumulate
        dim_t accu = 0;
        for (int rank=0; rank<m_mpiInfo->size; rank++) {
            const dim_t n = m_nodeDistribution[rank];
            m_nodeDistribution[rank] = accu;
            accu += n;
        }
        ESYS_ASSERT_MPI(accu == getNumDataPointsGlobal(),
                "something went wrong computing the DOF distribution!",
                m_mpiInfo);

        m_nodeDistribution[m_mpiInfo->size] = accu;
#endif
    } else {
        m_nodeDistribution[m_mpiInfo->size] = numDOF;
    }

    try {
        m_nodeId.resize(getNumNodes());
        m_dofId.resize(numDOF);
        m_elementId.resize(getNumElements());
    } catch (const std::length_error& le) {
        throw RipleyException("The system does not have sufficient memory for a domain of this size.");
    }

    // populate face element counts
    //left
    if (m_offset[0]==0)
        m_faceCount[0]=m_NE[1];
    else
        m_faceCount[0]=0;
    //right
    if (m_mpiInfo->rank%m_NX[0]==m_NX[0]-1)
        m_faceCount[1]=m_NE[1];
    else
        m_faceCount[1]=0;
    //bottom
    if (m_offset[1]==0)
        m_faceCount[2]=m_NE[0];
    else
        m_faceCount[2]=0;
    //top
    if (m_mpiInfo->rank/m_NX[0]==m_NX[1]-1)
        m_faceCount[3]=m_NE[0];
    else
        m_faceCount[3]=0;

    const dim_t NFE = getNumFaceElements();
    m_faceId.resize(NFE);

    const index_t left = getFirstInDim(0);
    const index_t bottom = getFirstInDim(1);
    const dim_t nDOF0 = getNumDOFInAxis(0);
    const dim_t nDOF1 = getNumDOFInAxis(1);
    const dim_t NE0 = m_NE[0];
    const dim_t NE1 = m_NE[1];

#define globalNodeId(x,y) \
    ((m_offset[0]+x)/nDOF0)*nDOF0*nDOF1+(m_offset[0]+x)%nDOF0 \
    + ((m_offset[1]+y)/nDOF1)*nDOF0*nDOF1*m_NX[0]+((m_offset[1]+y)%nDOF1)*nDOF0

    // set corner id's outside the parallel region
    m_nodeId[0] = globalNodeId(0, 0);
    m_nodeId[m_NN[0]-1] = globalNodeId(m_NN[0]-1, 0);
    m_nodeId[m_NN[0]*(m_NN[1]-1)] = globalNodeId(0, m_NN[1]-1);
    m_nodeId[m_NN[0]*m_NN[1]-1] = globalNodeId(m_NN[0]-1,m_NN[1]-1);
#undef globalNodeId

#pragma omp parallel
    {
        // populate degrees of freedom and own nodes (identical id)
#pragma omp for nowait
        for (dim_t i=0; i<nDOF1; i++) {
            for (dim_t j=0; j<nDOF0; j++) {
                const index_t nodeIdx=j+left+(i+bottom)*m_NN[0];
                const index_t dofIdx=j+i*nDOF0;
                m_dofId[dofIdx] = m_nodeId[nodeIdx]
                    = m_nodeDistribution[m_mpiInfo->rank]+dofIdx;
            }
        }

        // populate the rest of the nodes (shared with other ranks)
        if (m_faceCount[0]==0) { // left column
#pragma omp for nowait
            for (dim_t i=0; i<nDOF1; i++) {
                const index_t nodeIdx=(i+bottom)*m_NN[0];
                const index_t dofId=(i+1)*nDOF0-1;
                m_nodeId[nodeIdx]
                    = m_nodeDistribution[m_mpiInfo->rank-1]+dofId;
            }
        }
        if (m_faceCount[1]==0) { // right column
#pragma omp for nowait
            for (dim_t i=0; i<nDOF1; i++) {
                const index_t nodeIdx=(i+bottom+1)*m_NN[0]-1;
                const index_t dofId=i*nDOF0;
                m_nodeId[nodeIdx]
                    = m_nodeDistribution[m_mpiInfo->rank+1]+dofId;
            }
        }
        if (m_faceCount[2]==0) { // bottom row
#pragma omp for nowait
            for (dim_t i=0; i<nDOF0; i++) {
                const index_t nodeIdx=i+left;
                const index_t dofId=nDOF0*(nDOF1-1)+i;
                m_nodeId[nodeIdx]
                    = m_nodeDistribution[m_mpiInfo->rank-m_NX[0]]+dofId;
            }
        }
        if (m_faceCount[3]==0) { // top row
#pragma omp for nowait
            for (dim_t i=0; i<nDOF0; i++) {
                const index_t nodeIdx=m_NN[0]*(m_NN[1]-1)+i+left;
                const index_t dofId=i;
                m_nodeId[nodeIdx]
                    = m_nodeDistribution[m_mpiInfo->rank+m_NX[0]]+dofId;
            }
        }

        // populate element id's
#pragma omp for nowait
        for (dim_t i1=0; i1<NE1; i1++) {
            for (dim_t i0=0; i0<NE0; i0++) {
                m_elementId[i0+i1*NE0]=(m_offset[1]+i1)*m_gNE[0]+m_offset[0]+i0;
            }
        }

        // face elements
#pragma omp for
        for (dim_t k=0; k<NFE; k++)
            m_faceId[k]=k;
    } // end parallel section

    m_nodeTags.assign(getNumNodes(), 0);
    updateTagsInUse(Nodes);

    m_elementTags.assign(getNumElements(), 0);
    updateTagsInUse(Elements);

    // generate face offset vector and set face tags
    const index_t LEFT=1, RIGHT=2, BOTTOM=10, TOP=20;
    const index_t faceTag[] = { LEFT, RIGHT, BOTTOM, TOP };
    m_faceOffset.assign(4, -1);
    m_faceTags.clear();
    index_t offset=0;
    for (size_t i=0; i<4; i++) {
        if (m_faceCount[i]>0) {
            m_faceOffset[i]=offset;
            offset+=m_faceCount[i];
            m_faceTags.insert(m_faceTags.end(), m_faceCount[i], faceTag[i]);
        }
    }
    setTagMap("left", LEFT);
    setTagMap("right", RIGHT);
    setTagMap("bottom", BOTTOM);
    setTagMap("top", TOP);
    updateTagsInUse(FaceElements);

    populateDofMap();
}

paso::SystemMatrixPattern_ptr MultiRectangle::getPasoMatrixPattern(
                                                    bool reducedRowOrder,
                                                    bool reducedColOrder) const
{
    if (m_pattern.get())
        return m_pattern;

    // first call - create pattern, then return
    const dim_t numDOF = getNumDOF();
    const dim_t numShared = getNumNodes() - numDOF;
#pragma omp parallel for
    for (dim_t i = 0; i < numShared; i++) {
        sort(m_rowIndices[i].begin(), m_rowIndices[i].end());
    }

    // create main and couple blocks
    paso::Pattern_ptr mainPattern = createPasoPattern(getConnections(), numDOF);
    paso::Pattern_ptr colPattern = createPasoPattern(m_colIndices, numShared);
    paso::Pattern_ptr rowPattern = createPasoPattern(m_rowIndices, numDOF);

    // allocate paso distribution
    paso::Distribution_ptr distribution(new paso::Distribution(m_mpiInfo,
            const_cast<index_t*>(&m_nodeDistribution[0]), 1, 0));

    // finally create the system matrix pattern
    m_pattern.reset(new paso::SystemMatrixPattern(MATRIX_FORMAT_DEFAULT,
            distribution, distribution, mainPattern, colPattern, rowPattern,
            m_connector, m_connector));
    return m_pattern;
}

RankVector MultiRectangle::getOwnerVector(int fsType) const
{
    if (m_subdivisions != 1)
        throw RipleyException("Multiresolution domains only support ownership for the coarsest level");
    return Rectangle::getOwnerVector(fsType);
}

dim_t MultiRectangle::findNode(const double *coords) const
{
    const dim_t NOT_MINE = -1;
    //is the found element even owned by this rank
    // (inside owned or shared elements but will map to an owned element)
    for (int dim = 0; dim < m_numDim; dim++) {
        double min = m_origin[dim] + m_offset[dim]* m_dx[dim]
                - m_dx[dim]/2.; //allows for point outside mapping onto node
        double max = m_origin[dim] + (m_offset[dim] + m_NE[dim])*m_dx[dim]
                + m_dx[dim]/2.;
        if (min > coords[dim] || max < coords[dim]) {
            return NOT_MINE;
        }
    }
    // get distance from origin
    double x = coords[0] - m_origin[0];
    double y = coords[1] - m_origin[1];

    //check if the point is even inside the domain
    if (x < 0 || y < 0 || x > m_length[0] || y > m_length[1])
        return NOT_MINE;

    // distance in elements
    dim_t ex = (dim_t) floor((x + 0.01*m_dx[0]) / m_dx[0]);
    dim_t ey = (dim_t) floor((y + 0.01*m_dx[1]) / m_dx[1]);
    // set the min distance high enough to be outside the element plus a bit
    dim_t closest = NOT_MINE;
    double minDist = 1;
    for (int dim = 0; dim < m_numDim; dim++) {
        minDist += m_dx[dim]*m_dx[dim];
    }
    //find the closest node
    for (int dx = 0; dx < 1; dx++) {
        double xdist = (x - (ex + dx)*m_dx[0]);
        for (int dy = 0; dy < 1; dy++) {
            double ydist = (y - (ey + dy)*m_dx[1]);
            double total = xdist*xdist + ydist*ydist;
            if (total < minDist) {
                closest = INDEX2(ex+dx-m_offset[0], ey+dy-m_offset[1], m_NN[0]);
                minDist = total;
            }
        }
    }
    //if this happens, we've let a dirac point slip through, which is awful
    if (closest == NOT_MINE) {
        throw RipleyException("Unable to map appropriate dirac point to a node,"
                " implementation problem in MultiRectangle::findNode()");
    }
    return closest;
}

} // end of namespace ripley

