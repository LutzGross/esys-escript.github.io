
/*****************************************************************************
*
* Copyright (c) 2003-2020 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/

#include <ripley/MultiRectangle.h>
#include <ripley/blocktools.h>
#include <ripley/domainhelpers.h>

#include <escript/DataFactory.h>
#include <escript/FunctionSpaceFactory.h>
#include <escript/index.h>

#define FIRST_QUAD 0.21132486540518711775
#define SECOND_QUAD 0.78867513459481288225

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <limits>

using std::vector;
using std::string;

namespace ripley {

MultiRectangle::MultiRectangle(dim_t n0, dim_t n1, double x0, double y0,
                     double x1, double y1, int d0, int d1,
                     const vector<double>& points,
                     const vector<int>& tags,
                     const TagMap& tagnamestonums,
                     unsigned int subdivisions)
     : Rectangle(n0,n1, x0,y0, x1,y1, d0,d1, points, tags, tagnamestonums),
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

MultiRectangle::MultiRectangle(escript::JMPI jmpi, dim_t n0, dim_t n1, double x0, double y0,
                     double x1, double y1, int d0, int d1,
                     const vector<double>& points,
                     const vector<int>& tags,
                     const TagMap& tagnamestonums,
                     unsigned int subdivisions)
     : Rectangle(jmpi, n0,n1, x0,y0, x1,y1, d0,d1, points, tags, tagnamestonums),
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
    if (source.isComplex()!=target.isComplex())
    {
        throw RipleyException("Programmer Error: in and out parameters do not have the same complexity.");        
    }
    if (source.isComplex())
    {
        interpolateNodesToNodesFinerWorker(source, target, other, escript::DataTypes::cplx_t(0));
    }
    else
    {
        interpolateNodesToNodesFinerWorker(source, target, other, escript::DataTypes::real_t(0));      
    }
}

template <typename S>
void MultiRectangle::interpolateNodesToNodesFinerWorker(const escript::Data& source,
        escript::Data& target, const MultiRectangle& other, S sentinel) const
{
    const int scaling = other.getNumSubdivisionsPerElement()/m_subdivisions;
    const dim_t NN0 = m_NN[0], NN1 = m_NN[1], otherNN0 = other.getNumNodesPerDim()[0];
    const dim_t numComp = source.getDataPointSize();
    target.requireWrite();
#pragma omp parallel for
    for (dim_t ny = 0; ny < NN1 - 1; ny++) { //source nodes
        for (dim_t nx = 0; nx < NN0 - 1; nx++) {
            const S *x0y0 = source.getSampleDataRO(ny*NN0 + nx, sentinel);
            const S *x0y1 = source.getSampleDataRO((ny+1)*NN0 + nx, sentinel);
            const S *x1y0 = source.getSampleDataRO(ny*NN0 + nx + 1, sentinel);
            const S *x1y1 = source.getSampleDataRO((ny+1)*NN0 + nx + 1, sentinel);
            const S origin[2] = {getLocalCoordinate(nx, 0), getLocalCoordinate(ny, 1)};
            for (int sy = 0; sy < scaling + 1; sy++) { //target nodes
                const S y = (other.getLocalCoordinate(ny*scaling+sy, 1) - origin[1]) / m_dx[1];
                for (int sx = 0; sx < scaling + 1; sx++) {
                    const S x = (other.getLocalCoordinate(nx*scaling+sx, 0) - origin[0]) / m_dx[0];
                    S *out = target.getSampleDataRW(nx*scaling+sx + (ny*scaling+sy)*otherNN0, sentinel);
                    for (int comp = 0; comp < numComp; comp++) {
                        out[comp] = x0y0[comp]*(static_cast<S>(1)-x)*(static_cast<S>(1)-y) 
			   + x1y0[comp]*x*(static_cast<S>(1)-y) + x0y1[comp]*(static_cast<S>(1)-x)*y + x1y1[comp]*x*y;
                    }
                }
            }
        }
    }
}

void MultiRectangle::interpolateReducedToElementsFiner(const escript::Data& source,
        escript::Data& target, const MultiRectangle& other) const
{
    if (source.isComplex()!=target.isComplex())
    {
        throw RipleyException("Programmer Error: in and out parameters do not have the same complexity.");        
    }
    if (source.isComplex())
    {
        interpolateReducedToElementsFinerWorker(source, target, other, escript::DataTypes::cplx_t(0));
    }
    else
    {
        interpolateReducedToElementsFinerWorker(source, target, other, escript::DataTypes::real_t(0));      
    }
}

template <typename S>
void MultiRectangle::interpolateReducedToElementsFinerWorker(const escript::Data& source,
        escript::Data& target, const MultiRectangle& other, S sentinel) const
{
    const int scaling = other.getNumSubdivisionsPerElement()/m_subdivisions;
    const dim_t numComp = source.getDataPointSize();
    target.requireWrite();
    //for each of ours
#pragma omp parallel for
    for (dim_t ey = 0; ey < m_NE[1]; ey++) {
        for (dim_t ex = 0; ex < m_NE[0]; ex++) {
            const S *in = source.getSampleDataRO(ex + ey*m_NE[0], sentinel);
            //for each subelement
            for (dim_t sy = 0; sy < scaling; sy++) {
                const dim_t ty = ey*scaling + sy;
                for (dim_t sx = 0; sx < scaling; sx++) {
                    const dim_t tx = ex*scaling + sx;
                    S *out = target.getSampleDataRW(tx + ty*m_NE[0]*scaling, sentinel);
                    for (dim_t comp = 0; comp < numComp; comp++) {
                        const S quadvalue = in[comp];
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
    if (source.isComplex()!=target.isComplex())
    {
        throw RipleyException("Programmer Error: in and out parameters do not have the same complexity.");        
    }
    if (source.isComplex())
    {
        interpolateReducedToReducedFinerWorker(source, target, other, escript::DataTypes::cplx_t(0));
    }
    else
    {
        interpolateReducedToReducedFinerWorker(source, target, other, escript::DataTypes::real_t(0));      
    }
}

template <typename S>
void MultiRectangle::interpolateReducedToReducedFinerWorker(const escript::Data& source,
        escript::Data& target, const MultiRectangle& other, S sentinel) const
{
    const int scaling = other.getNumSubdivisionsPerElement()/m_subdivisions;
    const dim_t numComp = source.getDataPointSize();
    target.requireWrite();
    //for each of ours
#pragma omp parallel for
    for (dim_t ey = 0; ey < m_NE[1]; ey++) {
        for (dim_t ex = 0; ex < m_NE[0]; ex++) {
            const S *in = source.getSampleDataRO(ex + ey*m_NE[0], sentinel);
            //for each subelement
            for (dim_t sy = 0; sy < scaling; sy++) {
                const dim_t ty = ey*scaling + sy;
                for (dim_t sx = 0; sx < scaling; sx++) {
                    const dim_t tx = ex*scaling + sx;
                    S *out = target.getSampleDataRW(tx + ty*m_NE[0]*scaling, sentinel);
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
    if (source.isComplex()!=target.isComplex())
    {
        throw RipleyException("Programmer Error: in and out parameters do not have the same complexity.");        
    }
    if (source.isComplex())
    {
        interpolateNodesToElementsFinerWorker(source, target, other, escript::DataTypes::cplx_t(0));
    }
    else
    {
        interpolateNodesToElementsFinerWorker(source, target, other, escript::DataTypes::real_t(0));      
    }
}

template <typename S>
void MultiRectangle::interpolateNodesToElementsFinerWorker(const escript::Data& source,
        escript::Data& target, const MultiRectangle& other, S sentinel) const
{
    const int scaling = other.getNumSubdivisionsPerElement()/m_subdivisions;
    const dim_t NE0 = m_NE[0], NE1 = m_NE[1];
    const dim_t numComp = source.getDataPointSize();
    target.requireWrite();
#pragma omp parallel for
    for (dim_t ey = 0; ey < NE1; ey++) { //source nodes
        for (dim_t ex = 0; ex < NE0; ex++) {
            const S *x0y0 = source.getSampleDataRO(ey*(NE0+1) + ex, sentinel);
            const S *x0y1 = source.getSampleDataRO((ey+1)*(NE0+1) + ex, sentinel);
            const S *x1y0 = source.getSampleDataRO(ey*(NE0+1) + ex + 1, sentinel);
            const S *x1y1 = source.getSampleDataRO((ey+1)*(NE0+1) + ex + 1, sentinel);
            const S origin[2] = {getLocalCoordinate(ex, 0), getLocalCoordinate(ey, 1)};
            for (int sy = 0; sy < scaling; sy++) { //target elements
                for (int sx = 0; sx < scaling; sx++) {
                    const S x1 = (other.getLocalCoordinate(ex*scaling+sx, 0) - origin[0]) / m_dx[0] + FIRST_QUAD/scaling;
                    const S x2 = x1 + (SECOND_QUAD - FIRST_QUAD)/scaling;
                    const S y1 = (other.getLocalCoordinate(ey*scaling+sy, 1) - origin[1]) / m_dx[1] + FIRST_QUAD/scaling;
                    const S y2 = y1 + (SECOND_QUAD - FIRST_QUAD)/scaling;
		    const S mx1=static_cast<S>(1)-x1;
		    const S mx2=static_cast<S>(1)-x2;
		    const S my1=static_cast<S>(1)-y1;
		    const S my2=static_cast<S>(1)-y2;
                    S *out = target.getSampleDataRW(ex*scaling+sx + (ey*scaling+sy)*NE0*scaling, sentinel);
                    for (int comp = 0; comp < numComp; comp++) {
                        out[INDEX3(comp, 0, 0, numComp, 2)] = x0y0[comp]*(mx1)*(my1) + x1y0[comp]*x1*(my1) + x0y1[comp]*(mx1)*y1 + x1y1[comp]*x1*y1;
                        out[INDEX3(comp, 0, 1, numComp, 2)] = x0y0[comp]*(mx1)*(my2) + x1y0[comp]*x1*(my2) + x0y1[comp]*(mx1)*y2 + x1y1[comp]*x1*y2;
                        out[INDEX3(comp, 1, 0, numComp, 2)] = x0y0[comp]*(mx2)*(my1) + x1y0[comp]*x2*(my1) + x0y1[comp]*(mx2)*y1 + x1y1[comp]*x2*y1;
                        out[INDEX3(comp, 1, 1, numComp, 2)] = x0y0[comp]*(mx2)*(my2) + x1y0[comp]*x2*(my2) + x0y1[comp]*(mx2)*y2 + x1y1[comp]*x2*y2;
                    }
                }
            }
        }
    }
}

void MultiRectangle::interpolateElementsToElementsCoarser(const escript::Data& source,
        escript::Data& target, const MultiRectangle& other) const
{
    if (source.isComplex()!=target.isComplex())
    {
        throw RipleyException("Programmer Error: in and out parameters do not have the same complexity.");        
    }
    if (source.isComplex())
    {
        interpolateElementsToElementsCoarserWorker(source, target, other, escript::DataTypes::cplx_t(0));
    }
    else
    {
        interpolateElementsToElementsCoarserWorker(source, target, other, escript::DataTypes::real_t(0));      
    }
}

template <typename S>
void MultiRectangle::interpolateElementsToElementsCoarserWorker(const escript::Data& source,
        escript::Data& target, const MultiRectangle& other, S sentinel) const
{
    const int scaling = m_subdivisions/other.getNumSubdivisionsPerElement();
    const S scaling_volume = (1./scaling)*(1./scaling);
    const dim_t *theirNE = other.getNumElementsPerDim();
    const dim_t numComp = source.getDataPointSize();

    vector<S> points(scaling*2, 0);
    vector<S> first_lagrange(scaling*2, 1);
    vector<S> second_lagrange(scaling*2, 1);
    
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
            S *out = target.getSampleDataRW(tx + ty*theirNE[0], sentinel);
            //for each subelement
            for (dim_t sy = 0; sy < scaling; sy++) {
                const dim_t ey = ty*scaling + sy;
                for (dim_t sx = 0; sx < scaling; sx++) {
                    const dim_t ex = tx*scaling + sx;
                    const S *in = source.getSampleDataRO(ex + ey*m_NE[0], sentinel);
                    for (int quad = 0; quad < 4; quad++) {
                        int lx = sx*2 + quad%2;
                        int ly = sy*2 + quad/2;
                        for (dim_t comp = 0; comp < numComp; comp++) {
                            const S quadvalue = scaling_volume * in[comp + quad*numComp];
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
    if (source.isComplex()!=target.isComplex())
    {
        throw RipleyException("Programmer Error: in and out parameters do not have the same complexity.");        
    }
    if (source.isComplex())
    {
        interpolateElementsToElementsFinerWorker(source, target, other, escript::DataTypes::cplx_t(0));
    }
    else
    {
        interpolateElementsToElementsFinerWorker(source, target, other, escript::DataTypes::real_t(0));      
    }
}

template <typename S>
void MultiRectangle::interpolateElementsToElementsFinerWorker(const escript::Data& source,
        escript::Data& target, const MultiRectangle& other, S sentinel) const
{
    const int scaling = other.getNumSubdivisionsPerElement()/m_subdivisions;
    const dim_t numComp = source.getDataPointSize();

    vector<S> points(scaling*2, 0);
    vector<S> lagranges(scaling*4, 1);

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
            const S *in = source.getSampleDataRO(ex + ey*m_NE[0], sentinel);
            //for each subelement
            for (dim_t sy = 0; sy < scaling; sy++) {
                const dim_t ty = ey*scaling + sy;
                for (dim_t sx = 0; sx < scaling; sx++) {
                    const dim_t tx = ex*scaling + sx;
                    S *out = target.getSampleDataRW(tx + ty*m_NE[0]*scaling, sentinel);
                    for (int quad = 0; quad < 4; quad++) {
                        const int lx = scaling*2*(quad%2) + sx*2;
                        const int ly = scaling*2*(quad/2) + sy*2;
                        for (dim_t comp = 0; comp < numComp; comp++) {
                            const S quadvalue = in[comp + quad*numComp];
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

void MultiRectangle::readBinaryGridFromZipped(escript::Data& out, string filename,
                               const ReaderParameters& params) const
{
    if (m_subdivisions != 1)
        throw RipleyException("Non-parent MultiRectangles cannot read datafiles");
    Rectangle::readBinaryGridFromZipped(out, filename, params);
}

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
    // build node distribution vector first.
    // rank i owns m_nodeDistribution[i+1]-nodeDistribution[i] nodes.
    // Unlike regular ripley domains this is NOT constant for all ranks so
    // we do an Allgather (we could have also computed per rank but it's a bit
    // involved)
    m_nodeDistribution.assign(m_mpiInfo->size+1, 0);
    dim_t numDOF = getNumDOF();
    if (m_mpiInfo->size > 1) {
#if ESYS_MPI
        MPI_Allgather(&numDOF, 1, MPI_DIM_T, &m_nodeDistribution[0], 1,
                      MPI_DIM_T, m_mpiInfo->comm);

        // accumulate
        dim_t accu = 0;
        for (int rank = 0; rank < m_mpiInfo->size; rank++) {
            const dim_t n = m_nodeDistribution[rank];
            m_nodeDistribution[rank] = accu;
            accu += n;
        }
        ESYS_ASSERT(accu == getNumDataPointsGlobal(),
                "something went wrong computing the DOF distribution!");

        m_nodeDistribution[m_mpiInfo->size] = accu;
#endif
    } else {
        m_nodeDistribution[m_mpiInfo->size] = numDOF;
    }

    // degrees of freedom are numbered from left to right, bottom to top in
    // each rank, continuing on the next rank (ranks also go left-right,
    // bottom-top).
    // This means rank 0 has id 0...n0-1, rank 1 has id n0...n1-1 etc. which
    // helps when writing out data rank after rank.

    try {
        m_nodeId.assign(getNumNodes(), -1);
        m_dofMap.assign(getNumNodes(), -1);
        m_dofId.assign(numDOF, -1);
    } catch (const std::length_error& le) {
        throw RipleyException("The system does not have sufficient memory for a domain of this size.");
    }

    const index_t left = getFirstInDim(0);
    const index_t bottom = getFirstInDim(1);
    const dim_t nDOF0 = getNumDOFInAxis(0);
    const dim_t nDOF1 = getNumDOFInAxis(1);
    // populate node->DOF mapping, DOF IDs and own node IDs.
    // The rest of the node IDs are communicated further down.
#pragma omp parallel for
    for (dim_t i=0; i<nDOF1; i++) {
        for (dim_t j=0; j<nDOF0; j++) {
            const index_t nodeIdx = j+left + (i+bottom)*m_NN[0];
            const index_t dofIdx = j + i*nDOF0;
            m_dofMap[nodeIdx] = dofIdx;
            m_dofId[dofIdx] = m_nodeId[nodeIdx]
                = m_nodeDistribution[m_mpiInfo->rank] + dofIdx;
        }
    }

    // build list of shared components and neighbours
    m_colIndices.assign(numDOF, IndexVector());
    m_rowIndices.assign(getNumNodes() - numDOF, IndexVector());

    RankVector neighbour;
    IndexVector offsetInSharedSend(1,0);
    IndexVector offsetInSharedRecv(1,0);
    IndexVector sendShared, recvShared;
    const int x = m_mpiInfo->rank%m_NX[0];
    const int y = m_mpiInfo->rank/m_NX[0];
    // numShared will contain the number of shared DOFs after the following
    // blocks
    dim_t numShared = 0;
    // sharing bottom edge
    if (y > 0) {
        neighbour.push_back((y-1)*m_NX[0] + x);
        //joining edge, send and recv
        offsetInSharedSend.push_back(offsetInSharedSend.back()+nDOF0);
        offsetInSharedRecv.push_back(offsetInSharedRecv.back()+nDOF0*m_subdivisions);
        // add to send only
        for (dim_t i=0; i < nDOF0; i++) {
            sendShared.push_back(i);
        }
    
        for (unsigned sy = 0; sy < m_subdivisions; sy++) {
            for (index_t i = 0; i < nDOF0; i++, numShared++) {
                const index_t nodeIdx = left + i + sy*m_NN[0];
                const index_t dofIdx = i;
                recvShared.push_back(nodeIdx);
                m_dofMap[nodeIdx] = numDOF + numShared;
                if (i > 0)
                    doublyLink(m_colIndices, m_rowIndices, dofIdx - 1, numShared);
                doublyLink(m_colIndices, m_rowIndices, dofIdx, numShared);
                if (i < nDOF0 - 1)
                    doublyLink(m_colIndices, m_rowIndices, dofIdx + 1, numShared);
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
            for (index_t i = 0; i < nDOF0; i++) {
                sendShared.push_back(numDOF-nDOF0*(m_subdivisions - sy) + i);
            }
        }
        for (index_t i = 0; i < nDOF0; i++, numShared++) {
            const index_t nodeIdx = left + i + m_NN[0]*(m_NN[1]-1);
            const index_t dofIdx = numDOF - nDOF0 + i;
            recvShared.push_back(nodeIdx);
            m_dofMap[nodeIdx] = numDOF+numShared;
            if (i > 0)
                doublyLink(m_colIndices, m_rowIndices, dofIdx - 1, numShared);
            doublyLink(m_colIndices, m_rowIndices, dofIdx, numShared);
            if (i < nDOF0 - 1)
                doublyLink(m_colIndices, m_rowIndices, dofIdx + 1, numShared);
        }
    }
    // sharing left edge
    if (x > 0) {
        neighbour.push_back(y*m_NX[0] + x-1);
        offsetInSharedSend.push_back(offsetInSharedSend.back()+nDOF1);
        offsetInSharedRecv.push_back(offsetInSharedRecv.back()+nDOF1*m_subdivisions);
        for (index_t i = 0; i < nDOF1; i++) {
            const index_t dofIdx = i*nDOF0;
            sendShared.push_back(dofIdx);
            for (unsigned sx = 0; sx < m_subdivisions; sx++, numShared++) {
                const index_t nodeIdx = (bottom+i)*m_NN[0] + sx;
                recvShared.push_back(nodeIdx);
                m_dofMap[nodeIdx] = numDOF + numShared;
                if (i > 0)
                    doublyLink(m_colIndices, m_rowIndices, dofIdx - nDOF0, numShared);
                doublyLink(m_colIndices, m_rowIndices, dofIdx, numShared);
                if (i < nDOF1 - 1)
                    doublyLink(m_colIndices, m_rowIndices, dofIdx + nDOF0, numShared);
            }
        }
    }
    // sharing right edge
    if (x < m_NX[0] - 1) {
        neighbour.push_back(y*m_NX[0] + x+1);
        offsetInSharedSend.push_back(offsetInSharedSend.back()+nDOF1*m_subdivisions);
        offsetInSharedRecv.push_back(offsetInSharedRecv.back()+nDOF1);
        for (index_t i = 0; i < nDOF1; i++, numShared++) {
            for (unsigned sx = 0; sx < m_subdivisions - 1; sx++) {
                sendShared.push_back((i+1)*nDOF0-(m_subdivisions - sx));
            }
            const index_t nodeIdx = (bottom+1+i)*m_NN[0] - 1;
            const index_t dofIdx = (i+1)*nDOF0 - 1;
            sendShared.push_back(dofIdx);
            recvShared.push_back(nodeIdx);
            m_dofMap[nodeIdx] = numDOF + numShared;
            if (i > 0)
                doublyLink(m_colIndices, m_rowIndices, dofIdx - nDOF0, numShared);
            doublyLink(m_colIndices, m_rowIndices, dofIdx, numShared);
            if (i < nDOF1 - 1)
                doublyLink(m_colIndices, m_rowIndices, dofIdx + nDOF0, numShared);
        }
    }
    // sharing bottom-left block
    if (x > 0 && y > 0) {
        neighbour.push_back((y-1)*m_NX[0] + x-1);
        // sharing a node
        offsetInSharedSend.push_back(offsetInSharedSend.back()+1);
        offsetInSharedRecv.push_back(offsetInSharedRecv.back()+m_subdivisions*m_subdivisions);
        for (unsigned sy = 0; sy < m_subdivisions; sy++) {
            for (unsigned sx = 0; sx < m_subdivisions; sx++, numShared++) {
                const index_t nodeIdx = sx + sy*m_NN[0];
                m_dofMap[nodeIdx] = numDOF + numShared;
                recvShared.push_back(nodeIdx);
                doublyLink(m_colIndices, m_rowIndices, 0, numShared);
            }
        }
        sendShared.push_back(0);
    }
    // sharing top-left block
    if (x > 0 && y < m_NX[1]-1) {
        neighbour.push_back((y+1)*m_NX[0] + x-1);
        offsetInSharedSend.push_back(offsetInSharedSend.back()+m_subdivisions);
        offsetInSharedRecv.push_back(offsetInSharedRecv.back()+m_subdivisions);
        for (int s = 0; s < m_subdivisions; s++, numShared++) {
            const index_t nodeIdx = m_NN[0]*(m_NN[1]-1) + s;
            const index_t dofIdx = numDOF - (m_subdivisions - s)*nDOF0;
            sendShared.push_back(dofIdx);
            recvShared.push_back(nodeIdx);
            m_dofMap[nodeIdx] = numDOF + numShared;
            if (s > 0)
                doublyLink(m_colIndices, m_rowIndices, dofIdx - nDOF0, numShared);
            doublyLink(m_colIndices, m_rowIndices, dofIdx, numShared);
            if (s < m_subdivisions - 1)
                doublyLink(m_colIndices, m_rowIndices, dofIdx + nDOF0, numShared);            
        }
    }
    // sharing bottom-right block
    if (x < m_NX[0]-1 && y > 0) {
        neighbour.push_back((y-1)*m_NX[0] + x+1);
        offsetInSharedSend.push_back(offsetInSharedSend.back()+m_subdivisions);
        offsetInSharedRecv.push_back(offsetInSharedRecv.back()+m_subdivisions);
        for (int s = 0; s < m_subdivisions; s++, numShared++) {
            const index_t nodeIdx = (s+1)*m_NN[0] - 1;
            const dim_t ind = nDOF0 - (m_subdivisions - s);
            recvShared.push_back(nodeIdx);
            m_dofMap[nodeIdx] = numDOF + numShared;
            sendShared.push_back(ind);
            if (s > 0)
                doublyLink(m_colIndices, m_rowIndices, ind - 1, numShared);
            doublyLink(m_colIndices, m_rowIndices, ind, numShared);
            if (s < m_subdivisions - 1)
                doublyLink(m_colIndices, m_rowIndices, ind + 1, numShared);
        }
    }
    // sharing top-right block
    if (x < m_NX[0]-1 && y < m_NX[1]-1) {
        neighbour.push_back((y+1)*m_NX[0] + x+1);
        offsetInSharedSend.push_back(offsetInSharedSend.back()+m_subdivisions*m_subdivisions);
        offsetInSharedRecv.push_back(offsetInSharedRecv.back()+1);
        for (unsigned sy = 0; sy < m_subdivisions; sy++) {
            for (unsigned sx = 0; sx < m_subdivisions; sx++) {
                sendShared.push_back(numDOF-(m_subdivisions - sy - 1)*nDOF0 - (m_subdivisions - sx));
            }
        }
        const dim_t nodeIdx = m_NN[0]*m_NN[1] - 1;
        recvShared.push_back(nodeIdx);
        m_dofMap[nodeIdx] = numDOF+numShared;
        doublyLink(m_colIndices, m_rowIndices, numDOF-1, numShared);
        ++numShared;
    }

#ifdef ESYS_MPI
    if (m_mpiInfo->size > 1) {
        // now send off shared DOF IDs so nodeId will become a global node
        // labelling.
        const dim_t numSend = offsetInSharedSend.back();
        const dim_t numRecv = offsetInSharedRecv.back();
        IndexVector recvBuffer(numRecv);
        IndexVector sendBuffer(numSend);
        std::vector<MPI_Request> reqs(2*neighbour.size());
        std::vector<MPI_Status> stats(2*neighbour.size());

        // prepare the send buffer
#pragma omp parallel for
        for (index_t i = 0; i < numSend; ++i) {
            sendBuffer[i] = m_dofId[sendShared[i]];
        }
        for (index_t i = 0; i < neighbour.size(); i++) {
            MPI_Irecv(&recvBuffer[offsetInSharedRecv[i]],
                    offsetInSharedRecv[i+1] - offsetInSharedRecv[i],
                    MPI_DIM_T, neighbour[i],
                    m_mpiInfo->counter()+neighbour[i],
                    m_mpiInfo->comm, &reqs[2*i]);
            MPI_Issend(&sendBuffer[offsetInSharedSend[i]],
                    offsetInSharedSend[i+1] - offsetInSharedSend[i],
                    MPI_DIM_T, neighbour[i],
                    m_mpiInfo->counter()+m_mpiInfo->rank, m_mpiInfo->comm,
                    &reqs[2*i+1]);
        }
        m_mpiInfo->incCounter(m_mpiInfo->size);

        // do something else here...

        MPI_Waitall(2*neighbour.size(), &reqs[0], &stats[0]);

        // now populate rest of node IDs
#pragma omp parallel for
        for (index_t i=0; i < numRecv; i++) {
            const index_t nodeIdx = recvShared[i];
            m_nodeId[nodeIdx] = recvBuffer[m_dofMap[nodeIdx]-numDOF];
        }
    }
#endif // ESYS_MPI

#ifdef ESYS_HAVE_PASO
    createPasoConnector(neighbour, offsetInSharedSend, offsetInSharedRecv,
                        sendShared, recvShared);
#endif
}

void MultiRectangle::populateSampleIds()
{
    // label nodes and DOF first
    populateDofMap();

    m_elementId.assign(getNumElements(), -1);

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
    const dim_t NE0 = m_NE[0];
    const dim_t NE1 = m_NE[1];
    m_faceId.resize(NFE);

#pragma omp parallel
    {
        // populate element IDs
#pragma omp for nowait
        for (index_t i1 = 0; i1 < NE1; i1++) {
            for (index_t i0 = 0; i0 < NE0; i0++) {
                m_elementId[i0+i1*NE0]=(m_offset[1]+i1)*m_gNE[0]+m_offset[0]+i0;
            }
        }

        // face elements
#pragma omp for
        for (index_t k = 0; k < NFE; k++)
            m_faceId[k] = k;
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
    index_t offset = 0;
    for (size_t i=0; i<4; i++) {
        if (m_faceCount[i] > 0) {
            m_faceOffset[i]=offset;
            offset += m_faceCount[i];
            m_faceTags.insert(m_faceTags.end(), m_faceCount[i], faceTag[i]);
        }
    }
    setTagMap("left", LEFT);
    setTagMap("right", RIGHT);
    setTagMap("bottom", BOTTOM);
    setTagMap("top", TOP);
    updateTagsInUse(FaceElements);

}

#ifdef ESYS_HAVE_PASO
paso::SystemMatrixPattern_ptr MultiRectangle::getPasoMatrixPattern(
                                                    bool reducedRowOrder,
                                                    bool reducedColOrder) const
{
    if (!m_pattern) {
        // first call - create pattern, then return
        const dim_t numDOF = getNumDOF();
        const dim_t numShared = getNumNodes() - numDOF;
#pragma omp parallel for
        for (index_t i = 0; i < numShared; i++) {
            sort(m_rowIndices[i].begin(), m_rowIndices[i].end());
        }

        // create main and couple blocks
        paso::Pattern_ptr mainPattern = createPasoPattern(getConnections(), numDOF);
        paso::Pattern_ptr colPattern = createPasoPattern(m_colIndices, numShared);
        paso::Pattern_ptr rowPattern = createPasoPattern(m_rowIndices, numDOF);

        // allocate Paso distribution
        escript::Distribution_ptr distribution(new escript::Distribution(
                                               m_mpiInfo, m_nodeDistribution));

        // finally create the system matrix pattern
        m_pattern.reset(new paso::SystemMatrixPattern(MATRIX_FORMAT_DEFAULT,
                distribution, distribution, mainPattern, colPattern, rowPattern,
                getPasoConnector(), getPasoConnector()));
    }
    return m_pattern;
}
#endif

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

